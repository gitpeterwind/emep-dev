module Gravset_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Testing for ash gravitational settling
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
use Chemfields_ml,        only: xn_adv
use ChemChemicals_ml,     only: species,species_adv
use ChemSpecs_adv_ml
use ChemSpecs_shl_ml,     only: NSPEC_SHL
use ChemGroups_ml,        only: chemgroups
use DryDep_ml,            only: DDepMap
use GridValues_ml,        only: A_mid,B_mid,A_bnd,B_bnd
use MetFields_ml,         only: roa,th,ps
use ModelConstants_ml,    only: KMAX_MID,KMAX_BND,dt_advec
!use OwnDataTypes_ml,      only: depmap
use Par_ml,               only: MAXLIMAX,MAXLJMAX,limax,ljmax,me,li0,li1,lj0,lj1
use PhysicalConstants_ml, only: GRAV
use SmallUtils_ml,        only: find_index
use NetCDF_ml,            only: printCDF

implicit none
private

public :: gravset

real, public, allocatable, dimension(:,:,:,:) ::num_sed

contains

subroutine gravset()
  real                        :: slinnfac = 1.0,&
                                 density = 2.5E03 !3.0E03
  real                        :: ztempx,vt,zsedl,tempc,knut
  real                        :: Re,Re_f, X ,vt_old,z,ztemp
  real, save                  :: wil_hua
  real,dimension(7)           :: B_n = &
       (/-3.18657,0.992696,-0.00153193,-0.000987059,-0.000578878 ,&
       0.0000855176,-0.00000327815/)
  real,dimension(KMAX_MID)    :: zvis,p_mid,zlair, zflux,zdp1
  real,dimension(KMAX_BND)    :: p_full
  integer                     :: i,j,k,ispec,ash,volc_group
  integer,save                :: bins,volc
  integer                     :: v,b
  character(len=256)          :: name     ! Volcano (vent) name
  logical                     :: first_call = .true.,test_log = .false.
  real                        :: F = 0.8


! need to calculate
! zvis -> dynamic viscosity of air.. dependent on temperature

  type :: sediment
    integer :: spec
    real    :: diameter
  end type sediment
  type(sediment),save,allocatable,dimension(:) :: grav_sed


  if (first_call) then
    first_call = .false.

    write(*,*) "Gravset called!"

    wil_hua = F**(-0.828) + 2*SQRT(1.07-F)
    ash=find_index("ASH",chemgroups(:)%name)
    if(ash>0)then
       bins=size(chemgroups(ash)%specs)
       allocate(num_sed(bins,MAXLIMAX,MAXLJMAX,KMAX_MID))
       num_sed(:,:,:,:)=0.0
       allocate(grav_sed(bins))

       select case(bins)
       case(7)
         grav_sed(:)%spec = chemgroups(volc_group)%specs(:)-NSPEC_SHL
         grav_sed(:)%diameter = [0.1,0.3,1.0,3.0,10.0,30.0,100.0]*1e-6
       case(9)
         grav_sed(:)%spec = chemgroups(volc_group)%specs(:)-NSPEC_SHL
         grav_sed(:)%diameter = [4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,25.0]*1e-6
!      case(10)
!        grav_sed(:)%spec = chemgroups(volc_group)%specs(:)-NSPEC_SHL
!        grav_sed(:)%diameter = [2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,25.0]*1e-6
       case default
         call CheckStop("Unsupported number of ASH bins")
       end select
     end if
  end if !first_call

  do j = lj0,lj1
     do i = li0,li1
        do k = 1,KMAX_MID

           ! dynamic viscosity of air after Prup.Klett in [Pa s]
           tempc = th(i,j,k,1) - 273.15
           if (tempc >= 0.0 ) then
              zvis(k) = (1.718 + 0.0049*tempc)*1.E-5
           else
              zvis(k) = (1.718 + 0.0049*tempc - 1.2E-05*(tempc**2))*1.E-5
           endif

           ! mean free path of air (Prupp. Klett) in [10^-6 m]
           p_mid(k) = A_mid(k)+B_mid(k)*ps(i,j,1)
           zlair(k) = 0.066 *(1.01325E+5/p_mid(k))*(th(i,j,k,1)/293.15)*1.E-06

           ! air mass auxiliary  variable --> zdp1 [kg/(m^2 *s)]
           p_full(k) = A_bnd(k)+B_bnd(k)*ps(i,j,1)
        end do
        p_full(KMAX_BND) = A_bnd(KMAX_BND)+B_bnd(KMAX_BND)*ps(i,j,1)

        do k = 1,KMAX_MID
          zdp1(k)=(p_full(k+1) - p_full(k))/(GRAV*dt_advec) ! do outside of k-loop????
        end do

        do b = 1,bins
          do k = 1,KMAX_MID-1
            knut = 2*zlair(k)/grav_sed(b)%diameter
            ztemp = 2.*((grav_sed(b)%diameter/2)**2)*(density-roa(i,j,k,1))*GRAV/ &! roa [kg m-3]
                    (9.*zvis(k))![m/s]
          ! with Cunningham slip-flow correction
            vt = ztemp*slinnfac*     &
                (1.+ 1.257*knut+0.4*knut*EXP(-1.1/(knut))) ![m/s]
            Re = grav_sed(b)%diameter*vt/(zvis(k)/roa(i,j,k,1))
            vt_old = vt
            vt = vt/wil_hua

            num_sed(b,i,j,k)= vt
          ! calculation of sedimentation flux zflux[kg/(m^2 s)]=zsedl*zdp1
          ! definition of  zflux=vt*ztm1(:,:,jt)*zdens
          ! compute flux in terms of mixing ratio zsedl= zflux/zdp1 -->>zsedl [kg/kg]
          ! change of tracer tendency according to loss of tracer
          ! due to sedimentation from the box
          ! unit of zdp1 kg of air m-2 s-1

            ztempx = min(1.0,vt*roa(i,j,k,1)/zdp1(k)) ! 1, loss is limited to content of box
            zsedl = ztempx*xn_adv(grav_sed(b)%spec,i,j,k)        ! kg kg-1, loss in terms of mixing ratio ! blir likt uavhengig av vt
            zflux(k) = zsedl*zdp1(k)                     ! --> [kg m-2 s-1]

          ! loss of mass in layer
            xn_adv(grav_sed(b)%spec,i,j,k) = xn_adv(grav_sed(b)%spec,i,j,k) - zsedl !  kg kg-1
          end do

          ! teste å gjøre det i en ny k-loop så det ikke blir så effektivt

          do  k = 1,KMAX_MID-1
            ! "arrival" of sedimented mass in box below
            if (k <KMAX_MID) then
              xn_adv(grav_sed(b)%spec,i,j,k+1) =  xn_adv(grav_sed(b)%spec,i,j,k+1) +  zflux(k)/zdp1(k+1)
            end if
          end do
        end do
      end do
    end do
  end do

! Multilayer crossing is here no ralised!!!
! sedimentation velocity is in effect limited to z/delt
! sedimentation to the ground from first layer sflx --> [kg m-2 s-1]
! sflux = zflux(:,1) ! sedimenterer ikke på bakken, bare sender det til nederste laget
end subroutine gravset


end module Gravset_ml
