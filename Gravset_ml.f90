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

subroutine gravset



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
endtype sediment
type(sediment),save,allocatable,dimension(:,:) :: grav_sed

 


if (first_call) then
   first_call = .false.

   write(*,*) "Gravset called!"

   wil_hua = F**(-0.828) + 2*SQRT(1.07-F)
 !  write(*,*) "Wilson Huang: ",wil_hua
   volc = 0
   ash=find_index("ASH",chemgroups(:)%name)
   if(ash>0)then
      name="none"
      do i=1,size(chemgroups(ash)%ptr)
         if(species(chemgroups(ash)%ptr(i))%name(1:9)==name)cycle
         name=species(chemgroups(ash)%ptr(i))%name(1:9)
         !write(*,*) "name",name
         volc_group=find_index(name,chemgroups(:)%name)
        ! write(*,*) "volc_group ",volc_group
         
         if (volc_group .gt. 0) then
            volc = volc+1
       !     write(*,*) "inn her? ",volc
            bins = 9
         end if
      enddo
   
      allocate(num_sed(bins,MAXLIMAX,MAXLJMAX,KMAX_MID))
      num_sed(:,:,:,:)=0.0
      !write(*,*) "volcanoes: ",volc
      name="none"
      allocate(grav_sed(volc,bins))
      v = 0

      if (bins == 7 ) then 
         do i=1,size(chemgroups(ash)%specs)
            if(species(chemgroups(ash)%specs(i))%name(1:9)==name)cycle
            name=species(chemgroups(ash)%specs(i))%name(1:9)
            volc_group=find_index(name,chemgroups(:)%name)
            v = v+1
         !   write(*,*) "inn her? ",volc,v,i,trim(name),volc_group
            grav_sed(v,1) = sediment(chemgroups(volc_group)%specs(1)-NSPEC_SHL,1.00E-07)
            grav_sed(v,2) = sediment(chemgroups(volc_group)%specs(2)-NSPEC_SHL,3.00E-07)
            grav_sed(v,3) = sediment(chemgroups(volc_group)%specs(3)-NSPEC_SHL,1.00E-06)
            grav_sed(v,4) = sediment(chemgroups(volc_group)%specs(4)-NSPEC_SHL,3.00E-06)
            grav_sed(v,5) = sediment(chemgroups(volc_group)%specs(5)-NSPEC_SHL,10.00E-06)
            grav_sed(v,6) = sediment(chemgroups(volc_group)%specs(6)-NSPEC_SHL,30.00E-06)
            grav_sed(v,7) = sediment(chemgroups(volc_group)%specs(7)-NSPEC_SHL,100.00E-06)
        !    write(*,*) "specie: ",grav_sed(volc,:)%spec
        !    write(*,*) "diameter",grav_sed(volc,:)%diameter
         end do
      end if

      if (bins == 9 ) then 
         do i=1,size(chemgroups(ash)%ptr)
            if(species(chemgroups(ash)%ptr(i))%name(1:9)==name)cycle
            name=species(chemgroups(ash)%ptr(i))%name(1:9)
            volc_group=find_index(name,chemgroups(:)%name)
            grav_sed(volc,1) = sediment(chemgroups(volc_group)%ptr(1)-NSPEC_SHL,2.00E-06)
            grav_sed(volc,1) = sediment(chemgroups(volc_group)%ptr(1)-NSPEC_SHL,4.00E-06)
            grav_sed(volc,2) = sediment(chemgroups(volc_group)%ptr(2)-NSPEC_SHL,6.00E-06)
            grav_sed(volc,3) = sediment(chemgroups(volc_group)%ptr(3)-NSPEC_SHL,8.00E-06)
            grav_sed(volc,4) = sediment(chemgroups(volc_group)%ptr(4)-NSPEC_SHL,10.00E-06)
            grav_sed(volc,5) = sediment(chemgroups(volc_group)%ptr(5)-NSPEC_SHL,12.00E-06)
            grav_sed(volc,6) = sediment(chemgroups(volc_group)%ptr(6)-NSPEC_SHL,14.00E-06)
            grav_sed(volc,7) = sediment(chemgroups(volc_group)%ptr(7)-NSPEC_SHL,16.00E-06)
            grav_sed(volc,8) = sediment(chemgroups(volc_group)%ptr(8)-NSPEC_SHL,18.00E-06)
            grav_sed(volc,9) = sediment(chemgroups(volc_group)%ptr(9)-NSPEC_SHL,25.00E-06)

          !  grav_sed(volc,7) = sediment(chemgroups(volc_group)%ptr(7)-NSPEC_SHL,20.00E-06)
          !  grav_sed(volc,8) = sediment(chemgroups(volc_group)%ptr(8)-NSPEC_SHL,100.00E-06)
            ! grav_sed(volc,9) = sediment(chemgroups(volc_group)%ptr(9)-NSPEC_SHL,1000.00E-06)
         end do
      end if
       if (bins == 10 ) then 
         do i=1,size(chemgroups(ash)%ptr)
            if(species(chemgroups(ash)%ptr(i))%name(1:9)==name)cycle
            name=species(chemgroups(ash)%ptr(i))%name(1:9)
            volc_group=find_index(name,chemgroups(:)%name)
           ! do i=1,size(chemgroups(ash)%ptr)
            write(*,*) "chemgroups(volc_group) ",chemgroups(volc_group)%ptr(1),NSPEC_SHL,KMAX_MID
            grav_sed(volc,1) = sediment(chemgroups(volc_group)%ptr(1)-NSPEC_SHL,2.00E-06)
            grav_sed(volc,2) = sediment(chemgroups(volc_group)%ptr(1)-NSPEC_SHL,4.00E-06)
            grav_sed(volc,3) = sediment(chemgroups(volc_group)%ptr(2)-NSPEC_SHL,6.00E-06)
            grav_sed(volc,4) = sediment(chemgroups(volc_group)%ptr(3)-NSPEC_SHL,8.00E-06)
            grav_sed(volc,5) = sediment(chemgroups(volc_group)%ptr(4)-NSPEC_SHL,10.00E-06)
            grav_sed(volc,6) = sediment(chemgroups(volc_group)%ptr(5)-NSPEC_SHL,12.00E-06)
            grav_sed(volc,7) = sediment(chemgroups(volc_group)%ptr(6)-NSPEC_SHL,14.00E-06)
            grav_sed(volc,8) = sediment(chemgroups(volc_group)%ptr(7)-NSPEC_SHL,16.00E-06)
            grav_sed(volc,9) = sediment(chemgroups(volc_group)%ptr(8)-NSPEC_SHL,18.00E-06)
            grav_sed(volc,10) = sediment(chemgroups(volc_group)%ptr(9)-NSPEC_SHL,25.00E-06)
         end do
      end if
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
        
      do  v =1,volc ! loop trough all the ash sedimentation particles
         do b = 1,bins
            do k = 1,KMAX_MID-1
        
               knut = 2*zlair(k)/grav_sed(v,b)%diameter
               ztemp = 2.*((grav_sed(v,b)%diameter/2)**2)*(density-roa(i,j,k,1))*GRAV/ &! roa [kg m-3]
                    (9.*zvis(k))![m/s]
       ! with Cunningham slip-flow correction
               vt = ztemp*slinnfac*     &
                    (1.+ 1.257*knut+0.4*knut*EXP(-1.1/(knut))) ![m/s]
               Re = grav_sed(v,b)%diameter*vt/(zvis(k)/roa(i,j,k,1))
               vt_old = vt

           !    if (Re > 0.01) then ! Outside Stokes flow
           !       X = log(32*((grav_sed(v,b)%diameter/2)**3)*(density-roa(i,j,k,1))*roa(i,j,k,1)*GRAV/(3*(zvis(k))**2))
     
     !             Re_f = (1.+ 1.257*knut+0.4*knut*exp(-1.1/(knut)))*exp(B_n(1) + B_n(2)*X+B_n(3)*X**2 + B_n(4)*X**3 +&
     !                                                                       B_n(5)*X**4 + B_n(6)*X**5 + B_n(7)*X**6)
     !             vt = Re_f*(zvis(k)/roa(i,j,k,1))/ grav_sed(v,b)%diameter
  ! 
  !             end if
           

            vt = vt/wil_hua

               num_sed(b,i,j,k)= vt
          !  if (j .eq. 5 .and. i .eq.5 .and. b .eq.8 ) then 
          !     write(*,*) "Wilson Huang: ",wil_hua
          !     write(*,*) "num_sed ", vt,num_sed(b,i,j,k),b,k,me
          !  end if


         !   write(*,*) "vertical velocity: ",vt,b
           
               !  calculation of sedimentation flux zflux[kg/(m^2 s)]=zsedl*zdp1
               !  definition of  zflux=vt*ztm1(:,:,jt)*zdens 
               !  compute flux in terms of mixing ratio zsedl= zflux/zdp1 -->>zsedl [kg/kg]
               !  change of tracer tendency according to loss of tracer
               !  due to sedimentation from the box
               !  unit of zdp1 kg of air m-2 s-1
    
               ztempx = min(1.0,vt*roa(i,j,k,1)/zdp1(k)) ! 1, loss is limited to content of box
          
         
               zsedl = ztempx*xn_adv(grav_sed(v,b)%spec,i,j,k)        ! kg kg-1, loss in terms of mixing ratio ! blir likt uavhengig av vt
        
               zflux(k) = zsedl*zdp1(k)                     ! --> [kg m-2 s-1]
 
               
               !loss of mass in layer
               xn_adv(grav_sed(v,b)%spec,i,j,k) = xn_adv(grav_sed(v,b)%spec,i,j,k) - zsedl !  kg kg-1
         
            end do

            ! teste å gjøre det i en ny k-loop så det ikke blir så effektivt
               
            do  k = 1,KMAX_MID-1
               ! "arrival" of sedimented mass in box below
               if (k <KMAX_MID) then 
                  xn_adv(grav_sed(v,b)%spec,i,j,k+1) =  xn_adv(grav_sed(v,b)%spec,i,j,k+1) +  zflux(k)/zdp1(k+1)
               endif
            end do
         end do
      end do
   end do
end do

!end if
        ! Multilayer crossing is here no ralised!!!
           ! sedimentation velocity is in effect limited to z/delt

           ! sedimentation to the ground from first layer sflx --> [kg m-2 s-1]
           !sflux = zflux(:,1) ! sedimenterer ikke på bakken, bare sender det til nederste laget

          

end subroutine gravset

!end do
end module Gravset_ml
