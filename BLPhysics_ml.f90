module BLPhysics_ml
 ! Collection of boundary layer met routines. Will move more code from tiphys
 !  here in future. Try to keep 1-D or elemental to allow use in offline codes 
 ! (*No* routines in use, except for testing)
 ! Copiedfrom Unimod.rv18WTA  14th Feb 2010
 ! use GridValues_ml,  only : sigma_bnd
 use ModelConstants_ml,    only : KMAX_MID, KWINDTOP
 use PhysicalConstants_ml, only : GRAV
 implicit none
 private


! Queries  a
! - should we keep this min value. Very high for stable?
! - max value also seems high. We used to use 2500 m

 real, parameter, public :: PZPBL_MIN=200.   ! Amela's
 real, parameter, public :: PZPBL_MAX=5000.  ! Amela's


public :: RiB_Hmix
!public :: Grisogono
public :: risig1

contains

subroutine RiB_Hmix (windspeed, zm, theta, pzpbl)
  real, dimension(KWINDTOP:KMAX_MID), intent(in) :: windspeed
  real, dimension(KMAX_MID), intent(in) :: zm, theta
  real, intent(out) :: pzpbl
   integer :: k, km
   real :: Ri  ! Richardson number
   real, parameter :: Ric = 0.25  ! critical Ric


!     Mixing height calculation
!     Boundary layer height is calculated with the bulk Richardson number
!     method with critical value of Ric=0.25
! Query. Why k-1 below. Also, can never give first later as Hmix. Is
! this intenional? What if we have a really stable BL?
! zm is geopotential height. Is this what we want? I think so.....
! z_mid is also 


            do k=KMAX_MID,KWINDTOP+1,-1
               km = k + 1
              ! Ric0=0.115*((zm(i,km)-zm(i,k))*100.)**0.175
              ! Ric=max(0.25,Ric0)

               Ri = risig1(windspeed(k-1), &
                              theta(KMAX_MID), theta(k-1), &
                              zm(k-1))
               if(Ri >= Ric) then
                  pzpbl = zm(k)
                  exit
               endif
            enddo

            pzpbl=min(pzpbl,pzpbl_max )
   !QUERY         pzpbl=max(pzpbl, pzpbl_min)
end subroutine RiB_Hmix


!subroutine Grisogono (nr,KZ_MINIMUM,KZ_MAXIMUM,zm,exns,exnm)
!
!   !Amela Jericevic 2008
!   !this routine uses 'exponential' generalization of O_Brian 3rd order
!   ! polynomial function for Kz calculation; Grisogono and Oerlemans, 
!   ! 2002, TellusA, 54, 453-462.
!   !//////////////////////////////////////////////////////////////////
!   !
!      implicit none
!      integer, intent (in) :: nr
!      real, intent(in) :: zm(MAXLIMAX, MAXLJMAX, KMAX_MID)
!      real, intent(in) :: exns(MAXLIMAX,MAXLJMAX,KMAX_BND)
!      real, intent(in) :: exnm(MAXLIMAX,MAXLJMAX,KMAX_MID)
!      real, intent(in) :: KZ_MINIMUM, KZ_MAXIMUM
!
!      real, dimension (MAXLIMAX,MAXLJMAX) :: Kmax,&
!                                             h, &
!                                             const,&
!                                             help, &
!                                             z
!!AJ const can go out, help added
!      real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: ri
!      real :: dvdz,risig1,fac,fac2,dex12,ro,us_3,jedan,dva
!!DS: Do jedan, dva have any meaning? Odd names?
!!AJ: It is one-jedan and two-dva in Croatian used as auxiliary variables
!!AJ  for testing purposes
!!AJ now changed
!
!      integer :: i,j,k,n,m
!      real, parameter :: Ric=0.25
!!     real, parameter :: EPS=0.001  !DS is this same as KZ_MINIMUM?
!				    !AJ yes, I have changed it in the code it can go out now
!!DS: new constants suggested by Amela, 7th Jan 2009:
!      real, parameter :: KZ_CONST = 0.1  !DS is this same as KZ_MINIMUM?
!                                         !AJ-No it isn't the same
!
!!DS Not used:
!!DS      real, parameter :: orog_coef=0.17
!!DS      real, parameter :: oro_max=3400.
!!DS      real, parameter :: param=1.8
!!AJ      I would like to remove this part now, we do not need it! I 
!!c       have tested it and it works fine with 0.1 constatnt
!
!      !linear function for variable const(i,j)needed in Kmax 
!!DS       !determination is used depends on orography
!!DS       !orog_coef,oro_max and param are variables of the streight line
!
!      do j=1,ljmax
!         do i=1,limax
!
!      !DS Omit orography corrections for now.
!      !DS linear function for Kmax determination is used which depends on orography
!      !DS  const(i,j)=(param/oro_max)*abs(h_oro(i,j))/GRAV+orog_coef
!     !initialization of surface roughness
!
!     !AJ-have tested and this can be removed     ustar_nwp(i,j)=max(ustar_nwp(i,j), 0.1)
!
!     !DS Query - why use min value of 0.1 - this is quite high?
!     ! elsewhere in the code, we have:
!     !  ustar_nwp(i,j) = max( ustar_nwp(i,j), 1.0e-5 )
!     ! and a consistency between tau and ustar. Here you reset ustar
!     ! and this may affect other routines?
!
!!     Mixing height calculation
!!     Boundary layer height is calculated with the bulk Richardson number
!!      method with critical value of Ric=0.25
!
!            do k=KMAX_MID,2,-1
!               ri(i,j,k) = risig1(u(i,j,k-1,nr), v(i,j,k-1,nr), &
!                              th(i,j,KMAX_MID,nr), th(i,j,k-1,nr), &
!                              zm(i,j,k-1))
!               if(ri(i,j,k) >= ric) then
!                  pzpbl(i,j) = zm(i,j,k)
!                  exit
!               endif
!            enddo
!
!            pzpbl(i,j)=min(pzpbl(i,j),pzpbl_max )
!            pzpbl(i,j)=max(pzpbl(i,j), pzpbl_min)
!!AJ You have this minimum pbl height limit in earlier version.
!!AJ I kept it becouse it might influence dry deposition.
!
!! Determination of two input parameters needed in Grisogono aproach
!
!            !DS Kmax(i,j)=const(i,j)*vsm(i,j)*ustar_nwp(i,j)
!            !AJ Yes, this is better
!            Kmax(i,j)=KZ_CONST*pzpbl(i,j)*ustar_nwp(i,j)
!            h(i,j)=max(pzpbl(i,j)/3.,zm(i,j,KMAX_MID-1))
!
! !AJ changed 3 to 4 to , result of tests in 2001 on NO2, SO2 and SO4
! ! Grisogono approach is valid for boundary layer, so z(i,j) is defined as
! ! upper boundary 
! ! above the PBL K(z) calculated with Blackadar approach is used 
!
!!DS - query -what did you mean by nz? z perhaps?
!!AJ -yes z not nz, typing mistake
!
!!DS - faster to write (?):
!!AJ -OK
!               z(i,j)= 1.25 * pzpbl(i,j)
!
!!Calculation of the exponential profile
!
!            do k=KMAX_MID,2,-1
!!AJ               jedan=Kmax(i,j)*zm(i,j,k)*sqrt(exp(1.))/h(i,j)
!!AJ               dva=exp(-0.5*(zm(i,j,k)/h(i,j))**2)
!                  xksig(i,j,k)=Kmax(i,j)*zm(i,j,k)*sqrt(exp(1.))/h(i,j)&
!                              *exp(-0.5*(zm(i,j,k)/h(i,j))**2)+KZ_MINIMUM
!               if(zm(i,j,k)>z(i,j)) exit
!            enddo
!         enddo
!      enddo
!
!      !
!      ! storing of pbl height
!      ! Why use vsm at all, could use pzpbl(i,j) from start?
!      !
!      !AJ changed vsm to pzpbl
!      
!      !smoothing of pbl height
!
!        call smoosp(pzpbl,pzpbl_min,pzpbl_max)
!      !
!      ! transformation of xksig(i,j,k) into skh(i,j,k,nr)
!      ! from HF Kz(sigma)=Kz*ro**2*(GRAV/p*)
!      !
!      do k=2,KMAX_MID
!     !AJ added smoothing of Kz with Shapiro -same as with O_Brien 
!      do i=1,limax
!           do j=1,ljmax
!             help(i,j)=xksig(i,j,k)
!           enddo
!      enddo
!
!      call smoosp(help,KZ_MINIMUM,KZ_MAXIMUM)
!         do i=1,limax
!            do j=1,ljmax
!
!!I have exns=exf1,exnm=exf2,
!               fac = GRAV/(ps(i,j,nr) - PT)
!               fac2 = fac*fac
!               dex12 = th(i,j,k-1,nr)*(exnm(i,j,k)-exns(i,j,k)) &
!                     + th(i,j,k,nr)*(exns(i,j,k)-exnm(i,j,k-1))
!
!               ro = ((ps(i,j,nr)-PT)*sigma_bnd(k)+PT)*CP &
!                  *(exnm(i,j,k)-exnm(i,j,k-1))/(RGAS_KG*exns(i,j,k)*dex12)
!
!               skh(i,j,k,nr) = xksig(i,j,k)*ro*ro*fac2
!            enddo
!         enddo
!      enddo
!   end subroutine Grisogono
 !----------------------------------------------------------------------------
  function risig1(ws,th1,th2,z)    ! Amela
   !calculates the bulk Richardson number
      implicit none
      real, intent(in) :: ws   ! wind-speed
      real, intent(in) :: th1, th2
      real, intent(in) :: z
      real :: risig1
      real :: dvdz

      dvdz = ws+0.001
      risig1=(2.0*GRAV/(th1+th2))*(th2-th1)*z/dvdz
   end function risig1
 !----------------------------------------------------------------------------

end module BLPhysics_ml
