! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!> Provides Kz from various methods
!! Main complexity is those different methods, which we deal with using
!! key-word arguments, e.g. if KzMethod constant, we need one argument
!! but if KzMethod is power-law we need more. 

module Kz_ml
use PhysicalConstants_ml, only : KARMAN, GRAV
implicit none
private

  !> From J-P matlab code

  public :: def_kz_nsl
  public :: def_kz_inc
  public :: def_kz_pow
  private :: non_rect_hyperbola

  !> From EMEP model

  ! - Hmix routines
   public :: Venkatram_Hmix
   public :: Zilitinkevich_Hmix

  ! Kz routines

   public :: BrostWyngaardKz
   public :: JericevicKzE
   public :: O_BrienKz

   ! Misc
   public :: risig1

  
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function def_kz_nsl(z,nsl,ustar) result(K)
 !  neutral Kz in SL, constant above
  real, intent(in), dimension(:) :: z
  real, intent(in) ::  ustar
  integer, intent(in) :: nsl
  real, dimension(size(z)) :: K
  real :: kappa = 0.4
 
  K = kappa*ustar*z
  K = min(K,K(nsl))  !QUERY???

end function def_kz_nsl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> Leuning et al AFM 2000 Kz method
!! NB - stability correction (use of 1/L) not implemented

function def_Kz_inc(z,uStar,hVeg,hSL,invL,kappa,dPerh) result (Kz)
 !  neutral Kz in SL, constant above
  real, intent(in), dimension(:) :: z   !> heights at which Kz calculated
  real, intent(in) ::  &
     ustar, hVeg       &   !> Obvious?
    ,hsl               &   !> Height surface layer (m) and  1/L
    ,invL              &   !> 1/L (L=Obukhov length, m) NB: no stability corr. yet
    ,dPerh                 !> = d/hVeg
  real, intent(in) :: kappa !> von Karman constant

  real, dimension(size(z)) :: Kz   ! Resulting Kz
  real, dimension(size(z)) :: zPerh, fac, sigma_w, tau_L, Krsl
  real :: d
  integer :: nsl, nz, k

 !> Leuning equation params
  real, parameter :: ZH0 = 0.8, ZH1 = 0.25  &
    ,a0 = 0.850, b0 = 1.25, c0 = 0.98, d0 = -1, e = 0.2, f = 1.5 &
    ,a1 = 0.256, b1 = 0.40, c1 = 0.98, d1 =  1, s1 = 0.8  &
    ,a2 = 0.850, b2 = 0.41, c2 = 0.98, d2 = -1, s2 = 4


   d = dPerh*hVeg
   nsl = count( z < hsl+0.01 )
   nz = size(z)
   zPerh(:) = z(:)/hVeg
   print *, "ZPERH ", nsl, zPerh, z(nsl)

   !> From: Leuning et al., AFM, 2000 %%%
   !> sigma_w:

   where ( zPerh < ZH0)
     fac = e*exp(f*zPerh)
   elsewhere
     fac = non_rect_hyperbola(zPerh,a0,b0,c0,d0)
   end where

   sigma_w = fac*uStar

  !> tau_L:
   where ( zPerh < ZH1 )
     fac = non_rect_hyperbola(s2*zPerh,a2,b2,c2,d2)
   elsewhere
     fac = non_rect_hyperbola(zPerh-s1,a1,b1,c1,d1)
   end where

   tau_L = fac*hVeg/uStar

   Krsl = sigma_w**2 * tau_L     !! in and just above vegetation 
   Kz = kappa*uStar*(z-d)        !! M-O similarity theory

   do k = 1, nsl                 !! QUERY? Not sure of matlab code
     print *, "KLEU ", k, nsl,  Kz(k), Krsl(k) 
     Kz(k) = max(Kz(k),Krsl(k))             !! match
   end do

   Kz(nsl+1:nz) = maxval( Kz(1:nsl) )  !> constant above SL

end function def_Kz_inc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!>  non-rectangular hyperbola

elemental function non_rect_hyperbola(x,a,b,c,d) result(y)
  real, intent(in) :: x,a,b,c,d
  real :: y
  real :: dum

  dum = a*x + b;
  y = (dum + d*sqrt(dum**2 - 4*a*b*c*x))/(2*c);
end function non_rect_hyperbola



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elemental function def_kz_pow(z,za,Ka,n) result(K)
 ! power law K(z) = Ka*(z/za)^n
  real, intent(in) :: z,za,Ka,n
  real :: K

  K = Ka*(z/za)**n;

end function def_kz_pow
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 !----------------------------------------------------------------------------

function BrostWyngaardKz(z,h,ustar,invL,Kdef) result(Kz)
  real, intent(in) :: z    ! height
  real, intent(in) :: h    ! Boundary layer depth 
  real, intent(in) :: ustar!  u*
  real, intent(in) :: invL !  1/L  
  real, intent(in) :: Kdef !  1/L  
  real :: Kz

     if ( z < h ) then
        Kz = KARMAN * ustar * z * (1-z/h)**1.5 / (1+5*z*invL)
     else
       Kz =  Kdef
     end if

end function BrostWyngaardKz

 !----------------------------------------------------------------------------

elemental function JericevicKzE(z,h,ustar,Kdef) result(Kz)
  real, intent(in) :: z    ! height
  real, intent(in) :: h    ! Boundary layer depth 
  real, intent(in) :: ustar, Kdef !  u*, default Kz
  real :: Kz
  real :: Kmax, zmax


     if ( z < h ) then
        Kmax = 0.05 * h * ustar
        zmax = 0.21 * h
        Kz = 0.39 * ustar * z * exp( -0.5*(z/zmax)**2 )

     else ! open-source had this Kz=0.0 line. Not sure why
       Kz =  Kdef
     end if

end function JericevicKzE


!----------------------------------------------------------------------------
!>  From Venkatram, 1980, simple method!
!!
function Venkatram_Hmix (ustar) result(zi)
  real, intent(in) :: ustar
  real :: zi

     zi = 2.4e3 * ustar**1.5

end function Venkatram_Hmix
 !----------------------------------------------------------------------------

function Zilitinkevich_Hmix (ustar,invL,lat) result(zi)
 ! Not equator-proof yet!!
   real, intent(in) :: ustar, invL, lat
   real :: f, pi = 4.0*atan(1.0)
   real :: zi
   f = 1.46e-4 * sin(lat*pi/180.0)
   !zi = 0.2 *  ustar/f
   if( invL > f/ustar ) then
     zi = 0.4 * sqrt( ustar/(f*invL) )
   else
     zi = 0.2 *  ustar/f
   end if
end function Zilitinkevich_Hmix

 !----------------------------------------------------------------------------


 !************************************************************************!
   subroutine O_BrienKz( Hs, zi, zs_bnd, ustar, invL, &
        KzHs, KzHmix, Kz, debug_flag )       !
 !**********************************************************************!
 !..exchange coefficients for convective boundary layer:
 !..o'brien's profile formula:

    real,intent(in) :: zi, Hs   ! zi=Hmix, ESX: Hs = height surface layer
    real,intent(in), dimension(:) :: zs_bnd
    real,intent(in) :: ustar, invL
    real, intent(in) :: KzHs, KzHmix   !! Kz  at top of surf. and mixing layer
    real, intent(inout), dimension(:) :: Kz
    logical, intent(in) :: debug_flag

    !ESX real, parameter :: FHSL = 0.04 ! surface layer as fraction zi

    real :: ux3 , hsl, zimhs ,zimz ,zmhs
    real  ::  dKzdz ,Kzzi  !ESX  ,hs
    integer :: k

    !c..exchange parameter and its vertical derivative at z = Hs

      !ESX Kzhs=0.    ! Kz at hs
      dKzdz=0.   ! d Kz(Hs) /dz - derivateive at Hs (or KB)
      Kzzi=0.    ! Kz at top of BL, set to zero or KZ_MINIMUM?

      !...................................................................
      !..air density at ground level is always calculated diagnostically:

      ux3 = ustar*ustar*ustar

     !..........................
     !..unstable surface-layer.:

     !ESX Hs=FHSL*zi          !  height of surface layer

     !c..hsl=Hs/l where l is the monin-obhukov length
     !hsl = KARMAN*GRAV*Hs*fh*KAPPA /(ps*ux3)

     hsl = Hs*invL

     ! Garratt \Phi function - take from MicroMet

     !Kzhs = ustar*KARMAN*Hs*sqrt(1.0-16.0*hsl)  ! /Pr=1.00

     ! In ESX we want to match Kz at Hs. Could scale dKzdz also, but
     ! probably we can assume the gradient is okay.

     dKzdz = Kzhs*(1.-0.5*16.0*hsl/(1.0-16.0*hsl))/Hs

!ESXQUERY     Kz(KMAX_MID)=Kzhs  ! QUERY - should be at z_bnd(20)?
     if ( debug_flag ) write(*,"(a,f7.2,3es12.3)") "OBRIEN Kz20 ", &
        Hs, invL, Kzhs, dKzdz
     if ( 16.0*hsl >= 1.0 ) then
        write(*,"(a,f7.2,4es12.3)") "OBRIEN NEG ", Hs, hsl,  Kzhs, dKzdz
     end if

   !..exchange parameter at z = ziu

   do  k= 1, size(Kz) !! 2,KMAX_MID

         if( zs_bnd(k) >= zi) then

            Kzzi=Kz(k)  ! values above zi  stored

            if ( debug_flag ) write(*,"(a,i3,es12.3)") "OBRIEN Kzzi ", k,  Kzzi

         else 
            !.....................................................
            !..the obrien-profile for z<ziu                      .
            !.....................................................
            !
            if(zs_bnd(k) <= Hs) then
               Kz(k)=zs_bnd(k)*Kzhs/Hs
               if ( debug_flag ) &
                   write(*,"(a,i3,es12.3)") "OBRIEN Kzhs ", k, Kz(k)
            else !! if( zs_bnd(k) < zi) then
               zimhs = zi-Hs
               zimz  =zi-zs_bnd(k)
               zmhs  =zs_bnd(k)-Hs
               Kz(k) = Kzzi+(zimz/zimhs)*(zimz/zimhs)  &
                    *(Kzhs-Kzzi+zmhs*(dKzdz     &
                    + 2.*(Kzhs-Kzzi)/zimhs))
               if ( debug_flag ) &
                   write(*,"(a,i3,es12.3)") "OBRIEN Kz(k) ", k,  Kz(k)
            endif

         endif

   end do

  end subroutine O_BrienKz
 !----------------------------------------------------------------------------
  function risig1(ws,th1,th2,z) ! Amela Jericevic, not used in new version
   !calculates the bulk Richardson number
      implicit none
      real, intent(in) :: ws   ! wind-speed
      real, intent(in) :: th1, th2
      real, intent(in) :: z
      real :: risig1
      real :: dvdz

      dvdz = ws*ws+0.001
      risig1=(2.0*GRAV/(th1+th2))*(th2-th1)*z/dvdz
   end function risig1
 !----------------------------------------------------------------------------

end module Kz_ml
