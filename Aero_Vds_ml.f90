!=============================================================================
  module Aero_Vds_ml
!==============================================================================
  use PhysicalConstants_ml, only : FREEPATH, VISCO, BOLTZMANN, PI, GRAV, ROWATER
  use My_Aerosols_ml,       only : NSIZE, FINE_PM, COAR_PM, GIG_PM
  use ModelConstants_ml,    only : DEBUG_VDS
 
  ! DESCRIPTION
  ! Calculates laminar sub-layer resistance (rb) and gravitational settling 
  ! velocity (vs) for particles
  ! Finally: vd= vs+1/(ra+rb+ra*rb*vs),   where
  ! vs - gravitational settling velocity,
  ! ra - aerodynamic resistance, rb - viscous sub-layer resistance,
  ! rc - surface resistance (assumed zero for particles)
  !---------------------------------------------------------------------------

 implicit none
  private

!  public  :: Aero_Vds
     ! A variety of equations are presented. Two types of
     ! stability corrections are used in the literature,
     ! those based upon 1/L only, and those based upon
     ! zi/L, where zi is the PBL height. 
     !
     ! Using zi form gives max stab-fac of ca. 10 
     ! (range max 2.8 to 10.7 for zi = 200 - 2500)
     ! Using 300 form gives max stab-fac of ca. 6

     ! For EMEP, use:
     ! WT versions : zi/L and limit 1/L to -0.04, cf Pryor

  public  :: SettlingVelocity
  public  :: PetroffFit
  public  :: GPF_VdsZi  ! General, Gallager-Petroff fits (loosely!)
  public  :: GPF_Vds300 !  "  "
  public  :: WeselyWT
  public  :: Wesely300
  public  :: Gallagher1997
  public  :: Gallagher2002
  public  :: GallagherWT
  public  :: Nemitz2004
  public  :: RuijgrokDrySO4
  public  :: RuijgrokWetSO4
  public  :: Wesely1985
!TMP  public  :: McDonaldForest

!  public  :: Slinn - to be re-added from Aero_Rb, after mods

contains

   !------------------------------------------------------------------------
   function SettlingVelocity(tsK,roa) result(Vs)
    ! gravitational settling for 2 modes
    ! Equations Axx referred to here are from Appendix A,
    !  Binkowski+Shankar, JGR, 1995
    !  Note confusing notation in B+S, dp was used for diff. coeff.
    !  here we use Di

     real, intent(in) :: tsK, roa ! temp, air density, Rel.hum
!     integer, parameter :: NSIZE = 2  ! In My_Aerosols
     real  :: Vs(NSIZE)

   !QUERIES should FREEPATH, VISCO be constants? 

   ! here we use dp=0.33 equivalent to McDonald et al.2007, Table 3 
   ! data for accumulation mode, background.

     real, parameter, dimension(NSIZE) ::   &
                 diam   = (/ 0.33e-6, 4.0e-6, 8.5e-6 /),  &
                 sigma  = (/ 1.8, 2.0, 2.2 /),                    &
                 PMdens = (/ 1600.0, 2200.0, 2200.0 /) ! kg/m3
     real, parameter :: one2three = 1.0/3.0
     integer :: imod 
     real    :: lnsig2, dg, dg_dry, r_cm, rw, Rh, GF3, part_dens, & !slip, &
                knut, Di1, Di, vind, &
                stoke, schmidt, &  ! Stoke's and Schmidt numbers
                vsmo, vs1          ! Settling velocity


    do imod = 1, NSIZE

        lnsig2 = log(sigma(imod))**2

       !... mass median diameter -> geometric diameter 

        dg = exp (log(diam(imod)) - 3.* lnsig2 )

        knut = 2.0*FREEPATH/dg   ! Knut's number

!... slip correction coefficient  
!	slipmo= 1.+ knut*       &               ! for monodisperse
!                (1.257+0.4*exp(-1.1* /knut))
!NOT-USED?	slip =  1.+ 1.246*knut                  ! for polydisperse

!== monodisperse aerosols =====
!     Dimo =BOLTZMANN*tsK*slipmo/(3*PI*dg *VISCO*roa)        ! diffusion coefficient
!     vsmo =dg*dg *PMdens(imod) *GRAV*slipmo/(18.*VISCO*roa) ! gravitational settling

!== polydisperse aerosols (log-normal size distribution) =====

        Di1 =BOLTZMANN*tsK/(3*PI*dg *VISCO *roa)                    ! A30, dpg

        ! diffusion coefficient:
        Di = Di1*(exp(-2.5*lnsig2)+1.246*knut*exp(-4.*lnsig2))      ! A29, dpk 

        vs1= dg*dg * PMdens(imod) * GRAV / (18.0* VISCO*roa)           ! A32
        vs(imod) = vs1*(exp(8.0*lnsig2)+1.246*knut*exp(3.5*lnsig2)) ! A31, k=3

        if ( DEBUG_VDS ) write (6,'(a19,i3,f8.3)'), "** Settling Vd **",imod, vs(imod)*100.0

     end do !imod
   end function SettlingVelocity

! -------------------------------------------------------------------
   function PetroffFit(ustar,invL,SAI) result(Vds)
     ! "Simple" fitting of Petroff et al.
     ! Fig12b suggests  Vd = 0.3 cm/s for u* = 0.45, LAI=22
     ! ->  Vds = 0.007 * u*
     ! Fig.15 suggests that vds is approx prop to LAI for dp~0.4um
     ! We use SAI to keep some winter dep in decid forests
     ! To keep Vds above grassland, we use max(3.0, SAI)
     real, intent(in) :: ustar, invL,SAI
     real :: Vds

        Vds   = 0.007 * ustar * 0.1*max(SAI, 3.0) 

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-300.0 * max(-0.04,invL))**0.6667)
       end if
   end function PetroffFit
   !------------------------------------------------------------------------
   ! GallagherPetrof fits
   ! Two functions here, for different stability methods
     ! "Simple" fitting of Gallagher et al. (1997) and Petroff et al., which 
     ! roughly captures the differences between Speulderbos-type and typical
     ! forests, because of LAI. 
     ! Gallagher et al. had   Vds/u* = 0.0135 * Dp * stab function
     ! which gives 0.3 cm/s for neutral conditions, Dp=0.5
     !
     ! Fig12b suggests  Vd ~ 0.3 cm/s for u* = 0.45, LAI=22
     ! ->  Vds = 0.008 * u*
     ! Fig.15 suggests that vds is approx prop to LAI for dp~0.5um
     ! We use SAI to keep some winter dep in decid forests
     ! As Petroff started with a total LAI of 22, which is ca.
     ! 1-sided LAI=10, SAI=11, so we scale with SAI/11 = 0.09
     ! 
     ! We also limit the lowest Vds/u* to be 0.002, consistent with
     ! Wesely.

   function GPF_VdsZi(ustar,invL,SAI,zi) result(Vds)
     real, intent(in) :: ustar, invL,SAI, zi
     real :: Vds

        Vds   = ustar * max( 0.002, 0.008 * 0.1 * SAI )

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-0.3*zi*max(-0.04, invL))**0.6667)
       end if
   end function GPF_VdsZi
   !------------------------------------------------------------------------
   function GPF_Vds300(ustar,invL,SAI) result(Vds)
     real, intent(in) :: ustar, invL,SAI
     real :: Vds

        Vds   = ustar * max( 0.002, 0.008 * 0.1 * SAI )

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-300.0 * max(-0.04,invL))**0.6667)
       end if
   end function GPF_Vds300
   !------------------------------------------------------------------------
   function Wesely1985(ustar,invL, zi) result(Vds)
     real, intent(in) :: ustar, invL, zi
     real :: Vds

        Vds   = 0.002 * ustar  ! from grass

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-0.3*zi*max(-0.04, invL))**0.6667)
         !Alt: Vds = Vds * 0.0009 ( -zi*invL)**0.6667
       end if
   end function Wesely1985
   !------------------------------------------------------------------------
   function WeselyWT(ustar,invL, zi) result(Vds)
     ! Same as Wesely1985, but with 1/L limit
     real, intent(in) :: ustar, invL, zi
     real :: Vds

        Vds   = 0.002 * ustar  ! from grass

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-0.3*zi*max(-0.04,invL))**0.6667)
       end if
   end function WeselyWT
   !-- ----------------------------------------------------------------------
   function Wesely300(ustar,invL) result(Vds)
     ! Same as Wesely1985, but with 1/L limit and 300/L
     ! version of stability fac
     real, intent(in) :: ustar, invL
     real :: Vds

        Vds   = 0.002 * ustar  ! from grass

       if ( invL <  0.0 ) then
         Vds = Vds * (1.0+(-300.0 * max(-0.04,invL))**0.6667)
       end if
   end function Wesely300
   !-- ----------------------------------------------------------------------
   function RuijgrokDrySO4(ustar,RH,u_hveg) result(Vds)
     real, intent(in) :: u_hveg, RH, ustar ! RH in %
     real :: Vds

        Vds = 0.05*ustar**0.28 / u_hveg

        if ( RH > 80.0  ) then
           Vds = Vds * ( 1+ 0.18*exp( (RH-80)/20.0 ) )
        end if
   end function RuijgrokDrySO4
   !------------------------------------------------------------------------
   function RuijgrokWetSO4(u_hveg,RH,ustar) result(Vds)
     real, intent(in) :: u_hveg, RH, ustar ! RH in %
     real :: Vds

        Vds = 0.08*ustar**0.45 / u_hveg

        if ( RH > 80.0  ) then
           Vds = Vds * ( 1+ 0.37*exp( (RH-80)/20.0 ) )
        end if
   end function RuijgrokWetSO4
   !------------------------------------------------------------------------
   function Nemitz2004(dp,ustar,invL) result(Vds)
     real, intent(in) :: dp, ustar, invL
     real :: Vds

        Vds = 0.001*ustar 

        if ( invL < 0.0 ) then
           Vds = Vds *( 1+( -(960*dp-88.0)*invL )**0.6667)
        end if
   end function Nemitz2004 
   !------------------------------------------------------------------------
   function Gallagher1997(dp,ustar,invL) result(Vds)
     real, intent(in) :: dp, ustar, invL
     real :: Vds

        Vds = 0.0135 * ustar * dp 

        if ( invL < 0.0 ) then
           Vds = Vds * (1.0+(-300*invL)**0.6667 )
        end if
   end function Gallagher1997
   !------------------------------------------------------------------------
   function Gallagher2002(ustar,invL,z0) result(Vds)
     real, intent(in) :: ustar, invL, z0
     real :: Vds
     real :: k1, k2

     !if( log(z0) > 0.0 ) then ! z0 > ~0.04 m
     !if( log10(z0) > 0.0 ) then ! z0 > ~0.04 m
       !k1 = 0.001222 * log(z0) + 0.003906 
       k1 = 0.001222 * log10(z0) + 0.003906  ! Eqn (13)

     ! This equation has negative solutions. We set
     ! a small min value,  consistent wth 0.2 m/s from G02, Fig
     ! and u* = 0.5 m/s.
       k1 = min( 0.0004, k1)
     !else
     !  k1 = 0.001222
     !end if
     k2 = 0.0009

        if ( invL < 0.0 ) then
           Vds = ustar * &
          ( k1 + k2* (-300*invL)**0.6667 )
        else
           Vds = ustar * k1
        end if
   end function Gallagher2002
   !------------------------------------------------------------------------
   function GallagherWT(dp,ustar,invL, zi) result(Vds)
     ! Same as Gallagher1997, but with Wesely's zi form of
     ! stability correction, and 1/L limit
     real, intent(in) :: dp, ustar, invL,  zi
     real :: Vds

        Vds = 0.0135 * ustar * dp 

        if ( invL < 0.0 ) then
          Vds = Vds * (1.0+(-0.3*zi*max(-0.04,invL))**0.6667)
        end if
   end function GallagherWT
   !-- ----------------------------------------------------------------------
!TMP   function McDonaldForest(dpMed,ustar) result(Vds)
!TMP     real, intent(in) :: dpMed, ustar
!TMP     real :: Vds
!TMP     real :: X
!TMP     integer :: i
!TMP     real,dimension(7) :: a! = (/ 0.0000614,0.0012994,0.0023525, &
!TMP     a = (/ 0.0000614,0.0012994,0.0023525, &
!TMP                                -0.0647616,-0.0794396,1.2454391,3.8075140 /)
!TMP     Vds = a(7)
!TMP     X = log(dpMed)
!TMP     do i = 1, 6
!TMP         Vds = Vds + a(7-i) * X**i
!TMP     end do 
!TMP     Vds = 0.001 * exp(Vds)
!TMP
!TMP   end function McDonaldForest


     
   !-- ----------------------------------------------------------------------
! =================================================================

end module Aero_Vds_ml

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!=============================================================================
! From Aero_Rb, should be developed in Slinn function
! Need to reconsider snow, conv etc.
! 
!TMP! -------------------------------------------------------------------
!TMP
!TMP!// Stokes and Schmidt numbers:
!TMP
!TMP   ! == monodisperse ======
!TMP	! STmo=vsmo*ustar*ustar/VISCO/GRAV
!TMP	! SCmo=VISCO/dimo
!TMP   ! == polydisperse ======
!TMP	 schmidt = VISCO/Di              ! Schmidt number
!TMP	 stoke = vs(imod)*ustar*ustar/VISCO/GRAV ! Stoke number(based on depth 
!TMP                                              ! of laminar layer)
!TMP
!TMP !// collection efficiency  =======================
!TMP !      coleff=1./sc**(2./3.) + 1./10.**(3/stoke)
!TMP !=================================================
!TMP
!TMP         vind = max( 0.005, v50 )   ! wind at 50m height
!TMP
!TMP
!TMP
!TMP      if( LandType(lu)%is_water )    then !//===  WATER surface  ( Slinn & Slinn, 1980 ) =
!TMP
!TMP      	   coleff= ustar / (KARMAN * vind) *      &          ! polydisperse
!TMP                  (exp(-0.5*log(schmidt)) + exp(-3./stoke*log10) )   
!TMP
!TMP      elseif ( LandType(lu)%is_conif )  then !//===  CONIFEROUS ==============
!TMP
!TMP           stoke = vs(imod)*ustar/(AHAT*GRAV)   ! vegetation (Slinn, 1982)
!TMP           coleff= exp(-2./3.*log(schmidt)) + stoke/(1.+ stoke*stoke)  !  Slinn 
!TMP
!TMP      elseif ( LandType(lu)%is_veg ) then !//===  other VEGETATIVE surfaces ======
!TMP
!TMP          if ( snow > 0 .or. tsK <= 273.)   then   !... covered with snow or frozen 
!TMP
!TMP              coleff= exp(-2./3.*log(schmidt)) + exp(-3./stoke*log10)  ! polydisperse  
!TMP
!TMP          else   !... snowfree  
!TMP
!TMP               stoke = vs(imod)*ustar/(AHAT*GRAV)     ! Stoke for vegetation (Slinn, 1982)
!TMP               coleff= exp(-2./3.*log(schmidt)) + stoke/(1.+ stoke*stoke)  !  Slinn 
!TMP          endif
!TMP
!TMP      else    !//====  urban/desert/ice always ===============================
!TMP              !..   Slinn at al(1978), Seinfeld(1997), Binkowski 
!TMP
!TMP           coleff= exp(-2./3.*log(schmidt)) + exp(-3./stoke*log10)  ! polydisperse
!TMP
!TMP      endif        
!TMP
!TMP! .. laminar layer resistance .....................................
!TMP    !  rb= 1./ustar/colef                  !     Seinfeld
!TMP    !  rb= 0.4*vind/(ustar*ustar*colef)      !     Slinn
!TMP
!TMP!// ==  bounce-off for coarse particles (Slinn, 1982) ===============
!TMP   !... for fine aerosol
!TMP
!TMP       reb = 1.
!TMP 
!TMP   !... for coarse aerosol 
!TMP
!TMP      if(imod == NSIZE )  then
!TMP
!TMP        if( .not. (LandType(lu)%is_water .and. wetarea == 0.0))  & ! not on water/wet surface
!TMP
!TMP            reb = max (1.e-7, exp(-2. * sqrt(stoke)))           
!TMP      endif
!TMP
!TMP!//== enhanced dry.dep under convective conditions (Wesely at al,1985) ====
!TMP
!TMP      convfac = 0. 
!TMP     
!TMP!  only for (low) vegetation 
!TMP
!TMP      if (LandType(lu)%is_veg .and. snow == 0)  convfac = conv  
!TMP
!TMP!// == sub-laminar layer resistance ====================
!TMP
!TMP     rb (imod)  = 1. /(ustar*(1.+ 0.24 *convfac)) /(coleff*reb)
!TMP     rbw (imod) = 1. /(ustar*(1.+ 0.24 *convfac)) /coleff   ! no bounce-off on
!TMP                                                            ! wet surfaces 
!TMP
!TMP!..monodisperse: rb1= 1./(ustar*(1. + 0.24 * conv))/(colef1*reb) 
!TMP
!TMP       end do MODEloop
!=============================================================================
