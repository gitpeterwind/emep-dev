module Rsurface_ml

!================== Now under CVS control =================
! $Author: mifads $
! $Id: Rsurface_ml.f90,v 1.13 2003-03-24 14:07:09 mifads Exp $
! $Name: not supported by cvs2svn $
! =========================================================

!Changes from Version 1.0
!d1.1 - zenith angle range now to 89.9, not 85.0
!d1.3 - sinB now set to coszen, not 1/coszen
!d1.4 - restructered for use in full EMEP model.
!d1.7 - SO2 from CEH. Not LAI dependent, so a modified procedure
!       is needed for other gases.

use DepVariables_ml, only : forest, g_pot, g_temp,g_vpd,g_light,g_swp, &
                              b_inc     , albedo    ,   &
                              SAIadd, &
                              water, &     !ds rv1.8
                              g_max     , g_min     , g_lightfac,   &
                              g_temp_min, g_temp_opt, &
                              RgsS      , RgsO      , RextS, RextO, &
                              Rgs, Rinc, Gext, &
                              VPD_max   , VPD_min   , &
                              SWP_max   , PWP       , rootdepth 

use My_UKDep_ml, only : NDRYDEP_CALC, & ! no. of Wesely species  used
                         DRYDEP_CALC    ! array with species which are needed

                    
use PhysicalConstants_ml, only : KARMAN
use SoilWater_ml, only : SWP
use Wesely_ml, only  : Wesely_tab2, & ! Wesely Table 2 for 14 gases
                       WES_HNO3, & ! indices to identify HNO3
                       WES_NH3,  & ! indices to identify NH3  
                       WES_NO2 , & ! indices to identify NO2 
                       Rb_cor,  &! correction factor used in evaluating Rb
                       DRx       !  Ratio of diffusivities to ozone
implicit none
private

public   ::  Rsurface
!d1.7  public   ::  Conif_gpot
private  ::  g_stomatal    
private  ::  get_glight

logical, private, parameter :: DEBUG_RSURF = .false.
 
    
contains
! =======================================================================

  subroutine Rsurface(rh,lu,debug_flag, LAI,hveg,&
                      z0,ustar,Ts_C,vpd,SWP, &
                      psurf, precip, &                    !u7.lu
                      coszen, Idfuse, Idrctt, &       !u7.lu
                      snow, &                    !u7.lu
                      so2nh3ratio, &        !so2/nh3 ratio
                        g_sto,Rsur,Rsur_wet,Rb)
! =======================================================================
!
!     Description
!       calculates bulk surface resistance according to methodology
!       derived from EMEP MSC_W Note 6/00; the following pathways
!       apply for the surface resistance:
!
!        -- Rinc-- Rgs      In-canopy + soil/ground cover
!       |
!       |
!        -------- Rext      cuticular+other external surface
!       |
!       |
!        -------- Rsto      Stomatal
!
! with
!
!                           1
!  Rsur =           ____________________________
!                   LAI +  SAI  +       1
!                   ___    ___     ____________
!                   Rsto   Rext    Rinc + Rgs
!
! and we also calculate Rb here, since it (through u*) is land-use dependent 
!
!    Other gases
!
!  So far work has concentrated on ozone. However, other gases are 
!  considered using Wesely's idea of deriving resistances for ozone 
!  and SO2 first (e.g. RgsO, RgsS for Rgs) and then scaling using
!  effective Henry coefficients (H*) and reactivity coefficients (f0)
!  d1.7 - scaling applied to Gns, not individual resistances
!
!
! Structure of routine
!
!  1. Calculate:
!        Rlow             low-temperature correction
!        Rinc             in-canopy resistance
!        Rsur(HNO3)  
!        Gsto(O3)         stomatal conductance (if LAI > 0)
!
!       FOR EACH remaining gas (icmp is used as an index, since cmp is assumed 
!                               to  abbreviate "component".):
!  2. Calculate ground surface resistance, Rgs
!              if (LAI<0.1) go to 4  (for snow/ice/water...)
!  3. if (LAI>0.1)  calculate Gext
!  4. Calculate Rsur(icmp)
!       END
!
! =======================================================================

!......................................
! Input:

    real, intent(in) ::  rh             !hf d1.7 ddep
    integer, intent(in) :: lu           ! land-use index
    logical, intent(in) :: debug_flag   ! set true if output wanted 
                                        ! for this location
    real, intent(in) :: LAI             ! leaf area index (m2/m2)
    real, intent(in) :: hveg            ! height of vegetation (variable)
    real, intent(in) :: z0              ! vegetation roughness length (in m)
    real, intent(in) :: ustar           ! friction velocity (m/s)
    real, intent(in) :: Ts_C            ! surface temp. (degrees C)
    real, intent(in) :: vpd             ! vapour pressure deficit (kPa)
    real, intent(in) :: SWP             ! soil water potential (MPa)
    real, intent(in) ::  psurf
    real, intent(in) ::  precip         !acc precip at surface
    real, intent(in) ::  coszen
    real, intent(in) ::  Idfuse
    real, intent(in) ::  Idrctt
    integer, intent(in) :: snow         ! snow=1, non-snow=0
    real, intent(in) ::  so2nh3ratio    !so2/nh3 ratio

! Output:

   real, intent(out)             :: g_sto !  Stomatal conducatance (s/m)
   real,dimension(:),intent(out) :: Rsur ! bulk canopy surface resistance (s/m)
   real,dimension(:),intent(out) :: Rsur_wet !    " " for wet surfaces
   real,dimension(:),intent(out) :: Rb   ! quasi-laminar boundary layer
                                          ! resistance (s/m)


! Working values:
   
    integer :: icmp             ! gaseous species
    integer :: iwes             ! gaseous species, Wesely tables
    logical :: canopy,        & ! For SAI>0, .e.g grass, forest, also in winter
         leafy_canopy           ! For LAI>0, only when green
    real :: SAI                 ! Surface area index (m2/m2)
    real, parameter :: SMALLSAI= 0.05  ! arbitrary value but small enough
    real :: Hstar, f0           ! Wesely tabulated Henry's coeff.'s, reactivity
    real :: Rlow                ! adjustment for low temperatures (Wesely,
                                ! 1989, p.1296, left column) 
    real :: xRgsS, xRgsO        ! see  DepVariables_ml
    real :: xRgsS_wet, Rgs_wet  !  ??????
    real :: Ggs                 ! ground surface conductance, any gas
    real :: r_water             !hf CEH: non stomatal resistance for NH3 
    real :: Rhlim               !Rh limitation for wetness
    real :: so2nh3              !so2/nh3 after correction with 0.6
    real, parameter :: D_H2O = 0.21e-4  ! Diffusivity of H2O, m2/s
                                        ! From EMEP Notes
    real            :: D_i              ! Diffusivity of gas species, m2/s
!dep1_8  real            :: wRextS           !  TMP - crude wet Rext for S

!CEH resistances for SO2 in low NH3 conditions
     real, save :: CEHd = 180.0, CEHw = 100.0  !  dry, wet, m/s
     real :: cehfac   ! factor to interpolate between the 2,
                      ! crudely introduced by ds, Jan 2003....!


!dep1_7 : to interpolate CEH and O3 methods:
    real :: GnsO, GnsS_dry, GnsS_wet, Gns_dry, Gns_wet, GgsO


! START OF PROGRAMME: 


  !/** Do Rb first:
  !===========================================================================
 
  GASLOOP1: do icmp = 1, NDRYDEP_CALC
      iwes = DRYDEP_CALC(icmp)          ! index in Wesely table
                  
      !ds rv1.8 if   ( hveg  >=  0.0 ) then

    ! water, sea, rb is calculated as per Hicks and Liss (1976)
      if   ( water(lu) ) then

          D_i = D_H2O / Wesely_tab2(1,iwes)  ! CORR !

         !CORR 19/8/2002 - spotted by pw. Corrected diffc to D_i by ds...

          Rb(icmp) = log( z0 * KARMAN * ustar/ D_i )
          Rb(icmp) = Rb(icmp)/(ustar*KARMAN)

         ! CORR - Rb can be very large or even negative from this
         !        formula. We impose limits:

          Rb(icmp) = min( 1000.0, Rb(icmp) )    ! CORR - - gives min 0.001 m/s!
          Rb(icmp) = max(   10.0, Rb(icmp) )    ! CORR - - gives max 0.10 m/s!

      else 

          Rb(icmp) = 2.0 * Rb_cor(iwes)/(KARMAN*ustar)
      end if

      if ( DEBUG_RSURF .and. Rb(icmp) < 0.0  ) then
        print *, "ERROR RSURFB", lu, icmp, iwes, Rb_cor(iwes), Rb(icmp)
        return          
      end if

   end do GASLOOP1


  !/** Now begin Rsur calculations:
  !===========================================================================
  !/**  Adjustment for low temperatures (Wesely, 1989, p.1296, left column)

    Rlow = 1000.0*exp(-Ts_C - 4.0)
    Rlow = min(Rlow,9.9e4)    



!##############   1. Calculate In-Canopy Resistance, Rinc    ################

    canopy       = ( LAI+SAIadd(lu) > SMALLSAI ) ! - can include grass
    leafy_canopy = ( LAI     > SMALLSAI )        ! - can include grass

  !/** For canopies:

  !/** Calculate stomatal conductance if daytime and LAI > 0
    !d1.6 - set limit of coszen>0.001 which coresponds to zen = 89.94 degs.

   if( leafy_canopy  .and. coszen > 0.001 ) then  ! Daytime 

        call g_stomatal(debug_flag,lu,coszen,Idrctt,Idfuse, & 
                         Ts_C,psurf,LAI,vpd,SWP,g_sto)
   else
        g_sto = 0.0    ! Dark, no stomatal conductance...

   end if ! leafy canopy and daytime


  !/** Calculate Rinc, Gext 
  !NB-    previous if-test for snow replaced by more efficient
  !       multiplication, since snow=0 or 1.

   if(  canopy ) then   
         SAI  = LAI + SAIadd(lu)                  ! Accounts for bark/twigs
         Rinc = b_inc(lu)* SAI * hveg  / ustar    ! Lisa's eqn. 48 for u*=0.5
         Rinc = min( Rinc, 1000.0)

     ! TEST ds
         cehfac = 1.0
         so2nh3  = so2nh3ratio*0.6 !0.6=correction for local nh3
         if(so2nh3  < 2.0 ) cehfac = exp( -(2.0-so2nh3) )
     ! TEST ds
         xRgsS     = CEHd * cehfac  + Rlow  + snow * 2000.0
         xRgsS_wet = CEHw * cehfac  + Rlow  + snow * 2000.0

        ! for now, use CEH stuff for canopies, keep Ggs for non-canopy

         GnsS_dry = 1.0 /  xRgsS       ! For SO2, dry, low NH3 region
         GnsS_wet = 1.0 /  xRgsS_wet   ! For SO2, wet, low NH3 region

   else   ! No canopy present

        SAI  = 0.0     ! dep1.7
        Rinc = 0.0     ! dep1.8 -999.9
        !d1.7  Gext = -999.9  !  (not needed maybe, but...)

        !/ Here we preserve the values from the ukdep_gfac table
        !  giving higher deposition to water, less to deserts

        xRgsS     = RgsS(lu) + Rlow  + snow * 2000.0
        xRgsS_wet = xRgsS    ! Hard to know what's best here
      
   end if !  canopy

   !dep1_8 : For high RH we allow the canopy to take on wet characteristics.
   !         ds  introduced square function for RH>0.95, which has the value
   !         zero at RH=0.95 and 1 at RH=1.0:
   !hf Changed Rhlim according to Klemm. Not square anymore.

           cehfac = 0.0
           if (forest(lu)) then
              Rhlim=0.85     ! ds bug 85.
           else
              Rhlim=0.75     ! ds bug 75.
           endif
           if ( rh > Rhlim ) then
               cehfac = ( (rh-Rhlim)/(1.0-Rhlim) )
           end if


!####   2. Calculate Surface Resistance, Rsur, for HNO3 and  
!####      Ground Surface Resistance, Rgs, for the remaining 
!####      Gases of Interest                                

   !/ Ozone values....

     xRgsO  = RgsO(lu) + Rlow  + snow *    0.0    !u7.lu QUERY???
     GnsO   = SAI/RextO + 1.0/( xRgsO + Rinc ) ! (SAI=0 if no canopy)


!.........  Loop over all required gases   ................................

  GASLOOP2: do icmp = 1, NDRYDEP_CALC

     !-------------------------------------------------------------------------

     !  code obtained from Wesely during 1994 personal communication
     !  but changed (ds) to allow Vg(HNO3) to exceed Vg(SO2)

        if ( DRYDEP_CALC(icmp) == WES_HNO3 ) then
            !ds dep1_8 Rsur(icmp)  = max(10.0,Rlow)
            Rsur(icmp)  = max(1.0,Rlow)
            Rsur_wet(icmp)  = Rsur(icmp)
            cycle GASLOOP2
        end if

     !-------------------------------------------------------------------------
     ! Calculate the Wesely variables Hstar (solubility) and f0 (reactivity)

        iwes = DRYDEP_CALC(icmp)
        Hstar =Wesely_tab2(2,iwes)    !Extract H*'s 
        f0    =Wesely_tab2(5,iwes)    !Extract f0's
    
     !-------------------------------------------------------------------------

                          

     !   Use SAI to test for snow, ice, water, urban ...

       if ( canopy  ) then   

         ! ###   3. Calculate Cuticle conductance, Gext   ################
         ! ###      and  Ground surface conductance Ggs:

         ! Corrected for other species
         ! using Wesely's eqn. 7 approach. (We identify leaf surface resistance
         ! with Rext/SAI.)


    !dep1.8
    ! In earlier versions both O3 and SO2 had the same basic formulation, so
    ! that :
    !     Rsur(icmp) = 1.0/( LAI*DRx(iwes) *g_sto + SAI*Gext + Ggs )
    !
    ! However, adopting the CEH recommendation that we just use an Rns 
    ! (non-stomatal)
    ! independant of LAI or SAI, we need to somehow figure out what to do for
    ! gases in-between O3 and SO2. Instead of applying Wesely's Hstar, 
    ! f0 to Rext and
    ! Rgs separately, I suggest we apply it to the total Rns instead:

           !rv1.3b Gext  = 1.0e-5*Hstar/RextS + f0/RextO
           !dep1_7   Gext  = 1.0e-5*Hstar/wRextS + f0/RextO
  



         ! ##############   4. Calculate Rsur for canopies   ###############

!hf
           if ( DRYDEP_CALC(icmp) == WES_NH3 ) then
           !/** r_water   RIS eq. (24)

             if (Ts_C >0 ) then    ! Use "rh" - now in fraction 0..1.0

               r_water =10.0 * log10(Ts_C+2.0) * exp(100.0*(1.0-rh)/7.0)

               r_water = min( 200.0, r_water)  ! After discussion with Ron
               r_water = max(  10.0,r_water)

             else if ( Ts_C > -5 ) then

               r_water=200.0
             else

               r_water=1000.0
             end if !Ts_C

             Gns_dry = 1.0/r_water
             Gns_wet =  Gns_dry
!hf
           else  ! Not NH3 or NO2:

               Gns_dry = 1.0e-5*Hstar*GnsS_dry + f0 * GnsO   ! OLD SO2!
               Gns_wet = 1.0e-5*Hstar*GnsS_wet + f0 * GnsO 

               !.. and allow for partially wet surfaces at high RH
               Gns_dry = Gns_dry * (1.0-cehfac) + &
                         Gns_wet * cehfac
           end if  ! NH3 test

   !dep1_7 Rsur(icmp) = 1.0/( LAI*DRx(iwes) *g_sto + SAI*Gext + Ggs )


           Rsur(icmp)     = 1.0/( LAI*DRx(iwes) *g_sto + Gns_dry  )
           Rsur_wet(icmp) = 1.0/( LAI*DRx(iwes) *g_sto + Gns_wet  )


           if ( Rsur(icmp) < 1.0 ) then
                print *, "LOWRSUR", icmp, iwes, Rsur(icmp), lu
                print *, "LOWRSUR bits", LAI,DRx(iwes),g_sto, SAI,Gns_dry
                print *, "LOWRSUR H,f0", Hstar,f0
                print *, "LOWRSUR Rx",  RextS, RextO    
                print *, "LOWRSUR h,ustar,cosz",  hveg, ustar, coszen
           end if
              
       else   ! Non-Canopy modelling:

           Rgs     = 1.0/(1.0e-5*Hstar/xRgsS + f0/xRgsO)      ! Eqn. (9)
           Rgs_wet = 1.0/(1.0e-5*Hstar/xRgsS_wet + f0/xRgsO)  ! Eqn. (9)
           Rgs = min(Rgs,9999.9)  

           Rsur(icmp)     = Rgs
           Rsur_wet(icmp) = Rgs
           ! Ggs = 1.0/ Rgs

       end if  ! end of canopy tests 

       !/ --- for NO2 we increase Rsur with a factor of 2. This is a crude
       ! acknowledgement of the fact that NO2 should really be dealt with
       ! as a compensation-point problem, and measured fluxes are often
       ! from the ground (-ve Vg).

       !if ( DRYDEP_CALC(icmp) == WES_NO2 ) then
       !    Rsur(icmp)     = Rsur(icmp) * 2.0
       !    Rsur_wet(icmp) = Rsur_wet(icmp) * 2.0
       !end if


  end do GASLOOP2

 end subroutine Rsurface

!=======================================================================

    subroutine g_stomatal(debug_flag,lu, coszen,Idrctt,Idfuse, &
                         Ts_C,pres,LAI,vpd,SWP,g_sto)
!=======================================================================

!    Calculates stomatal conductance g_sto based upon methodology from 
!    EMEP MSC-W Note 6/00:
!
!    Gsto = [g_max * g_pot * g_light * g_temp * g_vpd * g_smd ]/41000.0
!
! Inputs:

  logical, intent(in) :: debug_flag   ! set true if output wanted 
  integer, intent(in) :: lu           ! land-use index (max = nlu)
  real, intent(in) :: coszen          ! cos of zenith angle
  real, intent(in) ::  Idfuse         ! Diffuse radn, u7.lu
  real, intent(in) ::  Idrctt         ! Direct  radn, u7.lu
  real, intent(in) :: Ts_C            ! surface temperature at height 
                                      ! 2m (deg. C)
  real, intent(in) :: pres            ! surface  pressure (Pa) 
  real, intent(in) :: LAI             ! leaf area index   (m2/m2)
  real, intent(in) :: vpd             ! vapour pressure deficit (kPa)
  real, intent(in) :: SWP             ! soil water potential (MPa)
!d1.4  integer, intent(in) :: SGS, EGS     ! start, end of growing season (day no.)

! Outputs:
 
  real, intent(out) :: g_sto         ! stomatal conductance


        
!..1 ) Calculate g_pot. Max value is 1.0.
!---------------------------------------
!u7.lu - these calculations only needed once per day - moved alongside
!        LAI calculations
!        Subroutine Conif_gpot still kept in this modukle though.
!... 


!..2 ) Calculate g_light 
!---------------------------------------
 ! (n.b. subroutine get_glight is defined below)

  call get_glight(coszen,Idrctt,Idfuse,g_lightfac(lu),LAI,albedo(lu),g_light)    
  

!..3) Calculate  g_temp
!---------------------------------------
  
  g_temp = (Ts_C - g_temp_opt(lu)) / &
           (g_temp_opt(lu) - g_temp_min(lu))
  g_temp = 1.0 - (g_temp*g_temp)
  g_temp = max(g_temp,g_min(lu) )


!..4) Calculate g_vpd
!---------------------------------------

 g_vpd = g_min(lu) + (1.0-g_min(lu)) * (VPD_min(lu)-vpd )/ &
                                    (VPD_min(lu)-VPD_max(lu) )
 g_vpd = min(g_vpd, 1.0)
 g_vpd = max(g_vpd, g_min(lu))


!..5) Calculate g_swp
!---------------------------------------

  !/  Use SWP_Mpa to get g_swp. We just need this updated
  !   once per day, but for simplicity we do it every time-step.

       g_swp = g_min(lu)+(1-g_min(lu))*(PWP(lu)-SWP)/(PWP(lu)-SWP_max(lu))
       g_swp = min(1.0,g_swp)
       !u7.lu g_swp = max(g_min(lu),g_swp)


!.. And finally,
!---------------------------------------
! (using factor 41000 from mmol O3/m2/s to s/m given in Jones, App. 3
!  for 20 deg.C )

   g_sto = (g_light * g_temp * g_vpd * g_swp )
   g_sto = max( g_sto,g_min(lu) )
   g_sto = ( g_max(lu) * g_pot * g_sto )/41000.0

  !if ( g_sto < 0 ) then
  !   print *, "GSTO NEG" , jday, g_sto, g_pot, g_light, g_temp, g_vpd
  ! end if


  end subroutine g_stomatal


! =====================================================================
!rv1.2     subroutine Conif_gpot(imm,g_pot)
!rv1.2 ! =====================================================================
!rv1.2 !   modifies g_pot (g_age) for effect of older needles, with the simple
!rv1.2 !   assumption that g_age(old) = 0.5.
!rv1.2 !
!rv1.2    !/ arguments
!rv1.2 
!rv1.2     integer, intent(in) :: imm    ! month
!rv1.2     real,   intent(inout) :: g_pot   ! Requires initial input of g_pot 
!rv1.2                                      ! (once obtained as output from g_stomatal)
!rv1.2 
!rv1.2    !/ Some parameters:
!rv1.2    !  Proportion of needles which are from current year:
!rv1.2     real, parameter, dimension(12) :: Pc = (/  &
!rv1.2                    0.53, 0.535, 0.54, 0.545, 0.1, 0.15,  &
!rv1.2                    0.27,  0.36, 0.42,  0.48, 0.5,  0.5  /)
!rv1.2 
!rv1.2     real, parameter :: G_POTOLD = 0.5  ! value of g_pot for old needles
!rv1.2 
!rv1.2 
!rv1.2 
!rv1.2 !needles from current year assumed to have g_pot as evaluated above;
!rv1.2 !needles from previous years assumed to have g_pot of 0.5
!rv1.2 !The sum of the g_pot's for the current year is added to the sum of the
!rv1.2 !g_pot's for previous years to obtain the overall g_pot for the landuse 
!rv1.2 !category temp. conif. forests.
!rv1.2 
!rv1.2     g_pot = Pc(imm)*g_pot + (1.0-Pc(imm))*G_POTOLD
!rv1.2 
!rv1.2   end subroutine Conif_gpot
! =====================================================================

!===========================================================================
    subroutine get_glight(coszen,Idrctt,Idfuse,g_lightfac,LAI,albedo,g_light)
!===========================================================================
!
!    Calculates g_light, using methodology as described in Emberson et 
!    al. (1998), eqns. 31-35, based upon sun/shade method of  
!    Norman (1979,1982)

!     input arguments:

    real, intent(in) :: coszen    ! cos(zen), zen=zenith angle
    real, intent(in) :: Idrctt    ! total direct solar radiation (W/m^2)
    real, intent(in) :: Idfuse    ! diffuse solar radiation (W/m^2)
    real, intent(in) :: g_lightfac ! land-use specific light factor
    real, intent(in) :: LAI       ! leaf area index (m^2/m^2), vegetation        
                                  ! specific
    
    real, intent(in) :: albedo    ! Fraction, 0-1

!     output arguments
  
    real, intent(out) :: g_light    

!   Some parameters

    real, parameter :: &
      PARfrac = 0.45,   & ! approximation to fraction (0.45 to 0.5) of total 
                          ! radiation in PAR waveband (400-700nm)
      Wm2_uE  = 4.57,   &  ! converts from W/m^2 to umol/m^2/s
      cosA    = 0.5       ! A = mean leaf inclination (60 deg.), where it is 
                          ! assumed that leaf inclination has
                          ! a spherical distribution

!     internal variables:

    
    real :: sinB      ! B = solar elevation angle  = complement of zenith angle
    
    real :: LAIsun    ! sun-fraction of LAI
    real :: LAIshade  ! shade-fraction of LAI
    real :: PARsun    ! sun-fraction of PAR
    real :: PARshade  ! shade-fraction of PAR
       
    real :: g_sun    ! sun-fraction contribution to g_light
    real :: g_shade  ! shade-fraction contribution to g_light


    !DEP1.2 error sinB = 1.0/coszen  ! = cos(zen)
    sinB = coszen  ! = cos(zen)


    LAIsun = (1.0 - exp(-0.5*LAI/sinB) ) * sinB/cosA
    
    LAIshade = LAI - LAIsun

! PAR flux densities evaluated using method of
! Norman (1982, p.79): 
! "conceptually, 0.07 represents a scattering coefficient"  

    PARshade = Idfuse * exp(-0.5*LAI**0.7) +  &
               0.07 * Idrctt  * (1.1-0.1*LAI)*exp(-sinB)   

    PARsun = Idrctt *cosA/sinB + PARshade

!.. Convert units, and to PAR fraction
!.. and multiply by albedo

    PARshade = PARshade * PARfrac * Wm2_uE * ( 1.0 - albedo )
    PARsun   = PARsun   * PARfrac * Wm2_uE * ( 1.0 - albedo )


    g_sun   = (1.0 - exp (-g_lightfac*PARsun  ) ) * (LAIsun  /LAI)
    g_shade = (1.0 - exp (-g_lightfac*PARshade) ) * (LAIshade/LAI)

    g_light = g_sun + g_shade

  end subroutine get_glight

!--------------------------------------------------------------------

end module Rsurface_ml

!COMMENTS
!In his treatment of surface resistance, Wesley does not recommend that NO_2 
! and HNO_3 be treated similarly to O_3 (see ibid. p.1296).

!The reason for deviating from Wesley here (with respect to NO_2) is that 
!whilst Wesley assumes that all NO_2 is deposited through stomata, it would 
!seem more likely that deposition on other surfaces of plants should also be 
!taken into consideration. 
