module Rsurface_ml

!================== Now under CVS control =================
! $Author: mifads $
! $Id: Rsurface_ml.f90,v 1.2 2002-09-21 14:45:04 mifads Exp $
! $Name: not supported by cvs2svn $
! =========================================================

!Changes from Version 1.0
!d1.1 - zenith angle range now to 89.9, not 85.0
!d1.3 - sinB now set to coszen, not 1/coszen
!d1.4 - restructered for use in full EMEP model.

use DepVariables_ml, only : g_pot, g_temp,g_vpd,g_light,g_swp,    & !u7.4
                              b_inc     , albedo    ,   &
                              SAIadd, &
                              g_max     , g_min     , g_lightfac,   &
                              g_temp_min, g_temp_opt, &
                              RgsS      , RgsO      , RextS, RextO, &
                              Rgs, Rinc, Gext, &
                              VPD_max   , VPD_min   , &
                              SWP_max   , PWP       , rootdepth 
!d1.4 use Io_ml,        only : IO_UKDEP, open_file
use My_UKDep_ml, only : NDRYDEP_CALC, & ! no. of Wesely species  used
                         DRYDEP_CALC!,  & ! array with species which are needed
                         !GASNAME    !array with names of gases considered

                    
use PhysicalConstants_ml, only : KARMAN
use SoilWater_ml, only : SWP
!d1.5.2use UKsetup_ml, only :   canopy, leafy_canopy, &   ! logical characteristics
!d1.5.2                         SAIadd          
use Wesely_ml, only  : Wesely_tab2, & ! Wesely Table 2 for 14 gases
                       WES_HNO3, & ! indices to identify HNO3
                       Rb_cor,  &! correction factor used in evaluating Rb
                       DRx       !  Ratio of diffusivities to ozone
implicit none
private

!u7.lu removed:
!u7.lu use Metdata_ml,   only : snow, psurf, t2, METSTEP, pr 
!u7.lu - pass in met params to preserve consistency with box-model:
!u7.lu use Met_ml,   only : snow, psurf, t2, pr 
!u7.lu use Radiation_ml, only : zen, coszen, Idfuse, Idrctt,SolBio

public   ::  Rsurface
public   ::  Conif_gpot
private  ::  g_stomatal    
private  ::  get_glight

logical, private, parameter :: DEBUG_RSURF = .false.
 
    
contains
! =======================================================================

  subroutine Rsurface(lu,debug_flag, LAI,hveg,&
                      z0,ustar,Ts_C,vpd,SWP, &
                      psurf, pr, &                    !u7.lu
                      coszen, Idfuse, Idrctt, &       !u7.lu
                      snow, &                    !u7.lu
                        g_sto,Rsur,Rb)
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
    real, intent(in) ::  pr
    real, intent(in) ::  coszen
    real, intent(in) ::  Idfuse
    real, intent(in) ::  Idrctt
    integer, intent(in) :: snow         ! snow=1, non-snow=0

! Output:

   real, intent(out)             :: g_sto !  Stomatal conducatance (s/m)
   real,dimension(:),intent(out) :: Rsur ! bulk canopy surface resistance (s/m)
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
    real :: Ggs                 ! ground surface conductance, any gas
!rv1.5.2   real :: Gext                ! external conductance
!CORR    real :: diffc               ! molecular gas diffusivity coefficient 
!CORR                                ! (extracted from Wesely_tab2)
  ! CORR: use D_H2O, D_i instead of diffc !
    real, parameter :: D_H2O = 0.21e-4  ! Diffusivity of H2O, m2/s
                                      ! From EMEP Notes
    real            :: D_i              ! Diffusivity of gas species, m2/s


! START OF PROGRAMME: 


  !/** Do Rb first:
  !===========================================================================
 
  GASLOOP1: do icmp = 1, NDRYDEP_CALC
      iwes = DRYDEP_CALC(icmp)          ! index in Wesely table
                  
      if   ( hveg  >=  0.0 ) then

          Rb(icmp) = 2.0 * Rb_cor(iwes)/(KARMAN*ustar)

      if ( DEBUG_RSURF .and. Rb(icmp) < 0.0  ) then
        print *, "ERROR RSURFB", lu, icmp, iwes, Rb_cor(iwes), Rb(icmp)
        return          
      end if

      else ! water, sea, rb is calculated as per Hicks and Liss (1976)

          ! CORR - old: ! 
          ! CORR diffc = Wesely_tab2(1,iwes) ! molecular diffusivity coeff.

          ! CORR - Wesely's table gives the ratio D_H2O/D_i
          ! CORR - therefore D_i = D_H2O / Wesely

          D_i = D_H2O / Wesely_tab2(1,iwes)  ! CORR !

!CORR 19/8/2002 - spotted by pw. Corrected diffc to D_i by ds...
! should double-check equations though.
          !CORR Rb(icmp) = log( z0 * KARMAN * ustar/ diffc )

          Rb(icmp) = log( z0 * KARMAN * ustar/ D_i )
          Rb(icmp) = Rb(icmp)/(ustar*KARMAN)

         ! CORR - Rb can be very large or even negative from this
         !        formula. We impose limits:

          Rb(icmp) = min( 1000.0, Rb(icmp) )    ! CORR - - gives min 0.001 m/s!
          Rb(icmp) = max(   10.0, Rb(icmp) )    ! CORR - - gives max 0.10 m/s!

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
    !u7.lu zen > 1.0e-15 .and. zen<=89.9 ) then  ! Daytime 

         call g_stomatal(debug_flag,lu,coszen,Idrctt,Idfuse, & 
                         Ts_C,psurf,LAI,vpd,SWP,g_sto)
   !d1.6 else
   !d1.6     g_sto = 0.0    ! Dark, no stomatal conductance...

   end if ! leafy canopy and daytime


  !/** Calculate Rinc, Gext 

   if(  canopy ) then   
         SAI  = LAI + SAIadd(lu)                  ! Accounts for bark/twigs
         Rinc = b_inc(lu)* SAI * hveg  / ustar    ! Lisa's eqn. 48 for u*=0.5
         Rinc = min( Rinc, 1000.0)

   else   ! No canopy present

        Rinc = -999.9
        Gext = -999.9  !  (not needed maybe, but...)
      
   end if !  canopy


!####   2. Calculate Surface Resistance, Rsur, for HNO3 and  
!####      Ground Surface Resistance, Rgs, for the remaining 
!####      Gases of Interest                                

!In his treatment of surface resistance, Wesley does not recommend that NO_2 
! and HNO_3 be treated similarly to O_3 (see ibid. p.1296).

!The reason for deviating from Wesley here (with respect to NO_2) is that 
!whilst Wesley assumes that all NO_2 is deposited through stomata, it would 
!seem more likely that deposition on other surfaces of plants should also be 
!taken into consideration. 

 !u7.lu - previous if-test for snow replaced by more efficient
 !        multiplication, since snow=0 or 1.

  xRgsS = RgsS(lu) + Rlow  + snow * 2000.0    !u7.lu
  xRgsO = RgsO(lu) + Rlow  + snow *    0.0    !u7.lu QUERY???


!.........  Loop over all required gases   ................................

  GASLOOP2: do icmp = 1, NDRYDEP_CALC

     !-------------------------------------------------------------------------

     !  code obtained from Wesely during 1994 personal communication

        if ( DRYDEP_CALC(icmp) == WES_HNO3 ) then
            Rsur(icmp)  = max(10.0,Rlow)
            cycle GASLOOP2
        end if

     !-------------------------------------------------------------------------
     ! Calculate the Wesely variables Hstar (solubility) and f0 (reactivity)

        iwes = DRYDEP_CALC(icmp)
        Hstar =Wesely_tab2(2,iwes)    !Extract H*'s 
        f0    =Wesely_tab2(5,iwes)    !Extract f0's
    
     !-------------------------------------------------------------------------

        Rgs = 1.0/(1.0e-5*Hstar/xRgsS + f0/xRgsO)      ! Eqn. (9)
        Rgs = min(Rgs,9999.9)  
                          

     !   Use SAI to test for snow, ice, water, urban ...

       if ( canopy  ) then   

         ! ###   3. Calculate Cuticle conductance, Gext   ################
         ! ###      and  Ground surface conductance Ggs:

         ! Corrected for other species
         ! using Wesely's eqn. 7 approach. (We identify leaf surface resistance
         ! with Rext/SAI.)

           Gext  = 1.0e-5*Hstar/RextS + f0/RextO
  
           Ggs = 1.0/( Rgs + Rinc )


         ! ##############   4. Calculate Rsur for canopies   ###############

           Rsur(icmp) = 1.0/( LAI*DRx(iwes) *g_sto + SAI*Gext + Ggs )

           if ( Rsur(icmp) < 1.0 ) then
             print *, "LOWRSUR", icmp, iwes, Rsur(icmp), lu
             print *, "LOWRSUR bits", LAI*DRx(iwes) *g_sto, SAI*Gext, Ggs
             print *, "LOWRSUR bit2", LAI,DRx(iwes),g_sto, SAI
             print *, "LOWRSUR H,f0", Hstar,f0
             print *, "LOWRSUR Rx",  RextS, RextO    
             print *, "LOWRSUR h,ustar,cosz",  hveg, ustar, coszen
             print *, "LOWRSUR gs",  Rgs, Rinc, RgsS(lu), Rlow, xRgsS
           end if
              
       else   ! Non-Canopy modelling:

           Rsur(icmp) = Rgs
           Ggs = 1.0/ Rgs

       end if  ! end of canopy tests 


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
    subroutine Conif_gpot(imm,g_pot)
! =====================================================================
!   modifies g_pot (g_age) for effect of older needles, with the simple
!   assumption that g_age(old) = 0.5.
!
   !/ arguments

    integer, intent(in) :: imm    ! month
    real,   intent(inout) :: g_pot   ! Requires initial input of g_pot 
                                     ! (once obtained as output from g_stomatal)

   !/ Some parameters:
   !  Proportion of needles which are from current year:
    real, parameter, dimension(12) :: Pc = (/  &
                   0.53, 0.535, 0.54, 0.545, 0.1, 0.15,  &
                   0.27,  0.36, 0.42,  0.48, 0.5,  0.5  /)

    real, parameter :: G_POTOLD = 0.5  ! value of g_pot for old needles



!needles from current year assumed to have g_pot as evaluated above;
!needles from previous years assumed to have g_pot of 0.5
!The sum of the g_pot's for the current year is added to the sum of the
!g_pot's for previous years to obtain the overall g_pot for the landuse 
!category temp. conif. forests.

    g_pot = Pc(imm)*g_pot + (1.0-Pc(imm))*G_POTOLD

  end subroutine Conif_gpot
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
