module DepVariables_ml

!DESCRIPTION
! Module containing main variables and parameters associated with the
! deposition module based upon Emberson et al. (EMEP Rep. 6/01) and
! Wesely et al. 

! So far, 17 land-use classes are defined. The classes temperate orchard and 
! Mediterranean orchard may be added later. 

implicit none
private

integer, public, parameter :: NLANDUSE = 17 !number of land-use classes
integer, public, parameter :: LU_CONIF =  1 ! Conif. forest (temp.)
integer, public, parameter :: LU_WATER =  15 !  Sea/lakes

integer, public, save :: lu            !land-use class index

! g_sto factors in UK dep method:

real, public, save :: &
        g_pot      &                  ! stomatal conductance age factor 
       ,g_temp     &                  ! stomatal conductance temperature factor
       ,g_vpd      &                  ! stomatal conductance vpd factor
       ,g_light    &                  ! stomatal conductance light factor
       ,g_swp                         ! stomatal conductance soil water factor


! In the Wesely method of dealing with many gases, the external resistance 
! equations for ozone and SO2 are used as reference points for the remaining
! gases, e.g. for  NH3, HNO3 and NO2.

    real, public, parameter :: &   ! external resistance for Ozone, Sulphur
          RextO =  2500.0   & ! prelim. value - gives Gext=0.2 cm/s for LAI=5
         ,RextS = 10000.0     ! prelim. value - high

! Here, "Gext=0.2cm/s" refers to the external conductance, G_ext, where 
! G_ext=LAI/R_ext. In many studies, it has been assumed 
! that G_ext should be low, particularly relative to stomatal conductance g_s.
! Results from a variety of experiments, however, have made the above 
! estimates  Rext0 and RextS plausible.  The above equation for G_ext has been
! designed on the basis of these experimental results. 

! Notice also that given the equations for the canopy resistance R_sur and the 
! deposition velocity V_g, V_g>=LAI/R_ext. The value of G_ext can therefore be
! interpreted as the minimum value for V_g.

!.. names
   character(len=20), public, dimension(NLANDUSE), save :: luname   

!.. common variables read from ukdep_biomass.dat   ...........................
! Note:
! SGS(Lat) = SGS50 + (Lat-50)*DSG
! The above equation is introduced as a means
! of estimating SGS(Lat), where "Lat" denotes latitude.
! SGS50 and DSGS are obtained from a curve relating 
! latitude to SGS.
! Similar remarks apply for EGS.

   real, public, save, dimension(NLANDUSE) :: &
          hveg_max    & ! max height of vegetation
       ,  b_inc       & ! RIVM factor for Rinc. If b_inc =0 no canopy mod'ng
       ,  albedo      & ! Albedo, 0-1
       ,  NH4_pl      & ! **sc: landuse dependent variable ranging over 
                        ! intercellular ammonium concentrations
       ,  SGS50       & ! Start of growing season (days) at 50 deg. N 
       ,  DSGS        & ! d(SGS)/d(Lat). 
       ,  EGS50       & ! end of growing season
       ,  DEGS        & ! d(EGS)/d(Lat)
       ,  LAImin      & ! min. value of LAI
       ,  LAImax      & ! max. value LAI
       ,  SLAIlen     & ! days from LAImin to LAImax at start of season
       ,  ELAIlen       ! days from LAImax to LAImin at end of season

!.. common variables read from ukdep_gfac1.dat   ...........................

   real, public, save, dimension(NLANDUSE) :: &
          g_pot_min         & ! min. value of g_pot
       ,  Sg_potlen        & ! length in days from g_pot=g_min to g_pot=1
       ,  Eg_potlen        & ! length in days from g_pot=1 to g_pot=g_potmin
       ,  g_max            & ! max. value conductance g_s
       ,  g_min            & ! min. value Conductance g_s
       ,  g_lightfac       & ! light coefficient 
       ,  g_temp_min       & ! temperature when g_temp starts
       ,  g_temp_opt       & ! temperature when g_temp max.   
       ,  g_temp_max         ! temperature when g_temp stops  

!.. .......and from ukdep_gfac2.dat   ........................................

   real, public, save, dimension(NLANDUSE) :: &
          RgsS           & ! ground surface resistance, Sulphur
       ,  RgsO           & ! ground surface resistance, Ozone   
       ,  VPD_max        & ! threshold VPD when relative g = g_min
       ,  VPD_min        & ! threshold VPD when relative g = 1
       ,  SWP_max        & ! threshold SWP when relative g = 1
       ,  PWP            & ! threshold SWP when relative g = g_min
                           ! and assumed equal to permanent wilting point
       ,  rootdepth        ! root depth (mm)
!..............................................................................
end module DepVariables_ml
