module DepVariables_ml
  use ModelConstants_ml, only : NLANDUSE

!DESCRIPTION
! Module containing main variables and parameters associated with the
! deposition module based upon Emberson et al. (EMEP Rep. 6/01) and
! Wesely et al. 

! So far, 19 land-use classes are defined. The classes temperate orchard and 
! Mediterranean orchard may be added later. 

implicit none

!DOSE integer, public, parameter :: NLANDUSE =   19 !number of land-use classes
integer, public, parameter :: LU_WATER =   14 !  Sea/lakes
integer, public, parameter :: LU_ICE   =   15 !  Ice  ! Pb210 - temporary
real,    public, parameter :: STUBBLE  = 0.01 ! Veg. ht. out of season
integer, public, save :: lc_water             !  Sea/lakes
integer, public, save :: lc_ice               !  Sea/lakes

!/-- assign UK 17 to categories - used in dry deposition outputs.
! Note that IAM_ stuff are "fakes", so not added in.
! See ~mifads/Unify/Data/Landuse.list for full list

integer, public, parameter,dimension(2) :: ECO_CONIF_FOREST =  (/ 1, 3 /)
integer, public, parameter,dimension(2) :: ECO_DECID_FOREST =  (/ 2, 4 /)
integer, public, parameter,dimension(3) :: ECO_CROP =  (/ 5, 6,7 /) 
! After discussions with Max Posch (CCE), we define seminat as basically
!   everything except awater, forest and crops
!   8=moorland, grass=9, medscrub=10,wetlands=11,tundra=12
integer, public, parameter,dimension(5) :: ECO_SEMINAT =  (/ 8, 9, 10, 11, 12 /)
integer, public, parameter,dimension(1) :: ECO_WETLAND =  (/ 11 /)
integer, public, parameter,dimension(1) :: ECO_WATER   =  (/ LU_WATER /)
integer, public, parameter  :: IAM_WHEAT =  17  !ds, JUN06
integer, public, parameter  :: IAM_BEECH =  18  !ds, JUN06
integer, public, parameter  :: IAM_MEDOAK =  19  !ds, JUN06


! g_sto factors in UK dep method:

real, public, save :: &
        f_phen     &                  ! stomatal conductance age factor 
       ,f_temp     &                  ! stomatal conductance temperature factor
       ,f_vpd      &                  ! stomatal conductance vpd factor
       ,f_light    &                  ! stomatal conductance light factor
       ,f_swp                         ! stomatal conductance soil water factor

 real, public, save  :: Rinc      ! in-canopy resistance (s/m)
 real, public, save  :: Rgs       ! ground surface resistance, any gas
 real, public, save  :: Gext      ! external surface (e.g. cuticular) conductance, any gas
 real, public, save  :: g_sto     ! stomatal conductance


! In the Wesely method of dealing with many gases, the external resistance 
! equations for ozone and SO2 are used as reference points for the remaining
! gases, e.g. for  NH3, HNO3 and NO2.

    real, public, parameter :: &   ! external resistance for Ozone, Sulphur
          RextO =  2500.0   & ! prelim. value - gives Gext=0.2 cm/s for LAI=5
         ,RextS = 10000.0     ! prelim. value - high - for dry conditions
        !,RextS =   200.0     ! more approriate to wet....

! Here, "Gext=0.2cm/s" refers to the external conductance, G_ext, where 
! G_ext=LAI/R_ext. In many studies, it has been assumed 
! that G_ext should be low, particularly relative to stomatal conductance g_s.
! Results from a variety of experiments, however, have made the above 
! estimates of Rext0 and RextS plausible.  The above equation for G_ext has 
! been designed on the basis of these experimental results. 

! Notice also that given the equations for the canopy resistance R_sur and the 
! deposition velocity V_g, V_g>=LAI/R_ext. The value of G_ext can therefore be
! interpreted as the minimum value for V_g.

!.. names
   character(len=20), public, dimension(NLANDUSE), save :: luname   

   logical, public, dimension(NLANDUSE), save :: &
              crops          &! true for veg which grows...
             ,bulk           &! true for land-classes without LAI
             ,water          &! true for  water, set with h < 0
             ,forest         &! Assumed when hveg_max > 5 m
             ,conif_forest   &! Assumed when hveg_max > 5 m and SGS<=1
             ,luflux_wanted  &! Set true if stomatal flux calcs possible
             ,vegetation     &! Assumed when hveg > 5 m and not urban
             ,urban           ! Assumed when hveg > 5 m && LAI < 0.0

  real, public,  dimension(NLANDUSE), save :: &
       SAIadd                  ! Additional surface area for bark, twigs

  real, public,  dimension(NLANDUSE), save :: &
      lai_flux,         & ! Fluxes to total LAI
      unit_flux,        & ! Fluxes per m2 of leaf area
      leaf_flux           ! Flag-leaf Fluxes per m2 of leaf area


!.. commen variables read from ukdep_biomass.dat   ...........................
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
       ,  NH4_pl      & ! landuse dependent variable ranging over 
                        ! intercellular ammonium concentrations
       ,  SGS50       & ! Start of growing season (days) at 50 deg. N 
       ,  DSGS        & ! d(SGS)/d(Lat). 
       ,  EGS50       & ! end of growing season
       ,  DEGS        & ! d(EGS)/d(Lat)
       ,  LAImin      & ! min. value of LAI
       ,  LAImax        ! max. value LAI

!.. common variables read from JUN06_gfac1.dat   ...........................

   integer, public, save, dimension(NLANDUSE) :: &
     !JUN06     SGS              & ! Start of growing season (days)
     !JUN06  ,  EGS              & ! End   of growing season (days)
         SLAIlen          & ! days from LAImin to LAImax at start of season
       , ELAIlen            ! days from LAImax to LAImin at end of season

   real, public, save, dimension(NLANDUSE) :: &
          f_phen_a         & ! min. value of f_phen
         ,f_phen_b         & ! min. value of f_phen
         ,f_phen_c         & ! min. value of f_phen
         ,f_phen_d         & ! min. value of f_phen
         ,f_phen_Slen      & ! Length of Startup  (days)
         ,f_phen_Elen      & ! Length of End period (days)
       ,  g_max            & ! max. value conductance g_s
       ,  f_min            & ! min. value Conductance g_s
       ,  f_lightfac       & ! light coefficient 
       ,  f_temp_min       & ! temperature when g_temp starts
       ,  f_temp_opt       & ! temperature when g_temp max.   
       ,  f_temp_max         ! temperature when g_temp stops  

!.. .......and from JUN06_gfac2.dat   ........................................

   real, public, save, dimension(NLANDUSE) :: &
          RgsS           & ! ground surface resistance, Sulphur
       ,  RgsO           & ! ground surface resistance, Ozone   
       ,  VPD_max        & ! threshold VPD when relative f = f_min
       ,  VPD_min        & ! threshold VPD when relative f = 1
       ,  SWP_max        & ! threshold SWP when relative f = 1
       ,  PWP            & ! threshold SWP when relative f = f_min
                           ! and assumed equal to permanent wilting point
       ,  rootdepth        ! root depth (mm)

!..............................................................................
end module DepVariables_ml
