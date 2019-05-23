module  My_Outputs_mod
! -----------------------------------------------------------------------
! Allows user to specify which species are output to various
! ascii and binary output files.
!
! Sites  - surface sites,     to sites.out
! Sondes - vertical profiles, to sondes.out
! Hourly - ascii output of selected species, selcted domain
! -----------------------------------------------------------------------

use CheckStop_mod,      only: CheckStop
use ChemDims_mod,       only: NSPEC_ADV, NSPEC_SHL
use ChemSpecs_mod
use Debug_module,       only: DEBUG   ! -> DEBUG%COLSRC and POLLEN
use SmallUtils_mod,     only: find_duplicates
use TimeDate_mod,       only: date
use Units_mod,          only: Init_Units

implicit none

logical, public, parameter :: out_binary = .false.
logical, public, parameter :: Ascii3D_WANTED = .false.

! Site outputs   (used in Sites_mod)
!==============================================================
! Specify the species to be output to the sites.out file
! For met params we have no simple index, so we use characters.
! These must be defined in Sites_mod.f90.

integer, private :: isite              ! To assign arrays, if needed
integer, public, parameter :: &
   NSITES_MAX =        99     & ! Max. no surface sites allowed
  ,FREQ_SITE  =         1     & ! Interval (hrs) between outputs
  ,NADV_SITE  = NSPEC_ADV     & ! No. advected species (1 up to NSPEC_ADV)
  ,NSHL_SITE  =        1      & ! Bosco OH NSPEC_SHL     & ! No. short-lived species
  ,NXTRA_SITE_MISC =    2     & ! No. Misc. met. params  ( e.g. T2, d_2d)
  ,NXTRA_SITE_D2D  =   20       ! Bosco = +5-4 No. Misc. met. params  ( e.g. T2, d_2d)
!Bosco  ,NXTRA_SITE_D2D  =  9+8       ! No. Misc. met. params  ( e.g. T2, d_2d)

integer, public, parameter, dimension(NADV_SITE) :: &
  SITE_ADV = [(isite, isite=1,NADV_SITE)]  ! Everything

integer, public, parameter, dimension(NSHL_SITE) :: &
  SITE_SHL = [ IXSHL_OH ]  ! All short-lived species
  !SITE_SHL = [(isite, isite=1,NSHL_SITE)]  ! All short-lived species

! Extra parameters - need to be coded in Sites_mod also. So far
! we can choose from hmix, T2, or th (pot. temp.) or d_2d fields.
!  d_2d fields can be accessed from Derived_mod by setting common index
!  "D2D" in SITE_XTRA and the actual field name (as defined in Derived_mod)
!  in SITE_XTRA_CODE (e.g. "D2_PM25 " or "D2_SIA") :

!** IMPORTANT!! Make sure the correspondence between selected for output
!** fields in SITE_XTRA and their names in SITE_XTRA_CODE

character(len=18), public, parameter, dimension(NXTRA_SITE_MISC) :: &
  SITE_XTRA_MISC=[character(len=18):: "th","T2"]

!integer, parameter :: MAX_NEXTRA_SITED2D=100
!character(len=24), public, save, dimension(MAX_NEXTRA_SITED2D) :: &
!   site_outputs_extraD2D = '-', sonde_outputs_extraD2D = '-'

!These variables must have been set in My_Derived for them to be used.
character(len=24), public, parameter, dimension(NXTRA_SITE_D2D) :: &
  SITE_XTRA_D2D=[character(len=24):: &
    "HMIX","PSURF", & ! Bosco skip: "ws_10m","rh2m",&
    "Emis_mgm2_BioNatC5H8","Emis_mgm2_BioNatBIOTERP",&
    "Emis_mgm2_BioNatNO","Emis_mgm2_nox",&
    'WDEP_PREC',&!''SNratio',&
    'met2d_uref','met2d_u10', 'met2d_v10','met2d_rh2m', &
    !'met2d_SMI1', 'met2d_SMI3',&
    'met2d_SMI_uppr', 'met2d_SMI_deep',&
    'met2d_ustar_nwp', 'met2d_LH_Wm2', 'met2d_SH_Wm2',&
    !BB 'SMI_deep','met2d_SMI_d','SMI_uppr','met2d_SMI_s',&
!Boscso Extra: +5
    !BB'USTAR_NWP', 
    'USTAR_DF','INVL_DF', &
    'met2d_PARdbh', 'met2d_PARdif' &
]

!End Boscso Extra:
  !OLD 'Idirect','Idiffuse' &
!   "SoilWater_deep","EVAP_CF","EVAP_DF",&
!   "EVAP_BF","EVAP_NF","WDEP_PREC",&
!   "RH_GR","CanopyO3_GR","VPD_GR","FstO3_GR",&
!   "RH_IAM_DF","CanopyO3_IAM_DF","VPD_IAM_DF","FstO3_IAM_DF",&
!   "COLUMN_CO_k20","COLUMN_C2H6_k20","COLUMN_HCHO_k20","COLUMN_CH4_k20",
!    "COLUMN_NO2_k20"]

!**** Aircraft outputs   (used in Polinat_mod)
!==============================================================
!   Specify the species to be output by Polinat for aircraft flight tracks

integer, public, parameter :: &
  NFLIGHT_MAX =    10         &   ! Max. no sondes allowed
 ,FREQ_FLIGHT =    12         &   ! Interval (hrs) between outputs
 ,NADV_FLIGHT =    1              ! No.  advected species

integer, public, parameter, dimension(NADV_FLIGHT) :: &
  FLIGHT_ADV =  (/ IXADV_O3 /)


!**** Sonde outputs   (used in Sites_mod)
!==============================================================
!     Specify the species to be output to the sondes.out file
!  We typically deal with fewer species for sonde output than
!  surface sites, so we use a different method to specify.
! For met params we have no simple index, so we use characters.
! These must be defined in Sites_mod.f90.

integer, public, parameter :: &
   NSONDES_MAX =    99        &   ! Max. no sondes allowed
  ,NLEVELS_SONDE =  20        &   ! No. k-levels (9 => 0--2500 m)
  ,FREQ_SONDE  =     1        &   ! Interval (hrs) between outputs
  ,NADV_SONDE  =    9         &   ! No.  advected species
  ,NSHL_SONDE  =    1         &   ! No. short-lived species
  ,NXTRA_SONDE =    4             ! No. Misc. met. params

integer, public, parameter, dimension(NADV_SONDE) :: &
  SONDE_ADV = [ IXADV_O3, IXADV_NO2, IXADV_NO, IXADV_PAN,  &
                IXADV_NO3_c, IXADV_NO3_f, IXADV_SO4, IXADV_NH4_f, IXADV_NH3 ]

integer, public, parameter, dimension(NSHL_SONDE) :: &
  SONDE_SHL = [ IXSHL_OH ] ! , IXSHL_OD, IXSHL_OP ]
character(len=10), public, parameter, dimension(NXTRA_SONDE) :: &
  SONDE_XTRA= [character(len=10):: &
   "NOy","z_mid","p_mid","th"]!,"Kz_m2s"]

!   can access d_3d fields through index here, by
!   setting "D3D" above and say D3_XKSIG12 here:

!*** wanted binary dates... specify days for which full binary
!    output is wanted. Replaces the hard-coding which was in wrtchem:
integer, public, parameter :: NBDATES = 3
type(date), public, save, dimension(NBDATES) :: wanted_dates_inst

!================================================================

public :: set_output_defs

contains

subroutine set_output_defs
  implicit none
  character(len=144) :: errmsg   ! Local error message
  character(len=*), parameter :: dtxt = 'My_Out:set:' ! debug txt

  call Init_Units()
  
  ! Safety checks for common mistakes:
  errmsg = find_duplicates(SITE_XTRA_D2D)
  call CheckStop ( errmsg /= 'no', dtxt//' Duplicate SITE_XTRA_D2D'//errmsg)

  !*** Wanted dates for instantaneous values output:
  !    specify months,days,hours for which full output is wanted.
  wanted_dates_inst(1) = date(-1,1,1,0,0)
  wanted_dates_inst(2) = date(-1,1,1,3,0)
  wanted_dates_inst(3) = date(-1,1,1,6,0)

endsubroutine set_output_defs
endmodule My_Outputs_mod
