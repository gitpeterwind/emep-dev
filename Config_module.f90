module Config_module
!----------------------------------------------------------------------------
! Specifies a number of constants used in the model, and reads namelist
! file to (re-)configure where possible.
! Note that physical constants (e.g. gravity, Cp, etc ( are specified in
! the module PhysicalConstants_mod.f90)
!----------------------------------------------------------------------------
use AeroConstants_mod,     only: AERO
use CheckStop_mod,         only: CheckStop
use ChemDims_mod,          only: NSPEC_ADV, NSPEC_SHL
use ChemSpecs_mod,         only: species, CM_schemes_ChemSpecs
use ChemGroups_mod,        only: chemgroups
use Debug_module,          only: DEBUG, DebugCell
use EmisDef_mod,           only: Emis_heights_sec_MAX, Emis_Nlevel_MAX, Emis_h, Emis_Plevels, Emis_h_pre
use Io_Nums_mod,           only: IO_NML, IO_LOG, IO_TMP
use OwnDataTypes_mod,      only: typ_ss, lf_sources, lf_country_group_type, uEMEP_type, Emis_id_type, &
                                 emis_in, EmisFile_id_type, Emis_sourceFile_id_type,&
                                 Sector_type,&
                                 TXTLEN_NAME, TXTLEN_FILE, TXTLEN_SHORT,&
                                 TXTLEN_DERIV, Emis_mask_type, &
                                 Deriv, typ_s1ind,typ_s5ind,O3cl_t,typ_s3,typ_s4
use TimeDate_mod,          only: date
use Precision_mod,         only: dp
use SmallUtils_mod,        only: find_index, key2str

implicit none
private

public :: Config_Constants
public :: WriteConfig_to_RunLog

!=============================================================================
! Experiment name:
!  EMEPSTD      Standard run & output
!  EMEP2010    EMEPSTD with Iceland Volcanic Eruption input
!  TFMM        EMEPSTD, but with INERIS_SNAP & TFMM hourly output
!  FORECAST    Forecast run, MACC-ENS hourly output & BC
!  EVA2010     FORECAST with MACC-EVA2010 hourly output & BC
!  EMERGENCY   FORECAST with ONLY Volcanic Eruption & Nuclear Accident.

!integer, public, parameter :: &
!  TXTLEN_NAME =  64, & !for performance, should be a multiple of 8
!  TXTLEN_FILE = 200    ! large enough for paths from namelists

CHARACTER(LEN=TXTLEN_NAME), public, save :: EXP_NAME="EMEPSTD"

! It is quite easy to end config files accidently through e.g. an extra
! & from a fortran comment. We place a marker at the end of the (first-called)
! config file, and test for this.
CHARACTER(LEN=TXTLEN_NAME), private, save :: LAST_CONFIG_LINE="NOTSET"
CHARACTER(LEN=TXTLEN_NAME), private, save :: LAST_CONFIG_LINE_DEFAULT

! EMEP daily measurements end at 6am, hence we typically adjust
! for that. For global though, zero would be more normal
  integer, save, public :: END_OF_EMEPDAY = 6 ! 

  type, private :: PBL_t
    ! Zi minimum value now generally calculated as z_mid(19), but we
    !   keep a fixed value for smoothing.  (old comment?)
    real :: ZiMIN = 50.0                     ! minimum mixing height
    real :: ZiMAX = 3000.0                   ! maximum mixing height
    character(len=10) :: HmixMethod = "JcRb_t2m"  ! Method used for Hmix 
      ! JcRb = Jericevic/Richardson number method
      ! JcRb_surfT is the new JcRb using pop T at skin
      ! JcRb_t2m is new JcRb using pop T at 2m
      ! "SbRb"= Seibert !"TIZi" = Original from Trond Iversen tiphysics
    real :: MIN_USTAR_LAND = 0.1 ! m/s - Defines stable BL height
    !
    ! Moved Feb 2021 from BLPhysics:
    logical :: NWP_Kz=.false.  ! hb 23.02.2010 Kz from meteo. NOT WORKING!
    logical :: USE_MIN_KZ =.false. ! "fix"
    character(len=9):: &
      KzMethod = "Mixed"  ! Set U, S separately, default, pre rv4-38
                      !   (defaulted to OBrien U + Jericevic S)
                      ! "TRONKz"   ! TROEN - U+S New default in testing :-)
                      ! "SILAMKz"  ! SILAM - U+S  
                      ! "JG"       ! Jericevic/Grisogono - U+S
    character(len=2) :: UnstableKzMethod = "OB" ! O'Brien
    character(len=2) :: StableKzMethod   = "JG" ! Jericevic/Grisogono
                                                !"BW"   ! Brost Wynngard
                                                !"Sb"   ! Seibert
  end type PBL_t
  type(PBL_t), public, save :: PBL = PBL_t()

  type, private :: EmBio_t
    character(len=10) :: GlobBvocMethod = 'GLC-CLM' ! can be MEGAN
    real :: IsopFac = 1.0                     ! for experiments
    real :: TerpFac = 1.0                     ! for experiments
  ! canopy light factor, 1/1.7=0.59, based on Lamb 1993 (cf MEGAN 0.57)
    real :: CLF     = 0.59                    ! canopy factor, leaf vs branch emissions
  end type EmBio_t
  type(EmBio_t), public, save :: EmBio = EmBio_t()

 ! NATBIO allows rcbio in CM_Reactions, but we access elements with
  ! the natbio indices here. These much match the indices used in rcbio
  ! We only use rcbio for isoprene and terpenes so far,  since
  ! soil NO, NH3 emissions etc are dealt with through rcemis.

  type, private :: natbio_t
    integer :: C5H8 = 1
    integer :: TERP = 2
    integer :: Nrcbio = 2  ! No. of rcbio defined in ChemFields/Biogenics_mod
    integer :: NO   = 3    ! used for EmisNat etc
    integer :: NH3  = 4
  end type natbio_t
  type(natbio_t), public, parameter :: NATBIO = natbio_t()

 ! We allow a flexible string which can switch between different
 ! experiments called by e.g. Solver. A but crude, but
 ! it makes sure the experiments are recorded in the config
 ! system

  character(len=100), save, public :: YieldModifications = 'VBS-T10' ! Default for EmChem16mt

  
  type, private :: LandCoverInputs_t
    character(len=TXTLEN_FILE), dimension(2) :: MapFile = 'NOTSET'  ! Usually PS European + global
    character(len=TXTLEN_FILE) :: LandDefs = 'DataDir/Inputs_LandDefs.csv'   !  LAI, h, etc (was Inputs_LandDefs
    character(len=TXTLEN_FILE) :: Do3seDefs = 'DataDir/Inputs_DO3SE.csv'  !  DO3SE inputs
  end type LandCoverInputs_t
  type(LandCoverInputs_t),target, public, save :: LandCoverInputs=LandCoverInputs_t()


! Namelist controlled:
! Some flags for model setup
!------------ NAMELIST VARIABLES - can be reset by emep_namelist.nml file

logical, private, parameter :: F = .false.
type, public :: emep_useconfig
  character(len=10) :: testname = "STD"
  logical :: &                   !
   ! emissions 
     FOREST_FIRES     = .true.  &!  Forest fire options
    ,EMIS             = .false. &! Uses ESX
    ,GRIDDED_EMIS_MONTHLY_FACTOR = .false. & ! .true. triggers ECLIPSE monthly factors
    ,DEGREEDAY_FACTORS = .true. &! will not be used if not found or global grid
    ,EMISSTACKS       = .false. &!
    ,BVOC             = .true.  &!triggers isoprene and terpene emissions
    ,SEASALT          = .true.  &!
    ,CONVECTION       = .false. &! false works best for Euro runs
    ,AIRCRAFT_EMIS    = .true.  &! Needs global file, see manual 
    ,LIGHTNING_EMIS   = .true.  &! 
    ,ROADDUST         = .false. &! TNO Road Dust routine. So far with simplified "climate-correction" factor
    ,DUST             = .false. &! Experimental
    ,EURO_SOILNOX     = .true.  &! ok, but diff for global + Euro runs
    ,GLOBAL_SOILNOX   = .false. &! Need to design better switch
    ,CAMS81_SOILNOX   = .false. &! Need to design better switch
    ,SOILNOX          = .true.  & ! DO NOT ALTER: Set after config
    ,OCEAN_DMS        = .false. &!set automatically true if found.
    ,OCEAN_NH3        = .false. &!set automatically true if found
    ,SOILNH3          = .false. &! DUMMY VALUES, DO NOT USE!
    ,ASH          = .true.  &! Ash from historical Volcanic Eruption
    ,PreADV       = .false. &! Column Emissions are preadvected when winds are very strong 
    ,NOCHEM       = .false. &! Turns of fchemistry for emergency runs
    ,POLLEN       = .false. &! EXPERIMENTAL. Only works if start Jan 1
    ,PROGRESS_FILES   = .false. &! write to file.msg for each output in file.nc
    ,SKIP_INCOMPLETE_OUTPUT = .false. & ! skip daily/montly/fullrun output for runs under 1/28/181 days
    ,NO_3DPRECIP      = .false. &! hack 3D precipitation off
    ,SURF_AREA        = .true.  &! For improved aerosol uptake
    ,MACEHEADFIX      = .true.  &! Correction to O3 BCs (Mace Head Obs.)
    ,MACEHEAD_AVG     = .false. &! Uses 10-year avg. Good for e.g. RCA runs.
    ,MINCONC          = .false. &! Experimental. To avoid problems with miniscule numbers
    ,FASTJ            = .false. & ! use FastJ_mod for computing rcphot
    ,AMINEAQ          = .false. & ! MKPS
!    ,ESX              = .false. &! Uses ESX
    ,PFT_MAPS         = .false. &! 
    ,uEMEP            = .false. &! make local fraction of pollutants
    ,LocalFractions   = .false. &! make local fraction of pollutants
    ! meteo related
    ,SOILWATER        = .false. &!
    ,EtaCOORDINATES   = .true.  &! default since October 2014
    ,WRF_MET_NAMES    = .false. &!to read directly WRF metdata
    ,ZREF             = .false. &! testing
    ,RH_FROM_NWP      = .true.  &! Use rh2m, not LE in Submet
    ,TIMEZONEMAP      = .false. & ! Uses new monthly_timezones_GLOBAL05 map
    ,EFFECTIVE_RESISTANCE = .true. ! Drydep method designed for shallow layer

 ! If USES%EMISTACKS, need to set:
  character(len=4)  :: PlumeMethod   = "none" !MKPS:"ASME","NILU","PVDI"

 ! N2O5 hydrolysis
 ! During 2015 the aersol surface area calculation was much improved, and this
 ! leads to the need for new n2o5 hydrolysis  methods. DO NOT USE EmepReimer,
 ! but one of :; 'SmixTen' , 'Smix', 'Gamma:0.002'

  character(len=20) :: n2o5HydrolysisMethod = 'Smix'

! Selection of method for Whitecap calculation for Seasalt
  character(len=15) :: WHITECAPS  = 'Callaghan'
! In development
   logical :: BIDIR       = .false.  !< FUTURE Bi-directional exchange
   character(len=20)      :: BiDirMethod = 'NOTSET'  ! FUTURE
   character(len=20)      :: MonthlyNH3  = 'NOTSET'  ! can be 'LOTOS'
end type emep_useconfig

type(emep_useconfig), public, save :: USES

logical,  public, save :: &
      FORCE_PFT_MAPS_FALSE = .false. !forces PFT_MAPS  = F, even if global grid

integer, parameter, public :: NSECTORS_ADD_MAX=  250  ! Max. total number of additional sector that can be read from config
type(Sector_type), public :: SECTORS_ADD(NSECTORS_ADD_MAX)
type(emis_in), public, dimension(50) :: emis_inputlist = emis_in()
type(Emis_sourceFile_id_type), public, save:: Emis_sourceFiles(20) !as read from config
type(Emis_mask_type), public, save :: EmisMask(10) !emission mask new format
!MaxNSECTORS to allow reading of SecEmisOutWanted before NSECTORS is defined
integer, public, parameter :: MaxNSECTORS = 100 
logical, public, save :: SecEmisOutWanted(MaxNSECTORS) = .false.
logical, public, save :: SecEmisTotalsWanted = .false.

logical, public, save :: EmisSplit_OUT = .false.

logical, public, save :: AOD_WANTED = .false.!set automatically to T, if AOD requested in output

logical, public, save  :: HourlyEmisOut = .false. !to output sector emissions hourly

!Note that we cannot define the settings as logical (T/F), because we need the state "NOTSET" also
character(len=TXTLEN_NAME), public, save :: EUROPEAN_settings = 'NOTSET'! The domain covers Europe
character(len=TXTLEN_NAME), public, save :: GLOBAL_settings = 'NOTSET'!The domain cover other regions

character(len=TXTLEN_FILE), public, save :: &
  EmisDir = '.',  &
  DataDir = '.',  &
  ZCMDIR = 'DataDir/ZCM_DIRS/ZCM_EmChem19',  & ! default EmChem19 - also used for EmChem19a etc
  OwnInputDir = '.',  &  ! user-defined location
  GRID = 'EECCA', & ! default grid
  meteo= 'DataDir/GRID/metdata_EC/YYYY/meteoYYYYMMDD.nc', & ! template for meteofile
  DegreeDayFactorsFile = 'MetDir/HDD18-GRID-YYYY.nc'        ! template for DegreeDayFactors.nc

integer, public, save :: startdate(4)=(/0,0,0,0/),enddate(4)=(/0,0,0,24/) ! start and end of the run
integer, public, save :: out_startdate(4)=(/-1,-1,-1,-1/) ! start of the output of data
integer, public, save :: spinup_enddate(4)=(/-1,-1,-1,-1/) ! end of spinup. Does not average concentration etc before that date

!-----------------------------------------------------------
! Convection factor - reduces convective fluxes (which can be
! too high in some NWPs)
real, public, save :: CONVECTION_FACTOR = 0.33   ! Pragmatic default
!-----------------------------------------------------------
logical, public, save ::             &
  TEGEN_DATA         = .true.        & ! Interpolate global data to make dust if  USES%DUST=.true.
 ,INERIS_SNAP1       = .false.       & !(EXP_NAME=="TFMM"), & ! Switches off decadal trend
 ,INERIS_SNAP2       = .false.       & !(EXP_NAME=="TFMM"), & ! Allows near-zero summer values
 ,ANALYSIS           = .false.       & ! EXPERIMENTAL: 3DVar data assimilation
 ,ZERO_ORDER_ADVEC   = .false.       & ! force zero order horizontal and vertical advection
 ,JUMPOVER29FEB      = .false.         ! When current date is 29th February, jump to next date.

type(uEMEP_type), public, save :: uEMEP ! The parameters steering uEMEP

integer, public, parameter :: MAXSRC=1000
type(lf_sources), public, save :: lf_src(MAXSRC)
integer, public, parameter :: Max_Country_list=100
character(len=10), public, save :: lf_country_list(Max_Country_list)='NOTSET'!new format "uEMEP" Local Fractions. List of countries
integer, public, parameter :: Max_Country_groups=30
type(lf_country_group_type), public, save :: lf_country_group(Max_Country_groups)
integer, public, parameter :: Max_Country_sectors=13
integer, public, save :: lf_country_sector_list(Max_Country_sectors)=-1!new format "uEMEP" Local Fractions. List of sectors for each country

integer, public, save :: &
  FREQ_HOURLY = 1  ! 3Dhourly netcdf special output frequency

! Soil NOx. Choose EURO for better spatial and temp res, but for
! global runs need global monthly. Variable USE_SOILNOX set from
! these below.
!
! Also, is scaling needed for EURO_SOILNOX?
! The Euro soil NO emissions are based upon average Nr-deposition calculated
!  for the 2000s, as given in the AnnualNdep.nc files. For future years a
!  new AnnualNdep.nc could be pre-calculated. A simpler but approximate
!  way is to scale with some other factor, e.g. the ratio of emissions over
!  some area (EMEP, or EU) in year YYYY divided by year 2005 values.
! Remember, soil-NO emissions are *very* uncertain.

  real, public, save :: EURO_SOILNOX_DEPSCALE = 1.0 !

!NB: *OCEAN*  are internal variables. Cannot be set manually.
  logical, public, save ::  FOUND_OCEAN_DMS = .false. !set automatically true if found

! Methane background:
!  -1 gives defaults in BoundaryConditions_mod
!  -26,-45 or -85 gives RC26, P45, 85 for iyr_trend
  real, public, save :: BGND_CH4 = -1  
! To skip rct value   (jAero work)
  integer, public, save, dimension(10) :: SKIP_RCT  = -1  ! -1 gives defaults
!

!ColumnsSource config
integer,  public, save ::  &
  NMAX_LOC = 7,   &! Max number of locations on processor/subdomain (increase to 24 for eEMEP)
  NMAX_EMS = 250   ! Max number of events def per location (increase to 6000 for eEMEP)
character(len=TXTLEN_FILE),  public, save :: &
  flocdef="columnsource_location.csv",  & ! see locdef
  femsdef="columnsource_emission.csv"     ! see emsdef
logical,  public, save ::          &
  need_topo    = .true.     ! do not use column emissions if topo file is not found

!Forest Fire config
character(len=4),  public , save:: BBMAP = 'FINN'
integer,  public, save ::    &
  BBverbose=1,        & ! debug verbosity 0,..,4
  persistence=1,    & ! persistence in days
  fire_year=-1        ! override current year
logical,  public, save ::    &
  BBneed_file=.true., & ! stop if don't find file
  BBneed_date=.true., & ! stop if don't find time record
  BBneed_poll=.true.    ! stop if don't find pollutant
character(len=TXTLEN_SHORT),  public, save :: BBMODE="DAILY_REC"
character(len=TXTLEN_FILE),  public, save :: &
  GFED_PATTERN = 'GFED_ForestFireEmis.nc',&
  FINN_PATTERN = 'FINN_ForestFireEmis_v15_YYYY.nc',&
  GFAS_PATTERN = 'GFAS_ForestFireEmis_YYYY.nc'

! Nest config
character(len=TXTLEN_SHORT),public, save ::  &
  NEST_MODE_READ='NONE',&  ! read  mode
  NEST_MODE_SAVE='NONE'    ! write mode
integer, public, save :: NEST_NHOURSAVE=3,NEST_NHOURREAD=1 ! write/read frequency
!if(NEST_NHOURREAD<NEST_NHOURSAVE) the data is interpolated in time

character(len=TXTLEN_FILE),public, target, save ::  &
  NEST_template_read_3D = 'EMEP_IN.nc',&       ! Different paths can be set here
  NEST_template_read_BC = 'EMEP_IN.nc',&       ! for each of the IO IC/BC files,
  NEST_template_write   = 'EMEP_OUT.nc',&      ! on namelist, if needed.
  NEST_template_dump    = 'EMEP_Dump.nc'       ! on namelist, if needed.
logical, public, save ::  &
  NEST_save_append = .false., & ! Append to an exixting NEST_template_write or stop the simulation.
  NEST_save_overwrite = .false. ! Overwrite an exixting NEST_template_write or stop the simulation.
! The run will stop if the file already exists, unless NEST_save_append=T or NEST_save_overwrite=T

integer,save, public ::   BC_DAYS=0   ! #days in one BC file, for use old BCs in a FORECAST
              ! 0 means "do not look for old files"

! Nested input/output on OUTDATE mode
integer,public,parameter  :: OUTDATE_NDUMP_MAX = 4  ! Number of nested output
integer, public, save     :: NEST_OUTDATE_NDUMP     = 0  ! Read by emepctm.f90
! on forecast run (1: start next forecast; 2-4: NMC statistics)
type(date), public :: NEST_outdate(OUTDATE_NDUMP_MAX)=date(-1,-1,-1,-1,-1)

character(len=TXTLEN_FILE),public, save :: NEST_MET_inner ='NOTSET' !path to metdata for inner grid
integer, save, public :: NEST_RUNDOMAIN_inner(4)=-1 ! RUNDOMAIN used for run in inner grid
! Limit output, e.g. for NMC statistics (3DVar)
character(len=TXTLEN_SHORT), public, save :: &
  NEST_WRITE_SPC(NSPEC_ADV)="", NEST_WRITE_GRP(size(chemgroups))=""

!coordinates of subdomain to write, relative to FULL domain (only used in write mode)
integer, public, save :: NEST_out_DOMAIN(4)=-1 ! =[istart,iend,jstart,jend]

logical,public, save ::  &               ! if IC/BC are in the same model/run
  NEST_native_grid_3D = .false.,&              ! grid, the expensive call to
  NEST_native_grid_BC = .false.,&              ! grid2grid_coeff in init_nest can be avoided
  NEST_omit_zero_write= .false.                ! skip const=0.0 variables

logical, public, save :: &
  USE_EXTERNAL_BIC = .false., & ! use external (non emepctm) BCs
  TOP_BC           = .false.    ! BCs include top level
character(len=TXTLEN_SHORT),public, save :: &
  EXTERNAL_BIC_NAME    = "DUMMY", EXTERNAL_BIC_VERSION = "use_any"
character(len=TXTLEN_FILE),public, target, save :: &
  filename_eta     = 'EMEP_IN_BC_eta.zaxis'

!Output_config variables

integer, public, parameter ::       &
  MAX_NUM_DERIV2D = 343,            &
  MAX_NUM_DDEP_ECOS = 9,            & ! Grid, Conif, etc.
  MAX_NUM_NEWMOS  = 30,             & !New system.
  ! Older system
  MAX_NUM_MOSCONCS  = 10,           & !careful here, we multiply by next:
  MAX_NUM_MOSLCS    = 10,           & !careful here, we multiply by prev:
  MAX_NUM_DDEP_WANTED = NSPEC_ADV,  & !plenty big
  MAX_NUM_WDEP_WANTED = NSPEC_ADV     !plenty big


integer, public, parameter :: &
   NSITES_MAX =        99     & ! Max. no surface sites allowed
  ,FREQ_SITE  =         1     & ! Interval (hrs) between outputs
  ,NSHL_SITE_MAX  =    10     & ! No. short-lived species
  ,NXTRA_SITE_MISC =    2     & ! No. Misc. met. params  ( e.g. T2, d_2d)
  ,NXTRA_SITE_D2D  =   18       ! No.  params from d_2d fields 
integer, public, parameter :: NSONDES_MAX = 99 ! Max. no sondes allowed

integer, private :: isite              ! To assign arrays, if needed

!**** Sonde outputs   (used in Sites_mod)
!==============================================================
!     Specify the species to be output to the sondes.out file
!  We typically deal with fewer species for sonde output than
!  surface sites, so we use a different method to specify.
! For met params we have no simple index, so we use characters.
! These must be defined in Sites_mod.f90.

integer, public, parameter :: &
    FREQ_SONDE  =     1        &   ! Interval (hrs) between outputs
   ,NXTRA_SONDE =    4             ! No. Misc. met. params


! Extra parameters - need to be coded in Sites_mod also. So far
! we can choose from hmix, T2, or th (pot. temp.) or d_2d fields.
!  d_2d fields can be accessed from Derived_mod by setting common index
!  "D2D" in SITE_XTRA and the actual field name (as defined in Derived_mod)
!  in SITE_XTRA_CODE (e.g. "D2_PM25 " or "D2_SIA") :

!** IMPORTANT!! Make sure the correspondence between selected for output
!** fields in SITE_XTRA and their names in SITE_XTRA_CODE

character(len=24), public, parameter, dimension(NXTRA_SITE_MISC) :: &
  SITE_XTRA_MISC=[character(len=18):: "th","T2"]

character(len=TXTLEN_SHORT), public :: SITE_SHL_names(NSPEC_SHL) = 'NOTSET'
character(len=TXTLEN_SHORT), public :: SONDE_SHL_names(NSPEC_SHL) = 'NOTSET'
character(len=TXTLEN_SHORT), public :: SONDE_ADV_names(NSPEC_ADV) = 'NOTSET'

!These variables must have been set in My_Derived for them to be used.
character(len=24), public, parameter, dimension(NXTRA_SITE_D2D) :: &
  SITE_XTRA_D2D=[character(len=24):: &
    "HMIX", & !Hmix is interpolated in time, unlike NWP version 
    "Emis_mgm2_BioNatC5H8","Emis_mgm2_BioNatTERP",&
    "Emis_mgm2_BioNatNO","Emis_mgm2_nox",&
    'WDEP_PREC',&!''SNratio',&
    'met2d_ps', 'met2d_uref','met2d_u10', & !u10 seems to be wind-speed
                                          !'met2d_v10','met2d_rh2m', &
    'met2d_SMI_uppr', 'met2d_SMI_deep',&
    'met2d_ustar_nwp', 'met2d_LH_Wm2', 'met2d_SH_Wm2',&
    'USTAR_DF','INVL_DF', &
    'met2d_PARdbh', 'met2d_PARdif' &
]
character(len=10), public, parameter, dimension(NXTRA_SONDE) :: &
  SONDE_XTRA= [character(len=10):: &
   "NOy","z_mid","p_mid","th"]!,"Kz_m2s"]



! Site/Sondes (under construction. DO NOT USE!)
integer, parameter :: MAX_NEXTRA_SITED2D=100
type, private :: sites_t
  integer :: freq_site = 1
  integer :: nmax = 99
  integer :: nadv = 0
  integer :: nshl = 0
  integer :: nd2d = 0
  integer :: nmisc = 0
  integer, allocatable, dimension(:) :: adv 
  integer, allocatable, dimension(:) :: shl
  integer, allocatable, dimension(:) :: d2d
  integer, allocatable, dimension(:) :: misc
end type sites_t
type(sites_t),save :: site_outputs, sonde_outputs
!character(len=24), public, save, dimension(MAX_NEXTRA_SITED2D) :: &
!   site_outputs_extraD2D = '-', sonde_outputs_extraD2D = '-'


type(Deriv), public, save, dimension(MAX_NUM_DERIV2D) :: OutputMisc= Deriv()
type(typ_s5ind), public, save, dimension(MAX_NUM_DERIV2D) :: &
  OutputConcs = typ_s5ind("-","-","-","-","-","-")

integer, parameter, private :: MAXNVO3  = 60
type(O3cl_t), public, save, dimension(MAXNVO3) :: OutputVegO3 = O3cl_t()

! Depositions
type(typ_s1ind), public, save, dimension(MAX_NUM_DDEP_ECOS) :: &
  DDEP_ECOS = typ_s1ind("-",'-') ! e.g. "Grid","YMD", 

type(typ_s3), public, save, dimension(MAX_NUM_DDEP_WANTED) :: &
  DDEP_WANTED = typ_s3('-','-','-'), & ! e.g. typ_s3("SO2",SPEC,"mgS"),
  SDEP_WANTED = typ_s3('-','-','-')    ! Stomatal deposition (for HTAP)

type(typ_s4), public, save, dimension(MAX_NUM_WDEP_WANTED) :: &
  WDEP_WANTED = typ_s4('-','-','-','-')

!- specify some species and land-covers we want to output
! dep. velocities for in netcdf files. Set in My_DryDep_mod.
! NewMosaic seems to mean new-style, to avoid  needing all combinations
! of MET & LC
type(typ_s5ind), public, save, dimension(MAX_NUM_NEWMOS) :: &
  NewMosaic = typ_s5ind('-','-','-','-','-','-')


! For met-data and canopy concs/fluxes ...
character(len=TXTLEN_DERIV), public, save, dimension(MAX_NUM_MOSCONCS) :: &
  MOSAIC_METCONCS = '-' ! = [character(len=TXTLEN_DERIV):: & 
     !,"VPD","FstO3","EVAP","Gsto" ,"USTAR","INVL"/)
! "USTAR","LAI","CanopyO3","FstO3"] ! SKIP CanopyO3
! "g_sto" needs more work - only set as L%g_sto

character(len=TXTLEN_DERIV), public, save, dimension(MAX_NUM_MOSLCS) :: &
  MET_LCS = '-' 
! [character(len=TXTLEN_DERIV)::  "DF","GR","BF","TC","IAM_DF","IAM_CR"]

character(len=10), public, save ::  Mosaic_timefmt='YM'  ! eg 'YMD'

!Machine_config variables
 character (len=TXTLEN_FILE), public :: DataPath(20) = 'NOTSET'
!
!Extra namelists
 character (len=TXTLEN_FILE), public :: ExtraConfigFile(20) = 'NOTSET'

!------------ END OF NAMELIST VARIABLES ------------------------------------!

! Some flags for model setup will be removed when code is sufficiently tested
! (for convection use foundconv in permanent code)
logical, public, parameter ::         &
  NO_CROPNH3DEP      = .true.,        & ! Stop NH3 deposition for growing crops
  EXTENDEDMASSBUDGET = .false.,       & ! extended massbudget outputs
  LANDIFY_MET        = .false.


! Boundary layer profiles. IN-TESTING
character(len=4), parameter, public :: &
  FluxPROFILE = "Iter"
! FluxPROFILE = "Ln95"  ! use Launiainen1995 EXPERIMENTAL. Fails in some areas

!The GEA emission data, which is used for EUCAARI runs on the HIRHAM domains
!have in several sea grid cells non-zero emissions in other sectors than SNAP8
!and there are also NH3 emission over sea areas. The former problem makes
!the code crash if the sea areas are defined  as sea (sea=T), so we treat
!them as land in the EUCAARI/HIRHAM runs (sea=F). This is a problem with GEA
!emission data only, not the HIRHAM domain! When e.g. interpolated EMEP emissions
!are used on the HIRHAM domain, this is not a problem.

logical, public, save :: &
  SEAFIX_GEA_NEEDED = .false. ! only if problems. Read from Model_config

!=============================================================================
!+ 1) Define first dimensions that might change quite often -  for different
!     run domains

!IN-TESTING (reset in NML if wanted)
!) Emissions. Standard =ascii emislist. CdfFractions possible for INERIS
!  and new cdf emission system in testing. Reset in config_ files
Logical , save, public :: &
  EMIS_OUT    = .false.     ! output emissions in separate files (memory demanding)

integer, public :: &  ! Full domain, set automatically from meteorology
  KMAX_MID, &         ! Number of points (levels) in vertical
  KMAX_BND, &         ! Number of level boundaries (=KMAX_MID+1)
  IIFULLDOM,JJFULLDOM ! Full grid

! Sub-domains, in fulldomain coordinates. Read on Model_config
! negative values means: rundomain size if not other specified.
integer, public, save, dimension(4) :: &
  RUNDOMAIN =-999,  & ! run sub-domain
  fullrun_DOMAIN=-999,   & ! fullrun (year) output sub-domain
  month_DOMAIN=-999,     & ! montly output sub-domain
  day_DOMAIN=-999,       & ! daily  output sub-domain
  hour_DOMAIN =-999        ! hourly output sub-domain

! 3D Output: all modell levels will be outputed by default
! see Init_My_Deriv and OutputSize_config for details
integer, public, save ::  &
  num_lev3d=0,lev3d(60)=0 ! numbers of levels,list of levels
logical, public, save ::  &
    lev3d_from_surface=.false. ! levels are to be read from surface up

integer, public, save ::  & ! Actual number of processors in longitude, latitude
  NPROCX, NPROCY, NPROC     ! and total. NPROCY must be 2 for GLOBAL runs.

CHARACTER(LEN=3), public, save :: &
  DOMAIN_DECOM_MODE=''      ! override parinit(Pole_singular) option (Par_mod)

!=============================================================================
!+ 2) Define  debug flags.

! We have one variable, to say if we are on master-processor
! or not: (kept here to avoid too many dependencies for box-model
! codes which don't need Par_mod. 
! See also DebugCell from Debug_mod

logical, public, save ::  MasterProc = .true.

!=============================================================================
! Some flags for model setup

! Debug flag DEBUG_XXX  applied in subroutine XXX
! logical, public, parameter :: PALEO_TEST = .false. 

!=============================================================================
! 3)  Source-receptor runs?
! We don't (generally) want daily outputs for SR runs, so in
! Derived_mod, we set all IOU_DAY false if SOURCE_RECPTOR = .true..

logical, public, save :: SOURCE_RECEPTOR = .false., VOLCANO_SR=.false.

! Compress NetCDF output? (nc4 feature, use -1 for netcdf3 output)
integer, public, save :: NETCDF_DEFLATE_LEVEL=4

!Hourly output in single file or monthly/daily files:
!NB: will not work well by default on Stallo per 14th Feb 2012 because of library bugs!
!Until this is fixed, you must compile with netcdf/4.1.3 and link and run with compiler 12.1.2
character(len=30), public, save :: &    ! ending depeding on date:
! HOURLYFILE_ending="YYYYMM.nc"    ! MM  -> month (01 .. 12)
! HOURLYFILE_ending="YYYYMMDD.nc"  ! DD  -> day of the month (00 .. 31)
! HOURLYFILE_ending="JJJ.nc"   ! JJJ -> the day of the year (001 .. 366)
     HOURLYFILE_ending=".nc"           ! default, just one file
!do not use  HOURLYFILE_ending="_hourExtra.nc"     ! keep the same for the whole run

! NH3 module as set up originally with U10 from met: kept for safety only.
! Will be replaced by sub.grid calculation of wind in future.
! Keep false until code re-implemented
logical, public, parameter :: NH3_U10 = .false.

!=============================================================================
!+ 4)  Define main model dimensions,  things that will
!       generally only change when switching Met-driver
integer, public, parameter ::  &
!TREEX  NLANDUSEMAX  = 19   &   ! Number of land use types in Inputs.Landuse file
  NLANDUSEMAX  = 40    &    ! Max num land use types in Inputs.Landuse file
, KTOP         = 1     &    ! K-value at top of domain
, KWINDTOP     = 5     &    ! Define extent needed for wind-speed array
, NMET         = 2     &    ! No. met fields in memory
, KCHEMTOP     = 2     &    ! chemistry not done for k=1
, KCLOUDTOP    = 8     &    ! limit of clouds (for MADE dj ??)
, KUPPER       = 6          ! limit of clouds (for wet dep.)

integer, public :: METSTEP = 3  ! time-step of met. (h). 3 hours default, but can be reset by metdata
real, public :: Zmix_ref = 50.0 !height at which concentration above different landuse are considered equal 

!> Namelist controlled: which veg do we want flux-outputs for
!! We will put the filename, and params (SGS, EGS, etc) in
!! the _Params array.
character(len=TXTLEN_SHORT), public, save, dimension(20) ::  &
   FLUX_VEGS=""    & ! e.g. WinterWheat
  ,FLUX_IGNORE=""  & ! e.g. Water, desert..
  ,VEG_2dGS=""
character(len=99), public, save, dimension(10) :: VEG_2dGS_Params=""
integer, public, save :: nFluxVegs = 0 ! reset in Landuse_mod

! To use external maps of plant functional types we need to
! map between EMEP codes and netcdf file codes
!NOT USED YET.
type(typ_ss), public, save, dimension(NLANDUSEMAX) :: &
  PFT_MAPPINGS=typ_ss('-','-')

real, public, save :: &
  dt_advec = -999.9,   & ! time-step for advection (s), grid resolution dependent
  dt_advec_inv  ! =1/dt_advec

! NTDAY:  Number of 2D O3 to be saved each day (for SOMO)
! 24/NTDAY is the time integration step for SOMO
! large value -> large memory use; too small value -> bad approx. for SOMO
! NB must be choosen:  24*3600/dt_advec <= NTDAY >=3 and
! preferably an integer fraction of 24*3600/dt_advec
integer, public, parameter :: NTDAY = 72

!/-- choose temperature range: from 148 K (-125C) ro 333K (+60C).
integer, parameter, public :: CHEMTMIN=148,CHEMTMAX=333

real, public, parameter :: &
  V_RAIN     = 5.           !approximate vertical speed of rain m/s

real, public, parameter :: &
  CW_THRESHOLD = 1.0E-7&!Cloudwater (kg/kg); above threshold allow possibility
                        ! for precipitations. Value could be adjusted.
, RH_THRESHOLD = 0.85  &!Relative humidity (fraction); above threshold allow
                        !possibility for precipitations.Value could be adjusted.
, CW2CC = 1.0E6         !Converts Cloudwater (kg/kg) into CloudCover in %
                        !Value could be adjusted.
!
!  additional parameters
!
integer, public, save   :: &
     nmax,  & ! Number of dt_advec steps in one METSTEP
     step_main, & ! Main time loop count
     iyr_trend ! Year specified for say BC changes

character(len=120), public, save :: runlabel1&!SHORT Allows explanatory text
                                  , runlabel2 !LONG  Read in from grun.pl

! Typically, we define as mainly sea when > 50% water, and
! likely_coastal when > 20%. See Landuse_mod
real, public, parameter, dimension(2) ::  SEA_LIMIT = (/ 0.2, 0.5 /)

real, public, parameter :: &
  MINCONC  = 1.0e-25   ! Experimental. To avoid problems with miniscule numbers
                       ! so that DEBUG runs have same values as normal. Uses if
                       ! USES%MINCONC = T, set in config

real, public, parameter :: &
  EPSIL  = 1.0e-30         &  ! small number
, PASCAL = 100.0           &  ! Conv. from hPa to Pa
, PPB    = 1.0e-9          &  ! parts per billion (mixing ratio)
, PPBINV = 1.0e+9          &
, PPT    = 1.0e-12         &  ! parts per trillion (mixing ratio)
, PPTINV = 1.0e+12         &
, PT     = 1.0e+4          &  ! Top of model region = 10000 Pa = 100 hPa
, TINY   = 1.0e-9             ! -1.E-9" is sometimes used  in order to avoid
                              ! different roundings on different machines.
real, public :: Pref   = 101325.0  ! Reference pressure in Pa used to define vertical levels


! Define output types.
!   Derived output types: types 1..6 (instantaneous,year,month,day,hour,hour_inst),
!                         refer to output variables defined in Derived_mod.
!   Hourly  output types: types 7..8 (hourly_out inst.,hourly_out_mean),
!                         refer to output variables defined in My_Outputs_mod.
integer, public, parameter ::  &
  IOU_INST=1,IOU_YEAR=2,IOU_MON=3,IOU_DAY=4,IOU_HOUR=5,IOU_HOUR_INST=6, & ! Derived output
  IOU_MAX_MAX=6                             ! Max values for of IOU (for array declarations)

character, public, parameter ::  & ! output shorthands, order should match IOU_*
  IOU_KEY(IOU_YEAR:IOU_HOUR_INST)=['Y','M','D','H','I']

character(len=*), public, parameter :: model="EMEP_MSC-W "
character(len=TXTLEN_FILE),target, public :: fileName_O3_Top = "NOTSET"
! Can use RCP values of CH4 for given iyr_trend.
character(len=TXTLEN_FILE),target, public :: fileName_CH4_ibcs = "NOTSET" ! eg ch4_rcp45

logical, public, parameter:: MANUAL_GRID=.false.!under developement.

!file names
type, public ::names
character(len=TXTLEN_FILE), pointer:: filename => null()
end type names
integer, public, parameter :: Size_InputFiles = 60
type(names), public, save :: InputFiles(Size_InputFiles)

!To add a new filename:
!1) add a line just here below, XXFile = '/default/Path/Default.name'
!2) add the XXFile in NAMELIST /Model_config/
!3) add a call associate_File(XXFile) near the end of Config_Constants
!4) In the routine using the file, add the XXFile under  "use Config_module"
!5) replace the name you used in the routine with XXFile
character(len=TXTLEN_FILE), target, save, public :: femisFile = 'DataDir/femis.dat'
character(len=TXTLEN_FILE), target, save, public :: Vertical_levelsFile = 'DataDir/Vertical_levels20_EC.txt'
character(len=TXTLEN_FILE), target, save, public :: EmisHeightsFile = 'DataDir/EmisHeights.txt'
character(len=TXTLEN_FILE), target, save, public :: SoilTypesFile = 'DataDir/SoilTypes_IFS.nc'
character(len=TXTLEN_FILE), target, save, public :: SurfacePressureFile = 'DataDir/SurfacePressure.nc'
character(len=TXTLEN_FILE), target, save, public :: AircraftEmis_FLFile = 'DataDir/AircraftEmis_FL.nc'
character(len=TXTLEN_FILE), target, save, public :: nox_emission_1996_2005File = 'DataDir/nox_emission_1996-2005.nc'
character(len=TXTLEN_FILE), target, save, public :: nox_emission_cams81File = 'DataDir/cams81_monthly_SoilEmissions_v2.2_GLOBAL05_2000_2018.nc'
!character(len=TXTLEN_FILE), target, save, public :: MonthlyFacFile = 'DataDir/inputs_emepdefaults_Jun2012/MonthlyFac.POLL'
!2021: added ECLIPSE6b-based factors for non-European areas
character(len=TXTLEN_FILE), target, save, public :: MonthlyFacFile = 'DataDir/Timefactors/MonthlyFacs_eclipse_V6b_snap_xJun2012/MonthlyFacs.POLL'
!POLL replaced by name of pollutant in Timefactors_mod
character(len=TXTLEN_FILE), target, save, public :: DailyFacFile = 'DataDir/inputs_emepdefaults_Jun2012/DailyFac.POLL'
character(len=TXTLEN_FILE), target, save, public :: HourlyFacFile = 'DataDir/inputs_emepdefaults_Jun2012/HourlyFacs.INERIS'
character(len=TXTLEN_FILE), target, save, public :: HourlyFacSpecialsFile = 'NOTSET'
! Chemical schemes have specific files:
!character(len=*), parameter :: ZCMDIR= 'DataDir/ZCM_CRI-R5-emep/'
character(len=TXTLEN_FILE), target, save, public :: &
  cmxbicDefaultFile          = 'ZCMDIR/CMX_BoundaryConditions.txt'   &
 ,cmxBiomassBurning_FINNv1p5 = 'ZCMDIR/CMX_BiomassBurning_FINNv1p5.txt' &
 ,cmxBiomassBurning_GFASv1   = 'ZCMDIR/CMX_BiomassBurning_GFASv1.txt' &
!POLL replaced by name of pollutant in EmisSplit
 ,SplitDefaultFile           = 'ZCMDIR/emissplit_run/emissplit.defaults.POLL' &
 ,SplitSpecialsFile          = 'ZCMDIR/emissplit_run/emissplit.specials.POLL'
character(len=TXTLEN_FILE), target, save, public :: RoadMapFile = 'DataDir/RoadMap.nc'
character(len=TXTLEN_FILE), target, save, public :: AVG_SMI_2005_2010File = 'DataDir/AVG_SMI_2005_2010.nc'
character(len=TXTLEN_FILE), target, save, public :: Soil_TegenFile = 'DataDir/Soil_Tegen.nc'
! default site/sond files use lat, lon and Kdown coords:
character(len=TXTLEN_FILE), target, save, public :: SitesFile = 'DataDir/sitesLLKD.dat'
character(len=TXTLEN_FILE), target, save, public :: SondesFile = 'DataDir/sondesLLKD.dat'
character(len=TXTLEN_FILE), target, save, public :: GLOBAL_LAInBVOCFile = 'DataDir/GLOBAL_LAInBVOC.nc'
character(len=TXTLEN_FILE), target, save, public :: EMEP_EuroBVOCFile = 'DataDir/LandInputs_Mar2011/EMEP_EuroBVOC.nc'
!SEASON replace by 'jan', 'apr', 'jul' or 'oct' in readdiss
character(len=TXTLEN_FILE), target, save, public :: jclearFile = 'DataDir/jclear.SEASON'
!SEASON replace by 'jan', 'apr', 'jul' or 'oct' in readdiss
character(len=TXTLEN_FILE), target, save, public :: jcl1kmFile = 'DataDir/jcl1.SEASON'
!SEASON replace by 'jan', 'apr', 'jul' or 'oct' in readdiss
character(len=TXTLEN_FILE), target, save, public :: jcl3kmFile = 'DataDir/jcl3.SEASON'
character(len=TXTLEN_FILE), target, save, public :: NdepFile = 'DataDir/AnnualNdep_PS50x_EECCA2005_2009.nc'
!MM replace by month in lightning()
character(len=TXTLEN_FILE), target, save, public :: lightningFile = 'DataDir/lt21-nox.datMM'
character(len=TXTLEN_FILE), target, save, public :: LoganO3File = 'DataDir/Logan_P.nc'
character(len=TXTLEN_FILE), target, save, public :: DustFile = 'DataDir/Dust2014_month.nc'
character(len=TXTLEN_FILE), target, save, public :: TopoFile = 'DataDir/GRID/topography.nc'
character(len=TXTLEN_FILE), target, save, public :: BiDirInputFile = 'NOTSET' ! FUTURE
character(len=TXTLEN_FILE), target, save, public :: Monthly_patternsFile = 'DataDir/ECLIPSEv5_monthly_patterns.nc'
character(len=TXTLEN_FILE), target, save, public :: Monthly_timezoneFile = 'DataDir/Timefactors/monthly_timezones_GLOBAL05.nc'

!----------------------------------------------------------------------------
contains
subroutine Config_Constants(iolog)
  integer, intent(in) :: iolog ! for Log file

  integer :: i, j, ispec, iostat
  logical,save :: first_call = .true.
  character(len=len(meteo)) ::  MetDir='./' ! path from meteo
  character(len=*), parameter ::  dtxt='Config_MC:'
  character(len=100 ) :: logtxt

  NAMELIST /Model_config/ &
    DegreeDayFactorsFile, meteo & !meteo template with full path
   ,END_OF_EMEPDAY &
   ,EXP_NAME &  ! e.g. EMEPSTD, FORECAST, TFMM, TodayTest, ....
   ,USES   & ! 
   ,AERO   & ! for aerosol equilibrium scheme
   ,PBL    & ! 
   ,EmBio  & ! 
   ,YieldModifications &  ! Allows dynamic change of chemical yields
   ,LandCoverInputs    &  ! for CLM, etc
   ,DEBUG  & !
   ,CONVECTION_FACTOR &
   ,EURO_SOILNOX_DEPSCALE &
   ,uEMEP & !old format . Avoid, will be removed in future versions
   ,lf_src & !new format "uEMEP" Local Fractions
   ,lf_country_list & !new format "uEMEP" Local Fractions. List of countries
   ,lf_country_group & !new format "uEMEP" Local Fractions. List of group of countries
   ,lf_country_sector_list & !new format "uEMEP" Local Fractions. List of sectors for each country
   ,INERIS_SNAP1, INERIS_SNAP2 &   ! Used for TFMM time-factors
   ,FREQ_HOURLY           &
   ,ANALYSIS, SOURCE_RECEPTOR, VOLCANO_SR &
   ,SEAFIX_GEA_NEEDED     & ! only if problems, see text above.
   ,BGND_CH4              & ! Can reset background CH4 values
   ,SKIP_RCT              & ! Can  skip some rct
   ,SECTORS_ADD           & ! additional definitions of sectors
   ,EMIS_OUT, emis_inputlist, EmisDir&
   ,EmisSplit_OUT         & ! Output of species emissions
   ,ZCMDIR                & ! location of emissplit and CMXfiles
   ,OwnInputDir           &  !
   ,Emis_sourceFiles      & ! new format
   ,EmisMask              & ! new format
   ,SecEmisOutWanted      & ! sector emissions to include in output
   ,SecEmisTotalsWanted    & ! give total per sectors and countries
   ,HourlyEmisOut         & ! to output emissions hourly
   ,FLUX_VEGS             & ! Allows user to add veg categories for eg IAM ouput
   ,FLUX_IGNORE           & ! Specify which landcovers don't need FLUX
   ,VEG_2dGS              & ! Allows 2d maps of growing seasons
   ,VEG_2dGS_Params       & ! Allows 2d maps of growing seasons
   ,PFT_MAPPINGS          & ! Allows use of external LAI maps
   ,NETCDF_DEFLATE_LEVEL,  RUNDOMAIN, DOMAIN_DECOM_MODE &
   ,JUMPOVER29FEB, HOURLYFILE_ending  &
   ,dt_advec              & ! can be set to override dt_advec
   ,METSTEP &
   ,ZERO_ORDER_ADVEC &! force zero order horizontal and vertical advection 
   ,EUROPEAN_settings & ! The domain covers Europe -> 
   ,GLOBAL_settings & ! The domain cover other regions too -> Convection
   ,fileName_O3_Top&
   ,fileName_CH4_ibcs&
   ,femisFile&
   ,Vertical_levelsFile&
   ,Emis_Plevels&
   ,Emis_h&
   ,SoilTypesFile&
   ,SurfacePressureFile&
   ,AircraftEmis_FLFile&
   ,nox_emission_1996_2005File&
   ,nox_emission_cams81File&
   ,MonthlyFacFile&
   ,DailyFacFile&
   ,HourlyFacFile&
   ,HourlyFacSpecialsFile&
   ,cmxbicDefaultFile&
   ,cmxBiomassBurning_FINNv1p5&
   ,cmxBiomassBurning_GFASv1&
   ,SplitDefaultFile&
   ,SplitSpecialsFile&
   ,RoadMapFile&
   ,AVG_SMI_2005_2010File&
   ,Soil_TegenFile&
   ,SitesFile&
   ,SondesFile&
   ,GLOBAL_LAInBVOCFile&
   ,EMEP_EuroBVOCFile&
   ,jclearFile&
   ,jcl1kmFile&
   ,jcl3kmFile&
   ,NdepFile&
   ,lightningFile&
   ,BiDirInputFile&
   ,LoganO3File&
   ,DustFile&
   ,TopoFile&
   ,Monthly_patternsFile&
   ,Monthly_timezoneFile&
   ,GRID,iyr_trend,runlabel1,runlabel2,startdate,enddate&
   ,NMAX_LOC,NMAX_EMS,flocdef,femsdef,need_topo&
   ,BBMODE,BBverbose,persistence,fire_year&
   ,BBneed_file,BBneed_date,BBneed_poll&
   ,BBMAP,GFED_PATTERN,FINN_PATTERN,GFAS_PATTERN&
   ,DataPath&
   ,ExtraConfigFile&
   ,NEST_MODE_READ,NEST_MODE_SAVE,NEST_NHOURREAD,NEST_NHOURSAVE &
   ,NEST_template_read_3D,NEST_template_read_BC,NEST_template_write&
   ,NEST_template_dump,BC_DAYS,NEST_save_append,NEST_save_overwrite&
   ,NEST_native_grid_3D,NEST_native_grid_BC,NEST_omit_zero_write,NEST_out_DOMAIN&
   ,NEST_MET_inner,NEST_RUNDOMAIN_inner&
   ,NEST_WRITE_SPC,NEST_WRITE_GRP,NEST_OUTDATE_NDUMP,NEST_outdate&
   ,USE_EXTERNAL_BIC,EXTERNAL_BIC_NAME,EXTERNAL_BIC_VERSION,TOP_BC,filename_eta&
   ,OutputMisc,OutputConcs,OutputVegO3&
   ,DDEP_ECOS, DDEP_WANTED, WDEP_WANTED,SDEP_WANTED&
   ,NewMosaic, MOSAIC_METCONCS, MET_LCS, Mosaic_timefmt&
   ,fullrun_DOMAIN,month_DOMAIN,day_DOMAIN&
   ,hour_DOMAIN, out_startdate, spinup_enddate&
   ,num_lev3d,lev3d,lev3d_from_surface&
   ,LAST_CONFIG_LINE &
   ,SITE_SHL_names,SONDE_SHL_names,SONDE_ADV_names

  LAST_CONFIG_LINE_DEFAULT = LAST_CONFIG_LINE !save default value
  DataPath(1) = '.'!default
  Emis_h = -1.0 ! to mark as not set
  Emis_Plevels = -1.0 ! to mark as not set
  open(IO_NML,file='config_emep.nml',delim='APOSTROPHE')
  read(IO_NML,NML=Model_config)
  ! do not close(IO_NML), other modules will be read namelist on this file
  if(MasterProc) write(*,*) dtxt//'DataPath',trim(DataPath(1))
  
  USES%SOILNOX = USES%EURO_SOILNOX .or. USES%GLOBAL_SOILNOx
  if(MasterProc) then
    write(logtxt,'(a,L2)') dtxt//'USES%SOILNOX ', USES%SOILNOX
    write(*,*) trim(logtxt), IOLOG
    write(IO_LOG,*) trim(logtxt)  ! Can't call PrintLog due to circularity
  end if
 
  USES%LocalFractions = USES%LocalFractions .or. USES%uEMEP !for backward compatibility

  ! Convert DEBUG%SPEC to index
  if(first_call)then
    if(MasterProc)&
         write(iolog,*) 'CHEM SCHEMES: ',trim(CM_schemes_ChemSpecs)

    ispec = find_index( DEBUG%SPEC, species(:)%name )
  ! print *, "debug%spec testing", ispec, trim(DEBUG%SPEC)
    call CheckStop(ispec<1,"debug%spec not found"//trim(DEBUG%SPEC))
    DEBUG%ISPEC = ispec
    first_call = .false.
  end if

  if(MasterProc)then
    write(*, * ) dtxt//"NAMELIST START "
    write(*,*)   dtxt//"LAST LINE after 1st config:"//trim(LAST_CONFIG_LINE)
!    write(iolog,*) dtxt//"NAMELIST IS "
!    write(iolog, NML=Model_config)
  end if

  do i=1,size(DataPath)
    if(DataPath(i)=="NOTSET")then
      if(MasterProc)then
        write(*,*)dtxt//'WARNING: Could not find valid DataDir. Tried:'
        do j=1,i-1
          write(*,*)trim(DataPath(j))
        end do
        stop
      end if
      exit
    end if
!   INQUIRE(...) does not behave consistently across intel/gfortran
    open(IO_TMP,file=trim(DataPath(i)),iostat=iostat,action='read')! does not work without action='read'
    if(iostat==0)then
      DataDir=trim(DataPath(i))
      if(MasterProc)write(*,*)dtxt//'DataDir set to',trim(DataDir)
      close(IO_TMP)
      exit
    end if
  end do

!before any conversion, read the additional namelists
  do i = 1, size(ExtraConfigFile)
     if(ExtraConfigFile(i)/="NOTSET")then
        !NB: replacements have not been made yet
        ExtraConfigFile(i) = key2str(ExtraConfigFile(i),'DataDir',DataDir)
        ExtraConfigFile(i) = key2str(ExtraConfigFile(i),'GRID',GRID)
        ExtraConfigFile(i) = key2str(ExtraConfigFile(i),'OwnInputDir',OwnInputDir)
        ExtraConfigFile(i) = key2str(ExtraConfigFile(i),'EmisDir',EmisDir)
        if(MasterProc) then
         write(*,*) dtxt//'Also reading namelist ',i,trim(ExtraConfigFile(i))
         write(*,*) dtxt//"LAST LINE:"//trim(LAST_CONFIG_LINE) ! for debugs
         write(iolog,*)'Also reading namelist ',i,trim(ExtraConfigFile(i))
        end if
        open(IO_tmp,file=trim(ExtraConfigFile(i)),delim='APOSTROPHE')
        read(IO_tmp,NML=Model_config)
        if(MasterProc) write(*,*) dtxt//'DataPath ExtraConf:', i, trim(DataPath(1))
        close(IO_tmp)
     endif
  enddo
  if(LAST_CONFIG_LINE==LAST_CONFIG_LINE_DEFAULT)then
     if(MasterProc) write(*,*) dtxt//"WARNING: LAST_CONFIG_LINE not modified"
     if(MasterProc) write(*,*) dtxt//"Probable syntax error in namelist!!!"
  else
     if(MasterProc) write(*,*) dtxt//"LAST LINE final:"//trim(LAST_CONFIG_LINE)
  endif
!EEEEEEEEEEEEEEEEEEEEEEEEE

  meteo = key2str(meteo,'DataDir',DataDir)
  meteo = key2str(meteo,'EmisDir',EmisDir)
  meteo = key2str(meteo,'GRID',GRID)
  MetDir= key2str(meteo,'meteoYYYYMMDD.nc','./')
  DegreeDayFactorsFile=key2str(DegreeDayFactorsFile,'MetDir',MetDir)
  DegreeDayFactorsFile=key2str(DegreeDayFactorsFile,'DataDir',DataDir)
  DegreeDayFactorsFile=key2str(DegreeDayFactorsFile,'GRID',GRID)
  DegreeDayFactorsFile=key2str(DegreeDayFactorsFile,'YYYY',startdate(1))
  if(MasterProc)then
    write(*,*)dtxt//'Defined DegreeDayFactorsFile as:'
    write(*,*)trim(DegreeDayFactorsFile)
  end if

 ! LandCoverInputs
  do i = 1, size(LandCoverInputs%MapFile(:))
    if ( LandCoverInputs%MapFile(i) /= 'NOTSET' ) then
       call associate_File(LandCoverInputs%MapFile(i))
    end if
  end do
  call associate_File(LandCoverInputs%LandDefs)
  call associate_File(LandCoverInputs%Do3seDefs)

  call associate_File(femisFile)
  call associate_File(Vertical_levelsFile)
  call associate_File(EmisHeightsFile)
  call associate_File(SoilTypesFile)
  call associate_File(SurfacePressureFile)
  call associate_File(AircraftEmis_FLFile)
  call associate_File(nox_emission_1996_2005File)
  call associate_File(nox_emission_cams81File)
  call associate_File(MonthlyFacFile)
  call associate_File(DailyFacFile)
  call associate_File(HourlyFacFile)
  call associate_File(HourlyFacSpecialsFile)
  call associate_File(cmxbicDefaultFile)
  call associate_File(cmxBiomassBurning_FINNv1p5)
  call associate_File(cmxBiomassBurning_GFASv1)
  call associate_File(SplitDefaultFile)
  call associate_File(SplitSpecialsFile)
  call associate_File(RoadMapFile)
  call associate_File(AVG_SMI_2005_2010File)
  call associate_File(Soil_TegenFile)
  call associate_File(SitesFile)
  call associate_File(SondesFile)
  call associate_File(GLOBAL_LAInBVOCFile)
  call associate_File(EMEP_EuroBVOCFile)
  call associate_File(jclearFile)
  call associate_File(jcl1kmFile)
  call associate_File(jcl3kmFile)
  call associate_File(NdepFile)
  call associate_File(lightningFile)
  call associate_File(BiDirInputFile)  ! FUTURE INPUT
  call associate_File(LoganO3File)
  call associate_File(DustFile)
  call associate_File(TopoFile)
  call associate_File(Monthly_patternsFile)
  call associate_File(Monthly_timezoneFile)
  call associate_File(fileName_O3_Top)
  call associate_File(fileName_CH4_ibcs)
  call associate_File(NEST_template_read_3D)
  call associate_File(NEST_template_read_BC)
  call associate_File(NEST_template_write)
  call associate_File(NEST_MET_inner)
  call associate_File(filename_eta)

!DS
  OwnInputDir= key2str(OwnInputDir,'DataDir',DataDir)

  do i = 1, size(Emis_sourceFiles)
     !part of a class cannot be a target (?) must therefore do this separately
     if(Emis_sourceFiles(i)%filename/='NOTSET')then
        Emis_sourceFiles(i)%filename = key2str(Emis_sourceFiles(i)%filename,'DataDir',DataDir)
        Emis_sourceFiles(i)%filename = key2str(Emis_sourceFiles(i)%filename,'EmisDir',EmisDir)
        Emis_sourceFiles(i)%filename = key2str(Emis_sourceFiles(i)%filename,'GRID',GRID)
        Emis_sourceFiles(i)%filename = &
          key2str(Emis_sourceFiles(i)%filename,'OwnInputDir',OwnInputDir)
     endif
  enddo
  do i = 1, size(EmisMask)
     if(EmisMask(i)%filename/='NOTSET')then
        EmisMask(i)%filename = key2str(EmisMask(i)%filename,'DataDir',DataDir)
        EmisMask(i)%filename = key2str(EmisMask(i)%filename,'EmisDir',EmisDir)
        EmisMask(i)%filename = key2str(EmisMask(i)%filename,'GRID',GRID)
        EmisMask(i)%filename = &
          key2str(EmisMask(i)%filename,'OwnInputDir',OwnInputDir)
     endif
  enddo
  do i = 1,size(InputFiles)
    if(associated(InputFiles(i)%filename))then
     InputFiles(i)%filename =key2str(InputFiles(i)%filename,'ZCMDIR',ZCMDIR)
     InputFiles(i)%filename =key2str(InputFiles(i)%filename,'DataDir',DataDir)
     InputFiles(i)%filename =key2str(InputFiles(i)%filename,'EmisDir',EmisDir)
     InputFiles(i)%filename =key2str(InputFiles(i)%filename,'GRID',GRID)
     InputFiles(i)%filename = &
            key2str(InputFiles(i)%filename,'OwnInputDir',OwnInputDir)
     InputFiles(i)%filename =key2str(InputFiles(i)%filename,'YYYY',startdate(1))
    endif
  enddo
  if(trim(fileName_O3_Top)/="NOTSET")then
     fileName_O3_Top = key2str(fileName_O3_Top,'YYYY',startdate(1))
     if(MasterProc)then
        write(*,*)dtxt//'Reading 3 hourly O3 at top from :'
        write(*,*)trim(fileName_O3_Top)
     end if
  endif
  if(trim(fileName_CH4_ibcs)/="NOTSET" .and. MasterProc)then
     write(*,*)dtxt//'Reading CH4 IBCs from:', iyr_trend, trim(fileName_CH4_ibcs)
  endif

end subroutine Config_Constants

! PRELIM. Just writes out USES so far. 
subroutine WriteConfig_to_RunLog(iolog)
  integer, intent(in) :: iolog ! for Log file
  NAMELIST /OutUSES/ USES
  if(MasterProc)then
     write(iolog,*) ' USES after 1st time-step'
     write(iolog,nml=OutUSES)
  endif
end subroutine WriteConfig_to_RunLog

subroutine associate_File(FileName)
  integer, save::ix=0
  character(len=*), target ::FileName
  ix = ix+1
  call CheckStop(ix > size(InputFiles) , "Config_module: Size_InputFiles too small")
  InputFiles(ix)%filename => FileName  
end subroutine associate_File

end module Config_module
!_____________________________________________________________________________
!!TSTEMX program testr
!!TSTEMX use Config_module
!!TSTEXM !FAILS due to mpi call Config_Constants()
!!TSTEMX end program testr
