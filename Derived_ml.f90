
!==============================================================================
module Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module performs the calculations associated with "derived" 2D and 3D,
  ! such as accumulated precipitation or sulphate, daily, monthly or yearly 
  ! averages, depositions. These fields are all typically output as binary 
  ! fields.
  !
  ! The definitions of the derived fields shoul have been specified in the
  ! user-defined My_Derived_ml.
  !
  ! These modules combine Steffen's rewriting of the output routines
  ! (IOU notation), elements of chemint_mach, Bud_ml, etc.
  ! Re-coded to use only 2 types of data (d_2 and d_3)
  ! and F90 types for Deriv by ds, Sept. 2001. 

  ! User-defined routines and treatments are often needed here. So far I have 
  ! added stuff for VOC, AOTs, accsu, taken from the MACHO chemint code. In 
  ! general
  ! such code should be added in such a way that it isn't activated if not
  ! needed. It then doesn't need to be commented out if not used.
  ! For "accsu" I haven't found a way to do this though.
  ! See the examples given.
  ! ds, 28/9/01
  !---------------------------------------------------------------------------

!current definitions:
!SOMO35: daily max is found over 00:00 to 24:00. (not emepday)
!SOMO35: accumulated over one year
!D2_MAXO3 :  daily max is found over an EMEPDAY
!D2_MAXO3 : accumulated into yearly_output from April to September
!AOTXXc: accumulated into yearly_output from May to July
!AOTXXf: accumulated into yearly_output from April to September
!D2_EUAOTXXWH: accumulated into yearly_output from May to July
!D2_EUAOTXXDF: accumulated into yearly_output from April to September
!D2_UNAOTXXWH: accumulated into yearly_output from May to July
!D2_UNAOTXXDF: accumulated into yearly_output from April to September
!D2_O3 is now yearly accumulated (is this correct ?)
!D2_FSTXXXX: unit correct? (accumulated but unit in /s)

use My_Derived_ml  ! Definitions of derived fields, NWDEP, etc., f_wdep, etc.
use Chemfields_ml, only : xn_adv, xn_shl, cfac,xn_bgn, PM_water
use GenSpec_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use GenSpec_shl_ml
use GenSpec_tot_ml
use GenChemicals_ml, only : species
!ds rv1_9_28: big-fix, i_glob hasn't been defined yet.
!ds use GridValues_ml,   only : i_glob, j_glob !ds  rv1_9_16 for debug
use Met_ml, only :   roa,pzpbl,xksig,ps,th
use ModelConstants_ml, &
                   only: KMAX_MID &   ! =>  z dimension
                        , atwS, atwN, ATWAIR  &
                        , PPBINV  &   !   1.0e9, for conversion of units 
                        , PPTINV  &   !   1.0e12, for conversion of units 
                        , MFAC    &   ! converts roa (kg/m3 to M, molec/cm3)
                        , AOT_HORIZON&! limit of daylight for AOT calcs
                        ,DEBUG_i, DEBUG_j & !ds rv_9_16
                        , current_date
use Par_ml,    only: MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     me, NPROC,         &   ! for gc_abort checks
                     gi0,gj0,ISMBEG,JSMBEG,&! rv1_9_28 for i_glob, j_glob
                     limax, ljmax           ! => used x, y area 
use PhysicalConstants_ml,  only : PI
use Radiation_ml,  only :  zen


implicit none
private

 !ds rv1_9_17  public  :: SumDerived           ! old chemint
 public  :: Init_Derived         !ds rv1_9_16 new deriv system 
 public  :: ResetDerived         ! Resets values to zero
 public  :: DerivedProds         ! Calculates any production terms
 private :: Define_Derived        !ds ** NEW rv1_9_13
 private :: Setups 

!ds 16/12/2003 - find indices of used derived stuff in total list:

   public :: find_one_index
   public :: find_index

 private :: Consistency_checks   !ds  - checks index numbers from My_Derived
 private :: Consistency_count    !ds  - checks index numbers from My_Derived
 private :: Setup_VOC            ! Defines VOC group
 public :: Derived              ! Calculations of sums, avgs etc.
 private :: voc_2dcalc           ! Calculates sum of VOC for 2d fields
 private :: voc_3dcalc           ! Calculates sum of VOC for 3d fields


  !ds START OF NEW SYSTEM - added 16/12/2003 ========================

    type, public:: Deriv
       integer  :: code     ! Identifier for DNMI/xfelt
       character(len=9) :: class ! Type of data, e.g. ADV or VOC
       logical  :: avg      ! True => average data (divide by nav at end),
                            !     else accumulate over run period
       integer  :: index    ! index in concentation array, or other
       real     :: scale    ! Scaling factor
       logical  :: rho      ! True when scale is ug (N or S)
       logical  :: inst     ! True when instantaneous values needed
       logical  :: year     ! True when yearly averages wanted
       logical  :: month    ! True when monthly averages wanted
       logical  :: day      ! True when daily averages wanted
       character(len=15) :: name ! Name of the variable (for netCDF output)
       character(len=10) :: unit ! Unit (writen in netCDF output)
    end type Deriv

   logical, private, parameter :: T = .true., F = .false. ! shorthands only

   integer, public, parameter ::  &
       NDEF_WDEP = 4       & ! Number of 2D Wet deposition fields defined
      ,NDEF_DDEP = 21      & ! Number of 2D dry deposition fields defined
      ,NDEF_DERIV_2D = 66  & ! Number of 2D derived fields defined
      ,NDEF_DERIV_3D = 17  & ! Number of 3D derived fields defined
      ,NTDAY = 72            ! Number of 2D O3 to be saved each day (for SOMO)

   integer, public, dimension(NWDEP),     save :: nused_wdep
   integer, public, dimension(NDDEP),     save :: nused_ddep
   integer, public, dimension(NDERIV_2D), save :: nused_2d
   integer, public, dimension(NDERIV_3D), save :: nused_3d

  ! We put definitions of all possible variables in def_xx (def_ddep, etc.)
  ! and copy the needed ones into f_xx. The data will go into d_xx
  ! (d_2d, d_3d, etc.)

    type(Deriv),private, dimension(NDEF_WDEP),     save :: def_wdep !wet dep
    type(Deriv),private, dimension(NDEF_DDEP),     save :: def_ddep !dry dep
    type(Deriv),private, dimension(NDEF_DERIV_2D), save :: def_2d   !other 2D 
    type(Deriv),private, dimension(NDEF_DERIV_3D), save :: def_3d   !other 3D

    type(Deriv),public, dimension(NWDEP),     save :: f_wdep
    type(Deriv),public, dimension(NDDEP),     save :: f_ddep
    type(Deriv),public, dimension(NDERIV_2D), save :: f_2d
    type(Deriv),public, dimension(NDERIV_3D), save :: f_3d


  ! And we move our definitionms of actual fields used here (was My_Derived):

  ! Define first the 4 possible former output types
  ! corresponding to instantaneous,year,month,day

   integer, public, parameter ::  &
        IOU_INST=1, IOU_YEAR=2, IOU_MON=3, IOU_DAY=4, IOU_HOUR=5

  ! The 2-d and 3-d fields use the above as a time-dimension. We define
  ! LENOUTxD according to how fine resolution we want on output. For 2d
  ! fields we use daily outputs. For the big 3d fields, monthly output
  ! is sufficient.

   integer, public, parameter ::  LENOUT2D = 4  ! Allows INST..DAY for 2d fields
   integer, public, parameter ::  LENOUT3D = 4  ! Allows INST..MON for 3d fields

    real, save,  public :: &
      wdep( NWDEP    ,MAXLIMAX, MAXLJMAX, LENOUT2D), &  !wet dep
      ddep( NDDEP    ,MAXLIMAX, MAXLJMAX, LENOUT2D), &  !dry dep
      d_2d( NDERIV_2D,MAXLIMAX, MAXLJMAX, LENOUT2D), &  !other deriv
      d_3d( NDERIV_3D,MAXLIMAX, MAXLJMAX, KMAX_MID, LENOUT3D )

!pw to be used for SOMO35
    real, save,  public :: &   !save O3 every hour during one day to find running max
     D2_O3_DAY( MAXLIMAX, MAXLJMAX, NTDAY) = 0.
 

  !ds END OF NEW SYSTEM - added 16/12/2003 ========================


  ! Counters to keep track of averaging  (old nav, nav8to16, nav_8to16, etc
  ! Initialise here to zero.

    integer, public, dimension(NWDEP,LENOUT2D),     save :: nav_wdep = 0
    integer, public, dimension(NDDEP,LENOUT2D),     save :: nav_ddep = 0
    integer, public, dimension(NDERIV_2D,LENOUT2D), save :: nav_2d   = 0
    integer, public, dimension(NDERIV_3D,LENOUT3D), save :: nav_3d   = 0

   ! Note - previous versions did not have the LENOUT2D dimension
   ! for wet and dry deposition. Why not?  Are annual or daily
   ! depositions never printed? Since I prefer to keep all 2d
   ! fields as similar as posisble, I have kept this dimension
   ! for now - ds


   !-- some variables for the VOC sum done for ozone models
   !   (have no effect in non-ozone models - leave in code)

   integer, private, save :: nvoc   ! No. VOCs 
   integer, private, dimension(NSPEC_ADV), save :: &
             voc_index, &     ! Index of VOC in xn_adv
             voc_carbon       ! Number of C atoms

   logical, private, parameter :: MY_DEBUG = .false.
   logical, private, save :: debug_flag
   integer, private, save :: i_debug=1, j_debug=1  !ds rv1_9_28 Initialised, 
                                                   ! reset if DEBUG 
   integer, private :: i,j,k,n, ivoc, index    ! Local loop variables
   integer, public, parameter:: startmonth_forest=4,endmonth_forest=9&
                                ,startmonth_crops=5,endmonth_crops=7

   contains

    !=========================================================================
    subroutine Init_Derived()
    !dssubroutine SumDerived(dt)
      !ds logical, save :: my_first_call = .true.
      !ds real, intent(in) :: dt  !  time-step used in intergrations

      !ds if ( my_first_call ) then
          write(*,*) "INITIALISE My DERIVED STUFF"
          !ds call Set_My_Derived()
          call Define_Derived()
          call Consistency_checks()  !ds - checks index numbers from My_Derived
          call Setups()
          !ds my_first_call = .false. 
      !ds end if

    end subroutine Init_Derived

   !=========================================================================
    subroutine Define_Derived()

   ! Set the parameters for the derived parameters, including the codes
   ! used by DNMI/xfelt and scaling factors. (The scaling factors may
   ! be changed later in Derived_ml.
   ! And, Initialise the fields to zero.
   
    real, save    :: ugS = atwS*PPBINV/ATWAIR
    real, save    :: ugN = atwN*PPBINV/ATWAIR
    real, save    :: ugSO4, ugHCHO, ugCH3CHO
    real, save    :: ugPMad, ugPMde, ugSS    !st advected and derived PM's & SeaS


  !ds rv1_9_28  - for debug  - now not affecting ModelConstants version
   integer, dimension(MAXLIMAX) :: i_glob
   integer, dimension(MAXLJMAX) :: j_glob


    !   same mol.wt=100 assumed for PPM25 and PPMco

     ugPMad = species(PM25)%molwt * PPBINV /ATWAIR 
     ugPMde = PPBINV /ATWAIR
     ugSS  = species( SSfi )%molwt * PPBINV /ATWAIR  !SeaS

     ugSO4 = species( SO4 )%molwt * PPBINV /ATWAIR
     ugHCHO   = species ( HCHO )%molwt * PPBINV /ATWAIR
     ugCH3CHO = species ( CH3CHO )%molwt * PPBINV /ATWAIR

!-- Deposition fields. Define all possible fields and their xfelt codes here:

      !code class  avg? ind scale rho Inst Yr Mn Day   name      unit  

def_wdep = (/&
 Deriv( 561, "PREC ", F, -1, 1.0,   F  , F  ,T ,T ,T ,"WDEP_PREC","mm"),&
 Deriv( 541, "WDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"WDEP_SOX","mg/m2"),&
 Deriv( 542, "WDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"WDEP_OXN","mg/m2"),&
 Deriv( 543, "WDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"WDEP_RDN","mg/m2")&
 /)

    ! Dry dep. !
    !--includes fields for ecosystem specific--- use index numbers 800+ 
    ! ecosystem codes: SW = sea/water, CF = conif forest, DF = decid forest, 
    !                SN = seminatural (grass/moorlande/tundra)

      !code class  avg? ind scale rho Inst Yr Mn Day   name      unit  

def_ddep = (/&
 Deriv( 521, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"DDEP_SOX","mg/m2")&
,Deriv( 522, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"DDEP_OXN","mg/m2")&
,Deriv( 523, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"DDEP_RDN","mg/m2")&
,Deriv( 824, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSSW","mg/m2")&
,Deriv( 825, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSCF","mg/m2")&
,Deriv( 826, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSDF","mg/m2")&
,Deriv( 827, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSCR","mg/m2")&
,Deriv( 828, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSSN","mg/m2")&
,Deriv( 829, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXSWE","mg/m2")&
!&
,Deriv( 830, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNSW","mg/m2")&
,Deriv( 831, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNCF","mg/m2")&
,Deriv( 832, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNDF","mg/m2")&
,Deriv( 833, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNCR","mg/m2")&
,Deriv( 834, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNSN","mg/m2")&
,Deriv( 835, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_OXNWE","mg/m2")&
!&
,Deriv( 836, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNSW","mg/m2")&
,Deriv( 837, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNCF","mg/m2")&
,Deriv( 838, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNDF","mg/m2")&
,Deriv( 839, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNCR","mg/m2")&
,Deriv( 840, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNSN","mg/m2")&
,Deriv( 841, "DDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,F ,"DDEP_RDNWE","mg/m2")&
 /)

!-- 2-D fields - the complex ones
! (multiplied with roa in layers?? ==>  rho "false" ) !ds - explain!

!       code class  avg? ind scale rho  Inst  Yr  Mn   Day  name      unit 

def_2d = (/&
!ds Added AOT20 - in case the health people decide on this!
 Deriv( 870, "AOT  ", F, 20, 1.0,   F  , F  ,  T , T ,  F,"D2_AOT20","ppb h")&
,Deriv( 630, "AOT  ", F, 30, 1.0,   F  , F  ,  T , T ,  F,"D2_AOT30","ppb h")&
,Deriv( 608, "AOT  ", F, 40, 1.0,   F  , F  ,  T , T ,  F,"D2_AOT40","ppb h")&
,Deriv( 609, "AOT  ", F, 60, 1.0,   F  , F  ,  T , T ,  F,"D2_AOT60","ppb h")&
,Deriv( 680, "AOT  ", F, 30, 1.0,   F  , F  ,  T , F ,  F,"D2_AOT30f","ppb h")&
,Deriv( 681, "AOT  ", F, 40, 1.0,   F  , F  ,  T , F ,  F,"D2_AOT40f","ppb h")&
,Deriv( 682, "AOT  ", F, 60, 1.0,   F  , F  ,  T , F ,  F,"D2_AOT60f","ppb h")&
,Deriv( 683, "AOT  ", F, 40, 1.0,   F  , F  ,  T , F ,  F,"D2_AOT40c","ppb h")&
!BUGGY ,Deriv( 611, "ACCSU", T, -1, ugSO4,   F  , F  ,  T , T ,  F,"D2_ACCSU","ug/m2")&
!&
! -- simple advected species. Note that some indices need to be set as dummys
!    in ACID, e.g. IXADV_O3
!&
,Deriv( 601, "ADV  ", T, IXADV_SO2, ugS, T, F , T , T , T ,"D2_SO2","ugS/m3")&
,Deriv( 620, "ADV  ", T, IXADV_SO4, ugS, T, F , T , T , T ,"D2_SO4","ugS/m3")&
,Deriv( 621, "ADV  ", T, IXADV_HNO3,ugN, T, F , T , T , T ,"D2_HNO3","ugN/m3")&
,Deriv( 604, "ADV  ", T, IXADV_PAN, ugN, T, F , T , T , T ,"D2_PAN","ugN/m3")&
,Deriv( 622, "ADV  ", T, IXADV_NH3, ugN, T, F , T , T , T ,"D2_NH3","ugN/m3")&
,Deriv( 623, "ADV  ", T, IXADV_NO , ugN, T, F , T , T , T ,"D2_NO","ugN/m3")&
,Deriv( 606, "ADV  ", T, IXADV_NO2, ugN, T, F , T , T , T ,"D2_NO2","ugN/m3")&
,Deriv( 600, "ADV  ", T,IXADV_aNH4, ugN, T, F , T , T , T ,"D2_aNH4","ugN/m3")&
!&
,Deriv( 607, "ADV  ",T,IXADV_O3 ,PPBINV, F, F , T, T , T ,"D2_O3","ppb")&
,Deriv( 612, "ADV  ",T,IXADV_CO ,PPBINV, F, F , T, T , T ,"D2_CO","ppb")&
,Deriv( 670, "ADV  ",T,IXADV_aNO3, ugN,  T, F , T, T , T ,"D2_aNO3","ugN/m3")&
,Deriv( 671, "ADV ", T,IXADV_pNO3, ugN,  T, F , T, T , T ,"D2_pNO3", "ugN/m3")&
!ds 31/3/2004 new:
,Deriv( 672,"NOX  ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_NOX","ugN/m3")&
,Deriv( 673,"NOZ  ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_NOZ","ugN/m3")&
,Deriv( 674,"OX   ", T,   -1  ,PPBINV , F , F,T,T,T,"D2_OX","ppb")&
!&
,Deriv( 615, "ADV  ",T,IXADV_PM25, ugPMad, T, F , T, T, T,"D2_PPM25","ug/m3")&
,Deriv( 616, "ADV  ",T,IXADV_PMco, ugPMad, T, F , T, T, T,"D2_PPMco","ug/m3")&
!SeaS
,Deriv( 659, "ADV  ",T,IXADV_SSfi, ugSS, T, F , T, T, T,"D2_SSfi","ug/m3")&
,Deriv( 660, "ADV  ",T,IXADV_SSco, ugSS, T, F , T, T, T,"D2_SSco","ug/m3")&
,Deriv(   8, "PS    ",T,  0 ,       1.0, F , T, T, T, T ,"D2_PS","m")&
,Deriv( 468, "HMIX  ",T,  0 ,       1.0, T , F, T, T, T ,"D2_HMIX","m")&
,Deriv( 469, "HMIX00",T,  0 ,       1.0, T , F, T, T, T ,"D2_HMIX00","m")&
,Deriv( 470, "HMIX12",T,  0 ,       1.0, T , F, T, T, T ,"D2_HMIX12","m")&
!&
!ds drydep
!   set as "external" parameters - ie set outside Derived subroutine
!   use index as lu, here 10=grass
!      code class   avg? ind scale rho  Inst Yr  Mn   Day    name      unit 
!
!Havn't worried about rho so far... does it matter?
!ds rv1_6_15 redef:
,Deriv( 850, "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_FSTDF00","nmol/m2/s")&
,Deriv( 851, "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_FSTDF08","nmol/m2/s")&
,Deriv( 852, "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_FSTDF16","nmol/m2/s")&
,Deriv( 853, "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_FSTWH00","nmol/m2/s")&
,Deriv( 854, "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_FSTWH20","nmol/m2/s")&
,Deriv( 855, "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_FSTWH30","nmol/m2/s")& !New IAM version
,Deriv( 856, "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_FSTWH40","nmol/m2/s")&
,Deriv( 857, "EXT  ", F, -1, 1. , F, F,T ,T ,T ,"D2_FSTWH60","nmol/m2/s")&

!      code class   avg? ind scale rho Inst Yr Mn  Day   name      unit 
,Deriv( 860, "EXT  ", T, -1, 1.   , F, F,T ,T ,T ,"D2_O3DF   ","ppb")&
,Deriv( 861, "EXT  ", T, -1, 1.   , F, F,T ,T ,T ,"D2_O3WH   ","ppb")&
!
! AOT30 and AOT40 for Wheart and Beech. May need daily here.
! Also, use field for EU definition (8-20 CET) and Mapping Manual/UNECE
! (daylight hours). 
! All of these use O3 at crop height, in contrast to the older AOT30, AOT40
! as defined above, and all allow daily output.
,Deriv( 862, "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_EUAOT30WH","ppb h")&
,Deriv( 863, "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_EUAOT40WH","ppb h")&
,Deriv( 864, "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_EUAOT30DF","ppb h")&
,Deriv( 865, "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_EUAOT40DF","ppb h")&
! UNECE:
,Deriv( 866, "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_UNAOT30WH","ppb h")&
,Deriv( 867, "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_UNAOT40WH","ppb h")&
,Deriv( 868, "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_UNAOT30DF","ppb h")&
,Deriv( 869, "EXT  ", F, -1, 1.   , F, F,T ,T ,T ,"D2_UNAOT40DF","ppb h")&
!
! --  time-averages - here 8-16
!
,Deriv( 613,"TADV ", T,IXADV_HCHO  ,ugHCHO,  T, F, T, T, T,"D2T_HCHO","ug/m3")&
,Deriv( 614,"TADV ", T,IXADV_CH3CHO,ugCH3CHO,T, F, T, T,T,"D2T_CH3CHO","ug/m3")&
,Deriv( 610,"VOC  ", T,  -1    ,PPBINV, F, F, T, T, T,"D2_VOC","ppb")&
!
! -- miscellaneous user-defined functions
!
!      code class   avg? ind scale rho Inst Yr Mn  Day   name      unit 
!ds ,Deriv( 602,"TSO4 ", T,   -1  ,ugS ,    T , F,T,T,T,"D2_SOX","ugS/m3")&
,Deriv( 603,"TOXN ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_OXN","ugN/m3")&
,Deriv( 605,"TRDN ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_REDN","ugN/m3")&
,Deriv( 624,"FRNIT", T,   -1  ,1.0 ,    F , F,T,T,T,"D2_FRNIT","(1)")&
,Deriv( 625,"MAXADV", F,IXADV_O3,PPBINV, F, F,T,T,T,"D2_MAXO3","ppb")&
,Deriv( 626,"MAXSHL", F,IXSHL_OH,1.0e13,F , F,T,F,T,"D2_MAXOH","?")&
!
,Deriv( 617, "tNO3 ", T, -1, ugN,    T, F, T, T, T,"D2_tNO3", "ugN/m3")&
,Deriv( 618, "SIA  ", T, -1, ugPMde, T, F, T, T, T,"D2_SIA" , "ug/m3")&
,Deriv( 619, "PMco ", T, -1, ugPMde, T, F, T, T, T,"D2_PMco", "ug/m3")&
,Deriv( 648, "PM25 ", T, -1, ugPMde, T, F, T, T, T,"D2_PM25", "ug/m3")&
,Deriv( 649, "PM10 ", T, -1, ugPMde, T, F, T, T, T,"D2_PM10", "ug/m3")&
,Deriv( 662, "H2O  ", T, -1,   1.0 , T, F, T, T, T,"D2_PM25_H2O ", "ug/m3")&   !water
,Deriv( 646, "SSalt", T, -1, ugSS,   T, F, T, T, T,"D2_SS  ", "ug/m3")& 
,Deriv( 627, "SOM", F, 35, 1.,   F, F, T, T, F,"D2_SOMO35", "ppb day")& 
,Deriv( 628, "SOM", F,  0, 1.,   F, F, T, T, F,"D2_SOMO0", "ppb day")& 
 /)

!-- 3-D fields

def_3d = (/ &
 Deriv(  18, "TH  ",T,  0 ,       1.0, F , T, T, T, T ,"D3_TH","m")&

, Deriv( 401, "ADV  ", T, IXADV_O3 , PPBINV, F, T, T, T, F ,"D3_O3","ppb")&
,Deriv( 402, "ADV  ", T, IXADV_SO2, PPBINV, F, T, T, T, F ,"D3_SO2","ppb")&
,Deriv( 403, "ADV  ", T, IXADV_PAN, PPBINV, F, T, T, T, F ,"D3_PAN","ppb")&
,Deriv( 404, "ADV  ", T, IXADV_HNO3,PPBINV, F, T, T, T, F ,"D3_HNO3","ppb")&
,Deriv( 405, "ADV  ", T, IXADV_aNO3,PPBINV, F, T, T, T, F ,"D3_aNO3","ppb")&
,Deriv( 406, "ADV  ", T, IXADV_NO2, PPBINV, F, T, T, T, F ,"D3_NO2","ppb")&
,Deriv( 407, "VOC  ", T,       -1 , PPBINV, F, T, T, T, F ,"D3_VOC","ppb")&
,Deriv( 408, "ADV  ", T, IXADV_aNH4,PPBINV, F, T, T, T, F ,"D3_aNH4","ppb")&
,Deriv( 409, "ADV  ", T, IXADV_SO4, PPBINV, F, T, T, T, F ,"D3_SO4","ppb")&
,Deriv( 410, "ADV  ", T, IXADV_H2O2,PPBINV, F, T, T, T, F ,"D3_H2O2","ppb")&
!&
!& rv1_9_28:
! Set Year true to allow debug - change later
,Deriv( 411, "SHL",   T, IXSHL_OH,  PPTINV, T, F, T, T, F ,"D3_OH","?")& ! ds from 1.0E16
,Deriv( 412, "ADV",   T, IXADV_CH3COO2,&
                                    PPTINV, F, F, T, T, F ,"D3_CH3COO2","?")&
,Deriv( 413, "MAX3DSHL", T,IXSHL_OH,PPTINV, T, F, T, T, F ,"D3_MAXOH","?")&   ! rho true for shl
,Deriv( 414, "MAX3DADV", T, IXADV_CH3COO2,&
                                    PPTINV, F, F, T, T, F ,"D3_MAXCH3COO2","?")&
,Deriv( 415, "PHNO3   ", T, IXSHL_PHNO3,1.0e8, F, F, T, T, F ,"D3_PHNO3","?")&
,Deriv( 416, "MAX3DADV", T, IXADV_O3,PPBINV,F, F, T, T, F ,"D3_MAXO3","?")&
/)

   !ds new, 25/3/2004:

     if ( SOURCE_RECEPTOR ) then  ! We assume that no daily outputs are wanted.
        def_2d(:)%day = .false.
        def_ddep(:)%day = .false.
        def_wdep(:)%day = .false.
     end if

   !Initialise to zero
   ! was previously in subroutine Set_My_Derived()

      wdep( :,:,:,:) = 0.0
      ddep( :,:,:,:) = 0.0
      d_2d( :,:,:,:) = 0.0
      d_3d( :,:,:,:,:) = 0.0


     ! Get indices of wanted fields in larger def_xx arrays:

      call find_index( DDEP_USED(:), def_ddep(:)%name, nused_ddep  ) 
      call find_index( WDEP_USED(:), def_wdep(:)%name, nused_wdep  ) 
      call find_index( D2_USED(:),   def_2d(:)%name, nused_2d  ) 

      if ( NDERIV_3D > 0 ) then  ! Special because of DUMMY possibility

         call find_index( D3_USED(1:NDERIV_3D),   def_3d(:)%name, nused_3d  ) 

        !ds was some strange bevaious for another array here!
         !do i = 1, 3
         ! print *, "AFTERNUSED ", i, nused_ddep(i), size(nused_ddep)
         !end do
      end if

     ! And set f_xx fields  from these indices and the fuller def_xx fields:

      f_wdep(:) = def_wdep( nused_wdep(:) )
      f_ddep(:) = def_ddep( nused_ddep(:) )
      f_2d(:)   = def_2d( nused_2d(:) )

      do n = 1, NDERIV_3D ! Special because of DUMMY possibility
         if( nused_3d(n) > 0) f_3d(n)   = def_3d( nused_3d(n) )
      end do

      debug_flag = .false.    !ds rv1_9_28
      if ( MY_DEBUG ) then

          !ds rv1_9_28: Need tod efine here since gridValues not yet set.

          i_glob = (/ (n + gi0 + ISMBEG - 2, n=1,MAXLIMAX) /)
          j_glob = (/ (n + gj0 + JSMBEG - 2, n=1,MAXLJMAX) /)

           do j = 1, ljmax
              do i = 1, limax
                  if ( i_glob(i)==DEBUG_i .and. j_glob(j)==DEBUG_j) then
                       debug_flag = .true.
                       i_debug = i
                       j_debug = j
                       write(*,*) "Derived DEBUG COORDS ", me, i_debug, j_debug
                  end if
              end do
            end do
      end if ! DEBUG

  end subroutine Define_Derived
 !=========================================================================

  function find_one_index(used, defined ) 
    character(len=*), intent(in) :: used
    character(len=*), dimension(:), intent(in) :: defined
    integer :: find_one_index
    integer :: n_match ! Count for safety
    integer :: nu, nf

    n_match  = 0

    if( used == "DUMMY") then
        write(6,*) "Dummy index found"
        find_one_index = -1
        return
    end if

    do nf = 1, size(defined)

         if ( used == defined(nf)  ) then
            find_one_index = nf
            n_match = n_match + 1
         end if
    end do ! nf

    if ( n_match /= 1 ) then
          print *, "ERROR: Find_one_index for:", used,  "Matches ", n_match
          call gc_abort(me,NPROC,"Find_one_index - no match!!")
    else
          if ( me == 0 ) then
             write(6,*) "FOUND_ONE_INDEX ", used, " => NF ", find_one_index
          end if
    end if ! ERROR check

  end function find_one_index

 !=========================================================================
  subroutine find_index(used, defined, nused) 

    !+
    ! Searches for the "used" characters in the list of defined species
    ! e.g. for contents of DDEP_USED in f_ddep%name
    ! Returns an integer array of these indices
    !+
    ! ds, 16/12/2003

    character(len=*), dimension(:), intent(in) :: used
    character(len=*), dimension(:), intent(in) :: defined
    integer, dimension(size(used)), intent(out) :: nused
    integer, dimension(size(used)) :: n_match ! Count for safety
    integer :: nu, nf

    n_match  = 0
    do nu = 1, size(used)

      if( used(nu) == "DUMMY") then
          write(6,*) "Dummy index found for nu", nu
          nused(nu) = -1
          n_match(nu) = n_match(nu) + 1
      end if

      do nf = 1, size(defined)
         if ( used(nu) == defined(nf) ) then
            nused(nu) = nf
            n_match(nu) = n_match(nu) + 1
         end if
      end do ! nf

      if ( n_match(nu) /= 1 ) then
          print *, "ERROR: Find_index2 Size U ", size(used), " DEF ", size(defined)
          print *, "ERROR: Find_index2 n", nu, " V ", used(nu), "Matches ", n_match(nu)
          call gc_abort(me,NPROC,"Find_index2 - no match!!")
      else
          if ( me == 0 ) then
             write(6,*) "FOUND_INDEX ", nu, used(nu), " => NF ", nused(nu)
          end if
      end if ! ERROR check
    end do ! nu

  end subroutine find_index

 !=========================================================================
  subroutine Consistency_checks()  

    !/** ds - checks index numbers from My_Derived to look for duplicates

       integer, parameter :: MAX_INDEX = 2000
       integer, dimension(MAX_INDEX) :: index_used
       integer :: i

       index_used = 0
       call  Consistency_count(MAX_INDEX, index_used, NWDEP, f_wdep) 
       call  Consistency_count(MAX_INDEX, index_used, NDDEP, f_ddep)
       call  Consistency_count(MAX_INDEX, index_used, NDERIV_2D, f_2d) 
       call  Consistency_count(MAX_INDEX, index_used, NDERIV_3D, f_3d) 

       if ( any( index_used > 1 ) ) then
           do i = 1, MAX_INDEX
             if( index_used(i) > 1 ) print *,  &
                   "Derived code problem, index: ",i, index_used(i)
           end do
           call gc_abort(me,NPROC,"Derived code problem!!")
       end if

     !ds 25/3/2004:

      if ((      SOURCE_RECEPTOR .and.  NDERIV_2D /= NSR_2D )   .or.  &
          ( .not.SOURCE_RECEPTOR .and.  NDERIV_2D == NSR_2D )) then
           print *, "Error! Wrong number of 2D params ", NDERIV_2D, NSR_2D
           call gc_abort(me,NPROC,"Derived code problem!!")
      end if
     end subroutine

    !=========================================================================
     subroutine Consistency_count(max,index_used, n,data)

    !/** ds adds up number of times each code from derived data array
    !    is used.

      integer, intent(in) :: max, n
      integer, dimension(:), intent(inout)  :: index_used
      type(Deriv), dimension(n), intent(in) :: data

      integer :: code, i

        do i = 1, n
           code = data(i)%code
           if ( code > max ) call gc_abort(me,NPROC,"My_Derived code >max!")
           index_used( code )  = index_used( code )  + 1
           !if( me == 0 ) print "(a6,i3,i4,2a15,i2)", "CODE ", i, code, 
           !     data(i)%name, data(i)%class, index_used(code)

           if( index_used(i) > 1 ) then
                print *, "Derived code problem, Consistency checks: me ",  me, &
                       " index: ",i, index_used(code)
                call gc_abort(me,NPROC,"Consistency code problem!!")
           end if

        end do
     end subroutine Consistency_count
    
    !=========================================================================
     subroutine Setups()

    !/** flexibility note. By making use of character-based tests such
    !    as for "VOC" below, we achieve code which can stay for both MADE and
    !    MACHO without having to define non-used indices. 
    !    Similarly, we avoid the previous "if NUM_ACCSU eq 1" type test,
    !    since the appropriate code will now only activate 

    !/ ** if voc wanted, set up voc_array. Works for all ozone chemistries
    !     (and code not called for MADE-type).

      if ( any(  f_2d(:)%class == "VOC" ) .or. &
           any(  f_3d(:)%class == "VOC" )  ) then
            call Setup_VOC()
            write(6,*) "Derived VOC setup returns ", nvoc, "vocs"
            write(6,"(a12,/,(20i3))")  "indices ", voc_index(1:nvoc)
            write(6,"(a12,/,(20i3))")  "carbons ", voc_carbon(1:nvoc)
      end if


    end subroutine Setups
    !=========================================================================

    subroutine Derived(dt,End_of_Day)   !ds rv1_9_28: End_of_Day added

    !/** DESCRIPTION
    !  Integration and averaging of chemical fields. Intended to be
    !  a more flexible version of the old chemint routine.
    !  Includes AOT40, AOT60 if present

      real, intent(in)    :: dt           !  time-step used in intergrations
      logical, intent(in) :: End_of_Day   !ds rv1_9_28: End_of_Day added

      character(len=len(f_2d%class)) :: typ  !  See defs of f_2d
      real :: thour                          ! Time of day (GMT)
      real :: timefrac                       ! dt as fraction of hour (3600/dt)
      real :: dayfrac                        ! fraction of day elapsed (in middle of dt)
      integer :: ntime                       ! 1...NTDAYS
      integer :: nhour                       ! hours of day (GMT) 
      real, dimension(MAXLIMAX,MAXLJMAX) :: density !  roa (kgair m-3 when 
                                                    ! scale in ug,  else 1

      real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: inv_air_density3D !ds rv1_9_28 
                ! Inverse of No. air mols/cm3 = 1/M 
                ! where M =  roa (kgair m-3) * MFAC  when ! scale in ug,  else 1
      logical :: accumulate_2dyear !flag to know when to accumulate d_2d (case "EXT")

      !ds rv1_9_17 integer :: ndef 

      timefrac = dt/3600.0
      thour = current_date%hour+current_date%seconds/3600.0  !ds rv1_9_28 moved here for 3D3D


     !/***** 2-D fields **************************

     do n = 1, NDERIV_2D

        accumulate_2dyear=.true.
        !rv1_9_17 ndef = nused_2d(n)
        typ = f_2d(n)%class


        if ( f_2d(n)%rho ) then
            forall ( i=1:limax, j=1:ljmax )
                density(i,j) = roa(i,j,KMAX_MID,1)
            end forall
        else
            density(:,:) = 1.0
        end if

        !/** user-defined time-averaging. Here we have defined TADV and TVOC
        !    so that 8-hour daytime averages will be calculated. 
        !    Just comment out if not wanted, or (better!) don't define any
        !    f_2d as TADV or TVOC

        if ( typ == "TADV" .or. typ == "TVOC" ) then
             !ds thour = current_date%hour+current_date%seconds/3600.0
             if(thour <= 8.0 .or. thour > 16.0 ) cycle  ! Start next species
        end if

       !hf hmix average at 00 and 12:

        if ( typ == "HMIX00" .or. typ == "XKSIG00" ) then
             !ds thour = current_date%hour+current_date%seconds/3600.0
             if(thour /= 0.0 ) cycle  ! Start next species
        end if

        if ( typ == "HMIX12" .or. typ == "XKSIG12" ) then
             !ds thour = current_date%hour+current_date%seconds/3600.0
             if(thour /= 12.0 ) cycle  ! Start next species
        end if

        !rv1_9_17 index = f_2d(ndef)%index
        index = f_2d(n)%index
        select case ( typ )

          case ( "PS" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = ps(i,j,1)*0.01
            end forall


          case ( "HMIX", "HMIX00", "HMIX12" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = pzpbl(i,j)  
            end forall

            if ( debug_flag ) then
             write(*,fmt="(a12,2i4,4f12.3)") "HMIX" , n , d_2d(n,i_debug,j_debug,IOU_INST)       
            end if

         ! Simple advected species:
          case ( "ADV", "TADV" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = xn_adv(index,i,j,KMAX_MID)  &
                                     * cfac(index,i,j) * density(i,j)  
            end forall

            if ( debug_flag ) then
             write(*,fmt="(a12,2i4,4f12.3)") "JUST ADV" , n, index  &
              ,d_2d(n,i_debug,j_debug,IOU_INST)*PPBINV &
              ,xn_adv(index,i_debug,j_debug,KMAX_MID)*PPBINV &
              ,density(i_debug,j_debug), cfac(index,i_debug,j_debug)
            end if

!water
          case ( "H2O" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = PM_water(i,j,KMAX_MID)  
            end forall
  

          case ( "MAXADV" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_DAY) = max( d_2d( n, i,j,IOU_DAY), &
                                xn_adv(index,i,j,KMAX_MID)  &
                                     * cfac(index,i,j) * density(i,j) )
            end forall

            if ( debug_flag ) then
             write(*,fmt="(a12,2i4,4f12.3)") "ADV MAX. ", n, index  &
                      , d_2d(n,i_debug,j_debug,IOU_DAY) * PPBINV      &
                      ,  xn_adv(index,i_debug,j_debug,KMAX_MID)* PPBINV  &
                      ,  density(i_debug,j_debug), cfac(index,i_debug,j_debug)

            end if

            !pw rv2_1
            !Monthly and yearly ARE averaged over days
            if(End_of_Day)then
              d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)
              nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
              if(    current_date%month >= 4 &
                 .or.current_date%month <= 9 )then
              d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY) 
              nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
              endif
            endif


          case ( "MAXSHL" )        ! Daily maxima - short-lived

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_DAY) = max( d_2d( n, i,j,IOU_DAY), &
                                xn_shl(index,i,j,KMAX_MID)  &
                                    / (density(i,j)*MFAC) )
                                   !u4  / (roa(:,:,KMAX_MID,1)*MFAC) )
            end forall


            if ( debug_flag ) then
               write(*, *) "SHL:MAX.,MFAC ", n, index  , MFAC
               write(*,fmt="(a12,2i4,4es12.3)") "SHL MAX. ", n, index  &
                      , d_2d(n,i_debug,j_debug,IOU_DAY) &
                      ,  xn_shl(index,i_debug,j_debug,KMAX_MID)  &
                      ,  density(i_debug,j_debug), MFAC
            end if

            !pw rv2_1
            !Monthly and yearly ARE averaged over days
            if(End_of_Day)then
              d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)
              nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
              if(    current_date%month >= 4 &
                 .or.current_date%month <= 9 )then
              d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY) 
              nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
              endif
            endif

          case ( "VOC", "TVOC" )

            call voc_2dcalc()

          case( "AOT" )

            call aot_calc( n, timefrac )

            if(     current_date%month<startmonth_forest&
                 .or.current_date%month>endmonth_forest)then
               if( f_2d(n)%name=="D2_AOT30f".or.& 
                   f_2d(n)%name=="D2_AOT40f".or.&
                   f_2d(n)%name=="D2_AOT60f")then
                   accumulate_2dyear=.false.
               endif

            endif
            if(     current_date%month<startmonth_crops&
               .or.current_date%month>endmonth_crops)then
               if( f_2d(n)%name=="D2_AOT30c".or.&
                   f_2d(n)%name=="D2_AOT40c".or.&
                   f_2d(n)%name=="D2_AOT60c")then
                   accumulate_2dyear=.false.
               endif
            endif


           case( "SOM" )


              dayfrac= (thour-(dt/7200.))/24. !must be < 1
              ntime=int(dayfrac*NTDAY )+1 !must be >=1 and <= NTDAY
              if(dayfrac<0)ntime=NTDAY !midnight

!last value  (not averaged): 
          D2_O3_DAY( : , : , ntime) =&
           xn_adv(IXADV_O3,:,:,KMAX_MID)*cfac(IXADV_O3,:,:)*PPBINV

              if(dayfrac<0)then !only at midnight: write on d_2d

                 
                 call som_calc( n )
!accumulate
                 d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY) 
 !                if(current_date%month>=4.and.current_date%month<=9)then
                 d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY) 
 !                endif
!overwritten anyway D2_O3_DAY = 0.
              endif



          case ( "EXT" )

             call setaccumulate_2dyear(n,accumulate_2dyear)
          ! Externally set for IOU_INST (in other routines); so no new work needed
            if ( debug_flag ) then
                 write(*,*) "EXTDer:Externally set d_2d should already have values"
                 write(*,*) "EXTDer: d_2d(i_debug,j_debug) for", typ, " is ", &
                           d_2d(n,i_debug,j_debug,IOU_INST)
            end if

          case  default

            if ( debug_flag ) then
                 write(*,*) "My_Deriv called for n=", n, "Type ",typ
                 write(*,*) "My_Deriv index", index
                 write(*,*) "My_Deriv avg? ", f_2d(n)%avg,  &
                                   " for nav ", nav_2d(n,IOU_INST)
                 write(*,*) "Deriv: Length of f_2d is ", len(f_2d%class)
                 write(*,*) "Deriv: f_2d class is ", f_2d(n)%class
             end if 

             call My_DerivFunc( d_2d(n,:,:,IOU_INST), n, typ, timefrac, density ) 

        end select


        !/** add to daily, monthly and yearly average, and increment counters
        !    /wdep, ddep not done here ???)
        !  Note that the MAXADV and MAXSHL and SOM needn't be summed here, but
        !  since the INST values are zero it doesn't harm, and the code is 
        !  shorter. These d_2d ( MAXADV, MAXSHL, SOM) are set elsewhere

        d_2d(n,:,:,IOU_DAY )  = d_2d(n,:,:,IOU_DAY )  + d_2d(n,:,:,IOU_INST) 
        if ( f_2d(n)%avg ) nav_2d(n,IOU_DAY) = nav_2d(n,IOU_DAY) + 1
        d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_INST) 
        if ( f_2d(n)%avg ) nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
        if(accumulate_2dyear)then
           d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_INST) 
           if ( f_2d(n)%avg ) nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
        endif
 
     end do   ! NDERIV_2D
     !/***** WET DEPOSITION **************************

       do n = 1, NWDEP  ! includes precip from Met_ml
        wdep(n,:,:,IOU_DAY )  = wdep(n,:,:,IOU_DAY )  + wdep(n,:,:,IOU_INST) 
        wdep(n,:,:,IOU_MON )  = wdep(n,:,:,IOU_MON )  + wdep(n,:,:,IOU_INST) 
        wdep(n,:,:,IOU_YEAR ) = wdep(n,:,:,IOU_YEAR ) + wdep(n,:,:,IOU_INST) 

        !ds if ( MY_DEBUG .and. me == 0 ) print *, "wet dep: ",wdep(n,3,3,IOU_DAY)

     end do  ! WET DEP.
     !/***** WET DEPOSITION **************************

     do n = 1, NDDEP

        ddep(n,:,:,IOU_DAY )  = ddep(n,:,:,IOU_DAY )  + ddep(n,:,:,IOU_INST) 
        ddep(n,:,:,IOU_MON )  = ddep(n,:,:,IOU_MON )  + ddep(n,:,:,IOU_INST) 
        ddep(n,:,:,IOU_YEAR ) = ddep(n,:,:,IOU_YEAR ) + ddep(n,:,:,IOU_INST) 

        !ds if ( MY_DEBUG .and. me == 0 ) print *, "dry dep: ",ddep(n,3,3,IOU_DAY)

     end do  ! DRY DEP.

     !/***** 3-D fields **************************

       if(debug_flag) then ! RUN through indices etc.
            write(*, "(a12,2i4,f12.3)") "3D3D TIME ",  me, NDERIV_3D, &
                     (current_date%hour+current_date%seconds/3600.0)
            !do n = 1, NDERIV_3D
            ! write(6,*) "3D3D CHECKING me, n,name,index,class ", me, &
            !        n, f_3d(n)%name,f_3d(n)%index, f_3d(n)%class
            !end do
        end if


     do n = 1, NDERIV_3D

        index = f_3d(n)%index

       if ( f_3d(n)%rho ) then
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                inv_air_density3D(i,j,k) = 1.0/( roa(i,j,k,1) * MFAC )
            end forall
        else
            !ds BUG air_density3D(i,j,k) = 1.0
            inv_air_density3D(:,:,:) = 1.0
        end if

        select case ( f_3d(n)%class )

         ! Simple advected species:
          case ( "ADV" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
            end forall

         case ( "BGN" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_bgn(index,i,j,k)
            end forall

         case ("XKSIG00", "XKSIG12" ) !hf hmix xksig

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xksig(i,j,k)
            end forall

         case ("TH  " ) !JEJ Pot. temp (needed for cross sections)

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = th(i,j,k,1)
            end forall


     ! ds rv1_9_28 changes -----------------------------------

         case ( "PHNO3" )   !ds-hf  rv1_9_28
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_shl(index,i,j,k)
            end forall

             if(debug_flag) write(*,"(a12,i4,2es12.3)") "3D3D PHNO3", n, &
                  xn_shl(index,i_debug,j_debug,KMAX_MID), &
                  d_3d(n,i_debug,j_debug,KMAX_MID,IOU_INST)

         case ( "MAX3DSHL" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )! Daily maxima - short-lived
              d_3d( n, i,j,k,IOU_INST) = max( d_3d( n, i,j,k,IOU_INST),&
                                      xn_shl(index,i,j,k) &
                                     * inv_air_density3D(i,j,k) )
            end forall

            if(debug_flag) write(*,"(a13,i4,f8.3,3es12.3)") "3D3D MAX3DSHL", n, thour, &
              xn_shl(index,i_debug,j_debug,KMAX_MID), &
              1.0/inv_air_density3D(i_debug,j_debug,KMAX_MID), &
              d_3d(n,i_debug,j_debug,KMAX_MID,IOU_INST)

          case ( "MAX3DADV" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) =  max( d_3d( n, i,j,k,IOU_INST),&
                                               xn_adv(index,i,j,k) )
            end forall

             if(debug_flag) write(*,"(a12,i4,f8.3,4es12.3)") "SET MAX3DADV", n, thour, &
                      xn_adv(index,i_debug,j_debug,KMAX_MID), &
                      d_3d(n,i_debug,j_debug,KMAX_MID,IOU_INST)

          case ( "SHL" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) =   xn_shl(index,i,j,k) * inv_air_density3D(i,j,k)
            end forall



          case ( "VOC" )

            call voc_3dcalc()

          case  default

           print *, "Derived 3D class NOT FOUND", n, index, &
                         f_3d(n)%name,f_3d(n)%class
           call gc_abort(me,NPROC,"3D FIELD - no match!!")



        end select
     
       !ds v1_9_28:


      !/** add to monthly and yearly average, and increment counters
       !    ( no daily averaging done for 3-D fields so far).


       !ds for the MAX3D possibilities, we store maximum value of the
       !   current day in the IOU_INST variables.
       !   These are then added into IOU_MON **only** at the end of each day.
       ! (NB there is an error made on 1st day used, since only 1st 6 hours
       !  are checked. Still, not much happens on 1st Jan.... ;-)

        if ( (f_3d(n)%class == "MAX3DSHL")  .or. &
             (f_3d(n)%class == "MAX3DADV") )then
           if (End_of_Day) then
              d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
                                     + d_3d(n,:,:,:,IOU_INST)
              d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
                                     + d_3d(n,:,:,:,IOU_INST) !men er d_3d reset paa samme tidspunkt?
              if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1 !only collected for end_of_day

              if( debug_flag ) then
                    write(*,fmt="(a20,a9,i4,f8.3,2es12.3)") "END_OF_DAY MAX3D", &
                      f_3d(n)%class, n, thour,  &
                      d_3d(n,i_debug,j_debug,KMAX_MID,IOU_MON ),&
                      d_3d(n,i_debug,j_debug,KMAX_MID,IOU_INST )
                    write(*,"(a20,i4,2x,6i6)") "END_OF_DAY NAV ", &
                      n, (nav_3d(n,i), i=1,LENOUT3D)
              end if

              d_3d(n,:,:,:,IOU_INST ) = 0.0  !! Reset d_3d  ! Was bug in Unimod.debug

           endif ! End_of_Day
        else
           d_3d(n,:,:,:,IOU_DAY ) = d_3d(n,:,:,:,IOU_DAY ) &
                + d_3d(n,:,:,:,IOU_INST)
           d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
                + d_3d(n,:,:,:,IOU_INST)
           d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
                + d_3d(n,:,:,:,IOU_INST)
           if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1
        endif



!      !/** add to monthly and yearly average, and increment counters
!      !    ( no daily averaging done for 3-D fields so far).
!
!       d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
!                              + d_3d(n,:,:,:,IOU_INST)
!       d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
!                              + d_3d(n,:,:,:,IOU_INST)
!
!       if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1
 
      end do
    end subroutine Derived
    !=========================================================================

    subroutine DerivedProds(text,dt)

    !/** DESCRIPTION
    !  Calculates chemical changes by comparing values before and  after 
    !  chemistry subroutine. Intended to be a more flexible version of the old 
    !  PRODO3  calculation

      character(len=*), intent(in) :: text  ! "Before" or "After"
      real,             intent(in) :: dt    ! timestep (s)

      real :: timefrac                      ! dt as fraction of hour (3600/dt)



      if (.not. any( f_3d%class == "PROD" ) ) return

      !bug: timefrac = 3600.0/dt
      timefrac = dt/3600.0
     !/***** 3-D fields **************************

     do n = 1, NDERIV_3D

        if ( f_3d(n)%class  == "PROD " ) then
           index = f_3d(n)%index

           select case ( text )

               case ( "Before" )   !! Initialise to xn_adv

                 forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                   d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
                 end forall

               case ( "After" )    !! Calculate change

                 forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                   d_3d( n, i,j,k,IOU_INST) = &
                      d_3d( n, i,j,k,IOU_INST) - xn_adv(index,i,j,k)
                 end forall

           end select
        end if
      end do
     
    end subroutine DerivedProds
    !=========================================================================

    subroutine ResetDerived(period)
      integer, intent(in) :: period   ! Either IOU_DAY or IOU_MON

       if ( period <= LENOUT2D ) then
           nav_wdep(:,period) = 0.0
           nav_ddep(:,period) = 0.0
           nav_2d  (:,period) = 0.0

           wdep(:,:,:,period) = 0.0
           ddep(:,:,:,period) = 0.0
           d_2d(:,:,:,period) = 0.0
       end if 


       if ( period <= LENOUT3D ) then
           nav_3d    (:,period) = 0.0
           d_3d(:,:,:,:,period) = 0.0
       end if

    end subroutine ResetDerived
 !=========================================================================

  subroutine Setup_VOC()
      !--------------------------------------------------------
      ! Searches through the advected species and colects the
      ! index and carbon content of nmhc species, as they were
      ! defined in GenOut_ml
      !
      ! Works for jej and ds chem
      !--------------------------------------------------------
       integer :: n
   

     !ds rv1_9_10: VOC now properly defined. Previous definition
     !             was for NMHC, not VOC.

      do n = 1, NSPEC_ADV
        !ds rv1_9_10 if ( species( NSPEC_SHL+n )%nmhc == 1 ) then

        if ( species( NSPEC_SHL+n )%carbons > 0 .and. &
             species( NSPEC_SHL+n )%name   /= "CO"  .and. &
             species( NSPEC_SHL+n )%name   /= "CH4" ) then

             nvoc = nvoc + 1
             voc_index(nvoc) = n
             voc_carbon(nvoc) = species( NSPEC_SHL+n )%carbons
        end if
      end do
  end subroutine Setup_VOC
 !=========================================================================

   subroutine voc_2dcalc()

    !/-- Sums up voc species using the indices defined earlier in Setup_VOCs

     ! We initialise d_2d first, the use a simple loop
     ! over voc. Some CPU could be saved by initialising
     ! with the 1st voc, then looping over 2, nvoc, but who cares...

      
      d_2d( n, 1:limax,1:ljmax,IOU_INST) =  0.0

      do ivoc = 1, nvoc

         index = voc_index(ivoc)           ! Gives which IXADV_ to use.
         forall ( i=1:limax, j=1:ljmax )
             d_2d( n, i,j,IOU_INST) = d_2d( n, i,j,IOU_INST)      &
                                    + xn_adv(index,i,j,KMAX_MID)  &
                                    * voc_carbon(ivoc) * cfac(index,i,j)
                               ! multiplied by nr. of C and "reduced to surface"
         end forall
      end do ! ivoc
   end subroutine voc_2dcalc

 !=========================================================================
   subroutine voc_3dcalc()

    !/-- as for voc_2dcalc

      d_3d( n, 1:limax,1:ljmax,1:KMAX_MID,IOU_INST) =  0.0

      do ivoc = 1, nvoc

         index = voc_index(ivoc)
         forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
             d_3d( n, i,j,k,IOU_INST) = d_3d( n, i,j,k,IOU_INST) + &
                     xn_adv(index,i,j,k)*voc_carbon(ivoc)
         end forall
      end do ! ivoc

   end subroutine voc_3dcalc
 !=========================================================================

  subroutine aot_calc( n, timefrac )

    !/-- Calcuates AOT values for input threshold. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )
    !    Only relevant in ozone models, so far....

    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    real,    intent(in) :: timefrac    ! Timestep as fraction of hour

    real    :: threshold               ! Threshold, e.g. 40 or 60 (ppb)
    integer :: izen                    ! integer of zenith angle
    real :: o3                         ! Ozone (ppb) - needed if AOTs

     threshold = f_2d(n)%index

      do i=1,limax
        do j=1,ljmax

           izen = max(1,int( zen(i,j) + 0.5))

           if ( izen < AOT_HORIZON ) then
                o3 = xn_adv(IXADV_O3,i,j,KMAX_MID) &
                     * cfac(IXADV_O3,i,j) * PPBINV 

               !! jej 26/8/2002 already in mixing ratio

                o3 = max( o3 - threshold , 0.0 )   ! Definition of AOTs

             ! d_2d values will be accumulated in Derived_ml

              d_2d(n, i,j,IOU_INST ) = o3 * timefrac  

           else ! rv1_9_19 bug-fix!!
               d_2d(n, i,j,IOU_INST ) = 0.0   
           end if
        end do
      end do
   end subroutine aot_calc

!=========================================================================

  subroutine som_calc( n )

!pw rv2_1

    !/-- Calculates SOM (8hours) values for input threshold.

    implicit none
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays

    real    :: threshold               ! Threshold, e.g. 35 (ppb)
    real :: o3                         ! Ozone (ppb) - needed if SOMs
    real :: sum8h
    integer, parameter :: N8h = (NTDAY*8)/24 !number of periods in 8 hours
    real, parameter :: N8h_inv=1./N8h
    integer :: nh


    threshold = f_2d(n)%index

      do i=1,limax
        do j=1,ljmax

           !find max running 8h sum O3
           sum8h=0.
           do nh=1,N8h
              sum8h = sum8h + D2_O3_DAY( i , j , nh)
           enddo
           o3=sum8h
           do nh=N8h+1,NTDAY
              sum8h =sum8h-D2_O3_DAY( i , j , nh-N8h)+D2_O3_DAY( i , j , nh)
              o3=max(o3,sum8h)
              if(n<0)write(*,*)o3 !pw fake for compiler!!
           enddo

           !divide by N8h to find 8h mean 
           o3=o3*N8h_inv

           o3 = max( o3 - threshold , 0.0 )   ! Definition of SOMs

             ! d_2d values will be accumulated in Derived_ml

           d_2d(n, i,j,IOU_DAY ) = o3  

        end do
      end do
   end subroutine som_calc

 !=========================================================================

   subroutine setaccumulate_2dyear(n,accumulate_2dyear)

! We don't want the yearly output to accumulate over the whole year
     integer, intent(in) :: n
      logical, intent(inout) :: accumulate_2dyear !flag to know when to accumulate d_2d (case "EXT")

      if( f_2d(n)%name=="D2_EUAOT30DF".or.&
          f_2d(n)%name=="D2_EUAOT40DF".or.&
          f_2d(n)%name=="D2_UNAOT30DF".or.&
          f_2d(n)%name=="D2_UNAOT40DF"    &
          )then
         if(   current_date%month<startmonth_forest&
              .or.current_date%month>endmonth_forest)then
            accumulate_2dyear=.false.
         endif   
      endif

       if(f_2d(n)%name=="D2_EUAOT30WH".or.&
          f_2d(n)%name=="D2_EUAOT40WH".or.&
          f_2d(n)%name=="D2_UNAOT30WH".or.&
          f_2d(n)%name=="D2_UNAOT40WH"    &
          )then
         if(   current_date%month<startmonth_crops&
              .or.current_date%month>endmonth_crops)then
            accumulate_2dyear=.false.

         endif   
      endif

    end subroutine setaccumulate_2dyear

end module Derived_ml


