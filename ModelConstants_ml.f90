! <ModelConstants_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2011 met.no
!*
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!
module ModelConstants_ml
 !+
 ! Specifies a number of constants used in the model. Note that
 ! physical constants (e.g. gravity, Cp, etc ( are specified in
 ! the module PhysicalConstants_ml.f90)
 !
 !----------------------------------------------------------------------------
  use PhysicalConstants_ml, only : AVOG

  implicit none
  private

!=============================================================================
! Some flags for model setup
! will be removed when code is sufficiently tested 
! (for convection use foundconv in permanent code)
 logical, public, parameter :: USE_CONVECTION   = .false.
 logical, public, parameter :: USE_SOILWATER    = .false.
 logical, public, parameter :: USE_FOREST_FIRES = .false.
 logical, public, parameter :: USE_SOIL_NOX     = .false.
 logical, public, parameter :: USE_BVOC_2010    = .false.
 logical, public, parameter :: USE_DUST         = .false.
 logical, public, parameter :: USE_PFT_MAPS     = .false.
 logical, public, parameter :: DO_SAHARA = .false.     ! Turn on/off BG Saharan Dust
! Biogenics. Use 3 even if no terpene chemistry - simplifies
! rest of code.
! iso = isoprene, mtp = monoterpenes from pools, mtl = monoterpenes
! with light dependence
  integer, public, parameter ::   NBVOC = 3
  character(len=4),public, save, dimension(NBVOC) :: &
       BVOC_USED = (/ "Eiso","Emt ","Emtl"/)
       !!BVOC_USED = (/ "isoprene","terpene "/)


  ! Next      ................................................................
  ! Next      ................................................................

!=============================================================================
!+ 1) Define first dimensions that might change quite often -  for different
!     run domains

 character(len=20), parameter, public :: DomainName = "EMEP-50kmEurope"
 !character(len=20), parameter, public :: DomainName = "EMEP-50kmEECCA"
 !character(len=22), parameter, public :: DomainName = "EMEPCWF-0.25degEurope"
 !character(len=22), parameter, public :: DomainName = "EMEPCWF-0.20degEurope"
 !character(len=20), parameter, public :: DomainName = "HIRHAM"

  logical, parameter, public :: IS_GLOBAL = .false.

  integer, public, parameter ::  &
  !  IIFULLDOM = 182, JJFULLDOM = 197  ! x,y-Dimensions of full HIRHAM domain 
   IIFULLDOM = 170, JJFULLDOM = 133  ! x,y-Dimensions of full EMEP domain
  ! IIFULLDOM = 132, JJFULLDOM = 159  ! x,y-Dimensions of full EECA domain
  ! IIFULLDOM = 360, JJFULLDOM = 180 ! .... full GLOBAL domain
  ! IIFULLDOM = 201, JJFULLDOM = 161 ! .... full GEMS 0.25 domain
  ! IIFULLDOM = 301, JJFULLDOM = 221 ! .... full GEMS 0.25 extended domain
  ! IIFULLDOM = 321, JJFULLDOM = 221 ! .... full MACC 0.20 domain

 !ds - I added these offsets, but now suspect I was thinking wrong.
 ! The difference between EMEP and EECCA is confusing...
 ! integer, public, parameter :: OFFSET_i= -35, OFFSET_j= -11 ! EECCA
  integer, public, parameter :: OFFSET_i= 0, OFFSET_j= 0 ! EMEP
  integer, public, parameter, dimension(4) ::  &
  !                x0   x1  y0   y1
  ! RUNDOMAIN = (/ 1, 182, 1, 197 /)      ! HIRHAM
  !RUNDOMAIN = (/  1, 132,  1, 159 /)     ! EECCA
  !RUNDOMAIN = (/ 1, 100, 1, 100 /)     ! EMEP domain in EECCA
  RUNDOMAIN = (/ 36, 167, 12, 122 /)     ! EMEP domain
  !RUNDOMAIN = (/ 56, 147, 12, 102 /)     ! EGU
  !RUNDOMAIN = (/ 75, 137, 32,  82 /)     ! EGU
  !RUNDOMAIN = (/  1, 360,  1, 180 /)     ! FULL GLOBAL
  !RUNDOMAIN = (/  1, 132,  1, 111 /)     ! EECCA, rep09
  !RUNDOMAIN = (/  1, 132,  1, 159 /)     ! EECCA, rep10
  !RUNDOMAIN = (/ 20, 167,  1, 122 /)     ! OSPAR/HELCOM domain
  !RUNDOMAIN = (/ 18, 169,  1, 124 /)     ! OSPAR/HELCOM domain+borders
  !RUNDOMAIN = (/  1, 201,  1, 161 /)     ! EMEP-CWF, GEMS 0.25 domain
  !RUNDOMAIN = (/  1, 301, 26, 221 /)     ! EMEP-CWF, GEMS 0.25 extended domain
  !RUNDOMAIN = (/  1, 321,  1, 221 /)     ! EMEP-CWF, MACC 0.20 domain
  !RUNDOMAIN = (/ 70+OFFSET_i, 90+OFFSET_i, 43+OFFSET_j,  63+OFFSET_j /)     ! (UK)
  !RUNDOMAIN = (/ 85+OFFSET_i, 120+OFFSET_i, 55+OFFSET_j,  70+OFFSET_j /)     ! (changeable)
  !RUNDOMAIN = (/ 75+OFFSET_i, 110+OFFSET_i, 45+OFFSET_j,  60+OFFSET_j /)     ! (gets Esk)
  !RUNDOMAIN = (/ 85+OFFSET_i, 120+OFFSET_i, 15+OFFSET_j,  40+OFFSET_j /)     ! (changeable)

  integer, public, parameter ::  &
    NPROCX      =   8        & ! Actual number of processors in longitude
  , NPROCY      =   8        & ! .. in latitude. NPROCY must be 2 for GLOBAL,
  , NPROC       = NPROCX * NPROCY

!=============================================================================
!+ 2) Define  debug flags.

  ! We have one variable, to say if we are on master-processor
  ! or not: (kept here to avoid too many dependencies for box-model
  ! codes which don't need Par_ml

  logical, public, save ::  MasterProc = .true.
  logical, public, save ::  DebugCell  = .false.

  ! For debugging, we often want to print out for  a specific location
  ! Set here:

 ! The coordinates given here only apply for the standard EMEP domain
 !integer, private, parameter :: DEBUG_ii= 79, DEBUG_jj= 56 ! Eskdalemuir
 !integer, private, parameter :: DEBUG_ii= 73, DEBUG_jj= 48 ! Mace Head
 !integer, private, parameter :: DEBUG_ii= 91, DEBUG_jj= 71 ! Rorvik
 !integer, private, parameter :: DEBUG_ii= 82, DEBUG_jj= 72 ! Voss, has some snow
 !integer, private, parameter :: DEBUG_ii=110, DEBUG_jj= 48 ! High Vg!
 !integer, private, parameter :: DEBUG_ii= 96, DEBUG_jj= 40 ! High VG_SO2_CF!
 !integer, private, parameter :: DEBUG_ii=111, DEBUG_jj= 54 ! High VG_PMCO_CF!
 !integer, private, parameter :: DEBUG_ii=101, DEBUG_jj= 51 ! Schauinsland
 !QUERY? integer, private, parameter :: DEBUG_ii= 87, DEBUG_jj= 20 ! Aveiro
 !QUERY? 
!integer, private, parameter :: DEBUG_ii= 88, DEBUG_jj= 21 ! Aveiro+i2
 !integer, private, parameter :: DEBUG_ii=103, DEBUG_jj= 50 ! Mid-Europe
 ! integer, private, parameter :: DEBUG_ii= 93, DEBUG_jj= 57 ! Elspeetsche (52d12',5d45') 92.83, 56.64
  integer, private, parameter :: DEBUG_ii= 92, DEBUG_jj= 56 ! Cabauw
 !integer, private, parameter :: DEBUG_ii= 97, DEBUG_jj= 62  ! Waldhof
 !integer, private, parameter :: DEBUG_ii=116, DEBUG_jj= 63 ! K-Puszta
 !integer, private, parameter :: DEBUG_ii=102, DEBUG_jj= 48 !  Payerne
 !integer, private, parameter :: DEBUG_ii=74, DEBUG_jj= 79 !  Payerne_HIRHAM
 !integer, private, parameter :: DEBUG_ii=85, DEBUG_jj= 50 !   Harwell
 !integer, private, parameter :: DEBUG_ii=85, DEBUG_jj= 15 !  biomass burnung, Aug 2003
 !integer, private, parameter :: DEBUG_ii=85, DEBUG_jj= 35 ! Sea, Bay of Biscay
 !integer, private, parameter :: DEBUG_ii=76, DEBUG_jj= 35 ! Sea,  North sea
 !integer, private, parameter :: DEBUG_ii=91, DEBUG_jj=67  ! Tange 
! integer, private, parameter :: DEBUG_ii=103, DEBUG_jj=32 ! Prades, SMDge 

 !integer, public, parameter :: DEBUG_i= 9, DEBUG_j= 201 ! MACC02

 integer, public, parameter :: &
    DEBUG_i= DEBUG_II+OFFSET_i, DEBUG_j= DEBUG_JJ+OFFSET_j

!=============================================================================
! Some flags for model setup

  ! Next      ................................................................
  ! Next      ................................................................
!=============================================================================
! Debug flag DEBUG_XXX  applied in subroutine XXX
 logical, public, parameter ::      &
     DEBUG_AQUEOUS        = .false. & !
    ,DEBUG_AOT            = .false. & !
    ,DEBUG_BCS            = .false. & !
    ,DEBUG_BIO            = .false. & !
    ,DEBUG_DERIVED        = .false. & !
       ,DEBUG_COLUMN      = .false. & ! extra option in Derived
    ,DEBUG_DO3SE          = .false. & !
    ,DEBUG_DRYRUN         = .false. & ! Skips chemistry and advection
    ,DEBUG_ECOSYSTEMS     = .false. & !
    ,DEBUG_FORESTFIRE     = .false. & !
    ,DEBUG_MET            = .false. & !
    ,DEBUG_BLM            = .false. & !    produces matrix of differnt Kz and Hmix
    ,DEBUG_Kz             = .false. & !
    ,DEBUG_MY_DERIVED     = .false. & !
    ,DEBUG_DRYDEP         = .false. & !
      ,DEBUG_VDS            = .false. & !
      ,DEBUG_MY_DRYDEP      = .false. & !
      ,DEBUG_CLOVER         = .false. & !
      ,DEBUG_STOFLUX        = .false. &
    ,DEBUG_EMISSIONS      = .false. &
    ,DEBUG_GETEMIS        = .false. &
    ,DEBUG_IOPROG         = .false. &
 !!! DEBUG_RUNCHEM is SPECIAL.. needed for indented debugs are to work
    ,DEBUG_MOSAICS        = .false.  & !
    ,DEBUG_RUNCHEM        = .false. &
      ,DEBUG_AEROSOL        = .false. & !
      ,DEBUG_MY_WETDEP      = .false. & !
      ,DEBUG_SOA            = .false. & !
      ,DEBUG_SOLVER         = .false. & !
      ,DEBUG_WETDEP         = .false. & !
    ,DEBUG_LANDDEFS       = .false. & !
    ,DEBUG_LANDUSE        = .false. & !
    ,DEBUG_LANDPFTS       = .false. &
    ,DEBUG_NETCDF         = .false. &
    ,DEBUG_NETCDF_RF      = .false. &  ! ReadField_CDF in NetCDF_ml
    ,DEBUG_PHYCHEM        = .false.  & !
    ,DEBUG_NH3            = .false. & ! hb NH3Emis
    ,DEBUG_RSUR           = .false. & !
    ,DEBUG_SEASALT        = .false. & !
    ,DEBUG_SOILNO         = .false. & !
    ,DEBUG_SUBMET         = .false. &
    ,DEBUG_SETUP_1DCHEM   = .false. & !
    ,DEBUG_SETUP_1DBIO    = .false. & !
    ,DEBUG_SITES          = .false. & !
    ,DEBUG_SOILWATER      = .false. & !
    ,DEBUG_VOLC           = .false. & ! volcanoes
    ,DEBUG_NEST_ICBC      = .false.   ! IFS-MOZART BC


!=============================================================================
  ! 3)  Source-receptor runs?
  ! We don't (generally) want daily outputs for SR runs, so in
  ! Derived_ml, we set all IOU_DAY false if SOURCE_RECPTOR = .true..

    logical, public, parameter :: SOURCE_RECEPTOR = .false.

  ! Forecast run?
  ! only dayly and hourly output is required on FORECAST mode, so in Derived_ml,
  ! we set all other output types to false if FORECAST=.true..
    logical, public, parameter :: FORECAST = .false.

  ! NH3 module as set up originally with U10 from met:
  !dshb - kept for safety only. Will be replaced by
  !sub.grid calculation of wind in future.

    logical, public, parameter :: NH3_U10 = .false.

  ! Nesting modes:
  ! produces netcdf dump of concentrations if wanted, or initialises mode runs from
  ! such a file. Used in Nest_ml

  integer,parameter ::NEST_MODE=0   !0=donothing , 1=write , 2=read , 3=read and write
  !10=write at end of run, 11=read at start , 12=read at start and write at end (BIC)


!=============================================================================
!+ 4)  Define main model dimensions,  things that will
!       generally only change when switching Met-driver
  integer, public, parameter ::  &
    NLANDUSEMAX  = 23    &    ! Number of land use types in Inputs.Landuse file
  , METSTEP      = 3     &    ! time-step of met. (h)
  , KMAX_MID     = 20    &    ! Number of points (levels) in vertical
  , KMAX_BND     = KMAX_MID+1 & ! Number of points (levels) in vertical + 1
  , KTOP         = 1     &    ! K-value at top of domain
  , KWINDTOP     = 5     &    ! Define extent needed for wind-speed array
  , NMET         = 2     &    ! No. met fields in memory
  , KCHEMTOP     = 2     &    ! chemistry not done for k=1
  , KCLOUDTOP    = 8     &    ! limit of clouds (for MADE dj ??)
  , KUPPER       = 6     &    ! limit of clouds (for wet dep.)
! And for My_Derived_ml VG_SPECS & LocalVariables_ml Vg_ref etc. we need a limit
  , NVGOUT_MAX   = 10         ! Max. no species for My_Derived VG outputs


  integer, public, parameter :: &   ! Groups for DDEP and WDEP
    SOX_INDEX = -1, OXN_INDEX = -2, RDN_INDEX = -3

! EMEP measurements end at 6am, used in  daily averages
  integer, public, parameter :: END_OF_EMEPDAY  = 6

  real, public, save :: dt_advec     ! time-step for advection (s)
  real, public, save :: dt_advec_inv ! =1/dt_advec

  ! NTDAY:  Number of 2D O3 to be saved each day (for SOMO)
  ! 24/NTDAY is the time integration step for SOMO
  ! large value -> large memory use; too small value -> bad approx. for SOMO
  ! NB must be choosen:  24*3600/dt_advec <= NTDAY >=3 and
  ! preferably an integer fraction of 24*3600/dt_advec
  integer, public, parameter :: NTDAY = 72

 !/-- choose temperature range: from 148 K (-125C) ro 333K (+60C).
  integer, parameter, public :: CHEMTMIN=148,CHEMTMAX=333

  real, public, parameter :: &
    V_RAIN     = 5.          & !pw approximate vertical speed of rain m/
  , CLOUDTHRES = 1.0e-5        !pw when cloudwater is larger than
                               !CLOUDTHRES, there are clouds.
                               !THIS VALUE MUST BE CHECKED BEFORE USE!
  real, public, parameter :: &
       CW_THRESHOLD = 1.0E-7&!Cloudwater (kg/kg); above threshold allow possibility 
                             ! for precipitations. Value could be adjusted.
     , RH_THRESHOLD = 0.85   !Relative humidity (fraction); above threshold allow 
                             ! possibility for precipitations.Value could be adjusted.
!
!  additional parameters
!
  integer, public, save   :: nterm, nmax, nstep, nprint,  nass, nbound

  integer, public, save   :: iyr_trend ! Year specified for say BC changes

  integer, public, save , dimension(20)   :: identi   !! ????
  integer, public, parameter :: TXTLEN_NAME = 20
  character(len=120), public, save :: runlabel1&!SHORT Allows explanatory text
                                    , runlabel2 !LONG  Read in from grun.pl

  real, public, parameter :: &
    EPSIL  = 1.0e-30         &  ! small number
  , PASCAL = 100.0           &  ! Conv. from hPa to Pa
  , PPB    = 1.0e-9          &  ! parts per billion (mixing ratio)
  , PPBINV = 1.0e+9          &
  , PPT    = 1.0e-12         &  ! parts per trillion (mixing ratio)
  , PPTINV = 1.0e+12         &
  , PT     = 1.0e+4          &  ! Top of model region = 10000 Pa = 100 hPa
  , Pref   = 101325.0        &  ! Reference pressure in Pa
  , TINY   = 1.0e-9             ! -1.E-9" is sometimes used  in order to avoid 
                                ! different roundings on different machines.

  real, public, parameter :: &
    ATWAIR = 28.964          & ! Mol. weight of air (Jones, 1992)
  , atwS   = 32.             & ! Atomic weight of Sulphur
  , atwN   = 14.             & ! Atomic weight of Nitrogen
  , atwPM  = 100.

  ! MFAC replaces earlier use of CHEFAC and ATWAIR - to scale from
  ! density (roa, kg/m3) to  molecules/cm3
  ! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)
    real, public, parameter   :: MFAC = 0.001*AVOG/ATWAIR


  ! Define 4 output types corresponding to instantaneous,year,month,day

   integer, public, parameter ::  &
        IOU_INST=1, IOU_YEAR=2, IOU_MON=3, IOU_DAY=4, &
        IOU_HOUR=5, IOU_HOUR_MEAN=6

   character(len=10),  public ,parameter :: model="EMEP_MSC-W"

end module ModelConstants_ml
!_____________________________________________________________________________
