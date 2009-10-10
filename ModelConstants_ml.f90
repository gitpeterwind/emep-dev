! <ModelConstants_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
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
!+ 1) Define first dimensions that might change quite often -  for different
!     run domains or debug points:

  integer, public, parameter, dimension(4) ::  &
 !                    x0   x1  y0   y1
    RUNDOMAIN = (/  36, 167, 12, 122 /)     ! EMEP domain
  !    RUNDOMAIN = (/  56, 147, 12, 102 /)     ! EGU
  !    RUNDOMAIN = (/  1, 360, 1, 180 /)     ! FULL GLOBAL
  !  RUNDOMAIN = (/  1, 132, 1, 111 /)       ! EECCA, rep09

   !  RUNDOMAIN = (/  85, 120,  55,  70 /)     ! (changeable)
    !  RUNDOMAIN = (/  20, 167,  1, 122 /)     !  OSPAR/HELCOM domain
    !  RUNDOMAIN = (/  18, 169,  1, 124 /)     !  OSPAR/HELCOM domain+borders

  integer, public, parameter ::  &
    NPROCX      =   8       & ! Actual number of processors in longitude
  , NPROCY      =   4       & ! .. in latitude. Must be 2 for GLOBAL!!!
  , NPROC       = NPROCX * NPROCY

  ! ds Jan2009
  ! We have one variable, to say if we are on master-processor
  ! or not: (kept here to avoid too many dependencies for box-model
  ! codes which don't need Par_ml

  logical, public, save ::  MasterProc = .true.
  logical, public, save ::  DebugCell  = .false.

  ! For debugging, we often want to print out for  a specific location
  ! Set here:

 ! The coordinates given here only apply for the standard EMEP domain
 !integer, public, parameter :: DEBUG_i=79, DEBUG_j=56 ! Eskdalemuir
 !integer, public, parameter :: DEBUG_i=73, DEBUG_j=48 ! Mace Head
 !integer, public, parameter :: DEBUG_i=91, DEBUG_j=71 ! Rorvik
 !!integer, public, parameter :: DEBUG_i=82, DEBUG_j=72 !  Voss, has some snow
 !!!!integer, public, parameter :: DEBUG_i=110, DEBUG_j=48 !   High Vg!
 !integer, public, parameter :: DEBUG_i=96, DEBUG_j=40 !   High VG_SO2_CF!
 !integer, public, parameter :: DEBUG_i=111, DEBUG_j=54 !  High VG_PMCO_CF!
 !integer, public, parameter :: DEBUG_i=101, DEBUG_j=51 !  Schauinsland
 ! integer, public, parameter :: DEBUG_i=87, DEBUG_j=20 !  Aveiro
 !integer, public, parameter :: DEBUG_i=103, DEBUG_j=50 !  Mid-Europe
 !integer, public, parameter :: DEBUG_i=93, DEBUG_j=57 !  Elspeetsche (52d12',5d45') 92.83, 56.64
 !integer, public, parameter :: DEBUG_i=97, DEBUG_j=62 !  Waldhof
 integer, public, parameter :: DEBUG_i=116, DEBUG_j=63 !  K-Puszta

!=============================================================================
! Debug flag DEBUG_XXX  applied in subroutine XXX
 logical, public, parameter ::      &
     DEBUG_AQUEOUS        = .false.  & ! 
    ,DEBUG_BCS            = .false.  & ! 
    ,DEBUG_BIO            = .false.  & !
    ,DEBUG_DERIVED        = .false. &  !
    ,DEBUG_MY_DERIVED     = .false. &  !
    ,DEBUG_DRYDEP         = .false. &  !
      ,DEBUG_VDS          = .false. &  !
    ,DEBUG_MY_DRYDEP      = .false. &  !
    ,DEBUG_EMISSIONS      = .false.   &
    ,DEBUG_GETEMIS        = .false.   &
    ,DEBUG_IOPROG         = .false.   &
 !!! DEBUG_RUNCHEM is SPECIAL.. needed for indented debugs are to work
    ,DEBUG_RUNCHEM        = .false.   &
        ,DEBUG_SOA        = .false.   & !
        ,DEBUG_SUBMET     = .false.  &
    ,DEBUG_LANDDEFS       = .false. &  !
    ,DEBUG_RSUR           = .false. &  !
    ,DEBUG_SETUP_1DCHEM   = .false. &  !
    ,DEBUG_SETUP_1DBIO    = .false.    !

!=============================================================================
  ! Source-receptor runs?
  ! We don't (generally) want daily outputs for SR runs, so in
  ! Derived_ml, we set all IOU_DAY false if SOURCE_RECPTOR = .true..

    logical, public, parameter :: SOURCE_RECEPTOR = .false.


!=============================================================================
!+ 2)  Define domain-name,  something that will
!       generally only change when switching Met-driver or large domain

  !character(len=20), parameter, public :: DomainName = "EMEP-50kmEurope"
  character(len=20), parameter, public :: DomainName = "EMEP-50kmEECCA"

  logical, parameter, public :: IS_GLOBAL = .false.

!=============================================================================
!+ 3)  Define main model dimensions,  things that will
!       generally only change when switching Met-driver or large domain
  integer, public, parameter ::  &
!EMEP:
   IIFULLDOM   = 170   &    ! x-Dimensions of full domain
!EMEP:
   , JJFULLDOM   = 133   &    ! y-Dimensions of full domain
!EECA:   IIFULLDOM   = 132   &    ! x-Dimensions of full domain
!EECA: ,JJFULLDOM   = 159   &    ! y-Dimensions of full domain
!GLOBAL:   IIFULLDOM   = 360   &    ! x-Dimensions of full domain
!GLOBAL: , JJFULLDOM   = 180   &    ! y-Dimensions of full domain
  , NLANDUSEMAX = 23    &    ! Number of land use types in Inputs.Landuse file
!
  , METSTEP      = 3     &    ! time-step of met. (h)
  , KMAX_MID     = 20    &    ! Number of points (levels) in vertical
  , KMAX_BND     = KMAX_MID+1 & ! Number of points (levels) in vertical + 1
  , KTOP         = 1     &    ! K-value at top of domain
  , NMET         = 2     &    ! No. met fields in memory
  , KCHEMTOP     = 2     &    ! chemistry not done for k=1
  , KCLOUDTOP    = 8     &    ! limit of clouds (for MADE dj ??) 
  , KUPPER       = 6     &    ! limit of clouds (for wet dep.)
  , AOT_HORIZON  = 89         ! Limit of daylight zenith angle for AOTs
      

! EMEP measurements end at 6am, used in  daily averages
  integer, public, parameter :: END_OF_EMEPDAY  = 6  

  real, public :: dt_advec   ! time-step for advection (s)
  real, public :: dt_advec_inv  ! =1/dt_advec

  !NTDAY:  Number of 2D O3 to be saved each day (for SOMO)       
  ! 24/NTDAY is the time integration step for SOMO
  !large value-> large memory use; too small value ->bad approximation for SOMO
  !NB must be choosen:  24*3600/dt_advec <= NTDAY >=3 and 
  !preferably an integer fraction of 24*3600/dt_advec
  integer, public, parameter ::      NTDAY = 72  

 !/-- choose temperature range: from 148 K (-125C) ro 333K (+60C).

 integer, parameter, public :: CHEMTMIN=148,CHEMTMAX=333

  real, public, parameter  ::    &
      V_RAIN   = 5.              &  ! pw approximate vertical speed of rain m/
     ,CLOUDTHRES =  1.0e-5         !pw when cloudwater is larger than 
                                   !CLOUDTHRES, there are clouds. 
                                   !THIS VALUE MUST BE CHECKED BEFORE USE! 
!
!  additional parameters
!
  integer, public, save   :: nterm, nmax, nstep, nprint,  nass, nbound

  integer, public, save   :: iyr_trend ! Year specified for say BC changes

  integer, public, save , dimension(20)   :: identi   !! ????

  character(len=120), public, save :: runlabel1& !SHORT Allows explanatory text
                                     ,runlabel2 !LONG  Read in from grun.pl

  real, public, parameter  ::    &
       EPSIL=1.0e-30             &  ! small number
    ,  PASCAL=100.0              &  ! Conv. from hPa to Pa
    ,  PPB = 1.0e-9              &  ! parts per billion (mixing ratio)
    ,  PPBINV = 1.0e+9           &
    ,  PPT = 1.0e-12		 &  ! parts per trillion (mixing ratio)
    ,  PPTINV = 1.0e+12		 &
    ,  PT = 1.0e+4                  ! Top of model region = 10000 Pa = 100 hPa

   real, public, parameter ::  &
    ATWAIR = 28.964                   & ! Mol. weight of air (Jones, 1992)
  , atwS   = 32.                      & ! Atomic weight of Sulphur
  , atwN   = 14.                      & ! Atomic weight of Nitrogen
  , atwPM  = 100.                    

  ! MFAC replaces earlier use of CHEFAC and ATWAIR - to scale from
  ! density (roa, kg/m3) to  molecules/cm3
  ! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)

    real, public, parameter   :: MFAC = 0.001*AVOG/ATWAIR


end module ModelConstants_ml
!_____________________________________________________________________________
