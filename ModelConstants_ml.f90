module ModelConstants_ml
 !+
 ! Specifies a number of constants used in the model. Note that
 ! physical constants (e.g. gravity, Cp, etc ( are specified in 
 ! the module PhysicalConstants_ml.f90)
 !
 ! Dependancies - none.
 ! Subroutines - assign_nmax.
 !----------------------------------------------------------------------------
  use Dates_ml, only : date   ! type giving yy, mm, dd, s
  use PhysicalConstants_ml, only : AVOG   
  implicit none
  private


 !/-- subroutine:

  public :: assign_nmax

 !/-- constants

  integer, public, parameter :: METSTEP   = 3  !u2 time-step of met. (h)
  real, public, save :: dt_advec  = 1200.0   ! time-step (s)

 !/-- choose temperature range: from 148 K (-125C) ro 333K (+60C).

 integer, parameter, public :: &
                 CHEMTMIN=148,CHEMTMAX=333    ! Min and max temp for rates, 

!
!  additional parameters, formerly set in par_ml, eulcon.inc
!
  integer, public, parameter :: &
    KMAX_MID   = 20           & ! Number of points (levels) in vertical
  , KMAX_BND   = KMAX_MID+1   & ! Number of points (levels) in vertical + 1
  , KTOP    = 1            & ! K-value at top of domain
  , NMET    = 2            & ! No. met fields in memory
  , KEMISTOP = KMAX_MID-3    ! k-limit for normal emissions input   !u4
!
! some values derived from MACHO

  integer, public, parameter :: &
    KCHEMTOP = 2           &  ! chemistry not done for k=1
  , KCLOUDTOP= 8           &  ! limit of clouds (for MADE dj ??) 
  , KUPPER   = 6           &  ! limit of clouds (for wet dep.)
  , AOT_HORIZON = 89          ! Limit of daylight zenith angle for AOTs

  real, public, parameter  ::    &
       CLOUDTHRES =  1.0e-5         !pw when cloudwater is larger than 
                                    !CLOUDTHRES, there are clouds. 
                                    !THIS VALUE MUST BE CHECKED BEFORE USE! 
 
!
!  additional parameters, formerly set in defcon, eulcon.inc
!

  integer, public, save   :: nterm, nmax, nstep, nprint,  nass, nbound

  ! was set in readpar_mach_ml, but not used!
  !real,    public, parameter   :: 
    ! alfa  = 1.5e-01   ! -> DRY DEP. FACTOR FOR SO2,SO4
    ! betaS = 5.0e-02   !  -> SO4 FACTOR OF SO2 EMISSIONS
    ! betaN = 1.0e-01   ! -> NO2 FACTOR OF NO  EMISSIONS

  integer, public, save , dimension(20)   :: identi   !! ????

  type(date), public, save :: current_date
  !u7.lu integer,    public, save :: daynumber    ! Days from 1st Jan   

  integer, public, parameter :: NNLANDUSE  = 10 ! Number of land use types 
                                                ! for rivm  (tmp)
      
  real, public, parameter  ::    &
       EPSIL=1.0e-30             &  ! small number
    ,  PASCAL=100.0              &  ! Conv. from hPa to Pa
    ,  PPB = 1.0e-9              &  ! parts per billion (mixing ratio)
    ,  PPBINV = 1.0e+9           &
    ,  PPT = 1.0e-12		 &  ! parts per trillion (mixing ratio)
    ,  PPTINV = 1.0e+12		 &
    ,  PT = 1.0e+4                  ! Top of model region = 100 hPa

   real, public, parameter ::  &
    ATWAIR = 28.964                   & ! Mol. weight of air (Jones, 1992)
  , atwS   = 32.                      & ! Atomic weight of Sulphur
  , atwN   = 14.                        ! Atomic weight of Nitrogen

  ! MFAC replaces earlier use of CHEFAC and ATWAIR - to scale from
  ! density (roa, kg/m3) to  molecules/cm3
  ! (kg/m3 = 1000 g/m3 = 0.001 * Avog/Atw molecules/cm3)

    real, public, parameter   :: MFAC = 0.001*AVOG/ATWAIR


!===========================================================================
! N2O5 -> nitrate calculation. Some constants for
! calculation of volume fraction of sulphate aerosol, and rate of uptake
! From EMEP Status Report 2/98
!
!     volume fraction of sulphate:
!     V = (so4 + ammonium sulphate - in moleculescm-3) x atw sulphate
!         ---------------------------------------------------------
!                 A0 X specific density of aerosols (assumed 4 OR 2 ??? )
!
!    Or, shorter, V = S x M0/(AVOG*rho)
!
!    where S is conc. sulphate (molecule/cm3), M0 is molwt. (96)
!
!    We do not want to include  concentrations  yet, so:
!
!     Let V0 =  M0/(AVOG*rho) = 96.0/(AVOG*2.0)
!   
!        
!  Rate coefficient, simplified form of Dentener and Crutzen
!    k =  V * 3/4*alpha * vav /raero
!        alpha is sticking coeff. for N2O5   (=0.1)
!        raero is mean aerosol radius in accumulation mode (0.034 x 10^-6 m)
!
!    Collect constants in VOLFAC:
!    VOLFAC =  V0* 3/4*alpha/raero
!===========================================================================
  real, parameter, public  :: VOLFAC = 96.0/(AVOG*2.0) * 0.75 *0.1/0.034e-6

 contains

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	subroutine assign_nmax(me,metstep)

	implicit none

!	input
	integer me,metstep

!	local
	integer nhelp
     
!     Assigne number of time-steps for the inner time-loop (over 6 hours)
!     from dt_advec

	nhelp = nint(dt_advec)
	if(mod(nhelp,60).ne.0) then
	  if (me .eq. 0) then
	    write(6,*)
	    write(6,*)'**********************************************'
	    write(6,*)'Impossible dt_advec,  dt_advec = (dt_advec/60) must be an integer'
	    write(6,*)
	  endif
	endif

	nhelp = nhelp/60

	if(mod(60,nhelp).ne.0) then
	  if (me .eq. 0) then
	    write(6,*)
	    write(6,*)'**********************************************'
	    write(6,*)'Impossible dt_advec,60/(dt_advec/60) must be an integer'
	    write(6,*)
	  endif
	endif

	nmax = 60/(nhelp)*metstep

	if (me .eq. 0) then
	  write(6,*)
	  write(6,*)'**********************************************'
	  write(6,*)'nmax and dt_advec : ',nmax,dt_advec
	  write(6,*)'**********************************************'
	  write(6,*)
	endif

	end subroutine assign_nmax
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


end module ModelConstants_ml
!_____________________________________________________________________________
!REMOVED:
!g11   ,  CHEFAC=6.023e20           &  ! Avog x 10-3
!gv    ,  PINC=1000.0               &
!gv    ,  PBAS=-PINC                &
!g11   ,  CHEFAC=6.023e20           &  ! Avog x 10-3
!gv    ,  PINC=1000.0               &
!gv    ,  PBAS=-PINC                &
