!_____________________________________________________________________________
! Simplified version of EMEP module
! Some basic model variables
!
module ModelConstants
  implicit none
  private

  !integer, public, parameter :: dp = kind(0.0d0)  ! Double precision real kind

  ! Sentinel values
  !real, public, parameter :: UNDEF_D = -huge(0.0_dp)
  real,    public, parameter :: UNDEF_R = -huge(0.0)
  integer, public, parameter :: UNDEF_I = -huge(0)


!/-- choose temperature range: from 148 K (-125C) ro 333K (+60C).
 integer, parameter, public :: &
                 CHEMTMIN=148,CHEMTMAX=333    ! Min and max temp for rates, 

!  integer, public, parameter :: &
!    KMAX_MID   = 3           & ! Number of points (levels) in vertical
!  , KTOP       = 1           & ! K-value at top of domain
!  , KCHEMTOP   = 1           &  ! chemistry not done for k=1
!  , KCLOUDTOP  = 1           &  ! limit of clouds (for MADE dj ??) 
!  , KUPPER     = 2              ! limit of clouds (for wet dep.)

  real, public, save :: dt_advec = 1200.0  ! Time step for CTM advection.
                                           ! => Outer loop

  logical, public, parameter :: DebugCell=.true.
  logical, public, parameter :: MasterProc=.true.

! OLD STUFF FROM EMEP.  WILL REMOVE MOST
  logical, public, parameter :: &
     DEBUG_DRYRUN=.false.  &
    ,DEBUG_SOLVER=.true.   &
    ,DEBUG_RUNCHEM=.true.


end module ModelConstants
!_____________________________________________________________________________
