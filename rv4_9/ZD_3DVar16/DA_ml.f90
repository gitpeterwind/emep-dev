!*****************************************************************************!
!
! Shared data for DA modules.
!
!*****************************************************************************!

module DA_ml
  
  !use ModelConstants_ml, only : PPT, PPTINV

  implicit none
  
  
  ! --- in/out ----------------------------------------
  
  public
  
  
  ! --- const -----------------------------------------

  ! print debug messages ?
  logical, parameter            ::  DEBUG_DA     = .false.  ! general purpose debug messages
  logical, parameter            ::  DEBUG_DA_OBS = .false.  ! observation info
  logical, parameter            ::  DEBUG_DA_3DV = .false.  ! 3DVar module

#ifdef with_ajs
  ! run only 1 DA step (no adv/chem)
  logical                       ::  DEBUG_DA_1STEP  = .false.
  ! hourly output before/after DA step
  logical                       ::  DEBUG_DA_OUTPUT = .false.
#else
  ! run only 1 DA step (no adv/chem)
  logical, parameter            ::  DEBUG_DA_1STEP  = .false.
  ! hourly output before/after DA step
  logical, parameter            ::  DEBUG_DA_OUTPUT = .false.
#endif

  ! settings:
  character(len=*), parameter   ::  DA_NAMELIST = "namelist.nml"
  
  ! write formats:
  character(len=*), parameter   ::  DA_FMT_DEF  = "('3DVar@PPP YYYY-MM-DD hh: ',A,'.')"
  character(len=*), parameter   ::  NMC_FMT_DEF = "('B-NMC@PPP YYYY-MM-DD hh: ',A,'.')"

  
  ! --- var -----------------------------------------

  ! to be filled in DA_3DVar_ml based on covariance units:
  real                              :: FGSCALE
  real                              :: FGSCALE_INV

  ! short variables for model domain shape:
  !integer                           ::  nx, ny
  integer                           ::  nlev

  integer                           ::  nchem
  integer                           ::  nchemobs
  integer, allocatable              ::  ichemObs(:)
  integer, allocatable              ::  ichemInv(:)

  ! storage for time stamps:
  character(len=len(DA_FMT_DEF))    ::  da_fmt_msg=''
  ! message line:
  character(len=128)                ::  da_msg=''

  ! timers:
  real, save                        :: datim_before, datim_after

end module DA_ml
