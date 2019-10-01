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

  ! version string:
  character(len=*), parameter   ::  DA_VERSION = '3DVar17 (patch 2019-09-30)'

  ! print debug messages ?
  logical, parameter            ::  DEBUG_DA     = .false.  ! general purpose debug messages
  logical, parameter            ::  DEBUG_DA_OBS = .false.  ! observation info
  logical, parameter            ::  DEBUG_DA_3DV = .false.  ! 3DVar module

#ifdef with_ajs
  ! run only 1 DA step (no adv/chem) ;
  ! here define as variable, since reset using rcfile settings:
  logical                       ::  DEBUG_DA_1STEP  = .false.
#else
  ! run only 1 DA step (no adv/chem)
  logical, parameter            ::  DEBUG_DA_1STEP  = .false.
#endif

  ! write formats:
  character(len=*), parameter   ::  DA_FMT_DEF  = "('3DVar@PPP YYYY-MM-DD hh: ',A,'.')"
  character(len=*), parameter   ::  NMC_FMT_DEF = "('B-NMC@PPP YYYY-MM-DD hh: ',A,'.')"

  
  ! --- var -----------------------------------------

  ! short variables for model domain shape:
  integer                           ::  nlev

  ! storage for time stamps:
  character(len=len(DA_FMT_DEF))    ::  da_fmt_msg=''
  ! message line:
  character(len=128)                ::  da_msg=''

  ! timers:
  real, save                        :: datim_before, datim_after

end module DA_ml
