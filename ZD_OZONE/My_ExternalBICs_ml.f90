module My_ExternalBICs_ml
!External Boundary and Initial Conditions, if needed:
! Here = Dummy values - when BICs set using EMEP defaults
! Replace this module by eg MOZART_ExternalBICs_ml if using
! input from other (eg global) models.
 use Io_ml, only : PrintLog
 implicit none
 private

  public :: set_extbic

  logical, public, parameter :: EXTERNAL_BIC_SET  = .false.
  logical, public, parameter :: TOP_BC = .false.
  character(len=30),public, parameter :: EXTERNAL_BIC_NAME = "DUMMY"

  type, public :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
    character(len=24) :: varname=""
    logical           :: wanted=.false.,found=.false.
  end type icbc

  type, public :: adv_icbc             ! IC/BC Set, included intended ixadv
    integer           :: ixadv=-1
    type(icbc)        :: icbc
  end type adv_icbc


  type(adv_icbc), dimension(9:0), public, save :: EXTERNAL_BC

contains
  subroutine set_extbic()
     call PrintLog("No external BICs set")
  end subroutine set_extbic

end module My_ExternalBICs_ml
