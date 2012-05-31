module My_ExternalBICs_ml
!External Boundary and Initial Conditions, if needed:
! Here = Dummy values - when BICs set using EMEP defaults
! Replace this module by eg MOZART_ExternalBICs_ml if using
! input from other (eg global) models.
use ModelConstants_ml,     only: MasterProc, DEBUG=>DEBUG_NEST_ICBC
use Io_ml,                 only: PrintLog
use TimeDate_ExtraUtil_ml, only: date2string
implicit none

private
public :: set_extbic

logical, public, parameter :: &
  EXTERNAL_BIC_SET  = .false.,&
  TOP_BC = .false.

character(len=30),public, parameter :: &
  EXTERNAL_BIC_NAME = "DUMMY"

! YYYY, YY, MM, DD, hh will be replaced by numbers by the program.
! Search for date2string in set_extbic and uncomment lines, if necessary.
! For details, see detail2str in TimeDate_ExtraUtil_ml.f90
character(len=*),private, parameter :: &
  template_read_3D = 'EMEP_IN_IC.nc'         , &
  template_read_BC = 'EMEP_IN_BC_YYYYMMDD.nc', &
  template_write   = 'EMEP_OUT.nc'
character(len=len(template_read_3D)),public, save :: &
  filename_read_3D = template_read_3D
character(len=len(template_read_BC)),public, save :: &
  filename_read_BC = template_read_BC
character(len=len(template_write)),public, save :: &
  filename_write   = template_write

character(len=*),public, parameter :: &
  filename_eta     = 'EMEP_IN_BC_eta.zaxis'

type, public :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
  character(len=24) :: varname=""
  logical           :: wanted=.false.,found=.false.
end type icbc

type,public :: adv_icbc             ! IC/BC Set, included intended ixadv
  integer           :: ixadv=-1
  type(icbc)        :: icbc=icbc('',.false.,.false.)
end type adv_icbc

type(adv_icbc), dimension(9:0), public, save :: EXTERNAL_BC

contains
subroutine set_extbic(idate)
  implicit none
  integer,intent(in),optional :: idate(4)
  logical, save :: first_call = .false.

  ! Set filename from idate
  if(present(idate))then
    if(MasterProc.and.DEBUG)&
      write(*,"(A,':',1X,A,1X,'''',A,'''.')")"set_extbic DEBUG",&
        "External BICs filenames for",EXTERNAL_BIC_NAME
!   filename_read_3D=date2string(trim(template_read_3D),idate,&
!     debug=MasterProc.and.DEBUG)
    filename_read_BC=date2string(trim(template_read_BC),idate,&
      debug=MasterProc.and.DEBUG)
!   filename_write=date2string(trim(template_write),idate,&
!     debug=MasterProc.and.DEBUG)
  endif

  if(MasterProc.and..not.first_call) &
    call PrintLog("No external BICs set")
  first_call = .true.
end subroutine set_extbic

end module My_ExternalBICs_ml
