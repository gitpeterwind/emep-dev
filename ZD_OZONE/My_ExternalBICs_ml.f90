module My_ExternalBICs_ml
! External Boundary and Initial Conditions
! are set from the value depending on Experimnt Name (EXP_NAME)
! Nothing in this file needs to be used if EXTERNAL_BIC_SET  = .false.
use ModelConstants_ml,     only: MasterProc, DEBUG=>DEBUG_NEST_ICBC
use CheckStop_ml,          only: CheckStop
use Io_ml,                 only: PrintLog
use TimeDate_ExtraUtil_ml, only: date2string
implicit none

private
public :: update_bicname,set_extbic

logical, public, parameter :: &
  EXTERNAL_BIC_SET  = .false.,&
  TOP_BC = .false.

character(len=30),public, parameter :: &
  EXTERNAL_BIC_NAME = "DUMMY"

! i West/East bnd; j North/South bnd; k Top
integer,save, public :: iw=-1, ie=-1, js=-1, jn=-1, kt=-1 ! i West/East bnd; j North/South bnd; k Top

character(len=*),private, parameter :: &  
  template_read_3D = 'EMEP_IN.nc', &      ! a different path can be set here
  template_read_BC = 'EMEP_IN.nc', &      ! for each of the IO IC/BC files,
  template_write   = 'EMEP_OUT.nc'        ! if needed.
! template_read_3D = 'EMEP_IN_IC.nc'         , & ! YYYY, YY, MM, DD, hh strings
! template_read_BC = 'EMEP_IN_BC_YYYYMMDD.nc', & ! will be replaced by numbers
! template_write   = 'EMEP_OUT_YYYYMMDD.nc'      ! on set_extbic.

character(len=len(template_read_3D)),public, save :: &
  filename_read_3D = template_read_3D   ! overwritten in update_bicname
character(len=len(template_read_BC)),public, save :: &
  filename_read_BC = template_read_BC   ! overwritten in update_bicname
character(len=len(template_write)),public, save :: &
  filename_write   = template_write     ! overwritten in update_bicname

character(len=*),public, parameter :: &
  filename_eta     = 'EMEP_IN_BC_eta.zaxis'

type, public :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
  integer           :: ixadv=-1
  character(len=24) :: varname="none"
  real              :: frac=1.0
  logical           :: wanted=.not.EXTERNAL_BIC_SET,found=.false.!default is all components for non-external
endtype icbc

type(icbc), dimension(:), public, pointer :: &
  EXTERNAL_BC=>null()

character(len=*),private, parameter :: &
  DEBUG_FMT="(A,' DEBUG: ',A,' ''',A,'''.')"

contains

subroutine update_bicname(idate)
  integer,intent(in) :: idate(4)
!--- Set filename from idate: on every call
  if(MasterProc.and.DEBUG) write(*,DEBUG_FMT),"update_bicname", &
    "External BICs filenames for",EXTERNAL_BIC_NAME
  filename_read_3D=date2string(template_read_3D,idate,debug=DEBUG.and.MasterProc)
  filename_read_BC=date2string(template_read_BC,idate,debug=DEBUG.and.MasterProc)
  filename_write  =date2string(template_write  ,idate,debug=DEBUG.and.MasterProc)
endsubroutine update_bicname

subroutine set_extbic(idate)
  integer,intent(in) :: idate(4) ! Needed on other versions of My_ExternalBICs_ml
  logical, save :: first_call=.true.

  if(first_call) return
  call PrintLog("No external BICs set",MasterProc)
  first_call = .false.
endsubroutine set_extbic

endmodule My_ExternalBICs_ml
