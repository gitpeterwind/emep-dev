module My_ExternalBICs_ml
!External Boundary and Initial Conditions in FORECAST mode:
! IFS-MOZART BCs for FORECAST mode
use ModelConstants_ml,     only: MasterProc, FORECAST, DEBUG=>DEBUG_NEST_ICBC
use CheckStop_ml,          only: CheckStop
use Io_ml,                 only: PrintLog
use TimeDate_ExtraUtil_ml, only: date2string
use ChemSpecs_adv_ml,      only: IXADV_O3,IXADV_NO,IXADV_NO2,IXADV_PAN,&
                                 IXADV_HNO3,IXADV_CO,IXADV_C2H6,IXADV_HCHO,&
                                 IXADV_CH3CHO
implicit none
private
public :: set_extbic

logical, public, save :: &
  EXTERNAL_BIC_SET  = .false.

logical, public, parameter :: &
  TOP_BC            = .true.

character(len=*),public, parameter :: &
  EXTERNAL_BIC_NAME = "IFS-MOZART"

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
! filename_eta='/global/work/mifapw/emep/Data/MACC02/Boundary_conditions/mozart_eta.zaxis'

type, public :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
  character(len=24) :: varname=""
  logical           :: wanted=.false.,found=.false.
end type icbc

type, public :: adv_icbc             ! IC/BC Set, included intended ixadv
  integer           :: ixadv=-1
  type(icbc)        :: icbc=icbc('',.false.,.false.)
end type adv_icbc

type(adv_icbc), dimension(:), public, pointer :: &
  EXTERNAL_BC=>null()

! BC from IFS-MOZ exp id: f7kn (before 2011-12-01)
type(adv_icbc), dimension(9), private, target :: &
  IFS_MOZ_f7kn=(/adv_icbc(IXADV_O3    ,icbc('O3_VMR_inst'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NO    ,icbc('NO_VMR_inst'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NO2   ,icbc('NO2_VMR_inst'   ,.true.,.false.)), &
                 adv_icbc(IXADV_PAN   ,icbc('PAN_VMR_inst'   ,.true.,.false.)), &
                 adv_icbc(IXADV_HNO3  ,icbc('HNO3_VMR_inst'  ,.true.,.false.)), &
                 adv_icbc(IXADV_CO    ,icbc('CO_VMR_inst'    ,.true.,.false.)), &
                 adv_icbc(IXADV_C2H6  ,icbc('C2H6_VMR_inst'  ,.true.,.false.)), &
                 adv_icbc(IXADV_HCHO  ,icbc('CH2O_VMR_inst'  ,.true.,.false.)), &
                 adv_icbc(IXADV_CH3CHO,icbc('CH3CHO_VMR_inst',.true.,.false.))/)

! BC from IFS-MOZ exp id: fkya (new since 2011-12-01)
type(adv_icbc), dimension(9), private, target :: &
  IFS_MOZ_fkya=(/adv_icbc(IXADV_O3    ,icbc('O3'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NO    ,icbc('NO'    ,.true.,.false.)), &
                 adv_icbc(IXADV_NO2   ,icbc('NO2'   ,.true.,.false.)), &
                 adv_icbc(IXADV_PAN   ,icbc('PAN'   ,.true.,.false.)), &
                 adv_icbc(IXADV_HNO3  ,icbc('HNO3'  ,.true.,.false.)), &
                 adv_icbc(IXADV_CO    ,icbc('CO'    ,.true.,.false.)), &
                 adv_icbc(IXADV_C2H6  ,icbc('C2H6'  ,.true.,.false.)), &
                 adv_icbc(IXADV_HCHO  ,icbc('CH2O'  ,.true.,.false.)), &
                 adv_icbc(IXADV_CH3CHO,icbc('CH3CHO',.true.,.false.))/)

contains
subroutine set_extbic(idate)
  implicit none
  integer,intent(in),optional :: idate(4)

  integer :: ydmh=0
  character(len=16) :: bctype_name='IFS_MOZ_????'

  ! BC used only on FORECAST mode
  call CheckStop(.not.EXTERNAL_BIC_SET.and..not.FORECAST,&
    "set_extbic: External BCs from "//EXTERNAL_BIC_NAME//&
    " are only set on FORECAST mode")
  if(.not.FORECAST)return

  ! A date is needed to set this BCs
  call CheckStop(.not.EXTERNAL_BIC_SET.and..not.present(idate),&
    "set_extbic: External BCs from "//EXTERNAL_BIC_NAME//&
    " are DATE depending")
  if(.not.present(idate))return

  ! Set filename from idate: on every call
  if(MasterProc.and.DEBUG)&
    write(*,"(A,':',1X,A,1X,'''',A,'''.')")"set_extbic DEBUG",&
      "External BICs filenames for",EXTERNAL_BIC_NAME
! filename_read_3D=date2string(trim(template_read_3D),idate,&
!   debug=MasterProc.and.DEBUG)
  filename_read_BC=date2string(trim(template_read_BC),idate,&
    debug=MasterProc.and.DEBUG)
! filename_write=date2string(trim(template_write),idate,&
!   debug=MasterProc.and.DEBUG)

  ! Set BC type  from idate: on first call only
  if(EXTERNAL_BIC_SET) return
  ydmh=idate(1)*1000000+idate(2)*10000+idate(3)*100+idate(4)
  select case (ydmh)
    case(:2011113023)         ! Untill 2011-11-30 23:00
      EXTERNAL_BC=>IFS_MOZ_f7kn
      bctype_name='IFS_MOZ_f7kn'
    case(2011120100:)         ! from   2011-12-01 00:00
      EXTERNAL_BC=>IFS_MOZ_fkya
      bctype_name='IFS_MOZ_fkya'
  endselect
  if(MasterProc.and.DEBUG)&
    write(*,"(A,':',1X,A,1X,'''',A,'''.')")"set_extbic DEBUG",&
    date2string("BCs for YYYY-MM-DD hh type",idate),trim(bctype_name)

  if(MasterProc) &
    call PrintLog("External BICs set for "//EXTERNAL_BIC_NAME)
  EXTERNAL_BIC_SET  = .true.
end subroutine set_extbic

end module My_ExternalBICs_ml
