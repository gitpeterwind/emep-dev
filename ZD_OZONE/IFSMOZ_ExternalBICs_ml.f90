module My_ExternalBICs_ml
!External Boundary and Initial Conditions in FORECAST mode:
! IFS-MOZART BCs for FORECAST mode
use ModelConstants_ml,     only: MasterProc, EXP_NAME, FORECAST, DEBUG=>DEBUG_NEST_ICBC
use CheckStop_ml,          only: CheckStop
use Io_ml,                 only: PrintLog
use TimeDate_ExtraUtil_ml, only: date2string
use ChemSpecs_adv_ml,      only: IXADV_O3,IXADV_NO,IXADV_NO2,IXADV_PAN,&
            IXADV_HNO3,IXADV_CO,IXADV_C2H6,IXADV_HCHO,IXADV_CH3CHO,IXADV_H2O2,&
            IXADV_C5H8,IXADV_C3H6,IXADV_NC4H10,IXADV_OXYL,IXADV_CH4,IXADV_SO2,&
            IXADV_SEASALT_F,IXADV_SEASALT_C,IXADV_DUST_SAH_F,IXADV_DUST_SAH_C,&
            IXADV_FFIRE_OC,IXADV_FFIRE_BC,IXADV_SO4
implicit none
private
public :: set_extbic

logical, public, save :: &
  EXTERNAL_BIC_SET  = .false.

logical, public, parameter :: &
  TOP_BC            = .true.

character(len=*),public, parameter :: &
  EXTERNAL_BIC_NAME = "IFS-MOZART"

! i West/East bnd; j North/South bnd; k Top
integer,save, public :: iw=-1, ie=-1, js=-1, jn=-1, kt=-1 ! i West/East bnd; j North/South bnd; k Top

! YYYY, YY, MM, DD, hh will be replaced by numbers by the program.
! Search for date2string in set_extbic and uncomment lines, if necessary.
! For details, see detail2str in TimeDate_ExtraUtil_ml.f90
character(len=*),private, parameter :: &
  template_read_3D = 'EMEP_IN_IC.nc'         , &
  template_read_BC = 'EMEP_IN_BC_YYYYMMDD.nc', &
  template_write   = 'EMEP_OUT.nc'
! template_write   = 'EMEP_OUT_YYYYMMDD.nc'
character(len=len(template_read_3D)),public, save :: &
  filename_read_3D = template_read_3D
character(len=len(template_read_BC)),public, save :: &
  filename_read_BC = template_read_BC
character(len=len(template_write)),public, save :: &
  filename_write   = template_write

character(len=*),public, parameter :: &
  filename_eta     = 'EMEP_IN_BC_eta.zaxis'

type, public :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
  integer           :: ixadv=-1
  character(len=24) :: varname="none"
  logical           :: wanted=.false.,found=.false.
endtype icbc

type(icbc), dimension(:), public, pointer :: &
  EXTERNAL_BC=>null()

logical, private, parameter :: T=.true.,F=.false.

! BC from IFS-MOZ exp id: f7kn (before 2011-12-01)
type(icbc), dimension(9), private, target :: &
  IFS_MOZ_f7kn=(/icbc(IXADV_O3    ,'O3_VMR_inst'    ,T,F), &
                 icbc(IXADV_NO    ,'NO_VMR_inst'    ,T,F), &
                 icbc(IXADV_NO2   ,'NO2_VMR_inst'   ,T,F), &
                 icbc(IXADV_PAN   ,'PAN_VMR_inst'   ,T,F), &
                 icbc(IXADV_HNO3  ,'HNO3_VMR_inst'  ,T,F), &
                 icbc(IXADV_CO    ,'CO_VMR_inst'    ,T,F), &
                 icbc(IXADV_C2H6  ,'C2H6_VMR_inst'  ,T,F), &
                 icbc(IXADV_HCHO  ,'CH2O_VMR_inst'  ,T,F), &
                 icbc(IXADV_CH3CHO,'CH3CHO_VMR_inst',T,F)/)

! BC from IFS-MOZ exp id: fkya (new since 2011-12-01)
type(icbc), dimension(9), private, target :: &
  IFS_MOZ_fkya=(/icbc(IXADV_O3    ,'O3'    ,T,F), &
                 icbc(IXADV_NO    ,'NO'    ,T,F), &
                 icbc(IXADV_NO2   ,'NO2'   ,T,F), &
                 icbc(IXADV_PAN   ,'PAN'   ,T,F), &
                 icbc(IXADV_HNO3  ,'HNO3'  ,T,F), &
                 icbc(IXADV_CO    ,'CO'    ,T,F), &
                 icbc(IXADV_C2H6  ,'C2H6'  ,T,F), &
                 icbc(IXADV_HCHO  ,'CH2O'  ,T,F), &
                 icbc(IXADV_CH3CHO,'CH3CHO',T,F)/)

! BC for EVA2010: IFS-MOZ (MACC-GRG) & IFS-AER (MACC-AER)
type(icbc), dimension(27), private, target :: &
  IFS_EVA_2010=(/icbc(IXADV_O3        ,'O3'          ,T,F), & ! same as fkya
                 icbc(IXADV_NO        ,'NO'          ,T,F), & ! same as fkya
                 icbc(IXADV_NO2       ,'NO2'         ,T,F), & ! same as fkya
                 icbc(IXADV_PAN       ,'PAN'         ,T,F), & ! same as fkya
                 icbc(IXADV_HNO3      ,'HNO3'        ,T,F), & ! same as fkya
                 icbc(IXADV_CO        ,'CO'          ,T,F), & ! same as fkya
                 icbc(IXADV_C2H6      ,'C2H6'        ,T,F), & ! same as fkya
                 icbc(IXADV_HCHO      ,'CH2O'        ,T,F), & ! same as fkya
                 icbc(IXADV_CH3CHO    ,'CH3CHO'      ,T,F), & ! same as fkya
                 icbc(IXADV_CH4       ,'CH4'         ,T,F), &
                 icbc(IXADV_SO2       ,'SO2'         ,T,F), &
                 icbc(IXADV_H2O2      ,'H2O2'        ,T,F), &
                 icbc(IXADV_C5H8      ,'ISOP'        ,T,F), &
                 icbc(-1              ,'OH'          ,F,F), &
                 icbc(IXADV_NO2       ,'HO2NO2'      ,T,F), &
                 icbc(IXADV_C3H6      ,'BIGENE'      ,T,F), &
                 icbc(IXADV_NC4H10    ,'BIGALK'      ,T,F), &
                 icbc(IXADV_OXYL      ,'TOLUENE'     ,T,F), &
                 icbc(IXADV_SEASALT_F ,'SeaSalt_f'   ,T,F), &
                 icbc(IXADV_SEASALT_C ,'SeaSalt_c'   ,T,F), &
                 icbc(-1              ,'SeaSalt_g'   ,F,F), &
                 icbc(IXADV_DUST_SAH_F,'DesertDust_f',T,F), &
                 icbc(IXADV_DUST_SAH_C,'DesertDust_c',T,F), &
                 icbc(-1              ,'DesertDust_g',F,F), &
                 icbc(IXADV_FFIRE_OC  ,'OC'          ,T,F), &
                 icbc(IXADV_FFIRE_BC  ,'BC'          ,T,F), &
                 icbc(IXADV_SO4       ,'SO4'         ,T,F)/)

contains
subroutine set_extbic(idate)
  implicit none
  integer,intent(in),optional :: idate(4)

  character(len=*), parameter :: &
    DEBUG_FMT="('set_extbic DEBUG: ',A,' ''',A,'''.')"
  character(len=16) :: bctype_name='IFS_MOZ_????'
  integer :: ydmh=0

  ! A date is needed to set this BCs
  call CheckStop(.not.EXTERNAL_BIC_SET.and..not.present(idate),&
    "set_extbic: External BCs from "//EXTERNAL_BIC_NAME//&
    " are DATE depending")
  if(.not.present(idate))return

!--- Set filename from idate: on every call
  if(MasterProc.and.DEBUG) write(*,DEBUG_FMT) &
    "External BICs filenames for",EXTERNAL_BIC_NAME
! filename_read_3D=date2string(template_read_3D,idate,debug=MasterProc.and.DEBUG)
  filename_read_BC=date2string(template_read_BC,idate,debug=MasterProc.and.DEBUG)
! filename_write  =date2string(template_write  ,idate,debug=MasterProc.and.DEBUG)

!--- Set BC type  from idate: on first call only
  if(EXTERNAL_BIC_SET) return
  select case (EXP_NAME)
  case("FORECAST")
   !ydmh=idate(1)*1000000+idate(2)*10000+idate(3)*100+idate(4)
    ydmh=dot_product(idate,nint((/1e6,1e4,1e2,1e0/)))
    select case (ydmh)
    case(:2011113023)         ! Untill 2011-11-30 23:00
      EXTERNAL_BC=>IFS_MOZ_f7kn
      bctype_name='IFS_MOZ_f7kn'
    case(2011120100:)         ! from   2011-12-01 00:00
      EXTERNAL_BC=>IFS_MOZ_fkya
      bctype_name='IFS_MOZ_fkya'
    endselect
  case("EVA2010")             ! EVA2010: GRG & AER
    EXTERNAL_BC=>IFS_EVA_2010
    bctype_name='IFS_EVA_2010'
  case default                ! BC used only on FORECAST mode
    call CheckStop("set_extbic: External BCs from "//EXTERNAL_BIC_NAME//&
      " are only set on FORECAST mode")
  endselect
  if(MasterProc.and.DEBUG) write(*,DEBUG_FMT) &
    date2string("BCs for YYYY-MM-DD hh type",idate),trim(bctype_name)

  if(MasterProc) &
    call PrintLog("External BICs set for "//EXTERNAL_BIC_NAME)
  EXTERNAL_BIC_SET  = .true.
endsubroutine set_extbic

endmodule My_ExternalBICs_ml
