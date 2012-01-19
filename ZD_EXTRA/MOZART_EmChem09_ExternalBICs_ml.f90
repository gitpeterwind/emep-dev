module My_ExternalBICs_ml
!EXTERNAL Boundary and Initial Conditions:
! Here =  IFS-MOZART; used for FORECAST
! Replaces the default My_ExternalBICs_ml
! Note: ONLY for EmChem09-type chemical mechanisms

 use ChemSpecs_adv_ml
 use Io_ml, only : PrintLog
 implicit none
 private

  public :: set_extbic

  logical, public, parameter :: EXTERNAL_BIC_SET  = .true.
  logical, public, parameter :: TOP_BC = .true.
  character(len=30),public, parameter :: EXTERNAL_BIC_NAME = "MOZART"

  type, public :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
    character(len=24) :: varname=""
    logical           :: wanted=.false.,found=.false.
  end type icbc

  type, public :: adv_icbc             ! IC/BC Set, included intended ixadv
    integer           :: ixadv=-1
    type(icbc)        :: icbc
  end type adv_icbc


  type(adv_icbc), dimension(9), public, save :: EXTERNAL_BC

contains

  subroutine set_extbic()

     call PrintLog("Setting MOZART external BICs")
      EXTERNAL_BC=(/adv_icbc(IXADV_O3    ,icbc('O3_VMR_inst'    ,.true.,.false.)), &
                adv_icbc(IXADV_NO    ,icbc('NO_VMR_inst'    ,.true.,.false.)), &
                adv_icbc(IXADV_NO2   ,icbc('NO2_VMR_inst'   ,.true.,.false.)), &
                adv_icbc(IXADV_PAN   ,icbc('PAN_VMR_inst'   ,.true.,.false.)), &
                adv_icbc(IXADV_HNO3  ,icbc('HNO3_VMR_inst'  ,.true.,.false.)), &
                adv_icbc(IXADV_CO    ,icbc('CO_VMR_inst'    ,.true.,.false.)), &
                adv_icbc(IXADV_C2H6  ,icbc('C2H6_VMR_inst'  ,.true.,.false.)), &
                adv_icbc(IXADV_HCHO  ,icbc('CH2O_VMR_inst'  ,.true.,.false.)), &
                adv_icbc(IXADV_CH3CHO,icbc('CH3CHO_VMR_inst',.true.,.false.))/)

  end subroutine set_extbic

end module My_ExternalBICs_ml
