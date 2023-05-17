module BiDir_module
  ! DUMMY 
  implicit none
  private

  character(len=*), parameter :: BiDir_module_status='TOBEDONE'

  ! Main:
  !public  :: BiDirXconcs
  !public  :: BiDirFluxes

  !public  :: BiDirResistances
  !public :: BiDirXwaterEuro
  !public :: BiDirXwaterOrig

  type, public :: BiDir_t

    logical :: EuroXwater = .false.
    logical :: OrigXwater = .false.
    logical :: skipForestDisp   = .false. 
    character(len=20) :: Method = 'NOTSET'
    ! allow for long file names
    character(len=500):: InputFile = 'NOTSET'
    character(len=500) :: InputDir  = 'NOTSET'
  end type BiDir_t
  type(BiDir_t), public, save :: BiDir= BiDir_t()

end module BiDir_module
