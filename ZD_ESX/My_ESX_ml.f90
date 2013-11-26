module ESX_ml

  use esx_Variables, only: Config_esx

  implicit none
  private

  public :: Init_ESX
  public :: ESX

contains

  subroutine Init_ESX()

    integer :: config_io

    open (newunit=config_io, file="config_esx.nml")
    call Config_esx(config_io, writelog=.false.)
    close (config_io)
  end subroutine Init_ESX

  subroutine ESX()
  end subroutine ESX

end module ESX_ml
