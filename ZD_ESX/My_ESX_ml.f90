module ESX_ml

  use esx_Variables, only: Config_esx, esx
  use esx_Zchem, only: init_Zchem
  use esx_Zgrid, only: init_ZGrid
  use esx_Zveg, only: Config_Zveg

  implicit none
  private

  public :: Init_ESX
  public :: Run_ESX

contains

  subroutine Init_ESX()

    integer :: config_io

    open (newunit=config_io, file="config_esx.nml")
    call Config_esx(config_io, writelog=.true.)
    call init_Zgrid()
    call Config_Zveg(config_io, writelog=.true.)
    call init_Zchem(esx%nz, config_io, debug_level=esx%debug_Zchem)
    close (config_io)
  end subroutine Init_ESX

  subroutine Run_ESX()
  end subroutine Run_ESX

end module ESX_ml
