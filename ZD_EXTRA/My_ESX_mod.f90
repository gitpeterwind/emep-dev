!> Dummy implementation of ESX_ml for models that don't use ESX.
module ESX_ml

  implicit none
  private

  public :: Init_ESX
  public :: Run_ESX

contains

  subroutine Init_ESX()
  end subroutine Init_ESX

  subroutine Run_ESX()
  end subroutine Run_ESX

end module ESX_ml

