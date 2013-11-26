!> Dummy implementation of ESX_ml for models that don't use ESX.
module ESX_ml

  implicit none
  private

  public :: Init_ESX
  public :: ESX

contains

  subroutine Init_ESX()
  end subroutine Init_ESX

  subroutine ESX()
  end subroutine ESX

end module ESX_ml

