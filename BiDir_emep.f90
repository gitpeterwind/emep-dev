module BiDir_emep
  ! DUMMY
  ! Will act as interface between emep-ctm and BiDir_module

  implicit none
  private

  public :: BiDir_ijInit
  public :: BiDir_ijRGs
  public :: BiDir_ijFluxes
  public :: BiDir_ijFinish 

contains
 subroutine BiDir_ijInit(i,j,NH3_ix)
     integer, intent(in) :: i,j,NH3_ix
 end subroutine BiDir_ijInit
 subroutine BiDir_ijRGs(ncall,iL,Rsur,Gsto)
     integer, intent(in):: ncall, iL
     real, intent(out) :: Rsur
     real, intent(in)  :: Gsto
 end subroutine BiDir_ijRGs
 subroutine BiDir_ijFluxes(i,j,iL,Vg_ref,Vg_eff,Rb,Rsur,Gsto)
     integer, intent(in) :: i,j,iL
     real, intent(in) :: Vg_ref,Vg_eff,Rb,Rsur,Gsto
 end subroutine BiDir_ijFluxes
 subroutine BiDir_ijFinish(i,j,gradfac,sumLand,DepLoss)
    integer, intent(in) :: i,j
    real, intent(inout) :: gradfac
    real, intent(in)    :: sumLand, DepLoss
 end subroutine BiDir_ijFinish

end module BiDir_emep
