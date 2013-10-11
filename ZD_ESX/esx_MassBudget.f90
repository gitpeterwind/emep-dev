! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!> MODULE
!! Keeps track of mass
!! especially changes in diffusion solver terms

module esx_MassBudget
  use ChemSpecs, only : NSPEC_TOT
  use esx_Variables, only : ESX_MAXNZ  ! Max no. layers in canopy
  implicit none
  private

  public  :: print_mass

  integer, public :: GRNDDEP=1,  TOPLOSS=2, CANDEP=3, TOT=4

 ! Will need to reconsider dimensions in future. The following is quick-n-dirty

  type, public :: massbudget_t
     real                       :: init = 0.0  ! budget terms, GRNDDEP.. TOT
     real, dimension(5)         :: b    = 0.0  ! budget terms, accumulated GRNDDEP.. TOT
     real, dimension(4)         :: loss = 0.0  ! budget terms, instant. loss, 1-3
     real, dimension(ESX_MAXNZ) :: zloss = 0.0 ! canopy losses
  end type massbudget_t

  type(massbudget_t), public, save, dimension(NSPEC_TOT) :: mass = massbudget_t()

contains

subroutine print_mass(ispec,idt,nz, c, txt)
 integer, intent(in) :: ispec  ! chemical index
 integer, intent(in) :: idt    ! time-step counter
 integer, intent(in) :: nz     ! number of layers
 real, dimension(:), intent(in) :: c
 character(len=*), optional, intent(in) :: txt ! Usually species name

 ! Local
 character(len=20):: label
 real :: pcntdiff


 label = "OUTMASS"
 if ( present(txt) )  label = "OUTMASS "//trim(txt)//" "

 associate ( m=> mass(ispec) )

   pcntdiff = -99.
   if( m%init > 1.0e-10 ) pcntdiff = 100.0*(sum(m%b(:)) - m%init ) / m%init

   print "(a20,i3,a,i5,12(a,es11.3))", label, ispec, &
      " ns ", idt, " Vd:", m%b(GRNDDEP), " Ve:", m%b(TOPLOSS), " can:", m%b(CANDEP), &
      " tot:", m%b(TOT),&
      " %diff: ", pcntdiff, " ctop: ", c(nz)
      !" sum: ",sum(m%b(1:5)), " init:", mass(ind,0)%b(5)
 end associate

end subroutine print_mass


end module esx_MassBudget

