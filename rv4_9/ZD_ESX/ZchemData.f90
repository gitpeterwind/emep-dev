!> MODULE  ZchemData
!! Storage of column arrays for chemical scheme
!!

 module ZchemData
  use AllocInits, only : AllocInit
  use ChemSpecs,  only : nspec => NSPEC_TOT
  use DefPhotolysis, only : MAXRCPHOT
  implicit none
  private

  public :: Alloc1Dchem       !> Sets up arrays of dimension nspec x nz


  !> 1-D chemical variables

  real, public, save, allocatable, dimension(:,:) :: & ! dims: nspec x nz
      rct               & !> chemical rate coeffs.  molec/cm3/s
     ,rcemis            & !> Emissions rate coeff.  molec/cm3/s
     ,DChem             & !> chemical tendency, molec/cm3/s
     ,rcphot              !> Photolysis rates

  ! We want pointers to this later:
  real, public, save, target,  allocatable, dimension(:,:) :: & ! dims: nspec x nz
      xChem               !> concentrations, molec/cm3

  real, public, save, allocatable, dimension(:) :: xSO4, xNO3, xNH4 ! for Riemer
 
  contains
  !--------------------------------------------------------------------------!
   subroutine Alloc1Dchem(nz, nrct, debug_level)
    integer, intent(in) :: nz, nrct, debug_level

    !> We need to re-allocate if changes in array size, e.g. moving from
    !! EMEP to ESX

      if( debug_level > 0 ) print *, "ALLOCATE xCHEM?? ", nspec, nz, size(xChem, 1)

      if ( .not. allocated( xChem) ) then
          if( debug_level > 0 ) print *, "ALLOCATES xCHEM ", nspec, nz
      end if

      call AllocInit( xChem,  0.0, nspec, nz, "set1D:xChem")
      call AllocInit( Dchem,  0.0, nspec, nz, "set1D:Dchem")
      call AllocInit( rcemis, 0.0, nspec, nz, "set1D:rcemis")
      call AllocInit( rct,    0.0,  nrct, nz, "set1D:rct")
      call AllocInit( rcphot, 0.0,  MAXRCPHOT, nz, "set1D:rct")

      call AllocInit( xSO4,    0.0, nz, "set1D:xSO4")
      call AllocInit( xNO3,    0.0, nz, "set1D:xNO3")
      call AllocInit( xNH4,    0.0, nz, "set1D:xNH4")

  end subroutine Alloc1Dchem

  !--------------------------------------------------------------------------!
 end module ZchemData

