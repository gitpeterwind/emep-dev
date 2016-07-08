module DO3SE_Resistance_ml

  use DO3SE_ModelConstants_ml, only: UNDEF, MAX_LAYERS

  implicit none
  private

  !> Resistance model.
  type, public :: ResistanceModel_t
    real :: Ra = UNDEF    !< Aerodynamic resistance (s m-1) between "decoupled"
                          !! height and the top of the canopy
    real :: Rb = UNDEF    !< Quasi-laminar boundary layer resistance (s m-1)
    real, dimension(MAX_LAYERS) :: &
      Rinc = UNDEF, &     !< In-canopy aerodynamic resistance (s m-1)
      Rext = UNDEF, &     !< External plant cuticle resistance (s m-1)
      Rsto = UNDEF, &     !< Stomatal resistance (s m-1)
      Rsur = UNDEF        !< Surface resistance (s m-1) from each layer downwards
    real :: Rgs = UNDEF   !< Ground surface resistance (s m-1)
    real :: nL = 0        !< Number of layers used
  end type

  public :: Rinc
  public :: Rinc_prototype
  public :: Rext
  public :: ML_Rsur

contains

  !> Estimate in-canopy aerodynamic resistance (Rinc, s m-1).
  !!
  !! This is the older single-layer DO3SE method.
  !!
  !! TODO: to use in a multilayer model, what does h represent?  Height above 
  !!       ground, or thickness of layer?
  pure real function Rinc(SAI, h, ustar)
    real, intent(in) :: SAI   !< Stand area index (m2 m-2)
    real, intent(in) :: h     !< Vegetation height (m)
    real, intent(in) :: ustar !< Friction velocity (m s-1)

    real, parameter :: Rinc_b = 14    ! Rinc coefficient

    Rinc = Rinc_b * SAI * h/ustar
  end function Rinc

  !> Estimate in-canopy aerodynamic resistance (Rinc, s m-1).
  !!
  !! This is the experimental method developed for the Keenley grassland
  !! multilayer model.
  !! 
  !! TODO: decide if we should keep this method
  pure real function Rinc_prototype(SAI, ustar)
    real, intent(in) :: SAI   !< Stand area index (m2 m-2)
    real, intent(in) :: ustar !< Friction velocity (m s-1)

    real, parameter :: Rinc_b = 14    ! Rinc coefficient

    Rinc_prototype = Rinc_b * SAI * Rinc_b/ustar
  end function Rinc_prototype

  !> Estimate external plant cuticle resistance (Rext, s m-1).
  pure real function Rext(SAI)
    real, intent(in) :: SAI   !< Stand area index (m2 m-2)

    real, parameter :: Rext_base = 2500

    Rext = Rext_base / SAI
  end function Rext

  !> Calculate multi-layer Rsur.
  !!
  !! Use per-layer Rsto, Rext and Rinc values to generate layered Rsur values.
  !! Counting the first layer as the top layer, Rsur for a particular layer is
  !! the total resistance for that layer and below, according to the diagram
  !! below.
  !!
  !!    |       Rsto(1)
  !!    |-----<
  !!    |       Rext(1)
  !!    |
  !!    | Rinc(1)
  !!    |
  !!    |       Rsto(2)
  !!    |-----<
  !!    |       Rext(2)
  !!    |
  !!    | Rinc(2)
  !!    |
  !!    | Rgs
  pure function ML_Rsur(Rsto, Rext, Rinc, Rgs) result(Rsur)
    real, dimension(:), intent(in) :: Rsto
    real, dimension(size(Rsto)), intent(in) :: Rext
    real, dimension(size(Rsto)), intent(in) :: Rinc
    real, intent(in) :: Rgs

    real, dimension(size(Rsto)) :: Rsur

    real, dimension(size(Rsto)+1) :: Rsur_tmp
    integer :: i

    Rsur_tmp(size(Rsto)+1) = Rgs
    do i = size(Rsto), 1, -1
      Rsur_tmp(i) = 1.0/(1.0/Rsto(i) + 1.0/Rext(i) + 1.0/(Rinc(i) + Rsur_tmp(i+1)))
    end do
    Rsur = Rsur_tmp(1:size(Rsto))
  end function ML_Rsur

end module DO3SE_Resistance_ml
