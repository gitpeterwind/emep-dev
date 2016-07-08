module DO3SE_Types_ml

  use DO3SE_ModelConstants_ml, only: UNDEF, IUNDEF, MAX_LAYERS
  use DO3SE_Util_ml
#include "DO3SE_Util_ml.h"

  implicit none
  private

  integer :: i

  !> Location properties.
  type, public :: Location_t
    real :: lat = UNDEF     !< Latitude (degrees North)
    real :: lon = UNDEF     !< Longitude (degrees East)
    real :: elev = UNDEF    !< Elevation (m above sea level)
    real :: albedo = UNDEF  !< Surface albedo (fraction)
    real :: Rsoil = UNDEF   !< Soil resistance (s m-1)
  end type Location_t

  !> Land cover season.
  type, public :: Season_t
    integer :: SGS = IUNDEF
    integer :: EGS = IUNDEF
    integer :: Astart = IUNDEF
    integer :: Aend = IUNDEF

    !> LAI method:
    !!    - "input all":      LAI supplied for every layer
    !!    - "input total":    Total LAI supplied as first layer's value,
    !!                        redistributed according to fLAI_layer
    !!    - "day PLF total":  Total LAI estimated from day of year according to
    !!                        a piecewise linear function and redistributed
    !!                        according to fLAI_layer.
    character(len=16) :: LAI_method = "day PLF total"
    ! LAI PLF parameters
    real :: LAI_a = UNDEF       !< LAI value at SGS (m2 m-2)
    real :: LAI_b = UNDEF       !< LAI value at SGS + LAI_1 (m2 m-2)
    real :: LAI_c = UNDEF       !< LAI value at EGS - LAI_2 (m2 m-2)
    real :: LAI_d = UNDEF       !< LAI value at EGS (m2 m-2)
    integer :: LAI_1 = IUNDEF   !< Time from LAI_a to LAI_b (days)
    integer :: LAI_2 = IUNDEF   !< Time from LAI_c to LAI_d (days)
    !> Distribution of LAI (and SAI) between layers
    real, dimension(MAX_LAYERS) :: fLAI_layer = UNDEF
    ! This method would be preferred but triggers several ifort bugs
    !real, dimension(MAX_LAYERS) :: fLAI_layer = (/1.0, (UNDEF, i = 2, MAX_LAYERS)/)

    !> SAI method:
    !!    - "input all":    SAI supplied for every layer
    !!    - "input total":  Total SAI supplied as first layer's value,
    !!                      redistributed according to fLAI_layer
    !!    - "LAI":          SAI = LAI
    !!    - "forest":       Forest method: SAI = LAI + 1
    !!    - "wheat":        Based on wheat lifecycle, requires LAI PLF parameters
    character(len=16) :: SAI_method = "LAI"
  end type Season_t

  !> Species parameters.
  type, public :: Species_t
    real :: h = UNDEF           !< Typical canopy height (m)
    real :: Lm = UNDEF          !< Leaf dimension (m)
  end type Species_t

  public :: check_Season

contains

  subroutine check_Season(season, nL)
    type(Season_t), intent(inout) :: season
    integer, intent(in) :: nL   ! Number of vegetetion layers being used

    ! Sanity-check/fix LAI distribution between layers
    if (nL < MAX_LAYERS) then
      ! Zero unused layers
      season%fLAI_layer(nL:MAX_LAYERS) = 0.0
    end if
    if (nL == 1) then
      ! Automatically set fLAI_layer for single-layer case
      season%fLAI_layer(1) = 1.0
    else
      ! Assert that all layers have been defined in multi-layer case
      ASSERT(all(season%fLAI_layer /= UNDEF))
      ASSERT(sum(season%fLAI_layer) > 0.0)
      ! Normalise layer fractions
      season%fLAI_layer = season%fLAI_layer / sum(season%fLAI_layer)
    end if
  end subroutine check_Season

end module DO3SE_Types_ml
