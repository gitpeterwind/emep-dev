module DO3SE_GstoTypes_ml

  use DO3SE_Util_ml
#include "DO3SE_Util_ml.h"
  use DO3SE_ModelConstants_ml, only: UNDEF, IUNDEF

  implicit none
  private

  !> Multiplicative stomatal conductance parameters.
  type, public :: GstoParams_t
    real :: fmin = UNDEF        !< Minimum stomatal conductance (fraction)
    real :: gmax = UNDEF        !< Maximum stomatal conductance (mmol O3 m-2 PLA s-1)
    real :: gmorph = 1.0        !< Sun/shade morphology factor (fraction)

    real :: f_phen = 1.0        !< Phenology-related effect on gsto (fraction)
    real :: leaf_f_phen = 1.0   !< Phenology-related effect on leaf gsto (fraction)
    real :: f_light = 1.0       !< Irradiance effect on gsto (fraction)
    real :: leaf_f_light = 1.0  !< Irradiance effect on leaf gsto (fraction)
    real :: f_temp = 1.0        !< Temperature effect on gsto (fraction)
    real :: f_VPD = 1.0         !< VPD effect on gsto (fraction)
    real :: f_SW = 1.0          !< Soil water effect on gsto (fraction)
    real :: f_O3 = 1.0          !< O3 effect on gsto (fraction)
  end type GstoParams_t

  !> Configuration for multiplicative stomatal conductance functions.
  type, public :: GstoConfig_t
    real :: fmin = UNDEF    !< Minimum stomatal conductance (fraction)
    real :: gmax = UNDEF    !< Maximum stomatal conductance (mmol O3 m-2 PLA s-1)
    real :: gmorph = 1.0    !< Sun/shade morphology factor (fraction)

    !> Stomatal conductance method:
    !!    - "multiplicative":   Use DO3SE multiplicative model
    !!    - "photosynthesis":   Use Farquar-based photosynthesis model (hybrid
    !!                          with some multiplicative components)
    character(len=16) :: method = "multiplicative"

    !> f_phen method:
    !!    - "disabled":         f_phen supplied (or left at default value of 1.0)
    !!    - "simple day PLF":   "single hump" method using 3 values (_a, _c and
    !!                          _e) and 2 slope periods (_1 and _4)
    !!    - "complex day PLF":  "double hump" method using 5 values (_a to _e)
    !!                          and 4 slopes: _1 and _4 for start and end, _2
    !!                          and _3 for the middle section, starting at _limA
    !!                          and ending at _limB
    character(len=16) :: f_phen_method = "simple"
    integer :: f_phen_limA = IUNDEF   !< Start of soil water limitation
    integer :: f_phen_limB = IUNDEF   !< End of soil water limitation
    real :: f_phen_a = UNDEF          !< f_phen at SGS
    real :: f_phen_b = UNDEF
    real :: f_phen_c = UNDEF
    real :: f_phen_d = UNDEF
    real :: f_phen_e = UNDEF          !< f_phen at EGS
    integer :: f_phen_1 = IUNDEF      !< Time from f_phen_a to f_phen_b (days)
    integer :: f_phen_2 = IUNDEF      !< Time from f_phen_b to f_phen_c (days)
    integer :: f_phen_3 = IUNDEF      !< Time from f_phen_c to f_phen_d (days)
    integer :: f_phen_4 = IUNDEF      !< Time from f_phen_d to f_phen_e (days)

    !> leaf_f_phen method:
    !!    - "disabled": leaf_f_phen supplied (or left at default value of 1.0)
    !!    - "f_phen":   leaf_f_phen = f_phen
    !!    - "day PLF:   "single hump" PLF between Astart and Aend
    character(len=16) :: leaf_f_phen_method = "f_phen"
    real :: leaf_f_phen_a = UNDEF       !< f_phen at Astart
    real :: leaf_f_phen_b = UNDEF       !< f_phen at mid-season peak
    real :: leaf_f_phen_c = UNDEF       !< f_phen at Aend
    integer :: leaf_f_phen_1 = IUNDEF   !< Time from _a to _b (days)
    integer :: leaf_f_phen_2 = IUNDEF   !< Time from _b to _c (days)

    !> f_light method:
    !!    - "disabled":   f_light and leaf_f_light supplied (or left at default
    !!                    value of 1.0)
    !!    - "enabled":    f_light and leaf_f_light calculated
    character(len=16) :: f_light_method = "enabled"
    real :: f_lightfac = 0.006  !< Single leaf f_light coefficient
    real :: cosA = 0.5          !< cos(A), A = mean leaf inclination (0.5 = 60 degrees)

    !> f_temp method:
    !!    - "disabled":       f_temp supplied (or left at default value of 1.0)
    !!    - "default":        Normal bell-shaped function over T_min -> T_opt -> T_max
    !!    - "square high":    Same as "default", but straight lines from
    !!                        (T_opt, 1.0) -> (T_max, 1.0) -> (T_max, 0.0)
    character(len=16) :: f_temp_method = "default"
    real :: T_min = UNDEF     !< Minimum temperature (degrees C)
    real :: T_opt = UNDEF     !< Optimum temperature, for max. gsto (degrees C)
    real :: T_max = UNDEF     !< Maximum temperature (degrees C)

    !> f_VPD method:
    !!    - "disabled":   f_VPD supplied (or left at default value of 1.0)
    !!    - "linear":     Linear f_VPD relationship between VPD_max and VPD_min
    !!    - "log":        Simple, unparameterised, logarithmic relationship
    character(len=16) :: f_VPD_method = "linear"
    real :: VPD_max = UNDEF     !< VPD for maximum gsto (kPa)
    real :: VPD_min = UNDEF     !< VPD for minimum gsto (kPa)

    !> f_SW method:
    !!    - "disabled":     f_SW supplied (or left at default value of 1.0)
    !!    - "fSWP exp":     Use fSWP exponential curve (see fSWP_exp_curve)
    !!    - "fSWP linear":  Use linear fSWP function (see SWP_min and SWP_max)
    !!    - "fLWP exp":     Use fSWP exponential curve, but with LWP instead of SWP
    !!    - "fPAW":         Use fPAW relationship
    character(len=16) :: f_SW_method = "disabled"

    ! fSWP linear parameters:
    real :: SWP_min = UNDEF   !< SWP for minimum gsto (MPa)
    real :: SWP_max = UNDEF   !< SWP for maximum gsto (MPa)

    !> fSWP exponential curve:
    !!    - "custom":         min(1, max(fmin, a * (-SWP)**b))
    !!    - "temperate":      a = 0.355, b = -0.706
    !!    - "mediterranean":  a = 0.619, b = -1.024
    character(len=16) :: fSWP_exp_curve = "temperate"
    real :: fSWP_exp_a = UNDEF
    real :: fSWP_exp_b = UNDEF

  end type GstoConfig_t

  public :: check_GstoConfig

contains

  subroutine check_GstoConfig(gc)
    type(GstoConfig_t), intent(inout) :: gc

    ASSERT_DEFINED(gc%fmin)
    ASSERT_DEFINED(gc%gmax)

    ! Disable multiplicative gsto components that photosynthesis gsto replaces
    if (gc%method == "photosynthesis") then
      gc%f_light_method = "disabled"
      gc%f_temp_method = "disabled"
      gc%f_VPD_method = "disabled"
    end if

    ! Ensure necessary f_phen parameters are present
    if (gc%f_phen_method /= "disabled") then
      ASSERT_DEFINED(gc%f_phen_a)
      ASSERT_DEFINED(gc%f_phen_c)
      ASSERT_DEFINED(gc%f_phen_e)
      ASSERT_DEFINED(gc%f_phen_1)
      ASSERT_DEFINED(gc%f_phen_4)
      if (gc%f_phen_method == "complex day PLF") then
        ASSERT_DEFINED(gc%f_phen_b)
        ASSERT_DEFINED(gc%f_phen_d)
        ASSERT_DEFINED(gc%f_phen_limA)
        ASSERT_DEFINED(gc%f_phen_limB)
        ASSERT_DEFINED(gc%f_phen_2)
        ASSERT_DEFINED(gc%f_phen_3)
      end if
    end if

    ! Ensure necessary leaf_f_phen parameters are present
    if (gc%leaf_f_phen_method == "day PLF") then
      ASSERT_DEFINED(gc%leaf_f_phen_a)
      ASSERT_DEFINED(gc%leaf_f_phen_b)
      ASSERT_DEFINED(gc%leaf_f_phen_c)
      ASSERT_DEFINED(gc%leaf_f_phen_1)
      ASSERT_DEFINED(gc%leaf_f_phen_2)
    end if

    ! Ensure necessary f_temp parameters are present
    if (gc%f_temp_method /= "disabled") then
      ASSERT_DEFINED(gc%T_min)
      ASSERT_DEFINED(gc%T_opt)
      ASSERT_DEFINED(gc%T_max)
    end if

    ! Ensure necessary f_VPD parameters are present
    if (gc%f_VPD_method == "linear") then
      ASSERT_DEFINED(gc%VPD_max)
      ASSERT_DEFINED(gc%VPD_min)
    end if

    ! Parameterise fSWP exponential curve
    select case (gc%fSWP_exp_curve)
    case ("custom")
      ! Nothing to do
    case ("temperate")
      gc%fSWP_exp_a = 0.355
      gc%fSWP_exp_b = -0.706
    case ("mediterranean")
      gc%fSWP_exp_a = 0.619
      gc%fSWP_exp_b = -1.024
    case default
      UNKNOWN_STRING(gc%fSWP_exp_curve)
    end select

    ! Check that appropriate f_SW parameters are defined
    select case (gc%f_SW_method)
    case ("fSWP exp", "fLWP exp")
      ASSERT_DEFINED(gc%fSWP_exp_a)
      ASSERT_DEFINED(gc%fSWP_exp_b)
    case ("fSWP linear")
      ASSERT_DEFINED(gc%SWP_min)
      ASSERT_DEFINED(gc%SWP_max)
    end select
  end subroutine check_GstoConfig

end module DO3SE_GstoTypes_ml
