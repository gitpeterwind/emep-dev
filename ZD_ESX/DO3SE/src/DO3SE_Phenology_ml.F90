module DO3SE_Phenology_ml

  use DO3SE_Types_ml, only: Season_t
  use DO3SE_Util_ml

  implicit none
  private

  public :: LAI_day_PLF
  public :: SAI_wheat

contains

  !> Estimate LAI based on a piecewise linear function of the day of year.
  !! The function wraps around to interpolate between EGS and SGS when the
  !! LAI does not conform to just a simple "bump" in the middle of the
  !! year.
  pure real function LAI_day_PLF(s, dd)
    type(Season_t), intent(in) :: s   !< Season configuration
    integer, intent(in) :: dd         !< Day of year

    LAI_day_PLF = PLF_value(reshape((/ real :: &
      (s%EGS - 365), s%LAI_d, &
      s%SGS, s%LAI_a, &
      (s%SGS + s%LAI_1), s%LAI_b, &
      (s%EGS - s%LAI_2), s%LAI_c, &
      s%EGS, s%LAI_d, &
      (365 + s%SGS), s%LAI_a &
      /),  (/2, 6/)), real(dd))
  end function LAI_day_PLF


  !> Wheat stand area index based on the growing season.
  pure real function SAI_wheat(s, dd, LAI)
    type(Season_t), intent(in) :: s   !< Season configuration
    integer, intent(in) :: dd         !< Day of year
    real, intent(in) :: LAI           !< Leaf area index (m2 m-2)

    if (dd < s%SGS .or. dd > s%EGS) then
      SAI_wheat = LAI
    else if (dd < s%SGS + s%LAI_1) then
      ! TODO: why not = LAI * (10.0/7.0) ?
      SAI_wheat = LAI + ((5.0/3.5) - 1) * LAI
    else
      SAI_wheat = LAI + 1.5
    end if
  end function SAI_wheat

end module DO3SE_Phenology_ml
