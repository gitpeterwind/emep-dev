module DO3SE_Util_ml

  use iso_fortran_env, only: ERROR_UNIT

  use DO3SE_ModelConstants_ml, only: UNDEF, IUNDEF

  implicit none
  private

  public :: PLF_value

  interface is_def
    module procedure is_def_real
    module procedure is_def_integer
  end interface is_def
  public :: is_def

  abstract interface
    subroutine conditional_error(check, message)
      logical, intent(in) :: check
      character(len=*), intent(in) :: message
    end subroutine conditional_error
  end interface

  procedure(conditional_error), pointer, public :: assert => null()
  procedure(conditional_error), pointer, public :: assert_not => null()
  public :: inverted_assert
  public :: inverted_assert_not
  public :: DO3SE_assert

contains

  !> Calculate the value of a piecewise linear function at a particular x value.
  !!
  !! Given an array of points (2xN) which describe a piecewise linear function,
  !! this function gets the y value for the given x value.  Values before/after
  !! the range of the function are set to the first/last y value respectively.
  !! Within the range of the function, zero-width segments are ignored so that
  !! discontinuous functions can be defined - be aware that the the value at the
  !! shared x value will come from the first of the two segments that share it.
  pure function PLF_value(points, x) result(y)
    real, dimension(:,:), intent(in) :: points
    real, intent(in) :: x
    real :: y

    ! TODO: sanity-check points: should be size 2 in first dimension, and
    !       x values should never decrease.

    integer :: i, n
    real :: ax, ay, bx, by

    n = ubound(points, 2)

    if (x < points(1,1)) then
      y = points(2,1)
    else if (x > points(1,n)) then
      y = points(2,n)
    else
      bx = points(1,1)
      by = points(2,1)

      do i = 2, n
        ax = bx
        ay = by
        bx = points(1,i)
        by = points(2,i)

        ! Skip zero-width pieces (this should be equivalent to an
        ! equality check, but checking floating point equality is evil
        ! and the compiler warns about it)
        if (abs(ax - bx) < epsilon(ax)) then
          continue
        end if

        if (x <= bx) then
          y = ay + (by - ay) * ((x - ax) / (bx - ax))
          exit
        end if
      end do
    end if
  end function PLF_value


  pure logical function is_def_real(x)
    real, intent(in) :: x

    is_def_real = x /= UNDEF
  end function is_def_real

  pure logical function is_def_integer(x)
    integer, intent(in) :: x

    is_def_integer = x /= IUNDEF
  end function is_def_integer


  subroutine inverted_assert(check, message)
    logical, intent(in) :: check
    character(len=*), intent(in) :: message

    call assert(.not. check, message)
  end subroutine inverted_assert

  subroutine inverted_assert_not(check, message)
    logical, intent(in) :: check
    character(len=*), intent(in) :: message

    call assert_not(.not. check, message)
  end subroutine inverted_assert_not

  subroutine DO3SE_assert(check, message)
    logical, intent(in) :: check
    character(len=*), intent(in) :: message

    if (.not. check) then
      write (ERROR_UNIT, '(a)') "ERROR: "//message
      stop
    end if
  end subroutine DO3SE_assert

end module DO3SE_Util_ml
