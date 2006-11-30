module Dates_ml
!-----------------------------------------------------------------------------
!+ defines numbers associated with dates, e.g. days per month, leap-years. 
!  - aslo "date" type and date "increment" function 
!-----------------------------------------------------------------------------
 implicit none
 private

 !/ subroutines:
 public  :: Init_nmdays   ! sets number of days per month, year, calls is_leap
 public  :: dayno        ! Calculates day number of input date

 !/ functions:
 public  :: is_leap      ! returns 1 for leap_year
 public  :: add_dates    ! adds two dates

 !--  Date type
 type, public :: date
  integer :: year
  integer :: month
  integer :: day
  integer :: hour
  integer :: seconds
 end type date

 interface add_dates
   module procedure add_2dates       ! To add two date-types
   module procedure add_dates_s      ! To add an integer No. of seconds to date
 end interface

 ! and define some useful increments:
  type(date), public, parameter :: ONE_YEAR   = date( 1, 0, 0, 0, 0)
  type(date), public, parameter :: ONE_MONTH  = date( 0, 1, 0, 0, 0)
  type(date), public, parameter :: ONE_DAY    = date( 0, 0, 1, 0, 0)
  type(date), public, parameter :: ONE_HOUR   = date( 0, 0, 0, 1, 0)
  type(date), public, parameter :: ONE_SECOND = date( 0, 0, 0, 0, 1)



 integer, public, dimension(12), save ::  nmdays       ! No. days per month
 integer, public,                save ::  nydays       ! No. days per year
 integer, public, dimension(12), save ::  mday1  ! Day no. of 1st day of month
 integer, public,                save ::  daynumber    ! Day no. (1st jan=1).


 integer, public, parameter :: NMONTHS  = 12   !  No. months per year (was NM)


!.. set start dates for years 1980-2001 ( Y2K compliant!) 

  integer, public, parameter :: MON=1, TUE=2, WED=3, THU=4,  &
                                FRI=5, SAT=6, SUN=7

  integer, public, parameter :: FIRST_SDYEAR = 1980, LAST_SDYEAR = 2006
  integer, public, parameter, dimension(FIRST_SDYEAR:LAST_SDYEAR) ::  &
    STARTDAY= (/ TUE,THU,FRI,SAT,SUN,TUE,WED,THU,FRI,SUN, &
    !...          80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 
                 MON,TUE,WED,FRI,SAT,SUN,MON,WED,THU,FRI,&
    !...          90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 
                  SAT,MON,TUE,WED,THU,SAT,SUN  /)
    !...          00, 01, 02, 03, 04, 05 , 06


  logical, private, save :: Init_done = .false.
contains
!---------------------------------------------------------------------
 subroutine Init_nmdays(indate)
  !/ ndm is number of days per month, for a non-leap year. This is modified
  !  for leap years by the subroutine is_leap(), called on start-up??

 type(date), intent(in) ::  indate
 !/..local
 integer             ::  yy, mm

  yy = indate%year
  if ( yy < 1900 .and.  yy > 50 ) yy = 1900 + indate%year
  if ( yy <   50 )                yy = 2000 + indate%year


   nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/) 
   nmdays(2) = 28  + is_leap(yy)
   nydays    = 365 + is_leap(yy)

   mday1(1) = 1  ! Jan 1st.
   do mm = 2, 12
      mday1(mm) = mday1(mm-1) + nmdays(mm-1) 
   end do

   Init_done = .true.     !/** so that other functions know !

 end subroutine Init_nmdays
!---------------------------------------------------------------------
!---------------------------------------------------------------------
 function is_leap(year) result (is)
 !-- a function to evaluate if year is leap or not, returns
 ! 1 is leap, 0 if not.

  integer, intent(in) :: year
  integer             :: iyear, is
  iyear = year

    if (iyear < 100 ) then
         ! print *,"Error in is_leap, needs 4 digit year!"
       iyear = year + 1900
    endif
    is = 0
    if (modulo(iyear,4)    ==  0 ) is = 1
    if (modulo(iyear,100)  ==  0 ) is = 0
    if (modulo(iyear,400)  ==  0 ) is = 1

  end function is_leap
!---------------------------------------------------------------------
!---------------------------------------------------------------------
  subroutine dayno(month,day,nday)
!---------------------------------------------------------------------
! this subroutine calculates the dayno of a particular date, starting
! with 1st jan as day number one . Leap years have already been 
! accounted for when nmdays was set.

!/ arguments
      integer, intent(in)  ::  month,day
      integer, intent(out) ::  nday      ! Day No.

!/ Local
      integer :: i, imon   

! the day number
      nday = 0
      imon = month - 1

      do i=1, imon
         nday = nday + nmdays(i)
      end do  ! i

      nday = nday + day

      end subroutine dayno

    !---------------------------------------------------------------------
      function add_2dates(date_1,date_2) result(newdate)
    !---------------------------------------------------------------------
    ! A function to add two dates, and return the new date. usually
    ! used to add some time-step to the current date. Only designed
    ! for date_2 values less than or equal to one month.
    !---------------------------------------------------------------------
        type(date), intent(in) ::  date_1, date_2
        type(date)             ::  newdate
        integer :: h_add, d_add, y_add

        if ( .not. Init_done .or.  date_2%month > 1 &
                             .or.  date_2%year > 0 ) then
            print *, " ERROR: Invalid date_2 in add_dates"
            newdate = date( -999, -999, -999, -999 , -999 )
            return 
        end if 

        h_add = 0
        d_add = 0
        y_add = 0

        newdate%seconds = date_1%seconds + date_2%seconds
        if ( newdate%seconds >= 3600 ) then
             newdate%seconds = newdate%seconds - 3600
             h_add = 1
        endif

        newdate%hour = date_1%hour + date_2%hour + h_add
        if ( newdate%hour >= 24 ) then
             newdate%hour = newdate%hour-24
             d_add = 1
        endif

       !/ we need to be careful now about ending up with the right month to 
       !/ get the nmdays value. We add the months first here.

        newdate%day   = date_1%day   + date_2%day   + d_add
        newdate%month = date_1%month + date_2%month

        if ( newdate%month > 12 ) then
            newdate%month = newdate%month - 12
            y_add = 1
        endif

        if ( newdate%day > nmdays( newdate%month )  ) then
             newdate%day = newdate%day - nmdays( newdate%month )
             newdate%month = newdate%month + 1
        endif

        if ( newdate%month > 12 ) then   ! pwhk again, test that newdate%month is not >12
             newdate%month = newdate%month - 12
             y_add = y_add + 1
        endif

        ! Not elegant, but we need to repeat to catch the possible
        ! addition of 31days+31days:

        if ( newdate%day > nmdays( newdate%month )  ) then
             newdate%day = newdate%day - nmdays( newdate%month )
             newdate%month = newdate%month + 1
        endif

        ! phew, got through that mess. Just the easy bits now...

        if ( newdate%month > 12 ) then   ! again, test
             newdate%month = newdate%month - 12
             y_add = y_add + 1
        endif

        newdate%year = date_1%year + date_2%year + y_add
    end function add_2dates

    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
      function add_dates_s(date_1,s) result(newdate)
    !---------------------------------------------------------------------
    ! A function to add a time in seconds to a date
    ! and return the new date. usually
    ! used to add some time-step to the current date. Only designed
    ! for s values less or equal to one hour (s<=3600).
    !---------------------------------------------------------------------
        type(date), intent(in) ::  date_1
        integer, intent(in)    ::  s          ! No. seconds
        type(date)             ::  newdate

        if ( .not. Init_done .or. s>3600) then
            print *, " ERROR: Invalid in add_dates_s",s
            newdate = date( -999, -999, -999, -999 , -999 )
            return 
        end if 
        if ( s < 0 ) then
            print *, " ERROR add_dates_s: Can't add negative time in seconds"
            newdate = date( -999, -999, -999, -999 , -999 )
            return 
        end if 

        newdate%seconds = date_1%seconds + s
        newdate%hour = date_1%hour
        newdate%day = date_1%day
        newdate%month = date_1%month
        newdate%year = date_1%year
        if ( newdate%seconds < 3600 ) return

! we reached next hour

        newdate%seconds = newdate%seconds - 3600
        newdate%hour = newdate%hour + 1
        if ( newdate%hour < 24 ) return

! we reached next day

        newdate%hour = 0
        newdate%day   = newdate%day + 1
        if ( newdate%day <= nmdays( newdate%month )  ) return

! we reached next month

        newdate%day = 1
        newdate%month = newdate%month + 1
        if ( newdate%month <= 12 ) return

! we reached next year

        newdate%month = 1
        newdate%year = newdate%year + 1

    end function add_dates_s

end module DATES_ML

