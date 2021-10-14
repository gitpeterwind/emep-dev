!###############################################################################
!
! CSO_DateTime  -  module to manipulate date structures
!
! TYPES
!
!  A structure is provided to store a date:
!
!   ! declare date types:
!   type(T_CSO_DateTime)     ::  t0, t1, t, dt
!   type(T_CSO_TimeDelta) ::  dt
!
!  with fields:
!
!    character(len=32)   ::  calendar     ! see 'CALENDARS'
!
!    integer             ::  year, month, day, hour, minute, sec, mili
!
!    integer             ::  zone       ! minutes; add this to obtain GMT
!
!
! CALENDARS
!
!   A number of different calendar types is supported:
!
!     'proleptic_gregorian', 
!       'gregorian', 
!       'standard'             : Gregorian calendar, some years have a Februari 29
!     '366_day'                : every year has a Februari 29
!     '365_day'                : a year has never a  Februari 29
!     '360_day'                : every month has 30 days
!
!     'incr'   : incremental time step: year=0, month=0, day >= 0
!
!   The 'incr' type is a special calendar which has no year
!    or month but might have any number of days.
!   Note that day==1 has the interpretation of 24 hours for an 'incr',
!   but means 'first' or 0 hours for one of the regular calendars.
!
!   Use the calendar '360_day' if only operations on years and months are required.
!
!
! CREATING DATE STRUCTURES
!
!   To initialize a new date structure, a few routines are available.
!
!   Use routine 'CSO_DateTime' to initialize some fields and to fill
!   the rest with zero's. If no calendar is specified,
!   the default value 'gregorian' is used (see also DEFAULTS).
!
!     t = CSO_DateTime( calendar='gregorian', year=2000, month=1, ... )
!
!   Use routine 'CSO_TimeDelta' to create a new increment:
!
!     dt = CSO_TimeDelta( year=2000, month=1 )
!
!   Fill the time from the system clock in a date structure:
!
!     t = go_SystemDate()
!
!   Read from line:
!     call subroutine goReadFromLine( '2001-02-03 04:05:06;xxx', t, status, &
!                                         sep=';', calendar='default' )
!
!
! FIELD MANIPULATION
!
!   Use 'Set' to fill some specific fields of a date structure.
!   Special arrays:
!     time4 = (/year,month,day,hour/)
!     time5 = (/year,month,day,hour,minute/)
!     time6 = (/year,month,day,hour,minute,sec/)
!   Example:
!
!     call Set( t [,year=2000] [,month=1] [,day=2] ... &
!                 [,time4=time4] [,time5=time5] [,time6=time6])
!
!   Use 'Get' to obtain some specific fields of a date structure.
!
!     call Get( t [,year=year] [,month=month] ... &
!                 [,time4=time4] [,time5=time5] [,time6=time6] )
!
!   Check contents of a date structure:
!
!     call Check( t )
!
!   Normalize hours to {0,..,23}, minutes to {0,..,59}, etc:
!
!     call Normalize( t )
!
!   Return minimum or maximum of two dates:
!
!     tmin = min( t1, t2 )
!     tmax = max( t1, t2 )
!
!   Retrun begin|end of year|month|day specified by t :
!
!     t1 = Get_Begin_Of( t, 'year'  )    ! 00:00 at 01 Jan of this year
!     t1 = Get_Begin_Of( t, 'month' )    ! 00:00 at first day of this month
!     t1 = Get_Begin_Of( t, 'day'   )    ! 00:00 at this day
!
!     t2 = Get_End_Of( t1, 'year'  )    ! 00:00 at 01 Jan of next year
!     t2 = Get_End_Of( t1, 'month' )    ! 00:00 at first day of next month
!     t2 = Get_End_Of( t1, 'day'   )    ! 00:00 at next day
!
!
! INQUIRY FUNCTIONS
!
!   A few inquiry functions are provided.
!
!   The logical functions 'LeapYear' and 'LeapDay' tell you
!   if the year has a Februari 29 and if is Februari 29:
!
!     l = LeapYear( t )
!     l = LeapDay( t )
!
!   Two integer functions are provided to count the total number
!   of days in a month or a year:
!
!     i = Days_in_Month( t )
!     i = Days_in_Year( t )
!
!   An integer function is provided to return the day number,
!   counting from 1 (Januari 1) to 360, 365, or 366 (last of December):
!
!     i = DayNumber( t )
!
!   The logical function 'Midnight' is true for times 00:00 :
!
!     l = Midnight( t )
!
!   The logical function 'Precisely' is true for times
!   that are precisely 1.5 hour/minute/etc:
!
!     l = Precisely( t,  1.5, 'hour' )    ! every 1.5 hour
!     l = Precisely( t, 24.0, 'hour' )    ! midnight
!
!
! OPERATORS
!
!   Operators '+' and '-' are redefined to perform operations
!   between two date structures.
!   Both should be of same calendar type, or one should be
!   an increment:
!
!     t = t1 + t2
!     t = t1 - t2
!
!   Operators '*' and '/' are redefined for multiplication with
!   or division by a real or an integer:
!
!     t = t1 + dt * 2
!     t = t1 + dt * 3.1415
!     t = t1 + dt / 3.1415
!
!
! LOGICAL OPERATORS
!
!   Operators '==', '/=', '<', '<=', '>', '>=' are defined to
!   compare two dates.
!
!
! SUMMATION ROUTINES
!
!   The total number in a certain unit is returned by 'rTotal'
!   (real value) or 'iTotal' (integer, error if not possible).
!   Currently supported units are 'year', 'month', 'day',
!   'hour', 'minute', 'sec', and 'mili'. If the total number is
!   not wel defined for a certain date, for example the
!   total number of years of today, an error message is produced.
!
!     r = rTotal( t, 'year'|'month'|... )
!     i = iTotal( t, 'year'|'month'|... )
!
!
! INTERPOLATION
!
!   For t in [t1,t2], return real coefficients alfa1 and alf2 such that:
!       t  =  alfa1 * t1  +  alfa2 * t2
!   Usefull for linear interpolation:
!       f(t)  ~  alfa1 * f(t1)  +  alfa2 * f(t2)
!
!     call InterpolFractions( t, t1, t2, alfa1, alfa2, status )
!
!   Return interval [tt(1),tt(2)] around t that matches with time resolution;
!   resolution specified by a step and unit:
!      3, 'hour'   # 00:00, 03:00, 06:00, ...
!
!     call Get_Surrounding_Interval( t, 3, 'hour', tt, status )
!
!
! UNIT CONVERSION
!
!  Expand time value from number to T_CSO_DateTime using the units description;
!  also 'months' is allowed as step:
!
!    call ExpandTime( 2.5, 'hours since 1900-01-01 00:00:0.0', t, status )
!
!  Compare 't' of type 'T_CSO_DateTime' with a numeric 'value' given its units
!  ("days since ...") and the calendar ("366_day") ;
!  return status: -1 if time and value do not match, 0 if match :
!
!    call Compare_Date_Num( t, value, 'days since 0000-01-00 00:00', '366_day', status )
!    if (status<0) stop 'time does not match with value'
!
!
! OUTPUT
!
!   To obtain a pretty formatted print of the value of a date,
!   the 'Pretty' routine is provided. Output differs based on
!   the calendar type.
!
!     print *, 't  = '//trim(Pretty(t))
!     print *, 'dt = '//trim(Pretty(dt))
!
!   Some compilers have problems with this kind of statements.
!   Therefore, also a routine is provided:
!
!     call PrintDate( 't = ', t )
!
!   Also printing to the 'gol' buffer from GO_Print is supported:
!
!     call wrtcsol( 'time       : ', t ); call csoPr
!     call wrtcsol( 'time step  : ', dt ); call csoPr
!     call wrtcsol( 'time range : ', t1, ' to ', t2 ); call csoPr
!
!
! DEFAULTS
!
!   For setting some default values, the subroutine 'CSO_DateTimeDefaults'
!   is available. All arguments are optional:
!
!     call CSO_DateTimeDefaults( [calendar='gregorian'] )
!
!
!### macro's ##################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call csoErr
!
#define IF_NOT_OK_STOP if (status/=0) then; TRACEBACK; stop; end if
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!###############################################################################


module CSO_DateTimes

  use CSO_Logging, only : csol, csoErr, csoPr

  implicit none

  ! --- in/out ---------------------------

  private

  public      ::  T_CSO_DateTime, T_CSO_TimeDelta

  public      ::  goDateDefaults

  public      ::  CSO_DateTime, CSO_TimeDelta, AnyDate, SystemDate
  public      ::  Get_Begin_Of, Get_End_Of

  public      ::  Set, Get
  public      ::  Check
  public      ::  Normalize

  public      ::  LeapYear
  public      ::  LeapDay
  public      ::  days_in_month, calc_days_in_month
  public      ::  days_in_year
  public      ::  DayNumber, calc_DayNumber
  public      ::  Midnight
  public      ::  Precisely

  public      ::  operator(+)
  public      ::  operator(-)
  public      ::  operator(*)
  public      ::  operator(/)

  public      ::  IsAnyDate
  public      ::  operator(==)
  public      ::  operator(/=)
  public      ::  operator(>)
  public      ::  operator(<)
  public      ::  operator(>=)
  public      ::  operator(<=)

  public      ::  min, max

  public      ::  rTotal, iTotal, dTotal

  public      ::  InterpolFractions

  public      ::  Get_Surrounding_Interval

  public      ::  CSO_ReadFromLine
  public      ::  ExpandTime
  public      ::  Compare_Date_Num

  public      ::  Pretty
  public      ::  wrtcsol
  public      ::  PrintDate

  ! --- const -----------------------------------

  character(len=*), parameter  ::  mname = 'CSO_DateTimes'

  ! --- types -------------------------------------

  ! Strucure with fields to store year, month, day,
  ! hour and minute.
  ! Operators for assignment (=), adding (+),
  ! and comparission (==,<,>,>= and <=)
  ! have been defined for operations between
  ! instances of this type.

  type T_CSO_DateTime
    ! type of calendar: 'gregorian', 'standard', '365_day', '360_day'
    character(len=32)   ::  calendar
    ! year, month etc:
    integer             ::  year, month, day, hour, minute, second, milisecond
    ! difference with Coordinated Universal Time (UTC)
    integer             ::  zone   ! minutes
    ! error status
    integer             ::  status ! = 1
  end type T_CSO_DateTime


  type T_CSO_TimeDelta
    ! days, hours, etc:
    integer             ::  day, hour, minute, second, milisecond
    ! error status
    integer             ::  status ! = 1
  end type T_CSO_TimeDelta

  !! testing ...
  !type DateTime_Interval
  !  ! boundaries:
  !  type(T_CSO_DateTime)         ::  left, right
  !  ! open ? default is closed:
  !  logical             ::  left_open, right_open
  !end type DateTime_Interval


  ! --- var --------------------------------

  ! default calendar type
  character(len=32)      ::  default_calendar = 'gregorian'


  ! --- interface ---------------------------

  interface Check
    module procedure date_Check
    module procedure incrdate_Check
  end interface

  ! *

  interface LeapYear
    module procedure calc_LeapYear
    module procedure date_LeapYear
  end interface

  interface LeapDay
    module procedure date_LeapDay
  end interface

  interface days_in_month
    module procedure date_days_in_month
  end interface

  interface days_in_year
    module procedure date_days_in_year
    module procedure calc_days_in_year_greg
  end interface

  interface DayNumber
    module procedure date_DayNumber
  end interface

  interface Midnight
    module procedure date_Midnight
  end interface

  interface Precisely
    module procedure date_Precisely
  end interface

  ! *

  interface Set
    module procedure date_Set
    module procedure incrdate_Set
  end interface

  interface Get
    module procedure date_Get
    module procedure incrdate_Get
  end interface

  ! * operators

  interface Normalize
    module procedure date_Normalize
    module procedure incrdate_Normalize
  end interface

  interface operator(+)
    module procedure t_plus_t
    module procedure t_plus_dt
    module procedure dt_plus_dt
  end interface

  interface operator(-)
    module procedure t_min_t
    module procedure t_min_dt
    module procedure dt_min_dt
  end interface

  interface operator(*)
    module procedure dt_times_r4
    module procedure r4_times_dt
    module procedure dt_times_r8
    module procedure r8_times_dt
    module procedure dt_times_i
    module procedure  i_times_dt
  end interface

  interface operator(/)
    module procedure dt_div_r
    module procedure dt_div_i
  end interface

  ! * logical operators

  interface IsAnyDate
    module procedure date_IsAnyDate
  end interface

  interface operator(==)
    module procedure date_eq_date
  end interface

  interface operator(/=)
    module procedure date_ne_date
  end interface

  interface operator(>)
    module procedure date_gt_date
  end interface

  interface operator(<)
    module procedure date_lt_date
  end interface

  interface operator(>=)
    module procedure date_ge_date
  end interface

  interface operator(<=)
    module procedure date_le_date
  end interface

  ! *

  interface min
    module procedure date_min
  end interface

  ! *

  interface max
    module procedure date_max
  end interface

  ! *
  
  interface iTotal
    module procedure date_iTotal
    module procedure incr_iTotal
  end interface
     
  interface rTotal
    module procedure date_rTotal
    module procedure incr_rTotal
  end interface

  interface dTotal
    module procedure date_dTotal
    module procedure incr_dTotal
  end interface

  ! *

  interface InterpolFractions
    module procedure date_InterpolFractions
  end interface

  ! *
  
  interface CSO_ReadFromLine
    module procedure CSO_ReadFromLine_t
  end interface

  interface Pretty
    module procedure date_Pretty
    module procedure incrdate_Pretty
  end interface

  interface wrtcsol
    module procedure wrtcsol_t
    module procedure wrtcsol_dt
    module procedure wrtcsol_t_dt
    module procedure wrtcsol_tt
    module procedure wrtcsol_t1_t2
    module procedure wrtcsol_t1_t2_t3
  end interface



contains


  ! ****************************************************
  ! ***
  ! *** set defaults
  ! ***
  ! ****************************************************


  subroutine goDateDefaults( calendar )

    ! --- in/out --------------------------------

    character(len=*), intent(in), optional    ::  calendar

    ! --- begin ----------------------------------

    if ( present(calendar) ) default_calendar = calendar

  end subroutine goDateDefaults


  ! ****************************************************
  ! ***
  ! *** check
  ! ***
  ! ****************************************************

  !
  ! Check fields of a date:
  !   range etc
  !

  subroutine date_Check( t, status )

    ! --- in/out ----------------------------------

    type(T_CSO_DateTime), intent(in)        ::  t
    integer, intent(out)           ::  status

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Check'

    ! --- begin -----------------------------------

    ! already error status ? then leave immediatelly:
    if ( t%status /= 0 ) then
      write (csol,'("found error status in date")'); call csoErr
      write (csol,'("  year,month,day : ",3i6)') t%year, t%month, t%day; call csoErr
      write (csol,'("  hour,minute,sec,mili : ",4i6)') t%hour, t%minute, t%second, t%milisecond; call csoErr
      TRACEBACK; status=1; return
    end if

    ! calendar specific
    select case ( t%calendar )
      case ( 'none' )
        ! always ok ...
        status = 0
        return
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day', '360_day' )
        ! check month
        if ( t%month<1 .or. t%month>12 ) then
          call wrtcsol( 'strange month in ', t ); call csoErr
          TRACEBACK; status=1; return
        end if
        ! check day
        if ( t%day<1 .or. t%day>days_in_month(t) ) then
          call wrtcsol( 'strange day in ', t ); call csoErr
          TRACEBACK; status=1; return
        end if
        ! zone should be zero:
        if ( t%zone /= 0 ) then
          call wrtcsol( 'expecting zero zone in date ', t ); call csoErr
          TRACEBACK; status=1; return
        end if
      case default
        write (csol,'("unknown calendar type: `",a,"`")') t%calendar; call csoErr
        write (csol,'("  year etc : ",6i5)') t%year, t%month, t%day, t%hour, t%minute, t%second; call csoErr
        TRACEBACK; status=1; return
    end select

    ! check minutes
    if ( t%minute<0 .or. t%minute>59 ) then
      call wrtcsol( 'found strange minutes in ', t ); call csoErr
      TRACEBACK; status=1; return
    end if

    ! check seconds
    if ( t%second<0 .or. t%second>59 ) then
      call wrtcsol( 'found strange seconds in ', t ); call csoErr
      TRACEBACK; status=1; return
    end if

    ! check mili
    if ( t%milisecond<0 .or. t%milisecond>999 ) then
      call wrtcsol( 'found strange mili seconds in ', t ); call csoErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine date_Check


  ! ***


  subroutine incrdate_Check( dt, status )

    ! --- in/out ----------------------------------

    type(T_CSO_TimeDelta), intent(in)    ::  dt
    integer, intent(out)           ::  status

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incrdate_Check'

    ! --- begin -----------------------------------

    ! already error status ? then leave immediatelly:
    if ( dt%status /= 0 ) then
      write (csol,'("found error status in incrdate")'); call csoErr
      write (csol,'("  day, hour,minu,sec,mili : ",5i6)') dt%day, dt%hour, dt%minute, dt%second, dt%milisecond; call csoErr
      TRACEBACK; status=1; return
    end if

    ! every value is allowed for increments ...

    ! ok
    status = 0

  end subroutine incrdate_Check



  ! ****************************************************
  ! ***
  ! *** computation
  ! ***
  ! ****************************************************

  ! Does this year have a 29 feb ?

  logical function calc_LeapYear( year )

    ! --- in/out -------------------------------

    integer, intent(in)           ::  year

    ! --- begin --------------------------------

    calc_LeapYear = ( (mod(year,4)==0) .and. .not.(mod(year,100)==0) ) &
                                        .or. (mod(year,400)==0)

  end function calc_LeapYear


  ! ***


  ! days per month

  integer function calc_days_in_month( calendar, year, month )

    ! --- in/out ---------------------------

    character(len=*), intent(in)   ::  calendar
    integer, intent(in)            ::  year, month

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_days_in_month'

    ! --- const -----------------------------

    ! days in a month                      1  2  3  4  5  6  7  8  9 10 11 12
    integer, parameter :: days365(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! normal
    integer, parameter :: days366(12) = (/31,29,31,30,31,30,31,31,30,31,30,31/) ! leap year
    integer, parameter :: days360(12) = (/30,30,30,30,30,30,30,30,30,30,30,30/) ! fixed month

    ! --- begin ----------------

    select case ( calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard' )
        if ( calc_LeapYear(year) ) then
          calc_days_in_month = days366(month)
        else
          calc_days_in_month = days365(month)
        end if
      case ( '366_day' )
        calc_days_in_month = days366(month)
      case ( '365_day' )
        calc_days_in_month = days365(month)
      case ( '360_day' )
        calc_days_in_month = days360(month)
      case ( 'none' )
        calc_days_in_month = 0
      case default
        calc_days_in_month = -1
        write (csol,'("unknown calendar type: ",a)') calendar; call csoErr
        TRACEBACK; stop
    end select

  end function calc_days_in_month


  ! ***


  ! days per year

  integer function calc_days_in_year( calendar, year )

    ! --- in/out ----------------

    character(len=*), intent(in)   ::  calendar
    integer, intent(in)            ::  year

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_days_in_year'

    ! --- begin ----------------

    select case ( calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard' )
        if ( calc_LeapYear(year) ) then
          calc_days_in_year = 366
        else
          calc_days_in_year = 365
        end if
      case ( '366_day' )
        calc_days_in_year = 366
      case ( '365_day' )
        calc_days_in_year = 365
      case ( '360_day' )
        calc_days_in_year = 360
      case ( 'none' )
        calc_days_in_year = 0
      case default
        write (csol,'("unknown calendar type: ",a)') calendar; call csoErr
        TRACEBACK; stop
    end select

  end function calc_days_in_year


  ! by default the gregorian calendar:

  integer function calc_days_in_year_greg( year )

    ! --- in/out ----------------

    integer, intent(in)            ::  year

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_days_in_year_greg'

    ! --- begin ----------------

    ! extract for Gregorian calendar:
    calc_days_in_year_greg = calc_days_in_year( 'gregorian', year )

  end function calc_days_in_year_greg


  ! ***


  ! Returns the number of the day spedified by the date iy/im/id.
  ! The existence of a february 29 is checked.
  !
  ! ndays( 1995,  1,  1 ) =   1
  ! ndays( 1995, 12, 31 ) = 365
  ! ndays( 1996, 12, 31 ) = 366        29 feb every 4 year ...
  ! ndays( 1900, 12, 31 ) = 365         except every 100 year ...
  ! ndays( 2000, 12, 31 ) = 366          except every 400 year ...

  integer function calc_DayNumber( calendar, year, month, day )

    ! --- in/out ----------------------------

    character(len=*), intent(in)   ::  calendar
    integer, intent(in)            ::  year, month, day

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/calc_DayNumber'

    ! --- local -----------------------------

    integer            ::  imonth

    ! --- begin ----------------------------

    select case ( calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day', '360_day' )
        calc_DayNumber = day
        do imonth = 1, month-1
          calc_DayNumber = calc_DayNumber + calc_days_in_month(calendar,year,imonth)
        end do
      case ( 'none' )
        calc_DayNumber = 0
      case default
        write (csol,'("unknown calendar type: ",a)') calendar; call csoErr
        TRACEBACK; stop
    end select

  end function calc_DayNumber


  ! **********************************************


  logical function date_LeapYear( t )

    ! --- in/out -------------------------------

    type(T_CSO_DateTime), intent(in)         ::  t

    ! --- begin --------------------------------

    ! calendar specific
    select case ( t%calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard' )
        ! see above ...
        date_LeapYear = calc_LeapYear( t%year )
      case default
        ! no leap years ...
        date_LeapYear = .false.
    end select

  end function date_LeapYear


  ! ***


  logical function date_LeapDay( t )

    ! --- in/out -------------------------------

    type(T_CSO_DateTime), intent(in)         ::  t

    ! --- begin --------------------------------

    ! check day ...
    date_LeapDay = (t%month == 2) .and. (t%day == 29)

  end function date_LeapDay


  ! ***


  integer function date_days_in_month( t )

    ! --- in/out ----------------

    type(T_CSO_DateTime), intent(in)         ::  t

    ! --- begin ----------------

    date_days_in_month = calc_days_in_month( t%calendar, t%year, t%month )

  end function date_days_in_month


  ! ***


  integer function date_days_in_year( t )

    ! --- in/out ----------------

    type(T_CSO_DateTime), intent(in)         ::  t

    ! --- begin ----------------

    date_days_in_year = calc_days_in_year( t%calendar, t%year )

  end function date_days_in_year


  ! ***


  integer function date_DayNumber( t )

    ! --- in/out ----------------

    type(T_CSO_DateTime), intent(in)        ::  t

    ! --- begin ----------------

    date_DayNumber = calc_DayNumber( t%calendar, t%year, t%month, t%day )

  end function date_DayNumber


  ! ***


  logical function date_Midnight( t )

    ! --- in/out -------------------------------

    type(T_CSO_DateTime), intent(in)         ::  t

    ! --- begin --------------------------------

    date_Midnight = all( (/t%hour,t%minute,t%second,t%milisecond/) == 0 )

  end function date_Midnight


  ! ***


  logical function date_Precisely( t, val, unit )

    ! --- in/out -------------------------------

    type(T_CSO_DateTime), intent(in)         ::  t
    real, intent(in)                ::  val
    character(len=*), intent(in)    ::  unit

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Precisely'

    ! --- local --------------------------------

    real      ::  rtot

    ! --- begin --------------------------------

    ! only some units supported yet:
    select case ( unit )

      ! every .. minute
      case ( 'min', 'minu' )
        ! check ...
        if ( (val <= 0.0) .or. (val > 60.0) ) then
          write (csol,'("strange value for minu : ",f8.2)') val; call csoErr
          TRACEBACK; stop
        end if
        ! total minute fraction in this day:
        rtot = rTotal( CSO_TimeDelta(0,t%hour,t%minute,t%second,t%milisecond), 'min' )
        ! modulo given value ?
        date_Precisely = modulo(rtot,val) == 0.0

      ! every .. hour
      case ( 'hour' )
        ! check ...
        if ( (val <= 0.0) .or. (val > 24.0) ) then
          write (csol,'("strange value for hour : ",f8.2)') val; call csoErr
          TRACEBACK; stop
        end if
        ! total hour fraction in this day:
        rtot = rTotal( CSO_TimeDelta(0,t%hour,t%minute,t%second,t%milisecond), 'hour' )
        ! modulo given value ?
        date_Precisely = modulo(rtot,val) == 0.0

      ! error ...
      case default
        write (csol,'("unsupported unit : ",a)') unit; call csoErr
        TRACEBACK; stop

    end select

  end function date_Precisely



  ! ****************************************************
  ! ***
  ! *** Set/Get fields in a structure T_CSO_DateTime
  ! ***
  ! ****************************************************

  !
  ! Fill fields in a 'T_CSO_DateTime' structure.
  !

  subroutine date_Set( date, year, month, day, hour, min, sec, mili, &
                             zone, calendar, time4, time5, time6 )

    ! --- in/out ------------------------------------

    type(T_CSO_DateTime), intent(inout)              ::  date
    integer, intent(in), optional            ::  year
    integer, intent(in), optional            ::  month
    integer, intent(in), optional            ::  day
    integer, intent(in), optional            ::  hour
    integer, intent(in), optional            ::  min
    integer, intent(in), optional            ::  sec
    integer, intent(in), optional            ::  mili
    integer, intent(in), optional            ::  zone
    character(len=*), intent(in), optional   ::  calendar
    integer, intent(in), optional            ::  time4(4)
    integer, intent(in), optional            ::  time5(5)
    integer, intent(in), optional            ::  time6(6)

    ! --- local ----------------------------------

    if ( present(calendar) ) date%calendar  = calendar

    if ( present(time4) ) then
      date%year  = time4(1)
      date%month = time4(2)
      date%day   = time4(3)
      date%hour  = time4(4)
    end if

    if ( present(time5) ) then
      date%year  = time5(1)
      date%month = time5(2)
      date%day   = time5(3)
      date%hour  = time5(4)
      date%minute   = time5(5)
    end if

    if ( present(time6) ) then
      date%year  = time6(1)
      date%month = time6(2)
      date%day   = time6(3)
      date%hour  = time6(4)
      date%minute   = time6(5)
      date%second   = time6(6)
    end if

    if ( present(year ) ) date%year  = year
    if ( present(month) ) date%month = month
    if ( present(day  ) ) date%day   = day
    if ( present(zone ) ) date%zone  = zone
    if ( present(hour ) ) date%hour  = hour
    if ( present(min  ) ) date%minute   = min
    if ( present(sec  ) ) date%second   = sec
    if ( present(mili ) ) date%milisecond  = mili

  end subroutine date_Set


  ! *

  subroutine incrdate_Set( date, day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    type(T_CSO_TimeDelta), intent(inout)             ::  date
    integer, intent(in), optional            ::  day
    integer, intent(in), optional            ::  hour
    integer, intent(in), optional            ::  min
    integer, intent(in), optional            ::  sec
    integer, intent(in), optional            ::  mili

    ! --- local ----------------------------------

    if ( present(day  ) ) date%day   = day
    if ( present(hour ) ) date%hour  = hour
    if ( present(min  ) ) date%minute   = min
    if ( present(sec  ) ) date%second   = sec
    if ( present(mili ) ) date%milisecond  = mili

  end subroutine incrdate_Set


  !
  ! Obtain fields from a 'T_CSO_DateTime' structure.
  !

  subroutine date_Get( date, &
                       year, month, day, hour, min, sec, mili, &
                       zone, calendar, time4, time5, time6 )

    ! --- in/out ------------------------------------

    type(T_CSO_DateTime), intent(in)                  ::  date
    integer, intent(out), optional           ::  year, month, day
    integer, intent(out), optional           ::  hour, min, sec, mili
    integer, intent(out), optional           ::  zone
    integer, intent(out), optional           ::  time4(4)
    integer, intent(out), optional           ::  time5(5)
    integer, intent(out), optional           ::  time6(6)
    character(len=*), intent(out), optional  ::  calendar

    ! --- local ----------------------------------

    if ( present(calendar) ) calendar = date%calendar

    if ( present(year)  ) year  = date%year
    if ( present(month) ) month = date%month
    if ( present(day  ) ) day   = date%day
    if ( present(zone ) ) zone  = date%zone
    if ( present(hour ) ) hour  = date%hour
    if ( present(min  ) ) min   = date%minute
    if ( present(sec  ) ) sec   = date%second
    if ( present(mili ) ) mili  = date%milisecond

    if ( present(time4) ) time4 = (/ date%year, date%month, date%day, date%hour /)

    if ( present(time5) ) time5 = (/ date%year, date%month, date%day, &
                                     date%hour, date%minute   /)

    if ( present(time6) ) time6 = (/ date%year, date%month, date%day, &
                                     date%hour, date%minute  , date%second /)

  end subroutine date_Get

  ! *

  subroutine incrdate_Get( date, day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    type(T_CSO_TimeDelta), intent(in)              ::  date
    integer, intent(out), optional           ::  day
    integer, intent(out), optional           ::  hour, min, sec, mili

    ! --- local ----------------------------------

    if ( present(day  ) ) day   = date%day
    if ( present(hour ) ) hour  = date%hour
    if ( present(min  ) ) min   = date%minute
    if ( present(sec  ) ) sec   = date%second
    if ( present(mili ) ) mili  = date%milisecond

  end subroutine incrdate_Get


  ! ****************************************************
  ! ***
  ! *** Return a date structure
  ! ***
  ! ****************************************************

  !
  ! Set fields to zero or fill some of them.
  !

  type(T_CSO_DateTime) function CSO_DateTime( year, month, day, hour, min, sec, mili, &
                                     zone, calendar, time4, time5, time6 )

    ! --- in/out ------------------------------------

    integer, intent(in), optional            ::  year, month, day
    integer, intent(in), optional            ::  hour, min, sec, mili
    integer, intent(in), optional            ::  zone
    character(len=*), intent(in), optional   ::  calendar
    integer, intent(in), optional            ::  time4(4)
    integer, intent(in), optional            ::  time5(5)
    integer, intent(in), optional            ::  time6(6)

    ! --- local ----------------------------------

    ! set default calendar type:
    CSO_DateTime%calendar = default_calendar

    ! Fields are zero by default:
    CSO_DateTime%year  = 0
    CSO_DateTime%month = 0
    CSO_DateTime%day   = 0
    CSO_DateTime%zone  = 0
    CSO_DateTime%hour  = 0
    CSO_DateTime%minute   = 0
    CSO_DateTime%second   = 0
    CSO_DateTime%milisecond  = 0

    ! Optionally, change some of them:
    if ( present(year    ) ) call Set( CSO_DateTime, year=year )
    if ( present(month   ) ) call Set( CSO_DateTime, month=month )
    if ( present(day     ) ) call Set( CSO_DateTime, day=day )
    if ( present(hour    ) ) call Set( CSO_DateTime, hour=hour )
    if ( present(min     ) ) call Set( CSO_DateTime, min=min )
    if ( present(sec     ) ) call Set( CSO_DateTime, sec=sec )
    if ( present(mili    ) ) call Set( CSO_DateTime, mili=mili )
    if ( present(zone    ) ) call Set( CSO_DateTime, zone=zone )
    if ( present(calendar) ) call Set( CSO_DateTime, calendar=calendar )
    if ( present(time4   ) ) call Set( CSO_DateTime, time4=time4 )
    if ( present(time5   ) ) call Set( CSO_DateTime, time5=time5 )
    if ( present(time6   ) ) call Set( CSO_DateTime, time6=time6 )

    ! normalize too small/too large values:
    !if ( CSO_DateTime%year /= 0000 ) call Normalize( CSO_DateTime )
    if ( CSO_DateTime%calendar /= 'none' ) call Normalize( CSO_DateTime )

    ! data filled, thus probably no error ...
    CSO_DateTime%status = 0

  end function CSO_DateTime


  ! ***


  type(T_CSO_TimeDelta) function CSO_TimeDelta( day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    integer, intent(in), optional            ::  day
    integer, intent(in), optional            ::  hour, min, sec, mili

    ! --- local ----------------------------------

    ! Fields are zero by default:
    CSO_TimeDelta%day   = 0
    CSO_TimeDelta%hour  = 0
    CSO_TimeDelta%minute   = 0
    CSO_TimeDelta%second   = 0
    CSO_TimeDelta%milisecond  = 0

    ! Optionally, change some of them:
    if ( present(day     ) ) call Set( CSO_TimeDelta, day=day )
    if ( present(hour    ) ) call Set( CSO_TimeDelta, hour=hour )
    if ( present(min     ) ) call Set( CSO_TimeDelta, min=min )
    if ( present(sec     ) ) call Set( CSO_TimeDelta, sec=sec )
    if ( present(mili    ) ) call Set( CSO_TimeDelta, mili=mili )

    ! normalize too small/too large values:
    call Normalize( CSO_TimeDelta )

    ! data filled, thus probably no error ...
    CSO_TimeDelta%status = 0

  end function CSO_TimeDelta


  ! ***


  !
  ! Set fields to zero, special calendar
  !

  type(T_CSO_DateTime) function AnyDate()

    ! --- local ----------------------------------

    ! Set some fields, other are automatically zero:
    AnyDate = CSO_DateTime( calendar='none' )

  end function AnyDate


  ! ***


  ! Fill with system time

  type(T_CSO_DateTime) function SystemDate()

    ! --- in/out ------------------------------

    ! none ...

    ! --- local ------------------------------

    integer           ::  values(8)

    ! --- begin ------------------------------

    !
    ! Optional character output of Date_and_Time:
    !
    !   date   '20020812'
    !   time   '211757.314'
    !   zone   '+0200'
    !

    ! obtain system date and time:
    call Date_and_Time( values=values )

    ! fill fields in structure:
    call Set( SystemDate, calendar='gregorian', &
                year=values(1), month=values(2), day=values(3), &
                zone=values(4), hour=values(5), &
                min=values(6), sec=values(7), mili=values(8) )

    ! no error ...
    SystemDate%status = 0

  end function SystemDate


  ! ***


  ! set to begin of : 'year', 'month', 'day'

  type(T_CSO_DateTime) function Get_Begin_Of( t, what )

    ! --- in/out ------------------------------

    type(T_CSO_DateTime), intent(in)           ::  t
    character(len=*), intent(in)      ::  what

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/Get_Begin_Of'

    ! --- local ------------------------------

    ! --- begin ------------------------------

    ! end of what ?
    select case ( what )

      case ( 'year' )

        ! begin of year:
        Get_Begin_Of = CSO_DateTime( t%year, 01, 01, 00, 00, 00 )

      case ( 'month' )

        ! begin of month:
        Get_Begin_Of = CSO_DateTime( t%year, t%month, 01, 00, 00, 00 )

      case ( 'day' )

        ! begin of day:
        Get_Begin_Of = CSO_DateTime( t%year, t%month, t%day, 00, 00, 00 )

      case default

        write (csol,'("end of `",a,"` not supported yet")') what; call csoErr
        TRACEBACK; stop

    end select

    ! no error ...
    Get_Begin_Of%status = 0

  end function Get_Begin_Of


  ! ***


  ! set to end of : 'month', ...

  type(T_CSO_DateTime) function Get_End_Of( t, what )

    ! --- in/out ------------------------------

    type(T_CSO_DateTime), intent(in)           ::  t
    character(len=*), intent(in)      ::  what

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/Get_End_Of'

    ! --- local ------------------------------

    ! --- begin ------------------------------

    ! first set to begin:
    Get_End_Of = Get_Begin_Of( t, what )

    ! end of what ?
    select case ( what )
      case ( 'year' )
        ! add number of days:
        Get_End_Of = Get_End_Of + CSO_TimeDelta(day=Days_in_Year(t))
      case ( 'month' )
        ! add number of days:
        Get_End_Of = Get_End_Of + CSO_TimeDelta(day=Days_in_Month(t))
      case ( 'day' )
        ! add single day:
        Get_End_Of = Get_End_Of + CSO_TimeDelta(day=1)
      case default
        write (csol,'("end of `",a,"` not supported yet")') what; call csoErr
        TRACEBACK; stop
    end select

    ! no error ...
    Get_End_Of%status = 0

  end function Get_End_Of


  ! ************************************************
  ! ***
  ! *** operators
  ! ***
  ! ************************************************


  subroutine date_Normalize( t )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(inout)       ::  t

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Normalize'

    ! --- begin ---------------------------------

    ! mili seconds
    do
      if ( t%milisecond >= 0 ) exit
      t%second = t%second - 1
      t%milisecond = t%milisecond + 1000
    end do
    do
      if ( t%milisecond <= 999 ) exit
      t%milisecond = t%milisecond - 1000
      t%second = t%second + 1
    end do

    ! seconds
    do
      if ( t%second >= 0 ) exit
      t%minute = t%minute - 1
      t%second = t%second + 60
    end do
    do
      if ( t%second <= 59 ) exit
      t%second = t%second - 60
      t%minute = t%minute + 1
    end do

    ! minutes
    do
      if ( t%minute >= 0 ) exit
      t%hour = t%hour - 1
      t%minute = t%minute + 60
    end do
    do
      if ( t%minute <= 59 ) exit
      t%minute = t%minute - 60
      t%hour = t%hour + 1
    end do

    ! hours
    do
      if ( t%hour >= 0 ) exit
      t%day = t%day - 1
      t%hour = t%hour + 24
    end do
    do
      if ( t%hour <= 23 ) exit
      t%hour = t%hour - 24
      t%day = t%day + 1
    end do

    ! days, months, year
    select case ( t%calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day', '360_day' )
        do
          if ( t%day >= 1 ) exit
          t%month = t%month - 1
          do
            if ( t%month >= 1 ) exit
            t%year = t%year - 1
            t%month = t%month + 12
          end do
          t%day = t%day + days_in_month(t)
        end do
        do
          if ( t%day <= days_in_month(t) ) exit
          t%day = t%day - days_in_month(t)
          t%month = t%month + 1
          do
            if ( t%month <= 12 ) exit
            t%month = t%month - 12
            t%year = t%year + 1
          end do
        end do
      case default
        write (csol,'("unsupported calendar type: ",a)') t%calendar; call csoErr
        TRACEBACK; stop
    end select

  end subroutine date_Normalize


  subroutine incrdate_Normalize( dt )

    ! --- in/out --------------------------------

    type(T_CSO_TimeDelta), intent(inout)    :: dt

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/incrdate_Normalize'

    ! --- begin ---------------------------------

    ! mili seconds
    do
      if ( dt%milisecond >= 0 ) exit
      dt%second = dt%second - 1
      dt%milisecond = dt%milisecond + 1000
    end do
    do
      if ( dt%milisecond <= 999 ) exit
      dt%milisecond = dt%milisecond - 1000
      dt%second = dt%second + 1
    end do

    ! seconds
    do
      if ( dt%second >= 0 ) exit
      dt%minute = dt%minute - 1
      dt%second = dt%second + 60
    end do
    do
      if ( dt%second <= 59 ) exit
      dt%second = dt%second - 60
      dt%minute = dt%minute + 1
    end do

    ! minutes
    do
      if ( dt%minute >= 0 ) exit
      dt%hour = dt%hour - 1
      dt%minute = dt%minute + 60
    end do
    do
      if ( dt%minute <= 59 ) exit
      dt%minute = dt%minute - 60
      dt%hour = dt%hour + 1
    end do

    ! hours
    do
      if ( dt%hour >= 0 ) exit
      dt%day = dt%day - 1
      dt%hour = dt%hour + 24
    end do
    do
      if ( dt%hour <= 23 ) exit
      dt%hour = dt%hour - 24
      dt%day = dt%day + 1
    end do

  end subroutine incrdate_Normalize


  ! *** date = t1 + t2 ************************


  !
  !        t1             +   t2             ->   t1+t2
  !
  !        greg               incr                greg
  !        366                incr                366
  !        365                incr                365
  !
  !        360                360                 360
  !        360                incr                360
  !
  !        incr               greg                greg
  !        incr               366                 366
  !        incr               365                 365
  !        incr               360                 360
  !        incr               incr                incr
  !

  type(T_CSO_DateTime) function t_plus_t( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_plus_t'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t1, status )
    IF_NOT_OK_STOP
    call Check( t2, status )
    IF_NOT_OK_STOP

    ! any date ? return any date ..
    if ( (t1%calendar == 'none') .or. (t2%calendar == 'none') ) then
      t_plus_t = AnyDate()
      return
    end if

    ! calendars should be the same:
    if ( t1%calendar /= t2%calendar ) then
      write (csol,'("calendars should be the same : ")'); call csoPr
      write (csol,'("  t1 : ",a)') trim(t1%calendar); call csoPr
      write (csol,'("  t2 : ",a)') trim(t2%calendar); call csoPr
      TRACEBACK; stop
    end if

    ! add all fields;
    t_plus_t = CSO_DateTime( calendar=t1%calendar, &
                        year  = t1%year  + t2%year  , &
                        month = t1%month + t2%month , &
                        day   = t1%day   + t2%day   , &
                        hour  = t1%hour  + t2%hour  , &
                        zone  = t1%zone  + t2%zone  , &
                        min   = t1%minute   + t2%minute   , &
                        sec   = t1%second   + t2%second   , &
                        mili  = t1%milisecond  + t2%milisecond     )

  end function t_plus_t


  ! *


  type(T_CSO_DateTime) function t_plus_dt( t, dt )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t
    type(T_CSO_TimeDelta), intent(in)   ::  dt

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_plus_dt'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t, status )
    IF_NOT_OK_STOP
    call Check( dt, status )
    IF_NOT_OK_STOP

    ! any date ? return any date ..
    if ( t%calendar == 'none' ) then
      t_plus_dt = AnyDate()
      return
    end if

    ! add fields; normalization is applied in routine:
    t_plus_dt = CSO_DateTime( calendar = t%calendar, &
                         year     = t%year             , &
                         month    = t%month            , &
                         day      = t%day   + dt%day   , &
                         hour     = t%hour  + dt%hour  , &
                         zone     = t%zone             , &
                         min      = t%minute   + dt%minute   , &
                         sec      = t%second   + dt%second   , &
                         mili     = t%milisecond  + dt%milisecond     )

  end function t_plus_dt


  ! *


  type(T_CSO_TimeDelta) function dt_plus_dt( dt1, dt2 )

    ! --- in/out --------------------------------

    type(T_CSO_TimeDelta), intent(in)   ::  dt1
    type(T_CSO_TimeDelta), intent(in)   ::  dt2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_plus_dt'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( dt1, status )
    IF_NOT_OK_STOP
    call Check( dt2, status )
    IF_NOT_OK_STOP

    ! add fields:
    dt_plus_dt = CSO_TimeDelta( day  = dt1%day  + dt2%day   , &
                           hour = dt1%hour + dt2%hour  , &
                           min  = dt1%minute  + dt2%minute   , &
                           sec  = dt1%second  + dt2%second   , &
                           mili = dt1%milisecond + dt2%milisecond     )

  end function dt_plus_dt


  ! *** date = t1 - t2

  !
  !        t1        ->  t2        ->   t1-t2      action
  !
  !        greg          greg           incr       difference
  !        greg          incr           greg       minus
  !
  !        366           366            incr       difference
  !        366           incr           366        minus
  !
  !        365           365            incr       difference
  !        365           incr           365        minus
  !
  !        360           360            360        difference
  !        360           incr           360        minus
  !
  !        incr          incr           incr       difference
  !


  type(T_CSO_TimeDelta) function t_min_t( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_min_t'

    ! --- local ---------------------------------

    ! --- begin ---------------------------------

    ! difference between two dates;
    ! algorithm implemented for positive differences only:

    ! difference should be positive:
    if ( t1 > t2 ) then

      ! call work routine:
      t_min_t = t_min_t_work( t1, t2 )

    else

      ! negative
      t_min_t = t_min_t_work( t2, t1 ) * (-1)

    end if

  end function t_min_t


  ! *

  type(T_CSO_TimeDelta) function t_min_t_work( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_min_t_work'

    ! --- local ---------------------------------

    integer             ::  status
    integer             ::  ndays
    type(T_CSO_DateTime)         ::  t

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t1, status )
    IF_NOT_OK_STOP
    call Check( t2, status )
    IF_NOT_OK_STOP

    ! any dates ? something wrong ...
    if ( (t1%calendar == 'none') .or. (t2%calendar == 'none') ) then
      write (csol,'("do not know how to compute difference between `any` dates ...")')
      TRACEBACK; stop
    end if

    ! calendars should be the same:
    if ( t1%calendar /= t2%calendar ) then
      write (csol,'("calendars should be the same : ")'); call csoPr
      write (csol,'("  t1 : ",a)') trim(t1%calendar); call csoPr
      write (csol,'("  t2 : ",a)') trim(t2%calendar); call csoPr
      TRACEBACK; stop
    end if

    ! difference between two dates;
    ! algorithm implemented for positive differences only:

    ! difference should be positive:
    if ( t1 < t2 ) then
      write (csol,'("expect t1 to exceed t2 :")'); call csoErr
      call wrtcsol( '  t1 : ', t1 ); call csoErr
      call wrtcsol( '  t2 : ', t2 ); call csoErr
      TRACEBACK; stop
    end if

    ! determine number of days between t1 and t2:
    t = t1
    ndays = daynumber(t) - 1
    do
      if ( t%year==t2%year ) exit
      t%year = t%year - 1
      ndays = ndays + days_in_year(t)
    end do
    ndays = ndays - (daynumber(t2)-1)

    ! store result:
    t_min_t_work = CSO_TimeDelta( day  = ndays, &
                             hour = t1%hour - t2%hour, &
                             min  = t1%minute  - t2%minute , &
                             sec  = t1%second  - t2%second , &
                             mili = t1%milisecond - t2%milisecond           )

  end function t_min_t_work


  ! *


  type(T_CSO_DateTime) function t_min_dt( t, dt )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t
    type(T_CSO_TimeDelta), intent(in)   ::  dt

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_min_dt'

    ! --- local ---------------------------------

    integer             ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t, status )
    IF_NOT_OK_STOP
    call Check( dt, status )
    IF_NOT_OK_STOP

    ! any date ? return any date ..
    if ( t%calendar == 'none' ) then
      t_min_dt = AnyDate()
      return
    end if

    ! result is of same type as t
    ! ~ first guess:
    t_min_dt = CSO_DateTime( calendar = t%calendar      , &
                        year     = t%year          , &
                        month    = t%month         , &
                        day      = t%day  -dt%day  , &
                        hour     = t%hour -dt%hour , &
                        zone     = t%zone          , &
                        min      = t%minute  -dt%minute  , &
                        sec      = t%second  -dt%second  , &
                        mili     = t%milisecond -dt%milisecond         )
    ! ~ normalize from negative or too large values into regular values:
    call Normalize( t_min_dt )

  end function t_min_dt


  ! *


  type(T_CSO_TimeDelta) function dt_min_dt( dt1, dt2 )

    ! --- in/out --------------------------------

    type(T_CSO_TimeDelta), intent(in)       ::  dt1
    type(T_CSO_TimeDelta), intent(in)       ::  dt2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_min_dt'

    ! --- local ---------------------------------

    integer             ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( dt1, status )
    IF_NOT_OK_STOP
    call Check( dt2, status )
    IF_NOT_OK_STOP

    ! fill result:
    dt_min_dt = CSO_TimeDelta( day  = dt1%day  - dt2%day , &
                          hour = dt1%hour - dt2%hour, &
                          min  = dt1%minute  - dt2%minute , &
                          sec  = dt1%second  - dt2%second , &
                          mili = dt1%milisecond - dt2%milisecond    )

  end function dt_min_dt


  ! *** date = t * r ************************************************


  ! multiply time with a real factor;
  ! use round for fractions

  type(T_CSO_TimeDelta) function dt_times_r4( dt, r )

    ! --- in/out --------------------------------

    type(T_CSO_TimeDelta), intent(in)   ::  dt
    real(4), intent(in)           ::  r

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_times_r4'

    ! --- local -----------------------------------

    real          ::  rday, rhour, rmin, rsec, rmili
    integer       ::  status

    ! --- begin ---------------------------------

    call Check( dt, status )
    IF_NOT_OK_STOP

    ! multiply each part, at remainder of previous:
    rday  = dt%day  * r
    rhour = dt%hour * r + (rday -nint(rday ))*24
    rmin  = dt%minute  * r + (rhour-nint(rhour))*60
    rsec  = dt%second  * r + (rmin -nint(rmin ))*60
    rmili = dt%milisecond * r + (rsec -nint(rsec ))*1000

    ! fill:
    dt_times_r4 = CSO_TimeDelta( day  = nint(rday ), &
                            hour = nint(rhour), &
                            min  = nint(rmin ), &
                            sec  = nint(rsec ), &
                            mili = nint(rmili)      )

  end function dt_times_r4


  ! *


  type(T_CSO_TimeDelta) function r4_times_dt( r, dt )

    ! --- in/out --------------------------------

    real(4), intent(in)           ::  r
    type(T_CSO_TimeDelta), intent(in)   ::  dt

    ! --- begin ---------------------------------

    r4_times_dt = dt * r

  end function r4_times_dt
  
  
  ! *
  

  type(T_CSO_TimeDelta) function dt_times_r8( dt, r )

    ! --- in/out --------------------------------

    type(T_CSO_TimeDelta), intent(in)   ::  dt
    real(8), intent(in)           ::  r

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_times_r8'

    ! --- local -----------------------------------

    real(8)       ::  rday, rhour, rmin, rsec, rmili
    integer       ::  status

    ! --- begin ---------------------------------

    call Check( dt, status )
    IF_NOT_OK_STOP

    ! multiply each part, at remainder of previous:
    rday  = dt%day  * r
    rhour = dt%hour * r + (rday -nint(rday ))*24
    rmin  = dt%minute  * r + (rhour-nint(rhour))*60
    rsec  = dt%second  * r + (rmin -nint(rmin ))*60
    rmili = dt%milisecond * r + (rsec -nint(rsec ))*1000

    ! fill:
    dt_times_r8 = CSO_TimeDelta( day  = nint(rday ), &
                            hour = nint(rhour), &
                            min  = nint(rmin ), &
                            sec  = nint(rsec ), &
                            mili = nint(rmili)      )

  end function dt_times_r8


  ! *


  type(T_CSO_TimeDelta) function r8_times_dt( r, dt )

    ! --- in/out --------------------------------

    real(8), intent(in)           ::  r
    type(T_CSO_TimeDelta), intent(in)   ::  dt

    ! --- begin ---------------------------------

    r8_times_dt = dt * r

  end function r8_times_dt


  ! *


  type(T_CSO_TimeDelta) function dt_times_i( dt, i )

    ! --- in/out --------------------------------

    type(T_CSO_TimeDelta), intent(in)   ::  dt
    integer, intent(in)           ::  i

    ! --- begin ---------------------------------

    dt_times_i = dt * (i*1.0)

  end function dt_times_i


  ! *


  type(T_CSO_TimeDelta) function i_times_dt( i, dt )

    ! --- in/out --------------------------------

    integer, intent(in)           ::  i
    type(T_CSO_TimeDelta), intent(in)   ::  dt

    ! --- begin ---------------------------------

    i_times_dt = dt * i

  end function i_times_dt


  ! *** dt = dt / r ************************************************


  type(T_CSO_TimeDelta) function dt_div_r( dt, r )

    ! --- in/out --------------------------------

    type(T_CSO_TimeDelta), intent(in)   ::  dt
    real, intent(in)              ::  r

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_div_r'

    ! --- local ---------------------------------

    integer    ::  status
    real       ::  rat
    integer    ::  intg
    real       ::  frac

    ! --- begin ---------------------------------

    call Check( dt, status )
    IF_NOT_OK_STOP

    ! days:
    rat = dt%day / r
    intg = floor( rat )
    frac = rat - intg
    dt_div_r = CSO_TimeDelta( day=intg )

    ! hours:
    rat = dt%hour / r + frac*24
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, hour=intg )

    ! mins:
    rat = dt%minute / r + frac*60
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, min=intg )

    ! seconds:
    rat = dt%second / r + frac*60
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, sec=intg )

    ! miliseconds:
    rat = dt%milisecond / r + frac*1000
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, mili=intg )

  end function dt_div_r

  ! *

  type(T_CSO_TimeDelta) function dt_div_i( dt, i )

    ! --- in/out --------------------------------

    type(T_CSO_TimeDelta), intent(in)   ::  dt
    integer, intent(in)           ::  i

    ! --- begin ---------------------------------

    dt_div_i = dt / (i*1.0)

  end function dt_div_i


  ! ************************************************
  ! ***
  ! *** logical operators
  ! ***
  ! ************************************************


  logical function date_IsAnyDate( t )

    ! --- in/out -------------------------------

    type(T_CSO_DateTime), intent(in)        ::  t

    ! --- begin --------------------------------

    date_IsAnyDate = t%calendar == 'none'

  end function date_IsAnyDate


  ! ***  date1 == date2

  logical function date_eq_date( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_eq_date'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    call Check( t1, status )
    IF_NOT_OK_STOP
    call Check( t2, status )
    IF_NOT_OK_STOP

    ! any date ? always equal
    if ( (t1%calendar == 'none') .or. (t2%calendar == 'none') ) then
      date_eq_date = .true.
      return
    end if

    ! compare values
    date_eq_date = &
        ( t1%year  == t2%year  ) .and. &
        ( t1%month == t2%month ) .and. &
        ( t1%day   == t2%day   ) .and. &
        ( t1%zone  == t2%zone  ) .and. &
        ( t1%hour  == t2%hour  ) .and. &
        ( t1%minute   == t2%minute   ) .and. &
        ( t1%second   == t2%second   ) .and. &
        ( t1%milisecond  == t2%milisecond  )

  end function date_eq_date


  ! ***  date1 /= date2

  logical function date_ne_date( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_ne_date'

    ! --- begin ---------------------------------

    date_ne_date = .not. ( t1 == t2 )

  end function date_ne_date


  ! ***  date1 > date2


  logical function date_gt_date( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_gt_date'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    call Check( t1, status )
    IF_NOT_OK_STOP
    call Check( t2, status )
    IF_NOT_OK_STOP

    ! any date ? always true
    if ( (t1%calendar == 'none') .or. (t2%calendar == 'none') ) then
      date_gt_date = .true.
      return
    end if

    if ( t1%year > t2%year ) then
      date_gt_date = .true.
      return
    else if ( t1%year < t2%year ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%month > t2%month ) then
      date_gt_date = .true.
      return
    else if ( t1%month < t2%month ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%day > t2%day ) then
      date_gt_date = .true.
      return
    else if ( t1%day < t2%day ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%hour > t2%hour ) then
      date_gt_date = .true.
      return
    else if ( t1%hour < t2%hour ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%minute > t2%minute ) then
      date_gt_date = .true.
      return
    else if ( t1%minute < t2%minute ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%second > t2%second ) then
      date_gt_date = .true.
      return
    else if ( t1%second < t2%second ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%milisecond > t2%milisecond ) then
      date_gt_date = .true.
      return
    else if ( t1%milisecond < t2%milisecond ) then
      date_gt_date = .false.
      return
    end if

    ! all fields are equal ...
    date_gt_date = .false.

  end function date_gt_date


  ! ***  date1 < date2


  logical function date_lt_date( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- begin ---------------------------------

    date_lt_date = (.not.( ( t1 == t2 ) .or. ( t1 > t2 ) )) .or. IsAnyDate(t1) .or. IsAnyDate(t2)

  end function date_lt_date


  ! ***  date1 >= date2 ************************


  logical function date_ge_date( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- begin ---------------------------------

    date_ge_date = ( t1 == t2 ) .or. ( t1 > t2 ) .or. IsAnyDate(t1) .or. IsAnyDate(t2)

  end function date_ge_date


  ! ***  date1 <= date2 ************************


  logical function date_le_date( t1, t2 )

    ! --- in/out --------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t1
    type(T_CSO_DateTime), intent(in)       ::  t2

    ! --- begin ---------------------------------

    date_le_date = (.not. ( t1 > t2 )) .or. IsAnyDate(t1) .or. IsAnyDate(t2)

  end function date_le_date


  ! ***********************************************
  ! ***
  ! *** min/max
  ! ***
  ! ***********************************************


  function date_min( t1, t2 )

    ! --- in/out ---------------------------------

    type(T_CSO_DateTime)                   ::  date_min
    type(T_CSO_DateTime), intent(in)       ::  t1, t2

    ! --- begin ----------------------------------

    if ( t1 <= t2 ) then
      date_min = t1
    else
      date_min = t2
    end if

  end function date_min


  ! ***


  function date_max( t1, t2 )

    ! --- in/out ---------------------------------

    type(T_CSO_DateTime)                   ::  date_max
    type(T_CSO_DateTime), intent(in)       ::  t1, t2

    ! --- begin ----------------------------------

    if ( t1 >= t2 ) then
      date_max = t1
    else
      date_max = t2
    end if

  end function date_max


  ! ***********************************************
  ! ***
  ! *** totals
  ! ***
  ! ***********************************************


  real function date_rTotal( t, unit )

    ! --- in/out ----------------------------

    type(T_CSO_DateTime), intent(in)          ::  t
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_rTotal'

    ! --- local -----------------------------------

    integer            ::  status
    real               ::  nday
    integer            ::  iyear

    ! --- begin -----------------------------

    call Check( t, status )
    IF_NOT_OK_STOP

    ! not all arguments are possible ...
    select case ( t%calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day' )
        select case ( unit )
          case ( 'year' )
            if ( any( (/t%month,t%day,t%hour,t%minute,t%second,t%milisecond/) /= 0 ) ) then
              write (csol,'("do not know how to count total:")'); call csoErr
              write (csol,'("  unit : ",a)') unit; call csoErr
              call wrtcsol( '  t    : ', t ); call csoErr
              TRACEBACK; stop
            end if
          case ( 'month' )
            if ( any( (/t%day,t%hour,t%minute,t%second,t%milisecond/) /= 0 ) ) then
              write (csol,'("do not know how to count total:")'); call csoErr
              write (csol,'("  unit : ",a)') unit; call csoErr
              call wrtcsol( '  t    : ', t ); call csoErr
              TRACEBACK; stop
            end if
        end select
      case ( 'incr' )
        select case ( unit )
          case ( 'year', 'month' )
            write (csol,'("do not know how to count total in incremental date:")') unit; call csoErr
            write (csol,'("  unit : ",a)') unit; call csoErr
            call wrtcsol( '  t    : ', t ); call csoErr
            TRACEBACK; stop
        end select
    end select

    ! precount total number of days for some of the units:
    select case ( unit )
      case ( 'day', 'hour', 'min', 'sec', 'mili' )
        nday = 0.0
        do iyear = 1, t%year-1
          nday  = nday  + calc_days_in_year(t%calendar,iyear)
        end do
        nday  = nday + DayNumber( t ) - 1
      case default
        write (csol,'("unsupported unit: ",a)') trim(unit); call csoErr
        TRACEBACK; stop
    end select

    ! count time units:
    select case ( unit )
      case ( 'year' )
        ! set 'nday' to a reference length of the year;
        ! if this length is not constant during the years, the
        ! values of t%month etc have been checked to be zero:
        nday = days_in_year(t) * 1.0
        ! count fractional years:
        date_rTotal = t%year                                           + &
                      t%month / 12.0                                   + &
                      t%day   / nday                                   + &
                      t%hour  / nday / 24.0                            + &
                      t%minute   / nday / 24.0 / 60.0                     + &
                      t%second   / nday / 24.0 / 60.0 / 60.0              + &
                      t%milisecond  / nday / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'month' )
        ! set 'nday' to a reference length of the month;
        ! if this length is not constant during the years, the
        ! values of t%day etc been checked to be zero:
        nday = days_in_month(t) * 1.0
        ! count fractional months:
        date_rTotal = t%year  * 12.0                                   + &
                      t%month                                          + &
                      t%day   / nday                                   + &
                      t%hour  / nday / 24.0                            + &
                      t%minute   / nday / 24.0 / 60.0                     + &
                      t%second   / nday / 24.0 / 60.0 / 60.0              + &
                      t%milisecond  / nday / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        date_rTotal = nday                                             + &
                      t%hour         / 24.0                            + &
                      t%minute          / 24.0 / 60.0                     + &
                      t%second          / 24.0 / 60.0 / 60.0              + &
                      t%milisecond         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        date_rTotal = nday           * 24.0                            + &
                      t%hour                                           + &
                      t%minute                 / 60.0                     + &
                      t%second                 / 60.0 / 60.0              + &
                      t%milisecond                / 60.0 / 60.0 / 1000.0
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        date_rTotal = nday           * 24.0 * 60.0                     + &
                      t%hour                * 60.0                     + &
                      t%minute                                            + &
                      t%second                        / 60.0              + &
                      t%milisecond                       / 60.0 / 1000.0
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        date_rTotal = nday           * 24.0 * 60.0 * 60.0              + &
                      t%hour                * 60.0 * 60.0              + &
                      t%minute                        * 60.0              + &
                      t%second                                            + &
                      t%milisecond                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        date_rTotal = nday           * 24.0 * 60.0 * 60.0 * 1000.0     + &
                      t%hour                * 60.0 * 60.0 * 1000.0     + &
                      t%minute                        * 60.0 * 1000.0     + &
                      t%second                               * 1000.0     + &
                      t%milisecond
      case default
        write (csol,'("do not know how to count time in unit : ",a)') trim(unit); call csoErr
        TRACEBACK; stop
    end select

  end function date_rTotal


  ! ***


  real function incr_rTotal( dt, unit )

    ! --- in/out ----------------------------

    type(T_CSO_TimeDelta), intent(in)      ::  dt
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incr_rTotal'

    ! --- local -----------------------------------

    integer            ::  status

    ! --- begin -----------------------------

    call Check( dt, status )
    IF_NOT_OK_STOP

    ! count time units:
    select case ( unit )
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        incr_rTotal = dt%day                                            + &
                      dt%hour         / 24.0                            + &
                      dt%minute          / 24.0 / 60.0                     + &
                      dt%second          / 24.0 / 60.0 / 60.0              + &
                      dt%milisecond         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        incr_rTotal = dt%day          * 24.0                            + &
                      dt%hour                                           + &
                      dt%minute                 / 60.0                     + &
                      dt%second                 / 60.0 / 60.0              + &
                      dt%milisecond                / 60.0 / 60.0 / 1000.0
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        incr_rTotal = dt%day          * 24.0 * 60.0                     + &
                      dt%hour                * 60.0                     + &
                      dt%minute                                            + &
                      dt%second                        / 60.0              + &
                      dt%milisecond                       / 60.0 / 1000.0
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        incr_rTotal = dt%day          * 24.0 * 60.0 * 60.0              + &
                      dt%hour                * 60.0 * 60.0              + &
                      dt%minute                        * 60.0              + &
                      dt%second                                            + &
                      dt%milisecond                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        incr_rTotal = dt%day          * 24.0 * 60.0 * 60.0 * 1000.0     + &
                      dt%hour                * 60.0 * 60.0 * 1000.0     + &
                      dt%minute                        * 60.0 * 1000.0     + &
                      dt%second                               * 1000.0     + &
                      dt%milisecond
      case default
        write (csol,'("do not know how to count time in unit : ",a)') trim(unit); call csoErr
        TRACEBACK; stop
    end select

  end function incr_rTotal


  ! ***


  integer function date_iTotal( t, unit )

    ! --- in/out ----------------------------

    type(T_CSO_DateTime), intent(in)          ::  t
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_iTotal'

    ! --- local -----------------------------

    integer            ::  status
    real               ::  rtot
    integer            ::  itot

    ! --- begin -----------------------------

    call Check( t, status )
    IF_NOT_OK_STOP

    ! determine total some as a real value:
    rtot = rTotal( t, unit )

    ! round to integer value:
    itot = nint(rtot)

    ! result should be pure integer ....
    if ( itot*1.0 == rtot ) then
      date_iTotal = itot
    else
      write (csol,'("date does not contain integer total:")'); call csoErr
      write (csol,'("  unit : ",a)') trim(unit); call csoErr
      call wrtcsol( '  t    : ', t ); call csoErr
      TRACEBACK; stop
    end if

  end function date_iTotal


  ! ***


  real(8) function date_dTotal( t, unit )

    ! --- in/out ----------------------------

    type(T_CSO_DateTime), intent(in)          ::  t
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_dTotal'

    ! --- local -----------------------------------

    integer            ::  status
    real(8)            ::  nday
    integer            ::  iyear

    ! --- begin -----------------------------

    call Check( t, status )
    IF_NOT_OK_STOP

    ! not all arguments are possible ...
    select case ( t%calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day' )
        select case ( unit )
          case ( 'year' )
            if ( any( (/t%month,t%day,t%hour,t%minute,t%second,t%milisecond/) /= 0 ) ) then
              write (csol,'("do not know how to count total:")'); call csoErr
              write (csol,'("  unit : ",a)') unit; call csoErr
              call wrtcsol( '  t    : ', t ); call csoErr
              TRACEBACK; stop
            end if
          case ( 'month' )
            if ( any( (/t%day,t%hour,t%minute,t%second,t%milisecond/) /= 0 ) ) then
              write (csol,'("do not know how to count total:")'); call csoErr
              write (csol,'("  unit : ",a)') unit; call csoErr
              call wrtcsol( '  t    : ', t ); call csoErr
              TRACEBACK; stop
            end if
        end select
      case ( 'incr' )
        select case ( unit )
          case ( 'year', 'month' )
            write (csol,'("do not know how to count total in incremental date:")'); call csoErr
            write (csol,'("  unit : ",a)') unit; call csoErr
            call wrtcsol( '  t    : ', t ); call csoErr
            TRACEBACK; stop
        end select
    end select

    ! precount total number of days for some of the units:
    select case ( unit )
      case ( 'day', 'hour', 'min', 'sec', 'mili' )
        nday = 0.0
        do iyear = 1, t%year-1
          nday  = nday  + calc_days_in_year(t%calendar,iyear)
        end do
        nday  = nday + DayNumber( t ) - 1
      case default
        write (csol,'("unsupported unit: ",a)') trim(unit); call csoErr
        TRACEBACK; stop
    end select

    ! count time units:
    select case ( unit )
      case ( 'year' )
        ! set 'nday' to a reference length of the year;
        ! if this length is not constant during the years, the
        ! values of t%month etc have been checked to be zero:
        nday = days_in_year(t) * 1.0
        ! count fractional years:
        date_dTotal = t%year                                           + &
                      t%month / 12.0                                   + &
                      t%day   / nday                                   + &
                      t%hour  / nday / 24.0                            + &
                      t%minute   / nday / 24.0 / 60.0                     + &
                      t%second   / nday / 24.0 / 60.0 / 60.0              + &
                      t%milisecond  / nday / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'month' )
        ! set 'nday' to a reference length of the month;
        ! if this length is not constant during the years, the
        ! values of t%day etc been checked to be zero:
        nday = days_in_month(t) * 1.0
        ! count fractional months:
        date_dTotal = t%year  * 12.0                                   + &
                      t%month                                          + &
                      t%day   / nday                                   + &
                      t%hour  / nday / 24.0                            + &
                      t%minute   / nday / 24.0 / 60.0                     + &
                      t%second   / nday / 24.0 / 60.0 / 60.0              + &
                      t%milisecond  / nday / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        date_dTotal = nday                                             + &
                      t%hour         / 24.0                            + &
                      t%minute          / 24.0 / 60.0                     + &
                      t%second          / 24.0 / 60.0 / 60.0              + &
                      t%milisecond         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        date_dTotal = nday           * 24.0                            + &
                      t%hour                                           + &
                      t%minute                 / 60.0                     + &
                      t%second                 / 60.0 / 60.0              + &
                      t%milisecond                / 60.0 / 60.0 / 1000.0
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        date_dTotal = nday           * 24.0 * 60.0                     + &
                      t%hour                * 60.0                     + &
                      t%minute                                            + &
                      t%second                        / 60.0              + &
                      t%milisecond                       / 60.0 / 1000.0
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        date_dTotal = nday           * 24.0 * 60.0 * 60.0              + &
                      t%hour                * 60.0 * 60.0              + &
                      t%minute                        * 60.0              + &
                      t%second                                            + &
                      t%milisecond                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        date_dTotal = nday           * 24.0 * 60.0 * 60.0 * 1000.0     + &
                      t%hour                * 60.0 * 60.0 * 1000.0     + &
                      t%minute                        * 60.0 * 1000.0     + &
                      t%second                               * 1000.0     + &
                      t%milisecond
      case default
        write (csol,'("do not know how to count time in unit : ",a)') trim(unit); call csoErr
        TRACEBACK; stop
    end select

  end function date_dTotal


  ! ***


  real(8) function incr_dTotal( dt, unit )

    ! --- in/out ----------------------------

    type(T_CSO_TimeDelta), intent(in)      ::  dt
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incr_rTotal'

    ! --- local -----------------------------------

    integer            ::  status

    ! --- begin -----------------------------

    call Check( dt, status )
    IF_NOT_OK_STOP

    ! count time units:
    select case ( unit )
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        incr_dTotal = dt%day                                            + &
                      dt%hour         / 24.0                            + &
                      dt%minute          / 24.0 / 60.0                     + &
                      dt%second          / 24.0 / 60.0 / 60.0              + &
                      dt%milisecond         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        incr_dTotal = dt%day          * 24.0                            + &
                      dt%hour                                           + &
                      dt%minute                 / 60.0                     + &
                      dt%second                 / 60.0 / 60.0              + &
                      dt%milisecond                / 60.0 / 60.0 / 1000.0
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        incr_dTotal = dt%day          * 24.0 * 60.0                     + &
                      dt%hour                * 60.0                     + &
                      dt%minute                                            + &
                      dt%second                        / 60.0              + &
                      dt%milisecond                       / 60.0 / 1000.0
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        incr_dTotal = dt%day          * 24.0 * 60.0 * 60.0              + &
                      dt%hour                * 60.0 * 60.0              + &
                      dt%minute                        * 60.0              + &
                      dt%second                                            + &
                      dt%milisecond                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        incr_dTotal = dt%day          * 24.0 * 60.0 * 60.0 * 1000.0     + &
                      dt%hour                * 60.0 * 60.0 * 1000.0     + &
                      dt%minute                        * 60.0 * 1000.0     + &
                      dt%second                               * 1000.0     + &
                      dt%milisecond
      case default
        write (csol,'("do not know how to count time in unit : ",a)') trim(unit); call csoErr
        TRACEBACK; stop
    end select

  end function incr_dTotal

  ! ***  
  
  integer function incr_iTotal( dt, unit )

    ! --- in/out ----------------------------

    type(T_CSO_TimeDelta), intent(in)      ::  dt
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incr_iTotal'

    ! --- local -----------------------------------

    integer            ::  status

    ! --- begin -----------------------------

    call Check( dt, status )
    IF_NOT_OK_STOP

    ! count time units:
    select case ( unit )
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        incr_iTotal = dt%day
        
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        incr_iTotal = dt%day          * 24                      + &
                      dt%hour
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        incr_iTotal = dt%day          * 24 * 60                 + &
                      dt%hour              * 60                 + &
                      dt%minute
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        incr_iTotal = dt%day          * 24 * 60 * 60            + &
                      dt%hour              * 60 * 60            + &
                      dt%minute                    * 60            + &
                      dt%second
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        incr_iTotal = dt%day          * 24 * 60 * 60 * 1000     + &
                      dt%hour              * 60 * 60 * 1000     + &
                      dt%minute                    * 60 * 1000     + &
                      dt%second                         * 1000     + &
                      dt%milisecond
      case default
        write (csol,'("do not know how to count time in unit : ",a)') trim(unit); call csoErr
        TRACEBACK; stop
    end select

  end function incr_iTotal



  ! ***********************************************
  ! ***
  ! *** interpolation
  ! ***
  ! ***********************************************

  !
  ! Return coeff such that
  !   t = alfa1 * t1 + alfa2 * t2
  !

  subroutine date_InterpolFractions( t, t1, t2, alfa1, alfa2, status )

    ! --- in/out -----------------------------

    type(T_CSO_DateTime), intent(in)    ::  t
    type(T_CSO_DateTime), intent(in)    ::  t1
    type(T_CSO_DateTime), intent(in)    ::  t2
    real, intent(out)          ::  alfa1
    real, intent(out)          ::  alfa2
    integer, intent(out)       ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_InterpolFractions'

    ! --- local ------------------------------

    real    ::  ds, ds1

    ! --- begin ------------------------------

    ! check ...
    if ( t1 > t2 ) then
      write (csol,'("required interval [t1,t2] :")'); call csoErr
      call wrtcsol( '  t1 = ', t1 ); call csoErr
      call wrtcsol( '  t2 = ', t2 ); call csoErr
      TRACEBACK; status=1; return
    end if

    ! check ...
    if ( (t < t1) .or. (t > t2) ) then
      write (csol,'("t not in [t1,t2] :")'); call csoErr
      call wrtcsol( '  t  = ', t  ); call csoErr
      call wrtcsol( '  t1 = ', t1 ); call csoErr
      call wrtcsol( '  t2 = ', t2 ); call csoErr
      TRACEBACK; status=1; return
    end if

    ! compute differences in seconds:
    ds  = rTotal( t2 - t1, 'sec' )
    ds1 = rTotal( t  - t1, 'sec' )

    ! return fractions
    if ( abs(ds) < tiny(ds) ) then
      alfa2 = 0.5
    else
      alfa2 = ds1 / ds
    end if
    alfa1 = 1.0 - alfa2

    ! ok
    status = 0

  end subroutine date_InterpolFractions


  ! ***********************************************
  ! ***
  ! *** intervals
  ! ***
  ! ***********************************************

  !
  ! Return interval [tt(1),tt(2)] around t that matches with time resolution;
  ! resolution specified by a step and unit:
  !    3, 'hour'   # 00:00, 03:00, 06:00, ...
  !   24, 'hour from 12'   # 12:00, 36:00, ...
  ! Interval boundaries are the same if t is exactly on a step.
  !

  subroutine Get_Surrounding_Interval( t, tref, interpolation, tt, status )
  
    use CSO_String, only : CSO_ReadFromLine, CSO_SplitString

    ! --- in/out -----------------------------

    type(T_CSO_DateTime), intent(in)       ::  t
    type(T_CSO_DateTime), intent(in)       ::  tref
    character(len=*), intent(in)  ::  interpolation
    type(T_CSO_DateTime), intent(out)      ::  tt(2)
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Get_Surrounding_Interval'

    ! --- local ------------------------------

    character(len=64)   ::  line
    integer             ::  ninterp
    character(len=64)   ::  interp(2)
    character(len=32)   ::  units
    character(len=32)   ::  from, after
    integer             ::  step
    integer             ::  offset
    integer             ::  tchange
    real                ::  r
    type(T_CSO_DateTime)         ::  t0

    ! --- begin ------------------------------
    
    if ( index(trim(interpolation), ';') > 0 ) then
      
      ! different units over time period
      ! split:
      call CSO_SplitString( interpolation, ninterp, interp, status, sep=';' )
      IF_NOT_OK_RETURN(status=1)
      
      if ( ninterp > 2 ) then
        write (csol, '("More than 2 time inervals found: not supported")' ) ; call csoErr
        TRACEBACK;status=1;return
      end if
      ! Read second part:
      ! format: * hour after *
      ! split at after
      call CSO_ReadFromLine( interp(2), step, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      call CSO_ReadFromLine( interp(2), units, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      call CSO_ReadFromLine( interp(2), after, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      call CSO_ReadFromLine( interp(2), tchange, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      offset = 0
      
      ! Check if time is before offset
      if ( itotal( t - tref, 'hour' ) < tchange ) then
        ! read first inerpolation part
        call CSO_ReadFromLine( interp(1), step, status, sep=' ' )
        IF_NOT_OK_RETURN(status=1)
        call CSO_ReadFromLine( interp(1), units, status, sep=' ' )
        IF_NOT_OK_RETURN(status=1)
        offset = 0
      endif
    else ! read first inerpolation part
      line = trim(interpolation)
      call CSO_ReadFromLine( line, step, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      call CSO_ReadFromLine( line, units, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      offset = 0
    endif
      
    ! which units ?
    select case ( units )

      ! x hourly
      case ( 'hour' )
        ! begin of day:
        t0 = Get_Begin_Of( t, 'day' )
        ! add offset:
        t0 = t0 + CSO_TimeDelta(hour=offset)
        ! back one day?
        if ( t0%hour > t%hour ) t0 = t0 - CSO_TimeDelta(hour=24)
        ! hour fraction:
        r = rTotal( t - t0, 'hour' )
        ! start of interval:
        tt(1) = t0 + CSO_TimeDelta( hour=int(floor(r/step))*step )
        ! end:
        if ( t == tt(1) ) then
          tt(2) = tt(1)
        else
          tt(2) = tt(1) + CSO_TimeDelta(hour=step)
        end if

      ! unknown ...
      case default
        write (csol,'("unsupported units `",a,"`")') trim(units); call csoErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine Get_Surrounding_Interval


  ! ***********************************************
  ! ***
  ! *** units
  ! ***
  ! ***********************************************

  ! extract time:
  !   1900-01-01 00:00:0.0

  subroutine CSO_ReadFromLine_t( input, t, status, sep, calendar )

    use CSO_String, only : CSO_ReadFromLine
    use CSO_String, only : CSO_Translate

    ! --- in/out ---------------------------------

    character(len=*), intent(inout)           ::  input
    type(T_CSO_DateTime), intent(out)                  ::  t
    integer, intent(out)                      ::  status
    character(len=1), intent(in), optional    ::  sep
    character(len=*), intent(in), optional    ::  calendar

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/CSO_ReadFromLine_t'

    ! --- local ----------------------------------

    character(len=32)     ::  line
    integer               ::  year, month, day, hour, minu
    real                  ::  seco

    ! --- begin ----------------------------------

    ! extract:
    call CSO_ReadFromLine( input, line, status, sep=sep )
    IF_NOT_OK_RETURN(status=1)
    
    ! replace seperators:
    call CSO_Translate( line, '-/:_T', ' ', status )
    IF_NOT_OK_RETURN(status=1)

    ! extract ref time values:
    call CSO_ReadFromLine( line, year, status, sep=' ' )
    IF_NOT_OK_RETURN(status=1)
    call CSO_ReadFromLine( line, month, status, sep=' ' )
    IF_NOT_OK_RETURN(status=1)
    call CSO_ReadFromLine( line, day, status, sep=' ' )
    IF_NOT_OK_RETURN(status=1)
    if ( len_trim(line) > 0 ) then
      call CSO_ReadFromLine( line, hour, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
    else
      hour = 0
    end if
    if ( len_trim(line) > 0 ) then
      call CSO_ReadFromLine( line, minu, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
    else
      minu = 0
    end if
    if ( len_trim(line) > 0 ) then
      call CSO_ReadFromLine( line, seco, status )
      IF_NOT_OK_RETURN(status=1)
    else
      seco = 0.0
    end if

    ! store:
    t = CSO_DateTime( year=year, month=month, day=day, &
                    hour=hour, min=minu, sec=nint(seco), &
                    calendar=calendar )

    ! ok
    status = 0

  end subroutine CSO_ReadFromLine_t
  
  ! *

  !
  ! Expand time value using units:
  !
  !   hours since 1900-01-01 00:00:0.0
  !
  ! Also step 'months' are accepted.
  !

  subroutine ExpandTime( value, units, calendar, t, status )

    use CSO_String, only : CSO_ReadFromLine

    ! --- in/out ---------------------------------

    real, intent(in)                ::  value
    character(len=*), intent(in)    ::  units
    character(len=*), intent(in)    ::  calendar
    type(T_CSO_DateTime), intent(out)        ::  t
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/ExpandTime'

    ! --- local ----------------------------------

    character(len=512)    ::  line
    character(len=32)     ::  steps
    character(len=32)     ::  dummy
    type(T_CSO_DateTime)  ::  tref
    integer               ::  nyear, nmonth, nday
    real                  ::  rmonth

    ! --- begin ----------------------------------

    ! copy:
    line = trim(units)

    ! extract step units:
    call CSO_ReadFromLine( line, steps, status, sep=' ' )
    IF_NOT_OK_RETURN(status=1)

    ! extract 'since'
    call CSO_ReadFromLine( line, dummy, status, sep=' ' )
    IF_NOT_OK_RETURN(status=1)
    ! check ...
    if ( trim(dummy) /= 'since' ) then
      write (csol,'("second field should be `since`, not `",a,"`")') trim(dummy); call csoPr
      TRACEBACK; status=1; return
    end if

    ! extract ref time values:
    call CSO_ReadFromLine( line, tref, status, calendar=calendar )
    IF_NOT_OK_RETURN(status=1)

    ! expand:
    select case ( trim(steps) )
    
      case ( 'month', 'months' )
        ! number of years and full months:
        nyear  = int(value/12)
        rmonth = value - nyear*12
        nmonth = int(rmonth)
        ! init result:
        t = CSO_DateTime( year=tref%year+nyear, month=tref%month+nmonth, day=1 )
        ! month:
        nday = Days_in_Month( t )
        ! remaining month fraction:
        rmonth = rmonth - nmonth
        ! add rest:
        t = t + CSO_TimeDelta( hour=int(rmonth*nday*24) )

      case ( 'day', 'days' )
        t = tref + CSO_TimeDelta(day=1) * value
      case ( 'hour', 'hours' )
        t = tref + CSO_TimeDelta(hour=1) * value
      case ( 'minute', 'minutes' )
        t = tref + CSO_TimeDelta(min=1) * value
      case ( 'second', 'seconds' )
        t = tref + CSO_TimeDelta(sec=1) * value

      case default
        write (csol,'("unsupported steps `",a,"`")') trim(steps); call csoPr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine ExpandTime


  ! ***


  !
  ! Compare 't' with a numeric value given its units
  !  ("days since ...") and the calendar ("366_day")
  !
  ! Return status:
  !   <0    : time and value do not match
  !    0    : time and value match
  !   >0    : some error
  !

  subroutine Compare_Date_Num( t, num, units, calendar, status )

    ! --- in/out ---------------------------------

    type(T_CSO_DateTime), intent(in)       ::  t
    real, intent(in)              ::  num
    character(len=*), intent(in)  ::  units
    character(len=*), intent(in)  ::  calendar
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Compare_Date_Num'

    ! --- local ----------------------------------

    type(T_CSO_DateTime)           ::  t_rec

    ! --- begin ----------------------------------

    ! check ...
    if ( t%calendar /= calendar ) then
      write (csol,'("calendar of t is `",a,"` instead of `",a,"`")') trim(t%calendar), trim(calendar); call csoErr
      TRACEBACK; status=1; return
    end if

    ! expand:
    call ExpandTime( num, units, calendar, t_rec, status )
    IF_NOT_OK_RETURN(status=1)
    ! compare:
    if ( t_rec /= t ) then
      status = -1; return
    end if

    ! ok
    status = 0

  end subroutine Compare_Date_Num


  ! ***********************************************
  ! ***
  ! *** print
  ! ***
  ! ***********************************************


  character(len=36) function date_Pretty( t )

    ! --- in/out -------------------------

    type(T_CSO_DateTime), intent(in)       ::  t

    ! --- const --------------------------

    character(len=3), parameter  ::  month_name(12) = &
                                      (/ 'jan','feb','mar','apr','may','jun', &
                                         'jul','aug','sep','oct','nov','dec' /)

    ! --- local --------------------------

    character(len=36)         ::  s

    ! --- begin --------------------------

    select case ( t%calendar )
      !case ( 'proleptic_gregorian', 'gregorian', 'standard' )
      !  if ( t%zone < 0 ) then
      !    zone_sign = '-'
      !  else
      !    zone_sign = '+'
      !  end if
      !  zone_abs  = abs(t%zone)
      !  zone_hour = floor(zone_abs/60.0)
      !  zone_min  = zone_abs - zone_hour*60
      !  write (s,'(i2,":",i2.2,":",i2.2,":",i3.3,       &
      !            & " ",i2.2," ",a3," ",i4.4,           &
      !            & " (GMT",a1,i2.2,":",i2.2,")")')     &
      !    t%hour, t%minute, t%second, t%milisecond, &
      !    t%day, month_name(t%month), t%year, &
      !    zone_sign, zone_hour, zone_min
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day', '360_day', 'none' )
        write (s,'(i4.4,"-",i2.2,"-",i2.2," ",i2,":",i2.2,":",i2.2)') &
          t%year, t%month, t%day, t%hour, t%minute, t%second
      case default
        s = 'no-calendar'
    end select

    date_Pretty = s

  end function date_Pretty


  ! *


  character(len=36) function incrdate_Pretty( dt )

    ! --- in/out -------------------------

    type(T_CSO_TimeDelta), intent(in)       ::  dt

    ! --- local --------------------------

    character(len=36)       ::  s

    ! --- begin --------------------------

    write (s,'(i5," days ",i2,":",i2.2,":",i2.2,":",i3.3)') &
      dt%day, dt%hour, dt%minute, dt%second, dt%milisecond

    incrdate_Pretty = s

  end function incrdate_Pretty


  ! *


  subroutine wrtcsol_t( msg, t )

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(T_CSO_DateTime), intent(in)         ::  t

    ! --- local ---------------------------------

    character(len=36)     ::  s

    ! --- begin -----------------------------------

    s = date_Pretty( t )
    write (csol,'(a,a)') msg, trim(s)

  end subroutine wrtcsol_t


  ! *


  subroutine wrtcsol_dt( msg, dt )

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(T_CSO_TimeDelta), intent(in)     ::  dt

    ! --- local ---------------------------------

    character(len=36)     ::  s

    ! --- begin -----------------------------------

    s = incrdate_Pretty( dt )
    write (csol,'(a,a)') msg, trim(s)

  end subroutine wrtcsol_dt


  ! *


  subroutine wrtcsol_t_dt( msg, t, msg2, dt )

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(T_CSO_DateTime), intent(in)         ::  t
    character(len=*), intent(in)    ::  msg2
    type(T_CSO_TimeDelta), intent(in)     ::  dt

    ! --- local ---------------------------------

    character(len=36)     ::  s
    character(len=36)     ::  s2

    ! --- begin -----------------------------------

    s  = date_Pretty( t  )
    s2 = incrdate_Pretty( dt )
    write (csol,'(a,a,a,a)') msg, trim(s), msg2, trim(s2)

  end subroutine wrtcsol_t_dt


  ! *


  subroutine wrtcsol_tt( msg, tt )

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(T_CSO_DateTime), intent(in)         ::  tt(2)

    ! --- local ---------------------------------

    character(len=36)     ::  s1
    character(len=36)     ::  s2

    ! --- begin -----------------------------------

    s1 = date_Pretty( tt(1) )
    s2 = date_Pretty( tt(2) )
    write (csol,'(a,"[",a,",",a,"]")') msg, trim(s1), trim(s2)

  end subroutine wrtcsol_tt


  ! *


  subroutine wrtcsol_t1_t2( msg, t, msg2, t2 )

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(T_CSO_DateTime), intent(in)         ::  t
    character(len=*), intent(in)    ::  msg2
    type(T_CSO_DateTime), intent(in)         ::  t2

    ! --- local ---------------------------------

    character(len=36)     ::  s
    character(len=36)     ::  s2

    ! --- begin -----------------------------------

    s  = date_Pretty( t  )
    s2 = date_Pretty( t2 )
    write (csol,'(a,a,a,a)') msg, trim(s), msg2, trim(s2)

  end subroutine wrtcsol_t1_t2


  ! *


  subroutine wrtcsol_t1_t2_t3( msg, t, msg2, t2, msg3, t3 )

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(T_CSO_DateTime), intent(in)         ::  t
    character(len=*), intent(in)    ::  msg2
    type(T_CSO_DateTime), intent(in)         ::  t2
    character(len=*), intent(in)    ::  msg3
    type(T_CSO_DateTime), intent(in)         ::  t3

    ! --- local ---------------------------------

    character(len=36)     ::  s
    character(len=36)     ::  s2
    character(len=36)     ::  s3

    ! --- begin -----------------------------------

    s  = date_Pretty( t  )
    s2 = date_Pretty( t2 )
    s3 = date_Pretty( t3 )
    write (csol,'(a,a,a,a,a,a)') msg, trim(s), msg2, trim(s2), msg3, trim(s3)

  end subroutine wrtcsol_t1_t2_t3


  ! *


  subroutine PrintDate( msg, t )

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(T_CSO_DateTime), intent(in)         ::  t

    ! --- local ---------------------------------

    character(len=36)     ::  s

    ! --- begin -----------------------------------

    s = date_Pretty( t )
    write (*,'(a,a)') msg, trim(s)

  end subroutine PrintDate


end module CSO_DateTimes


