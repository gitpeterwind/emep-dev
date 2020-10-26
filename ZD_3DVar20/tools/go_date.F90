!###############################################################################
!
! go_date  -  module to manipulate date structures
!
! TYPES
!
!  A structure is provided to store a date:
!
!   ! declare date types:
!   type(TDate)     ::  t0, t1, t, dt
!   type(TIncrDate) ::  dt
!
!  with fields:
!
!    character(len=32)   ::  calendar     ! see 'CALENDARS'
!
!    integer             ::  year, month, day, hour, min, sec, mili
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
!   Use routine 'NewDate' to initialize some fields and to fill
!   the rest with zero's. If no calendar is specified,
!   the default value 'gregorian' is used (see also DEFAULTS).
!
!     t = NewDate( calendar='gregorian', year=2000, month=1, ... )
!
!   Use routine 'IncrDate' to create a new increment:
!
!     dt = IncrDate( year=2000, month=1 )
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
!     time5 = (/year,month,day,hour,min/)
!     time6 = (/year,month,day,hour,min,sec/)
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
!   'hour', 'min', 'sec', and 'mili'. If the total number is
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
!  Expand time value from number to TDate using the units description;
!  also 'months' is allowed as step:
!
!    call ExpandTime( 2.5, 'hours since 1900-01-01 00:00:0.0', t, status )
!
!  Extract ref time and time step from units:
!
!    call Extract_Ref_and_Step( 'hours since 1900-01-01 00:00:0.0', t, dt, status,
!                                 calendar='gregorian' )
!
!  Compare 't' of type 'TDate' with a numeric 'value' given its units
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
!     call wrtgol( 'time       : ', t ); call goPr
!     call wrtgol( 'time step  : ', dt ); call goPr
!     call wrtgol( 'time range : ', t1, ' to ', t2 ); call goPr
!
!
! DEFAULTS
!
!   For setting some default values, the subroutine 'go_DateDefaults'
!   is available. All arguments are optional:
!
!     call go_DateDefaults( [calendar='standard'] )
!
!
!### macro's ##################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!###############################################################################


module GO_Date

  use GO_Print, only : gol, goErr, goPr

  implicit none

  ! --- in/out ---------------------------

  private

  public      ::  TDate, TIncrDate

  public      ::  goDateDefaults

  public      ::  NewDate, IncrDate, AnyDate, SystemDate
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

  public      ::  goReadFromLine
  public      ::  Extract_Ref_and_Step
  public      ::  ExpandTime
  public      ::  Compare_Date_Num

  public      ::  Pretty
  public      ::  wrtgol
  public      ::  PrintDate

  ! --- const -----------------------------------

  character(len=*), parameter  ::  mname = 'GO_Date'

  ! --- types -------------------------------------

  ! Strucure with fields to store year, month, day,
  ! hour and minute.
  ! Operators for assignment (=), adding (+),
  ! and comparission (==,<,>,>= and <=)
  ! have been defined for operations between
  ! instances of this type.

  type TDate
    ! type of calendar: 'gregorian', 'standard', '365_day', '360_day'
    character(len=32)   ::  calendar
    ! year, month etc:
    integer             ::  year, month, day, hour, min, sec, mili
    ! difference with Coordinated Universal Time (UTC)
    integer             ::  zone   ! minutes
    ! error status
    integer             ::  status ! = 1
  end type TDate


  type TIncrDate
    ! days, hours, etc:
    integer             ::  day, hour, min, sec, mili
    ! error status
    integer             ::  status ! = 1
  end type TIncrDate

  !! testing ...
  !type DateTime_Interval
  !  ! boundaries:
  !  type(TDate)         ::  left, right
  !  ! open ? default is closed:
  !  logical             ::  left_open, right_open
  !end type DateTime_Interval


  ! --- var --------------------------------

  ! default calendar type
  character(len=32)      ::  default_calendar = 'standard'


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

  ! *

  interface NewDate
    module procedure date_NewDate
  end interface

  interface AnyDate
    module procedure date_AnyDate
  end interface

  interface IncrDate
    module procedure incrdate_IncrDate
  end interface

  interface SystemDate
    module procedure date_SystemDate
  end interface

  interface Get_Begin_Of
    module procedure date_Get_Begin_Of
  end interface

  interface Get_End_Of
    module procedure date_Get_End_Of
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
  
  interface goReadFromLine
    module procedure goReadFromLine_t
  end interface

  interface Pretty
    module procedure date_Pretty
    module procedure incrdate_Pretty
  end interface

  interface wrtgol
    module procedure wrtgol_t
    module procedure wrtgol_dt
    module procedure wrtgol_t_dt
    module procedure wrtgol_tt
    module procedure wrtgol_t1_t2
    module procedure wrtgol_t1_t2_t3
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

    type(TDate), intent(in)        ::  t
    integer, intent(out)           ::  status

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Check'

    ! --- begin -----------------------------------

    ! already error status ? then leave immediatelly:
    if ( t%status /= 0 ) then
      write (gol,'("found error status in date")'); call goErr
      write (gol,'("  year,month,day : ",3i6)') t%year, t%month, t%day; call goErr
      write (gol,'("  hour,minu,sec,mili : ",4i6)') t%hour, t%min, t%sec, t%mili; call goErr
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
          call wrtgol( 'strange month in ', t ); call goErr
          TRACEBACK; status=1; return
        end if
        ! check day
        if ( t%day<1 .or. t%day>days_in_month(t) ) then
          call wrtgol( 'strange day in ', t ); call goErr
          TRACEBACK; status=1; return
        end if
        ! zone should be zero:
        if ( t%zone /= 0 ) then
          call wrtgol( 'expecting zero zone in date ', t ); call goErr
          TRACEBACK; status=1; return
        end if
      case default
        write (gol,'("unknown calendar type: `",a,"`")') t%calendar; call goErr
        write (gol,'("  year etc : ",6i5)') t%year, t%month, t%day, t%hour, t%min, t%sec; call goErr
        TRACEBACK; status=1; return
    end select

    ! check minutes
    if ( t%min<0 .or. t%min>59 ) then
      call wrtgol( 'found strange minutes in ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! check seconds
    if ( t%sec<0 .or. t%sec>59 ) then
      call wrtgol( 'found strange seconds in ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! check mili
    if ( t%mili<0 .or. t%mili>999 ) then
      call wrtgol( 'found strange mili seconds in ', t ); call goErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine date_Check


  ! ***


  subroutine incrdate_Check( dt, status )

    ! --- in/out ----------------------------------

    type(TIncrDate), intent(in)    ::  dt
    integer, intent(out)           ::  status

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incrdate_Check'

    ! --- begin -----------------------------------

    ! already error status ? then leave immediatelly:
    if ( dt%status /= 0 ) then
      write (gol,'("found error status in incrdate")'); call goErr
      write (gol,'("  day, hour,minu,sec,mili : ",5i6)') dt%day, dt%hour, dt%min, dt%sec, dt%mili; call goErr
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
        write (gol,'("unknown calendar type: ",a)') calendar; call goErr
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
        write (gol,'("unknown calendar type: ",a)') calendar; call goErr
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
        write (gol,'("unknown calendar type: ",a)') calendar; call goErr
        TRACEBACK; stop
    end select

  end function calc_DayNumber


  ! **********************************************


  logical function date_LeapYear( t )

    ! --- in/out -------------------------------

    type(TDate), intent(in)         ::  t

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

    type(TDate), intent(in)         ::  t

    ! --- begin --------------------------------

    ! check day ...
    date_LeapDay = (t%month == 2) .and. (t%day == 29)

  end function date_LeapDay


  ! ***


  integer function date_days_in_month( t )

    ! --- in/out ----------------

    type(TDate), intent(in)         ::  t

    ! --- begin ----------------

    date_days_in_month = calc_days_in_month( t%calendar, t%year, t%month )

  end function date_days_in_month


  ! ***


  integer function date_days_in_year( t )

    ! --- in/out ----------------

    type(TDate), intent(in)         ::  t

    ! --- begin ----------------

    date_days_in_year = calc_days_in_year( t%calendar, t%year )

  end function date_days_in_year


  ! ***


  integer function date_DayNumber( t )

    ! --- in/out ----------------

    type(TDate), intent(in)        ::  t

    ! --- begin ----------------

    date_DayNumber = calc_DayNumber( t%calendar, t%year, t%month, t%day )

  end function date_DayNumber


  ! ***


  logical function date_Midnight( t )

    ! --- in/out -------------------------------

    type(TDate), intent(in)         ::  t

    ! --- begin --------------------------------

    date_Midnight = all( (/t%hour,t%min,t%sec,t%mili/) == 0 )

  end function date_Midnight


  ! ***


  logical function date_Precisely( t, val, unit )

    ! --- in/out -------------------------------

    type(TDate), intent(in)         ::  t
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
          write (gol,'("strange value for minu : ",f8.2)') val; call goErr
          TRACEBACK; stop
        end if
        ! total minute fraction in this day:
        rtot = rTotal( IncrDate(0,t%hour,t%min,t%sec,t%mili), 'min' )
        ! modulo given value ?
        date_Precisely = modulo(rtot,val) == 0.0

      ! every .. hour
      case ( 'hour' )
        ! check ...
        if ( (val <= 0.0) .or. (val > 24.0) ) then
          write (gol,'("strange value for hour : ",f8.2)') val; call goErr
          TRACEBACK; stop
        end if
        ! total hour fraction in this day:
        rtot = rTotal( IncrDate(0,t%hour,t%min,t%sec,t%mili), 'hour' )
        ! modulo given value ?
        date_Precisely = modulo(rtot,val) == 0.0

      ! error ...
      case default
        write (gol,'("unsupported unit : ",a)') unit; call goErr
        TRACEBACK; stop

    end select

  end function date_Precisely



  ! ****************************************************
  ! ***
  ! *** Set/Get fields in a structure TDate
  ! ***
  ! ****************************************************

  !
  ! Fill fields in a 'TDate' structure.
  !

  subroutine date_Set( date, year, month, day, hour, min, sec, mili, &
                             zone, calendar, time4, time5, time6 )

    ! --- in/out ------------------------------------

    type(TDate), intent(inout)              ::  date
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
      date%min   = time5(5)
    end if

    if ( present(time6) ) then
      date%year  = time6(1)
      date%month = time6(2)
      date%day   = time6(3)
      date%hour  = time6(4)
      date%min   = time6(5)
      date%sec   = time6(6)
    end if

    if ( present(year ) ) date%year  = year
    if ( present(month) ) date%month = month
    if ( present(day  ) ) date%day   = day
    if ( present(zone ) ) date%zone  = zone
    if ( present(hour ) ) date%hour  = hour
    if ( present(min  ) ) date%min   = min
    if ( present(sec  ) ) date%sec   = sec
    if ( present(mili ) ) date%mili  = mili

  end subroutine date_Set


  ! *

  subroutine incrdate_Set( date, day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    type(TIncrDate), intent(inout)             ::  date
    integer, intent(in), optional            ::  day
    integer, intent(in), optional            ::  hour
    integer, intent(in), optional            ::  min
    integer, intent(in), optional            ::  sec
    integer, intent(in), optional            ::  mili

    ! --- local ----------------------------------

    if ( present(day  ) ) date%day   = day
    if ( present(hour ) ) date%hour  = hour
    if ( present(min  ) ) date%min   = min
    if ( present(sec  ) ) date%sec   = sec
    if ( present(mili ) ) date%mili  = mili

  end subroutine incrdate_Set


  !
  ! Obtain fields from a 'TDate' structure.
  !

  subroutine date_Get( date, &
                       year, month, day, hour, min, sec, mili, &
                       zone, calendar, time4, time5, time6 )

    ! --- in/out ------------------------------------

    type(TDate), intent(in)                  ::  date
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
    if ( present(min  ) ) min   = date%min
    if ( present(sec  ) ) sec   = date%sec
    if ( present(mili ) ) mili  = date%mili

    if ( present(time4) ) time4 = (/ date%year, date%month, date%day, date%hour /)

    if ( present(time5) ) time5 = (/ date%year, date%month, date%day, &
                                     date%hour, date%min   /)

    if ( present(time6) ) time6 = (/ date%year, date%month, date%day, &
                                     date%hour, date%min  , date%sec /)

  end subroutine date_Get

  ! *

  subroutine incrdate_Get( date, day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    type(TIncrDate), intent(in)              ::  date
    integer, intent(out), optional           ::  day
    integer, intent(out), optional           ::  hour, min, sec, mili

    ! --- local ----------------------------------

    if ( present(day  ) ) day   = date%day
    if ( present(hour ) ) hour  = date%hour
    if ( present(min  ) ) min   = date%min
    if ( present(sec  ) ) sec   = date%sec
    if ( present(mili ) ) mili  = date%mili

  end subroutine incrdate_Get


  ! ****************************************************
  ! ***
  ! *** Return a date structure
  ! ***
  ! ****************************************************

  !
  ! Set fields to zero or fill some of them.
  !

  type(TDate) function date_NewDate( year, month, day, hour, min, sec, mili, &
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
    date_NewDate%calendar = default_calendar

    ! Fields are zero by default:
    date_NewDate%year  = 0
    date_NewDate%month = 0
    date_NewDate%day   = 0
    date_NewDate%zone  = 0
    date_NewDate%hour  = 0
    date_NewDate%min   = 0
    date_NewDate%sec   = 0
    date_NewDate%mili  = 0

    ! Optionally, change some of them:
    if ( present(year    ) ) call Set( date_NewDate, year=year )
    if ( present(month   ) ) call Set( date_NewDate, month=month )
    if ( present(day     ) ) call Set( date_NewDate, day=day )
    if ( present(hour    ) ) call Set( date_NewDate, hour=hour )
    if ( present(min     ) ) call Set( date_NewDate, min=min )
    if ( present(sec     ) ) call Set( date_NewDate, sec=sec )
    if ( present(mili    ) ) call Set( date_NewDate, mili=mili )
    if ( present(zone    ) ) call Set( date_NewDate, zone=zone )
    if ( present(calendar) ) call Set( date_NewDate, calendar=calendar )
    if ( present(time4   ) ) call Set( date_NewDate, time4=time4 )
    if ( present(time5   ) ) call Set( date_NewDate, time5=time5 )
    if ( present(time6   ) ) call Set( date_NewDate, time6=time6 )

    ! normalize too small/too large values:
    !if ( date_NewDate%year /= 0000 ) call Normalize( date_NewDate )
    if ( date_NewDate%calendar /= 'none' ) call Normalize( date_NewDate )

    ! data filled, thus probably no error ...
    date_NewDate%status = 0

  end function date_NewDate


  ! ***


  type(TIncrDate) function incrdate_IncrDate( day, hour, min, sec, mili )

    ! --- in/out ------------------------------------

    integer, intent(in), optional            ::  day
    integer, intent(in), optional            ::  hour, min, sec, mili

    ! --- local ----------------------------------

    ! Fields are zero by default:
    incrdate_IncrDate%day   = 0
    incrdate_IncrDate%hour  = 0
    incrdate_IncrDate%min   = 0
    incrdate_IncrDate%sec   = 0
    incrdate_IncrDate%mili  = 0

    ! Optionally, change some of them:
    if ( present(day     ) ) call Set( incrdate_IncrDate, day=day )
    if ( present(hour    ) ) call Set( incrdate_IncrDate, hour=hour )
    if ( present(min     ) ) call Set( incrdate_IncrDate, min=min )
    if ( present(sec     ) ) call Set( incrdate_IncrDate, sec=sec )
    if ( present(mili    ) ) call Set( incrdate_IncrDate, mili=mili )

    ! normalize too small/too large values:
    call Normalize( incrdate_IncrDate )

    ! data filled, thus probably no error ...
    incrdate_IncrDate%status = 0

  end function incrdate_IncrDate


  ! ***


  !
  ! Set fields to zero, special calendar
  !

  type(TDate) function date_AnyDate()

    ! --- local ----------------------------------

    ! Set some fields, other are automatically zero:
    date_AnyDate = NewDate( calendar='none' )

  end function date_AnyDate


  ! ***


  ! Fill with system time

  type(TDate) function date_SystemDate()

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
    call Set( date_SystemDate, calendar='gregorian', &
                year=values(1), month=values(2), day=values(3), &
                zone=values(4), hour=values(5), &
                min=values(6), sec=values(7), mili=values(8) )

    ! no error ...
    date_SystemDate%status = 0

  end function date_SystemDate


  ! ***


  ! set to begin of : 'year', 'month', 'day'

  type(TDate) function date_Get_Begin_Of( t, what )

    ! --- in/out ------------------------------

    type(TDate), intent(in)           ::  t
    character(len=*), intent(in)      ::  what

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Get_Begin_Of'

    ! --- local ------------------------------

    ! --- begin ------------------------------

    ! end of what ?
    select case ( what )

      case ( 'year' )

        ! begin of year:
        date_Get_Begin_Of = NewDate( t%year, 01, 01, 00, 00, 00 )

      case ( 'month' )

        ! begin of month:
        date_Get_Begin_Of = NewDate( t%year, t%month, 01, 00, 00, 00 )

      case ( 'day' )

        ! begin of day:
        date_Get_Begin_Of = NewDate( t%year, t%month, t%day, 00, 00, 00 )

      case default

        write (gol,'("end of `",a,"` not supported yet")') what; call goErr
        TRACEBACK; stop

    end select

    ! no error ...
    date_Get_Begin_Of%status = 0

  end function date_Get_Begin_Of


  ! ***


  ! set to end of : 'month', ...

  type(TDate) function date_Get_End_Of( t, what )

    ! --- in/out ------------------------------

    type(TDate), intent(in)           ::  t
    character(len=*), intent(in)      ::  what

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Get_End_Of'

    ! --- local ------------------------------

    ! --- begin ------------------------------

    ! first set to begin:
    date_Get_End_Of = Get_Begin_Of( t, what )

    ! end of what ?
    select case ( what )
      case ( 'year' )
        ! add number of days:
        date_Get_End_Of = date_Get_End_Of + IncrDate(day=Days_in_Year(t))
      case ( 'month' )
        ! add number of days:
        date_Get_End_Of = date_Get_End_Of + IncrDate(day=Days_in_Month(t))
      case ( 'day' )
        ! add single day:
        date_Get_End_Of = date_Get_End_Of + IncrDate(day=1)
      case default
        write (gol,'("end of `",a,"` not supported yet")') what; call goErr
        TRACEBACK; stop
    end select

    ! no error ...
    date_Get_End_Of%status = 0

  end function date_Get_End_Of


  ! ************************************************
  ! ***
  ! *** operators
  ! ***
  ! ************************************************


  subroutine date_Normalize( t )

    ! --- in/out --------------------------------

    type(TDate), intent(inout)       ::  t

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/date_Normalize'

    ! --- begin ---------------------------------

    ! mili seconds
    do
      if ( t%mili >= 0 ) exit
      t%sec = t%sec - 1
      t%mili = t%mili + 1000
    end do
    do
      if ( t%mili <= 999 ) exit
      t%mili = t%mili - 1000
      t%sec = t%sec + 1
    end do

    ! seconds
    do
      if ( t%sec >= 0 ) exit
      t%min = t%min - 1
      t%sec = t%sec + 60
    end do
    do
      if ( t%sec <= 59 ) exit
      t%sec = t%sec - 60
      t%min = t%min + 1
    end do

    ! minutes
    do
      if ( t%min >= 0 ) exit
      t%hour = t%hour - 1
      t%min = t%min + 60
    end do
    do
      if ( t%min <= 59 ) exit
      t%min = t%min - 60
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
        write (gol,'("unsupported calendar type: ",a)') t%calendar; call goErr
        TRACEBACK; stop
    end select

  end subroutine date_Normalize


  subroutine incrdate_Normalize( dt )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(inout)    :: dt

    ! --- const -----------------------------

    character(len=*), parameter  ::  rname = mname//'/incrdate_Normalize'

    ! --- begin ---------------------------------

    ! mili seconds
    do
      if ( dt%mili >= 0 ) exit
      dt%sec = dt%sec - 1
      dt%mili = dt%mili + 1000
    end do
    do
      if ( dt%mili <= 999 ) exit
      dt%mili = dt%mili - 1000
      dt%sec = dt%sec + 1
    end do

    ! seconds
    do
      if ( dt%sec >= 0 ) exit
      dt%min = dt%min - 1
      dt%sec = dt%sec + 60
    end do
    do
      if ( dt%sec <= 59 ) exit
      dt%sec = dt%sec - 60
      dt%min = dt%min + 1
    end do

    ! minutes
    do
      if ( dt%min >= 0 ) exit
      dt%hour = dt%hour - 1
      dt%min = dt%min + 60
    end do
    do
      if ( dt%min <= 59 ) exit
      dt%min = dt%min - 60
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

  type(TDate) function t_plus_t( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_plus_t'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t1, status )
    IF_NOTOK_STOP
    call Check( t2, status )
    IF_NOTOK_STOP

    ! any date ? return any date ..
    if ( (t1%calendar == 'none') .or. (t2%calendar == 'none') ) then
      t_plus_t = AnyDate()
      return
    end if

    ! calendars should be the same:
    if ( t1%calendar /= t2%calendar ) then
      write (gol,'("calendars should be the same : ")'); call goPr
      write (gol,'("  t1 : ",a)') trim(t1%calendar); call goPr
      write (gol,'("  t2 : ",a)') trim(t2%calendar); call goPr
      TRACEBACK; stop
    end if

    ! add all fields;
    t_plus_t = NewDate( calendar=t1%calendar, &
                        year  = t1%year  + t2%year  , &
                        month = t1%month + t2%month , &
                        day   = t1%day   + t2%day   , &
                        hour  = t1%hour  + t2%hour  , &
                        zone  = t1%zone  + t2%zone  , &
                        min   = t1%min   + t2%min   , &
                        sec   = t1%sec   + t2%sec   , &
                        mili  = t1%mili  + t2%mili     )

  end function t_plus_t


  ! *


  type(TDate) function t_plus_dt( t, dt )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t
    type(TIncrDate), intent(in)   ::  dt

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_plus_dt'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t, status )
    IF_NOTOK_STOP
    call Check( dt, status )
    IF_NOTOK_STOP

    ! any date ? return any date ..
    if ( t%calendar == 'none' ) then
      t_plus_dt = AnyDate()
      return
    end if

    ! add fields; normalization is applied in routine:
    t_plus_dt = NewDate( calendar = t%calendar, &
                         year     = t%year             , &
                         month    = t%month            , &
                         day      = t%day   + dt%day   , &
                         hour     = t%hour  + dt%hour  , &
                         zone     = t%zone             , &
                         min      = t%min   + dt%min   , &
                         sec      = t%sec   + dt%sec   , &
                         mili     = t%mili  + dt%mili     )

  end function t_plus_dt


  ! *


  type(TIncrDate) function dt_plus_dt( dt1, dt2 )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt1
    type(TIncrDate), intent(in)   ::  dt2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_plus_dt'

    ! --- local  --------------------------------

    integer                ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( dt1, status )
    IF_NOTOK_STOP
    call Check( dt2, status )
    IF_NOTOK_STOP

    ! add fields:
    dt_plus_dt = IncrDate( day  = dt1%day  + dt2%day   , &
                           hour = dt1%hour + dt2%hour  , &
                           min  = dt1%min  + dt2%min   , &
                           sec  = dt1%sec  + dt2%sec   , &
                           mili = dt1%mili + dt2%mili     )

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


  type(TIncrDate) function t_min_t( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

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

  type(TIncrDate) function t_min_t_work( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_min_t_work'

    ! --- local ---------------------------------

    integer             ::  status
    integer             ::  ndays
    type(TDate)         ::  t

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t1, status )
    IF_NOTOK_STOP
    call Check( t2, status )
    IF_NOTOK_STOP

    ! any dates ? something wrong ...
    if ( (t1%calendar == 'none') .or. (t2%calendar == 'none') ) then
      write (gol,'("do not know how to compute difference between `any` dates ...")')
      TRACEBACK; stop
    end if

    ! calendars should be the same:
    if ( t1%calendar /= t2%calendar ) then
      write (gol,'("calendars should be the same : ")'); call goPr
      write (gol,'("  t1 : ",a)') trim(t1%calendar); call goPr
      write (gol,'("  t2 : ",a)') trim(t2%calendar); call goPr
      TRACEBACK; stop
    end if

    ! difference between two dates;
    ! algorithm implemented for positive differences only:

    ! difference should be positive:
    if ( t1 < t2 ) then
      write (gol,'("expect t1 to exceed t2 :")'); call goErr
      call wrtgol( '  t1 : ', t1 ); call goErr
      call wrtgol( '  t2 : ', t2 ); call goErr
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
    t_min_t_work = IncrDate( day  = ndays, &
                             hour = t1%hour - t2%hour, &
                             min  = t1%min  - t2%min , &
                             sec  = t1%sec  - t2%sec , &
                             mili = t1%mili - t2%mili           )

  end function t_min_t_work


  ! *


  type(TDate) function t_min_dt( t, dt )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t
    type(TIncrDate), intent(in)   ::  dt

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/t_min_dt'

    ! --- local ---------------------------------

    integer             ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( t, status )
    IF_NOTOK_STOP
    call Check( dt, status )
    IF_NOTOK_STOP

    ! any date ? return any date ..
    if ( t%calendar == 'none' ) then
      t_min_dt = AnyDate()
      return
    end if

    ! result is of same type as t
    ! ~ first guess:
    t_min_dt = NewDate( calendar = t%calendar      , &
                        year     = t%year          , &
                        month    = t%month         , &
                        day      = t%day  -dt%day  , &
                        hour     = t%hour -dt%hour , &
                        zone     = t%zone          , &
                        min      = t%min  -dt%min  , &
                        sec      = t%sec  -dt%sec  , &
                        mili     = t%mili -dt%mili         )
    ! ~ normalize from negative or too large values into regular values:
    call Normalize( t_min_dt )

  end function t_min_dt


  ! *


  type(TIncrDate) function dt_min_dt( dt1, dt2 )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)       ::  dt1
    type(TIncrDate), intent(in)       ::  dt2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_min_dt'

    ! --- local ---------------------------------

    integer             ::  status

    ! --- begin ---------------------------------

    ! check arguments
    call Check( dt1, status )
    IF_NOTOK_STOP
    call Check( dt2, status )
    IF_NOTOK_STOP

    ! fill result:
    dt_min_dt = IncrDate( day  = dt1%day  - dt2%day , &
                          hour = dt1%hour - dt2%hour, &
                          min  = dt1%min  - dt2%min , &
                          sec  = dt1%sec  - dt2%sec , &
                          mili = dt1%mili - dt2%mili    )

  end function dt_min_dt


  ! *** date = t * r ************************************************


  ! multiply time with a real factor;
  ! use round for fractions

  type(TIncrDate) function dt_times_r4( dt, r )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
    real(4), intent(in)           ::  r

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_times_r4'

    ! --- local -----------------------------------

    real          ::  rday, rhour, rmin, rsec, rmili
    integer       ::  status

    ! --- begin ---------------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

    ! multiply each part, at remainder of previous:
    rday  = dt%day  * r
    rhour = dt%hour * r + (rday -nint(rday ))*24
    rmin  = dt%min  * r + (rhour-nint(rhour))*60
    rsec  = dt%sec  * r + (rmin -nint(rmin ))*60
    rmili = dt%mili * r + (rsec -nint(rsec ))*1000

    ! fill:
    dt_times_r4 = IncrDate( day  = nint(rday ), &
                            hour = nint(rhour), &
                            min  = nint(rmin ), &
                            sec  = nint(rsec ), &
                            mili = nint(rmili)      )

  end function dt_times_r4


  ! *


  type(TIncrDate) function r4_times_dt( r, dt )

    ! --- in/out --------------------------------

    real(4), intent(in)           ::  r
    type(TIncrDate), intent(in)   ::  dt

    ! --- begin ---------------------------------

    r4_times_dt = dt * r

  end function r4_times_dt
  
  
  ! *
  

  type(TIncrDate) function dt_times_r8( dt, r )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
    real(8), intent(in)           ::  r

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/dt_times_r8'

    ! --- local -----------------------------------

    real          ::  rday, rhour, rmin, rsec, rmili
    integer       ::  status

    ! --- begin ---------------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

    ! multiply each part, at remainder of previous:
    rday  = dt%day  * r
    rhour = dt%hour * r + (rday -nint(rday ))*24
    rmin  = dt%min  * r + (rhour-nint(rhour))*60
    rsec  = dt%sec  * r + (rmin -nint(rmin ))*60
    rmili = dt%mili * r + (rsec -nint(rsec ))*1000

    ! fill:
    dt_times_r8 = IncrDate( day  = nint(rday ), &
                            hour = nint(rhour), &
                            min  = nint(rmin ), &
                            sec  = nint(rsec ), &
                            mili = nint(rmili)      )

  end function dt_times_r8


  ! *


  type(TIncrDate) function r8_times_dt( r, dt )

    ! --- in/out --------------------------------

    real(8), intent(in)           ::  r
    type(TIncrDate), intent(in)   ::  dt

    ! --- begin ---------------------------------

    r8_times_dt = dt * r

  end function r8_times_dt


  ! *


  type(TIncrDate) function dt_times_i( dt, i )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
    integer, intent(in)           ::  i

    ! --- begin ---------------------------------

    dt_times_i = dt * (i*1.0)

  end function dt_times_i


  ! *


  type(TIncrDate) function i_times_dt( i, dt )

    ! --- in/out --------------------------------

    integer, intent(in)           ::  i
    type(TIncrDate), intent(in)   ::  dt

    ! --- begin ---------------------------------

    i_times_dt = dt * i

  end function i_times_dt


  ! *** dt = dt / r ************************************************


  type(TIncrDate) function dt_div_r( dt, r )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
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
    IF_NOTOK_STOP

    ! days:
    rat = dt%day / r
    intg = floor( rat )
    frac = rat - intg
    dt_div_r = IncrDate( day=intg )

    ! hours:
    rat = dt%hour / r + frac*24
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, hour=intg )

    ! mins:
    rat = dt%min / r + frac*60
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, min=intg )

    ! seconds:
    rat = dt%sec / r + frac*60
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, sec=intg )

    ! miliseconds:
    rat = dt%mili / r + frac*1000
    intg = floor( rat )
    frac = rat - intg
    call Set( dt_div_r, mili=intg )

  end function dt_div_r

  ! *

  type(TIncrDate) function dt_div_i( dt, i )

    ! --- in/out --------------------------------

    type(TIncrDate), intent(in)   ::  dt
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

    type(TDate), intent(in)        ::  t

    ! --- begin --------------------------------

    date_IsAnyDate = t%calendar == 'none'

  end function date_IsAnyDate


  ! ***  date1 == date2

  logical function date_eq_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_eq_date'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    call Check( t1, status )
    IF_NOTOK_STOP
    call Check( t2, status )
    IF_NOTOK_STOP

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
        ( t1%min   == t2%min   ) .and. &
        ( t1%sec   == t2%sec   ) .and. &
        ( t1%mili  == t2%mili  )

  end function date_eq_date


  ! ***  date1 /= date2

  logical function date_ne_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_ne_date'

    ! --- begin ---------------------------------

    date_ne_date = .not. ( t1 == t2 )

  end function date_ne_date


  ! ***  date1 > date2


  logical function date_gt_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_gt_date'

    ! --- local -----------------------------------

    integer          ::  status

    ! --- begin ---------------------------------

    call Check( t1, status )
    IF_NOTOK_STOP
    call Check( t2, status )
    IF_NOTOK_STOP

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

    if ( t1%min > t2%min ) then
      date_gt_date = .true.
      return
    else if ( t1%min < t2%min ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%sec > t2%sec ) then
      date_gt_date = .true.
      return
    else if ( t1%sec < t2%sec ) then
      date_gt_date = .false.
      return
    end if

    if ( t1%mili > t2%mili ) then
      date_gt_date = .true.
      return
    else if ( t1%mili < t2%mili ) then
      date_gt_date = .false.
      return
    end if

    ! all fields are equal ...
    date_gt_date = .false.

  end function date_gt_date


  ! ***  date1 < date2


  logical function date_lt_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- begin ---------------------------------

    date_lt_date = (.not.( ( t1 == t2 ) .or. ( t1 > t2 ) )) .or. IsAnyDate(t1) .or. IsAnyDate(t2)

  end function date_lt_date


  ! ***  date1 >= date2 ************************


  logical function date_ge_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

    ! --- begin ---------------------------------

    date_ge_date = ( t1 == t2 ) .or. ( t1 > t2 ) .or. IsAnyDate(t1) .or. IsAnyDate(t2)

  end function date_ge_date


  ! ***  date1 <= date2 ************************


  logical function date_le_date( t1, t2 )

    ! --- in/out --------------------------------

    type(TDate), intent(in)       ::  t1
    type(TDate), intent(in)       ::  t2

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

    type(TDate)                   ::  date_min
    type(TDate), intent(in)       ::  t1, t2

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

    type(TDate)                   ::  date_max
    type(TDate), intent(in)       ::  t1, t2

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

    type(TDate), intent(in)          ::  t
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_rTotal'

    ! --- local -----------------------------------

    integer            ::  status
    real               ::  nday
    integer            ::  iyear

    ! --- begin -----------------------------

    call Check( t, status )
    IF_NOTOK_STOP

    ! not all arguments are possible ...
    select case ( t%calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day' )
        select case ( unit )
          case ( 'year' )
            if ( any( (/t%month,t%day,t%hour,t%min,t%sec,t%mili/) /= 0 ) ) then
              write (gol,'("do not know how to count total:")'); call goErr
              write (gol,'("  unit : ",a)') unit; call goErr
              call wrtgol( '  t    : ', t ); call goErr
              TRACEBACK; stop
            end if
          case ( 'month' )
            if ( any( (/t%day,t%hour,t%min,t%sec,t%mili/) /= 0 ) ) then
              write (gol,'("do not know how to count total:")'); call goErr
              write (gol,'("  unit : ",a)') unit; call goErr
              call wrtgol( '  t    : ', t ); call goErr
              TRACEBACK; stop
            end if
        end select
      case ( 'incr' )
        select case ( unit )
          case ( 'year', 'month' )
            write (gol,'("do not know how to count total in incremental date:")') unit; call goErr
            write (gol,'("  unit : ",a)') unit; call goErr
            call wrtgol( '  t    : ', t ); call goErr
            TRACEBACK; stop
        end select
    end select

    ! precount total number of days for some of the units:
    select case ( unit )
      case ( 'day', 'days', 'hour', 'hours', 'min', 'minutes', 'sec', 'seconds', 'mili' )
        nday = 0.0
        do iyear = 1, t%year-1
          nday  = nday  + calc_days_in_year(t%calendar,iyear)
        end do
        nday  = nday + DayNumber( t ) - 1
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
                      t%min   / nday / 24.0 / 60.0                     + &
                      t%sec   / nday / 24.0 / 60.0 / 60.0              + &
                      t%mili  / nday / 24.0 / 60.0 / 60.0 / 1000.0
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
                      t%min   / nday / 24.0 / 60.0                     + &
                      t%sec   / nday / 24.0 / 60.0 / 60.0              + &
                      t%mili  / nday / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'day', 'days' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        date_rTotal = nday                                             + &
                      t%hour         / 24.0                            + &
                      t%min          / 24.0 / 60.0                     + &
                      t%sec          / 24.0 / 60.0 / 60.0              + &
                      t%mili         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour', 'hours' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        date_rTotal = nday           * 24.0                            + &
                      t%hour                                           + &
                      t%min                 / 60.0                     + &
                      t%sec                 / 60.0 / 60.0              + &
                      t%mili                / 60.0 / 60.0 / 1000.0
      case ( 'min', 'minutes' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        date_rTotal = nday           * 24.0 * 60.0                     + &
                      t%hour                * 60.0                     + &
                      t%min                                            + &
                      t%sec                        / 60.0              + &
                      t%mili                       / 60.0 / 1000.0
      case ( 'sec', 'seconds' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        date_rTotal = nday           * 24.0 * 60.0 * 60.0              + &
                      t%hour                * 60.0 * 60.0              + &
                      t%min                        * 60.0              + &
                      t%sec                                            + &
                      t%mili                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        date_rTotal = nday           * 24.0 * 60.0 * 60.0 * 1000.0     + &
                      t%hour                * 60.0 * 60.0 * 1000.0     + &
                      t%min                        * 60.0 * 1000.0     + &
                      t%sec                               * 1000.0     + &
                      t%mili
      case default
        write (gol,'("do not know how to count time in unit : ",a)') trim(unit); call goErr
        TRACEBACK; stop
    end select

  end function date_rTotal


  ! ***


  real function incr_rTotal( dt, unit )

    ! --- in/out ----------------------------

    type(TIncrDate), intent(in)      ::  dt
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incr_rTotal'

    ! --- local -----------------------------------

    integer            ::  status

    ! --- begin -----------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

    ! count time units:
    select case ( unit )
      case ( 'day', 'days' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        incr_rTotal = dt%day                                            + &
                      dt%hour         / 24.0                            + &
                      dt%min          / 24.0 / 60.0                     + &
                      dt%sec          / 24.0 / 60.0 / 60.0              + &
                      dt%mili         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour', 'hours' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        incr_rTotal = dt%day          * 24.0                            + &
                      dt%hour                                           + &
                      dt%min                 / 60.0                     + &
                      dt%sec                 / 60.0 / 60.0              + &
                      dt%mili                / 60.0 / 60.0 / 1000.0
      case ( 'min', 'minutes' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        incr_rTotal = dt%day          * 24.0 * 60.0                     + &
                      dt%hour                * 60.0                     + &
                      dt%min                                            + &
                      dt%sec                        / 60.0              + &
                      dt%mili                       / 60.0 / 1000.0
      case ( 'sec', 'seconds' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        incr_rTotal = dt%day          * 24.0 * 60.0 * 60.0              + &
                      dt%hour                * 60.0 * 60.0              + &
                      dt%min                        * 60.0              + &
                      dt%sec                                            + &
                      dt%mili                              / 1000.0
      case ( 'mili', 'miliseconds' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        incr_rTotal = dt%day          * 24.0 * 60.0 * 60.0 * 1000.0     + &
                      dt%hour                * 60.0 * 60.0 * 1000.0     + &
                      dt%min                        * 60.0 * 1000.0     + &
                      dt%sec                               * 1000.0     + &
                      dt%mili
      case default
        write (gol,'("do not know how to count time in unit : ",a)') trim(unit); call goErr
        TRACEBACK; stop
    end select

  end function incr_rTotal


  ! ***


  integer function date_iTotal( t, unit )

    ! --- in/out ----------------------------

    type(TDate), intent(in)          ::  t
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_iTotal'

    ! --- local -----------------------------

    integer            ::  status
    real               ::  rtot
    integer            ::  itot

    ! --- begin -----------------------------

    call Check( t, status )
    IF_NOTOK_STOP

    ! determine total some as a real value:
    rtot = rTotal( t, unit )

    ! round to integer value:
    itot = nint(rtot)

    ! result should be pure integer ....
    if ( itot*1.0 == rtot ) then
      date_iTotal = itot
    else
      write (gol,'("date does not contain integer total:")'); call goErr
      write (gol,'("  unit : ",a)') trim(unit); call goErr
      call wrtgol( '  t    : ', t ); call goErr
      TRACEBACK; stop
    end if

  end function date_iTotal


  ! ***


  real(8) function date_dTotal( t, unit )

    ! --- in/out ----------------------------

    type(TDate), intent(in)          ::  t
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/date_dTotal'

    ! --- local -----------------------------------

    integer            ::  status
    real(8)            ::  nday
    integer            ::  iyear

    ! --- begin -----------------------------

    call Check( t, status )
    IF_NOTOK_STOP

    ! not all arguments are possible ...
    select case ( t%calendar )
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day' )
        select case ( unit )
          case ( 'year' )
            if ( any( (/t%month,t%day,t%hour,t%min,t%sec,t%mili/) /= 0 ) ) then
              write (gol,'("do not know how to count total:")'); call goErr
              write (gol,'("  unit : ",a)') unit; call goErr
              call wrtgol( '  t    : ', t ); call goErr
              TRACEBACK; stop
            end if
          case ( 'month' )
            if ( any( (/t%day,t%hour,t%min,t%sec,t%mili/) /= 0 ) ) then
              write (gol,'("do not know how to count total:")'); call goErr
              write (gol,'("  unit : ",a)') unit; call goErr
              call wrtgol( '  t    : ', t ); call goErr
              TRACEBACK; stop
            end if
        end select
      case ( 'incr' )
        select case ( unit )
          case ( 'year', 'month' )
            write (gol,'("do not know how to count total in incremental date:")') unit; call goErr
            write (gol,'("  unit : ",a)') unit; call goErr
            call wrtgol( '  t    : ', t ); call goErr
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
                      t%min   / nday / 24.0 / 60.0                     + &
                      t%sec   / nday / 24.0 / 60.0 / 60.0              + &
                      t%mili  / nday / 24.0 / 60.0 / 60.0 / 1000.0
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
                      t%min   / nday / 24.0 / 60.0                     + &
                      t%sec   / nday / 24.0 / 60.0 / 60.0              + &
                      t%mili  / nday / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        date_dTotal = nday                                             + &
                      t%hour         / 24.0                            + &
                      t%min          / 24.0 / 60.0                     + &
                      t%sec          / 24.0 / 60.0 / 60.0              + &
                      t%mili         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        date_dTotal = nday           * 24.0                            + &
                      t%hour                                           + &
                      t%min                 / 60.0                     + &
                      t%sec                 / 60.0 / 60.0              + &
                      t%mili                / 60.0 / 60.0 / 1000.0
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        date_dTotal = nday           * 24.0 * 60.0                     + &
                      t%hour                * 60.0                     + &
                      t%min                                            + &
                      t%sec                        / 60.0              + &
                      t%mili                       / 60.0 / 1000.0
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        date_dTotal = nday           * 24.0 * 60.0 * 60.0              + &
                      t%hour                * 60.0 * 60.0              + &
                      t%min                        * 60.0              + &
                      t%sec                                            + &
                      t%mili                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        date_dTotal = nday           * 24.0 * 60.0 * 60.0 * 1000.0     + &
                      t%hour                * 60.0 * 60.0 * 1000.0     + &
                      t%min                        * 60.0 * 1000.0     + &
                      t%sec                               * 1000.0     + &
                      t%mili
      case default
        write (gol,'("do not know how to count time in unit : ",a)') trim(unit); call goErr
        TRACEBACK; stop
    end select

  end function date_dTotal


  ! ***


  real(8) function incr_dTotal( dt, unit )

    ! --- in/out ----------------------------

    type(TIncrDate), intent(in)      ::  dt
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incr_rTotal'

    ! --- local -----------------------------------

    integer            ::  status

    ! --- begin -----------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

    ! count time units:
    select case ( unit )
      case ( 'day' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional months:
        incr_dTotal = dt%day                                            + &
                      dt%hour         / 24.0                            + &
                      dt%min          / 24.0 / 60.0                     + &
                      dt%sec          / 24.0 / 60.0 / 60.0              + &
                      dt%mili         / 24.0 / 60.0 / 60.0 / 1000.0
      case ( 'hour' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional hours:
        incr_dTotal = dt%day          * 24.0                            + &
                      dt%hour                                           + &
                      dt%min                 / 60.0                     + &
                      dt%sec                 / 60.0 / 60.0              + &
                      dt%mili                / 60.0 / 60.0 / 1000.0
      case ( 'min' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional minutes:
        incr_dTotal = dt%day          * 24.0 * 60.0                     + &
                      dt%hour                * 60.0                     + &
                      dt%min                                            + &
                      dt%sec                        / 60.0              + &
                      dt%mili                       / 60.0 / 1000.0
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        incr_dTotal = dt%day          * 24.0 * 60.0 * 60.0              + &
                      dt%hour                * 60.0 * 60.0              + &
                      dt%min                        * 60.0              + &
                      dt%sec                                            + &
                      dt%mili                              / 1000.0
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        incr_dTotal = dt%day          * 24.0 * 60.0 * 60.0 * 1000.0     + &
                      dt%hour                * 60.0 * 60.0 * 1000.0     + &
                      dt%min                        * 60.0 * 1000.0     + &
                      dt%sec                               * 1000.0     + &
                      dt%mili
      case default
        write (gol,'("do not know how to count time in unit : ",a)') trim(unit); call goErr
        TRACEBACK; stop
    end select

  end function incr_dTotal

  ! ***  
  
  integer function incr_iTotal( dt, unit )

    ! --- in/out ----------------------------

    type(TIncrDate), intent(in)      ::  dt
    character(len=*), intent(in)     ::  unit

    ! --- const -----------------------------------

    character(len=*), parameter  ::  rname = mname//'/incr_iTotal'

    ! --- local -----------------------------------

    integer            ::  status

    ! --- begin -----------------------------

    call Check( dt, status )
    IF_NOTOK_STOP

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
                      dt%min
      case ( 'sec' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional seconds:
        incr_iTotal = dt%day          * 24 * 60 * 60            + &
                      dt%hour              * 60 * 60            + &
                      dt%min                    * 60            + &
                      dt%sec
      case ( 'mili' )
        ! 'nday' has been set to the total number of days from 0 to t;
        ! count fractional mili seconds:
        incr_iTotal = dt%day          * 24 * 60 * 60 * 1000     + &
                      dt%hour              * 60 * 60 * 1000     + &
                      dt%min                    * 60 * 1000     + &
                      dt%sec                         * 1000     + &
                      dt%mili
      case default
        write (gol,'("do not know how to count time in unit : ",a)') trim(unit); call goErr
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

    type(TDate), intent(in)    ::  t
    type(TDate), intent(in)    ::  t1
    type(TDate), intent(in)    ::  t2
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
      write (gol,'("required interval [t1,t2] :")'); call goErr
      call wrtgol( '  t1 = ', t1 ); call goErr
      call wrtgol( '  t2 = ', t2 ); call goErr
      TRACEBACK; status=1; return
    end if

    ! check ...
    if ( (t < t1) .or. (t > t2) ) then
      write (gol,'("t not in [t1,t2] :")'); call goErr
      call wrtgol( '  t  = ', t  ); call goErr
      call wrtgol( '  t1 = ', t1 ); call goErr
      call wrtgol( '  t2 = ', t2 ); call goErr
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
  
    use GO_String, only : goReadFromLine, goSplitString

    ! --- in/out -----------------------------

    type(TDate), intent(in)       ::  t
    type(TDate), intent(in)       ::  tref
    character(len=*), intent(in)  ::  interpolation
    type(TDate), intent(out)      ::  tt(2)
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
    type(TDate)         ::  t0

    ! --- begin ------------------------------
    
    if ( index(trim(interpolation), ';') > 0 ) then
      
      ! different units over time period
      ! split:
      call goSplitString( interpolation, ninterp, interp, status, sep=';' )
      IF_NOTOK_RETURN(status=1)
      
      if ( ninterp > 2 ) then
        write( gol, '("More than 2 time inervals found: not supported")' ) ; call GoErr
        TRACEBACK;status=1;return
      end if
      ! Read second part:
      ! format: * hour after *
      ! split at after
      call goReadFromLine( interp(2), step, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( interp(2), units, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( interp(2), after, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( interp(2), tchange, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      offset = 0
      
      ! Check if time is before offset
      if ( itotal( t - tref, 'hour' ) < tchange ) then
        ! read first inerpolation part
        call goReadFromLine( interp(1), step, status, sep=' ' )
        IF_NOTOK_RETURN(status=1)
        call goReadFromLine( interp(1), units, status, sep=' ' )
        IF_NOTOK_RETURN(status=1)
        offset = 0
      endif
    else ! read first inerpolation part
      line = trim(interpolation)
      call goReadFromLine( line, step, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      call goReadFromLine( line, units, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
      offset = 0
    endif
      
    ! which units ?
    select case ( units )

      ! x hourly
      case ( 'hour' )
        ! begin of day:
        t0 = Get_Begin_Of( t, 'day' )
        ! add offset:
        t0 = t0 + IncrDate(hour=offset)
        ! back one day?
        if ( t0%hour > t%hour ) t0 = t0 - IncrDate(hour=24)
        ! hour fraction:
        r = rTotal( t - t0, 'hour' )
        ! start of interval:
        tt(1) = t0 + IncrDate( hour=int(floor(r/step))*step )
        ! end:
        if ( t == tt(1) ) then
          tt(2) = tt(1)
        else
          tt(2) = tt(1) + IncrDate(hour=step)
        end if

      ! unknown ...
      case default
        write (gol,'("unsupported units `",a,"`")') trim(units); call goErr
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

  subroutine goReadFromLine_t( input, t, status, sep, calendar )

    use GO_String, only : goReadFromLine
    use GO_String, only : goTranslate

    ! --- in/out ---------------------------------

    character(len=*), intent(inout)           ::  input
    type(TDate), intent(out)                  ::  t
    integer, intent(out)                      ::  status
    character(len=1), intent(in), optional    ::  sep
    character(len=*), intent(in), optional    ::  calendar

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/goReadFromLine_t'

    ! --- local ----------------------------------

    character(len=32)     ::  line
    integer               ::  year, month, day, hour, minu
    real                  ::  seco

    ! --- begin ----------------------------------

    ! extract:
    call goReadFromLine( input, line, status, sep=sep )
    IF_NOTOK_RETURN(status=1)
    
    ! replace seperators:
    call goTranslate( line, '-/:_TZ', ' ', status )
    IF_NOTOK_RETURN(status=1)

    ! extract ref time values:
    call goReadFromLine( line, year, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( line, month, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    call goReadFromLine( line, day, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    if ( len_trim(line) > 0 ) then
      call goReadFromLine( line, hour, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
    else
      hour = 0
    end if
    if ( len_trim(line) > 0 ) then
      call goReadFromLine( line, minu, status, sep=' ' )
      IF_NOTOK_RETURN(status=1)
    else
      minu = 0
    end if
    if ( len_trim(line) > 0 ) then
      call goReadFromLine( line, seco, status )
      IF_NOTOK_RETURN(status=1)
    else
      seco = 0.0
    end if

    ! store:
    t = NewDate( year=year, month=month, day=day, &
                    hour=hour, min=minu, sec=nint(seco), &
                    calendar=calendar )

    ! ok
    status = 0

  end subroutine goReadFromLine_t
  
  ! *

  ! extract step and reftime from units:
  !   hours since 1900-01-01 00:00:0.0

  subroutine Extract_Ref_and_Step( units, tref, tstep, status, calendar, stepunits )

    use GO_String, only : goReadFromLine

    ! --- in/out ---------------------------------

    character(len=*), intent(in)    ::  units
    type(TDate), intent(out)        ::  tref
    type(TIncrDate), intent(out)    ::  tstep
    integer, intent(out)            ::  status
    character(len=*), intent(in), optional    ::  calendar
    character(len=*), intent(out), optional   ::  stepunits

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Extract_Ref_and_Step'

    ! --- local ----------------------------------

    character(len=512)    ::  line
    character(len=32)     ::  steps
    character(len=32)     ::  dummy
    integer               ::  year, month, day, hour, minu
    real                  ::  seco

    ! --- begin ----------------------------------

    ! copy:
    line = trim(units)

    ! extract step units:
    call goReadFromLine( line, steps, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! convert ...
    select case ( trim(steps) )
      case ( 'days' )
        tstep = IncrDate(day=1)
      case ( 'hours' )
        tstep = IncrDate(hour=1)
      case ( 'seconds' )
        tstep = IncrDate(sec=1)
      case default
        write (gol,'("unsupported steps `",a,"`")') trim(steps); call goPr
        TRACEBACK; status=1; return
    end select
    
    ! return argument?
    if ( present(stepunits) ) stepunits = trim(steps)

    ! extract 'since'
    call goReadFromLine( line, dummy, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    if ( trim(dummy) /= 'since' ) then
      write (gol,'("second field should be `since`, not `",a,"`")') trim(dummy); call goPr
      TRACEBACK; status=1; return
    end if

    ! extract ref time values:
    call goReadFromLine( line, tref, status, calendar=calendar )
    IF_NOTOK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine Extract_Ref_and_Step
  
  ! *

  !
  ! Expand time value using units:
  !
  !   hours since 1900-01-01 00:00:0.0
  !
  ! Also step 'months' are accepted.
  !

  subroutine ExpandTime( value, units, calendar, t, status )

    use GO_String, only : goReadFromLine

    ! --- in/out ---------------------------------

    real, intent(in)                ::  value
    character(len=*), intent(in)    ::  units
    character(len=*), intent(in)    ::  calendar
    type(TDate), intent(out)        ::  t
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/ExpandTime'

    ! --- local ----------------------------------

    character(len=512)    ::  line
    character(len=32)     ::  steps
    character(len=32)     ::  dummy
    type(TDate)           ::  tref
    integer               ::  nyear, nmonth, nday
    real                  ::  rmonth

    ! --- begin ----------------------------------

    ! copy:
    line = trim(units)

    ! extract step units:
    call goReadFromLine( line, steps, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)

    ! extract 'since'
    call goReadFromLine( line, dummy, status, sep=' ' )
    IF_NOTOK_RETURN(status=1)
    ! check ...
    if ( trim(dummy) /= 'since' ) then
      write (gol,'("second field should be `since`, not `",a,"`")') trim(dummy); call goPr
      TRACEBACK; status=1; return
    end if

    ! extract ref time values:
    call goReadFromLine( line, tref, status, calendar=calendar )
    IF_NOTOK_RETURN(status=1)

    ! expand:
    select case ( trim(steps) )
    
      case ( 'month', 'months' )
        ! number of years and full months:
        nyear  = int(value/12)
        rmonth = value - nyear*12
        nmonth = int(rmonth)
        ! init result:
        t = NewDate( year=tref%year+nyear, month=tref%month+nmonth, day=1 )
        ! month:
        nday = Days_in_Month( t )
        ! remaining month fraction:
        rmonth = rmonth - nmonth
        ! add rest:
        t = t + IncrDate( hour=int(rmonth*nday*24) )

      case ( 'day', 'days' )
        t = tref + IncrDate(day=1) * value
      case ( 'hour', 'hours' )
        t = tref + IncrDate(hour=1) * value
      case ( 'minute', 'minutes' )
        t = tref + IncrDate(min=1) * value
      case ( 'second', 'seconds' )
        t = tref + IncrDate(sec=1) * value

      case default
        write (gol,'("unsupported steps `",a,"`")') trim(steps); call goPr
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

    type(TDate), intent(in)       ::  t
    real, intent(in)              ::  num
    character(len=*), intent(in)  ::  units
    character(len=*), intent(in)  ::  calendar
    integer, intent(out)          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Compare_Date_Num'

    ! --- local ----------------------------------

    type(TDate)           ::  t_rec

    ! --- begin ----------------------------------

    ! check ...
    if ( t%calendar /= calendar ) then
      write (gol,'("calendar of t is `",a,"` instead of `",a,"`")') trim(t%calendar), trim(calendar); call goErr
      TRACEBACK; status=1; return
    end if

    ! expand:
    call ExpandTime( num, units, calendar, t_rec, status )
    IF_NOTOK_RETURN(status=1)
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

    type(TDate), intent(in)       ::  t

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
      !    t%hour, t%min, t%sec, t%mili, &
      !    t%day, month_name(t%month), t%year, &
      !    zone_sign, zone_hour, zone_min
      case ( 'proleptic_gregorian', 'gregorian', 'standard', '366_day', '365_day', '360_day', 'none' )
        write (s,'(i4.4,"-",i2.2,"-",i2.2," ",i2,":",i2.2,":",i2.2)') &
          t%year, t%month, t%day, t%hour, t%min, t%sec
      case default
        s = 'no-calendar'
    end select

    date_Pretty = s

  end function date_Pretty


  ! *


  character(len=36) function incrdate_Pretty( dt )

    ! --- in/out -------------------------

    type(TIncrDate), intent(in)       ::  dt

    ! --- local --------------------------

    character(len=36)       ::  s

    ! --- begin --------------------------

    write (s,'(i5," days ",i2,":",i2.2,":",i2.2,":",i3.3)') &
      dt%day, dt%hour, dt%min, dt%sec, dt%mili

    incrdate_Pretty = s

  end function incrdate_Pretty


  ! *


  subroutine wrtgol_t( msg, t )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  t

    ! --- local ---------------------------------

    character(len=36)     ::  s

    ! --- begin -----------------------------------

    s = date_Pretty( t )
    write (gol,'(a,a)') msg, trim(s)

  end subroutine wrtgol_t


  ! *


  subroutine wrtgol_dt( msg, dt )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TIncrDate), intent(in)     ::  dt

    ! --- local ---------------------------------

    character(len=36)     ::  s

    ! --- begin -----------------------------------

    s = incrdate_Pretty( dt )
    write (gol,'(a,a)') msg, trim(s)

  end subroutine wrtgol_dt


  ! *


  subroutine wrtgol_t_dt( msg, t, msg2, dt )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  msg2
    type(TINcrDate), intent(in)     ::  dt

    ! --- local ---------------------------------

    character(len=36)     ::  s
    character(len=36)     ::  s2

    ! --- begin -----------------------------------

    s  = date_Pretty( t  )
    s2 = incrdate_Pretty( dt )
    write (gol,'(a,a,a,a)') msg, trim(s), msg2, trim(s2)

  end subroutine wrtgol_t_dt


  ! *


  subroutine wrtgol_tt( msg, tt )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  tt(2)

    ! --- local ---------------------------------

    character(len=36)     ::  s1
    character(len=36)     ::  s2

    ! --- begin -----------------------------------

    s1 = date_Pretty( tt(1) )
    s2 = date_Pretty( tt(2) )
    write (gol,'(a,"[",a,",",a,"]")') msg, trim(s1), trim(s2)

  end subroutine wrtgol_tt


  ! *


  subroutine wrtgol_t1_t2( msg, t, msg2, t2 )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  msg2
    type(TDate), intent(in)         ::  t2

    ! --- local ---------------------------------

    character(len=36)     ::  s
    character(len=36)     ::  s2

    ! --- begin -----------------------------------

    s  = date_Pretty( t  )
    s2 = date_Pretty( t2 )
    write (gol,'(a,a,a,a)') msg, trim(s), msg2, trim(s2)

  end subroutine wrtgol_t1_t2


  ! *


  subroutine wrtgol_t1_t2_t3( msg, t, msg2, t2, msg3, t3 )

    use go_print, only : gol

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  t
    character(len=*), intent(in)    ::  msg2
    type(TDate), intent(in)         ::  t2
    character(len=*), intent(in)    ::  msg3
    type(TDate), intent(in)         ::  t3

    ! --- local ---------------------------------

    character(len=36)     ::  s
    character(len=36)     ::  s2
    character(len=36)     ::  s3

    ! --- begin -----------------------------------

    s  = date_Pretty( t  )
    s2 = date_Pretty( t2 )
    s3 = date_Pretty( t3 )
    write (gol,'(a,a,a,a,a,a)') msg, trim(s), msg2, trim(s2), msg3, trim(s3)

  end subroutine wrtgol_t1_t2_t3


  ! *


  subroutine PrintDate( msg, t )

    ! --- in/out -----------------------------------

    character(len=*), intent(in)    ::  msg
    type(TDate), intent(in)         ::  t

    ! --- local ---------------------------------

    character(len=36)     ::  s

    ! --- begin -----------------------------------

    s = date_Pretty( t )
    write (*,'(a,a)') msg, trim(s)

  end subroutine PrintDate


end module GO_Date


