!###############################################################################
!
! go print : tools for standard output
!
! Example:
!
!   ! use logger instance from this module:
!   use GO_Print, only : lgr
!
!   ! messages printed by root only:
!   call lgr%Init( status, apply=myid==root, trace=.false. )
!   if (status/=0) stop
!
!   ! change destination of messages;
!   ! if no file name is provided, messages will be written to standard output:
!   call lgr%LogFile( status [,file='messages.log'] )
!   if (status/=0) stop
!
!   ! set routine label:
!   call goLabel( 'mymod/myroutine' )
!
!   ! write single message:
!   !   This is number  3
!   write (csol,'("This is number ",i2)') 3; call csoPr
!
!   ! write error message and traceback using the
!   ! previous defined routine label:
!   !   ERROR - Something wrong.
!   !   ERROR in mymod/myroutine
!   write (csol,'("Something wrong.")'); call csoErr
!   call csoErr
!
!   ! close label
!   call goLabel()
!
!   ! done
!   call lgr%Done( status )
!   if (status/=0) stop
!
!###############################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call csoErr
!
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!###############################################################################


module CSO_Logging

  implicit none

  ! --- in/out ---------------------------------

  private

  public  ::  csol
  public  ::  ncsol
  public  ::  csoPr, csoErr
!  public  ::  goBug
!  public  ::  goLabel

  public  ::  lgr


  ! --- const ---------------------------------

  character(len=*), parameter  ::  mname = 'CSO_Logging'

  ! standard file units
  integer, parameter    ::  csoStdErr = 0
  integer, parameter    ::  csoStdIn  = 5
  integer, parameter    ::  csoStdOut = 6


  ! buffer size:
  integer, parameter           ::  len_csol = 1024

!  ! max label stack height:
!  integer, parameter   ::  mstack = 400
!
!  ! white space for indents:
!  integer, parameter   ::  dindent = 2


  ! --- type -----------------------------------

  ! storage for logger:
  type TgoLogger

!    ! stack with labels:
!    character(len=64)    ::  labels(0:mstack)
!    integer              ::  istack = 0

    ! initialized ?
    ! some errors might be printed before initialization ...
    logical              ::  pr_initialized = .false.

    ! destination file unit, by default standard output:
    integer              ::  pr_fu = csoStdOut

    ! flags etc
    logical              ::  pr_apply
!    logical              ::  pr_trace
!
!    ! number of indents:
!    integer              ::  indent = 0

    ! writ to file ?
    logical              ::  pr_file
    character(len=256)   ::  pr_file_name
    logical              ::  pr_file_echo
  contains
    procedure   ::  Init      => Logger_Init
    procedure   ::  Done      => Logger_Done
    procedure   ::  Set       => Logger_Set
    procedure   ::  LogFile   => Logger_LogFile
!    procedure   ::  PrFU      => Logger_PrFU
!    procedure   ::  ErrFU     => Logger_ErrFU
  end type TgoLogger


  ! --- var ------------------------------------

  ! buffer for standard output
  character(len=len_csol)   ::  csol
  ! actual length, sometimes needed by error message routines:
  integer                   ::  ncsol

  ! current logger ;
  ! the 'save' attribute is needed by some compilers
  ! since 'pr_initialized' is given a value in the
  ! type definition above:
  type(TgoLogger), save     ::  lgr


contains



  ! ***************************************************************************
  ! ***
  ! *** printing
  ! ***
  ! ***************************************************************************


  subroutine csoPr

    ! --- local --------------------------------

!    integer             ::  nind

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/csoPr'

    ! --- begin --------------------------------

    ! not initialized yet ? then print to standard output:
    if ( .not. lgr%pr_initialized ) then
      write (*,'(a)') trim(csol)
      return
    end if

    ! print go line ?
    if ( lgr%pr_apply ) then

!      ! number of spaces to indent:
!      nind = max( 0, lgr%indent )
!
!      ! indention?
!      if ( nind > 0 ) then
!        ! write indention and line:
!        write (lgr%pr_fu,'(a,a)') repeat(' ',nind), trim(csol)
!      else
         ! line only:
         write (lgr%pr_fu,'(a)') trim(csol)
!      end if
!
!      ! show:
!      call flush(lgr%pr_fu)
!
!      ! to file and echo ?
!      if ( lgr%pr_file .and. lgr%pr_file_echo ) then
!        ! indention?
!        if ( nind > 0 ) then
!          ! write indention and line:
!          write (csoStdOut,'(a,a)') repeat(' ',nind), trim(csol)
!        else
!          ! line only:
!          write (csoStdOut,'(a)') trim(csol)
!        end if
!      end if

    end if

    ! clear output line:
    csol = ''

  end subroutine csoPr


  ! ***

  ! Print error message.
  ! Now printed to standard output, in future to standard error ?
  ! Make csol empty before leaving.
  ! If still empty in next call, this is a trace back
  !   (print error label, one label back)


  subroutine csoErr

    ! --- local -------------------------------

!    integer     ::  ilab

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/csoErr'

    ! --- local ----------------------------

    logical                    ::  save_pr_apply
    character(len=len(csol))   ::  csol2

    ! --- begin --------------------------------

    ! store original apply flag:
    save_pr_apply = lgr%pr_apply
    ! always print error messages:
    lgr%pr_apply = .true.

!    ! message in buffer ?
!    if ( len_trim(csol) > 0 ) then

      ! error message:
      csol2 = trim(csol)
      write (csol,'("ERROR - ",a)') trim(csol2); call csoPr

!    else
!
!      ! label index:
!      ilab = min( lgr%istack, mstack )
!
!      ! write error message:
!      write (csol,'("ERROR in ",a)') trim(lgr%labels(ilab)); call csoPr
!
!      ! one level back:
!      call goLabel()
!
!    end if

    ! restore apply flag:
    lgr%pr_apply = save_pr_apply

  end subroutine csoErr


!  ! ***
!
!
!  subroutine goBug
!
!    ! --- local ----------------------------
!
!    logical   ::  save_pr_apply
!
!    ! --- begin --------------------------------
!
!    ! store original apply flag:
!    save_pr_apply = lgr%pr_apply
!    ! always print bug messages:
!    lgr%pr_apply = .true.
!
!    ! write message
!    write (csol,'("BUG - ",a)') trim(csol); call csoPr
!
!    ! restore apply flag:
!    lgr%pr_apply = save_pr_apply
!
!  end subroutine goBug
!
!
!  ! ***************************************************************************
!  ! ***
!  ! *** routine labels
!  ! ***
!  ! ***************************************************************************
!
!
!  subroutine goLabel( label )
!
!    ! --- in/out -------------------------------
!
!    character(len=*), intent(in), optional  ::  label
!
!    ! --- const ----------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/goLabel'
!
!    ! --- begin --------------------------------
!
!    ! add new label to stack ?
!    if ( present(label) ) then
!      lgr%istack = lgr%istack + 1
!      if ( lgr%istack > mstack ) then
!        write (csol,'("BUG - stack too small; please increase mstack in go_print")'); call csoPr
!      else
!        lgr%labels(lgr%istack) = label
!      end if
!      if (lgr%pr_trace) then
!        write (csol,'("<",a,">")') trim(lgr%labels(lgr%istack)); call csoPr
!      end if
!      lgr%indent = lgr%indent + dindent
!    else
!      lgr%indent = lgr%indent - dindent
!      if (lgr%pr_trace) then
!        write (csol,'("(",a,")")') trim(lgr%labels(lgr%istack)); call csoPr
!      end if
!      lgr%istack = max( 0, lgr%istack - 1 )
!    end if
!
!  end subroutine goLabel



  ! =============================================================================
  ! ===
  ! === logger type
  ! ===
  ! =============================================================================


  subroutine Logger_Init( self, status )

    ! --- in/out ----------------------------

    class(TgoLogger), intent(out)            ::  self
    integer, intent(out)                     ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Logger_Init'

    ! --- local -----------------------------

    ! --- begin -----------------------------

    ! print or not ?
    self%pr_apply = .true.

!    ! trace labels ?
!    self%pr_trace = .false.

    ! do not write to file yet:
    self%pr_file      = .false.
    self%pr_file_name = ''
    self%pr_file_echo = .false.

    ! write to standard output:
    self%pr_fu = csoStdOut

    ! now the module is initialized ...
    self%pr_initialized = .true.

    ! ok
    status = 0

  end subroutine Logger_Init


  ! ***


  subroutine Logger_Done( self, status )

    ! --- in/out ----------------------------

    class(TgoLogger), intent(inout)   ::  self
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Logger_Done'

    ! --- begin -----------------------------

    ! close logfile if nesseary (switch to standard output)
    call self%Logfile( status )
    if (status/=0) then; write (*,'("ERROR - in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; status=1; return; end if

    ! ok
    status = 0

  end subroutine Logger_Done


  ! ***


!  subroutine Logger_Set( self, status, apply, trace )
  subroutine Logger_Set( self, status, apply, unit )

    ! --- in/out ----------------------------

    class(TgoLogger), intent(inout)          ::  self
    integer, intent(out)                     ::  status
    logical, intent(in), optional            ::  apply
    integer, intent(in), optional            ::  unit
!    logical, intent(in), optional            ::  trace

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Logger_Set'

    ! --- local -----------------------------

    ! --- begin -----------------------------

    ! print or not ?
    if ( present(apply) ) self%pr_apply = apply
    
    ! set file unit?
    if ( present(unit) ) self%pr_fu = unit

!    ! trace labels ?
!    if ( present(trace) ) self%pr_trace = trace

    ! ok
    status = 0

  end subroutine Logger_Set


  ! ***


!  subroutine Logger_Logfile( self, status, file, echo )
  subroutine Logger_Logfile( self, status, file )

    ! --- in/out ----------------------------

    class(TgoLogger), intent(inout)          ::  self
    integer, intent(out)                     ::  status
    character(len=*), intent(in), optional   ::  file
!    logical, intent(in), optional            ::  echo

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Logger_Logfile'

    ! --- local -----------------------------

    logical   ::  opened

    ! --- begin -----------------------------

    ! already a log file open ?
    if ( self%pr_file ) then
      ! close:
      close( self%pr_fu, iostat=status )
      if ( status/=0 ) then
        write (*,'("ERROR - closing output file:")')
        write (*,'("ERROR -   unit : ",i6)') self%pr_fu
        write (*,'("ERROR -   file : ",a)') trim(self%pr_file_name)
        write (*,'("ERROR - in ",a," (",a,i6,")")') rname, __FILE__, __LINE__; status=1; return
      end if
      ! do not write to files anymore:
      self%pr_file = .false.
      ! write to standard output from now on:
      self%pr_fu = csoStdOut
    end if

    ! open a new file ?
    if ( present(file) ) then

      ! store name of (new) logfile:
      self%pr_file_name = trim(file)

      ! info ...
      write (csol,'("switch to logfile : ",a)') trim(self%pr_file_name); call csoPr

      ! select free file unit:
      self%pr_fu = 789
      do
        inquire( self%pr_fu, opened=opened )
        if ( .not. opened ) exit
        self%pr_fu = self%pr_fu + 1
      end do
      ! open requested output file:
      open( unit=self%pr_fu, file=trim(self%pr_file_name), status='replace', iostat=status )
      if ( status/=0 ) then
        write (csol,'("opening file for output:")'); call csoPr
        write (csol,'("  unit : ",i6)') self%pr_fu; call csoPr
        write (csol,'("  file : ",a)') trim(self%pr_file_name); call csoPr
        TRACEBACK; status=1; return
      end if

      ! write to file from now on:
      self%pr_file = .true.

!      ! echo to std.out. too ?
!      if ( present(echo) ) then
!        self%pr_file_echo = echo
!      else
!        self%pr_file_echo = .false.
!      end if

    end if

    ! ok
    status = 0

  end subroutine Logger_Logfile


  ! ***


!  integer function Logger_PrFU( self )
!
!    ! --- in/out ----------------------------
!
!    class(TgoLogger), intent(in)          ::  self
!
!    ! --- begin -----------------------------
!
!    ! copy file unit:
!    Logger_PrFU = lgr%pr_fu
!
!  end function Logger_PrFU
!
!  ! *
!
!  integer function Logger_ErrFU( self )
!
!    ! --- in/out ----------------------------
!
!    class(TgoLogger), intent(in)          ::  self
!
!    ! --- begin -----------------------------
!
!    ! copy file unit:
!    Logger_ErrFU = lgr%pr_fu
!
!  end function Logger_ErrFU


end module CSO_Logging



! #############################################################################
! ###
! ### test program
! ###
! #############################################################################
!
!
!module testmod
!
!  implicit none
!
!  public
!
!contains
!
!  subroutine subr( i, status )
!
!    use go_print, only : goLabel, csol, csoPr, csoErr
!
!    ! --- in/out ----------------------------------------
!
!    integer, intent(in)           ::  i
!    integer, intent(out)          ::  status
!
!    ! --- begin -----------------------------------------
!
!    call goLabel( 'subr' )
!
!    write (csol,'("welcome to subr !")'); call csoPr
!
!    select case ( i )
!
!      case ( 0 )
!        write (csol,'("testing i : ",i2)') i; call csoPr
!
!      case ( 1 )
!        call subr2( 0, status )
!        if (status/=0) then; call csoErr; status=1; return; end if
!
!      case ( 2 )
!        call subr2( 1, status )
!        if (status/=0) then; call csoErr; status=1; return; end if
!
!      case default
!        write (csol,'("unsupported i : ",i2)') i; call csoErr
!        call csoErr; status=1; return
!
!    end select
!
!    call goLabel(); status=0
!
!  end subroutine subr
!
!
!  ! ***
!
!
!  subroutine subr2( i, status )
!
!    use go_print, only : goLabel, csol, csoPr, csoErr
!
!    ! --- in/out ----------------------------------------
!
!    integer, intent(in)           ::  i
!    integer, intent(out)          ::  status
!
!    ! --- begin -----------------------------------------
!
!    call goLabel('subr2')
!
!    write (csol,'("testing subr2")'); call csoPr
!
!    select case ( i )
!      case ( 0 )
!      case default
!        write (csol,'("wrong i : ",i2)') i; call csoErr
!        call csoErr; status=1; return
!    end select
!
!    call goLabel; status=0
!
!  end subroutine subr2
!
!
!
!end module testmod
!
!
! ################################################################
!
!
!program test
!
!  use go_print
!  use testmod
!
!  ! --- local -----------------------------------------
!
!  integer           ::  status
!
!  ! --- begin ------------------------------------------
!
!  call GO_Print_Init( status, trace=.false. )
!  call goLabel('test prog')
!
!  write (csol,'("begin of program")'); call csoPr
!
!  call Subr( 2, status )
!  if (status/=0) then; call csoErr; call exit(1); end if
!
!  write (csol,'("end of program")'); call csoPr
!
!  call goLabel()
!  call GO_Print_Done( status )
!
!end program test
!
