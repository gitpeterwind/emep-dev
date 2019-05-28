!###############################################################################
!
! NAME
!   GO_System  -  machine and/or compiler specific stuff
!
!
! DESCRIPTION
!
!   The module GO_System provides some basic constants for the
!   current compiler. In addition, some interfaces are defined
!   to routines for system calls, setting of exit statuses etc, 
!   which are non-standard Fortran, but often provided by the
!   vendor of the compiler. 
!   Since both constants and system routines differ from compiler
!   to compiler, this GO module is available in a number of copies,
!   each valid for a single compiler. If for some compiler a 
!   certain constant or system routine could not be filled,
!   a dummy value is used or a warning is issued.
!
!   The following system routines are defined:
!
!    o call goSystem( command, status )
!        Perform a system command, return exit status.
!
!    o call goExit( status )
!        Stop execution, set the exit status.
!
!    o call goArgCount( narg, status )
!        Count number of command line arguments.
!
!    o call goGetArg( nr, value, status )
!        Returns command line argument 'nr' in character string 'value'.
!
!    o call goSleep( nsec, status )
!        Wait for some seconds.
!
!    o call goMem( 'current memory usage (Mb) : ', status )
!
!
! GFORTRAN
!
!   Online manual:
!
!     gcc.gnu.org/onlinedocs/
!       (choose version)
!
!   Macro's defined from version 4.3.5 onwards
!   (section "Preprocessing and conditional compilation" in the manual):
!
!     __GFORTRAN__
!     __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__ 
!
!   According to the manual, this should work to list all defined macro's,
!   but when tested it only seemd to work with gcc:
!
!     gfortran -dM test.F90
!
!
! INTEL FORTRAN COMPILER
!
!   Online manual:
!
!       www.intel.com
!         Sitemap, Software, Find Product : Intel Compilers
!           Select 'Product Documentation' from the side menu
!             Intel Fortran Compiler 11.1 User and Reference Guides
!             (http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/index.htm)
!               Language Reference
!
!   Macro's defined:
!
!     __INTEL_COMPILER     ! evaluates to version number
!
!
! IBM XLF COMPILER
!
!   Online manuals are there but hard to locate; therefore some hints:
!
!     Language Reference - XL Fortran for AIX, V12.1
!       Service and utility procedures
!         General service and utility procedures
!
!   For the macro definitions, see:
!
!     XL C/C++ Compiler Reference Version 10.1
!       Compiler Predefined Macros
!         Macros indicating the XL C/C++ compiler product
!
!   Use '__IBMC__' rather than __xlc__ since the first evaluates
!   to a single integer number which can be tested with:
!     #if __IBMC__ == 1010
!   while the later evaluates to string like '10.1.0.4' .
!              
!   Compilation on ECMWF systems fails during linking because 'Exit_'
!   and 'Sleep_' could not be found. Although these are the standard names
!   according to the Compiler Reference, the ECMWF implementation only
!   recoqnizes 'Exit' and 'Sleep' (thus without underscores).
!   Adding the flag '-qnoextname' solves this, but induces a failure in
!   linking the HDF4 library. Therefore, a macro '__ecmwf__' should be
!   defined to distuinguish between XLF implementations at ECMWF and
!   at other institutes.
!
!
!###############################################################################
!
#define TRACEBACK write (gol,'("in ",a," (line",i5,")")') __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!###############################################################################

module GO_System

  use GO_Print, only : gol, goPr, goErr
  
  implicit none

  ! --- in/out ------------------------------
  
  private
  
  public   ::  goSystem
  public   ::  goExit
  public   ::  goArgCount, goGetArg
  public   ::  goSleep
  public   ::  goMemoryUsage
  public   ::  goMem
  
  public   ::  pathsep
  
  
  ! --- const ---------------------------------
  
  ! module name
  character(len=*), parameter  ::  mname = 'GO_System'
  
  ! path seperation:
  character(len=1), parameter  ::  pathsep = '/'
    

contains



  ! ############################################################################
  ! ###
  ! ### goSystem
  ! ###
  ! ############################################################################


  ! Execute a system command, return exit status.
  
  subroutine goSystem( command, status )
  
#ifdef __INTEL_COMPILER
    use IFPort, only : System
    use IFPort, only : iErrNo, E2BIG, ENOENT, ENOEXEC, ENOMEM
#endif
  
    ! --- in/out -----------------------------------------------
    
    character(len=*), intent(in)       ::  command
    integer, intent(inout)             ::  status
    
    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goSystem'
    
    ! --- local --------------------------------------------------
    
#ifdef __INTEL_COMPILER
    integer(4)     ::  stat
    integer(4)     ::  errno
#endif

    ! --- begin --------------------------------------------------
    
#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !   
    !   Intel Fortran Compiler 11.1 User and Reference Guides
    !     Language Reference
    !       A to Z Reference
    !         System
    !
    stat = System( command )
    
    ! trap errors in ifort system command
    if ( stat == -1 ) then
      write (gol,'("error in call to IFort Portability command `system`:")'); call goErr
      errno = iErrNo()
      select case ( errno )
        case ( E2BIG )
          write (gol,'("  ",a)') 'The argument list is too long.'; call goErr
        case ( ENOENT )
          write (gol,'("  ",a)') 'The command interpreter cannot be found.'; call goErr
        case ( ENOEXEC )
          write (gol,'("  ",a)') 'The command interpreter file has an invalid format and is not executable.'; call goErr
        case ( ENOMEM )
          write (gol,'("  ",a)') 'Not enough system resources are available to execute the command.'; call goErr
        case default
          write (gol,'("  unknown iErrNo ",i)') errno; call goErr
      end select
      TRACEBACK; status=stat; return
    end if
    
    ! if the shell command exit status is 'n',
    ! then the number returned by 'system' is 256 * n
    status = stat / 256

#else
#ifdef __GFORTRAN__ 

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GNU Fortran Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   gcc.gnu.org/onlinedocs/
    !     GCC 4.3.5 GNU Fortran 95 Manual
    !       6. Intrinsic Prodedures
    !         204. SYSTEM - Execute a shell command
    !
    call System( command, status )

#else
#ifdef __IBMC__ 

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! IBM XL Fortran Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !
    ! From manual:
    !
    !   Language Reference - XL Fortran for AIX, V12.1
    !     Intrinsic Procedures
    !
    call System( command, status )

#else

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !write (gol,'("could not evaluate system command : ",a)') trim(command); call goErr
    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
    !TRACEBACK; status=1; return

    ! try this, often works:
    call System( command, status )

#endif
#endif
#endif
  
  end subroutine goSystem
  
  
  ! ############################################################################
  ! ###
  ! ### goExit
  ! ###
  ! ############################################################################


  ! Stop execution, set exit status.
  
  subroutine goExit( status )
  
#ifdef __IBMC__
    use XLFUtility, only : Exit_
#endif

    ! --- in/out --------------------------------------------
    
    integer, intent(in)    ::  status
    
    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goExit'
    
    ! --- begin --------------------------------------------
    
#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !   
    !   Intel Fortran Compiler 11.1 User and Reference Guides
    !     Language Reference
    !       A to Z Reference
    !         Exit Subroutine
    !
    call Exit( status )

#else
#ifdef __GFORTRAN__
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GNU Fortran compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   gcc.gnu.org/onlinedocs/
    !     GCC 4.3.5 GNU Fortran 95 Manual
    !       6. Intrinsic Prodedures
    !         6.66. EXIT - Exit the program with status.
    !
    call Exit( status )
    
#else
#ifdef __IBMC__
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! IBM XLF compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   Language Reference - XL Fortran for AIX, V12.1
    !     Service and utility procedures
    !       General service and utility procedures
    !
#ifdef __ecmwf__
    call Exit( status )
#else
    call Exit_( status )
#endif
    
#else

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
    
    !! this is an emergency, so for one time, the Fortran stop is allowed ...
    !stop 'Fortran STOP in GO_System/goExit'

    ! try this, often works:
    call Exit( status )

#endif
#endif
#endif
  
  end subroutine goExit


  ! ############################################################################
  ! ###
  ! ### goArgCount
  ! ###
  ! ############################################################################


  ! Return number of command line arguments
  
  subroutine goArgCount( narg, status )
  
    ! --- in/out --------------------------------------------
    
    integer, intent(out)    ::  narg
    integer, intent(out)    ::  status

    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goArgCount'
    
    ! --- begin -------------------------------------------------
    
#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !   
    !   Intel Fortran Compiler 11.1 User and Reference Guides
    !     Language Reference
    !       A to Z Reference
    !         Command_Argument_Count
    !
    narg = Command_Argument_Count()
    
#else
#ifdef __GFORTRAN__

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GNU Fortran Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   gcc.gnu.org/onlinedocs/
    !     GCC 4.3.5 GNU Fortran 95 Manual
    !       6. Intrinsic Prodedures
    !         6.42. COMMAND_ARGUMENT_COUNT - Get number of command line arguments
    !
    narg = Command_Argument_Count()
    
#else
#ifdef __IBMC__
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! IBM XLF compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   Language Reference - XL Fortran for AIX, V12.1
    !     Intrinsic Procedures
    !
    narg = Command_Argument_Count()
    
#else

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !! always assign something ...
    !narg = -1

    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
    !TRACEBACK; status=1; return
    
    ! try this, often works:
    narg = iArgC()

#endif
#endif
#endif

    ! ok
    status = 0
  
  end subroutine goArgCount
  
  
  ! ############################################################################
  ! ###
  ! ### goGetArg
  ! ###
  ! ############################################################################



  ! Return a command line argument
  
  subroutine goGetArg( pos, value, status )

    ! --- in/out --------------------------------------------------
    
    integer, intent(in)             ::  pos
    character(len=*), intent(out)   ::  value
    integer, intent(inout)          ::  status
    
    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goGetArg'
    
    ! --- local -----------------------------------------------------
       
    ! --- begin -----------------------------------------------------

#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !   
    !   Intel Fortran Compiler 11.1 User and Reference Guides
    !     Language Reference
    !       A to Z Reference
    !         Get_Command_Argument
    !
    ! Following the F2003 standard:
    call Get_Command_Argument( pos, value=value, status=status )
    
#else
#ifdef __GFORTRAN__

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GNU Fortran Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   gcc.gnu.org/onlinedocs/
    !
    !     GCC 4.3.5 GNU Fortran 95 Manual
    !       6. Intrinsic Prodedures
    !         6.87. GET_COMMAND_ARGUMENT - Get command line arguments
    !                 call Get_Command_Argument( pos, value )
    !
    !     GCC 4.4.4 GNU Fortran 95 Manual
    !       7. Intrinsic Prodedures
    !         7.89. GET_COMMAND_ARGUMENT - Get command line arguments
    !                 call Get_Command_Argument( pos [,value, length, status] )
    !
#if __GNUC____GNUC_MINOR____GNUC_PATCHLEVEL__ == 435
    ! up to version 4.3.5 :
    call Get_Command_Argument( pos, value )
    ! no status returned ...
    status = 0
#else
    ! Following the F2003 standard from v4.4.4 onwards:
    call Get_Command_Argument( pos, value=value, status=status )
#endif

#else
#ifdef __IBMC__
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! IBM XLF compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   Language Reference - XL Fortran for AIX, V12.1
    !     Intrinsic Procedures
    !
    ! Following the F2003 standard:
    call Get_Command_Argument( pos, value=value, status=status )
    
#else

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !! use arguments to avoid compilation warnings:
    !status = pos
    !! always assign something ...
    !value = '?'
    
    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
    !TRACEBACK; status=1; return

    ! try this, often works:
    call GetArg( pos, value )
    ! no status returned ...
    status = 0

#endif
#endif
#endif
  
  end subroutine goGetArg
  
  
  ! ############################################################################
  ! ###
  ! ### goSleep
  ! ###
  ! ############################################################################


  ! wait some seconds ...
  
  subroutine goSleep( nsec, status )

#ifdef __INTEL_COMPILER
    use IFPort, only : Sleep
#endif
#ifdef __IBMC__
    use XLFUtility, only : Sleep_
#endif

    ! --- in/out --------------------------------------------------
    
    integer, intent(in)             ::  nsec
    integer, intent(out)            ::  status
    
    ! --- const ------------------------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goSleep'
    
    ! --- begin -----------------------------------------------------

#ifdef __INTEL_COMPILER

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Intel Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !   
    !   Intel Fortran Compiler 11.1 User and Reference Guides
    !     Language Reference
    !       A to Z Reference
    !         Sleep
    !
    call Sleep( nsec )

#else
#ifdef __GFORTRAN__

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! GNU Fortran Compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   gcc.gnu.org/onlinedocs/
    !     GCC 4.3.5 GNU Fortran 95 Manual
    !       6. Intrinsic Prodedures
    !         6.195. SLEEP - Sleep for the specified number of seconds
    !   
    call Sleep( nsec )

#else
#ifdef __IBMC__
    
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! IBM XLF compiler
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !
    ! From manual:
    !
    !   Language Reference - XL Fortran for AIX, V12.1
    !     Service and utility procedures
    !       General service and utility procedures
    !
#ifdef __ecmwf__
    call Sleep( nsec )
#else
    call Sleep_( nsec )
#endif
    
#else

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! error ...
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    !write (gol,'("subroutine not implemented for this compiler")'); call goErr
    !TRACEBACK; status=1; return
    
    ! try this, often works:
    call Sleep( nsec )

#endif
#endif
#endif

    ! ok
    status = 0
  
  end subroutine goSleep
  
  
  ! ***


  !
  ! Extact current memory usage from output of 'pmap <pid>' system call.
  !
  ! Alternative is to grep 'VmRSS' from the '/proc/<pid>/status' file, see:
  !   http://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90
  !
  ! However, RSS does not seem to be very useful, see the different memory counters in:
  !   https://utcc.utoronto.ca/~cks/space/blog/linux/LinuxMemoryStats
  ! Better use VSZ ?
  ! - seems in 'VmSize' in '/proc/<pid>/status' file.
  ! - first number of pmap total on vilje?
  !

  subroutine goMemoryUsage( usage, status, outline )

#ifdef __INTEL_COMPILER
    use IFPort, only : GetPid
    use IFPort, only : System
#endif
    use GO_Fu, only : goGetFU
  
    ! --- in/out ---------------------------------
    
    integer(8), intent(out)       ::  usage   ! bytes
    integer, intent(out)          ::  status
    character(len=*), optional    ::  outline
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goMemoryUsage'
    
    ! --- local ----------------------------------
    
    integer               ::  pid
    character(len=256)    ::  outfile
    character(len=256)    ::  command
    character(len=32)     ::  label, total
    character(len=1024)   ::  line
    integer               ::  fu
    integer               ::  l
    
    ! --- begin ----------------------------------
    
    !! wait ..
    !call Sleep( 2 )
    
    ! id of current processes:
    pid = GetPid()
    
    ! in future:
    !call Execute_Command_Line( 'pwd', exitsat=status, cmdmsg=msg )
    
    ! target file for commend output:
    write (outfile,'("/tmp/go-",i0,".out")') pid

    ! command to list total memory usage:
    write (command,'("pmap ",i0," | grep -i total > ",a)') pid, trim(outfile)
    ! run:
    status = System( command )
    if ( status /= 0 ) then
      write (gol,'("from system call:")'); call goErr
      write (gol,'("  command     : ",a)') trim(command); call goErr
      write (gol,'("  output log  : ",a)') trim(outfile); call goErr
      write (gol,'("  exit status : ",i0)') status; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! expected output:
    !    " total         12345k"
    ! get free file unit:
    call goGetFU( fu, status )
    IF_NOTOK_RETURN(status=1)
    ! open:
    open( fu, file=trim(outfile), status='old', iostat=status )
    if ( status /= 0 ) then
      write (gol,'("could not open output log: ",a)') trim(outfile); call goErr
      TRACEBACK; status=1; return
    end if
    ! read:
    read (fu,'(a)',iostat=status) line
    if ( status /= 0 ) then
      write (gol,'("could not read output line from log file: ",a)') trim(outfile); call goErr
      TRACEBACK; status=1; return
    end if
    ! done:
    close( fu, iostat=status )
    if ( status /= 0 ) then
      write (gol,'("closing log file: ",a)') trim(outfile); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! copy?
    if ( present(outline) ) outline = trim(line)
    
    ! read:
    read (line,*,iostat=status) label, total
    if ( status /= 0 ) then
      write (gol,'("could not read output from line: ",a)') trim(line); call goErr
      TRACEBACK; status=1; return
    end if

    ! extract and convert:
    l = len_trim(total)
    if ( total(l:l) == 'K' ) then
      read (total(1:l-1),*) usage
      usage = usage * 1000
    else
      write (gol,'("unsupported form `",a,"` of total")') trim(total); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine goMemoryUsage
  

  ! ***


  subroutine goMem( label, status )

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)    ::  label
    integer, intent(out)            ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goMem'
    
    ! --- local ----------------------------------
    
    integer(8)            ::  usage
    character(len=1024)   ::  outline
    
    ! --- begin ----------------------------------
    
    ! current memory usage in bytes:
    call goMemoryUsage( usage, status, outline=outline )
    IF_NOTOK_RETURN(status=1)
    
    ! show result in Mb:
    !write (gol,'(a," ",f0.3," Mb")') label, usage*1e-6; call goPr
    write (gol,'(a," ",a)') label, trim(outline); call goPr
    
    ! ok
    status = 0
    
  end subroutine goMem


end module GO_System


!! ######################################################################
!! ###
!! ### test
!! ###
!! ######################################################################
!!
!!  gfortran -o test.x go_fu.F90 go_print.F90 go_system.F90 && ./test.x a bc ; echo $?
!!
!!  ifort -o test.x go_fu.F90 go_print.F90 go_system.F90 && ./test.x a bc ; echo $?
!!
!!  xlf -qnoextname -o test.x go_fu.F90 go_print.F90 go_system.F90 && ./test.x a bc ; echo $?
!!
!!  xlf -qnoextname -o test.x -WF,-D__ecmwf__ go_fu.F90 go_print.F90 go_system.F90 && ./test.x a bc ; echo $?
!!
!#define IF_NOTOK_STOP if (status/=0) then; TRACEBACK; stop; end if
!!
!program test
!
!  use GO_Print, only : gol, goPr, goErr
!  use GO_System, only : goArgCount, goGetArg, goExit, goSleep, goSystem
!  
!  implicit none
!  
!  character(len=*), parameter  ::  rname = 'test'
!  
!  integer             ::  status
!  integer             ::  n, i
!  character(len=32)   ::  val
!  character(len=32)   ::  command
!  
!  print *, 'begin'
!
!  print *, ''
!  print *, 'Test value of predefined macros:'
!#ifdef __GFORTRAN__
!  print *, '  __GFORTRAN__          : ', __GFORTRAN__
!  print *, '  __GNUC__              : ', __GNUC__
!  print *, '  __GNUC_MINOR__        : ', __GNUC_MINOR__
!  print *, '  __GNUC_PATCHLEVEL__   : ', __GNUC_PATCHLEVEL__
!#else
!  print *, '  __GFORTRAN__          : undefined'
!#endif
!#ifdef __INTEL_COMPILER
!  print *, '  __INTEL_COMPILER      : ', __INTEL_COMPILER
!#else
!  print *, '  __INTEL_COMPILER      : undefined'
!#endif
!#ifdef __IBMC__
!  print *, '  __IBMC__              : ', __IBMC__
!#else
!  print *, '  __IBMC__              : undefined'
!#endif
!
!  print *, ''
!  print *, 'number of arguments ...'
!  call goArgCount( n, status )
!  IF_NOTOK_STOP
!  print *, n
!  
!  print *, 'get arguments ...'
!  do i = 1, n
!    print *, '  argument ', i, ' ...'
!    call goGetArg( i, val, status )
!    IF_NOTOK_STOP
!    print *, '"'//trim(val)//'"'
!  end do
!  
!  command = '/bin/ls -l'
!  print *, ''
!  print *, 'call system command : ', trim(command)
!  call goSystem( command, status )
!  IF_NOTOK_STOP
!  
!  n = 2
!  print *, ''
!  print *, 'wait ', n, ' seconds ...'
!  call goSleep( n, status )
!  IF_NOTOK_STOP
!  
!  status = 23
!  print *, ''
!  print *, 'exit with status ', status, ' ...'
!  call goExit( status )
!  
!  print *, 'end'
!
!end program test

