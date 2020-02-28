!#################################################################
!
! Usage:
!
!   call goDirname( '/storage/data/file.txt', dirname, status )
!      Return directory part of filename.
!
!   call goCreateDir( '/storage/data', status )
!      Create directory if not present yet.
!
!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#################################################################

module GO_Path

  use GO_Print, only : gol, goPr, goErr

  implicit none

  ! --- in/out -----------------------------

  private

  public  ::  pathsep
  public  ::  goDirname


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'GO_Path'
  
  ! path seperation character:
  character(len=*), parameter  ::  pathsep = '/'

  

contains


  !**********************************************************************

  
  subroutine goDirname( filename, dirname, status )
    
    ! --- in/out ----------------------------

    character(len=*), intent(in)      ::  filename
    character(len=*), intent(out)     ::  dirname
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goDirname'

    ! --- local -----------------------------

    integer                     ::  pos

    ! --- begin -----------------------------
    
    ! postion of last seperation character:
    pos = index( filename, pathsep, back=.true.)
    ! any path seperator?
    if ( pos > 1 ) then
      ! copy part:
      dirname = filename(1:pos-1)
    else
      ! empty:
      dirname = ''
    end if
    
    ! ok
    status = 0
    
  end subroutine goDirname
  
  ! *
  
  subroutine goCreateDir( dirname, status )
    
    ! --- in/out ----------------------------

    character(len=*), intent(in)      ::  dirname
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goCreateDir'

    ! --- local -----------------------------

    logical                         :: exists
    character(len=1024)             :: command
    integer                         :: exitstat
    integer                         :: cmdstat
    character(len=1024)             :: cmdmsg
    
    ! --- begin -----------------------------
    
    ! check if directory exists already:
    inquire( file=trim(dirname)//pathsep//'.', exist=exists )
    ! not yet?
    if ( .not. exists ) then
      ! create including parent directories:
      command = 'mkdir -p '//trim(dirname)
      call Execute_Command_Line( command, exitstat=exitstat, &
                                  cmdstat=cmdstat, cmdmsg=cmdmsg )
      if ( cmdstat /= 0 ) then
        write (gol,'("could not execute command:")'); call goErr
        write (gol,'("  ",a)') trim(command); call goErr
        TRACEBACK; status=1; return
      else if ( exitstat /= 0 ) then
        write (gol,'("non-zero exit status from command:")'); call goErr
        write (gol,'("  ",a)') trim(command); call goErr
        write (gol,'("exit status     : ",i0)') exitstat; call goErr
        write (gol,'("command message : ",a)') trim(cmdmsg); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0
    
  end subroutine goCreateDir
  
  ! *
  
  subroutine goCheckDir( filename, status )
    
    ! --- in/out ----------------------------

    character(len=*), intent(in)      ::  filename
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goCheckDir'

    ! --- local -----------------------------

    character(len=1024)     ::  dirname
    
    ! --- begin -----------------------------
    
    ! directory part:
    call goDirname( filename, dirname, status )
    IF_NOTOK_RETURN(status=1)
    ! defined?
    if ( len_trim(dirname) > 0 ) then
      ! create if necessary:
      call goCreateDir( dirname, status )
      IF_NOTOK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine goCheckDir


end module GO_Path
