!#################################################################
!
! File tools.
!  
!#################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!#################################################################

module CSO_File

  use CSO_Logging, only : csol, csoPr, csoErr

  implicit none

  ! --- in/out ---------------------

  private
  
  public  ::  CSO_GetFU
  public  ::  CSO_GetDirname
  public  ::  T_CSO_TextFile


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'CSO_File'

  ! standard file units
  integer, parameter    ::  goStdErr = 0
  integer, parameter    ::  goStdIn  = 5
  integer, parameter    ::  goStdOut = 6
  
  ! range of file units that might be used by this program:
  integer, parameter    ::  goFuRange(2) = (/200,999/)


  ! --- types ---------------------------------

  type T_CSO_TextFile
    character(len=200)      ::  name
    ! file unit:
    integer                 ::  fu
    ! comment ?
    logical                 ::  commented
    character(len=1)        ::  comment
  contains
    procedure ::  Init        =>  TextFile_Init
    procedure ::  Done        =>  TextFile_Done
    procedure ::  ReadLine    =>  TextFile_ReadLine
  end type T_CSO_TextFile



contains

 

  ! ==============================================================
  ! ===
  ! === file unit
  ! ===
  ! ==============================================================

  
  ! Return the first free available file unit number.

  subroutine CSO_GetFU( fu, status )
  
    ! --- in/out --------------------------

    integer, intent(out)      ::  fu
    integer, intent(out)      ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_GetFU'
    
    ! --- local --------------------------

    integer               ::  i
    character(len=256)    ::  fname
    logical               ::  opened

    ! --- local ---------------------------
    
    ! start with lowest possible unit:
    fu = goFuRange(1) - 1
    
    ! loop until unopned unit is found:
    do

      ! try next file unit:
      fu = fu + 1

      ! too large ?
      if ( fu > goFuRange(2) ) then
        write (csol,'("unable to select free file unit within allowed range ...")'); call csoErr
        write (csol,'("close some files or increase goFuRange in module GO_FU")'); call csoErr
        write (csol,'("current goFuRange : ",i6," .. ",i6)') goFuRange; call csoErr
        write (csol,'("open files:")'); call csoErr
        do i = goFuRange(1), goFuRange(2)
          inquire( unit=i, name=fname )
          write (csol,'(i6," : ",a)') i, trim(fname); call csoErr
        end do
        TRACEBACK; status=1; return
      end if
      
      ! skip ?
      if ( fu==goStdIn  ) cycle
      if ( fu==goStdOut ) cycle
      if ( fu==goStdErr ) cycle

      ! free available unit ? then ok
      inquire( unit=fu, opened=opened )
      if ( .not. opened ) exit
      
    end do

    ! ok
    status = 0

  end subroutine CSO_GetFU

 

  ! ==============================================================
  ! ===
  ! === paths
  ! ===
  ! ==============================================================

  
  ! return directory part of filename

  subroutine CSO_GetDirname( filename, dirname, status )
  
    ! --- in/out --------------------------

    character(len=*), intent(in)    ::  filename
    character(len=*), intent(out)   ::  dirname
    integer, intent(out)            ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_GetFU'
    
    ! --- local --------------------------

    integer               ::  i

    ! --- local ---------------------------
    
    ! search backward for path sep:
    i = index( filename, '/', back=.true. )
    ! found?
    if ( i > 1 ) then
      ! part up to sep:
      dirname = filename(1:i-1)
    else 
      ! empty:
      dirname = ''
    end if
    
    ! ok
    status = 0

  end subroutine CSO_GetDirname


  ! ==============================================================
  ! ===
  ! === text file
  ! ===
  ! ==============================================================


  !
  ! call Init( file, filename, iostat, [,status='unknown'|'old'|'new'] [,comment='\%'] )
  !
  ! Replaces the intrinsic 'open' command, but uses a
  ! a structure of type T_CSO_TextFile instead of a file unit number. \\
  ! Arguments passed are the same as for 'open'.\\
  ! In addition, a text file can be opened as a commented
  ! text file; with the 'ReadLine' command one is able to read
  ! lines from the file while skipping the lines starting
  ! with the specified comment.
  !

  subroutine TextFile_Init( self, filename, iostat, status, comment )

    ! --- in/out ------------------------

    class(T_CSO_TextFile), intent(out)            ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  iostat
    
    character(len=*), intent(in), optional    ::  status
    character(len=1), intent(in), optional    ::  comment

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/file_Init'
    
    ! --- local ----------------------------

    logical             ::  exist    
    character(len=10)   ::  statusX

    ! --- begin ----------------------------

    ! file exist ?
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (csol,'("commented text file not found:")'); call csoErr
      write (csol,'("  file name : ",a)') trim(filename); call csoErr
      TRACEBACK; iostat=1; return
    end if

    ! check file status : 'old', 'new', 'unknown'
    if (present(status)) then
      statusX = status
    else
      statusX = 'unknown'
    end if

    ! store filename:
    self%name = filename

    ! select free file unit:
    Call CSO_GetFU( self%fu, iostat )
    if (iostat/=0) then; TRACEBACK; iostat=1; return; end if

    ! open file:
    open( unit=self%fu, file=trim(filename), iostat=iostat, &
                                 status=statusX, form='formatted' )
    if ( iostat /= 0 ) then
      write (csol,'("from file open :")'); call csoErr
      write (csol,'("  file name : ",a)') trim(filename); call csoErr
      TRACEBACK; iostat=1; return; call csoErr
    end if
    
    ! check on comment lines ?
    if ( present(comment) ) then
      self%commented = .true.
      self%comment = comment
    else
      self%commented = .false.
      self%comment = 'x'
    end if

    ! ok
    iostat = 0

  end subroutine TextFile_Init


  ! ***

  subroutine TextFile_Done( self, status )

    ! --- in/out -----------------

    class(T_CSO_TextFile), intent(inout)   ::  self
    integer, intent(out)               ::  status

    ! --- const ----------------------
    
    character(len=*), parameter  ::  rname = mname//'/TextFile_Done'
    
    ! --- begin ------------------------

    ! close file:
    close( unit=self%fu, iostat=status )
    if ( status /= 0 ) then
      write (csol,'("from closing file:")'); call csoErr
      write (csol,'("  ",a)') trim(self%name); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine TextFile_Done


  ! ***


  !
  ! call ReadLine( file, s )
  !
  ! Reads the next line from a commented text file,
  ! but skips all lines starting with the 'comment'
  ! specified with the 'Init' command.
  ! Empty lines are skipped too.
  !

  subroutine TextFile_ReadLine( self, s, status  )

    ! --- in/out -------------------------

    class(T_CSO_TextFile), intent(inout)    ::  self
    character(len=*), intent(out)       ::  s
    integer, intent(out)                ::  status

    ! --- const --------------------------
    
    character(len=*), parameter  ::  rname = mname//'/TextFile_ReadLine'
    
    ! --- local --------------------------

    character(len=10)        ::  fmt

    ! --- begin --------------------------

    ! format (a100) etc:
    write (fmt,'("(a",i6.6,")")') len(s)

    ! loop until:
    !  o uncommented line has been read in s
    !  o eof is reached
    do
     
      ! read next line:
      read (self%fu,fmt,iostat=status) s
      if ( status < 0 ) then  ! eof
        s = ''
        status=-1; return
      else if ( status > 0 ) then
        write (csol,'("reading line from file:")'); call csoErr
        write (csol,'("  ",a)') trim(self%name); call csoErr
        TRACEBACK; status=1; return
      end if

      ! remove leading space:
      s = adjustl( s )

      ! empty ?
      if ( len_trim(s) == 0 ) cycle

      ! check for comment ?
      if ( self%commented .and. (scan(s,self%comment)==1) ) cycle
      
      ! s filled; leave loop
      exit
      
    end do
    
    ! ok
    status = 0

  end subroutine TextFile_ReadLine


end module CSO_File
