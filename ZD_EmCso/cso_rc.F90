!#################################################################
!
! call ReadRc( rcfile, 'test.flag', l, status [,default=.false.] )
!
! return status :
!   <0  : key not found, value set to default
!    0  : key found and value read without errors
!   >0  : some errors
!
! Search for extended keys:
!
!   call ReadRc( rcfile, 'test', (/'*  ','all','b  '/), flag, status, default=.true. )
!
! will search for (dots are inserted automatically):
!
!     test.*       :  F 
!     test.all     :  F 
!     test.b       :  T 
!
! The last found key overwrites all previous values.
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

module CSO_Rc

  use CSO_Logging, only : csol, csoPr, csoErr

  implicit none

  ! --- in/out ---------------------

  private
  
  public  ::  T_CSO_RcFile


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'CSO_Rc'

  ! maximum line length in rc file:
  integer, parameter     ::  buflen = 2000

  ! --- types ---------------------------------

  ! key/value pair:
  type T_KeyAndValue
    ! key and value:
    character(len=:), allocatable   ::  key
    character(len=:), allocatable   ::  value
  contains
    procedure ::  Init        =>  KeyAndValue_Init
    procedure ::  Done        =>  KeyAndValue_Done
  end type T_KeyAndValue

  ! file:
  type T_CSO_RcFile
    ! input file:
    character(len=:), allocatable       ::  fname
    ! key/value pairs:
    integer                             ::  n
    type(T_KeyAndValue), allocatable    ::  element(:)
    !
  contains
    procedure ::  Init        =>  RcFile_Init
    procedure ::  Done        =>  RcFile_Done
    procedure ::  FindIndex   =>  RcFile_FindIndex
    !
    procedure ::                  RcFile_Get_i
    procedure ::                  RcFile_Get_i1
    procedure ::                  RcFile_Get_r
    procedure ::                  RcFile_Get_l
    procedure ::                  RcFile_Get_s
    procedure ::                  RcFile_Get_datetime
    generic   ::  Get         =>  RcFile_Get_i , &
                                  RcFile_Get_i1, &
                                  RcFile_Get_r , &
                                  RcFile_Get_l , &
                                  RcFile_Get_s , &
                                  RcFile_Get_datetime
    !
  end type T_CSO_RcFile



contains


  ! ================================================================
  ! ===
  ! === key and value pairs
  ! ===
  ! ================================================================


  subroutine KeyAndValue_Init( self, key, value, status )

   ! --- in/out ---------------------------

    class(T_KeyAndValue), intent(out)   ::  self
    character(len=*), intent(in)        ::  key
    character(len=*), intent(in)        ::  value
    integer, intent(out)                ::  status
    
    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/KeyAndValue_Init'
    
    ! --- local --------------------------

    ! --- begin ---------------------------
    
    ! storage:
    allocate( character(len=len_trim(key)) :: self%key, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    self%key = trim(key)

    ! storage:
    allocate( character(len=len_trim(value)) :: self%value, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    self%value = trim(value)

    ! ok
    status = 0

  end subroutine KeyAndValue_Init


  ! ***


  subroutine KeyAndValue_Done( self,  status )

    ! --- in/out ---------------------------

    class(T_KeyAndValue), intent(inout)     ::  self
    integer, intent(out)                    ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/KeyAndValue_Done'
    
    ! --- local ---------------------------

    ! --- begin ---------------------------
    
    ! clear:
    deallocate( self%key, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%value, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine KeyAndValue_Done




  ! ================================================================
  ! ===
  ! === RcFile
  ! ===
  ! ================================================================


  subroutine RcFile_Init( self, filename, status )

    use CSO_String, only : CSO_SplitLine, CSO_Tab2Space
    use CSO_File  , only : T_CSO_TextFile

   ! --- in/out ---------------------------

    class(T_CSO_RcFile), intent(out)  ::  self
    character(len=*), intent(in)      ::  filename
    integer, intent(out)              ::  status
    
    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Init'
    
    ! --- local --------------------------
    
    logical                 ::  exist
    type(T_CSO_TextFile)    ::  file
    character(len=buflen)   ::  line
    integer                 ::  i
    integer                 ::  j
    character(len=buflen)   ::  key
    character(len=buflen)   ::  value

    ! --- begin ---------------------------

    ! file not present ?
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (csol,'("rcfile not found :")'); call csoErr
      write (csol,'("  ",a)') trim(filename); call csoErr
      TRACEBACK; status=1; return
    end if

    ! storage:
    allocate( character(len=len_trim(filename)) :: self%fname, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    self%fname = trim(filename)

    ! first pass: count number of key/value lines:
    self%n = 0
    ! open commented text file:
    call file%Init( trim(self%fname), status, status='old', comment='!' )
    IF_NOT_OK_RETURN(status=1)
    ! scan all lines 
    do
      ! read next non empty, non comment line:
      call file%ReadLine( line, status )
      if (status<0) exit  ! end of file
      IF_NOT_OK_RETURN(status=1)
      ! increase counter:
      self%n = self%n + 1
    end do
    ! close:
    call file%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! any values ?
    if ( self%n > 0 ) then
    
      ! storage:
      allocate( self%element(self%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! open commented text file:
      call file%Init( trim(self%fname), status, status='old', comment='!' )
      IF_NOT_OK_RETURN(status=1)
      ! loop over lines: 
      do i = 1, self%n
      
        ! read next non empty, non comment line:
        call file%ReadLine( line, status )
        IF_NOT_OK_RETURN(status=1)
      
        ! replace tabs:
        call CSO_Tab2Space( line )

        ! split at colon:
        call CSO_SplitLine( line, key, ':', value, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! check on doubles?
        if ( i > 1 ) then
          ! loop over current elements:
          do j = 1, i-1
            ! match?
            if ( trim(self%element(j)%key) == trim(key) ) then
              write (csol,'("found key `",a,"` at least twice in file:")') trim(key); call csoErr
              write (csol,'("  ",a)') trim(self%fname); call csoErr
              TRACEBACK; status=1; return
            end if  ! same?
          end do  ! previous elements
        end if ! check on doubles
        
        ! init key/value pair:
        call self%element(i)%Init( key, value, status )
        IF_NOT_OK_RETURN(status=1)
        
      end do
      ! close:
      call file%Done( status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! n > 0
    
    ! ok
    status = 0

  end subroutine RcFile_Init


  ! ***


  subroutine RcFile_Done( self,  status )

    ! --- in/out ---------------------------

    class(T_CSO_RcFile), intent(inout)    ::  self
    integer, intent(out)              ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Done'
    
    ! --- local ---------------------------
    
    integer     ::  i

    ! --- begin ---------------------------
    
    ! any elements?
    if ( self%n > 0 ) then
      ! loop over elements:
      do i = 1, self%n
        ! done:
        call self%element(i)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do ! elements
      ! clear:
      deallocate( self%element, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! n > 0

    ! clear:
    deallocate( self%fname, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine RcFile_Done


  ! ***


  ! 
  ! Return index of element with specified key,
  ! or negative number if not found.
  !
  
  subroutine RcFile_FindIndex( self, key, ind, status )

    ! --- in/out ---------------------------

    class(T_CSO_RcFile), intent(in)       ::  self
    character(len=*)                  ::  key
    integer, intent(out)              ::  ind
    integer, intent(out)              ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_FindIndex'
    
    ! --- local ---------------------------
    
    integer     ::  i

    ! --- begin ---------------------------
    
    ! init return value as not found:
    ind = -999
    
    ! any elements?
    if ( self%n > 0 ) then
      ! loop over elements:
      do i = 1, self%n
        ! compare:
        if ( trim(key) == self%element(i)%key ) then
          ! fill return value:
          ind = i
          ! leave:
          exit
        end if
      end do ! elements
    end if ! any elements
    
    ! ok
    status = 0

  end subroutine RcFile_FindIndex
  
  
  ! ***


  subroutine RcFile_Get_i( self, key, i, status, default )

    ! --- in/out ----------------------------

    class(T_CSO_RcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    integer, intent(out)                        ::  i
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_i'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOT_OK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        i = default
        status = -1 ; return
      else
        write (csol,'("key not found and no default specified ...")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      read (self%element(ind)%value,*,iostat=status) i
      if ( status /= 0 ) then
        write (csol,'("while reading integer:")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        write (csol,'("  value  : ",a)') trim(self%element(ind)%value); call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_i

  ! *

  subroutine RcFile_Get_i1( self, key, i, status, default )

    ! --- in/out ----------------------------

    class(T_CSO_RcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    integer, intent(out)                        ::  i(:)
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_i1'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOT_OK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        i = default
        status = -1 ; return
      else
        write (csol,'("key not found and no default specified ...")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      read (self%element(ind)%value,*,iostat=status) i
      if ( status /= 0 ) then
        write (csol,'("while reading integer:")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        write (csol,'("  value  : ",a)') trim(self%element(ind)%value); call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_i1


  ! ***


  subroutine RcFile_Get_r( self, key, r, status, default )

    ! --- in/out ----------------------------

    class(T_CSO_RcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    real, intent(out)                           ::  r
    integer, intent(out)                        ::  status
    
    real, intent(in), optional                  ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_r'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOT_OK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        r = default
        status = -1 ; return
      else
        write (csol,'("key not found and no default specified ...")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      read (self%element(ind)%value,*,iostat=status) r
      if ( status /= 0 ) then
        write (csol,'("while reading real :")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        write (csol,'("  value  : ",a)') trim(self%element(ind)%value); call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_r
  

  ! ***
  
  
  subroutine RcFile_Get_l( self, key, l, status, default )

    ! --- in/out ----------------------------

    class(T_CSO_RcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    logical, intent(out)                        ::  l
    integer, intent(out)                        ::  status
    
    logical, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_l'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOT_OK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with warning:
      if ( present(default) ) then
        l = default
        status = -1 ; return
      else
        write (csol,'("key not found and no default specified ...")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      read (self%element(ind)%value,*,iostat=status) l
      if ( status /= 0 ) then
        write (csol,'("while reading logical :")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        write (csol,'("  value  : ",a)') trim(self%element(ind)%value); call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_l


  ! ***


  subroutine RcFile_Get_s( self, key, s, status, default )
  
    ! --- in/out ----------------------------

    class(T_CSO_RcFile), intent(in)                 ::  self
    character(len=*), intent(in)                ::  key
    character(len=*), intent(out)               ::  s
    integer, intent(out)                        ::  status
    
    character(len=*), intent(in), optional      ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_s'

    ! --- local -----------------------------

    integer     ::  ind

    ! --- begin -----------------------------

    ! search index for specified key:
    call self%FindIndex( key, ind, status )
    IF_NOT_OK_RETURN(status=1)
    ! not found?
    if ( ind < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        s = trim(default)
        status = -1 ; return
      else
        write (csol,'("key not found and no default specified ...")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! key found; set value:
      s = trim(self%element(ind)%value)
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_s


  ! ***


  subroutine RcFile_Get_datetime( self, key, t, status, default )
  
    use CSO_DateTimes, only : T_CSO_DateTime, CSO_ReadFromLine
  
    ! --- in/out ----------------------------

    class(T_CSO_RcFile), intent(in)             ::  self
    character(len=*), intent(in)                ::  key
    type(T_CSO_DateTime), intent(out)           ::  t
    integer, intent(out)                        ::  status
    
    type(T_CSO_DateTime), intent(in), optional  ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/RcFile_Get_datetime'

    ! --- local -----------------------------

    character(len=64)     ::  line

    ! --- begin -----------------------------
    
    ! read character key, status <0 if default is used:
    call self%Get( key, line, status, default='None' )
    IF_ERROR_RETURN(status=1)
    if ( status < 0 ) then
      ! default?
      if ( present(default) ) then
        ! copy:
        t = default
        ! return with warning status:
        status = -1 ; return
      else
        write (csol,'("key not found and no default specified ...")'); call csoErr
        write (csol,'("  rcfile : ",a)') trim(self%fname); call csoErr
        write (csol,'("  key    : ",a)') trim(key); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! extract date/time from line:
      call CSO_ReadFromLine( line, t, status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0

  end subroutine RcFile_Get_datetime


end module CSO_Rc
