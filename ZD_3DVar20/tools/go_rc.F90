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
#define TRACEBACK write (gol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#################################################################

module GO_Rc

  implicit none

  ! --- in/out ---------------------

  private

  public  ::  TrcFile
  public  ::  Init, Done
  public  ::  ReadRc


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'GO_Rc'

  ! maximum line length in rc file:
  integer, parameter     ::  buflen = 1024

  ! --- types ---------------------------------

  type TrcFile
    character(len=80)      ::  fname
  end type TrcFile


  ! --- interfaces -------------------------------------

  interface Init
    module procedure rcfile_Init
  end interface

  interface Done
    module procedure rcfile_Done
  end interface

  interface ReadRc
    module procedure ReadRc_i
    module procedure ReadRcs_i
    module procedure ReadRc_i1
    module procedure ReadRc_r
    module procedure ReadRcs_r
    module procedure ReadRc_r1
    module procedure ReadRc_l
    module procedure ReadRcs_l
    module procedure ReadRc_s
    module procedure ReadRcs_s
  end interface


contains


  ! ================================================================
  ! ===
  ! === init, done
  ! ===
  ! ================================================================


  subroutine rcfile_Init( rcfile, fname, status )

    use GO_Print, only : gol, goErr

    ! --- in/out ---------------------------

    type(TrcFile), intent(out)    ::  rcfile
    character(len=*), intent(in)  ::  fname
    integer, intent(out)          ::  status
    
    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/rcfile_Init'
    
    ! --- local --------------------------
    
    logical          ::  exist

    ! --- begin ---------------------------

    ! file not present ?
    inquire( file=trim(fname), exist=exist )
    if ( .not. exist ) then
      write (gol,'("rcfile not found :")'); call goErr
      write (gol,'("  ",a)') trim(fname); call goErr
      TRACEBACK; status=1; return
    end if

    ! store file name:
    rcfile%fname = trim(fname)
    
    ! ok
    status = 0

  end subroutine rcfile_Init


  ! ***


  subroutine rcfile_Done( rcfile,  status )

    ! --- in/out ---------------------------

    type(TrcFile), intent(inout)    ::  rcfile
    integer, intent(out)            ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/rcfile_Done'
    
    ! --- begin ---------------------------

    ! dummy ...
    rcfile%fname = ''
    
    ! ok
    status = 0

  end subroutine rcfile_Done


  ! ================================================================
  ! ===
  ! === general read
  ! ===
  ! ================================================================


  ! Searches the file <filenameResource> for the string
  !  "<key> : "
  ! and save all characters behind the equal sign in <buffer>.
  ! The Resource file may contain comment lines starting with a "!"

  subroutine ReadRcItem( rcfile, key, buffer, status )

    use GO_Print , only : gol, goErr
    use GO_String, only : goSplitLine, goTab2Space
    use GO_File  , only : TTextFile, Init, Done, ReadLine

    ! --- in/out ----------------------

    type(TrcFile), intent(in)         ::  rcfile
    character(len=*), intent(in)      ::  key
    character(len=*), intent(out)     ::  buffer
    integer, intent(out)              ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRcItem'
    
    ! --- local -----------------------

    type(TTextFile)         ::  file
    Integer                 ::  nfound
    character(len=buflen)   ::  s, skey, sdata

    ! --- begin --------------------------

    ! open commented text file:
    call Init( file, rcfile%fname, status, status='old', comment='!' )
    IF_NOTOK_RETURN(status=1)

    ! no matching lines found yet ...    
    nfound = 0
    
    ! scan all lines 
    do

      ! read next non empty, non comment line:
      call ReadLine( file, s, status )
      if (status<0) exit  ! end of file
      IF_NOTOK_RETURN(status=1)
      
      ! Andy Jacobson, 10 Apr 2006.  Allows tabs in rc file.
      call goTab2Space( s )  

      ! split at colon:
      call goSplitLine( s, skey, ':', sdata, status )
      IF_NOTOK_RETURN(status=1)
      
      ! starts with requested key, and no extra text between key and colon ? then found!
      if ( (index(skey,key)==1) .and. (len_trim(key)==len_trim(skey))) then
        buffer = sdata
        nfound = nfound + 1
      end if

    end do

    ! close:
    call Done( file, status )
    IF_NOTOK_RETURN(status=1)
    
    ! not found ? warning status
    if ( nfound == 0 ) then
      status=-1; return
    end if

    ! multiple matches ?
    if ( nfound > 1 ) then
      write (gol,'("found more than one matching keys in rcfile:")'); call goErr
      write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
      write (gol,'("  key    : ",a)') trim(key); call goErr
      write (gol,'("  found  : ",i4," times")') nfound
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0

  end subroutine ReadRcItem
  
  
  ! ================================================================
  ! ===
  ! === integer
  ! ===
  ! ================================================================


  subroutine ReadRc_i( rcfile, key, i, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    integer, intent(out)                        ::  i
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRc_i'

    ! --- local -----------------------------

    character(len=buflen)     ::  buffer

    ! --- begin -----------------------------

    ! search key line in rcfile:
    call ReadRcItem( rcfile, key, buffer, status )
    if ( status < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        i = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else if ( status == 0 ) then
      ! key found; set value:
      read (buffer,*,iostat=status) i
      if ( status /= 0 ) then
        write (gol,'("while reading integer:")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(buffer); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! some error ...
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine ReadRc_i


  ! ***

  
  subroutine ReadRcs_i( rcfile, key, keys, i, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    character(len=*), intent(in)                ::  keys(:)
    integer, intent(out)                        ::  i
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRcs_i'

    ! --- local -----------------------------
    
    logical                 ::  found
    integer                 ::  ikey
    integer                 ::  i_curr

    ! --- begin -----------------------------
    
    ! pessimistic assumption ...
    found = .false. 
    
    ! loop over all key extensions:
    do ikey = 1, size(keys)
    
      ! try to read key; 
      ! provide default to return without error if key is not found:
      call ReadRc( rcfile, trim(key)//'.'//trim(keys(ikey)), i_curr, status, default=0 )
      if ( status < 0 ) then
        ! not found; try next
        cycle
      else if ( status == 0 ) then
        ! found and value read:
        found = .true.
        i = i_curr
      else
        ! error ...
        TRACEBACK; status=1; return
      end if
      
    end do   ! loop over keys
    
    ! not found ?
    if ( .not. found ) then
      ! default provided ?
      if ( present(default) ) then
        ! set to default:
        i = default
      else
        ! error ...
        write (gol,'("key(s) not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        do ikey = 1, size(keys)
          write (gol,'("  key    : ",a,".",a)') trim(key), trim(keys(ikey)); call goErr
        end do
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine ReadRcs_i


  ! ***


  subroutine ReadRc_i1( rcfile, key, i, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    integer, intent(out)                        ::  i(:)
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRc_i1'

    ! --- local -----------------------------

    character(len=buflen)     ::  buffer

    ! --- begin -----------------------------

    ! search key line in rcfile:
    call ReadRcItem( rcfile, key, buffer, status )
    if ( status < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        i = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else if ( status == 0 ) then
      ! key found; set value:
      read (buffer,*,iostat=status) i
      if ( status /= 0 ) then
        write (gol,'("while reading integer:")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(buffer); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! some error ...
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine ReadRc_i1


  ! ================================================================
  ! ===
  ! === real
  ! ===
  ! ================================================================


  subroutine ReadRc_r( rcfile, key, r, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    real, intent(out)                           ::  r
    integer, intent(out)                        ::  status
    
    real, intent(in), optional                  ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRc_r'

    ! --- local -----------------------------

    character(len=buflen)     ::  buffer

    ! --- begin -----------------------------

    ! search key line in rcfile:
    call ReadRcItem( rcfile, key, buffer, status )
    if ( status < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        r = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else if ( status == 0 ) then
      ! key found; set value:
      read (buffer,*,iostat=status) r
      if ( status /= 0 ) then
        write (gol,'("while reading real :")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(buffer); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! some error ...
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine ReadRc_r
  

  ! ***

  
  subroutine ReadRcs_r( rcfile, key, keys, r, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    character(len=*), intent(in)                ::  keys(:)
    real, intent(out)                           ::  r
    integer, intent(out)                        ::  status
    
    real, intent(in), optional                  ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRcs_r'

    ! --- local -----------------------------
    
    logical                 ::  found
    integer                 ::  ikey
    real                    ::  r_curr

    ! --- begin -----------------------------
    
    ! pessimistic assumption ...
    found = .false. 
    
    ! loop over all key extensions:
    do ikey = 1, size(keys)
    
      ! try to read key; 
      ! provide default to return without error if key is not found:
      call ReadRc( rcfile, trim(key)//'.'//trim(keys(ikey)), r_curr, status, default=0.0 )
      if ( status < 0 ) then
        ! not found; try next
        cycle
      else if ( status == 0 ) then
        ! found and value read:
        found = .true.
        r = r_curr
      else
        ! error ...
        TRACEBACK; status=1; return
      end if
      
    end do   ! loop over keys
    
    ! not found ?
    if ( .not. found ) then
      ! default provided ?
      if ( present(default) ) then
        ! set to default:
        r = default
      else
        ! error ...
        write (gol,'("key(s) not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        do ikey = 1, size(keys)
          write (gol,'("  key    : ",a,".",a)') trim(key), trim(keys(ikey)); call goErr
        end do
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine ReadRcs_r


  ! ***


  subroutine ReadRc_r1( rcfile, key, r, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    real, intent(out)                           ::  r(:)
    integer, intent(out)                        ::  status
    
    real, intent(in), optional                  ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRc_r1'

    ! --- local -----------------------------

    character(len=buflen)     ::  buffer

    ! --- begin -----------------------------

    ! search key line in rcfile:
    call ReadRcItem( rcfile, key, buffer, status )
    if ( status < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        r = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else if ( status == 0 ) then
      ! key found; set value:
      read (buffer,*,iostat=status) r
      if ( status /= 0 ) then
        write (gol,'("while reading real :")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(buffer); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! some error ...
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine ReadRc_r1


  ! ================================================================
  ! ===
  ! === logical
  ! ===
  ! ================================================================
  
  
  subroutine ReadRc_l( rcfile, key, l, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    logical, intent(out)                        ::  l
    integer, intent(out)                        ::  status
    
    logical, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRc_l'

    ! --- local -----------------------------

    character(len=buflen)     ::  buffer

    ! --- begin -----------------------------

    ! search key line in rcfile:
    call ReadRcItem( rcfile, key, buffer, status )
    if ( status < 0 ) then
      ! not found; set to default or leave with warning:
      if ( present(default) ) then
        l = default
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else if ( status == 0 ) then
      ! key found; set value:
      read (buffer,*,iostat=status) l
      if ( status /= 0 ) then
        write (gol,'("while reading logical :")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        write (gol,'("  value  : ",a)') trim(buffer); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! some error ...
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine ReadRc_l


  ! ***

  
  subroutine ReadRcs_l( rcfile, key, keys, l, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    character(len=*), intent(in)                ::  keys(:)
    logical, intent(out)                        ::  l
    integer, intent(out)                        ::  status
    
    logical, intent(in), optional               ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRcs_l'

    ! --- local -----------------------------
    
    logical                 ::  found
    integer                 ::  ikey
    logical                 ::  l_curr

    ! --- begin -----------------------------
    
    ! pessimistic assumption ...
    found = .false. 
    
    ! loop over all key extensions:
    do ikey = 1, size(keys)
    
      ! try to read key; 
      ! provide default to return without error if key is not found:
      call ReadRc( rcfile, trim(key)//'.'//trim(keys(ikey)), l_curr, status, default=.false. )
      if ( status < 0 ) then
        ! not found; try next
        cycle
      else if ( status == 0 ) then
        ! found and value read:
        found = .true.
        l = l_curr
      else
        ! error ...
        TRACEBACK; status=1; return
      end if
      
    end do   ! loop over keys
    
    ! not found ?
    if ( .not. found ) then
      ! default provided ?
      if ( present(default) ) then
        ! set to default and leave with warning:
        l = default
        status = -1 ; return
     else
        ! error ...
        write (gol,'("key(s) not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        do ikey = 1, size(keys)
          write (gol,'("  key    : ",a,".",a)') trim(key), trim(keys(ikey)); call goErr
        end do
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine ReadRcs_l


  ! ================================================================
  ! ===
  ! === character string
  ! ===
  ! ================================================================


  subroutine ReadRc_s( rcfile, key, s, status, default )
  
    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    character(len=*), intent(out)               ::  s
    integer, intent(out)                        ::  status
    
    character(len=*), intent(in), optional      ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRc_s'

    ! --- local -----------------------------

    character(len=buflen)     ::  buffer

    ! --- begin -----------------------------

    ! search key line in rcfile:
    call ReadRcItem( rcfile, key, buffer, status )
    if ( status < 0 ) then
      ! not found; set to default or leave with error:
      if ( present(default) ) then
        s = trim(default)
        status = -1 ; return
      else
        write (gol,'("key not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        write (gol,'("  key    : ",a)') trim(key); call goErr
        TRACEBACK; status=1; return
      end if
    else if ( status == 0 ) then
      ! key found; set value:
      s = trim(buffer)
    else
      ! some error ...
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0

  end subroutine ReadRc_s


  ! ***

  
  subroutine ReadRcs_s( rcfile, key, keys, s, status, default )

    use GO_Print, only : gol, goErr

    ! --- in/out ----------------------------

    type(TrcFile), intent(in)                   ::  rcfile
    character(len=*), intent(in)                ::  key
    character(len=*), intent(in)                ::  keys(:)
    character(len=*), intent(out)               ::  s
    integer, intent(out)                        ::  status
    
    character(len=*), intent(in), optional      ::  default
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadRcs_l'

    ! --- local -----------------------------
    
    logical                 ::  found
    integer                 ::  ikey
    character(len=buflen)   ::  s_curr

    ! --- begin -----------------------------
    
    ! pessimistic assumption ...
    found = .false. 
    
    ! loop over all key extensions:
    do ikey = 1, size(keys)
    
      ! try to read key; 
      ! provide default to return without error if key is not found:
      call ReadRc( rcfile, trim(key)//'.'//trim(keys(ikey)), s_curr, status, default='-' )
      if ( status < 0 ) then
        ! not found; try next
        cycle
      else if ( status == 0 ) then
        ! found and value read:
        found = .true.
        s = trim(s_curr)
      else
        ! error ...
        TRACEBACK; status=1; return
      end if
      
    end do   ! loop over keys
    
    ! not found ?
    if ( .not. found ) then
      ! default provided ?
      if ( present(default) ) then
        ! set to default:
        s = default
        ! warning status
        status=-1; return
      else
        ! error ...
        write (gol,'("key(s) not found and no default specified ...")'); call goErr
        write (gol,'("  rcfile : ",a)') trim(rcfile%fname); call goErr
        do ikey = 1, size(keys)
          write (gol,'("  key    : ",a,".",a)') trim(key), trim(keys(ikey)); call goErr
        end do
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ok
    status = 0

  end subroutine ReadRcs_s


end module GO_Rc
