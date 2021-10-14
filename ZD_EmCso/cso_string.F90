!###############################################################################
!
! NAME
!   CSO_String - general objects for character strings
!
! PROCEDURES
!
!   call CSO_SplitLine( 'ab#cd', s1, '#', s2, status )
!
!     Splits a string like 'ab#cd' at the first '#', and returns
!     the leading part in s1, and the rest in s2.
!     One or both of s1 and s2 might be empty.
!
!   call CSO_ReadFromLine( line, x, status [,sep=','] [,default=value] )
!    
!     Splits the string "line" at the first komma
!     (or the character specified by the optional argument "sep"), 
!     fills the integer|real|logical|character variable "x" with the
!     leading part, and returns the remainder in "line".
!     If the leading part is empty, the default is returned if presend
!     or otherwise an error is raised.
!     
!   call CSO_SplitString( 'aa c', n, values, status [,sep=' '] )
!
!     Split the input string at white spaces and return fields:
!       n      : integer, number of fields extracted
!       values : character array to store fields;
!                error status returned if size or lengths are not sufficient
!   
!   call CSO_MatchValue( 'aa' , (/'aa ','bbb','c  ','ddd'/), ind, status [,quiet=.false.] )
!
!     Compare character value with values in character list,
!     return index of matching element.
!     Negative status if not found.
!     Print error messages unless quiet is true.
!
!   call CSO_MatchValues( 'aa c' , (/'aa ','bbb','c  ','ddd'/), n, values, inds, status [,quiet=.false.] )
!   call CSO_MatchValues( '0 1 2', 1, 10                      , n, values      , status [,quiet=.false.] )
!
!     Read the values from the input line and compare with a
!     list of character values or a range of integers. Output:
!       n       : integer, number of values in list (and output arrays)
!       values  : array with found values, same type as input list
!       inds    : integer indices in list with possible values
!
!
!   bb = 'default'
!   call CSO_VarValue( 'aa=1;bb=xyz;cc=U123', ';', 'bb', '=', bb, status )
!
!     Read value from a line with multiple <name><is><value> triples,
!     seperated by the specified character.
!     If multiple matching values are found, the last one is returned.
!     Return status:
!       <0  :  variable not found, val remains the same
!        0  :  variable(s) found, val reset; 
!       >0  :  error
!
!   s = goNum2Str( i [,fmt='(i6)'] )
!
!     Returns a 6-character string with the representation of the
!     integer value i in the first characters.
!     Use 
!       trim(gonum2str(i))
!     to obtain a string of smallest size.
!
!   s2 = CSO_UpCase( trim(s1) )
!   s2 = CSO_LoCase( trim(s1) )
!
!     Convert to upper or lower case
!
!   call CSO_Tab2Space( s )
!
!     Replaces each tab-character in s by a space.
!
!   call CSO_Replace( line, key, s, status )
!   call CSO_Replace( line, key, fmt, i, status )
!   call CSO_Replace( line, key, fmt, r, status )
!
!     Replace all instances of the key in the line by the 
!     character replacement s, or by a formatted integer or real value.
!
!   call CSO_Translate( line, chars, repl, status )
!
!     Replace all instances of the characters in 'chars' by 'repl' .
!     Example:
!       s = '2000-01-02 03:04:05'
!       print *, 'before translation : "'//trim(s)//'"'
!       call CSO_Translate( s, '-:', ' ', status )
!       print *, 'after  translation : "'//trim(s)//'"'
!     will show:
!       before translation : "2000-01-02 03:04:05"
!       after  translation : "2000 01 02 03 04 05"
!    
!   call CSO_CharToString( chars, line, status )
!     Replaces 1D or nD array of characters to string
!     or array of strings (one dim less than chars).
!
!### macro's #####################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call csoErr
!
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!#################################################################

module CSO_String

  use CSO_Logging, only : csol, csoPr, csoErr

  implicit none

  ! --- in/out -----------------------------

  private

  public  ::  CSO_SplitLine
  public  ::  CSO_ReadFromLine
  public  ::  CSO_SplitString
  public  ::  CSO_MatchValue
  public  ::  CSO_MatchValues
  public  ::  CSO_VarValue
  public  ::  CSO_Num2str
  public  ::  CSO_IsNum
  public  ::  CSO_UpCase, CSO_LoCase
  public  ::  CSO_WriteKeyNum
  public  ::  CSO_Tab2Space
  public  ::  CSO_Replace, CSO_Translate
  public  ::  CSO_CharToString
  public  ::  CSO_NtsTrim


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'CSO_String'

  
  ! --- interfaces -------------------------------------

  interface CSO_ReadFromLine
    module procedure ReadFromLine_i
    module procedure ReadFromLine_r
    module procedure ReadFromLine_l
    module procedure ReadFromLine_s
  end interface

  interface CSO_MatchValue
    module procedure MatchValue_s
    module procedure MatchValue_r
    module procedure MatchValue_i
  end interface

  interface CSO_SplitString
    module procedure CSO_SplitString_s
    module procedure CSO_SplitString_r
    module procedure CSO_SplitString_i
  end interface

  interface CSO_MatchValues
    module procedure MatchValues_s
    module procedure MatchValues_i
  end interface

  interface CSO_VarValue
    module procedure CSO_VarValue_s
    module procedure CSO_VarValue_i
    module procedure CSO_VarValue_r
    module procedure CSO_VarValue_l
  end interface

  interface CSO_Num2str
    module procedure num2str_i
    module procedure num2str_r
  end interface

  interface CSO_UpCase
    module procedure UpCase
  end interface

  interface CSO_LoCase
    module procedure LoCase
  end interface

  interface CSO_WriteKeyNum
    module procedure WriteKeyNum
  end interface

  interface CSO_Replace
    module procedure CSO_Replace_s
    module procedure CSO_Replace_i
  end interface

  interface CSO_CharToString
    module procedure CSO_CharToString_1d
    module procedure CSO_CharToString_2d
  end interface

contains


  !**********************************************************************


  subroutine CSO_SplitLine( line, s1, c, s2, status )

    ! --- in/out ----------------------------

    character(len=*), intent(in)      ::  line
    character(len=*), intent(out)     ::  s1
    character(len=1), intent(in)      ::  c
    character(len=*), intent(out)     ::  s2
    integer, intent(out)              ::  status

    ! --- local -----------------------------

    integer                     ::  l, pos
    character(len=len(line))    ::  s

    ! --- begin -----------------------------

    s = line
    l = len_trim(s)

    pos = scan(s,c)
    if ( (pos<1) .or. (pos>l) ) then
      ! s='abcd'  -> s1='abcd', s2=''
      !call AdjustLeft( s1, s(1:l) )
      s1 = AdjustL( s(1:l) )
      s2 = ''
    else if (pos==1) then
      ! s=',' or s=',abcd'  ->  s1='', s2='' or 'abcd'
      s1 = ''
      if (l==1) then
        ! s=','
        s2 = ''
      else
        !call AdjustLeft( s2, s(pos+1:l) )
        s2 = AdjustL( s(pos+1:l) )
      end if
    else
      ! s='ab,' or s='ab,cd'
      !call AdjustLeft( s1, s(1:pos-1) )
      s1 = AdjustL( s(1:pos-1) )
      if (pos==l) then
        ! s='ab,'
        s2 = ''
      else
        ! s='ab,cd'
        !call AdjustLeft( s2, s(pos+1:l) )
        s2 = AdjustL( s(pos+1:l) )
      end if
    end if
    
    ! ok
    status = 0

  end subroutine CSO_SplitLine


  ! ***


!  subroutine AdjustLeft( t, s )
!
!    ! --- in/out ----------------------
!
!    character(len=*), intent(out)   ::  t
!    character(len=*), intent(in)    ::  s
!
!    ! --- local -----------------------
!
!    integer     ::  is,ls, lt
!
!    ! --- local -----------------------
!
!    lt = len(t)
!
!    ls = len_trim(s)
!    if (ls==0) then
!      t = ''
!    else
!      is = 0
!      do
!        is = is + 1
!        if (s(is:is)/=' ') exit
!        if (is==ls) exit
!      end do
!      if (ls-is+1 > lt) then
!        print *, 'AdjustLeft : error : target is to small ', &
!           '(',lt,') to contain "'//s//'".'
!        stop
!      end if
!      t = s(is:ls)
!    end if
!
! end subroutine AdjustLeft


 ! *****************************************************


  subroutine ReadFromLine_i( s, i, status, sep, default )

    ! --- in/out --------------------------

    character(len=*), intent(inout)         ::  s
    integer, intent(inout)                  ::  i
    integer, intent(out)                    ::  status
    
    character(len=1), intent(in), optional  ::  sep
    integer, intent(in), optional           ::  default

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadFromLine_i'
    
    ! --- local ----------------------------

    character(len=len(s))     ::  s1, s2
    character(len=1)          ::  thesep

    ! --- begin ----------------------------
    
    ! default seperation character provided as argument:
    thesep = ','
    if (present(sep)) thesep = sep

    ! split at seperation character:
    call CSO_SplitLine( s, s1, thesep, s2, status )
    IF_ERROR_RETURN(status=1)

    ! empty leading part ?
    if ( len_trim(s1) == 0 ) then
      ! default provided ?
      if ( present(default) ) then
        i = default
      else
        write (csol,'("found empty leading part while no default specified:")'); call csoErr
        write (csol,'("  line   : `",a,"`")') trim(s); call csoErr
        write (csol,'("  sep    : `",a,"`")') trim(thesep); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! read from leading part:
      read (s1,*,iostat=status) i
      if ( status /= 0 ) then
        write (csol,'("while reading integer out of `",a,"`")') trim(s); call csoErr
        TRACEBACK; status=1; return
      end if
    end if
   
    ! return remainder:
    s = s2
    
    ! ok
    status = 0

  end subroutine ReadFromLine_i


  ! ***


  subroutine ReadFromLine_r( s, r, status, sep, default )

    ! --- in/out --------------------------

    character(len=*), intent(inout)         ::  s
    real, intent(out)                       ::  r
    integer, intent(out)                    ::  status

    character(len=1), intent(in), optional  ::  sep
    real, intent(in), optional              ::  default

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadFromLine_r'
    
    ! --- local ----------------------------

    character(len=len(s))     ::  s1, s2
    character(len=1)          ::  thesep

    ! --- begin ----------------------------

    ! default seperation character provided as argument:
    thesep = ','
    if (present(sep)) thesep = sep

    ! split at seperation character:
    call CSO_SplitLine( s, s1, thesep, s2, status )
    IF_ERROR_RETURN(status=1)

    ! empty leading part ?
    if ( len_trim(s1) == 0 ) then
      ! default provided ?
      if ( present(default) ) then
        r = default
      else
        write (csol,'("found empty leading part while no default specified:")'); call csoErr
        write (csol,'("  line   : `",a,"`")') trim(s); call csoErr
        write (csol,'("  sep    : `",a,"`")') trim(thesep); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! read from leading part:
      read (s1,*,iostat=status) r
      if ( status /= 0 ) then
        write (csol,'("error while reading real out `",a,"`")') trim(s); call csoErr
        TRACEBACK; status=1; return
      end if
    end if

    ! return remainder:
    s = s2

    ! ok
    status = 0

  end subroutine ReadFromLine_r


  ! ***


  subroutine ReadFromLine_l( s, l, status, sep, default )

    ! --- in/out --------------------------

    character(len=*), intent(inout)         ::  s
    logical, intent(out)                    ::  l
    integer, intent(out)                    ::  status

    character(len=1), intent(in), optional  ::  sep
    logical, intent(in), optional           ::  default

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadFromLine_l'
    
    ! --- local ----------------------------

    character(len=len(s))     ::  s1, s2
    character(len=1)          ::  thesep

    ! --- begin ----------------------------

    ! default seperation character provided as argument:
    thesep = ','
    if (present(sep)) thesep = sep

    ! split at seperation character:
    call CSO_SplitLine( s, s1, thesep, s2, status )
    IF_ERROR_RETURN(status=1)

    ! empty leading part ?
    if ( len_trim(s1) == 0 ) then
      ! default provided ?
      if ( present(default) ) then
        l = default
      else
        write (csol,'("found empty leading part while no default specified:")'); call csoErr
        write (csol,'("  line   : `",a,"`")') trim(s); call csoErr
        write (csol,'("  sep    : `",a,"`")') trim(thesep); call csoErr
        TRACEBACK; status=1; return
      end if
    else
      ! read from leading part:
      read (s1,*,iostat=status) l
      if ( status /= 0 ) then
        write (csol,'("while reading logical out `",a,"`")') trim(s); call csoErr
        TRACEBACK; status=1; return
      end if
    end if

    ! return remainder:
    s = s2

    ! ok
    status = 0

  end subroutine ReadFromLine_l


  ! ***


  subroutine ReadFromLine_s( s, ss, status, sep )

    ! --- in/out --------------------------

    character(len=*), intent(inout)         ::  s
    character(len=*), intent(out)           ::  ss
    integer, intent(out)                    ::  status

    character(len=1), intent(in), optional  ::  sep

    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ReadFromLine_s'
    
    ! --- local ----------------------------

    character(len=len(s))     ::  s1, s2
    character(len=1)          ::  thesep
    integer                   ::  l, ll

    ! --- begin ----------------------------

    ! default seperation character provided as argument:
    thesep = ','
    if (present(sep)) thesep = sep

    ! split at seperation character:
    call CSO_SplitLine( s, s1, thesep, s2, status )
    IF_ERROR_RETURN(status=1)
    
    ! check storage:
    l = len_trim(s1)
    ll = len(ss)
    if ( ll < l ) then
      write (csol,'("size of output string not sufficient:")'); call csoErr
      write (csol,'("  first part of input : ",a )') trim(s1) ; call csoErr
      write (csol,'("  output length       : ",i4)') ll       ; call csoErr
      TRACEBACK; status=1; return
    end if
    ! store:
    ss = trim(s1)

    ! return remainder:
    s = s2

    ! ok
    status = 0

  end subroutine ReadFromLine_s


  ! *****************************************************


  subroutine CSO_SplitString_s( line, n, values, status, sep )
  
    ! --- in/out --------------------------------
    
    character(len=*), intent(in)            ::  line
    integer, intent(out)                    ::  n
    character(len=*), intent(out)           ::  values(:)
    integer, intent(out)                    ::  status
    character(len=1), intent(in), optional  ::  sep
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/goSplitString_s'
    
    ! --- local ---------------------------------
    
    character(len=1)            ::  the_sep
    character(len=len(line))    ::  line_curr
    character(len=len(line))    ::  val

    ! --- begin ---------------------------------
    
    ! seperation character:
    the_sep = ' '
    if ( present(sep) ) the_sep = sep
    
    ! copy input:
    line_curr = line
    
    ! no parts extracted yet:
    n = 0

    ! loop until all elements in line_curr are processed:
    do
      ! empty ? then finished:
      if ( len_trim(line_curr) == 0 ) exit
      ! next number:
      n = n + 1
      ! storage problem ?
      if ( n > size(values) ) then
        write (csol,'("output array is too small:")'); call csoErr
        write (csol,'("  input line   : ",a )') trim(line); call csoErr
        write (csol,'("  size(values) : ",i4)') size(values); call csoErr
        TRACEBACK; status=1; return
      end if
      ! extract leading name:
      call CSO_ReadFromLine( line_curr, val, status, sep=the_sep )
      IF_NOT_OK_RETURN(status=1)
      ! store value in output:
      values(n) = adjustl(val)
    end do
    
    ! ok
    status = 0
    
  end subroutine CSO_SplitString_s
  
  
  ! ***


  subroutine CSO_SplitString_r( line, n, values, status, sep )
  
    ! --- in/out --------------------------------
    
    character(len=*), intent(in)            ::  line
    integer, intent(out)                    ::  n
    real, intent(out)                       ::  values(:)
    integer, intent(out)                    ::  status
    character(len=1), intent(in), optional  ::  sep
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/goSplitString_r'
    
    ! --- local ---------------------------------
    
    character(len=1)            ::  the_sep
    character(len=len(line))    ::  line_curr
    real                        ::  val

    ! --- begin ---------------------------------
    
    ! seperation character:
    the_sep = ' '
    if ( present(sep) ) the_sep = sep
    
    ! copy input:
    line_curr = line
    
    ! no parts extracted yet:
    n = 0

    ! loop until all elements in line_curr are processed:
    do
      ! empty ? then finished:
      if ( len_trim(line_curr) == 0 ) exit
      ! next number:
      n = n + 1
      ! storage problem ?
      if ( n > size(values) ) then
        write (csol,'("output array is too small:")'); call csoErr
        write (csol,'("  input line   : ",a )') trim(line); call csoErr
        write (csol,'("  size(values) : ",i4)') size(values); call csoErr
        TRACEBACK; status=1; return
      end if
      ! extract leading name:
      call CSO_ReadFromLine( line_curr, val, status, sep=the_sep )
      IF_NOT_OK_RETURN(status=1)
      ! store value in output:
      values(n) = val
    end do
    
    ! ok
    status = 0
    
  end subroutine CSO_SplitString_r

  ! ***


  subroutine CSO_SplitString_i( line, n, values, status, sep )
  
    ! --- in/out --------------------------------
    
    character(len=*), intent(in)            ::  line
    integer, intent(out)                    ::  n
    integer, intent(out)                    ::  values(:)
    integer, intent(out)                    ::  status
    character(len=1), intent(in), optional  ::  sep
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/goSplitString_r'
    
    ! --- local ---------------------------------
    
    character(len=1)            ::  the_sep
    character(len=len(line))    ::  line_curr
    integer                     ::  val

    ! --- begin ---------------------------------
    
    ! seperation character:
    the_sep = ' '
    if ( present(sep) ) the_sep = sep
    
    ! copy input:
    line_curr = line
    
    ! no parts extracted yet:
    n = 0

    ! loop until all elements in line_curr are processed:
    do
      ! empty ? then finished:
      if ( len_trim(line_curr) == 0 ) exit
      ! next number:
      n = n + 1
      ! storage problem ?
      if ( n > size(values) ) then
        write (csol,'("output array is too small:")'); call csoErr
        write (csol,'("  input line   : ",a )') trim(line); call csoErr
        write (csol,'("  size(values) : ",i4)') size(values); call csoErr
        TRACEBACK; status=1; return
      end if
      ! extract leading name:
      call CSO_ReadFromLine( line_curr, val, status, sep=the_sep )
      IF_NOT_OK_RETURN(status=1)
      ! store value in output:
      values(n) = val
    end do
    
    ! ok
    status = 0
    
  end subroutine CSO_SplitString_i
  
  
  ! *******************************************************************
  

  subroutine MatchValue_s( val, list, ind, status, quiet )
  
    ! --- in/out --------------------------------
    
    character(len=*), intent(in)          ::  val
    character(len=*), intent(in)          ::  list(:)
    integer, intent(out)                  ::  ind
    integer, intent(out)                  ::  status
    
    logical, optional                     ::  quiet
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MatchValue_s'
    
    ! --- local ---------------------------------
    
    integer                     ::  nlist
    integer                     ::  i
    logical                     ::  verbose

    ! --- begin ---------------------------------
    
    ! shut up ?
    verbose = .true.
    if ( present(quiet) ) verbose = .not. quiet
    
    ! number of items in value list:
    nlist = size(list)
    
    ! search for this name in the global list:
    ind = -1
    do i = 1, nlist
      ! case indendent match ?
      if ( CSO_UpCase(trim(val)) == CSO_UpCase(trim(list(i))) ) then
        ! store index:
        ind = i
        ! do not search any further:
        exit
      end if
    end do
    ! not found ?
    if ( ind < 0 ) then
      if ( verbose ) then
        write (csol,'("name not supported:")'); call csoErr
        write (csol,'("  value           : ",a )') trim(val); call csoErr
        write (csol,'("  possible values : ")'); call csoErr
        do i = 1, nlist
          if ( len_trim(list(i)) > 0 ) then
            write (csol,'("    ",i4," ",a)') i, trim(list(i)); call csoErr
          end if
        end do
        TRACEBACK
      end if
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine MatchValue_s
  
  
  ! ***
  

  subroutine MatchValue_r( val, list, ind, status, quiet )
  
    ! --- in/out --------------------------------
    
    real, intent(in)                      ::  val
    real, intent(in)                      ::  list(:)
    integer, intent(out)                  ::  ind
    integer, intent(out)                  ::  status
    
    logical, optional                     ::  quiet
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MatchValue_r'
    
    ! --- local ---------------------------------
    
    integer                     ::  nlist
    integer                     ::  i
    logical                     ::  verbose

    ! --- begin ---------------------------------
    
    ! shut up ?
    verbose = .true.
    if ( present(quiet) ) verbose = .not. quiet
    
    ! number of items in value list:
    nlist = size(list)
    
    ! search for this name in the global list:
    ind = -1
    do i = 1, nlist
      ! case indendent match ?
      if ( val == list(i) ) then
        ! store index:
        ind = i
        ! do not search any further:
        exit
      end if
    end do
    ! not found ?
    if ( ind < 0 ) then
      if ( verbose ) then
        write (csol,'("name not supported:")'); call csoErr
        write (csol,'("  value           : ",i6 )') val; call csoErr
        write (csol,'("  possible values : ")'); call csoErr
        do i = 1, nlist
          write (csol,'("    ",i4," ",i6)') i, list(i); call csoErr
        end do
        TRACEBACK
      end if
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine MatchValue_r

  
  ! ***
  

  subroutine MatchValue_i( val, list, ind, status, quiet )
  
    ! --- in/out --------------------------------
    
    integer, intent(in)                   ::  val
    integer, intent(in)                   ::  list(:)
    integer, intent(out)                  ::  ind
    integer, intent(out)                  ::  status
    
    logical, optional                     ::  quiet
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MatchValue_i'
    
    ! --- local ---------------------------------
    
    integer                     ::  nlist
    integer                     ::  i
    logical                     ::  verbose

    ! --- begin ---------------------------------
    
    ! shut up ?
    verbose = .true.
    if ( present(quiet) ) verbose = .not. quiet
    
    ! number of items in value list:
    nlist = size(list)
    
    ! search for this name in the global list:
    ind = -1
    do i = 1, nlist
      ! case indendent match ?
      if ( val == list(i) ) then
        ! store index:
        ind = i
        ! do not search any further:
        exit
      end if
    end do
    ! not found ?
    if ( ind < 0 ) then
      if ( verbose ) then
        write (csol,'("name not supported:")'); call csoErr
        write (csol,'("  value           : ",i6 )') val; call csoErr
        write (csol,'("  possible values : ")'); call csoErr
        do i = 1, nlist
          write (csol,'("    ",i4," ",i6)') i, list(i); call csoErr
        end do
        TRACEBACK
      end if
      status=-1; return
    end if

    ! ok
    status = 0
    
  end subroutine MatchValue_i


  ! ***
  

  !
  ! Split the line into parts,
  ! try to match each part with the elements of a value list.
  !
  ! A warning status <0 is returned if a list part could not be matched;
  ! an informative messages is displayed too, unless quiet=.true.
  !
  ! If skipext is set (for example to '.'), then elements of the line
  ! (for example 'vector.u') are compared without the extension after '.',
  ! thus the values are scanned for 'vector' in this example.
  !

  subroutine MatchValues_s( line, list, &
                              n, values, inds, &
                              status, quiet, skipext )
  
    ! --- in/out --------------------------------
    
    character(len=*), intent(in)          ::  line
    character(len=*), intent(in)          ::  list(:)
    integer, intent(out)                  ::  n
    character(len=*), intent(out)         ::  values(:)
    integer, intent(out)                  ::  inds(:)
    integer, intent(out)                  ::  status
    
    logical, optional                     ::  quiet
    character(len=1), optional            ::  skipext
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MatchValues_s'
    
    ! --- local ---------------------------------
    
    integer                     ::  nlist
    character(len=len(line))    ::  line_curr
    character(len=32)           ::  val
    integer                     ::  ind
    logical                     ::  verbose
    integer                     ::  k

    ! --- begin ---------------------------------
    
    ! shut up ?
    verbose = .true.
    if ( present(quiet) ) verbose = .not. quiet
    
    ! nuber of items in value list:
    nlist = size(list)
    
    ! copy input:
    line_curr = line
    
    ! no matching list yet:
    n = 0

    ! loop until all elements in line_curr are processed:
    do
      ! empty ? then finished:
      if ( len_trim(line_curr) == 0 ) exit
      ! next number:
      n = n + 1
      ! storage problem ?
      if ( (n > size(values)) .or. (n > size(inds)) ) then
        write (csol,'("output array is too small:")'); call csoErr
        write (csol,'("  input line  : ",a )') trim(line); call csoErr
        write (csol,'("  size(values) : ",i4)') size(values); call csoErr
        write (csol,'("  size(inds  ) : ",i4)') size(inds  ); call csoErr
        TRACEBACK; status=1; return
      end if
      ! extract leading name:
      call CSO_ReadFromLine( line_curr, val, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      ! store value in output:
      values(n) = val
      ! skip extension for comparison with list?
      if ( present(skipext) ) then
        ! location of extension seperation character,
        ! or zero if not present:
        k = index( trim(val), skipext )
        ! present?
        if ( k > 0 ) then
          ! first part only:
          val = val(1:k-1)
        end if
      end if
      ! search for this name in the global list:
      call CSO_MatchValue( val, list, ind, status, quiet )
      ! not found ?
      if ( status /= 0 ) then
        if ( verbose ) then
          write (csol,'("unable to match value with list:")'); call csoErr
          write (csol,'("  line            : ",a )') trim(line); call csoErr
          write (csol,'("  line element    : ",i3)') n; call csoErr
          write (csol,'("  line value      : ",a )') trim(val); call csoErr
          TRACEBACK
        end if
        status=1; return
      end if
      ! store:
      inds(n) = ind
    end do
    
    ! empty ?
    if ( n < 1 ) then
      write (csol,'("no values extracted from line :")'); call csoErr
      write (csol,'("  ",a)') trim(line); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine MatchValues_s
  
  
  ! ***
  

  subroutine MatchValues_i( line, i1, i2, &
                              n, values, &
                              status, quiet )
  
    ! --- in/out --------------------------------
    
    character(len=*), intent(in)          ::  line
    integer, intent(in)                   ::  i1, i2
    integer, intent(out)                  ::  n
    integer, intent(out)                  ::  values(:)
    integer, intent(out)                  ::  status
    
    logical, optional                     ::  quiet
    
    ! --- const ----------------------------
    
    character(len=*), parameter   ::  rname = mname//'/MatchValues_i'
    
    ! --- local ---------------------------------
    
    character(len=len(line))    ::  line_curr
    integer                     ::  val
    logical                     ::  verbose

    ! --- begin ---------------------------------
    
    ! shut up ?
    verbose = .true.
    if ( present(quiet) ) verbose = .not. quiet
    
    ! copy input:
    line_curr = line
    
    ! no matching list yet:
    n = 0

    ! loop until all elements in line_curr are processed:
    do
      ! empty ? then finished:
      if ( len_trim(line_curr) == 0 ) exit
      ! next number:
      n = n + 1
      ! storage problem ?
      if ( n > size(values) ) then
        write (csol,'("output arrays are too small:")'); call csoErr
        write (csol,'("  input line  : ",a )') trim(line); call csoErr
        write (csol,'("  size(values) : ",i4)') size(values); call csoErr
        TRACEBACK; status=1; return
      end if
      ! extract leading name:
      call CSO_ReadFromLine( line_curr, val, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
      ! store value in output:
      values(n) = val
      ! out of range ?
      if ( (val < i1) .or. (val > i2) ) then
        if ( verbose ) then
          write (csol,'("value not in range:")'); call csoErr
          write (csol,'("  list            : ",a )') trim(line); call csoErr
          write (csol,'("  list element    : ",i3)') n; call csoErr
          write (csol,'("  list value      : ",i3)') val; call csoErr
          write (csol,'("  possible range  : ",i3," .. ",i3)') i1, i2; call csoErr
          TRACEBACK
        end if
        status=1; return
      end if
    end do
    
    ! empty ?
    if ( n < 1 ) then
      write (csol,'("no values extracted from line :")'); call csoErr
      write (csol,'("  ",a)') trim(line); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine MatchValues_i
  

  ! *****************************************************
  
  !
  !  Read value from line:
  !
  !    bb = 'default'
  !    call CSO_VarValue( 'aa=1;bb=xyz;cc=U123', ';', 'bb', '=', bb, status )
  !
  !  Return status:
  !    <0  :  variable not found, val remains the same
  !             print error message only if verbose is .true.
  !     0  :  variable found, val reset; 
  !             for multiple matches, last value is returned
  !    >0  :  error
  !
  
  subroutine CSO_VarValue_s( line, sep, var, is, val, status, verbose )
  
    use CSO_Logging, only : csol, csoPr, csoErr

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)     ::  line
    character(len=1), intent(in)     ::  sep
    character(len=*), intent(in)     ::  var
    character(len=1), intent(in)     ::  is
    character(len=*), intent(inout)  ::  val
    integer, intent(out)             ::  status
    
    logical, intent(in), optional    ::  verbose
    
    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goVarValue_s'
    
    ! --- local ----------------------------------
    
    logical                     ::  verb
    character(len=len(line))    ::  line2
    character(len=len(line))    ::  varval
    character(len=64)           ::  var2
    character(len=256)          ::  val2
    logical                     ::  found
    
    ! --- begin ----------------------------------
    
    ! shout?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! copy of input line:
    line2 = line
    
    ! by default not found:
    found = .false.

    ! loop over var=val keys :
    do
      ! no keys left ? then leave
      if ( len_trim(line2) == 0 ) exit
      ! remove leading var=value from line2 :
      call CSO_ReadFromLine( line2, varval, status, sep=sep )
      IF_ERROR_RETURN(status=1)
      ! split in var and value:
      call CSO_SplitLine( varval, var2, is, val2, status )
      IF_ERROR_RETURN(status=1)
      ! keys match ?
      if ( trim(var2) == trim(var) ) then
        ! store in output (might overwrite previously stored value):
        val = trim(val2)
        ! set flag:
        found = .true.
        ! continue, maybye the key is found twice ...
      end if
    end do

    ! set return status:
    if ( found ) then
      ! ok
      status = 0
    else
      ! shout?
      if ( verb ) then
        write (csol,'("no value defined for `",a,"` in (sep `",a,"`, assign `",a,"`) : ")') &
                           trim(var), sep, is; call csoErr
        write (csol,'("  ",a)') trim(line); call csoErr
      end if
      ! warning status
      status = -1
    end if
    
  end subroutine CSO_VarValue_s


  ! ***
      

  subroutine CSO_VarValue_i( line, sep, var, is, val, status, verbose )
  
    use CSO_Logging, only : csol, csoPr, csoErr

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)     ::  line
    character(len=1), intent(in)     ::  sep
    character(len=*), intent(in)     ::  var
    character(len=1), intent(in)     ::  is
    integer, intent(inout)           ::  val
    integer, intent(out)             ::  status
    
    logical, intent(in), optional    ::  verbose
    
    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goVarValue_i'
    
    ! --- local ----------------------------------
    
    logical                     ::  verb
    character(len=len(line))    ::  line2
    character(len=len(line))    ::  varval
    character(len=64)           ::  var2
    character(len=256)          ::  val2
    logical                     ::  found
    
    ! --- begin ----------------------------------
    
    ! shout?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! copy of input line:
    line2 = line
    
    ! by default not found:
    found = .false.

    ! loop over var=val keys :
    do
      ! no keys left ? then leave
      if ( len_trim(line2) == 0 ) exit
      ! remove leading var=value from line2 :
      call CSO_ReadFromLine( line2, varval, status, sep=sep )
      IF_ERROR_RETURN(status=1)
      ! split in var and value:
      call CSO_SplitLine( varval, var2, is, val2, status )
      IF_ERROR_RETURN(status=1)
      ! keys match ?
      if ( trim(var2) == trim(var) ) then
        ! store in output (might overwrite previously stored value):
        read (val2,'(i6)') val
        ! set flag:
        found = .true.
        ! continue, maybye the key is found twice ...
      end if
    end do

    ! set return status:
    if ( found ) then
      status = 0
    else
      ! shout?
      if ( verb ) then
        write (csol,'("no value defined for `",a,"` in (sep `",a,"`, assign `",a,"`) : ",a)') &
                           trim(var), sep, is; call csoErr
        write (csol,'("  ",a)') trim(line); call csoErr
      end if
      ! warning status
      status = -1
    end if
    
  end subroutine CSO_VarValue_i
    

  ! ***
      

  subroutine CSO_VarValue_r( line, sep, var, is, val, status, verbose )
  
    use CSO_Logging, only : csol, csoPr, csoErr

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)     ::  line
    character(len=1), intent(in)     ::  sep
    character(len=*), intent(in)     ::  var
    character(len=1), intent(in)     ::  is
    real, intent(inout)              ::  val
    integer, intent(out)             ::  status
    
    logical, intent(in), optional    ::  verbose
    
    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goVarValue_r'
    
    ! --- local ----------------------------------
    
    logical                     ::  verb
    character(len=len(line))    ::  line2
    character(len=len(line))    ::  varval
    character(len=64)           ::  var2
    character(len=256)          ::  val2
    logical                     ::  found
    
    ! --- begin ----------------------------------
    
    ! shout?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! copy of input line:
    line2 = line

    ! by default not found:
    found = .false.

    ! loop over var=val keys :
    do
      ! no keys left ? then leave
      if ( len_trim(line2) == 0 ) exit
      ! remove leading var=value from line2 :
      call CSO_ReadFromLine( line2, varval, status, sep=sep )
      IF_ERROR_RETURN(status=1)
      ! split in var and value:
      call CSO_SplitLine( varval, var2, is, val2, status )
      IF_ERROR_RETURN(status=1)
      ! keys match ?
      if ( trim(var2) == trim(var) ) then
        ! store in output (might overwrite previously stored value):
        read (val2,*) val
        ! set flag:
        found = .true.
        ! continue, maybye the key is found twice ...
      end if
    end do

    ! set return status:
    if ( found ) then
      status = 0
    else
      ! shout?
      if ( verb ) then
        write (csol,'("no value defined for `",a,"` in (sep `",a,"`, assign `",a,"`) : ",a)') &
                           trim(var), sep, is; call csoErr
        write (csol,'("  ",a)') trim(line); call csoErr
      end if
      ! warning status
      status = -1
    end if
    
  end subroutine CSO_VarValue_r


  ! ***
      

  subroutine CSO_VarValue_l( line, sep, var, is, val, status, verbose )
  
    use CSO_Logging, only : csol, csoPr, csoErr

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)     ::  line
    character(len=1), intent(in)     ::  sep
    character(len=*), intent(in)     ::  var
    character(len=1), intent(in)     ::  is
    logical, intent(inout)           ::  val
    integer, intent(out)             ::  status
    
    logical, intent(in), optional    ::  verbose
    
    ! --- const ------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goVarValue_l'
    
    ! --- local ----------------------------------
    
    logical                     ::  verb
    character(len=len(line))    ::  line2
    character(len=len(line))    ::  varval
    character(len=64)           ::  var2
    character(len=256)          ::  val2
    logical                     ::  found
    
    ! --- begin ----------------------------------
    
    ! shout?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! copy of input line:
    line2 = line
    
    ! by default not found:
    found = .false.

    ! loop over var=val keys :
    do
      ! no keys left ? then leave
      if ( len_trim(line2) == 0 ) exit
      ! remove leading var=value from line2 :
      call CSO_ReadFromLine( line2, varval, status, sep=sep )
      IF_ERROR_RETURN(status=1)
      ! split in var and value:
      call CSO_SplitLine( varval, var2, is, val2, status )
      IF_ERROR_RETURN(status=1)
      ! keys match ?
      if ( trim(var2) == trim(var) ) then
        ! store in output (might overwrite previously stored value):
        read (val2,'(l1)') val
        ! set flag:
        found = .true.
        ! continue, maybye the key is found twice ...
      end if
    end do

    ! set return status:
    if ( found ) then
      status = 0
    else
      ! shout?
      if ( verb ) then
        write (csol,'("no value defined for `",a,"` in (sep `",a,"`, assign `",a,"`) : ",a)') &
                           trim(var), sep, is; call csoErr
        write (csol,'("  ",a)') trim(line); call csoErr
      end if
      ! warning status
      status = -1
    end if
    
  end subroutine CSO_VarValue_l
    

  ! *****************************************************

  !---
  ! NAME
  !   gonum2str - prints number into character string
  !
  ! INTERFACE
  !   character(len=20) function CSO_num2str( x, fmt )
  !     integer [or real]  , intent(in)            ::  x
  !     character(len=*), intent(in), optional  ::  fmt
  !
  ! ARGUMENTS
  !   x
  !     Number to be converted.
  !   fmt
  !     Optional format, following the formats provided
  !     to the 'write' command.
  !     Default values:
  !
  !       type  x             fmt         example  (- is space)
  !       ------------------  ----------  ---------------------
  !       integer             '(i6)'      123---
  !
  ! CHANGES
  !   01/09/1999  Arjo Segers
  !---

  character(len=6) function num2str_i( i, fmt )

    ! --- in/out ----------------------

    integer, intent(in)                     ::  i
    character(len=*), intent(in), optional  ::  fmt

    ! --- local -----------------------

    character(len=6)   ::  s

    ! --- begin -----------------------

    if (present(fmt)) then
      write (s,fmt=fmt) i
    else
      write (s,'(i6)') i
    end if
    num2str_i=adjustl(s)

  end function num2str_i

  character(len=12) function num2str_r( r, fmt )

    ! --- in/out ----------------------

    real, intent(in)                     ::  r
    character(len=*), intent(in), optional  ::  fmt

    ! --- local -----------------------

    character(len=12)   ::  s

    ! --- begin -----------------------

    if (present(fmt)) then
      write (s,fmt=fmt) r
    else
      write (s,'(g10.3)') r
    end if
    num2str_r=adjustl(s)

  end function num2str_r
  
  
  ! ***
  
  function CSO_IsNum( s )
  
    ! --- in/out -----------------

    character(len=*), intent(in)    ::  s
    logical                         ::  CSO_IsNum

    ! --- local ------------------

    integer         ::  i

    ! --- begin ------------------

    ! init result:
    CSO_IsNum = .true.
    ! loop over characters:
    do i = 1, len_trim(s)
      ! no numeric value if strange character:
      if ( index('0123456789.eEdD+-',s(i:i)) < 1 ) then
        ! reset:
        CSO_IsNum = .false.
        ! leave:
        exit
      end if  ! strange char
    end do ! chars
    
  end function CSO_IsNum


  ! *** UpCase, LoCase ***


  function UpCase( s )

    ! --- in/out -----------------

    character(len=*), intent(in)    ::  s
    character(len=len(s))           ::  UpCase

    ! --- local ------------------

    integer         ::  i

    ! --- begin ------------------

    do i = 1, len_trim(s)
      select case (s(i:i))
        case ('a') ;  UpCase(i:i) = 'A'
        case ('b') ;  UpCase(i:i) = 'B'
        case ('c') ;  UpCase(i:i) = 'C'
        case ('d') ;  UpCase(i:i) = 'D'
        case ('e') ;  UpCase(i:i) = 'E'
        case ('f') ;  UpCase(i:i) = 'F'
        case ('g') ;  UpCase(i:i) = 'G'
        case ('h') ;  UpCase(i:i) = 'H'
        case ('i') ;  UpCase(i:i) = 'I'
        case ('j') ;  UpCase(i:i) = 'J'
        case ('k') ;  UpCase(i:i) = 'K'
        case ('l') ;  UpCase(i:i) = 'L'
        case ('m') ;  UpCase(i:i) = 'M'
        case ('n') ;  UpCase(i:i) = 'N'
        case ('o') ;  UpCase(i:i) = 'O'
        case ('p') ;  UpCase(i:i) = 'P'
        case ('q') ;  UpCase(i:i) = 'Q'
        case ('r') ;  UpCase(i:i) = 'R'
        case ('s') ;  UpCase(i:i) = 'S'
        case ('t') ;  UpCase(i:i) = 'T'
        case ('u') ;  UpCase(i:i) = 'U'
        case ('v') ;  UpCase(i:i) = 'V'
        case ('w') ;  UpCase(i:i) = 'W'
        case ('x') ;  UpCase(i:i) = 'X'
        case ('y') ;  UpCase(i:i) = 'Y'
        case ('z') ;  UpCase(i:i) = 'Z'
        case default
          UpCase(i:i) = s(i:i)
      end select
    end do

  end function UpCase


  ! ***


  function LoCase( s )

    ! --- in/out -----------------

    character(len=*), intent(in)   ::  s
    character(len=len(s))          ::  LoCase

    ! --- local ------------------

    integer         ::  i

    ! --- begin ------------------

    do i = 1, len_trim(s)
      select case (s(i:i))
        case ('A') ;  LoCase(i:i) = 'a'
        case ('B') ;  LoCase(i:i) = 'b'
        case ('C') ;  LoCase(i:i) = 'c'
        case ('D') ;  LoCase(i:i) = 'd'
        case ('E') ;  LoCase(i:i) = 'e'
        case ('F') ;  LoCase(i:i) = 'f'
        case ('G') ;  LoCase(i:i) = 'g'
        case ('H') ;  LoCase(i:i) = 'h'
        case ('I') ;  LoCase(i:i) = 'i'
        case ('J') ;  LoCase(i:i) = 'j'
        case ('K') ;  LoCase(i:i) = 'k'
        case ('L') ;  LoCase(i:i) = 'l'
        case ('M') ;  LoCase(i:i) = 'm'
        case ('N') ;  LoCase(i:i) = 'n'
        case ('O') ;  LoCase(i:i) = 'o'
        case ('P') ;  LoCase(i:i) = 'p'
        case ('Q') ;  LoCase(i:i) = 'q'
        case ('R') ;  LoCase(i:i) = 'r'
        case ('S') ;  LoCase(i:i) = 's'
        case ('T') ;  LoCase(i:i) = 't'
        case ('U') ;  LoCase(i:i) = 'u'
        case ('V') ;  LoCase(i:i) = 'v'
        case ('W') ;  LoCase(i:i) = 'w'
        case ('X') ;  LoCase(i:i) = 'x'
        case ('Y') ;  LoCase(i:i) = 'y'
        case ('Z') ;  LoCase(i:i) = 'z'
        case default
          LoCase(i:i) = s(i:i)
      end select
    end do

  end function LoCase
  
  
  ! ***
  
  
  subroutine WriteKeyNum( res, key, num )
  
    ! --- in/out ------------------------------
    
    character(len=*), intent(out)    ::  res
    character(len=*), intent(in)     ::  key
    integer, intent(in)              ::  num
    
    ! --- local -------------------------------
    
    integer    ::  anum
    
    ! --- begin -------------------------------
    
    anum = abs(num)
    
    if ( anum <= 9 ) then
      write (res,'(a,i1)') trim(key), anum
    else if ( anum <= 99 ) then
      write (res,'(a,i2)') trim(key), anum
    else if ( anum <= 999 ) then
      write (res,'(a,i3)') trim(key), anum
    else if ( anum <= 9999 ) then
      write (res,'(a,i4)') trim(key), anum
    else if ( anum <= 99999 ) then
      write (res,'(a,i5)') trim(key), anum
    else
      write (res,'(a,i6)') trim(key), anum
    end if
    
  end subroutine WriteKeyNum
    

  ! ***
      

  subroutine CSO_Tab2Space( s )

    ! --- in/out -----------------

    character(len=*), intent(inout)    ::  s

    ! --- local ------------------

    integer         ::  pos

    ! --- begin ------------------

    do 
      pos = scan(s,char(9))
      if ( pos == 0 ) exit
      s(pos:pos) = ' '
    end do

  end subroutine CSO_Tab2Space


  ! ***
      
      
  subroutine CSO_Replace_s( s, key, repl, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(inout)     ::  s
    character(len=*), intent(in)        ::  key
    character(len=*), intent(in)        ::  repl
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goReplace_s'
    
    ! --- local ----------------------------------
    
    integer                 ::  n_in, l_in
    integer                 ::  n_out
    integer                 ::  ind
    character(len=len(s))   ::  s_in
    
    ! --- begin ----------------------------------
    
    ! copy input:
    s_in = s
    
    ! empty target:
    s = ''
    n_out = 0
    
    ! number of characters in s:
    l_in = len_trim(s_in)
    
    ! characters copied from s_in to s:
    n_in = 0
    
    ! loop over all matches of key:
    do
      !print *, '---- n_in, l_in : ',n_in, l_in
      ! past end ?
      if ( n_in > l_in ) exit
      !print *, '  -- remaining : "'//s_in(n_in+1:l_in)//'" , key : "'//key//'"'
      ! search key in remaining part of input :
      if ( len(key) < 1 ) then
        ind = 0
      else
        ind = index(s_in(n_in+1:l_in),key)
      end if
      !print *, '  -- index     : ', ind
      ! not found ?
      if ( ind < 1 ) then
        ! add remaining part:
        s = s(1:n_out)//s_in(n_in+1:l_in)
        n_out = n_out + l_in-n_in
        ! leave:
        exit
      end if
      ! add first part:
      s = s(1:n_out)//s_in(n_in+1:n_in+ind-1)
      n_out = n_out + ind-1
      n_in  = n_in + ind-1
      !print *, '  -- out       : ',n_out,'"'//s(1:n_out)//'"'
      ! add replacement:
      s = s(1:n_out)//repl
      n_out = n_out + len(repl)
      n_in  = n_in + len(key)
      !print *, '  -- out       : ',n_out,'"'//s(1:n_out)//'"'
    end do
      
    ! ok
    status = 0
    
  end subroutine CSO_Replace_s


  ! ***
      
      
  subroutine CSO_Replace_i( s, key, fmt, i, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(inout)     ::  s
    character(len=*), intent(in)        ::  key
    character(len=*), intent(in)        ::  fmt
    integer, intent(in)                 ::  i
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goReplace_i'
    
    ! --- local ----------------------------------
    
    character(len=8)    ::  repl    
    
    ! --- begin ----------------------------------
    
    ! fill replacement:
    write (repl,fmt) i
    
    ! replace value:
    call CSO_Replace( s, key, trim(repl), status )
    IF_NOT_OK_RETURN(status=1)    
      
    ! ok
    status = 0
    
  end subroutine CSO_Replace_i
    

  ! ***
      
  ! replace in s all characters in 'chars' by the value in 'repl' .
  
  subroutine CSO_Translate( s, chars, repl, status )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(inout)     ::  s
    character(len=*), intent(in)        ::  chars
    character(len=*), intent(in)        ::  repl
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goTranslate'
    
    ! --- local ----------------------------------
    
    integer           ::  i
    
    ! --- begin ----------------------------------
    
    ! loop over characters to be replaced:
    do i = 1, len(chars)
      ! replace character:
      call CSO_Replace( s, chars(i:i), repl, status )
      IF_NOT_OK_RETURN(status=1)
    end do
      
    ! ok
    status = 0
    
  end subroutine CSO_Translate


  ! ***
  
  
  subroutine CSO_CharToString_1d( chars, line, status )
  
    ! --- in/out ---------------------------------
    
    character(len=1), intent(in)        ::  chars(:)
    character(len=*), intent(out)       ::  line
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goCharToString_1d'
    
    ! --- local ----------------------------------
    
    integer           ::  i
    
    ! --- begin ----------------------------------
    
    ! init result:
    line = ''
    ! loop over character:
    do i = 1, size(chars)
      ! stop at first null-character:
      if ( iachar(chars(i) ) == 0 ) exit
      ! check ..
      if ( i > len(line) ) then
        write (csol,'("found more characters than can be stored in line : ",i6)') len(line); call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      line(i:i) = chars(i)
    end do
    
    ! ok
    status = 0
    
  end subroutine CSO_CharToString_1d
  
  ! *
  
  subroutine CSO_CharToString_2d( chars, lines, status )
  
    ! --- in/out ---------------------------------
    
    character(len=1), intent(in)        ::  chars(:,:)
    character(len=*), intent(out)       ::  lines(:)
    integer, intent(out)                ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goCharToString_2d'
    
    ! --- local ----------------------------------
    
    integer           ::  j
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(chars,2) /= size(lines) ) then
      write (csol,'("shape of character array (",2i6,") differs from lines (",i6,")")') shape(chars), shape(lines); call csoErr
      TRACEBACK; status=1; return
    end if
  
    ! loop over lines:
    do j = 1, size(lines)
      ! apply to line:
      call CSO_CharToString_1d( chars(:,j), lines(j), status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! ok
    status = 0
    
  end subroutine CSO_CharToString_2d


  ! ********************************************************************
  ! ***
  ! *** null-terminated-string
  ! ***
  ! ********************************************************************
      

  ! truncate null-terminated-string at first null character: 
  
  function CSO_NtsTrim( s )
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)        ::  s
    character(len=len(s))               ::  CSO_NtsTrim
    
    ! --- local ----------------------------------
    
    integer           ::  n
    
    ! --- begin ----------------------------------
  
    ! key is a zero-terminated string, get the length:
    n = index(s,char(0))
    ! present?
    if ( n == 1 ) then
      ! empty:
      CSO_NtsTrim = ''
    else if ( n > 1 ) then
      ! extract:
      CSO_NtsTrim = s(1:n-1)
    else
      ! just trim:
      CSO_NtsTrim = trim(s)
    end if
    
  end function CSO_NtsTrim
  


end module CSO_String


! ######################################################################
! ###
! ### test
! ###
! ######################################################################
!
!program test
!
!  use CSO_String
!  
!  character(len=32)   ::  s
!  integer             ::  status
!  
!  print *, 'not found ...'
!  s='abcd' ; call CSO_Replace( s, 'q', 'x', status ) ; print *, '"'//trim(s)//'"'
!
!  print *, 'replace 1 character ...'
!  s='abcd' ; call CSO_Replace( s, 'a', 'x', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, 'b', 'x', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, 'c', 'x', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, 'd', 'x', status ) ; print *, '"'//trim(s)//'"'
!    
!  print *, 'empty arguments ...'
!  s=''     ; call CSO_Replace( s, 'a', 'x', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, '' , 'x', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, 'a', '' , status ) ; print *, '"'//trim(s)//'"'
!
!  print *, 'replace 1 by 2 characters ...'
!  s='abcd' ; call CSO_Replace( s, 'a', 'XY', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, 'b', 'XY', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, 'd', 'XY', status ) ; print *, '"'//trim(s)//'"'
!
!  print *, 'replace 2 characters by 1 ...'
!  s='abcd' ; call CSO_Replace( s, 'ab', 'x', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, 'bc', 'x', status ) ; print *, '"'//trim(s)//'"'
!  s='abcd' ; call CSO_Replace( s, 'cd', 'x', status ) ; print *, '"'//trim(s)//'"'
!
!  print *, 'replace all ...'
!  s='abcdabcda' ; call CSO_Replace( s, 'a', '_', status ) ; print *, '"'//trim(s)//'"'
!  s='abcdabcda' ; call CSO_Replace( s, 'a', '' , status ) ; print *, '"'//trim(s)//'"'
!
!end program test
!
!
