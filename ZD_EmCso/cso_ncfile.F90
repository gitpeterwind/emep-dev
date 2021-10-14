!###############################################################################
!
! CSO_File - tools for writing netcdf dims/variables/attributes
!
!###############################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; csol=nf90_strerror(status); call csoErr; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!###############################################################################


module CSO_NcFile

  use CSO_Logging     , only : csol, csoPr, csoErr
  use NetCDF          , only : NF90_StrError, NF90_NOERR

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private

  public    ::  T_NcDims
  public    ::  T_NcAttrs
  public    ::  T_NcFile
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_NcFile'
  
  
  ! --- types ------------------------------
  
  
  !
  ! * dims
  !
  
  type ::  T_NcDim
    ! dimension name:
    character(:), allocatable     ::  name
    ! dimension length (might be local size ..)
    integer                       ::  length
    ! offset in global grid, 
    ! if non-zero the length is the local size:
    integer                       ::  offset
    ! global size:
    integer                       ::  glb_length
    ! netcdf id:
    integer                       ::  dimid
    !
  contains
    procedure ::  Init         =>  NcDim_Init
    procedure ::  Done         =>  NcDim_Done
    procedure ::  Def          =>  NcDim_Def
  end type T_NcDim
  
  ! *
  
  type ::  P_NcDim
    type(T_NcDim), pointer       ::  p
  end type P_NcDim
  
  ! *
  
  type ::  T_NcDims
    integer                       ::  n
    type(P_NcDim), pointer        ::  values(:)
  contains
    procedure ::  Init         =>  NcDims_Init
    procedure ::  Done         =>  NcDims_Done
    procedure ::  InqID        =>  NcDims_InqID
    procedure ::  Append       =>  NcDims_Append
    procedure ::  Def          =>  NcDims_Def
    procedure ::  Get          =>  NcDims_Get
    procedure ::                   NcDims_GetDim_name
    procedure ::                   NcDims_GetDim_id
    generic   ::  GetDim       =>  NcDims_GetDim_name, &
                                   NcDims_GetDim_id
    procedure ::                   NcDims_SetDim_name
    procedure ::                   NcDims_SetDim_id
    generic   ::  SetDim       =>  NcDims_SetDim_name, &
                                   NcDims_SetDim_id
  end type T_NcDims
  
  
  !
  ! * attrs
  !
  
  type ::  T_NcAttr
    ! attribute name:
    character(:), allocatable     ::  name
    ! attribute type:
    character(:), allocatable     ::  atype
    ! values:
    integer                       ::  ivalue
    real                          ::  rvalue
    character(:), allocatable     ::  cvalue
    !
  contains
    procedure ::                   NcAttr_Init_i
    procedure ::                   NcAttr_Init_r
    procedure ::                   NcAttr_Init_c
    generic   ::  Init         =>  NcAttr_Init_i, &
                                   NcAttr_Init_r, &
                                   NcAttr_Init_c
    procedure ::  InitCopy     =>  NcAttr_InitCopy
    procedure ::  Done         =>  NcAttr_Done
    procedure ::  NcGet        =>  NcAttr_NcGet
    procedure ::  NcPut        =>  NcAttr_NcPut
  end type T_NcAttr
  
  ! *
  
  type ::  P_NcAttr
    type(T_NcAttr), pointer       ::  p
  end type P_NcAttr
  
  ! *
  
  type ::  T_NcAttrs
    integer                       ::  n
    type(P_NcAttr), pointer       ::  values(:)
  contains
    procedure ::  Init         =>  NcAttrs_Init
    procedure ::  InitCopy     =>  NcAttrs_InitCopy
    procedure ::  Done         =>  NcAttrs_Done
    procedure ::                   NcAttrs_Append_i
    procedure ::                   NcAttrs_Append_r
    procedure ::                   NcAttrs_Append_c
    generic   ::  Append       =>  NcAttrs_Append_i, &
                                   NcAttrs_Append_r, &
                                   NcAttrs_Append_c
    procedure ::  NcGet        =>  NcAttrs_NcGet
    procedure ::  NcPut        =>  NcAttrs_NcPut
    procedure ::  GetValue     =>  NcAttrs_GetValue
  end type T_NcAttrs
  
  
  !
  ! *** var
  !
  
  type ::  T_NcVar
    ! variable name:
    character(:), allocatable       ::  name
    ! shape:
    integer                         ::  ndim
    ! dimension names:
    character(len=32), allocatable  ::  dims(:)  ! (ndim)
    ! data type:
    character(len=16)               ::  dtype
    ! variable attributes:
    type(T_NcAttrs)                 ::  attrs
    ! netcdf id:
    integer                         ::  varid
    !
  contains
    procedure ::  Init         =>  NcVar_Init
    procedure ::  Done         =>  NcVar_Done
    procedure ::  GetDims      =>  NcVar_GetDims
    procedure ::  Def          =>  NcVar_Def
    procedure ::                   NcVar_Put_1d_r
    procedure ::                   NcVar_Put_2d_r
    generic   ::  Put          =>  NcVar_Put_1d_r, &
                                   NcVar_Put_2d_r
  end type T_NcVar
  
  ! *
  
  type ::  P_NcVar
    type(T_NcVar), pointer       ::  p
  end type P_NcVar
  
  ! *
  
  type ::  T_NcVars
    integer                       ::  n
    type(P_NcVar), pointer        ::  values(:)
  contains
    procedure ::  Init         =>  NcVars_Init
    procedure ::  Done         =>  NcVars_Done
    procedure ::  Append       =>  NcVars_Append
    procedure ::  GetIndex     =>  NcVars_GetIndex
    procedure ::                   NcVars_Set_Attr_i
    procedure ::                   NcVars_Set_Attr_r
    procedure ::                   NcVars_Set_Attr_c
    generic   ::  Set_Attr     =>  NcVars_Set_Attr_i, &
                                   NcVars_Set_Attr_r, &
                                   NcVars_Set_Attr_c
    procedure ::  Def          =>  NcVars_Def
  end type T_NcVars
  
  
  !
  ! *** file
  !
  
  type ::  T_NcFile
    ! file name:
    character(:), allocatable       ::  filename
    ! file:
    integer                         ::  ncid
    ! dims:
    type(T_NcDims)                  ::  dims
    ! variables:
    type(T_NcVars)                  ::  vars
    ! global attributes:
    type(T_NcAttrs)                 ::  attrs
    !
  contains
    procedure ::  Init         =>  NcFile_Init
    procedure ::  Done         =>  NcFile_Done
    procedure ::  Def_Dim      =>  NcFile_Def_Dim
    procedure ::  Def_Var      =>  NcFile_Def_Var
    procedure ::                   NcFile_Set_Attr_i
    procedure ::                   NcFile_Set_Attr_r
    procedure ::                   NcFile_Set_Attr_c
    generic   ::  Set_Attr     =>  NcFile_Set_Attr_i, &
                                   NcFile_Set_Attr_r, &
                                   NcFile_Set_Attr_c
    procedure ::  EndDef       =>  NcFile_EndDef
    procedure ::  Put_Var      =>  NcFile_Put_Var_1d_r
    procedure ::  Put_Var2D    =>  NcFile_Put_Var2D_r
  end type T_NcFile


contains


  ! ====================================================================
  ! ===
  ! === NcDim
  ! ===
  ! ====================================================================


  subroutine NcDim_Init( self, name, length, status, &
                            offset )
                            
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcDim), intent(out)           ::  self
    character(len=*), intent(in)          ::  name
    integer, intent(in)                   ::  length
    integer, intent(out)                  ::  status

    integer, intent(in), optional         ::  offset

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDim_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%name   = trim(name)
    self%length = length
    
    ! parallel dimension?
    if ( present(offset) ) then
      ! store:
      self%offset = offset
      ! minimum for global lenth ...
      self%glb_length = self%offset + self%length
      ! global length is maximum of offset+length, filled on root only:
      call csoc%Reduce( 'max', self%glb_length, status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! no offset:
      self%offset = 0
      ! current size is global:
      self%glb_length = length
    end if
    
    ! ok
    status = 0
    
  end subroutine NcDim_Init
  
  ! *  

  subroutine NcDim_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDim), intent(inout)         ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDim_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    self%name = ''

    ! ok
    status = 0
    
  end subroutine NcDim_Done
  
  ! *  

  subroutine NcDim_Def( self, ncid, status )

    use NetCDF, only : NF90_Def_Dim

    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcDim), intent(inout)         ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDim_Def'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
      
    ! written on root...
    if ( csoc%root ) then
      ! define:
      status = NF90_Def_Dim( ncid, self%name, self%glb_length, self%dimid )
      if ( status /= NF90_NOERR ) then
        csol = nf90_strerror(status); call csoErr
        write (csol,'("from definition of dimension `",a,"` with length ",i0)') &
                        trim(self%name), self%glb_length; call csoErr
        TRACEBACK; status=1; return
      end if
    end if ! root

    ! ok
    status = 0
    
  end subroutine NcDim_Def


  ! ====================================================================
  ! ===
  ! === NcDims
  ! ===
  ! ====================================================================


  subroutine NcDims_Init( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(out)          ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! empty:
    self%n = 0
    nullify( self%values )
    
    ! ok
    status = 0
    
  end subroutine NcDims_Init
  
  
  ! *
  

  subroutine NcDims_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(inout)        ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_Done'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! defined?
    if ( self%n > 0 ) then
      ! loop:
      do i = 1, self%n
        ! done with dimension:
        call self%values(i)%p%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do ! i
      ! clear:
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if  ! n > 0

    ! ok
    status = 0
    
  end subroutine NcDims_Done
  
  
  ! *
  
  
  !
  ! Return integer id of dimension name.
  ! If not found, an error message is printed and error error status (>0) is returned,
  ! unless quiet=.true. which gives no message but a warning status (<0).
  !
  

  subroutine NcDims_InqID( self, name, id, status, quiet )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(in)             ::  self
    character(len=*), intent(in)            ::  name
    integer, intent(out)                    ::  id
    integer, intent(out)                    ::  status
    logical, intent(in), optional           ::  quiet

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_InqID'
    
    ! --- local ----------------------------------
    
    logical                     ::  verbose
    integer                     ::  i
    
    ! --- begin ----------------------------------
    
    ! shout?
    verbose = .true.
    if ( present(quiet) ) verbose = .not. quiet

    ! search:
    id = -999
    do i = 1, self%n
      if ( trim(self%values(i)%p%name) == trim(name) ) then
        id = i
        exit
      end if
    end do
    ! check ...
    if ( id < 0 ) then
      if ( verbose ) then
        write (csol,'("could not find name `",a,"` in dimensions:")') trim(name); call csoErr
        do i = 1, self%n
          write (csol,'(i6," ",a)') i, trim(self%values(i)%p%name); call csoErr
        end do
        TRACEBACK; status=1; return
      else
        ! warning status ...
        status=-1; return
      end if
    end if

    ! ok
    status = 0
    
  end subroutine NcDims_InqID
  
  
  ! *
  

  !
  ! Append new dimension if not defined yet;
  ! if already defined, length should be the same.
  !
  
  subroutine NcDims_Append( self, name, length, status, &
                                    offset )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(inout)        ::  self
    character(len=*), intent(in)          ::  name
    integer, intent(in)                   ::  length
    integer, intent(out)                  ::  status

    integer, intent(in), optional         ::  offset

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_Append'
    
    ! --- local ----------------------------------
    
    logical                     ::  found
    type(P_NcDim), pointer      ::  new_values(:)
    integer                     ::  i
    
    ! --- begin ----------------------------------
    
    ! init flag:
    found = .false.
    ! check if already present ...
    if ( self%n > 0 ) then
      ! loop over existing dims:
      do i = 1, self%n
        ! compare:
        if ( trim(self%values(i)%p%name) == trim(name) ) then
          ! check length:
          if ( length /= self%values(i)%p%length ) then
            ! ok if new length is undefined (negative), then keep current:
            if ( length >= 0 ) then
              write (csol,'("dimension `",a,"` with length ",i0," already defined with length ",i0)') &
                                trim(name), length, self%values(i)%p%length; call csoErr
              TRACEBACK; status=1; return
            end if
          end if ! different length
          ! set flag:
          found = .true.
        end if ! same name
      end do ! idim
    end if
    
    ! not present yet?
    if ( .not. found ) then

      ! new storage:
      allocate( new_values(self%n+1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! allready values present?
      if ( self%n > 0 ) then
        ! loop over current values:
        do i = 1, self%n
          ! assign value:
          new_values(i)%p => self%values(i)%p
        end do
        ! clear current storage:
        deallocate( self%values, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if ! n > 0
      ! assign value list:
      self%values => new_values    

      ! increase counter:
      self%n = self%n + 1
      ! new value:
      allocate( self%values(self%n)%p, stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! init:
      call self%values(self%n)%p%Init( name, length, status, offset=offset )
      IF_NOT_OK_RETURN(status=1)
      
    end if ! not found

    ! ok
    status = 0
    
  end subroutine NcDims_Append
  
  
  ! *
  

  subroutine NcDims_Get( self, status, n )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(in)             ::  self
    integer, intent(out)                    ::  status
    integer, intent(out), optional          ::  n

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! number of dims:
    if ( present(n) ) n = self%n
    
    ! ok
    status = 0
    
  end subroutine NcDims_Get
  
  
  ! *
  

  subroutine NcDims_GetDim_name( self, name, status, &
                                  length, dimid, offset, glb_length )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(in)             ::  self
    character(len=*), intent(in)            ::  name
    integer, intent(out)                    ::  status
    
    integer, intent(out), optional          ::  length
    integer, intent(out), optional          ::  dimid
    integer, intent(out), optional          ::  offset
    integer, intent(out), optional          ::  glb_length

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_GetDim_name'
    
    ! --- local ----------------------------------
    
    integer                     ::  id
    
    ! --- begin ----------------------------------
    
    ! inquire id by name:
    call self%InqID( name, id, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! get values:
    call self%GetDim( id, status, length=length, dimid=dimid, offset=offset, glb_length=glb_length )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcDims_GetDim_name
  
  
  ! *
  

  subroutine NcDims_GetDim_id( self, id, status, &
                                name, length, dimid, offset, glb_length )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(in)             ::  self
    integer, intent(in)                     ::  id
    integer, intent(out)                    ::  status
    
    character(len=*), intent(out), optional ::  name
    integer, intent(out), optional          ::  length
    integer, intent(out), optional          ::  dimid
    integer, intent(out), optional          ::  offset
    integer, intent(out), optional          ::  glb_length

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_GetDim_id'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! check ...
    if ( (id < 1) .or. (id > self%n) ) then
      write (csol,'("id should be in 1,..,",i0)') id, self%n; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! name?
    if ( present(name) ) name = trim(self%values(id)%p%name)

    ! length:
    if ( present(length) ) length = self%values(id)%p%length

    ! netcdf dimension id:
    if ( present(dimid) ) dimid = self%values(id)%p%dimid

    ! offset and global length:
    if ( present(offset    ) ) offset     = self%values(id)%p%offset
    if ( present(glb_length) ) glb_length = self%values(id)%p%glb_length
      
    ! ok
    status = 0
    
  end subroutine NcDims_GetDim_id
  
  
  ! *
  

  subroutine NcDims_SetDim_name( self, name, length, status, offset )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(inout)          ::  self
    character(len=*), intent(in)            ::  name
    integer, intent(in)                     ::  length
    integer, intent(out)                    ::  status
    
    integer, intent(in), optional           ::  offset

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_SetDim_name'
    
    ! --- local ----------------------------------
    
    integer                     ::  id
    
    ! --- begin ----------------------------------
    
    ! inquire id by name:
    call self%InqID( name, id, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! get values:
    call self%SetDim( id, length, status, offset=offset )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcDims_SetDim_name
  
  
  ! *
  

  subroutine NcDims_SetDim_id( self, id, length, status, offset )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(inout)          ::  self
    integer, intent(in)                     ::  id
    integer, intent(in)                     ::  length
    integer, intent(out)                    ::  status
    
    integer, intent(in), optional           ::  offset

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_SetDim_id'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! check ...
    if ( (id < 1) .or. (id > self%n) ) then
      write (csol,'("id should be in 1,..,",i0)') id, self%n; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! length:
    self%values(id)%p%length = length
    
    ! offset defined or default?
    self%values(id)%p%offset = 0
    if ( present(offset) ) self%values(id)%p%offset = offset
    ! global length:
    self%values(id)%p%glb_length = self%values(id)%p%offset + self%values(id)%p%length
      
    ! ok
    status = 0
    
  end subroutine NcDims_SetDim_id
  
  
  ! *
  

  subroutine NcDims_Def( self, ncid, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcDims), intent(in)           ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcDims_Def'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! loop over dimensions:
    do i = 1, self%n
      ! define dimension in file:
      call self%values(i)%p%Def( ncid, status )
      IF_NOT_OK_RETURN(status=1)
    end do

    ! ok
    status = 0
    
  end subroutine NcDims_Def


  ! ====================================================================
  ! ===
  ! === NcAttr
  ! ===
  ! ====================================================================


  subroutine NcAttr_Init_i( self, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttr), intent(out)          ::  self
    character(len=*), intent(in)          ::  name
    integer, intent(in)                   ::  value
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttr_Init_i'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store name:
    self%name = trim(name)
    ! set type:
    self%atype  = 'integer'
    ! store:
    self%ivalue = value
    
    ! ok
    status = 0
    
  end subroutine NcAttr_Init_i
  
  ! * 

  subroutine NcAttr_Init_r( self, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttr), intent(out)          ::  self
    character(len=*), intent(in)          ::  name
    real, intent(in)                      ::  value
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttr_Init_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store name:
    self%name = trim(name)
    ! set type:
    self%atype  = 'real'
    ! store:
    self%rvalue = value
    
    ! ok
    status = 0
    
  end subroutine NcAttr_Init_r
  
  ! *  

  subroutine NcAttr_Init_c( self, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttr), intent(out)          ::  self
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  value
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttr_Init_c'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store name:
    self%name = trim(name)
    ! set type:
    self%atype  = 'character'
    ! store:
    self%cvalue = trim(value)
    
    ! ok
    status = 0
    
  end subroutine NcAttr_Init_c
  
  ! *  

  subroutine NcAttr_InitCopy( self, attr, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttr), intent(out)          ::  self
    class(T_NcAttr), intent(in)           ::  attr
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttr_InitCopy'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! copy name:
    self%name = trim(attr%name)
    ! copy type:
    self%atype  = trim(attr%atype)

    ! switch:
    select case ( trim(self%atype) )
      case ( 'character' )
        self%cvalue = attr%cvalue
      case ( 'integer' )
        self%ivalue = attr%ivalue
      case ( 'real' )
        self%rvalue = attr%rvalue
      case default
        write (csol,'("unsupported attribute type `",a,"`")') trim(self%atype); call csoErr
        TRACEBACK; status=1; return
    end select
    
    ! ok
    status = 0
    
  end subroutine NcAttr_InitCopy
  
  ! * 

  subroutine NcAttr_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttr), intent(inout)        ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttr_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    self%name = ''

    ! switch:
    select case ( trim(self%atype) )
      case ( 'character' )
        self%cvalue = ''
      case ( 'integer' )
        self%ivalue = 0
      case ( 'real' )
        self%rvalue = 0.0
      case default
        write (csol,'("unsupported attribute type `",a,"`")') trim(self%atype); call csoErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
  end subroutine NcAttr_Done
  
  ! *  

  subroutine NcAttr_NcGet( self, ncid, varid, aname, status )

    use NetCDF, only : NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE, NF90_CHAR
    use NetCDF, only : NF90_Inquire_Attribute
    use NetCDF, only : NF90_Get_Att

    ! --- in/out ---------------------------------
    
    class(T_NcAttr), intent(inout)        ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(in)                   ::  varid
    character(len=*), intent(in)          ::  aname
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttr_NcGet'
    
    ! --- local ----------------------------------
    
    integer     ::  xtype
    integer     ::  nchar
    
    ! --- begin ----------------------------------
    
    ! store name:
    self%name = trim(aname)
    
    ! get type number:
    status = NF90_Inquire_Attribute( ncid, varid, aname, xtype=xtype )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! switch:
    select case ( xtype )
      !~ integers:
      case ( NF90_BYTE, NF90_SHORT, NF90_INT )
        ! fill type:
        self%atype = 'integer'
        ! read:
        status = NF90_Get_Att( ncid, varid, trim(self%name), self%ivalue )
        IF_NF90_NOT_OK_RETURN(status=1)
      !~ reals:
      case ( NF90_FLOAT, NF90_DOUBLE )
        ! fill type:
        self%atype = 'real'
        ! read:
        status = NF90_Get_Att( ncid, varid, trim(self%name), self%rvalue )
        IF_NF90_NOT_OK_RETURN(status=1)
      !~ chars:
      case ( NF90_CHAR )
        ! fill type:
        self%atype = 'character'
        ! get character length:
        status = NF90_Inquire_Attribute( ncid, varid, aname, len=nchar )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! storage:
        allocate( character(len=nchar) :: self%cvalue, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! read:
        status = NF90_Get_Att( ncid, varid, trim(self%name), self%cvalue )
        IF_NF90_NOT_OK_RETURN(status=1)
      !~
      case default
        write (csol,'("unsupported attribute type `",a,"`")') trim(self%atype); call csoErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
  end subroutine NcAttr_NcGet
  
  ! *  

  subroutine NcAttr_NcPut( self, ncid, varid, status )

    use NetCDF, only : NF90_Put_Att

    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttr), intent(inout)        ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(in)                   ::  varid
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttr_NcPut'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
      
    ! written on root...
    if ( csoc%root ) then
      ! switch:
      select case ( trim(self%atype) )
        case ( 'integer' )
          status = NF90_Put_Att( ncid, varid, trim(self%name), self%ivalue )
          IF_NF90_NOT_OK_RETURN(status=1)
        case ( 'real' )
          status = NF90_Put_Att( ncid, varid, trim(self%name), self%rvalue )
          IF_NF90_NOT_OK_RETURN(status=1)
        case ( 'character' )
          status = NF90_Put_Att( ncid, varid, trim(self%name), self%cvalue )
          IF_NF90_NOT_OK_RETURN(status=1)
        case default
          write (csol,'("unsupported attribute type `",a,"`")') trim(self%atype); call csoErr
          TRACEBACK; status=1; return
      end select
    end if ! root

    ! ok
    status = 0
    
  end subroutine NcAttr_NcPut


  ! ====================================================================
  ! ===
  ! === NcAttrs
  ! ===
  ! ====================================================================


  subroutine NcAttrs_Init( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(out)          ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! empty:
    self%n = 0
    nullify( self%values )
    
    ! ok
    status = 0
    
  end subroutine NcAttrs_Init
  
  
  ! *


  subroutine NcAttrs_InitCopy( self, attrs, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(out)         ::  self
    class(T_NcAttrs), intent(in)          ::  attrs
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_InitCopy'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! copy size:
    self%n = attrs%n
    ! any?
    if ( self%n > 0 ) then
      ! storage:
      allocate( self%values(self%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! loop:
      do i = 1, self%n
        ! new value:
        allocate( self%values(i)%p, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! copy:
        call self%values(i)%p%InitCopy( attrs%values(i)%p, status )
        IF_NOT_OK_RETURN(status=1)
      end do ! i
    end if ! n > 0
    
    ! ok
    status = 0
    
  end subroutine NcAttrs_InitCopy
  
  
  ! *
  

  subroutine NcAttrs_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(inout)        ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_Done'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! defined?
    if ( self%n > 0 ) then
      ! loop:
      do i = 1, self%n
        ! done with attrension:
        call self%values(i)%p%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do ! i
      ! clear:
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if  ! n > 0

    ! ok
    status = 0
    
  end subroutine NcAttrs_Done
  
  
  ! *
  

  subroutine NcAttrs_Append_Empty( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(inout)       ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_Append_Empty'
    
    ! --- local ----------------------------------
    
    type(P_NcAttr), pointer     ::  new_values(:)
    integer                     ::  i
    
    ! --- begin ----------------------------------
    
    ! new storage:
    allocate( new_values(self%n+1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! allready values present?
    if ( self%n > 0 ) then
      ! loop over current values:
      do i = 1, self%n
        ! assign value:
        new_values(i)%p => self%values(i)%p
      end do
      ! clear current storage:
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! n > 0
    ! assign value list:
    self%values => new_values    
    
    ! increase counter:
    self%n = self%n + 1
    ! new value:
    allocate( self%values(self%n)%p, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcAttrs_Append_Empty
    
  ! *
  
  subroutine NcAttrs_Append_i( self, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(inout)       ::  self
    character(len=*), intent(in)          ::  name
    integer, intent(in)                   ::  value
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_Append_i'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! extra element:
    call NcAttrs_Append_Empty( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! init:
    call self%values(self%n)%p%Init( name, value, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcAttrs_Append_i
    
  ! *
  
  subroutine NcAttrs_Append_r( self, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(inout)       ::  self
    character(len=*), intent(in)          ::  name
    real, intent(in)                      ::  value
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_Append_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! extra element:
    call NcAttrs_Append_Empty( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! init:
    call self%values(self%n)%p%Init( name, value, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcAttrs_Append_r
    
  ! *
  
  subroutine NcAttrs_Append_c( self, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(inout)       ::  self
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  value
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_Append_c'
    
    ! --- local ----------------------------------
    
    integer     ::  ivalue
    real        ::  rvalue
    
    ! --- begin ----------------------------------
    
    ! extra element:
    call NcAttrs_Append_Empty( self, status )
    IF_NOT_OK_RETURN(status=1)

    ! check on conversions:
    if ( len_trim(value) > 6 ) then
      ! convert to real?
      if ( value(1:6) == 'float:' ) then
        ! extract:
        read (value(7:),*,iostat=status) rvalue
        if ( status /= 0 ) then
          write (csol,'("could not extract float from attributed value `",a,"`")') trim(value); call csoErr
          TRACEBACK; status=1; return
        end if
        ! store:
        call self%values(self%n)%p%Init( name, rvalue, status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! store:
        call self%values(self%n)%p%Init( name, value, status )
        IF_NOT_OK_RETURN(status=1)
      end if
      !
    else if ( len_trim(value) > 4 ) then
      ! convert to integer?
      if (  value(1:4) == 'int:' ) then
        ! extract:
        read (value(5:),*,iostat=status) ivalue
        if ( status /= 0 ) then
          write (csol,'("could not extract integer from attributed value `",a,"`")') trim(value); call csoErr
          TRACEBACK; status=1; return
        end if
        ! store:
        call self%values(self%n)%p%Init( name, ivalue, status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! store:
        call self%values(self%n)%p%Init( name, value, status )
        IF_NOT_OK_RETURN(status=1)
      end if
      !
    else
      ! store:
      call self%values(self%n)%p%Init( name, value, status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine NcAttrs_Append_c
  
  
  ! *
  

  subroutine NcAttrs_NcGet( self, ncid, varid, status )

    use NetCDF, only : NF90_Inquire_Variable
    use NetCDF, only : NF90_Inq_AttName
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(inout)       ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(in)                   ::  varid
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_NcGet'
    
    ! --- local ----------------------------------
    
    integer               ::  natts
    integer               ::  attnum
    character(len=256)    ::  aname
    
    ! --- begin ----------------------------------
    
    ! number of attributes:
    status = NF90_Inquire_Variable( ncid, varid, natts=natts )
    IF_NF90_NOT_OK_RETURN(status=1)    
    ! loop over attributes:
    do attnum = 1, natts
      ! get name:
      status = NF90_Inq_AttName( ncid, varid, attnum, aname )
      IF_NOT_OK_RETURN(status=1)
      ! filter ...
      if ( trim(aname) == '_FillValue' ) cycle
      ! extra element:
      call NcAttrs_Append_Empty( self, status )
      IF_NOT_OK_RETURN(status=1)
      ! fill from file:
      call self%values(self%n)%p%NcGet( ncid, varid, aname, status )
      IF_NOT_OK_RETURN(status=1)
    end do

    ! ok
    status = 0
    
  end subroutine NcAttrs_NcGet
  
  
  ! *
  

  subroutine NcAttrs_NcPut( self, ncid, varid, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(in)          ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(in)                   ::  varid
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_NcPut'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! loop over attributes:
    do i = 1, self%n
      ! define attrension in file:
      call self%values(i)%p%NcPut( ncid, varid, status )
      IF_NOT_OK_RETURN(status=1)
    end do

    ! ok
    status = 0
    
  end subroutine NcAttrs_NcPut
  
  
  ! *
  

  subroutine NcAttrs_GetValue( self, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(in)          ::  self
    character(len=*), intent(in)          ::  name
    character(len=*), intent(out)         ::  value
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_GetValue'
    
    ! --- local ----------------------------------
    
    logical     ::  found
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! init flag:
    found = .false.
    ! loop over attributes:
    do i = 1, self%n
      ! compare:
      found = trim(self%values(i)%p%name) == trim(name)
      ! found?
      if ( found ) then
        ! check type ...
        if ( trim(self%values(i)%p%atype) == 'character' ) then
          value = trim(self%values(i)%p%cvalue)
        else
          write (csol,'("attriubute `",a,"` type `",a,"` could not be converted to type `",a,"`")') &
                  trim(name), trim(self%values(i)%p%atype), 'character'; call csoErr
          TRACEBACK; status=1; return
        end if
        ! leave:
        exit
      end if
    end do ! attributes
    
    ! not found?
    if ( .not. found ) then
      write (csol,'("attribute `",a,"` not found")') trim(name); call csoErr
      TRACEBACK; status=1; return
    end if

    ! ok
    status = 0
    
  end subroutine NcAttrs_GetValue


  ! ====================================================================
  ! ===
  ! === NcVar
  ! ===
  ! ====================================================================


  subroutine NcVar_Init( self, name, dtype, dims, status )
                            
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcVar), intent(out)           ::  self
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  dtype
    character(len=*), intent(in)          ::  dims(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVar_Init'
    
    ! --- local ----------------------------------
    
    integer     ::  k
    
    ! --- begin ----------------------------------
    
    ! store:
    self%name  = trim(name)
    self%dtype = dtype
    
    ! shape:
    self%ndim = size(dims)
    ! storage:
    allocate( self%dims(self%ndim), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do k = 1, self%ndim
      ! copy:
      self%dims(k) = trim(dims(k))
    end do
    
    ! init attributes:
    call self%attrs%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcVar_Init
  
  ! *  

  subroutine NcVar_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVar), intent(inout)         ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVar_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%dims, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! done with attributes:
    call self%attrs%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcVar_Done
  
  ! *  

  subroutine NcVar_GetDims( self, ncdims, status, &
                               offsets, glb_shape )

    ! --- in/out ---------------------------------
    
    class(T_NcVar), intent(inout)         ::  self
    type(T_NcDims), intent(in)            ::  ncdims
    integer, intent(out)                  ::  status
    
    integer, intent(out), optional        ::  offsets(:)
    integer, intent(out), optional        ::  glb_shape(:)

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVar_GetDims'
    
    ! --- local ----------------------------------
    
    integer                 ::  k
    
    ! --- begin ----------------------------------
    
    ! loop over dims:
    do k = 1, self%ndim
      ! offset?
      if ( present(offsets) ) then
        call ncdims%GetDim( self%dims(k), status, offset=offsets(k) )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! global shape?
      if ( present(glb_shape) ) then
        call ncdims%GetDim( self%dims(k), status, glb_length=glb_shape(k) )
        IF_NOT_OK_RETURN(status=1)
      end if
    end do ! dims

    ! ok
    status = 0
    
  end subroutine NcVar_GetDims
  
  ! *  

  subroutine NcVar_Def( self, ncid, ncdims, status )

    use NetCDF  , only : NF90_FLOAT
    use NetCDF  , only : NF90_Def_Var
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcVar), intent(inout)         ::  self
    integer, intent(in)                   ::  ncid
    type(T_NcDims), intent(in)            ::  ncdims
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVar_Def'
    
    ! --- local ----------------------------------
    
    integer                 ::  dtype
    integer, allocatable    ::  dimids(:)
    integer                 ::  k
    
    ! --- begin ----------------------------------
    
    ! written on root...
    if ( csoc%root ) then

      ! data type:
      select case ( trim(self%dtype) )
        case ( 'float', 'real' )
          dtype = NF90_FLOAT
        case default
          write (csol,'("unsupported dtype `",a,"`")') trim(self%dtype); call csoErr
          TRACEBACK; status=1; return
      end select
      
      ! storage:
      allocate( dimids(self%ndim), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! loop over dims:
      do k = 1, self%ndim
        ! netcdf dimension id:
        call ncdims%GetDim( self%dims(k), status, dimid=dimids(k) )
        IF_NOT_OK_RETURN(status=1)
      end do
      
      ! define:
      status = NF90_Def_Var( ncid, self%name, dtype, dimids, self%varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      
      ! add attributes:
      call self%attrs%NcPut( ncid, self%varid, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( dimids, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if ! root

    ! ok
    status = 0
    
  end subroutine NcVar_Def
  
  
  ! *
  

  subroutine NcVar_Put_1d_r( self, ncid, values, status )

    use NetCDF  , only : NF90_Put_Var
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcVar), intent(inout)         ::  self
    integer, intent(in)                   ::  ncid
    real, intent(in)                      ::  values(:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVar_Put_1d_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! written on root...
    if ( csoc%root ) then
      ! put:
      status = NF90_Put_Var( ncid, self%varid, values )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine NcVar_Put_1d_r
  
  
  ! *
  

  subroutine NcVar_Put_2d_r( self, ncid, values, status )

    use NetCDF  , only : NF90_Put_Var
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcVar), intent(inout)         ::  self
    integer, intent(in)                   ::  ncid
    real, intent(in)                      ::  values(:,:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVar_Put_2d_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! written on root...
    if ( csoc%root ) then
      ! put:
      status = NF90_Put_Var( ncid, self%varid, values )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine NcVar_Put_2d_r


  ! ====================================================================
  ! ===
  ! === NcVars
  ! ===
  ! ====================================================================


  subroutine NcVars_Init( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVars), intent(out)          ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVars_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! empty:
    self%n = 0
    nullify( self%values )
    
    ! ok
    status = 0
    
  end subroutine NcVars_Init
  
  
  ! *
  

  subroutine NcVars_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVars), intent(inout)        ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVars_Done'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! defined?
    if ( self%n > 0 ) then
      ! loop:
      do i = 1, self%n
        ! done with varension:
        call self%values(i)%p%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do ! i
      ! clear:
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if  ! n > 0

    ! ok
    status = 0
    
  end subroutine NcVars_Done
  
  
  ! *
  

  subroutine NcVars_Append( self, name, dtype, dims, status, &
                               ivar )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVars), intent(inout)        ::  self
    character(len=*), intent(in)          ::  name
    character(len=*), intent(in)          ::  dtype
    character(len=*), intent(in)          ::  dims(:)
    integer, intent(out)                  ::  status
    
    integer, intent(out), optional        ::  ivar

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVars_Append'
    
    ! --- local ----------------------------------
    
    type(P_NcVar), pointer      ::  new_values(:)
    integer                     ::  i
    
    ! --- begin ----------------------------------
    
    ! new storage:
    allocate( new_values(self%n+1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! allready values present?
    if ( self%n > 0 ) then
      ! loop over current values:
      do i = 1, self%n
        ! assign value:
        new_values(i)%p => self%values(i)%p
      end do
      ! clear current storage:
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! n > 0
    ! assign value list:
    self%values => new_values    
    
    ! increase counter:
    self%n = self%n + 1
    ! new value:
    allocate( self%values(self%n)%p, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! init:
    call self%values(self%n)%p%Init( name, dtype, dims, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! return index?
    if ( present(ivar) ) ivar = self%n

    ! ok
    status = 0
    
  end subroutine NcVars_Append
  
  
  ! *
  

  subroutine NcVars_GetIndex( self, indx, status, ivar, varname )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVars), intent(inout)          ::  self
    integer, intent(out)                    ::  indx
    integer, intent(out)                    ::  status
    integer, intent(in), optional           ::  ivar
    character(len=*), intent(in), optional  ::  varname

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVars_Get'
    
    ! --- local ----------------------------------
    
    integer                     ::  i
    
    ! --- begin ----------------------------------
    
    ! index:
    if ( present(ivar) ) then
      ! check ...
      if ( (ivar < 1) .or. (ivar > self%n) ) then
        write (csol,'("ivar (",i0,") should be in 1,..,",i0)') ivar, self%n; call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      indx = ivar
    else if ( present(varname) ) then
      ! search:
      indx = -999
      do i = 1, self%n
        if ( trim(self%values(i)%p%name) == trim(varname) ) then
          indx = i
          exit
        end if
      end do
      ! check ...
      if ( indx < 0 ) then
        write (csol,'("could not find name `",a,"` in varensions:")') trim(varname); call csoErr
        do i = 1, self%n
          write (csol,'(i6," ",a)') i, trim(self%values(i)%p%name); call csoErr
        end do
        TRACEBACK; status=1; return
      end if
    else
      write (csol,'("either specify `ivar` or `name` argument")'); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine NcVars_GetIndex
  
  
  ! *

  subroutine NcVars_Set_Attr_i( self, ivar, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVars), intent(inout)          ::  self
    integer, intent(in)                     ::  ivar
    character(len=*), intent(in)            ::  name
    integer, intent(in)                     ::  value
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVars_Set_Attr_i'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! check ...
    if ( (ivar < 1) .or. (ivar > self%n) ) then
      write (csol,'("ivar (",i0,") should be in 1,..,",i0)') ivar, self%n; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! add attribute:
    call self%values(ivar)%p%attrs%Append( name, value, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcVars_Set_Attr_i
  
  ! *

  subroutine NcVars_Set_Attr_r( self, ivar, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVars), intent(inout)          ::  self
    integer, intent(in)                     ::  ivar
    character(len=*), intent(in)            ::  name
    real, intent(in)                        ::  value
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVars_Set_Attr_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ivar < 1) .or. (ivar > self%n) ) then
      write (csol,'("ivar (",i0,") should be in 1,..,",i0)') ivar, self%n; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! add attribute:
    call self%values(ivar)%p%attrs%Append( name, value, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcVars_Set_Attr_r
  
  ! *

  subroutine NcVars_Set_Attr_c( self, ivar, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVars), intent(inout)          ::  self
    integer, intent(in)                     ::  ivar
    character(len=*), intent(in)            ::  name
    character(len=*), intent(in)            ::  value
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVars_Set_Attr_c'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ivar < 1) .or. (ivar > self%n) ) then
      write (csol,'("ivar (",i0,") should be in 1,..,",i0)') ivar, self%n; call csoErr
      TRACEBACK; status=1; return
    end if

    ! store:
    call self%values(ivar)%p%attrs%Append( name, value, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcVars_Set_Attr_c
  
  
  ! *
  

  subroutine NcVars_Def( self, ncid, ncdims, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcVars), intent(inout)        ::  self
    integer, intent(in)                   ::  ncid
    type(T_NcDims), intent(in)            ::  ncdims
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcVars_Def'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! loop over varensions:
    do i = 1, self%n
      ! define variables in file:
      call self%values(i)%p%Def( ncid, ncdims, status )
      IF_NOT_OK_RETURN(status=1)
    end do

    ! ok
    status = 0
    
  end subroutine NcVars_Def


  ! ====================================================================
  ! ===
  ! === NcFile
  ! ===
  ! ====================================================================


  subroutine NcFile_Init( self, filename, status )

    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(out)          ::  self
    character(len=*), intent(in)          ::  filename
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Init'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! init dimension list:
    call self%dims%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init variable list:
    call self%vars%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init attributes list:
    call self%attrs%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcFile_Init
  
  
  ! *
  

  subroutine NcFile_Done( self, status )

    use NetCDF, only : NF90_Close
  
    use CSO_Comm        , only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! written on root...
    if ( csoc%root ) then
      ! close:
      status = NF90_Close( self%ncid )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if  ! root
    
    ! done with dimensions:
    call self%dims%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! done with variables:
    call self%vars%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! done with attributes:
    call self%attrs%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    self%filename = ''
      
    ! ok
    status = 0
    
  end subroutine NcFile_Done
  
  ! *

  subroutine NcFile_Def_Dim( self, name, length, status, &
                                offset )
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    character(len=*), intent(in)          ::  name
    integer, intent(in)                   ::  length
    integer, intent(out)                  ::  status

    integer, intent(in), optional         ::  offset

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Def_Dim'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    call self%dims%Append( name, length, status, offset=offset )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcFile_Def_Dim
  
  ! *

  subroutine NcFile_Def_Var( self, name, dims, status, &
                               ivar )
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)            ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  dims(:)
    integer, intent(out)                      ::  status
    
    integer, intent(out), optional            ::  ivar

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Def_Var'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! define:
    call self%vars%Append( name, 'float', dims, status, ivar=ivar )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcFile_Def_Var
  
  ! *

  subroutine NcFile_Set_Attr_i( self, ivar, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)            ::  self
    integer, intent(in)                       ::  ivar
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  value
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Set_Attr_i'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! file or variable attribute?
    if ( ivar < 1 ) then
      ! add file attribute:
      call self%attrs%Append( name, value, status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! add variable attribute:
      call self%vars%Set_Attr( ivar, name, value, status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine NcFile_Set_Attr_i
  
  ! *

  subroutine NcFile_Set_Attr_r( self, ivar, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)            ::  self
    integer, intent(in)                       ::  ivar
    character(len=*), intent(in)              ::  name
    real, intent(in)                          ::  value
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Set_Attr_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! file or variable attribute?
    if ( ivar < 1 ) then
      ! add file attribute:
      call self%attrs%Append( name, value, status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! add variable attribute:
      call self%vars%Set_Attr( ivar, name, value, status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine NcFile_Set_Attr_r
  
  ! *

  subroutine NcFile_Set_Attr_c( self, ivar, name, value, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)            ::  self
    integer, intent(in)                       ::  ivar
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  value
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Set_Attr_c'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! file or variable attribute?
    if ( ivar < 1 ) then
      ! add file attribute:
      call self%attrs%Append( name, value, status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! add variable attribute:
      call self%vars%Set_Attr( ivar, name, value, status )
     IF_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine NcFile_Set_Attr_c
    
  ! *
  
  subroutine NcFile_EndDef( self, status )
  
    use NetCDF  , only : NF90_Create
    use NetCDF  , only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF  , only : NF90_GLOBAL
    use NetCDF  , only : NF90_EndDef
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Def_Dim'
    
    ! --- local ----------------------------------
    
    integer                   ::  cmode
    
    ! --- begin ----------------------------------
    
    ! written on root...
    if ( csoc%root ) then

      ! set creation mode flag:
      cmode = NF90_CLOBBER       ! overwrite existing files
      !cmode = NF90_NOCLOBBER     ! do not overwrite existing files

      ! create file:
      status = NF90_Create( self%filename, cmode, self%ncid )
      if ( status /= NF90_NOERR ) then
         write (csol,'("creating file :")'); call csoErr
         write (csol,'("  ",a)') trim(self%filename); call csoErr
         TRACEBACK; status=1; return
      end if
      
      ! define all dimensions:
      call self%dims%Def( self%ncid, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! define all variables:
      call self%vars%Def( self%ncid, self%dims, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! add file attributes:
      call self%attrs%NcPut( self%ncid, NF90_GLOBAL, status )
      IF_NOT_OK_RETURN(status=1)

      ! end defintion mode:
      status = NF90_EndDef( self%ncid )
      IF_NF90_NOT_OK_RETURN(status=1)

    end if  ! root

    ! ok
    status = 0
    
  end subroutine NcFile_EndDef
  
  
  ! *
  

  subroutine NcFile_Put_Var_1d_r( self, ivar, values, status, empty )
  
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(in)                   ::  ivar
    real, intent(in)                      ::  values(:)
    integer, intent(out)                  ::  status
    
    logical, intent(in), optional         ::  empty

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Put_Var_1d_r'
    
    ! --- local ----------------------------------
    
    integer             ::  nloc
    integer             ::  nglb
    real, allocatable   ::  glb(:)
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ivar < 1) .or. (ivar > self%vars%n) ) then
      write (csol,'("ivar (",i0,") should be in 1,..,",i0)') ivar, self%vars%n; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! number of local values:
    nloc = size(values)
    ! might need to ignore the local values:
    if ( present(empty) ) then
      if ( empty ) nloc = 0
    end if
    
    ! total number:
    call csoc%ParInfo( nloc, status, ntot=nglb )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( glb(nglb), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! gather:
    call csoc%GatherV( values, glb, status, nloc=nloc )
    IF_NOT_OK_RETURN(status=1)
    
    ! write:
    call self%vars%values(ivar)%p%Put( self%ncid, glb, status )
    IF_NOT_OK_RETURN(status=1)
        
    ! clear:
    deallocate( glb, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcFile_Put_Var_1d_r
  
  
  ! *
  

  subroutine NcFile_Put_Var2D_r( self, ivar, values, status )
  
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(in)                   ::  ivar
    real, intent(in)                      ::  values(:,:)
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Put_Var2D_r'
    
    ! --- local ----------------------------------
    
    type(T_NcVar), pointer    ::  varp
    integer                   ::  offs(2)
    integer                   ::  gshp(2)
    real, allocatable         ::  glb(:,:)
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ivar < 1) .or. (ivar > self%vars%n) ) then
      write (csol,'("ivar (",i0,") should be in 1,..,",i0)') ivar, self%vars%n; call csoErr
      TRACEBACK; status=1; return
    end if
    ! short:
    varp => self%vars%values(ivar)%p
    
    ! offset and global shape:
    call varp%GetDims( self%dims, status, offsets=offs, glb_shape=gshp )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for global array, only needed on root:
    if ( csoc%root ) then
      allocate( glb(gshp(1),gshp(2)), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      allocate( glb(1,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! gather:
    call csoc%Gather2D( values, offs, glb, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write:
    call varp%Put( self%ncid, glb, status )
    IF_NOT_OK_RETURN(status=1)
        
    ! clear:
    deallocate( glb, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcFile_Put_Var2D_r
  

end module CSO_NcFile
  
