!###############################################################################
!
! CSO_File - tools for writing netcdf dims/variables/attributes
!
!
! HISTORY
!
! 2022-09, Arjo Segers
!   Support input and output of packed variables.
! 2022-10, Arjo Segers
!   Do not copy packing attributes from input files.
! 2023-01, Arjo Segers
!   When packing variables with single value, 
!   use add_offset=value and scale_factor=1
!   to avoid division by zero for all-zero variables.
! 2023-01, Arjo Segers
!   Support integer(1) and character variables.
!
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
  
  public    ::  rwp__out
  public    ::  iwp__packed
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_NcFile'
  
  ! real kind for output variables:
  integer, parameter            ::  rwp__out    = 4
  ! integer kind for packed variables:
  integer, parameter            ::  iwp__packed = 2
  
  
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
    real(4)                       ::  rvalue
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
    ! number of values:
    integer                       ::  n
    ! list of values:
    type(P_NcAttr), pointer       ::  values(:)
    ! how to read?
    logical                       ::  read_on_root
    logical                       ::  read_by_me
    !
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
    ! read or write?
    character(len=1)                ::  rwmode
    ! file:
    integer                         ::  ncid
    ! netcdf flavour:    
    integer                         ::  formatNum
    ! dims:
    type(T_NcDims)                  ::  dims
    ! variables:
    type(T_NcVars)                  ::  vars
    ! global attributes:
    type(T_NcAttrs)                 ::  attrs
    !
  contains
    procedure ::  Init             =>  NcFile_Init
    procedure ::  Done             =>  NcFile_Done
    procedure ::  Inquire          =>  NcFile_Inquire
    procedure ::  Inq_VarID        =>  NcFile_Inq_VarID
    procedure ::  Inq_Variable     =>  NcFile_Inq_Variable
    procedure ::  Inq_VarPacking   =>  NcFile_Inq_VarPacking
    procedure ::  Inq_VarMissing   =>  NcFile_Inq_VarMissing
    procedure ::  Inq_VarUnits     =>  NcFile_Inq_VarUnits
    !
    procedure   ::                     NcFile_Get_Var_i_1d
    procedure   ::                     NcFile_Get_Var_i1_1d
    procedure   ::                     NcFile_Get_Var_c_2d
    procedure   ::                     NcFile_Get_Var_i_2d
    procedure   ::                     NcFile_Get_Var_i_3d
    procedure   ::                     NcFile_Get_Var_r_1d
    procedure   ::                     NcFile_Get_Var_r_2d
    procedure   ::                     NcFile_Get_Var_r_3d
    generic     ::  Get_Var        =>  NcFile_Get_Var_i_1d, &
                                       NcFile_Get_Var_i1_1d, &
                                       NcFile_Get_Var_c_2d, &
                                       NcFile_Get_Var_i_2d, &
                                       NcFile_Get_Var_i_3d, &
                                       NcFile_Get_Var_r_1d, &
                                       NcFile_Get_Var_r_2d, &
                                       NcFile_Get_Var_r_3d
    !
    procedure ::  Def_Dim          =>  NcFile_Def_Dim
    procedure ::  Def_Var          =>  NcFile_Def_Var
    procedure ::                       NcFile_Set_Attr_i
    procedure ::                       NcFile_Set_Attr_r
    procedure ::                       NcFile_Set_Attr_c
    generic   ::  Set_Attr         =>  NcFile_Set_Attr_i, &
                                       NcFile_Set_Attr_r, &
                                       NcFile_Set_Attr_c
    procedure ::  EndDef           =>  NcFile_EndDef
    procedure ::  GetPacking       =>  NcFile_GetPacking
    procedure ::                       NcFile_Put_Var_1d_r
    procedure ::                       NcFile_Put_Var_2d_c
    procedure ::                       NcFile_Put_Var_2d_r
    procedure ::                       NcFile_Put_Var_3d_r
    generic   ::  Put_Var          =>  NcFile_Put_Var_1d_r, &
                                       NcFile_Put_Var_2d_c, &
                                       NcFile_Put_Var_2d_r, &
                                       NcFile_Put_Var_3d_r
    procedure ::  Put_Var1D        =>  NcFile_Put_Var1D_r
    procedure ::  Put_Var2D        =>  NcFile_Put_Var2D_r
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
  
  ! Read attribute value and store.
  ! Arguments:
  !  ncid       : netcdf file id
  !  varid      : netcdf variable id
  !  aname      : attribute name
  ! Optional arguments:
  !  read_on_root   :  set to .true. to let only root read the attribute,
  !                    default .false. (all processors read)
  !  only_here      :  if "read_on_root=.true." a broadcast to the other processors is done,
  !                    unless this flag is .true.
  ! Return value:
  !  status     : non-zero in case of error

  subroutine NcAttr_NcGet( self, ncid, varid, aname, status, &
                             read_on_root, only_here )

    use NetCDF, only : NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE, NF90_CHAR
    use NetCDF, only : NF90_Inquire_Attribute
    use NetCDF, only : NF90_Get_Att
    
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------
    
    class(T_NcAttr), intent(inout)        ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(in)                   ::  varid
    character(len=*), intent(in)          ::  aname
    integer, intent(out)                  ::  status
    logical, intent(in), optional         ::  read_on_root
    logical, intent(in), optional         ::  only_here

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttr_NcGet'
    
    ! --- local ----------------------------------
    
    logical     ::  with_bcast
    logical     ::  do_read_on_root
    logical     ::  do_read_by_me
    integer     ::  xtype
    integer     ::  nchar
    
    ! --- begin ----------------------------------
    
    ! store name:
    self%name = trim(aname)
    
    ! read by root and scatter?
    do_read_on_root = .false.
    if ( present(read_on_root) ) do_read_on_root = read_on_root
    ! read here?
    if ( do_read_on_root ) then
      do_read_by_me = csoc%root
    else
      do_read_by_me = .true.
    end if
    
    ! by default perform broadcasts in case of reading on root,
    ! but skip if only this pe should read:
    with_bcast = .true.
    if ( present(only_here) ) with_bcast = .not. only_here
    
    ! read here?
    if ( do_read_by_me ) then
      ! get type number:
      status = NF90_Inquire_Attribute( ncid, varid, aname, xtype=xtype )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if ! read by me
    ! need to broadcast?
    if ( do_read_on_root .and. with_bcast ) then
      ! broadcast from root:
      call csoc%BCast( csoc%root_id, xtype, status )
      IF_NOT_OK_RETURN(status=1)
    end if  ! broadcast

    ! switch:
    select case ( xtype )
      !~ integers;
      !  also try to read types not supported by NF90:
      !    7 unsigned byte
      case ( NF90_BYTE, NF90_SHORT, NF90_INT, 7 )
        ! fill type:
        self%atype = 'integer'
        ! read here?
        if ( do_read_by_me ) then
          ! read:
          status = NF90_Get_Att( ncid, varid, trim(self%name), self%ivalue )
          IF_NF90_NOT_OK_RETURN(status=1)
        end if ! read by me
        ! need to broadcast?
        if ( do_read_on_root .and. with_bcast ) then
          ! broadcast from root:
          call csoc%BCast( csoc%root_id, self%ivalue , status )
          IF_NOT_OK_RETURN(status=1)
        end if  ! broadcast
      !~ reals:
      case ( NF90_FLOAT, NF90_DOUBLE )
        ! fill type:
        self%atype = 'real'
        ! read here?
        if ( do_read_by_me ) then
          ! read:
          status = NF90_Get_Att( ncid, varid, trim(self%name), self%rvalue )
          IF_NF90_NOT_OK_RETURN(status=1)
        end if ! read by me
        ! need to broadcast?
        if ( do_read_on_root .and. with_bcast ) then
          ! broadcast from root:
          call csoc%BCast( csoc%root_id, self%rvalue , status )
          IF_NOT_OK_RETURN(status=1)
        end if  ! broadcast
      !~ chars:
      case ( NF90_CHAR )
        ! fill type:
        self%atype = 'character'
        ! read here?
        if ( do_read_by_me ) then
          ! get character length:
          status = NF90_Inquire_Attribute( ncid, varid, aname, len=nchar )
          IF_NF90_NOT_OK_RETURN(status=1)
        end if ! read by me
        ! need to broadcast?
        if ( do_read_on_root .and. with_bcast ) then
          ! broadcast from root:
          call csoc%BCast( csoc%root_id, nchar, status )
          IF_NOT_OK_RETURN(status=1)
        end if  ! broadcast
        ! storage:
        allocate( character(len=nchar) :: self%cvalue, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! read here?
        if ( do_read_by_me ) then
          ! read:
          status = NF90_Get_Att( ncid, varid, trim(self%name), self%cvalue )
          IF_NF90_NOT_OK_RETURN(status=1)
        end if ! read by me
        ! need to broadcast?
        if ( do_read_on_root .and. with_bcast ) then
          ! broadcast from root:
          call csoc%BCast( csoc%root_id, self%cvalue , status )
          IF_NOT_OK_RETURN(status=1)
        end if  ! broadcast
      !~
      case default
        write (csol,'("attribute `",a,"` has unsupported type id ",i0,"")') &
                               trim(aname), xtype; call csoErr
        write (csol,'("maybe unsigned byte or int?")'); call csoErr
        write (csol,'("supported by NF90 library:")'); call csoErr
        write (csol,'(i8,"  NF90_BYTE  ")') NF90_BYTE  ; call csoErr
        write (csol,'(i8,"  NF90_CHAR  ")') NF90_CHAR  ; call csoErr
        write (csol,'(i8,"  NF90_SHORT ")') NF90_SHORT ; call csoErr
        write (csol,'(i8,"  NF90_INT   ")') NF90_INT   ; call csoErr
        write (csol,'(i8,"  NF90_FLOAT ")') NF90_FLOAT ; call csoErr
        write (csol,'(i8,"  NF90_DOUBLE")') NF90_DOUBLE; call csoErr
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


  subroutine NcAttrs_Init( self, status, read_on_root )
  
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(out)         ::  self
    integer, intent(out)                  ::  status
    logical, intent(in), optional         ::  read_on_root

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! empty:
    self%n = 0
    nullify( self%values )
    
    ! read by root and scatter?
    self%read_on_root = .false.
    if ( present(read_on_root) ) self%read_on_root = read_on_root
    ! read here?
    if ( self%read_on_root ) then
      self%read_by_me = csoc%root
    else
      self%read_by_me = .true.
    end if
    
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
  

  subroutine NcAttrs_NcGet( self, ncid, varid, status, only_here )

    use NetCDF, only : NF90_Inquire_Variable
    use NetCDF, only : NF90_Inq_AttName
    
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcAttrs), intent(inout)       ::  self
    integer, intent(in)                   ::  ncid
    integer, intent(in)                   ::  varid
    integer, intent(out)                  ::  status
    logical, intent(in), optional         ::  only_here

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcAttrs_NcGet'
    
    ! --- local ----------------------------------
    
    logical               ::  with_bcast
    integer               ::  natts
    integer               ::  attnum
    character(len=256)    ::  aname
    
    ! --- begin ----------------------------------

    ! by default perform broadcasts in case of reading on root,
    ! but skip if only this pe should read:
    with_bcast = .true.
    if ( present(only_here) ) with_bcast = .not. only_here

    ! read here?
    if ( self%read_by_me ) then
      ! number of attributes:
      status = NF90_Inquire_Variable( ncid, varid, natts=natts )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if ! read by me
    ! need to broadcast?
    if ( self%read_on_root .and. with_bcast ) then
      ! broadcast from root:
      call csoc%BCast( csoc%root_id, natts, status )
      IF_NOT_OK_RETURN(status=1)
    end if  ! broadcast  

    ! loop over attributes:
    do attnum = 1, natts
      ! read here?
      if ( self%read_by_me ) then
        ! get name:
        status = NF90_Inq_AttName( ncid, varid, attnum, aname )
        IF_NOT_OK_RETURN(status=1)
      end if ! read by me
      ! need to broadcast?
      if ( self%read_on_root .and. with_bcast ) then
        ! broadcast from root:
        call csoc%BCast( csoc%root_id, aname, status )
        IF_NOT_OK_RETURN(status=1)
      end if  ! broadcast  
      ! do not copy special attributes,
      ! otherwise they are copied into the output files too:
      if ( trim(aname) == '_FillValue'   ) cycle
      if ( trim(aname) == 'add_offset'   ) cycle
      if ( trim(aname) == 'scale_factor' ) cycle
      ! extra element:
      call NcAttrs_Append_Empty( self, status )
      IF_NOT_OK_RETURN(status=1)
      ! fill from file:
      call self%values(self%n)%p%NcGet( ncid, varid, aname, status, &
                                         read_on_root=self%read_on_root, only_here=only_here )
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


  !
  ! Initialize file with provided name for reading or writing.
  ! - rwmode='r'
  !     The filename should exist, and the file is opened.
  ! - rwmode='o'
  !     File is alrady opened; filename is ignored, netcdf id should
  !     be provided using the 'ncid' argument.
  !     This form is to be replaced by using T_NcFile everywhere.
  ! - rwmode='w'
  !     File is created.
  !     Storage for variables and attributes is initialized.
  !     Calls to 'Def_Dim' and 'Def_Var' fill the storage.
  !     A call to 'End_Def' will finally define the dimensions and variables;
  !     the variables need to be filled by calls to 'Put' routines.
  !

  subroutine NcFile_Init( self, filename, rwmode, status, ncid )
  
    use NetCDF, only : NF90_Open
    use NetCDF, only : NF90_Create
    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_NOCLOBBER, NF90_CLOBBER
    use NetCDF, only : NF90_CLASSIC_MODEL, NF90_NETCDF4
    use NetCDF, only : NF90_Inquire

    use CSO_Comm, only : csoc
    use CSO_File, only : CSO_CheckDir
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(out)    ::  self
    character(len=*), intent(in)    ::  filename
    character(len=1), intent(in)    ::  rwmode
    integer, intent(out)            ::  status
    integer, intent(in), optional   ::  ncid

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Init'
    
    ! --- local ----------------------------------
    
    logical      ::  exist    
    integer      ::  cmode
    
    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! store mode:
    self%rwmode = rwmode
    
    ! init dimension list:
    call self%dims%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init variable list:
    call self%vars%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init attributes list:
    call self%attrs%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! switch:
    select case ( self%rwmode )
      ! read:
      case ( 'r' )
      
        ! check ..
        if ( present(ncid) ) then
          write (csol,'("unsupported argument `ncid` for rwmode `",a,"`")') self%rwmode; call csoErr
          TRACEBACK; status=1; return
        end if
    
        ! read on root...
        if ( csoc%root ) then
    
          ! check ..
          inquire( file=trim(self%filename), exist=exist )
          if ( .not. exist ) then
            write (csol,'("WARNING - file not found : ",a)') trim(self%filename); call csoPr
            status=-1; return
          end if

          !! Export filename to log
          !write (csol,'("Attempting to read from file: ", a)') trim(filename); call csoPr

          ! open file for reading:
          status = NF90_Open( trim(filename), NF90_NOWRITE, self%ncid )
          IF_NF90_NOT_OK_RETURN(status=1)

          ! info ...
          status = NF90_Inquire( self%ncid, formatNum=self%formatNum )
          IF_NF90_NOT_OK_RETURN(status=1)
          
        end if ! root
        
      ! already open:
      case ( 'o' )

        ! check ..
        if ( .not. present(ncid) ) then
          write (csol,'("no argument `ncid` provided for rwmode `",a,"`")') self%rwmode; call csoErr
          TRACEBACK; status=1; return
        end if
    
        ! store:
        self%ncid = ncid
        
      ! write:
      case ( 'w' )

        ! check ..
        if ( present(ncid) ) then
          write (csol,'("unsupported argument `ncid` for rwmode `",a,"`")') self%rwmode; call csoErr
          TRACEBACK; status=1; return
        end if
    
        ! written on root...
        if ( csoc%root ) then
        
          ! create directory if necessary:
          call CSO_CheckDir( filename, status )
          IF_NOT_OK_RETURN(status=1)

          ! set creation mode flag:
          cmode = 0
          !~ overwite?
          cmode = cmode + NF90_CLOBBER       ! overwrite existing files
          !cmode = cmode + NF90_NOCLOBBER     ! do not overwrite existing files
          !~ netcdf flavour:
          !cmode = cmode + NF90_CLASSIC_MODEL
          cmode = cmode + NF90_NETCDF4

          ! create file:
          status = NF90_Create( self%filename, cmode, self%ncid )
          if ( status /= NF90_NOERR ) then
             write (csol,'("creating file :")'); call csoErr
             write (csol,'("  ",a)') trim(self%filename); call csoErr
             TRACEBACK; status=1; return
          end if
          
        end if ! root

      ! unknown
      case default
        write (csol,'("unsupported rwmode `",a,"`")') self%rwmode; call csoErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
  
  end subroutine NcFile_Init


  ! *
  

  subroutine NcFile_Done( self, status )

    use NetCDF, only : NF90_Close
  
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! switch:
    select case ( self%rwmode )

      ! read, open:
      case ( 'r', 'o' )
        ! nothing to be done
        
      ! write:
      case ( 'w' )

        ! written on root...
        if ( csoc%root ) then
          ! close:
          status = NF90_Close( self%ncid )
          IF_NF90_NOT_OK_RETURN(status=1)
        end if  ! root

      ! unknown
      case default
        write (csol,'("unsupported rwmode `",a,"`")') self%rwmode; call csoErr
        TRACEBACK; status=1; return
    end select

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


  ! ***
  
  
  subroutine NcFile_Inquire( self, status, ndim, nvar )

    use NetCDF, only : NF90_Inquire
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)     ::  self
    integer, intent(out)            ::  status
    integer, intent(out), optional  ::  ndim
    integer, intent(out), optional  ::  nvar

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Inquire'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! obtain info:
    status = NF90_Inquire( self%ncid, nDimensions=ndim, nVariables=nvar )
    IF_NF90_NOT_OK_RETURN(status=1)   
    
    ! ok
    status = 0
    
  end subroutine NcFile_Inquire
  
  ! *

  !
  ! Return id of variable identified by it's name or an
  ! attribute value. The description is a ';' seperated list 
  ! of 'key=value' pairs, where 'var_name' should be used
  ! for the variable name and all other names are assumed
  ! to be attributes:
  !
  !    standard_name=pressure;long_name=air pressure;var_name=p
  !
  ! Duplicates are allowed, for example:
  !
  !    long_name=air pressure;long_name=Pressure
  !
  ! If a match is found the variable id is returned.
  !
  
  subroutine NcFile_Inq_VarID( self, description, varid, status )

    use NetCDF, only : NF90_INQ_VarID
    use NetCDF, only : NF90_Get_Att
    
    use CSO_String, only : CSO_ReadFromLine, CSO_NtsTrim
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    integer, intent(out)              ::  varid
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Inq_VarID'
    
    ! --- local ----------------------------------
    
    logical             ::  found
    integer             ::  nvar, ivar
    character(len=1024) ::  line
    character(len=64)   ::  name
    character(len=256)  ::  value
    character(len=256)  ::  attr_value
    
    ! --- begin ----------------------------------
    
    ! set flag:
    found = .false.
    
    ! copy:
    line = description
    ! loop over parts:
    do
      ! empty?
      if ( len_trim(line) == 0 ) exit

      ! split part:
      call CSO_ReadFromLine( line, value, status, sep=';' )
      IF_NOT_OK_RETURN(status=1)
      ! split-off first part:
      call CSO_ReadFromLine( value, name, status, sep='=' )
      IF_NOT_OK_RETURN(status=1)

      ! switch:
      select case ( name )
        !~
        case ( 'var_name' )
          ! if variable is present with name equal to the value,
          ! then the variable id is returned and status is not error;
          ! otherwise (no variable with this name) an error status is returned:
          status = NF90_INQ_VarID( self%ncid, trim(value), varid )
          found = status == NF90_NOERR
        !~
        case default
          ! number of variables:
          call self%Inquire( status, nvar=nvar )
          IF_NOT_OK_RETURN(status=1)
          ! loop:
          do ivar = 1, nvar
            ! variable id is number:
            varid = ivar
            ! get name if present, error status if not exists:
            status = NF90_Get_Att( self%ncid, varid, trim(name), attr_value )
            if ( status /= NF90_NOERR ) cycle
            ! fix nul-terminated strings:
            attr_value = CSO_NtsTrim( attr_value )
            ! match?
            found = trim(attr_value) == trim(value)
            ! leave ?
            if ( found ) exit
          end do
        !~
      end select
      
      ! leave?
      if ( found ) exit
      
    end do  ! key=value pairs
      
    ! check ...
    if ( .not. found ) then
      write (csol,'("no variable found matching description: ",a)') trim(description); call csoErr
      write (csol,'("  file: ",a)') trim(self%filename); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine NcFile_Inq_VarID
    
  ! *
  
  !
  ! Return variable shape in 'shp(:)'.
  ! Size of 'shp' should be at least number of dimensions,
  ! first dimensions are padded with length 1 if the file
  ! has less dimensions than the requested shape.
  !

  subroutine NcFile_Inq_Variable( self, varid, status, &
                                    ndims, shp )
  
    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable

    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    integer, intent(in)               ::  varid
    integer, intent(out)              ::  status
    integer, intent(out), optional    ::  ndims
    integer, intent(out), optional    ::  shp(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Inq_Variable'
    
    ! --- local ----------------------------------

    integer                 ::  nshp
    integer                 ::  nd
    integer, allocatable    ::  dimids(:)
    integer                 ::  idim
    
    ! --- begin ----------------------------------
    
    ! return number of dims?
    if ( present(ndims) ) then
      ! number of dimensions:
      status = NF90_Inquire_Variable( self%ncid, varid, ndims=ndims )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if ! ndim
    
    ! return shape?
    if ( present(shp) ) then

      ! number of dims to be filled:
      nshp = size(shp)

      ! number of dimensions:
      status = NF90_Inquire_Variable( self%ncid, varid, ndims=nd )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! check ..
      if ( nshp < nd ) then
        write (csol,'("argument shp has size ",i0," but variable has ",i0," dimensions")') &
                size(shp), nd; call csoErr
        TRACEBACK; status=1; return
      end if

      ! storage:
      allocate( dimids(nd), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! get dimension id's:
      status = NF90_Inquire_Variable( self%ncid, varid, dimids=dimids )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! loop over dimensions:
      do idim = 1, size(shp)
        ! scalar?
        if ( idim <= nshp-nd ) then
          shp(idim) = 1
        else
          ! get size:
          status = NF90_Inquire_Dimension( self%ncid, dimids(idim-(nshp-nd)), len=shp(idim) )
          IF_NF90_NOT_OK_RETURN(status=1)
        end if
      end do

      ! clear:
      deallocate( dimids, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! shp
    
    ! ok
    status = 0

  end subroutine NcFile_Inq_Variable
  
    
  ! *
  
  ! Return packing parameters 'add_offset' and 'scale_factor'.
  ! Sometimes files are found that have only one of these,
  ! the missing parameter is then set to 0.0 or 1.0 respectively.
  ! Return status -1 if none of the parameters is found.
  !

  subroutine NcFile_Inq_VarPacking( self, varid, add_offset, scale_factor, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    integer, intent(in)               ::  varid
    real, intent(out)                 ::  add_offset, scale_factor
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Inq_VarPacking'
    
    ! --- local ----------------------------------

    logical     ::  any_present
    
    ! --- begin ----------------------------------
    
    ! by default nothing found:
    any_present = .false.
    
    ! try to get packing variable, ok if not present:
    status = NF90_Get_Att( self%ncid, varid, 'add_offset', add_offset )
    if ( status == NF90_NOERR ) then
      ! reset flag:
      any_present = .true.
    else if ( status == NF90_ENOTATT ) then
      ! default value:
      add_offset = 0.0
    else
      csol=NF90_StrError(status); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! try to get packing variable, ok if not present:
    status = NF90_Get_Att( self%ncid, varid, 'scale_factor', scale_factor )
    if ( status == NF90_NOERR ) then
      ! reset flag:
      any_present = .true.
    else if ( status == NF90_ENOTATT ) then
      ! default value:
      scale_factor = 1.0
    else
      csol=NF90_StrError(status); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! any found?
    if ( any_present ) then
      status = 0
    else
      status = -1
    end if
    
  end subroutine NcFile_Inq_VarPacking


  ! *
  
  
  ! Return missing value, status -1 if not found.

  subroutine NcFile_Inq_VarMissing( self, varid, missing_value, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    integer, intent(in)               ::  varid
    real, intent(out)                 ::  missing_value
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Inq_VarMissing'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! try to get first packing variable:
    status = NF90_Get_Att( self%ncid, varid, 'missing_value', missing_value )
    ! no error, missing value available
    if ( status == NF90_NOERR ) then
      ! missing value available,
      ! remain in output argument
      
    ! attribute not found ?
    else if ( status == NF90_ENOTATT ) then
      ! no missing value, set default values:
      missing_value = -999.9
      ! warning status:
      status = -1; return
    !
    else
      ! some error ...
      csol=NF90_StrError(status); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine NcFile_Inq_VarMissing
  
  
  ! *
  
  
  !
  ! Return units from attribute of variabled opened with 'varid'.
  ! If attribute is not present, try to extract from 'description'
  ! that was used to obtain the variable id:
  !    description='var_name=pressure;units=Pa'
  !

  subroutine NcFile_Inq_VarUnits( self, varid, description, units, status )

    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
    
    use CSO_String, only : CSO_VarValue, CSO_NtsTrim
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    integer, intent(in)               ::  varid
    character(len=*), intent(in)      ::  description
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Inq_VarUnits'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! get units from attribute:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    !~ attribute found and read:
    if ( status == NF90_NOERR ) then
      ! remove nul-characters:
      units = CSO_NtsTrim( units )
      ! empty?
      if ( len_trim(units) == 0 ) then
        ! try to read from description, warning status<0 if not defined:
        call CSO_VarValue( description, ';', 'units', '=', units, status )
        IF_ERROR_RETURN(status=1)
        ! not found?
        if ( status < 0 ) then
          write (csol,'("empty units attribute in variable, and no explicit specification in description either:")'); call csoErr
          write (csol,'("  filename             : ",a)') trim(self%filename); call csoErr
          write (csol,'("  variable description : ",a)') trim(description); call csoErr
          TRACEBACK; status=1; return
        end if
      end if
    !~ attribute not found ...
    else if ( status == NF90_ENOTATT ) then
      ! try to read from description, warning status<0 if not defined:
      call CSO_VarValue( description, ';', 'units', '=', units, status )
      IF_ERROR_RETURN(status=1)
      ! not found?
      if ( status < 0 ) then
        write (csol,'("no units attribute in variable, and no explicit specification in description either:")'); call csoErr
        write (csol,'("  filename             : ",a)') trim(self%filename); call csoErr
        write (csol,'("  variable description : ",a)') trim(description); call csoErr
        TRACEBACK; status=1; return
      end if
    !~ other error ...
    else
      ! some error ...
      csol=NF90_StrError(status); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine NcFile_Inq_VarUnits
 

  !  ***


  subroutine NcFile_Get_Var_i1_1d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    integer(1), intent(out)           ::  values(:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Get_Var_i1_1d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call NcFile_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = nint( add_offset + scale_factor * values )
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcFile_Get_Var_i1_1d
 

  !  ***


  subroutine NcFile_Get_Var_i_1d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    integer, intent(out)              ::  values(:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Get_Var_i_1d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call NcFile_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = nint( add_offset + scale_factor * values )
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcFile_Get_Var_i_1d
  
  ! *
  
  subroutine NcFile_Get_Var_c_2d( self, description, values, units, status, &
                                   start, count )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    character(len=1), intent(out)     ::  values(:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Get_Var_c_2d'
    
    ! --- local ----------------------------------
    
    integer                         ::  varid
    character(len=:), allocatable   ::  cvalues(:)
    integer                         ::  i, j
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call NcFile_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! reading 2D char array does not work,
    ! use instead 1D array of strings:
    allocate( character(len=size(values,1)) :: cvalues(size(values,2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! safety ..
    if ( present(start) .or. present(count) ) then
      write (csol,'("optional arguments `start` or `count` not supported yet for char arrays")'); call csoErr
      TRACEBACK; status=1; return
    end if

    ! read:
    status = NF90_Get_Var( self%ncid, varid, cvalues )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! copy:
    do j = 1, size(values,2)
      do i = 1, size(values,1)
        values(i,j) = cvalues(j)(i:i)
      end do ! i
    end do ! j
    
    ! clear:
    deallocate( cvalues, stat=status )
    IF_NOT_OK_RETURN(status=1)
                
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcFile_Get_Var_c_2d
  
  ! *
  
  subroutine NcFile_Get_Var_i_2d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    integer, intent(out)              ::  values(:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Get_Var_i_2d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call NcFile_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = nint( add_offset + scale_factor * values )
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcFile_Get_Var_i_2d
  
  
  ! *
  
  
  subroutine NcFile_Get_Var_i_3d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    integer, intent(out)              ::  values(:,:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Get_Var_i_3d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call NcFile_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = int( add_offset + scale_factor * values )
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcFile_Get_Var_i_3d


  ! *
  

  subroutine NcFile_Get_Var_r_1d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    real, intent(out)                 ::  values(:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Get_Var_r_1d'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    real                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call NcFile_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)

    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcFile_Get_Var_r_1d
  
  ! *
  
  subroutine NcFile_Get_Var_r_2d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    real, intent(out)                 ::  values(:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Get_Var_r_2d'
    
    ! --- local ----------------------------------
    
    logical                 ::  combine
    integer, allocatable    ::  xstart(:)
    integer, allocatable    ::  xcount(:)
    integer                 ::  x1, nx
    integer                 ::  varid
    integer, allocatable    ::  dimids(:)
    real                    ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call NcFile_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! check start:
    if ( any((/present(start),present(count)/)) .and. &
         (.not. all((/present(start),present(count)/))) ) then
      write (csol,'("specify both start and count")'); call csoErr
      TRACEBACK; status=1; return
    end if

    ! combine slabs?
    combine = .false.
    if ( present(start) ) combine = start(1) < 1
    ! switch:
    if ( combine ) then

      ! storage:
      allocate( xstart(size(start)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( xcount(size(count)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dimids(size(count)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy:
      xstart = start
      xcount = count

      ! start index:
      x1 = xstart(1)
    
      ! get dimension id's:
      status = NF90_Inquire_Variable( self%ncid, varid, dimids=dimids )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! get input dimension length:
      status = NF90_Inquire_Dimension( self%ncid, dimids(1), len=nx )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! set input start and count for first slab:
      xstart(1) = nx + x1
      xcount(1) = nx - xstart(1) + 1
      ! read first slab:
      status = NF90_Get_Var( self%ncid, varid, values(1:xcount(1),:), &
                               start=xstart, count=xcount )
      if ( status /= NF90_NOERR ) then
        csol=NF90_StrError(status); call csoErr
        write (csol,'("while reading:")'); call csoErr
        write (csol,'("  file name              : ",a)') trim(self%filename); call csoErr
        write (csol,'("  variable description   : ",a)') trim(description); call csoErr
        write (csol,*) ' start                  : ', xstart; call csoErr
        write (csol,*) ' count                  : ', xcount; call csoErr
        TRACEBACK; status=1; return
      end if

      ! set input start and count for second slab:
      xstart(1) = 1
      xcount(1) = size(values,1) + x1 - 1
      ! read first slab:
      status = NF90_Get_Var( self%ncid, varid, values(2-x1:size(values,1),:), &
                               start=xstart, count=xcount )
      if ( status /= NF90_NOERR ) then
        csol=NF90_StrError(status); call csoErr
        write (csol,'("while reading:")'); call csoErr
        write (csol,'("  file name              : ",a)') trim(self%filename); call csoErr
        write (csol,'("  variable description   : ",a)') trim(description); call csoErr
        write (csol,*) ' start                  : ', xstart; call csoErr
        write (csol,*) ' count                  : ', xcount; call csoErr
        TRACEBACK; status=1; return
      end if
      
      ! clear:
      deallocate( xstart, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xcount, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dimids, stat=status )
      IF_NOT_OK_RETURN(status=1)

    else

      ! read:
      status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
      if ( status /= NF90_NOERR ) then
        csol=NF90_StrError(status); call csoErr
        write (csol,'("while reading:")'); call csoErr
        write (csol,'("  file name              : ",a)') trim(self%filename); call csoErr
        write (csol,'("  variable description   : ",a)') trim(description); call csoErr
        if ( present(start) ) then
          write (csol,*) ' start                  : ', start; call csoErr
        end if
        if ( present(count) ) then
          write (csol,*) ' count                  : ', count; call csoErr
        end if
        TRACEBACK; status=1; return
      end if
      
    end if
    
    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcFile_Get_Var_r_2d
  
  ! *
  
  subroutine NcFile_Get_Var_r_3d( self, description, values, units, status, &
                                   start, count, missing_value )

    use NetCDF, only : NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_ENOTATT
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(in)       ::  self
    character(len=*), intent(in)      ::  description
    real, intent(out)                 ::  values(:,:,:)
    character(len=*), intent(out)     ::  units
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  start(:), count(:)
    real, intent(out), optional       ::  missing_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NcFile_Get_Var_r_3d'
    
    ! --- local ----------------------------------
    
    logical                 ::  combine
    integer, allocatable    ::  xstart(:)
    integer, allocatable    ::  xcount(:)
    integer                 ::  x1, nx
    integer                 ::  varid
    integer, allocatable    ::  dimids(:)
    real                    ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! get variable id:
    call NcFile_Inq_VarID( self, description, varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! check start:
    if ( any((/present(start),present(count)/)) .and. &
         (.not. all((/present(start),present(count)/))) ) then
      write (csol,'("specify both start and count")'); call csoErr
      TRACEBACK; status=1; return
    end if

    ! combine slabs?
    combine = .false.
    if ( present(start) ) combine = start(1) < 1
    ! switch:
    if ( combine ) then

      ! storage:
      allocate( xstart(size(start)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( xcount(size(count)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dimids(size(count)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy:
      xstart = start
      xcount = count

      ! start index:
      x1 = xstart(1)
    
      ! get dimension id's:
      status = NF90_Inquire_Variable( self%ncid, varid, dimids=dimids )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! get input dimension length:
      status = NF90_Inquire_Dimension( self%ncid, dimids(1), len=nx )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! set input start and count for first slab:
      xstart(1) = nx + x1
      xcount(1) = nx - xstart(1) + 1
      ! read first slab:
      status = NF90_Get_Var( self%ncid, varid, values(1:xcount(1),:,:), &
                               start=xstart, count=xcount )
      if ( status /= NF90_NOERR ) then
        csol=NF90_StrError(status); call csoErr
        write (csol,'("while reading:")'); call csoErr
        write (csol,'("  file name              : ",a)') trim(self%filename); call csoErr
        write (csol,'("  variable description   : ",a)') trim(description); call csoErr
        write (csol,*) ' start                  : ', xstart; call csoErr
        write (csol,*) ' count                  : ', xcount; call csoErr
        TRACEBACK; status=1; return
      end if

      ! set input start and count for second slab:
      xstart(1) = 1
      xcount(1) = size(values,1) + x1 - 1
      ! read first slab:
      status = NF90_Get_Var( self%ncid, varid, values(2-x1:size(values,1),:,:), &
                               start=xstart, count=xcount )
      if ( status /= NF90_NOERR ) then
        csol=NF90_StrError(status); call csoErr
        write (csol,'("while reading:")'); call csoErr
        write (csol,'("  file name              : ",a)') trim(self%filename); call csoErr
        write (csol,'("  variable description   : ",a)') trim(description); call csoErr
        write (csol,*) ' start                  : ', xstart; call csoErr
        write (csol,*) ' count                  : ', xcount; call csoErr
        TRACEBACK; status=1; return
      end if
      
      ! clear:
      deallocate( xstart, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xcount, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dimids, stat=status )
      IF_NOT_OK_RETURN(status=1)

    else
    
      ! read:
      status = NF90_Get_Var( self%ncid, varid, values, start=start, count=count )
      if ( status /= NF90_NOERR ) then
        csol=NF90_StrError(status); call csoErr
        write (csol,'("while reading:")'); call csoErr
        write (csol,'("  file name              : ",a)') trim(self%filename); call csoErr
        write (csol,'("  variable description   : ",a)') trim(description); call csoErr
        if ( present(start) ) then
          write (csol,*) ' start                  : ', start; call csoErr
        end if
        if ( present(count) ) then
          write (csol,*) ' count                  : ', count; call csoErr
        end if
        TRACEBACK; status=1; return
      end if
      
    end if

    ! packed?
    call self%Inq_VarPacking( varid, add_offset, scale_factor, status )
    IF_ERROR_RETURN(status=1)
    if ( status == 0 ) then
      ! unpack:
      values = add_offset + scale_factor * values
    end if
    
    ! Missing value?
    if ( present( missing_value ) ) then
      call self%Inq_VarMissing( varid, missing_value, status )
      IF_ERROR_RETURN(status=1)
    end if
    
    ! get units:
    call self%Inq_VarUnits( varid, description, units, status )
    IF_ERROR_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NcFile_Get_Var_r_3d
  
  
  ! ***
  
  
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
  
    use NetCDF  , only : NF90_GLOBAL
    use NetCDF  , only : NF90_EndDef

    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_EndDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! written on root...
    if ( csoc%root ) then

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
  
  
  !
  ! Compute packing parameters given range of values.
  ! When packed into short integers, 
  ! the actual values are recalculated using the formula:
  !
  !   values = add_offset + scale_factor * packed_values
  !
  ! The packing parameters are computed using:
  !
  !   scale_factor = (vmax - vmin)/( fmax - fmin)
  !   add_offset   = vmin - scale_factor * fmin
  !
  ! such that the original range is:
  !
  !   vmin = add_offset + scale_factor * fmin
  !        = vmin - (vmax - vmin)/( fmax - fmin) * fmin + (vmax - vmin)/( fmax - fmin) * fmin
  !        = vmin
  !
  !   vmax = add_offset + scale_factor * fmax
  !        = vmin - (vmax - vmin)/( fmax - fmin) * fmin + (vmax - vmin)/( fmax - fmin) * fmax
  !        = vmin + (vmax - vmin)/( fmax - fmin) * (fmax-fmin)
  !        = vmax
  !
  ! where ``[vmin,vmax]`` is the range of input values,
  ! and ``[fmin,fmax]`` the range of possible values of the packed data type.
  !
  ! If only a single value is present (vmin==vmax) then is used:
  !   scale_factor = 1.0
  !   add_offset   = vmax
  !
  ! Original values might have no-data elements equal to 'fill_value'.
  ! If this is defined, then also 'fill_value__packed' should have been defined;
  ! its value should be outside [fmin,vmin], and is here checked to be fmin-1.
  !

  subroutine NcFile_GetPacking( self, vmin, vmax, &
                                  add_offset, scale_factor, status, &
                                  fill_value, fill_value__packed )
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    real, intent(in)                      ::  vmin, vmax
    real, intent(out)                     ::  add_offset, scale_factor
    integer, intent(out)                  ::  status
    
    real, intent(in), optional                    ::  fill_value
    integer(iwp__packed), intent(in), optional    ::  fill_value__packed

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_GetPacking'
    
    ! --- local ----------------------------------
    
    integer(iwp__packed)                ::  fmin, fmax
    
    ! --- begin ----------------------------------
    
    ! range of packed values, -huge(1) is reserved for fill_value ..
    fmin = -huge(int(1,kind=iwp__packed)) + 1
    fmax =  huge(int(1,kind=iwp__packed))
    
    ! any range?
    if ( vmin < vmax ) then
      ! packing parameters:
      scale_factor = (vmax - vmin)/( real(fmax) - real(fmin) )
      add_offset   = vmin - scale_factor * fmin
    else
      ! single value only;
      ! set offset to that (single) value and scale factor to 1.0,
      ! packed values will all be 0.0:
      scale_factor = 1.0
      add_offset   = vmax
    end if  ! range or single value

    ! need to check on no-data values?
    if ( present(fill_value) ) then
      ! check ..
      if ( .not. present(fill_value__packed) ) then
        write (csol,'("argument fill_value present but fill_value__packed not")'); call csoErr
        TRACEBACK; status=1; return
      end if
      ! check ..
      if ( fill_value__packed /= fmin-1 ) then
        write (csol,*) 'fill_value__packed is ', fill_value__packed, ' while expected ', fmin-1; call csoErr
        TRACEBACK; status=1; return
      end if
    end if  ! check on fill_value
    
    ! ok
    status = 0
    
  end subroutine NcFile_GetPacking
  
  
  ! *
  
  
  !
  ! (Adhoc routine, does not use the internal variable list yet ...)
  !
  ! Write 1D real data to variable with netcdf id 'varid'.
  ! Eventually pack as short integers.
  !

  subroutine NcFile_Put_Var_1d_r( self, varid, values, status, &
                                    fill_value, fill_value__out, &
                                    packed, fill_value__packed )
  
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Put_Att

    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(in)                   ::  varid
    real, intent(in)                      ::  values(:)
    integer, intent(out)                  ::  status
    
    real, intent(in), optional                    ::  fill_value
    real(rwp__out), intent(in), optional          ::  fill_value__out
    logical, intent(in), optional                 ::  packed
    integer(iwp__packed), intent(in), optional    ::  fill_value__packed

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Put_Var_1d_r'
    
    ! --- local ----------------------------------
    
    real(rwp__out), allocatable         ::  values__out(:)
    logical                             ::  with_packed
    integer(iwp__packed), allocatable   ::  values__packed(:)
    real                                ::  vmin, vmax
    real                                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! flag:
    with_packed = .false.
    if ( present(packed) ) with_packed = packed
    
    ! packed?
    if ( with_packed ) then
    
      ! The actual values are recalculated using the formula:
      !   values = add_offset + scale_factor * packed_values
      ! The packing parameters are computed using::
      !   scale_factor = (vmax - vmin)/( fmax - fmin)
      !   add_offset   = vmin - scale_factor * fmin
      ! where ``[vmin,vmax]`` is the range of input values,
      ! and ``[fmin,fmax]`` the range of possible values of the packed data type.

      ! storage:
      allocate( values__packed(size(values)), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! value range:
      if ( present(fill_value) ) then
        ! range, only values other than fill value:
        vmin = minval( values, mask=values/=fill_value )
        vmax = maxval( values, mask=values/=fill_value )
      else
        ! range:
        vmin = minval( values )
        vmax = maxval( values )
      end if
      
      ! get packing parameters, check fill values
      call self%GetPacking( vmin, vmax, add_offset, scale_factor, status, &
                             fill_value=fill_value, &
                             fill_value__packed=fill_value__packed )
      IF_NOT_OK_RETURN(status=1)

      ! need to check on no-data values?
      if ( present(fill_value) ) then
        ! encode into packed variable, except for no-data values:
        where ( values == fill_value )
          values__packed = fill_value__packed
        elsewhere
          values__packed = int( ( values - add_offset )/scale_factor, kind=iwp__packed )
        end where
      else
        ! encode into packed variable:
        values__packed = int( ( values - add_offset )/scale_factor, kind=iwp__packed )
      end if  ! check on fill_value

      ! write packed values:
      status = NF90_Put_Var( self%ncid, varid, values__packed )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! write packing attributes:
      status = NF90_Put_Att( self%ncid, varid, 'add_offset', add_offset )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( self%ncid, varid, 'scale_factor', scale_factor )
      IF_NF90_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( values__packed, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    else

      ! storage:
      allocate( values__out(size(values)), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! convert to output precission:
      if ( any((/present(fill_value),present(fill_value__out)/)) ) then
        ! check ...
        if ( .not. all((/present(fill_value),present(fill_value__out)/)) ) then
          write (csol,'("either none or both arguments `fill_value` and `` should be present")'); call csoErr
          TRACEBACK; status=1; return
        end if
        ! copy to output precission with change of fill values:
        where ( values == fill_value )
          values__out = fill_value__out
        elsewhere
          values__out = real(values,kind=rwp__out)
        end where
      else
         ! copy:
         values__out = real(values,kind=rwp__out)
      end if
    
      ! write:
      status = NF90_Put_Var( self%ncid, varid, values__out )
      IF_NF90_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( values__out, stat=status )
      IF_NOT_OK_RETURN(status=1)
    
    end if
    
    ! ok
    status = 0
    
  end subroutine NcFile_Put_Var_1d_r
  
  
  ! *
  
  
  !
  ! (Adhoc routine, does not use the internal variable list yet ...)
  !
  ! Write 2D char data to variable with netcdf id 'varid'.
  !

  subroutine NcFile_Put_Var_2d_c( self, varid, values, status )
  
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Put_Att

    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(in)                   ::  varid
    character(len=*), intent(in)          ::  values(:,:)
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Put_Var_2d_c'
    
    ! --- local ----------------------------------

    character(len=:), allocatable   ::  cvalues(:)
    integer                         ::  i, j
    
    ! --- begin ----------------------------------

    ! writing 2D char array does not work,
    ! use instead 1D array of strings:
    allocate( character(len=size(values,1)) :: cvalues(size(values,2)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy:
    do j = 1, size(values,2)
      do i = 1, size(values,1)
        cvalues(j)(i:i) = values(i,j)
      end do ! i
    end do ! j
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, cvalues )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( cvalues, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine NcFile_Put_Var_2d_c
  
  
  ! *
  
  
  !
  ! (Adhoc routine, does not use the internal variable list yet ...)
  !
  ! Write 2D real data to variable with netcdf id 'varid'.
  ! Eventually pack as short integers.
  !

  subroutine NcFile_Put_Var_2d_r( self, varid, values, status, &
                                    fill_value, fill_value__out, &
                                    packed, fill_value__packed )
  
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Put_Att

    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(in)                   ::  varid
    real, intent(in)                      ::  values(:,:)
    integer, intent(out)                  ::  status
    
    real, intent(in), optional                    ::  fill_value
    real(rwp__out), intent(in), optional          ::  fill_value__out
    logical, intent(in), optional                 ::  packed
    integer(iwp__packed), intent(in), optional    ::  fill_value__packed

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Put_Var_2d_r'
    
    ! --- local ----------------------------------
    
    real(rwp__out), allocatable         ::  values__out(:,:)
    logical                             ::  with_packed
    integer(iwp__packed), allocatable   ::  values__packed(:,:)
    real                                ::  vmin, vmax
    real                                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! flag:
    with_packed = .false.
    if ( present(packed) ) with_packed = packed
    
    ! packed?
    if ( with_packed ) then
    
      ! The actual values are recalculated using the formula:
      !   values = add_offset + scale_factor * packed_values
      ! The packing parameters are computed using::
      !   scale_factor = (vmax - vmin)/( fmax - fmin)
      !   add_offset   = vmin - scale_factor * fmin
      ! where ``[vmin,vmax]`` is the range of input values,
      ! and ``[fmin,fmax]`` the range of possible values of the packed data type.

      ! storage:
      allocate( values__packed(size(values,1),size(values,2)), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! value range:
      if ( present(fill_value) ) then
        ! range, only values other than fill value:
        vmin = minval( values, mask=values/=fill_value )
        vmax = maxval( values, mask=values/=fill_value )
      else
        ! range:
        vmin = minval( values )
        vmax = maxval( values )
      end if
      
      ! get packing parameters, check fill values
      call self%GetPacking( vmin, vmax, add_offset, scale_factor, status, &
                             fill_value=fill_value, &
                             fill_value__packed=fill_value__packed )
      IF_NOT_OK_RETURN(status=1)

      ! need to check on no-data values?
      if ( present(fill_value) ) then
        ! encode into packed variable, except for no-data values:
        where ( values == fill_value )
          values__packed = fill_value__packed
        elsewhere
          values__packed = int( ( values - add_offset )/scale_factor, kind=iwp__packed )
        end where
      else
        ! encode into packed variable:
        values__packed = int( ( values - add_offset )/scale_factor, kind=iwp__packed )
      end if  ! check on fill_value

      ! write packed values:
      status = NF90_Put_Var( self%ncid, varid, values__packed )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! write packing attributes:
      status = NF90_Put_Att( self%ncid, varid, 'add_offset', add_offset )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( self%ncid, varid, 'scale_factor', scale_factor )
      IF_NF90_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( values__packed, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    else

      ! storage:
      allocate( values__out(size(values,1),size(values,2)), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! convert to output precission:
      if ( any((/present(fill_value),present(fill_value__out)/)) ) then
        ! check ...
        if ( .not. all((/present(fill_value),present(fill_value__out)/)) ) then
          write (csol,'("either none or both arguments `fill_value` and `` should be present")'); call csoErr
          TRACEBACK; status=1; return
        end if
        ! copy to output precission with change of fill values:
        where ( values == fill_value )
          values__out = fill_value__out
        elsewhere
          values__out = real(values,kind=rwp__out)
        end where
      else
         ! copy:
         values__out = real(values,kind=rwp__out)
      end if
    
      ! write:
      status = NF90_Put_Var( self%ncid, varid, values__out )
      IF_NF90_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( values__out, stat=status )
      IF_NOT_OK_RETURN(status=1)
    
    end if
    
    ! ok
    status = 0
    
  end subroutine NcFile_Put_Var_2d_r
  
  
  ! *
  
  
  !
  ! (Adhoc routine, does not use the internal variable list yet ...)
  !
  ! Write 2D real data to variable with netcdf id 'varid'.
  ! Eventually pack as short integers.
  !

  subroutine NcFile_Put_Var_3d_r( self, varid, values, status, &
                                    fill_value, fill_value__out, &
                                    packed, fill_value__packed )
  
    use NetCDF, only : NF90_Put_Var
    use NetCDF, only : NF90_Put_Att

    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(in)                   ::  varid
    real, intent(in)                      ::  values(:,:,:)
    integer, intent(out)                  ::  status
    
    real, intent(in), optional                    ::  fill_value
    real(rwp__out), intent(in), optional          ::  fill_value__out
    logical, intent(in), optional                 ::  packed
    integer(iwp__packed), intent(in), optional    ::  fill_value__packed

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Put_Var_3d_r'
    
    ! --- local ----------------------------------
    
    real(rwp__out), allocatable         ::  values__out(:,:,:)
    logical                             ::  with_packed
    integer(iwp__packed), allocatable   ::  values__packed(:,:,:)
    real                                ::  vmin, vmax
    real                                ::  add_offset, scale_factor
    
    ! --- begin ----------------------------------
    
    ! flag:
    with_packed = .false.
    if ( present(packed) ) with_packed = packed
    
    ! packed?
    if ( with_packed ) then
    
      ! The actual values are recalculated using the formula:
      !   values = add_offset + scale_factor * packed_values
      ! The packing parameters are computed using::
      !   scale_factor = (vmax - vmin)/( fmax - fmin)
      !   add_offset   = vmin - scale_factor * fmin
      ! where ``[vmin,vmax]`` is the range of input values,
      ! and ``[fmin,fmax]`` the range of possible values of the packed data type.

      ! storage:
      allocate( values__packed(size(values,1),size(values,2),size(values,3)), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! value range:
      if ( present(fill_value) ) then
        ! range, only values other than fill value:
        vmin = minval( values, mask=values/=fill_value )
        vmax = maxval( values, mask=values/=fill_value )
      else
        ! range:
        vmin = minval( values )
        vmax = maxval( values )
      end if
      
      ! get packing parameters, check fill values
      call self%GetPacking( vmin, vmax, add_offset, scale_factor, status, &
                             fill_value=fill_value, &
                             fill_value__packed=fill_value__packed )
      IF_NOT_OK_RETURN(status=1)

      ! need to check on no-data values?
      if ( present(fill_value) ) then
        ! encode into packed variable, except for no-data values:
        where ( values == fill_value )
          values__packed = fill_value__packed
        elsewhere
          values__packed = int( ( values - add_offset )/scale_factor, kind=iwp__packed )
        end where
      else
        ! encode into packed variable:
        values__packed = int( ( values - add_offset )/scale_factor, kind=iwp__packed )
      end if  ! check on fill_value

      ! write packed values:
      status = NF90_Put_Var( self%ncid, varid, values__packed )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! write packing attributes:
      status = NF90_Put_Att( self%ncid, varid, 'add_offset', add_offset )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( self%ncid, varid, 'scale_factor', scale_factor )
      IF_NF90_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( values__packed, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    else

      ! storage:
      allocate( values__out(size(values,1),size(values,2),size(values,3)), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! convert to output precission:
      if ( any((/present(fill_value),present(fill_value__out)/)) ) then
        ! check ...
        if ( .not. all((/present(fill_value),present(fill_value__out)/)) ) then
          write (csol,'("either none or both arguments `fill_value` and `` should be present")'); call csoErr
          TRACEBACK; status=1; return
        end if
        ! copy to output precission with change of fill values:
        where ( values == fill_value )
          values__out = fill_value__out
        elsewhere
          values__out = real(values,kind=rwp__out)
        end where
      else
         ! copy:
         values__out = real(values,kind=rwp__out)
      end if
    
      ! write:
      status = NF90_Put_Var( self%ncid, varid, values__out )
      IF_NF90_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( values__out, stat=status )
      IF_NOT_OK_RETURN(status=1)
    
    end if
    
    ! ok
    status = 0
    
  end subroutine NcFile_Put_Var_3d_r
  

  ! *
  
  
  !
  ! Gather values on root, write to variable.
  !

  subroutine NcFile_Put_Var1D_r( self, ivar, values, status, empty )
  
    use CSO_Comm, only : csoc
  
    ! --- in/out ---------------------------------
    
    class(T_NcFile), intent(inout)        ::  self
    integer, intent(in)                   ::  ivar
    real, intent(in)                      ::  values(:)
    integer, intent(out)                  ::  status
    
    logical, intent(in), optional         ::  empty

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/NcFile_Put_Var1D_r'
    
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
    if ( csoc%root ) then
      allocate( glb(nglb), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      allocate( glb(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
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
    
  end subroutine NcFile_Put_Var1D_r
  
  
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
  
