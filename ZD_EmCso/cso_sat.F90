!###############################################################################
!
! CSO_Sat - simulations of satellite retrievals
!
!
! HISTORY
!
!   2022-09, Arjo Segers
!     Support input and output of packed variables.
!
!   2023-01, Arjo Segers
!     Support files with "sample" coordinate instead of "pixel",
!     and without "corner" data (also implies no "area").
!
!   2023-01, Arjo Segers
!     Added option to define datatype of variable;
!     mainly used define dtype "i1" (integer(1)) or "c" (character)
!     that are used for ground observation data.
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
  
module CSO_Sat

  use NetCDF      , only : NF90_StrError, NF90_NOERR

  use CSO_Logging , only : csol, csoPr, csoErr
  use CSO_NcFile  , only : T_NcFile
  use CSO_Pixels  , only : T_PixelDatas
  use CSO_Pixels  , only : T_Track_0D, T_Track_1D
  use CSO_Mapping , only : T_Mapping
  use CSO_Exchange, only : T_Exchange
  use CSO_Swapping, only : T_Swapping
  use CSO_Domains , only : T_CSO_Selection1D

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public  ::  T_CSO_Sat_Data
  public  ::  T_CSO_Sat_State
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Sat'


  ! --- types --------------------------------
  
  !
  ! output selection
  !
  
  type T_OutputSelection
    ! selection key:
    character(len=32)                 ::  key
    ! selected variables:
    integer                           ::  n
    character(len=32), pointer        ::  vars(:)
    !
  contains
    procedure :: Init            => OutputSelection_Init
    procedure :: InitCopy        => OutputSelection_InitCopy
    procedure :: Done            => OutputSelection_Done
  end type T_OutputSelection
  
  ! 
  
  type T_OutputSelections
    ! selections:
    integer                                 ::  n
    type(T_OutputSelection), allocatable    ::  sels(:)  ! (n)
    !
  contains
    procedure :: Init            => OutputSelections_Init
    procedure :: InitCopy        => OutputSelections_InitCopy
    procedure :: Done            => OutputSelections_Done
    procedure :: GetSelection    => OutputSelections_GetSelection
  end type T_OutputSelections
  
  !
  ! satelite data:
  !  footprints
  !  timestamps
  !  kernels
  !
  type T_CSO_Sat_Data
    ! base for rcfile keys:
    character(len=32)                 ::  rcbase
    ! filename:
    character(len=1024)               ::  filename
    ! flags:
    logical                           ::  with_track
    !
    ! how to read?
    logical                           ::  read_on_root
    logical                           ::  read_by_me
    !
    !~ netcdf id's
    type(T_NcFile)                    ::  ncf
    !  deflation level for output:
    integer                           ::  deflate_level
    !  pack variables on output?
    logical                           ::  packed
    !
    !~ global pixel arrays,
    !  could be full track including pixels outside model domain
    !
    ! storage type:
    !   'pixel'
    !   'sample'
    character(len=32)                 ::  storetype
    !
    ! number of pixels in global doamin:
    integer                           ::  nglb
    ! number of corners in footprint:
    integer                           ::  ncorner
    ! number of layers in profiles:
    integer                           ::  nlayer
    integer                           ::  nretr
    !
    ! selection flags:
    logical, pointer                  ::  glb_select(:)  ! (nglb)
    ! track info (2D footprint)
    integer                           ::  ntx, nty         ! cross track, along track
    type(T_Track_0D)                  ::  glb_track_lon    ! (ntx,nty)
    type(T_Track_0D)                  ::  glb_track_lat    ! (ntx,nty)
    type(T_Track_1D)                  ::  glb_track_clons  ! (ncorner,ntx,nty)
    type(T_Track_1D)                  ::  glb_track_clats  ! (ncorner,ntx,nty)
    ! 
    !~ local pixels
    !  (overlap with local domain)
    !
    ! number of pixels loaded:
    integer                           ::  npix
    !
    ! index in global arrays:
    integer, allocatable              ::  iglb(:)   ! (npix)
    !
    ! local selection of footprints:
    real, pointer                     ::  loc_lons(:)     ! (npix)
    real, pointer                     ::  loc_lats(:)     ! (npix)
    real, pointer                     ::  loc_clons(:,:)  ! (ncorner,npix)
    real, pointer                     ::  loc_clats(:,:)  ! (ncorner,npix)
    !
    ! selection info:
    type(T_CSO_Selection1D)           ::  slc
    !
    ! mapping from grid cells to pixels:
    type(T_Mapping)                   ::  mapping
    ! put out the mapping arrays?
    logical                           ::  putout_mapping
    !
    ! exchange of pixels with other domains:
    type(T_Exchange), pointer         ::  exch(:)  ! (0:npes-1)
    ! flag:
    logical                           ::  exchange_is_setup
    !
    ! info on how to swap pixels to other domain decomposition:
    type(T_Swapping), pointer         ::  swp
    ! flag:
    logical                           ::  swap_is_setup
    !
    ! pixel data:
    type(T_PixelDatas)                ::  pd
    !
  contains
    procedure :: Init            => CSO_Sat_Data_Init
    procedure :: InitSwap        => CSO_Sat_Data_InitSwap
    procedure :: Done            => CSO_Sat_Data_Done
    procedure :: Get             => CSO_Sat_Data_Get
    procedure :: GetData         => CSO_Sat_Data_GetData
    procedure :: Read            => CSO_Sat_Data_Read
    procedure :: GetPixel        => CSO_Sat_Data_GetPixel
    procedure :: SetPixel        => CSO_Sat_Data_SetPixel
    procedure :: SetMapping      => CSO_Sat_Data_SetMapping
    procedure :: SetupExchange   => CSO_Sat_Data_SetupExchange
    procedure :: PutOut          => CSO_Sat_Data_PutOut
  end type T_CSO_Sat_Data
  
  !
  ! state simulations
  !
  type T_CSO_Sat_State
    ! annote:
    character(len=256)                ::  description
    ! copied from data:
    integer                           ::  npix
    ! pixel data:
    type(T_PixelDatas)                ::  pd
    ! number of user-defined variables:
    integer                           ::  nvar
    ! number of user-defined dimensions and their names:
    integer                           ::  nudim
    character(len=32), allocatable    ::  udimnames(:)   ! (ndim)
    ! ouptut selection:
    type(T_OutputSelections)          ::  outkeys
    !
  contains
    procedure :: Init            => CSO_Sat_State_Init
    procedure :: InitSwap        => CSO_Sat_State_InitSwap
    procedure :: Done            => CSO_Sat_State_Done
    procedure :: GetData         => CSO_Sat_State_GetData
    procedure :: GatherData      => CSO_Sat_State_GatherData
    procedure :: GetPixel        => CSO_Sat_State_GetPixel
    procedure :: SetPixel        => CSO_Sat_State_SetPixel
    procedure :: Get             => CSO_Sat_State_Get
    procedure :: GetDim          => CSO_Sat_State_GetDim
    procedure :: SetDim          => CSO_Sat_State_SetDim
    procedure :: EndDef          => CSO_Sat_State_EndDef
    procedure :: Exchange        => CSO_Sat_State_Exchange
    procedure :: ApplyFormulas   => CSO_Sat_State_ApplyFormulas
    procedure :: GetForcing      => CSO_Sat_State_GetForcing
    procedure :: PutOut          => CSO_Sat_State_PutOut
    procedure :: ReadForcing     => CSO_Sat_State_ReadForcing
  end type T_CSO_Sat_State


  
contains


  ! ====================================================================
  ! ===
  ! === Output selection
  ! ===
  ! ====================================================================


  subroutine OutputSelection_Init( self, rcf, rcbase, key, status )
  
    use CSO_Comm  , only : csoc
    use CSO_Rc    , only : T_CSO_RcFile
    use CSO_String, only : CSO_SplitString

    ! --- in/out ---------------------------------
    
    class(T_OutputSelection), intent(out)       ::  self
    type(T_CSO_RcFile), intent(in)              ::  rcF
    character(len=*), intent(in)                ::  rcbase
    character(len=*), intent(in)                ::  key
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/OutputSelection_Init'
    
    ! --- local ----------------------------------

    character(len=1024)       ::  line
    character(len=32)         ::  values(100)
    integer                   ::  i

    ! --- begin ----------------------------------
    
    ! store:
    self%key = trim(key)

    ! line with variable names:
    call rcF%Get( trim(rcbase)//'.outkey.'//trim(key)//'.vars', line, status )
    IF_NOT_OK_RETURN(status=1)
    ! split:
    call CSO_SplitString( line, self%n, values, status )
    IF_NOT_OK_RETURN(status=1)
    ! any?
    if ( self%n > 0 ) then
      ! storage:
      allocate( self%vars(self%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy:
      do i = 1, self%n
        self%vars(i) = trim(values(i))
      end do
    end if ! n > 0

    ! ok
    status = 0
    
  end subroutine OutputSelection_Init


  ! ***
  

  subroutine OutputSelection_InitCopy( self, outsel, status )
  
    use CSO_Comm  , only : csoc
    use CSO_Rc    , only : T_CSO_RcFile
    use CSO_String, only : CSO_SplitString

    ! --- in/out ---------------------------------
    
    class(T_OutputSelection), intent(out)       ::  self
    type(T_OutputSelection), intent(in)         ::  outsel
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/OutputSelection_InitCopy'
    
    ! --- local ----------------------------------

    integer                   ::  i

    ! --- begin ----------------------------------
    
    ! copy number:
    self%n = outsel%n
    ! any?
    if ( self%n > 0 ) then
      ! storage:
      allocate( self%vars(self%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy:
      do i = 1, self%n
        self%vars(i) = outsel%vars(i)
      end do
    end if ! n > 0

    ! ok
    status = 0
    
  end subroutine OutputSelection_InitCopy
  
  
  ! ***


  subroutine OutputSelection_Done( self, status )

    use CSO_Comm, only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_OutputSelection), intent(inout)   ::  self
    integer, intent(out)                      ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/OutputSelection_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! any?
    if ( self%n > 0 ) then
      ! storage:
      deallocate( self%vars, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! n > 0

    ! ok
    status = 0
    
  end subroutine OutputSelection_Done

  
  ! ***
  
  
  subroutine OutputSelections_Init( self, rcf, rcbase, status )
  
    use CSO_Comm  , only : csoc
    use CSO_Rc    , only : T_CSO_RcFile
    use CSO_String, only : CSO_SplitString

    ! --- in/out ---------------------------------
    
    class(T_OutputSelections), intent(out)      ::  self
    type(T_CSO_RcFile), intent(in)              ::  rcF
    character(len=*), intent(in)                ::  rcbase
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/OutputSelections_Init'
    
    ! --- local ----------------------------------

    character(len=1024)       ::  line
    character(len=32)         ::  values(100)
    integer                   ::  i

    ! --- begin ----------------------------------
    
    ! line with variable names, might not be defined ...
    call rcF%Get( trim(rcbase)//'.outkeys', line, status, default='' )
    IF_ERROR_RETURN(status=1)
    ! split:
    call CSO_SplitString( line, self%n, values, status )
    IF_NOT_OK_RETURN(status=1)
    ! any?
    if ( self%n > 0 ) then
      ! storage:
      allocate( self%sels(self%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! loop:
      do i = 1, self%n
        ! init selection:
        call self%sels(i)%Init( rcF, rcbase, values(i), status )
        IF_NOT_OK_RETURN(status=1)
      end do
    end if ! n > 0

    ! ok
    status = 0
    
  end subroutine OutputSelections_Init

  
  ! ***
  
  
  subroutine OutputSelections_InitCopy( self, out, status )
  
    use CSO_Comm  , only : csoc
    use CSO_Rc    , only : T_CSO_RcFile
    use CSO_String, only : CSO_SplitString

    ! --- in/out ---------------------------------
    
    class(T_OutputSelections), intent(out)      ::  self
    type(T_OutputSelections), intent(in)        ::  out
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/OutputSelections_InitCopy'
    
    ! --- local ----------------------------------

    integer                   ::  i

    ! --- begin ----------------------------------
    
    ! copy:
    self%n = out%n
    ! any?
    if ( self%n > 0 ) then
      ! storage:
      allocate( self%sels(self%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! loop:
      do i = 1, self%n
        ! init selection:
        call self%sels(i)%InitCopy( out%sels(i), status )
        IF_NOT_OK_RETURN(status=1)
      end do
    end if ! n > 0

    ! ok
    status = 0
    
  end subroutine OutputSelections_InitCopy
  
  
  ! ***


  subroutine OutputSelections_Done( self, status )

    use CSO_Comm, only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_OutputSelections), intent(inout)    ::  self
    integer, intent(out)                        ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/OutputSelections_Done'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! any?
    if ( self%n > 0 ) then
      ! loop:
      do i = 1, self%n
        ! clear:
        call self%sels(i)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do
      ! storage:
      deallocate( self%sels, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! n > 0

    ! ok
    status = 0
    
  end subroutine OutputSelections_Done
  
  
  ! ***

  !
  ! Return pointer to variable list for key.
  !

  subroutine OutputSelections_GetSelection( self, key, vars, status )

    use CSO_Comm, only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_OutputSelections), intent(in)       ::  self
    character(len=*), intent(in)                ::  key
    character(len=*), pointer                   ::  vars(:)
    integer, intent(out)                        ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/OutputSelections_GetSelection'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    
    ! --- begin ----------------------------------
    
    ! init result:
    nullify( vars ) 
    
    ! any?
    if ( self%n > 0 ) then
      ! loop:
      do i = 1, self%n
        ! match?
        if ( trim(self%sels(i)%key) == trim(key) ) then
          ! assign:
          vars => self%sels(i)%vars
          ! found:
          exit
        end if
      end do
      ! check:
      if ( .not. associated(vars) ) then
        write (csol,'("output selection key `",a,"` not defined:")') trim(key); call csoErr
        do i = 1, self%n
          write (csol,'(i6,"  ",a)') i, trim(self%sels(i)%key); call csoErr
        end do
        TRACEBACK; status=1; return
      end if
    end if ! n > 0

    ! ok
    status = 0
    
  end subroutine OutputSelections_GetSelection




  ! ====================================================================
  ! ===
  ! === Data
  ! ===
  ! ====================================================================


  subroutine CSO_Sat_Data_Init( self, rcF, rcbase, filename, status )
  
    use CSO_Rc              , only : T_CSO_RcFile

    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_Open
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(out)          ::  self
    type(T_CSO_RcFile), intent(in)              ::  rcF
    character(len=*), intent(in)                ::  rcbase
    character(len=*), intent(in)                ::  filename
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_Data_Init'
    
    ! --- local ----------------------------------
    
    logical                 ::  exist
    integer                 ::  dimid
    integer                 ::  pid
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (csol,'(a,": initialize satellite data ...")') rname; call csoPr

    ! store:
    self%rcbase   = trim(rcbase)
    self%filename = trim(filename)
    
    ! flags:
    call rcF%Get( trim(self%rcbase)//'.with_track'  , self%with_track, status )
    IF_NOT_OK_RETURN(status=1)
    call rcF%Get( trim(self%rcbase)//'.with_mapping', self%putout_mapping, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! on output: enable deflation?
    call rcF%Get( trim(self%rcbase)//'.output.deflate_level', self%deflate_level, status )
    IF_NOT_OK_RETURN(status=1)
    ! on output: pack floats to short int?
    call rcF%Get( trim(self%rcbase)//'.output.packed', self%packed, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read on root, scatter to domains?
    ! by default .true. now, flag used for testing ...
    call rcF%Get( trim(self%rcbase)//'.read_on_root', self%read_on_root, status, &
                    default=.true. )
    IF_ERROR_RETURN(status=1)
    
    ! read here?
    if ( self%read_on_root ) then
      ! only root:
      self%read_by_me = csoc%root
    else
      ! all:
      self%read_by_me = .true.
    end if
    
    ! read by this pe?
    if ( self%read_by_me ) then
    
      ! open file for reading:
      call self%ncf%Init( trim(self%filename), 'r', status )
      IF_NOT_OK_RETURN(status=1)

      ! pixel or station dimension:
      status = NF90_Inq_DimID( self%ncf%ncid, 'pixel', dimid )
      if ( status == NF90_NOERR ) then
        ! number of pixels
        status = NF90_Inquire_Dimension( self%ncf%ncid, dimid, len=self%nglb )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! set key ..
        self%storetype = 'pixel'
      else
        ! no pixels, probably a stations file ...
        status = NF90_Inq_DimID( self%ncf%ncid, 'sample', dimid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! number of stations
        status = NF90_Inquire_Dimension( self%ncf%ncid, dimid, len=self%nglb )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! set key ..
        self%storetype = 'sample'
      end if

      ! OPTIONAL: corner dimension:
      status = NF90_Inq_DimID( self%ncf%ncid, 'corner', dimid )
      if ( status == NF90_NOERR ) then
        ! number of corner values;
        status = NF90_Inquire_Dimension( self%ncf%ncid, dimid, len=self%ncorner )
        IF_NF90_NOT_OK_RETURN(status=1)
      else
        ! dummy:
        self%ncorner = -999
      end if

      ! layer dimension (used for kernel):
      status = NF90_Inq_DimID( self%ncf%ncid, 'layer', dimid )
      if ( status == NF90_NOERR ) then
        status = NF90_Inquire_Dimension( self%ncf%ncid, dimid, len=self%nlayer )
        IF_NF90_NOT_OK_RETURN(status=1)
      else
        self%nlayer = -999
      end if

      ! retrieval dimension:
      status = NF90_Inq_DimID( self%ncf%ncid, 'retr', dimid )
      if ( status == NF90_NOERR ) then
        status = NF90_Inquire_Dimension( self%ncf%ncid, dimid, len=self%nretr )
        IF_NF90_NOT_OK_RETURN(status=1)
      else
        self%nretr = -999
      end if
      
    end if  ! read by me
    ! need to broadcast?
    if ( self%read_on_root ) then
      ! broadcast from root:
      call csoc%BCast( csoc%root_id, self%nglb     , status )
      IF_NOT_OK_RETURN(status=1)
      call csoc%BCast( csoc%root_id, self%storetype, status )
      IF_NOT_OK_RETURN(status=1)
      call csoc%BCast( csoc%root_id, self%ncorner  , status )
      IF_NOT_OK_RETURN(status=1)
      call csoc%BCast( csoc%root_id, self%nlayer   , status )
      IF_NOT_OK_RETURN(status=1)
      call csoc%BCast( csoc%root_id, self%nretr    , status )
      IF_NOT_OK_RETURN(status=1)
    end if
      

    ! * read global footprints

    ! info ...
    write (csol,'(a,":   number of pixels in file: ",i0)') rname, self%nglb; call csoPr

    ! init pixel data:
    call self%pd%Init( status, read_on_root=self%read_on_root )
    IF_NOT_OK_RETURN(status=1)

    ! define pixel data, read from ncfile, store global arrays (all pixels) on each domain:
    call self%pd%NcInit( 'longitude', '', self%ncf%ncid, 'longitude', status, glb=.true. )
    IF_NOT_OK_RETURN(status=1)
    call self%pd%NcInit( 'latitude' , '', self%ncf%ncid, 'latitude' , status, glb=.true. )
    IF_NOT_OK_RETURN(status=1)
    ! corners defined?
    if ( self%ncorner > 0 ) then
      call self%pd%NcInit( 'longitude_bounds', 'corner', self%ncf%ncid, 'longitude_bounds', status, glb=.true. )
      IF_NOT_OK_RETURN(status=1)
      call self%pd%NcInit( 'latitude_bounds' , 'corner', self%ncf%ncid, 'latitude_bounds' , status, glb=.true. )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! storage for selection flags:
    allocate( self%glb_select(self%nglb), source=.false., stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! track defined? only used for convenience in output, read on root only:
    if ( self%with_track .and. csoc%root ) then

      ! track dimension:
      status = NF90_Inq_DimID( self%ncf%ncid, 'track_pixel', dimid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! number of pixels
      status = NF90_Inquire_Dimension( self%ncf%ncid, dimid, len=self%ntx )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! track dimension:
      status = NF90_Inq_DimID( self%ncf%ncid, 'track_scan', dimid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! number of pixels
      status = NF90_Inquire_Dimension( self%ncf%ncid, dimid, len=self%nty )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! info ...
      write (csol,'(a,":   shape of track: ",i0," x ",i0)') rname, self%ntx, self%nty; call csoPr

      ! read global arrays from file:
      call self%glb_track_lon%NcInit( self%ncf, 'track_longitude', status )
      IF_NOT_OK_RETURN(status=1)
      call self%glb_track_lat%NcInit( self%ncf, 'track_latitude', status )
      IF_NOT_OK_RETURN(status=1)
      call self%glb_track_clons%NcInit( self%ncf, 'track_longitude_bounds', status )
      IF_NOT_OK_RETURN(status=1)
      call self%glb_track_clats%NcInit( self%ncf, 'track_latitude_bounds', status )
      IF_NOT_OK_RETURN(status=1)

      ! read compressed index in track:
      call self%pd%NcInit( 'pixel', '', self%ncf%ncid, 'pixel', status, &
                              xtype='integer', glb=.true., only_here=.true. )
      IF_NOT_OK_RETURN(status=1)

    end if ! track

    ! * local pixels

    ! no pixels yet:
    self%npix     = 0
    
    ! no local selection:
    nullify( self%loc_lons )
    nullify( self%loc_lats )
    nullify( self%loc_clons )
    nullify( self%loc_clats )
    
    ! init mapping info:
    call self%mapping%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! exchange info needed?
    if ( csoc%npes > 1 ) then
      ! storage:
      allocate( self%exch(0:csoc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! loop:
      do pid = 0, csoc%npes-1
        ! init empty:
        call self%exch(pid)%Init( pid, status )
        IF_NOT_OK_RETURN(status=1)
      end do ! pid
      ! not setup yet:
      self%exchange_is_setup = .false.
    end if ! npes > 0
    
    ! no swapping yet:
    self%swap_is_setup = .false.

    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_Init
  
  
  ! ***


  subroutine CSO_Sat_Data_InitSwap( self, sdata, doms, doms_f, status )
  
    use CSO_Comm    , only : csoc
    use CSO_Domains , only : T_CSO_Domains
    use CSO_Swapping, only : T_Swapping

    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(out)          ::  self
    type(T_CSO_Sat_Data), intent(in)            ::  sdata
    type(T_CSO_Domains), intent(in)             ::  doms
    type(T_CSO_Domains), intent(in)             ::  doms_f
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_Data_InitSwap'
    
    ! --- local ----------------------------------
    
    type(T_Swapping)          ::  mswp
    integer                   ::  pid
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (csol,'(a,": swap satellite data ...")') rname; call csoPr

    ! copy settings info:
    self%rcbase   = sdata%rcbase
    self%filename = sdata%filename
    
    ! copy flags:
    self%with_track = sdata%with_track
    
    ! copy dimensions:
    self%nglb    = sdata%nglb
    self%ncorner = sdata%ncorner
    self%nlayer  = sdata%nlayer
    self%nretr   = sdata%nretr
    
    ! copy track info (needed on root only):
    if ( csoc%root .and. (self%nglb > 0) ) then
      ! copy dims:
      self%ntx = sdata%ntx
      self%nty = sdata%nty
      ! copy tracks:
      call self%glb_track_lon%InitCopy( sdata%glb_track_lon, status )
      IF_NOT_OK_RETURN(status=1)
      call self%glb_track_lat%InitCopy( sdata%glb_track_lat, status )
      IF_NOT_OK_RETURN(status=1)
      call self%glb_track_clons%InitCopy( sdata%glb_track_clons, status )
      IF_NOT_OK_RETURN(status=1)
      call self%glb_track_clats%InitCopy( sdata%glb_track_clats, status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! copy flags:
    self%putout_mapping = sdata%putout_mapping
    
    ! ~ swap pixel array
    
    ! any pixels on some domain?
    if ( sdata%slc%nsel_tot > 0 ) then
    
      ! reset flag:
      self%swap_is_setup = .true.
      ! storage:
      allocate( self%swp, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! info on swapping from current decomposition 'doms' to decomposition 'doms_f';
      ! use the mapping info to know with which of the new domains a pixel overplaps ;
      ! define for swapping of (npix) arrays:
      call self%swp%Init( 'pix', sdata%mapping%map_n, sdata%mapping%map_ii, sdata%mapping%map_jj, &
                              doms, doms_f, status )
      IF_NOT_OK_RETURN(status=1)
      ! idem for (nmap) arrays:
      call mswp%Init( 'map', sdata%mapping%map_n, sdata%mapping%map_ii, sdata%mapping%map_jj, &
                              doms, doms_f, status )
      IF_NOT_OK_RETURN(status=1)

      ! new number of local pixels:
      self%npix = self%swp%nrecv

      ! storage per local pixel for global index:
      allocate( self%iglb(max(1,self%npix)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! swap global indices:
      call self%swp%Swap( sdata%iglb, self%iglb, status )
      IF_NOT_OK_RETURN(status=1)
      
    else
      ! no pixels:
      self%npix = 0
    end if

    ! store selection info, collect all info on root for:
    ! - scatter after  reading from root (optional)
    ! - gather  before writing from root (always needed):
    call self%slc%Init( self%nglb, self%npix, self%iglb, .true., status )
    IF_NOT_OK_RETURN(status=1)
    
    ! any pixels on some domain?
    if ( sdata%slc%nsel_tot > 0 ) then

      ! swap pixel data to target decomposition:
      call self%pd%InitSwap( sdata%pd, self%swp, status )
      IF_NOT_OK_RETURN(status=1)

      ! swap mapping info, requires swapping of both pixel and mapping arrays:
      call self%mapping%InitSwap( sdata%mapping, self%swp, mswp, status )
      IF_NOT_OK_RETURN(status=1)

      ! clear:
      call mswp%Done( status )
      IF_NOT_OK_RETURN(status=1)

    end if ! any pixels on some domain
    
    ! exchange info needed?
    if ( csoc%npes > 1 ) then
      ! storage:
      allocate( self%exch(0:csoc%npes-1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! loop:
      do pid = 0, csoc%npes-1
        ! init empty:
        call self%exch(pid)%Init( pid, status )
        IF_NOT_OK_RETURN(status=1)
      end do ! pid
    end if ! npes > 0
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_InitSwap


  ! ***


  subroutine CSO_Sat_Data_Done( self, status )

    use CSO_Comm, only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(inout)    ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_Data_Done'
    
    ! --- local ----------------------------------
    
    integer     ::  pid
    
    ! --- begin ----------------------------------
    
    ! any pixels in global domain?
    if ( self%nglb > 0 ) then
    
      ! storage for selection flags:
      deallocate( self%glb_select, stat=status )
      IF_NOT_OK_RETURN(status=1)
    
      ! track defined?
      if ( self%with_track .and. csoc%root ) then
        call self%glb_track_lon%Done( status )
        IF_NOT_OK_RETURN(status=1)
        call self%glb_track_lat%Done( status )
        IF_NOT_OK_RETURN(status=1)
        call self%glb_track_clons%Done( status )
        IF_NOT_OK_RETURN(status=1)
        call self%glb_track_clats%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end if ! track

      ! global index:
      deallocate( self%iglb, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! output collection:
      call self%slc%Done( status )
      IF_NOT_OK_RETURN(status=1)
      
      ! local pixels?
      if ( self%npix > 0 ) then
        deallocate( self%loc_lons, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( self%loc_lats, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! corners defined?
        if ( self%ncorner > 0 ) then
          deallocate( self%loc_clons, stat=status )
          IF_NOT_OK_RETURN(status=1)
          deallocate( self%loc_clats, stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if
      end if
      
    end if  ! nglb > 0
    
    ! exchange info needed?
    if ( csoc%npes > 1 ) then
      ! loop:
      do pid = 0, csoc%npes-1
        ! init empty:
        call self%exch(pid)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do ! pid
      ! clear:
      deallocate( self%exch, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! npes > 0
    
    ! swap is used?
    if ( self%swap_is_setup ) then
      ! clear:
      call self%swp%Done( status )
      IF_NOT_OK_RETURN(status=1)
      ! clear:
      deallocate( self%swp, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! reset flag:
      self%swap_is_setup = .false.
    end if ! swap is used

    ! done with mapping info:
    call self%mapping%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! done with pixel data:
    call self%pd%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_Done
  
  
  ! ***


  subroutine CSO_Sat_Data_Read( self, rcF, rcbase, status )
  
    use NetCDF     , only : NF90_Close

    use CSO_Comm   , only : csoc
    use CSO_Domains, only : T_CSO_Selection1D
    use CSO_Rc     , only : T_CSO_RcFile 
    use CSO_String , only : CSO_ReadFromLine 
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(inout)        ::  self
    type(T_CSO_RcFile), intent(in)              ::  rcF
    character(len=*), intent(in)                ::  rcbase
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_Data_Read'
    
    ! --- local ----------------------------------
    
    integer, allocatable    ::  iglb_local(:)  ! (nglb)
    integer                 ::  iglb
    integer                 ::  dimid
    integer                 ::  ipix
    integer                 ::  id
    character(len=1024)     ::  line
    character(len=32)       ::  vname
    character(len=256)      ::  dnames
    character(len=256)      ::  source
    character(len=32)       ::  xtype
    real, pointer           ::  glb_lons(:)
    real, pointer           ::  glb_lats(:)
    real, pointer           ::  glb_clons(:,:)
    real, pointer           ::  glb_clats(:,:)
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (csol,'(a,": read local pixels ...")') rname; call csoPr

    ! * read variables

    ! temporary mapping at maximum size, init with no-data:
    allocate( iglb_local(self%nglb), source=-999, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over all pixels:
    do iglb = 1, self%nglb
      ! selected?
      if ( self%glb_select(iglb) ) then
        ! increase counter for local pixels:
        self%npix = self%npix + 1
        ! store global index corresponndig to local pixel:
        iglb_local(self%npix) = iglb
      end if
    end do ! glb pixels
    
    ! info ..
    write (csol,'(a,": selected ",i0," local pixels out of ",i0," globally")') rname, self%npix, self%nglb; call csoPr

    ! any local pixels?
    if ( self%npix > 0 ) then

      ! storage for mapping:
      allocate( self%iglb(self%npix), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy from array with global length:
      self%iglb = iglb_local(1:self%npix)
      
    else

      ! dummy ...
      allocate( self%iglb(1), source=-999, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if
    
    ! clear:
    deallocate( iglb_local, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! store selection info, collect all info on root for scatter and gather:
    call self%slc%Init( self%nglb, self%npix, self%iglb, .true., status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (csol,'(a,": define state variables ...")') rname; call csoPr

    ! line with variable names:
    call rcF%Get( trim(rcbase)//'.dvars', line, status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over elements:
    do
      ! empty?
      if ( len_trim(line) == 0 ) exit
      ! next part:
      call CSO_ReadFromLine( line, vname, status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)

      ! dimension names:
      call rcF%Get( trim(rcbase)//'.dvar.'//trim(vname)//'.dims', dnames, status )
      IF_NOT_OK_RETURN(status=1)
      ! source variable:
      call rcF%Get( trim(rcbase)//'.dvar.'//trim(vname)//'.source', source, status )
      IF_NOT_OK_RETURN(status=1)
      ! source variable:
      call rcF%Get( trim(rcbase)//'.dvar.'//trim(vname)//'.dtype', xtype, status, default='real' )
      IF_ERROR_RETURN(status=1)
      ! info ...
      write (csol,'(a,":  variable `",a,"` ...")') rname, trim(vname); call csoPr
      write (csol,'(a,":    dimensions : ",a)') rname, trim(dnames); call csoPr
      write (csol,'(a,":    source     : ",a)') rname, trim(source); call csoPr
      write (csol,'(a,":    dtype      : ",a)') rname, trim(xtype); call csoPr

      ! define variable, and read (local subset) from nc file:
      call self%pd%NcInit( trim(vname), trim(dnames), self%ncf%ncid, trim(source), status, &
                             slc=self%slc, xtype=xtype )
      IF_NOT_OK_RETURN(status=1)

    end do ! dvars
    
    ! file opened?
    if ( self%read_by_me ) then
      ! close:
      call self%ncf%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end if ! read by me

    ! *
    
    ! local selection of footprints:
    if ( self%npix > 0 ) then
      ! pointers:
      call self%pd%GetData( status, name='longitude', data0=glb_lons )
      IF_NOT_OK_RETURN(status=1)
      call self%pd%GetData( status, name='latitude', data0=glb_lats )
      IF_NOT_OK_RETURN(status=1)
      ! corners defined?
      if ( self%ncorner > 0 ) then
        call self%pd%GetData( status, name='longitude_bounds', data1=glb_clons )
        IF_NOT_OK_RETURN(status=1)
        call self%pd%GetData( status, name='latitude_bounds', data1=glb_clats )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! storage:
      allocate( self%loc_lons(self%npix), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( self%loc_lats(self%npix), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! corners defined?
      if ( self%ncorner > 0 ) then
        allocate( self%loc_clons(self%ncorner,self%npix), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%loc_clats(self%ncorner,self%npix), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if  ! ncorner > 0
      ! loop:
      do ipix = 1, self%npix
        ! global index:
        iglb = self%iglb(ipix)
        ! copy:
        self%loc_lons(ipix)    = glb_lons(iglb)
        self%loc_lats(ipix)    = glb_lats(iglb)
        ! corners defined?
        if ( self%ncorner > 0 ) then
          self%loc_clons(:,ipix) = glb_clons(:,iglb)
          self%loc_clats(:,ipix) = glb_clats(:,iglb)
        end if
      end do ! ipix
    end if

    ! specials ...
    select case ( trim(self%storetype) )
      !~ satellite pixels with footprint:
      case ( 'pixel' )
        ! storage for pixel area, filled with zeros to ensure that 
        ! after gathering the sum is correct:
        call self%pd%Def( 'area', '', id, status, source=0.0 )
        IF_NOT_OK_RETURN(status=1)
        ! attributes:
        call self%pd%SetAttr( id, 'long_name', 'pixel area', status )
        IF_NOT_OK_RETURN(status=1)
        call self%pd%SetAttr( id, 'units', 'm2', status )
        IF_NOT_OK_RETURN(status=1)
      !~ samples
      case ( 'sample' )
        ! nothing to be done
      !~
      case default
        write (csol,'("unsupported storetype `",a,"`")') trim(self%storetype); call csoErr
        TRACEBACK; status=1; return
    end select
    
    ! end definition phase, allocate arrays (if not done yet by NcInit),
    ! this is probably only the just defined 'area':
    call self%pd%EndDef( self%npix, status )
    IF_NOT_OK_RETURN(status=1)

    ! *

    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_Read
  
  
  ! ***


  !
  ! Get info on data set:
  !   nglb       : number of pixels in global domain
  !   glb_lon    : pointer to (nglb) array with pixel center
  !   glb_lat    : pointer to (nglb) array with pixel center
  !
  !   nout       : total number of selected pixels out of global arrays (put out)
  !

  subroutine CSO_Sat_Data_Get( self, status, &
                                 nglb, glb_lon, glb_lat, glb_clons, glb_clats, glb_select, &
                                 npix, npix_all, nlayer, nretr, &
                                 lons, lats, clons, clats, &
                                 nw, iw0, ii, jj, ww, areas, &
                                 storetype, nout )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(in)           ::  self
    integer, intent(out)                        ::  status

    integer, intent(out), optional              ::  nglb
    real, pointer, optional                     ::  glb_lon(:)  ! (nglb)
    real, pointer, optional                     ::  glb_lat(:)  ! (nglb)
    real, pointer, optional                     ::  glb_clons(:,:)  ! (ncorner,nglb)
    real, pointer, optional                     ::  glb_clats(:,:)  ! (ncorner,nglb)
    logical, pointer, optional                  ::  glb_select(:)  ! (nglb)

    integer, intent(out), optional              ::  npix
    integer, intent(out), optional              ::  npix_all
    integer, intent(out), optional              ::  nlayer
    integer, intent(out), optional              ::  nretr

    real, pointer, optional                     ::  lons(:)     ! (npix)
    real, pointer, optional                     ::  lats(:)     ! (npix)
    real, pointer, optional                     ::  clons(:,:)  ! (ncorner,npix)
    real, pointer, optional                     ::  clats(:,:)  ! (ncorner,npix)

    integer, pointer, optional                  ::  nw(:)      ! (npix)
    integer, pointer, optional                  ::  iw0(:)     ! (npix)
    integer, pointer, optional                  ::  ii(:)      ! (nmap)
    integer, pointer, optional                  ::  jj(:)      ! (nmap)
    real, pointer, optional                     ::  ww(:)      ! (nmap)
    real, pointer, optional                     ::  areas(:)   ! (npix)
    
    character(len=*), intent(out), optional     ::  storetype
    integer, intent(out), optional              ::  nout

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_Data_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! return global size?
    if ( present(nglb) ) nglb = self%nglb
    
    ! centers?
    if ( present(glb_lon) ) then
      call self%pd%GetData( status, name='longitude', data0=glb_lon )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( present(glb_lat) ) then
      call self%pd%GetData( status, name='latitude', data0=glb_lat )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! corners?
    if ( present(glb_clons) ) then
      call self%pd%GetData( status, name='longitude_bounds', data1=glb_clons )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( present(glb_clats) ) then
      call self%pd%GetData( status, name='latitude_bounds', data1=glb_clats )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! selection flags:
    if ( present(glb_select) ) glb_select => self%glb_select
    
    ! sizes:
    if ( present(npix  ) ) npix   = self%npix
    
    ! total over all domains:
    if ( present(npix_all) ) then
      ! local number:
      npix_all = self%npix
      ! sum over domains:
      call csoc%AllReduce( 'sum', npix_all, status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! number of apriori layers:
    if ( present(nlayer) ) then
      ! check if defined ...
      if ( self%nlayer < 0 ) then
        write (csol,'("number of apriori layers not defined (not present in data file?)")'); call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      nlayer = self%nlayer
    end if
    
    ! number of retrieval layers:
    if ( present(nretr) ) then
      ! check if defined ...
      if ( self%nretr < 0 ) then
        write (csol,'("number of retrieval layes not defined (not present in data file?)")'); call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      nretr = self%nretr
    end if
    
    ! local footprints:
    if ( present(lons ) ) lons  => self%loc_lons
    if ( present(lats ) ) lats  => self%loc_lats
    
    ! corner locations:
    if ( present(clons) ) then
      ! check:
      if ( self%ncorner <=0 ) then
        write (csol,'("could not pointer to `clons`, no corners defined")'); call csoPr
        TRACEBACK; status=1; return
      end if
      ! set pointer:
      clons => self%loc_clons
    end if
    if ( present(clats) ) then
      ! check:
      if ( self%ncorner <=0 ) then
        write (csol,'("could not pointer to `clats`, no corners defined")'); call csoPr
        TRACEBACK; status=1; return
      end if
      ! set pointer:
      clats => self%loc_clats
    end if
    
    ! mapping arrays:
    if ( present(nw ) ) nw  => self%mapping%map_n
    if ( present(iw0) ) iw0 => self%mapping%map_i0
    if ( present(ii ) ) ii  => self%mapping%map_ii
    if ( present(jj ) ) jj  => self%mapping%map_jj
    if ( present(ww ) ) ww  => self%mapping%map_ww
    ! pixel areas:
    if ( present(areas) ) then
      ! get pointer:
      call self%pd%GetData( status, name='area', data0=areas )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! copy:
    if ( present(storetype) ) storetype = self%storetype
    
    ! number of pixels put out (unique pixels)
    if ( present(nout) ) nout = self%slc%nout

    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_Get
  
  
  ! ***


  !
  ! Get info on pixel data specified by either "id" number or "name"
  !

  subroutine CSO_Sat_Data_GetData( self, status, id, name,&
                                     data0, data1, data2, units )
  

    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(in)           ::  self
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional             ::  id
    character(len=*), intent(in), optional    ::  name

    real, pointer, optional                   ::  data0(:)
    real, pointer, optional                   ::  data1(:,:)
    real, pointer, optional                   ::  data2(:,:,:)
    character(len=*), intent(out), optional   ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_Data_GetData'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! get units:
    call self%pd%GetData( status, id=id, name=name, &
                            data0=data0, data1=data1, data2=data2, &
                            units=units )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_GetData


  ! ***


  subroutine CSO_Sat_Data_GetPixel( self, ipix, status, &
                                            glbid, lon, lat, clons, clats, &
                                            name, value, profile, &
                                            nw, ii, jj, ww, pix_area )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(in)             ::  self
    integer, intent(in)                           ::  ipix
    integer, intent(out)                          ::  status
    
    integer, intent(out), optional                ::  glbid
    real, pointer, optional                       ::  lon
    real, pointer, optional                       ::  lat
    real, pointer, optional                       ::  clons(:)  ! (ncorner)
    real, pointer, optional                       ::  clats(:)  ! (ncorner)
    character(len=*), intent(in), optional        ::  name
    real, pointer, optional                       ::  value
    real, pointer, optional                       ::  profile(:)
    integer, intent(out), optional                ::  nw
    integer, pointer, optional                    ::  ii(:)   ! (nmap)
    integer, pointer, optional                    ::  jj(:)   ! (nmap)
    real, pointer, optional                       ::  ww(:)   ! (nmap)
    real, intent(out), optional                   ::  pix_area
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_Data_GetPixel'
    
    ! --- local ----------------------------------
    
    integer           ::  iglb
    real, pointer     ::  data0(:)      ! (npix)
    real, pointer     ::  data1(:,:)    ! (nv,npix)
    real, pointer     ::  data2(:,:,:)  ! (mv,nv,npix)
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ipix < 1) .or. (ipix > self%npix) ) then
      write (csol,'("pixel ",i0," out of range 1:",i0)') ipix, self%npix; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! global index:
    iglb = self%iglb(ipix)
    
    ! id?
    if ( present(glbid) ) glbid = iglb
    
    ! location; extract using global index!
    if ( present(lon) ) then
      !lon => self%glb_lon%data(iglb)
      ! get pointer:
      call self%pd%GetData( status, name='longitude', data0=data0 )
      IF_NOT_OK_RETURN(status=1)
      ! assign:
      lon => data0(iglb)
    end if
    if ( present(lat) ) then
      !lat => self%glb_lat%data(iglb)
      ! get pointer:
      call self%pd%GetData( status, name='latitude', data0=data0 )
      IF_NOT_OK_RETURN(status=1)
      ! assign:
      lat => data0(iglb)
    end if
    ! footprint:
    if ( present(clons) ) then
      !clons => self%glb_clons%data(:,iglb)
      ! get pointer:
      call self%pd%GetData( status, name='longitude_bounds', data1=data1 )
      IF_NOT_OK_RETURN(status=1)
      ! assign:
      clons => data1(:,iglb)
    end if
    if ( present(clats) ) then
      !clats => self%glb_clats%data(:,iglb)
      ! get pointer:
      call self%pd%GetData( status, name='latitude_bounds', data1=data1 )
      IF_NOT_OK_RETURN(status=1)
      ! assign:
      clats => data1(:,iglb)
    end if
    
    ! return value?
    if ( present(value) ) then
      ! get pointer:
      call self%pd%GetData( status, name=name, data0=data0 )
      IF_NOT_OK_RETURN(status=1)
      ! assign:
      value => data0(ipix)
    end if
    
    ! return profile?
    if ( present(profile) ) then
      ! get pointer:
      call self%pd%GetData( status, name=name, data1=data1 )
      IF_NOT_OK_RETURN(status=1)
      ! assign:
      profile => data1(:,ipix)
    end if
    
    ! get mapping info:
    call self%mapping%GetPixel( ipix, status, nw=nw, ii=ii, jj=jj, ww=ww )
    IF_NOT_OK_RETURN(status=1)
    
    ! pixel area:
    if ( present(pix_area) ) then
      !pix_area = self%area%data(ipix)
      ! get pointer:
      call self%pd%GetData( status, name='area', data0=data0 )
      IF_NOT_OK_RETURN(status=1)
      ! assign:
      pix_area = data0(ipix)
    end if
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_GetPixel


  ! ***


  subroutine CSO_Sat_Data_SetPixel( self, ipix, status, &
                                       pix_area, ii, jj, ww )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(inout)          ::  self
    integer, intent(in)                           ::  ipix
    integer, intent(out)                          ::  status
    
    real, intent(in), optional                    ::  pix_area
    integer, intent(in), optional                 ::  ii(:)   ! (nmap)
    integer, intent(in), optional                 ::  jj(:)   ! (nmap)
    real, intent(in), optional                    ::  ww(:)   ! (nmap)
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_Data_SetPixel'
    
    ! --- local ----------------------------------

    integer             ::  id
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ipix < 1) .or. (ipix > self%npix) ) then
      write (csol,'("pixel ",i0," out of range 1:",i0)') ipix, self%npix; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! mapping elements?
    if ( any((/present(pix_area),present(ii),present(jj),present(ww)/)) ) then
      ! check ..
      if ( .not. all((/present(pix_area),present(ii),present(jj),present(ww)/)) ) then
        write (csol,'("either none or all of area/ii/jj/ww arguments should be present")'); call csoErr
        TRACEBACK; status=1; return
      end if
    
      ! fill area:
      call self%pd%InqID( 'area', id, status, quiet=.true. )
      if ( status == 0 ) then
        call self%pd%SetPixel( id, ipix, pix_area, status )
        IF_NOT_OK_RETURN(status=1)
      else if ( status > 0 ) then
        TRACEBACK; status=1; return
      end if
      
      ! set mapping info:
      call self%mapping%SetPixel( ipix, self%npix, self%iglb(ipix), ii,  jj, ww, status )
      IF_NOT_OK_RETURN(status=1)
    
    end if  ! mapping elements
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_SetPixel


  ! ***

  
  !
  ! Store mapping data:
  ! - area(npix)  : pixel area [m2]
  ! - nw(npix)    : number of mapping elements per pixel
  ! - ii(nmap), jj(nmap), ww(nmap)  :  source indices and weights, size nmap=sum(nw)
  !

  subroutine CSO_Sat_Data_SetMapping( self, area, nw, ii, jj, ww, status )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(inout)          ::  self
    real, intent(in)                              ::  area(:) ! (npix)
    integer, intent(in)                           ::  nw(:)   ! (npix)
    integer, intent(in)                           ::  ii(:)   ! (nmap)
    integer, intent(in)                           ::  jj(:)   ! (nmap)
    real, intent(in)                              ::  ww(:)   ! (nmap)
    integer, intent(out)                          ::  status
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_Data_SetMapping'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! specials ...
    select case ( trim(self%storetype) )
      !~ satellite pixels with footprint:
      case ( 'pixel' )
        ! store area:
        call self%pd%SetData( status, name='area', data0=area )
        IF_NOT_OK_RETURN(status=1)
      !~ samples
      case ( 'sample' )
        ! nothing to be done
      !~
      case default
        write (csol,'("unsupported storetype `",a,"`")') trim(self%storetype); call csoErr
        TRACEBACK; status=1; return
    end select
    
    ! store mapping info:
    call self%mapping%SetData( self%npix, self%iglb, nw, ii,  jj, ww, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_SetMapping


  ! ***


  subroutine CSO_Sat_Data_PutOut( self, filename, status )

    use NetCDF , only : NF90_Def_Dim
    use NetCDF , only : NF90_EndDef
  
    use CSO_Comm  , only : csoc
    use CSO_NcFile, only : T_NcFile
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(inout)          ::  self
    character(len=*), intent(in)                  ::  filename
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_Data_PutOut'
    
    ! --- local ----------------------------------
    
    type(T_NcFile)                ::  ncf
    integer                       ::  dimid_tx, dimid_ty
    integer                       ::  dimid_pixel
    integer                       ::  dimid_corner
    integer                       ::  dimid_layer
    integer                       ::  dimid_layeri
    integer                       ::  dimid_retr
    character(len=32), pointer    ::  vars(:)
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (csol,'(a,": put out sat data to: ",a)') rname, trim(filename); call csoPr
    
    ! no selection:
    nullify( vars )

    ! any pixels to be put out?
    if ( self%slc%nout > 0 ) then
    
      ! info ..
      write (csol,'(a,":   put out ",i0," pixels ...")') rname, self%slc%nout; call csoPr
    
      ! ~ collect mappings

      ! put out mapping arrays?
      if ( self%putout_mapping ) then
        ! collect mappings from grid cells to pixels:
        ! - fill self%nmap_all, 
        ! - collect arrays, store in map_*_out arrays
        call self%mapping%Collect( self%npix, self%slc%iout_tot, self%slc%nout, self%slc%iout_glb, status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      !! testing ...
      !write (csol,'(a,": nmap = ",i0," ; nmap_all = ",i0)') rname, self%nmap, self%nmap_all; call csoPr
      
      ! ~ create files
    
      ! collect on root, on some systems this is much faster than parallel write ...
      ! (in future, switch to support both?)
      if ( csoc%root ) then

        ! create file:
        call ncf%Init( filename, 'w', status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ define

        ! define dimensions:
        status = NF90_Def_Dim( ncf%ncid, trim(self%storetype), self%slc%nout, dimid_pixel )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! track defined?
        if ( self%with_track ) then
          status = NF90_Def_Dim( ncf%ncid, 'track_pixel', self%ntx, dimid_tx )
          IF_NF90_NOT_OK_RETURN(status=1)
          status = NF90_Def_Dim( ncf%ncid, 'track_scan', self%nty, dimid_ty )
          IF_NF90_NOT_OK_RETURN(status=1)
        end if ! track

        ! define netcdf variables:
        call self%pd%NcDef( ncf%ncid, dimid_pixel, vars, status, &
                               deflate_level=self%deflate_level, &
                               packed=self%packed )
        IF_NOT_OK_RETURN(status=1)
        
        ! corners defined?
        if ( self%ncorner > 0 ) then
          ! adhoc, not all into pd yet ...
          call self%pd%ncdims%GetDim( 'corner', status, dimid=dimid_corner )
          IF_NOT_OK_RETURN(status=1)
        end if

        ! track defined?
        if ( self%with_track ) then
          !~ global track arrays:
          call self%glb_track_lon%NcDef( ncf, 'track_longitude', &
                                           (/dimid_tx,dimid_ty/), status, &
                                           deflate_level=self%deflate_level, &
                                           packed=self%packed )
          IF_NOT_OK_RETURN(status=1)
          call self%glb_track_lat%NcDef( ncf, 'track_latitude', &
                                           (/dimid_tx,dimid_ty/), status, &
                                           deflate_level=self%deflate_level, &
                                           packed=self%packed )
          IF_NOT_OK_RETURN(status=1)
          call self%glb_track_clons%NcDef( ncf, 'track_longitude_bounds', &
                                           (/dimid_corner,dimid_tx,dimid_ty/), status, &
                                           deflate_level=self%deflate_level, &
                                           packed=self%packed )
          IF_NOT_OK_RETURN(status=1)
          call self%glb_track_clats%NcDef( ncf, 'track_latitude_bounds' , &
                                           (/dimid_corner,dimid_tx,dimid_ty/), status, &
                                           deflate_level=self%deflate_level, packed=self%packed )
          IF_NOT_OK_RETURN(status=1)
        end if ! track
        ! put out mapping arrays?
        if ( self%putout_mapping ) then
          ! define dims and variables:
          call self%mapping%NcDef( ncf%ncid, dimid_pixel, status, &
                                    deflate_level=self%deflate_level, &
                                    packed=self%packed )
          IF_NOT_OK_RETURN(status=1)
        end if ! mapping

        ! end defintion mode:
        status = NF90_EndDef( ncf%ncid )
        IF_NF90_NOT_OK_RETURN(status=1)
        
        ! ~ write

        ! put selection of "glb" arrays:
        call self%pd%NcPutGlbSelect( ncf, self%slc%iout_glb, status )
        IF_NOT_OK_RETURN(status=1)

        ! track defined?
        if ( self%with_track ) then
          ! write track arrays available as global arrays:
          call self%glb_track_lon%NcPutGlb( ncf, status )
          IF_NOT_OK_RETURN(status=1)
          call self%glb_track_lat%NcPutGlb( ncf, status )
          IF_NOT_OK_RETURN(status=1)
          call self%glb_track_clons%NcPutGlb( ncf, status )
          IF_NOT_OK_RETURN(status=1)
          call self%glb_track_clats%NcPutGlb( ncf, status )
          IF_NOT_OK_RETURN(status=1)
        end if  ! track
        
        ! put out mapping arrays?
        if ( self%putout_mapping ) then
          ! write variables:
          call self%mapping%NcPut( ncf%ncid, status )
          IF_NOT_OK_RETURN(status=1)
        end if ! mapping
        
      end if ! root
      
      ! collect distributed arrays on root and write from there,
      ! values found on multiple domains are the same so no need to add them:
      call self%pd%NcPutGather( ncf, self%slc%iout_tot, .false., vars, status )
      IF_NOT_OK_RETURN(status=1)

      ! written on root...
      if ( csoc%root ) then
        ! close:
        call ncf%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end if  ! root
    
    else
    
      ! info ..
      write (csol,'(a,":   no pixels to be put out ...")') rname; call csoPr
      
    end if  ! nout > 0

    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_PutOut


  ! ***
  
  
  !
  ! Setup information on exchange of simulations between domains,
  ! for pixels with contributions from more than one domain
  ! (footprint overlaps only partly with each domain)
  !
  ! Available after "Read":
  !   self%iglb_all(1:npix_all)
  !     Global pixel index for each of the local pixel contributions.
  !


  subroutine CSO_Sat_Data_SetupExchange( self, status )

    use CSO_Comm        , only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_Data), intent(inout)          ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_Data_SetupExchange'
    
    ! --- local ----------------------------------
    
    integer                          ::  ipix
    integer                          ::  np              ! 1=pid, 2=ipix
    integer, allocatable             ::  idi_loc(:,:)    ! (np,npix)
    integer, allocatable             ::  idi_all(:,:)    ! (np,npix_all) on root
    integer                          ::  maxpid
    integer, allocatable             ::  qn(:)         ! (nglb)
    integer, allocatable             ::  qidi(:,:,:)   ! (np,maxpid,nglb)
    integer                          ::  ipix_all
    integer                          ::  iglb
    integer                          ::  pid1, pid2
    integer                          ::  pid
    integer                          ::  nex_max
    integer                          ::  nex
    integer, allocatable             ::  exchi(:,:)   ! (nex_max,2)
    integer                          ::  i1, i2
    integer                          ::  k
    integer                          ::  itag1, itag2

    ! --- begin ----------------------------------
    
    !! info ...
    !write (csol,'(a,": determine pixels overlapping with multiple domains ...")') rname; call csoPr
    
    ! only need for multiple pe's, and if not done yet ...
    if ( (csoc%npes > 1) .and. (.not. self%exchange_is_setup) ) then

      ! storage parameters:
      !    pid, ipix
      np = 2

      ! any local pixels?
      if ( self%npix > 0 ) then
        ! storage for processor id:
        allocate( idi_loc(np,self%npix), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill local processor:
        idi_loc(1,:) = csoc%id
        ! fill local pixel indices:
        do k = 1, self%npix
          idi_loc(2,k) = k
        end do
        !! storage for local pixel id:
        !allocate( ipix(self%npix), source=csoc%id, stat=status )
        !IF_NOT_OK_RETURN(status=1)
      else
        ! dummy ...
        allocate( idi_loc(np,1), source=-999, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! storage for processor id's per pixel contribution, needed on root only:
      if ( csoc%root ) then
        allocate( idi_all(np,self%slc%nsel_tot), source=-999, stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( idi_all(np,1), source=-999, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! gather on root:
      call csoc%GatherV( idi_loc, idi_all, status, nloc=self%npix )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill exchange info on root:
      if ( csoc%root ) then
        ! assume max 4 domains contribute to single pixel:
        maxpid = 4
        ! storage for number of pids that contribute to pixel:
        allocate( qn(self%nglb), source=0, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! storage for pids that contribute to pixel:
        allocate( qidi(np,maxpid,self%nglb), source=-999, stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over pixel contributions:
        do ipix_all = 1, self%slc%nsel_tot
          ! global index:
          iglb = self%slc%iglb_tot(ipix_all)
          ! increase counter:
          qn(iglb) = qn(iglb) + 1
          ! check ...
          if ( qn(iglb) > maxpid ) then
            write (csol,'("found more than ",i0," contributions to pixel")') maxpid; call csoErr
            TRACEBACK; status=1; return
          end if
          ! store pe:
          qidi(:,qn(iglb),iglb) = idi_all(:,ipix_all)
        end do

        !! testing ...
        !print *, 'xxx0 pixels on multiple domains:'
        !do iglb = 1, self%nglb
        !  if ( qn(iglb) > 1 ) then
        !    print *, '--x0  ', iglb, qn(iglb)
        !  end if
        !end do
        
        ! number of exchanges:
        nex_max = count( qn > 1 )
        !! info ...
        !write (csol,'("counted ",i0," exchanges between domains")') nex_max; call csoPr
        
        ! pixel indices to send and receive:
        allocate( exchi(nex_max,2), source=-999, stat=status )
        IF_NOT_OK_RETURN(status=1)

      end if ! root

      !! testing ...
      !print *, 'xxx1 pairs of pe ...'
      ! loop over pe pairs (zero based id's):
      do pid1 = 0, csoc%npes-2
        do pid2 = pid1+1, csoc%npes-1
        
          ! unique tags for this pair:
          itag1 = pid1 * csoc%npes + pid2
          itag2 = pid2 * csoc%npes + pid1

          ! define exchange info on root:
          if ( csoc%root ) then
            !! testing ...
            !print *, '--x1 pair ', pid1, pid2

            ! number of pixels to be exchanged:
            nex = 0
            ! loop over pixels:
            do iglb = 1, self%nglb
              ! should be on multiple pe ...
              if ( qn(iglb) > 1 ) then
                ! search for first pe:
                i1 = -999
                do k = 1, qn(iglb)
                  if ( qidi(1,k,iglb) == pid1 ) then
                    i1 = k
                    exit
                  end if
                end do
                ! search for second pe:
                i2 = -999
                do k = 1, qn(iglb)
                  if ( qidi(1,k,iglb) == pid2 ) then
                    i2 = k
                    exit
                  end if
                end do
                ! contributions from this pair?
                if ( (i1 > 0) .and. (i2 > 0) ) then
                  !! testing ...
                  !print *, '  -- iglb ', iglb
                  ! increase counter:
                  nex = nex + 1
                  ! fill pixel indices:
                  exchi(nex,1) = qidi(2,i1,iglb)
                  exchi(nex,2) = qidi(2,i2,iglb)
                end if
              end if ! qn > 1
            end do ! iglb
            !! testing ..
            !print *, '  -- nex = ', nex
            
            ! send/recv from here (root)?
            if ( csoc%id == pid1 ) then
              !! info ...
              !write (csol,'(a,":   exchange ",i0," contributions with pe ",i0)') rname, nex, pid2; call csoPr
              ! init storage for exchange with pid2:
              call self%exch(pid2)%Alloc( nex, status )
              IF_NOT_OK_RETURN(status=1)
              ! any ?
              if ( nex > 0 ) then
                ! copy:
                self%exch(pid2)%ipix = exchi(1:nex,1)
                !! testing ...
                !do k = 1, self%exch(pid2)%nex
                !  ipix = self%exch(pid2)%ipix(k)
                !  print *, '  -- exchange ', k, ' ipix ', ipix, '(global ', self%iglb(ipix), ')'
                !end do ! k
              end if
            else
              ! send to first of pair:
              !~ size:
              call csoc%SendRecv( nex, csoc%root_id, pid1, itag1, status )
              IF_NOT_OK_RETURN(status=1)
              !~ pixel indices:
              if ( nex > 0 ) then
                call csoc%SendRecv( exchi(1:nex,1), csoc%root_id, pid1, itag1, status )
                IF_NOT_OK_RETURN(status=1)
              end if
            end if
            
            ! send to second of pair ..
            !~ size:
            call csoc%SendRecv( nex, csoc%root_id, pid2, itag2, status )
            IF_NOT_OK_RETURN(status=1)
            !~ pixel indices:
            if ( nex > 0 ) then
              call csoc%SendRecv( exchi(1:nex,2), csoc%root_id, pid2, itag2, status )
              IF_NOT_OK_RETURN(status=1)
            end if
            
          else if ( csoc%id == pid1 ) then
          
            ! receive number of exchanges:
            call csoc%SendRecv( nex, csoc%root_id, pid1, itag1, status )
            IF_NOT_OK_RETURN(status=1)          
            ! info ...
            write (csol,'(a,":   exchange ",i0," contributions with pe ",i0)') rname, nex, pid2; call csoPr
            ! storage:
            call self%exch(pid2)%Alloc( nex, status )
            IF_NOT_OK_RETURN(status=1)
            ! receive content:
            if ( nex > 0 ) then
              call csoc%SendRecv( self%exch(pid2)%ipix, csoc%root_id, pid1, itag1, status )
              IF_NOT_OK_RETURN(status=1)
              !! testing ...
              !do k = 1, self%exch(pid1)%nex
              !  ipix = self%exch(pid1)%ipix(k)
              !  print *, '  -- exchange ', k, ' ipix ', ipix, '(global ', self%iglb(ipix), ')'
              !end do ! k
            end if
              
          else if ( csoc%id == pid2 ) then
          
            ! send to first of pair ..
            call csoc%SendRecv( nex, csoc%root_id, pid2, itag2, status )
            IF_NOT_OK_RETURN(status=1)          
            !! info ...
            !write (csol,'(a,":   exchange ",i0," contributions with pe ",i0)') rname, nex, pid1; call csoPr
            ! storage:
            call self%exch(pid1)%Alloc( nex, status )
            IF_NOT_OK_RETURN(status=1)
            ! receive content:
            if ( nex > 0 ) then
              call csoc%SendRecv( self%exch(pid1)%ipix, csoc%root_id, pid2, itag2, status )
              IF_NOT_OK_RETURN(status=1)
              !! testing ...
              !do k = 1, self%exch(pid1)%nex
              !  ipix = self%exch(pid1)%ipix(k)
              !  print *, '  -- exchange ', k, ' ipix ', ipix, '(global ', self%iglb(ipix), ')'
              !end do ! k
            end if
            
          end if ! root, or pid1 or pid2

        end do ! pid2
      end do ! pid1

      ! clear:
      if ( csoc%root ) then
        ! clear:
        deallocate( exchi, stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! clear:
        deallocate( qn, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( qidi, stat=status )
        IF_NOT_OK_RETURN(status=1)

      end if  ! root

      ! clear:
      deallocate( idi_loc, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( idi_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! npes > 1
    
    ! set flag:
    self%exchange_is_setup = .true.
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_Data_SetupExchange


  ! ***


!  !
!  ! Fill array with local pixel size with random numbers out of N(0,1).
!  ! An array with length 'nglb' will be filled with random numbers first;
!  ! this array will be the same on all domains if the provided random generator
!  ! was initialized with the same seed on all domains.
!  ! The local pixel selection is then extracted from this array.
!  !
!
!  subroutine CSO_Sat_Data_GetRandom( self, rnd, v, status )
!  
!    use Num, only : T_Random
!    
!    ! --- in/out ---------------------------------
!    
!    class(T_CSO_Sat_Data), intent(inout)    ::  self
!    type(T_Random), intent(inout)                 ::  rnd
!    real, intent(out)                             ::  v(:)  ! (npix)
!    integer, intent(out)                          ::  status
!  
!    ! --- const ----------------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_Data_GetRandom'
!    
!    ! --- local ----------------------------------
!    
!    real, allocatable     ::  vglb(:)  ! (nglb)
!    integer               ::  iglb
!    integer               ::  ipix
!    
!    ! --- begin ----------------------------------
!    
!    ! check ...
!    if ( size(v) < self%npix ) then
!      write (csol,'("size v is ",i0," while nix is ",i0)') size(v), ipix; call csoErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! storage:
!    allocate( vglb(self%nglb), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! loop over global number of pixels:
!    do iglb = 1, self%nglb
!      ! generate random number out of N(0,1) distribution:
!      call rnd%Get_Normal( vglb(iglb), status )
!      IF_NOT_OK_RETURN(status=1)
!    end do ! i
!    
!    ! loop over local pixels:
!    do ipix = 1, self%npix
!      ! current global index:
!      iglb = self%iglb(ipix)
!      ! copy:
!      v(ipix) = vglb(iglb)
!    end do ! ipix
!    
!    ! clear:
!    deallocate( vglb, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! ok
!    status = 0
!  
!  end subroutine CSO_Sat_Data_GetRandom


  ! ====================================================================
  ! ===
  ! === State
  ! ===
  ! ====================================================================


  subroutine CSO_Sat_State_Init( self, sdata, rcF, rcbase, status, &
                                   description )
  
    use CSO_Rc              , only : T_CSO_RcFile
    use CSO_String          , only : CSO_ReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(out)         ::  self
    type(T_CSO_Sat_Data), intent(in)            ::  sdata
    type(T_CSO_RcFile), intent(in)              ::  rcF
    character(len=*), intent(in)                ::  rcbase
    integer, intent(out)                        ::  status
    
    character(len=*), intent(in), optional      ::  description

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_State_Init'
    
    ! --- local ----------------------------------
    
    character(len=1024)       ::  line
    integer                   ::  ivar
    integer                   ::  iudim
    integer                   ::  mxdim
    integer                   ::  i
    character(len=32)         ::  vname
    character(len=256)        ::  xnames
    character(len=256)        ::  aname
    character(len=1024)       ::  avalue
    character(len=256)        ::  units
    
    ! --- begin ----------------------------------
    
    ! defaults:
    self%description = 'simulated state'
    if ( present(description) ) self%description = trim(description)
    
    ! copy:
    self%npix   = sdata%npix
    
    ! any pixels in global domain?
    if ( sdata%nglb > 0 ) then

      ! info ...
      write (csol,'(a,": define state variables ...")') rname; call csoPr

      ! init extra output:
      call self%pd%Init( status )
      IF_NOT_OK_RETURN(status=1)

      ! line with variable names:
      call rcF%Get( trim(rcbase)//'.vars', line, status )
      IF_NOT_OK_RETURN(status=1)
      ! loop over elements:
      do
        ! empty?
        if ( len_trim(line) == 0 ) exit
        ! next part:
        call CSO_ReadFromLine( line, vname, status, sep=' ' )
        IF_NOT_OK_RETURN(status=1)
        
        ! dimension names:
        call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.dims', xnames, status )
        IF_NOT_OK_RETURN(status=1)
        ! info ...
        write (csol,'(a,":  variable `",a,"` ...")') rname, trim(vname); call csoPr
        !! testing ...
        !write (csol,'(a,":    dimensions : ",a)') rname, trim(xnames); call csoPr

        ! define new variable, return index:
        call self%pd%Def( vname, xnames, ivar, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! attribute names:
        call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.attrs', xnames, status )
        IF_NOT_OK_RETURN(status=1)
        ! loop:
        do
          ! empty?
          if ( len_trim(xnames) == 0 ) exit
          ! next name:
          call CSO_ReadFromLine( xnames, aname, status, sep=' ' )
          IF_NOT_OK_RETURN(status=1)
          ! attribute value:
          call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.attr.'//trim(aname), avalue, status )
          IF_NOT_OK_RETURN(status=1)
          ! store:
          call self%pd%SetAttr( ivar, aname, avalue, status )
          IF_NOT_OK_RETURN(status=1)
        end do ! attributes
        
        ! formula defined?
        call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.formula', aname, status, default='' )
        IF_ERROR_RETURN(status=1)
        ! defined?
        if ( len_trim(aname) > 0 ) then
          ! read formula terms:
          call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.formula_terms', avalue, status )
          IF_NOT_OK_RETURN(status=1)
          !! info ...
          !write (csol,'(a,":    formula    : ",a)') rname, trim(aname); call csoPr
          !write (csol,'(a,":      terms    : ",a)') rname, trim(avalue); call csoPr
          ! store:
          call self%pd%SetFormula( ivar, aname, avalue, status )
          IF_NOT_OK_RETURN(status=1)
        end if ! formula
          

      end do ! elements

      ! extract counter for user defined variables:
      call self%pd%Get( status, n=self%nvar )
      IF_NOT_OK_RETURN(status=1)

      ! extract total number of dimensions:
      call self%pd%Get( status, ndim=mxdim )
      IF_NOT_OK_RETURN(status=1)
      ! maximum storage for user defined dim names:
      allocate( self%udimnames(mxdim), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! init counter:
      iudim = 0
      ! loop over all dims:
      do i = 1, mxdim
        ! dimension name:
        call self%pd%GetDim( i, status, name=vname )
        IF_NOT_OK_RETURN(status=1)
        ! default or user defined?
        select case ( trim(vname) )
          !~ number of apriori layers:
          case ( 'layer' )
            ! set length:
            call self%pd%SetDim( vname, sdata%nlayer, status )
            IF_NOT_OK_RETURN(status=1)
          !~ number of retrieval layers:
          case ( 'retr' )
            ! set length:
            call self%pd%SetDim( vname, sdata%nretr, status )
            IF_NOT_OK_RETURN(status=1)
          !~ user defined:
          case default
            ! increase counter:
            iudim = iudim + 1
            ! store:
            self%udimnames(iudim) = trim(vname)
        end select
      end do ! dims
      ! store number of user defined dims:
      self%nudim = iudim
      
      ! output selections:
      call self%outkeys%Init( rcF, rcbase, status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! any global pixels
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_Init


  ! ***

  
  !
  ! Init state (defined on domain "doms_f") as swapped version of "sstate" (defined on "doms"),
  ! use the info in data "sdata_f" (already defined on "doms_f")
  !

  subroutine CSO_Sat_State_InitSwap( self, sdata_f, sstate, status )
  
    use CSO_Comm    , only : csoc
!    use CSO_Domains , only : T_CSO_Domains
!    use CSO_Swapping, only : T_Swapping
  
!    use CSO_String          , only : CSO_ReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(out)         ::  self
    type(T_CSO_Sat_Data), intent(in)            ::  sdata_f
    type(T_CSO_Sat_State), intent(in)           ::  sstate
!    type(T_CSO_Domains), intent(in)             ::  doms
!    type(T_CSO_Domains), intent(in)             ::  doms_f
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_State_InitSwap'
    
    ! --- local ----------------------------------
    
!    type(T_Swapping)          ::  swp
!    character(len=1024)       ::  line
!    integer                   ::  ivar
!    integer                   ::  iudim
!    integer                   ::  mxdim
!    integer                   ::  i
!    character(len=32)         ::  vname
!    character(len=256)        ::  xnames
!    character(len=256)        ::  aname
!    character(len=1024)       ::  avalue
!    character(len=256)        ::  units
    
    ! --- begin ----------------------------------
    
    ! copy:
    self%description = trim(sstate%description)
    
    ! copy:
    self%npix   = sdata_f%npix
    
    ! any pixels in global domain?
    if ( sdata_f%nglb > 0 ) then

      ! info ...
      write (csol,'(a,": define state variables ...")') rname; call csoPr
    
!      ! init extra output:
!      call self%pd%Init( status )
!      IF_NOT_OK_RETURN(status=1)
!
!      ! line with variable names:
!      call rcF%Get( trim(rcbase)//'.vars', line, status )
!      IF_NOT_OK_RETURN(status=1)
!      ! loop over elements:
!      do
!        ! empty?
!        if ( len_trim(line) == 0 ) exit
!        ! next part:
!        call CSO_ReadFromLine( line, vname, status, sep=' ' )
!        IF_NOT_OK_RETURN(status=1)
!        
!        ! dimension names:
!        call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.dims', xnames, status )
!        IF_NOT_OK_RETURN(status=1)
!        ! info ...
!        write (csol,'(a,":  variable `",a,"` ...")') rname, trim(vname); call csoPr
!        !! testing ...
!        !write (csol,'(a,":    dimensions : ",a)') rname, trim(xnames); call csoPr
!
!        ! define new variable, return index:
!        call self%pd%Def( vname, xnames, ivar, status )
!        IF_NOT_OK_RETURN(status=1)
!        
!        ! attribute names:
!        call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.attrs', xnames, status )
!        IF_NOT_OK_RETURN(status=1)
!        ! loop:
!        do
!          ! empty?
!          if ( len_trim(xnames) == 0 ) exit
!          ! next name:
!          call CSO_ReadFromLine( xnames, aname, status, sep=' ' )
!          IF_NOT_OK_RETURN(status=1)
!          ! attribute value:
!          call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.attr.'//trim(aname), avalue, status )
!          IF_NOT_OK_RETURN(status=1)
!          ! store:
!          call self%pd%SetAttr( ivar, aname, avalue, status )
!          IF_NOT_OK_RETURN(status=1)
!        end do ! attributes
!        
!        ! formula defined?
!        call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.formula', aname, status, default='' )
!        IF_ERROR_RETURN(status=1)
!        ! defined?
!        if ( len_trim(aname) > 0 ) then
!          ! read formula terms:
!          call rcF%Get( trim(rcbase)//'.var.'//trim(vname)//'.formula_terms', avalue, status )
!          IF_NOT_OK_RETURN(status=1)
!          !! info ...
!          !write (csol,'(a,":    formula    : ",a)') rname, trim(aname); call csoPr
!          !write (csol,'(a,":      terms    : ",a)') rname, trim(avalue); call csoPr
!          ! store:
!          call self%pd%SetFormula( ivar, aname, avalue, status )
!          IF_NOT_OK_RETURN(status=1)
!        end if ! formula
!
!      end do ! elements

      ! any pixels on some domain?
      if ( sdata_f%slc%nsel_tot > 0 ) then
      
        ! check ..
        if ( .not. sdata_f%swap_is_setup ) then
          write (csol,'("swap of data is not setup yet; sdata structure not swapped yet?")'); call csoErr
          TRACEBACK; status=1; return
        end if

        ! swap pixel data to target decomposition,
        ! use the swapping info that was setup in "sdata_f":
        call self%pd%InitSwap( sstate%pd, sdata_f%swp, status )
        IF_NOT_OK_RETURN(status=1)
        
      end if ! any pixels on some domain

!      ! extract counter for user defined variables:
!      call self%pd%Get( status, n=self%nvar )
!      IF_NOT_OK_RETURN(status=1)
!
!      ! extract total number of dimensions:
!      call self%pd%Get( status, ndim=mxdim )
!      IF_NOT_OK_RETURN(status=1)
!      ! maximum storage for user defined dim names:
!      allocate( self%udimnames(mxdim), stat=status )
!      IF_NOT_OK_RETURN(status=1)
!      ! init counter:
!      iudim = 0
!      ! loop over all dims:
!      do i = 1, mxdim
!        ! dimension name:
!        call self%pd%GetDim( i, status, name=vname )
!        IF_NOT_OK_RETURN(status=1)
!        ! default or user defined?
!        select case ( trim(vname) )
!          !~ number of apriori layers:
!          case ( 'layer' )
!            ! set length:
!            call self%pd%SetDim( vname, sdata%nlayer, status )
!            IF_NOT_OK_RETURN(status=1)
!          !~ number of retrieval layers:
!          case ( 'retr' )
!            ! set length:
!            call self%pd%SetDim( vname, sdata%nretr, status )
!            IF_NOT_OK_RETURN(status=1)
!          !~ user defined:
!          case default
!            ! increase counter:
!            iudim = iudim + 1
!            ! store:
!            self%udimnames(iudim) = trim(vname)
!        end select
!      end do ! dims
!      ! store number of user defined dims:
!      self%nudim = iudim
!      
!      ! output selections:
!      call self%outkeys%Init( rcF, rcbase, status )
!      IF_NOT_OK_RETURN(status=1)
      
      ! output selections:
      call self%outkeys%InitCopy( sstate%outkeys, status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! any global pixels
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_InitSwap


  ! ***


  subroutine CSO_Sat_State_Done( self, sdata, status )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    type(T_CSO_Sat_Data), intent(in)              ::  sdata
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! any pixels in global domain?
    if ( sdata%nglb > 0 ) then
    
      ! done with extra output:
      call self%pd%Done( status )
      IF_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( self%udimnames, stat=status )
      IF_NOT_OK_RETURN(status=1)

    end if
      
    ! clear output selections:
    call self%outkeys%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_Done


  ! ***


  ! Return info:
  !   - nudim      : number of user-defined dimensions to be filled
  !   - nuvar      : number of user-defined variables (all, or for specified outkey)
  !   - uvarnames  : names of user-defined variables (all, or for specified outkey)
  !   - uvarunits  : units of user-defined variables (all, or for specified outkey)
  !   - outkey     : optional limitation to specific output set

  subroutine CSO_Sat_State_Get( self, status, nudim, nuvar, uvarnames, uvarunits, outkey )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(in)            ::  self
    integer, intent(out)                          ::  status
    
    integer, intent(out), optional                ::  nudim
    integer, intent(out), optional                ::  nuvar
    character(len=*), intent(out), optional       ::  uvarnames(:)
    character(len=*), intent(out), optional       ::  uvarunits(:)
    character(len=*), intent(in), optional        ::  outkey
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_Get'
    
    ! --- local ----------------------------------

    character(len=32), pointer  ::  vars(:)
    
    ! --- begin ----------------------------------
    
    ! user defined dims:
    if ( present(nudim) ) nudim = self%nudim
    
    ! optional: pointer to array with selected variables:
    if ( present(outkey) ) then
      ! assign pointer:
      call self%outkeys%GetSelection( outkey, vars, status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! undefined:
      nullify( vars )
    end if
    
    ! user defined variable names and units:
    call self%pd%Get( status, nuvar=nuvar, uvarnames=uvarnames, uvarunits=uvarunits, vars=vars )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_Get


  ! ***


  subroutine CSO_Sat_State_GetDim( self, id, status, name )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(in)            ::  self
    integer, intent(in)                           ::  id
    integer, intent(out)                          ::  status
    
    character(len=*), intent(out), optional       ::  name
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_GetDim'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (id < 1) .or. (id > self%nudim) ) then
      write (csol,'("dimension ",i0," out of range 1:",i0)') id, self%nudim; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! user defined dim?
    if ( present(name) ) name = trim(self%udimnames(id))
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_GetDim


  ! ***


  subroutine CSO_Sat_State_SetDim( self, name, length, status )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    character(len=*), intent(in)                  ::  name
    integer, intent(in)                           ::  length
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_SetDim'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! set length by name:
    call self%pd%SetDim( name, length, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_SetDim


  ! ***


  subroutine CSO_Sat_State_EndDef( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_EndDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! allocate storage for user defined output:
    call self%pd%EndDef( self%npix, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_EndDef
  
  
  ! ***


  !
  ! Get info on pixel data specified by either "id" number or "name"
  !

  subroutine CSO_Sat_State_GetData( self, status, id, name,&
                                     data0, data1, data2, units )
  

    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(in)          ::  self
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional             ::  id
    character(len=*), intent(in), optional    ::  name

    real, pointer, optional                   ::  data0(:)     ! (npix)
    real, pointer, optional                   ::  data1(:,:)   ! (:,npix)
    real, pointer, optional                   ::  data2(:,:,:) ! (:,:,npix)
    character(len=*), intent(out), optional   ::  units

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_State_GetData'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! get units:
    call self%pd%GetData( status, id=id, name=name, &
                            data0=data0, data1=data1, data2=data2, &
                            units=units )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_GetData
  
  
  ! ***


  !
  ! Get gathered data on root, specified by either "id" number or "name"
  !

  subroutine CSO_Sat_State_GatherData( self, sdata, status, id, name, gdata1 )
  

    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(in)          ::  self
    type(T_CSO_Sat_Data), intent(in)            ::  sdata
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional             ::  id
    character(len=*), intent(in), optional    ::  name

    real, intent(out), optional               ::  gdata1(:,:)   ! (:,nout)

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Sat_State_GatherData'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! get units:
    call self%pd%GatherData( sdata%slc%iout_tot, status, id=id, name=name, &
                                gdata1=gdata1 )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_GatherData


  ! ***


  subroutine CSO_Sat_State_GetPixel( self, sdata, ipix, status, &
                                            y )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(in)            ::  self
    type(T_CSO_Sat_Data), intent(in)              ::  sdata
    integer, intent(in)                           ::  ipix
    integer, intent(out)                          ::  status
    
    real, pointer, optional                       ::  y(:)    ! (nr)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_GetPixel'
    
    ! --- local ----------------------------------
    
    real, pointer     ::  ydata(:,:)
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ipix < 1) .or. (ipix > sdata%npix) ) then
      write (csol,'("pixel ",i0," out of range 1:",i0)') ipix, sdata%npix; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! simulatied retrieval:
    if ( present(y) ) then
      ! get pointer:
      call self%pd%GetData( status, name='y', data1=ydata )
      IF_NOT_OK_RETURN(status=1)
      ! assign:
      y => ydata(:,ipix)
    end if
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_GetPixel


  ! ***


  subroutine CSO_Sat_State_SetPixel( self, sdata, ipix, status, &
                                            zero, &
                                            ivar, uval, uprof )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    type(T_CSO_Sat_Data), intent(in)              ::  sdata
    integer, intent(in)                           ::  ipix
    integer, intent(out)                          ::  status
    
    logical, intent(in), optional                 ::  zero
    
    integer, intent(in), optional                 ::  ivar
    real, intent(in), optional                    ::  uval
    real, intent(in), optional                    ::  uprof(:)
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_SetPixel'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (ipix < 1) .or. (ipix > sdata%npix) ) then
      write (csol,'("pixel ",i0," out of range 1:",i0)') ipix, sdata%npix; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! zero flag present?
    if ( present(zero) ) then
      ! reset values to zero?
      if ( zero ) then
        ! user data:
        call self%pd%SetPixelZero( ipix, status )
        IF_NOT_OK_RETURN(status=1)
      end if
    end if
    
    ! user defined value? requires ivar ...
    if ( present(uval) ) then
      ! check ..
      if ( .not. present(ivar) ) then
        write (csol,'("missing argument `ivar` with argument `uval`")'); call csoErr
        TRACEBACK; status=1; return
      end if
      ! fill:
      call self%pd%SetPixel( ivar, ipix, uval, status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! user defined value? requires ivar ...
    if ( present(uprof) ) then
      ! check ..
      if ( .not. present(ivar) ) then
        write (csol,'("missing argument `ivar` with argument `uprof`")'); call csoErr
        TRACEBACK; status=1; return
      end if
      ! fill:
      call self%pd%SetPixel( ivar, ipix, uprof, status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_SetPixel


  ! ***

  
  ! Exchange pixel contributions between domains
  !   - outkey     : optional limitation to specific output set
  
  subroutine CSO_Sat_State_Exchange( self, sdata, status, &
                                       outkey )
  
    use CSO_Comm, only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    type(T_CSO_Sat_Data), intent(in)              ::  sdata
    integer, intent(out)                          ::  status
    character(len=*), intent(in), optional        ::  outkey

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_Exchange'
    
    ! --- local ----------------------------------

    character(len=32), pointer  ::  vars(:)

    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. sdata%exchange_is_setup ) then
      write (csol,'("exchange information not setup yet ...")'); call csoPr
      TRACEBACK; status=1; return
    end if

    ! exchange needed?
    if ( csoc%npes > 1 ) then
    
      ! optional: pointer to array with selected variables:
      if ( present(outkey) ) then
        ! assign pointer:
        call self%outkeys%GetSelection( outkey, vars, status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! undefined:
        nullify( vars )
      end if

      ! apply:
      call self%pd%Exchange( sdata%exch, status ) !, vars=vars )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_Exchange


  ! ***

  
  ! Exchange pixel contributions between domains
  !   - outkey     : optional limitation to specific output set
  
  subroutine CSO_Sat_State_ApplyFormulas( self, sdata, status, &
                                            outkey )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    type(T_CSO_Sat_Data), intent(in)              ::  sdata
    integer, intent(out)                          ::  status
    character(len=*), intent(in), optional        ::  outkey

    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_ApplyFormulas'
    
    ! --- local ----------------------------------

    character(len=32), pointer  ::  vars(:)

    ! --- begin ----------------------------------
    
    ! optional: pointer to array with selected variables:
    if ( present(outkey) ) then
      ! assign pointer:
      call self%outkeys%GetSelection( outkey, vars, status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! undefined:
      nullify( vars )
    end if

    ! apply:
    call self%pd%ApplyFormulas( sdata%pd, status ) !, vars=vars )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_ApplyFormulas


  ! ***

  
  !
  ! Apply transposed convolution for gradient forcing:
  !     hf =  K^T  f
  ! where:
  !     f is adjoint forcing in units [1/(y_units)]
  ! Eventually including conversion from 1/ppb to 1/(mol/m2):
  !     hf = K^T  f / airm * 1e9
  !

  subroutine CSO_Sat_State_GetForcing( self, sdata, ipix, hf, status, verbose )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    type(T_CSO_Sat_Data), intent(in)              ::  sdata
    integer, intent(in)                           ::  ipix
    real, intent(out)                             ::  hf(:)  ! (nlayer)
    integer, intent(out)                          ::  status
    
    logical, intent(in), optional                 ::  verbose
    
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_GetForcing'
    
    ! --- local ----------------------------------
    
    logical         ::  verb
    integer         ::  i, j
    real, pointer   ::    y_data(:,:)    ! (nr,npix)
    real, pointer   ::    K_data(:,:,:)  ! (nr,na,npix)
    real, pointer   :: airm_data(:,:)    ! (nr,npix)
    
    ! --- begin ----------------------------------
    
    write (csol,'("update routine needed")'); call csoErr
    TRACEBACK; status=1; return
    ! dummy ...
    hf = 0.0
    
!    ! verbose?
!    verb = .false.
!    if ( present(verbose) ) verb = verbose
!    
!    ! check ...
!    if ( (ipix < 1) .or. (ipix > sdata%npix) ) then
!      write (csol,'("pixel ",i0," out of range 1:",i0)') ipix, sdata%npix; call csoErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! check ...
!    if ( size(hf) /= sdata%nlayer ) then
!      write (csol,'("input argument `hf` (",i0,") should have size `layer` (",i0,")")') size(hf), sdata%nlayer; call csoErr
!      TRACEBACK; status=1; return
!    end if
!    
!    ! pointer to simulated retrieval:
!    call self%pd%GetData( status, name='y', data1=y_data )
!    IF_NOT_OK_RETURN(status=1)
!    ! pointer to kernel:
!    call sdata%pd%GetData( status, name='K', data2=K_data )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! conversion needed?
!    if ( len_trim(sdata%conversion) > 0 ) then    
!      !
!      ! switch:
!      select case ( trim(sdata%conversion ) )
!        !~ vcd to vmr (in forward direction!)
!        case ( 'mol m-2 -> 1e-9' )
!          ! pointer to airmass:
!          call sdata%pd%GetData( status, name='airm', data1=airm_data )
!          IF_NOT_OK_RETURN(status=1)
!          ! loop over apri layers:
!          do j = 1, sdata%nlayer
!            ! init sum:
!            hf(j) = 0.0
!            ! loop over retrieval layers:
!            do i = 1, sdata%nr
!              ! devide by airmass, convert to ppb:
!              !                                      1/ppb      /   (mol air)/m2      ppb
!              hf(j) = hf(j) + K_data(i,j,ipix) * y_data(i,ipix) / airm_data(i,ipix) * 1e9  ! 1/((mol tr)/m2)
!              !! testing ...
!              !if (verb) print *, 'kkk1 ', i, j, sdata%K%data(i,j,ipix), sdata%K%data(i,j,ipix) / sdata%airm%data(i,ipix) * 1e9
!            end do ! i
!          end do ! j
!        !~ unknown
!        case default
!          write (csol,'("unsupported coversion: ",a)') trim(sdata%conversion); call csoErr
!          TRACEBACK; status=1; return
!      end select
!      !
!    else
!      !
!      ! loop over apri layers:
!      do j = 1, sdata%nlayer
!        ! init sum:
!        hf(j) = 0.0
!        ! loop over retrieval layers:
!        do i = 1, sdata%nr
!          ! devide by airmass, convert to ppb:
!          !                                    1/y_units
!          hf(j) = hf(j) + K_data(i,j,ipix) * y_data(i,ipix)
!        end do ! i
!      end do ! j
!      !
!    end if ! conversion
    
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_GetForcing


  ! ***


  subroutine CSO_Sat_State_PutOut( self, sdata, filename, status, &
                                      outkey )

    use NetCDF , only : NF90_Create, NF90_Close
    use NetCDF , only : NF90_Def_Dim
    use NetCDF , only : NF90_EndDef
  
    use CSO_Comm  , only : csoc
    use CSO_NcFile, only : T_NcFile
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    class(T_CSO_Sat_Data), intent(in)             ::  sdata
    character(len=*), intent(in)                  ::  filename
    integer, intent(out)                          ::  status
    
    character(len=*), intent(in), optional        ::  outkey
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_PutOut'
    
    ! --- local ----------------------------------
    
    type(T_NcFile)              ::  ncf
    integer                     ::  dimid_pixel
    integer                     ::  dimid_retr
    integer                     ::  dimid_mod_layer
    integer                     ::  dimid_mod_layeri
    character(len=32), pointer  ::  vars(:)
    
    ! --- begin ----------------------------------
    
    ! info ..
    write (csol,'(a,": put out sat state to: ",a)') rname, trim(filename); call csoPr
    
    ! optional: pointer to array with selected variables:
    if ( present(outkey) ) then
      ! assign pointer:
      call self%outkeys%GetSelection( outkey, vars, status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! undefined:
      nullify( vars )
    end if

    ! any pixels to be put out?
    if ( sdata%slc%nout > 0 ) then
    
      ! collect on root, this is much faster than parallel write ...
      if ( csoc%root ) then

        ! create file:
        call ncf%Init( filename, 'w', status )
        IF_NOT_OK_RETURN(status=1)

        ! pixel dimension:
        status = NF90_Def_Dim( ncf%ncid, 'pixel', sdata%slc%nout, dimid_pixel )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! user output:
        call self%pd%NcDef( ncf%ncid, dimid_pixel, vars, status, &
                              deflate_level=sdata%deflate_level, &
                              packed=sdata%packed )
        IF_NOT_OK_RETURN(status=1)

        ! end defintion mode:
        status = NF90_EndDef( ncf%ncid )
        IF_NF90_NOT_OK_RETURN(status=1)
        
      end if ! root
      
      ! collect distributed arrays on root and write from there,
      !~ contributions from multiple domains should be added together:
      !call self%pd%NcPutGather( ncid, sdata%iout_all, .true., status )
      !IF_NOT_OK_RETURN(status=1)
      !~ contributions from different domains have been exchanged already,
      !  so no need to add contributions:
      call self%pd%NcPutGather( ncf, sdata%slc%iout_tot, .false., vars, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! written on root...
      if ( csoc%root ) then
        ! close:
        call ncf%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end if  ! root

    else
    
      ! info ..
      write (csol,'(a,":   no pixels to be put out ...")') rname; call csoPr
      
    end if

    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_PutOut


  ! ***

  !
  ! Read local pixels.
  !

  subroutine CSO_Sat_State_ReadForcing( self, sdata, filename, status )

    use NetCDF  , only : NF90_NOWRITE
    use NetCDF  , only : NF90_Open, NF90_Close
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------
    
    class(T_CSO_Sat_State), intent(inout)         ::  self
    class(T_CSO_Sat_Data), intent(in)             ::  sdata
    character(len=*), intent(in)                  ::  filename
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Sat_State_ReadForcing'
    
    ! --- local ----------------------------------
    
    integer                   ::  ncid
    logical                   ::  exist
    character(len=64)         ::  varname
    
    ! --- begin ----------------------------------

    ! info ..
    write (csol,'(a,": read forcing from: ",a)') rname, trim(filename); call csoPr
    
    ! read on root:
    if ( csoc%root ) then

      ! check ...
      inquire( file=trim(filename), exist=exist )
      if ( .not. exist ) then
        write (csol,'("file not found: ",a)') trim(filename); call csoErr
        TRACEBACK; status=1; return
      end if

      ! open file:
      status = NF90_Open( trim(filename), NF90_NOWRITE, ncid )
      IF_NF90_NOT_OK_RETURN(status=1)
      
    else
      ! dummy ..
      ncid = -999
    end if ! root

    ! variable with forcing field:
    varname = 'forcing'
    !! read array, scatter to other domains:
    !call self%y%NcGetScatter( ncid, varname, sdata%iout_all, status )
    !IF_NOT_OK_RETURN(status=1)
    write (csol,'("no pd%NcGetScatter yet")'); call csoErr
    TRACEBACK; status=1; return

    ! read on root:
    if ( csoc%root ) then
      ! close:
      status = NF90_Close( ncid )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if ! root
      
    ! ok
    status = 0
    
  end subroutine CSO_Sat_State_ReadForcing


end module CSO_Sat

