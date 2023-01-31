!#######################################################################
!
! CSO_Comm - communication through MPI.
!
!
! MACRO'S
!
!   _MPI       : should be defined to enable MPI library
!
!
! SEE ALSO
!   MPI Tutorial:
!     https://computing.llnl.gov/tutorials/mpi/
!
!
! HISTORY
!
!   2023-01, Arjo Segers
!     Support integer(1) and character variables.
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_MPI_NOT_OK_RETURN(action) if (status/=MPI_SUCCESS) then; errorcode=status; call MPI_Error_String(errorcode,csol,ncsol,status); call csoErr; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!#######################################################################

module CSO_Comm

  use CSO_Logging, only : csol, ncsol, csoPr, csoErr
#ifdef _MPI
  use MPI_F08, only : MPI_SUCCESS, MPI_Error_String
  use MPI_F08, only : MPI_Comm
  use MPI_F08, only : MPI_Info
#endif

  implicit none


  ! --- in/out -----------------------------------

  private

  public  ::  csoc


  ! --- const ------------------------------------

  character(len=*), parameter   ::  mname = 'CSO_Comm'

  ! character lengths:
  integer, parameter  ::  LEN_NAME = 32


  ! --- types-------------------------------------

  ! comunicator:
  type T_CSO_Comm
    !! own name:
    !character(len=LEN_NAME)         ::  name
    ! use external communicator?
    logical                          ::  with_external_comm
#ifdef _MPI
    ! MPI communicator;
    ! use 'comm%mpi_val' for the integer:
    type(MPI_Comm)                  ::  comm
    ! MPI info structure;
    ! use 'info%mpi_val' for the integer:
    type(MPI_Info)                  ::  info
#endif
    ! size and rank:
    integer                         ::  npes
    integer                         ::  id
    ! root?
    integer                         ::  root_id
    logical                         ::  root
    !! info for each of the applications, including own:
    !type(T_Appl), allocatable       ::  appl(:)   ! (0:npes-1)
    !
  contains
    procedure   ::  Init        =>  CSO_Comm_Init
    procedure   ::  Done        =>  CSO_Comm_Done
    !procedure   ::  SetName     =>  CSO_Comm_SetName
    !procedure   ::  GetID       =>  CSO_Comm_GetID
    procedure   ::  Abort       =>  CSO_Comm_Abort
    procedure   ::  Barrier     =>  CSO_Comm_Barrier
    !
#ifdef _MPI
    procedure   ::  GetDataType =>  CSO_Comm_GetDataType
    procedure   ::  GetOper     =>  CSO_Comm_GetOper
#endif
    !
    procedure                       CSO_Comm_SendRecv_i4
    procedure                       CSO_Comm_SendRecv_i4_1d
    procedure                       CSO_Comm_SendRecv_i4_2d
    procedure                       CSO_Comm_SendRecv_r4_2d
    procedure                       CSO_Comm_SendRecv_r8_2d
    generic     ::  SendRecv =>     CSO_Comm_SendRecv_i4, &
                                    CSO_Comm_SendRecv_i4_1d, &
                                    CSO_Comm_SendRecv_i4_2d, &
                                    CSO_Comm_SendRecv_r4_2d, &
                                    CSO_Comm_SendRecv_r8_2d
    !
    procedure                       CSO_Comm_Reduce_i
    generic     ::  Reduce      =>  CSO_Comm_Reduce_i
    !
    procedure                       CSO_Comm_AllReduce_i
    procedure                       CSO_Comm_AllReduce_r
    procedure                       CSO_Comm_AllReduce_l
    procedure                       CSO_Comm_AllReduce_r4_1d
    procedure                       CSO_Comm_AllReduce_r4_2d
    procedure                       CSO_Comm_AllReduce_r4_3d
    procedure                       CSO_Comm_AllReduce_r4_4d
    procedure                       CSO_Comm_AllReduce_r4_5d
    procedure                       CSO_Comm_AllReduce_InPlace_r4_5d
    generic     ::  AllReduce   =>  CSO_Comm_AllReduce_i, &
                                    CSO_Comm_AllReduce_r, &
                                    CSO_Comm_AllReduce_l, &
                                    CSO_Comm_AllReduce_r4_1d, &
                                    CSO_Comm_AllReduce_r4_2d, &
                                    CSO_Comm_AllReduce_r4_3d, &
                                    CSO_Comm_AllReduce_r4_4d, &
                                    CSO_Comm_AllReduce_r4_5d, &
                                    CSO_Comm_AllReduce_InPlace_r4_5d
    !
    procedure ::                    CSO_Comm_BCast_i
    procedure ::                    CSO_Comm_BCast_i1_1d
    procedure ::                    CSO_Comm_BCast_i_1d
    procedure ::                    CSO_Comm_BCast_r4
    procedure ::                    CSO_Comm_BCast_r4_1d
    procedure ::                    CSO_Comm_BCast_r4_2d
    procedure ::                    CSO_Comm_BCast_r4_3d
    procedure ::                    CSO_Comm_BCast_r4_4d
    procedure ::                    CSO_Comm_BCast_r8
    procedure ::                    CSO_Comm_BCast_r8_1d
    procedure ::                    CSO_Comm_BCast_r8_2d
    procedure ::                    CSO_Comm_BCast_r8_3d
    procedure ::                    CSO_Comm_BCast_r8_4d
    procedure ::                    CSO_Comm_BCast_c
    generic   ::  BCast         =>  CSO_Comm_BCast_i, &
                                    CSO_Comm_BCast_i1_1d, &
                                    CSO_Comm_BCast_i_1d, &
                                    CSO_Comm_BCast_r4, &
                                    CSO_Comm_BCast_r4_1d, &
                                    CSO_Comm_BCast_r4_2d, &
                                    CSO_Comm_BCast_r4_3d, &
                                    CSO_Comm_BCast_r4_4d, &
                                    CSO_Comm_BCast_r8, &
                                    CSO_Comm_BCast_r8_1d, &
                                    CSO_Comm_BCast_r8_2d, &
                                    CSO_Comm_BCast_r8_3d, &
                                    CSO_Comm_BCast_r8_4d, &
                                    CSO_Comm_BCast_c
    !
    procedure ::                    CSO_Comm_Gather_i
    procedure ::                    CSO_Comm_Gather_i_1d
    generic   ::  Gather        =>  CSO_Comm_Gather_i, &
                                    CSO_Comm_Gather_i_1d
    !
    procedure ::                    CSO_Comm_AllGather_i
    procedure ::                    CSO_Comm_AllGather_r
    procedure ::                    CSO_Comm_AllGather_i1
    procedure ::                    CSO_Comm_AllGather_r1
    generic   ::  AllGather     =>  CSO_Comm_AllGather_i, &
                                    CSO_Comm_AllGather_r, &
                                    CSO_Comm_AllGather_i1, &
                                    CSO_Comm_AllGather_r1
    !
    procedure                       CSO_Comm_GatherV_i1_1d
    procedure                       CSO_Comm_GatherV_i_1d
    procedure                       CSO_Comm_GatherV_i_2d
    procedure                       CSO_Comm_GatherV_r4_1d
    procedure                       CSO_Comm_GatherV_r8_1d
    procedure                       CSO_Comm_GatherV_c1_2d
    procedure                       CSO_Comm_GatherV_r4_2d
    procedure                       CSO_Comm_GatherV_r8_2d
    procedure                       CSO_Comm_GatherV_r4_3d
    procedure                       CSO_Comm_GatherV_r8_3d
    generic   ::  GatherV       =>  CSO_Comm_GatherV_i1_1d, &
                                    CSO_Comm_GatherV_i_1d, &
                                    CSO_Comm_GatherV_i_2d, &
                                    CSO_Comm_GatherV_r4_1d, &
                                    CSO_Comm_GatherV_r8_1d, &
                                    CSO_Comm_GatherV_c1_2d, &
                                    CSO_Comm_GatherV_r4_2d, &
                                    CSO_Comm_GatherV_r8_2d, &
                                    CSO_Comm_GatherV_r4_3d, &
                                    CSO_Comm_GatherV_r8_3d
    !
    procedure                       CSO_Comm_Gather2D_r4
    procedure                       CSO_Comm_Gather2D_r8
    generic   ::  Gather2D      =>  CSO_Comm_Gather2D_r4, &
                                    CSO_Comm_Gather2D_r8
    !
    procedure                       CSO_Comm_ScatterV_i1_1d
    procedure                       CSO_Comm_ScatterV_r4_1d
    procedure                       CSO_Comm_ScatterV_c1_2d
    procedure                       CSO_Comm_ScatterV_r4_2d
    procedure                       CSO_Comm_ScatterV_r4_3d
    procedure                       CSO_Comm_ScatterV_r8_1d
    procedure                       CSO_Comm_ScatterV_r8_2d
    procedure                       CSO_Comm_ScatterV_r8_3d
    generic   ::  ScatterV      =>  CSO_Comm_ScatterV_i1_1d, &
                                    CSO_Comm_ScatterV_r4_1d, &
                                    CSO_Comm_ScatterV_c1_2d, &
                                    CSO_Comm_ScatterV_r4_2d, &
                                    CSO_Comm_ScatterV_r4_3d, &
                                    CSO_Comm_ScatterV_r8_1d, &
                                    CSO_Comm_ScatterV_r8_2d, &
                                    CSO_Comm_ScatterV_r8_3d
    !
    procedure ::                    CSO_Comm_AllToAll_i4_1d
    generic   ::  AllToAll      =>  CSO_Comm_AllToAll_i4_1d
    !
    procedure ::                    CSO_Comm_AllToAllV_i4_1d
    procedure ::                    CSO_Comm_AllToAllV_r4_1d
    procedure ::                    CSO_Comm_AllToAllV_r8_1d
    procedure ::                    CSO_Comm_AllToAllV_r4_2d
    procedure ::                    CSO_Comm_AllToAllV_r8_2d
    generic   ::  AllToAllV     =>  CSO_Comm_AllToAllV_i4_1d, &
                                    CSO_Comm_AllToAllV_r4_1d, &
                                    CSO_Comm_AllToAllV_r8_1d, &
                                    CSO_Comm_AllToAllV_r4_2d, &
                                    CSO_Comm_AllToAllV_r8_2d
    !
    procedure ::  ParInfo       =>  CSO_Comm_ParInfo
    procedure ::  Displacements =>  CSO_Comm_Displacements
    !
  end type T_CSO_Comm


  ! --- var -------------------------------------

  ! mpi error code; used as argument for 'MPI_Error_String' to avoid
  ! warnings about same argument 'status' being used for both 'errorcode' and 'ierror':
  integer                     ::  errorcode

  ! communicator:
  type(T_CSO_Comm)            ::  csoc



contains


  ! ====================================================================
  ! ===
  ! === comm
  ! ===
  ! ====================================================================


  subroutine CSO_Comm_Init( self, status, comm )

#ifdef _MPI
    use MPI_F08, only : MPI_COMM_WORLD
    use MPI_F08, only : MPI_INFO_NULL
    use MPI_F08, only : MPI_Init
    use MPI_F08, only : MPI_Comm_Size, MPI_Comm_Rank
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(out)    ::  self
    integer, intent(out)              ::  status
    
#ifdef _MPI
    type(MPI_Comm), intent(in), optional   ::  comm
#else
    integer, intent(in), optional          ::  comm
#endif

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_Init'

    ! --- local ----------------------------------

    integer         ::  id

    ! --- begin ----------------------------------
    
    ! default flag:
    self%with_external_comm = .false.

#ifdef _MPI
    ! communicator provided?
    if ( present(comm) ) then
    
      ! set flag:
      self%with_external_comm = .true.
      ! copy:
      self%comm = comm
      
    else

      ! init MPI interface:
      call MPI_Init( ierror=status )
      if ( status /= MPI_SUCCESS ) then
        write (csol,'("could not initialize MPI environment; exit code : ",i6)') status; call csoErr
        TRACEBACK; status=1; return
      end if

      ! use global communicator:
      self%comm = MPI_COMM_WORLD
      
    end if

    ! dummy info:
    self%info = MPI_INFO_NULL

    ! size and rank:
    call MPI_Comm_Size( self%comm, self%npes, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    call MPI_Comm_Rank( self%comm, self%id, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

#else

    ! single pe:
    self%npes = 1
    self%id = 0

#endif

    !! storage for info on applications:
    !allocate( self%appl(0:self%npes-1), stat=status )
    !IF_NOT_OK_RETURN(status=1)

    !! empty name:
    !self%name = ''

    ! set root id:
    self%root_id = 0
    ! root?
    self%root = self%id == self%root_id

    ! ok
    status = 0

  end subroutine CSO_Comm_Init


  ! ***


  subroutine CSO_Comm_Done( self, status )

#ifdef _MPI
    use MPI_F08, only : MPI_Finalize
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

#ifdef _MPI
    ! only if not external ...
    if ( .not. self%with_external_comm ) then
      ! done with MPI interface:
      call MPI_Finalize( ierror=status )
      IF_MPI_NOT_OK_RETURN(status=1)
    end if
#endif

    !! clear:
    !deallocate( self%appl, stat=status )
    !IF_NOT_OK_RETURN(status=1)

    !! reset name:
    !self%name = ''

    ! ok
    status = 0

  end subroutine CSO_Comm_Done


!  ! ***
!
!
!  subroutine CSO_Comm_SetName( self, name, status )
!
!#ifdef _MPI
!    use MPI_F08, only : MPI_COMM_WORLD
!    use MPI_F08, only : MPI_CHARACTER
!#endif
!
!    ! --- in/out ---------------------------------
!
!    class(T_CSO_Comm), intent(out)       ::  self
!    character(len=*), intent(in)      ::  name
!    integer, intent(out)              ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_SetName'
!
!    ! --- local ----------------------------------
!
!    integer         ::  id
!
!    ! --- begin ----------------------------------
!
!    ! store name:
!    self%name = trim(name)
!
!    ! store own:
!    self%appl(self%id)%name = trim(name)
!#ifdef _MPI
!    ! exchange names; loop over all pes:
!    do id = 0, self%npes-1
!      ! send from this id to all other:
!      call MPI_BCast( self%appl(id)%name, len(self%appl(id)%name), MPI_CHARACTER, &
!                          id, self%comm, status )
!      IF_MPI_NOT_OK_RETURN(status=1)
!    end do
!#endif
!
!    ! info ...
!    write (csol,'(a," - connected to the following pes:")') trim(self%name); call csoPr
!    do id = 0, self%npes-1
!      if ( id == self%id ) cycle
!      write (csol,'(a," -   pe ",i3," with application `",a,"`")') &
!               trim(self%name), id, trim(self%appl(id)%name); call csoPr
!    end do
!
!    ! ok
!    status = 0
!
!  end subroutine CSO_Comm_SetName
!
!
!  ! ***
!
!
!  ! Return process ID of requested application.
!
!  subroutine CSO_Comm_GetID( self, name, id, status )
!
!    ! --- in/out ---------------------------------
!
!    class(T_CSO_Comm), intent(inout)    ::  self
!    character(len=*), intent(in)        ::  name
!    integer, intent(out)                ::  id
!    integer, intent(out)                ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GetID'
!
!    ! --- local ----------------------------------
!
!    integer           ::  k
!
!    ! --- begin ----------------------------------
!
!    ! init as dummy ..
!    id = -999
!    ! loop over applications:
!    do k = 0, self%npes-1
!      ! compare:
!      if ( trim(name) == trim(self%appl(k)%name) ) then
!        ! found:
!        id = k
!        ! leave:
!        exit
!      end if
!    end do   ! pe's
!    ! not found ?
!    if ( id < 0 ) then
!      write (csol,'("could not find application with name `",a,"`")') trim(name); call csoErr
!      write (csol,'("available applications:")'); call csoErr
!      do k = 0, self%npes-1
!        write (csol,'("  - ",a)') trim(self%appl(k)%name); call csoErr
!      end do
!      TRACEBACK; status=1; return
!    end if
!
!    ! ok
!    status = 0
!
!  end subroutine CSO_Comm_GetID


  ! ***


  subroutine CSO_Comm_Barrier( self, status )

#ifdef _MPI
    use MPI_F08, only : MPI_Barrier
#endif

    ! --- const ---------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_Barrier'

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    integer, intent(out)                ::  status

    ! --- local -----------------------------------

    ! --- begin -----------------------------------

#ifdef _MPI
    ! start:
    call MPI_Barrier( self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_Barrier


  ! ***


  subroutine CSO_Comm_Abort( self, errcode )

#ifdef _MPI
    use MPI_F08, only : MPI_Abort
#endif

    ! --- in/out ----------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    integer, intent(in)                 ::  errcode

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/CSO_Comm_Abort'

    ! --- local -----------------------------

    integer        ::  status

    ! --- begin -----------------------------

#ifdef _MPI
    ! abort, get status back since even this might fail ...
    call MPI_Abort( self%comm, errcode, ierror=status )
    if ( status /= MPI_SUCCESS ) then
      errorcode=status
      call MPI_Error_String( errorcode, csol, ncsol, status ); call csoErr
      TRACEBACK
      ! alternative abort via system exit (non-standard routine!);
      ! this might not kill all other processes:
      call Exit( status )
    end if
#else
    ! system exit (non-standard routine!)
    call Exit( errcode )
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_Abort


  ! ====================================================================


#ifdef _MPI

  !
  ! Return 'MPI_DataType' for specified type ('real', 'integer') and kind.
  ! Note that the result does not use the communicatator,
  ! but since this object is needed for each MPI call anayway it is an
  ! easy way to ship it ...
  !

  subroutine CSO_Comm_GetDataType( self, typ, knd, dtype, status )

    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_LOGICAL
    use MPI_F08, only : MPI_BYTE, MPI_SHORT, MPI_INTEGER
    use MPI_F08, only : MPI_REAL, MPI_DOUBLE_PRECISION
    use MPI_F08, only : MPI_CHAR

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    character(len=*), intent(in)        ::  typ
    integer, intent(in)                 ::  knd
    type(MPI_DataType), intent(out)     ::  dtype
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GetDataType'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! switch:
    select case ( trim(typ) )

      !~ logical:
      case ( 'logical' )
        ! set result:
        dtype = MPI_LOGICAL

      !~ integer:
      case ( 'integer' )
        ! switch:
        select case ( knd )
          case ( 1 )
            dtype = MPI_BYTE
          case ( 2 )
            dtype = MPI_SHORT
          case ( 4 )
            dtype = MPI_INTEGER
          case default
            write (csol,'("could not set MPI data type for ",a," kind ",i0)') trim(typ), knd; call csoErr
            TRACEBACK; status=1; return
        end select

      !~ reals:
      case ( 'real' )
        ! switch:
        select case ( knd )
          case ( 4 )
            dtype = MPI_REAL
          case ( 8 )
            dtype = MPI_DOUBLE_PRECISION
          case default
            write (csol,'("could not set MPI data type for ",a," kind ",i0)') trim(typ), knd; call csoErr
            TRACEBACK; status=1; return
        end select

      !~ characters:
      case ( 'char' )
        ! set result:
        dtype = MPI_CHAR

      !~
      case default
        write (csol,'("could not set MPI data type for ",a," variables")') trim(typ); call csoErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine CSO_Comm_GetDataType


  ! *


  !
  ! Translate collective operator description ('min','max')
  ! to an MPI_Op type.
  !

  subroutine CSO_Comm_GetOper( self, oper, op, status )

    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_MIN, MPI_MAX, MPI_SUM
    use MPI_F08, only : MPI_LAND, MPI_LOR

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)          ::  self
    character(len=*), intent(in)        ::  oper
    type(MPI_Op), intent(out)           ::  op
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GetOper'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! operator:
    select case ( oper )
      case ( 'min' ) ; op = MPI_MIN
      case ( 'max' ) ; op = MPI_MAX
      case ( 'sum' ) ; op = MPI_SUM
      case ( 'and' ) ; op = MPI_LAND
      case ( 'or'  ) ; op = MPI_LOR
      case default
        write (csol,'("unsupported oper `",a,"`")') trim(oper); call csoErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0

  end subroutine CSO_Comm_GetOper

#endif


  ! ********************************************************************
  ! ***
  ! *** send from source, receive at dest
  ! ***
  ! ********************************************************************


  subroutine CSO_Comm_SendRecv_i4( self, value, source, dest, tag, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_Send
    use MPI_F08, only : MPI_Recv
!    use MPI_F08, only : MPI_SendRecv
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(inout)              ::  value
    integer, intent(in)                 ::  source
    integer, intent(in)                 ::  dest
    integer, intent(in)                 ::  tag
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_SendRecv_i4'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Status)      ::  mpi_stat
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! only needed if dest is different from source ...
    if ( dest /= source ) then
      ! data type:
      call self%GetDataType( 'integer', kind(value), dtype, status )
      IF_NOT_OK_RETURN(status=1)
      ! sending or receiving pe?
      if ( self%id == source ) then

        !! testing ...
        !write (csol,'("  x    send ",i0," integers (tag ",i0,")")') 1, tag; call csoPr

        ! send from this (source) pe to destination pe:
        call MPI_Send( value, 1, dtype, &
                         dest, tag, self%comm, &
                         ierror=status )
        IF_MPI_NOT_OK_RETURN(status=1)

      else if ( self%id == dest ) then

        !! testing ...
        !write (csol,'("  x    recv ",i0," integers (tag ",i0,")")') 1, tag; call csoPr

        ! receive on this (destination) pe from source pe:
        call MPI_Recv( value, 1, dtype, &
                         source, tag, self%comm, &
                         mpi_stat, ierror=status )
        IF_MPI_NOT_OK_RETURN(status=1)

      end if  ! source or dest
    end if ! dest differs from source

    ! NEW, but blocks ...
    !! send and receive:
    !call MPI_SendRecv( value, 1, dtype, dest  , tag, &
    !                   value, 1, dtype, source, tag, &
    !                   self%comm, mpi_stat, ierror=status )
    !IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_SendRecv_i4
  
  
  ! *


  subroutine CSO_Comm_SendRecv_i4_1d( self, values, source, dest, tag, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_Send
    use MPI_F08, only : MPI_Recv
    !use MPI_F08, only : MPI_SendRecv
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(inout)              ::  values(:)
    integer, intent(in)                 ::  source
    integer, intent(in)                 ::  dest
    integer, intent(in)                 ::  tag
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_SendRecv_i4_1d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Status)      ::  mpi_stat
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! only needed if dest is different from source ...
    if ( dest /= source ) then
      ! data type:
      call self%GetDataType( 'integer', kind(values), dtype, status )
      IF_NOT_OK_RETURN(status=1)
      ! sending or receiving pe?
      if ( self%id == source ) then

        !! testing ...
        !write (csol,'("  x    send ",i0," integers (tag ",i0,")")') size(values), tag; call csoPr

        ! send from this (source) pe to destination pe:
        call MPI_Send( values, size(values), dtype, &
                         dest, tag, self%comm, &
                         ierror=status )
        IF_MPI_NOT_OK_RETURN(status=1)

      else if ( self%id == dest ) then

        !! testing ...
        !write (csol,'("  x    recv ",i0," integers (tag ",i0,")")') size(values), tag; call csoPr

        ! receive on this (destination) pe from source pe:
        call MPI_Recv( values, size(values), dtype, &
                         source, tag, self%comm, &
                         mpi_stat, ierror=status )
        IF_MPI_NOT_OK_RETURN(status=1)

      end if  ! source or dest
    end if ! dest differs from source

    ! NEW, but blocks ...
    !! send and receive:
    !call MPI_SendRecv( values, size(values), dtype, dest  , tag, &
    !                   values, size(values), dtype, source, tag, &
    !                   self%comm, mpi_stat, ierror=status )
    !IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_SendRecv_i4_1d
  
  
  ! *


  subroutine CSO_Comm_SendRecv_i4_2d( self, values, source, dest, tag, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_Send
    use MPI_F08, only : MPI_Recv
    !use MPI_F08, only : MPI_SendRecv
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(inout)              ::  values(:,:)
    integer, intent(in)                 ::  source
    integer, intent(in)                 ::  dest
    integer, intent(in)                 ::  tag
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_SendRecv_i4_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Status)      ::  mpi_stat
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! only needed if dest is different from source ...
    if ( dest /= source ) then
      ! data type:
      call self%GetDataType( 'integer', kind(values), dtype, status )
      IF_NOT_OK_RETURN(status=1)
      ! sending or receiving pe?
      if ( self%id == source ) then

        !! testing ...
        !write (csol,'("  x    send ",i0," integers (tag ",i0,")")') size(values), tag; call csoPr

        ! send from this (source) pe to destination pe:
        call MPI_Send( values, size(values), dtype, &
                         dest, tag, self%comm, &
                         ierror=status )
        IF_MPI_NOT_OK_RETURN(status=1)

      else if ( self%id == dest ) then

        !! testing ...
        !write (csol,'("  x    recv ",i0," integers (tag ",i0,")")') size(values), tag; call csoPr

        ! receive on this (destination) pe from source pe:
        call MPI_Recv( values, size(values), dtype, &
                         source, tag, self%comm, &
                         mpi_stat, ierror=status )
        IF_MPI_NOT_OK_RETURN(status=1)

      end if  ! source or dest
    end if ! dest differs from source

    ! NEW, but blocks ...
    !! send and receive:
    !call MPI_SendRecv( values, size(values), dtype, dest  , tag, &
    !                   values, size(values), dtype, source, tag, &
    !                   self%comm, mpi_stat, ierror=status )
    !IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_SendRecv_i4_2d
  
  
  ! *


  subroutine CSO_Comm_SendRecv_r4_2d( self, values, source, dest, tag, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_Send
    use MPI_F08, only : MPI_Recv
    !use MPI_F08, only : MPI_SendRecv
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    real(4), intent(inout)              ::  values(:,:)
    integer, intent(in)                 ::  source
    integer, intent(in)                 ::  dest
    integer, intent(in)                 ::  tag
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_SendRecv_r4_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Status)      ::  mpi_stat
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! only needed if dest is different from source ...
    if ( dest /= source ) then
      ! data type:
      call self%GetDataType( 'real', kind(values), dtype, status )
      IF_NOT_OK_RETURN(status=1)
      ! sending or receiving pe?
      if ( self%id == source ) then

        !! testing ...
        !write (csol,'("  x    send ",i0," reals (tag ",i0,")")') size(values), tag; call csoPr

        ! send from this (source) pe to destination pe:
        call MPI_Send( values, size(values), dtype, &
                         dest, tag, self%comm, &
                         status )
        IF_MPI_NOT_OK_RETURN(status=1)

      else if ( self%id == dest ) then

        !! testing ...
        !write (csol,'("  x    recv ",i0," reals (tag ",i0,")")') size(values), tag; call csoPr

        ! receive on this (destination) pe from source pe:
        call MPI_Recv( values, size(values), dtype, &
                         source, tag, self%comm, &
                         mpi_stat, ierror=status )
        IF_MPI_NOT_OK_RETURN(status=1)

      end if  ! source or dest
    end if ! dest differs from source

    ! NEW, but blocks ...
    !! send and receive:
    !call MPI_SendRecv( values, size(values), dtype, dest  , tag, &
    !                   values, size(values), dtype, source, tag, &
    !                   self%comm, mpi_stat, ierror=status )
    !IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_SendRecv_r4_2d
  
  
  ! *


  subroutine CSO_Comm_SendRecv_r8_2d( self, values, source, dest, tag, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Status
    use MPI_F08, only : MPI_Send
    use MPI_F08, only : MPI_Recv
    !use MPI_F08, only : MPI_SendRecv
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    real(8), intent(inout)              ::  values(:,:)
    integer, intent(in)                 ::  source
    integer, intent(in)                 ::  dest
    integer, intent(in)                 ::  tag
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_SendRecv_r8_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Status)      ::  mpi_stat
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! only needed if dest is different from source ...
    if ( dest /= source ) then
      ! data type:
      call self%GetDataType( 'real', kind(values), dtype, status )
      IF_NOT_OK_RETURN(status=1)
      ! sending or receiving pe?
      if ( self%id == source ) then

        !! testing ...
        !write (csol,'("  x    send ",i0," doubles (tag ",i0,")")') size(values), tag; call csoPr

        ! send from this (source) pe to destination pe:
        call MPI_Send( values, size(values), dtype, &
                         dest, tag, self%comm, &
                         status )
        IF_MPI_NOT_OK_RETURN(status=1)

      else if ( self%id == dest ) then

        !! testing ...
        !write (csol,'("  x    recv ",i0," doubles (tag ",i0,")")') size(values), tag; call csoPr

        ! receive on this (destination) pe from source pe:
        call MPI_Recv( values, size(values), dtype, &
                         source, tag, self%comm, &
                         mpi_stat, ierror=status )
        IF_MPI_NOT_OK_RETURN(status=1)

      end if  ! source or dest
    end if ! dest differs from source

    ! NEW, but blocks ...
    !! send and receive:
    !call MPI_SendRecv( values, size(values), dtype, dest  , tag, &
    !                   values, size(values), dtype, source, tag, &
    !                   self%comm, mpi_stat, ierror=status )
    !IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_SendRecv_r8_2d


  ! ********************************************************************
  ! ***
  ! *** reduce
  ! ***
  ! ********************************************************************


  !
  ! reduce values to root;
  ! 'oper' is one of: 'sum' | 'min' | 'max'
  !
  
  subroutine CSO_Comm_Reduce_i( self, oper, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_Reduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    character(len=*), intent(in)        ::  oper
    integer, intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_Reduce_i'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Op)          ::  op
    integer               ::  value_loc
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! operator:
    call self%GetOper( oper, op, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    value_loc = value
    ! reduce over all processors:
    call MPI_Reduce( value_loc, value, 1, dtype, op, self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_Reduce_i


  ! ********************************************************************
  ! ***
  ! *** reduce, and send result to all
  ! ***
  ! ********************************************************************


  subroutine CSO_Comm_AllReduce_i( self, oper, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    character(len=*), intent(in)        ::  oper
    integer, intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_i'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Op)          ::  op
    integer               ::  value_loc
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! operator:
    call self%GetOper( oper, op, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    value_loc = value
    ! reduce over all processors:
    call MPI_AllReduce( value_loc, value, 1, dtype, op, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_i

  ! *

  subroutine CSO_Comm_AllReduce_r( self, oper, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    character(len=*), intent(in)        ::  oper
    real, intent(inout)                 ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_r'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Op)          ::  op
    real                  ::  value_loc
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! operator:
    call self%GetOper( oper, op, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    value_loc = value
    ! reduce over all processors:
    call MPI_AllReduce( value_loc, value, 1, dtype, op, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_r

  ! *

  subroutine CSO_Comm_AllReduce_l( self, oper, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Op
    use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    character(len=*), intent(in)        ::  oper
    logical, intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_l'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
    type(MPI_Op)          ::  op
    logical               ::  value_loc
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'logical', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! operator:
    call self%GetOper( oper, op, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    value_loc = value
    ! reduce over all processors:
    call MPI_AllReduce( value_loc, value, 1, dtype, op, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_l


  ! .......................................................................
  ! ...
  ! ... allreduce (from lekf v3.0.000, should become similar as above ...
  ! ...
  ! .......................................................................


  subroutine CSO_Comm_AllReduce_r4_1d( self, input, output, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    real(4), intent(in)                 ::  input(:)
    real(4), intent(out)                ::  output(:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_r4_1d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(input), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! sum inputs from all pe, store in output on each pe:
    call MPI_AllReduce( input, output, size(input), dtype, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

#else

    ! single pe, just copy:
    output = input

#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_r4_1d


  ! ***


  subroutine CSO_Comm_AllReduce_r4_2d( self, input, output, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    real(4), intent(in)                 ::  input(:,:)
    real(4), intent(out)                ::  output(:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_r4_2d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(input), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! sum inputs from all pe, store in output on each pe:
    call MPI_AllReduce( input, output, size(input), dtype, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

#else

    ! single pe, just copy:
    output = input

#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_r4_2d


  ! ***


  subroutine CSO_Comm_AllReduce_r4_3d( self, input, output, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    real(4), intent(in)                 ::  input(:,:,:)
    real(4), intent(out)                ::  output(:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_r4_3d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(input), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! sum inputs from all pe, store in output on each pe:
    call MPI_AllReduce( input, output, size(input), dtype, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

#else

    ! single pe, just copy:
    output = input

#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_r4_3d


  ! ***


  subroutine CSO_Comm_AllReduce_r4_4d( self, input, output, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    real(4), intent(in)                 ::  input(:,:,:,:)
    real(4), intent(out)                ::  output(:,:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_r4_4d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(input), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! sum inputs from all pe, store in output on each pe:
    call MPI_AllReduce( input, output, size(input), dtype, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

#else

    ! single pe, just copy:
    output = input

#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_r4_4d


  ! ***


  subroutine CSO_Comm_AllReduce_r4_5d( self, input, output, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_SUM
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    real(4), intent(in)                 ::  input(:,:,:,:,:)
    real(4), intent(out)                ::  output(:,:,:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_r4_5d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(input), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! sum inputs from all pe, store in output on each pe:
    call MPI_AllReduce( input, output, size(input), dtype, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

#else

    ! single pe, just copy:
    output = input

#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_r4_5d


  ! ***


  !
  ! Reduce "in place" (input and output buffer are the same)
  !

  subroutine CSO_Comm_AllReduce_InPlace_r4_5d( self, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_SUM
    use MPI_F08, only : MPI_IN_PLACE
    !use MPI_F08, only : MPI_AllReduce
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(inout)    ::  self
    real(4), intent(inout)              ::  values(:,:,:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllReduce_InPlace_r4_5d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! sum output from all pe, return in output on each pe;
    ! use special parameter 'MPI_IN_PLACE' for the send buffer,
    ! receive buffer is the input/output array:
    call MPI_AllReduce( MPI_IN_PLACE, values, size(values), dtype, &
                           MPI_SUM, self%comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllReduce_InPlace_r4_5d


  ! ********************************************************************
  ! ***
  ! *** broadcast
  ! ***
  ! ********************************************************************


  subroutine CSO_Comm_BCast_i( self, rootid, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    integer, intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_i'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( value, 1, dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_i


  ! ***


  subroutine CSO_Comm_BCast_i1_1d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    integer(1), intent(inout)           ::  values(:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_i1_1d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_i1_1d


  ! ***


  subroutine CSO_Comm_BCast_i_1d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    integer, intent(inout)              ::  values(:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_i_1d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_i_1d


  ! ***


  subroutine CSO_Comm_BCast_r4( self, rootid, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r4'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( value, 1, dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r4


  ! ***


  subroutine CSO_Comm_BCast_r4_1d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r4_1d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r4_1d


  ! ***


  subroutine CSO_Comm_BCast_r4_2d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)          ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r4_2d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r4_2d


  ! ***


  subroutine CSO_Comm_BCast_r4_3d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r4_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r4_3d


  ! ***


  subroutine CSO_Comm_BCast_r4_4d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    real(4), intent(inout)              ::  values(:,:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r4_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r4_4d


  ! ***


  subroutine CSO_Comm_BCast_r8( self, rootid, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    real(8), intent(inout)              ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r8'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( value, 1, dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r8


  ! ***


  subroutine CSO_Comm_BCast_r8_1d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    real(8), intent(inout)              ::  values(:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r8_1d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r8_1d


  ! ***


  subroutine CSO_Comm_BCast_r8_2d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)          ::  self
    integer, intent(in)                 ::  rootid
    real(8), intent(inout)              ::  values(:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r8_2d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r8_2d


  ! ***


  subroutine CSO_Comm_BCast_r8_3d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    real(8), intent(inout)              ::  values(:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r8_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r8_3d


  ! ***


  subroutine CSO_Comm_BCast_r8_4d( self, rootid, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    real(8), intent(inout)              ::  values(:,:,:,:)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_r8_2d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( values, size(values), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_r8_4d


  ! ***


  subroutine CSO_Comm_BCast_c( self, rootid, value, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_BCast
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  rootid
    character(len=*), intent(inout)     ::  value
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_BCast_c'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'char', 1, dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! send values from root to all other pe's:
    call MPI_BCast( value, len(value), dtype, rootid, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_BCast_c


  ! ********************************************************************
  ! ***
  ! *** gather
  ! ***
  ! ********************************************************************

  !
  ! Fill on root an array with length npes
  ! with values collected from everywhere.
  !

  subroutine CSO_Comm_Gather_i( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Gather
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)     ::  self
    integer, intent(in)               ::  value
    integer, intent(out)              ::  values(:)  ! (1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_Gather_i'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! collect from all pe and broadcast result:
    call MPI_Gather( value, 1, dtype, &
                     values, 1, dtype, &
                     self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values = value
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_Gather_i

  ! *

  subroutine CSO_Comm_Gather_i_1d( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_Gather
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)     ::  self
    integer, intent(in)               ::  value(:)     ! (n)
    integer, intent(out)              ::  values(:,:)  ! (n,1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_Gather_i_1d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(values), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! collect from all pe and broadcast result:
    call MPI_Gather( value, size(value), dtype, &
                     values, size(value), dtype, &
                     self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values(:,1) = value
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_Gather_i_1d


  ! ********************************************************************
  ! ***
  ! *** gather and broadcast
  ! ***
  ! ********************************************************************

  !
  ! Fill on each processor an array with length npes
  ! with values broadcasted from everywhere.
  !

  subroutine CSO_Comm_AllGather_i( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllGather
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)     ::  self
    integer, intent(in)               ::  value
    integer, intent(out)              ::  values(:)  ! (1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllGather_i'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

    ! check ...
    if ( size(values) /= self%npes ) then
      write (csol,'("output array has size ",i0," while npes=",i0)') size(values), self%npes; call csoErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! collect from all pe and broadcast result:
    call MPI_AllGather( value, 1, dtype, &
                        values, 1, dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values = value
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllGather_i

  ! *

  subroutine CSO_Comm_AllGather_r( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllGather
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)     ::  self
    real(4), intent(in)               ::  value
    real(4), intent(out)              ::  values(:)  ! (1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllGather_r'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

    ! check ...
    if ( size(values) /= self%npes ) then
      write (csol,'("output array has size ",i0," while npes=",i0)') size(values), self%npes; call csoErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! collect from all pe and broadcast result:
    call MPI_AllGather( value, 1, dtype, &
                        values, 1, dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values = value
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllGather_r

  ! *

  subroutine CSO_Comm_AllGather_i1( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllGather
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)     ::  self
    integer, intent(in)               ::  value(:)     ! (n)
    integer, intent(out)              ::  values(:,:)  ! (n,1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllGather_i1'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

    ! check ...
    if ( (size(values,1) /= size(value)) .or. (size(values,2) /= self%npes) ) then
      write (csol,'("output array has shape (",i0,",",i0,") while input has shape (",i0,") and npes=",i0)') &
          size(value), shape(values), self%npes; call csoErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! collect from all pe and broadcast result:
    call MPI_AllGather( value, size(value), dtype, &
                        values, size(value), dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values(:,1) = value(:)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllGather_i1

  ! *

  subroutine CSO_Comm_AllGather_r1( self, value, values, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllGather
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)     ::  self
    real(4), intent(in)               ::  value(:)     ! (n)
    real(4), intent(out)              ::  values(:,:)  ! (n,1:npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllGather_r1'

    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

    ! check ...
    if ( (size(values,1) /= size(value)) .or. (size(values,2) /= self%npes) ) then
      write (csol,'("output array has shape (",i0,",",i0,") while input has shape (",i0,") and npes=",i0)') &
          size(value), shape(values), self%npes; call csoErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'real', kind(value), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! collect from all pe and broadcast result:
    call MPI_AllGather( value, size(value), dtype, &
                        values, size(value), dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    values(:,1) = value(:)
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllGather_r1



  ! ********************************************************************
  ! ***
  ! *** GatherV
  ! ***
  ! ********************************************************************
  
  !
  ! Collect send buffers in 1D array on root.
  ! If send is supposed to be empty, use optional nloc=0 to specify this.
  !

  subroutine CSO_Comm_GatherV_i1_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer(1), intent(in)              ::  send(:)   ! (max(1,nloc))
    integer(1), intent(out)             ::  recv(:)   ! (sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_i1_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
    integer                 ::  ntot
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send)
    end if
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'integer', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( self%root ) then
      if ( size(recv) /= ntot ) then
        write (csol,'("receive buffer has size ",i0," while ntot is ",i0)') size(recv), ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, n, dtype, &
                      recv, recvcounts, displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_i1_1d
  
  ! *

  subroutine CSO_Comm_GatherV_i_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  send(:)   ! (max(1,nloc))
    integer, intent(out)                ::  recv(:)   ! (sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_i_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
    integer                 ::  ntot
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send)
    end if
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'integer', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( self%root ) then
      if ( size(recv) /= ntot ) then
        write (csol,'("receive buffer has size ",i0," while ntot is ",i0)') size(recv), ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, n, dtype, &
                      recv, recvcounts, displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_i_1d
  
  ! *

  subroutine CSO_Comm_GatherV_i_2d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer, intent(in)                 ::  send(:,:)   ! (m,max(1,nloc))
    integer, intent(out)                ::  recv(:,:)   ! (m,sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_i_2d'

    ! --- local ----------------------------------
    
    integer                 ::  m
    integer                 ::  n
#ifdef _MPI
    integer                 ::  ntot
    type(MPI_DataType)      ::  dtype
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send,2)
    end if
    
    ! first dim:
    m = size(send,1)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'integer', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check receive buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(recv) /= (/m,ntot/) ) ) then
        write (csol,'("receive buffer has shape (",i0,",",i0,") while (m,ntot) is (",i0,",",i0,")")') &
                        shape(recv), m,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, m*n, dtype, &
                      recv, m*recvcounts, m*displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_i_2d
  
  ! *

  subroutine CSO_Comm_GatherV_r4_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 4

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:)   ! (max(1,nloc))
    real(wp), intent(out)               ::  recv(:)   ! (sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_r4_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send)
    end if
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( self%root ) then
      if ( size(recv) /= ntot ) then
        write (csol,'("receive buffer has size ",i0," while ntot is ",i0)') size(recv), ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, n, dtype, &
                      recv, recvcounts, displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_r4_1d
  
  ! *

  subroutine CSO_Comm_GatherV_r8_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 8

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:)   ! (max(1,nloc))
    real(wp), intent(out)               ::  recv(:)   ! (sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_r8_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send)
    end if
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( self%root ) then
      if ( size(recv) /= ntot ) then
        write (csol,'("receive buffer has size ",i0," while ntot is ",i0)') size(recv), ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, n, dtype, &
                      recv, recvcounts, displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_r8_1d
  
  ! *

  subroutine CSO_Comm_GatherV_c1_2d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    class(T_CSO_Comm), intent(in)       ::  self
    character(len=1), intent(in)        ::  send(:,:)   ! (m,max(1,nloc))
    character(len=1), intent(out)       ::  recv(:,:)   ! (m,sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_c1_2d'

    ! --- local ----------------------------------
    
    integer                 ::  m
    integer                 ::  n
#ifdef _MPI
    integer                 ::  ntot
    type(MPI_DataType)      ::  dtype
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send,2)
    end if
    
    ! first dim:
    m = size(send,1)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'char', 1, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check receive buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(recv) /= (/m,ntot/) ) ) then
        write (csol,'("receive buffer has shape (",i0,",",i0,") while (m,ntot) is (",i0,",",i0,")")') &
                        shape(recv), m,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, m*n, dtype, &
                      recv, m*recvcounts, m*displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_c1_2d
  
  ! *

  subroutine CSO_Comm_GatherV_r4_2d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 4

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:)   ! (m,max(1,nloc))
    real(wp), intent(out)               ::  recv(:,:)   ! (m,sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_r4_2d'

    ! --- local ----------------------------------
    
    integer                 ::  m
    integer                 ::  n
#ifdef _MPI
    integer                 ::  ntot
    type(MPI_DataType)      ::  dtype
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send,2)
    end if
    
    ! first dim:
    m = size(send,1)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check receive buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(recv) /= (/m,ntot/) ) ) then
        write (csol,'("receive buffer has shape (",i0,",",i0,") while (m,ntot) is (",i0,",",i0,")")') &
                        shape(recv), m,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, m*n, dtype, &
                      recv, m*recvcounts, m*displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_r4_2d
  
  ! *

  subroutine CSO_Comm_GatherV_r8_2d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 8

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:)   ! (m,max(1,nloc))
    real(wp), intent(out)               ::  recv(:,:)   ! (m,sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_r8_2d'

    ! --- local ----------------------------------
    
    integer                 ::  m
    integer                 ::  n
#ifdef _MPI
    integer                 ::  ntot
    type(MPI_DataType)      ::  dtype
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send,2)
    end if
    
    ! first dim:
    m = size(send,1)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check receive buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(recv) /= (/m,ntot/) ) ) then
        write (csol,'("receive buffer has shape (",i0,",",i0,") while (m,ntot) is (",i0,",",i0,")")') &
                        shape(recv), m,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, m*n, dtype, &
                      recv, m*recvcounts, m*displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_r8_2d
  
  ! *

  subroutine CSO_Comm_GatherV_r4_3d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 4

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:,:)   ! (p,m,max(1,nloc))
    real(wp), intent(out)               ::  recv(:,:,:)   ! (p,m,sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_r4_2d'

    ! --- local ----------------------------------
    
    integer                 ::  p, m
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send,3)
    end if
    
    ! first dims:
    p = size(send,1)
    m = size(send,2)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check receive buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(recv) /= (/p,m,ntot/) ) ) then
        write (csol,'("receive buffer has shape (",i0,2(",",i0),") while (m,ntot) is (",i0,2(",",i0),")")') &
                        shape(recv), p,m,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, p*m*n, dtype, &
                      recv, p*m*recvcounts, p*m*displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(:,:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_r4_3d
  
  ! *

  subroutine CSO_Comm_GatherV_r8_3d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 8

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:,:)   ! (p,m,max(1,nloc))
    real(wp), intent(out)               ::  recv(:,:,:)   ! (p,m,sum nloc)
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_GatherV_r8_2d'

    ! --- local ----------------------------------
    
    integer                 ::  p, m
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  recvcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(send,3)
    end if
    
    ! first dims:
    p = size(send,1)
    m = size(send,2)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check receive buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(recv) /= (/p,m,ntot/) ) ) then
        write (csol,'("receive buffer has shape (",i0,2(",",i0),") while (m,ntot) is (",i0,2(",",i0),")")') &
                        shape(recv), p,m,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_GatherV( send, p*m*n, dtype, &
                      recv, p*m*recvcounts, p*m*displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv = send(:,:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_GatherV_r8_3d
  
  
  !
  ! *** gather 2D distributed array
  !

  subroutine CSO_Comm_Gather2D_r4( self, send, offset, recv, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 4

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:)   ! (m,n)
    integer, intent(in)                 ::  offset(2)
    real(wp), intent(out)               ::  recv(:,:)   ! (mtot,ntot)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_Gather2D_r4'

    ! --- local ----------------------------------
    
#ifdef _MPI
    integer                 ::  m, n
    integer                 ::  mtot, ntot
    type(MPI_DataType)      ::  dtype
    integer, allocatable    ::  offs(:,:), shps(:,:)      ! (2,npes)
    integer, allocatable    ::  recvcounts(:), displs(:)  ! (npes)
    integer                 ::  r, rtot
    real(wp), allocatable   ::  rbuf(:)
    integer                 ::  k
#endif

    ! --- begin ----------------------------------
    
#ifdef _MPI

    ! local size:
    m = size(send,1)
    n = size(send,2)
    
    ! data type:
    call self%GetDataType( 'real', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! local send size:
    r = m*n
    ! collect numbers:
    call self%ParInfo( r, status, ntot=rtot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    ! temporary receive buffer:
    allocate( rbuf(rtot), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! collect values from all pe's on root:
    call MPI_GatherV( send, r, dtype, &
                      rbuf, recvcounts, displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! storage:
    allocate( offs(2,0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( shps(2,0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect:
    call self%Gather( offset, offs, status )
    IF_NOT_OK_RETURN(status=1)
    call self%Gather( (/m,n/), shps, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! total shape:
    mtot = maxval( offs(1,:) + shps(1,:) )
    ntot = maxval( offs(2,:) + shps(2,:) )
    
    ! check receive buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(recv) /= (/mtot,ntot/) ) ) then
        write (csol,'("receive buffer has shape (",i0,",",i0,") while (mtot,ntot) is (",i0,",",i0,")")') &
                        shape(recv), mtot,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
      ! loop over pe's:
      do k = 0, self%npes-1
        ! insert result:
        recv(offs(1,k)+1:offs(1,k)+shps(1,k),offs(2,k)+1:offs(2,k)+shps(2,k)) = &
              reshape( rbuf(displs(k)+1:displs(k)+recvcounts(k)), (/shps(1,k),shps(2,k)/) )
      end do
    end if  ! root

    ! clear:
    deallocate( offs, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( shps, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( rbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! check ...
    if ( any( shape(recv) /= shape(send) ) ) then
      write (csol,'("receive buffer has shape (",i0,",",i0,") while send is (",i0,",",i0,")")') &
                      shape(recv), shape(send); call csoErr
      TRACEBACK; status=1; return
    end if
    ! just copy ...
    recv = send

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_Gather2D_r4
  
  ! *

  subroutine CSO_Comm_Gather2D_r8( self, send, offset, recv, status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_GatherV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 8

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:)   ! (m,n)
    integer, intent(in)                 ::  offset(2)
    real(wp), intent(out)               ::  recv(:,:)   ! (mtot,ntot)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_Gather2D_r8'

    ! --- local ----------------------------------
    
#ifdef _MPI
    integer                 ::  m, n
    integer                 ::  mtot, ntot
    type(MPI_DataType)      ::  dtype
    integer, allocatable    ::  offs(:,:), shps(:,:)      ! (2,npes)
    integer, allocatable    ::  recvcounts(:), displs(:)  ! (npes)
    integer                 ::  r, rtot
    real(wp), allocatable   ::  rbuf(:)
    integer                 ::  k
#endif

    ! --- begin ----------------------------------
    
#ifdef _MPI

    ! local size:
    m = size(send,1)
    n = size(send,2)
    
    ! data type:
    call self%GetDataType( 'real', kind(send), dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( recvcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! local send size:
    r = m*n
    ! collect numbers:
    call self%ParInfo( r, status, ntot=rtot, recvcounts=recvcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    ! temporary receive buffer:
    allocate( rbuf(rtot), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! collect values from all pe's on root:
    call MPI_GatherV( send, r, dtype, &
                      rbuf, recvcounts, displs, dtype, &
                      self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! storage:
    allocate( offs(2,0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( shps(2,0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect:
    call self%Gather( offset, offs, status )
    IF_NOT_OK_RETURN(status=1)
    call self%Gather( (/m,n/), shps, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! total shape:
    mtot = maxval( offs(1,:) + shps(1,:) )
    ntot = maxval( offs(2,:) + shps(2,:) )
    
    ! check receive buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(recv) /= (/mtot,ntot/) ) ) then
        write (csol,'("receive buffer has shape (",i0,",",i0,") while (mtot,ntot) is (",i0,",",i0,")")') &
                        shape(recv), mtot,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
      ! loop over pe's:
      do k = 0, self%npes-1
        ! insert result:
        recv(offs(1,k)+1:offs(1,k)+shps(1,k),offs(2,k)+1:offs(2,k)+shps(2,k)) = &
              reshape( rbuf(displs(k)+1:displs(k)+recvcounts(k)), (/shps(1,k),shps(2,k)/) )
      end do
    end if  ! root

    ! clear:
    deallocate( offs, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( shps, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( rbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! check ...
    if ( any( shape(recv) /= shape(send) ) ) then
      write (csol,'("receive buffer has shape (",i0,",",i0,") while send is (",i0,",",i0,")")') &
                      shape(recv), shape(send); call csoErr
      TRACEBACK; status=1; return
    end if
    ! just copy ...
    recv = send

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_Gather2D_r8



  ! ********************************************************************
  ! ***
  ! *** scatter
  ! ***
  ! ********************************************************************
  
  subroutine CSO_Comm_ScatterV_i1_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_ScatterV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 1

    class(T_CSO_Comm), intent(in)       ::  self
    integer(wp), intent(in)             ::  send(:)   ! (sum nloc)
    integer(wp), intent(out)            ::  recv(:)   ! (max(1,nloc))
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ScatterV_i1_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  sendcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(recv)
    end if
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'integer', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( sendcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, sendcounts=sendcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( self%root ) then
      if ( size(send) /= ntot ) then
        write (csol,'("send buffer has size ",i0," while ntot is ",i0)') size(send), ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_ScatterV( send, sendcounts, displs, dtype, &
                       recv, n                 , dtype, &
                       self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv(1:n) = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_ScatterV_i1_1d
  
  ! *
  
  subroutine CSO_Comm_ScatterV_r4_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_ScatterV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 4

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:)   ! (sum nloc)
    real(wp), intent(out)               ::  recv(:)   ! (max(1,nloc))
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ScatterV_r4_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  sendcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(recv)
    end if
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( sendcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, sendcounts=sendcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( self%root ) then
      if ( size(send) /= ntot ) then
        write (csol,'("send buffer has size ",i0," while ntot is ",i0)') size(send), ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_ScatterV( send, sendcounts, displs, dtype, &
                       recv, n                 , dtype, &
                       self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv(1:n) = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_ScatterV_r4_1d
  
  ! *
  
  subroutine CSO_Comm_ScatterV_c1_2d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_ScatterV
#endif

    ! --- in/out ---------------------------------
    
    class(T_CSO_Comm), intent(in)       ::  self
    character(len=1), intent(in)        ::  send(:,:)   ! (m,sum nloc)
    character(len=1), intent(out)       ::  recv(:,:)   ! (m,max(1,nloc))
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ScatterV_c1_2d'

    ! --- local ----------------------------------
    
    integer                 ::  m1
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  sendcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(recv,2)
    end if
    
    ! first dim:
    m1 = size(recv,1)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'char', 1, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( sendcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, sendcounts=sendcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check send buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(send) /= (/m1,ntot/) ) ) then
        write (csol,'("send buffer has shape (",i0,",",i0,") while (m,ntot) is (",i0,",",i0,")")') &
                        shape(send), m1,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_ScatterV( send, m1*sendcounts, m1*displs, dtype, &
                       recv, m1*n                    , dtype, &
                       self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv(:,1:n) = send(:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_ScatterV_c1_2d
  
  ! *
  
  subroutine CSO_Comm_ScatterV_r4_2d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_ScatterV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 4

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:)   ! (m,sum nloc)
    real(wp), intent(out)               ::  recv(:,:)   ! (m,max(1,nloc))
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ScatterV_r4_2d'

    ! --- local ----------------------------------
    
    integer                 ::  m1
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  sendcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(recv,2)
    end if
    
    ! first dim:
    m1 = size(recv,1)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( sendcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, sendcounts=sendcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check send buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(send) /= (/m1,ntot/) ) ) then
        write (csol,'("send buffer has shape (",i0,",",i0,") while (m,ntot) is (",i0,",",i0,")")') &
                        shape(send), m1,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_ScatterV( send, m1*sendcounts, m1*displs, dtype, &
                       recv, m1*n                    , dtype, &
                       self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv(:,1:n) = send(:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_ScatterV_r4_2d
  
  ! *
  
  subroutine CSO_Comm_ScatterV_r4_3d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_ScatterV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 4

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:,:)   ! (m1,m2,sum nloc)
    real(wp), intent(out)               ::  recv(:,:,:)   ! (m1,m2,max(1,nloc))
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ScatterV_r4_3d'

    ! --- local ----------------------------------
    
    integer                 ::  m1, m2
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  sendcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(recv,2)
    end if
    
    ! first dims:
    m1 = size(recv,1)
    m2 = size(recv,2)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( sendcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, sendcounts=sendcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check send buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(send) /= (/m1,m2,ntot/) ) ) then
        write (csol,'("send buffer has shape (",i0,2(",",i0),") while (m,ntot) is (",i0,2(",",i0),")")') &
                        shape(send), m1,m2,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_ScatterV( send, m1*m2*sendcounts, m1*m2*displs, dtype, &
                       recv, m1*m2*n                       , dtype, &
                       self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv(:,:,1:n) = send(:,:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_ScatterV_r4_3d
  
  
  !
  ! *** r8
  !
  
  
  subroutine CSO_Comm_ScatterV_r8_1d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_ScatterV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 8

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:)   ! (sum nloc)
    real(wp), intent(out)               ::  recv(:)   ! (max(1,nloc))
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ScatterV_r8_1d'

    ! --- local ----------------------------------
    
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  sendcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(recv)
    end if
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( sendcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, sendcounts=sendcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( self%root ) then
      if ( size(send) /= ntot ) then
        write (csol,'("send buffer has size ",i0," while ntot is ",i0)') size(send), ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_ScatterV( send, sendcounts, displs, dtype, &
                       recv, n                 , dtype, &
                       self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv(1:n) = send(1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_ScatterV_r8_1d
  
  ! *
  
  subroutine CSO_Comm_ScatterV_r8_2d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_ScatterV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 8

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:)   ! (m1,sum nloc)
    real(wp), intent(out)               ::  recv(:,:)   ! (m1,max(1,nloc))
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ScatterV_r8_2d'

    ! --- local ----------------------------------
    
    integer                 ::  m1
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  sendcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(recv,2)
    end if
    
    ! first dim:
    m1 = size(recv,1)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( sendcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, sendcounts=sendcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check send buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(send) /= (/m1,ntot/) ) ) then
        write (csol,'("send buffer has shape (",i0,",",i0,") while (m1,ntot) is (",i0,",",i0,")")') &
                        shape(send), m1,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_ScatterV( send, m1*sendcounts, m1*displs, dtype, &
                       recv, m1*n                    , dtype, &
                       self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv(:,1:n) = send(:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_ScatterV_r8_2d
  
  ! *
  
  subroutine CSO_Comm_ScatterV_r8_3d( self, send, recv, status, &
                                      nloc )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_ScatterV
#endif

    ! --- in/out ---------------------------------
    
    integer, parameter                  ::  wp = 8

    class(T_CSO_Comm), intent(in)       ::  self
    real(wp), intent(in)                ::  send(:,:,:)   ! (m1,m2,sum nloc)
    real(wp), intent(out)               ::  recv(:,:,:)   ! (m1,m2,max(1,nloc))
    integer, intent(out)                ::  status
    
    integer, intent(in), optional       ::  nloc

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ScatterV_r8_3d'

    ! --- local ----------------------------------
    
    integer                 ::  m1, m2
    integer                 ::  n
#ifdef _MPI
    type(MPI_DataType)      ::  dtype
    integer                 ::  ntot
    integer, allocatable    ::  sendcounts(:)  ! (npes)
    integer, allocatable    ::  displs(:)  ! (npes)
#endif

    ! --- begin ----------------------------------
    
    ! local size, take from optional argument if present (value is probably zero ..)
    if ( present(nloc) ) then
      n = nloc
    else
      n = size(recv,2)
    end if
    
    ! first dims:
    m1 = size(recv,1)
    m2 = size(recv,2)
    
#ifdef _MPI

    ! data type:
    call self%GetDataType( 'real', wp, dtype, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( sendcounts(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( displs(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! collect numbers:
    call self%ParInfo( n, status, ntot=ntot, sendcounts=sendcounts, displs=displs )
    IF_NOT_OK_RETURN(status=1)
    
    ! check send buffer ...
    if ( self%root ) then
      ! check ...
      if ( any( shape(send) /= (/m1,m2,ntot/) ) ) then
        write (csol,'("send buffer has shape (",i0,2(",",i0),") while (m,ntot) is (",i0,2(",",i0),")")') &
                        shape(send), m1,m2,ntot; call csoErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! collect values from all pe's on root:
    call MPI_ScatterV( send, m1*m2*sendcounts, m1*m2*displs, dtype, &
                       recv, m1*m2*n                       , dtype, &
                       self%root_id, self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( displs, stat=status )
    IF_NOT_OK_RETURN(status=1)

#else

    ! just copy ...
    if ( n > 0 ) recv(:,:,1:n) = send(:,:,1:n)

#endif
    
    ! ok
    status = 0

  end subroutine CSO_Comm_ScatterV_r8_3d



  ! ********************************************************************
  ! ***
  ! *** AllToAll
  ! ***
  ! ********************************************************************


  subroutine CSO_Comm_AllToAll_i4_1d( self, sendbuf, sendcount, &
                                            recvbuf, recvcount, &
                                            status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllToAll
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer(4), intent(in)              ::  sendbuf(:)
    integer, intent(in)                 ::  sendcount
    integer(4), intent(out)             ::  recvbuf(:)
    integer, intent(in)                 ::  recvcount
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllToAll_i4_1d'

    ! --- local ----------------------------------

#ifdef _MPI
    type(MPI_DataType)    ::  dtype
#endif

    ! --- begin ----------------------------------

#ifdef _MPI
    ! data type:
    call self%GetDataType( 'integer', kind(sendbuf), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! swap:
    call MPI_AllToAll( sendbuf, sendcount, dtype, &
                       recvbuf, recvcount, dtype, &
                       self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllToAll_i4_1d



  ! ********************************************************************
  ! ***
  ! *** AllToAllV
  ! ***
  ! ********************************************************************

  subroutine CSO_Comm_AllToAllV_i4_1d( self, sendbuf, sendcounts, &
                                             recvbuf, recvcounts, &
                                             status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    integer(4), intent(in)              ::  sendbuf(:)      ! (:)
    integer, intent(in)                 ::  sendcounts(:)   ! (nproc)
    integer(4), intent(out)             ::  recvbuf(:)      ! (:)
    integer, intent(in)                 ::  recvcounts(:)   ! (nproc)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllToAllV_i4_1d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    integer, allocatable      ::  sdispls(:)
    integer, allocatable      ::  rdispls(:)
    type(MPI_DataType)        ::  dtype
#endif

    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! storage:
    allocate( sdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! compute displacements:
    call csoc%Displacements( sendcounts, sdispls, status )
    IF_NOT_OK_RETURN(status=1)
    call csoc%Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! data type:
    call self%GetDataType( 'integer', kind(sendbuf), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
                        recvbuf, recvcounts, rdispls, dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! storage:
    deallocate( sdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllToAllV_i4_1d
  
  ! *

  subroutine CSO_Comm_AllToAllV_r4_1d( self, sendbuf, sendcounts, &
                                             recvbuf, recvcounts, &
                                             status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    real(4), intent(in)                 ::  sendbuf(:)      ! (:)
    integer, intent(in)                 ::  sendcounts(:)   ! (nproc)
    real(4), intent(out)                ::  recvbuf(:)      ! (:)
    integer, intent(in)                 ::  recvcounts(:)   ! (nproc)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllToAllV_r4_1d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    integer, allocatable      ::  sdispls(:)
    integer, allocatable      ::  rdispls(:)
    type(MPI_DataType)        ::  dtype
#endif

    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! storage:
    allocate( sdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! compute displacements:
    call csoc%Displacements( sendcounts, sdispls, status )
    IF_NOT_OK_RETURN(status=1)
    call csoc%Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! data type:
    call self%GetDataType( 'real', kind(sendbuf), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
                        recvbuf, recvcounts, rdispls, dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! storage:
    deallocate( sdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllToAllV_r4_1d
  
  ! *

  subroutine CSO_Comm_AllToAllV_r8_1d( self, sendbuf, sendcounts, &
                                             recvbuf, recvcounts, &
                                             status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    real(8), intent(in)                 ::  sendbuf(:)      ! (:)
    integer, intent(in)                 ::  sendcounts(:)   ! (nproc)
    real(8), intent(out)                ::  recvbuf(:)      ! (:)
    integer, intent(in)                 ::  recvcounts(:)   ! (nproc)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllToAllV_r8_1d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    integer, allocatable      ::  sdispls(:)
    integer, allocatable      ::  rdispls(:)
    type(MPI_DataType)        ::  dtype
#endif

    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! storage:
    allocate( sdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! compute displacements:
    call csoc%Displacements( sendcounts, sdispls, status )
    IF_NOT_OK_RETURN(status=1)
    call csoc%Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! data type:
    call self%GetDataType( 'real', kind(sendbuf), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
                        recvbuf, recvcounts, rdispls, dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! storage:
    deallocate( sdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllToAllV_r8_1d
  
  ! *

  subroutine CSO_Comm_AllToAllV_r4_2d( self, sendbuf, sendcounts, &
                                             recvbuf, recvcounts, &
                                             status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    real(4), intent(in)                 ::  sendbuf(:,:)    ! (n,:)
    integer, intent(in)                 ::  sendcounts(:)   ! (nproc)
    real(4), intent(out)                ::  recvbuf(:,:)    ! (n,:)
    integer, intent(in)                 ::  recvcounts(:)   ! (nproc)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllToAllV_r4_2d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    integer, allocatable      ::  sdispls(:)
    integer, allocatable      ::  rdispls(:)
    integer                   ::  n
    type(MPI_DataType)        ::  dtype
#endif

    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! storage:
    allocate( sdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! compute displacements:
    call csoc%Displacements( sendcounts, sdispls, status )
    IF_NOT_OK_RETURN(status=1)
    call csoc%Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! size of first dim:
    n = size(sendbuf,1)
    
    ! data type:
    call self%GetDataType( 'real', kind(sendbuf), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! specific call:
    call MPI_AllToAllV( sendbuf, n*sendcounts, n*sdispls, dtype, &
                        recvbuf, n*recvcounts, n*rdispls, dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! storage:
    deallocate( sdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllToAllV_r4_2d
  
  ! *

  subroutine CSO_Comm_AllToAllV_r8_2d( self, sendbuf, sendcounts, &
                                             recvbuf, recvcounts, &
                                             status )

#ifdef _MPI
    use MPI_F08, only : MPI_DataType
    use MPI_F08, only : MPI_AllToAllV
#endif

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)       ::  self
    real(8), intent(in)                 ::  sendbuf(:,:)    ! (n,:)
    integer, intent(in)                 ::  sendcounts(:)   ! (nproc)
    real(8), intent(out)                ::  recvbuf(:,:)    ! (n,:)
    integer, intent(in)                 ::  recvcounts(:)   ! (nproc)
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_AllToAllV_r8_2d'

    ! --- local ----------------------------------
    
#ifdef _MPI
    integer, allocatable      ::  sdispls(:)
    integer, allocatable      ::  rdispls(:)
    integer                   ::  n
    type(MPI_DataType)        ::  dtype
#endif

    ! --- begin ----------------------------------
    
#ifdef _MPI
    ! storage:
    allocate( sdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(csoc%npes), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! compute displacements:
    call csoc%Displacements( sendcounts, sdispls, status )
    IF_NOT_OK_RETURN(status=1)
    call csoc%Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! size of first dim:
    n = size(sendbuf,1)
    
    ! data type:
    call self%GetDataType( 'real', kind(sendbuf), dtype, status )
    IF_NOT_OK_RETURN(status=1)
    ! specific call:
    call MPI_AllToAllV( sendbuf, n*sendcounts, n*sdispls, dtype, &
                        recvbuf, n*recvcounts, n*rdispls, dtype, &
                        self%comm, ierror=status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! storage:
    deallocate( sdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif

    ! ok
    status = 0

  end subroutine CSO_Comm_AllToAllV_r8_2d



  ! ********************************************************************
  ! ***
  ! *** i/o tools
  ! ***
  ! ********************************************************************

  !
  ! Data with different size n per processor.
  ! Returns:
  ! - total size (output dimension) 
  ! - start index in output array (netcdf put_var argument)
  !

  subroutine CSO_Comm_ParInfo( self, n, status, &
                                ntot, istart, recvcounts, sendcounts, displs, displ )

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)     ::  self
    integer, intent(in)               ::  n
    integer, intent(out)              ::  status
    
    integer, intent(out), optional    ::  ntot
    integer, intent(out), optional    ::  istart
    integer, intent(out), optional    ::  recvcounts(:)  ! (npes)
    integer, intent(out), optional    ::  sendcounts(:)  ! (npes)
    integer, intent(out), optional    ::  displs(:)  ! (npes)
    integer, intent(out), optional    ::  displ

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ParInfo'

    ! --- local ----------------------------------
    
    integer, allocatable    ::  ns(:)
    integer                 ::  i

    ! --- begin ----------------------------------
    
    ! storage:
    allocate( ns(0:self%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of data values for each processor, broadcast to all:
    call self%AllGather( n, ns, status )
    IF_NOT_OK_RETURN(status=1)

    ! total number of data values:
    if ( present(ntot) ) ntot = sum(ns)

    ! start in output array:
    if ( present(istart) ) then
      if ( self%id == 0 ) then
        istart = 1
      else
        istart = sum(ns(0:self%id-1)) + 1
      end if      
    end if
    
    ! receive and send counts:
    if ( present(recvcounts) ) recvcounts = ns
    if ( present(sendcounts) ) sendcounts = ns

    ! displacements:
    if ( present(displs) ) then
      do i = 1, self%npes
        if ( i == 1 ) then
          displs(i) = 0
        else
          displs(i) = sum(ns(0:i-2))
        end if
      end do
    end if
    ! current displacement:
    if ( present(displ) ) then
      if ( self%id == 0 ) then
        displ = 0
      else
        displ = sum(ns(0:self%id-1))
      end if
    end if

    ! clear:
    deallocate( ns, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine CSO_Comm_ParInfo
  
  
  ! *
  
  
  ! compute displacements from counts
  
  subroutine CSO_Comm_Displacements( self, counts, displs, status )

    ! --- in/out ---------------------------------

    class(T_CSO_Comm), intent(in)     ::  self
    integer, intent(in)               ::  counts(:)   ! (npes)
    integer, intent(out)              ::  displs(:)   ! (npes)
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Comm_ParInfo'

    ! --- local ----------------------------------
    
    integer         ::  k

    ! --- begin ----------------------------------

    ! check ...
    if ( size(counts) /= self%npes ) then
      write (csol,'("argument `counts` has size ",i0," while npes=",i0)') size(counts), self%npes; call csoErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( size(displs) /= self%npes ) then
      write (csol,'("argument `displs` has size ",i0," while npes=",i0)') size(displs), self%npes; call csoErr
      TRACEBACK; status=1; return
    end if

    ! fill:
    displs(1) = 0
    do k = 2, size(displs)
      displs(k) = displs(k-1) + counts(k-1)
    end do

    ! ok
    status = 0

  end subroutine CSO_Comm_Displacements


end module CSO_Comm

