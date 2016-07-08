!#######################################################################
!
! MPIF90  -  generic interfaces around MPI routines
!
! Supported data types (if implemented):
!   integer(4)
!   real(4)
!   real(8)
!
! Supported buffer shapes (if implemented):
!   ''          :  0D (scalar)
!   (:)         :  1D
!   (:,:)       :  2D
!   (:,:,:)     :  3D
!   (:,:,:,:)   :  4D
!
! Broadcast data from root to all other process.
!
!   subroutine MPIF90_BCast( buffer, root, comm, status )
!     <dtype>, intent(inout)        ::  buffer<shp>
!     integer, intent(in)           ::  root
!     integer, intent(in)           ::  comm
!     integer, intent(out)          ::  status
!   end subroutine MPIF90_BCast
!
! Send same number of elements from each process to all others,
! and receive elements from all processes in return:
!
!   subroutine MPIF90_AllToAll( sendbuf, sendcount, &
!                               recvbuf, recvcount, &
!                               comm, status )
!     <dtype>, intent(in)           ::  sendbuf<shp>
!     integer, intent(in)           ::  sendcount
!     <dtype>, intent(out)          ::  recvbuf<shp>
!     integer, intent(in)           ::  recvcount
!     integer, intent(in)           ::  comm
!     integer, intent(out)          ::  status
!   end subroutine MPIF90_AllToAll
!
! Send varying number of elements from each process to all others,
! and receive elements from all processes in return:
!
!   subroutine MPIF90_AllToAllV( sendbuf, sendcounts, sdispls, &
!                                recvbuf, recvcounts, rdispls, &
!                                comm, status )
!     <dtype>, intent(in)           ::  sendbuf<shp>
!     integer, intent(in)           ::  sendcounts(<nproc>)
!     <dtype>, intent(out)          ::  recvbuf<shp>
!     integer, intent(in)           ::  recvcounts(<nproc>)
!     integer, intent(in)           ::  comm
!     integer, intent(out)          ::  status
!   end subroutine MPIF90_AllToAllV
!
! Compute sum or other 'op' over all processes,
! broadcast to the result to all:
!
!   subroutine MPIF90_AllReduce( sendbuf, recvbuf, op, comm, status )
!     <dtype>, intent(inout)        ::  sendbuf<shp>
!     <dtype>, intent(inout)        ::  recvbuf<shp>
!     integer, intent(in)           ::  op
!     integer, intent(in)           ::  comm
!     integer, intent(out)          ::  status
!   end subroutine MPIF90_AllReduce
!
!
!#######################################################################
!
#define TRACEBACK write (*,'("in ",a," (line",i5,")")') __FILE__, __LINE__
#define IF_MPI_NOT_OK_RETURN(action) if (status/=MPI_SUCCESS) then; call MPI_Error_String(status,msg,lenmsg,status); write (*,'(a)') msg(1:lenmsg); TRACEBACK; action; return; end if
!
!#######################################################################

module MPIF90

#ifdef _MPI  
  use MPI, only : MPI_SUCCESS
  use MPI, only : MPI_Error_String
  use MPI, only : MPI_COMM_WORLD
  use MPI, only : MPI_SUM, MPI_MIN, MPI_MAX, MPI_LAND, MPI_LOR
#endif

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public    ::  MPI_COMM_WORLD
  public    ::  MPI_SUCCESS
  public    ::  MPI_SUM, MPI_MIN, MPI_MAX, MPI_LAND, MPI_LOR

  public    ::  MPIF90_GetType
  public    ::  MPIF90_Displacements
  public    ::  MPIF90_BCast
  public    ::  MPIF90_Gather
  public    ::  MPIF90_AllGather
  public    ::  MPIF90_GatherV
  public    ::  MPIF90_AllToAll
  public    ::  MPIF90_AllToAllV
  public    ::  MPIF90_AllReduce
  
  
  ! --- const ------------------------------------
  
#ifdef _MPI
  ! parameters used from MPI module
#else
  ! dummy values for non-mpi code:
  integer, parameter    ::  MPI_COMM_WORLD = -1
  integer, parameter    ::  MPI_SUM        = -1
  integer, parameter    ::  MPI_MIN        = -1
  integer, parameter    ::  MPI_MAx        = -1
  integer, parameter    ::  MPI_LAND       = -1
  integer, parameter    ::  MPI_LOR        = -1
#endif
  
  
  ! --- interfaces -------------------------------
  
  interface MPIF90_GetType
    module procedure MPIF90_GetType_i4
    module procedure MPIF90_GetType_r4
    module procedure MPIF90_GetType_r8
  end interface MPIF90_GetType
  
  interface MPIF90_BCast
    module procedure MPIF90_BCast_i4_0d
    module procedure MPIF90_BCast_i4_1d
    module procedure MPIF90_BCast_i4_4d
    module procedure MPIF90_BCast_r4_0d
    module procedure MPIF90_BCast_r4_1d
    module procedure MPIF90_BCast_r4_4d
    module procedure MPIF90_BCast_r8_0d
    module procedure MPIF90_BCast_r8_1d
    module procedure MPIF90_BCast_r8_4d
  end interface MPIF90_BCast
  
  interface MPIF90_Gather
    module procedure MPIF90_Gather_i4_1d
    module procedure MPIF90_Gather_r4_1d
    module procedure MPIF90_Gather_r8_1d
  end interface MPIF90_Gather
  
  interface MPIF90_AllGather
    module procedure MPIF90_AllGather_i4_1d
  end interface MPIF90_AllGather
  
  interface MPIF90_GatherV
    module procedure MPIF90_GatherV_s1_1d
    module procedure MPIF90_GatherV_i4_1d
    module procedure MPIF90_GatherV_r4_1d
    module procedure MPIF90_GatherV_r8_1d
  end interface MPIF90_GatherV
  
  interface MPIF90_AllToAll
    module procedure MPIF90_AllToAll_i4_1d
    module procedure MPIF90_AllToAll_r4_1d
    module procedure MPIF90_AllToAll_r8_1d
  end interface MPIF90_AllToAll
  
  interface MPIF90_AllToAllV
    module procedure MPIF90_AllToAllV_i4_2d
    module procedure MPIF90_AllToAllV_i4_3d
    module procedure MPIF90_AllToAllV_r4_2d
    module procedure MPIF90_AllToAllV_r4_3d
    module procedure MPIF90_AllToAllV_r8_2d
    module procedure MPIF90_AllToAllV_r8_3d
  end interface MPIF90_AllToAllV
  
  interface MPIF90_AllReduce
    module procedure MPIF90_AllReduce_l_0d
    module procedure MPIF90_AllReduce_i4_0d
    module procedure MPIF90_AllReduce_r4_0d
    module procedure MPIF90_AllReduce_r8_0d
  end interface MPIF90_AllReduce
  
  
  ! --- local -------------------------------
  
  ! used for error messages:
  character(len=1024)     ::  msg
  integer                 ::  lenmsg
  
  
contains


  ! ********************************************************************
  ! mpi types
  ! ********************************************************************
  
  subroutine MPIF90_GetType_i4( value, mpi_type, status )
    ! external:
    use MPI, only : MPI_INTEGER
    ! arguments:
    integer(4), intent(in)      ::  value
    integer, intent(out)        ::  mpi_type
    integer, intent(out)        ::  status
    ! fill:
    mpi_type = MPI_INTEGER
    ! ok
    status = 0
  end subroutine MPIF90_GetType_i4
  
  ! *
  
  subroutine MPIF90_GetType_r4( value, mpi_type, status )
    ! external:
    use MPI, only : MPI_REAL
    ! arguments:
    real(4), intent(in)         ::  value
    integer, intent(out)        ::  mpi_type
    integer, intent(out)        ::  status
    ! fill:
    mpi_type = MPI_REAL
    ! ok
    status = 0
  end subroutine MPIF90_GetType_r4
  
  ! *
  
  subroutine MPIF90_GetType_r8( value, mpi_type, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    ! arguments:
    real(8), intent(in)         ::  value
    integer, intent(out)        ::  mpi_type
    integer, intent(out)        ::  status
    ! fill:
    mpi_type = MPI_DOUBLE_PRECISION
    ! ok
    status = 0
  end subroutine MPIF90_GetType_r8


  ! ********************************************************************
  ! displacements
  ! ********************************************************************
  
  ! compute displacements from counts
  
  subroutine MPIF90_Displacements( counts, displs, status )
    ! arguments:
    integer, intent(in)     ::  counts(:)
    integer, intent(out)    ::  displs(:)
    integer, intent(out)    ::  status
    ! local:
    integer         ::  k
    ! check ...
    if ( size(displs) /= size(counts) ) then
      write (*,'("size of displs (",i0,") does not match with size of counts (",i0,")")') &
                    size(displs), size(counts)
      write (*,'("ERROR - in MPIF90_Displacements")'); status=1; return
    end if
    ! fill:
    displs(1) = 0
    do k = 2, size(displs)
      displs(k) = displs(k-1) + counts(k-1)
    end do
    ! ok
    status = 0
  end subroutine MPIF90_Displacements


  ! ********************************************************************
  ! BCast
  ! ********************************************************************

  subroutine MPIF90_BCast_i4_0d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_BCast
    ! arguments:
    integer(4), intent(inout)     ::  buffer
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, 1, MPI_INTEGER, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_i4_0d
  
  ! *

  subroutine MPIF90_BCast_i4_1d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_BCast
    ! arguments:
    integer(4), intent(inout)     ::  buffer(:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, size(buffer), MPI_INTEGER, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_i4_1d
  
  ! *

  subroutine MPIF90_BCast_i4_4d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_BCast
    ! arguments:
    integer(4), intent(inout)     ::  buffer(:,:,:,:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, size(buffer), MPI_INTEGER, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_i4_4d
  
  ! *

  subroutine MPIF90_BCast_r4_0d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_BCast
    ! arguments:
    real(4), intent(inout)        ::  buffer
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, 1, MPI_REAL, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_r4_0d
  
  ! *

  subroutine MPIF90_BCast_r4_1d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_BCast
    ! arguments:
    real(4), intent(inout)        ::  buffer(:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, size(buffer), MPI_REAL, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_r4_1d
  
  ! *

  subroutine MPIF90_BCast_r4_4d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_BCast
    ! arguments:
    real(4), intent(inout)        ::  buffer(:,:,:,:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, size(buffer), MPI_REAL, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_r4_4d
  
  ! *

  subroutine MPIF90_BCast_r8_0d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_BCast
    ! arguments:
    real(8), intent(inout)        ::  buffer
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, 1, MPI_DOUBLE_PRECISION, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_r8_0d
  
  ! *

  subroutine MPIF90_BCast_r8_1d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_BCast
    ! arguments:
    real(8), intent(inout)        ::  buffer(:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, size(buffer), MPI_DOUBLE_PRECISION, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_r8_1d
  
  ! *

  subroutine MPIF90_BCast_r8_4d( buffer, root, comm, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_BCast
    ! arguments:
    real(8), intent(inout)        ::  buffer(:,:,:,:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_BCast( buffer, size(buffer), MPI_DOUBLE_PRECISION, root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_BCast_r8_4d


  ! ********************************************************************
  ! Gather
  ! ********************************************************************

  subroutine MPIF90_Gather_i4_1d( sendbuf, sendcount, &
                                  recvbuf, recvcount, &
                                  root, comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_Gather
    ! arguments:
    integer(4), intent(in)        ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    integer(4), intent(out)       ::  recvbuf(:)
    integer, intent(in)           ::  recvcount
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_Gather( sendbuf, sendcount, MPI_INTEGER, &
                     recvbuf, recvcount, MPI_INTEGER, &
                     root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_Gather_i4_1d
  
  ! *

  subroutine MPIF90_Gather_r4_1d( sendbuf, sendcount, &
                                  recvbuf, recvcount, &
                                  root, comm, status )
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_Gather
    ! arguments:
    real(4), intent(in)           ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    real(4), intent(out)          ::  recvbuf(:)
    integer, intent(in)           ::  recvcount
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_Gather( sendbuf, sendcount, MPI_REAL, &
                     recvbuf, recvcount, MPI_REAL, &
                     root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_Gather_r4_1d
  
  ! *

  subroutine MPIF90_Gather_r8_1d( sendbuf, sendcount, &
                                  recvbuf, recvcount, &
                                  root, comm, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_Gather
    ! arguments:
    real(8), intent(in)           ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    real(8), intent(out)          ::  recvbuf(:)
    integer, intent(in)           ::  recvcount
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_Gather( sendbuf, sendcount, MPI_DOUBLE_PRECISION, &
                     recvbuf, recvcount, MPI_DOUBLE_PRECISION, &
                     root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_Gather_r8_1d


  ! ********************************************************************
  ! AllGather
  ! ********************************************************************
  
  subroutine MPIF90_AllGather_i4_0d( sendbuf, recvbuf, comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_AllGather
    ! arguments:
    integer(4), intent(in)        ::  sendbuf       ! (1)
    integer(4), intent(out)       ::  recvbuf(:)    ! (1*nproc)
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllGather( sendbuf, 1, MPI_INTEGER, &
                        recvbuf, 1, MPI_INTEGER, &
                        comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllGather_i4_0d
  
  ! *
  
  subroutine MPIF90_AllGather_i4_1d( sendbuf, recvbuf, comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_AllGather
    ! arguments:
    integer(4), intent(in)        ::  sendbuf(:)   ! (n)
    integer(4), intent(out)       ::  recvbuf(:)   ! (n*nproc)
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! local:
    integer                       ::  recvcount
    !  number of values received from any process:
    recvcount = size(sendbuf)
    ! specific call:
    call MPI_AllGather( sendbuf, size(sendbuf), MPI_INTEGER, &
                        recvbuf, recvcount    , MPI_INTEGER, &
                        comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllGather_i4_1d


  ! ********************************************************************
  ! GatherV
  ! ********************************************************************
  
  subroutine MPIF90_GatherV_s1_1d( sendbuf, sendcount, MAXLEN, &
                                   recvbuf, recvcounts, rdispls, &
                                   root, comm, status )
    ! external:
    use MPI, only : MPI_CHAR
    use MPI, only : MPI_GatherV
    ! arguments:
    character(len=MAXLEN), intent(in)   ::  sendbuf(:)
    integer, intent(in)                 ::  sendcount
    integer, intent(in)                 ::  MAXLEN
    character(len=MAXLEN), intent(out)  ::  recvbuf(:)
    integer, intent(in)                 ::  recvcounts(:)
    integer, intent(in)                 ::  rdispls(:)
    integer, intent(in)                 ::  root
    integer, intent(in)                 ::  comm
    integer, intent(out)                ::  status
    ! string length:
    ! specific call:
    call MPI_GatherV( sendbuf, sendcount*MAXLEN                 , MPI_CHAR, &
                      recvbuf, recvcounts*MAXLEN, rdispls*MAXLEN, MPI_CHAR, &
                      root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_GatherV_s1_1d
  
  ! *
  
  subroutine MPIF90_GatherV_i4_1d( sendbuf, sendcount, &
                                   recvbuf, recvcounts, rdispls, &
                                   root, comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_GatherV
    ! arguments:
    integer(4), intent(in)        ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    integer(4), intent(out)       ::  recvbuf(:)
    integer, intent(in)           ::  recvcounts(:)
    integer, intent(in)           ::  rdispls(:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_GatherV( sendbuf, sendcount          , MPI_INTEGER, &
                      recvbuf, recvcounts, rdispls, MPI_INTEGER, &
                      root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_GatherV_i4_1d
  
  ! *

  subroutine MPIF90_GatherV_r4_1d( sendbuf, sendcount, &
                                   recvbuf, recvcounts, rdispls, &
                                   root, comm, status )
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_GatherV
    ! arguments:
    real(4), intent(in)           ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    real(4), intent(out)          ::  recvbuf(:)
    integer, intent(in)           ::  recvcounts(:)
    integer, intent(in)           ::  rdispls(:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_GatherV( sendbuf, sendcount          , MPI_REAL, &
                      recvbuf, recvcounts, rdispls, MPI_REAL, &
                      root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_GatherV_r4_1d
  
  ! *

  subroutine MPIF90_GatherV_r8_1d( sendbuf, sendcount, &
                                   recvbuf, recvcounts, rdispls, &
                                   root, comm, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_GatherV
    ! arguments:
    real(8), intent(in)           ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    real(8), intent(out)          ::  recvbuf(:)
    integer, intent(in)           ::  recvcounts(:)
    integer, intent(in)           ::  rdispls(:)
    integer, intent(in)           ::  root
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_GatherV( sendbuf, sendcount          , MPI_DOUBLE_PRECISION, &
                      recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, &
                      root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_GatherV_r8_1d


  ! ********************************************************************
  ! AllToAll
  ! ********************************************************************

  subroutine MPIF90_AllToAll_i4_1d( sendbuf, sendcount, &
                                    recvbuf, recvcount, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_AllToAll
    ! arguments:
    integer(4), intent(in)        ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    integer(4), intent(out)       ::  recvbuf(:)
    integer, intent(in)           ::  recvcount
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAll( sendbuf, sendcount, MPI_INTEGER, &
                       recvbuf, recvcount, MPI_INTEGER, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAll_i4_1d
  
  ! *

  subroutine MPIF90_AllToAll_r4_1d( sendbuf, sendcount, &
                                    recvbuf, recvcount, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_AllToAll
    ! arguments:
    real(4), intent(in)           ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    real(4), intent(out)          ::  recvbuf(:)
    integer, intent(in)           ::  recvcount
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAll( sendbuf, sendcount, MPI_REAL, &
                       recvbuf, recvcount, MPI_REAL, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAll_r4_1d
  
  ! *

  subroutine MPIF90_AllToAll_r8_1d( sendbuf, sendcount, &
                                    recvbuf, recvcount, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_AllToAll
    ! arguments:
    real(8), intent(in)           ::  sendbuf(:)
    integer, intent(in)           ::  sendcount
    real(8), intent(out)          ::  recvbuf(:)
    integer, intent(in)           ::  recvcount
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAll( sendbuf, sendcount, MPI_DOUBLE_PRECISION, &
                       recvbuf, recvcount, MPI_DOUBLE_PRECISION, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAll_r8_1d


  ! ********************************************************************
  ! AllToAllV
  ! ********************************************************************

  subroutine MPIF90_AllToAllV_i4_2d( sendbuf, sendcounts, sdispls, &
                                    recvbuf, recvcounts, rdispls, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_AllToAllV
    ! arguments:
    integer(4), intent(in)        ::  sendbuf(:,:)
    integer, intent(in)           ::  sendcounts(2)
    integer, intent(in)           ::  sdispls(2)
    integer(4), intent(out)       ::  recvbuf(:,:)
    integer, intent(in)           ::  recvcounts(2)
    integer, intent(in)           ::  rdispls(2)
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_INTEGER, &
                       recvbuf, recvcounts, rdispls, MPI_INTEGER, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAllV_i4_2d
  
  ! *

  subroutine MPIF90_AllToAllV_r4_2d( sendbuf, sendcounts, sdispls, &
                                    recvbuf, recvcounts, rdispls, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_AllToAllV
    ! arguments:
    real(4), intent(in)           ::  sendbuf(:,:)
    integer, intent(in)           ::  sendcounts(2)
    integer, intent(in)           ::  sdispls(2)
    real(4), intent(out)          ::  recvbuf(:,:)
    integer, intent(in)           ::  recvcounts(2)
    integer, intent(in)           ::  rdispls(2)
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_REAL, &
                       recvbuf, recvcounts, rdispls, MPI_REAL, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAllV_r4_2d
  
  ! *

  subroutine MPIF90_AllToAllV_r8_2d( sendbuf, sendcounts, sdispls, &
                                    recvbuf, recvcounts, rdispls, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_AllToAllV
    ! arguments:
    real(8), intent(in)           ::  sendbuf(:,:)
    integer, intent(in)           ::  sendcounts(2)
    integer, intent(in)           ::  sdispls(2)
    real(8), intent(out)          ::  recvbuf(:,:)
    integer, intent(in)           ::  recvcounts(2)
    integer, intent(in)           ::  rdispls(2)
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_DOUBLE_PRECISION, &
                       recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAllV_r8_2d
  
  ! ***

  subroutine MPIF90_AllToAllV_i4_3d( sendbuf, sendcounts, sdispls, &
                                    recvbuf, recvcounts, rdispls, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_AllToAllV
    ! arguments:
    integer(4), intent(in)        ::  sendbuf(:,:,:)
    integer, intent(in)           ::  sendcounts(3)
    integer, intent(in)           ::  sdispls(3)
    integer(4), intent(out)       ::  recvbuf(:,:,:)
    integer, intent(in)           ::  recvcounts(3)
    integer, intent(in)           ::  rdispls(3)
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_INTEGER, &
                       recvbuf, recvcounts, rdispls, MPI_INTEGER, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAllV_i4_3d
  
  ! *

  subroutine MPIF90_AllToAllV_r4_3d( sendbuf, sendcounts, sdispls, &
                                    recvbuf, recvcounts, rdispls, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_AllToAllV
    ! arguments:
    real(4), intent(in)           ::  sendbuf(:,:,:)
    integer, intent(in)           ::  sendcounts(3)
    integer, intent(in)           ::  sdispls(3)
    real(4), intent(out)          ::  recvbuf(:,:,:)
    integer, intent(in)           ::  recvcounts(3)
    integer, intent(in)           ::  rdispls(3)
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_REAL, &
                       recvbuf, recvcounts, rdispls, MPI_REAL, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAllV_r4_3d
  
  ! *

  subroutine MPIF90_AllToAllV_r8_3d( sendbuf, sendcounts, sdispls, &
                                    recvbuf, recvcounts, rdispls, &
                                    comm, status )
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_AllToAllV
    ! arguments:
    real(8), intent(in)           ::  sendbuf(:,:,:)
    integer, intent(in)           ::  sendcounts(3)
    integer, intent(in)           ::  sdispls(3)
    real(8), intent(out)          ::  recvbuf(:,:,:)
    integer, intent(in)           ::  recvcounts(3)
    integer, intent(in)           ::  rdispls(3)
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
    ! specific call:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, MPI_DOUBLE_PRECISION, &
                       recvbuf, recvcounts, rdispls, MPI_DOUBLE_PRECISION, &
                       comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! ok
    status = 0
  end subroutine MPIF90_AllToAllV_r8_3d
  

  ! ********************************************************************
  ! AllReduce
  ! ********************************************************************

  subroutine MPIF90_AllReduce_l_0d( sendbuf, recvbuf, op, comm, status )
#ifdef _MPI
    ! external:
    use MPI, only : MPI_LOGICAL
    use MPI, only : MPI_AllReduce
#endif
    ! arguments:
    logical, intent(in)           ::  sendbuf
    logical, intent(out)          ::  recvbuf
    integer, intent(in)           ::  op
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
#ifdef _MPI
    ! specific call:
    call MPI_AllReduce( sendbuf, recvbuf, 1, MPI_LOGICAL, op, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif
    ! ok
    status = 0
  end subroutine MPIF90_AllReduce_l_0d
  
  ! *

  subroutine MPIF90_AllReduce_i4_0d( sendbuf, recvbuf, op, comm, status )
#ifdef _MPI
    ! external:
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_AllReduce
#endif
    ! arguments:
    integer(4), intent(in)        ::  sendbuf
    integer(4), intent(out)       ::  recvbuf
    integer, intent(in)           ::  op
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
#ifdef _MPI
    ! specific call:
    call MPI_AllReduce( sendbuf, recvbuf, 1, MPI_INTEGER, op, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif
    ! ok
    status = 0
  end subroutine MPIF90_AllReduce_i4_0d
  
  ! *

  subroutine MPIF90_AllReduce_r4_0d( sendbuf, recvbuf, op, comm, status )
#ifdef _MPI
    ! external:
    use MPI, only : MPI_REAL
    use MPI, only : MPI_AllReduce
#endif
    ! arguments:
    real(4), intent(in)           ::  sendbuf
    real(4), intent(out)          ::  recvbuf
    integer, intent(in)           ::  op
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
#ifdef _MPI
    ! specific call:
    call MPI_AllReduce( sendbuf, recvbuf, 1, MPI_REAL, op, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif
    ! ok
    status = 0
  end subroutine MPIF90_AllReduce_r4_0d
  
  ! *

  subroutine MPIF90_AllReduce_r8_0d( sendbuf, recvbuf, op, comm, status )
#ifdef _MPI
    ! external:
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_AllReduce
#endif
    ! arguments:
    real(8), intent(in)           ::  sendbuf
    real(8), intent(out)          ::  recvbuf
    integer, intent(in)           ::  op
    integer, intent(in)           ::  comm
    integer, intent(out)          ::  status
#ifdef _MPI
    ! specific call:
    call MPI_AllReduce( sendbuf, recvbuf, 1, MPI_DOUBLE_PRECISION, op, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    recvbuf = sendbuf
#endif
    ! ok
    status = 0
  end subroutine MPIF90_AllReduce_r8_0d


  ! ********************************************************************
  ! end
  ! ********************************************************************
  
  
end module MPIF90
