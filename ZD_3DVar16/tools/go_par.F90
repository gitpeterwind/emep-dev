!#######################################################################
!
! Parallel tools using MPI.
!
! Initialization and finalization
!
!  #ifdef _MPI
!    use MPI, only : MPI_Init, MPI_Finalize, MPI_COMM_WORLD
!    use MPI, only : MPI_COMM_WORLD
!    use MPI, only : MPI_SUCCESS
!  #endif
!    use GO_Par, only : GO_Par_Init, GO_Par_Done
!
!    integer              ::  status
!    character(len=1024)  ::  msg
!    integer              ::  lenmsg
!    integer              ::  comm
!
!  #ifdef _MPI
!    ! init communication:
!    call MPI_Init( status )
!    if ( status /= MPI_SUCCESS ) then
!      call MPI_Error_String( status, msg, lenmsg, status )
!      write (*,'("Error from MPI_Init: ",a)') trim(msg)
!    end if
!    ! communicator:
!    comm = MPI_COMM_WORLD
!  #else
!    ! dummy ...
!    comm = 0
!  #endif
!
!    ! init module, provide communicator:
!    call GO_Par_Init( MPI_COMM_WORLD, status )
!    if (status/=0) stop
!
!    ! ...
!
!    ! done with module:
!    call GO_Par_Done( status )
!    if (status/=0) stop
!
!  #ifdef _MPI
!    ! done with communication:
!    call MPI_Finalize( status )
!    if ( status /= MPI_SUCCESS ) then
!      call MPI_Error_String( status, msg, lenmsg, status )
!      write (*,'("Error from MPI_Done: ",a)') trim(msg)
!    end if
!  #endif
!
!
! Preprocessor macro's
!
!  _MPI   : should be defined to enable MPI code, e.g. 'f90 -D_MPI ...'
!
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (line",i5,")")') __FILE__, __LINE__; call goErr
!
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status >0) then; TRACEBACK; action; return; end if
!
#define IF_MPI_NOT_OK_RETURN(action) if (status/=MPI_SUCCESS) then; call MPI_Error_String(status,gol,ngol,status); call goErr; TRACEBACK; action; return; end if
!
#include "go.inc"
!
!#######################################################################


module GO_Par

  use GO_Print, only : ngol, gol, goPr, goErr
#ifdef _MPI
  use MPI, only : MPI_SUCCESS, MPI_Error_String
#endif

  implicit none


  ! --- in/out ------------------------

  private

  public  ::  nproc, me
 
  public  ::  T_Domains
 
  public  ::  GO_Par_Init, GO_Par_Done, GO_Par_Setup
  public  ::  GO_Par_Barrier
  public  ::  GO_Par_Abort


  ! --- const --------------------------

  character(len=*), parameter   ::  mname = 'GO_Par'


  ! --- local  --------------------------
  
  ! mpi variables:
  integer             ::  comm
  integer             ::  nproc, me
  ! flag:
  logical             ::  module_setup


  ! --- types -----------------------------------
  
  ! domain ranges from all process to facilitate decomposition
  type T_Domains
    ! rank:
    integer                ::  ndim
    ! global bounds:
    integer, allocatable   ::  glbo(:,:)  ! (ndim,0:nproc-1)
    integer, allocatable   ::  gubo(:,:)  ! (ndim,0:nproc-1)
    ! local bounds:
    integer, allocatable   ::  lbo(:,:)   ! (ndim,0:nproc-1)
    integer, allocatable   ::  ubo(:,:)   ! (ndim,0:nproc-1)
    ! offset:
    integer, allocatable   ::  off(:,:)   ! (ndim,0:nproc-1)
    ! shape:
    integer, allocatable   ::  shp(:,:)   ! (ndim,0:nproc-1)
    ! total number of local elements:
    integer, allocatable   ::  n(:)   ! (0:nproc-1)
    !
  contains
    procedure   ::  Init         => Domains_Init
    procedure   ::  InitLocal    => Domains_InitLocal
    procedure   ::  Done         => Domains_Done
    procedure   ::  Inside       => Domains_Inside
    procedure   ::  Find         => Domains_Find
    procedure   ::  Intersection => Domains_Intersection
    procedure   ::                  Domains_Swap_2d_r8
    procedure   ::                  Domains_Swap_3d_r8
    procedure   ::                  Domains_Swap_4d_r8
    generic     ::  Swap         => Domains_Swap_2d_r8, &
                                    Domains_Swap_3d_r8, &
                                    Domains_Swap_4d_r8
    procedure   ::                  Domains_AllGather_1d_i4
    procedure   ::                  Domains_AllGather_1d_r8
    procedure   ::                  Domains_AllGather_3d_c8
    generic     ::  AllGather    => Domains_AllGather_1d_i4, &
                                    Domains_AllGather_1d_r8, &
                                    Domains_AllGather_3d_c8
    procedure   ::  Extract      => Domains_Extract_1d_r8
  end type T_Domains


  
contains


  ! ***************************************************************************
  ! ***
  ! *** module init/done
  ! ***
  ! ***************************************************************************
  
  
  subroutine GO_Par_Init( status )

    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GO_Par_Init'
    
    ! --- local -----------------------------
    
    ! --- begin -----------------------------
    
    ! dummy init as single pe:
    comm = 0
    nproc = 1
    me = 0
    
    ! not setup yet ...
    module_setup = .false.
    
    ! ok
    status = 0
    
  end subroutine GO_Par_Init
  
  
  ! ***
  
  
  subroutine GO_Par_Done( status )

    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GO_Par_Done'
    
    ! --- begin -----------------------------
  
    ! ok
    status = 0
    
  end subroutine GO_Par_Done


  ! ***************************************************************************
  ! ***
  ! *** module setup (after initialization)
  ! ***
  ! ***************************************************************************
  
  
  subroutine GO_Par_Setup( mpi_comm, status )

#ifdef _MPI
    use MPI, only : MPI_Init
    use MPI, only : MPI_COMM_WORLD
    use MPI, only : MPI_Comm_Size, MPI_Comm_Rank
#endif
  
    ! --- in/out ----------------------------
    
    integer, intent(in)            ::  mpi_comm
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GO_Par_Setup'
    
    ! --- local -----------------------------
    
    ! --- begin -----------------------------
    
#ifdef _MPI
    ! communicator:
    comm = mpi_comm
    ! number of process:
    call MPI_Comm_Size( comm, nproc, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! own:
    call MPI_Comm_Rank( comm, me, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! info ...
    write (gol,'("could not setup GO_Par without MPI environment; define macro _MPI ...")'); call goErr
    TRACEBACK; status=1; return
#endif
    
    ! reset flag:
    module_setup = .true.
    
    ! ok
    status = 0
    
  end subroutine GO_Par_Setup
  
  


  ! ***************************************************************************
  ! ***
  ! *** tools
  ! ***
  ! ***************************************************************************
  
  subroutine GO_Par_Barrier( status )
  
#ifdef _MPI
    use MPI, only : MPI_Barrier
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/GO_Par_Barrier'
    
    ! --- in/out ---------------------------------
    
    integer, intent(out)        ::  status
    
    ! --- local -----------------------------------
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if

#ifdef _MPI
    ! start:
    call MPI_Barrier( comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#endif

    ! ok
    status = 0
    
  end subroutine GO_Par_Barrier
  
  
  ! ***
  
  
  subroutine GO_Par_Abort( ierr, status )

#ifdef _MPI
    use MPI, only : MPI_Abort
#endif
  
    ! --- in/out ----------------------------
    
    integer, intent(in)           ::  ierr
    integer, intent(out)          ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/GO_Par_Done'          
    
    ! --- begin -----------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if
  
#ifdef _MPI
    ! start:
    call MPI_Abort( comm, ierr, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! system exit (non-standard routine!)
    call Exit( 1 )
#endif
    
    ! ok
    status = 0
    
  end subroutine GO_Par_Abort


  ! ********************************************************************
  ! ***
  ! *** Domains 2D
  ! ***
  ! ********************************************************************

  
  !
  ! initialize local domain using global bounds
  !

  subroutine Domains_Init( self, glbo, gubo, status )

    use MPI, only : MPI_AllGather
    use MPI, only : MPI_INTEGER
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Init'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(out)     ::  self
    integer, intent(in)               ::  glbo(:)
    integer, intent(in)               ::  gubo(:)
    integer, intent(out)              ::  status

    ! --- local ----------------------------------
    
    integer    ::  i
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! count:
    self%ndim = size(glbo)
    ! check ...
    if ( size(gubo) /= self%ndim ) then
      write (gol,'("size of argument gubo (",i0,") does not match with size of glbo (",i0,")")') &
                      size(gubo), self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( self%glbo(self%ndim,0:nproc-1) )
    allocate( self%gubo(self%ndim,0:nproc-1) )
    allocate( self%off(self%ndim,0:nproc-1) )
    allocate( self%lbo(self%ndim,0:nproc-1) )
    allocate( self%ubo(self%ndim,0:nproc-1) )
    allocate( self%shp(self%ndim,0:nproc-1) )
    allocate( self%n(0:nproc-1) )
    
#ifdef _MPI
    ! exchange global bounds:
    call MPI_AllGather(      glbo, self%ndim, MPI_INTEGER, &
                        self%glbo, self%ndim, MPI_INTEGER, &
                        comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    call MPI_AllGather(      gubo, self%ndim, MPI_INTEGER, &
                        self%gubo, self%ndim, MPI_INTEGER, &
                        comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    self%glbo = glbo
    self%gubo = gubo
#endif
    
    ! local offset in global space:
    self%off = self%glbo - 1
    
    ! local shapes:
    self%shp = self%gubo - self%glbo + 1
    ! trap undefined per process:
    do i = 0, nproc-1
      ! any dimension undefined ? set all to zero:
      if ( any(self%shp(:,i) <= 0) ) self%shp(:,i) = 0
    end do
    
    ! total number of local elements,
    ! will be zero for processes with empty domain:
    self%n = product( self%shp, dim=1 )
    
    ! local bounds:
    self%lbo = 1
    self%ubo = self%shp
    
    ! ok
    status = 0

  end subroutine Domains_Init


  ! ***
  
  
  subroutine Domains_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(inout)      ::  self
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! clear:
    deallocate( self%glbo )
    deallocate( self%gubo )
    deallocate( self%lbo )
    deallocate( self%ubo )
    deallocate( self%off )
    deallocate( self%shp )
    deallocate( self%n   )
    
    ! ok
    status = 0

  end subroutine Domains_Done


  ! ***

  
  !
  ! initialize subdomains in 1D array;
  ! each pe calls this routine with the local size,
  ! global position is then determined by assuming
  ! that all sub-domains are in processor ordering
  !

  subroutine Domains_InitLocal( self, n, status )

#ifdef _MPI
    use MPI, only : MPI_AllGather
    use MPI, only : MPI_INTEGER
#endif
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Init'
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(out)     ::  self
    integer, intent(in)               ::  n
    integer, intent(out)              ::  status

    ! --- local ----------------------------------
    
    integer, allocatable   ::  counts(:)     ! (nproc)
    integer, allocatable   ::  displs(:)     ! (nproc)
    integer                ::  k
    
    ! --- begin ----------------------------------
    
    ! storage:
    allocate( counts(0:nproc-1) )
    allocate( displs(0:nproc-1) )
    
#ifdef _MPI
    ! exchange global bounds:
    call MPI_AllGather(      n, 1, MPI_INTEGER, &
                        counts, 1, MPI_INTEGER, &
                        comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    counts = n
#endif

    ! compute displacements:
    displs(0) = 0
    do k = 1, nproc-1
      displs(k) = displs(k-1) + counts(k-1)
    end do
    
    ! general init for 1D arrays, specify global bounds:
    call self%Init( (/displs(me)+1/), (/displs(me)+counts(me)/), status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    deallocate( counts )
    deallocate( displs )
    
    ! ok
    status = 0

  end subroutine Domains_InitLocal

  ! ***
  
  !
  ! Compute local intersection bounds and shape.
  ! Return status:
  !   -1  : no intersection
  !    0  : ok
  !   else error
  !

  subroutine Domains_Intersection( self, glbo, gubo, lbo, ubo, shp, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)         ::  self
    integer, intent(in)                  ::  glbo(:)
    integer, intent(in)                  ::  gubo(:)
    integer, intent(out)                 ::  lbo(:)
    integer, intent(out)                 ::  ubo(:)
    integer, intent(out)                 ::  shp(:)
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( any( (/ size(glbo), size(gubo), size(lbo), size(ubo), size(shp) /) /= self%ndim ) ) then
      write (gol,'("arguments should have size ndim = ",i0)') self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check on empty domains first:
    if ( any( self%shp(:,me) <= 0 ) ) then

      ! empty local domain:
      status = -1

    else if ( any( gubo < glbo ) ) then

      ! empty inquired domain:
      status = -1

    else

      ! local bounds of intersection:
      lbo = max( self%glbo(:,me), glbo ) - self%off(:,me)
      ubo = min( self%gubo(:,me), gubo ) - self%off(:,me)
      ! local shape:
      shp = ubo - lbo + 1

      ! set return status:
      if ( any(shp <= 0) ) then
        status = -1
      else
        status = 0
      end if
      
    end if
    
  end subroutine Domains_Intersection


  ! ***
  
  
  ! return status:
  !   -1  : no
  !    0  : yes
  !   else error
  
  subroutine Domains_Inside( self, ind, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)            ::  self
    integer, intent(in)                     ::  ind(:)  ! global indices
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Inside'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( size(ind) /= self%ndim ) then
      write (gol,'("argument ind has size ",i0," while ndim = ",i0)') size(ind), self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! compare:
    if ( all(self%glbo(:,me) <= ind) .and. all(ind <= self%gubo(:,me)) ) then
      ! in domain:
      status = 0
    else
      ! outside ..
      status = -1
    end if
    
  end subroutine Domains_Inside
  
  
  ! index of domain holding cell with provided global indices
  
  subroutine Domains_Find( self, ind, iproc, status )
    
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)            ::  self
    integer, intent(in)                     ::  ind(:)  ! global indices
    integer, intent(out)                    ::  iproc   ! 0:nproc-1
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Find'
    
    ! --- local ----------------------------------
    
    integer     ::  k
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( size(ind) /= self%ndim ) then
      write (gol,'("argument ind has size ",i0," while ndim = ",i0)') size(ind), self%ndim; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! dummy result:
    iproc = -999
    ! loop:
    do k = 0, nproc-1
      ! compare:
      if ( all(self%glbo(:,k) <= ind) .and. all(ind <= self%gubo(:,k)) ) then
        ! found !
        iproc = k
        ! leave:
        exit
      end if  ! location on proc domain
    end do  ! procs
    ! check ...
    if ( iproc < 0 ) then
      write (gol,*) 'could not find domain holding global indices : ', ind; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine Domains_Find


  ! ***************************************************************************
  ! ***
  ! *** swap decomposition
  ! ***
  ! ***************************************************************************


  subroutine Domains_Swap_2d_r8( self, input, doms, output, status )
  
#ifdef _MPI
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_2d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 2
    integer, parameter  ::  wpr = 8
#ifdef _MPI
    integer, parameter  ::  dtype = MPI_DOUBLE_PRECISION
#endif
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:,:)
    type(T_Domains), intent(in)   ::  doms
    real(wpr), intent(out)        ::  output(:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer, allocatable     ::  sendcounts(:)  ! (nproc)
    integer, allocatable     ::  sdispls(:)     ! (nproc)
    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable     ::  recvcounts(:)  ! (nproc)
    integer, allocatable     ::  rdispls(:)     ! (nproc)
    integer                  ::  pid
    integer                  ::  ilbo(ndim)
    integer                  ::  iubo(ndim)
    integer                  ::  ishp(ndim)
    integer                  ::  n
    integer                  ::  displ
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! local input slab ?
    if ( all( self%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,me) ) ) then
        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
                      shape(input), self%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! local output slab ?
    if ( all( doms%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(output) /= doms%shp(:,me) ) ) then
        write (gol,'("shape of output (",i0,3(",",i0),") while doms define (",i0,3(",",i0),")")') &
                      shape(output), doms%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty output domain, careful ...")'); call goPr
    end if
    
    ! ~ create send buffer
    
    ! storage:
    allocate( sendbuf(self%n(me)) )
    allocate( sendcounts(0:nproc-1) )
    allocate( sdispls   (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over target processes:
    do pid = 0, nproc-1
    
      ! intersection of local input domain with output domain on target process:
      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,1(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2), &
        !        displ+1, displ+n; call goPr
        ! copy slab into buffer:
        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
                                                  ilbo(2):iubo(2)), (/n/) )
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      sendcounts(pid) = n
      ! displacement:
      sdispls(pid) = displ
      
      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= self%n(me) ) then
      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)') displ, self%n(me); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(doms%n(me)) )
    allocate( recvcounts(0:nproc-1) )
    allocate( rdispls   (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, nproc-1
    
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      rdispls(pid) = displ

      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= doms%n(me) ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, doms%n(me); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
    !write (gol,*) '   sdispls      : ', sdispls; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
    !write (gol,*) '   rdispls      : ', rdispls; call goPr
    ! transfer data:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
                        recvbuf, recvcounts, rdispls, dtype, &
                        comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, nproc-1
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
      else if ( status == 0 ) then
        ! valid intersection found; extract:
        n     = recvcounts(pid)
        displ = rdispls(pid)
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,1(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2); call goPr
        ! extract output slab from buffer:
        output(ilbo(1):iubo(1),&
               ilbo(2):iubo(2)) = reshape( recvbuf(displ+1:displ+n), ishp )
      else
        TRACEBACK; status=1; return
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( rdispls )
    
    ! clear:
    deallocate( sendbuf )
    deallocate( sendcounts )
    deallocate( sdispls )
    
    ! ok
    status = 0
    
  end subroutine Domains_Swap_2d_r8


  ! ***

  subroutine Domains_Swap_3d_r8( self, input, doms, output, status )
  
#ifdef _MPI
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_3d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 3
    integer, parameter  ::  wpr = 8
#ifdef _MPI
    integer, parameter  ::  dtype = MPI_DOUBLE_PRECISION
#endif
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:,:,:)
    type(T_Domains), intent(in)   ::  doms
    real(wpr), intent(out)        ::  output(:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer, allocatable     ::  sendcounts(:)  ! (nproc)
    integer, allocatable     ::  sdispls(:)     ! (nproc)
    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable     ::  recvcounts(:)  ! (nproc)
    integer, allocatable     ::  rdispls(:)     ! (nproc)
    integer                  ::  pid
    integer                  ::  ilbo(ndim)
    integer                  ::  iubo(ndim)
    integer                  ::  ishp(ndim)
    integer                  ::  n
    integer                  ::  displ
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! local input slab ?
    if ( all( self%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,me) ) ) then
        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
                      shape(input), self%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! local output slab ?
    if ( all( doms%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(output) /= doms%shp(:,me) ) ) then
        write (gol,'("shape of output (",i0,3(",",i0),") while doms define (",i0,3(",",i0),")")') &
                      shape(output), doms%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty output domain, careful ...")'); call goPr
    end if
    
    ! ~ create send buffer
    
    ! storage:
    allocate( sendbuf(self%n(me)) )
    allocate( sendcounts(0:nproc-1) )
    allocate( sdispls   (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over target processes:
    do pid = 0, nproc-1
    
      ! intersection of local input domain with output domain on target process:
      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3), &
        !        displ+1, displ+n; call goPr
        ! copy slab into buffer:
        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
                                                  ilbo(2):iubo(2),&
                                                  ilbo(3):iubo(3)), (/n/) )
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      sendcounts(pid) = n
      ! displacement:
      sdispls(pid) = displ
      
      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= self%n(me) ) then
      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)') displ, self%n(me); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(doms%n(me)) )
    allocate( recvcounts(0:nproc-1) )
    allocate( rdispls   (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, nproc-1
    
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      rdispls(pid) = displ

      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= doms%n(me) ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, doms%n(me); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
    !write (gol,*) '   sdispls      : ', sdispls; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
    !write (gol,*) '   rdispls      : ', rdispls; call goPr
    ! transfer data:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
                        recvbuf, recvcounts, rdispls, dtype, &
                        comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, nproc-1
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
      else if ( status == 0 ) then
        ! valid intersection found; extract:
        n     = recvcounts(pid)
        displ = rdispls(pid)
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3); call goPr
        ! extract output slab from buffer:
        output(ilbo(1):iubo(1),&
               ilbo(2):iubo(2),&
               ilbo(3):iubo(3)) = reshape( recvbuf(displ+1:displ+n), ishp )
      else
        TRACEBACK; status=1; return
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( rdispls )
    
    ! clear:
    deallocate( sendbuf )
    deallocate( sendcounts )
    deallocate( sdispls )
    
    ! ok
    status = 0
    
  end subroutine Domains_Swap_3d_r8


  ! ***


  subroutine Domains_Swap_4d_r8( self, input, doms, output, status )
  
#ifdef _MPI
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_4d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 4
    integer, parameter  ::  wpr = 8
#ifdef _MPI
    integer, parameter  ::  dtype = MPI_DOUBLE_PRECISION
#endif
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:,:,:,:)
    type(T_Domains), intent(in)   ::  doms
    real(wpr), intent(out)        ::  output(:,:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer, allocatable     ::  sendcounts(:)  ! (nproc)
    integer, allocatable     ::  sdispls(:)     ! (nproc)
    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable     ::  recvcounts(:)  ! (nproc)
    integer, allocatable     ::  rdispls(:)     ! (nproc)
    integer                  ::  pid
    integer                  ::  ilbo(ndim)
    integer                  ::  iubo(ndim)
    integer                  ::  ishp(ndim)
    integer                  ::  n
    integer                  ::  displ
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! local input slab ?
    if ( all( self%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,me) ) ) then
        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
                      shape(input), self%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! local output slab ?
    if ( all( doms%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(output) /= doms%shp(:,me) ) ) then
        write (gol,'("shape of output (",i0,3(",",i0),") while doms define (",i0,3(",",i0),")")') &
                      shape(output), doms%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty output domain, careful ...")'); call goPr
    end if
    
    ! ~ create send buffer
    
    ! storage:
    allocate( sendbuf(self%n(me)) )
    allocate( sendcounts(0:nproc-1) )
    allocate( sdispls   (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over target processes:
    do pid = 0, nproc-1
    
      ! intersection of local input domain with output domain on target process:
      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,3(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3),ilbo(4),iubo(4), &
        !        displ+1, displ+n; call goPr
        ! copy slab into buffer:
        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
                                                  ilbo(2):iubo(2),&
                                                  ilbo(3):iubo(3),&
                                                  ilbo(4):iubo(4)), (/n/) )
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      sendcounts(pid) = n
      ! displacement:
      sdispls(pid) = displ
      
      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= self%n(me) ) then
      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)') displ, self%n(me); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(doms%n(me)) )
    allocate( recvcounts(0:nproc-1) )
    allocate( rdispls   (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, nproc-1
    
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
        n = 0
      else if ( status == 0 ) then
        ! valid intersection found; size:
        n = product(ishp)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        TRACEBACK; status=1; return
      end if
      
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      rdispls(pid) = displ

      ! increase displacement:
      displ = displ + n
    
    end do ! iproc
    ! check ...
    if ( displ /= doms%n(me) ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, doms%n(me); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
    !write (gol,*) '   sdispls      : ', sdispls; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
    !write (gol,*) '   rdispls      : ', rdispls; call goPr
    ! transfer data:
    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
                        recvbuf, recvcounts, rdispls, dtype, &
                        comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, nproc-1
      ! intersection of local output domain with input domain on source process:
      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
                               ilbo, iubo, ishp, status )
      if ( status == -1 ) then
        ! no intersection
      else if ( status == 0 ) then
        ! valid intersection found; extract:
        n     = recvcounts(pid)
        displ = rdispls(pid)
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,3(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3),ilbo(4),iubo(4); call goPr
        ! extract output slab from buffer:
        output(ilbo(1):iubo(1),&
               ilbo(2):iubo(2),&
               ilbo(3):iubo(3),&
               ilbo(4):iubo(4)) = reshape( recvbuf(displ+1:displ+n), ishp )
      else
        TRACEBACK; status=1; return
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( rdispls )
    
    ! clear:
    deallocate( sendbuf )
    deallocate( sendcounts )
    deallocate( sdispls )
    
    ! ok
    status = 0
    
  end subroutine Domains_Swap_4d_r8
  


  ! ***************************************************************************
  ! ***
  ! *** gather and broadcast to all
  ! ***
  ! ***************************************************************************


  subroutine Domains_AllGather_1d_i4( self, input, output, status )
  
#ifdef _MPI
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_AllGather_1d_i4'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 1
    integer, parameter  ::  wpr = 4
#ifdef _MPI
    integer, parameter  ::  dtype = MPI_INTEGER
#endif
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    integer(wpr), intent(in)      ::  input(:)
    integer(wpr), intent(out)     ::  output(:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    integer(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    integer(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! local input slab ?
    if ( all( self%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,me) ) ) then
        write (gol,'("shape of input (",i0,2(",",i0),") while expected (",i0,2(",",i0),")")') &
                      shape(input), self%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, nproc-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! check ...
    if ( any( shape(output) /= gshp ) ) then
      write (gol,'("shape of output (",i0,2(",",i0),") while assumed global shape (",i0,2(",",i0),")")') &
                    shape(output), gshp; call goErr
      TRACEBACK; status=1; return
    end if
    ! total:
    gn = product( gshp )
    ! offset, arguments are 1-based:
    goff = glbo - 1
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(me)
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n) )
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,me),self%ubo(1,me),self%lbo(2,me),self%ubo(2,me),&
      !           self%lbo(3,me),self%ubo(3,me); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,me):self%ubo(1,me)), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1) )
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(gn) )
    allocate( recvcounts(0:nproc-1) )
    allocate( displs    (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, nproc-1
      ! source defined ?
      if ( self%n(pid) > 0 ) then
        ! all values:
        n = self%n(pid)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        ! empty:
        n = 0
      end if
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      displs(pid) = displ
      ! increase displacement:
      displ = displ + n
    end do ! iproc
    ! check ...
    if ( displ /= gn ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_AllGatherV( sendbuf, sendcount, dtype, &
                         recvbuf, recvcounts, displs, dtype, &
                         comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, nproc-1
      ! extract:
      n     = recvcounts(pid)
      displ = displs(pid)
      ! any received from this process ?
      if ( n > 0 ) then
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
        !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
        ! copy shape, needed to avoid strange compilation error ... 
        lshp = self%shp(:,pid)
        ! extract output slab from buffer:
        output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1)) = &
            reshape( recvbuf(displ+1:displ+n), lshp )
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( displs )
    
    ! clear:
    deallocate( sendbuf )
    
    ! ok
    status = 0
    
  end subroutine Domains_AllGather_1d_i4
  
  
  ! ***


  subroutine Domains_AllGather_1d_r8( self, input, output, status )
  
#ifdef _MPI
    use MPI, only : MPI_DOUBLE_PRECISION
    use MPI, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_AllGather_1d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 1
    integer, parameter  ::  wpr = 8
#ifdef _MPI
    integer, parameter  ::  dtype = MPI_DOUBLE_PRECISION
#endif
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:)
    real(wpr), intent(out)        ::  output(:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    real(wpr), allocatable      ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    real(wpr), allocatable      ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! local input slab ?
    if ( all( self%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,me) ) ) then
        write (gol,'("shape of input (",i0,") while expected (",i0,")")') &
                      shape(input), self%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, nproc-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! check ...
    if ( any( shape(output) /= gshp ) ) then
      write (gol,'("shape of output (",i0,") while assumed global shape (",i0,")")') &
                    shape(output), gshp; call goErr
      TRACEBACK; status=1; return
    end if
    ! total:
    gn = product( gshp )
    ! offset, arguments are 1-based:
    goff = glbo - 1
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(me)
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n) )
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,me),self%ubo(1,me),self%lbo(2,me),self%ubo(2,me),&
      !           self%lbo(3,me),self%ubo(3,me); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,me):self%ubo(1,me)), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1) )
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(gn) )
    allocate( recvcounts(0:nproc-1) )
    allocate( displs    (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, nproc-1
      ! source defined ?
      if ( self%n(pid) > 0 ) then
        ! all values:
        n = self%n(pid)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        ! empty:
        n = 0
      end if
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      displs(pid) = displ
      ! increase displacement:
      displ = displ + n
    end do ! iproc
    ! check ...
    if ( displ /= gn ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_AllGatherV( sendbuf, sendcount, dtype, &
                         recvbuf, recvcounts, displs, dtype, &
                         comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, nproc-1
      ! extract:
      n     = recvcounts(pid)
      displ = displs(pid)
      ! any received from this process ?
      if ( n > 0 ) then
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
        !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
        ! copy shape, needed to avoid strange compilation error ... 
        lshp = self%shp(:,pid)
        ! extract output slab from buffer:
        output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1)) = &
            reshape( recvbuf(displ+1:displ+n), lshp )
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( displs )
    
    ! clear:
    deallocate( sendbuf )
    
    ! ok
    status = 0
    
  end subroutine Domains_AllGather_1d_r8


  ! ***
  

  subroutine Domains_AllGather_3d_c8( self, input, output, status )
  
#ifdef _MPI
    use MPI, only : MPI_DOUBLE_COMPLEX
    use MPI, only : MPI_AllToAllV
#endif

    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_AllGather_3d_c8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 3
    integer, parameter  ::  wpr = 8
#ifdef _MPI
    integer, parameter  ::  dtype = MPI_DOUBLE_COMPLEX
#endif
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    complex(wpr), intent(in)      ::  input(:,:,:)
    complex(wpr), intent(out)     ::  output(:,:,:)
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------
    
    integer                     ::  gn
    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  lshp(ndim)
    complex(wpr), allocatable   ::  sendbuf(:)     ! (input n)
    integer                     ::  sendcount
    complex(wpr), allocatable   ::  recvbuf(:)     ! (output n)
    integer, allocatable        ::  recvcounts(:)  ! (nproc)
    integer, allocatable        ::  displs(:)      ! (nproc)
    integer                     ::  pid
    integer                     ::  n
    integer                     ::  displ
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! local input slab ?
    if ( all( self%shp(:,me) > 0 ) ) then
      ! check ...
      if ( any( shape(input) /= self%shp(:,me) ) ) then
        write (gol,'("shape of input (",i0,2(",",i0),") while expected (",i0,2(",",i0),")")') &
                      shape(input), self%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if
    !else
    !  ! info ...
    !  write (gol,'("empty input domain, careful ...")'); call goPr
    end if

    ! global bounds:
    glbo = self%glbo(:,0)
    gubo = self%gubo(:,0)
    do pid = 1, nproc-1
      if ( self%n(pid) > 0 ) then
        glbo = min( glbo, self%glbo(:,pid) )
        gubo = max( gubo, self%gubo(:,pid) )
      end if
    end do
    ! global shape:
    gshp = gubo - glbo + 1
    ! check ...
    if ( any( shape(output) /= gshp ) ) then
      write (gol,'("shape of output (",i0,2(",",i0),") while assumed global shape (",i0,2(",",i0),")")') &
                    shape(output), gshp; call goErr
      TRACEBACK; status=1; return
    end if
    ! total:
    gn = product( gshp )
    ! offset, arguments are 1-based:
    goff = glbo - 1
    
    ! ~ create send buffer
    
    ! local size:
    n = self%n(me)
    ! store:
    sendcount = n
    ! defined ?
    if ( n > 0 ) then
      ! storage:
      allocate( sendbuf(n) )
      !! info ...
      !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf")') &
      !        n, self%lbo(1,me),self%ubo(1,me),self%lbo(2,me),self%ubo(2,me),&
      !           self%lbo(3,me),self%ubo(3,me); call goPr
      ! copy slab into buffer:
      sendbuf = reshape( input(self%lbo(1,me):self%ubo(1,me),&
                               self%lbo(2,me):self%ubo(2,me),&
                               self%lbo(3,me):self%ubo(3,me)), (/n/) )
    else
      ! dummy:
      allocate( sendbuf(1) )
    end if

    ! ~ prepare receive buffer
    
    ! storage:
    allocate( recvbuf(gn) )
    allocate( recvcounts(0:nproc-1) )
    allocate( displs    (0:nproc-1) )
    
    ! init displacement:
    displ = 0
    ! loop over source processes:
    do pid = 0, nproc-1
      ! source defined ?
      if ( self%n(pid) > 0 ) then
        ! all values:
        n = self%n(pid)
        !! info ...
        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
      else
        ! empty:
        n = 0
      end if
      ! number of values:
      recvcounts(pid) = n
      ! displacement:
      displs(pid) = displ
      ! increase displacement:
      displ = displ + n
    end do ! iproc
    ! check ...
    if ( displ /= gn ) then
      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') displ, gn; call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ 
    
#ifdef _MPI
    !! info ...
    !write (gol,'("exchange data ...")'); call goPr
    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
    !write (gol,*) '   sendcount   : ', sendcount; call goPr
    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
    !write (gol,*) '   recvcounts  : ', recvcounts; call goPr
    !write (gol,*) '   displs      : ', displs; call goPr
    ! transfer data:
    call MPI_AllGatherV( sendbuf, sendcount, dtype, &
                         recvbuf, recvcounts, displs, dtype, &
                         comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    !! info ...
    !write (gol,'("copy data ...")'); call goPr
    ! copy:
    recvbuf = sendbuf
#endif

    ! ~ unpack

    ! loop over source processes:
    do pid = 0, nproc-1
      ! extract:
      n     = recvcounts(pid)
      displ = displs(pid)
      ! any received from this process ?
      if ( n > 0 ) then
        !! info ...
        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
        !        n, displ+1, displ+n, self%glbo(1,pid),self%gubo(1,pid),self%glbo(2,pid),self%gubo(2,pid),&
        !                             self%glbo(3,pid),self%gubo(3,pid); call goPr
        ! copy shape, needed to avoid strange compilation error ... 
        lshp = self%shp(:,pid)
        ! extract output slab from buffer:
        output(self%glbo(1,pid)-goff(1):self%gubo(1,pid)-goff(1),&
               self%glbo(2,pid)-goff(2):self%gubo(2,pid)-goff(2),&
               self%glbo(3,pid)-goff(3):self%gubo(3,pid)-goff(3)) = &
            reshape( recvbuf(displ+1:displ+n), lshp )
      end if
    end do ! iproc

    ! clear:
    deallocate( recvbuf )
    deallocate( recvcounts )
    deallocate( displs )
    
    ! clear:
    deallocate( sendbuf )
    
    ! ok
    status = 0
    
  end subroutine Domains_AllGather_3d_c8
  

  ! ***************************************************************************
  ! ***
  ! *** extract local slab
  ! ***
  ! ***************************************************************************


  subroutine Domains_Extract_1d_r8( self, input, output, status )
  
    ! --- const ---------------------------------
    
    character(len=*), parameter   ::  rname = mname//'/Domains_Extract_1d_r8'
    
    ! rank and kind:
    integer, parameter  ::  ndim = 1
    integer, parameter  ::  wpr = 8
  
    ! --- in/out ---------------------------------
    
    class(T_Domains), intent(in)  ::  self
    real(wpr), intent(in)         ::  input(:)  ! global array
    real(wpr), intent(out)        ::  output(:)  ! local array
    integer, intent(out)          ::  status
    
    ! --- local -----------------------------------

    integer                     ::  glbo(ndim), gubo(ndim), gshp(ndim), goff(ndim)
    integer                     ::  pid
    
    ! --- begin -----------------------------------
    
    ! check ...
    if ( .not. module_setup ) then
      write (gol,'("module not setup yet ..")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! local slab has non-zero shape ?
    if ( all( self%shp(:,me) > 0 ) ) then

      ! global bounds:
      glbo = self%glbo(:,0)
      gubo = self%gubo(:,0)
      do pid = 1, nproc-1
        if ( self%n(pid) > 0 ) then
          glbo = min( glbo, self%glbo(:,pid) )
          gubo = max( gubo, self%gubo(:,pid) )
        end if
      end do
      ! global shape:
      gshp = gubo - glbo + 1
      ! check ...
      if ( any( shape(input) /= gshp ) ) then
        write (gol,'("shape of input (",i0,") while expected global shape (",i0,")")') &
                      shape(input), gshp; call goErr
        TRACEBACK; status=1; return
      end if
      ! offset, arguments are 1-based:
      goff = glbo - 1

      ! check ...
      if ( any( shape(output) /= self%shp(:,me) ) ) then
        write (gol,'("shape of output (",i0,") while expected local shape (",i0,")")') &
                      shape(output), self%shp(:,me); call goErr
        TRACEBACK; status=1; return
      end if

      !! testing ..
      !print *, me, ': extract ', self%glbo(1,me)-goff(1), self%gubo(1,me)-goff(1)
      
      ! copy local slab:
      output = input(self%glbo(1,me)-goff(1):self%gubo(1,me)-goff(1))

    end if

    ! ok
    status = 0
    
  end subroutine Domains_Extract_1d_r8


end module GO_Par
