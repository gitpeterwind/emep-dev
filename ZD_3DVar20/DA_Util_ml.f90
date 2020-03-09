!*****************************************************************************!
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!*****************************************************************************!

module DA_Util_ml

  use GO, only : ngol, gol, goPr, goErr

!!#ifdef with_ajs
!!  use GO
!  use GO, only : T_Domains
!  use GO, only : ngol, gol, goPr, goErr, GO_Print_Set
!  !se go,       only: GO_Init, GO_Done
!  !se go_par,   only: GO_Par_Setup, T_Domains
!  !se go_print, only: ngol, gol, goPr, goErr, GO_Print_Set
!  !se go_timer, only: GO_Timer_Def, GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
!  !se go_string,only: goReadFromLine, goMatchValue, goCharToString
!  !se go_date,  only: TDate, Get, AnyDate, IsAnyDate, TIncrDate, &
!  !                   Extract_Ref_and_Step, operator(-), rTotal, Num_to_Date
!!#else
!!  use Config_module    , only : nproc
!!  use Par_mod          , only : me
!!  use MPI              , only : MPI_SUCCESS
!!  use MPI              , only : MPI_INTEGER, MPI_DOUBLE_PRECISION
!!#endif


  ! --- in/out -------------------------------------

  private 

  public :: DA_Util_Init, DA_Util_Done
!  public :: goGetFU
!  public :: gol, ngol, goPr, goErr, GO_Print_Set
!  public :: T_Domains
!!#ifdef with_ajs
!!  public :: GO_Timer_Def, GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
!!#endif
  public :: norm


  ! --- const --------------------------------------

  character(len=*), parameter :: mname = 'DA_Util_ml'
  
!#ifndef with_ajs
!  ! range of file units that might be used by this program:
!  integer, parameter    ::  goFuRange(2) = (/200,999/)
!
!  ! go_print
!  !logical, parameter  :: GO_QUIET = .true.
!  logical, parameter  :: GO_QUIET = .false.
!#endif


  ! --- var --------------------------------------

!#ifndef with_ajs
!  ! go_print
!  character(len=1024) ::  gol = '' ! buffer for standard output 
!  integer             :: ngol = 0  ! call MPI_Error_String(status,gol,ngol,status)
!  logical             :: pr_apply = .not. GO_QUIET
!
!  ! go_par
!  integer             :: comm=0
!  logical             :: module_setup=.false.
!#endif


  ! --- types --------------------------------------

!#ifndef with_ajs
!  ! domain ranges from all process to facilitate decomposition
!  type T_Domains
!    ! rank:
!    integer                ::  ndim
!    ! global bounds:
!    integer, allocatable   ::  glbo(:,:)  ! (ndim,0:nproc-1)
!    integer, allocatable   ::  gubo(:,:)  ! (ndim,0:nproc-1)
!    ! local bounds:
!    integer, allocatable   ::  lbo(:,:)   ! (ndim,0:nproc-1)
!    integer, allocatable   ::  ubo(:,:)   ! (ndim,0:nproc-1)
!    ! offset:
!    integer, allocatable   ::  off(:,:)   ! (ndim,0:nproc-1)
!    ! shape:
!    integer, allocatable   ::  shp(:,:)   ! (ndim,0:nproc-1)
!    ! total number of local elements:
!    integer, allocatable   ::  n(:)   ! (0:nproc-1)
!    !
!  contains
!    procedure :: Init         => Domains_Init
!    procedure :: InitLocal    => Domains_InitLocal
!    procedure :: Done         => Domains_Done
!    procedure :: Inside       => Domains_Inside
!    procedure :: Find         => Domains_Find
!    procedure :: Intersection => Domains_Intersection
!    procedure ::                 Domains_Swap_2d_r8
!    procedure ::                 Domains_Swap_3d_r8
!    procedure ::                 Domains_Swap_4d_r8
!    generic   :: Swap         => Domains_Swap_2d_r8, &
!                                 Domains_Swap_3d_r8, &
!                                 Domains_Swap_4d_r8
!  ! procedure ::                 Domains_AllGather_1d_i4
!  ! procedure ::                 Domains_AllGather_1d_r8
!  ! procedure ::                 Domains_AllGather_3d_c8
!  ! generic   :: AllGather    => Domains_AllGather_1d_i4, &
!  !                              Domains_AllGather_1d_r8, &
!  !                              Domains_AllGather_3d_c8
!  ! procedure :: Extract      => Domains_Extract_1d_r8
!  end type T_Domains
!
!#endif


  ! --- interfaces -------------------------------

  interface norm
    module procedure norm_r8d1
    module procedure norm_r8d2
    module procedure norm_r8d3
    module procedure norm_r8d4
    module procedure norm_c8d1
    module procedure norm_c8d2
    module procedure norm_c8d3
    module procedure norm_c8d4
  end interface norm


contains


  !+------------------------------------------------------------------
  ! module init/done
  !+------------------------------------------------------------------
  
  subroutine DA_Util_Init( status )

    use GO            , only : GO_Init
    use GO            , only : GO_Print_Set
    use GO            , only : GO_Par_Setup
    use Config_module , only : masterProc
    use MPI_Groups_mod, only : MPI_COMM_CALC

    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/DA_Util_Init'
    
    ! --- local -----------------------------
    
    ! --- begin -----------------------------
  
    ! initialize GO tools:
    call GO_Init( status )
    if ( status /= 0 ) then
      write (*,'("in ",a," (line",i5,")")') __FILE__, __LINE__
      stop
    end if
    ! from now on, the gol/goPr/goPr logging could be used ...

    ! log from root only:
    call GO_Print_Set( status, apply=MasterProc, prompt_pe=.true. )
    IF_NOT_OK_RETURN(status=1)
  
    ! setup parallel tools in GO modules:
    call GO_Par_Setup( MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine DA_Util_Init
  
  
  ! ***
  
  
  subroutine DA_Util_Done( status )

    use GO, only : GO_Done
    
    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/DA_Util_Done'
    
    ! --- begin -----------------------------

    ! done with GO modules:
    call GO_Done( status )
    if ( status /= 0 ) then
      write (*,'("in ",a," (line",i5,")")') __FILE__, __LINE__
      stop
    end if
  
    ! ok
    status = 0
    
  end subroutine DA_Util_Done

!#ifndef with_ajs
!
!  !+------------------------------------------------------------------
!  ! go_fu
!  !+------------------------------------------------------------------
!  
!  ! Return the first free available file unit number.
!
!  subroutine goGetFU( fu, status )
!  
!    ! --- in/out --------------------------
!
!    integer, intent(out)      ::  fu
!    integer, intent(out)      ::  status
!
!    ! --- const ---------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/goGetFU'
!    
!    ! --- local --------------------------
!
!    integer               ::  i
!    character(len=256)    ::  fname
!    logical               ::  opened
!
!    ! --- local ---------------------------
!    
!    ! start with lowest possible unit:
!    fu = goFuRange(1) - 1
!    
!    ! loop until unopned unit is found:
!    do
!
!      ! try next file unit:
!      fu = fu + 1
!
!      ! too large ?
!      if ( fu > goFuRange(2) ) then
!        write (*,'("unable to select free file unit within allowed range ...")')
!        write (*,'("close some files or increase goFuRange in module GO_FU")')
!        write (*,'("current goFuRange : ",i6," .. ",i6)') goFuRange
!        write (*,'("open files:")')
!        do i = goFuRange(1), goFuRange(2)
!          inquire( unit=i, name=fname )
!          write (*,'(i6," : ",a)') i, trim(fname)
!        end do
!        write (*,'("in ",a)') rname; status=1; return
!      end if
!      
!      !! skip ?
!      !if ( fu==goStdIn  ) cycle
!      !if ( fu==goStdOut ) cycle
!      !if ( fu==goStdErr ) cycle
!
!      ! free available unit ? then ok
!      inquire( unit=fu, opened=opened )
!      if ( .not. opened ) exit
!      
!    end do
!
!    ! ok
!    status = 0
!
!  end subroutine goGetFU
!
!  !+------------------------------------------------------------------
!  ! go_print
!  !+------------------------------------------------------------------
!
!  subroutine goPr()
!    if ( pr_apply ) write (*,'(a)') trim(gol)
!    gol=''
!  end subroutine
!
!  subroutine goErr()
!    write (*,'("ERROR - ",a)') trim(gol)
!    gol=''
!  end subroutine
!
!  subroutine GO_Print_Set( status, apply, prompt_pe, pe, trace )
!    integer, intent(out)          :: status
!    logical, intent(in), optional :: apply
!    logical, intent(in), optional :: prompt_pe
!    integer, intent(in), optional :: pe
!    logical, intent(in), optional :: trace
!
!    if(present(apply)) pr_apply = apply.and..not.GO_QUIET
!    status=0
!  end subroutine
!
!  !+------------------------------------------------------------------
!  ! go_par
!  !+------------------------------------------------------------------
!
!  subroutine GO_Par_Setup( mpi_comm, status )
!
!    ! --- in/out ----------------------------    
!
!    integer, intent(in)            ::  mpi_comm
!    integer, intent(out)           ::  status
!
!    ! --- const ----------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/GO_Par_Setup'
!
!    ! --- begin -----------------------------
!
!#ifndef _MPI
!    ERROR(gol="could not setup GO_Par without MPI environment; define macro _MPI ...")
!#else
!    comm=mpi_comm       ! communicator
!    module_setup=.true. ! reset flag
!    status=0            ! ok
!#endif
!
!  end subroutine GO_Par_Setup
!
!  !+------------------------------------------------------------------
!  ! Domains 2D
!  !+------------------------------------------------------------------
!
!  subroutine Domains_Init( self, glbo, gubo, status )
!    ! initialize local domain using global bounds
!    ! --- const --------------------------------------
!    character(len=*), parameter  ::  rname = mname//'/Domains_Init'
!    ! --- in/out ---------------------------------
!    class(T_Domains), intent(out)     ::  self
!    integer, intent(in)               ::  glbo(:)
!    integer, intent(in)               ::  gubo(:)
!    integer, intent(out)              ::  status
!    ! --- local ----------------------------------
!    integer    ::  i
!    ! --- begin ----------------------------------
!
!    ! check ...
!    if ( .not. module_setup ) then
!      ERROR(gol="module not setup yet ..")
!    endif
!
!    ! count:
!    self%ndim = size(glbo)
!    ! check ...
!    if(size(gubo)/=self%ndim ) then
!      write(gol,'("size of argument gubo (",i0,") &
!        &does not match with size of glbo (",i0,")")')size(gubo), self%ndim
!      ERROR()
!    endif
!
!    ! storage:
!    allocate( self%glbo(self%ndim,0:nproc-1) )
!    allocate( self%gubo(self%ndim,0:nproc-1) )
!    allocate( self%off(self%ndim,0:nproc-1) )
!    allocate( self%lbo(self%ndim,0:nproc-1) )
!    allocate( self%ubo(self%ndim,0:nproc-1) )
!    allocate( self%shp(self%ndim,0:nproc-1) )
!    allocate( self%n(0:nproc-1) )
!
!    ! exchange global bounds:
!    call MPI_AllGather(      glbo, self%ndim, MPI_INTEGER, &
!                        self%glbo, self%ndim, MPI_INTEGER, &
!                        comm, status )
!    IF_NOT_OK_ERROR(MPI_SUCCESS,call MPI_Error_String(status,gol,ngol,status))
!
!    call MPI_AllGather(      gubo, self%ndim, MPI_INTEGER, &
!                        self%gubo, self%ndim, MPI_INTEGER, &
!                        comm, status )
!    IF_NOT_OK_ERROR(MPI_SUCCESS,call MPI_Error_String(status,gol,ngol,status))
!
!    ! local offset in global space:
!    self%off = self%glbo - 1
!
!    ! local shapes:
!    self%shp = self%gubo - self%glbo + 1
!    ! trap undefined per process:
!    do i = 0, nproc-1
!      ! any dimension undefined ? set all to zero:
!      if ( any(self%shp(:,i) <= 0) ) self%shp(:,i) = 0
!    enddo
!
!    ! total number of local elements,
!    ! will be zero for processes with empty domain:
!    self%n = product( self%shp, dim=1 )
!
!    ! local bounds:
!    self%lbo = 1
!    self%ubo = self%shp
!
!    ! ok
!    status = 0
!
!  end subroutine Domains_Init
!  
!  ! *
!
!  subroutine Domains_Done( self, status )
!  ! --- in/out ---------------------------------
!    class(T_Domains), intent(inout)      ::  self
!    integer, intent(out)                 ::  status
!  ! --- const --------------------------------------
!    character(len=*), parameter  ::  rname = mname//'/Domains_Done'
!  ! --- begin ----------------------------------
!
!    ! check ...
!    if( .not. module_setup ) then
!      ERROR(gol="module not setup yet ..")
!    endif
!
!    ! clear:
!    deallocate( self%glbo )
!    deallocate( self%gubo )
!    deallocate( self%lbo )
!    deallocate( self%ubo )
!    deallocate( self%off )
!    deallocate( self%shp )
!    deallocate( self%n   )
!
!    ! ok
!    status = 0
!  end subroutine Domains_Done
!
!  subroutine Domains_InitLocal( self, n, status )
!  ! initialize subdomains in 1D array;
!  ! each pe calls this routine with the local size,
!  ! global position is then determined by assuming
!  ! that all sub-domains are in processor ordering
!  ! --- const --------------------------------------
!    character(len=*), parameter  ::  rname = mname//'/Domains_Init'
!  ! --- in/out ---------------------------------
!    class(T_Domains), intent(out)     ::  self
!    integer, intent(in)               ::  n
!    integer, intent(out)              ::  status
!  ! --- local ----------------------------------
!    integer, allocatable   ::  counts(:)     ! (nproc)
!    integer, allocatable   ::  displs(:)     ! (nproc)
!    integer                ::  k
!  ! --- begin ----------------------------------
!
!    ! storage:
!    allocate( counts(0:nproc-1) )
!    allocate( displs(0:nproc-1) )
!
!    ! exchange global bounds:
!    call MPI_AllGather(      n, 1, MPI_INTEGER, &
!                        counts, 1, MPI_INTEGER, &
!                        comm, status )
!    IF_NOT_OK_ERROR(MPI_SUCCESS,call MPI_Error_String(status,gol,ngol,status))
!
!    ! compute displacements:
!    displs(0) = 0
!    do k = 1, nproc-1
!      displs(k) = displs(k-1) + counts(k-1)
!    enddo
!
!    ! general init for 1D arrays, specify global bounds:
!    call self%Init( [displs(me)+1], [displs(me)+counts(me)], status )
!    IF_NOT_OK_ERROR(0,)
!
!    ! storage:
!    deallocate( counts )
!    deallocate( displs )
!
!    ! ok
!    status = 0
!  end subroutine Domains_InitLocal
!
!
!  subroutine Domains_Intersection( self, glbo, gubo, lbo, ubo, shp, status )
!  ! Compute local intersection bounds and shape.
!  ! Return status:
!  !   -1  : no intersection
!  !    0  : ok
!  !   else error
!
!  ! --- in/out ---------------------------------
!    class(T_Domains), intent(in)         ::  self
!    integer, intent(in)                  ::  glbo(:)
!    integer, intent(in)                  ::  gubo(:)
!    integer, intent(out)                 ::  lbo(:)
!    integer, intent(out)                 ::  ubo(:)
!    integer, intent(out)                 ::  shp(:)
!    integer, intent(out)                 ::  status
!  ! --- const --------------------------------------
!    character(len=*), parameter  ::  rname = mname//'/Domains_Done'
!  ! --- begin ----------------------------------
!
!    ! check ...
!    if( .not. module_setup ) then
!      ERROR(gol="module not setup yet ..")
!    endif
!
!    ! check ..
!    if(any([size(glbo),size(gubo),size(lbo),size(ubo),size(shp)] /= self%ndim)) then
!      write(gol,'("arguments should have size ndim = ",i0)') self%ndim
!      ERROR()
!    endif
!
!    ! check on empty domains first:
!    if ( any( self%shp(:,me) <= 0 ) ) then
!
!      ! empty local domain:
!      status = -1
!
!    elseif ( any( gubo < glbo ) ) then
!
!      ! empty inquired domain:
!      status = -1
!
!    else
!
!      ! local bounds of intersection:
!      lbo = max( self%glbo(:,me), glbo ) - self%off(:,me)
!      ubo = min( self%gubo(:,me), gubo ) - self%off(:,me)
!      ! local shape:
!      shp = ubo - lbo + 1
!
!      ! set return status:
!      if( any(shp <= 0) ) then
!        status = -1
!      else
!        status = 0
!      endif
!    endif
!  end subroutine Domains_Intersection
!
!  subroutine Domains_Inside( self, ind, status )
!  ! return status:
!  !   -1  : no
!  !    0  : yes
!  !   else error
!
!  ! --- in/out ---------------------------------
!    class(T_Domains), intent(in)            ::  self
!    integer, intent(in)                     ::  ind(:)  ! global indices
!    integer, intent(out)                    ::  status
!  ! --- const --------------------------------------
!    character(len=*), parameter  ::  rname = mname//'/Domains_Inside'
!  ! --- begin ----------------------------------
!
!    ! check ...
!    if( .not. module_setup ) then
!      ERROR(gol="module not setup yet ..")
!    endif
!
!      ! check ..
!    if ( size(ind) /= self%ndim ) then
!      write(gol,'("argument ind has size ",i0," while ndim = ",i0)')&
!        size(ind), self%ndim
!      ERROR()
!    endif
!
!    ! compare:
!    if ( all(self%glbo(:,me) <= ind) .and. all(ind <= self%gubo(:,me)) ) then
!      ! in domain:
!      status = 0
!    else
!      ! outside ..
!      status = -1
!    endif
!  end subroutine Domains_Inside
!
!  subroutine Domains_Find( self, ind, iproc, status )
!  ! index of domain holding cell with provided global indices
!
!  ! --- in/out ---------------------------------
!    class(T_Domains), intent(in)            ::  self
!    integer, intent(in)                     ::  ind(:)  ! global indices
!    integer, intent(out)                    ::  iproc   ! 0:nproc-1
!    integer, intent(out)                    ::  status
!  ! --- const --------------------------------------
!    character(len=*), parameter  ::  rname = mname//'/Domains_Find'
!  ! --- local ----------------------------------
!    integer     ::  k
!  ! --- begin ----------------------------------
!
!    ! check ...
!    if( .not. module_setup ) then
!      ERROR(gol="module not setup yet ..")
!    endif
!
!      ! check ..
!    if ( size(ind) /= self%ndim ) then
!      write (gol,'("argument ind has size ",i0," while ndim = ",i0)') &
!        size(ind), self%ndim
!      ERROR()
!    endif
!
!    ! dummy result:
!    iproc = -999
!    ! loop:
!    do k = 0, nproc-1
!      ! compare:
!      if ( all(self%glbo(:,k) <= ind) .and. all(ind <= self%gubo(:,k)) ) then
!        ! found !
!        iproc = k
!        ! leave:
!        exit
!      endif  ! location on proc domain
!    enddo  ! procs
!
!    ! check ...
!    if ( iproc < 0 ) then
!      write (gol,*) 'could not find domain holding global indices : ', ind
!      ERROR()
!    endif
!
!    ! ok
!    status = 0
!  end subroutine Domains_Find
!  !+------------------------------------------------------------------
!  ! swap decomposition
!  !+------------------------------------------------------------------
!  subroutine Domains_Swap_2d_r8( self, input, doms, output, status )
!  ! --- const ---------------------------------
!    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_2d_r8'
!    ! rank and kind:
!    integer, parameter :: ndim = 2
!    integer, parameter :: wpr = 8
!    integer, parameter :: dtype = MPI_DOUBLE_PRECISION
!  ! --- in/out ---------------------------------
!    class(T_Domains), intent(in)  ::  self
!    real(wpr), intent(in)         ::  input(:,:)
!    type(T_Domains), intent(in)   ::  doms
!    real(wpr), intent(out)        ::  output(:,:)
!    integer, intent(out)          ::  status
!  ! --- local -----------------------------------
!    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
!    integer, allocatable     ::  sendcounts(:)  ! (nproc)
!    integer, allocatable     ::  sdispls(:)     ! (nproc)
!    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
!    integer, allocatable     ::  recvcounts(:)  ! (nproc)
!    integer, allocatable     ::  rdispls(:)     ! (nproc)
!    integer                  ::  pid
!    integer                  ::  ilbo(ndim)
!    integer                  ::  iubo(ndim)
!    integer                  ::  ishp(ndim)
!    integer                  ::  n
!    integer                  ::  displ
!  ! --- begin -----------------------------------
!
!    ! check ...
!    if( .not. module_setup ) then
!      ERROR(gol="module not setup yet ..")
!    endif
!
!    ! local input slab ?
!    if ( all( self%shp(:,me) > 0 ) ) then
!      ! check ...
!      if ( any( shape(input) /= self%shp(:,me) ) ) then
!        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
!                      shape(input), self%shp(:,me)
!        ERROR()
!      endif
!    endif
!
!    ! local output slab ?
!    if ( all( doms%shp(:,me) > 0 ) ) then
!      ! check ...
!      if ( any( shape(output) /= doms%shp(:,me) ) ) then
!        write (gol,'("shape of output (",i0,3(",",i0),")&
!          & while doms define (",i0,3(",",i0),")")') &
!          shape(output), doms%shp(:,me)
!        ERROR()
!      endif
!    endif
!
!    ! ~ create send buffer
!
!    ! storage:
!    allocate( sendbuf(self%n(me)) )
!    allocate( sendcounts(0:nproc-1) )
!    allocate( sdispls   (0:nproc-1) )
!
!    ! init displacement:
!    displ = 0
!    ! loop over target processes:
!    do pid = 0, nproc-1
!
!      ! intersection of local input domain with output domain on target process:
!      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if ( status == -1 ) then
!        ! no intersection
!        n = 0
!      elseif( status == 0 ) then
!        ! valid intersection found; size:
!        n = product(ishp)
!        !! info ...
!        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,1(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
!        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2), &
!        !        displ+1, displ+n; call goPr
!        ! copy slab into buffer:
!        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
!                                                  ilbo(2):iubo(2)), (/n/) )
!      else
!        ERROR()
!      endif
!
!      ! number of values:
!      sendcounts(pid) = n
!      ! displacement:
!      sdispls(pid) = displ
!
!      ! increase displacement:
!      displ = displ + n
!
!    enddo ! iproc
!
!    ! check ...
!    if ( displ /= self%n(me) ) then
!      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)')&
!        displ, self%n(me)
!      ERROR()
!    endif
!
!    ! ~ prepare receive buffer
!
!    ! storage:
!    allocate( recvbuf(doms%n(me)) )
!    allocate( recvcounts(0:nproc-1) )
!    allocate( rdispls   (0:nproc-1) )
!
!    ! init displacement:
!    displ = 0
!    ! loop over source processes:
!    do pid = 0, nproc-1
!
!      ! intersection of local output domain with input domain on source process:
!      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if ( status == -1 ) then
!        ! no intersection
!        n = 0
!      elseif ( status == 0 ) then
!        ! valid intersection found; size:
!        n = product(ishp)
!        !! info ...
!        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
!      else
!        TRACEBACK; status=1; return
!      endif
!
!      ! number of values:
!      recvcounts(pid) = n
!      ! displacement:
!      rdispls(pid) = displ
!
!      ! increase displacement:
!      displ = displ + n
!
!    enddo ! iproc
!
!    ! check ...
!    if ( displ /= doms%n(me) ) then
!      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)')&
!        displ, doms%n(me)
!      ERROR()
!    endif
!
!  !! info ...
!  !write (gol,'("exchange data ...")'); call goPr
!  !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
!  !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
!  !write (gol,*) '   sdispls      : ', sdispls; call goPr
!  !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
!  !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
!  !write (gol,*) '   rdispls      : ', rdispls; call goPr
!
!    ! transfer data:
!    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
!                        recvbuf, recvcounts, rdispls, dtype, &
!                        comm, status )
!    IF_NOT_OK_ERROR(MPI_SUCCESS,call MPI_Error_String(status,gol,ngol,status))
!
!    ! ~ unpack
!
!    ! loop over source processes:
!    do pid = 0, nproc-1
!      ! intersection of local output domain with input domain on source process:
!      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if ( status == -1 ) then
!        ! no intersection
!      elseif ( status == 0 ) then
!        ! valid intersection found; extract:
!        n     = recvcounts(pid)
!        displ = rdispls(pid)
!        !! info ...
!        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,1(",",i0,":",i0),")")') &
!        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2); call goPr
!        ! extract output slab from buffer:
!        output(ilbo(1):iubo(1),&
!               ilbo(2):iubo(2)) = reshape( recvbuf(displ+1:displ+n), ishp )
!      else
!        ERROR()
!      endif
!    enddo ! iproc
!
!    ! clear:
!    deallocate( recvbuf )
!    deallocate( recvcounts )
!    deallocate( rdispls )
!
!    ! clear:
!    deallocate( sendbuf )
!    deallocate( sendcounts )
!    deallocate( sdispls )
!
!    ! ok
!    status = 0
!  end subroutine Domains_Swap_2d_r8
!
!  subroutine Domains_Swap_3d_r8( self, input, doms, output, status )
!  ! --- const ---------------------------------
!    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_3d_r8'
!    ! rank and kind:
!    integer, parameter :: ndim = 3
!    integer, parameter :: wpr = 8
!    integer, parameter :: dtype = MPI_DOUBLE_PRECISION
!  ! --- in/out ---------------------------------
!    class(T_Domains), intent(in)  ::  self
!    real(wpr), intent(in)         ::  input(:,:,:)
!    type(T_Domains), intent(in)   ::  doms
!    real(wpr), intent(out)        ::  output(:,:,:)
!    integer, intent(out)          ::  status
!  ! --- local -----------------------------------
!    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
!    integer, allocatable     ::  sendcounts(:)  ! (nproc)
!    integer, allocatable     ::  sdispls(:)     ! (nproc)
!    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
!    integer, allocatable     ::  recvcounts(:)  ! (nproc)
!    integer, allocatable     ::  rdispls(:)     ! (nproc)
!    integer                  ::  pid
!    integer                  ::  ilbo(ndim)
!    integer                  ::  iubo(ndim)
!    integer                  ::  ishp(ndim)
!    integer                  ::  n
!    integer                  ::  displ
!  ! --- begin -----------------------------------
!
!    ! check ...
!    if( .not. module_setup ) then
!      ERROR(gol="module not setup yet ..")
!    endif
!
!    ! local input slab ?
!    if ( all( self%shp(:,me) > 0 ) ) then
!      ! check ...
!      if ( any( shape(input) /= self%shp(:,me) ) ) then
!        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
!                      shape(input), self%shp(:,me)
!        ERROR()
!      endif
!    endif
!
!    ! local output slab ?
!    if ( all( doms%shp(:,me) > 0 ) ) then
!      ! check ...
!      if ( any( shape(output) /= doms%shp(:,me) ) ) then
!        write (gol,'("shape of output (",i0,3(",",i0),") while doms define (",i0,3(",",i0),")")') &
!                      shape(output), doms%shp(:,me)
!        ERROR()
!      endif
!    endif
!
!    ! ~ create send buffer
!
!    ! storage:
!    allocate( sendbuf(self%n(me)) )
!    allocate( sendcounts(0:nproc-1) )
!    allocate( sdispls   (0:nproc-1) )
!
!    ! init displacement:
!    displ = 0
!    ! loop over target processes:
!    do pid = 0, nproc-1
!
!      ! intersection of local input domain with output domain on target process:
!      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if( status == -1 ) then
!        ! no intersection
!        n = 0
!      elseif ( status == 0 ) then
!        ! valid intersection found; size:
!        n = product(ishp)
!        !! info ...
!        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,2(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
!        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3), &
!        !        displ+1, displ+n; call goPr
!        ! copy slab into buffer:
!        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
!                                                  ilbo(2):iubo(2),&
!                                                  ilbo(3):iubo(3)), (/n/) )
!      else
!        ERROR()
!      endif
!
!      ! number of values:
!      sendcounts(pid) = n
!      ! displacement:
!      sdispls(pid) = displ
!
!      ! increase displacement:
!      displ = displ + n
!
!    enddo ! iproc
!
!    ! check ...
!    if ( displ /= self%n(me) ) then
!      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)') &
!        displ, self%n(me)
!      ERROR()
!    endif
!
!    ! ~ prepare receive buffer
!
!    ! storage:
!    allocate( recvbuf(doms%n(me)) )
!    allocate( recvcounts(0:nproc-1) )
!    allocate( rdispls   (0:nproc-1) )
!
!    ! init displacement:
!    displ = 0
!    ! loop over source processes:
!    do pid = 0, nproc-1
!
!      ! intersection of local output domain with input domain on source process:
!      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if ( status == -1 ) then
!        ! no intersection
!        n = 0
!      elseif ( status == 0 ) then
!        ! valid intersection found; size:
!        n = product(ishp)
!        !! info ...
!        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
!      else
!        ERROR()
!      endif
!
!      ! number of values:
!      recvcounts(pid) = n
!      ! displacement:
!      rdispls(pid) = displ
!
!      ! increase displacement:
!      displ = displ + n
!
!    end do ! iproc
!
!    ! check ...
!    if ( displ /= doms%n(me) ) then
!      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)') &
!        displ, doms%n(me)
!      ERROR()
!    endif
!
!    !! info ...
!    !write (gol,'("exchange data ...")'); call goPr
!    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
!    !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
!    !write (gol,*) '   sdispls      : ', sdispls; call goPr
!    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
!    !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
!    !write (gol,*) '   rdispls      : ', rdispls; call goPr
!    ! transfer data:
!    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
!                        recvbuf, recvcounts, rdispls, dtype, &
!                        comm, status )
!    IF_NOT_OK_ERROR(MPI_SUCCESS,call MPI_Error_String(status,gol,ngol,status))
!
!    ! ~ unpack
!
!    ! loop over source processes:
!    do pid = 0, nproc-1
!      ! intersection of local output domain with input domain on source process:
!      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if ( status == -1 ) then
!        ! no intersection
!      else if ( status == 0 ) then
!        ! valid intersection found; extract:
!        n     = recvcounts(pid)
!        displ = rdispls(pid)
!        !! info ...
!        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,2(",",i0,":",i0),")")') &
!        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3); call goPr
!        ! extract output slab from buffer:
!        output(ilbo(1):iubo(1),&
!               ilbo(2):iubo(2),&
!               ilbo(3):iubo(3)) = reshape( recvbuf(displ+1:displ+n), ishp )
!      else
!        ERROR()
!      endif
!    enddo ! iproc
!
!    ! clear:
!    deallocate( recvbuf )
!    deallocate( recvcounts )
!    deallocate( rdispls )
!
!    ! clear:
!    deallocate( sendbuf )
!    deallocate( sendcounts )
!    deallocate( sdispls )
!
!    ! ok
!    status = 0
!  end subroutine Domains_Swap_3d_r8
!
!  subroutine Domains_Swap_4d_r8( self, input, doms, output, status )
!  ! --- const ---------------------------------
!    character(len=*), parameter   ::  rname = mname//'/Domains_Swap_4d_r8'
!    ! rank and kind:
!    integer, parameter :: ndim = 4
!    integer, parameter :: wpr = 8
!    integer, parameter :: dtype = MPI_DOUBLE_PRECISION
!  ! --- in/out ---------------------------------
!    class(T_Domains), intent(in)  ::  self
!    real(wpr), intent(in)         ::  input(:,:,:,:)
!    type(T_Domains), intent(in)   ::  doms
!    real(wpr), intent(out)        ::  output(:,:,:,:)
!    integer, intent(out)          ::  status
!  ! --- local -----------------------------------
!    real(wpr), allocatable   ::  sendbuf(:)     ! (input n)
!    integer, allocatable     ::  sendcounts(:)  ! (nproc)
!    integer, allocatable     ::  sdispls(:)     ! (nproc)
!    real(wpr), allocatable   ::  recvbuf(:)     ! (output n)
!    integer, allocatable     ::  recvcounts(:)  ! (nproc)
!    integer, allocatable     ::  rdispls(:)     ! (nproc)
!    integer                  ::  pid
!    integer                  ::  ilbo(ndim)
!    integer                  ::  iubo(ndim)
!    integer                  ::  ishp(ndim)
!    integer                  ::  n
!    integer                  ::  displ
!  ! --- begin -----------------------------------
!
!    ! check ...
!    if( .not. module_setup ) then
!      ERROR(gol="module not setup yet ..")
!    endif
!
!    ! local input slab ?
!    if ( all( self%shp(:,me) > 0 ) ) then
!      ! check ...
!      if ( any( shape(input) /= self%shp(:,me) ) ) then
!        write (gol,'("shape of input (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
!                      shape(input), self%shp(:,me)
!        ERROR()
!      endif
!    endif
!
!    ! local output slab ?
!    if ( all( doms%shp(:,me) > 0 ) ) then
!      ! check ...
!      if ( any( shape(output) /= doms%shp(:,me) ) ) then
!        write (gol,'("shape of output (",i0,3(",",i0),") while doms define (",i0,3(",",i0),")")') &
!                      shape(output), doms%shp(:,me)
!        ERROR()
!      endif
!    endif
!
!    ! ~ create send buffer
!
!    ! storage:
!    allocate( sendbuf(self%n(me)) )
!    allocate( sendcounts(0:nproc-1) )
!    allocate( sdispls   (0:nproc-1) )
!
!    ! init displacement:
!    displ = 0
!    ! loop over target processes:
!    do pid = 0, nproc-1
!
!      ! intersection of local input domain with output domain on target process:
!      call self%Intersection( doms%glbo(:,pid), doms%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if ( status == -1 ) then
!        ! no intersection
!        n = 0
!      elseif ( status == 0 ) then
!        ! valid intersection found; size:
!        n = product(ishp)
!        !! info ...
!        !write (gol,'("collect ",i0," elements from input(",i0,":",i0,3(",",i0,":",i0),") into sendbuf(",i0,":",i0,")")') &
!        !        n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3),ilbo(4),iubo(4), &
!        !        displ+1, displ+n; call goPr
!        ! copy slab into buffer:
!        sendbuf(displ+1:displ+n) = reshape( input(ilbo(1):iubo(1),&
!                                                  ilbo(2):iubo(2),&
!                                                  ilbo(3):iubo(3),&
!                                                  ilbo(4):iubo(4)), (/n/) )
!      else
!        ERROR()
!      endif
!
!      ! number of values:
!      sendcounts(pid) = n
!      ! displacement:
!      sdispls(pid) = displ
!
!      ! increase displacement:
!      displ = displ + n
!
!    end do ! iproc
!
!    ! check ...
!    if ( displ /= self%n(me) ) then
!      write (gol,'("stored ",i0," values in sendbuf, but local input size is ",i0)')&
!         displ, self%n(me)
!      ERROR()
!    endif
!
!    ! ~ prepare receive buffer
!
!    ! storage:
!    allocate( recvbuf(doms%n(me)) )
!    allocate( recvcounts(0:nproc-1) )
!    allocate( rdispls   (0:nproc-1) )
!
!    ! init displacement:
!    displ = 0
!    ! loop over source processes:
!    do pid = 0, nproc-1
!
!      ! intersection of local output domain with input domain on source process:
!      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if ( status == -1 ) then
!        ! no intersection
!        n = 0
!      elseif ( status == 0 ) then
!        ! valid intersection found; size:
!        n = product(ishp)
!        !! info ...
!        !write (gol,'("prepare recvbuf to hold ",i0," elements from process ",i0)') n, pid; call goPr
!      else
!        ERROR()
!      endif
!
!      ! number of values:
!      recvcounts(pid) = n
!      ! displacement:
!      rdispls(pid) = displ
!
!      ! increase displacement:
!      displ = displ + n
!
!    enddo ! iproc
!
!    ! check ...
!    if ( displ /= doms%n(me) ) then
!      write (gol,'("prepared ",i0," values in recvbuf, but local output size is ",i0)')&
!         displ, doms%n(me)
!      ERROR()
!    endif
!
!    !! info ...
!    !write (gol,'("exchange data ...")'); call goPr
!    !write (gol,'("   size sendbuf : ",i0)') size(sendbuf); call goPr
!    !write (gol,*) '   sendcounts   : ', sendcounts; call goPr
!    !write (gol,*) '   sdispls      : ', sdispls; call goPr
!    !write (gol,'("   size recvbuf : ",i0)') size(recvbuf); call goPr
!    !write (gol,*) '   recvcounts   : ', recvcounts; call goPr
!    !write (gol,*) '   rdispls      : ', rdispls; call goPr
!    ! transfer data:
!    call MPI_AllToAllV( sendbuf, sendcounts, sdispls, dtype, &
!                        recvbuf, recvcounts, rdispls, dtype, &
!                        comm, status )
!    IF_NOT_OK_ERROR(MPI_SUCCESS,call MPI_Error_String(status,gol,ngol,status))
!
!    ! ~ unpack
!
!    ! loop over source processes:
!    do pid = 0, nproc-1
!      ! intersection of local output domain with input domain on source process:
!      call doms%Intersection( self%glbo(:,pid), self%gubo(:,pid), &
!                               ilbo, iubo, ishp, status )
!      if ( status == -1 ) then
!        ! no intersection
!      elseif ( status == 0 ) then
!        ! valid intersection found; extract:
!        n     = recvcounts(pid)
!        displ = rdispls(pid)
!        !! info ...
!        !write (gol,'("store ",i0," elements from recvbuf(",i0,":",i0,") into output(",i0,":",i0,3(",",i0,":",i0),")")') &
!        !        n, displ+1, displ+n, ilbo(1),iubo(1),ilbo(2),iubo(2),ilbo(3),iubo(3),ilbo(4),iubo(4); call goPr
!        ! extract output slab from buffer:
!        output(ilbo(1):iubo(1),&
!               ilbo(2):iubo(2),&
!               ilbo(3):iubo(3),&
!               ilbo(4):iubo(4)) = reshape( recvbuf(displ+1:displ+n), ishp )
!      else
!        ERROR()
!      endif
!    enddo ! iproc
!
!    ! clear:
!    deallocate( recvbuf )
!    deallocate( recvcounts )
!    deallocate( rdispls )
!
!    ! clear:
!    deallocate( sendbuf )
!    deallocate( sendcounts )
!    deallocate( sdispls )
!
!    ! ok
!    status = 0
!  end subroutine Domains_Swap_4d_r8
!
!#endif
!! with_ajs defined

  !+------------------------------------------------------------------
  ! norms
  !+------------------------------------------------------------------

  function norm_r8d1(a,squared) result(norm)
    implicit none
    intent(in)      :: a,squared
    optional        :: squared
    logical         :: squared,sq
    real(kind=8)    :: a(:),norm
    sq=.false.;if(present(squared))sq=squared
    norm=sum(a**2)
    if(.not.sq)norm=sqrt(norm)
  end function norm_r8d1

  function norm_c8d1(a,squared) result(norm)
    implicit none
    intent(in)      :: a,squared
    optional        :: squared
    logical         :: squared,sq
    complex(kind=8) :: a(:)
    real(kind=8)    :: norm
    sq=.false.;if(present(squared))sq=squared
    norm=sum(conjg(a)*a)
    if(.not.sq)norm=sqrt(norm)
  end function norm_c8d1
  
  function norm_r8d2(a,squared) result(norm)
    implicit none
    intent(in)      :: a,squared
    optional        :: squared
    logical         :: squared,sq
    real(kind=8)    :: a(:,:),norm
    sq=.false.;if(present(squared))sq=squared
    norm=sum(a**2)
    if(.not.sq)norm=sqrt(norm)
  end function norm_r8d2
  
  function norm_c8d2(a,squared) result(norm)
    implicit none
    intent(in)      :: a,squared
    optional        :: squared
    logical         :: squared,sq
    complex(kind=8) :: a(:,:)
    real(kind=8)    :: norm
    sq=.false.;if(present(squared))sq=squared
    norm=sum(conjg(a)*a)
    if(.not.sq)norm=sqrt(norm)
  end function norm_c8d2
  
  function norm_r8d3(a,squared) result(norm)
    implicit none
    intent(in)      :: a,squared
    optional        :: squared
    logical         :: squared,sq
    real(kind=8)    :: a(:,:,:),norm
    sq=.false.;if(present(squared))sq=squared
    norm=sum(a**2)
    if(.not.sq)norm=sqrt(norm)
  end function norm_r8d3
  
  function norm_c8d3(a,squared) result(norm)
    implicit none
    intent(in)      :: a,squared
    optional        :: squared
    logical         :: squared,sq
    complex(kind=8) :: a(:,:,:)
    real(kind=8)    :: norm
    sq=.false.;if(present(squared))sq=squared
    norm=sum(conjg(a)*a)
    if(.not.sq)norm=sqrt(norm)
  end function norm_c8d3
  
  function norm_r8d4(a,squared) result(norm)
    implicit none
    intent(in)      :: a,squared
    optional        :: squared
    logical         :: squared,sq
    real(kind=8)    :: a(:,:,:,:),norm
    sq=.false.;if(present(squared))sq=squared
    norm=sum(a**2)
    if(.not.sq)norm=sqrt(norm)
  end function norm_r8d4
  
  function norm_c8d4(a,squared) result(norm)
    implicit none
    intent(in)      :: a,squared
    optional        :: squared
    logical         :: squared,sq
    complex(kind=8) :: a(:,:,:,:)
    real(kind=8)    :: norm
    sq=.false.;if(present(squared))sq=squared
    norm=sum(conjg(a)*a)
    if(.not.sq)norm=sqrt(norm)
  end function norm_c8d4

end module DA_Util_ml
