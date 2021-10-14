!#######################################################################
!
! Data domain decomposition tools.
!
! Preprocessor macro's
!
!  _MPI   : should be defined to enable MPI code, e.g. 'f90 -D_MPI ...'
!
!
! Manuals
!
! - MPI 3.0 manual:
!     http://mpi-forum.org/docs/mpi-3.0/mpi30-report.pdf
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


module CSO_Domains

  use CSO_Logging, only : csol, ncsol, csoPr, csoErr
#ifdef _MPI
  use MPI_F08, only : MPI_SUCCESS, MPI_Error_String
#endif

  implicit none


  ! --- in/out ------------------------

  private

  public  ::  T_CSO_Domains


  ! --- const --------------------------

  character(len=*), parameter   ::  mname = 'CSO_Domains'
  
  
  ! --- types -----------------------------------
  
  ! domain ranges from all process to facilitate decomposition
  type T_CSO_Domains
    ! rank:
    integer                ::  ndim
    ! global bounds:
    integer, allocatable   ::  glbo(:,:)  ! (ndim,0:csoc%npes-1)
    integer, allocatable   ::  gubo(:,:)  ! (ndim,0:csoc%npes-1)
    ! local bounds:
    integer, allocatable   ::  lbo(:,:)   ! (ndim,0:csoc%npes-1)
    integer, allocatable   ::  ubo(:,:)   ! (ndim,0:csoc%npes-1)
    ! offset:
    integer, allocatable   ::  off(:,:)   ! (ndim,0:csoc%npes-1)
    ! shape:
    integer, allocatable   ::  shp(:,:)   ! (ndim,0:csoc%npes-1)
    ! total number of local elements:
    integer, allocatable   ::  n(:)   ! (0:csoc%npes-1)
    !
  contains
    procedure   ::  Init         => CSO_Domains_Init
    procedure   ::  Done         => CSO_Domains_Done
    procedure   ::  Get          => CSO_Domains_Get
    procedure   ::  Find         => CSO_Domains_Find
  end type T_CSO_Domains


  ! --- var -------------------------------------
  
  ! mpi error code; used as argument for 'MPI_Error_String' to avoid
  ! warnings about same argument 'status' being used for both 'errorcode' and 'ierror':
  integer                           ::  errorcode


  
contains


  ! ********************************************************************
  ! ***
  ! *** Domains
  ! ***
  ! ********************************************************************


  !
  ! initialize domain composition using the global bounds of the local domain
  !

  subroutine CSO_Domains_Init( self, glbo, gubo, status )

    use CSO_Comm, only : csoc
    
    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Init'
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Domains), intent(out)     ::  self
    integer, intent(in)                   ::  glbo(:)
    integer, intent(in)                   ::  gubo(:)
    integer, intent(out)                  ::  status

    ! --- local ----------------------------------
    
    integer    ::  i
    
    ! --- begin ----------------------------------
    
    ! count:
    self%ndim = size(glbo)
    ! check ...
    if ( size(gubo) /= self%ndim ) then
      write (csol,'("size of argument gubo (",i0,") does not match with size of glbo (",i0,")")') &
                      size(gubo), self%ndim; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( self%glbo(self%ndim,0:csoc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%gubo(self%ndim,0:csoc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%off(self%ndim,0:csoc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%lbo(self%ndim,0:csoc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ubo(self%ndim,0:csoc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%shp(self%ndim,0:csoc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%n(0:csoc%npes-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
#ifdef _MPI
    ! exchange global bounds:
    call csoc%AllGather( glbo, self%glbo, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    call csoc%AllGather( gubo, self%gubo, status )
    IF_MPI_NOT_OK_RETURN(status=1)
#else
    ! copy:
    self%glbo(:,0) = glbo
    self%gubo(:,0) = gubo
#endif
    
    ! local offset in global space:
    self%off = self%glbo - 1
    
    ! local shapes:
    self%shp = self%gubo - self%glbo + 1
    ! trap undefined per process:
    do i = 0, csoc%npes-1
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

  end subroutine CSO_Domains_Init


  ! ***
  
  
  subroutine CSO_Domains_Done( self, status )
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Domains), intent(inout)      ::  self
    integer, intent(out)                 ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%glbo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%gubo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%lbo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%ubo, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%off, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%shp, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%n  , stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine CSO_Domains_Done


  ! ***
  
  
  subroutine CSO_Domains_Get( self, status, shp, off, glbo, gubo )

    use CSO_Comm, only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Domains), intent(in)      ::  self
    integer, intent(out)                  ::  status

    integer, intent(out), optional        ::  shp(:)
    integer, intent(out), optional        ::  off(:)
    integer, intent(out), optional        ::  glbo(:), gubo(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/CSO_Domains_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! return local shape?
    if ( present(shp) ) then
      ! check ..
      if ( size(shp) /= self%ndim ) then
        write (csol,'("argument shp has size ",i0," while dim is ",i0)') size(shp), self%ndim; call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      shp = self%shp(:,csoc%id)
    end if
    
    ! return local offset?
    if ( present(off) ) then
      ! check ..
      if ( size(off) /= self%ndim ) then
        write (csol,'("argument off has size ",i0," while dim is ",i0)') size(off), self%ndim; call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      off = self%off(:,csoc%id)
    end if
    
    ! return global lower bounds?
    if ( present(glbo) ) then
      ! check ..
      if ( size(glbo) /= self%ndim ) then
        write (csol,'("argument glbo has size ",i0," while dim is ",i0)') size(off), self%ndim; call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      glbo = self%glbo(:,csoc%id)
    end if
    
    ! return global upper bounds?
    if ( present(gubo) ) then
      ! check ..
      if ( size(gubo) /= self%ndim ) then
        write (csol,'("argument gubo has size ",i0," while dim is ",i0)') size(off), self%ndim; call csoErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      gubo = self%gubo(:,csoc%id)
    end if
    
    ! ok
    status = 0

  end subroutine CSO_Domains_Get
  

  ! ***

  
  ! index of domain holding cell with provided global indices
  
  subroutine CSO_Domains_Find( self, ind, iproc, status, locind )
  
    use CSO_Comm, only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_CSO_Domains), intent(in)            ::  self
    integer, intent(in)                     ::  ind(:)  ! global indices
    integer, intent(out)                    ::  iproc   ! 0:csoc%npes-1
    integer, intent(out)                    ::  status

    integer, intent(out), optional          ::  locind(:)  ! local indices

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Domains_Find'
    
    ! --- local ----------------------------------
    
    integer     ::  k
    
    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(ind) /= self%ndim ) then
      write (csol,'("argument ind has size ",i0," while ndim = ",i0)') size(ind), self%ndim; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! dummy result:
    iproc = -999
    ! loop:
    do k = 0, csoc%npes-1
      ! compare:
      if ( all(self%glbo(:,k) <= ind) .and. all(ind <= self%gubo(:,k)) ) then
        ! found !
        iproc = k
        ! set output?
        if ( present(locind) ) then
          ! local indices are global minus offset:
          locind = ind - self%off(:,k)
        end if
        ! leave:
        exit
      end if  ! location on proc domain
    end do  ! procs
    ! check ...
    if ( iproc < 0 ) then
      write (csol,*) 'could not find domain holding global indices : ', ind; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine CSO_Domains_Find


end module CSO_Domains
