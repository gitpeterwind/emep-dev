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


module CSO_Domains

  use CSO_Logging, only : csol, ncsol, csoPr, csoErr
#ifdef _MPI
  use MPI_F08, only : MPI_SUCCESS, MPI_Error_String
#endif

  implicit none


  ! --- in/out ------------------------

  private

  public  ::  T_CSO_Domains
  public  ::  T_CSO_Selection1D


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
  
  
  ! * info on distributed selections;
  !   for example used to read global array on root,
  !   and scatter selections to other pe
  
  type T_CSO_Selection1D
    ! number of global values:
    integer                         ::  nglb
    ! number of locally selected values:
    integer                         ::  nsel
    ! selected indices in range 1,..,nglb:
    integer, allocatable            ::  iglb_sel(:)   ! (nsel)
    !
    ! scatter input from root?
    ! then root will collected info from the other domins in "*_all" and "*_tot" arrays:
    logical                         ::  with_root
    ! number of selected elements per pe
    ! (needed on root only)
    integer, allocatable            ::  nsel_all(:)   ! (npes)
    ! total number of selected elements
    ! (sum over ns)
    integer                         ::  nsel_tot
    ! storage for all selections
    integer, allocatable            ::  iglb_tot(:)   ! (ns_tot), indices in 1,..,nglb
    !
    ! number of unique pixels used at some domain,
    ! equal or less than nglb in case some pixels are outside domain:
    integer                         ::  nout
    ! mapping from 1:nglb to 1:nout,
    ! for example used to write glb_lon/etc arrays
    integer, allocatable            ::  iout_glb(:)   ! (nglb)
    ! mapping from 1:nsel_tot to 1:nout:
    integer, allocatable            ::  iout_tot(:)   ! (nsel_tot)
    !
  contains
    procedure   ::  Init        =>  CSO_Selection1D_Init
    procedure   ::  Done        =>  CSO_Selection1D_Done
    procedure   ::                  CSO_Selection1D_Copy_i1_1d
    procedure   ::                  CSO_Selection1D_Copy_r_1d
    procedure   ::                  CSO_Selection1D_Copy_c1_2d
    procedure   ::                  CSO_Selection1D_Copy_r_2d
    procedure   ::                  CSO_Selection1D_Copy_r_3d
    generic     ::  Copy        =>  CSO_Selection1D_Copy_i1_1d, &
                                    CSO_Selection1D_Copy_r_1d, &
                                    CSO_Selection1D_Copy_c1_2d, &
                                    CSO_Selection1D_Copy_r_2d, &
                                    CSO_Selection1D_Copy_r_3d
    procedure   ::                  CSO_Selection1D_ScatterV_i1_1d
    procedure   ::                  CSO_Selection1D_ScatterV_r_1d
    procedure   ::                  CSO_Selection1D_ScatterV_c1_2d
    procedure   ::                  CSO_Selection1D_ScatterV_r_2d
    procedure   ::                  CSO_Selection1D_ScatterV_r_3d
    generic     ::  ScatterV    =>  CSO_Selection1D_ScatterV_i1_1d, &
                                    CSO_Selection1D_ScatterV_r_1d, &
                                    CSO_Selection1D_ScatterV_c1_2d, &
                                    CSO_Selection1D_ScatterV_r_2d, &
                                    CSO_Selection1D_ScatterV_r_3d
  end type T_CSO_Selection1D


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


  ! ====================================================================
  ! ===
  ! === Selection
  ! ===
  ! ====================================================================


  subroutine CSO_Selection1D_Init( self, nglb, nsel, iglb_sel, with_root, status )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(out)  ::  self
    integer, intent(in)                    ::  nglb
    integer, intent(in)                    ::  nsel          ! could be 0
    integer, intent(in)                    ::  iglb_sel(:)   ! (nsel)
    logical                                ::  with_root
    integer, intent(out)                   ::  status
    
    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_Init'

    ! --- local ----------------------------------
    
    logical, allocatable    ::  used_glb(:)    ! (nglb)
    integer                 ::  k
    integer                 ::  iglb
    integer                 ::  iout

    ! --- begin ----------------------------------
    
    ! store:
    self%nglb = nglb
    ! check ...
    if ( self%nglb <= 0 ) then
      write (csol,'("number of global elements should be > 0, found: ",i0)') self%nglb; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! store:
    self%nsel = nsel
    ! check ...
    if ( self%nsel < 0 ) then
      write (csol,'("number of selected elements should be >= 0, found: ",i0)') self%nsel; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! any selected?
    if ( self%nsel > 0 ) then
      ! check ..
      if ( size(iglb_sel) < self%nsel ) then
        write (csol,'("size of array with selected indices (",i0,") should be at least ",i0)') &
                           size(iglb_sel), self%nsel; call csoErr
        TRACEBACK; status=1; return
      end if
      ! check ..
      if ( any(iglb_sel < 1) .or. any( iglb_sel > self%nglb) ) then
        write (csol,'("expected selection indices in range 1,..",i0,", found range: ",i0,"..",i0)') &
                           self%nglb, minval(iglb_sel), maxval(iglb_sel); call csoErr
        TRACEBACK; status=1; return
      end if
      ! storage for copy:
      allocate( self%iglb_sel(self%nsel), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy:
      self%iglb_sel = iglb_sel(1:self%nsel)
    else
      ! dummy:
      allocate( self%iglb_sel(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if  ! nsel > 0
    
    ! copy flag:
    self%with_root = with_root
    
    ! setup scatter?
    if ( self%with_root ) then

      ! storage:
      if ( csoc%root ) then
        allocate( self%nsel_all(csoc%npes), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        allocate( self%nsel_all(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! collect on root:
      call csoc%Gather( self%nsel, self%nsel_all, status )
      IF_NOT_OK_RETURN(status=1)

      ! total number of selected elements:
      self%nsel_tot = sum(self%nsel_all)
      ! broadcast counter:
      call csoc%BCast( csoc%root_id, self%nsel_tot, status )
      IF_NOT_OK_RETURN(status=1)

      ! any selections?
      if ( self%nsel_tot > 0 ) then

        ! storage for all selection indices in single array:
        if ( csoc%root ) then
          allocate( self%iglb_tot(self%nsel_tot), stat=status )
          IF_NOT_OK_RETURN(status=1)
        else
          allocate( self%iglb_tot(1), stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if
        ! collect on root:
        call csoc%GatherV( iglb_sel, self%iglb_tot, status, nloc=self%nsel )
        IF_NOT_OK_RETURN(status=1)

        ! on root, info on output
        if ( csoc%root ) then

          ! storage for flag per global pixel, disable by default:
          allocate( used_glb(self%nglb), source=.false., stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop over all locally used pixesl:
          do k = 1, self%nsel_tot
            ! global pixel index:
            iglb = self%iglb_tot(k)
            ! enable:
            used_glb(iglb) = .true.
          end do

          ! count number of global pixels somewhere used:
          self%nout = count( used_glb )

          ! storage for mapping from 1:nglb to 1:nout
          allocate( self%iout_glb(self%nglb), source=-999, stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! init counter:
          iout = 0
          ! loop over global pixels:
          do iglb = 1, self%nglb
            ! in use?
            if ( used_glb(iglb) ) then
              ! next:
              iout = iout + 1
              ! store mapping:
              self%iout_glb(iglb) = iout
            end if
          end do ! iglb

          ! clear:
          deallocate( used_glb, stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! storage for mapping from "tot" array to "out":
          allocate( self%iout_tot(self%nsel_tot), stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop over all locally used pixesl:
          do k = 1, self%nsel_tot
            ! global pixel index:
            iglb = self%iglb_tot(k)
            ! mapping:
            self%iout_tot(k) = self%iout_glb(iglb)
          end do ! k

        else
          ! dummy ...
          allocate( self%iout_glb(1), source=-999, stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! dummy ...
          allocate( self%iout_tot(1), source=-999, stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if ! root

        ! broadcast counter:
        call csoc%BCast( csoc%root_id, self%nout, status )
        IF_NOT_OK_RETURN(status=1)

      end if ! nsel_tot > 0
      
    end if ! with root
    
    ! ok
    status = 0

  end subroutine CSO_Selection1D_Init


  ! ***


  subroutine CSO_Selection1D_Done( self, status )

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(inout)    ::  self
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%iglb_sel, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! scatter from root?
    if ( self%with_root ) then

      ! clear:
      deallocate( self%nsel_all, stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! any selections?
      if ( self%nsel_tot > 0 ) then
        ! clear:
        deallocate( self%iglb_tot, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! clear:
        deallocate( self%iout_glb, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( self%iout_tot, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if  ! nsel_tot > 0
      
    end if ! with root

    ! ok
    status = 0

  end subroutine CSO_Selection1D_Done


  ! ***

  !
  ! Copy elements from global array.
  !

  subroutine CSO_Selection1D_Copy_i1_1d( self, glb, loc, status )

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    integer(1), intent(in)                   ::  glb(:)  ! (nglb)
    integer(1), intent(out)                  ::  loc(:)  ! (nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_Copy_i1_1d'

    ! --- local ----------------------------------
    
    integer     ::  isel

    ! --- begin ----------------------------------
    
    ! check ...
    if ( any( shape(glb) /= (/self%nglb/) ) ) then
      write (csol,'("shape of glb (",i0,") while expected (",i0,")")') &
                       shape(glb), self%nglb; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! any selected?
    if ( self%nsel > 0 ) then
      ! check ...
      if ( any( shape(loc) /= (/self%nsel/) ) ) then
        write (csol,'("shape of loc (",i0,") while expected (",i0,")")') &
                         shape(glb), self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! loop over target elements:
      do isel = 1, self%nsel
        ! copy:
        loc(isel) = glb(self%iglb_sel(isel))
      end do ! i
    end if ! nsel > 0

    ! ok
    status = 0

  end subroutine CSO_Selection1D_Copy_i1_1d
  
  ! *

  subroutine CSO_Selection1D_Copy_r_1d( self, glb, loc, status )

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    real, intent(in)                         ::  glb(:)  ! (nglb)
    real, intent(out)                        ::  loc(:)  ! (nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_Copy_r_1d'

    ! --- local ----------------------------------
    
    integer     ::  isel

    ! --- begin ----------------------------------
    
    ! check ...
    if ( any( shape(glb) /= (/self%nglb/) ) ) then
      write (csol,'("shape of glb (",i0,") while expected (",i0,")")') &
                       shape(glb), self%nglb; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! any selected?
    if ( self%nsel > 0 ) then
      ! check ...
      if ( any( shape(loc) /= (/self%nsel/) ) ) then
        write (csol,'("shape of loc (",i0,") while expected (",i0,")")') &
                         shape(glb), self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! loop over target elements:
      do isel = 1, self%nsel
        ! copy:
        loc(isel) = glb(self%iglb_sel(isel))
      end do ! i
    end if ! nsel > 0

    ! ok
    status = 0

  end subroutine CSO_Selection1D_Copy_r_1d
  
  ! *

  subroutine CSO_Selection1D_Copy_c1_2d( self, glb, loc, status )

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    character(len=1), intent(in)             ::  glb(:,:)  ! (n1,nglb)
    character(len=1), intent(out)            ::  loc(:,:)  ! (n1,nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_Copy_c1_2d'

    ! --- local ----------------------------------
    
    integer     ::  n1
    integer     ::  isel

    ! --- begin ----------------------------------
    
    ! shape:
    n1 = size(glb,1)
    
    ! check ...
    if ( any( shape(glb) /= (/n1,self%nglb/) ) ) then
      write (csol,'("shape of glb (",i0,",",i0,") while expected (",i0,",",i0,")")') &
                       shape(glb), n1, self%nglb; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! any selected?
    if ( self%nsel > 0 ) then
      ! check ...
      if ( any( shape(loc) /= (/n1,self%nsel/) ) ) then
        write (csol,'("shape of loc (",i0,",",i0,") while expected (",i0,",",i0,")")') &
                         shape(glb), n1, self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! loop over target elements:
      do isel = 1, self%nsel
        ! copy:
        loc(:,isel) = glb(:,self%iglb_sel(isel))
      end do ! i
    end if ! nsel > 0

    ! ok
    status = 0

  end subroutine CSO_Selection1D_Copy_c1_2d
  
  ! *

  subroutine CSO_Selection1D_Copy_r_2d( self, glb, loc, status )

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    real, intent(in)                         ::  glb(:,:)  ! (n1,nglb)
    real, intent(out)                        ::  loc(:,:)  ! (n1,nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_Copy_r_2d'

    ! --- local ----------------------------------
    
    integer     ::  n1
    integer     ::  isel

    ! --- begin ----------------------------------
    
    ! shape:
    n1 = size(glb,1)
    
    ! check ...
    if ( any( shape(glb) /= (/n1,self%nglb/) ) ) then
      write (csol,'("shape of glb (",i0,",",i0,") while expected (",i0,",",i0,")")') &
                       shape(glb), n1, self%nglb; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! any selected?
    if ( self%nsel > 0 ) then
      ! check ...
      if ( any( shape(loc) /= (/n1,self%nsel/) ) ) then
        write (csol,'("shape of loc (",i0,",",i0,") while expected (",i0,",",i0,")")') &
                         shape(glb), n1, self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! loop over target elements:
      do isel = 1, self%nsel
        ! copy:
        loc(:,isel) = glb(:,self%iglb_sel(isel))
      end do ! i
    end if ! nsel > 0

    ! ok
    status = 0

  end subroutine CSO_Selection1D_Copy_r_2d
  
  ! *

  subroutine CSO_Selection1D_Copy_r_3d( self, glb, loc, status )

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    real, intent(in)                         ::  glb(:,:,:)  ! (n1,n2,nglb)
    real, intent(out)                        ::  loc(:,:,:)  ! (n1,n2,nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_Copy_r_3d'

    ! --- local ----------------------------------
    
    integer     ::  n1, n2
    integer     ::  isel

    ! --- begin ----------------------------------
    
    ! shape:
    n1 = size(glb,1)
    n2 = size(glb,2)
    
    ! check ...
    if ( any( shape(glb) /= (/n1,n2,self%nglb/) ) ) then
      write (csol,'("shape of glb (",i0,2(",",i0),") while expected (",i0,2(",",i0),")")') &
                       shape(glb), n1, n2, self%nglb; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! any selected?
    if ( self%nsel > 0 ) then
      ! check ...
      if ( any( shape(loc) /= (/n1,n2,self%nsel/) ) ) then
        write (csol,'("shape of loc (",i0,2(",",i0),") while expected (",i0,2(",",i0),")")') &
                         shape(glb), n1, n2, self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! loop over target elements:
      do isel = 1, self%nsel
        ! copy:
        loc(:,:,isel) = glb(:,:,self%iglb_sel(isel))
      end do ! i
    end if ! nsel > 0

    ! ok
    status = 0

  end subroutine CSO_Selection1D_Copy_r_3d


  ! ***

  !
  ! Scatter selections from global array on root
  !

  subroutine CSO_Selection1D_ScatterV_i1_1d( self, glb, loc, status )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    integer(1), intent(in)                   ::  glb(:)  ! (nglb)
    integer(1), intent(out)                  ::  loc(:)  ! (nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_ScatterV_i1_1d'

    ! --- local ----------------------------------
    
    integer                     ::  itot
    integer(1), allocatable     ::  sendbuf(:)

    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. self%with_root ) then
      write (csol,'("could not scatter if info is not collected on root")'); call csoPr
      TRACEBACK; status=1; return
    end if
    
    ! data to be scattered:
    if ( csoc%root ) then
      ! check ..
      if ( any( shape(glb) /= (/self%nglb/) ) ) then
        write (csol,'("shape of glb (",i0,") while expected (",i0,")")') &
                         shape(glb), self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      allocate( sendbuf(self%nsel_tot), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill:
      do itot = 1, self%nsel_tot
        sendbuf(itot) = glb(self%iglb_tot(itot))
      end do
    else
      ! dummy:
      allocate( sendbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! scatter:
    call csoc%ScatterV( sendbuf, loc, status, nloc=self%nsel )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine CSO_Selection1D_ScatterV_i1_1d
  
  ! *

  subroutine CSO_Selection1D_ScatterV_r_1d( self, glb, loc, status )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    real, intent(in)                         ::  glb(:)  ! (nglb)
    real, intent(out)                        ::  loc(:)  ! (nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_ScatterV_r_1d'

    ! --- local ----------------------------------
    
    integer               ::  itot
    real, allocatable     ::  sendbuf(:)

    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. self%with_root ) then
      write (csol,'("could not scatter if info is not collected on root")'); call csoPr
      TRACEBACK; status=1; return
    end if
    
    ! data to be scattered:
    if ( csoc%root ) then
      ! check ..
      if ( any( shape(glb) /= (/self%nglb/) ) ) then
        write (csol,'("shape of glb (",i0,") while expected (",i0,")")') &
                         shape(glb), self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      allocate( sendbuf(self%nsel_tot), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill:
      do itot = 1, self%nsel_tot
        sendbuf(itot) = glb(self%iglb_tot(itot))
      end do
    else
      ! dummy:
      allocate( sendbuf(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! scatter:
    call csoc%ScatterV( sendbuf, loc, status, nloc=self%nsel )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine CSO_Selection1D_ScatterV_r_1d
  
  ! *

  subroutine CSO_Selection1D_ScatterV_c1_2d( self, glb, loc, status )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    character(len=1), intent(in)             ::  glb(:,:)  ! (n1,nglb) on root
    character(len=1), intent(out)            ::  loc(:,:)  ! (n1,nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_ScatterV_c1_2d'

    ! --- local ----------------------------------
    
    integer                           ::  n1
    integer                           ::  itot
    character(len=1), allocatable     ::  sendbuf(:,:)

    ! --- begin ----------------------------------
    
    ! shape:
    n1 = size(glb,1)
    
    ! check ...
    if ( .not. self%with_root ) then
      write (csol,'("could not scatter if info is not collected on root")'); call csoPr
      TRACEBACK; status=1; return
    end if
    
    ! data to be scattered:
    if ( csoc%root ) then
      ! check ..
      if ( any( shape(glb) /= (/n1,self%nglb/) ) ) then
        write (csol,'("shape of glb (",i0,",",i0,") while expected (",i0,",",i0,")")') &
                         shape(glb), n1, self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      allocate( sendbuf(n1,self%nsel_tot), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill:
      do itot = 1, self%nsel_tot
        sendbuf(:,itot) = glb(:,self%iglb_tot(itot))
      end do
    else
      ! dummy:
      allocate( sendbuf(1,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! scatter:
    call csoc%ScatterV( sendbuf, loc, status, nloc=self%nsel )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine CSO_Selection1D_ScatterV_c1_2d
  
  ! *

  subroutine CSO_Selection1D_ScatterV_r_2d( self, glb, loc, status )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    real, intent(in)                         ::  glb(:,:)  ! (n1,nglb) on root
    real, intent(out)                        ::  loc(:,:)  ! (n1,nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_ScatterV_r_2d'

    ! --- local ----------------------------------
    
    integer               ::  n1
    integer               ::  itot
    real, allocatable     ::  sendbuf(:,:)

    ! --- begin ----------------------------------
    
    ! shape:
    n1 = size(glb,1)
    
    ! check ...
    if ( .not. self%with_root ) then
      write (csol,'("could not scatter if info is not collected on root")'); call csoPr
      TRACEBACK; status=1; return
    end if
    
    ! data to be scattered:
    if ( csoc%root ) then
      ! check ..
      if ( any( shape(glb) /= (/n1,self%nglb/) ) ) then
        write (csol,'("shape of glb (",i0,",",i0,") while expected (",i0,",",i0,")")') &
                         shape(glb), n1, self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      allocate( sendbuf(n1,self%nsel_tot), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill:
      do itot = 1, self%nsel_tot
        sendbuf(:,itot) = glb(:,self%iglb_tot(itot))
      end do
    else
      ! dummy:
      allocate( sendbuf(1,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! scatter:
    call csoc%ScatterV( sendbuf, loc, status, nloc=self%nsel )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine CSO_Selection1D_ScatterV_r_2d
  
  ! *

  subroutine CSO_Selection1D_ScatterV_r_3d( self, glb, loc, status )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------

    class(T_CSO_Selection1D), intent(in)     ::  self
    real, intent(in)                         ::  glb(:,:,:)  ! (n1,n2,nglb) on root
    real, intent(out)                        ::  loc(:,:,:)  ! (n1,n2,nsel)
    integer, intent(out)                     ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/CSO_Selection1D_ScatterV_r_3d'

    ! --- local ----------------------------------
    
    integer               ::  n1, n2
    integer               ::  itot
    real, allocatable     ::  sendbuf(:,:,:)

    ! --- begin ----------------------------------
    
    ! shape:
    n1 = size(glb,1)
    n2 = size(glb,2)
    
    ! check ...
    if ( .not. self%with_root ) then
      write (csol,'("could not scatter if info is not collected on root")'); call csoPr
      TRACEBACK; status=1; return
    end if
    
    ! data to be scattered:
    if ( csoc%root ) then
      ! check ..
      if ( any( shape(glb) /= (/n1,n2,self%nglb/) ) ) then
        write (csol,'("shape of glb (",i0,2(",",i0),") while expected (",i0,2(",",i0),")")') &
                         shape(glb), n1, n2, self%nglb; call csoErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      allocate( sendbuf(n1,n2,self%nsel_tot), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill:
      do itot = 1, self%nsel_tot
        sendbuf(:,:,itot) = glb(:,:,self%iglb_tot(itot))
      end do
    else
      ! dummy:
      allocate( sendbuf(1,1,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! scatter:
    call csoc%ScatterV( sendbuf, loc, status, nloc=self%nsel )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine CSO_Selection1D_ScatterV_r_3d



end module CSO_Domains
