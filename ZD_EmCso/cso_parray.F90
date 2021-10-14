!###############################################################################
!
! NAME
!   CSO_PArray - pointer array tools
!
! PROCEDURES
!
!   call CSO_Parray_Reshape( x, n )
!     <type>, pointer                   ::  x(:)
!     integer, intent(in)               ::  n
!     integer, intent(out)              ::  status
!
!     Reshape size of pointer array with n elements, copy current values.
!
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

module CSO_PArray

  use CSO_Logging, only : csol, csoPr, csoErr

  implicit none

  ! --- in/out -----------------------------

  private

  public  ::  CSO_Parray_Init
  public  ::  CSO_Parray_Done
  public  ::  CSO_Parray_Reshape


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'CSO_PArray'

  
  ! --- interfaces -------------------------------------

  interface CSO_Parray_Init
    module procedure CSO_Parray_Init_i1
    module procedure CSO_Parray_Init_r1
    module procedure CSO_Parray_Init_r2
  end interface

  interface CSO_Parray_Done
    module procedure CSO_Parray_Done_i1
    module procedure CSO_Parray_Done_r1
    module procedure CSO_Parray_Done_r2
  end interface

  interface CSO_Parray_Reshape
    module procedure CSO_Parray_Reshape_i1
    module procedure CSO_Parray_Reshape_r1
    module procedure CSO_Parray_Reshape_r2
  end interface


contains


  !**********************************************************************
  !***
  !*** i1
  !***
  !**********************************************************************


  subroutine CSO_Parray_Init_i1( x, status )

    ! --- in/out ----------------------------

    integer, pointer                  ::  x(:)
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Init_i1'

    ! --- local -----------------------------

    ! --- begin -----------------------------
    
    ! not associated:
    nullify( x )
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Init_i1


  ! ***


  subroutine CSO_Parray_Done_i1( x, status )

    ! --- in/out ----------------------------

    integer, pointer                  ::  x(:)
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Done_i1'

    ! --- local -----------------------------

    ! --- begin -----------------------------
    
    ! associated?
    if ( associated(x) ) then
      ! clear:
      deallocate( x, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! safety ...
    nullify( x )
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Done_i1


  ! ***


  ! Increase size with n elements,
  ! copy current values (if any).
  
  subroutine CSO_Parray_Reshape_i1( x, n, status, source )

    ! --- in/out ----------------------------

    integer, pointer                  ::  x(:)
    integer, intent(in)               ::  n
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  source

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Reshape_i1'

    ! --- local -----------------------------

    integer                 ::  nx
    integer                 ::  nc
    integer, pointer        ::  new_x(:)

    ! --- begin -----------------------------
    
    ! not associated yet?
    if ( .not. associated(x) ) then
      ! new storage (cannot use optional argument source):
      allocate( x(n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill values?
      if ( present(source) ) x = source
    else
      ! current size:
      nx = size(x)
      ! different?
      if ( nx /= n ) then
        ! new storage (cannot use optional argument source):
        allocate( new_x(n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill values?
        if ( present(source) ) x = source
        ! copy current:
        nc = min(nx,n)
        new_x(1:nc) = x(1:nc)
        ! clear current:
        deallocate( x, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! re-assign:
        x => new_x
      end if ! Reshape
    end if ! new
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Reshape_i1


  !**********************************************************************
  !***
  !*** r1
  !***
  !**********************************************************************


  subroutine CSO_Parray_Init_r1( x, status )

    ! --- in/out ----------------------------

    real, pointer                     ::  x(:)
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Init_r1'

    ! --- local -----------------------------

    ! --- begin -----------------------------
    
    ! not associated:
    nullify( x )
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Init_r1


  ! ***


  subroutine CSO_Parray_Done_r1( x, status )

    ! --- in/out ----------------------------

    real, pointer                     ::  x(:)
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Done_r1'

    ! --- local -----------------------------

    ! --- begin -----------------------------
    
    ! associated?
    if ( associated(x) ) then
      ! clear:
      deallocate( x, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! safety ...
    nullify( x )
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Done_r1


  ! ***


  ! Increase size with n elements,
  ! copy current values (if any).

  subroutine CSO_Parray_Reshape_r1( x, n, status, source )

    ! --- in/out ----------------------------

    real, pointer                     ::  x(:)
    integer, intent(in)               ::  n
    integer, intent(out)              ::  status
    real, intent(in), optional        ::  source

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Reshape_r1'

    ! --- local -----------------------------

    integer                 ::  nx
    integer                 ::  nc
    real, pointer           ::  new_x(:)

    ! --- begin -----------------------------
    
    ! not associated yet?
    if ( .not. associated(x) ) then
      ! new storage (cannot use optional argument source):
      allocate( x(n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill values?
      if ( present(source) ) x = source
    else
      ! current size:
      nx = size(x)
      ! different?
      if ( nx /= n ) then
        ! new storage (cannot use optional argument source):
        allocate( new_x(n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill values?
        if ( present(source) ) x = source
        ! copy current:
        nc = min(nx,n)
        new_x(1:nc) = x(1:nc)
        ! clear current:
        deallocate( x, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! re-assign:
        x => new_x
      end if ! Reshape
    end if ! new
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Reshape_r1


  !**********************************************************************
  !***
  !*** r1
  !***
  !**********************************************************************


  subroutine CSO_Parray_Init_r2( x, status )

    ! --- in/out ----------------------------

    real, pointer                     ::  x(:,:)
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Init_r2'

    ! --- local -----------------------------

    ! --- begin -----------------------------
    
    ! not associated:
    nullify( x )
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Init_r2


  ! ***


  subroutine CSO_Parray_Done_r2( x, status )

    ! --- in/out ----------------------------

    real, pointer                     ::  x(:,:)
    integer, intent(out)              ::  status

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Done_r2'

    ! --- local -----------------------------

    ! --- begin -----------------------------
    
    ! associated?
    if ( associated(x) ) then
      ! clear:
      deallocate( x, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! safety ...
    nullify( x )
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Done_r2


  ! ***


  ! Increase size with n elements,
  ! copy current values (if any).

  subroutine CSO_Parray_Reshape_r2( x, n, status, source )

    ! --- in/out ----------------------------

    real, pointer                     ::  x(:,:)
    integer, intent(in)               ::  n(2)
    integer, intent(out)              ::  status
    real, intent(in), optional        ::  source

    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Parray_Reshape_r21'

    ! --- local -----------------------------

    integer                 ::  nx(2)
    integer                 ::  nc(2)
    real, pointer           ::  new_x(:,:)

    ! --- begin -----------------------------
    
    ! not associated yet?
    if ( .not. associated(x) ) then
      ! new storage (cannot use optional argument source):
      allocate( x(n(1),n(2)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill values?
      if ( present(source) ) x = source
    else
      ! current size:
      nx = shape(x)
      ! different?
      if ( any( nx /= n ) ) then
        ! new storage (cannot use optional argument source):
        allocate( new_x(n(1),n(2)), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill values?
        if ( present(source) ) x = source
        ! copy current:
        nc = min(nx,n)
        new_x(1:nc(1),1:nc(2)) = x(1:nc(1),1:nc(2))
        ! clear current:
        deallocate( x, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! re-assign:
        x => new_x
      end if ! Reshape
    end if ! new
    
    ! ok
    status = 0
    
  end subroutine CSO_Parray_Reshape_r2


end module CSO_Parray

