!###############################################################################
!
! Exchange pixels between domains
!
!###############################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!###############################################################################

module CSO_Exchange

  use CSO_Logging, only : csol, csoPr, csoErr

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private

  public  ::  T_Exchange
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Exchange'


  ! --- types --------------------------------
  
  ! Exchange of pixels with other domains
  type T_Exchange
    ! processor id of other pe:
    integer                           ::  other
    ! number of pixels to exchange:
    integer                           ::  nex
    ! local pixel indices involved in exchange:
    integer, allocatable              ::  ipix(:)
    ! temporary exchange buffer:
    real, pointer                     ::  sendbuf(:,:)  ! (nval,nex)
    real, pointer                     ::  recvbuf(:,:)  ! (nval,nex)
    !
  contains
    procedure :: Init            => Exchange_Init
    procedure :: Done            => Exchange_Done
    procedure :: Alloc           => Exchange_Alloc
    procedure :: Get             => Exchange_Get
    procedure :: AllocValues     => Exchange_AllocValues
  end type T_Exchange
  

  
contains


  ! ====================================================================
  ! ===
  ! === exchange info
  ! ===
  ! ====================================================================


  subroutine Exchange_Init( self, pid, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Exchange), intent(out)              ::  self
    integer, intent(in)                         ::  pid
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Exchange_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! exchange with pe:
    self%other = pid
    
    ! no content yet:
    self%nex = 0
    
    ! no buffer yet:
    nullify( self%sendbuf )
    nullify( self%recvbuf )
    
    ! ok
    status = 0
    
  end subroutine Exchange_Init


  ! ***


  subroutine Exchange_Done( self, status )

    ! --- in/out ---------------------------------
    
    class(T_Exchange), intent(inout)              ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Exchange_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! any content?
    if ( self%nex > 0 ) then
      ! clear:
      deallocate( self%ipix, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! nex > 0

    ! clear buffer?
    if ( associated(self%sendbuf) ) then
      deallocate( self%sendbuf, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( associated(self%recvbuf) ) then
      deallocate( self%recvbuf, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine Exchange_Done


  ! ***


  subroutine Exchange_Alloc( self, nex, status )

    ! --- in/out ---------------------------------
    
    class(T_Exchange), intent(inout)              ::  self
    integer, intent(in)                           ::  nex
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Exchange_Alloc'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( nex < 0 ) then
      write (csol,'("unsupported nex ",i0)') nex; call csoErr
      TRACEBACK; status=1; return
    end if
    ! store:
    self%nex = nex
    ! storage?
    if ( self%nex > 0 ) then
      allocate( self%ipix(self%nex), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! nex > 0
    
    ! ok
    status = 0
    
  end subroutine Exchange_Alloc


  ! ***


  subroutine Exchange_AllocValues( self, nval, status )

    ! --- in/out ---------------------------------
    
    class(T_Exchange), intent(inout)              ::  self
    integer, intent(in)                           ::  nval
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Exchange_AllocValues'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! storage:
    allocate( self%sendbuf(nval,self%nex), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%recvbuf(nval,self%nex), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Exchange_AllocValues


  ! ***


  subroutine Exchange_Get( self, status, nex, sendbuf, recvbuf )

    ! --- in/out ---------------------------------
    
    class(T_Exchange), intent(inout)              ::  self
    integer, intent(out)                          ::  status

    integer, intent(out), optional                ::  nex
    real, pointer, optional                       ::  sendbuf(:,:)
    real, pointer, optional                       ::  recvbuf(:,:)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Exchange_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! number of values to exchange:
    if ( present(nex) ) nex = self%nex
    
    ! pointer to buffer:
    if ( present(sendbuf) ) then
      if ( associated(self%sendbuf) ) then
        sendbuf => self%sendbuf
      else
        nullify( sendbuf )
      end if
    end if
    
    ! pointer to buffer:
    if ( present(recvbuf) ) then
      if ( associated(self%recvbuf) ) then
        recvbuf => self%recvbuf
      else
        nullify( recvbuf )
      end if
    end if

    ! ok
    status = 0
    
  end subroutine Exchange_Get



end module CSO_Exchange

