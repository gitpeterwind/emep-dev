!###############################################################################
!
! CSO_Profile - conversion of vertical profiles
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


!! testing ...
!module CSO_Logging
!  character(len=1024)   ::  csol
!contains
!  subroutine csoPr()
!    write (*,'(a)') trim(csol)
!  end subroutine csoPr
!  subroutine csoErr()
!    write (*,'("ERROR - ",a)') trim(csol)
!  end subroutine csoErr
!end module CSO_Logging


  
module CSO_Profile

  use CSO_Logging     , only : csol, csoPr, csoErr

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public    ::  T_ProfileMapping

  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Profile'
  
  
  ! --- types ------------------------------
  
  type T_ProfileMapping
    ! sizes:
    integer                   ::  nx
    integer                   ::  ny
    ! interval lengths:
    real, allocatable         ::  dx(:)  ! (nx)
    real, allocatable         ::  dy(:)  ! (ny)
    ! mapping:
    integer                   ::  nw
    integer, allocatable      ::  ii(:)  ! (nwmax)
    integer, allocatable      ::  jj(:)  ! (nwmax)
    real, allocatable         ::  ww(:)  ! (nwmax)
    !
  contains
    procedure ::  Init                   => ProfileMapping_Init
    procedure ::  Done                   => ProfileMapping_Done
    procedure ::  Setup                  => ProfileMapping_Setup
    procedure ::  Apply_Sum              => ProfileMapping_Apply_Sum
    procedure ::  Apply_Sum_Adj          => ProfileMapping_Apply_Sum_Adj
    procedure ::  Apply_WeightedSum      => ProfileMapping_Apply_WeightedSum
    procedure ::  Apply_WeightedSumAdj   => ProfileMapping_Apply_WeightedSumAdj
    procedure ::  Apply_WeightedMean     => ProfileMapping_Apply_WeightedMean
  end type T_ProfileMapping

  
contains


  !===============================================================================
  !
  ! The 'T_ProfileMapping' class is used to perform mappings
  ! from a profile f(1:nx) to g(1:ny).
  ! The coordinates x and y are defined as grid cells,
  ! with f and g constant in a cell.
  !
  ! As exaple, define the source profile:
  !
  !   ! source size:
  !   integer, parameter       ::  nx = 5
  !   ! interval boundaries:
  !   real                     ::  xx(0:nx) = (/ 0.0,   1.0,   2.0,   3.0,   4.0,   5.0 /)
  !   ! profile:
  !   real                     ::  f(nx)    = (/   100.0, 200.0, 300.0, 400.0, 500.0 /)
  !
  ! and target cells:
  !
  !   ! target size:
  !   integer, parameter       ::  ny = 3
  !   ! interval boundaries::
  !   real                     ::  yy(0:ny) = (/ 0.0, 1.5, 3.5, 5.0 /)
  !
  ! Define a mapping object using:
  !
  !   type(T_ProfileMapping)          ::  mapping
  !   integer                         ::  status
  !
  ! Initialize using the object using the input and output size:
  !
  !   call mapping%Init( nx, ny, status )
  !   if (status/=0) stop
  !
  ! Compute weights etc for mapping from xx to yy:
  !
  !   call mapping%Setup( xx, yy, status )
  !   if (status/=0) stop
  !
  ! Three routines are available to apply the weights.
  !
  ! * The weights can be used to compute partial sums at the target grid:
  !
  !     call mapping%Apply_Sum( f, g, status )
  !     if (status/=0) stop
  !
  !   In case f is in [kg], this will conserve the total mass:
  !
  !     g(j) = sum_i f(i) * w(i,j)
  !     [kg]         [kg]
  !
  !   with w(i,j) the faction of source interval i 
  !   that overlaps with target interval j.
  !
  ! * A weighted sum could be computed in case an integral is needed:
  !
  !     call mapping%Apply_WeightedSum( f, g, status )
  !     if (status/=0) stop
  !
  !   For example, if f is a mass mixing ratio in [(kg tracer)/(kg air)]
  !   and x and y are vertical airmass column density in [(kg air)/m2]
  !   the the result is a tracer column density:
  !
  !        g(j)     = sum_i  f(i)  *  dx(i) * w(i,j)
  !    [(kg tr)/m2]        [kg/kg]  [kg/m2]
  !
  !   with w(i,j) the faction of source interval i with length dx(i)
  !   that overlaps with target interval j.
  !   A coordinate x with this units can be computed from pressure [Pa]
  !   devided by the gravitiy constant g (~9.8 [m/s2]):
  !
  !           x         = pressure /    g
  !      [(kg air)/m2]      [Pa]   /  [m/s2]
  !
  ! * A weighted mean will also normalize with the target interval:
  !
  !     call mapping%Apply_WeightedMean( f, g, status )
  !     if (status/=0) stop
  !
  !   For example, if f is a mass mixing ratio in [(kg tracer)/(kg air)]
  !   and x and y are vertical airmass column density in [(kg air)/m2]
  !   the the result is a mass mixing ratio:
  !
  !        g(j)     = [ sum_i  f(i)  *  dx(i) * w(i,j) ] / dy(j)
  !      [kg/kg]             [kg/kg]  [kg/m2]             [kg/m2]
  !
  !   with w(i,j) the faction of source interval i with length dx(i)
  !   that overlaps with target interval j with length dy(j)
  !
  ! Clear the object using:
  !
  !   call mapping%Done( status )
  !   if (status/=0) stop
  !
  !===============================================================================


  !
  ! Initialize for mappings from xx(1:nx) to yy(1:ny).
  !
  
  subroutine ProfileMapping_Init( self, nx, ny, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ProfileMapping), intent(out)      ::  self
    integer, intent(in)                       ::  nx
    integer, intent(in)                       ::  ny
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ProfileMapping_Init'

    ! --- local ----------------------------------
    
    integer        ::    nwmax

    ! --- begin ----------------------------------
    
    ! store:
    self%nx = nx
    self%ny = ny
    
    ! storage:
    allocate( self%dx(self%nx), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%dy(self%ny), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! maximum number of mappings:
    nwmax = self%nx * self%ny
    ! storage:
    allocate( self%ii(nwmax), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj(nwmax), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ww(nwmax), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! no weights yet:
    self%nw = 0
    
    ! ok
    status = 0
    
  end subroutine ProfileMapping_Init
  
  ! *

  subroutine ProfileMapping_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ProfileMapping), intent(inout)    ::  self
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ProfileMapping_Done'

    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%dx, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%dy, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%ii, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%jj, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%ww, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine ProfileMapping_Done
  
  ! *
  
  !
  ! Compute weights for mapping 
  ! from profile x(1:nx) with boundaries xx(1:nx+1)
  !   to profile y(1:ny) with boundaries yy(1:ny+1)
  !

  subroutine ProfileMapping_Setup( self, xx, yy, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ProfileMapping), intent(inout)    ::  self
    real, intent(in)                          ::  xx(:)  ! (1:nx+1)
    real, intent(in)                          ::  yy(:)  ! (1:ny+1)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ProfileMapping_Setup'

    ! --- local ----------------------------------
    
    integer     ::  xdir, ydir
    integer     ::  i, j
    integer     ::  is, in
    integer     ::  js, jn
    integer     ::  i1, i2
    real        ::  y1, y2
    real        ::  w
    real        ::  wsum

    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(xx) /= self%nx+1 ) then
      write (csol,'("xx has size ",i0," instead of nx+1 for nx ",i0)') size(xx), self%nx; call csoErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    if ( size(yy) /= self%ny+1 ) then
      write (csol,'("yy has size ",i0," instead of ny+1 for ny ",i0)') size(yy), self%ny; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! source cell range and step in increasing direction:
    if ( xx(1) < xx(self%nx+1) ) then
      is   = 1
      in   = self%nx
      xdir = 1
    else
      is   = self%nx
      in   = 1
      xdir = -1
    end if
 
    ! target cell range and step in increasing direction:
    if ( yy(1) < yy(self%ny+1) ) then
      js   = 1
      jn   = self%ny
      ydir = 1
    else
      js   = self%ny
      jn   = 1
      ydir = -1
    end if
    
    ! interval lengths:
    self%dx = xdir * ( xx(2:self%nx+1) - xx(1:self%nx) )
    ! check ..
    if ( any(self%dx <= 0.0) ) then
      write (csol,'("found zero or negative interval length in x:")'); call csoErr
      do i = 1, self%nx
        write (csol,'("  cell ",i6," [",f16.6,",",f16.6,"] ",f16.6)') i, xx(i), xx(i+1), self%dx(i); call csoErr
      end do
      TRACEBACK; status=1; return
    end if
    
    ! interval lengths:
    self%dy = ydir * ( yy(2:self%ny+1) - yy(1:self%ny) )
    ! check ..
    if ( any(self%dy <= 0.0) ) then
      write (csol,'("found zero or negative interval length in y:")'); call csoErr
      do i = 1, self%ny
        write (csol,'("  cell ",i6," [",f16.6,",",f16.6,"] ",f16.6)') i, yy(i), yy(i+1), self%dy(i); call csoErr
      end do
      TRACEBACK; status=1; return
    end if
    
    ! reset weights:
    self%nw = 0
    ! no first source cell yet:
    i1 = -999
    ! loop over target elements:
    do j = js, jn, ydir
    
      ! target cell bounds:
      if ( ydir > 0 ) then
        y1 = yy(j)
        y2 = yy(j+1)
      else
        y1 = yy(j+1)
        y2 = yy(j)
      end if
      
      !! testing ...
      !print *, 'target cell ', j, ' range ', y1, y2
      
      ! search x cell holding y1 if not copied from i2:
      if ( i1 < 0 ) then
        ! loop over source cells until found:
        do i = is, in, xdir
          ! in cell bounds?
          if ( (xdir*xx(i) <= xdir*y1) .and. (xdir*y1 <= xdir*xx(i+1)) ) then
            i1 = i
            exit
          end if
        end do ! i
      end if ! search i1
      ! check ..
      if ( i1 < 0 ) then
        write (csol,'("could not find interval holding y1 = ",f16.6)') y1; call csoErr
        do i = 1, self%nx
          write (csol,'("  cell ",i6," [",f16.6,",",f16.6,"]")') i, xx(i), xx(i+1); call csoErr
        end do
        TRACEBACK; status=1; return
      end if
      
      ! search x cell holding y2:
      i2 = -999
      do i = i1, in, xdir
        !! testing ...
        !print *, '  search i2 in ', i, ';', xdir*xx(i), '<=', xdir*y2, '=', xdir*xx(i) <= xdir*y2, ';', &
        !                                    xdir*y2, '<=', xdir*xx(i+1), '=', xdir*y2 < xdir*xx(i+1)
        ! in cell bounds?
        if ( (xdir*xx(i) <= xdir*y2) .and. (xdir*y2 <= xdir*xx(i+1)) ) then
          i2 = i
          exit
        end if
      end do ! i
      ! check ..
      if ( i2 < 0 ) then
        write (csol,'("could not find interval holding y2 = ",f16.6)') y2; call csoErr
        do i = 1, self%nx
          write (csol,'("  cell ",i6," [",f16.6,",",f16.6,"]")') i, xx(i), xx(i+1); call csoErr
        end do
        TRACEBACK; status=1; return
      end if
      
      !! testing ...
      !print *, '  source range: ', i1, i2
      ! loop over source cells:
      do i = i1, i2, xdir
        ! fraction of source cell that overlaps with [y1,y2]:
        if ( xdir > 0 ) then
          w = ( min(xx(i+1),y2) - max(xx(i  ),y1) )/self%dx(i)
        else
          w = ( min(xx(i  ),y2) - max(xx(i+1),y1) )/self%dx(i)
        end if
        !! testing ..
        !print *, '    cell ', i, '[', xx(i), xx(i+1), ']', w
        ! increase counter:
        self%nw = self%nw + 1
        ! store:
        self%ii(self%nw) = i
        self%jj(self%nw) = j
        self%ww(self%nw) = w
      end do ! i
    
      ! next:
      i1 = i2
    end do ! j
    
    ! check ..
    wsum = sum(self%ww(1:self%nw))
    if ( abs(wsum - 1.0) <= 1.0e-4 ) then
      write (csol,'("sum of weights is ",es16.6," instead of 1.0")') wsum; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine ProfileMapping_Setup
  
  ! *
  
  !
  ! Compute partial sums:
  !
  !   g(j) = sum f(i) * w(i,j)
  !           i
  !
  ! with w(i,j) the faction of source interval i 
  ! that overlaps with target interval j.
  !

  subroutine ProfileMapping_Apply_Sum( self, f, g, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ProfileMapping), intent(inout)    ::  self
    real, intent(in)                          ::  f(:)  ! (1:nx)
    real, intent(out)                         ::  g(:)  ! (1:ny)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ProfileMapping_Apply_Sum'

    ! --- local ----------------------------------
    
    integer     ::  iw

    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(f) /= self%nx ) then
      write (csol,'("f has size ",i0," instead of nx ",i0)') size(f), self%nx; call csoErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    if ( size(g) /= self%ny ) then
      write (csol,'("g has size ",i0," instead of ny ",i0)') size(g), self%ny; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! init result:
    g = 0.0
    ! loop over mapping weights:
    do iw = 1, self%nw
      ! add contribution:
      g(self%jj(iw)) = g(self%jj(iw)) + f(self%ii(iw)) * self%ww(iw)
    end do ! iw
    
    ! ok
    status = 0
    
  end subroutine ProfileMapping_Apply_Sum
  
  ! *
  
  !
  ! Return adjoint of partial sums:
  !
  !   f(i) = sum g(j) * w(i,j)
  !           j
  !
  ! with w(i,j) the faction of source interval i 
  ! that overlaps with target interval j.
  !

  subroutine ProfileMapping_Apply_Sum_Adj( self, g, f, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ProfileMapping), intent(inout)    ::  self
    real, intent(in)                          ::  g(:)  ! (1:ny)
    real, intent(out)                         ::  f(:)  ! (1:nx)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ProfileMapping_Apply_Sum_Adj'

    ! --- local ----------------------------------
    
    integer     ::  iw

    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(f) /= self%nx ) then
      write (csol,'("f has size ",i0," instead of nx ",i0)') size(f), self%nx; call csoErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    if ( size(g) /= self%ny ) then
      write (csol,'("g has size ",i0," instead of ny ",i0)') size(g), self%ny; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! init result:
    f = 0.0
    ! loop over mapping weights:
    do iw = 1, self%nw
      ! add contribution:
      f(self%ii(iw)) = f(self%ii(iw)) + g(self%jj(iw)) * self%ww(iw)
    end do ! iw
    
    ! ok
    status = 0
    
  end subroutine ProfileMapping_Apply_Sum_Adj
  
  ! *
  
  !
  ! Compute partial sums weighted with interval lenths:
  !
  !   g(j) = sum f(i) * dx(i) * w(i,j)
  !           i
  !
  ! with w(i,j) the faction of source interval i with length dx(i)
  ! that overlaps with target interval j.
  !

  subroutine ProfileMapping_Apply_WeightedSum( self, f, g, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ProfileMapping), intent(inout)    ::  self
    real, intent(in)                          ::  f(:)  ! (1:nx)
    real, intent(out)                         ::  g(:)  ! (1:ny)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ProfileMapping_Apply_WeightedSum'

    ! --- local ----------------------------------
    
    integer     ::  iw

    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(f) /= self%nx ) then
      write (csol,'("f has size ",i0," instead of nx ",i0)') size(f), self%nx; call csoErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    if ( size(g) /= self%ny ) then
      write (csol,'("g has size ",i0," instead of ny ",i0)') size(g), self%ny; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! init result:
    g = 0.0
    ! loop over mapping weights:
    do iw = 1, self%nw
      ! add contribution:
      g(self%jj(iw)) = g(self%jj(iw)) + f(self%ii(iw)) * self%dx(self%ii(iw)) * self%ww(iw)
    end do ! iw
    
    ! ok
    status = 0
    
  end subroutine ProfileMapping_Apply_WeightedSum
  
  ! *
  
  !
  ! Adjoint of partial sums weighted with interval lenths:
  !
  !   forward:   g(j) = sum_i f(i) * dx(i) * w(i,j)
  !
  !   adjoint:   f(i) = sum_j g(j) * dx(i) * w(i,j)
  !
  ! with w(i,j) the faction of source interval i with length dx(i)
  ! that overlaps with target interval j.
  !

  subroutine ProfileMapping_Apply_WeightedSumAdj( self, g, f, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ProfileMapping), intent(inout)    ::  self
    real, intent(in)                          ::  g(:)  ! (1:ny)
    real, intent(out)                         ::  f(:)  ! (1:nx)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ProfileMapping_Apply_WeightedSumAdj'

    ! --- local ----------------------------------
    
    integer     ::  iw

    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(f) /= self%nx ) then
      write (csol,'("f has size ",i0," instead of nx ",i0)') size(f), self%nx; call csoErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    if ( size(g) /= self%ny ) then
      write (csol,'("g has size ",i0," instead of ny ",i0)') size(g), self%ny; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! init result:
    f = 0.0
    ! loop over mapping weights:
    do iw = 1, self%nw
      ! add contribution:
      f(self%ii(iw)) = f(self%ii(iw)) + g(self%jj(iw)) * self%dx(self%ii(iw)) * self%ww(iw)
    end do ! iw
    
    ! ok
    status = 0
    
  end subroutine ProfileMapping_Apply_WeightedSumAdj
  
  ! *
  
  !
  ! Compute partial sums normalized with interval lenths:
  !
  !   g(j) = [ sum f(i) * dx(i) * w(i,j) ] / dy(j)
  !             i
  !
  ! with w(i,j) the faction of source interval i with length dx(i)
  ! that overlaps with target interval j with length dy(j)
  !

  subroutine ProfileMapping_Apply_WeightedMean( self, f, g, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ProfileMapping), intent(inout)    ::  self
    real, intent(in)                          ::  f(:)  ! (1:nx)
    real, intent(out)                         ::  g(:)  ! (1:ny)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter   ::  rname = mname//'/ProfileMapping_Apply_WeightedMean'

    ! --- local ----------------------------------
    
    integer     ::  iw

    ! --- begin ----------------------------------
    
    ! check ..
    if ( size(f) /= self%nx ) then
      write (csol,'("f has size ",i0," instead of nx ",i0)') size(f), self%nx; call csoErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    if ( size(g) /= self%ny ) then
      write (csol,'("g has size ",i0," instead of ny ",i0)') size(g), self%ny; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! init result:
    g = 0.0
    ! loop over mapping weights:
    do iw = 1, self%nw
      ! add contribution:
      g(self%jj(iw)) = g(self%jj(iw)) + f(self%ii(iw)) * self%dx(self%ii(iw)) * self%ww(iw) / self%dy(self%jj(iw))
    end do ! iw
    
    ! ok
    status = 0
    
  end subroutine ProfileMapping_Apply_WeightedMean
                               
                               
end module CSO_Profile
