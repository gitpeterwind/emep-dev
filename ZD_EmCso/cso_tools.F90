!###############################################################################
!
! CSO_Tools
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


  
module CSO_Tools

  use CSO_Logging     , only : csol, csoPr, csoErr

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public    ::  TriangleArea
  public    ::  LonLatTriangleArea
  public    ::  LonLatRectangleArea  
  public    ::  GetPolygonPoints
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Tools'
  
  ! earth radius used in area compuations:
  real, parameter         ::  EarthRadius = 6.371e6     ! m
  
  ! something with circles ...
  real, parameter         ::  pi = 3.14159265358979
  ! conversion from degrees to radians:
  real, parameter         ::  deg2rad = pi/180.0     ! rad/deg

  
contains



  !
  ! Use Heron's formula for area of triangle given corners:
  !
  !   https://en.wikipedia.org/wiki/Heron%27s_formula
  !
  ! Formula below is the standard form:
  ! - sides of triangle: a, b, c
  ! - semi-perimeter:   s = 0.5 * ( a + b + c )
  ! - area:  A^2 =  s * (s-a) * (s-b) * (s-c) 
  ! 
  ! Alternative is to use the "numerical stable" form?
  ! - order: a >= b >= c
  ! - area: (4*A)^2 = (a+(b+c)) * (c - (a-b)) * (c+(a-b)) * (a+(b-c))
  !
  ! However, for problematic "needle" shaped triangles also the second
  ! formula could give slightly negative area^2. 
  ! Because it requires less if-statements the starndard form is used,
  ! and negative area's are truncated to zero.
  !
  
  subroutine TriangleArea( xx, yy, area, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(3), yy(3)
    real, intent(out)                    ::  area
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/TriangleArea'
    
    ! --- local ----------------------------------
    
    real                    ::  a, b, c
    real                    ::  s
    real                    ::  area2

    !real                    ::  qa, qb, qc

    ! --- begin ----------------------------------
    
    ! side lengths:
    a = sqrt( (xx(2) - xx(1))**2 + (yy(2) - yy(1))**2 )
    b = sqrt( (xx(3) - xx(2))**2 + (yy(3) - yy(2))**2 )
    c = sqrt( (xx(1) - xx(3))**2 + (yy(1) - yy(3))**2 )
    ! semi-perimeter:
    s = 0.5 * ( a + b + c )
    ! squared area:
    area2 = s * (s-a) * (s-b) * (s-c)
    ! check ...
    if ( area2 < 0.0 ) then
      ! truncate ...
      area = 0.0
      !! testing ...
      !write (csol,*) 'could not calculate area of triangle:'; call csoErr
      !write (csol,*) '  xx     = ', xx; call csoErr
      !write (csol,*) '  yy     = ', yy; call csoErr
      !write (csol,*) '  area^2 = ', area2; call csoErr
      !TRACEBACK; status=1; return
    else
      ! area:
      area = sqrt( area2 )
    end if
    
    !! sort such that qa >= qb >= qc:
    !!~ initial copy: qa, qb, qc
    !qa = a
    !qb = b
    !qc = c
    !! swap first pair if needed:
    !if ( qa < qb ) then
    !  qa = b
    !  qb = a
    !end if
    !! now:  qa >= qb, qc
    !! need to swap second pair?
    !if ( qb < qc ) then
    !  qc = qb
    !  qb = c
    !  ! now:  qa, qb >= qc
    !  ! need to swap firt pair again?
    !  if ( qa < qb ) then
    !    qb = qa
    !    qa = c
    !  end if
    !  ! now: qa >= qb >= qc
    !end if
    !! compute area:
    !TriangleArea = 0.25*sqrt( (qa+(qb+qc)) * (qc - (qa-qb)) * (qc+(qa-qb)) * (qa+(qb-qc)) )

    ! ok
    status = 0

  end subroutine TriangleArea
  
  !
  ! Area of triangle (xx(:),yy(:)) as (lon,lat) in degrees
  ! with R the earth radius:
  !
  !    area =  R**2 * int cos(y) dy dx
  !                   x,y
  !
  ! Approximate this by first computing the area in degrees**2
  ! and use the average latitude:
  !
  !    area ~  R**2 * int dy dx * cos(y_aver)
  !                   x,y
  !
  
  subroutine LonLatTriangleArea( xx, yy, A, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(3), yy(3)   ! degrees
    real, intent(out)                    ::  A
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LonLatTriangleArea'
    
    ! --- local ----------------------------------
    
    real                    ::  lambda
    real                    ::  Adeg2

    ! --- begin ----------------------------------
    
    ! average latitude:
    lambda = sum(yy)/size(yy) * deg2rad  ! rad
    ! area in deg2:
    call TriangleArea( xx, yy, Adeg2, status )
    IF_NOT_OK_RETURN(status=1)
    ! use triangle area:
    !                            deg2              rad2/deg2                      m2
    A = Adeg2 * deg2rad**2 * cos(lambda) * EarthRadius**2  ! m2

    ! ok
    status = 0

  end subroutine LonLatTriangleArea
  
  !
  ! Approximate area of rectangle in (lon,lat) in degrees using 2 triangles:
  !
  !      4 o---o 3
  !        | / | 
  !      1 o--o 2
  !
  ! Althouth the exact formula is not very difficult, this routine is used
  ! to ensure that area's are approximated in the same way ...
  !
  
  subroutine LonLatRectangleArea( xx, yy, A, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(4), yy(4)   ! degrees
    real, intent(out)                    ::  A
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LonLatRectangleArea'
    
    ! --- local ----------------------------------

    real    ::  A1, A2

    ! --- begin ----------------------------------
    
    ! triangle:
    call LonLatTriangleArea( xx(1:3), yy(1:3), A1, status )
    IF_NOT_OK_RETURN(status=1)
    ! triangle:
    call LonLatTriangleArea( (/xx(3),xx(4),xx(1)/), (/yy(3),yy(4),yy(1)/), A2, status )
    IF_NOT_OK_RETURN(status=1)
    ! combine:
    A = A1 + A2
    
    ! ok
    status = 0

  end subroutine LonLatRectangleArea



  
  !
  ! Arrays xx(:) and yy(:) define corners of a polygon (longitude and latitudes).
  ! If the polygon has 4 or more sides, first devide into triangles:
  !
  !           o-----o
  !          |\  / |
  !         |  o  |
  !        |/  \ |
  !       o-----o
  !
  ! Recursively divide triangles into 2 new triangles using the
  ! midle of the longest edage and the oposite tip as corners:
  !
  !           o        
  !         / | \      ;  do not remove whitespace after the slash ...
  !       /   |   \    ;  do not remove whitespace after the slash ...
  !     o-----*-----o
  !
  ! Return values:
  !   xxp(:), yyp(:) : centroids (mass middle points) of triangles
  !   wwp(:)         : triangle areas in m2
  !
  ! Level 0 returns one point.
  !

  recursive subroutine LonLatPolygonCentroids( xx, yy, level, maxlevel, xxp, yyp, wwp, np, status )
    
    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(:), yy(:)   ! [degree]
    integer, intent(in)                  ::  level
    integer, intent(in)                  ::  maxlevel
    real, intent(inout)                  ::  xxp(:), yyp(:)   ! (np+) centroids   [degree]
    real, intent(inout)                  ::  wwp(:)           ! (np+) triangle area's [m2]
    integer, intent(out)                 ::  np
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/LonLatPolygonCentroids'
    
    ! --- local ----------------------------------
    
    integer                 ::  nc
    integer                 ::  nmax
    integer                 ::  i
    real                    ::  xc, yc
    real                    ::  area
    real                    ::  x1, x2
    real                    ::  y1, y2
    real                    ::  dd(3)
    integer                 ::  ibase
    integer                 ::  ib1, ib2
    integer                 ::  itop
    real                    ::  xmid, ymid
    real, allocatable       ::  xxt(:), yyt(:)
    real, allocatable       ::  wwt(:)
    integer                 ::  nt

    ! --- begin ----------------------------------
    
    ! number of corners:
    nc = size(xx)
    
    ! centroid:
    xc = sum(xx)/nc
    yc = sum(yy)/nc
    
    ! reached end?
    if ( level == maxlevel ) then
    
      ! return center:
      xxp(1) = xc
      yyp(1) = yc
      ! single centroid requested (in polygon) ?
      if ( maxlevel == 0 ) then
        ! single centroid, assign polygon area;
        !~ init sum:
        wwp(1) = 0.0
        !~ loop over sides:
        do i = 1, nc
          ! start and end point:
          x1 = xx(i)
          y1 = yy(i)
          if ( i == nc ) then
            x2 = xx(1)
            y2 = yy(1)
          else
            x2 = xx(i+1)
            y2 = yy(i+1)
          end if
          ! new triangle:
          call LonLatTriangleArea( (/x1,x2,xc/), (/y1,y2,yc/), area, status )
          IF_NOT_OK_RETURN(status=1)
          ! add contribution:
          wwp(1) = wwp(1) + area
        end do ! sides
      else
        ! weight by area:
        call LonLatTriangleArea( xx, yy, wwp(1), status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! set counter:
      np = 1
      
      !! OUTPUT FOR TEST PROGRAM, DO NOT REMOVE!
      !print *, xc, yc, xx, yy
      
    else
    
      ! storage:
      allocate( xxt(size(xxp)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( yyt(size(yyp)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( wwt(size(wwp)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! no data yet:
      np = 0

      ! devide in triangles?
      if ( nc /= 3 ) then

        !           o-----o
        !          |\  / |
        !         |  o  |
        !        |/  \ |
        !       o-----o
      
        ! loop over sides:
        do i = 1, nc
          ! start and end point:
          x1 = xx(i)
          y1 = yy(i)
          if ( i == nc ) then
            x2 = xx(1)
            y2 = yy(1)
          else
            x2 = xx(i+1)
            y2 = yy(i+1)
          end if
          ! new triangle:
          call LonLatPolygonCentroids( (/x1,x2,xc/), (/y1,y2,yc/), level+1, maxlevel, &
                                       xxt, yyt, wwt, nt, status )
          IF_NOT_OK_RETURN(status=1)
          ! extend:
          xxp(np+1:np+nt) = xxt(1:nt)
          yyp(np+1:np+nt) = yyt(1:nt)
          wwp(np+1:np+nt) = wwt(1:nt)
          np = np+nt
        end do ! sides
      
      else
        
        !        top
        !         3
        !         o
        !     3 / | \ 2     ;   do not remove whitespace after the slash ...
        !     /   |   \     ;   do not remove whitespace after the slash ...
        !   o-----*-----o
        !   1     1     2
        !        mid
        !
        
        ! need to search for longest edge;
        ! compute squared edge lengths:
        dd(1) = (xx(2) - xx(1))**2 + (yy(2) - yy(1))**2
        dd(2) = (xx(3) - xx(2))**2 + (yy(3) - yy(2))**2
        dd(3) = (xx(1) - xx(3))**2 + (yy(1) - yy(3))**2
        ! edge with longest length is the base:
        ibase = maxloc( dd, 1 )
        ! indices:
        ib1  = ibase                   !  1  2  3
        ib2  = mod( ibase  , 3 ) + 1   !  2  3  1
        itop = mod( ibase+1, 3 ) + 1   !  3  1  2
        ! mid of base:
        xmid = ( xx(ib1) + xx(ib2) )/2.0
        ymid = ( yy(ib1) + yy(ib2) )/2.0
      
        ! new triangle:
        call LonLatPolygonCentroids( (/xx(ib1),xmid,xx(itop)/), (/yy(ib1),ymid,yy(itop)/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOT_OK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
        ! new triangle:
        call LonLatPolygonCentroids( (/xmid,xx(ib2),xx(itop)/), (/ymid,yy(ib2),yy(itop)/), level+1, maxlevel, &
                                      xxt, yyt, wwt, nt, status )
        IF_NOT_OK_RETURN(status=1)
        ! extend:
        xxp(np+1:np+nt) = xxt(1:nt)
        yyp(np+1:np+nt) = yyt(1:nt)
        wwp(np+1:np+nt) = wwt(1:nt)
        np = np+nt
        
      end if

      ! clear:
      deallocate( xxt, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( yyt, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( wwt, stat=status )
      IF_NOT_OK_RETURN(status=1)
    
    end if  ! end recursion?
    
    ! ok
    status = 0
    
  end subroutine LonLatPolygonCentroids
  
  ! *
  
  !
  ! Given polygon defined by corners [(xx(:),yy(:)],
  ! return arrays [xxp(:),yyp(:)] with locations of points distributed within the polygon;
  ! also return wwp(:) with area corresponding to each point.
  !
  ! The polygon is recursively devided in triangles, the points are their centroids.
  ! Each triangle is devided into 2 new triangles using the
  ! midle of the longest edage and the oposite tip as corners:
  !
  !           o        
  !         / | \      ;  do not remove whitespace after the slash ...
  !       /   |   \    ;  do not remove whitespace after the slash ...
  !     o-----*-----o
  !
  ! Recursion is repeated up to "levels" deep;
  ! for a 6-sided polygon the number of points is
  ! 1 (levels==0), 6 (levels=1), 12 (levels==2), 24 (levels=3), 48 (levels=4), 96 (levels=5, default)
  ! With "levels==0" a single point is returned.
  ! Output arrays are allocated at return.
  !
  
  subroutine GetPolygonPoints( xx, yy, xxp, yyp, wwp, status, levels )

    ! --- in/out ---------------------------------
    
    real, intent(in)                     ::  xx(:), yy(:)     ! [degree]
    real, allocatable, intent(inout)     ::  xxp(:), yyp(:)   ! [degree]
    real, allocatable, intent(inout)     ::  wwp(:)           ! [m2]
    integer, intent(out)                 ::  status
    integer, intent(in), optional        ::  levels

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GetPolygonPoints'
    
    ! --- local ----------------------------------
    
    integer                 ::  nc
    integer                 ::  ns
    integer                 ::  maxp
    integer                 ::  np
    integer                 ::  i, j
    integer                 ::  maxlevel

    ! --- begin ----------------------------------
    
    ! number of corners:
    nc = size(xx)
    
    ! default:
    maxlevel = 5
    if ( present(levels) ) maxlevel = levels

    ! number of centroids expected:
    do i = 0, maxlevel
      if ( i == 0 ) then
        maxp = 1
      else if ( i == 1 ) then
        maxp = nc
      else
        maxp = maxp * 2
      end if
    end do

    ! already allocated?
    if ( allocated(xxp) ) then
      ! check size:
      if ( size(xxp) /= maxp ) then
        ! clear:
        deallocate( xxp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( yyp, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( wwp, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
    end if
    ! not allocated ?
    if ( .not. allocated(xxp) ) then
      ! storage:
      allocate( xxp(maxp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( yyp(maxp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( wwp(maxp), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! recursive call:
    call LonLatPolygonCentroids( xx, yy, 0, maxlevel, xxp, yyp, wwp, np, status )
    IF_NOT_OK_RETURN(status=1)

    ! check ...
    if ( np /= maxp ) then
      write (csol,'("filled ",i0," centroids while expected ",i0)') np, maxp; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
    
  end subroutine GetPolygonPoints
  

end module CSO_Tools


! ***


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Testing generation of sample points in footprints
!!!
!!! - Uncomment in this file:
!!!   - the dummy "CSO_Logging" module in the top
!!!   - the "Test" program below
!!!
!!! - In subroutine "LonLatPolygonCentroids" enable the print statement:
!!!     print *, xc, yc, xx, yy
!!!
!!! - Compile with:
!!!     make test.x
!!!
!!! - Run with redirection of output:
!!!     ./test.x > test.out
!!!
!!! - Create a python script from the code at the bottom of this file
!!!   and try to make a plot of the footprint.
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!program Test
!
!  use CSO_Logging, only : csol, csoPr, csoErr  
!  use CSO_Tools  , only : GetPolygonPoints
!  
!  implicit None
!    
!  character(len=*), parameter   :: rname = 'Test'
!  
!  integer, parameter    ::  nc = 4
!  real, parameter       ::  xx(nc) = (/0.5,6.5,7.5,1.5/), yy(nc) = (/1.5,0.5,2.5,3.5/)
!  
!  !integer, parameter    ::  nc = 6
!  !real, parameter       ::  xx(nc) = (/0.5,2.5,4.5,5.5,3.5,1.5/), yy(nc) = (/2.5,0.5,1.5,3.5,5.5,4.5/) 
!  
!  integer, parameter    ::  maxlevel = 5
!
!  real, allocatable     ::  xxp(:), yyp(:), wwp(:)
!  integer               ::  status
!
!  ! get points:
!  call GetPolygonPoints( xx, yy, xxp, yyp, wwp, status, levels=maxlevel )
!  IF_NOT_OK_RETURN(status=1)
!
!end program Test
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! Test program for plotting.
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  #! /usr/bin/env python
!
!  # read output from 'test.x' program
!  with open('test.out','r') as f : lines = f.readlines()
!
!  # new:
!  fig = plt.figure( figsize=(8,4) )
!  ax = fig.add_axes([0.1,0.1,0.8,0.8], aspect=1 )
!
!  # rows with:
!  #   xc, yc, xx, yy
!  for line in lines :
!      # convert:
!      values = list(map(float,line.split()))
!      # extract:
!      xc = values[0]
!      yc = values[1]
!      nc = int((len(values)-2)/2)
!      xx = values[2:2+nc]
!      yy = values[2+nc:2+2*nc]
!      # add:
!      ax.plot( xx+[xx[0]], yy+[yy[0]], linestyle='-', color='k' )
!      ax.plot( xc, yc, marker='.', color='k' )
!  #endfor
!
!  # axes:
!  ax.set_xlim([0,8])
!  ax.set_ylim([0,4])
!
!  # show:
!  plt.show()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! End
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
