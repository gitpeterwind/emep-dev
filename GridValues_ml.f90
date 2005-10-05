
                          Module GridValues_ml

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Define parameters and variables associated with 3-D grid and its
!  geography.
!
! History: 
! March - changed folllwing Steffen's optimisation/correction of sigma_mid.
! January 2001 : Created by ds from old defconstants. I have made enough changes
! to get this into F90, and to make x,y inputs to the position subroutine,
! but the basic equations are untouched.
! October 2001 hf added call to ReadField (which now does global2local)
! Nov. 2001 - tidied up a bit (ds). Use statements moved to top of module
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 use Io_ml,                only : IO_ROUGH
 use ModelConstants_ml,    only : KMAX_BND, KMAX_MID  &! vertical extent
      ,DEBUG_i, DEBUG_j     ! global coordinate of debug-site !dsjun2005
 use Par_ml, only : &
        MAXLIMAX,MAXLJMAX   & ! max. possible i, j in this domain
      ,limax,ljmax          & ! actual max.   i, j in this domain
      ,ISMBEG,JSMBEG        & ! start of user-specified domain
      ,gi0,gj0              & ! global coordinates of domain lower l.h. corner
      ,NPROC,me &               ! Number of processors, local processor
!hf BC
      ,IILARDOM,JJLARDOM
 use PhysicalConstants_ml, only : GRAV, PI       ! gravity, pi
 implicit none
 private

 !-- contains subroutine:

 Public :: DefGrid     ! =>  GRIDWIDTH_M, map-factor stuff, calls other routines
 Public :: ij2lb       !pw grid to longitude latitude
 Public :: lb2ij       !pw longitude latitude to grid
 Public :: ij2ij       !pw grid1 to grid2

!Hf BC
 Public :: GlobalPosition     ! => 
 private :: Position   ! => lat(gb), long (gl)
!emep1.2beta private :: SetZ0      ! =>  z0 (from RIVM ???)
!pw moved to DryDep_ml private :: inpar    ! input of iclass
!ds moved instead to Met_ml


  !** 1) Public (saved) Variables from module:

  real, public, save :: &
          xp, yp  &  ! Coordinates of North pole (from infield)
        , fi      &  ! projections rotation angle around y axis (from infield)
        , AN      &  ! Distance on the map from pole to equator (No. of cells)
        ,GRIDWIDTH_M &! width of grid at 60N, in meters (old "h")(from infield)
        ,ref_latitude ! latitude at which projection is true (degrees)

  !/ Variables to define global coordinates of local i,j values. 

  integer, public, save, dimension(0:MAXLIMAX+1) :: i_glob !global coordinates
  integer, public, save, dimension(0:MAXLJMAX+1) :: j_glob !of local i,j


  real, public, save,  dimension(KMAX_BND) ::  &
           sigma_bnd ! sigma, layer boundary 

  real, public, save,  dimension(KMAX_MID) ::  &     ! su corrected!
           sigma_mid   ! sigma layer midpoint

  real, public, save,  dimension(KMAX_MID) ::  carea    ! for budgets ???

  real, public, save,  dimension(MAXLIMAX,MAXLJMAX) :: &
            gl   &               !longitude of EMEP grid center
           ,gb                  !latitude  of EMEP grid center

    real, public, save,  dimension(IILARDOM,JJLARDOM) :: &    !hf BC
            gb_glob,   &               !longitude of EMEP grid center
            gl_glob                    !longitude of EMEP grid center


  real, public, save :: gbacmax,gbacmin,glacmax,glacmin

! EMEP grid definitions (old and official)
  real, public, parameter :: xp_EMEP_official=8.&
                            ,yp_EMEP_official=110.&
                            ,fi_EMEP=-32.&
                            ,GRIDWIDTH_M_EMEP=50000.&
                            ,an_EMEP = 237.7316364 &! = 6.370e6*(1.0+0.5*sqrt(3.0))/50000.
                            ,xp_EMEP_old =  43.0&
                            ,yp_EMEP_old = 121.0

  !/** Map factor stuff:

  real, public, save, dimension(0:MAXLIMAX+1,0:MAXLJMAX+1) ::  &
                xm_i     &    ! map-factor
               ,xm_j     &    ! map-factor
               ,xm2    &    ! xm*xm
               ,xmd    &    ! 1/(xm*xm)    ! ds 1/xm
               ,rpol2       ! square of (distance from pole to i,j
                            ! divided by AN )

  real, public, save, dimension(0:MAXLJMAX+1,0:MAXLIMAX+1) ::  &
               xm2ji  &
              ,xmdji

 ! ds jun2005:
  integer, public, save :: &
              debug_li, debug_lj         ! Local Coordinates of debug-site
  logical, public, save :: debug_proc          ! Processor with debug-site

  !/** internal parameters

  logical, private, parameter ::  DEBUG_GRID = .false.  ! for debugging
  logical,public,parameter :: METEOfelt=.true. !.true. if uses "old" (not CDF) meteo input

contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine DefGrid()
  !-------------------------------------------------------------------! 
  !     defines map parameters and fields for the model.
  !-------------------------------------------------------------------! 

    integer ::  i, j, k, n
    real    ::  an2, x, y

  ! Earth radius = 6.37e6 m, gives gridwidth:

!pw u3    GRIDWIDTH_M = 6.370e6*(1.0+0.5*sqrt(3.0))/AN    ! = 50000.0 m


! NB! HIRLAM uses Earth radius = 6.371e6 m : 
! AN = No. grids from pole to equator
! AN = 6.371e6*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M = 237.768957  

  !/ Define global coordinates of local i,j values. We need to account for
  !  the fact that each parallel domain has its starting cordinate
  !  gi0, gj0, and the user may specify a set of lower-left starting
  !  coordinates for running the model, ISMBEG, JSMBEG
  !       i_glob(i)  = i + gi0 + ISMBEG - 2
  !       j_glob(j)  = j + gj0 + JSMBEG - 2

    i_glob = (/ (n + gi0 + ISMBEG - 2, n=0,MAXLIMAX+1) /) 
    j_glob = (/ (n + gj0 + JSMBEG - 2, n=0,MAXLJMAX+1) /) 

    if ( DEBUG_GRID ) then
        write(*,*) "DefGrid: me,ISMBEG,JSMBEG,gi0,gj0",me,ISMBEG,JSMBEG,gi0,gj0
        write(*,*) "DefGrid: me,i_glob",  me, i_glob(1), i_glob(MAXLIMAX+1)
        write(*,*) "DefGrid: me,j_glob",  me, j_glob(1), j_glob(MAXLJMAX+1)
        write(*,*) "DefGrid: me,MAXLIMAX, limax",me,   MAXLIMAX, limax
        write(*,*) "DefGrid: me,MAXLJMAX, ljmax",me,   MAXLJMAX, ljmax
    end if

    !ds jun2005:
    !ds -------------- Find debug coords  and processor ------------------

    do i = 1, MAXLIMAX
       do j = 1, MAXLJMAX
          if( i_glob(i) == DEBUG_i .and. j_glob(j) == DEBUG_j ) then
             debug_li = i
             debug_lj = j
             debug_proc = .true.
          end if
       end do
    end do
    if( debug_proc ) write(*,*) "GridValues debug:", me, debug_li, debug_lj
    !------------------------------------------------------------------


  !  map factor, and map factor squared.

 if( METEOfelt)then

    ref_latitude=60.
    AN = 6.370e6*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M    ! = 237.7316364 for GRIDWIDTH_M=50 km
    an2 = AN*AN

    do j=0,MAXLJMAX+1           ! ds - changed from ljmax+1
          y = j_glob(j) - yp     ! ds - changed from gj0+JSMBEG-2+j
          y = y*y
          do i=0,MAXLIMAX+1        ! ds - changed from limax+1
              x = i_glob(i) - xp   ! ds - changed
              rpol2(i,j) = (x*x + y)/an2

              xm_i(i,j) = 0.5*(1.0+sin(PI/3.0))*(1.0 + rpol2(i,j))
              xm_j(i,j) = xm_i(i,j)

              xm2(i,j) = xm_i(i,j)*xm_j(i,j)
              xmd(i,j) = 1.0/xm2(i,j)
              xm2ji(j,i) = xm2(i,j)
              xmdji(j,i) = xmd(i,j)
          end do
      end do 

   else
!mapping factor xm and ref_latitude are read from the meteo file

    AN = 6.370e6*(1.0+sin( ref_latitude*PI/180.))/GRIDWIDTH_M    ! = 237.7316364 for GRIDWIDTH_M=50 km and ref_latitude=60

    do j=0,MAXLJMAX+1        
       do i=0,MAXLIMAX+1     
          xm2(i,j) = xm_i(i,j)*xm_j(i,j)
          xmd(i,j) = 1.0/xm2(i,j)
          xm2ji(j,i) = xm2(i,j)
          xmdji(j,i) = xmd(i,j)
       enddo
    enddo

   endif
!     definition of the half-sigma levels (boundaries between layers) from the full levels. 

          sigma_bnd(KMAX_BND) = 1.
          do k = KMAX_MID,2,-1
             sigma_bnd(k) = 2.*sigma_mid(k) - sigma_bnd(k+1)
          enddo
          sigma_bnd(1) = 0.
          
   !
   !     some conversion coefficients needed for budget calculations
   !
    do  k=1,KMAX_MID
        carea(k) = (sigma_bnd(k+1) - sigma_bnd(k))/GRAV*GRIDWIDTH_M*GRIDWIDTH_M
        !write(6,*)'carea,sigma_bnd,h',carea(k),sigma_bnd(k),GRIDWIDTH_M
        !su  cflux(k) = (sigma_bnd(k+1) - sigma_bnd(k))/G*h
    end do

   ! set latitude, longitude

       call Position()

    if ( DEBUG_GRID ) then
        print *," After position" 
        print *, "DefGrid: ", ISMBEG, gi0
        print *, "DefGrid: i_glob",  i_glob(1), i_glob(MAXLIMAX+1)
        print *, "DefGrid:  min, max gl ",  minval(gl), maxval(gl)
        print *, "DefGrid:  min, max gb ",  minval(gb), maxval(gb)
        print *, "DefGrid:  gb - 1,1, MAXLIMAX, MAXLJMAX",  &
                              gb(1,1), gb(MAXLIMAX,MAXLJMAX)
    end if

  end subroutine DefGrid
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine Position()
  !-------------------------------------------------------------------! 
  !      calculates l(lat),b(long) (geographical coord.) 
  !      in every grid point defined by the i_glob, j_glob arrays. 
  !
  !      input:  xp,yp:   coord. of the polar point.
  !              AN:      number of grid-distances from pole to equator.
  !              fi:      rotational angle for the x,y grid.
  !              imax,jmax:   number of points in  x- og y- direction
  !              glmin:   gives min.value of geographical lenght
  !                       =>  glmin <= l <= glmin+360.  
  !                           (example glmin = -180. or 0.)
  !                       if "geopos","georek" is used
  !                       then glmin must be the lenght i(1,1) in the
  !                       geographical grid (gl1 to "geopos")
  !      output: gl(ii,jj): latitude glmin <= l <= glmin+360. 
  !              gb(ii,jj): longitude  -90. <= b <= +90. 
  !-------------------------------------------------------------------! 
  ! ds  - replaced pi by PI, glx by glmax, removed gln
  ! ds  - evaluate gl, gb over whole domain given by MAXLIMAX, MAXLJMAX
  !       to safeguard against possible use of non-defined gl,bb squares. 
  ! ds  - note, we could use rpol2(i,j) to save some computations here, 
  !       but for now we leave it. This stuff is only done once anyway

    real    :: glmin, glmax, om, om2, dy, dy2,rp,rb, rl, dx, dr !,fi read in Met_ml pw u3 
    integer :: i, j, info
      !                                                    
      !      call xylbf(43.,121.,237.,-32.,gl,gb,imax,jmax,-180.)
      !      subroutine xylbf(xp,yp,AN,fi,gl,gb,ii,jj,glmin)
      !su    xp,yp read in infield                
      !su    xp = 43.
      !su    yp = 121.

!    fi = -32.0   !read in Met_ml pw u3
    glmin = -180.0

    glmax = glmin + 360.0
    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees (om=Norwegian omvendt?)
    om2   = om * 2.0

    do j = 1, MAXLJMAX          ! ds - changed from ljmax
       dy  = yp - j_glob(j)     ! ds - changed from gj0+JSMBEG-2+j
       dy2 = dy*dy
       do i = 1, MAXLIMAX       ! ds - changed from limax
         dx = i_glob(i) - xp    ! ds - changed
         rp = sqrt(dx*dx+dy2)           ! => distance to pole
         rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
         rl = 0.0
         if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         gl(i,j)=rl                     !     longitude
         gb(i,j)=rb                     !     latitude

       end do ! i
    end do ! j

    !su - added test to find global-min and max lat/long values

    gbacmax = maxval(gb(:,:))
    gbacmin = minval(gb(:,:))
    glacmax = maxval(gl(:,:))
    glacmin = minval(gl(:,:))
    call gc_rmax(1, NPROC, info, gbacmax)
    call gc_rmin(1, NPROC, info, gbacmin)
    call gc_rmax(1, NPROC, info, glacmax)
    call gc_rmin(1, NPROC, info, glacmin)

    if(me==0) write(unit=6,fmt=*) "max/min for gb,gl", &
                     gbacmax,gbacmin,glacmax,glacmin

    if ( DEBUG_GRID ) then
       do j = 1, MAXLJMAX
         do i = 1, MAXLIMAX
           if ( i_glob(i) == 91 .and. j_glob(j) == 71 ) then ! Rorvik
             print "(a20,4i4,2f8.2,f7.3)", &
              "Position: i,j,i_glob,j_glob,gl,gb,rp: ",  &
              i,j, i_glob(i), j_glob(j), gl(i,j), gb(i,j),rp
           end if
         end do
       end do
    end if

  end subroutine Position
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine GlobalPosition

    integer i,j
    real :: dr,om,om2,rb,rl,rp,dx,dy,dy2,glmax,glmin

    glmin = -180.0
    glmax = glmin + 360.0

    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees (om=Norwegian omvendt?)
    om2   = om * 2.0

    do j = 1, JJLARDOM
       dy  = yp - j  
       dy2 = dy*dy
       do i = 1, IILARDOM
         dx = i - xp    
         rp = sqrt(dx*dx+dy2)           ! => distance to pole
         rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
         rl = 0.0
         if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         gl_glob(i,j)=rl                     !     longitude
         gb_glob(i,j)=rb                     !     latitude

       end do ! i
    end do ! j

end subroutine GlobalPosition


  subroutine lb2ij(imax,jmax,gl,gb,ir2,jr2,fi2,an2,xp2,yp2)
    !-------------------------------------------------------------------! 
    !      calculates coordinates ir2, jr2 (real values) from gl(lat),gb(long) 
    !
    !      input:  xp2,yp2:   coord. of the polar point in grid2
    !              an2:   number of grid-distances from pole to equator in grid2.
    !              fi2:      rotational angle for the grid2 (at i2=0).
    !              i1max,j1max: number of points (grid1) in  x- og y- direction
    !
    !
    !      output: i2(i1,j1): i coordinates in grid2 
    !              j2(i1,j1): j coordinates in grid2 
    !-------------------------------------------------------------------! 

!    use Par_ml,   only : MAXLIMAX, MAXLJMAX

    implicit none


    integer :: imax,jmax,i1, j1
    real    :: fi2,an2,xp2,yp2
    real    :: gl(imax,jmax),gb(imax,jmax)
    real    :: ir2(imax,jmax),jr2(imax,jmax)

    real, parameter :: PI=3.14159265358979323
    real    :: PId4,dr,dr2


    PId4    = PI/4.      
    dr2    = PI/180.0/2.      ! degrees to radians /2
    dr    = PI/180.0      ! degrees to radians 

    do j1 = 1, jmax
       do i1 = 1, imax

          ir2(i1,j1)=xp2+an2*tan(PId4-gb(i1,j1)*dr2)*sin(dr*(gl(i1,j1)-fi2))
          jr2(i1,j1)=yp2-an2*tan(PId4-gb(i1,j1)*dr2)*cos(dr*(gl(i1,j1)-fi2))
!          write(*,*)i1,j1,ir2(i1,j1),jr2(i1,j1),gl(i1,j1),gb(i1,j1)

       end do ! i
    end do ! j

    return
  end subroutine lb2ij

  subroutine ij2lb(imax,jmax,gl,gb,fi,an,xp,yp)
  !-------------------------------------------------------------------! 
  !      calculates l(lat),b(long) (geographical coord.) 
  !      in every grid point. 
  !
  !      input:  xp,yp:   coord. of the polar point.
  !              an:      number of grid-distances from pole to equator.
  !              fi:      rotational angle for the x,y grid (at i=0).
  !              imax,jmax:   number of points in  x- og y- direction
  !              glmin:   gives min.value of geographical lenght
  !                       =>  glmin <= l <= glmin+360.  
  !                           (example glmin = -180. or 0.)
  !                       if "geopos","georek" is used
  !                       then glmin must be the lenght i(1,1) in the
  !                       geographical grid (gl1 to "geopos")
  !      output: gl(ii,jj): longitude glmin <= l <= glmin+360. 
  !              gb(ii,jj): latitude  -90. <= b <= +90. 
  !-------------------------------------------------------------------! 



    implicit none


    integer :: i, j, imax, jmax
    real    :: gl(imax,jmax),gb(imax,jmax)
    real    :: fi, an, xp, yp 
    real    :: om, om2, glmin, glmax,dy, dy2,rp,rb, rl, dx, dr
    real, parameter :: PI=3.14159265358979323


!    fi = -32.0
    glmin = -180.0

    glmax = glmin + 360.0
    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees (om=Norwegian omvendt?)
    om2   = om * 2.0

    do j = 1, jmax          
       dy  = yp - j            
       dy2 = dy*dy
       do i = 1, imax       

         dx = i - xp    ! ds - changed
         rp = sqrt(dx*dx+dy2)           ! => distance to pole
         rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
         rl = 0.0
         if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         gl(i,j)=rl                     !     longitude
         gb(i,j)=rb                     !     latitude
!         write(*,*)i,j,gl(i,j),gb(i,j)
       end do ! i
    end do ! j

   return
  end subroutine ij2lb

  subroutine ij2ij(in_field,imaxin,jmaxin,out_field,imaxout,jmaxout, &
                   fiin,anin,xpin,ypin,fiout,anout,xpout,ypout)

!   Converts data (in_field) stored in coordinates (fiin,anin,xpin,ypin) 
!   into data (out_field) in coordinates (fiout,anout,xpout,ypout)
!   pw august 2002

    use Par_ml   ,      only :  me, NPROC

	implicit none
        integer, intent(in) :: imaxin,jmaxin,imaxout,jmaxout
        real, intent(in) :: fiin,anin,xpin,ypin,fiout,anout,xpout,ypout
        real, intent(in) :: in_field(imaxin,jmaxin)! Field to be transformed
        real, intent(out) :: out_field(imaxout,jmaxout)! Field to be transformed

        real, allocatable,dimension(:,:) :: x,y,gb,gl
        integer alloc_err,i,j,i2,j2
        logical :: interpolate
        real :: f11,f12,f21,f22

        interpolate = .true.
!        interpolate = .false.

       allocate(x(imaxout,jmaxout), stat=alloc_err)
       allocate(y(imaxout,jmaxout), stat=alloc_err)
       allocate(gb(imaxout,jmaxout), stat=alloc_err)
       allocate(gl(imaxout,jmaxout), stat=alloc_err)
       if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "ij2ij alloc failed")

! find longitude, latitude of wanted area
    call ij2lb(imaxout,jmaxout,gl,gb,fiout,anout,xpout,ypout)

! find corresponding coordinates (i,j) in in_field coordinates 
    call lb2ij(imaxout,jmaxout,gl,gb,x,y,fiin,anin,xpin,ypin)


    ! check if the corners of the domain are inside the area covered by the 
    ! in_grid: (In principle we should test for all i,j , but test the corners
    ! should be good enough in practice) 

    if(int(x(1,1)) < 1 .or. int(x(1,1))+1 > imaxin .or. &
         int(x(imaxout,1)) < 1 .or. int(x(imaxout,1))+1 > imaxin .or. &
         int(x(1,jmaxout)) < 1 .or. int(x(1,jmaxout))+1 > imaxin .or. &
         int(x(imaxout,jmaxout)) < 1 .or. &
          int(x(imaxout,jmaxout))+1 > imaxin .or. &
         int(y(1,1)) < 1 .or. int(y(1,1))+1 > jmaxin .or. &
         int(y(imaxout,1)) < 1 .or. int(y(imaxout,1))+1 > jmaxin .or. &
         int(y(1,jmaxout)) < 1 .or. int(y(1,jmaxout))+1 > jmaxin .or. &
         int(y(imaxout,jmaxout)) < 1 .or. &
         int(y(imaxout,jmaxout))+1 > jmaxin ) then
       write(*,*)'Did not find all the necessary data in in_field'
       write(*,*)'values needed: '
       write(*,*)x(1,1),y(1,1)
       write(*,*)x(imaxout,1),y(imaxout,1)
       write(*,*)x(1,jmaxout),y(1,jmaxout)
       write(*,*)x(imaxout,jmaxout),y(imaxout,jmaxout)
       write(*,*)'max values found: ',imaxin ,jmaxin
       call gc_abort(me,NPROC, "ij2ij: area to small")
    endif



!  interpolate fields if required
!
    if(interpolate)then
    do j = 1, jmaxout
       do i = 1,imaxout
          i2=int(x(i,j))
          j2=int(y(i,j))
          f11=(1.-(x(i,j)-i2))*(1.-(y(i,j)-j2))
          f12=(1.-(x(i,j)-i2))*((y(i,j)-j2))
          f21=((x(i,j)-i2))*(1.-(y(i,j)-j2))
          f22=((x(i,j)-i2))*((y(i,j)-j2))

          out_field(i,j) =  &
               f11 * in_field(i2,j2) +  & 
               f12 * in_field(i2,j2+1) +  & 
               f21 * in_field(i2+1,j2) +  & 
               f22 * in_field(i2+1,j2+1)

       enddo
    enddo
    else

    do j = 1, jmaxout
       do i = 1,imaxout
          out_field(i,j) =in_field(nint(x(i,j)),nint(y(i,j)))
       enddo
    enddo

    endif

       deallocate(x,stat=alloc_err)
       deallocate(y,stat=alloc_err)
       deallocate(gb,stat=alloc_err)
       deallocate(gl,stat=alloc_err)
       if ( alloc_err /= 0 ) call gc_abort(me,NPROC,"ij2ij de-alloc_err")
    
    end subroutine ij2ij

  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module GridValues_ml
!==============================================================================
