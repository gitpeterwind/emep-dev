
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

!u1  use My_Runmode_ml,        only: DEBUG
 use Io_ml,                only : IO_ROUGH
 use ModelConstants_ml,    only : KMAX_BND, KMAX_MID   ! vertical extent
 use ReadField_ml,         only : Readfield
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
!Hf BC
 Public :: GlobalPosition     ! => 
 private :: Position   ! => lat(gb), long (gl)
 private :: SetZ0      ! =>  z0 (from RIVM ???)
 private :: inpar    ! input of iclass


  !** 1) Public (saved) Variables from module:

  real, public, save :: &
          xp, yp  &  ! Coordinates of North pole (from infield)
        , fi      &  ! projections rotation angle around y axis (from infield)
        , AN      &  ! Distance on the map from pole to equator (No. of cells)
        ,GRIDWIDTH_M ! width of grid at 60N, in meters (old "h")(from infield)

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
           ,gb   &               !latitude  of EMEP grid center
          , z0              !ds: From LAM50 - replace a.s.a.p. !

!hf BC
    real, public, save,  dimension(IILARDOM,JJLARDOM) :: &
            gb_glob,   &               !longitude of EMEP grid center
            gl_glob                    !longitude of EMEP grid center


  real, public, save :: gbacmax,gbacmin,glacmax,glacmin

  integer, public, save, &
          dimension(MAXLIMAX,MAXLJMAX) ::  iclass  ! roughness class

  !/** Map factor stuff:

  real, public, save, dimension(0:MAXLIMAX+1,0:MAXLJMAX+1) ::  &
                xm     &    ! map-factor
               ,xm2    &    ! xm*xm
               ,xmd    &    ! 1/xm
               ,rpol2       ! square of (distance from pole to i,j
                            ! divided by AN )

  real, public, save, dimension(0:MAXLJMAX+1,0:MAXLIMAX+1) ::  &
               xm2ji  &
              ,xmdji

  !/** internal parameters

!pw u3  real, private, parameter :: AN = 237.731644  !No. grids from pole to equator?

  logical, private, parameter ::  DEBUG_GRID = .false.  ! for debugging

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

    AN = 6.370e6*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M    ! = 237.7316364 for GRIDWIDTH_M=50 km

! NB! HIRLAM uses Earth radius = 6.371e6 m : 
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
        print *, "DefGrid: me,ISMBEG,JSMBEG,gi0,gj0",me,ISMBEG,JSMBEG,gi0,gj0
        print *, "DefGrid: me,i_glob",  me, i_glob(1), i_glob(MAXLIMAX+1)
        print *, "DefGrid: me,j_glob",  me, j_glob(1), j_glob(MAXLJMAX+1)
        print *, "DefGrid: me,MAXLIMAX, limax",me,   MAXLIMAX, limax
        print *, "DefGrid: me,MAXLJMAX, ljmax",me,   MAXLJMAX, ljmax
    end if


  !  map factor, and map factor squared.  !ko proper EMEP grid definition

     an2 = AN*AN
     do j=0,MAXLJMAX+1           ! ds - changed from ljmax+1
          y = j_glob(j) - yp     ! ds - changed from gj0+JSMBEG-2+j
          y = y*y
          do i=0,MAXLIMAX+1        ! ds - changed from limax+1
              x = i_glob(i) - xp   ! ds - changed
              rpol2(i,j) = (x*x + y)/an2

!ko           xm(i,j) = 0.933/an2*(an2 + rpol2(i,j))
              xm(i,j) = 0.5*(1.0+sin(PI/3.0))*(1.0 + rpol2(i,j))
!ko
              xm2(i,j) = xm(i,j)*xm(i,j)
              xmd(i,j) = 1.0/xm2(i,j)
              xm2ji(j,i) = xm2(i,j)
              xmdji(j,i) = xmd(i,j)
          end do
      end do 


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

   ! read iclass
    call inpar

   ! set z0
    call SetZ0()

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
    real :: dr,om,om2,dy,rb,rl,rp,dx,dy,dy2,glmax,glmin

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
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine inpar

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !     This subroutine reads parameterfields from file
   !     reading surface roughness classes from file: rough.170
   !
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    character*20 fname
    integer i

    integer il,jl,ill,jll

    real, dimension(MAXLIMAX, MAXLJMAX) :: &
                 r_class  ! Roughness (real) in rough.170


      if (me ==  0) then
       write(fname,fmt='(''rough.170'')') 
       write(6,*) 'filename for landuse ',fname
     endif

    !6b Just need to read real numbers from rough.170 into r_class:

      call ReadField(IO_ROUGH,fname,r_class)

     ! And convert from real to integer field

      iclass(:,:)=nint(r_class(:,:)) 


  end subroutine inpar
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine SetZ0()

    !  Set the surface roughness equal to the lam50e values.
    !  (ds comment - all this rougness stuff has to be reviewed/replaced 
    !   with new  deposition modules, at least that from RIVM)

    integer              :: i, j, icl
    real, dimension(0:6) ::  class   ! Some surface roughnesses ??

    class =  (/ 1.0e-4,1.0e-3,3.0e-1,3.0e-1,3.0e-1,3.0e-1,1.0e-3 /)

    do j = 1,ljmax
      do i = 1,limax
         icl = iclass(i,j)
         z0(i,j) = class(icl)
      end do
    end do

  end subroutine SetZ0
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module GridValues_ml

!==============================================================================
!REMOVED
!REM    use par_ml , only : IILARDOM,JJLARDOM        &
!REM            ,MSG_INIT0            &
!
!hf start      open(IO_ROUGH,file='rough.170',err=62)
!
!      do i=1,IILARDOM*JJLARDOM
!        read(IO_ROUGH,*,err=62)ill,jll,r_cl
!          r_class(ill,jll)=nint(r_cl)
!      enddo
!      close(IO_ROUGH)
!    endif
!
!    call global2local_int(r_class,iclass,MSG_INIT0        &
!        ,IILARDOM,JJLARDOM,1,ISMBEG,JSMBEG)
!    goto 1000
!62    continue
!    write(6,*) 'error in opening land use file, landuse'
!
!hf end
!REM 1000    return

  !!6b real, dimension(MAXLIMAX,MAXLJMAX) ::  loc_class ! Real iclass
  !!6b                                         !(Local,real roughness field)
    !!6b real r_cl(IILARDOM,JJLARDOM)
    !!6b integer r_class(IILARDOM,JJLARDOM)
