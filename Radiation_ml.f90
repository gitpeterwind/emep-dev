module Radiation_ml
  !---------------------------------------------------------------------------
  !DESCRIPTION
  !* written by brian a. ridley at york univ., ca. 1978
  !*  based on equations given by w.h.smart "textbook of spherical astronomy,"
  !*                6th ed. cambridge u. press
  !* this subroutine calculates solar zenith and azimuth angles for a particular
  !* time and location.  must specify:
  !*
  !* input:
  !*       gb - latitude in decimal degrees 
  !*       gl - longitude in decimal degrees
  !*       idate - date at greenwich - specify year (19yy), month (mm), day (dd)
  !*              format is six-digit integer:  yymmdd
  !*       thour - erstatter gmt
  !*       gmt  - greenwich mean time - decimal military eg.
  !*               22.75 = 45 min after ten pm gmt
  !* output  
  !*       zeta
  !!
  !* Optimised by Steffen Unger, 2000?
  !* converted to F90 module, Dave Simpson, Nov. 2001
  !* ToDo: Check leap-years for centuries etc.
  !---------------------------------------------------------------------------

    use Par_ml   , only : MAXLIMAX,MAXLJMAX,li0,li1,lj0,lj1
    use Dates_ml , only: nmdays
    use GridValues_ml,         only: gl,gb
    use ModelConstants_ml ,    only: current_date
    use PhysicalConstants_ml , only: PI
    implicit none
    private

    real, public, dimension(MAXLIMAX, MAXLJMAX), save:: &
        zen    &  ! Zenith angle (degrees)
       ,coszen    ! cos of zenith angle

  public :: ZenAng   ! subroutine

contains
 !6c subroutine zenit_angle(zeta,thour)
 !<===========================================================================
  subroutine ZenAng(thour)

    real, intent(in) :: thour

!6c    real :: zeta(MAXLIMAX,MAXLJMAX)
    real :: dr
    integer :: ii,jj,i,iiyear,iyear,imth,iday,iiy,nyears
    integer :: leap,noleap,ijd,in,ij
    real :: d,jd,yr,ml,rml,w,wr,ec,epsi,yt,pepsi,cww
    real :: sw,ssw, eyt, feqt1, feqt2, feqt3, feqt4, feqt5, &
            feqt6, feqt7, feqt,eqt,ra,reqt,zpt
    real :: rdecl,sinrdecl,cosrdecl,rlt,lzgmt,lbgmt,csz,rra,tab

!* convert to radians

    dr = pi/180.

!* parse date

    iiyear = current_date%year
    iyear = 19*100 + iiyear
    imth = current_date%month
    iday = current_date%day

!* identify and correct leap years

    iiy = (iiyear/4)*4

!* count days from dec.31,1973 to jan 1

    nyears = iyear - 1974
    leap = (nyears+1)/4
    if(nyears <= -1) leap = (nyears-2)/4
    noleap = nyears - leap
    yr = 365.0*noleap + 366.0*leap

    ijd = 0
    in = imth - 1
!ds    if(in == 0) go to 40
!ds    do i=1,in
!ds      ijd = ijd + nmdays(i)
!ds    end do
!ds    ijd = ijd + iday
!ds    go to 50
!ds40  ijd = iday
!ds50  ij = iyear - 1973

    if(in == 0) then
        ijd = iday
    else
        do i=1,in
          ijd = ijd + nmdays(i)
        end do
        ijd = ijd + iday
    end if   
    ij = iyear - 1973

!* julian days current "ijd"  
!ds comment - actually Julian dates started with jday=1 at 12 noon
! in the year 13xx so the word is used incorrectly here.

    jd = ijd + yr
    d = jd + thour/24.0

!* calc geom mean longitude

    ml = 279.2801988 + .9856473354*d + 2.267e-13*d*d
    rml = ml*dr

!* calc equation of time in sec
!*  w = mean long of perigee
!*  e = eccentricity
!*  epsi = mean obliquity of ecliptic

    w = 282.4932328 + 4.70684e-5*d + 3.39e-13*d*d
    wr = w*dr
    ec = 1.6720041e-2 - 1.1444e-9*d - 9.4e-17*d*d
    epsi = 23.44266511 - 3.5626e-7*d - 1.23e-15*d*d
    pepsi = epsi*dr
    yt = (tan(pepsi/2.0))**2
    cww = cos(wr)
    sw = sin(wr)
    ssw = sin(2.0*wr)
    eyt = 2.*ec*yt
    feqt1 = sin(rml)*(-eyt*cww - 2.*ec*cww)
    feqt2 = cos(rml)*(2.*ec*sw - eyt*sw)
    feqt3 = sin(2.*rml)*(yt - (5.*ec**2/4.)*(cww**2-sw**2))
    feqt4 = cos(2.*rml)*(5.*ec**2*ssw/4.) 
    feqt5 = sin(3.*rml)*(eyt*cww)
    feqt6 = cos(3.*rml)*(-eyt*sw)
    feqt7 = -sin(4.*rml)*(.5*yt**2)
    feqt = feqt1 + feqt2 + feqt3 + feqt4 + feqt5 + feqt6 + feqt7
    eqt = feqt*13751.0

!* equation of time in hrs:

!    eqh = eqt/3600.

!* convert eq of time from sec to deg

    reqt = eqt/240.

!* calc right ascension in rads

    ra = ml - reqt
    rra = ra*dr

!* calc declination in rads, deg

    tab = 0.43360*sin(rra)
    rdecl = atan(tab)
    sinrdecl = sin(rdecl)
    cosrdecl = cos(rdecl)
!su    decl = rdecl/dr         

!* calc local hour angle

    do jj=lj0,lj1
      do ii=li0,li1
        rlt = gb(ii,jj)*dr
        lbgmt = 12.0 - eqt/3600. - gl(ii,jj)*24./360.
        lzgmt = 15.0*(thour - lbgmt)
        zpt = lzgmt*dr
        csz = sin(rlt)*sinrdecl + cos(rlt)*cosrdecl*cos(zpt)
!su        csz = sin(rlt)*sin(rdecl) + cos(rlt)*cos(rdecl)*cos(zpt)
    
        if(csz  >  1.) csz = 1.
        if(csz  < -1.) csz = -1.
        coszen(ii,jj) = csz
        zen(ii,jj) = acos(csz)/dr 
!6c        zeta(ii,jj) = csz
!        raz = acos(csz)
!        zeta(ii,jj) = raz/dr
      end do
    end do

 end subroutine ZenAng
 !6cend subroutine zenit_angle
 !<===========================================================================

end module Radiation_ml
