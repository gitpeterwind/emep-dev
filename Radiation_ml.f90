module Radiation_ml
  !---------------------------------------------------------------------------
  !DESCRIPTION
  !*ZenAng*  written by brian a. ridley at york univ., ca. 1978
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
  !DESCRIPTION
  !*solbio*
  !* from T. Pierce.
  !---------------------------------------------------------------------------
  !* ZenAng Optimised by Steffen Unger, 2000?
  !* converted to F90 module, Dave Simpson, Nov. 2001
  !' Added solbio, UK stuff, ds May 2002
  !* ToDo: Check leap-years for centuries etc.
  !---------------------------------------------------------------------------

    use Par_ml   , only : MAXLIMAX,MAXLJMAX,li0,li1,lj0,lj1
    use Dates_ml , only: nmdays
    use GridValues_ml,         only: gl,gb
    use ModelConstants_ml ,    only: current_date
    use PhysicalConstants_ml , only: PI,DEG2RAD,RAD2DEG
    implicit none
    private

 !  Subroutines:

  public :: ZenAng         ! => coszen=cos(zen), zen=zenith angle (degrees) 
  public :: UKZenAng         ! => coszen=cos(zen), zen=zenith angle (degrees) 
  public :: SolBio         ! => irradiance (W/m2)
  public :: SolBio2D       ! => 2D fields of direct and diffuse irradiance (W/m2)

 ! Public Variables:

    real, public, dimension(MAXLIMAX, MAXLJMAX), save:: &
        zen          &  ! Zenith angle (degrees)
       ,coszen       &  ! cos of zenith angle
       ,Idiffuse     &  ! diffuse solar radiation (W/m^2)        ! ds rv1_6_x
       ,Idirectt        ! total direct solar radiation (W/m^2)   ! ds rv1_6_x

    real, public, save :: daylength(0:90), solarnoon(-180:180)    ! ds_OH

    real, public, parameter :: &
      PARfrac = 0.45,   & ! approximation to fraction (0.45 to 0.5) of total 
                          ! radiation in PAR waveband (400-700nm)
      Wm2_uE  = 4.57      ! converts from W/m^2 to umol/m^2/s

 !u7.lu removed from Dep part
 ! real, public, save :: Idrctn      ! => irradiance (W/m^2), normal to beam
 !                      solar     & ! => irradiance (W/m^2)

  
  logical, private, parameter :: DEBUG = .false.

contains
 !<===========================================================================
  subroutine ZenAng(thour)

    real, intent(in) :: thour

    real :: dr
    integer :: ii,jj,i,iiyear,iyear,imth,iday,iiy,nyears
    integer :: leap,noleap,ijd,in,ij
    real :: d,jd,yr,ml,rml,w,wr,ec,epsi,yt,pepsi,cww
    real :: sw,ssw, eyt, feqt1, feqt2, feqt3, feqt4, feqt5, &
            feqt6, feqt7, feqt,eqt,ra,reqt,zpt
    real :: rdecl,sinrdecl,cosrdecl,rlt,lzgmt,lbgmt,csz,rra  !dsOH,tab
    real :: eqt_h, tan_decl, arg     !dsOH  (tan_decl replaces tab)
    integer :: ilat, ilong           !dsOH
    logical, parameter :: MY_DEBUG = .false.

!* convert to radians

    dr = pi/180.

!* parse date

    iiyear = current_date%year
    !dsOH BUG FIX!! iyear = 19*100 + iiyear
       iyear = iiyear
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

    eqt_h = eqt/3600.   !dsOH added back

!* convert eq of time from sec to deg

    reqt = eqt/240.

!* calc right ascension in rads

    ra = ml - reqt
    rra = ra*dr

!* calc declination in rads, deg

    !dsOH tab = 0.43360*sin(rra)
    tan_decl = 0.43360*sin(rra)

!* ds added: calculation of solar noon and daylength
  ! Caculate solar noon, length of day, following Jones, Appendix 7:

     do ilat = 0, 90

        arg =  -tan_decl * tan( DEG2RAD*ilat )

        if( arg <= -1.0 ) then 
           daylength(ilat) = 24.0  !! Polar summer
        else if( arg >= 1.0 ) then 
           daylength(ilat) = 0.0   !! Polar night
        else
           daylength(ilat) = acos( -tan_decl * tan( DEG2RAD*ilat ) ) * RAD2DEG * 2.0/15.0  ! eqn A7.2, Jones
        end if
     end do

     do ilong = -180, 180
        solarnoon(ilong) = 12.0 - eqt_h - ilong*24./360.
     end do


    rdecl = atan(tan_decl)
    sinrdecl = sin(rdecl)
    cosrdecl = cos(rdecl)
!su    decl = rdecl/dr         

    if ( MY_DEBUG ) then
        write(*,*) "DEBUG ZenAng: iyear ", iyear, " noleap? ", noleap, " leap? ", leap, " yr ", yr
        write(*,*) "DEBUG ZenAng: ijd ", ijd, " decl ", rdecl/dr , " eqt(mins) ", eqt/60.0
    end if
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

  ! ======================================================================
 subroutine UKZenAng(lon, lat, daynr, nydays, hr,coszen,zen)
  ! ======================================================================
  !  routine determines (approximate) cos(zen), where "zen" denotes the zenith 
  !  angle, (in accordance with Iversen and Nordeng (1987, p.28))
  !
  ! arguments:              
    real,    intent(in) ::  lon      ! longitude (degrees), east is positive
    real,    intent(in) ::  lat      ! latitude (degrees), north is positive
    integer, intent(in) ::  daynr    ! day nr. (1..366)
    integer, intent(in) ::  nydays   ! number of days per year (365 or 366)
    real,    intent(in) ::  hr       ! hour  (0-24, gmt)  ! ds - was integer
  
    
    real, intent(out) :: coszen     ! cos(zen) 
    real, intent(out) :: zen        ! zenith angle (degrees)                              
  !
  !                                       dnmi  29-9-95  Hugo Jakobsen

  ! ZenAng to be compared with zenith angle routine in UK model
  ! ======================================================================

     !/ Local....
     real    :: lonr, latr, arg, decl, tangle

     lonr=lon*(PI/180.0)          ! convert to radians
     latr=lat*(PI/180.0)          ! convert to radians

     arg = ((daynr - 80.0)/nydays) * 2.0 * PI  !matches monthly year angle
                                               !calculation in UK model     

     decl = 23.5 * sin(arg) * (PI/180.0)  !EMEP procedure for calculating the 
                                          !sun's declination (in radians): 
                                          !differs       
                                          !from the UK procedure (compare)
     
     tangle = lonr + (hr/12.0-1.0)*PI      !no time correction (change?)
     coszen =(sin(decl)*sin(latr)+cos(decl)*cos(latr)*cos(tangle))
     zen = acos(coszen)*180/PI
  end subroutine UKZenAng

  !=============================================================================
  !u7.lu subroutine SolBio(jday,coszen,zen,cl,pres, solar,Idrctn,Idfuse,Idrctt)
  !u7.lu removed some superfluous variables
  subroutine SolBio(jday,coszen,cl,pres, Idfuse,Idrctt)
  !=============================================================================
  !
  !
  !  function:
  !
  !     computes the radiation terms needed for stomatal calculations
  !     one term is computed: total solar radiation (W/m^2)
  !     methodology for this calculation taken from Iqbal, M., 1983,
  !     An introduction to solar radiation, Academic Press, New York, 
  !     pp. 202-210.
  !
  !  history:
  !
  !     modified for EMEP model
  !     Development of this routine was prompted by the need for a
  !     horizontal rather than an actinic flux calculation (which had
  !     been performed by Soleng). Furthermore, Soleng computed total
  !     radiation only out to the near-ir spectrum. This program
  !     is designed only for approximate radiation estimates to be used
  !     for stomatal calculations.
  !
  !     8/90    initial development of SolBio by T. Pierce
  !     9/95    modified by Hugo Jakobsen, 29/9-95 
  !     1-3/01  modified by Dave Simpson, with albedo as input
  !     and Idrct, Idfuse as output, and converted to F90
  !     19/9/01 - albedo removed, Solar thus gives total incoming
  !     ToDo
  !          - could consider alternative schemes for cloud attenuation,
  !      notably use of low, medium and high cloud (see e.g.  Stull (1988) 
  !          - check clearness number, pressure assumptions (too American?)
  !
  !  argument list description:
  !
      integer, intent(in)  :: jday   ! day no. from 1 Jan (Julian day)
      real,    intent(in)  :: coszen ! cos(zen) where zen is zenith angle

      real,    intent(in)  :: cl     ! opaque sky cover (fraction, from 0 to 1)
      
      real,    intent(in)  :: pres   ! surface air pressure 
     
   

  !     output arguments:

      real, intent(out) :: Idfuse  ! diffuse solar radiation (W/m^2)
      real, intent(out) :: Idrctt  ! total direct solar radiation (W/m^2)
  !
  !     internal arguments:
  !
  !        cn    - clearness number (defined as the ratio of normal
  !                incident radiation, under local mean water vapour,
  !                divided by the normal incident radiation, for
  !                water vapour in a basic atmosphere)  
  !              - currently, this value set equal to 1.0, but eventually
  !                may vary as a function of latitude and month pending further
  !                literature review.
  !        a     - solar constant at sea-level, varies by day (W/m^2)
  !        aday  - fixed values of a used in the table look up
  !        b     - inverse air mass, varies by day (atm^-1)
  !        bday  - fixed values of b used in the table look up
  !        pres0 - std sea-level pressure (101300 N/m^2)
  !        c     - constant which accounts for water vapour, varies by
  !                Julian day (unitless)
  !        cday  - fixed values of c used in the table look up
  !        iday  - fixed values of Julian day corresponding to aday,
  !                bday and cday
  !        dayinc - day increment used in interpolating between days

  !**********************************************************************
  !     declarations:


  real    ::  a, b, c, dayinc  ! As above
  real :: expa            ! exp. function with zenith angle and air mass
  real :: catten          ! cloud attenuation (unitless), from 0 to 1
  integer :: i
  integer, dimension(14), save ::   &
       nday = (/  1, 21, 52, 81,112,142,173, 203,234,265,295,326,356,366/)

  real, dimension(14), save ::   &
       aday = (/1203.0,1202.0,1187.0,1164.0,1130.0,1106.0,1092.0,   &
                1093.0,1107.0,1136.0,1136.0,1190.0,1204.0,1203.0/)  &
      ,bday = (/0.141,0.141,0.142,0.149,0.164,0.177,0.185,          &
                0.186,0.182,0.165,0.152,0.144,0.141,0.141/)         &
      ,cday = (/0.103,0.103,0.104,0.109,0.120,0.130,0.137,          &
                0.138,0.134,0.121,0.111,0.106,0.103,0.103/)

  real, save  :: cn = 1.0
  real, save  :: pres0 = 101300.0 

  real :: solar  !total solar radiation, diff.+direct (W/m^2)
  real :: Idrctn  ! direct normal solar radiation (W/m^2)



 if ( coszen  < 1.0e-15 ) then  ! night or too close to zen=90.0 
                                              ! for comfort!
     Idfuse = 0.0
     Idrctt = 0.0
     return  !  Nothing else to do!

  else !...compute radiation

     ! first, perform the table look up
      do i = 1, 14
        if (jday <=   nday(i)) then 
      exit    ! ds go to 20
      end if
      end  do
      if ( DEBUG .and. i < 1 .or. i > 14) then
        write(unit=6,fmt=*) "solbio: day index out of range"
      end if

      if (nday(i) == 1) then
        a = aday(1)
        b = bday(1)
        c = cday(1)
      else
        dayinc = real(jday-nday(i-1)) / real(nday(i)-nday(i-1))
        a = aday(i-1) + (aday(i)-aday(i-1))*dayinc
        b = bday(i-1) + (bday(i)-bday(i-1))*dayinc
        c = cday(i-1) + (cday(i)-cday(i-1))*dayinc
      end if

      expa = exp(-b*(pres/pres0)/coszen)
     

    !...compute cloud attenuation

      catten = 1.0 - 0.75*cl**3.4     !(source: Kasten & Czeplak (1980)) 

    !ds - catten applied first to Idrctn

      Idrctn  = cn*a*expa* catten
      Idfuse =  c*Idrctn
      Idrctt  = Idrctn*coszen

     !u7.lu solar  = Idrctt + Idfuse

    ! ds- removed
    !...compute absorbed radiation (W/m2)...
    !...the albedo is set constant equal to 0.23 ... life is hard #$@!
    ! real,    intent(in)  ::  albedo !  albedo (0..1, default = 0.23)
    !..   solar = solar *(1-albedo)

  end if ! daylight

  end subroutine SolBio

  !=============================================================================
  subroutine SolBio2D(jday,cl,pres)
  !=============================================================================
  !
  !
  !  function:
  !
  !     computes the radiation terms needed for stomatal calculations
  !     one term is computed: total solar radiation (W/m^2)
  !     methodology for this calculation taken from Iqbal, M., 1983,
  !     An introduction to solar radiation, Academic Press, New York, 
  !     pp. 202-210.
  !
  !  history:
  !
  !     modified for EMEP model
  !     Development of this routine was prompted by the need for a
  !     horizontal rather than an actinic flux calculation (which had
  !     been performed by Soleng). Furthermore, Soleng computed total
  !     radiation only out to the near-ir spectrum. This program
  !     is designed only for approximate radiation estimates to be used
  !     for stomatal calculations.
  !
  !     8/90    initial development of SolBio by T. Pierce
  !     9/95    modified by Hugo Jakobsen, 29/9-95 
  !     1-3/01  modified by Dave Simpson, with albedo as input
  !     and Idrct, Idfuse as output, and converted to F90
  !     19/9/01 - albedo removed, Solar thus gives total incoming
  !     ToDo
  !          - could consider alternative schemes for cloud attenuation,
  !      notably use of low, medium and high cloud (see e.g.  Stull (1988) 
  !          - check clearness number, pressure assumptions (too American?)
  !
  !  argument list description:
  !
      integer, intent(in)  :: jday   ! day no. from 1 Jan (Julian day)

      real,    intent(in), dimension(:,:) :: &
!ds               coszen   &! cos(zen) where zen is zenith angle
               cl       &! opaque sky cover (fraction, from 0 to 1)
              ,pres      ! surface air pressure 

  !     output arguments:
  !ds    real,    intent(out), dimension(:,:) :: &
  !ds                         Idfuse &! diffuse solar radiation (W/m^2)
  !ds                        ,Idrctt  ! total direct solar radiation (W/m^2)
  !
  !     internal arguments:
  !
  !        cn    - clearness number (defined as the ratio of normal
  !                incident radiation, under local mean water vapour,
  !                divided by the normal incident radiation, for
  !                water vapour in a basic atmosphere)  
  !              - currently, this value set equal to 1.0, but eventually
  !                may vary as a function of latitude and month pending further
  !                literature review.
  !        a     - solar constant at sea-level, varies by day (W/m^2)
  !        aday  - fixed values of a used in the table look up
  !        b     - inverse air mass, varies by day (atm^-1)
  !        bday  - fixed values of b used in the table look up
  !        pres0 - std sea-level pressure (101300 N/m^2)
  !        c     - constant which accounts for water vapour, varies by
  !                Julian day (unitless)
  !        cday  - fixed values of c used in the table look up
  !        iday  - fixed values of Julian day corresponding to aday,
  !                bday and cday
  !        dayinc - day increment used in interpolating between days

  !**********************************************************************
  !     declarations:


  real    ::  a, b, c, dayinc  ! As above
!ds  real :: expa            ! exp. function with zenith angle and air mass
!ds  real :: catten          ! cloud attenuation (unitless), from 0 to 1
  integer :: i
  integer, dimension(14), save ::   &
       nday = (/  1, 21, 52, 81,112,142,173, 203,234,265,295,326,356,366/)

  real, dimension(14), save ::   &
       aday = (/1203.0,1202.0,1187.0,1164.0,1130.0,1106.0,1092.0,   &
                1093.0,1107.0,1136.0,1136.0,1190.0,1204.0,1203.0/)  &
      ,bday = (/0.141,0.141,0.142,0.149,0.164,0.177,0.185,          &
                0.186,0.182,0.165,0.152,0.144,0.141,0.141/)         &
      ,cday = (/0.103,0.103,0.104,0.109,0.120,0.130,0.137,          &
                0.138,0.134,0.121,0.111,0.106,0.103,0.103/)

  real, save  :: cn = 1.0
  real, save  :: pres0 = 101300.0 

  real :: solar  !total solar radiation, diff.+direct (W/m^2)
!ds  real :: Idrctn  ! direct normal solar radiation (W/m^2)
!ds:
  real, dimension(size(coszen,1),size(coszen,2)) :: expa, catten, Idrctn

 ! first, perform the table look up
    do i = 1, 14
        if (jday <=   nday(i)) then 
            exit    ! ds go to 20
        end if
    end  do
    if ( DEBUG .and. i < 1 .or. i > 14) then
        write(unit=6,fmt=*) "solbio: day index out of range"
    end if

    if (nday(i) == 1) then
        a = aday(1)
        b = bday(1)
        c = cday(1)
    else
        dayinc = real(jday-nday(i-1)) / real(nday(i)-nday(i-1))
        a = aday(i-1) + (aday(i)-aday(i-1))*dayinc
        b = bday(i-1) + (bday(i)-bday(i-1))*dayinc
        c = cday(i-1) + (cday(i)-cday(i-1))*dayinc
    end if

 where ( coszen  < 1.0e-15 ) ! night or too close to zen=90.0 
                             ! for comfort!
     Idiffuse = 0.0
     Idirectt = 0.0

  else where !...compute radiation

      expa = exp(-b*(pres/pres0)/coszen)

    !...compute cloud attenuation

      catten = 1.0 - 0.75*cl**3.4     !(source: Kasten & Czeplak (1980)) 

    !ds - catten applied first to Idrctn

      Idrctn    = cn*a*expa* catten
      Idiffuse  = c*Idrctn
      Idirectt  = Idrctn*coszen

    ! ds- removed
    !...compute absorbed radiation (W/m2)...
    !...the albedo is set constant equal to 0.23 ... life is hard #$@!
    ! real,    intent(in)  ::  albedo !  albedo (0..1, default = 0.23)
    !..   solar = solar *(1-albedo)

  end where ! daylight

  end subroutine SolBio2D
!===============================================================
end module Radiation_ml
!===============================================================
