module DO3SE_Met_ml

  use DO3SE_PhysicalConstants_ml, only: DEG2RAD, seaP, C2K, PI, PARfrac, PAR_Wm2_to_photons
  use DO3SE_ModelConstants_ml, only: UNDEF
  use DO3SE_Types_ml
  use DO3SE_MetTypes_ml

  implicit none
  private

  public :: met_data_fixup
  public :: PAR_sun_shade

contains

  pure subroutine met_data_fixup(met, loc, dd, hr)
    type(MetData_t), intent(inout) :: met
    type(Location_t), intent(in) :: loc
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)

    ! Humidity
    met%esat = saturated_vapour_pressure(met%Ts_C)
    if (met%VPD == UNDEF .and. met%RH == UNDEF) then
      ! TODO: error if neither VPD or RH supplied
    else if (met%RH == UNDEF) then
      ! Calculate relative humidity from VPD
      met%eact = met%esat - met%VPD
      met%RH = met%eact / met%esat
    else if (met%VPD == UNDEF) then
      ! Calculate VPD from relative humidity
      met%eact = met%esat * met%RH
      met%VPD = met%esat - met%eact
    end if

    ! Irradiance
    if (met%PAR == UNDEF) then
      ! Calculate PAR from some other available source
      if (met%Idrctt /= UNDEF .and. met%Idfuse /= UNDEF) then
        ! PAR is sum of direct and diffuse PAR
        met%PAR = met%Idrctt + met%Idfuse
      else if (met%PPFD /= UNDEF) then
        ! PAR from PPFD
        met%PAR = met%PPFD / PAR_Wm2_to_photons
      else if (met%R /= UNDEF) then
        ! Estimate PAR from global radiation
        met%PAR = met%R * PARfrac
      else
        ! TODO: error if no source of PAR is available
      end if
    end if
    if (met%PPFD == UNDEF) then
      ! If PPFD is missing, convert from PAR
      met%PPFD = met%PAR * PAR_Wm2_to_photons
    end if
    if (met%Idrctt == UNDEF .or. met%Idfuse == UNDEF) then
      ! If direct + diffuse are missing, estimate from PAR
      call PAR_direct_diffuse(met%PAR, met%sinB, met%P, met%Idrctt, met%Idfuse)
    end if
    if (met%R == UNDEF) then
      ! If global radiation is absent, estimate from PAR
      met%R = met%PAR / PARfrac
    end if

    ! Solar elevation
    met%sinB = solar_elevation(loc%lat, loc%lon, dd, hr)

    ! Net radiation
    if (met%Rn == UNDEF) then
      met%Rn = net_radiation(loc%lat, loc%lon, loc%elev, loc%albedo, &
                             dd, hr, met%sinB, met%R, met%Ts_C, met%eact)
    end if
  end subroutine met_data_fixup


  !> Calculate saturated vapour pressure (kPa).
  pure real function saturated_vapour_pressure(Ts_C)
    real, intent(in) :: Ts_C    !< Surface air temperature (degrees C)

    saturated_vapour_pressure = 0.611 * exp(17.27 * Ts_C / (Ts_C + 237.3))
  end function saturated_vapour_pressure


  pure function solar_elevation(lat, lon, dd, hr) result(sinB)
    real, intent(in) :: lat     !< Latitude (degrees North)
    real, intent(in) :: lon     !< Longitude (degrees East)
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)

    real :: sinB      ! Output: sin() of solar elevation angle

    real :: t0, dec, h

    t0 = solar_noon(lon, dd)
    dec = solar_declination(dd)

    ! Hour-angle of the sun
    h = (15 * (hr - t0)) * DEG2RAD

    ! sin() of solar elevation angle
    sinB = sin(lat * DEG2RAD)*sin(dec) + cos(lat * DEG2RAD)*cos(dec)*cos(h)
    ! TODO: should this line be removed? what effect does it have?  Does any 
    !       use of sinB happen when sinB < 0?
    sinB = max(0.0, sinB)
  end function solar_elevation


  !> Calculate solar noon (hour of day).
  pure real function solar_noon(lon, dd)
    real, intent(in) :: lon     !< Longitude (degrees East)
    integer, intent(in) :: dd   !< Day of year

    real :: f, e, lonm, LC

    ! Solar noon correction for day of year
    f = (279.575 + (0.9856 * dd)) * DEG2RAD
    e = (-104.7*sin(f) + 596.2*sin(2*f) + 4.3*sin(3*f) - 12.7*sin(4*f) &
        - 429.3*cos(f) - 2.0*cos(2*f) + 19.3*cos(3*f)) / 3600
    ! Calculate the longitudinal meridian
    lonm = nint(lon / 15.0) * 15.0
    ! Solar noon, with day of year and longitudinal correction
    LC = (lon - lonm) / 15
    solar_noon = 12 - LC - e
  end function solar_noon


  !> Calculate solar declination (radians).
  pure real function solar_declination(dd)
    integer, intent(in) :: dd   !< Day of year

    solar_declination = (-23.4 * cos((360 * ((dd + 10) / 365.0))*DEG2RAD))*DEG2RAD
  end function solar_declination


  !> Estimate diffuse and direct PAR components.
  pure subroutine PAR_direct_diffuse(PAR, sinB, P, Idrctt, Idfuse)
    real, intent(in) :: PAR       !< Photosynthetically active radiation (W m-2)
    real, intent(in) :: sinB      !< sin() of solar elevation angle
    real, intent(in) :: P         !< Atmospheric pressure (kPa)

    real, intent(out) :: Idrctt   !< Direct PAR irradiance (W m-2)
    real, intent(out) :: Idfuse   !< Diffuse PAR irradiance (W m-2)

    real :: m, pPARdir, pPARdif, pPARtotal, ST, fPARdir, fPARdif

    m = 1.0 / sinB

    ! Potential direct PAR
    pPARdir = 600 * exp(-0.185 * (P / seaP) * m) * sinB
    ! Potential diffuse PAR
    pPARdif = 0.4 * (600 - pPARdir) * sinB
    ! Potential total PAR
    pPARtotal = pPARdir + pPARdif

    ! Sky transmissivity
    ! TODO: this differs from documentation
    ST = max(0.21, min(0.9, PAR / pPARtotal))

    ! Direct and diffuse fractions
    fPARdir = (pPARdir / pPARtotal) * (1.0 - ((0.9 - ST) / 0.7)**(2.0/3.0))
    fPARdif = 1 - fPARdir

    ! Apply calculated direct and diffuse fractions to PARtotal
    ! TODO: warning: f_light calculations expect PPFD(?)
    Idrctt = fPARdir * PAR
    Idfuse = fPARdif * PAR
  end subroutine PAR_direct_diffuse


  !> Estimate net radiation (MJ m-2 h-1).
  pure real function net_radiation(lat, lon, elev, albedo, dd, hr, sinB, R, Ts_C, eact)
    real, intent(in) :: lat     !< Latitude (degrees North)
    real, intent(in) :: lon     !< Longitude (degrees East)
    real, intent(in) :: elev    !< Elevation (m above sea level)
    real, intent(in) :: albedo  !< Surface albedo (fraction)
    integer, intent(in) :: dd   !< Day of year (1--365)
    integer, intent(in) :: hr   !< Hour of day (0--23)
    real, intent(in) :: sinB    !< sin() of solar elevation angle
    real, intent(in) :: R       !< Global radiation (W m-2)
    real, intent(in) :: Ts_C    !< Surface air temperature (degrees C)
    real, intent(in) :: eact    !< Actual vapour pressure (kPa)

    real, parameter :: Gsc = 0.082          ! Solar constant (MJ/m^2/min)
    real, parameter :: SBC = 4.903e-9 / 24  ! Stephan Boltzman constant

    real :: lat_rad, R_MJ, t0, h, h1, h2, dr, dec, Re, pR, Rnl, Rns

    if (sinB <= 0) then
      net_radiation = 0.0
    else
      ! Latitude in radians
      lat_rad = lat * DEG2RAD

      ! Convert global radiation W m-2 to MJ m-2 s-1
      R_MJ = R * 0.0036

      ! Hour-angle of the sun
      t0 = solar_noon(lon, dd)
      h = (15 * (hr - t0)) * DEG2RAD
      h1 = h - (PI/24)
      h2 = h + (PI/24)

      dr = 1 + (0.033 * cos(((2 * PI) / 365) * dd))
      dec = solar_declination(dd)
      ! External radiation (with fix to stop div by zero)
      ! TODO: fix this to be less hackish
      Re = max(0.00000000001, &
               ((12 * 60) / PI) * Gsc * dr * ((h2 - h1) * sin(lat_rad) * sin(dec) &
               + cos(lat_rad) * cos(dec) * (sin(h2) - sin(h1))))
      ! TODO: what was this for?
      !Re = max(0.0, ((12*60)/PI)*Gsc*dr*sinB)

      ! Calculate net longwave radiation
      pR = (0.75 + (2e-5 * elev)) * Re

      Rnl = max(0.0, (SBC*((Ts_C + C2K)**4)) * (0.34-(0.14*sqrt(eact))) &
                     * ((1.35*(min(1.0, R_MJ/pR)))-0.35))
      Rns = (1 - albedo) * R_MJ

      net_radiation = max(0.0, Rns - Rnl)
    end if
  end function net_radiation


  !> Estimate PAR received by sun and shade leaves within the canopy.
  pure subroutine PAR_sun_shade(Idrctt, Idfuse, sinB, cosA, LAI, PARsun, PARshade)
    real, intent(in) :: Idrctt        !< Direct PAR irradiance (W m-2)
    real, intent(in) :: Idfuse        !< Diffuse PAR irradiance (W m-2)
    real, intent(in) :: sinB          !< sin() of solar elevation angle
    real, intent(in) :: cosA          !< cos(A), A = mean leaf inclination (0.5 = 60 degrees)
    real, intent(in) :: LAI           !< Leaf area index (m2 m-2)

    real, intent(out) :: PARsun       !< PAR received by sunlit leaves (W m-2)
    real, intent(out) :: PARshade     !< PAR received by shaded leaves (W m-2)

    ! PAR flux densities evaluated using method of Norman (1982, p.79):
    ! "conceptually, 0.07 represents a scattering coefficient"
    PARshade = Idfuse * exp(-0.5 * LAI**0.8) + &
      0.07 * Idrctt * (1.1 - (0.1 * LAI)) * exp(-sinB)
    PARsun = Idrctt * 0.8 * (cosA / sinB) + PARshade
  end subroutine PAR_sun_shade

end module DO3SE_Met_ml
