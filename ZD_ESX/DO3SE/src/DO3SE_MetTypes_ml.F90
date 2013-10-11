module DO3SE_MetTypes_ml

  use DO3SE_ModelConstants_ml, only: UNDEF, MAX_LAYERS

  implicit none
  private

  type, public :: MetData_t
    real :: Ts_C = UNDEF      !< Surface air temperature (degrees C)
    real :: P = UNDEF         !< Atmospheric pressure (kPa)
    real :: precip = UNDEF    !< Precipitations (mm)
    real :: ustar = UNDEF     !< Friction velocity over target canopy (m s-1)
    real :: u = UNDEF         !< Measured windspeed (m s-1)
    real :: O3 = UNDEF        !< Measured O3 concentration (ppb)

    real :: sinB = UNDEF      !< sin() of solar elevation angle
    real :: Rn = UNDEF        !< Net radiation (MJ m-2 h-1)
    real :: R = UNDEF         !< Global radiation (W m-2)
    real :: PAR = UNDEF       !< Photosynthetically active radiation (W m-2)
    real :: PPFD = UNDEF      !< Photosynthetic photon flux density (umol m-2 s-1)
    real :: Idrctt = UNDEF    !< Direct PAR irradiance (W m-2)
    real :: Idfuse = UNDEF    !< Diffuse PAR irradiance (W m-2)

    real :: VPD = UNDEF       !< Vapour pressure deficit (kPa)
    real :: esat = UNDEF      !< Saturated vapour pressure (kPa)
    real :: eact = UNDEF      !< Actual vapour pressure (kPa)
    real :: RH = UNDEF        !< Relative humidity (fraction)

    real :: Tleaf_C           !< Leaf temperature (degrees C)
    real :: CO2 = UNDEF       !< CO2 concentration (ppm)
  end type MetData_t

  type, public :: MetConfig_t
  end type MetConfig_t

contains
end module DO3SE_MetTypes_ml
