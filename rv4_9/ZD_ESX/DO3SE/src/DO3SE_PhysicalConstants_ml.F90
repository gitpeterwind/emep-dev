module DO3SE_PhysicalConstants_ml

  implicit none
  public

  real, parameter :: PI = 3.141592653589793238
  real, parameter :: DEG2RAD = 0.017453292519943295

  ! Atmospheric pressure at sea level (kPa)
  real, parameter :: seaP = 101.325

  ! Degrees Celcius to degrees Kelvin
  real, parameter :: C2K = 273.15

  ! Approximate fraction of global radiation in PAR waveband
  real, parameter :: PARfrac = 0.45
  ! PAR conversion from W m-2 to umol photons m-2 s-1
  real, parameter :: PAR_Wm2_to_photons = 4.57

contains
end module DO3SE_PhysicalConstants_ml
