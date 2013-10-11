module DO3SE_SMD_ml

  use DO3SE_ModelConstants_ml, only: UNDEF
  use DO3SE_Resistance_ml, only: ResistanceModel_t
  use DO3SE_Util_ml
#include "DO3SE_Util_ml.h"

  implicit none
  private

  !> SMD soil texture parameters.
  type, public :: Soil_t
    real :: b = UNDEF         !< Texture dependent soil conductivity parameter
    real :: FC = UNDEF        !< Field capacity (m3 m-3)
    real :: SWP_AE = UNDEF    !< Water potential at air entry (MPa)
    real :: Ksat = UNDEF      !< Saturated soil conductance (s-2 MPa-1 mm-1)
  end type Soil_t

  ! Commonly used soil textures:
  !> Sandy loam
  type(Soil_t), parameter, public :: SOIL_SANDY_LOAM = &
    Soil_t(b = 3.31, FC = 0.16, SWP_AE = -0.00091, Ksat = 0.0009576)
  !> Silt loam
  type(Soil_t), parameter, public :: SOIL_SILT_LOAM = &
    Soil_t(b = 4.38, FC = 0.26, SWP_AE = -0.00158, Ksat = 0.0002178)
  !> Loam
  type(Soil_t), parameter, public :: SOIL_LOAM = &
    Soil_t(b = 6.58, FC = 0.29, SWP_AE = -0.00188, Ksat = 0.0002286)
  !> Clay loam (Ksat estimated)
  type(Soil_t), parameter, public :: SOIL_CLAY_LOAM = &
    Soil_t(b = 7.00, FC = 0.37, SWP_AE = -0.00588, Ksat = 0.00016)

  type, public :: SMDConfig_t
    !> Soil texture:
    !!    - "sandy loam"
    !!    - "silt loam"
    !!    - "loam"
    !!    - "clay loam"
    !!    - "" or "custom":   Soil texture parameters set individually
    character(len=16) :: soil_texture = "loam"
    !> Soil texture parameters
    type(Soil_t) :: soil = Soil_t()

    real :: root = 1.2        !< Root depth (m)
    real :: PWP = -4.0        !< "Permanent wilting point", minimum level for SWP (MPa)
    real :: ASW_FC = UNDEF    !< Available soil water at field capacity (m)
                              !! (calculated by check_SMDConfig)

    !> Soil water data source:
    !!    - "disabled":   No soil water data
    !!    - "input SWP":  Input soil water potential (MPa)
    !!    - "input SWC":  Input soil water volumetric content (m3 m-3)
    !!    - "P-M":        Use Penman-Monteith method to track soil water content
    character(len=16) :: source = "disabled"
    ! Penman-Monteith method parameters
    real :: PM_initial_SWC = UNDEF    !< Initial soil water content, defaults to soil%FC
  end type SMDConfig_t

  type, public :: SMDData_t
    real :: Sn = UNDEF        !< Soil water content (m3 m-3)
    real :: SWP = UNDEF       !< Soil water potential (MPa)
    real :: ASW = UNDEF       !< Available soil water (m)
    real :: SMD = UNDEF       !< Soil moisture deficit (m)
  end type SMDData_t

  type, public :: PM_State_t
    ! Latest evapotranspiration values
    real :: Ei = 0.0          !< Evaporation from canopy (m)
    real :: Et = 0.0          !< Plant transpiration (m)
    real :: Es = 0.0          !< Soil surface evaporation (m)
    real :: Eat = 0.0         !< Evapotranspiration (m)

    ! Accumulated values (processed and cleared daily)
    real :: Ei_acc = 0.0      !< Accumulated evaporation from canopy (m)
    real :: Et_acc = 0.0      !< Accumulated plant transpiration (m)
    real :: Es_acc = 0.0      !< Accumulated soil surface evaporation (m)
    real :: Eat_acc = 0.0     !< Accumulated evapotranspiration (m)
    real :: precip_acc = 0.0  !< Accumulated precipitation (m)

    ! Soil water content tracking
    real :: Sn = UNDEF        !< Volumetric soil water content (m3 m-3)
    real :: Sn_diff = 0.0     !< Latest change in Sn (m3 m-3)
  end type PM_State_t

  public :: check_SMDConfig
  public :: soil_moisture_from_SWC
  public :: soil_moisture_from_SWP
  public :: penman_monteith_hourly
  public :: penman_monteith_daily

contains

  !> Check that a Soil_t is valid.  The optional *name* argument can be used to
  !! select one of the preset values, passing name="" or name="custom" has the
  !! same effect as omitting the argument.
  subroutine check_Soil(soil, name)
    type(Soil_t), intent(inout) :: soil
    character(len=*), intent(in), optional :: name

    if (present(name)) then
      select case (name)
      case ("sandy loam")
        soil = SOIL_SANDY_LOAM
      case ("silt loam")
        soil = SOIL_SILT_LOAM
      case ("loam")
        soil = SOIL_LOAM
      case ("clay loam")
        soil = SOIL_CLAY_LOAM
      case ("", "custom")
        ! Parameters should have been set manually
      case default
        ERROR("unknown soil texture: "//trim(name))
      end select
    end if

    ! Check that parameters were set by the name or manually
    ASSERT_DEFINED(soil%b)
    ASSERT_DEFINED(soil%FC)
    ASSERT_DEFINED(soil%SWP_AE)
    ASSERT_DEFINED(soil%Ksat)
  end subroutine check_Soil

  !> Check that a SMDConfig_t is valid.
  !!
  !!    - Checks that the Soil_t is valid
  !!    - Applies named soil texture if necessary
  !!    - Calculates various soil-dependent constants
  subroutine check_SMDConfig(config)
    type(SMDConfig_t), intent(inout) :: config

    type(SMDData_t) :: tmp_data

    ! Check/set soil texture parameters
    call check_Soil(config%soil, config%soil_texture)

    ! Calculate available soil water at field capacity
    tmp_data = SMDData_t()
    call soil_moisture_from_SWC(config, tmp_data, config%soil%FC)
    config%ASW_FC = tmp_data%ASW

    ! Set up soil water data source parameters
    if (config%source == "P-M") then
      if (.not. is_def(config%PM_initial_SWC)) then
        config%PM_initial_SWC = config%soil%FC
      end if
    end if
  end subroutine check_SMDConfig

  !> Use soil water release curve to convert from soil water content (m3 m-3) to
  !! soil water potential (MPa).
  pure real function SWC_to_SWP(SWP_AE, b, SWC)
    real, intent(in) :: SWP_AE    !< Water potential at air entry (MPa)
    real, intent(in) :: b         !< Texture dependent soil conductivity parameter
    real, intent(in) :: SWC       !< Soil water content (m3 m-3)

    real, parameter :: SWC_sat = 0.4 ! Saturated soil water content for soil water release curve

    SWC_to_SWP = SWP_AE * ((SWC_sat / SWC)**b)
  end function SWC_to_SWP

  !> Use soil water release curve to convert from soil water potential (MPa) to
  !! soil water content (m3 m-3).
  pure real function SWP_to_SWC(SWP_AE, b, SWP)
    real, intent(in) :: SWP_AE    !< Water potential at air entry (MPa)
    real, intent(in) :: b         !< Texture dependent soil conductivity parameter
    real, intent(in) :: SWP       !< Soil water potential (MPa)

    real, parameter :: SWC_sat = 0.4 ! Saturated soil water content for soil water release curve

    SWP_to_SWC = 1.0 / (((SWP/SWP_AE)**(1.0/b)) / SWC_sat)
  end function SWP_to_SWC

  !> Fill soil moisture data from a soil water content value.
  pure subroutine soil_moisture_from_SWC(config, data, Sn)
    type(SMDConfig_t), intent(in) :: config
    type(SMDData_t), intent(inout) :: data
    real, intent(in) :: Sn      !< Soil water content (m3 m-3)

    real :: PWP_vol

    associate (soil => config%soil)
      ! Convert PWP to volumetric content to use as a minimum soil water content
      PWP_vol = SWP_to_SWC(soil%SWP_AE, soil%b, config%PWP)

      ! Constrain soil water content to be between field capacity and PWP
      data%Sn = max(PWP_vol, min(soil%FC, Sn))

      ! Calculate soil water potential (SWP)
      data%SWP = SWC_to_SWP(soil%SWP_AE, soil%b, data%Sn)

      ! Calculate available soil water (ASW)
      data%ASW = (data%Sn - PWP_vol) * config%root

      ! Calculate soil moisture deficit (SMD)
      data%SMD = (soil%FC - data%Sn) * config%root
    end associate
  end subroutine soil_moisture_from_SWC

  !> Fill soil moisture data from a soil water potential value.
  pure subroutine soil_moisture_from_SWP(config, data, SWP)
    type(SMDConfig_t), intent(in) :: config
    type(SMDData_t), intent(inout) :: data
    real, intent(in) :: SWP     !< Soil water potential (MPa)

    associate (soil => config%soil)
      call soil_moisture_from_SWC(config, data, SWP_to_SWC(soil%SWP_AE, soil%b, SWP))
    end associate
  end subroutine soil_moisture_from_SWP


  !> Hourly Penman-Monteith calculations for evaporation and transpiration.  
  !! Fills in values and updates accumulators in *state* for Ei, Et, Es and Eat.
  !!
  !! TODO: should we still be using Rsoil, instead of rm%Rgs?
  subroutine penman_monteith_hourly(Rn, P, Ts_C, VPD, rm, LAI, Es_blocked, state)
    real, intent(in) :: Rn              !< Net radiation (J)
    real, intent(in) :: P               !< Atmospheric pressure (P)
    real, intent(in) :: Ts_C            !< Air temperature (degrees C)
    real, intent(in) :: VPD             !< Vapour pressure deficit (Pa)
    type(ResistanceModel_t), intent(in) :: rm
    real, intent(in) :: LAI             !< Leaf area index (m2 m-2)
    logical, intent(in) :: Es_blocked   !< Is soil evaporation blocked?

    type(PM_State_t), intent(inout) :: state

    real, parameter :: Ts_K = 273.16    ! Conversion from degrees Celsius to Kelvin
    real, parameter :: Dratio = 0.663   ! Ratio between molecular diffusivity of O3 and H2O

    real :: esat, eact, Tvir, delta, lambda, psychro, Pair, Cair, G

    real :: Et_1, Et_2, Ei_3, Et_3
    real :: t, Es_Rn, Es_G, Es_1, Es_2, Es_3
    real :: SW_a, SW_s, SW_c, C_canopy, C_soil

    ! This model (probably) makes some one-layer assumptions, so don't allow 
    ! multi-layer resistance model.
    ASSERT(rm%nL == 1)


    associate (Ra => rm%Ra, &
               Rb_H2O => rm%Rb, &
               Rinc => rm%Rinc(1), &
               Rsto_c => rm%Rsto(1), &
               Rsoil => rm%Rgs, &
               Ei => state%Ei, &
               Et => state%Et, &
               Es => state%Es, &
               Eat => state%Eat)

      esat = 611 * exp(17.27 * Ts_C / (Ts_C + 237.3))
      eact = esat - VPD
      Tvir = (Ts_c+Ts_K)/(1-(0.378*(eact/P)))
      delta= ((4098*esat)/((Ts_C+237.3)**2))
      lambda = (2501000-(2361*Ts_C))                ! approx. 2.5e6
      psychro = 1628.6 * (P/lambda)                 ! 66.1 Pa/K at 100 kPa, 20C
      Pair = (0.003486*(P/Tvir))
      Cair = (0.622*((lambda*psychro)/P))

      G = 0.1 * Rn

      Et_1 = (delta * (Rn - G)) / lambda
      Et_2 = 3600 * Pair * Cair * VPD / Rb_H2O / lambda

      Ei_3 = delta + psychro
      Ei = (Et_1 + Et_2) / Ei_3 / 1000

      Et_3 = delta + psychro * (1 + (Rsto_c * Dratio) / Rb_H2O)
      Et = (Et_1 + Et_2) / Et_3 / 1000

      if (Es_blocked) then
        Es = 0
      else
        t = exp(-0.5 * LAI)
        Es_Rn = Rn * t
        Es_G = 0.1 * Es_Rn
        Es_1 = (delta * (Rn - G)) / lambda
        Es_2 = ((3600 * Pair * Cair * VPD) - (delta * Rinc * ((Rn - G) - (Es_Rn - Es_G)))) / (Rinc + Rb_H2O) / lambda
        Es_3 = delta + (psychro * (1.0 + (Rsoil / (Rb_H2O + Rinc))))
        Es = (Es_1 + Es_2) / Es_3 / 1000
      end if

      ! Calculate Eat from Et and Es (after Shuttleworth and Wallace, 1985)
      SW_a = (delta + psychro) * (Ra + Rb_H2O)
      SW_s = (delta + psychro) * Rinc + (psychro * Rsoil)
      SW_c = (delta + psychro) * 0 + (psychro + Rsto_c) ! Boundary layer 
      ! resistance = 0
      C_canopy = 1 / (1 + ((SW_c * SW_a) / (SW_s * (SW_c + SW_a))))
      C_soil = 1 / (1 + ((SW_s * SW_a) / (SW_c * (SW_s + SW_a))))
      if (Es <= 0) then
        Eat = Et
      else
        Eat = (C_canopy * Et) + (C_soil * Es)
      end if

    end associate

    ! Accumulate values
    state%Ei_acc = state%Ei_acc + state%Ei
    state%Et_acc = state%Et_acc + state%Et
    state%Es_acc = state%Es_acc + state%Es
    state%Eat_acc = state%Eat_acc + state%Eat
  end subroutine penman_monteith_hourly

  !> Daily soil water content update from accumulated Penman-Monteith values.  
  !! Sets state%Sn_diff, updates state%Sn and clears accumulators.  The soil 
  !! water content value isn't constrained in any way here, so this step must be 
  !! done by whatever calls this subroutine.
  pure subroutine penman_monteith_daily(state, LAI, root)
    type(PM_State_t), intent(inout) :: state
    real, intent(in) :: LAI     !< Leaf area index (m2 m-2)
    real, intent(in) :: root    !< Root depth (m)

    real :: P_input

    ! Estimate precipitation that recharges the soil water
    ! The first (0.0001*LAI) of precipitation is assumed to be intercepted, and then
    ! the evaporation (Ei) is applied to this amount.  Soil water cannot be lost through Ei.
    if (state%precip_acc > 0) then
      P_input = (state%precip_acc - (0.0001*LAI)) + max(0.0, (0.0001*LAI) - state%Ei_acc)
      P_input = max(0.0, P_input)
    else
      P_input = 0.0
    end if

    ! Total water change = incoming not evaporated - evapotranspiration from plant and soil
    ! Converted to volumetric change using root depth
    state%Sn_diff = (P_input - state%Eat_acc) / root
    state%Sn = state%Sn + state%Sn_diff

    ! Clear accumulated variables
    state%Ei_acc = 0.0
    state%Et_acc = 0.0
    state%Es_acc = 0.0
    state%Eat_acc = 0.0
    state%precip_acc = 0.0
  end subroutine penman_monteith_daily

end module DO3SE_SMD_ml
