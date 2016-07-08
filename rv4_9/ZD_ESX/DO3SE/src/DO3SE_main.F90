program DO3SE_main
  use DO3SE_ModelConstants_ml
  use DO3SE_Types_ml
  use DO3SE_Util_ml
#include "DO3SE_Util_ml.h"
  use DO3SE_MetTypes_ml
  use DO3SE_Met_ml
  use DO3SE_GstoTypes_ml
  use DO3SE_Gsto_ml
  use DO3SE_Photosynthesis_ml
  use DO3SE_Phenology_ml
  use DO3SE_Resistance_ml
  use DO3SE_SMD_ml

  implicit none

  character(len=256) :: config_file
  integer :: config_unit

  type :: Config_t
    type(Location_t) :: location = Location_t()
    type(SMDConfig_t) :: SMD = SMDConfig_t()
    integer :: nL = 1
  end type

  type :: LandCover_t
    character(len=128) :: name = ""
    type(Season_t) :: season = Season_t()
    type(GstoConfig_t) :: gsto = GstoConfig_t()
    type(PnGstoConfig_t) :: pn_gsto = PnGstoConfig_t()
    real :: Lm = UNDEF    !< Leaf dimension (m)
  end type

  type(Config_t) :: conf
  type(LandCover_t), dimension(MAX_LC) :: LCs

  type(MetData_t) :: met
  type(SMDData_t) :: SMD

  type(GstoParams_t), dimension(MAX_LAYERS,MAX_LC) :: gsto_params
  type(ResistanceModel_t) :: rmodel
  integer :: nL, nLC

  integer :: dd   !< Day of year
  integer :: hr   !< Hour of day (0--23)
  real :: h       !< Canopy height
  real :: PARsun            !< PAR received by sunlit leaves (W m-2)
  real :: PARshade          !< PAR received by shaded leaves (W m-2)
  real, dimension(MAX_LAYERS,MAX_LC) :: LAI = UNDEF
  real, dimension(MAX_LAYERS,MAX_LC) :: LAIsunfrac = UNDEF
  real, dimension(MAX_LAYERS,MAX_LC) :: SAI = UNDEF
  real, dimension(MAX_LAYERS,MAX_LC) :: leaf_gsto = UNDEF
  real, dimension(MAX_LAYERS,MAX_LC) :: mean_gsto = UNDEF
  real, dimension(MAX_LAYERS,MAX_LC) :: bulk_gsto = UNDEF

  integer :: iL, iLC
  logical :: success

  ! Set up error handling
  assert => DO3SE_assert
  assert_not => inverted_assert

  ! Read configuration
  call get_command_argument(1, config_file)
  open (newunit=config_unit, file=config_file, status="old", action="read", position="rewind")
  call assert(read_DO3SE_Config(config_unit, conf), "Failed to read DO3SE_Config")
  nL = conf%nL  ! TODO: replace this with a better mechanism
  rewind (config_unit)
  nLC = 0
  do while (nLC < MAX_LC)
    nLC = nLC + 1
    success = read_DO3SE_LandCover(config_unit, LCs(nLC))
    if (.not. success) then
      nLC = nLC - 1
      exit
    end if
  end do
  close (config_unit)

  ! Sanitise/fill in configuration
  call check_SMDConfig(conf%SMD)
  do iLC = 1, nLC
    associate (LC => LCs(iLC))
      call check_Season(LC%season, nL)
      call check_GstoConfig(LC%gsto)
      if (LC%gsto%method == "photosynthesis") then
        call check_PnGstoConfig(LC%pn_gsto)
      end if
    end associate
  end do

  ! Initialise modules

  do
    exit

    ! Prepare for reading input row

    ! Read input row

    ! Fixup met data
    call met_data_fixup(met, conf%location, dd, hr)
    call PAR_sun_shade(met%Idrctt, met%Idfuse, met%sinB, &
                       sum(LCs(1:nLC)%gsto%cosA)/nLC, sum(LAI(1:nL,1:nLC)), &
                       PARsun, PARshade)

    ! Do daily P-M if necessary
    ! Do LAI, SAI, canopy height

    ! Ensure we have LAI and SAI values
    do iLC = 1, nLC
      associate (season => LCs(iLC)%season)
        select case (season%LAI_method)
        case ("input all")
          ! nothing to do
        case ("input total")
          ! Redistribute top layer LAI across all layers
          LAI(1:nL,iLC) = season%fLAI_layer(1:nL) * LAI(1,iLC)
        case ("day PLF total")
          ! Estimate total LAI and then redistribute across layers
          LAI(1,iLC) = LAI_day_PLF(season, dd)
          LAI(1:nL,iLC) = season%fLAI_layer(1:nL) * LAI(1,iLC)
        case default
          UNKNOWN_STRING(season%LAI_method)
        end select

        select case (season%SAI_method)
        case ("input all")
          ! nothing to do
        case ("input total")
          ! Redistribute top layer SAI across all layers
          SAI(1:nL,iLC) = season%fLAI_layer(1:nL) * SAI(1,iLC)
        case ("LAI")
          ! SAI = LAI
          SAI(1:nL,iLC) = LAI(1:nL,iLC)
        case ("forest")
          ! Forest SAI = LAI + 1, this is applied by splitting the "+ 1" across
          ! the layers according to the LAI distribution.
          SAI(1:nL,iLC) = LAI(1:nL,iLC) + 1.0 * season%fLAI_layer(1:nL)
        case ("wheat")
          ! Wheat SAI based on growing season - only valid if LAI is similarly
          ! estimated
          if (season%LAI_method == "day PLF total") then
            ! Estimate SAI from total LAI and then redistribute across layers
            SAI(1,iLC) = SAI_wheat(season, dd, sum(LAI(1:nL,iLC)))
            SAI(1:nL,iLC) = season%fLAI_layer(1:nL) * SAI(1,iLC)
          else
            ERROR('season%SAI_method="wheat" requires season%LAI_method="day PLF total"')
          end if
        case default
          UNKNOWN_STRING(season%SAI_method)
        end select
      end associate
    end do

    ! Estimate sunlit LAI (used later in f_light)
    call MLMC_sunlit_LAI(LAI(1:nL,1:nLC), met%sinB, LAIsunfrac(1:nL,1:nLC))

    ! Calculate gsto
    do iL = 1, nL
      do iLC = 1, nLC
        associate (LC => LCs(iLC), &
                   season => LCs(iLC)%season, &
                   gc => LCs(iLC)%gsto, &
                   gp => gsto_params(iL,iLC))
          ! Initialise gsto parameters from configuration
          gp = GstoParams_t(fmin=gc%fmin, gmax=gc%gmax, gmorph=gc%gmorph)

          ! Calculate f_phen
          select case (gc%f_phen_method)
          case ("disabled")
            ! Nothing to do
          case ("simple day PLF")
            gp%f_phen = f_phen_simple_PLF(gc, season%SGS, season%EGS, dd)
          case ("complex day PLF")
            gp%f_phen = f_phen_complex_PLF(gc, season%SGS, season%EGS, dd)
          case default
            UNKNOWN_STRING(gc%f_phen_method)
          end select

          ! Calculate leaf_f_phen
          select case (gc%leaf_f_phen_method)
          case ("disabled")
            ! Nothing to do
          case ("f_phen")
            gp%leaf_f_phen = gp%f_phen
          case ("day PLF")
            gp%leaf_f_phen = leaf_f_phen_PLF(gc, season%Astart, season%Aend, dd)
          case default
            UNKNOWN_STRING(gc%leaf_f_phen_method)
          end select

          select case (gc%f_light_method)
          case ("disabled")
            ! Nothing to do
          case ("enabled")
            ! Calculate f_light and leaf_f_light
            ! TODO: attenuate PAR properly through the canopy
            gp%f_light = f_light(gc%f_lightfac, PARsun, PARshade, LAIsunfrac(iL,iLC))
            ! TODO: "grassland multilayer" model used leaf_flight = Flightsun, i.e.
            !       leaf_f_light(gc%f_lightfac, met%PARsun) - which version is right?
            gp%leaf_f_light = leaf_f_light(gc%f_lightfac, met%PAR)
          case default
            UNKNOWN_STRING(gc%f_light_method)
          end select

          ! Calculate f_temp
          select case (gc%f_temp_method)
          case ("disabled")
            ! Nothing to do
          case ("default")
            gp%f_temp = f_temp(met%Ts_C, gc%T_min, gc%T_opt, gc%T_max, gc%fmin)
          case ("square high")
            gp%f_temp = f_temp_square_high(met%Ts_C, gc%T_min, gc%T_opt, gc%T_max, gc%fmin)
          case default
            UNKNOWN_STRING(gc%f_temp_method)
          end select

          ! Calculate f_VPD
          select case (gc%f_VPD_method)
          case ("disabled")
            ! Nothing to do
          case ("linear")
            gp%f_VPD = f_VPD_linear(met%VPD, gc%VPD_max, gc%VPD_min, gc%fmin)
          case ("log")
            gp%f_VPD = f_VPD_log(met%VPD, gc%fmin)
          case default
            UNKNOWN_STRING(gc%f_VPD_method)
          end select

          ! Calculate f_SW
          select case (gc%f_SW_method)
          case ("disabled")
            ! Nothing to do
          case ("fSWP exp")
            gp%f_SW = f_SWP_exp(gc%fSWP_exp_a, gc%fSWP_exp_b, gc%fmin, SMD%SWP)
          case ("fSWP linear")
            gp%f_SW = f_SWP_linear(gc%SWP_min, gc%SWP_max, gc%fmin, SMD%SWP)
          ! TODO: implement LWP
          !case ("fLWP exp")
          !  gp%f_SW = f_SWP_exp(gc%fSWP_exp_a, gc%fSWP_exp_b, gc%fmin, SMD%LWP)
          case ("fPAW")
            gp%f_SW = f_PAW(conf%SMD%ASW_FC, gc%fmin, SMD%ASW)
          case default
            UNKNOWN_STRING(gc%f_SW_method)
          end select

          ! Calculate gsto
          select case (gc%method)
          case ("multiplicative")
            leaf_gsto(iL,iLC) = gsto_leaf(gp)
            mean_gsto(iL,iLC) = gsto_mean(gp)
          case ("photosynthesis")
            ! Hybrid photosynthesis + multiplicative model
            gp%gmax = gsto_pn(LC%pn_gsto, met%Ts_C, met%Tleaf_C, met%u, met%CO2, met%RH, met%PPFD, LC%Lm)
            leaf_gsto(iL,iLC) = gsto_leaf(gp)
            mean_gsto(iL,iLC) = gsto_mean(gp)
          case default
            UNKNOWN_STRING(gc%method)
          end select
          ! Scale mean gsto up to bulk gsto
          bulk_gsto(iL,iLC) = mean_gsto(iL,iLC) * LAI(iL,iLC)
        end associate
      end do
    end do

    ! Calculate resistance model
    rmodel = ResistanceModel_t()
    !rmodel%Ra = ...
    !rmodel%Rb = ...
    do iL = 1, nL
      rmodel%Rinc(iL) = Rinc_prototype(sum(SAI(iL,1:nLC)), met%ustar)
      rmodel%Rext(iL) = Rext(sum(SAI(iL,1:nLC)))
      !rmodel%Rsto(iL) = MC_Rsto(gsto_params(iL,1:nLC), LAI(iL,1:nLC))
    end do
    rmodel%Rsur(1:nL) = ML_Rsur(rmodel%Rsto(1:nL), rmodel%Rext(1:nL), &
                                rmodel%Rinc(1:nL), conf%location%Rsoil)

    ! Calculate O3 deposition
    ! Do hourly P-M
    ! Write output row
  end do
contains

  !> Read a land cover definition from a namelist file.  Returns .true. or
  !! .false. depending on whether or not the read succeeded.
  logical function read_DO3SE_LandCover(unit, LC)
    integer, intent(in) :: unit
    type(LandCover_t), intent(out), target :: LC

    character(len=128), pointer :: name
    type(Season_t), pointer :: season
    type(GstoConfig_t), pointer :: gsto
    type(PnGstoConfig_t), pointer :: pn_gsto
    real, pointer :: Lm
    namelist /DO3SE_LandCover/ LC, name, season, gsto, pn_gsto, Lm

    integer :: ios

    LC = LandCover_t()
    name => LC%name
    season => LC%season
    gsto => LC%gsto
    pn_gsto => LC%pn_gsto
    Lm => LC%Lm

    read (unit=unit, nml=DO3SE_LandCover, iostat=ios)
    ! Set true/false exit status based on iostat
    read_DO3SE_LandCover = ios == 0
  end function read_DO3SE_LandCover

  !> Read general configuration.
  logical function read_DO3SE_Config(unit, conf)
    integer, intent(in) :: unit
    type(Config_t), intent(out), target :: conf

    type(Location_t), pointer :: location
    type(SMDConfig_t), pointer :: SMD

    namelist /DO3SE_Config/ location, SMD

    integer :: ios

    conf = Config_t()
    location => conf%location
    SMD => conf%SMD

    read (unit=unit, nml=DO3SE_Config, iostat=ios)
    read_DO3SE_Config = ios == 0
  end function read_DO3SE_Config

  !> Calculate the sunlit LAI from total LAI and solar elevation angle.
  !!
  !! TODO: (2 * sinB) should be sinB/cosA
  pure real function sunlit_LAI(LAI, sinB)
    real, intent(in) :: LAI     !< Leaf area index (m^2/m^2)
    real, intent(in) :: sinB    !< sin() of solar elevation angle

    sunlit_LAI = ((1 - exp(-0.5 * LAI / sinB)) * (2 * sinB))
  end function sunlit_LAI

  !> Multi-layer multi-component sunlit LAI model.  LAI and LAIsunfrac must be
  !! the same dimension(nL,nLC).
  pure subroutine MLMC_sunlit_LAI(LAI, sinB, LAIsunfrac)
    real, dimension(:,:), intent(in) :: LAI  !< Leaf area index (m^2/m^2)
    real, intent(in) :: sinB                 !< sin() of solar elevation angle

    real, dimension(size(LAI,1),size(LAI,2)), intent(out) :: LAIsunfrac  !< Fraction of LAI that is sunlit

    real, dimension(0:size(LAI,1)) :: sunLAI_acc
    real :: sunLAI_layer
    integer :: iL

    sunLAI_acc(0) = 0.0
    do iL = 1, size(LAI, 1)
      ! How much of "canopy so far" is sunlit?
      sunLAI_acc(iL) = sunlit_LAI(sum(LAI(1:iL,:)), sinB)
      ! How much of that is in this layer?
      sunLAI_layer = sunLAI_acc(iL) - sunLAI_acc(iL - 1)
      ! Fraction of LAI which is sunlit
      LAIsunfrac(iL,:) = sunLAI_layer / sum(LAI(iL,:))
    end do
  end subroutine MLMC_sunlit_LAI

end program DO3SE_main
