module esx_DO3SE

  use CheckStop_ml, only: CheckStop_TF, CheckStop
  use DO3SE_GstoTypes_ml, only: GstoConfig_t, GstoParams_t, check_GstoConfig
  use DO3SE_Gsto_ml
  use DO3SE_Photosynthesis_ml, only: gsto_pn, &
                                     PnGstoConfig_t, &
                                     check_PnGstoConfig
  use DO3SE_SMD_ml, only: SMDConfig_t, &
                          SMDData_t, &
                          check_SMDConfig, &
                          soil_moisture_from_SWC
  use DO3SE_Util_ml, only: assert, &
                           assert_not, &
                           inverted_assert_not
  use DO3SE_PhysicalConstants_ml, only: PAR_Wm2_to_photons
  use esx_Variables, only: esx, Zmet, Zveg, Loc, NZMAX => ESX_MAXNZ
  use esx_Zveg, only: Veg

  implicit none
  private

  public :: config_DO3SE
  public :: DO3SE_gleaf

  type(GstoConfig_t), public, save :: conf
  type(PnGstoConfig_t), public, save :: pn_conf
  type(SMDConfig_t), public, save :: SMD_conf
  type(GstoParams_t), dimension(NZMAX), public, save :: gsto_params

contains

  subroutine config_DO3SE()
    integer :: config_io, io

    namelist /do3se_config/ conf, pn_conf, SMD_conf

    assert => inverted_assert_not
    assert_not => CheckStop_TF

    conf = GstoConfig_t()
    pn_conf = PnGstoConfig_t()
    SMD_conf = SMDConfig_t()
    open (newunit=config_io, file="config_esx.nml")
    rewind(config_io)
    read (config_io, nml=do3se_config)
    close(config_io)
    open(newunit=io,file=trim(esx%odir)//"/LogConfig.do3se")
    write(io,nml=do3se_config)
    close(io)
    if ( esx%debug_gleaf>0) write(*,nml=do3se_config)

    call check_GstoConfig(conf)
    call check_PnGstoConfig(pn_conf)
    call check_SMDConfig(SMD_conf)
  end subroutine config_DO3SE

  subroutine DO3SE_gleaf()
    integer :: iz
    real :: mmol2sm   ! mmol O3/m2/s to s/m
    type(SMDData_t) :: SMD

    ! Reset gsto parameters
    gsto_params = GstoParams_t(fmin=conf%fmin, gmax=conf%gmax, gmorph=conf%gmorph)

    ! Influences which don't vary with height
    select case (conf%f_phen_method)
    case ("disabled")
      ! Nothing to do
    case ("simple day PLF")
      gsto_params%f_phen = f_phen_simple_PLF(conf, int(Veg%SGS), int(Veg%EGS), esx%daynumber)
    case ("complex day PLF")
      gsto_params%f_phen = f_phen_complex_PLF(conf, int(Veg%SGS), int(Veg%EGS), esx%daynumber)
    case default
      call CheckStop("DO3SE_gleaf: unknown f_phen_method: "//trim(conf%f_phen_method))
    end select

    select case (conf%leaf_f_phen_method)
    case ("disabled")
      ! Nothing to do
    case ("f_phen")
      gsto_params%leaf_f_phen = gsto_params%f_phen
    case ("simple day PLF")
      gsto_params%leaf_f_phen = leaf_f_phen_PLF(conf, int(Veg%Astart), int(Veg%Aend), esx%daynumber)
    case default
      call CheckStop("DO3SE_gleaf: unknown leaf_f_phen_method: "//trim(conf%leaf_f_phen_method))
    end select

    ! Calculate f_SW
    SMD = SMDData_t()
    ! TODO: real soil moisture data
    call soil_moisture_from_SWC(SMD_conf, SMD, SMD_conf%soil%FC)
    select case (conf%f_SW_method)
    case ("disabled")
      ! Nothing to do
    case ("fSWP exp")
      gsto_params%f_SW = f_SWP_exp(conf%fSWP_exp_a, conf%fSWP_exp_b, conf%fmin, SMD%SWP)
    case ("fSWP linear")
      gsto_params%f_SW = f_SWP_linear(conf%SWP_min, conf%SWP_max, conf%fmin, SMD%SWP)
    ! TODO: implement LWP
    !case ("fLWP exp")
    !  gp%f_SW = f_SWP_exp(gc%fSWP_exp_a, gc%fSWP_exp_b, gc%fmin, SMD%LWP)
    case ("fPAW")
      gsto_params%f_SW = f_PAW(SMD_conf%ASW_FC, conf%fmin, SMD%ASW)
    case default
      call CheckStop("DO3SE_gleaf: unknown f_SW_method: "//trim(conf%f_SW_method))
    end select

    !gsto_p%f_O3 = ???

    ! Influences which do vary with height
    do iz = esx%nhVeg, 1, -1
      associate (gp => gsto_params(iz), &
                 met => Zmet(iz))

        select case (conf%f_light_method)
        case ("disabled")
          ! Nothing to do
        case ("enabled")
          ! TODO: needs radiation profiles for PARsun and PARshade, plus sunlit LAI fractions
          !gp%f_light = f_light(conf%f_lightfac, Zveg(iz)%PARsun, Zveg(iz)%PARshade, Zveg(iz)%LAIsunfrac)
          gp%leaf_f_light = leaf_f_light(conf%f_lightfac, Zveg(iz)%PARz)
        case default
          call CheckStop( "DO3SE_gleaf: unknown f_light_method: "//trim(conf%f_light_method))
        end select

        select case (conf%f_temp_method)
        case ("disabled")
          ! Nothing to do
        case ("default")
          gp%f_temp = f_temp(met%tzC, conf%T_min, conf%T_opt, conf%T_max, conf%fmin)
        case ("square high")
          gp%f_temp = f_temp_square_high(met%tzC, conf%T_min, conf%T_opt, conf%T_max, conf%fmin)
        case default
          call CheckStop("DO3SE_gleaf: unknown f_temp_method: "//trim(conf%f_temp_method))
        end select

        select case (conf%f_VPD_method)
        case ("disabled")
          ! Nothing to do
        case ("linear")
          gp%f_VPD = f_VPD_linear(met%VPD, conf%VPD_max, conf%VPD_min, conf%fmin)
        case ("log")
          gp%f_VPD = f_VPD_log(met%VPD, conf%fmin)
        case default
          call CheckStop("DO3SE_gleaf: unknown f_VPD_method: "//trim(conf%f_VPD_method))
        end select

        ! Calculate mean stomatal conductance
        select case (conf%method)
        case ("multiplicative")
          Zveg(iz)%gsto = gsto_mean(gp)
        case ("photosynthesis")
          gp%gmax = gsto_pn(pn_conf, met%tzC, met%tleafC, met%uz, &
                            met%CO2, met%rh, Zveg(iz)%PARz * PAR_Wm2_to_photons, Veg%Lm)
          Zveg(iz)%gsto = gsto_mean(gp)
        case default
          call CheckStop("DO3SE_gleaf: unknown gsto method: "//trim(conf%method))
        end select

        ! Convert to m/s units
        ! From mmol O3/m2/s to s/m given in Jones, App. 3, gives 41000 for 20 deg.C )
        ! (Assumes P=100 hPa?)

        mmol2sm = 8.3144e-8 * met%tzK    ! 0.001 * RT/P. Or should it be Tleaf?

        Zveg(iz)%gsto = Zveg(iz)%gsto * mmol2sm

        if ( esx%debug_gleaf > 0 ) then

           if( iz == esx%nhVeg )  then
              print "(3(1x,a))", "DEBUG_gleaf F-phens methods:",&
                  trim(conf%f_phen_method), trim(conf%leaf_f_phen_method)

              print "(a,a7,2a12,88a8)", "DEBUG_gleaf:","z", "PARz","gsto(cm/s)",&
               "fmin","gmax", "gmorph","fphen","lfphen","flight","lflight",&
               "ftemp","fVPD","fSW","fO3"
           end if

           print "(a,i7,2f12.3,88f8.2)", "DEBUG_gleaf:", iz, Zveg(iz)%PARz, &
              100.0*Zveg(iz)%gsto, gsto_params(iz) !, gsto_params%leaf_f_phen

        end if ! debug

      end associate
    end do
  end subroutine DO3SE_gleaf

end module esx_DO3SE
