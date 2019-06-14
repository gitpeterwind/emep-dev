! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem16zi VBS_EmChem16z APINENE_BVOC_emis BVOC_EmChem16z Aqueous_EmChem16x Aero2017nx Ash_EmChem16z ShipNOx16z FFireInert16z SeaSalt16z DustExtended16z
module ChemRates_mod

  use AeroConstants_mod  ! => AERO%PM etc, ...
  use AeroFunctions_mod  ! => UptakeRates, ...
  use ChemDims_mod       ! => NSPEC_TOT, NCHEMRATES, ....
  use ChemFunctions_mod  ! => IUPAC_troe, RiemerN2O5, ....
  use ChemSpecs_mod      ! => PINALD, .... for FgasJ08
  use DefPhotolysis_mod  ! => IDNO2 etc.
  use ZchemData_mod      ! => rct
  
  implicit none
  private
  character(len=*),parameter, public :: CM_schemes_ChemRates = " EmChem16zi VBS_EmChem16z APINENE_BVOC_emis BVOC_EmChem16z Aqueous_EmChem16x Aero2017nx Ash_EmChem16z ShipNOx16z FFireInert16z SeaSalt16z DustExtended16z"
  
  
  public :: setChemRates
  public :: setPhotolUsed
  
  ! Photolysis rates
  integer, save, public, dimension(NPHOTOLRATES) :: photol_used

contains
  
  subroutine setPhotolUsed()
    photol_used = (/ &
        IDMEK  &
      , IDCH3COY  &
      , IDNO2  &
      , IDRCOCHO  &
      , IDCH3O2H  &
      , IDHCHO_H2  &
      , IDCH3CHO  &
      , IDO3_O1D  &
      , IDO3_O3P  &
      , IDH2O2  &
      , IDNO3  &
      , IDHONO  &
      , IDHNO3  &
      , IDHO2NO2  &
      , IDHCHO_H  &
      , IDCHOCHO  &
    /)  
  end subroutine setPhotolUsed
  
  subroutine setChemRates()
    !integer, intent(in) :: debug_level
  
    rct(1,:) = ((5.681e-34*EXP(-2.6*(LOG(TEMP/300))))*O2)*M
    rct(2,:) = (2.15e-11*EXP(110.0*TINV))*N2
    rct(3,:) = (3.2e-11*EXP(67.0*TINV))*O2
    rct(4,:) = (2.14e-10)*H2O
    rct(5,:) = 1.4e-12*EXP(-1310.0*TINV)
    rct(6,:) = 1.4e-13*EXP(-2470.0*TINV)
    rct(7,:) = 1.7e-12*EXP(-940.0*TINV)
    rct(8,:) = 2.03e-16*EXP(-4.57*(LOG(300/TEMP)))*EXP(693.0*TINV)
    rct(9,:) = 1.8e-11*EXP(110.0*TINV)
    rct(10,:) = 3.45e-12*EXP(270.0*TINV)
    rct(11,:) = 4.5e-14*EXP(-1260.0*TINV)
    rct(12,:) = 4.8e-11*EXP(250.0*TINV)
    rct(13,:) = 2.9e-12*EXP(-160.0*TINV)
    rct(14,:) = 7.7e-12*EXP(-2100.0*TINV)
    rct(15,:) = (KMT3(2.4e-14,  &
              & 460.0,  &
              & 6.5e-34,  &
              & 1335.0,  &
              & 2.7e-17,  &
              & 2199.0))
    rct(16,:) = ((1.0+1.4e-21*H2O*EXP(2200.0*TINV)))*2.2e-13*EXP(600.0*TINV)
    rct(17,:) = (((1.0+1.4e-21*H2O*EXP(2200.0*TINV)))  &
              & *1.9e-33*EXP(980.0*TINV))*M
    rct(18,:) = 2.5e-12*EXP(260.0*TINV)
    rct(19,:) = (IUPAC_TROE(1.4e-31*EXP(3.1*(LOG(300/TEMP))),  &
              & 4.0e-12,  &
              & 0.4,  &
              & M,  &
              & 0.75-1.27*LOG10(0.4)))
    rct(20,:) = (IUPAC_TROE(4.10e-05*EXP(-10650.0*TINV),  &
              & 6.0e+15*EXP(-11170.0*TINV),  &
              & 0.4,  &
              & M,  &
              & 0.75-1.27*LOG10(0.4)))
    rct(21,:) = 3.2e-13*EXP(690.0*TINV)
    rct(22,:) = 1.85e-12*EXP(-1690.0*TINV)
    rct(23,:) = (1.44e-13*(1+(M/4.2e+19)))
    rct(24,:) = 2.3e-12*EXP(360.0*TINV)
    rct(25,:) = 7.4e-13*EXP(-520.0*TINV)*2
    rct(26,:) = 2*(1.03e-13*EXP(365.0*TINV)-7.4e-13*EXP(-520.0*TINV))
    rct(27,:) = 3.8e-13*EXP(780.0*TINV)
    rct(28,:) = 2.85e-12*EXP(-345.0*TINV)
    rct(29,:) = 5.3e-12*EXP(190.0*TINV)
    rct(30,:) = 5.4e-12*EXP(135.0*TINV)
    rct(31,:) = 2.0e-12*EXP(-2440.0*TINV)
    rct(32,:) = 6.9e-12*EXP(-1000.0*TINV)
    rct(33,:) = 2.55e-12*EXP(380.0*TINV)
    rct(34,:) = 6.4e-13*EXP(710.0*TINV)
    rct(35,:) = 4.7e-12*EXP(345.0*TINV)
    rct(36,:) = (1.4e-12*EXP(-1860.0*TINV))
    rct(37,:) = (7.5e-12*EXP(290.0*TINV))
    rct(38,:) = (3.14e-12*EXP(580.0*TINV))
    rct(39,:) = 9.8e-12*EXP(-425.0*TINV)
    rct(40,:) = (2.91e-13*EXP(1300.0*TINV))*0.625
    rct(41,:) = (2.7e-12*EXP(360.0*TINV))
    rct(42,:) = 1.53e-13*EXP(1300.0*TINV)
    rct(43,:) = 6.82e-15*EXP(-2500.0*TINV)
    rct(44,:) = 5.77e-15*EXP(-1880.0*TINV)
    rct(45,:) = (2.91e-13*EXP(1300.0*TINV))*0.52
    rct(46,:) = 4.60e-13*EXP(-1155*TINV)
    rct(47,:) = 2.3e-12*EXP(-190.0*TINV)
    rct(48,:) = 1.8e-12*EXP(340.0*TINV)
    rct(49,:) = (2.91e-13*EXP(1300.0*TINV))*0.859
    rct(50,:) = (2.3e-12)
    rct(51,:) = (2.91e-13*EXP(1300.0*TINV))*0.706
    rct(52,:) = (1.4e-12*EXP(-1860.0*TINV))*5.5
    rct(53,:) = 3.1e-12*EXP(340.0*TINV)
    rct(54,:) = 1.9e-12*EXP(575.0*TINV)
    rct(55,:) = 2.7e-11*EXP(390.0*TINV)
    rct(56,:) = 1.05e-14*EXP(-2000.0*TINV)
    rct(57,:) = 2.95e-12*EXP(-450.0*TINV)
    rct(58,:) = 1.30e-12*EXP(610.0*TINV)
    rct(59,:) = 4.00e-12*EXP(380.0*TINV)
    rct(60,:) = 4.26e-16*EXP(-1520.0*TINV)
    rct(61,:) = 7.0e-16*EXP(-2100.0*TINV)
    rct(62,:) = 1.60e-12*EXP(305.0*TINV)
    rct(63,:) = (IUPAC_TROE(3.28e-28*EXP(6.87*(LOG(300/TEMP))),  &
              & 1.125e-11*EXP(1.105*(LOG(300/TEMP))),  &
              & 0.3,  &
              & M,  &
              & 0.75-1.27*LOG10(0.3)))*0.107
    rct(64,:) = (IUPAC_TROE(1.10e-05*EXP(-10100.0*TINV),  &
              & 1.90e17*EXP(-14100.0*TINV),  &
              & 0.3,  &
              & M,  &
              & 0.75-1.27*LOG10(0.3)))
    rct(65,:) = 1.88e09*EXP(-7261*TINV)*0.052
    rct(66,:) = 1.45e12*EXP(-10688*TINV)
    rct(67,:) = (IUPAC_TROE(1.0e-31*EXP(1.6*(LOG(300/TEMP))),  &
              & 5.0e-11*EXP(0.3*(LOG(300/TEMP))),  &
              & 0.85,  &
              & M,  &
              & 0.75-1.27*LOG10(0.85)))
    rct(68,:) = (IUPAC_TROE(3.6e-30*EXP(4.1*(LOG(300/TEMP))),  &
              & 1.9e-12*EXP(-0.2*(LOG(300/TEMP))),  &
              & 0.35,  &
              & M,  &
              & 0.75-1.27*LOG10(0.35)))
    rct(69,:) = (IUPAC_TROE(1.3e-3*EXP(3.5*(LOG(300/TEMP)))*EXP(-11000.0*TINV),  &
              & 9.70e14*EXP(-0.1*(LOG(300/TEMP)))*EXP(-11080.0*TINV),  &
              & 0.35,  &
              & M,  &
              & 0.75-1.27*LOG10(0.35)))
    rct(70,:) = (IUPAC_TROE(3.2e-30*EXP(4.5*(LOG(300/TEMP))),  &
              & 3.0e-11,  &
              & 0.41,  &
              & M,  &
              & 0.75-1.27*LOG10(0.41)))
    rct(71,:) = (IUPAC_TROE(3.28e-28*EXP(6.87*(LOG(300/TEMP))),  &
              & 1.125e-11*EXP(1.105*(LOG(300/TEMP))),  &
              & 0.3,  &
              & M,  &
              & 0.75-1.27*LOG10(0.3)))
    rct(72,:) = (IUPAC_TROE(8.6e-29*EXP(3.1*(LOG(300/TEMP))),  &
              & 9.0e-12*EXP(0.85*(LOG(300/TEMP))),  &
              & 0.48,  &
              & M,  &
              & 0.75-1.27*LOG10(0.48)))
    rct(73,:) = (IUPAC_TROE(8.0e-27*EXP(3.5*(LOG(300/TEMP))),  &
              & 9.0e-9*TINV,  &
              & 0.5,  &
              & M,  &
              & 0.75-1.27*LOG10(0.5)))
    rct(74,:) = (IUPAC_TROE(7.4e-31*EXP(2.4*(LOG(300/TEMP))),  &
              & 3.3e-11*EXP(0.3*(LOG(300/TEMP))),  &
              & 0.81,  &
              & M,  &
              & 0.75-1.27*LOG10(0.81)))
    rct(75,:) = 1.20e-11*EXP(440.0*TINV)
    rct(76,:) = 8.05e-16*EXP(-640.0*TINV)
    rct(77,:) = 1.2e-12*EXP(490.0*TINV)
    rct(78,:) = 4.0e-12*FGAS(ASOC_UG1,  &
              & :)
    rct(79,:) = 4.0e-12*FGAS(ASOC_UG10,  &
              & :)
    rct(80,:) = 4.0e-12*FGAS(ASOC_UG1e2,  &
              & :)
    rct(81,:) = 4.0e-12*FGAS(ASOC_UG1e3,  &
              & :)
    rct(82,:) = 4.0e-12*FGAS(NON_C_ASOA_UG1,  &
              & :)
    rct(83,:) = 4.0e-12*FGAS(NON_C_ASOA_UG10,  &
              & :)
    rct(84,:) = 4.0e-12*FGAS(NON_C_ASOA_UG1e2,  &
              & :)
    rct(85,:) = 4.0e-12*FGAS(NON_C_ASOA_UG1e3,  &
              & :)
    rct(86,:) = 4.0e-12*FGAS(BSOC_UG1,  &
              & :)
    rct(87,:) = 4.0e-12*FGAS(BSOC_UG10,  &
              & :)
    rct(88,:) = 4.0e-12*FGAS(BSOC_UG1e2,  &
              & :)
    rct(89,:) = 4.0e-12*FGAS(BSOC_UG1e3,  &
              & :)
    rct(90,:) = 4.0e-12*FGAS(NON_C_BSOA_UG1,  &
              & :)
    rct(91,:) = 4.0e-12*FGAS(NON_C_BSOA_UG10,  &
              & :)
    rct(92,:) = 4.0e-12*FGAS(NON_C_BSOA_UG1e2,  &
              & :)
    rct(93,:) = 4.0e-12*FGAS(NON_C_BSOA_UG1e3,  &
              & :)
    rct(94,:) = EC_ageING_RATE()
    rct(95,:) = 1.34e-11*EXP(410.0*TINV)
    rct(96,:) = 7.5e-13*EXP(700.0*TINV)
    rct(97,:) = 3.8e-12*EXP(200.0*TINV)
    rct(98,:) = 2.6e-12*EXP(610.0*TINV)
    rct(99,:) = 8.5e-16*EXP(-1520.0*TINV)
    rct(100,:) = HYDROLYSISN2O5()
    rct(101,:) = UPTAKERATE(CNO3,  &
               & GAM=0.001,  &
               & S=S_M2M3(AERO%PM,  &
               & :))
    rct(102,:) = UPTAKERATE(CHNO3,  &
               & GAM=0.1,  &
               & S=S_M2M3(AERO%DU_C,  &
               & :))
    rct(103,:) = UPTAKERATE(CHNO3,  &
               & GAM=0.01,  &
               & S=S_M2M3(AERO%SS_C,  &
               & :))
    rct(104,:) = UPTAKERATE(CHO2,  &
               & GAM=0.2,  &
               & S=S_M2M3(AERO%PM,  &
               & :))
    rct(105,:) = 1.44e-13+M*3.43e-33
  
  end subroutine setChemRates

end module ChemRates_mod
