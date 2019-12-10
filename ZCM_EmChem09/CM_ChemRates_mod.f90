! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem09 VBS_acp2012 APINENE_BVOC_emis BVOC_EmChem19 Aqueous_EmChem16x Aero2017nx Ash ShipNOx FFireInert SeaSalt DustExtended Pollen
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
  character(len=*),parameter, public :: CM_schemes_ChemRates = " EmChem09 VBS_acp2012 APINENE_BVOC_emis BVOC_EmChem19 Aqueous_EmChem16x Aero2017nx Ash ShipNOx FFireInert SeaSalt DustExtended Pollen"
  
  
  public :: setChemRates
  public :: setPhotolUsed
  
  ! Photolysis rates
  integer, save, public, dimension(NPHOTOLRATES) :: photol_used

contains
  
  subroutine setPhotolUsed()
    photol_used = (/ &
        IDO3_O3P  &
      , IDO3_O1D  &
      , IDNO2  &
      , IDH2O2  &
      , IDHNO3  &
      , IDHCHO_H  &
      , IDHCHO_H2  &
      , IDCH3CHO  &
      , IDNO3_NO2  &
      , IDCH3O2H  &
      , IDCHOCHO  &
      , IDRCOCHO  &
      , IDMEK  &
    /)  
  end subroutine setPhotolUsed
  
  subroutine setChemRates()
    !integer, intent(in) :: debug_level
  
    rct(1,:) = (6.0e-34*O2+5.6e-34*N2)*O2*EXP(-2.6*(LOG(TEMP/300)))
    rct(2,:) = 1.8e-11*N2*EXP(107.0*TINV)
    rct(3,:) = 3.2e-11*O2*EXP(67.0*TINV)
    rct(4,:) = 2.2e-10*H2O
    rct(5,:) = 1.4e-12*EXP(-1310.0*TINV)
    rct(6,:) = 1.4e-13*EXP(-2470.0*TINV)
    rct(7,:) = 1.7e-12*EXP(-940.0*TINV)
    rct(8,:) = 2.03e-16*EXP(-4.57*(LOG(300/TEMP)))*EXP(693.0*TINV)
    rct(9,:) = 1.8e-11*EXP(110.0*TINV)
    rct(10,:) = 3.6e-12*EXP(270.0*TINV)
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
    rct(17,:) = ((1.0+1.4e-21*H2O*EXP(2200.0*TINV)))*1.9e-33*EXP(980.0*TINV)  &
              & *M
    rct(18,:) = 2.5e-12*EXP(-260.0*TINV)
    rct(19,:) = 1.85e-20*EXP(2.82*LOG(TEMP))*EXP(-987.0*TINV)
    rct(20,:) = 1.44e-13+M*3.43e-33
    rct(21,:) = 2.3e-12*EXP(360.0*TINV)
    rct(22,:) = 7.4e-13*EXP(-520.0*TINV)
    rct(23,:) = 1.03e-13*EXP(365.0*TINV)-7.4e-13*EXP(-520.0*TINV)
    rct(24,:) = 6.38e-18*(TEMP**2)*EXP(144.0*TINV)
    rct(25,:) = 3.8e-13*EXP(780.0*TINV)
    rct(26,:) = 5.3e-12*EXP(190.0*TINV)
    rct(27,:) = 1.25e-17*(TEMP**2)*EXP(615.0*TINV)
    rct(28,:) = 2.0e-12*EXP(-2440.0*TINV)
    rct(29,:) = 6.9e-12*EXP(-1000.0*TINV)
    rct(30,:) = 2.55e-12*EXP(380.0*TINV)
    rct(31,:) = 3.8e-13*EXP(900.0*TINV)
    rct(32,:) = (1.9e-12*EXP(190.0*TINV))
    rct(33,:) = 4.4e-12*EXP(365.0*TINV)
    rct(34,:) = 7.5e-12*EXP(290.0*TINV)
    rct(35,:) = 2.0e-12*EXP(500.0*TINV)
    rct(36,:) = 2.9e-12*EXP(500.0*TINV)
    rct(37,:) = 5.2e-13*EXP(980.0*TINV)
    rct(38,:) = 1.9e-12*EXP(190.0*TINV)
    rct(39,:) = 6.7e-18*(TEMP**2)*EXP(511.0*TINV)
    rct(40,:) = 2.03e-17*(TEMP**2)*EXP(78.0*TINV)
    rct(41,:) = (2.54e-12*EXP(360.0*TINV))
    rct(42,:) = 0.625*(2.91e-13*EXP(1300.0*TINV))
    rct(43,:) = 2.53e-18*(TEMP**2)*EXP(503.0*TINV)
    rct(44,:) = 9.1e-15*EXP(-2580.0*TINV)
    rct(45,:) = 5.5e-15*EXP(-1880.0*TINV)
    rct(46,:) = 0.52*(2.91e-13*EXP(1300.0*TINV))
    rct(47,:) = 0.859*(2.91e-13*EXP(1300.0*TINV))
    rct(48,:) = 0.706*(2.91e-13*EXP(1300.0*TINV))
    rct(49,:) = 6.6e-18*(TEMP**2)*EXP(820.0*TINV)
    rct(50,:) = 1.9e-12*EXP(575.0*TINV)
    rct(51,:) = 1.03e-14*EXP(-1995.0*TINV)
    rct(52,:) = 2.7e-11*EXP(390.0*TINV)
    rct(53,:) = 2.6e-12*EXP(610.0*TINV)
    rct(54,:) = 1.36e-15*EXP(-2112.0*TINV)
    rct(55,:) = 8.0e-12*EXP(380.0*TINV)
    rct(56,:) = 7.6e-12*EXP(180.0*TINV)
    rct(57,:) = (2.5e-12)
    rct(58,:) = 1.6e-12*EXP(305.0*TINV)
    rct(59,:) = 8.7e-12*EXP(290.0*TINV)
    rct(60,:) = 8.5e-16*EXP(-1520.0*TINV)
    rct(61,:) = 3.15e-12*EXP(-450.0*TINV)
    rct(62,:) = 4.3e-13*EXP(1040.0*TINV)
    rct(63,:) = (IUPAC_TROE(1.0e-31*EXP(1.6*(LOG(300/TEMP))),  &
              & 3.0e-11*EXP(-0.3*(LOG(300/TEMP))),  &
              & 0.85,  &
              & M,  &
              & 0.75-1.27*LOG10(0.85)))
    rct(64,:) = (IUPAC_TROE(3.6e-30*EXP(4.1*(LOG(300/TEMP))),  &
              & 1.9e-12*EXP(-0.2*(LOG(300/TEMP))),  &
              & 0.35,  &
              & M,  &
              & 0.75-1.27*LOG10(0.35)))
    rct(65,:) = (IUPAC_TROE(1.3e-3*EXP(3.5*(LOG(300/TEMP)))*EXP(-11000.0*TINV),  &
              & 9.70e14*EXP(-0.1*(LOG(300/TEMP)))*EXP(-11080.0*TINV),  &
              & 0.35,  &
              & M,  &
              & 0.75-1.27*LOG10(0.35)))
    rct(66,:) = (IUPAC_TROE(3.3e-30*EXP(3.0*(LOG(300/TEMP))),  &
              & 4.1e-11,  &
              & 0.40,  &
              & M,  &
              & 0.75-1.27*LOG10(0.4)))
    rct(67,:) = (IUPAC_TROE(2.7e-28*EXP(7.1*(LOG(300/TEMP))),  &
              & 1.2e-11*EXP(0.9*(LOG(300/TEMP))),  &
              & 0.3,  &
              & M,  &
              & 0.75-1.27*LOG10(0.3)))
    rct(68,:) = (IUPAC_TROE(4.9e-3*EXP(-12100.0*TINV),  &
              & 5.4e16*EXP(-13830.0*TINV),  &
              & 0.3,  &
              & M,  &
              & 0.75-1.27*LOG10(0.3)))
    rct(69,:) = (IUPAC_TROE(8.6e-29*EXP(3.1*(LOG(300/TEMP))),  &
              & 9.0e-12*EXP(0.85*(LOG(300/TEMP))),  &
              & 0.48,  &
              & M,  &
              & 0.75-1.27*LOG10(0.48)))
    rct(70,:) = (IUPAC_TROE(8.0e-27*EXP(3.5*(LOG(300/TEMP))),  &
              & 3.0e-11*300.0*TINV,  &
              & 0.5,  &
              & M,  &
              & 0.75-1.27*LOG10(0.5)))
    rct(71,:) = (IUPAC_TROE(7.4e-31*EXP(2.4*(LOG(300/TEMP))),  &
              & 3.3e-11*EXP(0.3*(LOG(300/TEMP))),  &
              & EXP(-TEMP/1420.0),  &
              & M,  &
              & 0.75+3.884e-4*TEMP))
    rct(72,:) = 6.3e-16*EXP(-580.0*TINV)
    rct(73,:) = 1.2e-11*EXP(444.0*TINV)
    rct(74,:) = 1.2e-12*EXP(490.0*TINV)
    rct(75,:) = 0.914*(2.91e-13*EXP(1300.0*TINV))
    rct(76,:) = 4.0e-12*FGAS(ASOC_UG1,  &
              & :)
    rct(77,:) = 4.0e-12*FGAS(ASOC_UG10,  &
              & :)
    rct(78,:) = 4.0e-12*FGAS(ASOC_UG1e2,  &
              & :)
    rct(79,:) = 4.0e-12*FGAS(ASOC_UG1e3,  &
              & :)
    rct(80,:) = 4.0e-12*FGAS(NON_C_ASOA_UG1,  &
              & :)
    rct(81,:) = 4.0e-12*FGAS(NON_C_ASOA_UG10,  &
              & :)
    rct(82,:) = 4.0e-12*FGAS(NON_C_ASOA_UG1e2,  &
              & :)
    rct(83,:) = 4.0e-12*FGAS(NON_C_ASOA_UG1e3,  &
              & :)
    rct(84,:) = 4.0e-12*FGAS(BSOC_UG1,  &
              & :)
    rct(85,:) = 4.0e-12*FGAS(BSOC_UG10,  &
              & :)
    rct(86,:) = 4.0e-12*FGAS(BSOC_UG1e2,  &
              & :)
    rct(87,:) = 4.0e-12*FGAS(BSOC_UG1e3,  &
              & :)
    rct(88,:) = 4.0e-12*FGAS(NON_C_BSOA_UG1,  &
              & :)
    rct(89,:) = 4.0e-12*FGAS(NON_C_BSOA_UG10,  &
              & :)
    rct(90,:) = 4.0e-12*FGAS(NON_C_BSOA_UG1e2,  &
              & :)
    rct(91,:) = 4.0e-12*FGAS(NON_C_BSOA_UG1e3,  &
              & :)
    rct(92,:) = EC_AGEING_RATE()
    rct(93,:) = 1.34e-11*EXP(410.0*TINV)
    rct(94,:) = 7.5e-13*EXP(700.0*TINV)
    rct(95,:) = 3.8e-12*EXP(200.0*TINV)
    rct(96,:) = 8.22e-16*EXP(-640.0*TINV)
    rct(97,:) = HYDROLYSISN2O5()
    rct(98,:) = UPTAKERATE(CNO3,  &
              & GAM=0.001,  &
              & S=S_M2M3(AERO%PM,  &
              & :))
    rct(99,:) = UPTAKERATE(CHNO3,  &
              & GAM=0.1,  &
              & S=S_M2M3(AERO%DU_C,  &
              & :))
    rct(100,:) = UPTAKERATE(CHNO3,  &
               & GAM=0.01,  &
               & S=S_M2M3(AERO%SS_C,  &
               & :))
    rct(101,:) = UPTAKERATE(CHO2,  &
               & GAM=0.2,  &
               & S=S_M2M3(AERO%PM,  &
               & :))
    rct(102,:) = (((3.2e-30*M*(TEMP/300)**(-4.5))*(3.0e-11))  &
               & *(10**(LOG10(0.41)/(1+(LOG10(((3.2e-30*M*(TEMP/300)**(-4.5))  &
               & /(3.0e-11)))/(0.75-1.27*(LOG10(0.41))))**2)))  &
               & /((3.2e-30*M*(TEMP/300)**(-4.5))+(3.0e-11)))
  
  end subroutine setChemRates

end module ChemRates_mod
