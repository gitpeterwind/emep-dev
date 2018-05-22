! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem16x ShipNOx BVOC_EmChem16x Aero2017nx VBS_acp2012 FFireInert SeaSalt Dust
module ChemGroups_mod

  use ChemSpecs_mod        ! => species indices
  use OwnDataTypes_mod  ! => typ_sp, typ_factors, typ_maps
  
  implicit none
  private
  character(len=*),parameter, public :: CM_schemes_ChemGroups = " EmChem16x ShipNOx BVOC_EmChem16x Aero2017nx VBS_acp2012 FFireInert SeaSalt Dust"
  
  
  ! Assignment of groups from GenIn_Species.csv
  public :: Init_ChemGroups
  
  type(typ_sp), dimension(124), public, save :: chemgroups
  type(typ_factors), dimension(2), public, save :: chemgroups_factors
  type(typ_maps), dimension(1), public, save :: chemgroups_maps
  
  integer, public, target, save, dimension (14) :: &
    RO2_GROUP = (/  &
      CH3O2, C2H5O2, SECC4H9O2, ISRO2, ETRO2, PRRO2, OXYO2, MEKO2,  &
      MALO2, MACRO2, CH3CO3, TERPO2, XMTO3_RO2, TERPPeroxy  &
    /)
  
  integer, public, target, save, dimension (2) :: &
    OX_GROUP = (/ O3, NO2 /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_OX_GROUP = (/ O3, NO2 /)
  
  integer, public, target, save, dimension (3) :: &
    NOX_GROUP = (/ NO, NO2, SHIPNOX /)
  
  integer, public, target, save, dimension (14) :: &
    OXN_GROUP = (/  &
      NO, NO2, PAN, MPAN, NO3, N2O5, HNO3, HONO, HO2NO2, ISON,  &
      NALD, NO3_f, NO3_c, SHIPNOX  &
    /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_NOX_GROUP = (/ NO2, SHIPNOX /)
  
  integer, public, target, save, dimension (9) :: &
    DDEP_OXN_GROUP = (/ NO2, PAN, MPAN, HNO3, HONO, HO2NO2, NO3_f, NO3_c, SHIPNOX /)
  
  integer, public, target, save, dimension (1) :: &
    DAOBS_GROUP = (/ NO2 /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_daObs_GROUP = (/ NO2 /)
  
  integer, public, target, save, dimension (2) :: &
    PAN_GROUP = (/ PAN, MPAN /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_PAN_GROUP = (/ PAN, MPAN /)
  
  integer, public, target, save, dimension (5) :: &
    WDEP_OXN_GROUP = (/ HNO3, HONO, HO2NO2, NO3_f, NO3_c /)
  
  integer, public, target, save, dimension (7) :: &
    BVOC_GROUP = (/ C5H8, APINENE, MVK, BIOTERP, BPINENE, XTERP, SQT_SOA_NV /)
  
  integer, public, target, save, dimension (14) :: &
    ROOH_GROUP = (/  &
      CH3OOH, C2H5OOH, BURO2H, ETRO2H, PRRO2H, OXYO2H, MEKO2H,  &
      MALO2H, MACROOH, ISRO2H, H2O2, CH3CO3H, HPALD, TERPOOH  &
    /)
  
  integer, public, target, save, dimension (4) :: &
    DDEP_ROOH_GROUP = (/ CH3OOH, C2H5OOH, H2O2, TERPOOH /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_ROOH_GROUP = (/ H2O2 /)
  
  integer, public, target, save, dimension (2) :: &
    SOX_GROUP = (/ SO2, SO4 /)
  
  integer, public, target, save, dimension (2) :: &
    WDEP_SOX_GROUP = (/ SO2, SO4 /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_SOX_GROUP = (/ SO2, SO4 /)
  
  integer, public, target, save, dimension (21) :: &
    PM10_GROUP = (/  &
      SO4, NO3_f, NO3_c, NH4_f, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, REMPPM25, REMPPM_c, OM25_p, VBS_TEST, FFIRE_BC,  &
      FFIRE_REMPPM25, SeaSalt_f, SeaSalt_c, Dust_f, Dust_c  &
    /)
  
  integer, public, target, save, dimension (19) :: &
    WDEP_PM10_GROUP = (/  &
      SO4, NO3_f, NO3_c, NH4_f, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, REMPPM25, REMPPM_c, FFIRE_BC, FFIRE_REMPPM25,  &
      SeaSalt_f, SeaSalt_c, Dust_f, Dust_c  &
    /)
  
  integer, public, target, save, dimension (19) :: &
    DDEP_PM10_GROUP = (/  &
      SO4, NO3_f, NO3_c, NH4_f, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, REMPPM25, REMPPM_c, FFIRE_BC, FFIRE_REMPPM25,  &
      SeaSalt_f, SeaSalt_c, Dust_f, Dust_c  &
    /)
  
  integer, public, target, save, dimension (14) :: &
    PMFINE_GROUP = (/  &
      SO4, NO3_f, NH4_f, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, REMPPM25, OM25_p, VBS_TEST,  &
      FFIRE_BC, FFIRE_REMPPM25, SeaSalt_f, Dust_f  &
    /)
  
  integer, public, target, save, dimension (12) :: &
    WDEP_PMFINE_GROUP = (/  &
      SO4, NO3_f, NH4_f, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, REMPPM25, FFIRE_BC,  &
      FFIRE_REMPPM25, SeaSalt_f, Dust_f  &
    /)
  
  integer, public, target, save, dimension (12) :: &
    DDEP_PMFINE_GROUP = (/  &
      SO4, NO3_f, NH4_f, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, REMPPM25, FFIRE_BC,  &
      FFIRE_REMPPM25, SeaSalt_f, Dust_f  &
    /)
  
  integer, public, target, save, dimension (4) :: &
    SIA_GROUP = (/ SO4, NO3_f, NO3_c, NH4_f /)
  
  integer, public, target, save, dimension (4) :: &
    WDEP_SIA_GROUP = (/ SO4, NO3_f, NO3_c, NH4_f /)
  
  integer, public, target, save, dimension (4) :: &
    DDEP_SIA_GROUP = (/ SO4, NO3_f, NO3_c, NH4_f /)
  
  integer, public, target, save, dimension (2) :: &
    RDN_GROUP = (/ NH3, NH4_f /)
  
  integer, public, target, save, dimension (2) :: &
    WDEP_RDN_GROUP = (/ NH3, NH4_f /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_RDN_GROUP = (/ NH3, NH4_f /)
  
  integer, public, target, save, dimension (2) :: &
    TNO3_GROUP = (/ NO3_f, NO3_c /)
  
  integer, public, target, save, dimension (2) :: &
    WDEP_TNO3_GROUP = (/ NO3_f, NO3_c /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_TNO3_GROUP = (/ NO3_f, NO3_c /)
  
  integer, public, target, save, dimension (5) :: &
    PMCO_GROUP = (/ NO3_c, POM_c_FFUEL, REMPPM_c, SeaSalt_c, Dust_c /)
  
  integer, public, target, save, dimension (5) :: &
    WDEP_PMCO_GROUP = (/ NO3_c, POM_c_FFUEL, REMPPM_c, SeaSalt_c, Dust_c /)
  
  integer, public, target, save, dimension (5) :: &
    DDEP_PMCO_GROUP = (/ NO3_c, POM_c_FFUEL, REMPPM_c, SeaSalt_c, Dust_c /)
  
  integer, public, target, save, dimension (1) :: &
    TMPX_GROUP = (/ SHIPNOX /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_TMPX_GROUP = (/ SHIPNOX /)
  
  integer, public, target, save, dimension (4) :: &
    MONOTERP_GROUP = (/ BIOTERP, BPINENE, XTERP, SQT_SOA_NV /)
  
  integer, public, target, save, dimension (26) :: &
    OM25_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, OM25_BGND, ASOC_ng100, ASOC_ug1,  &
      ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3, NON_C_ASOA_ng100,  &
      NON_C_ASOA_ug1, NON_C_ASOA_ug10, NON_C_ASOA_ug1e2,  &
      NON_C_ASOA_ug1e3, BSOC_ng100, BSOC_ug1, BSOC_ug10,  &
      BSOC_ug1e2, BSOC_ug1e3, NON_C_BSOA_ng100, NON_C_BSOA_ug1,  &
      NON_C_BSOA_ug10, NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3,  &
      FFFUEL_ng10, WOODOA_ng10, FFIRE_OM  &
    /)
  
  integer, public, target, save, dimension (25) :: &
    WDEP_OM25_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, ASOC_ng100, ASOC_ug1, ASOC_ug10,  &
      ASOC_ug1e2, ASOC_ug1e3, NON_C_ASOA_ng100, NON_C_ASOA_ug1,  &
      NON_C_ASOA_ug10, NON_C_ASOA_ug1e2, NON_C_ASOA_ug1e3,  &
      BSOC_ng100, BSOC_ug1, BSOC_ug10, BSOC_ug1e2, BSOC_ug1e3,  &
      NON_C_BSOA_ng100, NON_C_BSOA_ug1, NON_C_BSOA_ug10,  &
      NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3, FFFUEL_ng10,  &
      WOODOA_ng10, FFIRE_OM  &
    /)
  
  integer, public, target, save, dimension (23) :: &
    DDEP_OM25_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, ASOC_ng100, ASOC_ug1, ASOC_ug10,  &
      ASOC_ug1e2, ASOC_ug1e3, NON_C_ASOA_ng100, NON_C_ASOA_ug1,  &
      NON_C_ASOA_ug10, NON_C_ASOA_ug1e2, NON_C_ASOA_ug1e3,  &
      BSOC_ng100, BSOC_ug1, BSOC_ug10, BSOC_ug1e2, BSOC_ug1e3,  &
      NON_C_BSOA_ng100, NON_C_BSOA_ug1, NON_C_BSOA_ug10,  &
      NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3, FFIRE_OM  &
    /)
  
  integer, public, target, save, dimension (30) :: &
    PCM_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, ASOC_ng100, ASOC_ug1,  &
      ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3, NON_C_ASOA_ng100,  &
      NON_C_ASOA_ug1, NON_C_ASOA_ug10, NON_C_ASOA_ug1e2,  &
      NON_C_ASOA_ug1e3, BSOC_ng100, BSOC_ug1, BSOC_ug10,  &
      BSOC_ug1e2, BSOC_ug1e3, NON_C_BSOA_ng100, NON_C_BSOA_ug1,  &
      NON_C_BSOA_ug10, NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3,  &
      FFFUEL_ng10, WOODOA_ng10, FFIRE_OM, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (30) :: &
    WDEP_PCM_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, ASOC_ng100, ASOC_ug1,  &
      ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3, NON_C_ASOA_ng100,  &
      NON_C_ASOA_ug1, NON_C_ASOA_ug10, NON_C_ASOA_ug1e2,  &
      NON_C_ASOA_ug1e3, BSOC_ng100, BSOC_ug1, BSOC_ug10,  &
      BSOC_ug1e2, BSOC_ug1e3, NON_C_BSOA_ng100, NON_C_BSOA_ug1,  &
      NON_C_BSOA_ug10, NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3,  &
      FFFUEL_ng10, WOODOA_ng10, FFIRE_OM, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (28) :: &
    DDEP_PCM_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, ASOC_ng100, ASOC_ug1,  &
      ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3, NON_C_ASOA_ng100,  &
      NON_C_ASOA_ug1, NON_C_ASOA_ug10, NON_C_ASOA_ug1e2,  &
      NON_C_ASOA_ug1e3, BSOC_ng100, BSOC_ug1, BSOC_ug10,  &
      BSOC_ug1e2, BSOC_ug1e3, NON_C_BSOA_ng100, NON_C_BSOA_ug1,  &
      NON_C_BSOA_ug10, NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3,  &
      FFIRE_OM, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (7) :: &
    PPM25_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, REMPPM25  &
    /)
  
  integer, public, target, save, dimension (7) :: &
    WDEP_PPM25_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, REMPPM25  &
    /)
  
  integer, public, target, save, dimension (7) :: &
    DDEP_PPM25_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_f_FFUEL_new, EC_f_FFUEL_age, REMPPM25  &
    /)
  
  integer, public, target, save, dimension (11) :: &
    PPM10_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, REMPPM25, REMPPM_c  &
    /)
  
  integer, public, target, save, dimension (11) :: &
    WDEP_PPM10_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, REMPPM25, REMPPM_c  &
    /)
  
  integer, public, target, save, dimension (11) :: &
    DDEP_PPM10_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, REMPPM25, REMPPM_c  &
    /)
  
  integer, public, target, save, dimension (1) :: &
    NVWOODOC25_GROUP = (/ POM_f_WOOD /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_NVWOODOC25_GROUP = (/ POM_f_WOOD /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_NVWOODOC25_GROUP = (/ POM_f_WOOD /)
  
  integer, public, target, save, dimension (12) :: &
    NONVOLPCM_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, OM25_BGND, FFIRE_OM, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (11) :: &
    WDEP_NONVOLPCM_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, FFIRE_OM, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (11) :: &
    DDEP_NONVOLPCM_GROUP = (/  &
      POM_f_WOOD, POM_f_FFUEL, POM_c_FFUEL, EC_f_WOOD_new,  &
      EC_f_WOOD_age, EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age,  &
      EC_c_FFUEL, FFIRE_OM, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (4) :: &
    NVABSOM_GROUP = (/ POM_f_WOOD, POM_f_FFUEL, OM25_BGND, FFIRE_OM /)
  
  integer, public, target, save, dimension (3) :: &
    WDEP_NVABSOM_GROUP = (/ POM_f_WOOD, POM_f_FFUEL, FFIRE_OM /)
  
  integer, public, target, save, dimension (3) :: &
    DDEP_NVABSOM_GROUP = (/ POM_f_WOOD, POM_f_FFUEL, FFIRE_OM /)
  
  integer, public, target, save, dimension (1) :: &
    NVFFUELOC25_GROUP = (/ POM_f_FFUEL /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_NVFFUELOC25_GROUP = (/ POM_f_FFUEL /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_NVFFUELOC25_GROUP = (/ POM_f_FFUEL /)
  
  integer, public, target, save, dimension (4) :: &
    PPM_C_GROUP = (/ POM_c_FFUEL, EC_c_WOOD, EC_c_FFUEL, REMPPM_c /)
  
  integer, public, target, save, dimension (4) :: &
    WDEP_PPM_c_GROUP = (/ POM_c_FFUEL, EC_c_WOOD, EC_c_FFUEL, REMPPM_c /)
  
  integer, public, target, save, dimension (4) :: &
    DDEP_PPM_c_GROUP = (/ POM_c_FFUEL, EC_c_WOOD, EC_c_FFUEL, REMPPM_c /)
  
  integer, public, target, save, dimension (3) :: &
    PM10ANTHR_GROUP = (/ POM_c_FFUEL, EC_c_WOOD, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (3) :: &
    WDEP_PM10anthr_GROUP = (/ POM_c_FFUEL, EC_c_WOOD, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (3) :: &
    DDEP_PM10anthr_GROUP = (/ POM_c_FFUEL, EC_c_WOOD, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (1) :: &
    NVFFUELOC_COARSE_GROUP = (/ POM_c_FFUEL /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_NVFFUELOC_COARSE_GROUP = (/ POM_c_FFUEL /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_NVFFUELOC_COARSE_GROUP = (/ POM_c_FFUEL /)
  
  integer, public, target, save, dimension (1) :: &
    OMCOARSE_GROUP = (/ POM_c_FFUEL /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_OMCOARSE_GROUP = (/ POM_c_FFUEL /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_OMCOARSE_GROUP = (/ POM_c_FFUEL /)
  
  integer, public, target, save, dimension (5) :: &
    EC_F_GROUP = (/  &
      EC_f_WOOD_new, EC_f_WOOD_age, EC_f_FFUEL_new,  &
      EC_f_FFUEL_age, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (5) :: &
    WDEP_EC_f_GROUP = (/  &
      EC_f_WOOD_new, EC_f_WOOD_age, EC_f_FFUEL_new,  &
      EC_f_FFUEL_age, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (5) :: &
    DDEP_EC_f_GROUP = (/  &
      EC_f_WOOD_new, EC_f_WOOD_age, EC_f_FFUEL_new,  &
      EC_f_FFUEL_age, FFIRE_BC  &
    /)
  
  integer, public, target, save, dimension (3) :: &
    WOODEC_GROUP = (/ EC_f_WOOD_new, EC_f_WOOD_age, EC_c_WOOD /)
  
  integer, public, target, save, dimension (3) :: &
    WDEP_WOODEC_GROUP = (/ EC_f_WOOD_new, EC_f_WOOD_age, EC_c_WOOD /)
  
  integer, public, target, save, dimension (3) :: &
    DDEP_WOODEC_GROUP = (/ EC_f_WOOD_new, EC_f_WOOD_age, EC_c_WOOD /)
  
  integer, public, target, save, dimension (4) :: &
    ECFINE_GROUP = (/ EC_f_WOOD_new, EC_f_WOOD_age, EC_f_FFUEL_new, EC_f_FFUEL_age /)
  
  integer, public, target, save, dimension (4) :: &
    WDEP_ECFINE_GROUP = (/ EC_f_WOOD_new, EC_f_WOOD_age, EC_f_FFUEL_new, EC_f_FFUEL_age /)
  
  integer, public, target, save, dimension (4) :: &
    DDEP_ECFINE_GROUP = (/ EC_f_WOOD_new, EC_f_WOOD_age, EC_f_FFUEL_new, EC_f_FFUEL_age /)
  
  integer, public, target, save, dimension (2) :: &
    ECCOARSE_GROUP = (/ EC_c_WOOD, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (2) :: &
    WDEP_ECCOARSE_GROUP = (/ EC_c_WOOD, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_ECCOARSE_GROUP = (/ EC_c_WOOD, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (3) :: &
    FFUELEC_GROUP = (/ EC_f_FFUEL_new, EC_f_FFUEL_age, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (3) :: &
    WDEP_FFUELEC_GROUP = (/ EC_f_FFUEL_new, EC_f_FFUEL_age, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (3) :: &
    DDEP_FFUELEC_GROUP = (/ EC_f_FFUEL_new, EC_f_FFUEL_age, EC_c_FFUEL /)
  
  integer, public, target, save, dimension (3) :: &
    OC_GROUP = (/ OM25_p, VBS_TEST, FFIRE_OM /)
  
  integer, public, target, save, dimension (10) :: &
    ASOA_GROUP = (/  &
      ASOC_ng100, ASOC_ug1, ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3,  &
      NON_C_ASOA_ng100, NON_C_ASOA_ug1, NON_C_ASOA_ug10,  &
      NON_C_ASOA_ug1e2, NON_C_ASOA_ug1e3  &
    /)
  
  integer, public, target, save, dimension (10) :: &
    WDEP_ASOA_GROUP = (/  &
      ASOC_ng100, ASOC_ug1, ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3,  &
      NON_C_ASOA_ng100, NON_C_ASOA_ug1, NON_C_ASOA_ug10,  &
      NON_C_ASOA_ug1e2, NON_C_ASOA_ug1e3  &
    /)
  
  integer, public, target, save, dimension (10) :: &
    DDEP_ASOA_GROUP = (/  &
      ASOC_ng100, ASOC_ug1, ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3,  &
      NON_C_ASOA_ng100, NON_C_ASOA_ug1, NON_C_ASOA_ug10,  &
      NON_C_ASOA_ug1e2, NON_C_ASOA_ug1e3  &
    /)
  
  integer, public, target, save, dimension (10) :: &
    BSOA_GROUP = (/  &
      BSOC_ng100, BSOC_ug1, BSOC_ug10, BSOC_ug1e2, BSOC_ug1e3,  &
      NON_C_BSOA_ng100, NON_C_BSOA_ug1, NON_C_BSOA_ug10,  &
      NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3  &
    /)
  
  integer, public, target, save, dimension (10) :: &
    WDEP_BSOA_GROUP = (/  &
      BSOC_ng100, BSOC_ug1, BSOC_ug10, BSOC_ug1e2, BSOC_ug1e3,  &
      NON_C_BSOA_ng100, NON_C_BSOA_ug1, NON_C_BSOA_ug10,  &
      NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3  &
    /)
  
  integer, public, target, save, dimension (10) :: &
    DDEP_BSOA_GROUP = (/  &
      BSOC_ng100, BSOC_ug1, BSOC_ug10, BSOC_ug1e2, BSOC_ug1e3,  &
      NON_C_BSOA_ng100, NON_C_BSOA_ug1, NON_C_BSOA_ug10,  &
      NON_C_BSOA_ug1e2, NON_C_BSOA_ug1e3  &
    /)
  
  integer, public, target, save, dimension (1) :: &
    PFFUELOA25_GROUP = (/ FFFUEL_ng10 /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_PFFUELOA25_GROUP = (/ FFFUEL_ng10 /)
  
  integer, public, target, save, dimension (1) :: &
    SVFFUELOA25_GROUP = (/ FFFUEL_ng10 /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_SVFFUELOA25_GROUP = (/ FFFUEL_ng10 /)
  
  integer, public, target, save, dimension (1) :: &
    PWOODOA25_GROUP = (/ WOODOA_ng10 /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_PWOODOA25_GROUP = (/ WOODOA_ng10 /)
  
  integer, public, target, save, dimension (1) :: &
    SVWOODOA25_GROUP = (/ WOODOA_ng10 /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_SVWOODOA25_GROUP = (/ WOODOA_ng10 /)
  
  integer, public, target, save, dimension (1) :: &
    TRACER_GROUP = (/ FFIRE_CO /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_OC_GROUP = (/ FFIRE_OM /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_OC_GROUP = (/ FFIRE_OM /)
  
  integer, public, target, save, dimension (3) :: &
    PPM25_FIRE_GROUP = (/ FFIRE_OM, FFIRE_BC, FFIRE_REMPPM25 /)
  
  integer, public, target, save, dimension (3) :: &
    WDEP_PPM25_FIRE_GROUP = (/ FFIRE_OM, FFIRE_BC, FFIRE_REMPPM25 /)
  
  integer, public, target, save, dimension (3) :: &
    DDEP_PPM25_FIRE_GROUP = (/ FFIRE_OM, FFIRE_BC, FFIRE_REMPPM25 /)
  
  integer, public, target, save, dimension (1) :: &
    NVFFIREOC25_GROUP = (/ FFIRE_OM /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_NVFFIREOC25_GROUP = (/ FFIRE_OM /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_NVFFIREOC25_GROUP = (/ FFIRE_OM /)
  
  integer, public, target, save, dimension (1) :: &
    FFIREBC_GROUP = (/ FFIRE_BC /)
  
  integer, public, target, save, dimension (1) :: &
    WDEP_FFIREBC_GROUP = (/ FFIRE_BC /)
  
  integer, public, target, save, dimension (1) :: &
    DDEP_FFIREBC_GROUP = (/ FFIRE_BC /)
  
  integer, public, target, save, dimension (2) :: &
    SS_GROUP = (/ SeaSalt_f, SeaSalt_c /)
  
  integer, public, target, save, dimension (2) :: &
    WDEP_SS_GROUP = (/ SeaSalt_f, SeaSalt_c /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_SS_GROUP = (/ SeaSalt_f, SeaSalt_c /)
  
  integer, public, target, save, dimension (2) :: &
    DUST_GROUP = (/ Dust_f, Dust_c /)
  
  integer, public, target, save, dimension (2) :: &
    WDEP_DUST_GROUP = (/ Dust_f, Dust_c /)
  
  integer, public, target, save, dimension (2) :: &
    DDEP_DUST_GROUP = (/ Dust_f, Dust_c /)
  
  integer, public, target, save, dimension (22) :: &
    CSTAR_GROUP = (/  &
      ASOC_ng100, ASOC_ug1, ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3,  &
      NON_C_ASOA_ng100, NON_C_ASOA_ug1, NON_C_ASOA_ug10,  &
      NON_C_ASOA_ug1e2, NON_C_ASOA_ug1e3, BSOC_ng100, BSOC_ug1,  &
      BSOC_ug10, BSOC_ug1e2, BSOC_ug1e3, NON_C_BSOA_ng100,  &
      NON_C_BSOA_ug1, NON_C_BSOA_ug10, NON_C_BSOA_ug1e2,  &
      NON_C_BSOA_ug1e3, FFFUEL_ng10, WOODOA_ng10  &
    /)
  
  real, public, target, save, dimension (22) :: &
    CSTAR_GROUP_FACTORS = (/  &
      0.1, 1.0, 10.0, 1.0e2, 1.0e3, 0.1, 1.0, 10.0, 1.0e2, 1.0e3,  &
      0.1, 1.0, 10.0, 1.0e2, 1.0e3, 0.1, 1.0, 10.0, 1.0e2, 1.0e3,  &
      0.01, 0.01  &
    /)
  
  integer, public, target, save, dimension (22) :: &
    DELTAH_GROUP = (/  &
      ASOC_ng100, ASOC_ug1, ASOC_ug10, ASOC_ug1e2, ASOC_ug1e3,  &
      NON_C_ASOA_ng100, NON_C_ASOA_ug1, NON_C_ASOA_ug10,  &
      NON_C_ASOA_ug1e2, NON_C_ASOA_ug1e3, BSOC_ng100, BSOC_ug1,  &
      BSOC_ug10, BSOC_ug1e2, BSOC_ug1e3, NON_C_BSOA_ng100,  &
      NON_C_BSOA_ug1, NON_C_BSOA_ug10, NON_C_BSOA_ug1e2,  &
      NON_C_BSOA_ug1e3, FFFUEL_ng10, WOODOA_ng10  &
    /)
  
  real, public, target, save, dimension (22) :: &
    DELTAH_GROUP_FACTORS = (/  &
      30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0,  &
      30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0,  &
      112.0, 112.0  &
    /)
  
  integer, public, target, save, dimension (14) :: &
    EXTINC_GROUP = (/  &
      SO4, NO3_f, NO3_c, NH4_f, EC_f_WOOD_new, EC_f_WOOD_age,  &
      EC_c_WOOD, EC_f_FFUEL_new, EC_f_FFUEL_age, EC_c_FFUEL,  &
      REMPPM25, REMPPM_c, FFIRE_BC, FFIRE_REMPPM25  &
    /)
  
  character(len=TXTLEN_SHORT), public, target, save, dimension (14) :: &
    EXTINC_GROUP_MAPBACK = [ character(len=TXTLEN_SHORT) :: &
    "SO4", "NO3f", "NO3c", "NH4f", "ECn", "ECa", "EC", "ECn",  &
      "ECa", "EC", "DDf", "DDc", "EC", "DDf"  &
    ]

contains
  
  subroutine Init_ChemGroups()
    
    chemgroups(1)%name="RO2"
    chemgroups(1)%specs=>RO2_GROUP
    
    chemgroups(2)%name="OX"
    chemgroups(2)%specs=>OX_GROUP
    
    chemgroups(3)%name="DDEP_OX"
    chemgroups(3)%specs=>DDEP_OX_GROUP
    
    chemgroups(4)%name="NOX"
    chemgroups(4)%specs=>NOX_GROUP
    
    chemgroups(5)%name="OXN"
    chemgroups(5)%specs=>OXN_GROUP
    
    chemgroups(6)%name="DDEP_NOX"
    chemgroups(6)%specs=>DDEP_NOX_GROUP
    
    chemgroups(7)%name="DDEP_OXN"
    chemgroups(7)%specs=>DDEP_OXN_GROUP
    
    chemgroups(8)%name="DAOBS"
    chemgroups(8)%specs=>DAOBS_GROUP
    
    chemgroups(9)%name="DDEP_daObs"
    chemgroups(9)%specs=>DDEP_daObs_GROUP
    
    chemgroups(10)%name="PAN"
    chemgroups(10)%specs=>PAN_GROUP
    
    chemgroups(11)%name="DDEP_PAN"
    chemgroups(11)%specs=>DDEP_PAN_GROUP
    
    chemgroups(12)%name="WDEP_OXN"
    chemgroups(12)%specs=>WDEP_OXN_GROUP
    
    chemgroups(13)%name="BVOC"
    chemgroups(13)%specs=>BVOC_GROUP
    
    chemgroups(14)%name="ROOH"
    chemgroups(14)%specs=>ROOH_GROUP
    
    chemgroups(15)%name="DDEP_ROOH"
    chemgroups(15)%specs=>DDEP_ROOH_GROUP
    
    chemgroups(16)%name="WDEP_ROOH"
    chemgroups(16)%specs=>WDEP_ROOH_GROUP
    
    chemgroups(17)%name="SOX"
    chemgroups(17)%specs=>SOX_GROUP
    
    chemgroups(18)%name="WDEP_SOX"
    chemgroups(18)%specs=>WDEP_SOX_GROUP
    
    chemgroups(19)%name="DDEP_SOX"
    chemgroups(19)%specs=>DDEP_SOX_GROUP
    
    chemgroups(20)%name="PM10"
    chemgroups(20)%specs=>PM10_GROUP
    
    chemgroups(21)%name="WDEP_PM10"
    chemgroups(21)%specs=>WDEP_PM10_GROUP
    
    chemgroups(22)%name="DDEP_PM10"
    chemgroups(22)%specs=>DDEP_PM10_GROUP
    
    chemgroups(23)%name="PMFINE"
    chemgroups(23)%specs=>PMFINE_GROUP
    
    chemgroups(24)%name="WDEP_PMFINE"
    chemgroups(24)%specs=>WDEP_PMFINE_GROUP
    
    chemgroups(25)%name="DDEP_PMFINE"
    chemgroups(25)%specs=>DDEP_PMFINE_GROUP
    
    chemgroups(26)%name="SIA"
    chemgroups(26)%specs=>SIA_GROUP
    
    chemgroups(27)%name="WDEP_SIA"
    chemgroups(27)%specs=>WDEP_SIA_GROUP
    
    chemgroups(28)%name="DDEP_SIA"
    chemgroups(28)%specs=>DDEP_SIA_GROUP
    
    chemgroups(29)%name="RDN"
    chemgroups(29)%specs=>RDN_GROUP
    
    chemgroups(30)%name="WDEP_RDN"
    chemgroups(30)%specs=>WDEP_RDN_GROUP
    
    chemgroups(31)%name="DDEP_RDN"
    chemgroups(31)%specs=>DDEP_RDN_GROUP
    
    chemgroups(32)%name="TNO3"
    chemgroups(32)%specs=>TNO3_GROUP
    
    chemgroups(33)%name="WDEP_TNO3"
    chemgroups(33)%specs=>WDEP_TNO3_GROUP
    
    chemgroups(34)%name="DDEP_TNO3"
    chemgroups(34)%specs=>DDEP_TNO3_GROUP
    
    chemgroups(35)%name="PMCO"
    chemgroups(35)%specs=>PMCO_GROUP
    
    chemgroups(36)%name="WDEP_PMCO"
    chemgroups(36)%specs=>WDEP_PMCO_GROUP
    
    chemgroups(37)%name="DDEP_PMCO"
    chemgroups(37)%specs=>DDEP_PMCO_GROUP
    
    chemgroups(38)%name="TMPX"
    chemgroups(38)%specs=>TMPX_GROUP
    
    chemgroups(39)%name="DDEP_TMPX"
    chemgroups(39)%specs=>DDEP_TMPX_GROUP
    
    chemgroups(40)%name="MONOTERP"
    chemgroups(40)%specs=>MONOTERP_GROUP
    
    chemgroups(41)%name="OM25"
    chemgroups(41)%specs=>OM25_GROUP
    
    chemgroups(42)%name="WDEP_OM25"
    chemgroups(42)%specs=>WDEP_OM25_GROUP
    
    chemgroups(43)%name="DDEP_OM25"
    chemgroups(43)%specs=>DDEP_OM25_GROUP
    
    chemgroups(44)%name="PCM"
    chemgroups(44)%specs=>PCM_GROUP
    
    chemgroups(45)%name="WDEP_PCM"
    chemgroups(45)%specs=>WDEP_PCM_GROUP
    
    chemgroups(46)%name="DDEP_PCM"
    chemgroups(46)%specs=>DDEP_PCM_GROUP
    
    chemgroups(47)%name="PPM25"
    chemgroups(47)%specs=>PPM25_GROUP
    
    chemgroups(48)%name="WDEP_PPM25"
    chemgroups(48)%specs=>WDEP_PPM25_GROUP
    
    chemgroups(49)%name="DDEP_PPM25"
    chemgroups(49)%specs=>DDEP_PPM25_GROUP
    
    chemgroups(50)%name="PPM10"
    chemgroups(50)%specs=>PPM10_GROUP
    
    chemgroups(51)%name="WDEP_PPM10"
    chemgroups(51)%specs=>WDEP_PPM10_GROUP
    
    chemgroups(52)%name="DDEP_PPM10"
    chemgroups(52)%specs=>DDEP_PPM10_GROUP
    
    chemgroups(53)%name="NVWOODOC25"
    chemgroups(53)%specs=>NVWOODOC25_GROUP
    
    chemgroups(54)%name="WDEP_NVWOODOC25"
    chemgroups(54)%specs=>WDEP_NVWOODOC25_GROUP
    
    chemgroups(55)%name="DDEP_NVWOODOC25"
    chemgroups(55)%specs=>DDEP_NVWOODOC25_GROUP
    
    chemgroups(56)%name="NONVOLPCM"
    chemgroups(56)%specs=>NONVOLPCM_GROUP
    
    chemgroups(57)%name="WDEP_NONVOLPCM"
    chemgroups(57)%specs=>WDEP_NONVOLPCM_GROUP
    
    chemgroups(58)%name="DDEP_NONVOLPCM"
    chemgroups(58)%specs=>DDEP_NONVOLPCM_GROUP
    
    chemgroups(59)%name="NVABSOM"
    chemgroups(59)%specs=>NVABSOM_GROUP
    
    chemgroups(60)%name="WDEP_NVABSOM"
    chemgroups(60)%specs=>WDEP_NVABSOM_GROUP
    
    chemgroups(61)%name="DDEP_NVABSOM"
    chemgroups(61)%specs=>DDEP_NVABSOM_GROUP
    
    chemgroups(62)%name="NVFFUELOC25"
    chemgroups(62)%specs=>NVFFUELOC25_GROUP
    
    chemgroups(63)%name="WDEP_NVFFUELOC25"
    chemgroups(63)%specs=>WDEP_NVFFUELOC25_GROUP
    
    chemgroups(64)%name="DDEP_NVFFUELOC25"
    chemgroups(64)%specs=>DDEP_NVFFUELOC25_GROUP
    
    chemgroups(65)%name="PPM_C"
    chemgroups(65)%specs=>PPM_C_GROUP
    
    chemgroups(66)%name="WDEP_PPM_c"
    chemgroups(66)%specs=>WDEP_PPM_c_GROUP
    
    chemgroups(67)%name="DDEP_PPM_c"
    chemgroups(67)%specs=>DDEP_PPM_c_GROUP
    
    chemgroups(68)%name="PM10ANTHR"
    chemgroups(68)%specs=>PM10ANTHR_GROUP
    
    chemgroups(69)%name="WDEP_PM10anthr"
    chemgroups(69)%specs=>WDEP_PM10anthr_GROUP
    
    chemgroups(70)%name="DDEP_PM10anthr"
    chemgroups(70)%specs=>DDEP_PM10anthr_GROUP
    
    chemgroups(71)%name="NVFFUELOC_COARSE"
    chemgroups(71)%specs=>NVFFUELOC_COARSE_GROUP
    
    chemgroups(72)%name="WDEP_NVFFUELOC_COARSE"
    chemgroups(72)%specs=>WDEP_NVFFUELOC_COARSE_GROUP
    
    chemgroups(73)%name="DDEP_NVFFUELOC_COARSE"
    chemgroups(73)%specs=>DDEP_NVFFUELOC_COARSE_GROUP
    
    chemgroups(74)%name="OMCOARSE"
    chemgroups(74)%specs=>OMCOARSE_GROUP
    
    chemgroups(75)%name="WDEP_OMCOARSE"
    chemgroups(75)%specs=>WDEP_OMCOARSE_GROUP
    
    chemgroups(76)%name="DDEP_OMCOARSE"
    chemgroups(76)%specs=>DDEP_OMCOARSE_GROUP
    
    chemgroups(77)%name="EC_F"
    chemgroups(77)%specs=>EC_F_GROUP
    
    chemgroups(78)%name="WDEP_EC_f"
    chemgroups(78)%specs=>WDEP_EC_f_GROUP
    
    chemgroups(79)%name="DDEP_EC_f"
    chemgroups(79)%specs=>DDEP_EC_f_GROUP
    
    chemgroups(80)%name="WOODEC"
    chemgroups(80)%specs=>WOODEC_GROUP
    
    chemgroups(81)%name="WDEP_WOODEC"
    chemgroups(81)%specs=>WDEP_WOODEC_GROUP
    
    chemgroups(82)%name="DDEP_WOODEC"
    chemgroups(82)%specs=>DDEP_WOODEC_GROUP
    
    chemgroups(83)%name="ECFINE"
    chemgroups(83)%specs=>ECFINE_GROUP
    
    chemgroups(84)%name="WDEP_ECFINE"
    chemgroups(84)%specs=>WDEP_ECFINE_GROUP
    
    chemgroups(85)%name="DDEP_ECFINE"
    chemgroups(85)%specs=>DDEP_ECFINE_GROUP
    
    chemgroups(86)%name="ECCOARSE"
    chemgroups(86)%specs=>ECCOARSE_GROUP
    
    chemgroups(87)%name="WDEP_ECCOARSE"
    chemgroups(87)%specs=>WDEP_ECCOARSE_GROUP
    
    chemgroups(88)%name="DDEP_ECCOARSE"
    chemgroups(88)%specs=>DDEP_ECCOARSE_GROUP
    
    chemgroups(89)%name="FFUELEC"
    chemgroups(89)%specs=>FFUELEC_GROUP
    
    chemgroups(90)%name="WDEP_FFUELEC"
    chemgroups(90)%specs=>WDEP_FFUELEC_GROUP
    
    chemgroups(91)%name="DDEP_FFUELEC"
    chemgroups(91)%specs=>DDEP_FFUELEC_GROUP
    
    chemgroups(92)%name="OC"
    chemgroups(92)%specs=>OC_GROUP
    
    chemgroups(93)%name="ASOA"
    chemgroups(93)%specs=>ASOA_GROUP
    
    chemgroups(94)%name="WDEP_ASOA"
    chemgroups(94)%specs=>WDEP_ASOA_GROUP
    
    chemgroups(95)%name="DDEP_ASOA"
    chemgroups(95)%specs=>DDEP_ASOA_GROUP
    
    chemgroups(96)%name="BSOA"
    chemgroups(96)%specs=>BSOA_GROUP
    
    chemgroups(97)%name="WDEP_BSOA"
    chemgroups(97)%specs=>WDEP_BSOA_GROUP
    
    chemgroups(98)%name="DDEP_BSOA"
    chemgroups(98)%specs=>DDEP_BSOA_GROUP
    
    chemgroups(99)%name="PFFUELOA25"
    chemgroups(99)%specs=>PFFUELOA25_GROUP
    
    chemgroups(100)%name="WDEP_PFFUELOA25"
    chemgroups(100)%specs=>WDEP_PFFUELOA25_GROUP
    
    chemgroups(101)%name="SVFFUELOA25"
    chemgroups(101)%specs=>SVFFUELOA25_GROUP
    
    chemgroups(102)%name="WDEP_SVFFUELOA25"
    chemgroups(102)%specs=>WDEP_SVFFUELOA25_GROUP
    
    chemgroups(103)%name="PWOODOA25"
    chemgroups(103)%specs=>PWOODOA25_GROUP
    
    chemgroups(104)%name="WDEP_PWOODOA25"
    chemgroups(104)%specs=>WDEP_PWOODOA25_GROUP
    
    chemgroups(105)%name="SVWOODOA25"
    chemgroups(105)%specs=>SVWOODOA25_GROUP
    
    chemgroups(106)%name="WDEP_SVWOODOA25"
    chemgroups(106)%specs=>WDEP_SVWOODOA25_GROUP
    
    chemgroups(107)%name="TRACER"
    chemgroups(107)%specs=>TRACER_GROUP
    
    chemgroups(108)%name="WDEP_OC"
    chemgroups(108)%specs=>WDEP_OC_GROUP
    
    chemgroups(109)%name="DDEP_OC"
    chemgroups(109)%specs=>DDEP_OC_GROUP
    
    chemgroups(110)%name="PPM25_FIRE"
    chemgroups(110)%specs=>PPM25_FIRE_GROUP
    
    chemgroups(111)%name="WDEP_PPM25_FIRE"
    chemgroups(111)%specs=>WDEP_PPM25_FIRE_GROUP
    
    chemgroups(112)%name="DDEP_PPM25_FIRE"
    chemgroups(112)%specs=>DDEP_PPM25_FIRE_GROUP
    
    chemgroups(113)%name="NVFFIREOC25"
    chemgroups(113)%specs=>NVFFIREOC25_GROUP
    
    chemgroups(114)%name="WDEP_NVFFIREOC25"
    chemgroups(114)%specs=>WDEP_NVFFIREOC25_GROUP
    
    chemgroups(115)%name="DDEP_NVFFIREOC25"
    chemgroups(115)%specs=>DDEP_NVFFIREOC25_GROUP
    
    chemgroups(116)%name="FFIREBC"
    chemgroups(116)%specs=>FFIREBC_GROUP
    
    chemgroups(117)%name="WDEP_FFIREBC"
    chemgroups(117)%specs=>WDEP_FFIREBC_GROUP
    
    chemgroups(118)%name="DDEP_FFIREBC"
    chemgroups(118)%specs=>DDEP_FFIREBC_GROUP
    
    chemgroups(119)%name="SS"
    chemgroups(119)%specs=>SS_GROUP
    
    chemgroups(120)%name="WDEP_SS"
    chemgroups(120)%specs=>WDEP_SS_GROUP
    
    chemgroups(121)%name="DDEP_SS"
    chemgroups(121)%specs=>DDEP_SS_GROUP
    
    chemgroups(122)%name="DUST"
    chemgroups(122)%specs=>DUST_GROUP
    
    chemgroups(123)%name="WDEP_DUST"
    chemgroups(123)%specs=>WDEP_DUST_GROUP
    
    chemgroups(124)%name="DDEP_DUST"
    chemgroups(124)%specs=>DDEP_DUST_GROUP
    
    chemgroups_factors(1)%name="CSTAR"
    chemgroups_factors(1)%species=>CSTAR_GROUP
    chemgroups_factors(1)%factors=>CSTAR_GROUP_FACTORS
    
    chemgroups_factors(2)%name="DELTAH"
    chemgroups_factors(2)%species=>DELTAH_GROUP
    chemgroups_factors(2)%factors=>DELTAH_GROUP_FACTORS
    
    chemgroups_maps(1)%name="EXTINC"
    chemgroups_maps(1)%species=>EXTINC_GROUP
    chemgroups_maps(1)%maps=>EXTINC_GROUP_MAPBACK
  
  end subroutine Init_ChemGroups

end module ChemGroups_mod
