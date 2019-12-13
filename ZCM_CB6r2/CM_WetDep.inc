type(dep_t), public, dimension(NWETDEP_ADV), save :: CM_WDepMap = (/ &
!        species     surrogate      setrate
  dep_t("AACD      ", "H2O2      ", -999.0 ) &
, dep_t("EPOX      ", "HNO3      ", -999.0 ) &
, dep_t("FACD      ", "H2O2      ", -999.0 ) &
, dep_t("FORM      ", "ROOH      ", -999.0 ) &
, dep_t("GLY       ", "H2O2      ", -999.0 ) &
, dep_t("GLYD      ", "0p6       ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HONO      ", "H2O2      ", -999.0 ) &
, dep_t("HPLD      ", "H2O2      ", -999.0 ) &
, dep_t("INTR      ", "0p7       ", -999.0 ) &
, dep_t("ISPX      ", "H2O2      ", -999.0 ) &
, dep_t("MGLY      ", "0p3       ", -999.0 ) &
, dep_t("N2O5      ", "HNO3      ", -999.0 ) &
, dep_t("OPEN      ", "H2O2      ", -999.0 ) &
, dep_t("PACD      ", "ROOH      ", -999.0 ) &
, dep_t("PNA       ", "0p8       ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("XOPN      ", "H2O2      ", -999.0 ) &
, dep_t("NTR2      ", "HCHO      ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("POM_f_wood", "PMf       ", -999.0 ) &
, dep_t("POM_c_wood", "PMc       ", -999.0 ) &
, dep_t("POM_f_ffuel", "PMf       ", -999.0 ) &
, dep_t("POM_c_ffuel", "PMc       ", -999.0 ) &
, dep_t("EC_f_wood_new", "ECfn      ", -999.0 ) &
, dep_t("EC_f_wood_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_wood ", "PMc       ", -999.0 ) &
, dep_t("EC_f_ffuel_new", "ECfn      ", -999.0 ) &
, dep_t("EC_f_ffuel_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_ffuel", "PMc       ", -999.0 ) &
, dep_t("pSO4f     ", "PMf       ", -999.0 ) &
, dep_t("pSO4c     ", "PMc       ", -999.0 ) &
, dep_t("remPPM25  ", "PMf       ", -999.0 ) &
, dep_t("remPPM_c  ", "PMc       ", -999.0 ) &
, dep_t("ASOC_ng10 ", "H2O2      ", -999.0 ) &
, dep_t("ASOC_ng1e2", "H2O2      ", -999.0 ) &
, dep_t("ASOC_ug1  ", "H2O2      ", -999.0 ) &
, dep_t("ASOC_ug10 ", "H2O2      ", -999.0 ) &
, dep_t("ASOC_ug1e2", "H2O2      ", -999.0 ) &
, dep_t("ASOC_ug1e3", "H2O2      ", -999.0 ) &
, dep_t("NON_C_ASOA_ng10", "H2O2      ", -999.0 ) &
, dep_t("NON_C_ASOA_ng1e2", "H2O2      ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1", "H2O2      ", -999.0 ) &
, dep_t("NON_C_ASOA_ug10", "H2O2      ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1e2", "H2O2      ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1e3", "H2O2      ", -999.0 ) &
, dep_t("BSOC_ng10 ", "HNO3      ", -999.0 ) &
, dep_t("BSOC_ng1e2", "HNO3      ", -999.0 ) &
, dep_t("BSOC_ug1  ", "H2O2      ", -999.0 ) &
, dep_t("BSOC_ug10 ", "H2O2      ", -999.0 ) &
, dep_t("BSOC_ug1e2", "H2O2      ", -999.0 ) &
, dep_t("BSOC_ug1e3", "H2O2      ", -999.0 ) &
, dep_t("NON_C_BSOA_ng10", "HNO3      ", -999.0 ) &
, dep_t("NON_C_BSOA_ng1e2", "HNO3      ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1", "H2O2      ", -999.0 ) &
, dep_t("NON_C_BSOA_ug10", "H2O2      ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1e2", "H2O2      ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1e3", "H2O2      ", -999.0 ) &
, dep_t("ffire_OM  ", "PMf       ", -999.0 ) &
, dep_t("ffire_BC  ", "PMf       ", -999.0 ) &
, dep_t("ffire_remPPM25", "PMf       ", -999.0 ) &
, dep_t("ffire_c   ", "PMc       ", -999.0 ) &
, dep_t("SeaSalt_f ", "SSf       ", -999.0 ) &
, dep_t("SeaSalt_c ", "SSc       ", -999.0 ) &
, dep_t("Dust_road_f", "PMf       ", -999.0 ) &
, dep_t("Dust_road_c", "PMc       ", -999.0 ) &
, dep_t("Dust_wb_f ", "PMf       ", -999.0 ) &
, dep_t("Dust_wb_c ", "PMc       ", -999.0 ) &
, dep_t("Dust_sah_f", "PMf       ", -999.0 ) &
, dep_t("Dust_sah_c", "PMc       ", -999.0 ) &
, dep_t("Ash_f     ", "PMf       ", -999.0 ) &
, dep_t("Ash_c     ", "PMc       ", -999.0 ) &
, dep_t("SQT_SOA_NV", "PMf       ", -999.0 ) &
/)
