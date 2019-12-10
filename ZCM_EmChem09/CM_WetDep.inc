type(dep_t), public, dimension(NWETDEP_ADV), save :: CM_WDepMap = (/ &
!        species     surrogate      setrate
  dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HONO      ", "HNO3      ", -999.0 ) &
, dep_t("HCHO      ", "HCHO      ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMf       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("POM_f_wood", "PMf       ", -999.0 ) &
, dep_t("POM_f_ffuel", "PMf       ", -999.0 ) &
, dep_t("POM_c_ffuel", "PMc       ", -999.0 ) &
, dep_t("EC_f_wood_new", "ECfn      ", -999.0 ) &
, dep_t("EC_f_wood_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_wood ", "PMc       ", -999.0 ) &
, dep_t("EC_f_ffuel_new", "ECfn      ", -999.0 ) &
, dep_t("EC_f_ffuel_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_ffuel", "PMc       ", -999.0 ) &
, dep_t("REMPPM25  ", "PMf       ", -999.0 ) &
, dep_t("REMPPM_c  ", "PMc       ", -999.0 ) &
, dep_t("ASOC_ng100", "PMf       ", -999.0 ) &
, dep_t("ASOC_ug1  ", "PMf       ", -999.0 ) &
, dep_t("ASOC_ug10 ", "PMf       ", -999.0 ) &
, dep_t("ASOC_ug1e2", "PMf       ", -999.0 ) &
, dep_t("ASOC_ug1e3", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ng100", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ug10", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1e2", "PMf       ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1e3", "PMf       ", -999.0 ) &
, dep_t("BSOC_ng100", "PMf       ", -999.0 ) &
, dep_t("BSOC_ug1  ", "PMf       ", -999.0 ) &
, dep_t("BSOC_ug10 ", "PMf       ", -999.0 ) &
, dep_t("BSOC_ug1e2", "PMf       ", -999.0 ) &
, dep_t("BSOC_ug1e3", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ng100", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ug10", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1e2", "PMf       ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1e3", "PMf       ", -999.0 ) &
, dep_t("Fffuel_ng10", "PMf       ", -999.0 ) &
, dep_t("woodOA_ng10", "PMf       ", -999.0 ) &
, dep_t("SQT_SOA_nv", "PMf       ", -999.0 ) &
, dep_t("TERPOOH   ", "H2O2      ", -999.0 ) &
, dep_t("Ash_f     ", "PMf       ", -999.0 ) &
, dep_t("Ash_c     ", "PMc       ", -999.0 ) &
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
, dep_t("POLLEN_BIRCH", "POLLw     ", -999.0 ) &
, dep_t("POLLEN_OLIVE", "POLLw     ", -999.0 ) &
, dep_t("POLLEN_RWEED", "POLLw     ", -999.0 ) &
, dep_t("POLLEN_GRASS", "POLLw     ", -999.0 ) &
/)
