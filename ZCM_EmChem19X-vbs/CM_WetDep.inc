type(dep_t), public, dimension(NWETDEP_ADV), save :: CM_WDepMap = (/ &
!        species     surrogate      setrate
  dep_t("N2O5      ", "HNO3      ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("HONO      ", "H2O2      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HO2NO2    ", "0p8       ", -999.0 ) &
, dep_t("OXYCATECH ", "H2O2      ", -999.0 ) &
, dep_t("HCHO      ", "ROOH      ", -999.0 ) &
, dep_t("HOCH2CHO  ", "0p6       ", -999.0 ) &
, dep_t("HCOCH2OOH ", "0p3       ", -999.0 ) &
, dep_t("ACETOL    ", "HCHO      ", -999.0 ) &
, dep_t("GLYOX     ", "H2O2      ", -999.0 ) &
, dep_t("MGLYOX    ", "0p3       ", -999.0 ) &
, dep_t("EGLYOX    ", "ROOH      ", -999.0 ) &
, dep_t("CO2H3CHO  ", "H2O2      ", -999.0 ) &
, dep_t("CO23C3CHO ", "H2O2      ", -999.0 ) &
, dep_t("BIACETOH  ", "0p4       ", -999.0 ) &
, dep_t("C5DICARB  ", "H2O2      ", -999.0 ) &
, dep_t("HOCH2CO2H ", "H2O2      ", -999.0 ) &
, dep_t("CH3COCO2H ", "HNO3      ", -999.0 ) &
, dep_t("MACO2H    ", "HCHO      ", -999.0 ) &
, dep_t("ETRO2H    ", "H2O2      ", -999.0 ) &
, dep_t("MEKO2H    ", "1p3       ", -999.0 ) &
, dep_t("ISRO2H    ", "H2O2      ", -999.0 ) &
, dep_t("C5DICAROOH", "H2O2      ", -999.0 ) &
, dep_t("OCATEC1OOH", "H2O2      ", -999.0 ) &
, dep_t("HPALD     ", "H2O2      ", -999.0 ) &
, dep_t("MVKO2H    ", "H2O2      ", -999.0 ) &
, dep_t("MACROOH   ", "H2O2      ", -999.0 ) &
, dep_t("OXYBPEROOH", "H2O2      ", -999.0 ) &
, dep_t("CH3CO3H   ", "ROOH      ", -999.0 ) &
, dep_t("CO2H3CO3H ", "H2O2      ", -999.0 ) &
, dep_t("IEPOX     ", "HNO3      ", -999.0 ) &
, dep_t("MVKNO3    ", "H2O2      ", -999.0 ) &
, dep_t("OXNCATECH ", "H2O2      ", -999.0 ) &
, dep_t("NOXYOL1OOH", "1p2       ", -999.0 ) &
, dep_t("CO2N3CHO  ", "H2O2      ", -999.0 ) &
, dep_t("NALD      ", "ROOH      ", -999.0 ) &
, dep_t("NOA       ", "ROOH      ", -999.0 ) &
, dep_t("ISON      ", "0p7       ", -999.0 ) &
, dep_t("ISON2     ", "ROOH      ", -999.0 ) &
, dep_t("NC4MDCO2H ", "H2O2      ", -999.0 ) &
, dep_t("NOXYOLOOH ", "HNO3      ", -999.0 ) &
, dep_t("OXNCATCOOH", "HNO3      ", -999.0 ) &
, dep_t("C4PAN6    ", "1p1       ", -999.0 ) &
, dep_t("C5PAN18   ", "ROOH      ", -999.0 ) &
, dep_t("OXNNCATOOH", "HNO3      ", -999.0 ) &
, dep_t("TERPOOH   ", "H2O2      ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
, dep_t("ASOC_ng1e2", "ROOH      ", -999.0 ) &
, dep_t("ASOC_ug1  ", "ROOH      ", -999.0 ) &
, dep_t("ASOC_ug10 ", "ROOH      ", -999.0 ) &
, dep_t("ASOC_ug1e2", "ROOH      ", -999.0 ) &
, dep_t("ASOC_ug1e3", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ng1e2", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ug1", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ug10", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ug1e2", "ROOH      ", -999.0 ) &
, dep_t("non_C_ASOA_ug1e3", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ng1e2", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ug1  ", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ug10 ", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ug1e2", "ROOH      ", -999.0 ) &
, dep_t("BSOC_ug1e3", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ng1e2", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ug1", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ug10", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ug1e2", "ROOH      ", -999.0 ) &
, dep_t("non_C_BSOA_ug1e3", "ROOH      ", -999.0 ) &
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
, dep_t("SQT_SOA_NV", "PMf       ", -999.0 ) &
/)
