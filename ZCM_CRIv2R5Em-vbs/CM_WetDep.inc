type(dep_t), public, dimension(NWETDEP_ADV), save :: CM_WDepMap = (/ &
!        species     surrogate      setrate
  dep_t("N2O5      ", "HNO3      ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("HONO      ", "H2O2      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HO2NO2    ", "0p8       ", -999.0 ) &
, dep_t("CH3CO2H   ", "H2O2      ", -999.0 ) &
, dep_t("HCHO      ", "ROOH      ", -999.0 ) &
, dep_t("HCOOH     ", "H2O2      ", -999.0 ) &
, dep_t("HCOCO2H   ", "HCHO      ", -999.0 ) &
, dep_t("HCOCO3H   ", "H2O2      ", -999.0 ) &
, dep_t("GLYOX     ", "H2O2      ", -999.0 ) &
, dep_t("MGLYOX    ", "0p3       ", -999.0 ) &
, dep_t("ACETOL    ", "HCHO      ", -999.0 ) &
, dep_t("RN19NO3   ", "0p8       ", -999.0 ) &
, dep_t("HOC2H4NO3 ", "0p6       ", -999.0 ) &
, dep_t("RN9NO3    ", "HCHO      ", -999.0 ) &
, dep_t("RN12NO3   ", "HCHO      ", -999.0 ) &
, dep_t("RN15NO3   ", "0p6       ", -999.0 ) &
, dep_t("RN19OOH   ", "H2O2      ", -999.0 ) &
, dep_t("CH3CO3H   ", "ROOH      ", -999.0 ) &
, dep_t("RN8OOH    ", "0p4       ", -999.0 ) &
, dep_t("RN11OOH   ", "1p3       ", -999.0 ) &
, dep_t("RN14OOH   ", "H2O2      ", -999.0 ) &
, dep_t("RN17OOH   ", "HNO3      ", -999.0 ) &
, dep_t("HOC2H4OOH ", "H2O2      ", -999.0 ) &
, dep_t("RN9OOH    ", "H2O2      ", -999.0 ) &
, dep_t("RN12OOH   ", "H2O2      ", -999.0 ) &
, dep_t("RN15OOH   ", "HNO3      ", -999.0 ) &
, dep_t("NRN6OOH   ", "0p6       ", -999.0 ) &
, dep_t("NRN9OOH   ", "0p4       ", -999.0 ) &
, dep_t("NRN12OOH  ", "0p3       ", -999.0 ) &
, dep_t("RU13OOH   ", "H2O2      ", -999.0 ) &
, dep_t("CARB13    ", "H2O2      ", -999.0 ) &
, dep_t("CARB12    ", "0p4       ", -999.0 ) &
, dep_t("HOCH2CHO  ", "0p6       ", -999.0 ) &
, dep_t("HOCH2CO3H ", "ROOH      ", -999.0 ) &
, dep_t("CCARB12   ", "H2O2      ", -999.0 ) &
, dep_t("UCARB12   ", "H2O2      ", -999.0 ) &
, dep_t("RU14NO3   ", "0p7       ", -999.0 ) &
, dep_t("RU14OOH   ", "H2O2      ", -999.0 ) &
, dep_t("RU12OOH   ", "H2O2      ", -999.0 ) &
, dep_t("RU10OOH   ", "H2O2      ", -999.0 ) &
, dep_t("NUCARB12  ", "ROOH      ", -999.0 ) &
, dep_t("NRU14OOH  ", "0p2       ", -999.0 ) &
, dep_t("NRU12OOH  ", "H2O2      ", -999.0 ) &
, dep_t("NOA       ", "ROOH      ", -999.0 ) &
, dep_t("RA13NO3   ", "H2O2      ", -999.0 ) &
, dep_t("RA13OOH   ", "H2O2      ", -999.0 ) &
, dep_t("UDCARB8   ", "H2O2      ", -999.0 ) &
, dep_t("AROH14    ", "ROOH      ", -999.0 ) &
, dep_t("RA16NO3   ", "H2O2      ", -999.0 ) &
, dep_t("RA16OOH   ", "H2O2      ", -999.0 ) &
, dep_t("UDCARB11  ", "H2O2      ", -999.0 ) &
, dep_t("UDCARB14  ", "0p7       ", -999.0 ) &
, dep_t("RA19NO3   ", "H2O2      ", -999.0 ) &
, dep_t("RA19OOH   ", "H2O2      ", -999.0 ) &
, dep_t("RTN28NO3  ", "0p8       ", -999.0 ) &
, dep_t("RTN28OOH  ", "H2O2      ", -999.0 ) &
, dep_t("TNCARB26  ", "HCHO      ", -999.0 ) &
, dep_t("RTN26OOH  ", "H2O2      ", -999.0 ) &
, dep_t("NRTN28OOH ", "0p3       ", -999.0 ) &
, dep_t("RTN26PAN  ", "ROOH      ", -999.0 ) &
, dep_t("RTN25OOH  ", "1p3       ", -999.0 ) &
, dep_t("RTN24OOH  ", "H2O2      ", -999.0 ) &
, dep_t("RTN23OOH  ", "HNO3      ", -999.0 ) &
, dep_t("RTN14OOH  ", "HNO3      ", -999.0 ) &
, dep_t("RTN10OOH  ", "H2O2      ", -999.0 ) &
, dep_t("TNCARB10  ", "H2O2      ", -999.0 ) &
, dep_t("RTN25NO3  ", "H2O2      ", -999.0 ) &
, dep_t("TNCARB15  ", "H2O2      ", -999.0 ) &
, dep_t("RCOOH25   ", "H2O2      ", -999.0 ) &
, dep_t("RTX28NO3  ", "0p8       ", -999.0 ) &
, dep_t("RTX28OOH  ", "H2O2      ", -999.0 ) &
, dep_t("RTX24NO3  ", "HCHO      ", -999.0 ) &
, dep_t("RTX24OOH  ", "H2O2      ", -999.0 ) &
, dep_t("TXCARB22  ", "0p7       ", -999.0 ) &
, dep_t("RTX22NO3  ", "H2O2      ", -999.0 ) &
, dep_t("RTX22OOH  ", "H2O2      ", -999.0 ) &
, dep_t("NRTX28OOH ", "0p3       ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("HPUCARB12 ", "H2O2      ", -999.0 ) &
, dep_t("HUCARB9   ", "0p8       ", -999.0 ) &
, dep_t("IEPOX     ", "HNO3      ", -999.0 ) &
, dep_t("HMML      ", "0p5       ", -999.0 ) &
, dep_t("DHPR12OOH ", "HNO3      ", -999.0 ) &
, dep_t("DHPCARB9  ", "H2O2      ", -999.0 ) &
, dep_t("DHCARB9   ", "ROOH      ", -999.0 ) &
, dep_t("CO2C3CHO  ", "0p3       ", -999.0 ) &
, dep_t("CO2C3CO3H ", "H2O2      ", -999.0 ) &
, dep_t("RU13NO3   ", "H2O2      ", -999.0 ) &
, dep_t("RU12NO3   ", "H2O2      ", -999.0 ) &
, dep_t("RU10NO3   ", "H2O2      ", -999.0 ) &
, dep_t("CO2C3PAN  ", "0p2       ", -999.0 ) &
, dep_t("HCC7CO    ", "H2O2      ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
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
