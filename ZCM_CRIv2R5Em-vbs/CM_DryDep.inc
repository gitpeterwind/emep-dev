type(dep_t), public, dimension(NDRYDEP_ADV), save :: CM_DDepMap = (/ &
!        species     surrogate      setrate
  dep_t("NO2       ", "NO2       ", -999.0 ) &
, dep_t("O3        ", "O3        ", -999.0 ) &
, dep_t("O3_REG    ", "O3        ", -999.0 ) &
, dep_t("N2O5      ", "N2O5      ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("HONO      ", "HNO2      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HO2NO2    ", "HO2NO2    ", -999.0 ) &
, dep_t("CH3CHO    ", "ALD       ", -999.0 ) &
, dep_t("C2H5CHO   ", "MEK       ", -999.0 ) &
, dep_t("CH3CO2H   ", "CH3CO2H   ", -999.0 ) &
, dep_t("CH3COCH3  ", "MEK       ", -999.0 ) &
, dep_t("CARB14    ", "MEK       ", -999.0 ) &
, dep_t("HCHO      ", "HCHO      ", -999.0 ) &
, dep_t("HCOOH     ", "ORA       ", -999.0 ) &
, dep_t("HCOCO2H   ", "ACETOL    ", -999.0 ) &
, dep_t("HCOCO3H   ", "HCOCO3H   ", -999.0 ) &
, dep_t("MEK       ", "MEK       ", -999.0 ) &
, dep_t("MVK       ", "MEK       ", -999.0 ) &
, dep_t("MACR      ", "MEK       ", -999.0 ) &
, dep_t("GLYOX     ", "GLYOX     ", -999.0 ) &
, dep_t("MGLYOX    ", "MGLYOX    ", -999.0 ) &
, dep_t("ACETOL    ", "ACETOL    ", -999.0 ) &
, dep_t("CH3NO3    ", "CH3NO3    ", -999.0 ) &
, dep_t("C2H5NO3   ", "CH3NO3    ", -999.0 ) &
, dep_t("RN10NO3   ", "VLSOLNO3  ", -999.0 ) &
, dep_t("IC3H7NO3  ", "VLSOLNO3  ", -999.0 ) &
, dep_t("RN13NO3   ", "VLSOLNO3  ", -999.0 ) &
, dep_t("RN16NO3   ", "VLSOLNO3  ", -999.0 ) &
, dep_t("RN19NO3   ", "MDNO3OH   ", -999.0 ) &
, dep_t("HOC2H4NO3 ", "SMNO3OH   ", -999.0 ) &
, dep_t("RN9NO3    ", "MDSOLNO3  ", -999.0 ) &
, dep_t("RN12NO3   ", "MDSOLNO3  ", -999.0 ) &
, dep_t("RN15NO3   ", "MDNO3OH   ", -999.0 ) &
, dep_t("CH3OOH    ", "MEOOH     ", -999.0 ) &
, dep_t("C2H5OOH   ", "MEOOH     ", -999.0 ) &
, dep_t("RN10OOH   ", "C3H7OOH   ", -999.0 ) &
, dep_t("IC3H7OOH  ", "C3H7OOH   ", -999.0 ) &
, dep_t("RN13OOH   ", "C3H7OOH   ", -999.0 ) &
, dep_t("RN16OOH   ", "MPAN      ", -999.0 ) &
, dep_t("RN19OOH   ", "LHISOLOOH ", -999.0 ) &
, dep_t("CH3CO3H   ", "PAA       ", -999.0 ) &
, dep_t("RN8OOH    ", "HYPERACET ", -999.0 ) &
, dep_t("RN11OOH   ", "MDSOLOOH  ", -999.0 ) &
, dep_t("RN14OOH   ", "MDSOLOOH  ", -999.0 ) &
, dep_t("RN17OOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("HOC2H4OOH ", "SHISOLOOH ", -999.0 ) &
, dep_t("RN9OOH    ", "SHISOLOOH ", -999.0 ) &
, dep_t("RN12OOH   ", "RN12OOH   ", -999.0 ) &
, dep_t("RN15OOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("NRN6OOH   ", "SMNO3OH   ", -999.0 ) &
, dep_t("NRN9OOH   ", "MNO3OOH   ", -999.0 ) &
, dep_t("NRN12OOH  ", "MNO3OOH   ", -999.0 ) &
, dep_t("RU13OOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("PAN       ", "PAN       ", -999.0 ) &
, dep_t("CARB10    ", "MEK       ", -999.0 ) &
, dep_t("CARB13    ", "CH3CO2H   ", -999.0 ) &
, dep_t("CARB16    ", "MEK       ", -999.0 ) &
, dep_t("CARB9     ", "MEK       ", -999.0 ) &
, dep_t("CARB12    ", "CARB12    ", -999.0 ) &
, dep_t("C2H5CO3H  ", "C3H7OOH   ", -999.0 ) &
, dep_t("PPN       ", "PAN       ", -999.0 ) &
, dep_t("HOCH2CHO  ", "HOCH2CHO  ", -999.0 ) &
, dep_t("HOCH2CO3H ", "HOCH2CO3H ", -999.0 ) &
, dep_t("PHAN      ", "PAN       ", -999.0 ) &
, dep_t("CCARB12   ", "HISOLF0   ", -999.0 ) &
, dep_t("UCARB12   ", "DICARB    ", -999.0 ) &
, dep_t("UCARB13   ", "MEK       ", -999.0 ) &
, dep_t("RU14NO3   ", "MDNO3OH   ", -999.0 ) &
, dep_t("RU14OOH   ", "SHISOLOOH ", -999.0 ) &
, dep_t("RU12OOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("RU10OOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("MPAN      ", "MPAN      ", -999.0 ) &
, dep_t("NUCARB12  ", "LOSOLNO3  ", -999.0 ) &
, dep_t("NRU14OOH  ", "MDSOLNO3  ", -999.0 ) &
, dep_t("NRU12OOH  ", "VHISOLNO3 ", -999.0 ) &
, dep_t("NOA       ", "VLSOLNO3  ", -999.0 ) &
, dep_t("RA13NO3   ", "VHISOLNO3 ", -999.0 ) &
, dep_t("RA13OOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("UDCARB8   ", "DICARB    ", -999.0 ) &
, dep_t("AROH14    ", "PHENOL    ", -999.0 ) &
, dep_t("RA16NO3   ", "HISOLNO3  ", -999.0 ) &
, dep_t("RA16OOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("UDCARB11  ", "DICARB    ", -999.0 ) &
, dep_t("UDCARB14  ", "MCARB     ", -999.0 ) &
, dep_t("RA19NO3   ", "HISOLNO3  ", -999.0 ) &
, dep_t("RA19OOH   ", "VHISOLOOH ", -999.0 ) &
, dep_t("RTN28NO3  ", "C10H17NO4 ", -999.0 ) &
, dep_t("RTN28OOH  ", "LHISOLOOH ", -999.0 ) &
, dep_t("TNCARB26  ", "PINONALDEHYDE", -999.0 ) &
, dep_t("RTN26OOH  ", "PERPINONIC", -999.0 ) &
, dep_t("NRTN28OOH ", "C10NO3OOH ", -999.0 ) &
, dep_t("RTN26PAN  ", "C10PAN2   ", -999.0 ) &
, dep_t("RTN25OOH  ", "C96OOH    ", -999.0 ) &
, dep_t("RTN24OOH  ", "VHISOLOOH ", -999.0 ) &
, dep_t("RTN23OOH  ", "VHISOLOOH ", -999.0 ) &
, dep_t("RTN14OOH  ", "VHISOLOOH ", -999.0 ) &
, dep_t("RTN10OOH  ", "VHISOLOOH ", -999.0 ) &
, dep_t("TNCARB10  ", "CO23C4CHO ", -999.0 ) &
, dep_t("RTN25NO3  ", "VHISOLNO3 ", -999.0 ) &
, dep_t("TNCARB15  ", "HISOLF0   ", -999.0 ) &
, dep_t("RCOOH25   ", "PINONIC   ", -999.0 ) &
, dep_t("RTX28NO3  ", "C10H17NO4 ", -999.0 ) &
, dep_t("RTX28OOH  ", "LHISOLOOH ", -999.0 ) &
, dep_t("RTX24NO3  ", "LOSOLNO3  ", -999.0 ) &
, dep_t("RTX24OOH  ", "NOPINAOOH ", -999.0 ) &
, dep_t("TXCARB22  ", "MCARB     ", -999.0 ) &
, dep_t("RTX22NO3  ", "VHISOLNO3 ", -999.0 ) &
, dep_t("RTX22OOH  ", "VHISOLOOH ", -999.0 ) &
, dep_t("NRTX28OOH ", "C10NO3OOH ", -999.0 ) &
, dep_t("ANHY      ", "ANHY      ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("HPUCARB12 ", "MDSOLOOH  ", -999.0 ) &
, dep_t("HUCARB9   ", "MCARB     ", -999.0 ) &
, dep_t("IEPOX     ", "HISOLF0   ", -999.0 ) &
, dep_t("HMML      ", "CARB12    ", -999.0 ) &
, dep_t("DHPR12OOH ", "VHISOLOOH ", -999.0 ) &
, dep_t("DHPCARB9  ", "VHISOLOOH ", -999.0 ) &
, dep_t("DHCARB9   ", "MACROH    ", -999.0 ) &
, dep_t("CO2C3CHO  ", "MGLYOX    ", -999.0 ) &
, dep_t("CO2C3CO3H ", "SHISOLOOH ", -999.0 ) &
, dep_t("RU13NO3   ", "VHISOLNO3 ", -999.0 ) &
, dep_t("RU12NO3   ", "VHISOLNO3 ", -999.0 ) &
, dep_t("RU10NO3   ", "VHISOLNO3 ", -999.0 ) &
, dep_t("CO2C3PAN  ", "CO2C3PAN  ", -999.0 ) &
, dep_t("HCC7CO    ", "HCC7CO    ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH4_f     ", "PMfNH4    ", -999.0 ) &
, dep_t("NO3_f     ", "PMfNO3    ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("ASOC_ng1e2", "ALD       ", -999.0 ) &
, dep_t("ASOC_ug1  ", "ALD       ", -999.0 ) &
, dep_t("ASOC_ug10 ", "ALD       ", -999.0 ) &
, dep_t("ASOC_ug1e2", "ALD       ", -999.0 ) &
, dep_t("ASOC_ug1e3", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ng1e2", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ug1", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ug10", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ug1e2", "ALD       ", -999.0 ) &
, dep_t("non_C_ASOA_ug1e3", "ALD       ", -999.0 ) &
, dep_t("BSOC_ng1e2", "ALD       ", -999.0 ) &
, dep_t("BSOC_ug1  ", "ALD       ", -999.0 ) &
, dep_t("BSOC_ug10 ", "ALD       ", -999.0 ) &
, dep_t("BSOC_ug1e2", "ALD       ", -999.0 ) &
, dep_t("BSOC_ug1e3", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ng1e2", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ug1", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ug10", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ug1e2", "ALD       ", -999.0 ) &
, dep_t("non_C_BSOA_ug1e3", "ALD       ", -999.0 ) &
, dep_t("shipNOx   ", "NO2       ", -999.0 ) &
, dep_t("ffire_OM  ", "PMf       ", -999.0 ) &
, dep_t("ffire_BC  ", "PMf       ", -999.0 ) &
, dep_t("ffire_remPPM25", "PMf       ", -999.0 ) &
, dep_t("ffire_c   ", "PMc       ", -999.0 ) &
, dep_t("SeaSalt_f ", "SSf       ", -999.0 ) &
, dep_t("SeaSalt_c ", "SSc       ", -999.0 ) &
, dep_t("Dust_road_f", "DUf       ", -999.0 ) &
, dep_t("Dust_road_c", "DUc       ", -999.0 ) &
, dep_t("Dust_wb_f ", "DUf       ", -999.0 ) &
, dep_t("Dust_wb_c ", "DUc       ", -999.0 ) &
, dep_t("Dust_sah_f", "DUf       ", -999.0 ) &
, dep_t("Dust_sah_c", "DUc       ", -999.0 ) &
, dep_t("Ash_f     ", "PMf       ", -999.0 ) &
, dep_t("Ash_c     ", "PMc       ", -999.0 ) &
, dep_t("POM_f_wood", "PMf       ", -999.0 ) &
, dep_t("POM_c_wood", "PMc       ", -999.0 ) &
, dep_t("POM_f_ffuel", "PMf       ", -999.0 ) &
, dep_t("POM_c_ffuel", "PMc       ", -999.0 ) &
, dep_t("EC_f_wood_new", "PMf       ", -999.0 ) &
, dep_t("EC_f_wood_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_wood ", "PMc       ", -999.0 ) &
, dep_t("EC_f_ffuel_new", "PMf       ", -999.0 ) &
, dep_t("EC_f_ffuel_age", "PMf       ", -999.0 ) &
, dep_t("EC_c_ffuel", "PMc       ", -999.0 ) &
, dep_t("pSO4f     ", "PMf       ", -999.0 ) &
, dep_t("pSO4c     ", "PMc       ", -999.0 ) &
, dep_t("remPPM25  ", "PMf       ", -999.0 ) &
, dep_t("remPPM_c  ", "PMc       ", -999.0 ) &
, dep_t("SQT_SOA_NV", "PMf       ", -999.0 ) &
/)
