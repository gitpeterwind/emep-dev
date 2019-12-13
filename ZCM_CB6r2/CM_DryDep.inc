type(dep_t), public, dimension(NDRYDEP_ADV), save :: CM_DDepMap = (/ &
!        species     surrogate      setrate
  dep_t("AACD      ", "CH3CO2H   ", -999.0 ) &
, dep_t("ACET      ", "MEK       ", -999.0 ) &
, dep_t("ALD2      ", "ALD       ", -999.0 ) &
, dep_t("ALDX      ", "MEK       ", -999.0 ) &
, dep_t("EPOX      ", "HISOLF0   ", -999.0 ) &
, dep_t("FACD      ", "ORA       ", -999.0 ) &
, dep_t("FORM      ", "HCHO      ", -999.0 ) &
, dep_t("GLY       ", "GLYOX     ", -999.0 ) &
, dep_t("GLYD      ", "HOCH2CHO  ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HONO      ", "HNO2      ", -999.0 ) &
, dep_t("HPLD      ", "MDSOLOOH  ", -999.0 ) &
, dep_t("INTR      ", "MDNO3OH   ", -999.0 ) &
, dep_t("ISPD      ", "MEK       ", -999.0 ) &
, dep_t("ISPX      ", "SHISOLOOH ", -999.0 ) &
, dep_t("CH3OOH    ", "MEOOH     ", -999.0 ) &
, dep_t("MGLY      ", "MGLYOX    ", -999.0 ) &
, dep_t("N2O5      ", "N2O5      ", -999.0 ) &
, dep_t("NO2       ", "NO2       ", -999.0 ) &
, dep_t("O3        ", "O3        ", -999.0 ) &
, dep_t("OPAN      ", "PAN       ", -999.0 ) &
, dep_t("OPEN      ", "DICARB    ", -999.0 ) &
, dep_t("PACD      ", "PAA       ", -999.0 ) &
, dep_t("PAN       ", "PAN       ", -999.0 ) &
, dep_t("PANX      ", "PAN       ", -999.0 ) &
, dep_t("PNA       ", "HO2NO2    ", -999.0 ) &
, dep_t("ROOH      ", "C3H7OOH   ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("XOPN      ", "DICARB    ", -999.0 ) &
, dep_t("NTR1      ", "VLSOLNO3  ", -999.0 ) &
, dep_t("NTR2      ", "MDSOLNO3  ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH4_f     ", "PMfNH4    ", -999.0 ) &
, dep_t("NO3_f     ", "PMfNO3    ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
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
, dep_t("ASOC_ng10 ", "LVASOA    ", -999.0 ) &
, dep_t("ASOC_ng1e2", "HCC7CO    ", -999.0 ) &
, dep_t("ASOC_ug1  ", "HCC7CO    ", -999.0 ) &
, dep_t("ASOC_ug10 ", "SVASOA    ", -999.0 ) &
, dep_t("ASOC_ug1e2", "SVASOA    ", -999.0 ) &
, dep_t("ASOC_ug1e3", "SVASOA    ", -999.0 ) &
, dep_t("NON_C_ASOA_ng10", "LVASOA    ", -999.0 ) &
, dep_t("NON_C_ASOA_ng1e2", "HCC7CO    ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1", "HCC7CO    ", -999.0 ) &
, dep_t("NON_C_ASOA_ug10", "SVASOA    ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1e2", "SVASOA    ", -999.0 ) &
, dep_t("NON_C_ASOA_ug1e3", "SVASOA    ", -999.0 ) &
, dep_t("BSOC_ng10 ", "LSVBSOA   ", -999.0 ) &
, dep_t("BSOC_ng1e2", "LSVBSOA   ", -999.0 ) &
, dep_t("BSOC_ug1  ", "LSVBSOA   ", -999.0 ) &
, dep_t("BSOC_ug10 ", "LSVBSOA   ", -999.0 ) &
, dep_t("BSOC_ug1e2", "SVBSOA    ", -999.0 ) &
, dep_t("BSOC_ug1e3", "PINONIC   ", -999.0 ) &
, dep_t("NON_C_BSOA_ng10", "LSVBSOA   ", -999.0 ) &
, dep_t("NON_C_BSOA_ng1e2", "LSVBSOA   ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1", "LSVBSOA   ", -999.0 ) &
, dep_t("NON_C_BSOA_ug10", "LSVBSOA   ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1e2", "SVBSOA    ", -999.0 ) &
, dep_t("NON_C_BSOA_ug1e3", "PINONIC   ", -999.0 ) &
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
, dep_t("SQT_SOA_NV", "PMf       ", -999.0 ) &
/)
