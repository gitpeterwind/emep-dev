type(dep_t), public, dimension(NWETDEP_ADV), save :: CM_WDepMap = (/ &
!        species     surrogate      setrate
  dep_t("N2O5      ", "HNO3      ", -999.0 ) &
, dep_t("H2O2      ", "H2O2      ", -999.0 ) &
, dep_t("HONO      ", "H2O2      ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HO2NO2    ", "0p8       ", -999.0 ) &
, dep_t("HCHO      ", "ROOH      ", -999.0 ) &
, dep_t("ACETOL    ", "HCHO      ", -999.0 ) &
, dep_t("GLYOX     ", "H2O2      ", -999.0 ) &
, dep_t("MGLYOX    ", "0p3       ", -999.0 ) &
, dep_t("C5DICARB  ", "H2O2      ", -999.0 ) &
, dep_t("ETRO2H    ", "H2O2      ", -999.0 ) &
, dep_t("MEKO2H    ", "1p3       ", -999.0 ) &
, dep_t("ISRO2H    ", "H2O2      ", -999.0 ) &
, dep_t("C5DICAROOH", "H2O2      ", -999.0 ) &
, dep_t("HPALD     ", "H2O2      ", -999.0 ) &
, dep_t("MACROOH   ", "H2O2      ", -999.0 ) &
, dep_t("OXYO2H    ", "H2O2      ", -999.0 ) &
, dep_t("CH3CO3H   ", "ROOH      ", -999.0 ) &
, dep_t("PACALD    ", "H2O2      ", -999.0 ) &
, dep_t("IEPOX     ", "HNO3      ", -999.0 ) &
, dep_t("NALD      ", "ROOH      ", -999.0 ) &
, dep_t("ISON      ", "0p7       ", -999.0 ) &
, dep_t("TERPOOH   ", "H2O2      ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "PMf       ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
, dep_t("EC_f_new  ", "ECfn      ", -999.0 ) &
, dep_t("EC_f_age  ", "PMf       ", -999.0 ) &
/)
