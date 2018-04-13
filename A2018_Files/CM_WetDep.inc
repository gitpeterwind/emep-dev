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
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
/)
