type(dep_t), public, dimension(NWETDEP_ADV), save :: CM_WDepMap = (/ &
!        species     surrogate      setrate
  dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "SO4       ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMf       ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMf       ", -999.0 ) &
, dep_t("HCHO      ", "HCHO      ", -999.0 ) &
, dep_t("PPM25_FIRE", "PMf       ", -999.0 ) &
, dep_t("NPP_Cs137 ", "PMf       ", -999.0 ) &
, dep_t("NPP_I131  ", "PMf       ", -999.0 ) &
, dep_t("NPP_Sr90  ", "PMf       ", -999.0 ) &
, dep_t("ASH_1     ", "PMf       ", -999.0 ) &
, dep_t("ASH_2     ", "PMf       ", -999.0 ) &
, dep_t("ASH_3     ", "PMf       ", -999.0 ) &
, dep_t("ASH_4     ", "PMf       ", -999.0 ) &
, dep_t("ASH_5     ", "PMc       ", -999.0 ) &
, dep_t("ASH_6     ", "PMc       ", -999.0 ) &
, dep_t("ASH_7     ", "PMc       ", -999.0 ) &
/)
