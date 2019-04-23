type(dep_t), public, dimension(NDRYDEP_ADV), save :: CM_DDepMap = (/ &
!        species     surrogate      setrate
  dep_t("O3        ", "O3        ", -999.0 ) &
, dep_t("NO2       ", "NO2       ", -999.0 ) &
, dep_t("SHIPNOX   ", "NO2       ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("CH3CHO    ", "ALD       ", -999.0 ) &
, dep_t("PAN       ", "PAN       ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "PMfS      ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMfN      ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMfN      ", -999.0 ) &
, dep_t("HCHO      ", "HCHO      ", -999.0 ) &
, dep_t("PPM25_FIRE", "PMfS      ", -999.0 ) &
, dep_t("V1702A02B_1", "PMfS      ", -999.0 ) &
, dep_t("V1702A02B_2", "PMfS      ", -999.0 ) &
, dep_t("V1702A02B_3", "PMfS      ", -999.0 ) &
, dep_t("V1702A02B_4", "PMfS      ", -999.0 ) &
, dep_t("V1702A02B_5", "PMc       ", -999.0 ) &
, dep_t("V1702A02B_6", "PMc       ", -999.0 ) &
, dep_t("V1702A02B_7", "PMc       ", -999.0 ) &
, dep_t("SeaSalt_f ", "SSf       ", -999.0 ) &
, dep_t("SeaSalt_c ", "SSc       ", -999.0 ) &
/)
