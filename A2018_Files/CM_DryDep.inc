type(dep_t), public, dimension(NDRYDEP_ADV), save :: CM_DDepMap = (/ &
!        species     surrogate      setrate
  dep_t("O3        ", "O3        ", -999.0 ) &
, dep_t("NO2       ", "NO2       ", -999.0 ) &
, dep_t("PAN       ", "PAN       ", -999.0 ) &
, dep_t("MPAN      ", "PAN       ", -999.0 ) &
, dep_t("HNO3      ", "HNO3      ", -999.0 ) &
, dep_t("HONO      ", "HNO2      ", -999.0 ) &
, dep_t("MAL       ", "ALD       ", -999.0 ) &
, dep_t("HCHO      ", "HCHO      ", -999.0 ) &
, dep_t("CH3CHO    ", "ALD       ", -999.0 ) &
, dep_t("CH3OOH    ", "ROOH      ", -999.0 ) &
, dep_t("C2H5OOH   ", "ROOH      ", -999.0 ) &
, dep_t("H2O2      ", "SO2       ", -999.0 ) &
, dep_t("ACETOL    ", "ALD       ", -999.0 ) &
, dep_t("SO2       ", "SO2       ", -999.0 ) &
, dep_t("SO4       ", "PMfS      ", -999.0 ) &
, dep_t("NH3       ", "NH3       ", -999.0 ) &
, dep_t("NO3_f     ", "PMfNO3    ", -999.0 ) &
, dep_t("NO3_c     ", "PMc       ", -999.0 ) &
, dep_t("NH4_f     ", "PMfNH4    ", -999.0 ) &
/)
