  integer, public, parameter ::  NDRYDEP_ADV  = 25
  type(depmap), public, dimension(NDRYDEP_ADV), parameter:: DDepMap= (/ &
       depmap( IXADV_O3, CDDEP_O3, -1) & 
      , depmap( IXADV_NO2, CDDEP_NO2, -1) & 
      , depmap( IXADV_PAN, CDDEP_PAN, -1) & 
      , depmap( IXADV_MPAN, CDDEP_PAN, -1) & 
      , depmap( IXADV_HNO3, CDDEP_HNO3, -1) & 
      , depmap( IXADV_HONO, CDDEP_HNO2, -1) & 
      , depmap( IXADV_MAL, CDDEP_ALD, -1) & 
      , depmap( IXADV_HCHO, CDDEP_HCHO, -1) & 
      , depmap( IXADV_CH3CHO, CDDEP_ALD, -1) & 
      , depmap( IXADV_CH3O2H, CDDEP_ROOH, -1) & 
      , depmap( IXADV_C2H5OOH, CDDEP_ROOH, -1) & 
      , depmap( IXADV_H2O2, CDDEP_SO2, -1) & 
      , depmap( IXADV_ACETOL, CDDEP_ALD, -1) & 
      , depmap( IXADV_SO2, CDDEP_SO2, -1) & 
      , depmap( IXADV_SO4, CDDEP_PMf, -1) & 
      , depmap( IXADV_pNO3_c, CDDEP_PMc, -1) & 
      , depmap( IXADV_NH3, CDDEP_NH3, -1) & 
      , depmap( IXADV_aNH4, CDDEP_PMf, -1) & 
      , depmap( IXADV_pNO3_f, CDDEP_PMf, -1) & 
      , depmap( IXADV_PPM25, CDDEP_PMf, -1) & 
      , depmap( IXADV_PPM25_FIRE, CDDEP_PMf, -1) & 
      , depmap( IXADV_PPMCO, CDDEP_PMc, -1) & 
      , depmap( IXADV_SeaSalt_f, CDDEP_PMf, -1) & 
      , depmap( IXADV_SeaSalt_c, CDDEP_PMc, -1) & 
      , depmap( IXADV_Pb210, CDDEP_PMf, -1) & 
  /)
