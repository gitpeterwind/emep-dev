  integer, public, parameter ::  NWETDEP_ADV  = 16
  type(depmap), public, dimension(NWETDEP_ADV), parameter:: WDepMap= (/ &
       depmap( HNO3, CWDEP_HNO3, -1) & 
      , depmap( HONO, CWDEP_HNO3, -1) & 
      , depmap( HCHO, CWDEP_HCHO, -1) & 
      , depmap( H2O2, CWDEP_H2O2, -1) & 
      , depmap( SO2, CWDEP_SO2, -1) & 
      , depmap( SO4, CWDEP_SO4, -1) & 
      , depmap( pNO3_c, CWDEP_PMc, -1) & 
      , depmap( NH3, CWDEP_NH3, -1) & 
      , depmap( aNH4, CWDEP_PMf, -1) & 
      , depmap( pNO3_f, CWDEP_PMf, -1) & 
      , depmap( PPM25, CWDEP_PMf, -1) & 
      , depmap( PPM25_FIRE, CWDEP_PMf, -1) & 
      , depmap( PPMCO, CWDEP_PMc, -1) & 
      , depmap( SeaSalt_f, CWDEP_PMf, -1) & 
      , depmap( SeaSalt_c, CWDEP_PMc, -1) & 
      , depmap( Pb210, CWDEP_PMf, -1) & 
  /)
