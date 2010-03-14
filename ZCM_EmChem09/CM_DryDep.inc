  integer, public, parameter ::  NDRYDEP_ADV  = 24
  type(depmap), public, dimension(NDRYDEP_ADV), parameter:: Dep= (/ &
       depmap( IXADV_HNO3, CDEP_HNO3, -1) & 
      , depmap( IXADV_HONO, CDEP_HNO2, -1) & 
      , depmap( IXADV_PAN, CDEP_PAN, -1) & 
      , depmap( IXADV_NO2, CDEP_NO2, -1) & 
      , depmap( IXADV_SO2, CDEP_SO2, -1) & 
      , depmap( IXADV_SO4, CDEP_FIN, -1) & 
      , depmap( IXADV_NH3, CDEP_NH3, -1) & 
      , depmap( IXADV_ANH4, CDEP_FIN, -1) & 
      , depmap( IXADV_ANO3, CDEP_FIN, -1) & 
      , depmap( IXADV_O3, CDEP_O3, -1) & 
      , depmap( IXADV_H2O2, CDEP_SO2, -1) & 
      , depmap( IXADV_MPAN, CDEP_PAN, -1) & 
      , depmap( IXADV_HCHO, CDEP_HCHO, -1) & 
      , depmap( IXADV_CH3CHO, CDEP_ALD, -1) & 
      , depmap( IXADV_CH3O2H, CDEP_OP, -1) & 
      , depmap( IXADV_C2H5OOH, CDEP_OP, -1) & 
      , depmap( IXADV_ACETOL, CDEP_ALD, -1) & 
      , depmap( IXADV_PNO3, CDEP_COA, -1) & 
      , depmap( IXADV_MAL, CDEP_ALD, -1) & 
      , depmap( IXADV_PPM25, CDEP_FIN, -1) & 
      , depmap( IXADV_PPMCO, CDEP_COA, -1) & 
      , depmap( IXADV_SSFI, CDEP_FIN, -1) & 
      , depmap( IXADV_SSCO, CDEP_COA, -1) & 
      , depmap( IXADV_PB210, CDEP_FIN, -1) & 
  /)
