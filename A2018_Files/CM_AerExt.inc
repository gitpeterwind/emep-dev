integer, public, parameter :: NUM_EXT = 24 - 10 - 10   !DS
type(ExtEffMap), public, dimension(NUM_EXT), parameter :: ExtMap=(/&
!DS  ExtEffMap(ASH_C           ,CEXT_DDc ),&
!DS  ExtEffMap(ASH_F           ,CEXT_DDf ),&
!DS  ExtEffMap(DUST_ROAD_C     ,CEXT_DDc ),&
!DS  ExtEffMap(DUST_ROAD_F     ,CEXT_DDf ),&
!DS  ExtEffMap(DUST_SAH_C      ,CEXT_DDc ),&
!DS  ExtEffMap(DUST_SAH_F      ,CEXT_DDf ),&
!DS  ExtEffMap(DUST_WB_C       ,CEXT_DDc ),&
!DS  ExtEffMap(DUST_WB_F       ,CEXT_DDf ),&
!DS  ExtEffMap(EC_F_FFUEL_AGE  ,CEXT_ECa ),&
!DS  ExtEffMap(EC_F_FFUEL_NEW  ,CEXT_ECn ),&
!DS  ExtEffMap(EC_F_WOOD_AGE   ,CEXT_ECa ),&
!DS  ExtEffMap(EC_F_WOOD_NEW   ,CEXT_ECn ),&
!DS  ExtEffMap(FFIRE_BC        ,CEXT_EC  ),&
!DS  ExtEffMap(FFIRE_OM        ,CEXT_OC  ),&
!DS  ExtEffMap(FFIRE_REMPPM25  ,CEXT_DDf ),&
  ExtEffMap(NH4_F           ,CEXT_NH4f),&
  ExtEffMap(NO3_C           ,CEXT_NO3c),&
  ExtEffMap(NO3_F           ,CEXT_NO3f),&
!DS  ExtEffMap(OM25_P          ,CEXT_OC  ),&
!DS  ExtEffMap(REMPPM25        ,CEXT_DDf ),&
!DS  ExtEffMap(REMPPM_C        ,CEXT_DDc ),&
!DS  ExtEffMap(SEASALT_C       ,CEXT_SSc ),&
!DS  ExtEffMap(SEASALT_F       ,CEXT_SSf ),&
  ExtEffMap(SO4             ,CEXT_SO4 )/)
