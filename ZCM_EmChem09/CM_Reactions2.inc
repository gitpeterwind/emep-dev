
!-> SO4 

      P =   &
          2e-12*AQRCK(ICLOHSO2,K) * xnew(OH ) * xnew(SO2 )  &
        + AQRCK(ICLRC1,K) * xnew(SO2 ) * xnew(H2O2 )  &
        + AQRCK(ICLRC2,K) * xnew(SO2 ) * xnew(O3 )  &
        + AQRCK(ICLRC3,K) * xnew(SO2 )  &
        + RCEMIS(SO4,K) 
      ! L = 0.0
 

      xnew(SO4)=  max(0.0,  xold(SO4) + dt2 * P)

!-> PNO3_C 

      P =   &
          rcmisc(19,k) * xnew(HNO3 ) 
      ! L = 0.0
 

      xnew(PNO3_C)=  max(0.0,  xold(PNO3_C) + dt2 * P)

!-> NH3 

      P =   &
          RCEMIS(NH3,K)  &
        + RCNH3(K) 
      ! L = 0.0
 

      xnew(NH3)=  max(0.0,  xold(NH3) + dt2 * P)

!-> ANH4 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for ANH4! xnew(ANH4)=  max(0.0, xold(ANH4))

!-> PNO3_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for PNO3_F! xnew(PNO3_F)=  max(0.0, xold(PNO3_F))

!-> PPM25 

      P =   &
          RCEMIS(PPM25,K) 
      ! L = 0.0
 

      xnew(PPM25)=  max(0.0,  xold(PPM25) + dt2 * P)

!-> PPM25_FIRE 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for PPM25_FIRE! xnew(PPM25_FIRE)=  max(0.0, xold(PPM25_FIRE))

!-> PPMCO 

      P =   &
          RCEMIS(PPMCO,K) 
      ! L = 0.0
 

      xnew(PPMCO)=  max(0.0,  xold(PPMCO) + dt2 * P)

!-> SEASALT_F 

      P =   &
          RCSS(QSSFI,K) 
      ! L = 0.0
 

      xnew(SEASALT_F)=  max(0.0,  xold(SEASALT_F) + dt2 * P)

!-> SEASALT_C 

      P =   &
          RCSS(QSSCO,K) 
      ! L = 0.0
 

      xnew(SEASALT_C)=  max(0.0,  xold(SEASALT_C) + dt2 * P)

!-> RN222 

      P =   &
          RC_RN222(K) 

      L =  &
          2.1e-6 

      xnew(RN222)=  max(0.0, ( xold(RN222) + dt2 * P)) /(1.0 + dt2*L )

!-> PB210 

      P =   &
          2.1e-6 * xnew(RN222 ) 

      L =  &
          rcmisc(18,k) 

      xnew(PB210)=  max(0.0, ( xold(PB210) + dt2 * P)) /(1.0 + dt2*L )
