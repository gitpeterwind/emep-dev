
!-> OD 

      P =   &
          rcphot(IDBO3,k) * xnew(O3 ) 

      L =  &
          rct(2,k) 

      xnew(OD)=  ( xold(OD) + dt2 * P) /(1.0 + dt2*L )

!-> OP 

      P =   &
          rcphot(IDNO2,k) * xnew(NO2 ) 

      L =  &
          rct(1,k) 

      xnew(OP)=  ( xold(OP) + dt2 * P) /(1.0 + dt2*L )

!-> OH 

      P =   &
          2.*rct(2,k) * xnew(OD )  &
        + rct(6,k) * xnew(NO ) * xnew(HO2 ) 

      L =  &
          rct(4,k)* xnew(VOC )  &
        + rct(5,k)* xnew(CH3CHO )  &
        + 1.11e-11* xnew(NO2 ) 

      xnew(OH)=  ( xold(OH) + dt2 * P) /(1.0 + dt2*L )

!-> HO2 

      P =   &
          rct(7,k) * xnew(NO ) * xnew(RO2 )  &
        + rcphot(IDCH3CHO,k) * xnew(CH3CHO ) 

      L =  &
          rct(6,k)* xnew(NO ) 

      xnew(HO2)=  ( xold(HO2) + dt2 * P) /(1.0 + dt2*L )

!-> RO2 

      P =   &
          rct(4,k) * xnew(OH ) * xnew(VOC )  &
        + rct(8,k) * xnew(CH3COO2 ) * xnew(NO )  &
        + rcphot(IDCH3CHO,k) * xnew(CH3CHO ) 

      L =  &
          rct(7,k)* xnew(NO ) 

      xnew(RO2)=  ( xold(RO2) + dt2 * P) /(1.0 + dt2*L )

!-> TRACER 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for TRACER! xnew(TRACER)=  max(0.0, xold(TRACER))

!-> O3 

      P =   &
          rct(1,k) * xnew(OP ) 

      L =  &
          rct(3,k)* xnew(NO )  &
        + rcphot(IDBO3,k) 

      xnew(O3)=  ( xold(O3) + dt2 * P) /(1.0 + dt2*L )

!-> NO 

      P =   &
          rcphot(IDNO2,k) * xnew(NO2 )  &
        + rcemis(NO,k) 

      L =  &
          rct(3,k)* xnew(O3 )  &
        + rct(6,k)* xnew(HO2 )  &
        + rct(7,k)* xnew(RO2 )  &
        + rct(8,k)* xnew(CH3COO2 ) 

      xnew(NO)=  ( xold(NO) + dt2 * P) /(1.0 + dt2*L )

!-> NO2 

      P =   &
          rct(3,k) * xnew(O3 ) * xnew(NO )  &
        + rct(6,k) * xnew(NO ) * xnew(HO2 )  &
        + rct(7,k) * xnew(NO ) * xnew(RO2 )  &
        + rct(8,k) * xnew(CH3COO2 ) * xnew(NO )  &
        + rct(9,k) * xnew(PAN ) 

      L =  &
          1.11e-11* xnew(OH )  &
        + 4.7e-12* xnew(CH3COO2 )  &
        + rcphot(IDNO2,k) 

      xnew(NO2)=  ( xold(NO2) + dt2 * P) /(1.0 + dt2*L )

!-> HNO3 

      P =   &
          1.11e-11 * xnew(NO2 ) * xnew(OH ) 
      ! L = 0.0
 

      xnew(HNO3)=  xold(HNO3) + dt2 * P

!-> VOC 

      P =   &
          rcemis(VOC,k) 

      L =  &
          rct(4,k)* xnew(OH ) 

      xnew(VOC)=  ( xold(VOC) + dt2 * P) /(1.0 + dt2*L )

!-> CH3CHO 

      P =   &
          rct(7,k) * xnew(NO ) * xnew(RO2 ) 

      L =  &
          rct(5,k)* xnew(OH )  &
        + rcphot(IDCH3CHO,k) 

      xnew(CH3CHO)=  ( xold(CH3CHO) + dt2 * P) /(1.0 + dt2*L )

!-> PAN 

      P =   &
          4.7e-12 * xnew(CH3COO2 ) * xnew(NO2 ) 

      L =  &
          rct(9,k) 

      xnew(PAN)=  ( xold(PAN) + dt2 * P) /(1.0 + dt2*L )

!-> CH3COO2 

      P =   &
          rct(5,k) * xnew(OH ) * xnew(CH3CHO )  &
        + rct(9,k) * xnew(PAN ) 

      L =  &
          rct(8,k)* xnew(NO )  &
        + 4.7e-12* xnew(NO2 ) 

      xnew(CH3COO2)=  ( xold(CH3COO2) + dt2 * P) /(1.0 + dt2*L )

!-> SO2 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for SO2! xnew(SO2)=  max(0.0, xold(SO2))

!-> SO4 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for SO4! xnew(SO4)=  max(0.0, xold(SO4))

!-> NH3 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NH3! xnew(NH3)=  max(0.0, xold(NH3))

!-> NO3_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NO3_F! xnew(NO3_F)=  max(0.0, xold(NO3_F))

!-> NH4_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NH4_F! xnew(NH4_F)=  max(0.0, xold(NH4_F))
