
!-> OD 

      P =   &
          rcphot(IDBO3,K) * xnew(O3 ) 

      L =  &
          rcmisc(2,k)  &
        + rcmisc(3,k)  &
        + rcmisc(4,k) 

      xnew(OD)=  max(0.0, ( xold(OD) + dt2 * P)) /(1.0 + dt2*L )

!-> OP 

      P =   &
          rcmisc(2,k) * xnew(OD )  &
        + rcmisc(3,k) * xnew(OD )  &
        + 0.3*rct(32,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.2*rct(37,k) * xnew(MVK ) * xnew(O3 )  &
        + rcphot(IDAO3,K) * xnew(O3 )  &
        + rcphot(IDNO2,K) * xnew(NO2 )  &
        + rcphot(IDNO3,K) * xnew(NO3 ) 

      L =  &
          rcmisc(1,k)  &
        + rcmisc(11,k)* xnew(NO ) 

      xnew(OP)=  max(0.0, ( xold(OP) + dt2 * P)) /(1.0 + dt2*L )

!-> OH 

      P =   &
          2.*rcmisc(4,k) * xnew(OD )  &
        + rct(4,k) * xnew(O3 ) * xnew(HO2 )  &
        + rct(6,k) * xnew(NO ) * xnew(HO2 )  &
        + rct(18,k) * xnew(CH3O2H ) * xnew(OH )  &
        + rct(23,k) * xnew(C2H5OOH ) * xnew(OH )  &
        + rct(23,k) * xnew(BURO2H ) * xnew(OH )  &
        + rct(23,k) * xnew(ETRO2H ) * xnew(OH )  &
        + 0.15*rct(31,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(23,k) * xnew(PRRO2H ) * xnew(OH )  &
        + 0.55*rct(32,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 2.0e-11 * xnew(ISRO2H ) * xnew(OH )  &
        + 0.08*rct(37,k) * xnew(MVK ) * xnew(O3 )  &
        + 2*rcphot(IDH2O2,K) * xnew(H2O2 )  &
        + rcphot(IDHNO3,K) * xnew(HNO3 )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3O2H )  &
        + rcphot(IDCH3O2H,K) * xnew(C2H5OOH )  &
        + rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MEKO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3COO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H ) 

      L =  &
          rct(3,k)* xnew(O3 )  &
        + rct(8,k)* xnew(HO2 )  &
        + rct(9,k)* xnew(H2O2 )  &
        + rct(10,k)* xnew(H2 )  &
        + rct(11,k)* xnew(HNO3 )  &
        + 2e-12*AQRCK(ICLOHSO2,K)* xnew(SO2 )  &
        + rct(12,k)* xnew(CH4 )  &
        + rcmisc(7,k)* xnew(CO )  &
        + rct(16,k)* xnew(CH3OH )  &
        + rct(18,k)* xnew(CH3O2H )  &
        + rct(19,k)* xnew(CH3O2H )  &
        + rct(20,k)* xnew(HCHO )  &
        + rct(21,k)* xnew(C2H6 )  &
        + rct(23,k)* xnew(C2H5OOH )  &
        + rct(19,k)* xnew(C2H5OOH )  &
        + rct(24,k)* xnew(CH3CHO )  &
        + rct(19,k)* xnew(CH3COO2H )  &
        + rct(28,k)* xnew(C2H5OH )  &
        + rct(29,k)* xnew(NC4H10 )  &
        + 1.15e-12* xnew(MEK )  &
        + 4.8e-12* xnew(MEKO2H )  &
        + rct(19,k)* xnew(BURO2H )  &
        + rct(23,k)* xnew(BURO2H )  &
        + rct(23,k)* xnew(ETRO2H )  &
        + rct(19,k)* xnew(ETRO2H )  &
        + rct(23,k)* xnew(PRRO2H )  &
        + rct(19,k)* xnew(PRRO2H )  &
        + 1.37e-11* xnew(OXYL )  &
        + 1.7e-11* xnew(OXYO2H )  &
        + 2.0e-11* xnew(MAL )  &
        + 2.4e-11* xnew(MALO2H )  &
        + 1.1e-11* xnew(GLYOX )  &
        + 1.5e-11* xnew(MGLYOX )  &
        + rct(33,k)* xnew(C5H8 )  &
        + rct(34,k)* xnew(MVK )  &
        + 2.0e-11* xnew(ISRO2H )  &
        + rct(35,k)* xnew(MACR )  &
        + 3.35e-11* xnew(ISNI )  &
        + 3.2e-11* xnew(CH2CO2HCH3 )  &
        + 2.0e-11* xnew(ISONO3H )  &
        + 2.2e-11* xnew(MVKO2H )  &
        + 3.7e-11* xnew(ISNIRH )  &
        + 3.7e-11* xnew(MARO2H )  &
        + rcmisc(14,k)* xnew(NO2 )  &
        + rcmisc(17,k)* xnew(C2H4 )  &
        + rcmisc(18,k)* xnew(C3H6 ) 

      xnew(OH)=  max(0.0, ( xold(OH) + dt2 * P)) /(1.0 + dt2*L )

!-> HO2 

      P =   &
          rct(3,k) * xnew(O3 ) * xnew(OH )  &
        + rct(9,k) * xnew(OH ) * xnew(H2O2 )  &
        + rct(10,k) * xnew(OH ) * xnew(H2 )  &
        + 2e-12*AQRCK(ICLOHSO2,K) * xnew(OH ) * xnew(SO2 )  &
        + rcmisc(7,k) * xnew(OH ) * xnew(CO )  &
        + rct(13,k) * xnew(CH3O2 ) * xnew(NO )  &
        + 2.*rct(14,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rct(16,k) * xnew(OH ) * xnew(CH3OH )  &
        + rct(20,k) * xnew(OH ) * xnew(HCHO )  &
        + 5.8e-16 * xnew(NO3 ) * xnew(HCHO )  &
        + 8.7e-12 * xnew(C2H5O2 ) * xnew(NO )  &
        + 0.5*1.1e-11 * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + rct(28,k) * xnew(OH ) * xnew(C2H5OH )  &
        + 0.65*rct(13,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(13,k) * xnew(ETRO2 ) * xnew(NO )  &
        + 0.12*rct(30,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.28*rct(31,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(13,k) * xnew(NO ) * xnew(PRRO2 )  &
        + rct(13,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(13,k) * xnew(MALO2 ) * xnew(NO )  &
        + 1.1e-11 * xnew(OH ) * xnew(GLYOX )  &
        + 0.06*rct(32,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.78*rct(13,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.95*rct(13,k) * xnew(MVKO2 ) * xnew(NO )  &
        + rct(13,k) * xnew(CH2CCH3 ) * xnew(NO )  &
        + 0.06*rct(37,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.05*rct(13,k) * xnew(ISNIR ) * xnew(NO )  &
        + 0.8*rct(13,k) * xnew(ISONO3 ) * xnew(NO )  &
        + 2*rcphot(IDACH2O,K) * xnew(HCHO )  &
        + rcphot(IDCH3CHO,K) * xnew(CH3CHO )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3O2H )  &
        + 0.5*rcphot(IDHCOHCO,K) * xnew(GLYOX )  &
        + rcphot(IDRCOHCO,K) * xnew(MGLYOX )  &
        + rcphot(IDCH3O2H,K) * xnew(C2H5OOH )  &
        + rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + 0.65*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(PRRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H ) 

      L =  &
          rct(4,k)* xnew(O3 )  &
        + rct(6,k)* xnew(NO )  &
        + rct(8,k)* xnew(OH )  &
        + rcmisc(5,k)* xnew(HO2 )  &
        + rcmisc(5,k)* xnew(HO2 )  &
        + rcmisc(6,k)* xnew(HO2 )  &
        + rcmisc(6,k)* xnew(HO2 )  &
        + rct(17,k)* xnew(CH3O2 )  &
        + rct(22,k)* xnew(C2H5O2 )  &
        + rct(26,k)* xnew(CH3COO2 )  &
        + rct(27,k)* xnew(CH3COO2 )  &
        + 1.0e-11* xnew(SECC4H9O2 )  &
        + 1.0e-11* xnew(MEKO2 )  &
        + 1.0e-11* xnew(ETRO2 )  &
        + 1.0e-11* xnew(PRRO2 )  &
        + 1.0e-11* xnew(OXYO2 )  &
        + 1.0e-11* xnew(MALO2 )  &
        + 1.0e-11* xnew(ISRO2 )  &
        + 1.0e-11* xnew(ISONO3 )  &
        + 1.0e-11* xnew(MVKO2 )  &
        + 1.0e-11* xnew(MACRO2 )  &
        + 1.0e-11* xnew(CH2CCH3 )  &
        + 1.0e-11* xnew(ISNIR ) 

      xnew(HO2)=  max(0.0, ( xold(HO2) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3O2 

      P =   &
          rct(12,k) * xnew(OH ) * xnew(CH4 )  &
        + rct(19,k) * xnew(CH3O2H ) * xnew(OH )  &
        + 2.0e-11 * xnew(CH3COO2 ) * xnew(NO )  &
        + 0.5*1.1e-11 * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + rct(25,k) * xnew(CH3COO2 ) * xnew(CH3COO2 )  &
        + rct(25,k) * xnew(CH3COO2 ) * xnew(CH3COO2 )  &
        + 0.31*rct(31,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rcphot(IDCH3CHO,K) * xnew(CH3CHO )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3COO2H ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + rct(14,k)* xnew(CH3O2 )  &
        + rct(14,k)* xnew(CH3O2 )  &
        + rct(15,k)* xnew(CH3O2 )  &
        + rct(15,k)* xnew(CH3O2 )  &
        + rct(17,k)* xnew(HO2 )  &
        + 1.1e-11* xnew(CH3COO2 ) 

      xnew(CH3O2)=  max(0.0, ( xold(CH3O2) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H5O2 

      P =   &
          rct(21,k) * xnew(OH ) * xnew(C2H6 )  &
        + rct(19,k) * xnew(C2H5OOH ) * xnew(OH )  &
        + 0.35*rct(13,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + 0.35*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3COX,K) * xnew(MEK ) 

      L =  &
          8.7e-12* xnew(NO )  &
        + rct(22,k)* xnew(HO2 ) 

      xnew(C2H5O2)=  max(0.0, ( xold(C2H5O2) + dt2 * P)) /(1.0 + dt2*L )

!-> SECC4H9O2 

      P =   &
          rct(29,k) * xnew(OH ) * xnew(NC4H10 )  &
        + rct(19,k) * xnew(BURO2H ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(SECC4H9O2)=  max(0.0, ( xold(SECC4H9O2) + dt2 * P)) /(1.0 + dt2*L )

!-> ISRO2 

      P =   &
          rct(33,k) * xnew(C5H8 ) * xnew(OH )  &
        + 0.12*rct(13,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 2.0e-11 * xnew(ISRO2H ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(ISRO2)=  max(0.0, ( xold(ISRO2) + dt2 * P)) /(1.0 + dt2*L )

!-> ETRO2 

      P =   &
          rct(19,k) * xnew(ETRO2H ) * xnew(OH )  &
        + rcmisc(17,k) * xnew(C2H4 ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(ETRO2)=  max(0.0, ( xold(ETRO2) + dt2 * P)) /(1.0 + dt2*L )

!-> PRRO2 

      P =   &
          rct(19,k) * xnew(PRRO2H ) * xnew(OH )  &
        + rcmisc(18,k) * xnew(OH ) * xnew(C3H6 ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(PRRO2)=  max(0.0, ( xold(PRRO2) + dt2 * P)) /(1.0 + dt2*L )

!-> OXYO2 

      P =   &
          1.37e-11 * xnew(OXYL ) * xnew(OH )  &
        + 1.7e-11 * xnew(OXYO2H ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(OXYO2)=  max(0.0, ( xold(OXYO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MEKO2 

      P =   &
          1.15e-12 * xnew(OH ) * xnew(MEK )  &
        + 4.8e-12 * xnew(MEKO2H ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(MEKO2)=  max(0.0, ( xold(MEKO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MALO2 

      P =   &
          2.0e-11 * xnew(MAL ) * xnew(OH )  &
        + 2.4e-11 * xnew(MALO2H ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(MALO2)=  max(0.0, ( xold(MALO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MVKO2 

      P =   &
          rct(34,k) * xnew(MVK ) * xnew(OH )  &
        + 2.2e-11 * xnew(MVKO2H ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(MVKO2)=  max(0.0, ( xold(MVKO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MACRO2 

      P =   &
          0.5*rct(35,k) * xnew(MACR ) * xnew(OH )  &
        + rct(36,k) * xnew(MPAN )  &
        + 3.7e-11 * xnew(MARO2H ) * xnew(OH ) 

      L =  &
          1.0e-11* xnew(NO2 )  &
        + 2.0e-11* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(MACRO2)=  max(0.0, ( xold(MACRO2) + dt2 * P)) /(1.0 + dt2*L )

!-> O3 

      P =   &
          rcmisc(1,k) * xnew(OP )  &
        + rct(27,k) * xnew(CH3COO2 ) * xnew(HO2 ) 

      L =  &
          rct(1,k)* xnew(NO )  &
        + rct(2,k)* xnew(NO2 )  &
        + rct(3,k)* xnew(OH )  &
        + rct(4,k)* xnew(HO2 )  &
        + AQRCK(ICLRC2,K)* xnew(SO2 )  &
        + rct(30,k)* xnew(C2H4 )  &
        + rct(31,k)* xnew(C3H6 )  &
        + rct(32,k)* xnew(C5H8 )  &
        + 8.0e-18* xnew(ISRO2H )  &
        + rct(37,k)* xnew(MVK )  &
        + rcphot(IDAO3,K)  &
        + rcphot(IDBO3,K) 

      xnew(O3)=  max(0.0, ( xold(O3) + dt2 * P)) /(1.0 + dt2*L )

!-> NO 

      P =   &
          rct(7,k) * xnew(NO2 ) * xnew(NO3 )  &
        + rcphot(IDNO2,K) * xnew(NO2 )  &
        + RCEMIS(QRCNO,K) 

      L =  &
          rct(1,k)* xnew(O3 )  &
        + rct(5,k)* xnew(NO3 )  &
        + rct(6,k)* xnew(HO2 )  &
        + rct(13,k)* xnew(CH3O2 )  &
        + 8.7e-12* xnew(C2H5O2 )  &
        + 2.0e-11* xnew(CH3COO2 )  &
        + rct(13,k)* xnew(SECC4H9O2 )  &
        + rct(13,k)* xnew(MEKO2 )  &
        + rct(13,k)* xnew(ETRO2 )  &
        + rct(13,k)* xnew(PRRO2 )  &
        + rct(13,k)* xnew(OXYO2 )  &
        + rct(13,k)* xnew(MALO2 )  &
        + rct(13,k)* xnew(ISRO2 )  &
        + rct(13,k)* xnew(MVKO2 )  &
        + 2.0e-11* xnew(MACRO2 )  &
        + rct(13,k)* xnew(CH2CCH3 )  &
        + rct(13,k)* xnew(ISNIR )  &
        + rct(13,k)* xnew(ISONO3 )  &
        + rcmisc(11,k)* xnew(OP ) 

      xnew(NO)=  max(0.0, ( xold(NO) + dt2 * P)) /(1.0 + dt2*L )

!-> NO2 

      P =   &
          rct(1,k) * xnew(O3 ) * xnew(NO )  &
        + rct(5,k) * xnew(NO ) * xnew(NO3 )  &
        + rct(5,k) * xnew(NO ) * xnew(NO3 )  &
        + rct(6,k) * xnew(NO ) * xnew(HO2 )  &
        + rct(7,k) * xnew(NO2 ) * xnew(NO3 )  &
        + rct(13,k) * xnew(CH3O2 ) * xnew(NO )  &
        + 8.7e-12 * xnew(C2H5O2 ) * xnew(NO )  &
        + 2.0e-11 * xnew(CH3COO2 ) * xnew(NO )  &
        + rct(13,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(13,k) * xnew(MEKO2 ) * xnew(NO )  &
        + rct(13,k) * xnew(ETRO2 ) * xnew(NO )  &
        + rct(13,k) * xnew(NO ) * xnew(PRRO2 )  &
        + rct(13,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(13,k) * xnew(MALO2 ) * xnew(NO )  &
        + 0.86*rct(13,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.95*rct(13,k) * xnew(MVKO2 ) * xnew(NO )  &
        + rct(36,k) * xnew(MPAN )  &
        + 2.0e-11 * xnew(MACRO2 ) * xnew(NO )  &
        + rct(13,k) * xnew(CH2CCH3 ) * xnew(NO )  &
        + 1.9*rct(13,k) * xnew(ISNIR ) * xnew(NO )  &
        + 1.1*rct(13,k) * xnew(ISONO3 ) * xnew(NO )  &
        + rcmisc(11,k) * xnew(OP ) * xnew(NO )  &
        + rcmisc(13,k) * xnew(N2O5 )  &
        + rcmisc(16,k) * xnew(PAN )  &
        + rcphot(IDHNO3,K) * xnew(HNO3 )  &
        + rcphot(IDNO3,K) * xnew(NO3 )  &
        + RCEMIS(QRCNO2,K) 

      L =  &
          rct(2,k)* xnew(O3 )  &
        + rct(7,k)* xnew(NO3 )  &
        + 1.0e-11* xnew(MACRO2 )  &
        + rcmisc(12,k)* xnew(NO3 )  &
        + rcmisc(14,k)* xnew(OH )  &
        + rcmisc(15,k)* xnew(CH3COO2 )  &
        + rcphot(IDNO2,K) 

      xnew(NO2)=  max(0.0, ( xold(NO2) + dt2 * P)) /(1.0 + dt2*L )

!-> PAN 

      P =   &
          rcmisc(15,k) * xnew(CH3COO2 ) * xnew(NO2 ) 

      L =  &
          rcmisc(16,k) 

      xnew(PAN)=  max(0.0, ( xold(PAN) + dt2 * P)) /(1.0 + dt2*L )

!-> MPAN 

      P =   &
          1.0e-11 * xnew(MACRO2 ) * xnew(NO2 ) 

      L =  &
          rct(36,k) 

      xnew(MPAN)=  max(0.0, ( xold(MPAN) + dt2 * P)) /(1.0 + dt2*L )

!-> NO3 

      P =   &
          rct(2,k) * xnew(O3 ) * xnew(NO2 )  &
        + rct(11,k) * xnew(OH ) * xnew(HNO3 )  &
        + rcmisc(13,k) * xnew(N2O5 ) 

      L =  &
          rct(5,k)* xnew(NO )  &
        + rct(7,k)* xnew(NO2 )  &
        + 5.8e-16* xnew(HCHO )  &
        + 7.8e-13* xnew(C5H8 )  &
        + rcmisc(12,k)* xnew(NO2 )  &
        + rcphot(IDNO3,K) 

      xnew(NO3)=  max(0.0, ( xold(NO3) + dt2 * P)) /(1.0 + dt2*L )

!-> N2O5 

      P =   &
          rcmisc(12,k) * xnew(NO2 ) * xnew(NO3 ) 

      L =  &
          rcmisc(8,k)  &
        + rcmisc(13,k) 

      xnew(N2O5)=  max(0.0, ( xold(N2O5) + dt2 * P)) /(1.0 + dt2*L )

!-> ISONO3 

      P =   &
          7.8e-13 * xnew(C5H8 ) * xnew(NO3 )  &
        + 2.0e-11 * xnew(ISONO3H ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(ISONO3)=  max(0.0, ( xold(ISONO3) + dt2 * P)) /(1.0 + dt2*L )

!-> HNO3 

      P =   &
          5.8e-16 * xnew(NO3 ) * xnew(HCHO )  &
        + 2.*rcmisc(8,k) * xnew(N2O5 )  &
        + rcmisc(14,k) * xnew(NO2 ) * xnew(OH ) 

      L =  &
          rct(11,k)* xnew(OH )  &
        + rcmisc(10,k)  &
        + rcphot(IDHNO3,K) 

      xnew(HNO3)=  max(0.0, ( xold(HNO3) + dt2 * P)) /(1.0 + dt2*L )

!-> CH2CCH3 

      P =   &
          2.0e-11 * xnew(MACRO2 ) * xnew(NO )  &
        + 3.2e-11 * xnew(CH2CO2HCH3 ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(CH2CCH3)=  max(0.0, ( xold(CH2CCH3) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3COO2 

      P =   &
          rct(24,k) * xnew(OH ) * xnew(CH3CHO )  &
        + rct(19,k) * xnew(CH3COO2H ) * xnew(OH )  &
        + rct(13,k) * xnew(MEKO2 ) * xnew(NO )  &
        + 1.5e-11 * xnew(OH ) * xnew(MGLYOX )  &
        + 0.684*rct(13,k) * xnew(MVKO2 ) * xnew(NO )  &
        + rcmisc(16,k) * xnew(PAN )  &
        + rcphot(IDRCOHCO,K) * xnew(MGLYOX )  &
        + rcphot(IDCH3O2H,K) * xnew(MEKO2H )  &
        + rcphot(IDCH3COX,K) * xnew(MEK ) 

      L =  &
          2.0e-11* xnew(NO )  &
        + 1.1e-11* xnew(CH3O2 )  &
        + rct(25,k)* xnew(CH3COO2 )  &
        + rct(25,k)* xnew(CH3COO2 )  &
        + rct(26,k)* xnew(HO2 )  &
        + rct(27,k)* xnew(HO2 )  &
        + rcmisc(15,k)* xnew(NO2 ) 

      xnew(CH3COO2)=  max(0.0, ( xold(CH3COO2) + dt2 * P)) /(1.0 + dt2*L )

!-> MACR 

      P =   &
          0.67*rct(32,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.32*rct(13,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.1*rct(13,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          rct(35,k)* xnew(OH ) 

      xnew(MACR)=  max(0.0, ( xold(MACR) + dt2 * P)) /(1.0 + dt2*L )

!-> ISNI 

      P =   &
          0.14*rct(13,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.05*rct(13,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.85*rct(13,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          3.35e-11* xnew(OH ) 

      xnew(ISNI)=  max(0.0, ( xold(ISNI) + dt2 * P)) /(1.0 + dt2*L )

!-> ISNIR 

      P =   &
          3.35e-11 * xnew(ISNI ) * xnew(OH )  &
        + 3.7e-11 * xnew(ISNIRH ) * xnew(OH ) 

      L =  &
          rct(13,k)* xnew(NO )  &
        + 1.0e-11* xnew(HO2 ) 

      xnew(ISNIR)=  max(0.0, ( xold(ISNIR) + dt2 * P)) /(1.0 + dt2*L )

!-> GLYOX 

      P =   &
          rct(13,k) * xnew(MALO2 ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H ) 

      L =  &
          1.1e-11* xnew(OH )  &
        + rcphot(IDHCOHCO,K) 

      xnew(GLYOX)=  max(0.0, ( xold(GLYOX) + dt2 * P)) /(1.0 + dt2*L )

!-> MGLYOX 

      P =   &
          rct(13,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rct(13,k) * xnew(MALO2 ) * xnew(NO )  &
        + 0.266*rct(13,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.82*rct(37,k) * xnew(MVK ) * xnew(O3 )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MALO2H ) 

      L =  &
          1.5e-11* xnew(OH )  &
        + rcphot(IDRCOHCO,K) 

      xnew(MGLYOX)=  max(0.0, ( xold(MGLYOX) + dt2 * P)) /(1.0 + dt2*L )

!-> MAL 

      P =   &
          rct(13,k) * xnew(OXYO2 ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(OXYO2H ) 

      L =  &
          2.0e-11* xnew(OH ) 

      xnew(MAL)=  max(0.0, ( xold(MAL) + dt2 * P)) /(1.0 + dt2*L )

!-> MEK 

      P =   &
          0.65*rct(13,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(23,k) * xnew(BURO2H ) * xnew(OH )  &
        + rct(23,k) * xnew(PRRO2H ) * xnew(OH )  &
        + rct(13,k) * xnew(CH2CCH3 ) * xnew(NO )  &
        + 0.95*rct(13,k) * xnew(ISNIR ) * xnew(NO )  &
        + 0.65*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + RCEMIS(QRCMEK,K) 

      L =  &
          1.15e-12* xnew(OH )  &
        + rcphot(IDCH3COX,K) 

      xnew(MEK)=  max(0.0, ( xold(MEK) + dt2 * P)) /(1.0 + dt2*L )

!-> MVK 

      P =   &
          0.26*rct(32,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.42*rct(13,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.05*rct(13,k) * xnew(ISONO3 ) * xnew(NO ) 

      L =  &
          rct(34,k)* xnew(OH )  &
        + rct(37,k)* xnew(O3 ) 

      xnew(MVK)=  max(0.0, ( xold(MVK) + dt2 * P)) /(1.0 + dt2*L )

!-> HCHO 

      P =   &
          rct(13,k) * xnew(CH3O2 ) * xnew(NO )  &
        + 2.*rct(14,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rct(15,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + rct(16,k) * xnew(OH ) * xnew(CH3OH )  &
        + rct(18,k) * xnew(CH3O2H ) * xnew(OH )  &
        + 1.1e-11 * xnew(CH3O2 ) * xnew(CH3COO2 )  &
        + 2.*rct(13,k) * xnew(ETRO2 ) * xnew(NO )  &
        + rct(30,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.5*rct(31,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(13,k) * xnew(NO ) * xnew(PRRO2 )  &
        + 0.8*rct(32,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.74*rct(13,k) * xnew(ISRO2 ) * xnew(NO )  &
        + 0.266*rct(13,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.7*8.0e-18 * xnew(ISRO2H ) * xnew(O3 )  &
        + 0.8*rct(37,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.15*rct(13,k) * xnew(ISONO3 ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(CH3O2H )  &
        + 0.1*rcphot(IDHCOHCO,K) * xnew(GLYOX )  &
        + 1.56*rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(PRRO2H )  &
        + RCEMIS(QRCHCHO,K) 

      L =  &
          rct(20,k)* xnew(OH )  &
        + 5.8e-16* xnew(NO3 )  &
        + rcphot(IDACH2O,K)  &
        + rcphot(IDBCH2O,K) 

      xnew(HCHO)=  max(0.0, ( xold(HCHO) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3CHO 

      P =   &
          8.7e-12 * xnew(C2H5O2 ) * xnew(NO )  &
        + rct(23,k) * xnew(C2H5OOH ) * xnew(OH )  &
        + rct(28,k) * xnew(OH ) * xnew(C2H5OH )  &
        + 0.35*rct(13,k) * xnew(NO ) * xnew(SECC4H9O2 )  &
        + rct(13,k) * xnew(MEKO2 ) * xnew(NO )  &
        + rct(23,k) * xnew(ETRO2H ) * xnew(OH )  &
        + 0.5*rct(31,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rct(13,k) * xnew(NO ) * xnew(PRRO2 )  &
        + 0.684*rct(13,k) * xnew(MVKO2 ) * xnew(NO )  &
        + 0.04*rct(37,k) * xnew(MVK ) * xnew(O3 )  &
        + 0.95*rct(13,k) * xnew(ISNIR ) * xnew(NO )  &
        + rcphot(IDCH3O2H,K) * xnew(C2H5OOH )  &
        + 0.22*rcphot(IDCH3O2H,K) * xnew(ETRO2H )  &
        + 0.35*rcphot(IDCH3O2H,K) * xnew(BURO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(PRRO2H )  &
        + rcphot(IDCH3O2H,K) * xnew(MEKO2H )  &
        + RCEMIS(QRCCH3CHO,K) 

      L =  &
          rct(24,k)* xnew(OH )  &
        + rcphot(IDCH3CHO,K) 

      xnew(CH3CHO)=  max(0.0, ( xold(CH3CHO) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H6 

      P =   &
          RCEMIS(QRCC2H6,K) 

      L =  &
          rct(21,k)* xnew(OH ) 

      xnew(C2H6)=  max(0.0, ( xold(C2H6) + dt2 * P)) /(1.0 + dt2*L )

!-> NC4H10 

      P =   &
          RCEMIS(QRCNC4H10,K) 

      L =  &
          rct(29,k)* xnew(OH ) 

      xnew(NC4H10)=  max(0.0, ( xold(NC4H10) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H4 

      P =   &
          RCEMIS(QRCC2H4,K) 

      L =  &
          rct(30,k)* xnew(O3 )  &
        + rcmisc(17,k)* xnew(OH ) 

      xnew(C2H4)=  max(0.0, ( xold(C2H4) + dt2 * P)) /(1.0 + dt2*L )

!-> C3H6 

      P =   &
          0.07*rct(32,k) * xnew(C5H8 ) * xnew(O3 )  &
        + RCEMIS(QRCC3H6,K) 

      L =  &
          rct(31,k)* xnew(O3 )  &
        + rcmisc(18,k)* xnew(OH ) 

      xnew(C3H6)=  max(0.0, ( xold(C3H6) + dt2 * P)) /(1.0 + dt2*L )

!-> OXYL 

      P =   &
          RCEMIS(QRCOXYL,K)  &
        + RCEMIS(QRCOXYL,K) 

      L =  &
          1.37e-11* xnew(OH ) 

      xnew(OXYL)=  max(0.0, ( xold(OXYL) + dt2 * P)) /(1.0 + dt2*L )

!-> C5H8 

      P =   &
          RCEMIS(QRCC5H8,K) 

      L =  &
          rct(32,k)* xnew(O3 )  &
        + rct(33,k)* xnew(OH )  &
        + 7.8e-13* xnew(NO3 ) 

      xnew(C5H8)=  max(0.0, ( xold(C5H8) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3O2H 

      P =   &
          rct(17,k) * xnew(HO2 ) * xnew(CH3O2 ) 

      L =  &
          rct(18,k)* xnew(OH )  &
        + rct(19,k)* xnew(OH )  &
        + 1.0e-5  &
        + rcphot(IDCH3O2H,K) 

      xnew(CH3O2H)=  max(0.0, ( xold(CH3O2H) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H5OOH 

      P =   &
          rct(22,k) * xnew(C2H5O2 ) * xnew(HO2 ) 

      L =  &
          rct(23,k)* xnew(OH )  &
        + rct(19,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(C2H5OOH)=  max(0.0, ( xold(C2H5OOH) + dt2 * P)) /(1.0 + dt2*L )

!-> BURO2H 

      P =   &
          1.0e-11 * xnew(SECC4H9O2 ) * xnew(HO2 ) 

      L =  &
          rct(19,k)* xnew(OH )  &
        + rct(23,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(BURO2H)=  max(0.0, ( xold(BURO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> ETRO2H 

      P =   &
          1.0e-11 * xnew(ETRO2 ) * xnew(HO2 ) 

      L =  &
          rct(23,k)* xnew(OH )  &
        + rct(19,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(ETRO2H)=  max(0.0, ( xold(ETRO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> PRRO2H 

      P =   &
          1.0e-11 * xnew(PRRO2 ) * xnew(HO2 ) 

      L =  &
          rct(23,k)* xnew(OH )  &
        + rct(19,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(PRRO2H)=  max(0.0, ( xold(PRRO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> OXYO2H 

      P =   &
          1.0e-11 * xnew(OXYO2 ) * xnew(HO2 ) 

      L =  &
          1.7e-11* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(OXYO2H)=  max(0.0, ( xold(OXYO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> MEKO2H 

      P =   &
          1.0e-11 * xnew(MEKO2 ) * xnew(HO2 ) 

      L =  &
          4.8e-12* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(MEKO2H)=  max(0.0, ( xold(MEKO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> MALO2H 

      P =   &
          1.0e-11 * xnew(MALO2 ) * xnew(HO2 ) 

      L =  &
          2.4e-11* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(MALO2H)=  max(0.0, ( xold(MALO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> MVKO2H 

      P =   &
          1.0e-11 * xnew(MVKO2 ) * xnew(HO2 ) 

      L =  &
          2.2e-11* xnew(OH ) 

      xnew(MVKO2H)=  max(0.0, ( xold(MVKO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> MARO2H 

      P =   &
          1.0e-11 * xnew(MACRO2 ) * xnew(HO2 ) 

      L =  &
          3.7e-11* xnew(OH ) 

      xnew(MARO2H)=  max(0.0, ( xold(MARO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> ISRO2H 

      P =   &
          1.0e-11 * xnew(ISRO2 ) * xnew(HO2 ) 

      L =  &
          2.0e-11* xnew(OH )  &
        + 8.0e-18* xnew(O3 ) 

      xnew(ISRO2H)=  max(0.0, ( xold(ISRO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> H2O2 

      P =   &
          rcmisc(5,k) * xnew(HO2 ) * xnew(HO2 )  &
        + rcmisc(6,k) * xnew(HO2 ) * xnew(HO2 ) 

      L =  &
          rct(9,k)* xnew(OH )  &
        + AQRCK(ICLRC1,K)* xnew(SO2 )  &
        + 1.0e-5  &
        + rcphot(IDH2O2,K) 

      xnew(H2O2)=  max(0.0, ( xold(H2O2) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3COO2H 

      P =   &
          rct(26,k) * xnew(CH3COO2 ) * xnew(HO2 ) 

      L =  &
          rct(19,k)* xnew(OH )  &
        + rcphot(IDCH3O2H,K) 

      xnew(CH3COO2H)=  max(0.0, ( xold(CH3COO2H) + dt2 * P)) /(1.0 + dt2*L )

!-> CH2CO2HCH3 

      P =   &
          1.0e-11 * xnew(CH2CCH3 ) * xnew(HO2 ) 

      L =  &
          3.2e-11* xnew(OH ) 

      xnew(CH2CO2HCH3)=  max(0.0, ( xold(CH2CO2HCH3) + dt2 * P)) /(1.0 + dt2*L )

!-> ISONO3H 

      P =   &
          1.0e-11 * xnew(ISONO3 ) * xnew(HO2 ) 

      L =  &
          2.0e-11* xnew(OH ) 

      xnew(ISONO3H)=  max(0.0, ( xold(ISONO3H) + dt2 * P)) /(1.0 + dt2*L )

!-> ISNIRH 

      P =   &
          1.0e-11 * xnew(ISNIR ) * xnew(HO2 ) 

      L =  &
          3.7e-11* xnew(OH ) 

      xnew(ISNIRH)=  max(0.0, ( xold(ISNIRH) + dt2 * P)) /(1.0 + dt2*L )

!-> CH3OH 

      P =   &
          rct(15,k) * xnew(CH3O2 ) * xnew(CH3O2 )  &
        + RCEMIS(QRCCH3OH,K) 

      L =  &
          rct(16,k)* xnew(OH ) 

      xnew(CH3OH)=  max(0.0, ( xold(CH3OH) + dt2 * P)) /(1.0 + dt2*L )

!-> C2H5OH 

      P =   &
          RCEMIS(QRCC2H5OH,K) 

      L =  &
          rct(28,k)* xnew(OH ) 

      xnew(C2H5OH)=  max(0.0, ( xold(C2H5OH) + dt2 * P)) /(1.0 + dt2*L )

!-> H2 

      P =   &
          0.13*rct(30,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.07*rct(31,k) * xnew(O3 ) * xnew(C3H6 )  &
        + rcphot(IDBCH2O,K) * xnew(HCHO ) 

      L =  &
          rct(10,k)* xnew(OH ) 

      xnew(H2)=  max(0.0, ( xold(H2) + dt2 * P)) /(1.0 + dt2*L )

!-> CO 

      P =   &
          rct(20,k) * xnew(OH ) * xnew(HCHO )  &
        + 5.8e-16 * xnew(NO3 ) * xnew(HCHO )  &
        + 0.44*rct(30,k) * xnew(C2H4 ) * xnew(O3 )  &
        + 0.4*rct(31,k) * xnew(O3 ) * xnew(C3H6 )  &
        + 1.1e-11 * xnew(OH ) * xnew(GLYOX )  &
        + 1.1e-11 * xnew(OH ) * xnew(GLYOX )  &
        + 1.5e-11 * xnew(OH ) * xnew(MGLYOX )  &
        + 0.05*rct(32,k) * xnew(C5H8 ) * xnew(O3 )  &
        + 0.05*rct(37,k) * xnew(MVK ) * xnew(O3 )  &
        + rcphot(IDACH2O,K) * xnew(HCHO )  &
        + rcphot(IDBCH2O,K) * xnew(HCHO )  &
        + rcphot(IDCH3CHO,K) * xnew(CH3CHO )  &
        + 1.9*rcphot(IDHCOHCO,K) * xnew(GLYOX )  &
        + rcphot(IDRCOHCO,K) * xnew(MGLYOX )  &
        + RCEMIS(QRCCO,K) 

      L =  &
          rcmisc(7,k)* xnew(OH ) 

      xnew(CO)=  max(0.0, ( xold(CO) + dt2 * P)) /(1.0 + dt2*L )

!-> CH4 

      P =   &
          0.07*rct(31,k) * xnew(O3 ) * xnew(C3H6 ) 

      L =  &
          rct(12,k)* xnew(OH ) 

      xnew(CH4)=  max(0.0, ( xold(CH4) + dt2 * P)) /(1.0 + dt2*L )

!-> SO2 

      P =   &
          (1-0.05)*RCEMIS(QRCSO2,K) 

      L =  &
          2e-12*AQRCK(ICLOHSO2,K)* xnew(OH )  &
        + AQRCK(ICLRC1,K)* xnew(H2O2 )  &
        + AQRCK(ICLRC2,K)* xnew(O3 )  &
        + AQRCK(ICLRC3,K) 

      xnew(SO2)=  max(0.0, ( xold(SO2) + dt2 * P)) /(1.0 + dt2*L )
