!-> OD

  P = rcphot(IDO3_O1D,k) * xnew(O3)
  
  L = rct(2,k)  &
    + rct(3,k)  &
    + rct(4,k)
  
  xnew(OD) = (xold(OD) + dt2 * P) / (1.0 + dt2 * L)


!-> OP

  P = rct(2,k) * xnew(OD)  &
    + rct(3,k) * xnew(OD)  &
    + 0.3* rct(51,k) * xnew(C5H8) * xnew(O3)  &
    + 0.2* rct(60,k) * xnew(MVK) * xnew(O3)  &
    + rcphot(IDO3_O3P,k) * xnew(O3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + rcphot(IDNO3_NO2,k) * xnew(NO3)
  
  L = rct(1,k)  &
    + rct(63,k) * xnew(NO)
  
  xnew(OP) = (xold(OP) + dt2 * P) / (1.0 + dt2 * L)


!-> OH

  P = 2.* rct(4,k) * xnew(OD)  &
    + rct(8,k) * xnew(O3) * xnew(HO2)  &
    + rct(10,k) * xnew(NO) * xnew(HO2)  &
    + 0.4* rct(26,k) * xnew(CH3OOH) * xnew(OH)  &
    + 8.01e-12 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.44* rct(37,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 2.15e-11 * xnew(BURO2H) * xnew(OH)  &
    + 1.38e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.13* rct(44,k) * xnew(C2H4) * xnew(O3)  &
    + 0.36* rct(45,k) * xnew(O3) * xnew(C3H6)  &
    + 2.44e-11 * xnew(PRRO2H) * xnew(OH)  &
    + 0.55* rct(51,k) * xnew(C5H8) * xnew(O3)  &
    + 7.5e-11 * xnew(ISRO2H) * xnew(OH)  &
    + 0.82* rct(54,k) * xnew(MACR) * xnew(O3)  &
    + 0.08* rct(60,k) * xnew(MVK) * xnew(O3)  &
    + 2* rcphot(IDH2O2,k) * xnew(H2O2)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)
  
  !-> OH cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MALO2H)
  
  L = rct(7,k) * xnew(O3)  &
    + rct(12,k) * xnew(HO2)  &
    + rct(13,k) * xnew(H2O2)  &
    + rct(14,k) * xnew(H2)  &
    + rct(15,k) * xnew(HNO3)  &
    + rct(18,k) * xnew(HONO)  &
    + rct(19,k) * xnew(CH4)  &
    + 1.44e-13 * xnew(CO)  &
    + rct(20,k) * xnew(CO)  &
    + rct(24,k) * xnew(CH3OH)  &
    + rct(26,k) * xnew(CH3OOH)  &
    + rct(27,k) * xnew(HCHO)  &
    + rct(29,k) * xnew(C2H6)  &
    + 8.01e-12 * xnew(C2H5OOH)  &
    + rct(32,k) * xnew(C2H5OOH)  &
    + rct(33,k) * xnew(CH3CHO)  &
    + rct(38,k) * xnew(CH3CO3H)  &
    + rct(39,k) * xnew(C2H5OH)  &
    + rct(40,k) * xnew(NC4H10)  &
    + rct(43,k) * xnew(MEK)
  
  !-> OH cont.
    L = L +   rct(38,k) * xnew(MEKO2H)  &
    + rct(32,k) * xnew(BURO2H)  &
    + 2.15e-11 * xnew(BURO2H)  &
    + 1.38e-11 * xnew(ETRO2H)  &
    + rct(32,k) * xnew(ETRO2H)  &
    + 2.44e-11 * xnew(PRRO2H)  &
    + rct(32,k) * xnew(PRRO2H)  &
    + 1.36e-11 * xnew(OXYL)  &
    + 4.2e-11 * xnew(OXYO2H)  &
    + 5.58e-11 * xnew(MAL)  &
    + rct(38,k) * xnew(MALO2H)  &
    + rct(49,k) * xnew(GLYOX)  &
    + rct(50,k) * xnew(MGLYOX)  &
    + rct(52,k) * xnew(C5H8)  &
    + rct(53,k) * xnew(MVK)  &
    + 7.5e-11 * xnew(ISRO2H)  &
    + rct(55,k) * xnew(MACR)  &
    + 2.82e-11 * xnew(MACROOH)  &
    + rct(58,k) * xnew(ACETOL)  &
    + 5.96e-11 * xnew(ISNI)
  
  !-> OH cont.
    L = L +   1.87e-11 * xnew(MACO3H)  &
    + 1.51e-11 * xnew(MACO2H)  &
    + 2.0e-11 * xnew(ISONO3H)  &
    + 2.2e-11 * xnew(MVKO2H)  &
    + 3.7e-11 * xnew(ISNIRH)  &
    + 2.9e-11 * xnew(MPAN)  &
    + rct(66,k) * xnew(NO2)  &
    + rct(69,k) * xnew(C2H4)  &
    + rct(70,k) * xnew(C3H6)  &
    + rct(72,k) * xnew(NO)
  
  xnew(OH) = (xold(OH) + dt2 * P) / (1.0 + dt2 * L)


!-> RO2POOL

  ! P = 0.0
  
  ! L = 0.0
  
  xnew(RO2POOL) = sum(xnew(RO2_GROUP))


!-> HO2

  P = rct(7,k) * xnew(O3) * xnew(OH)  &
    + rct(13,k) * xnew(OH) * xnew(H2O2)  &
    + rct(14,k) * xnew(OH) * xnew(H2)  &
    + 1.44e-13 * xnew(OH) * xnew(CO)  &
    + rct(20,k) * xnew(OH) * xnew(CO)  &
    + rct(21,k) * xnew(CH3O2) * xnew(NO)  &
    + 2.* rct(22,k) * xnew(CH3O2) * xnew(CH3O2)  &
    + rct(24,k) * xnew(OH) * xnew(CH3OH)  &
    + rct(27,k) * xnew(OH) * xnew(HCHO)  &
    + rct(28,k) * xnew(NO3) * xnew(HCHO)  &
    + rct(30,k) * xnew(C2H5O2) * xnew(NO)  &
    + 0.9* rct(35,k) * xnew(CH3O2) * xnew(CH3CO3)  &
    + rct(39,k) * xnew(OH) * xnew(C2H5OH)  &
    + 0.65* rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + rct(41,k) * xnew(ETRO2) * xnew(NO)  &
    + 0.13* rct(44,k) * xnew(C2H4) * xnew(O3)  &
    + 0.28* rct(45,k) * xnew(O3) * xnew(C3H6)  &
    + rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + rct(41,k) * xnew(MALO2) * xnew(NO)
  
  !-> HO2 cont.
    P = P +   rct(49,k) * xnew(OH) * xnew(GLYOX)  &
    + 0.06* rct(51,k) * xnew(C5H8) * xnew(O3)  &
    + 0.78* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.95* rct(41,k) * xnew(MVKO2) * xnew(NO)  &
    + 0.41* rct(54,k) * xnew(MACR) * xnew(O3)  &
    + 0.95* rct(56,k) * xnew(MACRO2) * xnew(NO)  &
    + rct(57,k) * xnew(MACRO2) * xnew(NO3)  &
    + rct(58,k) * xnew(ACETOL) * xnew(OH)  &
    + 0.06* rct(60,k) * xnew(MVK) * xnew(O3)  &
    + 0.05* rct(41,k) * xnew(ISNIR) * xnew(NO)  &
    + 0.8* rct(41,k) * xnew(ISONO3) * xnew(NO)  &
    + 2* rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + 2* rcphot(IDCHOCHO_2CHO,k) * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + 0.65* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)
  
  !-> HO2 cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MALO2H)
  
  L = rct(8,k) * xnew(O3)  &
    + rct(10,k) * xnew(NO)  &
    + rct(12,k) * xnew(OH)  &
    + rct(16,k) * xnew(HO2)  &
    + rct(16,k) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2)  &
    + rct(25,k) * xnew(CH3O2)  &
    + rct(31,k) * xnew(C2H5O2)  &
    + rct(37,k) * xnew(CH3CO3)  &
    + rct(42,k) * xnew(SECC4H9O2)  &
    + rct(42,k) * xnew(MEKO2)  &
    + 1.2e-11 * xnew(ETRO2)  &
    + rct(46,k) * xnew(PRRO2)  &
    + rct(47,k) * xnew(OXYO2)  &
    + rct(48,k) * xnew(MALO2)  &
    + rct(48,k) * xnew(ISRO2)  &
    + rct(42,k) * xnew(MACRO2)  &
    + rct(48,k) * xnew(ISONO3)  &
    + rct(42,k) * xnew(MVKO2)
  
  !-> HO2 cont.
    L = L +   rct(62,k) * xnew(MACO3)  &
    + rct(48,k) * xnew(ISNIR)
  
  xnew(HO2) = (xold(HO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3O2

  P = rct(19,k) * xnew(OH) * xnew(CH4)  &
    + 0.6* rct(26,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.05* rct(33,k) * xnew(OH) * xnew(CH3CHO)  &
    + rct(34,k) * xnew(CH3CO3) * xnew(NO)  &
    + 0.9* rct(35,k) * xnew(CH3O2) * xnew(CH3CO3)  &
    + rct(36,k) * xnew(CH3CO3) * xnew(CH3CO3)  &
    + rct(36,k) * xnew(CH3CO3) * xnew(CH3CO3)  &
    + 0.44* rct(37,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.28* rct(45,k) * xnew(O3) * xnew(C3H6)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)
  
  L = rct(21,k) * xnew(NO)  &
    + rct(22,k) * xnew(CH3O2)  &
    + rct(22,k) * xnew(CH3O2)  &
    + rct(23,k) * xnew(CH3O2)  &
    + rct(23,k) * xnew(CH3O2)  &
    + rct(25,k) * xnew(HO2)  &
    + rct(35,k) * xnew(CH3CO3)
  
  xnew(CH3O2) = (xold(CH3O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5O2

  P = rct(29,k) * xnew(OH) * xnew(C2H6)  &
    + rct(32,k) * xnew(C2H5OOH) * xnew(OH)  &
    + 0.35* rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + 0.35* rcphot(IDCH3O2H,k) * xnew(BURO2H)
  
  L = rct(30,k) * xnew(NO)  &
    + rct(31,k) * xnew(HO2)
  
  xnew(C2H5O2) = (xold(C2H5O2) + dt2 * P) / (1.0 + dt2 * L)


!-> SECC4H9O2

  P = rct(40,k) * xnew(OH) * xnew(NC4H10)  &
    + rct(32,k) * xnew(BURO2H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(42,k) * xnew(HO2)
  
  xnew(SECC4H9O2) = (xold(SECC4H9O2) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2

  P = rct(52,k) * xnew(C5H8) * xnew(OH)  &
    + 0.12* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 7.5e-11 * xnew(ISRO2H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(48,k) * xnew(HO2)
  
  xnew(ISRO2) = (xold(ISRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2

  P = rct(32,k) * xnew(ETRO2H) * xnew(OH)  &
    + rct(69,k) * xnew(C2H4) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + 1.2e-11 * xnew(HO2)
  
  xnew(ETRO2) = (xold(ETRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> PRRO2

  P = rct(32,k) * xnew(PRRO2H) * xnew(OH)  &
    + rct(70,k) * xnew(OH) * xnew(C3H6)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(46,k) * xnew(HO2)
  
  xnew(PRRO2) = (xold(PRRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYO2

  P = 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + 4.2e-11 * xnew(OXYO2H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(47,k) * xnew(HO2)
  
  xnew(OXYO2) = (xold(OXYO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2

  P = rct(43,k) * xnew(OH) * xnew(MEK)  &
    + rct(38,k) * xnew(MEKO2H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(42,k) * xnew(HO2)
  
  xnew(MEKO2) = (xold(MEKO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MALO2

  P = 5.58e-11 * xnew(MAL) * xnew(OH)  &
    + rct(38,k) * xnew(MALO2H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(48,k) * xnew(HO2)
  
  xnew(MALO2) = (xold(MALO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MVKO2

  P = rct(53,k) * xnew(MVK) * xnew(OH)  &
    + 2.2e-11 * xnew(MVKO2H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(42,k) * xnew(HO2)
  
  xnew(MVKO2) = (xold(MVKO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACRO2

  P = 0.5* rct(55,k) * xnew(MACR) * xnew(OH)  &
    + 2.82e-11 * xnew(MACROOH) * xnew(OH)
  
  L = rct(56,k) * xnew(NO)  &
    + rct(57,k) * xnew(NO3)  &
    + rct(42,k) * xnew(HO2)
  
  xnew(MACRO2) = (xold(MACRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACO3

  P = 0.5* rct(55,k) * xnew(MACR) * xnew(OH)  &
    + 1.87e-11 * xnew(MACO3H) * xnew(OH)  &
    + rct(71,k) * xnew(MPAN)
  
  L = rct(59,k) * xnew(NO)  &
    + rct(62,k) * xnew(HO2)  &
    + rct(67,k) * xnew(NO2)
  
  xnew(MACO3) = (xold(MACO3) + dt2 * P) / (1.0 + dt2 * L)


!-> NMVOC

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NMVOC! xnew(NMVOC) = max(0.0, xold(NMVOC))


!-> O3

  P = rct(1,k) * xnew(OP)  &
    + 0.15* rct(37,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.29* rct(62,k) * xnew(MACO3) * xnew(HO2)
  
  L = rct(5,k) * xnew(NO)  &
    + rct(6,k) * xnew(NO2)  &
    + rct(7,k) * xnew(OH)  &
    + rct(8,k) * xnew(HO2)  &
    + rct(44,k) * xnew(C2H4)  &
    + rct(45,k) * xnew(C3H6)  &
    + rct(51,k) * xnew(C5H8)  &
    + rct(54,k) * xnew(MACR)  &
    + rct(60,k) * xnew(MVK)  &
    + rcphot(IDO3_O1D,k)  &
    + rcphot(IDO3_O3P,k)
  
  xnew(O3) = (xold(O3) + dt2 * P) / (1.0 + dt2 * L)


!-> NO

  P = rct(11,k) * xnew(NO2) * xnew(NO3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + rcphot(IDNO3_NO,k) * xnew(NO3)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcemis(NO,k)
  
  L = rct(5,k) * xnew(O3)  &
    + rct(9,k) * xnew(NO3)  &
    + rct(10,k) * xnew(HO2)  &
    + rct(21,k) * xnew(CH3O2)  &
    + rct(30,k) * xnew(C2H5O2)  &
    + rct(34,k) * xnew(CH3CO3)  &
    + rct(41,k) * xnew(SECC4H9O2)  &
    + rct(41,k) * xnew(MEKO2)  &
    + rct(41,k) * xnew(ETRO2)  &
    + rct(41,k) * xnew(PRRO2)  &
    + rct(41,k) * xnew(OXYO2)  &
    + rct(41,k) * xnew(MALO2)  &
    + rct(41,k) * xnew(ISRO2)  &
    + rct(41,k) * xnew(MVKO2)  &
    + rct(56,k) * xnew(MACRO2)  &
    + rct(59,k) * xnew(MACO3)  &
    + rct(41,k) * xnew(ISNIR)  &
    + rct(41,k) * xnew(ISONO3)  &
    + rct(63,k) * xnew(OP)  &
    + rct(72,k) * xnew(OH)
  
  xnew(NO) = (xold(NO) + dt2 * P) / (1.0 + dt2 * L)


!-> NO2

  P = rct(5,k) * xnew(O3) * xnew(NO)  &
    + rct(9,k) * xnew(NO) * xnew(NO3)  &
    + rct(9,k) * xnew(NO) * xnew(NO3)  &
    + rct(10,k) * xnew(NO) * xnew(HO2)  &
    + rct(11,k) * xnew(NO2) * xnew(NO3)  &
    + rct(18,k) * xnew(OH) * xnew(HONO)  &
    + rct(21,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(30,k) * xnew(C2H5O2) * xnew(NO)  &
    + rct(34,k) * xnew(CH3CO3) * xnew(NO)  &
    + rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + rct(41,k) * xnew(MEKO2) * xnew(NO)  &
    + rct(41,k) * xnew(ETRO2) * xnew(NO)  &
    + rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + rct(41,k) * xnew(MALO2) * xnew(NO)  &
    + 0.86* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.95* rct(41,k) * xnew(MVKO2) * xnew(NO)  &
    + 0.95* rct(56,k) * xnew(MACRO2) * xnew(NO)  &
    + rct(57,k) * xnew(MACRO2) * xnew(NO3)  &
    + rct(59,k) * xnew(MACO3) * xnew(NO)
  
  !-> NO2 cont.
    P = P +   1.9* rct(41,k) * xnew(ISNIR) * xnew(NO)  &
    + 1.1* rct(41,k) * xnew(ISONO3) * xnew(NO)  &
    + 2.9e-11 * xnew(MPAN) * xnew(OH)  &
    + rct(63,k) * xnew(OP) * xnew(NO)  &
    + rct(65,k) * xnew(N2O5)  &
    + rct(68,k) * xnew(PAN)  &
    + rct(71,k) * xnew(MPAN)  &
    + rcphot(IDNO3_NO2,k) * xnew(NO3)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + rcemis(NO2,k)
  
  L = rct(6,k) * xnew(O3)  &
    + rct(11,k) * xnew(NO3)  &
    + rct(64,k) * xnew(NO3)  &
    + rct(66,k) * xnew(OH)  &
    + rct(67,k) * xnew(CH3CO3)  &
    + rct(67,k) * xnew(MACO3)  &
    + rcphot(IDNO2,k)
  
  xnew(NO2) = (xold(NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> PAN

  P = rct(67,k) * xnew(CH3CO3) * xnew(NO2)
  
  L = rct(68,k)
  
  xnew(PAN) = (xold(PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> MPAN

  P = rct(67,k) * xnew(MACO3) * xnew(NO2)
  
  L = 2.9e-11 * xnew(OH)  &
    + rct(71,k)
  
  xnew(MPAN) = (xold(MPAN) + dt2 * P) / (1.0 + dt2 * L)


!-> NO3

  P = rct(6,k) * xnew(O3) * xnew(NO2)  &
    + rct(15,k) * xnew(OH) * xnew(HNO3)  &
    + rct(65,k) * xnew(N2O5)
  
  L = rct(9,k) * xnew(NO)  &
    + rct(11,k) * xnew(NO2)  &
    + rct(28,k) * xnew(HCHO)  &
    + rct(57,k) * xnew(MACRO2)  &
    + rct(61,k) * xnew(C5H8)  &
    + rct(64,k) * xnew(NO2)  &
    + rcphot(IDNO3_NO,k)  &
    + rcphot(IDNO3_NO2,k)
  
  xnew(NO3) = (xold(NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> N2O5

  P = rct(64,k) * xnew(NO2) * xnew(NO3)
  
  L = rct(65,k)
  
  xnew(N2O5) = (xold(N2O5) + dt2 * P) / (1.0 + dt2 * L)


!-> ISONO3

  P = rct(61,k) * xnew(C5H8) * xnew(NO3)  &
    + 2.0e-11 * xnew(ISONO3H) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(48,k) * xnew(HO2)
  
  xnew(ISONO3) = (xold(ISONO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HNO3

  P = rct(28,k) * xnew(NO3) * xnew(HCHO)  &
    + rct(66,k) * xnew(NO2) * xnew(OH)
  
  L = rct(15,k) * xnew(OH)  &
    + rcphot(IDHNO3,k)
  
  xnew(HNO3) = (xold(HNO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HONO

  P = rct(72,k) * xnew(OH) * xnew(NO)
  
  L = rct(18,k) * xnew(OH)  &
    + rcphot(IDHONO,k)
  
  xnew(HONO) = (xold(HONO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3

  P = 0.95* rct(33,k) * xnew(OH) * xnew(CH3CHO)  &
    + rct(38,k) * xnew(CH3CO3H) * xnew(OH)  &
    + rct(41,k) * xnew(MEKO2) * xnew(NO)  &
    + rct(50,k) * xnew(OH) * xnew(MGLYOX)  &
    + 0.684* rct(41,k) * xnew(MVKO2) * xnew(NO)  &
    + 0.41* rct(54,k) * xnew(MACR) * xnew(O3)  &
    + rct(59,k) * xnew(MACO3) * xnew(NO)  &
    + 1.51e-11 * xnew(MACO2H) * xnew(OH)  &
    + rct(68,k) * xnew(PAN)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)
  
  L = rct(34,k) * xnew(NO)  &
    + rct(35,k) * xnew(CH3O2)  &
    + rct(36,k) * xnew(CH3CO3)  &
    + rct(36,k) * xnew(CH3CO3)  &
    + rct(37,k) * xnew(HO2)  &
    + rct(67,k) * xnew(NO2)
  
  xnew(CH3CO3) = (xold(CH3CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> MACR

  P = 0.67* rct(51,k) * xnew(C5H8) * xnew(O3)  &
    + 0.32* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.1* rct(41,k) * xnew(ISONO3) * xnew(NO)
  
  L = rct(54,k) * xnew(O3)  &
    + rct(55,k) * xnew(OH)
  
  xnew(MACR) = (xold(MACR) + dt2 * P) / (1.0 + dt2 * L)


!-> ISNI

  P = 0.14* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.05* rct(41,k) * xnew(MVKO2) * xnew(NO)  &
    + 0.05* rct(56,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.05* rct(41,k) * xnew(ISNIR) * xnew(NO)  &
    + 0.85* rct(41,k) * xnew(ISONO3) * xnew(NO)
  
  L = 5.96e-11 * xnew(OH)
  
  xnew(ISNI) = (xold(ISNI) + dt2 * P) / (1.0 + dt2 * L)


!-> ISNIR

  P = 5.96e-11 * xnew(ISNI) * xnew(OH)  &
    + 3.7e-11 * xnew(ISNIRH) * xnew(OH)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(48,k) * xnew(HO2)
  
  xnew(ISNIR) = (xold(ISNIR) + dt2 * P) / (1.0 + dt2 * L)


!-> GLYOX

  P = rct(41,k) * xnew(MALO2) * xnew(NO)  &
    + rcphot(IDCH3O2H,k) * xnew(MALO2H)  &
    + rcemis(GLYOX,k)
  
  L = rct(49,k) * xnew(OH)  &
    + rcphot(IDCHOCHO_2CO,k)  &
    + rcphot(IDCHOCHO_HCHO,k)  &
    + rcphot(IDCHOCHO_2CHO,k)
  
  xnew(GLYOX) = (xold(GLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MGLYOX

  P = rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + rct(41,k) * xnew(MALO2) * xnew(NO)  &
    + 0.266* rct(41,k) * xnew(MVKO2) * xnew(NO)  &
    + 0.59* rct(54,k) * xnew(MACR) * xnew(O3)  &
    + rct(58,k) * xnew(ACETOL) * xnew(OH)  &
    + 0.82* rct(60,k) * xnew(MVK) * xnew(O3)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MALO2H)  &
    + rcemis(MGLYOX,k)
  
  L = rct(50,k) * xnew(OH)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(MGLYOX) = (xold(MGLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MAL

  P = rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYO2H)
  
  L = 5.58e-11 * xnew(OH)
  
  xnew(MAL) = (xold(MAL) + dt2 * P) / (1.0 + dt2 * L)


!-> MEK

  P = 0.65* rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + 2.15e-11 * xnew(BURO2H) * xnew(OH)  &
    + 0.65* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcemis(MEK,k)
  
  L = rct(43,k) * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(MEK) = (xold(MEK) + dt2 * P) / (1.0 + dt2 * L)


!-> MVK

  P = 0.26* rct(51,k) * xnew(C5H8) * xnew(O3)  &
    + 0.42* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.05* rct(41,k) * xnew(ISONO3) * xnew(NO)
  
  L = rct(53,k) * xnew(OH)  &
    + rct(60,k) * xnew(O3)
  
  xnew(MVK) = (xold(MVK) + dt2 * P) / (1.0 + dt2 * L)


!-> HCHO

  P = rct(21,k) * xnew(CH3O2) * xnew(NO)  &
    + 2.* rct(22,k) * xnew(CH3O2) * xnew(CH3O2)  &
    + rct(23,k) * xnew(CH3O2) * xnew(CH3O2)  &
    + rct(24,k) * xnew(OH) * xnew(CH3OH)  &
    + 0.1* rct(25,k) * xnew(HO2) * xnew(CH3O2)  &
    + 0.4* rct(26,k) * xnew(CH3OOH) * xnew(OH)  &
    + rct(35,k) * xnew(CH3O2) * xnew(CH3CO3)  &
    + 2.* rct(41,k) * xnew(ETRO2) * xnew(NO)  &
    + 1.14* rct(44,k) * xnew(C2H4) * xnew(O3)  &
    + 0.545* rct(45,k) * xnew(O3) * xnew(C3H6)  &
    + rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.8* rct(51,k) * xnew(C5H8) * xnew(O3)  &
    + 0.74* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.266* rct(41,k) * xnew(MVKO2) * xnew(NO)  &
    + 0.534* rct(54,k) * xnew(MACR) * xnew(O3)  &
    + rct(57,k) * xnew(MACRO2) * xnew(NO3)  &
    + rct(59,k) * xnew(MACO3) * xnew(NO)  &
    + 0.8* rct(60,k) * xnew(MVK) * xnew(O3)  &
    + 0.15* rct(41,k) * xnew(ISONO3) * xnew(NO)  &
    + 1.51e-11 * xnew(MACO2H) * xnew(OH)
  
  !-> HCHO cont.
    P = P +   rcphot(IDCHOCHO_HCHO,k) * xnew(GLYOX)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + 1.56* rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)  &
    + rcemis(HCHO,k)
  
  L = rct(27,k) * xnew(OH)  &
    + rct(28,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k)  &
    + rcphot(IDHCHO_H2,k)
  
  xnew(HCHO) = (xold(HCHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CHO

  P = rct(30,k) * xnew(C2H5O2) * xnew(NO)  &
    + 8.01e-12 * xnew(C2H5OOH) * xnew(OH)  &
    + rct(39,k) * xnew(OH) * xnew(C2H5OH)  &
    + 0.35* rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + rct(41,k) * xnew(MEKO2) * xnew(NO)  &
    + 1.38e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.545* rct(45,k) * xnew(O3) * xnew(C3H6)  &
    + rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.684* rct(41,k) * xnew(MVKO2) * xnew(NO)  &
    + 0.04* rct(60,k) * xnew(MVK) * xnew(O3)  &
    + 0.95* rct(41,k) * xnew(ISNIR) * xnew(NO)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + 0.22* rcphot(IDCH3O2H,k) * xnew(ETRO2H)  &
    + 0.35* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(PRRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcemis(CH3CHO,k)
  
  L = rct(33,k) * xnew(OH)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(CH3CHO) = (xold(CH3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H6

  P = rcemis(C2H6,k)
  
  L = rct(29,k) * xnew(OH)
  
  xnew(C2H6) = (xold(C2H6) + dt2 * P) / (1.0 + dt2 * L)


!-> NC4H10

  P = rcemis(NC4H10,k)
  
  L = rct(40,k) * xnew(OH)
  
  xnew(NC4H10) = (xold(NC4H10) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H4

  P = rcemis(C2H4,k)
  
  L = rct(44,k) * xnew(O3)  &
    + rct(69,k) * xnew(OH)
  
  xnew(C2H4) = (xold(C2H4) + dt2 * P) / (1.0 + dt2 * L)


!-> C3H6

  P = 0.07* rct(51,k) * xnew(C5H8) * xnew(O3)  &
    + rcemis(C3H6,k)
  
  L = rct(45,k) * xnew(O3)  &
    + rct(70,k) * xnew(OH)
  
  xnew(C3H6) = (xold(C3H6) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYL

  P = rcemis(OXYL,k)
  
  L = 1.36e-11 * xnew(OH)
  
  xnew(OXYL) = (xold(OXYL) + dt2 * P) / (1.0 + dt2 * L)


!-> C5H8

  P = rcemis(C5H8,k)  &
    + RCBIO(NATBIO%C5H8,K)
  
  L = rct(51,k) * xnew(O3)  &
    + rct(52,k) * xnew(OH)  &
    + rct(61,k) * xnew(NO3)
  
  xnew(C5H8) = (xold(C5H8) + dt2 * P) / (1.0 + dt2 * L)


!-> APINENE

  ! P = 0.0
  
  L = rct(73,k) * xnew(O3)  &
    + rct(74,k) * xnew(OH)  &
    + rct(75,k) * xnew(NO3)
  
  xnew(APINENE) = xold(APINENE) / (1.0 + dt2 * L)


!-> BVOCTEST

  P = RCBIO(NATBIO%TERP,K)
  
  ! L = 0.0
  
  xnew(BVOCTEST) = xold(BVOCTEST) + dt2 * P


!-> CH3OOH

  P = 0.9* rct(25,k) * xnew(HO2) * xnew(CH3O2)
  
  L = rct(26,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(CH3OOH) = (xold(CH3OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OOH

  P = rct(31,k) * xnew(C2H5O2) * xnew(HO2)
  
  L = 8.01e-12 * xnew(OH)  &
    + rct(32,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(C2H5OOH) = (xold(C2H5OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> BURO2H

  P = 0.95* rct(42,k) * xnew(SECC4H9O2) * xnew(HO2)
  
  L = rct(32,k) * xnew(OH)  &
    + 2.15e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(BURO2H) = (xold(BURO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2H

  P = 1.2e-11 * xnew(ETRO2) * xnew(HO2)
  
  L = 1.38e-11 * xnew(OH)  &
    + rct(32,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ETRO2H) = (xold(ETRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> PRRO2H

  P = 0.795* rct(46,k) * xnew(PRRO2) * xnew(HO2)
  
  L = 2.44e-11 * xnew(OH)  &
    + rct(32,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(PRRO2H) = (xold(PRRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYO2H

  P = 0.33* rct(47,k) * xnew(OXYO2) * xnew(HO2)
  
  L = 4.2e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(OXYO2H) = (xold(OXYO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2H

  P = rct(42,k) * xnew(MEKO2) * xnew(HO2)
  
  L = rct(38,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(MEKO2H) = (xold(MEKO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MALO2H

  P = rct(48,k) * xnew(MALO2) * xnew(HO2)
  
  L = rct(38,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(MALO2H) = (xold(MALO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MVKO2H

  P = rct(42,k) * xnew(MVKO2) * xnew(HO2)
  
  L = 2.2e-11 * xnew(OH)
  
  xnew(MVKO2H) = (xold(MVKO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MACROOH

  P = rct(42,k) * xnew(MACRO2) * xnew(HO2)
  
  L = 2.82e-11 * xnew(OH)
  
  xnew(MACROOH) = (xold(MACROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> MACO3H

  P = 0.71* rct(62,k) * xnew(MACO3) * xnew(HO2)
  
  L = 1.87e-11 * xnew(OH)
  
  xnew(MACO3H) = (xold(MACO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> MACO2H

  P = 0.29* rct(62,k) * xnew(MACO3) * xnew(HO2)
  
  L = 1.51e-11 * xnew(OH)
  
  xnew(MACO2H) = (xold(MACO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2H

  P = 0.97* rct(48,k) * xnew(ISRO2) * xnew(HO2)
  
  L = 7.5e-11 * xnew(OH)
  
  xnew(ISRO2H) = (xold(ISRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> H2O2

  P = rct(16,k) * xnew(HO2) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2) * xnew(HO2)  &
    + 0.14* rct(44,k) * xnew(C2H4) * xnew(O3)  &
    + 0.09* rct(45,k) * xnew(O3) * xnew(C3H6)  &
    + 0.124* rct(54,k) * xnew(MACR) * xnew(O3)
  
  L = rct(13,k) * xnew(OH)  &
    + rcphot(IDH2O2,k)
  
  xnew(H2O2) = (xold(H2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3H

  P = 0.41* rct(37,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = rct(38,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(CH3CO3H) = (xold(CH3CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> ISONO3H

  P = rct(48,k) * xnew(ISONO3) * xnew(HO2)
  
  L = 2.0e-11 * xnew(OH)
  
  xnew(ISONO3H) = (xold(ISONO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> ISNIRH

  P = rct(48,k) * xnew(ISNIR) * xnew(HO2)
  
  L = 3.7e-11 * xnew(OH)
  
  xnew(ISNIRH) = (xold(ISNIRH) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OH

  P = rct(23,k) * xnew(CH3O2) * xnew(CH3O2)  &
    + rcemis(CH3OH,k)
  
  L = rct(24,k) * xnew(OH)
  
  xnew(CH3OH) = (xold(CH3OH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OH

  P = rcemis(C2H5OH,k)
  
  L = rct(39,k) * xnew(OH)
  
  xnew(C2H5OH) = (xold(C2H5OH) + dt2 * P) / (1.0 + dt2 * L)


!-> ACETOL

  P = 2.44e-11 * xnew(PRRO2H) * xnew(OH)  &
    + 0.95* rct(56,k) * xnew(MACRO2) * xnew(NO)  &
    + rct(57,k) * xnew(MACRO2) * xnew(NO3)  &
    + 0.95* rct(41,k) * xnew(ISNIR) * xnew(NO)  &
    + 2.9e-11 * xnew(MPAN) * xnew(OH)
  
  L = rct(58,k) * xnew(OH)
  
  xnew(ACETOL) = (xold(ACETOL) + dt2 * P) / (1.0 + dt2 * L)


!-> H2

  P = rcphot(IDHCHO_H2,k) * xnew(HCHO)
  
  L = rct(14,k) * xnew(OH)
  
  xnew(H2) = (xold(H2) + dt2 * P) / (1.0 + dt2 * L)


!-> CO

  P = rct(27,k) * xnew(OH) * xnew(HCHO)  &
    + rct(28,k) * xnew(NO3) * xnew(HCHO)  &
    + 0.05* rct(33,k) * xnew(OH) * xnew(CH3CHO)  &
    + 0.63* rct(44,k) * xnew(C2H4) * xnew(O3)  &
    + 0.56* rct(45,k) * xnew(O3) * xnew(C3H6)  &
    + rct(49,k) * xnew(OH) * xnew(GLYOX)  &
    + rct(49,k) * xnew(OH) * xnew(GLYOX)  &
    + rct(50,k) * xnew(OH) * xnew(MGLYOX)  &
    + 0.05* rct(51,k) * xnew(C5H8) * xnew(O3)  &
    + 0.82* rct(54,k) * xnew(MACR) * xnew(O3)  &
    + 0.95* rct(56,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.05* rct(60,k) * xnew(MVK) * xnew(O3)  &
    + 2.9e-11 * xnew(MPAN) * xnew(OH)  &
    + rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDHCHO_H2,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + 2* rcphot(IDCHOCHO_2CO,k) * xnew(GLYOX)  &
    + rcphot(IDCHOCHO_HCHO,k) * xnew(GLYOX)  &
    + 2* rcphot(IDCHOCHO_2CHO,k) * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)
  
  !-> CO cont.
    P = P +   rcemis(CO,k)
  
  L = 1.44e-13 * xnew(OH)  &
    + rct(20,k) * xnew(OH)
  
  xnew(CO) = (xold(CO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH4

  P = 0.1* rct(45,k) * xnew(O3) * xnew(C3H6)
  
  L = rct(19,k) * xnew(OH)
  
  xnew(CH4) = (xold(CH4) + dt2 * P) / (1.0 + dt2 * L)


!-> SO2

  P = rcemis(SO2,k)
  
  L = 1.0e-5
  
  xnew(SO2) = (xold(SO2) + dt2 * P) / (1.0 + dt2 * L)


!-> HYPROPO2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HYPROPO2! xnew(HYPROPO2) = max(0.0, xold(HYPROPO2))


!-> IPROPOLO2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for IPROPOLO2! xnew(IPROPOLO2) = max(0.0, xold(IPROPOLO2))


!-> newC5H8

  P = RCBIO(NATBIO%C5H8,K)
  
  ! L = 0.0
  
  xnew(newC5H8) = xold(newC5H8) + dt2 * P


!-> newBVOCTEST

  P = RCBIO(NATBIO%TERP,K)
  
  ! L = 0.0
  
  xnew(newBVOCTEST) = xold(newBVOCTEST) + dt2 * P


!-> POM_f_WOOD

  P = rcemis(POM_F_WOOD,k)
  
  ! L = 0.0
  
  xnew(POM_f_WOOD) = xold(POM_f_WOOD) + dt2 * P


!-> POM_f_FFUEL

  P = rcemis(POM_F_FFUEL,k)
  
  ! L = 0.0
  
  xnew(POM_f_FFUEL) = xold(POM_f_FFUEL) + dt2 * P


!-> POM_c_FFUEL

  P = rcemis(POM_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(POM_c_FFUEL) = xold(POM_c_FFUEL) + dt2 * P


!-> EC_f_WOOD_new

  P = rcemis(EC_F_WOOD_NEW,k)
  
  L = rct(93,k)
  
  xnew(EC_f_WOOD_new) = (xold(EC_f_WOOD_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_WOOD_age

  P = rcemis(EC_F_WOOD_AGE,k)  &
    + rct(93,k) * xnew(EC_f_WOOD_new)
  
  ! L = 0.0
  
  xnew(EC_f_WOOD_age) = xold(EC_f_WOOD_age) + dt2 * P


!-> EC_c_WOOD

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for EC_c_WOOD! xnew(EC_c_WOOD) = max(0.0, xold(EC_c_WOOD))


!-> EC_f_FFUEL_new

  P = rcemis(EC_F_FFUEL_NEW,k)
  
  L = rct(93,k)
  
  xnew(EC_f_FFUEL_new) = (xold(EC_f_FFUEL_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_FFUEL_age

  P = rcemis(EC_F_FFUEL_AGE,k)  &
    + rct(93,k) * xnew(EC_f_FFUEL_new)
  
  ! L = 0.0
  
  xnew(EC_f_FFUEL_age) = xold(EC_f_FFUEL_age) + dt2 * P


!-> EC_c_FFUEL

  P = rcemis(EC_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(EC_c_FFUEL) = xold(EC_c_FFUEL) + dt2 * P


!-> REMPPM25

  P = rcemis(REMPPM25,k)
  
  ! L = 0.0
  
  xnew(REMPPM25) = xold(REMPPM25) + dt2 * P


!-> REMPPM_c

  P = rcemis(REMPPM_C,k)
  
  ! L = 0.0
  
  xnew(REMPPM_c) = xold(REMPPM_c) + dt2 * P


!-> FFIRE_OM

  P = rcemis(FFIRE_OM,k)
  
  ! L = 0.0
  
  xnew(FFIRE_OM) = xold(FFIRE_OM) + dt2 * P


!-> FFIRE_BC

  P = rcemis(FFIRE_BC,k)
  
  ! L = 0.0
  
  xnew(FFIRE_BC) = xold(FFIRE_BC) + dt2 * P


!-> FFIRE_REMPPM25

  P = rcemis(FFIRE_REMPPM25,k)
  
  ! L = 0.0
  
  xnew(FFIRE_REMPPM25) = xold(FFIRE_REMPPM25) + dt2 * P


!-> OM25_BGND

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_BGND! xnew(OM25_BGND) = max(0.0, xold(OM25_BGND))


!-> OM25_p

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_p! xnew(OM25_p) = max(0.0, xold(OM25_p))


!-> TERPPeroxy

  P = rct(73,k) * xnew(O3) * xnew(APINENE)  &
    + rct(74,k) * xnew(OH) * xnew(APINENE)  &
    + rct(75,k) * xnew(NO3) * xnew(APINENE)
  
  L = rct(41,k) * xnew(NO)  &
    + rct(76,k) * xnew(HO2)
  
  xnew(TERPPeroxy) = (xold(TERPPeroxy) + dt2 * P) / (1.0 + dt2 * L)


!-> VBS_TEST

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for VBS_TEST! xnew(VBS_TEST) = max(0.0, xold(VBS_TEST))


!-> ASOC_ng100

  P = rct(77,k) * xnew(OH) * xnew(ASOC_ug1)
  
  ! L = 0.0
  
  xnew(ASOC_ng100) = xold(ASOC_ng100) + dt2 * P


!-> ASOC_ug1

  P = 0.00206* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.010294* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 0.008413* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 0.315476* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + rct(78,k) * xnew(OH) * xnew(ASOC_ug10)
  
  L = rct(77,k) * xnew(OH)
  
  xnew(ASOC_ug1) = (xold(ASOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug10

  P = 0.108* rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + 0.2132* rct(42,k) * xnew(SECC4H9O2) * xnew(HO2)  &
    + 0.01029* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.018529* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 0.820238* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 1.261905* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + rct(79,k) * xnew(OH) * xnew(ASOC_ug1e2)
  
  L = rct(78,k) * xnew(OH)
  
  xnew(ASOC_ug10) = (xold(ASOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e2

  P = 0.078235* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.123529* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 1.261905* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 1.577381* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + rct(80,k) * xnew(OH) * xnew(ASOC_ug1e3)
  
  L = rct(79,k) * xnew(OH)
  
  xnew(ASOC_ug1e2) = (xold(ASOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e3

  P = 0.30882* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.463235* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 1.829762* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.208333* rct(47,k) * xnew(OXYO2) * xnew(HO2)
  
  L = rct(80,k) * xnew(OH)
  
  xnew(ASOC_ug1e3) = (xold(ASOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ng100

  P = 0.9* rct(77,k) * xnew(OH) * xnew(ASOC_ug1)  &
    + 1.075* rct(81,k) * xnew(OH) * xnew(NON_C_ASOA_ug1)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ng100) = xold(NON_C_ASOA_ng100) + dt2 * P


!-> NON_C_ASOA_ug1

  P = 0.0173* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.08647* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 0.11105* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 4.16429* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + 0.9* rct(78,k) * xnew(OH) * xnew(ASOC_ug10)  &
    + 1.075* rct(82,k) * xnew(OH) * xnew(NON_C_ASOA_ug10)
  
  L = rct(81,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1) = (xold(NON_C_ASOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug10

  P = 0.90753* rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + 1.7912* rct(42,k) * xnew(SECC4H9O2) * xnew(HO2)  &
    + 0.08647* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.15565* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 10.8271* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 16.6571* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + 0.9* rct(79,k) * xnew(OH) * xnew(ASOC_ug1e2)  &
    + 1.075* rct(83,k) * xnew(OH) * xnew(NON_C_ASOA_ug1e2)
  
  L = rct(82,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug10) = (xold(NON_C_ASOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug1e2

  P = 0.65718* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 1.03765* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 16.6571* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 20.82143* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + 0.9* rct(80,k) * xnew(OH) * xnew(ASOC_ug1e3)  &
    + 1.075* rct(84,k) * xnew(OH) * xnew(NON_C_ASOA_ug1e3)
  
  L = rct(83,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1e2) = (xold(NON_C_ASOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug1e3

  P = 2.59412* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 3.89118* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 24.15286* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 29.15* rct(47,k) * xnew(OXYO2) * xnew(HO2)
  
  L = rct(84,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1e3) = (xold(NON_C_ASOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ng100

  P = rct(85,k) * xnew(OH) * xnew(BSOC_ug1)
  
  ! L = 0.0
  
  xnew(BSOC_ng100) = xold(BSOC_ng100) + dt2 * P


!-> BSOC_ug1

  P = 0.002833* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.0255* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 0.08* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 0.715333* rct(76,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + rct(86,k) * xnew(OH) * xnew(BSOC_ug10)
  
  L = rct(85,k) * xnew(OH)
  
  xnew(BSOC_ug1) = (xold(BSOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug10

  P = 0.065167* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.085* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 0.813333* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 0.612* rct(76,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + rct(87,k) * xnew(OH) * xnew(BSOC_ug1e2)
  
  L = rct(86,k) * xnew(OH)
  
  xnew(BSOC_ug10) = (xold(BSOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e2

  P = 0.0425* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.0425* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 1.34* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 2.391333* rct(76,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + rct(88,k) * xnew(OH) * xnew(BSOC_ug1e3)
  
  L = rct(87,k) * xnew(OH)
  
  xnew(BSOC_ug1e2) = (xold(BSOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e3

  P = 3.333333* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 4.05* rct(76,k) * xnew(HO2) * xnew(TERPPeroxy)
  
  L = rct(88,k) * xnew(OH)
  
  xnew(BSOC_ug1e3) = (xold(BSOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ng100

  P = 0.9* rct(85,k) * xnew(OH) * xnew(BSOC_ug1)  &
    + 1.075* rct(89,k) * xnew(OH) * xnew(NON_C_BSOA_ug1)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ng100) = xold(NON_C_BSOA_ng100) + dt2 * P


!-> NON_C_BSOA_ug1

  P = 0.034* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.306* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 0.672* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 6.009* rct(76,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + 0.9* rct(86,k) * xnew(OH) * xnew(BSOC_ug10)  &
    + 1.075* rct(90,k) * xnew(OH) * xnew(NON_C_BSOA_ug10)
  
  L = rct(89,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1) = (xold(NON_C_BSOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug10

  P = 0.782* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 1.02* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 6.832* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 5.1408* rct(76,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + 0.9* rct(87,k) * xnew(OH) * xnew(BSOC_ug1e2)  &
    + 1.075* rct(91,k) * xnew(OH) * xnew(NON_C_BSOA_ug1e2)
  
  L = rct(90,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug10) = (xold(NON_C_BSOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug1e2

  P = 0.51* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.51* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 11.256* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 20.0872* rct(76,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + 0.9* rct(88,k) * xnew(OH) * xnew(BSOC_ug1e3)  &
    + 1.075* rct(92,k) * xnew(OH) * xnew(NON_C_BSOA_ug1e3)
  
  L = rct(91,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1e2) = (xold(NON_C_BSOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug1e3

  P = 28.* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 34.02* rct(76,k) * xnew(HO2) * xnew(TERPPeroxy)
  
  L = rct(92,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1e3) = (xold(NON_C_BSOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> FFFUEL_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for FFFUEL_ng10! xnew(FFFUEL_ng10) = max(0.0, xold(FFFUEL_ng10))


!-> WOODOA_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for WOODOA_ng10! xnew(WOODOA_ng10) = max(0.0, xold(WOODOA_ng10))


!-> FFIREOA_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for FFIREOA_ng10! xnew(FFIREOA_ng10) = max(0.0, xold(FFIREOA_ng10))


