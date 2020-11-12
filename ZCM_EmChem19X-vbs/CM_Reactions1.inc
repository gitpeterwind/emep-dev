!-> OD

  P = rcphot(IDO3_O1D,k) * xnew(O3)
  
  L = rct(2,k)  &
    + rct(3,k)  &
    + rct(4,k)
  
  xnew(OD) = (xold(OD) + dt2 * P) / (1.0 + dt2 * L)


!-> OP

  P = rct(2,k) * xnew(OD)  &
    + rct(3,k) * xnew(OD)  &
    + rcphot(IDO3_O3P,k) * xnew(O3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + 0.873* rcphot(IDNO3,k) * xnew(NO3)
  
  L = rct(1,k)  &
    + rct(10,k) * xnew(NO)
  
  xnew(OP) = (xold(OP) + dt2 * P) / (1.0 + dt2 * L)


!-> OH

  P = 2* rct(4,k) * xnew(OD)  &
    + rct(14,k) * xnew(O3) * xnew(HO2)  &
    + rct(16,k) * xnew(NO) * xnew(HO2)  &
    + 4.0e-12 * xnew(HO2) * xnew(NO3)  &
    + 0.4* rct(34,k) * xnew(CH3OOH) * xnew(OH)  &
    + 0.69* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.454* rct(45,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.89* 3.3e-11 * xnew(HCOCH2OOH) * xnew(OH)  &
    + 0.847* 2.4e-11 * xnew(BURO2H) * xnew(OH)  &
    + 0.826* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + 0.75* rct(45,k) * xnew(C2H5CO3) * xnew(HO2)  &
    + 0.79* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.75* rct(45,k) * xnew(HOCH2CO3) * xnew(HO2)  &
    + 0.17* rct(57,k) * xnew(C2H4) * xnew(O3)  &
    + 0.3646* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + 0.584* rct(68,k) * xnew(GLYOX) * xnew(OH)  &
    + 0.25* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.93* 6.0e-11 * xnew(ISRO2H) * xnew(OH)  &
    + 0.2392* rct(74,k) * xnew(MACR) * xnew(O3)  &
    + rct(75,k) * xnew(MACRO2)
  
  !-> OH cont.
    P = P +   0.44* rct(45,k) * xnew(MACO3) * xnew(HO2)  &
    + 3.77e-11 * xnew(MACROOH) * xnew(OH)  &
    + 0.36* rct(80,k) * xnew(MVK) * xnew(O3)  &
    + 4.5e-11 * xnew(MVKO2H) * xnew(OH)  &
    + 0.89* 2.4e-17 * xnew(ISON2) * xnew(O3)  &
    + rct(45,k) * xnew(NC4CO3) * xnew(HO2)  &
    + rct(84,k) * xnew(ISRO2)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + 0.44* rct(45,k) * xnew(CO2H3CO3) * xnew(HO2)  &
    + 0.7* rct(88,k) * xnew(APINENE) * xnew(O3)  &
    + 2* rcphot(IDH2O2,k) * xnew(H2O2)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.333* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + 0.58* rcphot(IDCH3COY,k) * xnew(CO23C3CHO)  &
    + rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)
  
  !-> OH cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(CH3CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + rcphot(IDMEK,k) * xnew(MEKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(HCOCH2OOH)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HCOCH2OOH)  &
    + rcphot(IDHCHO_H,k) *1.87 * xnew(MACROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(MVKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(OXNCATCOOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OXNNCATOOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OCATEC1OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(DM123OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NOXYOLOOH)  &
    + rcphot(IDMEK,k) *0.78 * xnew(NOXYOLOOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NOXYOL1OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYBPEROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(CO2H3CO3H)
  
  L = rct(5,k) * xnew(HO2)  &
    + rct(6,k) * xnew(H2O2)  &
    + rct(7,k) * xnew(H2)  &
    + rct(13,k) * xnew(O3)  &
    + rct(20,k) * xnew(NO2)  &
    + rct(21,k) * xnew(HNO3)  &
    + rct(22,k) * xnew(NO)  &
    + rct(23,k) * xnew(HONO)  &
    + 2.0e-11 * xnew(NO3)  &
    + rct(26,k) * xnew(HO2NO2)  &
    + rct(27,k) * xnew(CH4)  &
    + rct(28,k) * xnew(CO)  &
    + rct(33,k) * xnew(CH3OH)  &
    + rct(34,k) * xnew(CH3OOH)  &
    + rct(35,k) * xnew(HCHO)  &
    + rct(37,k) * xnew(C2H6)  &
    + 1.17e-11 * xnew(C2H5OOH)  &
    + rct(40,k) * xnew(CH3CHO)  &
    + 1.1e-11 * xnew(CH3CO3H)  &
    + 3.3e-11 * xnew(HCOCH2OOH)
  
  !-> OH cont.
    L = L +   3.22e-12 * xnew(C2H5OH)  &
    + rct(48,k) * xnew(NC4H10)  &
    + 2.4e-11 * xnew(BURO2H)  &
    + 8.6e-13 * xnew(SC4H9NO3)  &
    + 1.1e-12 * xnew(MEK)  &
    + 2.1e-11 * xnew(MEKO2H)  &
    + 2.3e-13 * xnew(BIACET)  &
    + 1.27e-12 * xnew(PPN)  &
    + 1.39e-11 * xnew(EGLYOX)  &
    + rct(53,k) * xnew(C2H4)  &
    + 1.7475e-11 * xnew(ETRO2H)  &
    + 1.00e-11 * xnew(HOCH2CHO)  &
    + 2.73e-12 * xnew(HOCH2CO2H)  &
    + 1.12e-12 * xnew(PHAN)  &
    + rct(58,k) * xnew(C3H6)  &
    + 6.70e-13 * xnew(NOA)  &
    + rct(62,k) * xnew(BENZENE)  &
    + rct(63,k) * xnew(TOLUENE)  &
    + 1.36e-11 * xnew(OXYL)  &
    + 6.1e-11 * xnew(OXYLOL)
  
  !-> OH cont.
    L = L +   2.05e-10 * xnew(OXYCATECH)  &
    + 6.83e-12 * xnew(OXNCATECH)  &
    + rct(65,k) * xnew(NC4MDCO2H)  &
    + rct(65,k) * xnew(OXNCATCOOH)  &
    + rct(65,k) * xnew(OXNNCATOOH)  &
    + 2.05e-10 * xnew(OCATEC1OOH)  &
    + 1.53e-12 * xnew(OXY1OHNO2)  &
    + 1.53e-12 * xnew(NOXYOL1OOH)  &
    + 4.65e-11 * xnew(DM123OOH)  &
    + 1.43e-10 * xnew(NOXYOLOOH)  &
    + 7.62e-11 * xnew(OXYBPEROOH)  &
    + 4.45e-11 * xnew(BZFUONE)  &
    + 6.2e-11 * xnew(C5DICARB)  &
    + rct(65,k) * xnew(C5DICAROOH)  &
    + rct(68,k) * xnew(GLYOX)  &
    + rct(69,k) * xnew(MGLYOX)  &
    + rct(70,k) * xnew(C5H8)  &
    + 6.0e-11 * xnew(ISRO2H)  &
    + 1.16e-11 * xnew(IEPOX)  &
    + 3.25e-11 * xnew(ISON)
  
  !-> OH cont.
    L = L +   1.55e-12 * xnew(NALD)  &
    + rct(73,k) * xnew(MACR)  &
    + 2.9e-11 * xnew(MPAN)  &
    + 3.77e-11 * xnew(MACROOH)  &
    + rct(78,k) * xnew(ACETOL)  &
    + rct(79,k) * xnew(MVK)  &
    + 4.5e-11 * xnew(MVKO2H)  &
    + 2.96e-12 * xnew(BIACETOH)  &
    + 1.23e-11 * xnew(CO23C3CHO)  &
    + 2.45e-11 * xnew(CO2H3CHO)  &
    + 1.33e-12 * xnew(MVKNO3)  &
    + 4.16e-11 * xnew(ISON2)  &
    + 2.16e-11 * xnew(C5PAN18)  &
    + 7.61e-11 * xnew(HPALD)  &
    + 7.34e-12 * xnew(CO2H3CO3H)  &
    + 3.74e-12 * xnew(C4PAN6)  &
    + 1.51e-11 * xnew(MACO2H)  &
    + rct(85,k) * xnew(APINENE)  &
    + rct(87,k) * xnew(TERPOOH)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(SO2)
  
  !-> OH cont.
    L = L +   rct(20,k) * xnew(shipNOx)
  
  xnew(OH) = (xold(OH) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2

  P = rct(6,k) * xnew(OH) * xnew(H2O2)  &
    + rct(7,k) * xnew(OH) * xnew(H2)  &
    + rct(13,k) * xnew(O3) * xnew(OH)  &
    + 2.0e-11 * xnew(OH) * xnew(NO3)  &
    + rct(25,k) * xnew(HO2NO2)  &
    + rct(28,k) * xnew(CO) * xnew(OH)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(29,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(30,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + rct(33,k) * xnew(CH3OH) * xnew(OH)  &
    + rct(35,k) * xnew(HCHO) * xnew(OH)  &
    + rct(36,k) * xnew(HCHO) * xnew(NO3)  &
    + rct(38,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + 0.6* 3.10e-13 * xnew(RO2POOL) * xnew(C2H5O2)  &
    + rct(46,k) * xnew(HCOCH2O2) * xnew(NO)  &
    + 2.3e-12 * xnew(HCOCH2O2) * xnew(NO3)  &
    + 0.6* 2.0e-12 * xnew(RO2POOL) * xnew(HCOCH2O2)  &
    + 0.95* 3.22e-12 * xnew(C2H5OH) * xnew(OH)  &
    + 0.654* rct(46,k) * xnew(C4H9O2) * xnew(NO)
  
  !-> HO2 cont.
    P = P +   0.716* 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + rct(55,k) * xnew(HOCH2CH2O)  &
    + rct(56,k) * xnew(HOCH2CH2O)  &
    + 0.2* 1.00e-11 * xnew(HOCH2CHO) * xnew(OH)  &
    + 0.75* rct(45,k) * xnew(HOCH2CO3) * xnew(HO2)  &
    + rct(42,k) * xnew(HOCH2CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(HOCH2CO3) * xnew(NO3)  &
    + 0.7* 1.00e-11 * xnew(RO2POOL) * xnew(HOCH2CO3)  &
    + 2.73e-12 * xnew(HOCH2CO2H) * xnew(OH)  &
    + 0.27* rct(57,k) * xnew(C2H4) * xnew(O3)  &
    + 0.038* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + rct(46,k) * xnew(PRRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(PRRO2) * xnew(NO3)  &
    + 0.6* 8.8e-13 * xnew(RO2POOL) * xnew(PRRO2)  &
    + 0.87* rct(46,k) * xnew(NRN9O2) * xnew(NO)  &
    + 0.87* 2.3e-12 * xnew(NRN9O2) * xnew(NO3)  &
    + rct(46,k) * xnew(IBUTOLBO2) * xnew(NO)  &
    + 2.3e-12 * xnew(IBUTOLBO2) * xnew(NO3)  &
    + 0.654* rct(62,k) * xnew(BENZENE) * xnew(OH)  &
    + 0.29* rct(63,k) * xnew(TOLUENE) * xnew(OH)
  
  !-> HO2 cont.
    P = P +   0.16* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + 0.83* 6.1e-11 * xnew(OXYLOL) * xnew(OH)  &
    + rct(46,k) * xnew(OXNCATECO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXNCATECO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(OXNCATECO2)  &
    + 0.862* rct(46,k) * xnew(OXYBIPERO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYBIPERO2) * xnew(NO3)  &
    + rct(46,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(66,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.416* rct(68,k) * xnew(GLYOX) * xnew(OH)  &
    + 0.1725* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.9* rct(46,k) * xnew(ISRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3) * xnew(ISRO2)  &
    + 1.45e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.6* 3.25e-11 * xnew(ISON) * xnew(OH)  &
    + 0.1408* rct(74,k) * xnew(MACR) * xnew(O3)  &
    + rct(46,k) * xnew(MACRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MACRO2) * xnew(NO3)  &
    + rct(78,k) * xnew(ACETOL) * xnew(OH)  &
    + 0.28* rct(80,k) * xnew(MVK) * xnew(O3)
  
  !-> HO2 cont.
    P = P +   0.28944* rct(46,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.3* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.18* 1.22e-12 * xnew(RO2POOL) * xnew(HMVKO2)  &
    + 2.96e-12 * xnew(OH) * xnew(BIACETOH)  &
    + 0.67* 1.33e-12 * xnew(MVKNO3) * xnew(OH)  &
    + rct(46,k) * xnew(NISOPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NISOPO2) * xnew(NO3)  &
    + 0.6* 1.3e-12 * xnew(RO2POOL) * xnew(NISOPO2)  &
    + 0.445* 2.4e-17 * xnew(ISON2) * xnew(O3)  &
    + rct(45,k) * xnew(NC4CO3) * xnew(HO2)  &
    + rct(42,k) * xnew(NC4CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(NC4CO3) * xnew(NO3)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(NC4CO3)  &
    + rct(83,k) * xnew(ISRO2)  &
    + 0.44* rct(45,k) * xnew(CO2H3CO3) * xnew(HO2)  &
    + rct(42,k) * xnew(CO2H3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CO2H3CO3) * xnew(NO3)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(CO2H3CO3)  &
    + rct(46,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)
  
  !-> HO2 cont.
    P = P +   rct(88,k) * xnew(APINENE) * xnew(O3)  &
    + 0.667* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + 2* rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + 2* rcphot(IDHCHO_H,k) *0.6 * xnew(HOCH2CHO)  &
    + 2* rcphot(IDHCHO_H,k) *0.6 * xnew(CO2H3CHO)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + 0.5* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(EGLYOX)  &
    + 0.35* rcphot(IDRCOCHO,k) * xnew(CH3COCO2H)  &
    + 0.41* rcphot(IDCH3COY,k) * xnew(CO23C3CHO)  &
    + rcphot(IDRCOCHO,k) * xnew(CO23C3CHO)  &
    + rcphot(IDHCHO_H2,k) *0.04767 * xnew(MACR)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + rcphot(IDNO2,k) *0.2 * xnew(C5DICARB)  &
    + 0.68* rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)
  
  !-> HO2 cont.
    P = P +   0.716* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(HCOCH2OOH)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HCOCH2OOH)  &
    + rcphot(IDHCHO_H,k) *1.87 * xnew(MACROOH)  &
    + 0.3* rcphot(IDCH3O2H,k) * xnew(MVKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(OXNCATCOOH)  &
    + rcphot(IDMEK,k) *0.78 * xnew(NOXYOLOOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OXYBPEROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(CO2H3CO3H)  &
    + 0.652* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rcphot(IDHCHO_H,k) *0.1814 * xnew(ISON)  &
    + rcphot(IDMEK,k) *133.34 * xnew(ISON2)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)
  
  L = rct(5,k) * xnew(OH)  &
    + rct(8,k) * xnew(HO2)  &
    + rct(8,k) * xnew(HO2)  &
    + rct(9,k) * xnew(HO2)  &
    + rct(9,k) * xnew(HO2)  &
    + rct(14,k) * xnew(O3)  &
    + rct(16,k) * xnew(NO)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(24,k) * xnew(NO2)  &
    + rct(32,k) * xnew(CH3O2)  &
    + rct(39,k) * xnew(C2H5O2)  &
    + rct(45,k) * xnew(CH3CO3)  &
    + rct(47,k) * xnew(HCOCH2O2)  &
    + rct(49,k) * xnew(C4H9O2)  &
    + rct(49,k) * xnew(MEKO2)  &
    + rct(45,k) * xnew(C2H5CO3)  &
    + rct(54,k) * xnew(ETRO2)  &
    + rct(45,k) * xnew(HOCH2CO3)  &
    + rct(60,k) * xnew(PRRO2)  &
    + rct(60,k) * xnew(NRN9O2)
  
  !-> HO2 cont.
    L = L +   rct(49,k) * xnew(IBUTOLBO2)  &
    + rct(64,k) * xnew(OXNCATECO2)  &
    + rct(64,k) * xnew(OXNNCATCO2)  &
    + rct(64,k) * xnew(OXCATEC1O2)  &
    + rct(64,k) * xnew(NOXYOL1O2)  &
    + rct(64,k) * xnew(DM123O2)  &
    + rct(64,k) * xnew(NOXYOLO2)  &
    + rct(64,k) * xnew(OXYBIPERO2)  &
    + rct(67,k) * xnew(C5DICARBO2)  &
    + rct(67,k) * xnew(ISRO2)  &
    + rct(49,k) * xnew(MACRO2)  &
    + rct(45,k) * xnew(MACO3)  &
    + rct(49,k) * xnew(HMVKO2)  &
    + rct(67,k) * xnew(NISOPO2)  &
    + rct(45,k) * xnew(NC4CO3)  &
    + rct(45,k) * xnew(CO2H3CO3)  &
    + rct(86,k) * xnew(TERPO2)  &
    + rct(110,k)
  
  xnew(HO2) = (xold(HO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3O2

  P = rct(27,k) * xnew(CH4) * xnew(OH)  &
    + 0.6* rct(34,k) * xnew(CH3OOH) * xnew(OH)  &
    + rct(42,k) * xnew(CH3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CH3CO3) * xnew(NO3)  &
    + 0.454* rct(45,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.7* 1.0e-11 * xnew(RO2POOL) * xnew(CH3CO3)  &
    + 0.286* rct(45,k) * xnew(MACO3) * xnew(HO2)  &
    + 0.455* 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)  &
    + 0.65* rct(76,k) * xnew(MACO3) * xnew(NO)  &
    + 0.65* 4.0e-12 * xnew(MACO3) * xnew(NO3)  &
    + 0.65* 1.51e-11 * xnew(MACO2H) * xnew(OH)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + 0.325* rcphot(IDHCHO_H2,k) *0.04767 * xnew(MACR)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)
  
  L = 1.2e-12 * xnew(NO3)  &
    + rct(29,k) * xnew(NO)  &
    + rct(30,k) * xnew(RO2POOL)  &
    + rct(31,k) * xnew(RO2POOL)  &
    + rct(32,k) * xnew(HO2)
  
  xnew(CH3O2) = (xold(CH3O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5O2

  P = rct(37,k) * xnew(C2H6) * xnew(OH)  &
    + 0.31* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.264* rct(46,k) * xnew(C4H9O2) * xnew(NO)  &
    + 0.284* 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + 0.75* rct(45,k) * xnew(C2H5CO3) * xnew(HO2)  &
    + rct(50,k) * xnew(C2H5CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(C2H5CO3) * xnew(NO3)  &
    + 0.7* 1.0e-11 * xnew(RO2POOL) * xnew(C2H5CO3)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + 0.284* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + 0.348* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)
  
  L = rct(38,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(39,k) * xnew(HO2)  &
    + 3.10e-13 * xnew(RO2POOL)
  
  xnew(C2H5O2) = (xold(C2H5O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C4H9O2

  P = rct(48,k) * xnew(NC4H10) * xnew(OH)  &
    + 0.153* 2.4e-11 * xnew(BURO2H) * xnew(OH)
  
  L = rct(49,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)
  
  xnew(C4H9O2) = (xold(C4H9O2) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2

  P = rct(70,k) * xnew(C5H8) * xnew(OH)  &
    + 0.07* 6.0e-11 * xnew(ISRO2H) * xnew(OH)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(67,k) * xnew(HO2)  &
    + 1.45e-12 * xnew(RO2POOL)  &
    + rct(83,k)  &
    + rct(84,k)
  
  xnew(ISRO2) = (xold(ISRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> HCOCH2O2

  P = 0.05* rct(40,k) * xnew(CH3CHO) * xnew(OH)  &
    + 0.11* 3.3e-11 * xnew(HCOCH2OOH) * xnew(OH)  &
    + 0.3* rct(59,k) * xnew(C3H6) * xnew(O3)
  
  L = rct(47,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + 2.0e-12 * xnew(RO2POOL)
  
  xnew(HCOCH2O2) = (xold(HCOCH2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2

  P = 0.05* 3.22e-12 * xnew(C2H5OH) * xnew(OH)  &
    + rct(53,k) * xnew(C2H4) * xnew(OH)  &
    + 0.21* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(54,k) * xnew(HO2)  &
    + 1.9e-12 * xnew(RO2POOL)
  
  xnew(ETRO2) = (xold(ETRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> PRRO2

  P = rct(58,k) * xnew(C3H6) * xnew(OH)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(60,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(PRRO2) = (xold(PRRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> HMVKO2

  P = 0.74* 1.16e-11 * xnew(IEPOX) * xnew(OH)  &
    + rct(79,k) * xnew(MVK) * xnew(OH)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(49,k) * xnew(HO2)  &
    + 1.22e-12 * xnew(RO2POOL)
  
  xnew(HMVKO2) = (xold(HMVKO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYBIPERO2

  P = 0.346* rct(62,k) * xnew(BENZENE) * xnew(OH)  &
    + 0.71* rct(63,k) * xnew(TOLUENE) * xnew(OH)  &
    + 0.84* 1.36e-11 * xnew(OXYL) * xnew(OH)  &
    + 7.62e-11 * xnew(OXYBPEROOH) * xnew(OH)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(64,k) * xnew(HO2)
  
  xnew(OXYBIPERO2) = (xold(OXYBIPERO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACO3

  P = 0.45* rct(73,k) * xnew(MACR) * xnew(OH)  &
    + 3.4e-15 * xnew(MACR) * xnew(NO3)  &
    + rct(77,k) * xnew(MPAN)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.04767 * xnew(MACR)
  
  L = rct(45,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)  &
    + rct(76,k) * xnew(NO)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(43,k) * xnew(NO2)
  
  xnew(MACO3) = (xold(MACO3) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2

  P = 1.1e-12 * xnew(MEK) * xnew(OH)  &
    + 0.174* 2.1e-11 * xnew(MEKO2H) * xnew(OH)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(49,k) * xnew(HO2)  &
    + 1.04e-12 * xnew(RO2POOL)
  
  xnew(MEKO2) = (xold(MEKO2) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICARBO2

  P = 6.2e-11 * xnew(C5DICARB) * xnew(OH)  &
    + rct(65,k) * xnew(C5DICAROOH) * xnew(OH)  &
    + rcphot(IDNO2,k) *0.2 * xnew(C5DICARB)
  
  L = rct(46,k) * xnew(NO)  &
    + rct(66,k) * xnew(NO3)  &
    + rct(67,k) * xnew(HO2)
  
  xnew(C5DICARBO2) = (xold(C5DICARBO2) + dt2 * P) / (1.0 + dt2 * L)


!-> DM123O2

  P = 2.86e-13 * xnew(OXY1O) * xnew(O3)  &
    + 4.65e-11 * xnew(DM123OOH) * xnew(OH)
  
  L = rct(64,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)
  
  xnew(DM123O2) = (xold(DM123O2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXCATEC1O2

  P = 2.86e-13 * xnew(OXCATEC1O) * xnew(O3)  &
    + 2.05e-10 * xnew(OCATEC1OOH) * xnew(OH)
  
  L = rct(64,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(OXCATEC1O2) = (xold(OXCATEC1O2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH2O2CO3H

  P = 0.5* 1.1e-11 * xnew(CH3CO3H) * xnew(OH)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)
  
  xnew(CH2O2CO3H) = (xold(CH2O2CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> MACRO2

  P = 0.15* 1.16e-11 * xnew(IEPOX) * xnew(OH)  &
    + 0.55* rct(73,k) * xnew(MACR) * xnew(OH)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(49,k) * xnew(HO2)  &
    + rct(75,k)
  
  xnew(MACRO2) = (xold(MACRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> HOCH2CO3

  P = rct(46,k) * xnew(CH2O2CO3H) * xnew(NO)  &
    + 2.3e-12 * xnew(CH2O2CO3H) * xnew(NO3)  &
    + 0.8* 1.00e-11 * xnew(HOCH2CHO) * xnew(OH)  &
    + rct(41,k) * xnew(HOCH2CHO) * xnew(NO3)  &
    + rct(44,k) * xnew(PHAN)  &
    + rcphot(IDCH3COY,k) * xnew(BIACETOH)
  
  L = rct(45,k) * xnew(HO2)  &
    + rct(42,k) * xnew(NO)  &
    + rct(43,k) * xnew(NO2)  &
    + 4.0e-12 * xnew(NO3)  &
    + 1.00e-11 * xnew(RO2POOL)
  
  xnew(HOCH2CO3) = (xold(HOCH2CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3

  P = 0.95* rct(40,k) * xnew(CH3CHO) * xnew(OH)  &
    + rct(41,k) * xnew(CH3CHO) * xnew(NO3)  &
    + rct(44,k) * xnew(PAN)  &
    + 0.5* 1.1e-11 * xnew(CH3CO3H) * xnew(OH)  &
    + 0.854* rct(46,k) * xnew(MEKO2) * xnew(NO)  &
    + 0.854* 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + 0.6405* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + 8.3337* rcphot(IDMEK,k) * xnew(NOA)  &
    + rct(69,k) * xnew(MGLYOX) * xnew(OH)  &
    + 0.1535* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.0984* rct(74,k) * xnew(MACR) * xnew(O3)  &
    + 0.154* rct(45,k) * xnew(MACO3) * xnew(HO2)  &
    + 0.245* 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)  &
    + 0.35* rct(76,k) * xnew(MACO3) * xnew(NO)  &
    + 0.35* 4.0e-12 * xnew(MACO3) * xnew(NO3)  &
    + 8.2e-18 * xnew(MPAN) * xnew(O3)  &
    + 0.28* rct(80,k) * xnew(MVK) * xnew(O3)  &
    + 0.67536* rct(46,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.7* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)
  
  !-> CH3CO3 cont.
    P = P +   0.42* 1.22e-12 * xnew(RO2POOL) * xnew(HMVKO2)  &
    + 1.23e-11 * xnew(CO23C3CHO) * xnew(OH)  &
    + rct(81,k) * xnew(CO23C3CHO) * xnew(NO3)  &
    + 0.35* 1.51e-11 * xnew(MACO2H) * xnew(OH)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + 0.35* rcphot(IDRCOCHO,k) * xnew(CH3COCO2H)  &
    + 2* rcphot(IDCH3COY,k) * xnew(BIACET)  &
    + rcphot(IDCH3COY,k) * xnew(BIACETOH)  &
    + rcphot(IDCH3COY,k) * xnew(CO23C3CHO)  &
    + rcphot(IDRCOCHO,k) * xnew(CO23C3CHO)  &
    + 0.175* rcphot(IDHCHO_H2,k) *0.04767 * xnew(MACR)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + 0.32* rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + 0.854* rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + 0.854* rcphot(IDMEK,k) * xnew(MEKO2H)  &
    + 0.7* rcphot(IDCH3O2H,k) * xnew(MVKO2H)  &
    + rcphot(IDMEK,k) *13. * xnew(MVKNO3)  &
    + rcphot(IDMEK,k) *83. * xnew(CO2N3CHO)
  
  L = rct(42,k) * xnew(NO)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(43,k) * xnew(NO2)  &
    + rct(45,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(CH3CO3) = (xold(CH3CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5CO3

  P = 0.146* rct(46,k) * xnew(MEKO2) * xnew(NO)  &
    + 0.146* 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + 0.1095* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + rct(51,k) * xnew(PPN)  &
    + rct(52,k) * xnew(EGLYOX) * xnew(NO3)  &
    + 1.39e-11 * xnew(EGLYOX) * xnew(OH)  &
    + rcphot(IDRCOCHO,k) * xnew(EGLYOX)  &
    + 0.146* rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + 0.146* rcphot(IDMEK,k) * xnew(MEKO2H)
  
  L = rct(45,k) * xnew(HO2)  &
    + rct(50,k) * xnew(NO)  &
    + rct(43,k) * xnew(NO2)  &
    + 4.0e-12 * xnew(NO3)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(C2H5CO3) = (xold(C2H5CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> CO2H3CO3

  P = 2.45e-11 * xnew(CO2H3CHO) * xnew(OH)  &
    + rct(81,k) * xnew(CO2H3CHO) * xnew(NO3)  &
    + 7.34e-12 * xnew(CO2H3CO3H) * xnew(OH)  &
    + rct(44,k) * xnew(C4PAN6)
  
  L = rct(45,k) * xnew(HO2)  &
    + rct(42,k) * xnew(NO)  &
    + 4.0e-12 * xnew(NO3)  &
    + 1.0e-11 * xnew(RO2POOL)  &
    + rct(43,k) * xnew(NO2)
  
  xnew(CO2H3CO3) = (xold(CO2H3CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> IBUTOLBO2

  P = 3.0e-13 * xnew(BZFUONE) * xnew(NO3)  &
    + 4.45e-11 * xnew(BZFUONE) * xnew(OH)
  
  L = rct(49,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)
  
  xnew(IBUTOLBO2) = (xold(IBUTOLBO2) + dt2 * P) / (1.0 + dt2 * L)


!-> HOCH2CH2O

  P = rct(46,k) * xnew(ETRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(ETRO2) * xnew(NO3)  &
    + 0.6* 1.9e-12 * xnew(RO2POOL) * xnew(ETRO2)  &
    + rcphot(IDCH3O2H,k) * xnew(ETRO2H)
  
  L = rct(55,k)  &
    + rct(56,k)
  
  xnew(HOCH2CH2O) = (xold(HOCH2CH2O) + dt2 * P) / (1.0 + dt2 * L)


!-> OXY1O

  P = 0.39* 2.1e-11 * xnew(OXYLOL) * xnew(NO3)  &
    + 0.17* 6.1e-11 * xnew(OXYLOL) * xnew(OH)  &
    + rct(46,k) * xnew(DM123O2) * xnew(NO)  &
    + 2.3e-12 * xnew(DM123O2) * xnew(NO3)  &
    + rcphot(IDCH3O2H,k) * xnew(DM123OOH)
  
  L = 2.08e-12 * xnew(NO2)  &
    + 2.86e-13 * xnew(O3)
  
  xnew(OXY1O) = (xold(OXY1O) + dt2 * P) / (1.0 + dt2 * L)


!-> OXCATEC1O

  P = 2.01e-10 * xnew(OXYCATECH) * xnew(NO3)  &
    + 2.05e-10 * xnew(OXYCATECH) * xnew(OH)  &
    + rct(46,k) * xnew(OXCATEC1O2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXCATEC1O2) * xnew(NO3)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(OXCATEC1O2)  &
    + rcphot(IDCH3O2H,k) * xnew(OCATEC1OOH)
  
  L = 2.08e-12 * xnew(NO2)  &
    + 2.86e-13 * xnew(O3)
  
  xnew(OXCATEC1O) = (xold(OXCATEC1O) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPO2

  P = rct(85,k) * xnew(APINENE) * xnew(OH)  &
    + rct(87,k) * xnew(TERPOOH) * xnew(OH)  &
    + rct(89,k) * xnew(APINENE) * xnew(NO3)
  
  L = rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + rct(86,k) * xnew(HO2)
  
  xnew(TERPO2) = (xold(TERPO2) + dt2 * P) / (1.0 + dt2 * L)


!-> RO2POOL

  ! P = 0.0
  
  ! L = 0.0
  
  xnew(RO2POOL) = sum(xnew(RO2_GROUP))


!-> NRN9O2

  P = rct(61,k) * xnew(C3H6) * xnew(NO3)
  
  L = rct(60,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)
  
  xnew(NRN9O2) = (xold(NRN9O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NISOPO2

  P = rct(72,k) * xnew(C5H8) * xnew(NO3)
  
  L = rct(67,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + 1.3e-12 * xnew(RO2POOL)
  
  xnew(NISOPO2) = (xold(NISOPO2) + dt2 * P) / (1.0 + dt2 * L)


!-> NC4CO3

  P = 4.16e-11 * xnew(ISON2) * xnew(OH)  &
    + rct(82,k) * xnew(ISON2) * xnew(NO3)  &
    + rct(44,k) * xnew(C5PAN18)
  
  L = rct(45,k) * xnew(HO2)  &
    + rct(42,k) * xnew(NO)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(43,k) * xnew(NO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(NC4CO3) = (xold(NC4CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> NOXYOLO2

  P = 0.61* 2.1e-11 * xnew(OXYLOL) * xnew(NO3)  &
    + 1.43e-10 * xnew(NOXYOLOOH) * xnew(OH)
  
  L = rct(64,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL)
  
  xnew(NOXYOLO2) = (xold(NOXYOLO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXNCATECO2

  P = 6.83e-12 * xnew(OXNCATECH) * xnew(OH)  &
    + rct(65,k) * xnew(OXNCATCOOH) * xnew(OH)
  
  L = rct(64,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL)
  
  xnew(OXNCATECO2) = (xold(OXNCATECO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXNNCATCO2

  P = 5.03e-12 * xnew(OXNCATECH) * xnew(NO3)  &
    + rct(65,k) * xnew(OXNNCATOOH) * xnew(OH)
  
  L = rct(64,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL)
  
  xnew(OXNNCATCO2) = (xold(OXNNCATCO2) + dt2 * P) / (1.0 + dt2 * L)


!-> NOXYOL1O2

  P = 2.86e-13 * xnew(NOXYOL1O) * xnew(O3)  &
    + 1.53e-12 * xnew(NOXYOL1OOH) * xnew(OH)
  
  L = rct(64,k) * xnew(HO2)  &
    + rct(46,k) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3)
  
  xnew(NOXYOL1O2) = (xold(NOXYOL1O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NOXYOL1O

  P = 3.13e-13 * xnew(OXY1OHNO2) * xnew(NO3)  &
    + 1.53e-12 * xnew(OXY1OHNO2) * xnew(OH)  &
    + rct(46,k) * xnew(NOXYOL1O2) * xnew(NO)  &
    + 2.3e-12 * xnew(NOXYOL1O2) * xnew(NO3)  &
    + rcphot(IDCH3O2H,k) * xnew(NOXYOL1OOH)
  
  L = 2.08e-12 * xnew(NO2)  &
    + 2.86e-13 * xnew(O3)
  
  xnew(NOXYOL1O) = (xold(NOXYOL1O) + dt2 * P) / (1.0 + dt2 * L)


!-> NMVOC

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NMVOC! xnew(NMVOC) = max(0.0, xold(NMVOC))


!-> O3

  P = rct(1,k) * xnew(OP)  &
    + 0.162* rct(45,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.25* rct(45,k) * xnew(C2H5CO3) * xnew(HO2)  &
    + 0.25* rct(45,k) * xnew(HOCH2CO3) * xnew(HO2)  &
    + 0.56* rct(45,k) * xnew(MACO3) * xnew(HO2)
  
  L = rct(11,k) * xnew(NO)  &
    + rct(12,k) * xnew(NO2)  &
    + rct(13,k) * xnew(OH)  &
    + rct(14,k) * xnew(HO2)  &
    + rct(57,k) * xnew(C2H4)  &
    + rct(59,k) * xnew(C3H6)  &
    + 2.86e-13 * xnew(OXCATEC1O)  &
    + 2.86e-13 * xnew(OXY1O)  &
    + 2.86e-13 * xnew(NOXYOL1O)  &
    + rct(71,k) * xnew(C5H8)  &
    + rct(74,k) * xnew(MACR)  &
    + 8.2e-18 * xnew(MPAN)  &
    + rct(80,k) * xnew(MVK)  &
    + 2.4e-17 * xnew(ISON2)  &
    + rct(88,k) * xnew(APINENE)  &
    + rcphot(IDO3_O1D,k)  &
    + rcphot(IDO3_O3P,k)  &
    + AQRCK(ICLRC2,K) * xnew(SO2)
  
  xnew(O3) = (xold(O3) + dt2 * P) / (1.0 + dt2 * L)


!-> NO

  P = rct(17,k) * xnew(NO2) * xnew(NO3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + 0.127* rcphot(IDNO3,k) * xnew(NO3)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcemis(NO,k)
  
  L = rct(10,k) * xnew(OP)  &
    + rct(11,k) * xnew(O3)  &
    + rct(15,k) * xnew(NO3)  &
    + rct(16,k) * xnew(HO2)  &
    + rct(22,k) * xnew(OH)  &
    + rct(29,k) * xnew(CH3O2)  &
    + rct(38,k) * xnew(C2H5O2)  &
    + rct(42,k) * xnew(CH3CO3)  &
    + rct(46,k) * xnew(CH2O2CO3H)  &
    + rct(46,k) * xnew(HCOCH2O2)  &
    + rct(46,k) * xnew(C4H9O2)  &
    + rct(46,k) * xnew(MEKO2)  &
    + rct(50,k) * xnew(C2H5CO3)  &
    + rct(46,k) * xnew(ETRO2)  &
    + rct(42,k) * xnew(HOCH2CO3)  &
    + rct(46,k) * xnew(PRRO2)  &
    + rct(46,k) * xnew(NRN9O2)  &
    + rct(46,k) * xnew(IBUTOLBO2)  &
    + rct(46,k) * xnew(OXNCATECO2)  &
    + rct(46,k) * xnew(OXNNCATCO2)
  
  !-> NO cont.
    L = L +   rct(46,k) * xnew(OXCATEC1O2)  &
    + rct(46,k) * xnew(NOXYOL1O2)  &
    + rct(46,k) * xnew(DM123O2)  &
    + rct(46,k) * xnew(NOXYOLO2)  &
    + rct(46,k) * xnew(OXYBIPERO2)  &
    + rct(46,k) * xnew(C5DICARBO2)  &
    + rct(46,k) * xnew(ISRO2)  &
    + rct(46,k) * xnew(MACRO2)  &
    + rct(76,k) * xnew(MACO3)  &
    + rct(46,k) * xnew(HMVKO2)  &
    + rct(46,k) * xnew(NISOPO2)  &
    + rct(42,k) * xnew(NC4CO3)  &
    + rct(42,k) * xnew(CO2H3CO3)  &
    + rct(46,k) * xnew(TERPO2)
  
  xnew(NO) = (xold(NO) + dt2 * P) / (1.0 + dt2 * L)


!-> NO2

  P = rct(10,k) * xnew(OP) * xnew(NO)  &
    + rct(11,k) * xnew(O3) * xnew(NO)  &
    + 2* rct(15,k) * xnew(NO) * xnew(NO3)  &
    + rct(16,k) * xnew(NO) * xnew(HO2)  &
    + rct(17,k) * xnew(NO2) * xnew(NO3)  &
    + rct(19,k) * xnew(N2O5)  &
    + rct(23,k) * xnew(OH) * xnew(HONO)  &
    + 2.0e-11 * xnew(OH) * xnew(NO3)  &
    + 4.0e-12 * xnew(HO2) * xnew(NO3)  &
    + rct(25,k) * xnew(HO2NO2)  &
    + rct(26,k) * xnew(OH) * xnew(HO2NO2)  &
    + 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(29,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(38,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + rct(42,k) * xnew(CH3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CH3CO3) * xnew(NO3)  &
    + rct(44,k) * xnew(PAN)  &
    + rct(46,k) * xnew(CH2O2CO3H) * xnew(NO)  &
    + 2.3e-12 * xnew(CH2O2CO3H) * xnew(NO3)
  
  !-> NO2 cont.
    P = P +   rct(46,k) * xnew(HCOCH2O2) * xnew(NO)  &
    + 2.3e-12 * xnew(HCOCH2O2) * xnew(NO3)  &
    + 0.917* rct(46,k) * xnew(C4H9O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + 8.6e-13 * xnew(SC4H9NO3) * xnew(OH)  &
    + rct(46,k) * xnew(MEKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + rct(50,k) * xnew(C2H5CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(C2H5CO3) * xnew(NO3)  &
    + 1.27e-12 * xnew(PPN) * xnew(OH)  &
    + rct(51,k) * xnew(PPN)  &
    + rct(46,k) * xnew(ETRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(ETRO2) * xnew(NO3)  &
    + rct(42,k) * xnew(HOCH2CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(HOCH2CO3) * xnew(NO3)  &
    + 1.12e-12 * xnew(PHAN) * xnew(OH)  &
    + rct(44,k) * xnew(PHAN)  &
    + rct(46,k) * xnew(PRRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(PRRO2) * xnew(NO3)  &
    + 1.13* rct(46,k) * xnew(NRN9O2) * xnew(NO)
  
  !-> NO2 cont.
    P = P +   1.13* 2.3e-12 * xnew(NRN9O2) * xnew(NO3)  &
    + 8.3337* rcphot(IDMEK,k) * xnew(NOA)  &
    + 6.70e-13 * xnew(NOA) * xnew(OH)  &
    + rct(46,k) * xnew(IBUTOLBO2) * xnew(NO)  &
    + 2.3e-12 * xnew(IBUTOLBO2) * xnew(NO3)  &
    + rct(46,k) * xnew(OXNCATECO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXNCATECO2) * xnew(NO3)  &
    + rct(65,k) * xnew(NC4MDCO2H) * xnew(OH)  &
    + 2* rct(46,k) * xnew(OXNNCATCO2) * xnew(NO)  &
    + 2* 2.3e-12 * xnew(OXNNCATCO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(OXNNCATCO2)  &
    + rct(46,k) * xnew(OXCATEC1O2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXCATEC1O2) * xnew(NO3)  &
    + rct(46,k) * xnew(NOXYOL1O2) * xnew(NO)  &
    + 2.3e-12 * xnew(NOXYOL1O2) * xnew(NO3)  &
    + rct(46,k) * xnew(DM123O2) * xnew(NO)  &
    + 2.3e-12 * xnew(DM123O2) * xnew(NO3)  &
    + 2* rct(46,k) * xnew(NOXYOLO2) * xnew(NO)  &
    + 2* 2.3e-12 * xnew(NOXYOLO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(NOXYOLO2)
  
  !-> NO2 cont.
    P = P +   0.862* rct(46,k) * xnew(OXYBIPERO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXYBIPERO2) * xnew(NO3)  &
    + 3.0e-13 * xnew(BZFUONE) * xnew(NO3)  &
    + rct(46,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(66,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.9* rct(46,k) * xnew(ISRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3) * xnew(ISRO2)  &
    + 0.4* 3.25e-11 * xnew(ISON) * xnew(OH)  &
    + 1.55e-12 * xnew(NALD) * xnew(OH)  &
    + rct(46,k) * xnew(MACRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(MACRO2) * xnew(NO3)  &
    + rct(76,k) * xnew(MACO3) * xnew(NO)  &
    + 4.0e-12 * xnew(MACO3) * xnew(NO3)  &
    + rct(77,k) * xnew(MPAN)  &
    + 0.9648* rct(46,k) * xnew(HMVKO2) * xnew(NO)  &
    + 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.33* 1.33e-12 * xnew(MVKNO3) * xnew(OH)  &
    + rct(46,k) * xnew(NISOPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NISOPO2) * xnew(NO3)  &
    + 0.445* 2.4e-17 * xnew(ISON2) * xnew(O3)
  
  !-> NO2 cont.
    P = P +   rct(42,k) * xnew(NC4CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(NC4CO3) * xnew(NO3)  &
    + rct(44,k) * xnew(C5PAN18)  &
    + 2.16e-11 * xnew(C5PAN18) * xnew(OH)  &
    + rct(42,k) * xnew(CO2H3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CO2H3CO3) * xnew(NO3)  &
    + 3.74e-12 * xnew(C4PAN6) * xnew(OH)  &
    + rct(44,k) * xnew(C4PAN6)  &
    + rct(46,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rct(89,k) * xnew(APINENE) * xnew(NO3)  &
    + 0.873* rcphot(IDNO3,k) * xnew(NO3)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.667* rcphot(IDHO2NO2,k) * xnew(HO2NO2)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rcphot(IDCH3O2H,k) * xnew(OXNNCATOOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NOXYOLOOH)  &
    + rcphot(IDMEK,k) *0.78 * xnew(NOXYOLOOH)  &
    + rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rcphot(IDMEK,k) *13. * xnew(MVKNO3)
  
  !-> NO2 cont.
    P = P +   rcphot(IDMEK,k) *83. * xnew(CO2N3CHO)  &
    + rcphot(IDHCHO_H,k) *0.1814 * xnew(ISON)  &
    + rcphot(IDMEK,k) *133.34 * xnew(ISON2)  &
    + rcemis(NO2,k)
  
  L = rct(12,k) * xnew(O3)  &
    + rct(17,k) * xnew(NO3)  &
    + rct(18,k) * xnew(NO3)  &
    + rct(20,k) * xnew(OH)  &
    + rct(24,k) * xnew(HO2)  &
    + rct(43,k) * xnew(CH3CO3)  &
    + rct(43,k) * xnew(C2H5CO3)  &
    + rct(43,k) * xnew(HOCH2CO3)  &
    + 2.08e-12 * xnew(OXCATEC1O)  &
    + 2.08e-12 * xnew(OXY1O)  &
    + 2.08e-12 * xnew(NOXYOL1O)  &
    + rct(43,k) * xnew(MACO3)  &
    + rct(43,k) * xnew(NC4CO3)  &
    + rct(43,k) * xnew(CO2H3CO3)  &
    + rcphot(IDNO2,k)
  
  xnew(NO2) = (xold(NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> NO3

  P = rct(12,k) * xnew(O3) * xnew(NO2)  &
    + rct(19,k) * xnew(N2O5)  &
    + rct(21,k) * xnew(OH) * xnew(HNO3)  &
    + 2.9e-11 * xnew(MPAN) * xnew(OH)  &
    + 8.2e-18 * xnew(MPAN) * xnew(O3)  &
    + 0.333* rcphot(IDHO2NO2,k) * xnew(HO2NO2)
  
  L = rct(15,k) * xnew(NO)  &
    + rct(17,k) * xnew(NO2)  &
    + rct(18,k) * xnew(NO2)  &
    + 2.0e-11 * xnew(OH)  &
    + 4.0e-12 * xnew(HO2)  &
    + 1.2e-12 * xnew(CH3O2)  &
    + rct(36,k) * xnew(HCHO)  &
    + 2.3e-12 * xnew(C2H5O2)  &
    + rct(41,k) * xnew(CH3CHO)  &
    + 4.0e-12 * xnew(CH3CO3)  &
    + 2.3e-12 * xnew(CH2O2CO3H)  &
    + 2.3e-12 * xnew(HCOCH2O2)  &
    + 2.3e-12 * xnew(C4H9O2)  &
    + 2.3e-12 * xnew(MEKO2)  &
    + 4.0e-12 * xnew(C2H5CO3)  &
    + rct(52,k) * xnew(EGLYOX)  &
    + 2.3e-12 * xnew(ETRO2)  &
    + rct(41,k) * xnew(HOCH2CHO)  &
    + 4.0e-12 * xnew(HOCH2CO3)  &
    + 2.3e-12 * xnew(PRRO2)
  
  !-> NO3 cont.
    L = L +   rct(61,k) * xnew(C3H6)  &
    + 2.3e-12 * xnew(NRN9O2)  &
    + 2.3e-12 * xnew(IBUTOLBO2)  &
    + 2.1e-11 * xnew(OXYLOL)  &
    + 2.01e-10 * xnew(OXYCATECH)  &
    + 5.03e-12 * xnew(OXNCATECH)  &
    + 2.3e-12 * xnew(OXNCATECO2)  &
    + 2.3e-12 * xnew(OXNNCATCO2)  &
    + 2.3e-12 * xnew(OXCATEC1O2)  &
    + 3.13e-13 * xnew(OXY1OHNO2)  &
    + 2.3e-12 * xnew(NOXYOL1O2)  &
    + 2.3e-12 * xnew(DM123O2)  &
    + 2.3e-12 * xnew(NOXYOLO2)  &
    + 2.3e-12 * xnew(OXYBIPERO2)  &
    + 3.0e-13 * xnew(BZFUONE)  &
    + rct(66,k) * xnew(C5DICARBO2)  &
    + rct(72,k) * xnew(C5H8)  &
    + 2.3e-12 * xnew(ISRO2)  &
    + 3.4e-15 * xnew(MACR)  &
    + 2.3e-12 * xnew(MACRO2)
  
  !-> NO3 cont.
    L = L +   4.0e-12 * xnew(MACO3)  &
    + 2.3e-12 * xnew(HMVKO2)  &
    + rct(81,k) * xnew(CO23C3CHO)  &
    + rct(81,k) * xnew(CO2H3CHO)  &
    + 2.3e-12 * xnew(NISOPO2)  &
    + rct(82,k) * xnew(ISON2)  &
    + 4.0e-12 * xnew(NC4CO3)  &
    + 4.0e-12 * xnew(CO2H3CO3)  &
    + 2.3e-12 * xnew(TERPO2)  &
    + rct(89,k) * xnew(APINENE)  &
    + rcphot(IDNO3,k)  &
    + rct(107,k)
  
  xnew(NO3) = (xold(NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> N2O5

  P = rct(18,k) * xnew(NO2) * xnew(NO3)
  
  L = rct(19,k)  &
    + rct(106,k)
  
  xnew(N2O5) = (xold(N2O5) + dt2 * P) / (1.0 + dt2 * L)


!-> H2

  P = 0.18* rct(57,k) * xnew(C2H4) * xnew(O3)  &
    + 0.0684* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + rcphot(IDHCHO_H2,k) * xnew(HCHO)
  
  L = rct(7,k) * xnew(OH)
  
  xnew(H2) = (xold(H2) + dt2 * P) / (1.0 + dt2 * L)


!-> H2O2

  P = rct(8,k) * xnew(HO2) * xnew(HO2)  &
    + rct(9,k) * xnew(HO2) * xnew(HO2)  &
    + 0.1575* rct(57,k) * xnew(C2H4) * xnew(O3)  &
    + 0.169* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + 0.28* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.1437* rct(74,k) * xnew(MACR) * xnew(O3)  &
    + 0.09* rct(80,k) * xnew(MVK) * xnew(O3)  &
    + 0.0756* 2.4e-17 * xnew(ISON2) * xnew(O3)  &
    + 0.5* rct(110,k) * xnew(HO2)
  
  L = rct(6,k) * xnew(OH)  &
    + rcphot(IDH2O2,k)  &
    + AQRCK(ICLRC1,K) * xnew(SO2)  &
    + 1.0e-5
  
  xnew(H2O2) = (xold(H2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> HONO

  P = rct(22,k) * xnew(OH) * xnew(NO)
  
  L = rct(23,k) * xnew(OH)  &
    + rcphot(IDHONO,k)
  
  xnew(HONO) = (xold(HONO) + dt2 * P) / (1.0 + dt2 * L)


!-> HNO3

  P = rct(20,k) * xnew(NO2) * xnew(OH)  &
    + rct(36,k) * xnew(HCHO) * xnew(NO3)  &
    + rct(41,k) * xnew(CH3CHO) * xnew(NO3)  &
    + rct(52,k) * xnew(EGLYOX) * xnew(NO3)  &
    + rct(41,k) * xnew(HOCH2CHO) * xnew(NO3)  &
    + 0.39* 2.1e-11 * xnew(OXYLOL) * xnew(NO3)  &
    + 2.01e-10 * xnew(OXYCATECH) * xnew(NO3)  &
    + 3.13e-13 * xnew(OXY1OHNO2) * xnew(NO3)  &
    + 3.4e-15 * xnew(MACR) * xnew(NO3)  &
    + rct(81,k) * xnew(CO23C3CHO) * xnew(NO3)  &
    + rct(81,k) * xnew(CO2H3CHO) * xnew(NO3)  &
    + rct(82,k) * xnew(ISON2) * xnew(NO3)  &
    + 2* rct(106,k) * xnew(N2O5)  &
    + rct(107,k) * xnew(NO3)  &
    + rct(20,k) * xnew(shipNOx) * xnew(OH)  &
    + 3.2e-5 * xnew(shipNOx)
  
  L = rct(21,k) * xnew(OH)  &
    + rcphot(IDHNO3,k)  &
    + rct(108,k)  &
    + rct(109,k)
  
  xnew(HNO3) = (xold(HNO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2NO2

  P = rct(24,k) * xnew(HO2) * xnew(NO2)
  
  L = rct(25,k)  &
    + rct(26,k) * xnew(OH)  &
    + rcphot(IDHO2NO2,k)
  
  xnew(HO2NO2) = (xold(HO2NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CO

  P = rct(35,k) * xnew(HCHO) * xnew(OH)  &
    + rct(36,k) * xnew(HCHO) * xnew(NO3)  &
    + rct(46,k) * xnew(HCOCH2O2) * xnew(NO)  &
    + 2.3e-12 * xnew(HCOCH2O2) * xnew(NO3)  &
    + 0.6* 2.0e-12 * xnew(RO2POOL) * xnew(HCOCH2O2)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + 1.27e-12 * xnew(PPN) * xnew(OH)  &
    + rct(52,k) * xnew(EGLYOX) * xnew(NO3)  &
    + 1.39e-11 * xnew(EGLYOX) * xnew(OH)  &
    + 1.12e-12 * xnew(PHAN) * xnew(OH)  &
    + 0.35* rct(57,k) * xnew(C2H4) * xnew(O3)  &
    + 0.183* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + 1.417* rct(68,k) * xnew(GLYOX) * xnew(OH)  &
    + rct(69,k) * xnew(MGLYOX) * xnew(OH)  &
    + 0.364* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 1.55e-12 * xnew(NALD) * xnew(OH)  &
    + 0.6528* rct(74,k) * xnew(MACR) * xnew(O3)  &
    + 0.85* rct(46,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.85* 2.3e-12 * xnew(MACRO2) * xnew(NO3)  &
    + rct(75,k) * xnew(MACRO2)
  
  !-> CO cont.
    P = P +   0.286* rct(45,k) * xnew(MACO3) * xnew(HO2)  &
    + 0.455* 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)  &
    + 0.65* rct(76,k) * xnew(MACO3) * xnew(NO)  &
    + 0.65* 4.0e-12 * xnew(MACO3) * xnew(NO3)  &
    + 2.9e-11 * xnew(MPAN) * xnew(OH)  &
    + 3.77e-11 * xnew(MACROOH) * xnew(OH)  &
    + 0.56* rct(80,k) * xnew(MVK) * xnew(O3)  &
    + 2* 1.23e-11 * xnew(CO23C3CHO) * xnew(OH)  &
    + 2* rct(81,k) * xnew(CO23C3CHO) * xnew(NO3)  &
    + 0.89* 2.4e-17 * xnew(ISON2) * xnew(O3)  &
    + rct(45,k) * xnew(NC4CO3) * xnew(HO2)  &
    + rct(42,k) * xnew(NC4CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(NC4CO3) * xnew(NO3)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(NC4CO3)  &
    + 2* 2.16e-11 * xnew(C5PAN18) * xnew(OH)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + 3.74e-12 * xnew(C4PAN6) * xnew(OH)  &
    + 0.65* 1.51e-11 * xnew(MACO2H) * xnew(OH)  &
    + rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDHCHO_H2,k) * xnew(HCHO)
  
  !-> CO cont.
    P = P +   rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HOCH2CHO)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(CO2H3CHO)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + 1.9* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(EGLYOX)  &
    + 0.05* rcphot(IDRCOCHO,k) * xnew(CH3COCO2H)  &
    + 1.41* rcphot(IDCH3COY,k) * xnew(CO23C3CHO)  &
    + 2* rcphot(IDRCOCHO,k) * xnew(CO23C3CHO)  &
    + 0.825* rcphot(IDHCHO_H2,k) *0.04767 * xnew(MACR)  &
    + rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + 1.68* rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + rcphot(IDCH3O2H,k) * xnew(HCOCH2OOH)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HCOCH2OOH)  &
    + rcphot(IDHCHO_H,k) *1.87 * xnew(MACROOH)  &
    + 2* rcphot(IDMEK,k) *133.34 * xnew(ISON2)  &
    + rcemis(CO,k)
  
  L = rct(28,k) * xnew(OH)
  
  xnew(CO) = (xold(CO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH4

  P = 0.04* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + 0.11* 1.16e-11 * xnew(IEPOX) * xnew(OH)
  
  L = rct(27,k) * xnew(OH)
  
  xnew(CH4) = (xold(CH4) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H6

  P = rct(46,k) * xnew(IBUTOLBO2) * xnew(NO)  &
    + 2.3e-12 * xnew(IBUTOLBO2) * xnew(NO3)  &
    + 0.19* rct(46,k) * xnew(TERPO2) * xnew(NO)  &
    + 0.19* 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rcemis(C2H6,k)
  
  L = rct(37,k) * xnew(OH)
  
  xnew(C2H6) = (xold(C2H6) + dt2 * P) / (1.0 + dt2 * L)


!-> NC4H10

  P = rcemis(NC4H10,k)
  
  L = rct(48,k) * xnew(OH)
  
  xnew(NC4H10) = (xold(NC4H10) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H4

  P = rcemis(C2H4,k)
  
  L = rct(53,k) * xnew(OH)  &
    + rct(57,k) * xnew(O3)
  
  xnew(C2H4) = (xold(C2H4) + dt2 * P) / (1.0 + dt2 * L)


!-> C3H6

  P = rct(46,k) * xnew(NOXYOLO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NOXYOLO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(NOXYOLO2)  &
    + 0.0475* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + rcphot(IDCH3O2H,k) * xnew(NOXYOLOOH)  &
    + rcphot(IDMEK,k) *0.78 * xnew(NOXYOLOOH)  &
    + rcemis(C3H6,k)
  
  L = rct(58,k) * xnew(OH)  &
    + rct(59,k) * xnew(O3)  &
    + rct(61,k) * xnew(NO3)
  
  xnew(C3H6) = (xold(C3H6) + dt2 * P) / (1.0 + dt2 * L)


!-> BENZENE

  P = rcemis(BENZENE,k)
  
  L = rct(62,k) * xnew(OH)
  
  xnew(BENZENE) = (xold(BENZENE) + dt2 * P) / (1.0 + dt2 * L)


!-> TOLUENE

  P = rcemis(TOLUENE,k)
  
  L = rct(63,k) * xnew(OH)
  
  xnew(TOLUENE) = (xold(TOLUENE) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYL

  P = rcemis(OXYL,k)
  
  L = 1.36e-11 * xnew(OH)
  
  xnew(OXYL) = (xold(OXYL) + dt2 * P) / (1.0 + dt2 * L)


!-> C5H8

  P = rcemis(C5H8,k)  &
    + RCBIO(1,K)
  
  L = rct(70,k) * xnew(OH)  &
    + rct(71,k) * xnew(O3)  &
    + rct(72,k) * xnew(NO3)
  
  xnew(C5H8) = (xold(C5H8) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OH

  P = 0.5* rct(31,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.05* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + rct(46,k) * xnew(IBUTOLBO2) * xnew(NO)  &
    + 2.3e-12 * xnew(IBUTOLBO2) * xnew(NO3)  &
    + 0.89* 1.16e-11 * xnew(IEPOX) * xnew(OH)  &
    + rcemis(CH3OH,k)
  
  L = rct(33,k) * xnew(OH)
  
  xnew(CH3OH) = (xold(CH3OH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OH

  P = 0.2* 3.10e-13 * xnew(RO2POOL) * xnew(C2H5O2)  &
    + rct(49,k) * xnew(IBUTOLBO2) * xnew(HO2)  &
    + rcemis(C2H5OH,k)
  
  L = 3.22e-12 * xnew(OH)
  
  xnew(C2H5OH) = (xold(C2H5OH) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYLOL

  P = 0.39* rct(62,k) * xnew(BENZENE) * xnew(OH)  &
    + 0.16* rct(63,k) * xnew(TOLUENE) * xnew(OH)  &
    + 0.16* 1.36e-11 * xnew(OXYL) * xnew(OH)
  
  L = 2.1e-11 * xnew(NO3)  &
    + 6.1e-11 * xnew(OH)
  
  xnew(OXYLOL) = (xold(OXYLOL) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYCATECH

  P = 0.83* 6.1e-11 * xnew(OXYLOL) * xnew(OH)
  
  L = 2.01e-10 * xnew(NO3)  &
    + 2.05e-10 * xnew(OH)
  
  xnew(OXYCATECH) = (xold(OXYCATECH) + dt2 * P) / (1.0 + dt2 * L)


!-> HCHO

  P = 1.2e-12 * xnew(CH3O2) * xnew(NO3)  &
    + rct(29,k) * xnew(CH3O2) * xnew(NO)  &
    + rct(30,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.5* rct(31,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.1* rct(32,k) * xnew(CH3O2) * xnew(HO2)  &
    + rct(33,k) * xnew(CH3OH) * xnew(OH)  &
    + 0.4* rct(34,k) * xnew(CH3OOH) * xnew(OH)  &
    + rct(46,k) * xnew(HCOCH2O2) * xnew(NO)  &
    + 2.3e-12 * xnew(HCOCH2O2) * xnew(NO3)  &
    + 0.6* 2.0e-12 * xnew(RO2POOL) * xnew(HCOCH2O2)  &
    + 0.146* rct(46,k) * xnew(MEKO2) * xnew(NO)  &
    + 0.146* 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + 0.1095* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 2.3e-13 * xnew(BIACET) * xnew(OH)  &
    + 2* rct(55,k) * xnew(HOCH2CH2O)  &
    + 0.75* rct(45,k) * xnew(HOCH2CO3) * xnew(HO2)  &
    + rct(42,k) * xnew(HOCH2CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(HOCH2CO3) * xnew(NO3)  &
    + 0.7* 1.00e-11 * xnew(RO2POOL) * xnew(HOCH2CO3)  &
    + 2.73e-12 * xnew(HOCH2CO2H) * xnew(OH)
  
  !-> HCHO cont.
    P = P +   1.12e-12 * xnew(PHAN) * xnew(OH)  &
    + 1.1575* rct(57,k) * xnew(C2H4) * xnew(O3)  &
    + 0.68* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + rct(46,k) * xnew(PRRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(PRRO2) * xnew(NO3)  &
    + 0.6* 8.8e-13 * xnew(RO2POOL) * xnew(PRRO2)  &
    + 0.13* rct(46,k) * xnew(NRN9O2) * xnew(NO)  &
    + 0.13* 2.3e-12 * xnew(NRN9O2) * xnew(NO3)  &
    + 8.3337* rcphot(IDMEK,k) * xnew(NOA)  &
    + rct(46,k) * xnew(IBUTOLBO2) * xnew(NO)  &
    + 2.3e-12 * xnew(IBUTOLBO2) * xnew(NO3)  &
    + 0.7775* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.9* rct(46,k) * xnew(ISRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NO3) * xnew(ISRO2)  &
    + 1.45e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.4* 3.25e-11 * xnew(ISON) * xnew(OH)  &
    + 1.55e-12 * xnew(NALD) * xnew(OH)  &
    + 0.2421* rct(74,k) * xnew(MACR) * xnew(O3)  &
    + 0.15* rct(46,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.15* 2.3e-12 * xnew(MACRO2) * xnew(NO3)
  
  !-> HCHO cont.
    P = P +   0.44* rct(45,k) * xnew(MACO3) * xnew(HO2)  &
    + 0.7* 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)  &
    + rct(76,k) * xnew(MACO3) * xnew(NO)  &
    + 4.0e-12 * xnew(MACO3) * xnew(NO3)  &
    + 8.2e-18 * xnew(MPAN) * xnew(O3)  &
    + 0.645* rct(80,k) * xnew(MVK) * xnew(O3)  &
    + 0.28944* rct(46,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.3* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.18* 1.22e-12 * xnew(RO2POOL) * xnew(HMVKO2)  &
    + rct(84,k) * xnew(ISRO2)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + 1.51e-11 * xnew(MACO2H) * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HOCH2CHO)  &
    + rcphot(IDCH3CHO,k) * xnew(NALD)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + 0.1* rcphot(IDCHOCHO,k) * xnew(GLYOX)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.04767 * xnew(MACR)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + 0.146* rcphot(IDCH3O2H,k) * xnew(MEKO2H)
  
  !-> HCHO cont.
    P = P +   0.146* rcphot(IDMEK,k) * xnew(MEKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(HCOCH2OOH)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HCOCH2OOH)  &
    + 0.3* rcphot(IDCH3O2H,k) * xnew(MVKO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + rcphot(IDHCHO_H,k) *0.1814 * xnew(ISON)  &
    + rcemis(HCHO,k)
  
  L = rct(35,k) * xnew(OH)  &
    + rct(36,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k)  &
    + rcphot(IDHCHO_H2,k)
  
  xnew(HCHO) = (xold(HCHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CHO

  P = rct(38,k) * xnew(C2H5O2) * xnew(NO)  &
    + 2.3e-12 * xnew(C2H5O2) * xnew(NO3)  &
    + 0.8* 3.10e-13 * xnew(RO2POOL) * xnew(C2H5O2)  &
    + 0.69* 1.17e-11 * xnew(C2H5OOH) * xnew(OH)  &
    + 0.95* 3.22e-12 * xnew(C2H5OH) * xnew(OH)  &
    + 0.264* rct(46,k) * xnew(C4H9O2) * xnew(NO)  &
    + 0.284* 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + 0.854* rct(46,k) * xnew(MEKO2) * xnew(NO)  &
    + 0.854* 2.3e-12 * xnew(MEKO2) * xnew(NO3)  &
    + 0.6405* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 1.27e-12 * xnew(PPN) * xnew(OH)  &
    + 0.489* rct(59,k) * xnew(C3H6) * xnew(O3)  &
    + rct(46,k) * xnew(PRRO2) * xnew(NO)  &
    + 2.3e-12 * xnew(PRRO2) * xnew(NO3)  &
    + 0.6* 8.8e-13 * xnew(RO2POOL) * xnew(PRRO2)  &
    + 0.13* rct(46,k) * xnew(NRN9O2) * xnew(NO)  &
    + 0.13* 2.3e-12 * xnew(NRN9O2) * xnew(NO3)  &
    + 0.1* rct(80,k) * xnew(MVK) * xnew(O3)  &
    + 0.6* rcphot(IDRCOCHO,k) * xnew(CH3COCO2H)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)
  
  !-> CH3CHO cont.
    P = P +   0.284* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + 0.854* rcphot(IDCH3O2H,k) * xnew(MEKO2H)  &
    + 0.854* rcphot(IDMEK,k) * xnew(MEKO2H)  &
    + 0.348* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rcemis(CH3CHO,k)
  
  L = rct(40,k) * xnew(OH)  &
    + rct(41,k) * xnew(NO3)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(CH3CHO) = (xold(CH3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> MACR

  P = 0.363* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.3* rct(46,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.333* 2.3e-12 * xnew(NO3) * xnew(ISRO2)  &
    + 0.3335* 1.45e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.2* 3.25e-11 * xnew(ISON) * xnew(OH)  &
    + rct(84,k) * xnew(ISRO2)  &
    + rct(46,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rct(88,k) * xnew(APINENE) * xnew(O3)  &
    + 0.16* rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + 0.333* rcphot(IDHCHO_H,k) *0.1814 * xnew(ISON)  &
    + rcemis(MACR,k)
  
  L = rct(73,k) * xnew(OH)  &
    + rct(74,k) * xnew(O3)  &
    + 3.4e-15 * xnew(NO3)  &
    + rcphot(IDHCHO_H2,k) *0.04767
  
  xnew(MACR) = (xold(MACR) + dt2 * P) / (1.0 + dt2 * L)


!-> MEK

  P = 0.847* 2.4e-11 * xnew(BURO2H) * xnew(OH)  &
    + 0.654* rct(46,k) * xnew(C4H9O2) * xnew(NO)  &
    + 0.716* 2.3e-12 * xnew(C4H9O2) * xnew(NO3)  &
    + 8.6e-13 * xnew(SC4H9NO3) * xnew(OH)  &
    + 0.716* rcphot(IDCH3O2H,k) * xnew(BURO2H)  &
    + 0.652* rcphot(IDMEK,k) *0.78 * xnew(SC4H9NO3)  &
    + rcemis(MEK,k)
  
  L = 1.1e-12 * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(MEK) = (xold(MEK) + dt2 * P) / (1.0 + dt2 * L)


!-> MVK

  P = 0.312* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.6* rct(46,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.667* 2.3e-12 * xnew(NO3) * xnew(ISRO2)  &
    + 0.6665* 1.45e-12 * xnew(RO2POOL) * xnew(ISRO2)  &
    + 0.2* 3.25e-11 * xnew(ISON) * xnew(OH)  &
    + rct(46,k) * xnew(TERPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(TERPO2) * xnew(NO3)  &
    + rct(88,k) * xnew(APINENE) * xnew(O3)  &
    + 0.84* rcphot(IDCH3O2H,k) * xnew(ISRO2H)  &
    + 0.667* rcphot(IDHCHO_H,k) *0.1814 * xnew(ISON)
  
  L = rct(79,k) * xnew(OH)  &
    + rct(80,k) * xnew(O3)  &
    + rcphot(IDHCHO_H2,k) *0.078
  
  xnew(MVK) = (xold(MVK) + dt2 * P) / (1.0 + dt2 * L)


!-> HOCH2CHO

  P = 0.2* 2.0e-12 * xnew(RO2POOL) * xnew(HCOCH2O2)  &
    + 0.4* 1.9e-12 * xnew(RO2POOL) * xnew(ETRO2)  &
    + rct(56,k) * xnew(HOCH2CH2O)  &
    + 0.79* 1.7475e-11 * xnew(ETRO2H) * xnew(OH)  &
    + 0.67536* rct(46,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.7* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.42* 1.22e-12 * xnew(RO2POOL) * xnew(HMVKO2)  &
    + 0.32* rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + 0.7* rcphot(IDCH3O2H,k) * xnew(MVKO2H)  &
    + rcphot(IDMEK,k) *13. * xnew(MVKNO3)
  
  L = 1.00e-11 * xnew(OH)  &
    + rct(41,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k) *0.6
  
  xnew(HOCH2CHO) = (xold(HOCH2CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> HCOCH2OOH

  P = rct(47,k) * xnew(HCOCH2O2) * xnew(HO2)
  
  L = 3.3e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDHCHO_H,k) *0.6
  
  xnew(HCOCH2OOH) = (xold(HCOCH2OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> ACETOL

  P = rct(60,k) * xnew(PRRO2) * xnew(HO2)  &
    + 0.4* 8.8e-13 * xnew(RO2POOL) * xnew(PRRO2)  &
    + 0.6* 3.25e-11 * xnew(ISON) * xnew(OH)  &
    + 0.85* rct(46,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.85* 2.3e-12 * xnew(MACRO2) * xnew(NO3)  &
    + rct(75,k) * xnew(MACRO2)  &
    + 2.9e-11 * xnew(MPAN) * xnew(OH)  &
    + 3.77e-11 * xnew(MACROOH) * xnew(OH)  &
    + 0.68* rcphot(IDHCHO_H2,k) *12. * xnew(HPALD)  &
    + rcphot(IDHCHO_H,k) *1.87 * xnew(MACROOH)  &
    + rcphot(IDMEK,k) *133.34 * xnew(ISON2)
  
  L = rct(78,k) * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(ACETOL) = (xold(ACETOL) + dt2 * P) / (1.0 + dt2 * L)


!-> GLYOX

  P = 0.2* 2.0e-12 * xnew(RO2POOL) * xnew(HCOCH2O2)  &
    + 0.89* 3.3e-11 * xnew(HCOCH2OOH) * xnew(OH)  &
    + 0.2* 1.00e-11 * xnew(HOCH2CHO) * xnew(OH)  &
    + rct(46,k) * xnew(NOXYOLO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NOXYOLO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(NOXYOLO2)  &
    + rct(46,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(66,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.555* 2.4e-17 * xnew(ISON2) * xnew(O3)  &
    + rcphot(IDCH3O2H,k) * xnew(NOXYOLOOH)  &
    + rcphot(IDMEK,k) *0.78 * xnew(NOXYOLOOH)  &
    + rcphot(IDMEK,k) *83. * xnew(CO2N3CHO)  &
    + rcemis(GLYOX,k)
  
  L = rct(68,k) * xnew(OH)  &
    + rcphot(IDCHOCHO,k)
  
  xnew(GLYOX) = (xold(GLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MGLYOX

  P = 6.70e-13 * xnew(NOA) * xnew(OH)  &
    + 0.28015* rct(46,k) * xnew(OXYBIPERO2) * xnew(NO)  &
    + 0.325* 2.3e-12 * xnew(OXYBIPERO2) * xnew(NO3)  &
    + rct(46,k) * xnew(C5DICARBO2) * xnew(NO)  &
    + rct(66,k) * xnew(C5DICARBO2) * xnew(NO3)  &
    + 0.9016* rct(74,k) * xnew(MACR) * xnew(O3)  &
    + 0.15* rct(46,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.15* 2.3e-12 * xnew(MACRO2) * xnew(NO3)  &
    + rct(78,k) * xnew(ACETOL) * xnew(OH)  &
    + 0.545* rct(80,k) * xnew(MVK) * xnew(O3)  &
    + 0.28944* rct(46,k) * xnew(HMVKO2) * xnew(NO)  &
    + 0.3* 2.3e-12 * xnew(HMVKO2) * xnew(NO3)  &
    + 0.18* 1.22e-12 * xnew(RO2POOL) * xnew(HMVKO2)  &
    + 0.445* 2.4e-17 * xnew(ISON2) * xnew(O3)  &
    + 7.61e-11 * xnew(OH) * xnew(HPALD)  &
    + 0.44* rct(45,k) * xnew(CO2H3CO3) * xnew(HO2)  &
    + rct(42,k) * xnew(CO2H3CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(CO2H3CO3) * xnew(NO3)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(CO2H3CO3)  &
    + 3.74e-12 * xnew(C4PAN6) * xnew(OH)
  
  !-> MGLYOX cont.
    P = P +   rcphot(IDHCHO_H,k) *0.6 * xnew(CO2H3CHO)  &
    + 0.3* rcphot(IDCH3O2H,k) * xnew(MVKO2H)  &
    + 0.325* rcphot(IDCH3O2H,k) * xnew(OXYBPEROOH)  &
    + rcphot(IDCH3O2H,k) * xnew(CO2H3CO3H)  &
    + rcemis(MGLYOX,k)
  
  L = rct(69,k) * xnew(OH)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(MGLYOX) = (xold(MGLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> EGLYOX

  P = 0.0365* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.1206* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + 0.28015* rct(46,k) * xnew(OXYBIPERO2) * xnew(NO)  &
    + 0.325* 2.3e-12 * xnew(OXYBIPERO2) * xnew(NO3)  &
    + 0.325* rcphot(IDCH3O2H,k) * xnew(OXYBPEROOH)
  
  L = rct(52,k) * xnew(NO3)  &
    + 1.39e-11 * xnew(OH)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(EGLYOX) = (xold(EGLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> CO2H3CHO

  P = 0.3* 4.5e-11 * xnew(MVKO2H) * xnew(OH)  &
    + 0.12* 1.22e-12 * xnew(RO2POOL) * xnew(HMVKO2)
  
  L = 2.45e-11 * xnew(OH)  &
    + rct(81,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k) *0.6
  
  xnew(CO2H3CHO) = (xold(CO2H3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CO23C3CHO

  P = 2.96e-12 * xnew(OH) * xnew(BIACETOH)
  
  L = 1.23e-11 * xnew(OH)  &
    + rct(81,k) * xnew(NO3)  &
    + rcphot(IDCH3COY,k)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(CO23C3CHO) = (xold(CO23C3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> BIACET

  P = 0.2135* 1.04e-12 * xnew(RO2POOL) * xnew(MEKO2)  &
    + 0.7054* 2.1e-11 * xnew(MEKO2H) * xnew(OH)  &
    + 0.3017* rct(46,k) * xnew(OXYBIPERO2) * xnew(NO)  &
    + 0.35* 2.3e-12 * xnew(OXYBIPERO2) * xnew(NO3)  &
    + 0.35* rcphot(IDCH3O2H,k) * xnew(OXYBPEROOH)
  
  L = 2.3e-13 * xnew(OH)  &
    + rcphot(IDCH3COY,k)
  
  xnew(BIACET) = (xold(BIACET) + dt2 * P) / (1.0 + dt2 * L)


!-> BIACETOH

  P = 0.11* 1.16e-11 * xnew(IEPOX) * xnew(OH)  &
    + 0.7* 4.5e-11 * xnew(MVKO2H) * xnew(OH)  &
    + 0.28* 1.22e-12 * xnew(RO2POOL) * xnew(HMVKO2)  &
    + 0.33* 1.33e-12 * xnew(MVKNO3) * xnew(OH)
  
  L = 2.96e-12 * xnew(OH)  &
    + rcphot(IDCH3COY,k)
  
  xnew(BIACETOH) = (xold(BIACETOH) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICARB

  P = 0.28015* rct(46,k) * xnew(OXYBIPERO2) * xnew(NO)  &
    + 0.325* 2.3e-12 * xnew(OXYBIPERO2) * xnew(NO3)  &
    + 0.325* rcphot(IDCH3O2H,k) * xnew(OXYBPEROOH)
  
  L = 6.2e-11 * xnew(OH)  &
    + rcphot(IDNO2,k) *0.2
  
  xnew(C5DICARB) = (xold(C5DICARB) + dt2 * P) / (1.0 + dt2 * L)


!-> BZFUONE

  P = 0.58185* rct(46,k) * xnew(OXYBIPERO2) * xnew(NO)  &
    + 0.675* 2.3e-12 * xnew(OXYBIPERO2) * xnew(NO3)  &
    + 0.675* rcphot(IDCH3O2H,k) * xnew(OXYBPEROOH)
  
  L = 3.0e-13 * xnew(NO3)  &
    + 4.45e-11 * xnew(OH)
  
  xnew(BZFUONE) = (xold(BZFUONE) + dt2 * P) / (1.0 + dt2 * L)


!-> HOCH2CO2H

  P = 0.25* rct(45,k) * xnew(HOCH2CO3) * xnew(HO2)  &
    + 0.3* 1.00e-11 * xnew(RO2POOL) * xnew(HOCH2CO3)
  
  L = 2.73e-12 * xnew(OH)
  
  xnew(HOCH2CO2H) = (xold(HOCH2CO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3COCO2H

  P = rct(46,k) * xnew(OXNCATECO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXNCATECO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(OXNCATECO2)  &
    + rct(46,k) * xnew(OXNNCATCO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXNNCATCO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(OXNNCATCO2)  &
    + rct(46,k) * xnew(NOXYOLO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NOXYOLO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(NOXYOLO2)  &
    + 0.075* rct(80,k) * xnew(MVK) * xnew(O3)  &
    + rcphot(IDCH3O2H,k) * xnew(OXNCATCOOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OXNNCATOOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NOXYOLOOH)  &
    + rcphot(IDMEK,k) *0.78 * xnew(NOXYOLOOH)
  
  L = rcphot(IDRCOCHO,k)
  
  xnew(CH3COCO2H) = (xold(CH3COCO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MACO2H

  P = 0.105* rct(71,k) * xnew(C5H8) * xnew(O3)  &
    + 0.56* rct(45,k) * xnew(MACO3) * xnew(HO2)  &
    + 0.3* 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)
  
  L = 1.51e-11 * xnew(OH)
  
  xnew(MACO2H) = (xold(MACO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OOH

  P = 0.9* rct(32,k) * xnew(CH3O2) * xnew(HO2)
  
  L = rct(34,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + 1.0e-5
  
  xnew(CH3OOH) = (xold(CH3OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OOH

  P = rct(39,k) * xnew(C2H5O2) * xnew(HO2)  &
    + rct(49,k) * xnew(IBUTOLBO2) * xnew(HO2)
  
  L = 1.17e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(C2H5OOH) = (xold(C2H5OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> BURO2H

  P = rct(49,k) * xnew(C4H9O2) * xnew(HO2)
  
  L = 2.4e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(BURO2H) = (xold(BURO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> DM123OOH

  P = rct(64,k) * xnew(DM123O2) * xnew(HO2)
  
  L = 4.65e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(DM123OOH) = (xold(DM123OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> ETRO2H

  P = rct(54,k) * xnew(ETRO2) * xnew(HO2)
  
  L = 1.7475e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ETRO2H) = (xold(ETRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MEKO2H

  P = rct(49,k) * xnew(MEKO2) * xnew(HO2)
  
  L = 2.1e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDMEK,k)
  
  xnew(MEKO2H) = (xold(MEKO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> ISRO2H

  P = rct(67,k) * xnew(ISRO2) * xnew(HO2)
  
  L = 6.0e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ISRO2H) = (xold(ISRO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> C5DICAROOH

  P = rct(67,k) * xnew(C5DICARBO2) * xnew(HO2)
  
  L = rct(65,k) * xnew(OH)
  
  xnew(C5DICAROOH) = (xold(C5DICAROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> OCATEC1OOH

  P = rct(64,k) * xnew(OXCATEC1O2) * xnew(HO2)
  
  L = 2.05e-10 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(OCATEC1OOH) = (xold(OCATEC1OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> HPALD

  P = rct(83,k) * xnew(ISRO2)
  
  L = 7.61e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H2,k) *12.
  
  xnew(HPALD) = (xold(HPALD) + dt2 * P) / (1.0 + dt2 * L)


!-> MVKO2H

  P = rct(49,k) * xnew(HMVKO2) * xnew(HO2)
  
  L = 4.5e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(MVKO2H) = (xold(MVKO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> MACROOH

  P = rct(49,k) * xnew(MACRO2) * xnew(HO2)
  
  L = 3.77e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *1.87
  
  xnew(MACROOH) = (xold(MACROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYBPEROOH

  P = rct(64,k) * xnew(OXYBIPERO2) * xnew(HO2)
  
  L = 7.62e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(OXYBPEROOH) = (xold(OXYBPEROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3H

  P = 0.384* rct(45,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = 1.1e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(CH3CO3H) = (xold(CH3CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> CO2H3CO3H

  P = 0.56* rct(45,k) * xnew(CO2H3CO3) * xnew(HO2)
  
  L = 7.34e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(CO2H3CO3H) = (xold(CO2H3CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> IEPOX

  P = 0.93* 6.0e-11 * xnew(ISRO2H) * xnew(OH)
  
  L = 1.16e-11 * xnew(OH)
  
  xnew(IEPOX) = (xold(IEPOX) + dt2 * P) / (1.0 + dt2 * L)


!-> SC4H9NO3

  P = 0.083* rct(46,k) * xnew(C4H9O2) * xnew(NO)
  
  L = 8.6e-13 * xnew(OH)  &
    + rcphot(IDMEK,k) *0.78
  
  xnew(SC4H9NO3) = (xold(SC4H9NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> MVKNO3

  P = 0.0352* rct(46,k) * xnew(HMVKO2) * xnew(NO)
  
  L = 1.33e-12 * xnew(OH)  &
    + rcphot(IDMEK,k) *13.
  
  xnew(MVKNO3) = (xold(MVKNO3) + dt2 * P) / (1.0 + dt2 * L)


!-> OXY1OHNO2

  P = 2.08e-12 * xnew(OXY1O) * xnew(NO2)
  
  L = 3.13e-13 * xnew(NO3)  &
    + 1.53e-12 * xnew(OH)
  
  xnew(OXY1OHNO2) = (xold(OXY1OHNO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OXNCATECH

  P = 2.08e-12 * xnew(OXCATEC1O) * xnew(NO2)
  
  L = 5.03e-12 * xnew(NO3)  &
    + 6.83e-12 * xnew(OH)
  
  xnew(OXNCATECH) = (xold(OXNCATECH) + dt2 * P) / (1.0 + dt2 * L)


!-> NOXYOL1OOH

  P = rct(64,k) * xnew(NOXYOL1O2) * xnew(HO2)
  
  L = 1.53e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(NOXYOL1OOH) = (xold(NOXYOL1OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> CO2N3CHO

  P = 0.67* 1.33e-12 * xnew(MVKNO3) * xnew(OH)
  
  L = rcphot(IDMEK,k) *83.
  
  xnew(CO2N3CHO) = (xold(CO2N3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> NALD

  P = 0.6* 3.25e-11 * xnew(ISON) * xnew(OH)
  
  L = 1.55e-12 * xnew(OH)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(NALD) = (xold(NALD) + dt2 * P) / (1.0 + dt2 * L)


!-> NOA

  P = rct(60,k) * xnew(NRN9O2) * xnew(HO2)  &
    + 0.87* rct(46,k) * xnew(NRN9O2) * xnew(NO)  &
    + 0.87* 2.3e-12 * xnew(NRN9O2) * xnew(NO3)  &
    + 0.555* 2.4e-17 * xnew(ISON2) * xnew(O3)  &
    + rct(45,k) * xnew(NC4CO3) * xnew(HO2)  &
    + rct(42,k) * xnew(NC4CO3) * xnew(NO)  &
    + 4.0e-12 * xnew(NC4CO3) * xnew(NO3)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(NC4CO3)  &
    + 2.16e-11 * xnew(C5PAN18) * xnew(OH)
  
  L = 8.3337* rcphot(IDMEK,k)  &
    + 6.70e-13 * xnew(OH)
  
  xnew(NOA) = (xold(NOA) + dt2 * P) / (1.0 + dt2 * L)


!-> ISON

  P = 0.1* rct(46,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.2* 1.3e-12 * xnew(RO2POOL) * xnew(NISOPO2)
  
  L = 3.25e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *0.1814
  
  xnew(ISON) = (xold(ISON) + dt2 * P) / (1.0 + dt2 * L)


!-> ISON2

  P = rct(67,k) * xnew(NISOPO2) * xnew(HO2)  &
    + rct(46,k) * xnew(NISOPO2) * xnew(NO)  &
    + 2.3e-12 * xnew(NISOPO2) * xnew(NO3)  &
    + 0.8* 1.3e-12 * xnew(RO2POOL) * xnew(NISOPO2)
  
  L = 4.16e-11 * xnew(OH)  &
    + 2.4e-17 * xnew(O3)  &
    + rct(82,k) * xnew(NO3)  &
    + rcphot(IDMEK,k) *133.34
  
  xnew(ISON2) = (xold(ISON2) + dt2 * P) / (1.0 + dt2 * L)


!-> NC4MDCO2H

  P = rct(46,k) * xnew(OXNCATECO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXNCATECO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(OXNCATECO2)  &
    + rct(46,k) * xnew(OXNNCATCO2) * xnew(NO)  &
    + 2.3e-12 * xnew(OXNNCATCO2) * xnew(NO3)  &
    + 8.0e-13 * xnew(RO2POOL) * xnew(OXNNCATCO2)  &
    + rcphot(IDCH3O2H,k) * xnew(OXNCATCOOH)  &
    + rcphot(IDCH3O2H,k) * xnew(OXNNCATOOH)
  
  L = rct(65,k) * xnew(OH)
  
  xnew(NC4MDCO2H) = (xold(NC4MDCO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> NOXYOLOOH

  P = rct(64,k) * xnew(NOXYOLO2) * xnew(HO2)  &
    + 0.138* rct(46,k) * xnew(OXYBIPERO2) * xnew(NO)
  
  L = 1.43e-10 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDMEK,k) *0.78
  
  xnew(NOXYOLOOH) = (xold(NOXYOLOOH) + dt2 * P) / (1.0 + dt2 * L)


!-> OXNCATCOOH

  P = rct(64,k) * xnew(OXNCATECO2) * xnew(HO2)
  
  L = rct(65,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(OXNCATCOOH) = (xold(OXNCATCOOH) + dt2 * P) / (1.0 + dt2 * L)


!-> PAN

  P = rct(43,k) * xnew(CH3CO3) * xnew(NO2)
  
  L = rct(44,k)
  
  xnew(PAN) = (xold(PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> PHAN

  P = rct(43,k) * xnew(HOCH2CO3) * xnew(NO2)
  
  L = 1.12e-12 * xnew(OH)  &
    + rct(44,k)
  
  xnew(PHAN) = (xold(PHAN) + dt2 * P) / (1.0 + dt2 * L)


!-> PPN

  P = rct(43,k) * xnew(C2H5CO3) * xnew(NO2)
  
  L = 1.27e-12 * xnew(OH)  &
    + rct(51,k)
  
  xnew(PPN) = (xold(PPN) + dt2 * P) / (1.0 + dt2 * L)


!-> MPAN

  P = rct(43,k) * xnew(MACO3) * xnew(NO2)
  
  L = rct(77,k)  &
    + 2.9e-11 * xnew(OH)  &
    + 8.2e-18 * xnew(O3)
  
  xnew(MPAN) = (xold(MPAN) + dt2 * P) / (1.0 + dt2 * L)


!-> C4PAN6

  P = rct(43,k) * xnew(CO2H3CO3) * xnew(NO2)
  
  L = 3.74e-12 * xnew(OH)  &
    + rct(44,k)
  
  xnew(C4PAN6) = (xold(C4PAN6) + dt2 * P) / (1.0 + dt2 * L)


!-> C5PAN18

  P = rct(43,k) * xnew(NC4CO3) * xnew(NO2)
  
  L = rct(44,k)  &
    + 2.16e-11 * xnew(OH)
  
  xnew(C5PAN18) = (xold(C5PAN18) + dt2 * P) / (1.0 + dt2 * L)


!-> OXNNCATOOH

  P = rct(64,k) * xnew(OXNNCATCO2) * xnew(HO2)
  
  L = rct(65,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(OXNNCATOOH) = (xold(OXNNCATOOH) + dt2 * P) / (1.0 + dt2 * L)


!-> APINENE

  P = RCBIO(2,K)
  
  L = rct(85,k) * xnew(OH)  &
    + rct(88,k) * xnew(O3)  &
    + rct(89,k) * xnew(NO3)
  
  xnew(APINENE) = (xold(APINENE) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPOOH

  P = rct(86,k) * xnew(TERPO2) * xnew(HO2)
  
  L = rct(87,k) * xnew(OH)
  
  xnew(TERPOOH) = (xold(TERPOOH) + dt2 * P) / (1.0 + dt2 * L)


!-> SO2

  P = rcemis(SO2,k)
  
  L = AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH)  &
    + AQRCK(ICLRC1,K) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(O3)  &
    + AQRCK(ICLRC3,K)
  
  xnew(SO2) = (xold(SO2) + dt2 * P) / (1.0 + dt2 * L)


!-> shipNOx

  P = rcemis(SHIPNOX,k)
  
  L = rct(20,k) * xnew(OH)  &
    + 3.2e-5
  
  xnew(shipNOx) = (xold(shipNOx) + dt2 * P) / (1.0 + dt2 * L)


!-> Dust_road_f

  P = rcemis(DUST_ROAD_F,K)
  
  ! L = 0.0
  
  xnew(Dust_road_f) = xold(Dust_road_f) + dt2 * P


!-> Dust_road_c

  P = rcemis(DUST_ROAD_C,K)
  
  ! L = 0.0
  
  xnew(Dust_road_c) = xold(Dust_road_c) + dt2 * P


!-> Dust_wb_f

  P = rcemis(DUST_WB_F,K)
  
  ! L = 0.0
  
  xnew(Dust_wb_f) = xold(Dust_wb_f) + dt2 * P


!-> Dust_wb_c

  P = rcemis(DUST_WB_C,K)
  
  ! L = 0.0
  
  xnew(Dust_wb_c) = xold(Dust_wb_c) + dt2 * P


!-> Dust_sah_f

  P = rcemis(DUST_SAH_F,K)
  
  ! L = 0.0
  
  xnew(Dust_sah_f) = xold(Dust_sah_f) + dt2 * P


!-> Dust_sah_c

  P = rcemis(DUST_SAH_C,K)
  
  ! L = 0.0
  
  xnew(Dust_sah_c) = xold(Dust_sah_c) + dt2 * P


!-> SQT_SOA_NV

  P = rcemis(SQT_SOA_NV,k)  &
    + 0.00383* RCBIO(2,K)
  
  ! L = 0.0
  
  xnew(SQT_SOA_NV) = xold(SQT_SOA_NV) + dt2 * P


