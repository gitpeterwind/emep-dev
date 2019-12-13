!-> BZO2

  P = 0.352* rct(121,k) * xnew(BENZ) * xnew(OH)
  
  L = rct(56,k) * xnew(NO)  &
    + rct(46,k) * xnew(C2O3)  &
    + rct(122,k) * xnew(HO2)  &
    + rct(52,k) * xnew(RO2)
  
  xnew(BZO2) = (xold(BZO2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2O3

  P = rct(44,k) * xnew(PAN)  &
    + 0.6* rcphot(IDHO2NO2,k) *0.14 * xnew(PAN)  &
    + rct(46,k) * xnew(C2O3) * xnew(RO2)  &
    + rct(58,k) * xnew(PACD) * xnew(OH)  &
    + rct(69,k) * xnew(ALD2) * xnew(OP)  &
    + rct(70,k) * xnew(ALD2) * xnew(OH)  &
    + rct(71,k) * xnew(ALD2) * xnew(NO3)  &
    + 0.8* rct(75,k) * xnew(GLYD) * xnew(OH)  &
    + rct(71,k) * xnew(GLYD) * xnew(NO3)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLY)  &
    + rct(71,k) * xnew(MGLY) * xnew(NO3)  &
    + rct(77,k) * xnew(MGLY) * xnew(OH)  &
    + 0.5* rcphot(IDACETON,k) *1.5 * xnew(KET)  &
    + 0.62* rcphot(IDACETON,k) * xnew(ACET)  &
    + rct(84,k) * xnew(ACET) * xnew(OH)  &
    + 0.269* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.143* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.208* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.3* rcphot(IDNO2,k) *0.087 * xnew(XOPN)  &
    + 0.6* rct(132,k) * xnew(XOPN) * xnew(O3)
  
  !-> C2O3 cont.
    P = P +   0.12* rct(135,k) * xnew(OPEN) * xnew(O3)
  
  L = rct(42,k) * xnew(NO)  &
    + rct(43,k) * xnew(NO2)  &
    + rct(45,k) * xnew(HO2)  &
    + rct(46,k) * xnew(RO2)  &
    + rct(47,k) * xnew(C2O3)  &
    + rct(47,k) * xnew(C2O3)  &
    + rct(47,k) * xnew(CXO3)  &
    + rct(55,k) * xnew(MEO2)  &
    + rct(46,k) * xnew(XO2H)  &
    + rct(46,k) * xnew(XO2)  &
    + rct(46,k) * xnew(XO2N)  &
    + rct(46,k) * xnew(ISO2)  &
    + rct(46,k) * xnew(EPX2)  &
    + rct(46,k) * xnew(BZO2)  &
    + rct(46,k) * xnew(TO2)  &
    + rct(46,k) * xnew(XLO2)  &
    + rct(47,k) * xnew(OPO3)
  
  xnew(C2O3) = (xold(C2O3) + dt2 * P) / (1.0 + dt2 * L)


!-> CRO

  P = 0.2* rct(125,k) * xnew(CRES) * xnew(OH)  &
    + 0.3* rct(126,k) * xnew(CRES) * xnew(NO3)  &
    + 0.5* rct(129,k) * xnew(CRON) * xnew(OH)  &
    + 0.5* rct(130,k) * xnew(CRON) * xnew(NO3)  &
    + 0.5* rct(136,k) * xnew(CAT1) * xnew(OH)  &
    + rct(137,k) * xnew(CAT1) * xnew(NO3)
  
  L = rct(127,k) * xnew(NO2)  &
    + rct(128,k) * xnew(HO2)
  
  xnew(CRO) = (xold(CRO) + dt2 * P) / (1.0 + dt2 * L)


!-> CXO3

  P = rct(44,k) * xnew(PANX)  &
    + 0.6* rcphot(IDHO2NO2,k) *0.14 * xnew(PANX)  &
    + rct(72,k) * xnew(ALDX) * xnew(OP)  &
    + rct(73,k) * xnew(ALDX) * xnew(OH)  &
    + rct(74,k) * xnew(ALDX) * xnew(NO3)  &
    + 0.5* rcphot(IDACETON,k) *1.5 * xnew(KET)  &
    + 0.25* rct(104,k) * xnew(C5H8) * xnew(OP)  &
    + 0.2* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.717* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + 0.39* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + 0.2* rct(138,k) * xnew(OPO3) * xnew(NO)
  
  L = rct(47,k) * xnew(C2O3)  &
    + rct(48,k) * xnew(NO)  &
    + rct(43,k) * xnew(NO2)  &
    + rct(45,k) * xnew(HO2)  &
    + rct(46,k) * xnew(RO2)  &
    + rct(49,k) * xnew(CXO3)  &
    + rct(49,k) * xnew(CXO3)
  
  xnew(CXO3) = (xold(CXO3) + dt2 * P) / (1.0 + dt2 * L)


!-> EPX2

  P = rct(115,k) * xnew(EPOX) * xnew(OH)
  
  L = rct(106,k) * xnew(HO2)  &
    + rct(105,k) * xnew(NO)  &
    + rct(46,k) * xnew(C2O3)  &
    + rct(52,k) * xnew(RO2)
  
  xnew(EPX2) = (xold(EPX2) + dt2 * P) / (1.0 + dt2 * L)


!-> HCO3

  P = rct(65,k) * xnew(FORM) * xnew(HO2)
  
  L = rct(66,k)  &
    + rct(67,k) * xnew(NO)  &
    + rct(68,k) * xnew(HO2)
  
  xnew(HCO3) = (xold(HCO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2

  P = rct(9,k) * xnew(O3) * xnew(OH)  &
    + rct(11,k) * xnew(OH) * xnew(OP)  &
    + rct(18,k) * xnew(H2O2) * xnew(OH)  &
    + rct(19,k) * xnew(H2O2) * xnew(OP)  &
    + rct(26,k) * xnew(NO3) * xnew(OH)  &
    + rct(40,k) * xnew(PNA)  &
    + 0.59* rcphot(IDHO2NO2,k) * xnew(PNA)  &
    + rct(51,k) * xnew(RO2) * xnew(HO2)  &
    + rct(53,k) * xnew(MEO2) * xnew(NO)  &
    + 0.9* rct(55,k) * xnew(MEO2) * xnew(C2O3)  &
    + 0.37* rct(52,k) * xnew(MEO2) * xnew(RO2)  &
    + rct(56,k) * xnew(XO2H) * xnew(NO)  &
    + 0.8* rct(46,k) * xnew(XO2H) * xnew(C2O3)  &
    + 0.6* rct(52,k) * xnew(XO2H) * xnew(RO2)  &
    + 0.8* rct(46,k) * xnew(XO2N) * xnew(C2O3)  &
    + rcphot(IDCH3O2H,k) * xnew(ROOH)  &
    + rct(60,k) * xnew(FACD) * xnew(OH)  &
    + rct(62,k) * xnew(FORM) * xnew(OH)  &
    + 2* rcphot(IDHCHO_H,k) * xnew(FORM)  &
    + rct(63,k) * xnew(FORM) * xnew(OP)
  
  !-> HO2 cont.
    P = P +   rct(64,k) * xnew(FORM) * xnew(NO3)  &
    + rct(66,k) * xnew(HCO3)  &
    + rct(67,k) * xnew(HCO3) * xnew(NO)  &
    + 0.2* rct(68,k) * xnew(HCO3) * xnew(HO2)  &
    + rcphot(IDCH3CHO,k) * xnew(ALD2)  &
    + rcphot(IDMEK,k) *5.4 * xnew(ALDX)  &
    + 0.2* rct(75,k) * xnew(GLYD) * xnew(OH)  &
    + 1.4* rcphot(IDACETON,k) *10. * xnew(GLYD)  &
    + rct(76,k) * xnew(GLY) * xnew(OH)  &
    + 2* rcphot(IDCHOCHO,k) * xnew(GLY)  &
    + rct(71,k) * xnew(GLY) * xnew(NO3)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLY)  &
    + rct(78,k) * xnew(H2) * xnew(OH)  &
    + rct(79,k) * xnew(CO) * xnew(OH)  &
    + rct(82,k) * xnew(MEOH) * xnew(OH)  &
    + 0.9* rct(83,k) * xnew(ETOH) * xnew(OH)  &
    + rct(88,k) * xnew(ROR)  &
    + 0.3* rct(90,k) * xnew(ETHY) * xnew(OH)  &
    + rct(91,k) * xnew(ETH) * xnew(OP)  &
    + 0.16* rct(93,k) * xnew(ETH) * xnew(O3)
  
  !-> HO2 cont.
    P = P +   0.1* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.08* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.25* rct(104,k) * xnew(C5H8) * xnew(OP)  &
    + 0.818* rct(105,k) * xnew(ISO2) * xnew(NO)  &
    + 0.12* rct(106,k) * xnew(ISO2) * xnew(HO2)  &
    + 0.728* rct(46,k) * xnew(ISO2) * xnew(C2O3)  &
    + 0.728* rct(52,k) * xnew(ISO2) * xnew(RO2)  &
    + rct(107,k) * xnew(ISO2)  &
    + 0.066* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.137* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.398* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.76* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.825* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 0.825* rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 0.66* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 0.825* rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + 0.53* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + 0.918* rct(56,k) * xnew(BZO2) * xnew(NO)  &
    + rct(46,k) * xnew(BZO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(BZO2) * xnew(RO2)
  
  !-> HO2 cont.
    P = P +   0.18* rct(123,k) * xnew(TOL) * xnew(OH)  &
    + 0.86* rct(56,k) * xnew(TO2) * xnew(NO)  &
    + rct(46,k) * xnew(TO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(TO2) * xnew(RO2)  &
    + 0.155* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + 0.86* rct(56,k) * xnew(XLO2) * xnew(NO)  &
    + rct(46,k) * xnew(XLO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(XLO2) * xnew(RO2)  &
    + rct(125,k) * xnew(CRES) * xnew(OH)  &
    + rcphot(IDHONO,k) *0.076 * xnew(CRON)  &
    + 0.7* rcphot(IDNO2,k) *0.087 * xnew(XOPN)  &
    + rcphot(IDNO2,k) *0.087 * xnew(OPEN)  &
    + 0.56* rct(135,k) * xnew(OPEN) * xnew(O3)  &
    + 0.2* rct(136,k) * xnew(CAT1) * xnew(OH)  &
    + 0.8* rct(138,k) * xnew(OPO3) * xnew(NO)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)
  
  L = rct(10,k) * xnew(O3)  &
    + rct(12,k) * xnew(OP)  &
    + rct(15,k) * xnew(OH)  &
    + rct(16,k) * xnew(HO2)  &
    + rct(16,k) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2)  &
    + rct(21,k) * xnew(NO)  &
    + rct(27,k) * xnew(NO3)  &
    + rct(39,k) * xnew(NO2)  &
    + rct(45,k) * xnew(C2O3)  &
    + rct(45,k) * xnew(CXO3)  &
    + rct(51,k) * xnew(RO2)  &
    + rct(54,k) * xnew(MEO2)  &
    + rct(57,k) * xnew(XO2H)  &
    + rct(57,k) * xnew(XO2)  &
    + rct(57,k) * xnew(XO2N)  &
    + rct(65,k) * xnew(FORM)  &
    + rct(68,k) * xnew(HCO3)  &
    + rct(106,k) * xnew(ISO2)
  
  !-> HO2 cont.
    L = L +   rct(106,k) * xnew(EPX2)  &
    + rct(122,k) * xnew(BZO2)  &
    + rct(122,k) * xnew(TO2)  &
    + rct(122,k) * xnew(XLO2)  &
    + rct(128,k) * xnew(CRO)  &
    + rct(45,k) * xnew(OPO3)  &
    + rct(145,k)
  
  xnew(HO2) = (xold(HO2) + dt2 * P) / (1.0 + dt2 * L)


!-> ISO2

  P = rct(103,k) * xnew(C5H8) * xnew(OH)  &
    + 0.067* rct(113,k) * xnew(ISPX) * xnew(OH)
  
  L = rct(105,k) * xnew(NO)  &
    + rct(106,k) * xnew(HO2)  &
    + rct(46,k) * xnew(C2O3)  &
    + rct(52,k) * xnew(RO2)  &
    + rct(107,k)
  
  xnew(ISO2) = (xold(ISO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MEO2

  P = rct(42,k) * xnew(C2O3) * xnew(NO)  &
    + 0.4* rcphot(IDHO2NO2,k) *0.14 * xnew(PAN)  &
    + 0.44* rct(45,k) * xnew(C2O3) * xnew(HO2)  &
    + 2* rct(47,k) * xnew(C2O3) * xnew(C2O3)  &
    + rct(47,k) * xnew(C2O3) * xnew(CXO3)  &
    + 0.9* rct(55,k) * xnew(MEO2) * xnew(C2O3)  &
    + 0.8* rct(46,k) * xnew(XO2H) * xnew(C2O3)  &
    + 0.8* rct(46,k) * xnew(XO2) * xnew(C2O3)  &
    + 0.8* rct(46,k) * xnew(XO2N) * xnew(C2O3)  &
    + 0.6* rct(58,k) * xnew(CH3OOH) * xnew(OH)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rct(61,k) * xnew(AACD) * xnew(OH)  &
    + rcphot(IDCH3CHO,k) * xnew(ALD2)  &
    + rct(80,k) * xnew(CH4) * xnew(OH)  &
    + 0.5* rcphot(IDACETON,k) *1.5 * xnew(KET)  &
    + 1.38* rcphot(IDACETON,k) * xnew(ACET)  &
    + 0.8* rct(46,k) * xnew(ISO2) * xnew(C2O3)  &
    + 0.115* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.34* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.8* rct(46,k) * xnew(EPX2) * xnew(C2O3)
  
  !-> MEO2 cont.
    P = P +   rct(46,k) * xnew(BZO2) * xnew(C2O3)  &
    + rct(46,k) * xnew(TO2) * xnew(C2O3)  &
    + rct(46,k) * xnew(XLO2) * xnew(C2O3)  &
    + rct(47,k) * xnew(OPO3) * xnew(C2O3)  &
    + rct(80,k) * xnew(ECH4) * xnew(OH)
  
  L = rct(53,k) * xnew(NO)  &
    + rct(54,k) * xnew(HO2)  &
    + rct(55,k) * xnew(C2O3)  &
    + rct(52,k) * xnew(RO2)
  
  xnew(MEO2) = (xold(MEO2) + dt2 * P) / (1.0 + dt2 * L)


!-> OP

  P = rcphot(IDNO2,k) * xnew(NO2)  &
    + rcphot(IDO3_O3P,k) * xnew(O3)  &
    + rct(7,k) * xnew(OD)  &
    + rct(13,k) * xnew(OH) * xnew(OH)  &
    + 0.887* rcphot(IDNO3,k) * xnew(NO3)
  
  L = rct(1,k)  &
    + rct(3,k) * xnew(NO)  &
    + rct(4,k) * xnew(NO2)  &
    + rct(5,k) * xnew(NO2)  &
    + rct(6,k) * xnew(O3)  &
    + rct(11,k) * xnew(OH)  &
    + rct(12,k) * xnew(HO2)  &
    + rct(19,k) * xnew(H2O2)  &
    + rct(25,k) * xnew(NO3)  &
    + rct(63,k) * xnew(FORM)  &
    + rct(69,k) * xnew(ALD2)  &
    + rct(72,k) * xnew(ALDX)  &
    + rct(91,k) * xnew(ETH)  &
    + rct(95,k) * xnew(OLE)  &
    + rct(99,k) * xnew(IOLE)  &
    + rct(104,k) * xnew(C5H8)  &
    + rct(117,k) * xnew(TERP)
  
  xnew(OP) = (xold(OP) + dt2 * P) / (1.0 + dt2 * L)


!-> OD

  P = rcphot(IDO3_O1D,k) * xnew(O3)
  
  L = rct(7,k)  &
    + rct(8,k)
  
  xnew(OD) = (xold(OD) + dt2 * P) / (1.0 + dt2 * L)


!-> OH

  P = 2* rct(8,k) * xnew(OD)  &
    + rct(10,k) * xnew(O3) * xnew(HO2)  &
    + rct(12,k) * xnew(HO2) * xnew(OP)  &
    + 2* rcphot(IDH2O2,k) * xnew(H2O2)  &
    + rct(19,k) * xnew(H2O2) * xnew(OP)  &
    + rct(21,k) * xnew(HO2) * xnew(NO)  &
    + rct(27,k) * xnew(NO3) * xnew(HO2)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + 0.41* rcphot(IDHO2NO2,k) * xnew(PNA)  &
    + 0.44* rct(45,k) * xnew(C2O3) * xnew(HO2)  &
    + 0.44* rct(45,k) * xnew(CXO3) * xnew(HO2)  &
    + 0.4* rct(58,k) * xnew(CH3OOH) * xnew(OH)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + 0.4* rct(58,k) * xnew(ROOH) * xnew(OH)  &
    + rcphot(IDCH3O2H,k) * xnew(ROOH)  &
    + rct(63,k) * xnew(FORM) * xnew(OP)  &
    + 0.2* rct(68,k) * xnew(HCO3) * xnew(HO2)  &
    + rct(69,k) * xnew(ALD2) * xnew(OP)  &
    + rct(72,k) * xnew(ALDX) * xnew(OP)
  
  !-> OH cont.
    P = P +   0.19* rcphot(IDACETON,k) *10. * xnew(GLYD)  &
    + 0.7* rct(90,k) * xnew(ETHY) * xnew(OH)  &
    + 0.3* rct(91,k) * xnew(ETH) * xnew(OP)  &
    + 0.16* rct(93,k) * xnew(ETH) * xnew(O3)  &
    + 0.1* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.334* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.5* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.12* rct(106,k) * xnew(ISO2) * xnew(HO2)  &
    + 0.266* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.461* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.933* rct(113,k) * xnew(ISPX) * xnew(OH)  &
    + rcphot(IDCHOCHO,k) *7.2 * xnew(HPLD)  &
    + 1.125* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 0.125* rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 0.1* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 0.125* rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + 0.57* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + 0.118* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + 0.1* rct(123,k) * xnew(TOL) * xnew(OH)  &
    + 0.244* rct(124,k) * xnew(XYL) * xnew(OH)
  
  !-> OH cont.
    P = P +   0.5* rct(132,k) * xnew(XOPN) * xnew(O3)  &
    + 0.5* rct(135,k) * xnew(OPEN) * xnew(O3)  &
    + 0.44* rct(45,k) * xnew(OPO3) * xnew(HO2)
  
  L = rct(9,k) * xnew(O3)  &
    + rct(11,k) * xnew(OP)  &
    + rct(13,k) * xnew(OH)  &
    + rct(13,k) * xnew(OH)  &
    + rct(14,k) * xnew(OH)  &
    + rct(14,k) * xnew(OH)  &
    + rct(15,k) * xnew(HO2)  &
    + rct(18,k) * xnew(H2O2)  &
    + rct(26,k) * xnew(NO3)  &
    + rct(33,k) * xnew(NO)  &
    + rct(36,k) * xnew(HONO)  &
    + rct(37,k) * xnew(NO2)  &
    + rct(38,k) * xnew(HNO3)  &
    + rct(41,k) * xnew(PNA)  &
    + rct(58,k) * xnew(CH3OOH)  &
    + rct(58,k) * xnew(ROOH)  &
    + rct(59,k) * xnew(NTR1)  &
    + rct(60,k) * xnew(FACD)  &
    + rct(61,k) * xnew(AACD)  &
    + rct(58,k) * xnew(PACD)
  
  !-> OH cont.
    L = L +   rct(62,k) * xnew(FORM)  &
    + rct(70,k) * xnew(ALD2)  &
    + rct(73,k) * xnew(ALDX)  &
    + rct(75,k) * xnew(GLYD)  &
    + rct(76,k) * xnew(GLY)  &
    + rct(77,k) * xnew(MGLY)  &
    + rct(78,k) * xnew(H2)  &
    + rct(79,k) * xnew(CO)  &
    + rct(80,k) * xnew(CH4)  &
    + rct(81,k) * xnew(ETHA)  &
    + rct(82,k) * xnew(MEOH)  &
    + rct(83,k) * xnew(ETOH)  &
    + rct(84,k) * xnew(ACET)  &
    + rct(85,k) * xnew(PRPA)  &
    + rct(86,k) * xnew(PAR)  &
    + rct(90,k) * xnew(ETHY)  &
    + rct(92,k) * xnew(ETH)  &
    + rct(96,k) * xnew(OLE)  &
    + rct(100,k) * xnew(IOLE)  &
    + rct(103,k) * xnew(C5H8)
  
  !-> OH cont.
    L = L +   rct(110,k) * xnew(ISPD)  &
    + rct(113,k) * xnew(ISPX)  &
    + rct(115,k) * xnew(EPOX)  &
    + rct(116,k) * xnew(INTR)  &
    + rct(118,k) * xnew(TERP)  &
    + rct(121,k) * xnew(BENZ)  &
    + rct(123,k) * xnew(TOL)  &
    + rct(124,k) * xnew(XYL)  &
    + rct(125,k) * xnew(CRES)  &
    + rct(129,k) * xnew(CRON)  &
    + rct(131,k) * xnew(XOPN)  &
    + rct(134,k) * xnew(OPEN)  &
    + rct(136,k) * xnew(CAT1)  &
    + rct(117,k) * xnew(OPAN)  &
    + rct(133,k) * xnew(PANX)  &
    + rct(80,k) * xnew(ECH4)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(SO2)  &
    + rct(146,k) * xnew(shipNOx)
  
  xnew(OH) = (xold(OH) + dt2 * P) / (1.0 + dt2 * L)


!-> OPO3

  P = 0.457* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.48* rct(126,k) * xnew(CRES) * xnew(NO3)  &
    + rcphot(IDNO2,k) *0.087 * xnew(OPEN)  &
    + 0.6* rct(134,k) * xnew(OPEN) * xnew(OH)  &
    + rct(130,k) * xnew(OPEN) * xnew(NO3)  &
    + rct(44,k) * xnew(OPAN)
  
  L = rct(138,k) * xnew(NO)  &
    + rct(43,k) * xnew(NO2)  &
    + rct(45,k) * xnew(HO2)  &
    + rct(47,k) * xnew(C2O3)  &
    + rct(46,k) * xnew(RO2)
  
  xnew(OPO3) = (xold(OPO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RO2

  P = rct(42,k) * xnew(C2O3) * xnew(NO)  &
    + 0.4* rcphot(IDHO2NO2,k) *0.14 * xnew(PAN)  &
    + 0.44* rct(45,k) * xnew(C2O3) * xnew(HO2)  &
    + 2* rct(47,k) * xnew(C2O3) * xnew(C2O3)  &
    + 2* rct(47,k) * xnew(C2O3) * xnew(CXO3)  &
    + rct(48,k) * xnew(CXO3) * xnew(NO)  &
    + 0.4* rcphot(IDHO2NO2,k) *0.14 * xnew(PANX)  &
    + 0.44* rct(45,k) * xnew(CXO3) * xnew(HO2)  &
    + 0.8* rct(46,k) * xnew(CXO3) * xnew(RO2)  &
    + 2* rct(49,k) * xnew(CXO3) * xnew(CXO3)  &
    + 0.9* rct(55,k) * xnew(MEO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(MEO2) * xnew(RO2)  &
    + 0.8* rct(46,k) * xnew(XO2H) * xnew(C2O3)  &
    + rct(52,k) * xnew(XO2H) * xnew(RO2)  &
    + 0.8* rct(46,k) * xnew(XO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(XO2) * xnew(RO2)  &
    + 0.8* rct(46,k) * xnew(XO2N) * xnew(C2O3)  &
    + rct(52,k) * xnew(XO2N) * xnew(RO2)  &
    + 0.6* rct(58,k) * xnew(CH3OOH) * xnew(OH)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)
  
  !-> RO2 cont.
    P = P +   0.6* rct(58,k) * xnew(ROOH) * xnew(OH)  &
    + rct(61,k) * xnew(AACD) * xnew(OH)  &
    + rcphot(IDCH3CHO,k) * xnew(ALD2)  &
    + rcphot(IDMEK,k) *5.4 * xnew(ALDX)  &
    + 0.11* rcphot(IDACETON,k) *10. * xnew(GLYD)  &
    + 0.2* rct(76,k) * xnew(GLY) * xnew(OH)  &
    + 0.5* rct(71,k) * xnew(GLY) * xnew(NO3)  &
    + rct(71,k) * xnew(MGLY) * xnew(NO3)  &
    + rct(80,k) * xnew(CH4) * xnew(OH)  &
    + rct(81,k) * xnew(ETHA) * xnew(OH)  &
    + 0.1* rct(83,k) * xnew(ETOH) * xnew(OH)  &
    + rcphot(IDACETON,k) *1.5 * xnew(KET)  &
    + 1.38* rcphot(IDACETON,k) * xnew(ACET)  &
    + rct(84,k) * xnew(ACET) * xnew(OH)  &
    + rct(85,k) * xnew(PRPA) * xnew(OH)  &
    + rct(86,k) * xnew(PAR) * xnew(OH)  &
    + 0.98* rct(87,k) * xnew(ROR)  &
    + 0.7* rct(91,k) * xnew(ETH) * xnew(OP)  &
    + rct(92,k) * xnew(ETH) * xnew(OH)  &
    + rct(94,k) * xnew(ETH) * xnew(NO3)
  
  !-> RO2 cont.
    P = P +   0.21* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 1.195* rct(96,k) * xnew(OLE) * xnew(OH)  &
    + 0.15* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 0.1* rct(99,k) * xnew(IOLE) * xnew(OP)  &
    + rct(100,k) * xnew(IOLE) * xnew(OH)  &
    + 0.3* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + rct(102,k) * xnew(IOLE) * xnew(NO3)  &
    + rct(103,k) * xnew(C5H8) * xnew(OH)  &
    + 0.25* rct(104,k) * xnew(C5H8) * xnew(OP)  &
    + 0.082* rct(105,k) * xnew(ISO2) * xnew(NO)  &
    + 0.872* rct(46,k) * xnew(ISO2) * xnew(C2O3)  &
    + 0.072* rct(52,k) * xnew(ISO2) * xnew(RO2)  &
    + 0.2* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + rct(109,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.658* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.284* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + 0.84* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.067* rct(113,k) * xnew(ISPX) * xnew(OH)  &
    + rct(115,k) * xnew(EPOX) * xnew(OH)
  
  !-> RO2 cont.
    P = P +   0.8* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + rct(116,k) * xnew(INTR) * xnew(OH)  &
    + 1.5* rct(118,k) * xnew(TERP) * xnew(OH)  &
    + 0.94* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + 1.28* rct(120,k) * xnew(TERP) * xnew(NO3)  &
    + 0.352* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + rct(46,k) * xnew(BZO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(BZO2) * xnew(RO2)  &
    + 0.72* rct(123,k) * xnew(TOL) * xnew(OH)  &
    + rct(46,k) * xnew(TO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(TO2) * xnew(RO2)  &
    + 0.602* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + rct(46,k) * xnew(XLO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(XLO2) * xnew(RO2)  &
    + 0.02* rct(125,k) * xnew(CRES) * xnew(OH)  &
    + 0.7* rct(126,k) * xnew(CRES) * xnew(NO3)  &
    + 2* rct(131,k) * xnew(XOPN) * xnew(OH)  &
    + 0.3* rct(132,k) * xnew(XOPN) * xnew(O3)  &
    + rct(133,k) * xnew(XOPN) * xnew(NO3)
  
  !-> RO2 cont.
    P = P +   0.4* rct(134,k) * xnew(OPEN) * xnew(OH)  &
    + 0.44* rct(45,k) * xnew(OPO3) * xnew(HO2)  &
    + 2* rct(47,k) * xnew(OPO3) * xnew(C2O3)  &
    + 1.8* rct(46,k) * xnew(OPO3) * xnew(RO2)  &
    + rct(80,k) * xnew(ECH4) * xnew(OH)
  
  L = rct(46,k) * xnew(C2O3)  &
    + rct(46,k) * xnew(CXO3)  &
    + rct(50,k) * xnew(NO)  &
    + rct(51,k) * xnew(HO2)  &
    + rct(52,k) * xnew(RO2)  &
    + rct(52,k) * xnew(RO2)  &
    + rct(52,k) * xnew(MEO2)  &
    + rct(52,k) * xnew(XO2H)  &
    + rct(52,k) * xnew(XO2)  &
    + rct(52,k) * xnew(XO2N)  &
    + rct(52,k) * xnew(ISO2)  &
    + rct(52,k) * xnew(EPX2)  &
    + rct(52,k) * xnew(BZO2)  &
    + rct(52,k) * xnew(TO2)  &
    + rct(52,k) * xnew(XLO2)  &
    + rct(46,k) * xnew(OPO3)
  
  xnew(RO2) = (xold(RO2) + dt2 * P) / (1.0 + dt2 * L)


!-> ROR

  P = 0.76* rct(86,k) * xnew(PAR) * xnew(OH)  &
    + 0.02* rct(87,k) * xnew(ROR)
  
  L = rct(87,k)  &
    + rct(88,k)  &
    + rct(89,k) * xnew(NO2)
  
  xnew(ROR) = (xold(ROR) + dt2 * P) / (1.0 + dt2 * L)


!-> TO2

  P = 0.65* rct(123,k) * xnew(TOL) * xnew(OH)
  
  L = rct(56,k) * xnew(NO)  &
    + rct(46,k) * xnew(C2O3)  &
    + rct(122,k) * xnew(HO2)  &
    + rct(52,k) * xnew(RO2)
  
  xnew(TO2) = (xold(TO2) + dt2 * P) / (1.0 + dt2 * L)


!-> XLO2

  P = 0.544* rct(124,k) * xnew(XYL) * xnew(OH)
  
  L = rct(56,k) * xnew(NO)  &
    + rct(122,k) * xnew(HO2)  &
    + rct(46,k) * xnew(C2O3)  &
    + rct(52,k) * xnew(RO2)
  
  xnew(XLO2) = (xold(XLO2) + dt2 * P) / (1.0 + dt2 * L)


!-> XO2

  P = 0.2* rct(76,k) * xnew(GLY) * xnew(OH)  &
    + 0.5* rct(71,k) * xnew(GLY) * xnew(NO3)  &
    + rct(71,k) * xnew(MGLY) * xnew(NO3)  &
    + rct(84,k) * xnew(ACET) * xnew(OH)  &
    + 0.76* rct(86,k) * xnew(PAR) * xnew(OH)  &
    + 0.5* rct(94,k) * xnew(ETH) * xnew(NO3)  &
    + 0.195* rct(96,k) * xnew(OLE) * xnew(OH)  &
    + 0.48* rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 0.48* rct(102,k) * xnew(IOLE) * xnew(NO3)  &
    + 0.25* rct(104,k) * xnew(C5H8) * xnew(OP)  &
    + 0.2* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.33* rct(109,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.521* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.142* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + 0.16* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.63* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + 0.5* rct(118,k) * xnew(TERP) * xnew(OH)  &
    + 0.69* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + 0.75* rct(120,k) * xnew(TERP) * xnew(NO3)  &
    + 0.48* rct(126,k) * xnew(CRES) * xnew(NO3)
  
  !-> XO2 cont.
    P = P +   0.45* rct(133,k) * xnew(XOPN) * xnew(NO3)  &
    + rct(47,k) * xnew(OPO3) * xnew(C2O3)
  
  L = rct(56,k) * xnew(NO)  &
    + rct(57,k) * xnew(HO2)  &
    + rct(46,k) * xnew(C2O3)  &
    + rct(52,k) * xnew(RO2)
  
  xnew(XO2) = (xold(XO2) + dt2 * P) / (1.0 + dt2 * L)


!-> XO2H

  P = rct(47,k) * xnew(C2O3) * xnew(CXO3)  &
    + rct(48,k) * xnew(CXO3) * xnew(NO)  &
    + 0.4* rcphot(IDHO2NO2,k) *0.14 * xnew(PANX)  &
    + 0.44* rct(45,k) * xnew(CXO3) * xnew(HO2)  &
    + 0.8* rct(46,k) * xnew(CXO3) * xnew(RO2)  &
    + 2* rct(49,k) * xnew(CXO3) * xnew(CXO3)  &
    + 0.54* rct(58,k) * xnew(ROOH) * xnew(OH)  &
    + rcphot(IDMEK,k) *5.4 * xnew(ALDX)  &
    + 0.11* rcphot(IDACETON,k) *10. * xnew(GLYD)  &
    + 0.991* rct(81,k) * xnew(ETHA) * xnew(OH)  &
    + 0.1* rct(83,k) * xnew(ETOH) * xnew(OH)  &
    + 0.5* rcphot(IDACETON,k) *1.5 * xnew(KET)  &
    + 0.97* rct(85,k) * xnew(PRPA) * xnew(OH)  &
    + 0.11* rct(86,k) * xnew(PAR) * xnew(OH)  &
    + 0.94* rct(87,k) * xnew(ROR)  &
    + 0.7* rct(91,k) * xnew(ETH) * xnew(OP)  &
    + rct(92,k) * xnew(ETH) * xnew(OH)  &
    + 0.5* rct(94,k) * xnew(ETH) * xnew(NO3)  &
    + 0.2* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.976* rct(96,k) * xnew(OLE) * xnew(OH)
  
  !-> XO2H cont.
    P = P +   0.15* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.48* rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 0.1* rct(99,k) * xnew(IOLE) * xnew(OP)  &
    + rct(100,k) * xnew(IOLE) * xnew(OH)  &
    + 0.3* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.48* rct(102,k) * xnew(IOLE) * xnew(NO3)  &
    + 0.082* rct(105,k) * xnew(ISO2) * xnew(NO)  &
    + 0.072* rct(46,k) * xnew(ISO2) * xnew(C2O3)  &
    + 0.072* rct(52,k) * xnew(ISO2) * xnew(RO2)  &
    + 0.64* rct(109,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.142* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + 0.34* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.37* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + 0.75* rct(118,k) * xnew(TERP) * xnew(OH)  &
    + 0.07* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + 0.28* rct(120,k) * xnew(TERP) * xnew(NO3)  &
    + 0.07* rct(123,k) * xnew(TOL) * xnew(OH)  &
    + 0.058* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + 0.12* rct(126,k) * xnew(CRES) * xnew(NO3)  &
    + rcphot(IDNO2,k) *0.087 * xnew(XOPN)
  
  !-> XO2H cont.
    P = P +   2* rct(131,k) * xnew(XOPN) * xnew(OH)  &
    + 0.3* rct(132,k) * xnew(XOPN) * xnew(O3)  &
    + 0.45* rct(133,k) * xnew(XOPN) * xnew(NO3)  &
    + 0.4* rct(134,k) * xnew(OPEN) * xnew(OH)  &
    + 0.44* rct(45,k) * xnew(OPO3) * xnew(HO2)  &
    + 0.8* rct(46,k) * xnew(OPO3) * xnew(RO2)
  
  L = rct(56,k) * xnew(NO)  &
    + rct(57,k) * xnew(HO2)  &
    + rct(46,k) * xnew(C2O3)  &
    + rct(52,k) * xnew(RO2)
  
  xnew(XO2H) = (xold(XO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> XO2N

  P = 0.06* rct(58,k) * xnew(ROOH) * xnew(OH)  &
    + 0.009* rct(81,k) * xnew(ETHA) * xnew(OH)  &
    + 0.03* rct(85,k) * xnew(PRPA) * xnew(OH)  &
    + 0.13* rct(86,k) * xnew(PAR) * xnew(OH)  &
    + 0.04* rct(87,k) * xnew(ROR)  &
    + 0.01* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.024* rct(96,k) * xnew(OLE) * xnew(OH)  &
    + 0.04* rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 0.04* rct(102,k) * xnew(IOLE) * xnew(NO3)  &
    + 0.03* rct(109,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.022* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.25* rct(118,k) * xnew(TERP) * xnew(OH)  &
    + 0.18* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + 0.25* rct(120,k) * xnew(TERP) * xnew(NO3)  &
    + 0.02* rct(125,k) * xnew(CRES) * xnew(OH)  &
    + 0.1* rct(126,k) * xnew(CRES) * xnew(NO3)  &
    + 0.1* rct(133,k) * xnew(XOPN) * xnew(NO3)
  
  L = rct(56,k) * xnew(NO)  &
    + rct(57,k) * xnew(HO2)  &
    + rct(46,k) * xnew(C2O3)  &
    + rct(52,k) * xnew(RO2)
  
  xnew(XO2N) = (xold(XO2N) + dt2 * P) / (1.0 + dt2 * L)


!-> NO3

  P = rct(5,k) * xnew(OP) * xnew(NO2)  &
    + rct(22,k) * xnew(NO2) * xnew(O3)  &
    + rct(31,k) * xnew(N2O5)  &
    + rcphot(IDN2O5,k) * xnew(N2O5)  &
    + rct(38,k) * xnew(HNO3) * xnew(OH)  &
    + 0.41* rcphot(IDHO2NO2,k) * xnew(PNA)  &
    + 0.4* rcphot(IDHO2NO2,k) *0.14 * xnew(PAN)  &
    + 0.4* rcphot(IDHO2NO2,k) *0.14 * xnew(PANX)  &
    + 0.185* rct(116,k) * xnew(INTR) * xnew(OH)
  
  L = rcphot(IDNO3,k)  &
    + rct(23,k) * xnew(NO)  &
    + rct(24,k) * xnew(NO2)  &
    + rct(25,k) * xnew(OP)  &
    + rct(26,k) * xnew(OH)  &
    + rct(27,k) * xnew(HO2)  &
    + rct(28,k) * xnew(O3)  &
    + rct(29,k) * xnew(NO3)  &
    + rct(29,k) * xnew(NO3)  &
    + rct(30,k) * xnew(NO2)  &
    + rct(64,k) * xnew(FORM)  &
    + rct(71,k) * xnew(ALD2)  &
    + rct(74,k) * xnew(ALDX)  &
    + rct(71,k) * xnew(GLYD)  &
    + rct(71,k) * xnew(GLY)  &
    + rct(71,k) * xnew(MGLY)  &
    + rct(94,k) * xnew(ETH)  &
    + rct(98,k) * xnew(OLE)  &
    + rct(102,k) * xnew(IOLE)  &
    + rct(109,k) * xnew(C5H8)
  
  !-> NO3 cont.
    L = L +   rct(112,k) * xnew(ISPD)  &
    + rct(114,k) * xnew(HPLD)  &
    + rct(120,k) * xnew(TERP)  &
    + rct(126,k) * xnew(CRES)  &
    + rct(130,k) * xnew(CRON)  &
    + rct(133,k) * xnew(XOPN)  &
    + rct(130,k) * xnew(OPEN)  &
    + rct(137,k) * xnew(CAT1)  &
    + rct(142,k)
  
  xnew(NO3) = (xold(NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> AACD

  P = 0.15* rct(45,k) * xnew(C2O3) * xnew(HO2)  &
    + 0.15* rct(45,k) * xnew(CXO3) * xnew(HO2)  &
    + 0.1* rct(55,k) * xnew(MEO2) * xnew(C2O3)  &
    + 0.2* rct(46,k) * xnew(XO2H) * xnew(C2O3)  &
    + 0.2* rct(46,k) * xnew(XO2) * xnew(C2O3)  &
    + 0.2* rct(46,k) * xnew(XO2N) * xnew(C2O3)  &
    + 0.13* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.08* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.2* rct(46,k) * xnew(ISO2) * xnew(C2O3)  &
    + 0.2* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 0.15* rct(45,k) * xnew(OPO3) * xnew(HO2)  &
    + 0.2* rct(46,k) * xnew(OPO3) * xnew(RO2)  &
    + rcemis(AACD,k)
  
  L = rct(61,k) * xnew(OH)
  
  xnew(AACD) = (xold(AACD) + dt2 * P) / (1.0 + dt2 * L)


!-> ACET

  P = 0.71* rct(85,k) * xnew(PRPA) * xnew(OH)  &
    + 0.42* rct(87,k) * xnew(ROR)  &
    + 0.137* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.17* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.17* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + rcemis(ACET,k)
  
  L = rcphot(IDACETON,k)  &
    + rct(84,k) * xnew(OH)
  
  xnew(ACET) = (xold(ACET) + dt2 * P) / (1.0 + dt2 * L)


!-> ALD2

  P = rct(47,k) * xnew(C2O3) * xnew(CXO3)  &
    + rct(48,k) * xnew(CXO3) * xnew(NO)  &
    + 0.4* rcphot(IDHO2NO2,k) *0.14 * xnew(PANX)  &
    + 0.44* rct(45,k) * xnew(CXO3) * xnew(HO2)  &
    + 0.8* rct(46,k) * xnew(CXO3) * xnew(RO2)  &
    + 2* rct(49,k) * xnew(CXO3) * xnew(CXO3)  &
    + rcphot(IDMEK,k) *5.4 * xnew(ALDX)  &
    + 0.991* rct(81,k) * xnew(ETHA) * xnew(OH)  &
    + 0.95* rct(83,k) * xnew(ETOH) * xnew(OH)  &
    + 0.5* rcphot(IDACETON,k) *1.5 * xnew(KET)  &
    + 0.74* rct(87,k) * xnew(ROR)  &
    + 0.2* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.488* rct(96,k) * xnew(OLE) * xnew(OH)  &
    + 0.295* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.25* rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 1.24* rct(99,k) * xnew(IOLE) * xnew(OP)  &
    + 1.3* rct(100,k) * xnew(IOLE) * xnew(OH)  &
    + 0.732* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.5* rct(102,k) * xnew(IOLE) * xnew(NO3)  &
    + 0.04* rct(111,k) * xnew(ISPD) * xnew(O3)
  
  !-> ALD2 cont.
    P = P +   0.1* rct(132,k) * xnew(XOPN) * xnew(O3)  &
    + 0.02* rct(135,k) * xnew(OPEN) * xnew(O3)  &
    + rct(133,k) * xnew(PANX) * xnew(OH)  &
    + rcemis(ALD2,k)
  
  L = rct(69,k) * xnew(OP)  &
    + rct(70,k) * xnew(OH)  &
    + rct(71,k) * xnew(NO3)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(ALD2) = (xold(ALD2) + dt2 * P) / (1.0 + dt2 * L)


!-> ALDX

  P = 0.26* rct(85,k) * xnew(PRPA) * xnew(OH)  &
    + 0.11* rct(86,k) * xnew(PAR) * xnew(OH)  &
    + 0.37* rct(87,k) * xnew(ROR)  &
    + 0.3* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.488* rct(96,k) * xnew(OLE) * xnew(OH)  &
    + 0.27* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.375* rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 0.66* rct(99,k) * xnew(IOLE) * xnew(OP)  &
    + 0.7* rct(100,k) * xnew(IOLE) * xnew(OH)  &
    + 0.442* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.625* rct(102,k) * xnew(IOLE) * xnew(NO3)  &
    + 0.15* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.029* rct(113,k) * xnew(ISPX) * xnew(OH)  &
    + 0.078* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + 0.15* rct(117,k) * xnew(TERP) * xnew(OP)  &
    + 0.47* rct(118,k) * xnew(TERP) * xnew(OH)  &
    + 0.21* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + 0.47* rct(120,k) * xnew(TERP) * xnew(NO3)  &
    + 0.44* rct(45,k) * xnew(OPO3) * xnew(HO2)  &
    + rct(47,k) * xnew(OPO3) * xnew(C2O3)
  
  !-> ALDX cont.
    P = P +   0.8* rct(46,k) * xnew(OPO3) * xnew(RO2)  &
    + rcemis(ALDX,k)
  
  L = rct(72,k) * xnew(OP)  &
    + rct(73,k) * xnew(OH)  &
    + rct(74,k) * xnew(NO3)  &
    + rcphot(IDMEK,k) *5.4
  
  xnew(ALDX) = (xold(ALDX) + dt2 * P) / (1.0 + dt2 * L)


!-> BENZ

  P = rcemis(BENZ,k)
  
  L = rct(121,k) * xnew(OH)
  
  xnew(BENZ) = (xold(BENZ) + dt2 * P) / (1.0 + dt2 * L)


!-> CAT1

  P = 0.732* rct(125,k) * xnew(CRES) * xnew(OH)
  
  L = rct(136,k) * xnew(OH)  &
    + rct(137,k) * xnew(NO3)
  
  xnew(CAT1) = (xold(CAT1) + dt2 * P) / (1.0 + dt2 * L)


!-> CO

  P = rct(62,k) * xnew(FORM) * xnew(OH)  &
    + rcphot(IDHCHO_H,k) * xnew(FORM)  &
    + rcphot(IDHCHO_H2,k) * xnew(FORM)  &
    + rct(63,k) * xnew(FORM) * xnew(OP)  &
    + rct(64,k) * xnew(FORM) * xnew(NO3)  &
    + rcphot(IDCH3CHO,k) * xnew(ALD2)  &
    + rcphot(IDMEK,k) *5.4 * xnew(ALDX)  &
    + 0.89* rcphot(IDACETON,k) *10. * xnew(GLYD)  &
    + 1.8* rct(76,k) * xnew(GLY) * xnew(OH)  &
    + 2* rcphot(IDCHOCHO,k) * xnew(GLY)  &
    + 1.5* rct(71,k) * xnew(GLY) * xnew(NO3)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLY)  &
    + rct(77,k) * xnew(MGLY) * xnew(OH)  &
    + 0.38* rcphot(IDACETON,k) * xnew(ACET)  &
    + 0.3* rct(90,k) * xnew(ETHY) * xnew(OH)  &
    + rct(91,k) * xnew(ETH) * xnew(OP)  &
    + 0.51* rct(93,k) * xnew(ETH) * xnew(O3)  &
    + 0.2* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.378* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.1* rct(99,k) * xnew(IOLE) * xnew(OP)
  
  !-> CO cont.
    P = P +   0.245* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.066* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.137* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.543* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.251* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 0.251* rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 0.2* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 0.251* rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + 0.001* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + 0.7* rcphot(IDNO2,k) *0.087 * xnew(XOPN)  &
    + 0.5* rct(132,k) * xnew(XOPN) * xnew(O3)  &
    + rcphot(IDNO2,k) *0.087 * xnew(OPEN)  &
    + 1.98* rct(135,k) * xnew(OPEN) * xnew(O3)  &
    + 0.5* rct(138,k) * xnew(OPO3) * xnew(NO)  &
    + rct(117,k) * xnew(OPAN) * xnew(OH)  &
    + rcemis(CO,k)
  
  L = rct(79,k) * xnew(OH)
  
  xnew(CO) = (xold(CO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH4

  P = rcemis(CH4,k)
  
  L = rct(80,k) * xnew(OH)
  
  xnew(CH4) = (xold(CH4) + dt2 * P) / (1.0 + dt2 * L)


!-> CRES

  P = 0.53* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + 0.18* rct(123,k) * xnew(TOL) * xnew(OH)  &
    + 0.155* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + rct(128,k) * xnew(CRO) * xnew(HO2)
  
  L = rct(125,k) * xnew(OH)  &
    + rct(126,k) * xnew(NO3)
  
  xnew(CRES) = (xold(CRES) + dt2 * P) / (1.0 + dt2 * L)


!-> CRON

  P = rct(127,k) * xnew(CRO) * xnew(NO2)
  
  L = rct(129,k) * xnew(OH)  &
    + rct(130,k) * xnew(NO3)  &
    + rcphot(IDHONO,k) *0.076
  
  xnew(CRON) = (xold(CRON) + dt2 * P) / (1.0 + dt2 * L)


!-> EPOX

  P = 0.904* rct(113,k) * xnew(ISPX) * xnew(OH)
  
  L = rct(115,k) * xnew(OH)
  
  xnew(EPOX) = (xold(EPOX) + dt2 * P) / (1.0 + dt2 * L)


!-> ETH

  P = rcemis(ETH,k)
  
  L = rct(91,k) * xnew(OP)  &
    + rct(92,k) * xnew(OH)  &
    + rct(93,k) * xnew(O3)  &
    + rct(94,k) * xnew(NO3)
  
  xnew(ETH) = (xold(ETH) + dt2 * P) / (1.0 + dt2 * L)


!-> ETHA

  P = rcemis(ETHA,k)
  
  L = rct(81,k) * xnew(OH)
  
  xnew(ETHA) = (xold(ETHA) + dt2 * P) / (1.0 + dt2 * L)


!-> ETHY

  P = rcemis(ETHY,k)
  
  L = rct(90,k) * xnew(OH)
  
  xnew(ETHY) = (xold(ETHY) + dt2 * P) / (1.0 + dt2 * L)


!-> ETOH

  P = rcemis(ETOH,k)
  
  L = rct(83,k) * xnew(OH)
  
  xnew(ETOH) = (xold(ETOH) + dt2 * P) / (1.0 + dt2 * L)


!-> FACD

  P = rct(67,k) * xnew(HCO3) * xnew(NO)  &
    + 0.5* rct(68,k) * xnew(HCO3) * xnew(HO2)  &
    + 0.3* rct(90,k) * xnew(ETHY) * xnew(OH)  &
    + 0.37* rct(93,k) * xnew(ETH) * xnew(O3)  &
    + 0.09* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.15* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.074* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 0.185* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + rcemis(FACD,k)
  
  L = rct(60,k) * xnew(OH)
  
  xnew(FACD) = (xold(FACD) + dt2 * P) / (1.0 + dt2 * L)


!-> FORM

  P = rct(53,k) * xnew(MEO2) * xnew(NO)  &
    + 0.1* rct(54,k) * xnew(MEO2) * xnew(HO2)  &
    + rct(55,k) * xnew(MEO2) * xnew(C2O3)  &
    + 0.685* rct(52,k) * xnew(MEO2) * xnew(RO2)  &
    + 0.4* rct(58,k) * xnew(CH3OOH) * xnew(OH)  &
    + rct(66,k) * xnew(HCO3)  &
    + 0.74* rcphot(IDACETON,k) *10. * xnew(GLYD)  &
    + rct(82,k) * xnew(MEOH) * xnew(OH)  &
    + 0.078* rct(83,k) * xnew(ETOH) * xnew(OH)  &
    + rct(84,k) * xnew(ACET) * xnew(OH)  &
    + rct(91,k) * xnew(ETH) * xnew(OP)  &
    + 1.56* rct(92,k) * xnew(ETH) * xnew(OH)  &
    + rct(93,k) * xnew(ETH) * xnew(O3)  &
    + 1.125* rct(94,k) * xnew(ETH) * xnew(NO3)  &
    + 0.2* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.781* rct(96,k) * xnew(OLE) * xnew(OH)  &
    + 0.555* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.5* rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 0.128* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.5* rct(104,k) * xnew(C5H8) * xnew(OP)
  
  !-> FORM cont.
    P = P +   0.673* rct(105,k) * xnew(ISO2) * xnew(NO)  &
    + 0.12* rct(106,k) * xnew(ISO2) * xnew(HO2)  &
    + 0.598* rct(46,k) * xnew(ISO2) * xnew(C2O3)  &
    + 0.598* rct(52,k) * xnew(ISO2) * xnew(RO2)  &
    + 0.6* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.35* rct(109,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.231* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.26* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.375* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 0.375* rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 0.3* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 0.375* rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + 0.592* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + 0.28* rct(118,k) * xnew(TERP) * xnew(OH)  &
    + 0.24* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + rcphot(IDHONO,k) *0.076 * xnew(CRON)  &
    + 0.08* rct(135,k) * xnew(OPEN) * xnew(O3)  &
    + 0.14* rct(136,k) * xnew(CAT1) * xnew(OH)  &
    + rcemis(FORM,k)
  
  L = rct(62,k) * xnew(OH)  &
    + rcphot(IDHCHO_H,k)  &
    + rcphot(IDHCHO_H2,k)  &
    + rct(63,k) * xnew(OP)  &
    + rct(64,k) * xnew(NO3)  &
    + rct(65,k) * xnew(HO2)
  
  xnew(FORM) = (xold(FORM) + dt2 * P) / (1.0 + dt2 * L)


!-> GLY

  P = 0.2* rct(75,k) * xnew(GLYD) * xnew(OH)  &
    + 0.11* rcphot(IDACETON,k) *10. * xnew(GLYD)  &
    + 0.7* rct(90,k) * xnew(ETHY) * xnew(OH)  &
    + 0.075* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.24* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.17* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.275* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 0.275* rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 0.22* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 0.275* rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + 0.918* rct(56,k) * xnew(BZO2) * xnew(NO)  &
    + rct(46,k) * xnew(BZO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(BZO2) * xnew(RO2)  &
    + 0.417* rct(56,k) * xnew(TO2) * xnew(NO)  &
    + 0.48* rct(46,k) * xnew(TO2) * xnew(C2O3)  &
    + 0.48* rct(52,k) * xnew(TO2) * xnew(RO2)  &
    + 0.221* rct(56,k) * xnew(XLO2) * xnew(NO)  &
    + 0.26* rct(46,k) * xnew(XLO2) * xnew(C2O3)  &
    + 0.26* rct(52,k) * xnew(XLO2) * xnew(RO2)  &
    + 0.025* rct(125,k) * xnew(CRES) * xnew(OH)
  
  !-> GLY cont.
    P = P +   0.24* rct(126,k) * xnew(CRES) * xnew(NO3)  &
    + 0.4* rcphot(IDNO2,k) *0.087 * xnew(XOPN)  &
    + 0.4* rct(131,k) * xnew(XOPN) * xnew(OH)  &
    + 0.4* rct(134,k) * xnew(OPEN) * xnew(OH)  &
    + 1.4* rct(135,k) * xnew(OPEN) * xnew(O3)  &
    + 0.5* rct(138,k) * xnew(OPO3) * xnew(NO)  &
    + 0.5* rct(117,k) * xnew(OPAN) * xnew(OH)  &
    + rcemis(GLY,k)
  
  L = rct(76,k) * xnew(OH)  &
    + rcphot(IDCHOCHO,k)  &
    + rct(71,k) * xnew(NO3)
  
  xnew(GLY) = (xold(GLY) + dt2 * P) / (1.0 + dt2 * L)


!-> GLYD

  P = 0.011* rct(83,k) * xnew(ETOH) * xnew(OH)  &
    + 0.22* rct(92,k) * xnew(ETH) * xnew(OH)  &
    + 0.269* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.113* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + 0.128* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.275* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 0.275* rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 0.22* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 0.275* rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + 0.331* rct(116,k) * xnew(INTR) * xnew(OH)
  
  L = rct(75,k) * xnew(OH)  &
    + rcphot(IDACETON,k) *10.  &
    + rct(71,k) * xnew(NO3)
  
  xnew(GLYD) = (xold(GLYD) + dt2 * P) / (1.0 + dt2 * L)


!-> H2O2

  P = rct(14,k) * xnew(OH) * xnew(OH)  &
    + rct(16,k) * xnew(HO2) * xnew(HO2)  &
    + rct(17,k) * xnew(HO2) * xnew(HO2)  &
    + 0.04* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.08* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.5* rct(145,k) * xnew(HO2)
  
  L = rcphot(IDH2O2,k)  &
    + rct(18,k) * xnew(OH)  &
    + rct(19,k) * xnew(OP)  &
    + AQRCK(ICLRC1,K) * xnew(SO2)  &
    + 1.0e-5
  
  xnew(H2O2) = (xold(H2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> HNO3

  P = 2* rct(32,k) * xnew(N2O5)  &
    + rct(37,k) * xnew(NO2) * xnew(OH)  &
    + rct(64,k) * xnew(FORM) * xnew(NO3)  &
    + rct(71,k) * xnew(ALD2) * xnew(NO3)  &
    + rct(74,k) * xnew(ALDX) * xnew(NO3)  &
    + rct(71,k) * xnew(GLYD) * xnew(NO3)  &
    + rct(71,k) * xnew(GLY) * xnew(NO3)  &
    + rct(71,k) * xnew(MGLY) * xnew(NO3)  &
    + 0.717* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + rct(114,k) * xnew(HPLD) * xnew(NO3)  &
    + rct(126,k) * xnew(CRES) * xnew(NO3)  &
    + rct(130,k) * xnew(CRON) * xnew(NO3)  &
    + rct(130,k) * xnew(OPEN) * xnew(NO3)  &
    + rct(137,k) * xnew(CAT1) * xnew(NO3)  &
    + rct(139,k) * xnew(NTR2)  &
    + 2* rct(141,k) * xnew(N2O5)  &
    + rct(142,k) * xnew(NO3)  &
    + rct(146,k) * xnew(shipNOx) * xnew(OH)  &
    + 3.2e-5 * xnew(shipNOx)
  
  L = rct(38,k) * xnew(OH)  &
    + rcphot(IDHNO3,k)  &
    + rct(143,k)  &
    + rct(144,k)
  
  xnew(HNO3) = (xold(HNO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HONO

  P = rct(33,k) * xnew(NO) * xnew(OH)  &
    + 2* rct(34,k) * xnew(NO) * xnew(NO2)  &
    + rcphot(IDHONO,k) *0.076 * xnew(CRON)
  
  L = rct(35,k) * xnew(HONO)  &
    + rct(35,k) * xnew(HONO)  &
    + rcphot(IDHONO,k)  &
    + rct(36,k) * xnew(OH)
  
  xnew(HONO) = (xold(HONO) + dt2 * P) / (1.0 + dt2 * L)


!-> HPLD

  P = rct(107,k) * xnew(ISO2)
  
  L = rcphot(IDCHOCHO,k) *7.2  &
    + rct(114,k) * xnew(NO3)
  
  xnew(HPLD) = (xold(HPLD) + dt2 * P) / (1.0 + dt2 * L)


!-> INTR

  P = 0.1* rct(105,k) * xnew(ISO2) * xnew(NO)  &
    + 0.104* rct(116,k) * xnew(INTR) * xnew(OH)
  
  L = rct(116,k) * xnew(OH)
  
  xnew(INTR) = (xold(INTR) + dt2 * P) / (1.0 + dt2 * L)


!-> IOLE

  P = 0.029* rct(113,k) * xnew(ISPX) * xnew(OH)  &
    + rcemis(IOLE,k)
  
  L = rct(99,k) * xnew(OP)  &
    + rct(100,k) * xnew(OH)  &
    + rct(101,k) * xnew(O3)  &
    + rct(102,k) * xnew(NO3)
  
  xnew(IOLE) = (xold(IOLE) + dt2 * P) / (1.0 + dt2 * L)


!-> C5H8

  P = rcemis(C5H8,k)  &
    + RCBIO(1,K)
  
  L = rct(103,k) * xnew(OH)  &
    + rct(104,k) * xnew(OP)  &
    + rct(108,k) * xnew(O3)  &
    + rct(109,k) * xnew(NO3)
  
  xnew(C5H8) = (xold(C5H8) + dt2 * P) / (1.0 + dt2 * L)


!-> ISPD

  P = 0.75* rct(104,k) * xnew(C5H8) * xnew(OP)  &
    + 0.9* rct(105,k) * xnew(ISO2) * xnew(NO)  &
    + 0.12* rct(106,k) * xnew(ISO2) * xnew(HO2)  &
    + rct(46,k) * xnew(ISO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(ISO2) * xnew(RO2)  &
    + 0.65* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.35* rct(109,k) * xnew(C5H8) * xnew(NO3)  &
    + rcphot(IDCHOCHO,k) *7.2 * xnew(HPLD)  &
    + rct(114,k) * xnew(HPLD) * xnew(NO3)
  
  L = rct(110,k) * xnew(OH)  &
    + rct(111,k) * xnew(O3)  &
    + rct(112,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H2,k) *0.62
  
  xnew(ISPD) = (xold(ISPD) + dt2 * P) / (1.0 + dt2 * L)


!-> ISPX

  P = 0.88* rct(106,k) * xnew(ISO2) * xnew(HO2)
  
  L = rct(113,k) * xnew(OH)
  
  xnew(ISPX) = (xold(ISPX) + dt2 * P) / (1.0 + dt2 * L)


!-> KET

  P = 0.2* rct(87,k) * xnew(ROR)  &
    + rct(88,k) * xnew(ROR)  &
    + rcemis(KET,k)
  
  L = rcphot(IDACETON,k) *1.5
  
  xnew(KET) = (xold(KET) + dt2 * P) / (1.0 + dt2 * L)


!-> MEOH

  P = 0.315* rct(52,k) * xnew(MEO2) * xnew(RO2)  &
    + 0.15* rcphot(IDACETON,k) *10. * xnew(GLYD)  &
    + rcemis(MEOH,k)
  
  L = rct(82,k) * xnew(OH)
  
  xnew(MEOH) = (xold(MEOH) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OOH

  P = 0.9* rct(54,k) * xnew(MEO2) * xnew(HO2)  &
    + 0.5* rct(68,k) * xnew(HCO3) * xnew(HO2)
  
  L = rct(58,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + 1.0e-5
  
  xnew(CH3OOH) = (xold(CH3OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> MGLY

  P = 0.075* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + 0.06* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + 0.115* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.531* rct(111,k) * xnew(ISPD) * xnew(O3)  &
    + 0.113* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + 0.275* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 0.275* rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 0.22* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 0.275* rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + 0.443* rct(56,k) * xnew(TO2) * xnew(NO)  &
    + 0.52* rct(46,k) * xnew(TO2) * xnew(C2O3)  &
    + 0.52* rct(52,k) * xnew(TO2) * xnew(RO2)  &
    + 0.675* rct(56,k) * xnew(XLO2) * xnew(NO)  &
    + 0.77* rct(46,k) * xnew(XLO2) * xnew(C2O3)  &
    + 0.77* rct(52,k) * xnew(XLO2) * xnew(RO2)  &
    + 0.24* rct(126,k) * xnew(CRES) * xnew(NO3)  &
    + rct(131,k) * xnew(XOPN) * xnew(OH)  &
    + 1.2* rct(132,k) * xnew(XOPN) * xnew(O3)  &
    + 0.25* rct(133,k) * xnew(XOPN) * xnew(NO3)  &
    + 0.24* rct(135,k) * xnew(OPEN) * xnew(O3)
  
  !-> MGLY cont.
    P = P +   rcemis(MGLY,k)
  
  L = rcphot(IDRCOCHO,k)  &
    + rct(71,k) * xnew(NO3)  &
    + rct(77,k) * xnew(OH)
  
  xnew(MGLY) = (xold(MGLY) + dt2 * P) / (1.0 + dt2 * L)


!-> N2O5

  P = rct(30,k) * xnew(NO3) * xnew(NO2)
  
  L = rct(31,k)  &
    + rcphot(IDN2O5,k)  &
    + rct(32,k)  &
    + rct(141,k)
  
  xnew(N2O5) = (xold(N2O5) + dt2 * P) / (1.0 + dt2 * L)


!-> NO

  P = rcphot(IDNO2,k) * xnew(NO2)  &
    + rct(4,k) * xnew(OP) * xnew(NO2)  &
    + 0.113* rcphot(IDNO3,k) * xnew(NO3)  &
    + rct(24,k) * xnew(NO3) * xnew(NO2)  &
    + rct(35,k) * xnew(HONO) * xnew(HONO)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rct(50,k) * xnew(RO2) * xnew(NO)  &
    + rcemis(NO,k)
  
  L = rct(2,k) * xnew(O3)  &
    + rct(3,k) * xnew(OP)  &
    + rct(20,k) * xnew(NO)  &
    + rct(20,k) * xnew(NO)  &
    + rct(21,k) * xnew(HO2)  &
    + rct(23,k) * xnew(NO3)  &
    + rct(33,k) * xnew(OH)  &
    + rct(34,k) * xnew(NO2)  &
    + rct(42,k) * xnew(C2O3)  &
    + rct(48,k) * xnew(CXO3)  &
    + rct(50,k) * xnew(RO2)  &
    + rct(53,k) * xnew(MEO2)  &
    + rct(56,k) * xnew(XO2H)  &
    + rct(56,k) * xnew(XO2)  &
    + rct(56,k) * xnew(XO2N)  &
    + rct(67,k) * xnew(HCO3)  &
    + rct(105,k) * xnew(ISO2)  &
    + rct(105,k) * xnew(EPX2)  &
    + rct(56,k) * xnew(BZO2)  &
    + rct(56,k) * xnew(TO2)
  
  !-> NO cont.
    L = L +   rct(56,k) * xnew(XLO2)  &
    + rct(138,k) * xnew(OPO3)
  
  xnew(NO) = (xold(NO) + dt2 * P) / (1.0 + dt2 * L)


!-> NO2

  P = rct(2,k) * xnew(O3) * xnew(NO)  &
    + rct(3,k) * xnew(OP) * xnew(NO)  &
    + 2* rct(20,k) * xnew(NO) * xnew(NO)  &
    + rct(21,k) * xnew(HO2) * xnew(NO)  &
    + 0.887* rcphot(IDNO3,k) * xnew(NO3)  &
    + 2* rct(23,k) * xnew(NO3) * xnew(NO)  &
    + rct(24,k) * xnew(NO3) * xnew(NO2)  &
    + rct(25,k) * xnew(NO3) * xnew(OP)  &
    + rct(26,k) * xnew(NO3) * xnew(OH)  &
    + rct(27,k) * xnew(NO3) * xnew(HO2)  &
    + rct(28,k) * xnew(NO3) * xnew(O3)  &
    + 2* rct(29,k) * xnew(NO3) * xnew(NO3)  &
    + rct(31,k) * xnew(N2O5)  &
    + rcphot(IDN2O5,k) * xnew(N2O5)  &
    + rct(35,k) * xnew(HONO) * xnew(HONO)  &
    + rct(36,k) * xnew(HONO) * xnew(OH)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + rct(40,k) * xnew(PNA)  &
    + 0.59* rcphot(IDHO2NO2,k) * xnew(PNA)  &
    + rct(41,k) * xnew(PNA) * xnew(OH)
  
  !-> NO2 cont.
    P = P +   rct(42,k) * xnew(C2O3) * xnew(NO)  &
    + rct(44,k) * xnew(PAN)  &
    + 0.6* rcphot(IDHO2NO2,k) *0.14 * xnew(PAN)  &
    + rct(48,k) * xnew(CXO3) * xnew(NO)  &
    + rct(44,k) * xnew(PANX)  &
    + 0.6* rcphot(IDHO2NO2,k) *0.14 * xnew(PANX)  &
    + rct(53,k) * xnew(MEO2) * xnew(NO)  &
    + rct(56,k) * xnew(XO2H) * xnew(NO)  &
    + rct(56,k) * xnew(XO2) * xnew(NO)  &
    + rct(67,k) * xnew(HCO3) * xnew(NO)  &
    + 0.5* rct(94,k) * xnew(ETH) * xnew(NO3)  &
    + 0.5* rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 0.5* rct(102,k) * xnew(IOLE) * xnew(NO3)  &
    + 0.9* rct(105,k) * xnew(ISO2) * xnew(NO)  &
    + 0.35* rct(109,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.142* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 0.444* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + 0.47* rct(120,k) * xnew(TERP) * xnew(NO3)  &
    + 0.918* rct(56,k) * xnew(BZO2) * xnew(NO)
  
  !-> NO2 cont.
    P = P +   0.86* rct(56,k) * xnew(TO2) * xnew(NO)  &
    + 0.86* rct(56,k) * xnew(XLO2) * xnew(NO)  &
    + 0.5* rct(133,k) * xnew(XOPN) * xnew(NO3)  &
    + rct(138,k) * xnew(OPO3) * xnew(NO)  &
    + rct(44,k) * xnew(OPAN)  &
    + 0.5* rct(117,k) * xnew(OPAN) * xnew(OH)  &
    + rct(133,k) * xnew(PANX) * xnew(OH)  &
    + rcemis(NO2,k)
  
  L = rcphot(IDNO2,k)  &
    + rct(4,k) * xnew(OP)  &
    + rct(5,k) * xnew(OP)  &
    + rct(22,k) * xnew(O3)  &
    + rct(24,k) * xnew(NO3)  &
    + rct(30,k) * xnew(NO3)  &
    + rct(34,k) * xnew(NO)  &
    + rct(37,k) * xnew(OH)  &
    + rct(39,k) * xnew(HO2)  &
    + rct(43,k) * xnew(C2O3)  &
    + rct(43,k) * xnew(CXO3)  &
    + rct(89,k) * xnew(ROR)  &
    + rct(127,k) * xnew(CRO)  &
    + rct(43,k) * xnew(OPO3)
  
  xnew(NO2) = (xold(NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> O3

  P = rct(1,k) * xnew(OP)  &
    + 0.15* rct(45,k) * xnew(C2O3) * xnew(HO2)  &
    + 0.15* rct(45,k) * xnew(CXO3) * xnew(HO2)  &
    + 0.15* rct(45,k) * xnew(OPO3) * xnew(HO2)
  
  L = rct(2,k) * xnew(NO)  &
    + rct(6,k) * xnew(OP)  &
    + rcphot(IDO3_O3P,k)  &
    + rcphot(IDO3_O1D,k)  &
    + rct(9,k) * xnew(OH)  &
    + rct(10,k) * xnew(HO2)  &
    + rct(22,k) * xnew(NO2)  &
    + rct(28,k) * xnew(NO3)  &
    + rct(93,k) * xnew(ETH)  &
    + rct(97,k) * xnew(OLE)  &
    + rct(101,k) * xnew(IOLE)  &
    + rct(108,k) * xnew(C5H8)  &
    + rct(111,k) * xnew(ISPD)  &
    + rct(119,k) * xnew(TERP)  &
    + rct(132,k) * xnew(XOPN)  &
    + rct(135,k) * xnew(OPEN)  &
    + AQRCK(ICLRC2,K) * xnew(SO2)
  
  xnew(O3) = (xold(O3) + dt2 * P) / (1.0 + dt2 * L)


!-> OLE

  P = 0.24* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 0.098* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + rcemis(OLE,k)
  
  L = rct(95,k) * xnew(OP)  &
    + rct(96,k) * xnew(OH)  &
    + rct(97,k) * xnew(O3)  &
    + rct(98,k) * xnew(NO3)
  
  xnew(OLE) = (xold(OLE) + dt2 * P) / (1.0 + dt2 * L)


!-> OPAN

  P = rct(43,k) * xnew(OPO3) * xnew(NO2)
  
  L = rct(44,k)  &
    + rct(117,k) * xnew(OH)
  
  xnew(OPAN) = (xold(OPAN) + dt2 * P) / (1.0 + dt2 * L)


!-> OPEN

  P = 0.118* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + 0.918* rct(56,k) * xnew(BZO2) * xnew(NO)  &
    + rct(46,k) * xnew(BZO2) * xnew(C2O3)  &
    + rct(52,k) * xnew(BZO2) * xnew(RO2)  &
    + 0.1* rct(123,k) * xnew(TOL) * xnew(OH)  &
    + 0.66* rct(56,k) * xnew(TO2) * xnew(NO)  &
    + 0.77* rct(46,k) * xnew(TO2) * xnew(C2O3)  &
    + 0.77* rct(52,k) * xnew(TO2) * xnew(RO2)  &
    + 0.3* rct(56,k) * xnew(XLO2) * xnew(NO)  &
    + 0.35* rct(46,k) * xnew(XLO2) * xnew(C2O3)  &
    + 0.35* rct(52,k) * xnew(XLO2) * xnew(RO2)  &
    + 0.025* rct(125,k) * xnew(CRES) * xnew(OH)  &
    + rcphot(IDHONO,k) *0.076 * xnew(CRON)  &
    + 0.25* rct(133,k) * xnew(XOPN) * xnew(NO3)
  
  L = rcphot(IDNO2,k) *0.087  &
    + rct(134,k) * xnew(OH)  &
    + rct(135,k) * xnew(O3)  &
    + rct(130,k) * xnew(NO3)
  
  xnew(OPEN) = (xold(OPEN) + dt2 * P) / (1.0 + dt2 * L)


!-> PACD

  P = 0.41* rct(45,k) * xnew(C2O3) * xnew(HO2)  &
    + 0.41* rct(45,k) * xnew(CXO3) * xnew(HO2)  &
    + 0.41* rct(45,k) * xnew(OPO3) * xnew(HO2)
  
  L = rct(58,k) * xnew(OH)
  
  xnew(PACD) = (xold(PACD) + dt2 * P) / (1.0 + dt2 * L)


!-> PAN

  P = rct(43,k) * xnew(C2O3) * xnew(NO2)
  
  L = rct(44,k)  &
    + rcphot(IDHO2NO2,k) *0.14
  
  xnew(PAN) = (xold(PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> PANX

  P = rct(43,k) * xnew(CXO3) * xnew(NO2)
  
  L = rct(44,k)  &
    + rcphot(IDHO2NO2,k) *0.14  &
    + rct(133,k) * xnew(OH)
  
  xnew(PANX) = (xold(PANX) + dt2 * P) / (1.0 + dt2 * L)


!-> PAR

  P = 0.26* rct(85,k) * xnew(PRPA) * xnew(OH)  &
    + 0.2* rct(95,k) * xnew(OLE) * xnew(OP)  &
    + 0.1* rct(99,k) * xnew(IOLE) * xnew(OP)  &
    + 0.29* rct(101,k) * xnew(IOLE) * xnew(O3)  &
    + rct(102,k) * xnew(IOLE) * xnew(NO3)  &
    + 0.25* rct(104,k) * xnew(C5H8) * xnew(OP)  &
    + 0.35* rct(108,k) * xnew(C5H8) * xnew(O3)  &
    + 0.117* rct(110,k) * xnew(ISPD) * xnew(OH)  &
    + 0.717* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + 0.24* rcphot(IDHCHO_H2,k) *0.62 * xnew(ISPD)  &
    + 2.175* rct(106,k) * xnew(EPX2) * xnew(HO2)  &
    + 2.175* rct(105,k) * xnew(EPX2) * xnew(NO)  &
    + 1.74* rct(46,k) * xnew(EPX2) * xnew(C2O3)  &
    + 2.175* rct(52,k) * xnew(EPX2) * xnew(RO2)  &
    + 2.7* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + 5.12* rct(117,k) * xnew(TERP) * xnew(OP)  &
    + 1.66* rct(118,k) * xnew(TERP) * xnew(OH)  &
    + 7* rct(119,k) * xnew(TERP) * xnew(O3)  &
    + rcemis(PAR,k)
  
  L = rct(86,k) * xnew(OH)  &
    + 1.0e6 * xnew(PARLOSS)
  
  xnew(PAR) = (xold(PAR) + dt2 * P) / (1.0 + dt2 * L)


!-> PNA

  P = rct(39,k) * xnew(HO2) * xnew(NO2)
  
  L = rct(40,k)  &
    + rcphot(IDHO2NO2,k)  &
    + rct(41,k) * xnew(OH)
  
  xnew(PNA) = (xold(PNA) + dt2 * P) / (1.0 + dt2 * L)


!-> PRPA

  P = rcemis(PRPA,k)
  
  L = rct(85,k) * xnew(OH)
  
  xnew(PRPA) = (xold(PRPA) + dt2 * P) / (1.0 + dt2 * L)


!-> ROOH

  P = rct(57,k) * xnew(XO2H) * xnew(HO2)  &
    + rct(57,k) * xnew(XO2) * xnew(HO2)  &
    + rct(57,k) * xnew(XO2N) * xnew(HO2)
  
  L = rct(58,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(ROOH) = (xold(ROOH) + dt2 * P) / (1.0 + dt2 * L)


!-> SO2

  P = rcemis(SO2,k)
  
  L = AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH)  &
    + AQRCK(ICLRC1,K) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(O3)  &
    + AQRCK(ICLRC3,K)
  
  xnew(SO2) = (xold(SO2) + dt2 * P) / (1.0 + dt2 * L)


!-> TERP

  P = rcemis(TERP,k)  &
    + RCBIO(2,K)
  
  L = rct(117,k) * xnew(OP)  &
    + rct(118,k) * xnew(OH)  &
    + rct(119,k) * xnew(O3)  &
    + rct(120,k) * xnew(NO3)
  
  xnew(TERP) = (xold(TERP) + dt2 * P) / (1.0 + dt2 * L)


!-> TOL

  P = rcemis(TOL,k)
  
  L = rct(123,k) * xnew(OH)
  
  xnew(TOL) = (xold(TOL) + dt2 * P) / (1.0 + dt2 * L)


!-> XOPN

  P = 0.2* rct(56,k) * xnew(TO2) * xnew(NO)  &
    + 0.23* rct(46,k) * xnew(TO2) * xnew(C2O3)  &
    + 0.23* rct(52,k) * xnew(TO2) * xnew(RO2)  &
    + 0.244* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + 0.56* rct(56,k) * xnew(XLO2) * xnew(NO)  &
    + 0.65* rct(46,k) * xnew(XLO2) * xnew(C2O3)  &
    + 0.65* rct(52,k) * xnew(XLO2) * xnew(RO2)
  
  L = rcphot(IDNO2,k) *0.087  &
    + rct(131,k) * xnew(OH)  &
    + rct(132,k) * xnew(O3)  &
    + rct(133,k) * xnew(NO3)
  
  xnew(XOPN) = (xold(XOPN) + dt2 * P) / (1.0 + dt2 * L)


!-> XYL

  P = rcemis(XYL,k)
  
  L = rct(124,k) * xnew(OH)
  
  xnew(XYL) = (xold(XYL) + dt2 * P) / (1.0 + dt2 * L)


!-> NTR1

  P = 0.5* rct(56,k) * xnew(XO2N) * xnew(NO)  &
    + rct(89,k) * xnew(ROR) * xnew(NO2)  &
    + 0.5* rct(94,k) * xnew(ETH) * xnew(NO3)  &
    + 0.5* rct(98,k) * xnew(OLE) * xnew(NO3)  &
    + 0.5* rct(102,k) * xnew(IOLE) * xnew(NO3)
  
  L = rct(59,k) * xnew(OH)
  
  xnew(NTR1) = (xold(NTR1) + dt2 * P) / (1.0 + dt2 * L)


!-> NTR2

  P = 0.5* rct(56,k) * xnew(XO2N) * xnew(NO)  &
    + rct(59,k) * xnew(NTR1) * xnew(OH)  &
    + 0.65* rct(109,k) * xnew(C5H8) * xnew(NO3)  &
    + 0.142* rct(112,k) * xnew(ISPD) * xnew(NO3)  &
    + 0.266* rct(116,k) * xnew(INTR) * xnew(OH)  &
    + 0.53* rct(120,k) * xnew(TERP) * xnew(NO3)  &
    + 0.082* rct(56,k) * xnew(BZO2) * xnew(NO)  &
    + 0.14* rct(56,k) * xnew(TO2) * xnew(NO)  &
    + 0.14* rct(56,k) * xnew(XLO2) * xnew(NO)  &
    + rct(129,k) * xnew(CRON) * xnew(OH)  &
    + rct(130,k) * xnew(CRON) * xnew(NO3)  &
    + 0.5* rct(133,k) * xnew(XOPN) * xnew(NO3)  &
    + 0.5* rct(117,k) * xnew(OPAN) * xnew(OH)
  
  L = rct(139,k)
  
  xnew(NTR2) = (xold(NTR2) + dt2 * P) / (1.0 + dt2 * L)


!-> ECH4

  ! P = 0.0
  
  L = rct(80,k) * xnew(OH)
  
  xnew(ECH4) = xold(ECH4) / (1.0 + dt2 * L)


!-> CL2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CL2! xnew(CL2) = max(0.0, xold(CL2))


!-> CL

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CL! xnew(CL) = max(0.0, xold(CL))


!-> CLO

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CLO! xnew(CLO) = max(0.0, xold(CLO))


!-> HCL

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HCL! xnew(HCL) = max(0.0, xold(HCL))


!-> HOCL

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HOCL! xnew(HOCL) = max(0.0, xold(HOCL))


!-> CLN2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CLN2! xnew(CLN2) = max(0.0, xold(CLN2))


!-> CLN3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CLN3! xnew(CLN3) = max(0.0, xold(CLN3))


!-> FMCL

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for FMCL! xnew(FMCL) = max(0.0, xold(FMCL))


!-> BR2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for BR2! xnew(BR2) = max(0.0, xold(BR2))


!-> BR

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for BR! xnew(BR) = max(0.0, xold(BR))


!-> BRO

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for BRO! xnew(BRO) = max(0.0, xold(BRO))


!-> HBR

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HBR! xnew(HBR) = max(0.0, xold(HBR))


!-> HOBR

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HOBR! xnew(HOBR) = max(0.0, xold(HOBR))


!-> BRN2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for BRN2! xnew(BRN2) = max(0.0, xold(BRN2))


!-> BRN3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for BRN3! xnew(BRN3) = max(0.0, xold(BRN3))


!-> FMBR

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for FMBR! xnew(FMBR) = max(0.0, xold(FMBR))


!-> I2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for I2! xnew(I2) = max(0.0, xold(I2))


!-> Iodine

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for Iodine! xnew(Iodine) = max(0.0, xold(Iodine))


!-> IO

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for IO! xnew(IO) = max(0.0, xold(IO))


!-> OIO

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OIO! xnew(OIO) = max(0.0, xold(OIO))


!-> I2O2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for I2O2! xnew(I2O2) = max(0.0, xold(I2O2))


!-> IXOY

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for IXOY! xnew(IXOY) = max(0.0, xold(IXOY))


!-> HI

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HI! xnew(HI) = max(0.0, xold(HI))


!-> HOI

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HOI! xnew(HOI) = max(0.0, xold(HOI))


!-> HIO3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HIO3! xnew(HIO3) = max(0.0, xold(HIO3))


!-> INO2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for INO2! xnew(INO2) = max(0.0, xold(INO2))


!-> INO3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for INO3! xnew(INO3) = max(0.0, xold(INO3))


!-> CH3I

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CH3I! xnew(CH3I) = max(0.0, xold(CH3I))


!-> MI2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MI2! xnew(MI2) = max(0.0, xold(MI2))


!-> MIB

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MIB! xnew(MIB) = max(0.0, xold(MIB))


!-> MIC

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MIC! xnew(MIC) = max(0.0, xold(MIC))


!-> MBC

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MBC! xnew(MBC) = max(0.0, xold(MBC))


!-> MB2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MB2! xnew(MB2) = max(0.0, xold(MB2))


!-> MBC2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MBC2! xnew(MBC2) = max(0.0, xold(MBC2))


!-> MB2C

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MB2C! xnew(MB2C) = max(0.0, xold(MB2C))


!-> MB3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MB3! xnew(MB3) = max(0.0, xold(MB3))


!-> SSCL

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for SSCL! xnew(SSCL) = max(0.0, xold(SSCL))


!-> SSBR

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for SSBR! xnew(SSBR) = max(0.0, xold(SSBR))


!-> SSN3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for SSN3! xnew(SSN3) = max(0.0, xold(SSN3))


!-> IALK

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for IALK! xnew(IALK) = max(0.0, xold(IALK))


!-> H2

  P = rcphot(IDHCHO_H2,k) * xnew(FORM)
  
  L = rct(78,k) * xnew(OH)
  
  xnew(H2) = (xold(H2) + dt2 * P) / (1.0 + dt2 * L)


!-> NH3

  P = rcemis(NH3,k)
  
  ! L = 0.0
  
  xnew(NH3) = xold(NH3) + dt2 * P


!-> SO4

  P = rcemis(SO4,k)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)  &
    + AQRCK(ICLRC1,K) * xnew(SO2) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(SO2) * xnew(O3)  &
    + AQRCK(ICLRC3,K) * xnew(SO2)
  
  ! L = 0.0
  
  xnew(SO4) = xold(SO4) + dt2 * P


!-> NH4_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NH4_f! xnew(NH4_f) = max(0.0, xold(NH4_f))


!-> NO3_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NO3_f! xnew(NO3_f) = max(0.0, xold(NO3_f))


!-> NO3_c

  P = rct(143,k) * xnew(HNO3)  &
    + rct(144,k) * xnew(HNO3)
  
  ! L = 0.0
  
  xnew(NO3_c) = xold(NO3_c) + dt2 * P


!-> PARLOSS

  P = 2.5* rcphot(IDACETON,k) *1.5 * xnew(KET)  &
    + 0.11* rct(86,k) * xnew(PAR) * xnew(OH)  &
    + 2.7* rct(87,k) * xnew(ROR)  &
    + 0.73* rct(96,k) * xnew(OLE) * xnew(OH)  &
    + 0.79* rct(97,k) * xnew(OLE) * xnew(O3)  &
    + rct(98,k) * xnew(OLE) * xnew(NO3)
  
  L = 1.0e6 * xnew(PAR)
  
  xnew(PARLOSS) = (xold(PARLOSS) + dt2 * P) / (1.0 + dt2 * L)


!-> TERPPeroxy

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for TERPPeroxy! xnew(TERPPeroxy) = max(0.0, xold(TERPPeroxy))


!-> VBS_TEST

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for VBS_TEST! xnew(VBS_TEST) = max(0.0, xold(VBS_TEST))


!-> shipNOx

  P = rcemis(SHIPNOX,k)
  
  L = rct(146,k) * xnew(OH)  &
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

  P = 0.00383* RCBIO(2,K)
  
  ! L = 0.0
  
  xnew(SQT_SOA_NV) = xold(SQT_SOA_NV) + dt2 * P


