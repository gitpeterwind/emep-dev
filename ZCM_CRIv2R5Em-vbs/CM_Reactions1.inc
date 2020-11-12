!-> OD

  P = rcphot(IDO3_O1D,k) * xnew(O3)
  
  L = rct(7,k)  &
    + rct(8,k)  &
    + rct(16,k)
  
  xnew(OD) = (xold(OD) + dt2 * P) / (1.0 + dt2 * L)


!-> OP

  P = rct(7,k) * xnew(OD)  &
    + rct(8,k) * xnew(OD)  &
    + rcphot(IDO3_O3P,k) * xnew(O3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + rcphot(IDNO3_NO2,k) * xnew(NO3)
  
  L = rct(1,k)  &
    + rct(2,k)  &
    + rct(3,k) * xnew(O3)  &
    + rct(4,k) * xnew(NO)  &
    + rct(5,k) * xnew(NO2)  &
    + rct(6,k) * xnew(NO2)
  
  xnew(OP) = (xold(OP) + dt2 * P) / (1.0 + dt2 * L)


!-> OH

  P = 2* rct(16,k) * xnew(OD)  &
    + rct(21,k) * xnew(HO2) * xnew(O3)  &
    + rct(27,k) * xnew(HO2) * xnew(NO)  &
    + 4.0e-12 * xnew(HO2) * xnew(NO3)  &
    + 0.13* rct(49,k) * xnew(O3) * xnew(C2H4)  &
    + 0.36* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 0.57* rct(54,k) * xnew(O3) * xnew(TBUT2ENE)  &
    + 0.636* rct(55,k) * xnew(OH) * xnew(C2H2)  &
    + 0.25* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.8* rct(65,k) * xnew(APINENE) * xnew(O3)  &
    + 0.35* rct(67,k) * xnew(BPINENE) * xnew(O3)  &
    + 0.239* rct(70,k) * xnew(O3) * xnew(MACR)  &
    + 0.36* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + 0.75* 2.4e-17 * xnew(O3) * xnew(UCARB12)  &
    + 0.508* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 5.2e-11 * xnew(OH) * xnew(HPUCARB12)  &
    + rct(77,k) * xnew(NO3) * xnew(HPUCARB12)  &
    + 2* 2.4e-17 * xnew(O3) * xnew(HPUCARB12)  &
    + 3.64e-11 * xnew(OH) * xnew(DHPCARB9)  &
    + 0.4* rct(79,k) * xnew(OH) * xnew(CH3OOH)
  
  !-> OH cont.
    P = P +   1.36e-11 * xnew(OH) * xnew(C2H5OOH)  &
    + 2.78e-11 * xnew(OH) * xnew(IC3H7OOH)  &
    + 2.13e-11 * xnew(OH) * xnew(HOC2H4OOH)  &
    + 1.2e-11 * xnew(OH) * xnew(RN8OOH)  &
    + 2.5e-11 * xnew(OH) * xnew(RN9OOH)  &
    + 1.89e-11 * xnew(OH) * xnew(RN10OOH)  &
    + 2.5e-11 * xnew(OH) * xnew(RN11OOH)  &
    + 3.25e-11 * xnew(OH) * xnew(RN12OOH)  &
    + 3.57e-11 * xnew(OH) * xnew(RN13OOH)  &
    + 3.2e-11 * xnew(OH) * xnew(RN14OOH)  &
    + 3.74e-11 * xnew(OH) * xnew(RN15OOH)  &
    + 4.21e-11 * xnew(OH) * xnew(RN16OOH)  &
    + 1.27e-10 * xnew(OH) * xnew(RN19OOH)  &
    + 3.84e-11 * xnew(OH) * xnew(RU10OOH)  &
    + 0.94* 7.14e-11 * xnew(OH) * xnew(RU14OOH)  &
    + 9.77e-11 * xnew(OH) * xnew(RA13OOH)  &
    + 9.64e-11 * xnew(OH) * xnew(RA16OOH)  &
    + 2.38e-11 * xnew(OH) * xnew(RTN28OOH)  &
    + 7.5e-11 * xnew(OH) * xnew(RTX22OOH)  &
    + 8.59e-11 * xnew(OH) * xnew(RTX24OOH)
  
  !-> OH cont.
    P = P +   5.64e-11 * xnew(OH) * xnew(DHPR12OOH)  &
    + 5.22e-12 * xnew(OH) * xnew(NRN6OOH)  &
    + 6.5e-12 * xnew(OH) * xnew(NRN9OOH)  &
    + 7.15e-12 * xnew(OH) * xnew(NRN12OOH)  &
    + 1.03e-10 * xnew(OH) * xnew(NRU14OOH)  &
    + 9.5e-12 * xnew(OH) * xnew(NRTN28OOH)  &
    + 0.7* 4.33e-12 * xnew(OH) * xnew(HMML)  &
    + rct(92,k) * xnew(BRU14O2)  &
    + rct(95,k) * xnew(HCOCO)  &
    + rct(97,k) * xnew(DHPR12O2)  &
    + rct(98,k) * xnew(RU12O2)  &
    + rct(99,k) * xnew(MACRO2)  &
    + rct(109,k) * xnew(DHPR12O2) * xnew(NO)  &
    + rct(113,k) * xnew(DHPR12O2) * xnew(NO3)  &
    + 0.454* rct(129,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(C2H5CO3) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(MACO3) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(HOCH2CO3) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(HCOCO3) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(CO2C3CO3) * xnew(HO2)
  
  !-> OH cont.
    P = P +   0.44* rct(129,k) * xnew(RTN26O2) * xnew(HO2)  &
    + 7.6e-13 * xnew(RO2POOL) * xnew(DHPR12O2)  &
    + 2* rcphot(IDH2O2,k) * xnew(H2O2)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + rcphot(IDHCHO_H2,k) *6 * xnew(HUCARB9)  &
    + 1.5* rcphot(IDHCHO_H2,k) *12 * xnew(HPUCARB12)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(DHPCARB9)  &
    + 2* rcphot(IDCH3O2H,k) *2 * xnew(DHPCARB9)  &
    + 0.5* rcphot(IDMEK,k) *132.76 * xnew(NUCARB12)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(IC3H7OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(HOC2H4OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN8OOH)  &
    + rcphot(IDMEK,k) * xnew(RN8OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN9OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN10OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN11OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN12OOH)
  
  !-> OH cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(RN13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN14OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN15OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN16OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN17OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN19OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RU10OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RU12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RU13OOH)  &
    + rcphot(IDMEK,k) * xnew(RU13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA16OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA19OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN10OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN14OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN23OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN24OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN25OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN26OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN28OOH)
  
  !-> OH cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(RTX22OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX24OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX28OOH)  &
    + rcphot(IDHCHO_H,k) *1.67 * xnew(DHPR12OOH)  &
    + 2* rcphot(IDCH3O2H,k) *3 * xnew(DHPR12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRN6OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRN9OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRN12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRU12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRU14OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRTN28OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRTX28OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(HOCH2CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(HCOCO3H)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HCOCO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(CO2C3CO3H)  &
    + rcphot(IDMEK,k) * xnew(CO2C3CO3H)
  
  L = rct(17,k) * xnew(O3)  &
    + rct(18,k) * xnew(H2)  &
    + rct(19,k) * xnew(CO)  &
    + rct(20,k) * xnew(H2O2)  &
    + rct(22,k) * xnew(HO2)  &
    + rct(25,k) * xnew(NO)  &
    + rct(26,k) * xnew(NO2)  &
    + 2.0e-11 * xnew(NO3)  &
    + rct(30,k) * xnew(HO2NO2)  &
    + rct(31,k) * xnew(HONO)  &
    + rct(32,k) * xnew(HNO3)  &
    + rct(33,k) * xnew(CH3OH)  &
    + rct(34,k) * xnew(C2H5OH)  &
    + rct(35,k) * xnew(NPROPOL)  &
    + rct(36,k) * xnew(IPROPOL)  &
    + rct(37,k) * xnew(HCHO)  &
    + rct(39,k) * xnew(CH3CHO)  &
    + rct(41,k) * xnew(C2H5CHO)  &
    + rct(43,k) * xnew(CH4)  &
    + rct(44,k) * xnew(C2H6)
  
  !-> OH cont.
    L = L +   rct(45,k) * xnew(C3H8)  &
    + rct(46,k) * xnew(NC4H10)  &
    + rct(47,k) * xnew(C2H4)  &
    + rct(50,k) * xnew(C3H6)  &
    + rct(53,k) * xnew(TBUT2ENE)  &
    + rct(55,k) * xnew(C2H2)  &
    + rct(56,k) * xnew(BENZENE)  &
    + rct(57,k) * xnew(TOLUENE)  &
    + 1.36e-11 * xnew(OXYL)  &
    + 2.63e-11 * xnew(AROH14)  &
    + 9.0e-13 * xnew(ARNOH14)  &
    + 4.65e-11 * xnew(AROH17)  &
    + 1.53e-12 * xnew(ARNOH17)  &
    + rct(58,k) * xnew(C5H8)  &
    + rct(61,k) * xnew(CH3COCH3)  &
    + rct(62,k) * xnew(MEK)  &
    + rct(63,k) * xnew(APINENE)  &
    + rct(66,k) * xnew(BPINENE)  &
    + 4.5e-13 * xnew(HCOOH)  &
    + rct(68,k) * xnew(CH3CO2H)
  
  !-> OH cont.
    L = L +   1.23e-11 * xnew(HCOCO2H)  &
    + 6.65e-12 * xnew(RCOOH25)  &
    + rct(69,k) * xnew(MACR)  &
    + rct(71,k) * xnew(MVK)  &
    + rct(73,k) * xnew(GLYOX)  &
    + rct(74,k) * xnew(MGLYOX)  &
    + rct(75,k) * xnew(ACETOL)  &
    + 8.0e-12 * xnew(HOCH2CHO)  &
    + 7.15e-11 * xnew(CO2C3CHO)  &
    + 2.4e-13 * xnew(CARB9)  &
    + 5.86e-12 * xnew(CARB10)  &
    + 1.38e-12 * xnew(CARB12)  &
    + 1.65e-11 * xnew(CARB13)  &
    + 1.87e-11 * xnew(CARB14)  &
    + 1.25e-11 * xnew(CARB16)  &
    + 1.19e-10 * xnew(HCC7CO)  &
    + 6.42e-11 * xnew(UCARB12)  &
    + 4.19e-11 * xnew(UCARB13)  &
    + 5.78e-11 * xnew(HUCARB9)  &
    + 5.2e-11 * xnew(HPUCARB12)
  
  !-> OH cont.
    L = L +   3.42e-11 * xnew(DHCARB9)  &
    + 3.64e-11 * xnew(DHPCARB9)  &
    + 4.79e-12 * xnew(CCARB12)  &
    + 5.2e-11 * xnew(UDCARB8)  &
    + 5.58e-11 * xnew(UDCARB11)  &
    + 7.0e-11 * xnew(UDCARB14)  &
    + 1.0e-10 * xnew(TNCARB10)  &
    + 1.0e-12 * xnew(TNCARB15)  &
    + 4.2e-11 * xnew(TNCARB26)  &
    + 4.55e-12 * xnew(TXCARB22)  &
    + 1.55e-11 * xnew(TXCARB24)  &
    + 4.16e-11 * xnew(NUCARB12)  &
    + rct(79,k) * xnew(CH3OOH)  &
    + 1.36e-11 * xnew(C2H5OOH)  &
    + 2.78e-11 * xnew(IC3H7OOH)  &
    + 2.13e-11 * xnew(HOC2H4OOH)  &
    + 1.2e-11 * xnew(RN8OOH)  &
    + 2.5e-11 * xnew(RN9OOH)  &
    + 1.89e-11 * xnew(RN10OOH)  &
    + 2.5e-11 * xnew(RN11OOH)
  
  !-> OH cont.
    L = L +   3.25e-11 * xnew(RN12OOH)  &
    + 3.57e-11 * xnew(RN13OOH)  &
    + 3.2e-11 * xnew(RN14OOH)  &
    + 3.74e-11 * xnew(RN15OOH)  &
    + 4.21e-11 * xnew(RN16OOH)  &
    + 7.06e-11 * xnew(RN17OOH)  &
    + 1.27e-10 * xnew(RN19OOH)  &
    + 3.84e-11 * xnew(RU10OOH)  &
    + 3.5e-11 * xnew(RU12OOH)  &
    + 6.71e-11 * xnew(RU13OOH)  &
    + 7.14e-11 * xnew(RU14OOH)  &
    + 9.77e-11 * xnew(RA13OOH)  &
    + 9.64e-11 * xnew(RA16OOH)  &
    + 7.62e-11 * xnew(RA19OOH)  &
    + 4.23e-12 * xnew(RTN10OOH)  &
    + 8.69e-11 * xnew(RTN14OOH)  &
    + 2.05e-11 * xnew(RTN23OOH)  &
    + 1.05e-11 * xnew(RTN24OOH)  &
    + 1.66e-11 * xnew(RTN25OOH)  &
    + 1.2e-11 * xnew(RTN26OOH)
  
  !-> OH cont.
    L = L +   2.38e-11 * xnew(RTN28OOH)  &
    + 7.5e-11 * xnew(RTX22OOH)  &
    + 8.59e-11 * xnew(RTX24OOH)  &
    + 2.0e-11 * xnew(RTX28OOH)  &
    + 5.64e-11 * xnew(DHPR12OOH)  &
    + 5.22e-12 * xnew(NRN6OOH)  &
    + 6.5e-12 * xnew(NRN9OOH)  &
    + 7.15e-12 * xnew(NRN12OOH)  &
    + 2.65e-11 * xnew(NRU12OOH)  &
    + 1.03e-10 * xnew(NRU14OOH)  &
    + 9.5e-12 * xnew(NRTN28OOH)  &
    + 9.58e-12 * xnew(NRTX28OOH)  &
    + 3.7e-12 * xnew(CH3CO3H)  &
    + 4.42e-12 * xnew(C2H5CO3H)  &
    + 6.19e-12 * xnew(HOCH2CO3H)  &
    + 1.58e-11 * xnew(HCOCO3H)  &
    + 4.18e-12 * xnew(CO2C3CO3H)  &
    + 1.16e-11 * xnew(IEPOX)  &
    + 1.5e-12 * xnew(ANHY)  &
    + 4.33e-12 * xnew(HMML)
  
  !-> OH cont.
    L = L +   rct(80,k) * xnew(CH3NO3)  &
    + rct(81,k) * xnew(C2H5NO3)  &
    + rct(82,k) * xnew(IC3H7NO3)  &
    + 6.7e-13 * xnew(NOA)  &
    + 1.09e-12 * xnew(HOC2H4NO3)  &
    + 1.31e-12 * xnew(RN9NO3)  &
    + 5.8e-13 * xnew(RN10NO3)  &
    + 1.79e-12 * xnew(RN12NO3)  &
    + 9.2e-13 * xnew(RN13NO3)  &
    + 1.03e-11 * xnew(RN15NO3)  &
    + 1.85e-12 * xnew(RN16NO3)  &
    + 9.6e-11 * xnew(RN19NO3)  &
    + 5.26e-13 * xnew(RU10NO3)  &
    + 2.5e-12 * xnew(RU12NO3)  &
    + 1.6e-12 * xnew(RU13NO3)  &
    + 3.0e-11 * xnew(RU14NO3)  &
    + 7.3e-11 * xnew(RA13NO3)  &
    + 7.16e-11 * xnew(RA16NO3)  &
    + 8.31e-11 * xnew(RA19NO3)  &
    + 2.88e-12 * xnew(RTN25NO3)
  
  !-> OH cont.
    L = L +   4.35e-12 * xnew(RTN28NO3)  &
    + 4.74e-12 * xnew(RTX22NO3)  &
    + 6.48e-12 * xnew(RTX24NO3)  &
    + 3.53e-12 * xnew(RTX28NO3)  &
    + 3.0e-14 * xnew(PAN)  &
    + 1.12e-12 * xnew(PHAN)  &
    + 1.27e-12 * xnew(PPN)  &
    + 5.93e-13 * xnew(CO2C3PAN)  &
    + 2.9e-11 * xnew(MPAN)  &
    + 3.66e-12 * xnew(RTN26PAN)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(SO2)  &
    + rct(26,k) * xnew(shipNOx)
  
  xnew(OH) = (xold(OH) + dt2 * P) / (1.0 + dt2 * L)


!-> CRU14

  P = 0.538* rct(58,k) * xnew(OH) * xnew(C5H8)  &
    + rct(84,k) * xnew(CRU14O2)  &
    + rct(86,k) * xnew(BRU14O2)
  
  L = rct(83,k)  &
    + rct(85,k)
  
  xnew(CRU14) = (xold(CRU14) + dt2 * P) / (1.0 + dt2 * L)


!-> TRU14

  P = 0.4* rct(58,k) * xnew(OH) * xnew(C5H8)  &
    + rct(88,k) * xnew(TRU14O2)  &
    + rct(90,k) * xnew(BRU14O2)
  
  L = rct(87,k)  &
    + rct(89,k)
  
  xnew(TRU14) = (xold(TRU14) + dt2 * P) / (1.0 + dt2 * L)


!-> HCOCO

  P = rct(73,k) * xnew(OH) * xnew(GLYOX)  &
    + 4.0e-16 * xnew(NO3) * xnew(GLYOX)
  
  L = rct(93,k)  &
    + rct(94,k)  &
    + rct(95,k)  &
    + rct(96,k)
  
  xnew(HCOCO) = (xold(HCOCO) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2

  P = rct(17,k) * xnew(OH) * xnew(O3)  &
    + rct(18,k) * xnew(OH) * xnew(H2)  &
    + rct(19,k) * xnew(OH) * xnew(CO)  &
    + rct(20,k) * xnew(OH) * xnew(H2O2)  &
    + 2.0e-11 * xnew(OH) * xnew(NO3)  &
    + rct(29,k) * xnew(HO2NO2)  &
    + rct(33,k) * xnew(OH) * xnew(CH3OH)  &
    + 0.95* rct(34,k) * xnew(OH) * xnew(C2H5OH)  &
    + 0.49* rct(35,k) * xnew(OH) * xnew(NPROPOL)  &
    + 0.861* rct(36,k) * xnew(OH) * xnew(IPROPOL)  &
    + rct(37,k) * xnew(OH) * xnew(HCHO)  &
    + rct(38,k) * xnew(NO3) * xnew(HCHO)  &
    + 0.13* rct(49,k) * xnew(O3) * xnew(C2H4)  &
    + 0.28* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 0.125* rct(54,k) * xnew(O3) * xnew(TBUT2ENE)  &
    + 0.364* rct(55,k) * xnew(OH) * xnew(C2H2)  &
    + 0.53* rct(56,k) * xnew(OH) * xnew(BENZENE)  &
    + 0.18* rct(57,k) * xnew(OH) * xnew(TOLUENE)  &
    + 0.062* rct(58,k) * xnew(OH) * xnew(C5H8)  &
    + 0.173* rct(60,k) * xnew(O3) * xnew(C5H8)
  
  !-> HO2 cont.
    P = P +   4.5e-13 * xnew(HCOOH) * xnew(OH)  &
    + 1.23e-11 * xnew(OH) * xnew(HCOCO2H)  &
    + 0.141* rct(70,k) * xnew(O3) * xnew(MACR)  &
    + 0.18* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + rct(75,k) * xnew(ACETOL) * xnew(OH)  &
    + 0.2* 8.0e-12 * xnew(OH) * xnew(HOCH2CHO)  &
    + 5.86e-12 * xnew(OH) * xnew(CARB10)  &
    + 0.5* 2.4e-17 * xnew(O3) * xnew(UCARB12)  &
    + 0.389* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 5.78e-11 * xnew(OH) * xnew(HUCARB9)  &
    + 3.42e-11 * xnew(OH) * xnew(DHCARB9)  &
    + 0.5* 5.2e-11 * xnew(OH) * xnew(UDCARB8)  &
    + 3.5e-11 * xnew(OH) * xnew(RU12OOH)  &
    + 6.71e-11 * xnew(OH) * xnew(RU13OOH)  &
    + 2.65e-11 * xnew(OH) * xnew(NRU12OOH)  &
    + 1.6e-12 * xnew(OH) * xnew(RU13NO3)  &
    + 0.66* 3.0e-11 * xnew(OH) * xnew(RU14NO3)  &
    + 0.5* rct(91,k) * xnew(CRU14O2)  &
    + rct(93,k) * xnew(HCOCO)  &
    + rct(94,k) * xnew(HCOCO)
  
  !-> HO2 cont.
    P = P +   rct(100,k) * xnew(HCOCO3) * xnew(NO2)  &
    + 0.999* rct(105,k) * xnew(CH3O2) * xnew(NO)  &
    + 0.967* rct(106,k) * xnew(C2H5O2) * xnew(NO)  &
    + 0.98* rct(107,k) * xnew(RN10O2) * xnew(NO)  &
    + 0.965* rct(108,k) * xnew(IC3H7O2) * xnew(NO)  &
    + 0.995* rct(109,k) * xnew(HOCH2CH2O2) * xnew(NO)  &
    + 0.987* rct(109,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.979* rct(109,k) * xnew(RN9O2) * xnew(NO)  &
    + 0.959* rct(109,k) * xnew(RN12O2) * xnew(NO)  &
    + 0.552* rct(109,k) * xnew(RN13O2) * xnew(NO)  &
    + 0.936* rct(109,k) * xnew(RN15O2) * xnew(NO)  &
    + 0.722* rct(109,k) * xnew(RN19O2) * xnew(NO)  &
    + 0.975* rct(109,k) * xnew(RN15AO2) * xnew(NO)  &
    + 0.295* rct(109,k) * xnew(RU10O2) * xnew(NO)  &
    + 0.795* rct(109,k) * xnew(RU12O2) * xnew(NO)  &
    + 0.949* rct(109,k) * xnew(RU13O2) * xnew(NO)  &
    + 0.908* rct(109,k) * xnew(CRU14O2) * xnew(NO)  &
    + 0.908* rct(109,k) * xnew(TRU14O2) * xnew(NO)  &
    + 0.889* rct(109,k) * xnew(BRU14O2) * xnew(NO)  &
    + 0.918* rct(109,k) * xnew(RA13O2) * xnew(NO)
  
  !-> HO2 cont.
    P = P +   0.889* rct(109,k) * xnew(RA16O2) * xnew(NO)  &
    + 0.862* rct(109,k) * xnew(RA19AO2) * xnew(NO)  &
    + 0.862* rct(109,k) * xnew(RA19CO2) * xnew(NO)  &
    + rct(109,k) * xnew(RTN14O2) * xnew(NO)  &
    + 0.702* rct(109,k) * xnew(RTN28O2) * xnew(NO)  &
    + 0.506* rct(109,k) * xnew(RTX24O2) * xnew(NO)  &
    + 0.702* rct(109,k) * xnew(RTX28O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRU12O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRU14O2) * xnew(NO)  &
    + rct(110,k) * xnew(HOCH2CO3) * xnew(NO)  &
    + rct(110,k) * xnew(HCOCO3) * xnew(NO)  &
    + rct(112,k) * xnew(CH3O2) * xnew(NO3)  &
    + rct(113,k) * xnew(C2H5O2) * xnew(NO3)  &
    + rct(113,k) * xnew(IC3H7O2) * xnew(NO3)  &
    + rct(113,k) * xnew(HOCH2CH2O2) * xnew(NO3)  &
    + rct(113,k) * xnew(MACRO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN9O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN10O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN12O2) * xnew(NO3)  &
    + 0.602* rct(113,k) * xnew(RN13O2) * xnew(NO3)
  
  !-> HO2 cont.
    P = P +   rct(113,k) * xnew(RN15O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN19O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN15AO2) * xnew(NO3)  &
    + 0.3* rct(113,k) * xnew(RU10O2) * xnew(NO3)  &
    + 0.83* rct(113,k) * xnew(RU12O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RU13O2) * xnew(NO3)  &
    + rct(113,k) * xnew(CRU14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(TRU14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(BRU14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA13O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA16O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA19AO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA19CO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN28O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTX24O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTX28O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRU12O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRU14O2) * xnew(NO3)  &
    + rct(114,k) * xnew(HOCH2CO3) * xnew(NO3)
  
  !-> HO2 cont.
    P = P +   rct(115,k) * xnew(HCOCO3) * xnew(NO3)  &
    + 0.44* rct(129,k) * xnew(HOCH2CO3) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(HCOCO3) * xnew(HO2)  &
    + rct(130,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.6* 3.1e-13 * xnew(RO2POOL) * xnew(C2H5O2)  &
    + 0.6* 4.0e-14 * xnew(RO2POOL) * xnew(IC3H7O2)  &
    + 2.0e-12 * xnew(RO2POOL) * xnew(HOCH2CH2O2)  &
    + 3.6e-13 * xnew(RO2POOL) * xnew(MACRO2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN9O2)  &
    + 0.6* 6.0e-13 * xnew(RO2POOL) * xnew(RN10O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN12O2)  &
    + 0.602* 2.5e-13 * xnew(RO2POOL) * xnew(RN13O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN15O2)  &
    + 0.75* 2.5e-13 * xnew(RO2POOL) * xnew(RN19O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN15AO2)  &
    + 0.3* 1.22e-12 * xnew(RO2POOL) * xnew(RU10O2)  &
    + 0.83* 4.2e-13 * xnew(RO2POOL) * xnew(RU12O2)  &
    + 6.0e-13 * xnew(RO2POOL) * xnew(RU13O2)  &
    + 2.4e-12 * xnew(RO2POOL) * xnew(CRU14O2)  &
    + 2.4e-12 * xnew(RO2POOL) * xnew(TRU14O2)
  
  !-> HO2 cont.
    P = P +   1.26e-12 * xnew(RO2POOL) * xnew(BRU14O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA13O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA16O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA19AO2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA19CO2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RTN14O2)  &
    + 2.85e-13 * xnew(RO2POOL) * xnew(RTN28O2)  &
    + 2.5e-13 * xnew(RO2POOL) * xnew(RTX24O2)  &
    + 2.0e-12 * xnew(RO2POOL) * xnew(RTX28O2)  &
    + 9.6e-13 * xnew(RO2POOL) * xnew(NRU12O2)  &
    + 1.3e-12 * xnew(RO2POOL) * xnew(NRU14O2)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(HOCH2CO3)  &
    + 0.7* 1.0e-11 * xnew(RO2POOL) * xnew(HCOCO3)  &
    + 2* rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDMEK,k) *5.2 * xnew(C2H5CHO)  &
    + 2* rcphot(IDRCOCHO,k) * xnew(HCOCO2H)  &
    + rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + 2* rcphot(IDCH3O2H,k) *4.79 * xnew(GLYOX)
  
  !-> HO2 cont.
    P = P +   rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + 2* rcphot(IDHCHO_H,k) *0.6 * xnew(HOCH2CHO)  &
    + rcphot(IDMEK,k) * xnew(CARB10)  &
    + rcphot(IDMEK,k) *3.00 * xnew(CARB13)  &
    + rcphot(IDMEK,k) *3.35 * xnew(CARB16)  &
    + 1.25* rcphot(IDHCHO_H2,k) *0.0474 * xnew(UCARB12)  &
    + rcphot(IDHCHO_H2,k) *6 * xnew(HUCARB9)  &
    + 0.425* rcphot(IDHCHO_H2,k) *12 * xnew(HPUCARB12)  &
    + 2* rcphot(IDHCHO_H,k) *1.67 * xnew(DHCARB9)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(DHPCARB9)  &
    + (rcphot(IDNO2,k) *(1.28e-02)) * xnew(UDCARB8)  &
    + (2* rcphot(IDNO2,k) *(7.20e-03)) * xnew(UDCARB8)  &
    + (rcphot(IDNO2,k) *(9.00e-03)) * xnew(UDCARB11)  &
    + (rcphot(IDNO2,k) *(1.10e-02)) * xnew(UDCARB11)  &
    + (rcphot(IDNO2,k) *(9.00e-03)) * xnew(UDCARB14)  &
    + (rcphot(IDNO2,k) *(1.10e-02)) * xnew(UDCARB14)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(TNCARB26)  &
    + 0.5* rcphot(IDMEK,k) *132.76 * xnew(NUCARB12)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)
  
  !-> HO2 cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(C2H5OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(IC3H7OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(HOC2H4OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN9OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN10OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN12OOH)  &
    + 0.602* rcphot(IDCH3O2H,k) * xnew(RN13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN15OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN19OOH)  &
    + 0.383* rcphot(IDCH3O2H,k) * xnew(RU10OOH)  &
    + 0.83* rcphot(IDCH3O2H,k) * xnew(RU12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RU13OOH)  &
    + rcphot(IDMEK,k) * xnew(RU13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA16OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA19OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN14OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN28OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX24OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX28OOH)
  
  !-> HO2 cont.
    P = P +   rcphot(IDHCHO_H,k) *1.67 * xnew(DHPR12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRU12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRU14OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(HOCH2CO3H)  &
    + rcphot(IDCH3O2H,k) * xnew(HCOCO3H)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HCOCO3H)  &
    + rcphot(IDMEK,k) *0.293 * xnew(CH3NO3)  &
    + rcphot(IDHNO3,k) *1.932 * xnew(C2H5NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(IC3H7NO3)  &
    + rcphot(IDHNO3,k) *2.62 * xnew(RN10NO3)  &
    + 0.602* rcphot(IDMEK,k) *0.778 * xnew(RN13NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RN19NO3)  &
    + rcphot(IDMEK,k) *83. * xnew(RU12NO3)  &
    + rcphot(IDMEK,k) * xnew(RU13NO3)  &
    + rcphot(IDHNO3,k) *2.62 * xnew(RU14NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RA13NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RA16NO3)  &
    + rcphot(IDHCHO_H,k) *0.234 * xnew(RA19NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RTX24NO3)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)
  
  L = rct(21,k) * xnew(O3)  &
    + rct(22,k) * xnew(OH)  &
    + rct(23,k) * xnew(HO2)  &
    + rct(23,k) * xnew(HO2)  &
    + rct(24,k) * xnew(HO2)  &
    + rct(24,k) * xnew(HO2)  &
    + rct(27,k) * xnew(NO)  &
    + rct(28,k) * xnew(NO2)  &
    + 4.0e-12 * xnew(NO3)  &
    + rct(116,k) * xnew(CH3O2)  &
    + rct(117,k) * xnew(C2H5O2)  &
    + rct(118,k) * xnew(IC3H7O2)  &
    + rct(119,k) * xnew(HOCH2CH2O2)  &
    + rct(120,k) * xnew(MACRO2)  &
    + rct(121,k) * xnew(RN8O2)  &
    + rct(121,k) * xnew(RN9O2)  &
    + rct(121,k) * xnew(RN10O2)  &
    + rct(120,k) * xnew(RN11O2)  &
    + rct(120,k) * xnew(RN12O2)  &
    + rct(120,k) * xnew(RN13O2)
  
  !-> HO2 cont.
    L = L +   rct(122,k) * xnew(RN14O2)  &
    + rct(122,k) * xnew(RN15O2)  &
    + rct(122,k) * xnew(RN16O2)  &
    + rct(123,k) * xnew(RN17O2)  &
    + rct(123,k) * xnew(RN19O2)  &
    + rct(120,k) * xnew(RN13AO2)  &
    + rct(122,k) * xnew(RN15AO2)  &
    + rct(122,k) * xnew(RN16AO2)  &
    + rct(120,k) * xnew(RU10O2)  &
    + rct(122,k) * xnew(RU12O2)  &
    + rct(122,k) * xnew(RU13O2)  &
    + rct(122,k) * xnew(CRU14O2)  &
    + rct(122,k) * xnew(TRU14O2)  &
    + rct(122,k) * xnew(BRU14O2)  &
    + rct(123,k) * xnew(RA13O2)  &
    + rct(124,k) * xnew(RA16O2)  &
    + rct(125,k) * xnew(RA19AO2)  &
    + rct(125,k) * xnew(RA19CO2)  &
    + rct(122,k) * xnew(RTN10O2)  &
    + rct(123,k) * xnew(RTN14O2)
  
  !-> HO2 cont.
    L = L +   rct(126,k) * xnew(RTN23O2)  &
    + rct(126,k) * xnew(RTN24O2)  &
    + rct(126,k) * xnew(RTN25O2)  &
    + rct(127,k) * xnew(RTN28O2)  &
    + rct(127,k) * xnew(RTX28O2)  &
    + rct(126,k) * xnew(RTX24O2)  &
    + rct(126,k) * xnew(RTX22O2)  &
    + rct(122,k) * xnew(DHPR12O2)  &
    + rct(128,k) * xnew(NRN6O2)  &
    + rct(121,k) * xnew(NRN9O2)  &
    + rct(120,k) * xnew(NRN12O2)  &
    + rct(122,k) * xnew(NRU12O2)  &
    + rct(122,k) * xnew(NRU14O2)  &
    + rct(127,k) * xnew(NRTN28O2)  &
    + rct(127,k) * xnew(NRTX28O2)  &
    + rct(129,k) * xnew(CH3CO3)  &
    + rct(129,k) * xnew(C2H5CO3)  &
    + rct(129,k) * xnew(MACO3)  &
    + rct(129,k) * xnew(HOCH2CO3)  &
    + rct(129,k) * xnew(HCOCO3)
  
  !-> HO2 cont.
    L = L +   rct(129,k) * xnew(CO2C3CO3)  &
    + rct(129,k) * xnew(RTN26O2)  &
    + rct(159,k)
  
  xnew(HO2) = (xold(HO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3O2

  P = rct(43,k) * xnew(OH) * xnew(CH4)  &
    + 0.28* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 0.695* rct(54,k) * xnew(O3) * xnew(TBUT2ENE)  &
    + 0.019* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + rct(68,k) * xnew(CH3CO2H) * xnew(OH)  &
    + 0.45* 5.58e-11 * xnew(OH) * xnew(UDCARB11)  &
    + 0.6* rct(79,k) * xnew(OH) * xnew(CH3OOH)  &
    + rct(110,k) * xnew(CH3CO3) * xnew(NO)  &
    + 0.65* rct(111,k) * xnew(MACO3) * xnew(NO)  &
    + rct(114,k) * xnew(CH3CO3) * xnew(NO3)  &
    + 0.65* rct(115,k) * xnew(MACO3) * xnew(NO3)  &
    + 0.454* rct(129,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.286* rct(129,k) * xnew(MACO3) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(CH3CO3)  &
    + 0.65* 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDACETON,k) * xnew(CH3COCH3)  &
    + 0.325* rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + (rcphot(IDNO2,k) *(9.00e-03)) * xnew(UDCARB11)  &
    + 0.097* rcphot(IDCH3O2H,k) * xnew(RU10OOH)
  
  !-> CH3O2 cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(CH3CO3H)
  
  L = rct(105,k) * xnew(NO)  &
    + rct(112,k) * xnew(NO3)  &
    + rct(116,k) * xnew(HO2)  &
    + rct(130,k) * xnew(RO2POOL)  &
    + rct(131,k) * xnew(RO2POOL)
  
  xnew(CH3O2) = (xold(CH3O2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5O2

  P = rct(44,k) * xnew(OH) * xnew(C2H6)  &
    + 0.5* 5.2e-11 * xnew(OH) * xnew(UDCARB8)  &
    + 0.45* 7.0e-11 * xnew(OH) * xnew(UDCARB14)  &
    + 0.365* rct(109,k) * xnew(RN13O2) * xnew(NO)  &
    + rct(110,k) * xnew(C2H5CO3) * xnew(NO)  &
    + 0.398* rct(113,k) * xnew(RN13O2) * xnew(NO3)  &
    + rct(114,k) * xnew(C2H5CO3) * xnew(NO3)  &
    + 0.44* rct(129,k) * xnew(C2H5CO3) * xnew(HO2)  &
    + 0.398* 2.5e-13 * xnew(RO2POOL) * xnew(RN13O2)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(C2H5CO3)  &
    + rcphot(IDMEK,k) *5.2 * xnew(C2H5CHO)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + (rcphot(IDNO2,k) *(1.28e-02)) * xnew(UDCARB8)  &
    + (rcphot(IDNO2,k) *(9.00e-03)) * xnew(UDCARB14)  &
    + 0.398* rcphot(IDCH3O2H,k) * xnew(RN13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5CO3H)  &
    + 0.398* rcphot(IDMEK,k) *0.778 * xnew(RN13NO3)
  
  L = rct(106,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(117,k) * xnew(HO2)  &
    + 3.1e-13 * xnew(RO2POOL)
  
  xnew(C2H5O2) = (xold(C2H5O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN10O2

  P = 0.264* rct(45,k) * xnew(OH) * xnew(C3H8)  &
    + 0.55* 5.58e-11 * xnew(OH) * xnew(UDCARB11)  &
    + rcphot(IDMEK,k) *4.74 * xnew(CARB14)  &
    + (rcphot(IDNO2,k) *(1.10e-02)) * xnew(UDCARB11)  &
    + rcphot(IDCH3O2H,k) * xnew(RN11OOH)
  
  L = rct(107,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(121,k) * xnew(HO2)  &
    + 6.0e-13 * xnew(RO2POOL)
  
  xnew(RN10O2) = (xold(RN10O2) + dt2 * P) / (1.0 + dt2 * L)


!-> IC3H7O2

  P = 0.736* rct(45,k) * xnew(OH) * xnew(C3H8)
  
  L = rct(108,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(118,k) * xnew(HO2)  &
    + 4.0e-14 * xnew(RO2POOL)
  
  xnew(IC3H7O2) = (xold(IC3H7O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN13O2

  P = rct(46,k) * xnew(OH) * xnew(NC4H10)  &
    + 1.65e-11 * xnew(OH) * xnew(CARB13)  &
    + 0.55* 7.0e-11 * xnew(OH) * xnew(UDCARB14)  &
    + 0.7* rct(109,k) * xnew(RTX22O2) * xnew(NO)  &
    + rct(113,k) * xnew(RTX22O2) * xnew(NO3)  &
    + 2.5e-13 * xnew(RO2POOL) * xnew(RTX22O2)  &
    + (rcphot(IDNO2,k) *(1.10e-02)) * xnew(UDCARB14)  &
    + rcphot(IDCH3O2H,k) * xnew(RN14OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX22OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(120,k) * xnew(HO2)  &
    + 2.5e-13 * xnew(RO2POOL)
  
  xnew(RN13O2) = (xold(RN13O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN16O2

  P = 1.25e-11 * xnew(OH) * xnew(CARB16)  &
    + rcphot(IDCH3O2H,k) * xnew(RN17OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 2.5e-13 * xnew(RO2POOL)
  
  xnew(RN16O2) = (xold(RN16O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN19O2

  P = 0.065* rct(109,k) * xnew(RTN28O2) * xnew(NO)  &
    + 0.065* rct(109,k) * xnew(RTX28O2) * xnew(NO)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(123,k) * xnew(HO2)  &
    + 2.5e-13 * xnew(RO2POOL)
  
  xnew(RN19O2) = (xold(RN19O2) + dt2 * P) / (1.0 + dt2 * L)


!-> HOCH2CH2O2

  P = 0.05* rct(34,k) * xnew(OH) * xnew(C2H5OH)  &
    + rct(47,k) * xnew(OH) * xnew(C2H4)  &
    + 1.5e-12 * xnew(OH) * xnew(ANHY)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(119,k) * xnew(HO2)  &
    + 2.0e-12 * xnew(RO2POOL)
  
  xnew(HOCH2CH2O2) = (xold(HOCH2CH2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN9O2

  P = 0.51* rct(35,k) * xnew(OH) * xnew(NPROPOL)  &
    + 0.139* rct(36,k) * xnew(OH) * xnew(IPROPOL)  &
    + rct(50,k) * xnew(OH) * xnew(C3H6)  &
    + 2.4e-13 * xnew(OH) * xnew(CARB9)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(121,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RN9O2) = (xold(RN9O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN12O2

  P = rct(53,k) * xnew(OH) * xnew(TBUT2ENE)  &
    + 1.38e-12 * xnew(OH) * xnew(CARB12)  &
    + 4.79e-12 * xnew(OH) * xnew(CCARB12)  &
    + rct(109,k) * xnew(RN13AO2) * xnew(NO)  &
    + rct(113,k) * xnew(RN13AO2) * xnew(NO3)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN13AO2)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(120,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RN12O2) = (xold(RN12O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN15O2

  P = rct(109,k) * xnew(RN16AO2) * xnew(NO)  &
    + rct(113,k) * xnew(RN16AO2) * xnew(NO3)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN16AO2)  &
    + rcphot(IDHNO3,k) *2.62 * xnew(RN16NO3)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RN15O2) = (xold(RN15O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NRN6O2

  P = rct(48,k) * xnew(NO3) * xnew(C2H4)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(128,k) * xnew(HO2)  &
    + 6.0e-13 * xnew(RO2POOL)
  
  xnew(NRN6O2) = (xold(NRN6O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NRN9O2

  P = rct(51,k) * xnew(NO3) * xnew(C3H6)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(121,k) * xnew(HO2)  &
    + 2.3e-13 * xnew(RO2POOL)
  
  xnew(NRN9O2) = (xold(NRN9O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NRN12O2

  P = 3.9e-13 * xnew(NO3) * xnew(TBUT2ENE)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(120,k) * xnew(HO2)  &
    + 2.5e-13 * xnew(RO2POOL)
  
  xnew(NRN12O2) = (xold(NRN12O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN11O2

  P = rct(62,k) * xnew(MEK) * xnew(OH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(120,k) * xnew(HO2)  &
    + 1.4e-12 * xnew(RO2POOL)
  
  xnew(RN11O2) = (xold(RN11O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN14O2

  P = 1.87e-11 * xnew(OH) * xnew(CARB14)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 1.4e-12 * xnew(RO2POOL)
  
  xnew(RN14O2) = (xold(RN14O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN8O2

  P = rct(61,k) * xnew(OH) * xnew(CH3COCH3)  &
    + 0.077* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + rct(109,k) * xnew(RTN10O2) * xnew(NO)  &
    + rct(110,k) * xnew(CO2C3CO3) * xnew(NO)  &
    + rct(113,k) * xnew(RTN10O2) * xnew(NO3)  &
    + rct(115,k) * xnew(CO2C3CO3) * xnew(NO3)  &
    + 0.44* rct(129,k) * xnew(CO2C3CO3) * xnew(HO2)  &
    + 2.0e-12 * xnew(RO2POOL) * xnew(RTN10O2)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(CO2C3CO3)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(CO2C3CHO)  &
    + rcphot(IDCH3COY,k) * xnew(CARB12)  &
    + rcphot(IDMEK,k) *3.00 * xnew(CARB13)  &
    + rcphot(IDMEK,k) *3.35 * xnew(CARB16)  &
    + rcphot(IDMEK,k) * xnew(UCARB13)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN10OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(CO2C3CO3H)  &
    + rcphot(IDMEK,k) * xnew(CO2C3CO3H)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(121,k) * xnew(HO2)  &
    + 1.4e-12 * xnew(RO2POOL)
  
  xnew(RN8O2) = (xold(RN8O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN17O2

  P = 1.19e-10 * xnew(OH) * xnew(HCC7CO)  &
    + 7.06e-11 * xnew(OH) * xnew(RN17OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(123,k) * xnew(HO2)  &
    + 1.4e-12 * xnew(RO2POOL)
  
  xnew(RN17O2) = (xold(RN17O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN13AO2

  P = 0.337* rct(109,k) * xnew(RTX24O2) * xnew(NO)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(120,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RN13AO2) = (xold(RN13AO2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN16AO2

  P = rct(109,k) * xnew(RN17O2) * xnew(NO)  &
    + rct(113,k) * xnew(RN17O2) * xnew(NO3)  &
    + 1.4e-12 * xnew(RO2POOL) * xnew(RN17O2)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RN16AO2) = (xold(RN16AO2) + dt2 * P) / (1.0 + dt2 * L)


!-> RN15AO2

  P = 1.0e-12 * xnew(OH) * xnew(TNCARB15)  &
    + 0.877* rct(109,k) * xnew(RN16O2) * xnew(NO)  &
    + rct(113,k) * xnew(RN16O2) * xnew(NO3)  &
    + 2.5e-13 * xnew(RO2POOL) * xnew(RN16O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RN16OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RN15AO2) = (xold(RN15AO2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5CO3

  P = rct(41,k) * xnew(OH) * xnew(C2H5CHO)  &
    + rct(42,k) * xnew(NO3) * xnew(C2H5CHO)  &
    + 4.42e-12 * xnew(OH) * xnew(C2H5CO3H)  &
    + rct(102,k) * xnew(PPN)  &
    + rct(109,k) * xnew(RN14O2) * xnew(NO)  &
    + rct(113,k) * xnew(RN14O2) * xnew(NO3)  &
    + 1.4e-12 * xnew(RO2POOL) * xnew(RN14O2)
  
  L = rct(100,k) * xnew(NO2)  &
    + rct(110,k) * xnew(NO)  &
    + rct(114,k) * xnew(NO3)  &
    + rct(129,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(C2H5CO3) = (xold(C2H5CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HOCH2CO3

  P = 0.8* 8.0e-12 * xnew(OH) * xnew(HOCH2CHO)  &
    + rct(40,k) * xnew(NO3) * xnew(HOCH2CHO)  &
    + 6.19e-12 * xnew(OH) * xnew(HOCH2CO3H)  &
    + rct(101,k) * xnew(PHAN)  &
    + 0.17* rct(109,k) * xnew(RU12O2) * xnew(NO)  &
    + 0.17* rct(113,k) * xnew(RU12O2) * xnew(NO3)  &
    + 0.17* 4.2e-13 * xnew(RO2POOL) * xnew(RU12O2)  &
    + 0.17* rcphot(IDCH3O2H,k) * xnew(RU12OOH)
  
  L = rct(100,k) * xnew(NO2)  &
    + rct(110,k) * xnew(NO)  &
    + rct(114,k) * xnew(NO3)  &
    + rct(129,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(HOCH2CO3) = (xold(HOCH2CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> CO2C3CO3

  P = 7.15e-11 * xnew(OH) * xnew(CO2C3CHO)  &
    + rct(76,k) * xnew(CO2C3CHO) * xnew(NO3)  &
    + 4.18e-12 * xnew(OH) * xnew(CO2C3CO3H)  &
    + rct(101,k) * xnew(CO2C3PAN)
  
  L = rct(100,k) * xnew(NO2)  &
    + rct(110,k) * xnew(NO)  &
    + rct(115,k) * xnew(NO3)  &
    + rct(129,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(CO2C3CO3) = (xold(CO2C3CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RU12O2

  P = 6.42e-11 * xnew(OH) * xnew(UCARB12)  &
    + rct(77,k) * xnew(NO3) * xnew(UCARB12)  &
    + 1.16e-11 * xnew(OH) * xnew(IEPOX)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.0474 * xnew(UCARB12)
  
  L = rct(98,k)  &
    + rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 4.2e-13 * xnew(RO2POOL)
  
  xnew(RU12O2) = (xold(RU12O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RU10O2

  P = rct(71,k) * xnew(OH) * xnew(MVK)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(120,k) * xnew(HO2)  &
    + 1.22e-12 * xnew(RO2POOL)
  
  xnew(RU10O2) = (xold(RU10O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NRU14O2

  P = rct(59,k) * xnew(NO3) * xnew(C5H8)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 1.3e-12 * xnew(RO2POOL)
  
  xnew(NRU14O2) = (xold(NRU14O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NRU12O2

  P = 4.16e-11 * xnew(OH) * xnew(NUCARB12)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 9.6e-13 * xnew(RO2POOL)
  
  xnew(NRU12O2) = (xold(NRU12O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RA13O2

  P = 0.47* rct(56,k) * xnew(OH) * xnew(BENZENE)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(123,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RA13O2) = (xold(RA13O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RA16O2

  P = 0.82* rct(57,k) * xnew(OH) * xnew(TOLUENE)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(124,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RA16O2) = (xold(RA16O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RA19AO2

  P = 0.7* 1.36e-11 * xnew(OH) * xnew(OXYL)  &
    + 7.62e-11 * xnew(OH) * xnew(RA19OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(125,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RA19AO2) = (xold(RA19AO2) + dt2 * P) / (1.0 + dt2 * L)


!-> RA19CO2

  P = 0.3* 1.36e-11 * xnew(OH) * xnew(OXYL)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(125,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RA19CO2) = (xold(RA19CO2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN28O2

  P = rct(63,k) * xnew(APINENE) * xnew(OH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(127,k) * xnew(HO2)  &
    + 2.85e-13 * xnew(RO2POOL)
  
  xnew(RTN28O2) = (xold(RTN28O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN26O2

  P = 0.8* rct(65,k) * xnew(APINENE) * xnew(O3)  &
    + 4.2e-11 * xnew(OH) * xnew(TNCARB26)  &
    + 3.8e-14 * xnew(NO3) * xnew(TNCARB26)  &
    + 1.2e-11 * xnew(OH) * xnew(RTN26OOH)  &
    + rct(101,k) * xnew(RTN26PAN)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(TNCARB26)
  
  L = rct(104,k) * xnew(NO2)  &
    + rct(110,k) * xnew(NO)  &
    + rct(114,k) * xnew(NO3)  &
    + rct(129,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(RTN26O2) = (xold(RTN26O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NRTN28O2

  P = rct(64,k) * xnew(APINENE) * xnew(NO3)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(127,k) * xnew(HO2)  &
    + 1.0e-13 * xnew(RO2POOL)
  
  xnew(NRTN28O2) = (xold(NRTN28O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN25O2

  P = 6.65e-12 * xnew(OH) * xnew(RCOOH25)  &
    + 1.66e-11 * xnew(OH) * xnew(RTN25OOH)  &
    + rct(110,k) * xnew(RTN26O2) * xnew(NO)  &
    + rct(114,k) * xnew(RTN26O2) * xnew(NO3)  &
    + 0.44* rct(129,k) * xnew(RTN26O2) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(RTN26O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN26OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(126,k) * xnew(HO2)  &
    + 1.3e-12 * xnew(RO2POOL)
  
  xnew(RTN25O2) = (xold(RTN25O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN24O2

  P = 1.05e-11 * xnew(OH) * xnew(RTN24OOH)  &
    + 0.84* rct(109,k) * xnew(RTN25O2) * xnew(NO)  &
    + rct(113,k) * xnew(RTN25O2) * xnew(NO3)  &
    + 1.3e-12 * xnew(RO2POOL) * xnew(RTN25O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN25OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(126,k) * xnew(HO2)  &
    + 6.7e-15 * xnew(RO2POOL)
  
  xnew(RTN24O2) = (xold(RTN24O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN23O2

  P = 2.05e-11 * xnew(OH) * xnew(RTN23OOH)  &
    + rct(109,k) * xnew(RTN24O2) * xnew(NO)  &
    + rct(113,k) * xnew(RTN24O2) * xnew(NO3)  &
    + 6.7e-15 * xnew(RO2POOL) * xnew(RTN24O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN24OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(126,k) * xnew(HO2)  &
    + 6.7e-15 * xnew(RO2POOL)
  
  xnew(RTN23O2) = (xold(RTN23O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN14O2

  P = 8.69e-11 * xnew(OH) * xnew(RTN14OOH)  &
    + rct(109,k) * xnew(RTN23O2) * xnew(NO)  &
    + rct(113,k) * xnew(RTN23O2) * xnew(NO3)  &
    + 6.7e-15 * xnew(RO2POOL) * xnew(RTN23O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN23OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(123,k) * xnew(HO2)  &
    + 8.8e-13 * xnew(RO2POOL)
  
  xnew(RTN14O2) = (xold(RTN14O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN10O2

  P = 1.0e-10 * xnew(OH) * xnew(TNCARB10)  &
    + rct(78,k) * xnew(NO3) * xnew(TNCARB10)  &
    + 4.23e-12 * xnew(OH) * xnew(RTN10OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 2.0e-12 * xnew(RO2POOL)
  
  xnew(RTN10O2) = (xold(RTN10O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX28O2

  P = rct(66,k) * xnew(BPINENE) * xnew(OH)  &
    + 2.0e-11 * xnew(OH) * xnew(RTX28OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(127,k) * xnew(HO2)  &
    + 2.0e-12 * xnew(RO2POOL)
  
  xnew(RTX28O2) = (xold(RTX28O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX24O2

  P = 0.35* rct(67,k) * xnew(BPINENE) * xnew(O3)  &
    + 1.55e-11 * xnew(OH) * xnew(TXCARB24)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(126,k) * xnew(HO2)  &
    + 2.5e-13 * xnew(RO2POOL)
  
  xnew(RTX24O2) = (xold(RTX24O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX22O2

  P = 4.55e-12 * xnew(OH) * xnew(TXCARB22)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(126,k) * xnew(HO2)  &
    + 2.5e-13 * xnew(RO2POOL)
  
  xnew(RTX22O2) = (xold(RTX22O2) + dt2 * P) / (1.0 + dt2 * L)


!-> NRTX28O2

  P = 2.51e-12 * xnew(BPINENE) * xnew(NO3)  &
    + 9.58e-12 * xnew(OH) * xnew(NRTX28OOH)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(127,k) * xnew(HO2)  &
    + 9.2e-14 * xnew(RO2POOL)
  
  xnew(NRTX28O2) = (xold(NRTX28O2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACRO2

  P = 0.55* rct(69,k) * xnew(OH) * xnew(MACR)
  
  L = rct(99,k)  &
    + rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(120,k) * xnew(HO2)  &
    + 3.6e-13 * xnew(RO2POOL)
  
  xnew(MACRO2) = (xold(MACRO2) + dt2 * P) / (1.0 + dt2 * L)


!-> MACO3

  P = 0.45* rct(69,k) * xnew(OH) * xnew(MACR)  &
    + 3.4e-15 * xnew(NO3) * xnew(MACR)  &
    + rct(103,k) * xnew(MPAN)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)
  
  L = rct(100,k) * xnew(NO2)  &
    + rct(111,k) * xnew(NO)  &
    + rct(115,k) * xnew(NO3)  &
    + rct(129,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(MACO3) = (xold(MACO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HCOCO3

  P = 1.58e-11 * xnew(OH) * xnew(HCOCO3H)  &
    + rct(96,k) * xnew(HCOCO)
  
  L = rct(100,k) * xnew(NO2)  &
    + rct(110,k) * xnew(NO)  &
    + rct(115,k) * xnew(NO3)  &
    + rct(129,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(HCOCO3) = (xold(HCOCO3) + dt2 * P) / (1.0 + dt2 * L)


!-> CRU14O2

  P = rct(83,k) * xnew(CRU14)
  
  L = rct(84,k)  &
    + rct(91,k)  &
    + rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 2.4e-12 * xnew(RO2POOL)
  
  xnew(CRU14O2) = (xold(CRU14O2) + dt2 * P) / (1.0 + dt2 * L)


!-> BRU14O2

  P = 0.06* 7.14e-11 * xnew(OH) * xnew(RU14OOH)  &
    + rct(85,k) * xnew(CRU14)  &
    + rct(89,k) * xnew(TRU14)
  
  L = rct(86,k)  &
    + rct(90,k)  &
    + rct(92,k)  &
    + rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 1.26e-12 * xnew(RO2POOL)
  
  xnew(BRU14O2) = (xold(BRU14O2) + dt2 * P) / (1.0 + dt2 * L)


!-> TRU14O2

  P = rct(87,k) * xnew(TRU14)
  
  L = rct(88,k)  &
    + rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 2.4e-12 * xnew(RO2POOL)
  
  xnew(TRU14O2) = (xold(TRU14O2) + dt2 * P) / (1.0 + dt2 * L)


!-> DHPR12O2

  P = 0.5* rct(91,k) * xnew(CRU14O2)
  
  L = rct(97,k)  &
    + rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 7.6e-13 * xnew(RO2POOL)
  
  xnew(DHPR12O2) = (xold(DHPR12O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RU13O2

  P = 4.19e-11 * xnew(OH) * xnew(UCARB13)
  
  L = rct(109,k) * xnew(NO)  &
    + rct(113,k) * xnew(NO3)  &
    + rct(122,k) * xnew(HO2)  &
    + 6.0e-13 * xnew(RO2POOL)
  
  xnew(RU13O2) = (xold(RU13O2) + dt2 * P) / (1.0 + dt2 * L)


!-> RO2POOL

  ! P = 0.0
  
  ! L = 0.0
  
  xnew(RO2POOL) = sum(xnew(RO2_GROUP))


!-> NO

  P = rct(5,k) * xnew(OP) * xnew(NO2)  &
    + rct(13,k) * xnew(NO2) * xnew(NO3)  &
    + rcphot(IDNO2,k) * xnew(NO2)  &
    + rcphot(IDNO3_NO,k) * xnew(NO3)  &
    + rcphot(IDHONO,k) * xnew(HONO)  &
    + rcemis(NO,k)
  
  L = rct(4,k) * xnew(OP)  &
    + rct(9,k) * xnew(O3)  &
    + rct(11,k) * xnew(NO)  &
    + rct(11,k) * xnew(NO)  &
    + rct(12,k) * xnew(NO3)  &
    + rct(25,k) * xnew(OH)  &
    + rct(27,k) * xnew(HO2)  &
    + rct(105,k) * xnew(CH3O2)  &
    + rct(106,k) * xnew(C2H5O2)  &
    + rct(107,k) * xnew(RN10O2)  &
    + rct(108,k) * xnew(IC3H7O2)  &
    + rct(109,k) * xnew(HOCH2CH2O2)  &
    + rct(109,k) * xnew(MACRO2)  &
    + rct(109,k) * xnew(RN8O2)  &
    + rct(109,k) * xnew(RN9O2)  &
    + rct(109,k) * xnew(RN11O2)  &
    + rct(109,k) * xnew(RN12O2)  &
    + rct(109,k) * xnew(RN13O2)  &
    + rct(109,k) * xnew(RN14O2)  &
    + rct(109,k) * xnew(RN15O2)
  
  !-> NO cont.
    L = L +   rct(109,k) * xnew(RN16O2)  &
    + rct(109,k) * xnew(RN17O2)  &
    + rct(109,k) * xnew(RN19O2)  &
    + rct(109,k) * xnew(RN13AO2)  &
    + rct(109,k) * xnew(RN15AO2)  &
    + rct(109,k) * xnew(RN16AO2)  &
    + rct(109,k) * xnew(RU10O2)  &
    + rct(109,k) * xnew(RU12O2)  &
    + rct(109,k) * xnew(RU13O2)  &
    + rct(109,k) * xnew(CRU14O2)  &
    + rct(109,k) * xnew(TRU14O2)  &
    + rct(109,k) * xnew(BRU14O2)  &
    + rct(109,k) * xnew(RA13O2)  &
    + rct(109,k) * xnew(RA16O2)  &
    + rct(109,k) * xnew(RA19AO2)  &
    + rct(109,k) * xnew(RA19CO2)  &
    + rct(109,k) * xnew(RTN10O2)  &
    + rct(109,k) * xnew(RTN14O2)  &
    + rct(109,k) * xnew(RTN23O2)  &
    + rct(109,k) * xnew(RTN24O2)
  
  !-> NO cont.
    L = L +   rct(109,k) * xnew(RTN25O2)  &
    + rct(109,k) * xnew(RTN28O2)  &
    + rct(109,k) * xnew(RTX22O2)  &
    + rct(109,k) * xnew(RTX24O2)  &
    + rct(109,k) * xnew(RTX28O2)  &
    + rct(109,k) * xnew(DHPR12O2)  &
    + rct(109,k) * xnew(NRN6O2)  &
    + rct(109,k) * xnew(NRN9O2)  &
    + rct(109,k) * xnew(NRN12O2)  &
    + rct(109,k) * xnew(NRU12O2)  &
    + rct(109,k) * xnew(NRU14O2)  &
    + rct(109,k) * xnew(NRTN28O2)  &
    + rct(109,k) * xnew(NRTX28O2)  &
    + rct(110,k) * xnew(CH3CO3)  &
    + rct(110,k) * xnew(C2H5CO3)  &
    + rct(111,k) * xnew(MACO3)  &
    + rct(110,k) * xnew(HOCH2CO3)  &
    + rct(110,k) * xnew(HCOCO3)  &
    + rct(110,k) * xnew(CO2C3CO3)  &
    + rct(110,k) * xnew(RTN26O2)
  
  xnew(NO) = (xold(NO) + dt2 * P) / (1.0 + dt2 * L)


!-> NO2

  P = rct(4,k) * xnew(OP) * xnew(NO)  &
    + rct(9,k) * xnew(NO) * xnew(O3)  &
    + 2* rct(11,k) * xnew(NO) * xnew(NO)  &
    + 2* rct(12,k) * xnew(NO) * xnew(NO3)  &
    + rct(13,k) * xnew(NO2) * xnew(NO3)  &
    + rct(15,k) * xnew(N2O5)  &
    + 2.0e-11 * xnew(OH) * xnew(NO3)  &
    + rct(27,k) * xnew(HO2) * xnew(NO)  &
    + rct(29,k) * xnew(HO2NO2)  &
    + rct(30,k) * xnew(OH) * xnew(HO2NO2)  &
    + 4.0e-12 * xnew(HO2) * xnew(NO3)  &
    + rct(31,k) * xnew(OH) * xnew(HONO)  &
    + 9.0e-13 * xnew(OH) * xnew(ARNOH14)  &
    + 9.0e-14 * xnew(NO3) * xnew(ARNOH14)  &
    + 1.53e-12 * xnew(OH) * xnew(ARNOH17)  &
    + 3.13e-13 * xnew(NO3) * xnew(ARNOH17)  &
    + 5.22e-12 * xnew(OH) * xnew(NRN6OOH)  &
    + 6.5e-12 * xnew(OH) * xnew(NRN9OOH)  &
    + 7.15e-12 * xnew(OH) * xnew(NRN12OOH)  &
    + 9.5e-12 * xnew(OH) * xnew(NRTN28OOH)
  
  !-> NO2 cont.
    P = P +   rct(80,k) * xnew(OH) * xnew(CH3NO3)  &
    + rct(81,k) * xnew(OH) * xnew(C2H5NO3)  &
    + rct(82,k) * xnew(OH) * xnew(IC3H7NO3)  &
    + 6.7e-13 * xnew(OH) * xnew(NOA)  &
    + 1.09e-12 * xnew(OH) * xnew(HOC2H4NO3)  &
    + 1.31e-12 * xnew(OH) * xnew(RN9NO3)  &
    + 5.8e-13 * xnew(OH) * xnew(RN10NO3)  &
    + 1.79e-12 * xnew(OH) * xnew(RN12NO3)  &
    + 9.2e-13 * xnew(OH) * xnew(RN13NO3)  &
    + 1.03e-11 * xnew(OH) * xnew(RN15NO3)  &
    + 1.85e-12 * xnew(OH) * xnew(RN16NO3)  &
    + 9.6e-11 * xnew(OH) * xnew(RN19NO3)  &
    + 5.26e-13 * xnew(OH) * xnew(RU10NO3)  &
    + 2.5e-12 * xnew(OH) * xnew(RU12NO3)  &
    + 0.34* 3.0e-11 * xnew(OH) * xnew(RU14NO3)  &
    + 7.3e-11 * xnew(OH) * xnew(RA13NO3)  &
    + 7.16e-11 * xnew(OH) * xnew(RA16NO3)  &
    + 8.31e-11 * xnew(OH) * xnew(RA19NO3)  &
    + 2.88e-12 * xnew(OH) * xnew(RTN25NO3)  &
    + 4.35e-12 * xnew(OH) * xnew(RTN28NO3)
  
  !-> NO2 cont.
    P = P +   4.74e-12 * xnew(OH) * xnew(RTX22NO3)  &
    + 6.48e-12 * xnew(OH) * xnew(RTX24NO3)  &
    + 3.53e-12 * xnew(OH) * xnew(RTX28NO3)  &
    + 3.0e-14 * xnew(OH) * xnew(PAN)  &
    + 1.12e-12 * xnew(OH) * xnew(PHAN)  &
    + 1.27e-12 * xnew(OH) * xnew(PPN)  &
    + 5.93e-13 * xnew(OH) * xnew(CO2C3PAN)  &
    + 3.66e-12 * xnew(OH) * xnew(RTN26PAN)  &
    + rct(101,k) * xnew(PAN)  &
    + rct(102,k) * xnew(PPN)  &
    + rct(103,k) * xnew(MPAN)  &
    + rct(101,k) * xnew(PHAN)  &
    + rct(101,k) * xnew(CO2C3PAN)  &
    + rct(101,k) * xnew(RTN26PAN)  &
    + 0.999* rct(105,k) * xnew(CH3O2) * xnew(NO)  &
    + 0.967* rct(106,k) * xnew(C2H5O2) * xnew(NO)  &
    + 0.98* rct(107,k) * xnew(RN10O2) * xnew(NO)  &
    + 0.965* rct(108,k) * xnew(IC3H7O2) * xnew(NO)  &
    + 0.995* rct(109,k) * xnew(HOCH2CH2O2) * xnew(NO)  &
    + 0.987* rct(109,k) * xnew(MACRO2) * xnew(NO)
  
  !-> NO2 cont.
    P = P +   rct(109,k) * xnew(RN8O2) * xnew(NO)  &
    + 0.979* rct(109,k) * xnew(RN9O2) * xnew(NO)  &
    + rct(109,k) * xnew(RN11O2) * xnew(NO)  &
    + 0.959* rct(109,k) * xnew(RN12O2) * xnew(NO)  &
    + 0.917* rct(109,k) * xnew(RN13O2) * xnew(NO)  &
    + rct(109,k) * xnew(RN14O2) * xnew(NO)  &
    + 0.936* rct(109,k) * xnew(RN15O2) * xnew(NO)  &
    + 0.877* rct(109,k) * xnew(RN16O2) * xnew(NO)  &
    + rct(109,k) * xnew(RN17O2) * xnew(NO)  &
    + 0.722* rct(109,k) * xnew(RN19O2) * xnew(NO)  &
    + rct(109,k) * xnew(RN13AO2) * xnew(NO)  &
    + 0.975* rct(109,k) * xnew(RN15AO2) * xnew(NO)  &
    + rct(109,k) * xnew(RN16AO2) * xnew(NO)  &
    + 0.965* rct(109,k) * xnew(RU10O2) * xnew(NO)  &
    + 0.965* rct(109,k) * xnew(RU12O2) * xnew(NO)  &
    + 0.949* rct(109,k) * xnew(RU13O2) * xnew(NO)  &
    + 0.908* rct(109,k) * xnew(CRU14O2) * xnew(NO)  &
    + 0.908* rct(109,k) * xnew(TRU14O2) * xnew(NO)  &
    + 0.889* rct(109,k) * xnew(BRU14O2) * xnew(NO)  &
    + 0.918* rct(109,k) * xnew(RA13O2) * xnew(NO)
  
  !-> NO2 cont.
    P = P +   0.889* rct(109,k) * xnew(RA16O2) * xnew(NO)  &
    + 0.862* rct(109,k) * xnew(RA19AO2) * xnew(NO)  &
    + 0.862* rct(109,k) * xnew(RA19CO2) * xnew(NO)  &
    + rct(109,k) * xnew(RTN10O2) * xnew(NO)  &
    + rct(109,k) * xnew(RTN14O2) * xnew(NO)  &
    + rct(109,k) * xnew(RTN23O2) * xnew(NO)  &
    + rct(109,k) * xnew(RTN24O2) * xnew(NO)  &
    + 0.84* rct(109,k) * xnew(RTN25O2) * xnew(NO)  &
    + 0.767* rct(109,k) * xnew(RTN28O2) * xnew(NO)  &
    + 0.7* rct(109,k) * xnew(RTX22O2) * xnew(NO)  &
    + 0.843* rct(109,k) * xnew(RTX24O2) * xnew(NO)  &
    + 0.767* rct(109,k) * xnew(RTX28O2) * xnew(NO)  &
    + rct(109,k) * xnew(DHPR12O2) * xnew(NO)  &
    + 2* rct(109,k) * xnew(NRN6O2) * xnew(NO)  &
    + 2* rct(109,k) * xnew(NRN9O2) * xnew(NO)  &
    + 2* rct(109,k) * xnew(NRN12O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRU12O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRU14O2) * xnew(NO)  &
    + 2* rct(109,k) * xnew(NRTN28O2) * xnew(NO)  &
    + 2* rct(109,k) * xnew(NRTX28O2) * xnew(NO)
  
  !-> NO2 cont.
    P = P +   rct(110,k) * xnew(CH3CO3) * xnew(NO)  &
    + rct(110,k) * xnew(C2H5CO3) * xnew(NO)  &
    + rct(111,k) * xnew(MACO3) * xnew(NO)  &
    + rct(110,k) * xnew(HOCH2CO3) * xnew(NO)  &
    + rct(110,k) * xnew(HCOCO3) * xnew(NO)  &
    + rct(110,k) * xnew(CO2C3CO3) * xnew(NO)  &
    + rct(110,k) * xnew(RTN26O2) * xnew(NO)  &
    + rct(112,k) * xnew(CH3O2) * xnew(NO3)  &
    + rct(113,k) * xnew(C2H5O2) * xnew(NO3)  &
    + rct(113,k) * xnew(IC3H7O2) * xnew(NO3)  &
    + rct(113,k) * xnew(HOCH2CH2O2) * xnew(NO3)  &
    + rct(113,k) * xnew(MACRO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN8O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN9O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN10O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN11O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN12O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN13O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN15O2) * xnew(NO3)
  
  !-> NO2 cont.
    P = P +   rct(113,k) * xnew(RN16O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN17O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN19O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN13AO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN15AO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN16AO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RU10O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RU12O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RU13O2) * xnew(NO3)  &
    + rct(113,k) * xnew(CRU14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(TRU14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(BRU14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA13O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA16O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA19AO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA19CO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN10O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN23O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN24O2) * xnew(NO3)
  
  !-> NO2 cont.
    P = P +   rct(113,k) * xnew(RTN25O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN28O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTX22O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTX24O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTX28O2) * xnew(NO3)  &
    + rct(113,k) * xnew(DHPR12O2) * xnew(NO3)  &
    + 2* rct(113,k) * xnew(NRN6O2) * xnew(NO3)  &
    + 2* rct(113,k) * xnew(NRN9O2) * xnew(NO3)  &
    + 2* rct(113,k) * xnew(NRN12O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRU12O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRU14O2) * xnew(NO3)  &
    + 2* rct(113,k) * xnew(NRTN28O2) * xnew(NO3)  &
    + 2* rct(113,k) * xnew(NRTX28O2) * xnew(NO3)  &
    + rct(114,k) * xnew(CH3CO3) * xnew(NO3)  &
    + rct(114,k) * xnew(C2H5CO3) * xnew(NO3)  &
    + rct(115,k) * xnew(MACO3) * xnew(NO3)  &
    + rct(114,k) * xnew(HOCH2CO3) * xnew(NO3)  &
    + rct(115,k) * xnew(HCOCO3) * xnew(NO3)  &
    + rct(115,k) * xnew(CO2C3CO3) * xnew(NO3)  &
    + rct(114,k) * xnew(RTN26O2) * xnew(NO3)
  
  !-> NO2 cont.
    P = P +   6.0e-13 * xnew(RO2POOL) * xnew(NRN6O2)  &
    + 2.3e-13 * xnew(RO2POOL) * xnew(NRN9O2)  &
    + 2.5e-13 * xnew(RO2POOL) * xnew(NRN12O2)  &
    + 1.0e-13 * xnew(RO2POOL) * xnew(NRTN28O2)  &
    + 9.2e-14 * xnew(RO2POOL) * xnew(NRTX28O2)  &
    + rcphot(IDNO3_NO2,k) * xnew(NO3)  &
    + rcphot(IDHNO3,k) * xnew(HNO3)  &
    + rcphot(IDMEK,k) *132.76 * xnew(NUCARB12)  &
    + rcphot(IDCH3O2H,k) * xnew(NRN6OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRN9OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRN12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRTN28OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRTX28OOH)  &
    + rcphot(IDMEK,k) *0.293 * xnew(CH3NO3)  &
    + rcphot(IDHNO3,k) *1.932 * xnew(C2H5NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(IC3H7NO3)  &
    + rcphot(IDMEK,k) *8.3 * xnew(NOA)  &
    + rcphot(IDHNO3,k) *2.62 * xnew(RN10NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RN13NO3)  &
    + rcphot(IDHNO3,k) *2.62 * xnew(RN16NO3)
  
  !-> NO2 cont.
    P = P +   rcphot(IDMEK,k) *0.778 * xnew(RN19NO3)  &
    + rcphot(IDMEK,k) *13.3 * xnew(RU10NO3)  &
    + rcphot(IDMEK,k) *83. * xnew(RU12NO3)  &
    + rcphot(IDMEK,k) * xnew(RU13NO3)  &
    + rcphot(IDHNO3,k) *2.62 * xnew(RU14NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RA13NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RA16NO3)  &
    + rcphot(IDHCHO_H,k) *0.234 * xnew(RA19NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RTX24NO3)  &
    + rcemis(NO2,k)
  
  L = rct(5,k) * xnew(OP)  &
    + rct(6,k) * xnew(OP)  &
    + rct(10,k) * xnew(O3)  &
    + rct(13,k) * xnew(NO3)  &
    + rct(14,k) * xnew(NO3)  &
    + rct(26,k) * xnew(OH)  &
    + rct(28,k) * xnew(HO2)  &
    + 2.08e-12 * xnew(RAROH14)  &
    + 2.08e-12 * xnew(RAROH17)  &
    + rct(100,k) * xnew(CH3CO3)  &
    + rct(100,k) * xnew(C2H5CO3)  &
    + rct(100,k) * xnew(MACO3)  &
    + rct(100,k) * xnew(HOCH2CO3)  &
    + rct(100,k) * xnew(HCOCO3)  &
    + rct(100,k) * xnew(CO2C3CO3)  &
    + rct(104,k) * xnew(RTN26O2)  &
    + rcphot(IDNO2,k)
  
  xnew(NO2) = (xold(NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> CO

  P = rct(37,k) * xnew(OH) * xnew(HCHO)  &
    + rct(38,k) * xnew(NO3) * xnew(HCHO)  &
    + 0.63* rct(49,k) * xnew(O3) * xnew(C2H4)  &
    + 0.56* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 0.57* rct(54,k) * xnew(O3) * xnew(TBUT2ENE)  &
    + 0.364* rct(55,k) * xnew(OH) * xnew(C2H2)  &
    + 0.364* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.2* rct(67,k) * xnew(BPINENE) * xnew(O3)  &
    + 1.23e-11 * xnew(OH) * xnew(HCOCO2H)  &
    + 0.652* rct(70,k) * xnew(O3) * xnew(MACR)  &
    + 0.56* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + rct(74,k) * xnew(OH) * xnew(MGLYOX)  &
    + 5.0e-16 * xnew(NO3) * xnew(MGLYOX)  &
    + 0.75* 2.4e-17 * xnew(O3) * xnew(UCARB12)  &
    + 0.505* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 5.78e-11 * xnew(OH) * xnew(HUCARB9)  &
    + 5.2e-11 * xnew(OH) * xnew(HPUCARB12)  &
    + rct(77,k) * xnew(NO3) * xnew(HPUCARB12)  &
    + 3.64e-11 * xnew(OH) * xnew(DHPCARB9)  &
    + 3.84e-11 * xnew(OH) * xnew(RU10OOH)
  
  !-> CO cont.
    P = P +   3.5e-11 * xnew(OH) * xnew(RU12OOH)  &
    + 5.64e-11 * xnew(OH) * xnew(DHPR12OOH)  &
    + 5.26e-13 * xnew(OH) * xnew(RU10NO3)  &
    + 3.0e-14 * xnew(OH) * xnew(PAN)  &
    + 1.12e-12 * xnew(OH) * xnew(PHAN)  &
    + 1.27e-12 * xnew(OH) * xnew(PPN)  &
    + 5.93e-13 * xnew(OH) * xnew(CO2C3PAN)  &
    + 0.22* 2.9e-11 * xnew(OH) * xnew(MPAN)  &
    + 2* rct(93,k) * xnew(HCOCO)  &
    + 2* rct(94,k) * xnew(HCOCO)  &
    + rct(95,k) * xnew(HCOCO)  &
    + rct(97,k) * xnew(DHPR12O2)  &
    + rct(98,k) * xnew(RU12O2)  &
    + rct(99,k) * xnew(MACRO2)  &
    + rct(100,k) * xnew(HCOCO3) * xnew(NO2)  &
    + 0.987* rct(109,k) * xnew(MACRO2) * xnew(NO)  &
    + rct(109,k) * xnew(RTN10O2) * xnew(NO)  &
    + 0.5* rct(109,k) * xnew(NRU12O2) * xnew(NO)  &
    + 0.65* rct(111,k) * xnew(MACO3) * xnew(NO)  &
    + rct(110,k) * xnew(HCOCO3) * xnew(NO)
  
  !-> CO cont.
    P = P +   rct(113,k) * xnew(MACRO2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN10O2) * xnew(NO3)  &
    + 0.5* rct(113,k) * xnew(NRU12O2) * xnew(NO3)  &
    + 0.65* rct(115,k) * xnew(MACO3) * xnew(NO3)  &
    + rct(115,k) * xnew(HCOCO3) * xnew(NO3)  &
    + 0.286* rct(129,k) * xnew(MACO3) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(HCOCO3) * xnew(HO2)  &
    + 3.6e-13 * xnew(RO2POOL) * xnew(MACRO2)  &
    + 2.0e-12 * xnew(RO2POOL) * xnew(RTN10O2)  &
    + 0.5* 9.6e-13 * xnew(RO2POOL) * xnew(NRU12O2)  &
    + 0.65* 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)  &
    + 0.7* 1.0e-11 * xnew(RO2POOL) * xnew(HCOCO3)  &
    + rcphot(IDHCHO_H,k) * xnew(HCHO)  &
    + rcphot(IDHCHO_H2,k) * xnew(HCHO)  &
    + rcphot(IDCH3CHO,k) * xnew(CH3CHO)  &
    + rcphot(IDMEK,k) *5.2 * xnew(C2H5CHO)  &
    + rcphot(IDRCOCHO,k) * xnew(HCOCO2H)  &
    + 0.825* rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + 2* rcphot(IDRCOCHO,k) *0.46 * xnew(GLYOX)
  
  !-> CO cont.
    P = P +   rcphot(IDRCOCHO,k) *0.07 * xnew(GLYOX)  &
    + 2* rcphot(IDCH3O2H,k) *4.79 * xnew(GLYOX)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HOCH2CHO)  &
    + 0.75* rcphot(IDHCHO_H2,k) *0.0474 * xnew(UCARB12)  &
    + rcphot(IDHCHO_H2,k) *6 * xnew(HUCARB9)  &
    + 1.425* rcphot(IDHCHO_H2,k) *12 * xnew(HPUCARB12)  &
    + rcphot(IDHCHO_H,k) *1.67 * xnew(DHCARB9)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(DHPCARB9)  &
    + rcphot(IDCH3COY,k) *0.5 * xnew(TNCARB10)  &
    + 1.5* rcphot(IDMEK,k) *132.76 * xnew(NUCARB12)  &
    + 0.254* rcphot(IDCH3O2H,k) * xnew(RU10OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN10OOH)  &
    + rcphot(IDHCHO_H,k) *1.67 * xnew(DHPR12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(HCOCO3H)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HCOCO3H)  &
    + rcemis(CO,k)
  
  L = rct(19,k) * xnew(OH)
  
  xnew(CO) = (xold(CO) + dt2 * P) / (1.0 + dt2 * L)


!-> O3

  P = rct(1,k) * xnew(OP)  &
    + rct(2,k) * xnew(OP)  &
    + 0.162* rct(129,k) * xnew(CH3CO3) * xnew(HO2)  &
    + 0.15* rct(129,k) * xnew(HCOCO3) * xnew(HO2)
  
  L = rct(3,k) * xnew(OP)  &
    + rct(9,k) * xnew(NO)  &
    + rct(10,k) * xnew(NO2)  &
    + rct(17,k) * xnew(OH)  &
    + rct(21,k) * xnew(HO2)  &
    + rct(49,k) * xnew(C2H4)  &
    + rct(52,k) * xnew(C3H6)  &
    + rct(54,k) * xnew(TBUT2ENE)  &
    + rct(60,k) * xnew(C5H8)  &
    + rct(65,k) * xnew(APINENE)  &
    + rct(67,k) * xnew(BPINENE)  &
    + rct(70,k) * xnew(MACR)  &
    + rct(72,k) * xnew(MVK)  &
    + 2.4e-17 * xnew(UCARB12)  &
    + 1.19e-17 * xnew(UCARB13)  &
    + 2.4e-17 * xnew(HPUCARB12)  &
    + rcphot(IDO3_O1D,k)  &
    + rcphot(IDO3_O3P,k)  &
    + AQRCK(ICLRC2,K) * xnew(SO2)
  
  xnew(O3) = (xold(O3) + dt2 * P) / (1.0 + dt2 * L)


!-> O3_REG

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for O3_REG! xnew(O3_REG) = max(0.0, xold(O3_REG))


!-> NO3

  P = rct(6,k) * xnew(OP) * xnew(NO2)  &
    + rct(10,k) * xnew(NO2) * xnew(O3)  &
    + rct(15,k) * xnew(N2O5)  &
    + rct(32,k) * xnew(OH) * xnew(HNO3)  &
    + 2.9e-11 * xnew(OH) * xnew(MPAN)  &
    + rct(100,k) * xnew(HCOCO3) * xnew(NO2)
  
  L = rct(12,k) * xnew(NO)  &
    + rct(13,k) * xnew(NO2)  &
    + rct(14,k) * xnew(NO2)  &
    + 2.0e-11 * xnew(OH)  &
    + 4.0e-12 * xnew(HO2)  &
    + rct(38,k) * xnew(HCHO)  &
    + rct(40,k) * xnew(CH3CHO)  &
    + rct(42,k) * xnew(C2H5CHO)  &
    + rct(48,k) * xnew(C2H4)  &
    + rct(51,k) * xnew(C3H6)  &
    + 3.9e-13 * xnew(TBUT2ENE)  &
    + 3.78e-12 * xnew(AROH14)  &
    + 9.0e-14 * xnew(ARNOH14)  &
    + 1.25e-11 * xnew(AROH17)  &
    + 3.13e-13 * xnew(ARNOH17)  &
    + rct(59,k) * xnew(C5H8)  &
    + rct(64,k) * xnew(APINENE)  &
    + 2.51e-12 * xnew(BPINENE)  &
    + 3.4e-15 * xnew(MACR)  &
    + 4.0e-16 * xnew(GLYOX)
  
  !-> NO3 cont.
    L = L +   5.0e-16 * xnew(MGLYOX)  &
    + rct(40,k) * xnew(HOCH2CHO)  &
    + rct(76,k) * xnew(CO2C3CHO)  &
    + rct(77,k) * xnew(UCARB12)  &
    + rct(77,k) * xnew(HPUCARB12)  &
    + rct(78,k) * xnew(TNCARB10)  &
    + 3.8e-14 * xnew(TNCARB26)  &
    + rct(112,k) * xnew(CH3O2)  &
    + rct(113,k) * xnew(C2H5O2)  &
    + rct(113,k) * xnew(IC3H7O2)  &
    + rct(113,k) * xnew(HOCH2CH2O2)  &
    + rct(113,k) * xnew(MACRO2)  &
    + rct(113,k) * xnew(RN8O2)  &
    + rct(113,k) * xnew(RN9O2)  &
    + rct(113,k) * xnew(RN10O2)  &
    + rct(113,k) * xnew(RN11O2)  &
    + rct(113,k) * xnew(RN12O2)  &
    + rct(113,k) * xnew(RN13O2)  &
    + rct(113,k) * xnew(RN14O2)  &
    + rct(113,k) * xnew(RN15O2)
  
  !-> NO3 cont.
    L = L +   rct(113,k) * xnew(RN16O2)  &
    + rct(113,k) * xnew(RN17O2)  &
    + rct(113,k) * xnew(RN19O2)  &
    + rct(113,k) * xnew(RN13AO2)  &
    + rct(113,k) * xnew(RN15AO2)  &
    + rct(113,k) * xnew(RN16AO2)  &
    + rct(113,k) * xnew(RU10O2)  &
    + rct(113,k) * xnew(RU12O2)  &
    + rct(113,k) * xnew(RU13O2)  &
    + rct(113,k) * xnew(CRU14O2)  &
    + rct(113,k) * xnew(TRU14O2)  &
    + rct(113,k) * xnew(BRU14O2)  &
    + rct(113,k) * xnew(RA13O2)  &
    + rct(113,k) * xnew(RA16O2)  &
    + rct(113,k) * xnew(RA19AO2)  &
    + rct(113,k) * xnew(RA19CO2)  &
    + rct(113,k) * xnew(RTN10O2)  &
    + rct(113,k) * xnew(RTN14O2)  &
    + rct(113,k) * xnew(RTN23O2)  &
    + rct(113,k) * xnew(RTN24O2)
  
  !-> NO3 cont.
    L = L +   rct(113,k) * xnew(RTN25O2)  &
    + rct(113,k) * xnew(RTN28O2)  &
    + rct(113,k) * xnew(RTX22O2)  &
    + rct(113,k) * xnew(RTX24O2)  &
    + rct(113,k) * xnew(RTX28O2)  &
    + rct(113,k) * xnew(DHPR12O2)  &
    + rct(113,k) * xnew(NRN6O2)  &
    + rct(113,k) * xnew(NRN9O2)  &
    + rct(113,k) * xnew(NRN12O2)  &
    + rct(113,k) * xnew(NRU12O2)  &
    + rct(113,k) * xnew(NRU14O2)  &
    + rct(113,k) * xnew(NRTN28O2)  &
    + rct(113,k) * xnew(NRTX28O2)  &
    + rct(114,k) * xnew(CH3CO3)  &
    + rct(114,k) * xnew(C2H5CO3)  &
    + rct(115,k) * xnew(MACO3)  &
    + rct(114,k) * xnew(HOCH2CO3)  &
    + rct(115,k) * xnew(HCOCO3)  &
    + rct(115,k) * xnew(CO2C3CO3)  &
    + rct(114,k) * xnew(RTN26O2)
  
  !-> NO3 cont.
    L = L +   rcphot(IDNO3_NO,k)  &
    + rcphot(IDNO3_NO2,k)  &
    + rct(156,k)
  
  xnew(NO3) = (xold(NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> N2O5

  P = rct(14,k) * xnew(NO2) * xnew(NO3)
  
  L = rct(15,k)  &
    + rct(155,k)
  
  xnew(N2O5) = (xold(N2O5) + dt2 * P) / (1.0 + dt2 * L)


!-> H2

  P = rcphot(IDHCHO_H2,k) * xnew(HCHO)  &
    + rcphot(IDRCOCHO,k) *0.46 * xnew(GLYOX)
  
  L = rct(18,k) * xnew(OH)
  
  xnew(H2) = (xold(H2) + dt2 * P) / (1.0 + dt2 * L)


!-> H2O2

  P = rct(23,k) * xnew(HO2) * xnew(HO2)  &
    + rct(24,k) * xnew(HO2) * xnew(HO2)  &
    + 0.139* rct(49,k) * xnew(O3) * xnew(C2H4)  &
    + 0.09* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 0.0675* rct(54,k) * xnew(O3) * xnew(TBUT2ENE)  &
    + 0.385* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.075* rct(65,k) * xnew(APINENE) * xnew(O3)  &
    + 0.2* rct(67,k) * xnew(BPINENE) * xnew(O3)  &
    + 0.142* rct(70,k) * xnew(O3) * xnew(MACR)  &
    + 0.165* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + 0.165* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 0.5* rct(159,k) * xnew(HO2)
  
  L = rct(20,k) * xnew(OH)  &
    + rcphot(IDH2O2,k)  &
    + AQRCK(ICLRC1,K) * xnew(SO2)  &
    + 1.0e-5
  
  xnew(H2O2) = (xold(H2O2) + dt2 * P) / (1.0 + dt2 * L)


!-> HONO

  P = rct(25,k) * xnew(OH) * xnew(NO)
  
  L = rct(31,k) * xnew(OH)  &
    + rcphot(IDHONO,k)
  
  xnew(HONO) = (xold(HONO) + dt2 * P) / (1.0 + dt2 * L)


!-> HNO3

  P = rct(26,k) * xnew(OH) * xnew(NO2)  &
    + rct(38,k) * xnew(NO3) * xnew(HCHO)  &
    + rct(40,k) * xnew(NO3) * xnew(CH3CHO)  &
    + rct(42,k) * xnew(NO3) * xnew(C2H5CHO)  &
    + 3.78e-12 * xnew(NO3) * xnew(AROH14)  &
    + 9.0e-14 * xnew(NO3) * xnew(ARNOH14)  &
    + 1.25e-11 * xnew(NO3) * xnew(AROH17)  &
    + 3.13e-13 * xnew(NO3) * xnew(ARNOH17)  &
    + 3.4e-15 * xnew(NO3) * xnew(MACR)  &
    + 4.0e-16 * xnew(NO3) * xnew(GLYOX)  &
    + 5.0e-16 * xnew(NO3) * xnew(MGLYOX)  &
    + rct(40,k) * xnew(NO3) * xnew(HOCH2CHO)  &
    + rct(76,k) * xnew(CO2C3CHO) * xnew(NO3)  &
    + rct(77,k) * xnew(NO3) * xnew(UCARB12)  &
    + rct(77,k) * xnew(NO3) * xnew(HPUCARB12)  &
    + rct(78,k) * xnew(NO3) * xnew(TNCARB10)  &
    + 3.8e-14 * xnew(NO3) * xnew(TNCARB26)  &
    + 2* rct(155,k) * xnew(N2O5)  &
    + rct(156,k) * xnew(NO3)  &
    + rct(26,k) * xnew(shipNOx) * xnew(OH)
  
  !-> HNO3 cont.
    P = P +   3.2e-5 * xnew(shipNOx)
  
  L = rct(32,k) * xnew(OH)  &
    + rcphot(IDHNO3,k)  &
    + rct(157,k)  &
    + rct(158,k)
  
  xnew(HNO3) = (xold(HNO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2NO2

  P = rct(28,k) * xnew(HO2) * xnew(NO2)
  
  L = rct(29,k)  &
    + rct(30,k) * xnew(OH)
  
  xnew(HO2NO2) = (xold(HO2NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H6

  P = rcemis(C2H6,k)
  
  L = rct(44,k) * xnew(OH)
  
  xnew(C2H6) = (xold(C2H6) + dt2 * P) / (1.0 + dt2 * L)


!-> C3H8

  P = rcemis(C3H8,k)
  
  L = rct(45,k) * xnew(OH)
  
  xnew(C3H8) = (xold(C3H8) + dt2 * P) / (1.0 + dt2 * L)


!-> NC4H10

  P = rcemis(NC4H10,k)
  
  L = rct(46,k) * xnew(OH)
  
  xnew(NC4H10) = (xold(NC4H10) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H4

  P = rcemis(C2H4,k)
  
  L = rct(47,k) * xnew(OH)  &
    + rct(48,k) * xnew(NO3)  &
    + rct(49,k) * xnew(O3)
  
  xnew(C2H4) = (xold(C2H4) + dt2 * P) / (1.0 + dt2 * L)


!-> C3H6

  P = 0.5* rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + rcemis(C3H6,k)
  
  L = rct(50,k) * xnew(OH)  &
    + rct(51,k) * xnew(NO3)  &
    + rct(52,k) * xnew(O3)
  
  xnew(C3H6) = (xold(C3H6) + dt2 * P) / (1.0 + dt2 * L)


!-> TBUT2ENE

  P = rcemis(TBUT2ENE,k)
  
  L = rct(53,k) * xnew(OH)  &
    + 3.9e-13 * xnew(NO3)  &
    + rct(54,k) * xnew(O3)
  
  xnew(TBUT2ENE) = (xold(TBUT2ENE) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CHO

  P = 0.95* rct(34,k) * xnew(OH) * xnew(C2H5OH)  &
    + 0.545* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 1.0675* rct(54,k) * xnew(O3) * xnew(TBUT2ENE)  &
    + 0.047* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.1* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + 1.36e-11 * xnew(OH) * xnew(C2H5OOH)  &
    + 6.5e-12 * xnew(OH) * xnew(NRN9OOH)  &
    + 2* 7.15e-12 * xnew(OH) * xnew(NRN12OOH)  &
    + rct(81,k) * xnew(OH) * xnew(C2H5NO3)  &
    + 1.27e-12 * xnew(OH) * xnew(PPN)  &
    + 0.967* rct(106,k) * xnew(C2H5O2) * xnew(NO)  &
    + 0.979* rct(109,k) * xnew(RN9O2) * xnew(NO)  &
    + rct(109,k) * xnew(RN11O2) * xnew(NO)  &
    + 1.918* rct(109,k) * xnew(RN12O2) * xnew(NO)  &
    + 0.365* rct(109,k) * xnew(RN13O2) * xnew(NO)  &
    + rct(109,k) * xnew(RN14O2) * xnew(NO)  &
    + 0.936* rct(109,k) * xnew(RN15O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRN9O2) * xnew(NO)  &
    + 2* rct(109,k) * xnew(NRN12O2) * xnew(NO)  &
    + rct(113,k) * xnew(C2H5O2) * xnew(NO3)
  
  !-> CH3CHO cont.
    P = P +   rct(113,k) * xnew(RN9O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN11O2) * xnew(NO3)  &
    + 2* rct(113,k) * xnew(RN12O2) * xnew(NO3)  &
    + 0.398* rct(113,k) * xnew(RN13O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN15O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRN9O2) * xnew(NO3)  &
    + 2* rct(113,k) * xnew(NRN12O2) * xnew(NO3)  &
    + 0.8* 3.1e-13 * xnew(RO2POOL) * xnew(C2H5O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN9O2)  &
    + 1.4e-12 * xnew(RO2POOL) * xnew(RN11O2)  &
    + 2* 8.8e-13 * xnew(RO2POOL) * xnew(RN12O2)  &
    + 0.398* 2.5e-13 * xnew(RO2POOL) * xnew(RN13O2)  &
    + 1.4e-12 * xnew(RO2POOL) * xnew(RN14O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN15O2)  &
    + 2.3e-13 * xnew(RO2POOL) * xnew(NRN9O2)  &
    + 2* 2.5e-13 * xnew(RO2POOL) * xnew(NRN12O2)  &
    + rcphot(IDMEK,k) * xnew(CARB10)  &
    + rcphot(IDMEK,k) *3.00 * xnew(CARB13)  &
    + rcphot(IDCH3O2H,k) * xnew(C2H5OOH)
  
  !-> CH3CHO cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(RN9OOH)  &
    + 2* rcphot(IDCH3O2H,k) * xnew(RN12OOH)  &
    + 0.398* rcphot(IDCH3O2H,k) * xnew(RN13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN15OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRN9OOH)  &
    + 2* rcphot(IDCH3O2H,k) * xnew(NRN12OOH)  &
    + rcphot(IDHNO3,k) *1.932 * xnew(C2H5NO3)  &
    + 0.398* rcphot(IDMEK,k) *0.778 * xnew(RN13NO3)  &
    + rcemis(CH3CHO,k)
  
  L = rct(39,k) * xnew(OH)  &
    + rct(40,k) * xnew(NO3)  &
    + rcphot(IDCH3CHO,k)
  
  xnew(CH3CHO) = (xold(CH3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5CHO

  P = 0.49* rct(35,k) * xnew(OH) * xnew(NPROPOL)  &
    + 1.89e-11 * xnew(OH) * xnew(RN10OOH)  &
    + 5.8e-13 * xnew(OH) * xnew(RN10NO3)  &
    + 0.98* rct(107,k) * xnew(RN10O2) * xnew(NO)  &
    + 0.936* rct(109,k) * xnew(RN15O2) * xnew(NO)  &
    + rct(113,k) * xnew(RN10O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN15O2) * xnew(NO3)  &
    + 0.8* 6.0e-13 * xnew(RO2POOL) * xnew(RN10O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN15O2)  &
    + rcphot(IDMEK,k) *3.35 * xnew(CARB16)  &
    + rcphot(IDCH3O2H,k) * xnew(RN10OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN15OOH)  &
    + rcphot(IDHNO3,k) *2.62 * xnew(RN10NO3)  &
    + rcemis(C2H5CHO,k)
  
  L = rct(41,k) * xnew(OH)  &
    + rct(42,k) * xnew(NO3)  &
    + rcphot(IDMEK,k) *5.2
  
  xnew(C2H5CHO) = (xold(C2H5CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO2H

  P = 0.075* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 0.1125* rct(54,k) * xnew(O3) * xnew(TBUT2ENE)  &
    + 0.162* rct(129,k) * xnew(CH3CO3) * xnew(HO2)  &
    + rcemis(CH3CO2H,k)
  
  L = rct(68,k) * xnew(OH)
  
  xnew(CH3CO2H) = (xold(CH3CO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3COCH3

  P = 0.861* rct(36,k) * xnew(OH) * xnew(IPROPOL)  &
    + 2.78e-11 * xnew(OH) * xnew(IC3H7OOH)  &
    + 7.5e-11 * xnew(OH) * xnew(RTX22OOH)  &
    + rct(82,k) * xnew(OH) * xnew(IC3H7NO3)  &
    + 2.88e-12 * xnew(OH) * xnew(RTN25NO3)  &
    + 4.74e-12 * xnew(OH) * xnew(RTX22NO3)  &
    + 3.66e-12 * xnew(OH) * xnew(RTN26PAN)  &
    + 0.965* rct(108,k) * xnew(IC3H7O2) * xnew(NO)  &
    + rct(109,k) * xnew(RTN23O2) * xnew(NO)  &
    + 0.065* rct(109,k) * xnew(RTN28O2) * xnew(NO)  &
    + 0.7* rct(109,k) * xnew(RTX22O2) * xnew(NO)  &
    + 0.337* rct(109,k) * xnew(RTX24O2) * xnew(NO)  &
    + 0.065* rct(109,k) * xnew(RTX28O2) * xnew(NO)  &
    + rct(113,k) * xnew(IC3H7O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN23O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTX22O2) * xnew(NO3)  &
    + 0.8* 4.0e-14 * xnew(RO2POOL) * xnew(IC3H7O2)  &
    + 6.7e-15 * xnew(RO2POOL) * xnew(RTN23O2)  &
    + 2.5e-13 * xnew(RO2POOL) * xnew(RTX22O2)  &
    + rcphot(IDCH3O2H,k) * xnew(IC3H7OOH)
  
  !-> CH3COCH3 cont.
    P = P +   rcphot(IDCH3O2H,k) * xnew(RTN23OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX22OOH)  &
    + rcphot(IDMEK,k) *0.778 * xnew(IC3H7NO3)  &
    + rcemis(CH3COCH3,k)
  
  L = rct(61,k) * xnew(OH)  &
    + rcphot(IDACETON,k)
  
  xnew(CH3COCH3) = (xold(CH3COCH3) + dt2 * P) / (1.0 + dt2 * L)


!-> CARB14

  P = 4.21e-11 * xnew(OH) * xnew(RN16OOH)  &
    + 1.85e-12 * xnew(OH) * xnew(RN16NO3)
  
  L = 1.87e-11 * xnew(OH)  &
    + rcphot(IDMEK,k) *4.74
  
  xnew(CARB14) = (xold(CARB14) + dt2 * P) / (1.0 + dt2 * L)


!-> C5H8

  P = rcemis(C5H8,k)  &
    + RCBIO(1,K)
  
  L = rct(58,k) * xnew(OH)  &
    + rct(59,k) * xnew(NO3)  &
    + rct(60,k) * xnew(O3)
  
  xnew(C5H8) = (xold(C5H8) + dt2 * P) / (1.0 + dt2 * L)


!-> BENZENE

  P = rcemis(BENZENE,k)
  
  L = rct(56,k) * xnew(OH)
  
  xnew(BENZENE) = (xold(BENZENE) + dt2 * P) / (1.0 + dt2 * L)


!-> TOLUENE

  P = rcemis(TOLUENE,k)
  
  L = rct(57,k) * xnew(OH)
  
  xnew(TOLUENE) = (xold(TOLUENE) + dt2 * P) / (1.0 + dt2 * L)


!-> OXYL

  P = rcemis(OXYL,k)
  
  L = 1.36e-11 * xnew(OH)
  
  xnew(OXYL) = (xold(OXYL) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H2

  P = rcemis(C2H2,k)
  
  L = rct(55,k) * xnew(OH)
  
  xnew(C2H2) = (xold(C2H2) + dt2 * P) / (1.0 + dt2 * L)


!-> HCHO

  P = rct(33,k) * xnew(OH) * xnew(CH3OH)  &
    + 1.139* rct(49,k) * xnew(O3) * xnew(C2H4)  &
    + 0.545* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 0.825* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.45* rct(67,k) * xnew(BPINENE) * xnew(O3)  &
    + 0.242* rct(70,k) * xnew(O3) * xnew(MACR)  &
    + 0.645* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + 0.636* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 0.4* rct(79,k) * xnew(OH) * xnew(CH3OOH)  &
    + 2* 6.71e-11 * xnew(OH) * xnew(RU13OOH)  &
    + 2* 5.22e-12 * xnew(OH) * xnew(NRN6OOH)  &
    + 6.5e-12 * xnew(OH) * xnew(NRN9OOH)  &
    + rct(80,k) * xnew(OH) * xnew(CH3NO3)  &
    + 2* 1.6e-12 * xnew(OH) * xnew(RU13NO3)  &
    + 0.48* 3.0e-11 * xnew(OH) * xnew(RU14NO3)  &
    + 3.53e-12 * xnew(OH) * xnew(RTX28NO3)  &
    + 3.0e-14 * xnew(OH) * xnew(PAN)  &
    + 1.12e-12 * xnew(OH) * xnew(PHAN)  &
    + rct(92,k) * xnew(BRU14O2)  &
    + 0.999* rct(105,k) * xnew(CH3O2) * xnew(NO)
  
  !-> HCHO cont.
    P = P +   1.544* rct(109,k) * xnew(HOCH2CH2O2) * xnew(NO)  &
    + rct(109,k) * xnew(RN8O2) * xnew(NO)  &
    + 0.979* rct(109,k) * xnew(RN9O2) * xnew(NO)  &
    + 0.295* rct(109,k) * xnew(RU10O2) * xnew(NO)  &
    + 0.949* rct(109,k) * xnew(RU13O2) * xnew(NO)  &
    + 0.889* rct(109,k) * xnew(BRU14O2) * xnew(NO)  &
    + rct(109,k) * xnew(RTN14O2) * xnew(NO)  &
    + 0.337* rct(109,k) * xnew(RTX24O2) * xnew(NO)  &
    + 0.702* rct(109,k) * xnew(RTX28O2) * xnew(NO)  &
    + 2* rct(109,k) * xnew(NRN6O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRN9O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRTX28O2) * xnew(NO)  &
    + rct(111,k) * xnew(MACO3) * xnew(NO)  &
    + rct(110,k) * xnew(HOCH2CO3) * xnew(NO)  &
    + rct(112,k) * xnew(CH3O2) * xnew(NO3)  &
    + 1.552* rct(113,k) * xnew(HOCH2CH2O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN8O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN9O2) * xnew(NO3)  &
    + 0.3* rct(113,k) * xnew(RU10O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RU13O2) * xnew(NO3)
  
  !-> HCHO cont.
    P = P +   rct(113,k) * xnew(BRU14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTN14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RTX28O2) * xnew(NO3)  &
    + 2* rct(113,k) * xnew(NRN6O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRN9O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRTX28O2) * xnew(NO3)  &
    + rct(115,k) * xnew(MACO3) * xnew(NO3)  &
    + rct(114,k) * xnew(HOCH2CO3) * xnew(NO3)  &
    + 0.1* rct(116,k) * xnew(CH3O2) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(MACO3) * xnew(HO2)  &
    + 0.44* rct(129,k) * xnew(HOCH2CO3) * xnew(HO2)  &
    + rct(130,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 0.5* rct(131,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + 1.552* 2.0e-12 * xnew(RO2POOL) * xnew(HOCH2CH2O2)  &
    + 1.4e-12 * xnew(RO2POOL) * xnew(RN8O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN9O2)  &
    + 0.3* 1.22e-12 * xnew(RO2POOL) * xnew(RU10O2)  &
    + 6.0e-13 * xnew(RO2POOL) * xnew(RU13O2)  &
    + 1.26e-12 * xnew(RO2POOL) * xnew(BRU14O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RTN14O2)
  
  !-> HCHO cont.
    P = P +   2.0e-12 * xnew(RO2POOL) * xnew(RTX28O2)  &
    + 2* 6.0e-13 * xnew(RO2POOL) * xnew(NRN6O2)  &
    + 2.3e-13 * xnew(RO2POOL) * xnew(NRN9O2)  &
    + 9.2e-14 * xnew(RO2POOL) * xnew(NRTX28O2)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)  &
    + 1.0e-11 * xnew(RO2POOL) * xnew(HOCH2CO3)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + rcphot(IDRCOCHO,k) *0.07 * xnew(GLYOX)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(HOCH2CHO)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(CO2C3CHO)  &
    + rcphot(IDCH3O2H,k) *2 * xnew(DHPCARB9)  &
    + rcphot(IDCH3O2H,k) * xnew(CH3OOH)  &
    + 2* rcphot(IDCH3O2H,k) * xnew(HOC2H4OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN8OOH)  &
    + rcphot(IDMEK,k) * xnew(RN8OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RN9OOH)  &
    + 0.376* rcphot(IDCH3O2H,k) * xnew(RU10OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RU13OOH)
  
  !-> HCHO cont.
    P = P +   rcphot(IDMEK,k) * xnew(RU13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN14OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX28OOH)  &
    + 2* rcphot(IDCH3O2H,k) * xnew(NRN6OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRN9OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRTX28OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(HOCH2CO3H)  &
    + rcphot(IDMEK,k) *0.293 * xnew(CH3NO3)  &
    + rcphot(IDMEK,k) *8.3 * xnew(NOA)  &
    + rcphot(IDMEK,k) * xnew(RU13NO3)  &
    + rcphot(IDHNO3,k) *2.62 * xnew(RU14NO3)  &
    + rcemis(HCHO,k)
  
  L = rct(37,k) * xnew(OH)  &
    + rct(38,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k)  &
    + rcphot(IDHCHO_H2,k)
  
  xnew(HCHO) = (xold(HCHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3

  P = rct(39,k) * xnew(OH) * xnew(CH3CHO)  &
    + rct(40,k) * xnew(NO3) * xnew(CH3CHO)  &
    + 0.154* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.098* rct(70,k) * xnew(O3) * xnew(MACR)  &
    + 0.28* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + rct(74,k) * xnew(OH) * xnew(MGLYOX)  &
    + 5.0e-16 * xnew(NO3) * xnew(MGLYOX)  &
    + 0.25* 2.4e-17 * xnew(O3) * xnew(UCARB12)  &
    + 0.158* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 3.7e-12 * xnew(OH) * xnew(CH3CO3H)  &
    + 0.3* 4.33e-12 * xnew(OH) * xnew(HMML)  &
    + rct(101,k) * xnew(PAN)  &
    + rct(109,k) * xnew(RN8O2) * xnew(NO)  &
    + rct(109,k) * xnew(RN11O2) * xnew(NO)  &
    + 0.67* rct(109,k) * xnew(RU10O2) * xnew(NO)  &
    + 0.35* rct(111,k) * xnew(MACO3) * xnew(NO)  &
    + rct(113,k) * xnew(RN8O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RN11O2) * xnew(NO3)  &
    + 0.7* rct(113,k) * xnew(RU10O2) * xnew(NO3)  &
    + 0.35* rct(115,k) * xnew(MACO3) * xnew(NO3)
  
  !-> CH3CO3 cont.
    P = P +   0.154* rct(129,k) * xnew(MACO3) * xnew(HO2)  &
    + 1.4e-12 * xnew(RO2POOL) * xnew(RN8O2)  &
    + 1.4e-12 * xnew(RO2POOL) * xnew(RN11O2)  &
    + 0.7* 1.22e-12 * xnew(RO2POOL) * xnew(RU10O2)  &
    + 0.35* 1.0e-11 * xnew(RO2POOL) * xnew(MACO3)  &
    + rcphot(IDACETON,k) * xnew(CH3COCH3)  &
    + rcphot(IDMEK,k) * xnew(MEK)  &
    + 0.175* rcphot(IDHCHO_H2,k) *0.0474 * xnew(MACR)  &
    + 0.5* rcphot(IDHCHO_H2,k) *0.078 * xnew(MVK)  &
    + rcphot(IDRCOCHO,k) * xnew(MGLYOX)  &
    + rcphot(IDMEK,k) * xnew(ACETOL)  &
    + 2* rcphot(IDCH3COY,k) * xnew(CARB9)  &
    + rcphot(IDMEK,k) * xnew(CARB10)  &
    + rcphot(IDCH3COY,k) * xnew(CARB12)  &
    + rcphot(IDMEK,k) *4.74 * xnew(CARB14)  &
    + 0.25* rcphot(IDHCHO_H2,k) *0.0474 * xnew(UCARB12)  &
    + rcphot(IDMEK,k) * xnew(UCARB13)  &
    + 0.075* rcphot(IDHCHO_H2,k) *12 * xnew(HPUCARB12)  &
    + 2* rcphot(IDCH3COY,k) *0.5 * xnew(TNCARB10)  &
    + rcphot(IDCH3O2H,k) * xnew(RN8OOH)
  
  !-> CH3CO3 cont.
    P = P +   rcphot(IDMEK,k) * xnew(RN8OOH)  &
    + 0.520* rcphot(IDCH3O2H,k) * xnew(RU10OOH)  &
    + rcphot(IDMEK,k) *8.3 * xnew(NOA)  &
    + rcphot(IDMEK,k) *13.3 * xnew(RU10NO3)
  
  L = rct(100,k) * xnew(NO2)  &
    + rct(110,k) * xnew(NO)  &
    + rct(114,k) * xnew(NO3)  &
    + rct(129,k) * xnew(HO2)  &
    + 1.0e-11 * xnew(RO2POOL)
  
  xnew(CH3CO3) = (xold(CH3CO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HCOOH

  P = 0.231* rct(49,k) * xnew(O3) * xnew(C2H4)  &
    + 0.075* rct(52,k) * xnew(O3) * xnew(C3H6)  &
    + 0.364* rct(55,k) * xnew(OH) * xnew(C2H2)  &
    + 0.175* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.204* rct(70,k) * xnew(O3) * xnew(MACR)  &
    + 0.075* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + 0.044* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 0.3* 4.33e-12 * xnew(OH) * xnew(HMML)  &
    + rcemis(HCOOH,k)
  
  L = 4.5e-13 * xnew(OH)
  
  xnew(HCOOH) = (xold(HCOOH) + dt2 * P) / (1.0 + dt2 * L)


!-> HCOCO2H

  P = 0.15* rct(129,k) * xnew(HCOCO3) * xnew(HO2)  &
    + 0.3* 1.0e-11 * xnew(RO2POOL) * xnew(HCOCO3)
  
  L = 1.23e-11 * xnew(OH)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(HCOCO2H) = (xold(HCOCO2H) + dt2 * P) / (1.0 + dt2 * L)


!-> HCOCO3H

  P = 0.41* rct(129,k) * xnew(HCOCO3) * xnew(HO2)
  
  L = 1.58e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDHCHO_H,k) *0.6
  
  xnew(HCOCO3H) = (xold(HCOCO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> MEK

  P = 3.57e-11 * xnew(OH) * xnew(RN13OOH)  &
    + 9.2e-13 * xnew(OH) * xnew(RN13NO3)  &
    + 0.552* rct(109,k) * xnew(RN13O2) * xnew(NO)  &
    + 0.602* rct(113,k) * xnew(RN13O2) * xnew(NO3)  &
    + 0.602* 2.5e-13 * xnew(RO2POOL) * xnew(RN13O2)  &
    + 0.602* rcphot(IDCH3O2H,k) * xnew(RN13OOH)  &
    + 0.602* rcphot(IDMEK,k) *0.778 * xnew(RN13NO3)  &
    + rcemis(MEK,k)
  
  L = rct(62,k) * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(MEK) = (xold(MEK) + dt2 * P) / (1.0 + dt2 * L)


!-> MVK

  P = 0.312* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.6* rct(92,k) * xnew(BRU14O2)  &
    + 0.593* rct(109,k) * xnew(BRU14O2) * xnew(NO)  &
    + 0.667* rct(113,k) * xnew(BRU14O2) * xnew(NO3)  &
    + 0.333* 1.26e-12 * xnew(RO2POOL) * xnew(BRU14O2)  &
    + 0.333* rcphot(IDHNO3,k) *2.62 * xnew(RU14NO3)
  
  L = rct(71,k) * xnew(OH)  &
    + rct(72,k) * xnew(O3)  &
    + rcphot(IDHCHO_H2,k) *0.078
  
  xnew(MVK) = (xold(MVK) + dt2 * P) / (1.0 + dt2 * L)


!-> MACR

  P = 0.468* rct(60,k) * xnew(O3) * xnew(C5H8)  &
    + 0.4* rct(92,k) * xnew(BRU14O2)  &
    + 0.296* rct(109,k) * xnew(BRU14O2) * xnew(NO)  &
    + 0.333* rct(113,k) * xnew(BRU14O2) * xnew(NO3)  &
    + 0.667* 1.26e-12 * xnew(RO2POOL) * xnew(BRU14O2)  &
    + 0.667* rcphot(IDHNO3,k) *2.62 * xnew(RU14NO3)  &
    + rcemis(MACR,k)
  
  L = rct(69,k) * xnew(OH)  &
    + 3.4e-15 * xnew(NO3)  &
    + rct(70,k) * xnew(O3)  &
    + rcphot(IDHCHO_H2,k) *0.0474
  
  xnew(MACR) = (xold(MACR) + dt2 * P) / (1.0 + dt2 * L)


!-> GLYOX

  P = 0.636* rct(55,k) * xnew(OH) * xnew(C2H2)  &
    + 0.2* 8.0e-12 * xnew(OH) * xnew(HOCH2CHO)  &
    + 0.5* 2.4e-17 * xnew(O3) * xnew(UCARB12)  &
    + 0.158* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 2.4e-17 * xnew(O3) * xnew(HPUCARB12)  &
    + 9.77e-11 * xnew(OH) * xnew(RA13OOH)  &
    + 9.64e-11 * xnew(OH) * xnew(RA16OOH)  &
    + 2.65e-11 * xnew(OH) * xnew(NRU12OOH)  &
    + 2.5e-12 * xnew(OH) * xnew(RU12NO3)  &
    + 7.3e-11 * xnew(OH) * xnew(RA13NO3)  &
    + 7.16e-11 * xnew(OH) * xnew(RA16NO3)  &
    + 0.388* rct(109,k) * xnew(RU12O2) * xnew(NO)  &
    + 0.918* rct(109,k) * xnew(RA13O2) * xnew(NO)  &
    + 0.622* rct(109,k) * xnew(RA16O2) * xnew(NO)  &
    + 0.862* rct(109,k) * xnew(RA19AO2) * xnew(NO)  &
    + rct(109,k) * xnew(DHPR12O2) * xnew(NO)  &
    + 0.5* rct(109,k) * xnew(NRU12O2) * xnew(NO)  &
    + 0.415* rct(113,k) * xnew(RU12O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA13O2) * xnew(NO3)  &
    + 0.7* rct(113,k) * xnew(RA16O2) * xnew(NO3)
  
  !-> GLYOX cont.
    P = P +   rct(113,k) * xnew(RA19AO2) * xnew(NO3)  &
    + rct(113,k) * xnew(DHPR12O2) * xnew(NO3)  &
    + 0.5* rct(113,k) * xnew(NRU12O2) * xnew(NO3)  &
    + 0.415* 4.2e-13 * xnew(RO2POOL) * xnew(RU12O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA13O2)  &
    + 0.7* 8.8e-13 * xnew(RO2POOL) * xnew(RA16O2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA19AO2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA19CO2)  &
    + 7.6e-13 * xnew(RO2POOL) * xnew(DHPR12O2)  &
    + 0.5* 9.6e-13 * xnew(RO2POOL) * xnew(NRU12O2)  &
    + 0.415* rcphot(IDCH3O2H,k) * xnew(RU12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA13OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA16OOH)  &
    + rcphot(IDCH3O2H,k) *3 * xnew(DHPR12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRU12OOH)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RA13NO3)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RA16NO3)  &
    + rcemis(GLYOX,k)
  
  L = rct(73,k) * xnew(OH)  &
    + 4.0e-16 * xnew(NO3)  &
    + rcphot(IDRCOCHO,k) *0.46  &
    + rcphot(IDRCOCHO,k) *0.07  &
    + rcphot(IDCH3O2H,k) *4.79
  
  xnew(GLYOX) = (xold(GLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> MGLYOX

  P = 0.902* rct(70,k) * xnew(O3) * xnew(MACR)  &
    + 0.62* rct(72,k) * xnew(O3) * xnew(MVK)  &
    + rct(75,k) * xnew(ACETOL) * xnew(OH)  &
    + 0.5* 2.4e-17 * xnew(O3) * xnew(UCARB12)  &
    + 0.158* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 5.78e-11 * xnew(OH) * xnew(HUCARB9)  &
    + 2.4e-17 * xnew(O3) * xnew(HPUCARB12)  &
    + 3.42e-11 * xnew(OH) * xnew(DHCARB9)  &
    + 1.2e-11 * xnew(OH) * xnew(RN8OOH)  &
    + 0.7* 4.33e-12 * xnew(OH) * xnew(HMML)  &
    + 6.7e-13 * xnew(OH) * xnew(NOA)  &
    + 8.31e-11 * xnew(OH) * xnew(RA19NO3)  &
    + 5.93e-13 * xnew(OH) * xnew(CO2C3PAN)  &
    + 0.295* rct(109,k) * xnew(RU10O2) * xnew(NO)  &
    + 0.407* rct(109,k) * xnew(RU12O2) * xnew(NO)  &
    + 0.267* rct(109,k) * xnew(RA16O2) * xnew(NO)  &
    + 0.3* rct(113,k) * xnew(RU10O2) * xnew(NO3)  &
    + 0.415* rct(113,k) * xnew(RU12O2) * xnew(NO3)  &
    + 0.3* rct(113,k) * xnew(RA16O2) * xnew(NO3)  &
    + 0.3* 1.22e-12 * xnew(RO2POOL) * xnew(RU10O2)
  
  !-> MGLYOX cont.
    P = P +   0.415* 4.2e-13 * xnew(RO2POOL) * xnew(RU12O2)  &
    + 0.3* 8.8e-13 * xnew(RO2POOL) * xnew(RA16O2)  &
    + rcphot(IDHCHO_H2,k) *6 * xnew(HUCARB9)  &
    + rcphot(IDCH3O2H,k) *2 * xnew(DHPCARB9)  &
    + 0.226* rcphot(IDCH3O2H,k) * xnew(RU10OOH)  &
    + 0.415* rcphot(IDCH3O2H,k) * xnew(RU12OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA19OOH)  &
    + rcphot(IDMEK,k) *83. * xnew(RU12NO3)  &
    + rcphot(IDHCHO_H,k) *0.234 * xnew(RA19NO3)  &
    + rcemis(MGLYOX,k)
  
  L = rct(74,k) * xnew(OH)  &
    + 5.0e-16 * xnew(NO3)  &
    + rcphot(IDRCOCHO,k)
  
  xnew(MGLYOX) = (xold(MGLYOX) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OH

  P = 0.5* rct(131,k) * xnew(RO2POOL) * xnew(CH3O2)  &
    + rcemis(CH3OH,k)
  
  L = rct(33,k) * xnew(OH)
  
  xnew(CH3OH) = (xold(CH3OH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OH

  P = 0.2* 3.1e-13 * xnew(RO2POOL) * xnew(C2H5O2)  &
    + rcemis(C2H5OH,k)
  
  L = rct(34,k) * xnew(OH)
  
  xnew(C2H5OH) = (xold(C2H5OH) + dt2 * P) / (1.0 + dt2 * L)


!-> NPROPOL

  P = 0.2* 6.0e-13 * xnew(RO2POOL) * xnew(RN10O2)  &
    + rcemis(NPROPOL,k)
  
  L = rct(35,k) * xnew(OH)
  
  xnew(NPROPOL) = (xold(NPROPOL) + dt2 * P) / (1.0 + dt2 * L)


!-> IPROPOL

  P = 0.2* 4.0e-14 * xnew(RO2POOL) * xnew(IC3H7O2)  &
    + rcemis(IPROPOL,k)
  
  L = rct(36,k) * xnew(OH)
  
  xnew(IPROPOL) = (xold(IPROPOL) + dt2 * P) / (1.0 + dt2 * L)


!-> ACETOL

  P = 0.25* 2.4e-17 * xnew(O3) * xnew(UCARB12)  &
    + 2.5e-11 * xnew(OH) * xnew(RN9OOH)  &
    + 3.84e-11 * xnew(OH) * xnew(RU10OOH)  &
    + 1.31e-12 * xnew(OH) * xnew(RN9NO3)  &
    + 5.26e-13 * xnew(OH) * xnew(RU10NO3)  &
    + 2.5e-12 * xnew(OH) * xnew(RU12NO3)  &
    + 0.22* 2.9e-11 * xnew(OH) * xnew(MPAN)  &
    + rct(99,k) * xnew(MACRO2)  &
    + 0.987* rct(109,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.558* rct(109,k) * xnew(RU12O2) * xnew(NO)  &
    + rct(113,k) * xnew(MACRO2) * xnew(NO3)  &
    + 0.585* rct(113,k) * xnew(RU12O2) * xnew(NO3)  &
    + 3.6e-13 * xnew(RO2POOL) * xnew(MACRO2)  &
    + 0.585* 4.2e-13 * xnew(RO2POOL) * xnew(RU12O2)  &
    + 0.25* rcphot(IDHCHO_H2,k) *0.0474 * xnew(UCARB12)  &
    + 0.425* rcphot(IDHCHO_H2,k) *12 * xnew(HPUCARB12)  &
    + rcphot(IDHCHO_H,k) *1.67 * xnew(DHCARB9)  &
    + 0.5* rcphot(IDMEK,k) *132.76 * xnew(NUCARB12)  &
    + 0.157* rcphot(IDCH3O2H,k) * xnew(RU10OOH)  &
    + 0.585* rcphot(IDCH3O2H,k) * xnew(RU12OOH)
  
  !-> ACETOL cont.
    P = P +   rcemis(ACETOL,k)
  
  L = rct(75,k) * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(ACETOL) = (xold(ACETOL) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3NO3

  P = 0.001* rct(105,k) * xnew(CH3O2) * xnew(NO)
  
  L = rct(80,k) * xnew(OH)  &
    + rcphot(IDMEK,k) *0.293
  
  xnew(CH3NO3) = (xold(CH3NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5NO3

  P = 0.033* rct(106,k) * xnew(C2H5O2) * xnew(NO)
  
  L = rct(81,k) * xnew(OH)  &
    + rcphot(IDHNO3,k) *1.932
  
  xnew(C2H5NO3) = (xold(C2H5NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RN10NO3

  P = 0.02* rct(107,k) * xnew(RN10O2) * xnew(NO)
  
  L = 5.8e-13 * xnew(OH)  &
    + rcphot(IDHNO3,k) *2.62
  
  xnew(RN10NO3) = (xold(RN10NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> IC3H7NO3

  P = 0.035* rct(108,k) * xnew(IC3H7O2) * xnew(NO)
  
  L = rct(82,k) * xnew(OH)  &
    + rcphot(IDMEK,k) *0.778
  
  xnew(IC3H7NO3) = (xold(IC3H7NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RN13NO3

  P = 0.083* rct(109,k) * xnew(RN13O2) * xnew(NO)
  
  L = 9.2e-13 * xnew(OH)  &
    + rcphot(IDMEK,k) *0.778
  
  xnew(RN13NO3) = (xold(RN13NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RN16NO3

  P = 0.123* rct(109,k) * xnew(RN16O2) * xnew(NO)
  
  L = 1.85e-12 * xnew(OH)  &
    + rcphot(IDHNO3,k) *2.62
  
  xnew(RN16NO3) = (xold(RN16NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RN19NO3

  P = 0.278* rct(109,k) * xnew(RN19O2) * xnew(NO)
  
  L = 9.6e-11 * xnew(OH)  &
    + rcphot(IDMEK,k) *0.778
  
  xnew(RN19NO3) = (xold(RN19NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> HOC2H4NO3

  P = 0.005* rct(109,k) * xnew(HOCH2CH2O2) * xnew(NO)
  
  L = 1.09e-12 * xnew(OH)
  
  xnew(HOC2H4NO3) = (xold(HOC2H4NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RN9NO3

  P = 0.021* rct(109,k) * xnew(RN9O2) * xnew(NO)
  
  L = 1.31e-12 * xnew(OH)
  
  xnew(RN9NO3) = (xold(RN9NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RN12NO3

  P = 0.041* rct(109,k) * xnew(RN12O2) * xnew(NO)
  
  L = 1.79e-12 * xnew(OH)
  
  xnew(RN12NO3) = (xold(RN12NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RN15NO3

  P = 0.064* rct(109,k) * xnew(RN15O2) * xnew(NO)  &
    + 0.025* rct(109,k) * xnew(RN15AO2) * xnew(NO)
  
  L = 1.03e-11 * xnew(OH)
  
  xnew(RN15NO3) = (xold(RN15NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3OOH

  P = 0.9* rct(116,k) * xnew(CH3O2) * xnew(HO2)
  
  L = rct(79,k) * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + 1.0e-5
  
  xnew(CH3OOH) = (xold(CH3OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5OOH

  P = rct(117,k) * xnew(C2H5O2) * xnew(HO2)
  
  L = 1.36e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(C2H5OOH) = (xold(C2H5OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN10OOH

  P = rct(121,k) * xnew(RN10O2) * xnew(HO2)
  
  L = 1.89e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN10OOH) = (xold(RN10OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> IC3H7OOH

  P = rct(118,k) * xnew(IC3H7O2) * xnew(HO2)
  
  L = 2.78e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(IC3H7OOH) = (xold(IC3H7OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN13OOH

  P = rct(120,k) * xnew(RN13O2) * xnew(HO2)  &
    + rct(120,k) * xnew(RN13AO2) * xnew(HO2)
  
  L = 3.57e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN13OOH) = (xold(RN13OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN16OOH

  P = rct(122,k) * xnew(RN16O2) * xnew(HO2)  &
    + rct(122,k) * xnew(RN16AO2) * xnew(HO2)
  
  L = 4.21e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN16OOH) = (xold(RN16OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN19OOH

  P = rct(123,k) * xnew(RN19O2) * xnew(HO2)
  
  L = 1.27e-10 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN19OOH) = (xold(RN19OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CO3H

  P = 0.384* rct(129,k) * xnew(CH3CO3) * xnew(HO2)
  
  L = 3.7e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(CH3CO3H) = (xold(CH3CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> RN8OOH

  P = 3.64e-11 * xnew(OH) * xnew(DHPCARB9)  &
    + 6.71e-11 * xnew(OH) * xnew(RU13OOH)  &
    + rct(109,k) * xnew(DHPR12O2) * xnew(NO)  &
    + rct(113,k) * xnew(DHPR12O2) * xnew(NO3)  &
    + rct(121,k) * xnew(RN8O2) * xnew(HO2)  &
    + 7.6e-13 * xnew(RO2POOL) * xnew(DHPR12O2)  &
    + rcphot(IDHCHO_H,k) *0.6 * xnew(DHPCARB9)  &
    + rcphot(IDCH3O2H,k) *3 * xnew(DHPR12OOH)
  
  L = 1.2e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDMEK,k)
  
  xnew(RN8OOH) = (xold(RN8OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN11OOH

  P = rct(120,k) * xnew(RN11O2) * xnew(HO2)
  
  L = 2.5e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN11OOH) = (xold(RN11OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN14OOH

  P = rct(122,k) * xnew(RN14O2) * xnew(HO2)
  
  L = 3.2e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN14OOH) = (xold(RN14OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN17OOH

  P = rct(123,k) * xnew(RN17O2) * xnew(HO2)
  
  L = 7.06e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN17OOH) = (xold(RN17OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> HOC2H4OOH

  P = rct(119,k) * xnew(HOCH2CH2O2) * xnew(HO2)
  
  L = 2.13e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(HOC2H4OOH) = (xold(HOC2H4OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN9OOH

  P = rct(121,k) * xnew(RN9O2) * xnew(HO2)
  
  L = 2.5e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN9OOH) = (xold(RN9OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN12OOH

  P = rct(120,k) * xnew(RN12O2) * xnew(HO2)
  
  L = 3.25e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN12OOH) = (xold(RN12OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RN15OOH

  P = rct(122,k) * xnew(RN15O2) * xnew(HO2)  &
    + rct(122,k) * xnew(RN15AO2) * xnew(HO2)
  
  L = 3.74e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RN15OOH) = (xold(RN15OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> NRN6OOH

  P = rct(128,k) * xnew(NRN6O2) * xnew(HO2)
  
  L = 5.22e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(NRN6OOH) = (xold(NRN6OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> NRN9OOH

  P = rct(121,k) * xnew(NRN9O2) * xnew(HO2)
  
  L = 6.5e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(NRN9OOH) = (xold(NRN9OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> NRN12OOH

  P = rct(120,k) * xnew(NRN12O2) * xnew(HO2)
  
  L = 7.15e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(NRN12OOH) = (xold(NRN12OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RU13OOH

  P = rct(122,k) * xnew(RU13O2) * xnew(HO2)
  
  L = 6.71e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDMEK,k)
  
  xnew(RU13OOH) = (xold(RU13OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> PAN

  P = rct(100,k) * xnew(CH3CO3) * xnew(NO2)
  
  L = 3.0e-14 * xnew(OH)  &
    + rct(101,k)
  
  xnew(PAN) = (xold(PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> APINENE

  P = rcemis(APINENE,k)  &
    + RCBIO(2,K)
  
  L = rct(63,k) * xnew(OH)  &
    + rct(64,k) * xnew(NO3)  &
    + rct(65,k) * xnew(O3)
  
  xnew(APINENE) = (xold(APINENE) + dt2 * P) / (1.0 + dt2 * L)


!-> BPINENE

  P = rcemis(BPINENE,k)
  
  L = rct(66,k) * xnew(OH)  &
    + 2.51e-12 * xnew(NO3)  &
    + rct(67,k) * xnew(O3)
  
  xnew(BPINENE) = (xold(BPINENE) + dt2 * P) / (1.0 + dt2 * L)


!-> CARB10

  P = 3.25e-11 * xnew(OH) * xnew(RN12OOH)  &
    + 1.79e-12 * xnew(OH) * xnew(RN12NO3)
  
  L = 5.86e-12 * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(CARB10) = (xold(CARB10) + dt2 * P) / (1.0 + dt2 * L)


!-> CARB13

  P = 9.0e-13 * xnew(OH) * xnew(ARNOH14)  &
    + 9.0e-14 * xnew(NO3) * xnew(ARNOH14)  &
    + 3.74e-11 * xnew(OH) * xnew(RN15OOH)  &
    + 1.03e-11 * xnew(OH) * xnew(RN15NO3)  &
    + 0.975* rct(109,k) * xnew(RN15AO2) * xnew(NO)  &
    + rct(113,k) * xnew(RN15AO2) * xnew(NO3)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RN15AO2)
  
  L = 1.65e-11 * xnew(OH)  &
    + rcphot(IDMEK,k) *3.00
  
  xnew(CARB13) = (xold(CARB13) + dt2 * P) / (1.0 + dt2 * L)


!-> CARB16

  P = 1.53e-12 * xnew(OH) * xnew(ARNOH17)  &
    + 3.13e-13 * xnew(NO3) * xnew(ARNOH17)  &
    + 3.66e-12 * xnew(OH) * xnew(RTN26PAN)
  
  L = 1.25e-11 * xnew(OH)  &
    + rcphot(IDMEK,k) *3.35
  
  xnew(CARB16) = (xold(CARB16) + dt2 * P) / (1.0 + dt2 * L)


!-> CARB9

  P = 5.86e-12 * xnew(OH) * xnew(CARB10)  &
    + 2.5e-11 * xnew(OH) * xnew(RN11OOH)  &
    + 0.862* rct(109,k) * xnew(RA19CO2) * xnew(NO)  &
    + rct(113,k) * xnew(RA19CO2) * xnew(NO3)
  
  L = 2.4e-13 * xnew(OH)  &
    + rcphot(IDCH3COY,k)
  
  xnew(CARB9) = (xold(CARB9) + dt2 * P) / (1.0 + dt2 * L)


!-> CARB12

  P = 3.2e-11 * xnew(OH) * xnew(RN14OOH)
  
  L = 1.38e-12 * xnew(OH)  &
    + rcphot(IDCH3COY,k)
  
  xnew(CARB12) = (xold(CARB12) + dt2 * P) / (1.0 + dt2 * L)


!-> C2H5CO3H

  P = 0.56* rct(129,k) * xnew(C2H5CO3) * xnew(HO2)
  
  L = 4.42e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(C2H5CO3H) = (xold(C2H5CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> PPN

  P = rct(100,k) * xnew(C2H5CO3) * xnew(NO2)
  
  L = 1.27e-12 * xnew(OH)  &
    + rct(102,k)
  
  xnew(PPN) = (xold(PPN) + dt2 * P) / (1.0 + dt2 * L)


!-> HOCH2CHO

  P = 0.25* 2.4e-17 * xnew(O3) * xnew(UCARB12)  &
    + 2.13e-11 * xnew(OH) * xnew(HOC2H4OOH)  &
    + 1.09e-12 * xnew(OH) * xnew(HOC2H4NO3)  &
    + 0.223* rct(109,k) * xnew(HOCH2CH2O2) * xnew(NO)  &
    + 0.67* rct(109,k) * xnew(RU10O2) * xnew(NO)  &
    + 0.407* rct(109,k) * xnew(RU12O2) * xnew(NO)  &
    + 0.224* rct(113,k) * xnew(HOCH2CH2O2) * xnew(NO3)  &
    + 0.7* rct(113,k) * xnew(RU10O2) * xnew(NO3)  &
    + 0.415* rct(113,k) * xnew(RU12O2) * xnew(NO3)  &
    + 0.224* 2.0e-12 * xnew(RO2POOL) * xnew(HOCH2CH2O2)  &
    + 0.7* 1.22e-12 * xnew(RO2POOL) * xnew(RU10O2)  &
    + 0.415* 4.2e-13 * xnew(RO2POOL) * xnew(RU12O2)  &
    + 0.25* rcphot(IDHCHO_H2,k) *0.0474 * xnew(UCARB12)  &
    + 0.075* rcphot(IDHCHO_H2,k) *12 * xnew(HPUCARB12)  &
    + 0.467* rcphot(IDCH3O2H,k) * xnew(RU10OOH)  &
    + 0.415* rcphot(IDCH3O2H,k) * xnew(RU12OOH)  &
    + rcphot(IDMEK,k) *13.3 * xnew(RU10NO3)  &
    + rcphot(IDMEK,k) *83. * xnew(RU12NO3)  &
    + rcemis(HOCH2CHO,k)
  
  L = 8.0e-12 * xnew(OH)  &
    + rct(40,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k) *0.6
  
  xnew(HOCH2CHO) = (xold(HOCH2CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> HOCH2CO3H

  P = 0.56* rct(129,k) * xnew(HOCH2CO3) * xnew(HO2)
  
  L = 6.19e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(HOCH2CO3H) = (xold(HOCH2CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> PHAN

  P = rct(100,k) * xnew(HOCH2CO3) * xnew(NO2)
  
  L = 1.12e-12 * xnew(OH)  &
    + rct(101,k)
  
  xnew(PHAN) = (xold(PHAN) + dt2 * P) / (1.0 + dt2 * L)


!-> CCARB12

  P = 7.5e-11 * xnew(OH) * xnew(RTX22OOH)  &
    + 4.74e-12 * xnew(OH) * xnew(RTX22NO3)
  
  L = 4.79e-12 * xnew(OH)
  
  xnew(CCARB12) = (xold(CCARB12) + dt2 * P) / (1.0 + dt2 * L)


!-> UCARB12

  P = 0.09* 7.14e-11 * xnew(OH) * xnew(RU14OOH)  &
    + 0.34* 3.0e-11 * xnew(OH) * xnew(RU14NO3)  &
    + 0.908* rct(109,k) * xnew(CRU14O2) * xnew(NO)  &
    + 0.908* rct(109,k) * xnew(TRU14O2) * xnew(NO)  &
    + rct(113,k) * xnew(CRU14O2) * xnew(NO3)  &
    + rct(113,k) * xnew(TRU14O2) * xnew(NO3)  &
    + 2.4e-12 * xnew(RO2POOL) * xnew(CRU14O2)  &
    + 2.4e-12 * xnew(RO2POOL) * xnew(TRU14O2)
  
  L = 6.42e-11 * xnew(OH)  &
    + rct(77,k) * xnew(NO3)  &
    + 2.4e-17 * xnew(O3)  &
    + rcphot(IDHCHO_H2,k) *0.0474
  
  xnew(UCARB12) = (xold(UCARB12) + dt2 * P) / (1.0 + dt2 * L)


!-> UCARB13

  P = 0.062* rct(58,k) * xnew(OH) * xnew(C5H8)
  
  L = 4.19e-11 * xnew(OH)  &
    + 1.19e-17 * xnew(O3)  &
    + rcphot(IDMEK,k)
  
  xnew(UCARB13) = (xold(UCARB13) + dt2 * P) / (1.0 + dt2 * L)


!-> RU14NO3

  P = 0.092* rct(109,k) * xnew(CRU14O2) * xnew(NO)  &
    + 0.092* rct(109,k) * xnew(TRU14O2) * xnew(NO)  &
    + 0.111* rct(109,k) * xnew(BRU14O2) * xnew(NO)
  
  L = 3.0e-11 * xnew(OH)  &
    + rcphot(IDHNO3,k) *2.62
  
  xnew(RU14NO3) = (xold(RU14NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RU14OOH

  P = rct(122,k) * xnew(CRU14O2) * xnew(HO2)  &
    + rct(122,k) * xnew(TRU14O2) * xnew(HO2)  &
    + rct(122,k) * xnew(BRU14O2) * xnew(HO2)
  
  L = 7.14e-11 * xnew(OH)
  
  xnew(RU14OOH) = (xold(RU14OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RU12OOH

  P = rct(122,k) * xnew(RU12O2) * xnew(HO2)
  
  L = 3.5e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RU12OOH) = (xold(RU12OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RU10OOH

  P = 3.5e-11 * xnew(OH) * xnew(RU12OOH)  &
    + rct(120,k) * xnew(MACRO2) * xnew(HO2)  &
    + rct(120,k) * xnew(RU10O2) * xnew(HO2)  &
    + 0.56* rct(129,k) * xnew(MACO3) * xnew(HO2)
  
  L = 3.84e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RU10OOH) = (xold(RU10OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> MPAN

  P = rct(100,k) * xnew(MACO3) * xnew(NO2)
  
  L = 2.9e-11 * xnew(OH)  &
    + rct(103,k)
  
  xnew(MPAN) = (xold(MPAN) + dt2 * P) / (1.0 + dt2 * L)


!-> NUCARB12

  P = 1.03e-10 * xnew(OH) * xnew(NRU14OOH)  &
    + rct(109,k) * xnew(NRU14O2) * xnew(NO)  &
    + rct(113,k) * xnew(NRU14O2) * xnew(NO3)  &
    + 1.3e-12 * xnew(RO2POOL) * xnew(NRU14O2)  &
    + rcphot(IDCH3O2H,k) * xnew(NRU14OOH)
  
  L = 4.16e-11 * xnew(OH)  &
    + rcphot(IDMEK,k) *132.76
  
  xnew(NUCARB12) = (xold(NUCARB12) + dt2 * P) / (1.0 + dt2 * L)


!-> NRU14OOH

  P = rct(122,k) * xnew(NRU14O2) * xnew(HO2)
  
  L = 1.03e-10 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(NRU14OOH) = (xold(NRU14OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> NRU12OOH

  P = rct(122,k) * xnew(NRU12O2) * xnew(HO2)
  
  L = 2.65e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(NRU12OOH) = (xold(NRU12OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> NOA

  P = 2.65e-11 * xnew(OH) * xnew(NRU12OOH)  &
    + 1.6e-12 * xnew(OH) * xnew(RU13NO3)  &
    + rct(109,k) * xnew(NRU12O2) * xnew(NO)  &
    + rct(113,k) * xnew(NRU12O2) * xnew(NO3)  &
    + 9.6e-13 * xnew(RO2POOL) * xnew(NRU12O2)  &
    + rcphot(IDCH3O2H,k) * xnew(NRU12OOH)
  
  L = 6.7e-13 * xnew(OH)  &
    + rcphot(IDMEK,k) *8.3
  
  xnew(NOA) = (xold(NOA) + dt2 * P) / (1.0 + dt2 * L)


!-> RA13NO3

  P = 0.082* rct(109,k) * xnew(RA13O2) * xnew(NO)
  
  L = 7.3e-11 * xnew(OH)  &
    + rcphot(IDMEK,k) *0.778
  
  xnew(RA13NO3) = (xold(RA13NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RA13OOH

  P = rct(123,k) * xnew(RA13O2) * xnew(HO2)
  
  L = 9.77e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RA13OOH) = (xold(RA13OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> UDCARB8

  P = 9.77e-11 * xnew(OH) * xnew(RA13OOH)  &
    + 7.3e-11 * xnew(OH) * xnew(RA13NO3)  &
    + 0.918* rct(109,k) * xnew(RA13O2) * xnew(NO)  &
    + 0.267* rct(109,k) * xnew(RA16O2) * xnew(NO)  &
    + 0.862* rct(109,k) * xnew(RA19CO2) * xnew(NO)  &
    + rct(113,k) * xnew(RA13O2) * xnew(NO3)  &
    + 0.3* rct(113,k) * xnew(RA16O2) * xnew(NO3)  &
    + rct(113,k) * xnew(RA19CO2) * xnew(NO3)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA13O2)  &
    + 0.3* 8.8e-13 * xnew(RO2POOL) * xnew(RA16O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RA13OOH)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RA13NO3)
  
  L = 5.2e-11 * xnew(OH)  &
    + (rcphot(IDNO2,k) *(1.28e-02))  &
    + (rcphot(IDNO2,k) *(7.20e-03))
  
  xnew(UDCARB8) = (xold(UDCARB8) + dt2 * P) / (1.0 + dt2 * L)


!-> AROH14

  P = 0.53* rct(56,k) * xnew(OH) * xnew(BENZENE)
  
  L = 2.63e-11 * xnew(OH)  &
    + 3.78e-12 * xnew(NO3)
  
  xnew(AROH14) = (xold(AROH14) + dt2 * P) / (1.0 + dt2 * L)


!-> RAROH14

  P = 2.63e-11 * xnew(OH) * xnew(AROH14)  &
    + 3.78e-12 * xnew(NO3) * xnew(AROH14)
  
  L = 2.08e-12 * xnew(NO2)
  
  xnew(RAROH14) = (xold(RAROH14) + dt2 * P) / (1.0 + dt2 * L)


!-> ARNOH14

  P = 2.08e-12 * xnew(RAROH14) * xnew(NO2)
  
  L = 9.0e-13 * xnew(OH)  &
    + 9.0e-14 * xnew(NO3)
  
  xnew(ARNOH14) = (xold(ARNOH14) + dt2 * P) / (1.0 + dt2 * L)


!-> RA16NO3

  P = 0.111* rct(109,k) * xnew(RA16O2) * xnew(NO)
  
  L = 7.16e-11 * xnew(OH)  &
    + rcphot(IDMEK,k) *0.778
  
  xnew(RA16NO3) = (xold(RA16NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RA16OOH

  P = rct(124,k) * xnew(RA16O2) * xnew(HO2)
  
  L = 9.64e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RA16OOH) = (xold(RA16OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> UDCARB11

  P = 9.64e-11 * xnew(OH) * xnew(RA16OOH)  &
    + 7.16e-11 * xnew(OH) * xnew(RA16NO3)  &
    + 8.31e-11 * xnew(OH) * xnew(RA19NO3)  &
    + 0.622* rct(109,k) * xnew(RA16O2) * xnew(NO)  &
    + 0.7* rct(113,k) * xnew(RA16O2) * xnew(NO3)  &
    + 0.7* 8.8e-13 * xnew(RO2POOL) * xnew(RA16O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RA16OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(RA19OOH)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RA16NO3)  &
    + rcphot(IDHCHO_H,k) *0.234 * xnew(RA19NO3)
  
  L = 5.58e-11 * xnew(OH)  &
    + (rcphot(IDNO2,k) *(9.00e-03))  &
    + (rcphot(IDNO2,k) *(1.10e-02))
  
  xnew(UDCARB11) = (xold(UDCARB11) + dt2 * P) / (1.0 + dt2 * L)


!-> AROH17

  P = 0.18* rct(57,k) * xnew(OH) * xnew(TOLUENE)
  
  L = 4.65e-11 * xnew(OH)  &
    + 1.25e-11 * xnew(NO3)
  
  xnew(AROH17) = (xold(AROH17) + dt2 * P) / (1.0 + dt2 * L)


!-> RAROH17

  P = 4.65e-11 * xnew(OH) * xnew(AROH17)  &
    + 1.25e-11 * xnew(NO3) * xnew(AROH17)
  
  L = 2.08e-12 * xnew(NO2)
  
  xnew(RAROH17) = (xold(RAROH17) + dt2 * P) / (1.0 + dt2 * L)


!-> ARNOH17

  P = 2.08e-12 * xnew(RAROH17) * xnew(NO2)
  
  L = 1.53e-12 * xnew(OH)  &
    + 3.13e-13 * xnew(NO3)
  
  xnew(ARNOH17) = (xold(ARNOH17) + dt2 * P) / (1.0 + dt2 * L)


!-> UDCARB14

  P = 0.862* rct(109,k) * xnew(RA19AO2) * xnew(NO)  &
    + rct(113,k) * xnew(RA19AO2) * xnew(NO3)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA19AO2)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RA19CO2)
  
  L = 7.0e-11 * xnew(OH)  &
    + (rcphot(IDNO2,k) *(9.00e-03))  &
    + (rcphot(IDNO2,k) *(1.10e-02))
  
  xnew(UDCARB14) = (xold(UDCARB14) + dt2 * P) / (1.0 + dt2 * L)


!-> RA19NO3

  P = 0.138* rct(109,k) * xnew(RA19AO2) * xnew(NO)  &
    + 0.138* rct(109,k) * xnew(RA19CO2) * xnew(NO)
  
  L = 8.31e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *0.234
  
  xnew(RA19NO3) = (xold(RA19NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RA19OOH

  P = rct(125,k) * xnew(RA19AO2) * xnew(HO2)  &
    + rct(125,k) * xnew(RA19CO2) * xnew(HO2)
  
  L = 7.62e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RA19OOH) = (xold(RA19OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN28NO3

  P = 0.233* rct(109,k) * xnew(RTN28O2) * xnew(NO)
  
  L = 4.35e-12 * xnew(OH)
  
  xnew(RTN28NO3) = (xold(RTN28NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN28OOH

  P = rct(127,k) * xnew(RTN28O2) * xnew(HO2)
  
  L = 2.38e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTN28OOH) = (xold(RTN28OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> TNCARB26

  P = 0.075* rct(65,k) * xnew(APINENE) * xnew(O3)  &
    + 2.38e-11 * xnew(OH) * xnew(RTN28OOH)  &
    + 9.5e-12 * xnew(OH) * xnew(NRTN28OOH)  &
    + 4.35e-12 * xnew(OH) * xnew(RTN28NO3)  &
    + 0.702* rct(109,k) * xnew(RTN28O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRTN28O2) * xnew(NO)  &
    + rct(113,k) * xnew(RTN28O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRTN28O2) * xnew(NO3)  &
    + 2.85e-13 * xnew(RO2POOL) * xnew(RTN28O2)  &
    + 1.0e-13 * xnew(RO2POOL) * xnew(NRTN28O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN28OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRTN28OOH)
  
  L = 4.2e-11 * xnew(OH)  &
    + 3.8e-14 * xnew(NO3)  &
    + rcphot(IDHCHO_H,k) *0.6
  
  xnew(TNCARB26) = (xold(TNCARB26) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN26OOH

  P = 0.56* rct(129,k) * xnew(RTN26O2) * xnew(HO2)
  
  L = 1.2e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTN26OOH) = (xold(RTN26OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> NRTN28OOH

  P = rct(127,k) * xnew(NRTN28O2) * xnew(HO2)
  
  L = 9.5e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(NRTN28OOH) = (xold(NRTN28OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN26PAN

  P = rct(104,k) * xnew(RTN26O2) * xnew(NO2)
  
  L = 3.66e-12 * xnew(OH)  &
    + rct(101,k)
  
  xnew(RTN26PAN) = (xold(RTN26PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN25OOH

  P = rct(126,k) * xnew(RTN25O2) * xnew(HO2)
  
  L = 1.66e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTN25OOH) = (xold(RTN25OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN24OOH

  P = rct(126,k) * xnew(RTN24O2) * xnew(HO2)
  
  L = 1.05e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTN24OOH) = (xold(RTN24OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN23OOH

  P = rct(126,k) * xnew(RTN23O2) * xnew(HO2)
  
  L = 2.05e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTN23OOH) = (xold(RTN23OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN14OOH

  P = rct(123,k) * xnew(RTN14O2) * xnew(HO2)
  
  L = 8.69e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTN14OOH) = (xold(RTN14OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN10OOH

  P = rct(122,k) * xnew(RTN10O2) * xnew(HO2)
  
  L = 4.23e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTN10OOH) = (xold(RTN10OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> TNCARB10

  P = rct(109,k) * xnew(RTN14O2) * xnew(NO)  &
    + rct(113,k) * xnew(RTN14O2) * xnew(NO3)  &
    + 8.8e-13 * xnew(RO2POOL) * xnew(RTN14O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RTN14OOH)
  
  L = 1.0e-10 * xnew(OH)  &
    + rct(78,k) * xnew(NO3)  &
    + rcphot(IDCH3COY,k) *0.5
  
  xnew(TNCARB10) = (xold(TNCARB10) + dt2 * P) / (1.0 + dt2 * L)


!-> RTN25NO3

  P = 0.16* rct(109,k) * xnew(RTN25O2) * xnew(NO)
  
  L = 2.88e-12 * xnew(OH)
  
  xnew(RTN25NO3) = (xold(RTN25NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> TNCARB15

  P = 2.88e-12 * xnew(OH) * xnew(RTN25NO3)
  
  L = 1.0e-12 * xnew(OH)
  
  xnew(TNCARB15) = (xold(TNCARB15) + dt2 * P) / (1.0 + dt2 * L)


!-> RCOOH25

  P = 0.125* rct(65,k) * xnew(APINENE) * xnew(O3)
  
  L = 6.65e-12 * xnew(OH)
  
  xnew(RCOOH25) = (xold(RCOOH25) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX28NO3

  P = 0.233* rct(109,k) * xnew(RTX28O2) * xnew(NO)
  
  L = 3.53e-12 * xnew(OH)
  
  xnew(RTX28NO3) = (xold(RTX28NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX28OOH

  P = rct(127,k) * xnew(RTX28O2) * xnew(HO2)
  
  L = 2.0e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTX28OOH) = (xold(RTX28OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> TXCARB24

  P = 0.4* rct(67,k) * xnew(BPINENE) * xnew(O3)  &
    + 3.53e-12 * xnew(OH) * xnew(RTX28NO3)  &
    + 0.702* rct(109,k) * xnew(RTX28O2) * xnew(NO)  &
    + rct(109,k) * xnew(NRTX28O2) * xnew(NO)  &
    + rct(113,k) * xnew(RTX28O2) * xnew(NO3)  &
    + rct(113,k) * xnew(NRTX28O2) * xnew(NO3)  &
    + 2.0e-12 * xnew(RO2POOL) * xnew(RTX28O2)  &
    + 9.2e-14 * xnew(RO2POOL) * xnew(NRTX28O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX28OOH)  &
    + rcphot(IDCH3O2H,k) * xnew(NRTX28OOH)
  
  L = 1.55e-11 * xnew(OH)
  
  xnew(TXCARB24) = (xold(TXCARB24) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX24NO3

  P = 0.157* rct(109,k) * xnew(RTX24O2) * xnew(NO)
  
  L = 6.48e-12 * xnew(OH)  &
    + rcphot(IDMEK,k) *0.778
  
  xnew(RTX24NO3) = (xold(RTX24NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX24OOH

  P = rct(126,k) * xnew(RTX24O2) * xnew(HO2)
  
  L = 8.59e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTX24OOH) = (xold(RTX24OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> TXCARB22

  P = 0.25* rct(67,k) * xnew(BPINENE) * xnew(O3)  &
    + 8.59e-11 * xnew(OH) * xnew(RTX24OOH)  &
    + 6.48e-12 * xnew(OH) * xnew(RTX24NO3)  &
    + 0.506* rct(109,k) * xnew(RTX24O2) * xnew(NO)  &
    + rct(113,k) * xnew(RTX24O2) * xnew(NO3)  &
    + 2.5e-13 * xnew(RO2POOL) * xnew(RTX24O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RTX24OOH)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RTX24NO3)
  
  L = 4.55e-12 * xnew(OH)
  
  xnew(TXCARB22) = (xold(TXCARB22) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX22NO3

  P = 0.3* rct(109,k) * xnew(RTX22O2) * xnew(NO)
  
  L = 4.74e-12 * xnew(OH)
  
  xnew(RTX22NO3) = (xold(RTX22NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RTX22OOH

  P = rct(126,k) * xnew(RTX22O2) * xnew(HO2)
  
  L = 7.5e-11 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(RTX22OOH) = (xold(RTX22OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> NRTX28OOH

  P = rct(127,k) * xnew(NRTX28O2) * xnew(HO2)
  
  L = 9.58e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)
  
  xnew(NRTX28OOH) = (xold(NRTX28OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> ANHY

  P = 0.5* 5.2e-11 * xnew(OH) * xnew(UDCARB8)  &
    + 0.45* 5.58e-11 * xnew(OH) * xnew(UDCARB11)  &
    + 0.45* 7.0e-11 * xnew(OH) * xnew(UDCARB14)  &
    + (rcphot(IDNO2,k) *(7.20e-03)) * xnew(UDCARB8)  &
    + (rcphot(IDNO2,k) *(9.00e-03)) * xnew(UDCARB11)  &
    + (rcphot(IDNO2,k) *(9.00e-03)) * xnew(UDCARB14)
  
  L = 1.5e-12 * xnew(OH)
  
  xnew(ANHY) = (xold(ANHY) + dt2 * P) / (1.0 + dt2 * L)


!-> SO2

  P = rcemis(SO2,k)
  
  L = AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH)  &
    + AQRCK(ICLRC1,K) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(O3)  &
    + AQRCK(ICLRC3,K)
  
  xnew(SO2) = (xold(SO2) + dt2 * P) / (1.0 + dt2 * L)


!-> HPUCARB12

  P = 0.5* rct(91,k) * xnew(CRU14O2)
  
  L = 5.2e-11 * xnew(OH)  &
    + rct(77,k) * xnew(NO3)  &
    + 2.4e-17 * xnew(O3)  &
    + rcphot(IDHCHO_H2,k) *12
  
  xnew(HPUCARB12) = (xold(HPUCARB12) + dt2 * P) / (1.0 + dt2 * L)


!-> HUCARB9

  P = 5.2e-11 * xnew(OH) * xnew(HPUCARB12)  &
    + rct(77,k) * xnew(NO3) * xnew(HPUCARB12)  &
    + 0.5* rcphot(IDHCHO_H2,k) *12 * xnew(HPUCARB12)  &
    + 0.5* rcphot(IDMEK,k) *132.76 * xnew(NUCARB12)
  
  L = 5.78e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H2,k) *6
  
  xnew(HUCARB9) = (xold(HUCARB9) + dt2 * P) / (1.0 + dt2 * L)


!-> IEPOX

  P = 0.85* 7.14e-11 * xnew(OH) * xnew(RU14OOH)
  
  L = 1.16e-11 * xnew(OH)
  
  xnew(IEPOX) = (xold(IEPOX) + dt2 * P) / (1.0 + dt2 * L)


!-> HMML

  P = 0.78* 2.9e-11 * xnew(OH) * xnew(MPAN)
  
  L = 4.33e-12 * xnew(OH)
  
  xnew(HMML) = (xold(HMML) + dt2 * P) / (1.0 + dt2 * L)


!-> DHPR12OOH

  P = rct(122,k) * xnew(DHPR12O2) * xnew(HO2)
  
  L = 5.64e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *1.67  &
    + rcphot(IDCH3O2H,k) *3
  
  xnew(DHPR12OOH) = (xold(DHPR12OOH) + dt2 * P) / (1.0 + dt2 * L)


!-> DHPCARB9

  P = 5.64e-11 * xnew(OH) * xnew(DHPR12OOH)  &
    + rct(97,k) * xnew(DHPR12O2)  &
    + rcphot(IDHCHO_H,k) *1.67 * xnew(DHPR12OOH)
  
  L = 3.64e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *0.6  &
    + rcphot(IDCH3O2H,k) *2
  
  xnew(DHPCARB9) = (xold(DHPCARB9) + dt2 * P) / (1.0 + dt2 * L)


!-> DHCARB9

  P = rct(98,k) * xnew(RU12O2)
  
  L = 3.42e-11 * xnew(OH)  &
    + rcphot(IDHCHO_H,k) *1.67
  
  xnew(DHCARB9) = (xold(DHCARB9) + dt2 * P) / (1.0 + dt2 * L)


!-> CO2C3CHO

  P = 0.607* 1.19e-17 * xnew(O3) * xnew(UCARB13)  &
    + 0.949* rct(109,k) * xnew(RU13O2) * xnew(NO)  &
    + rct(113,k) * xnew(RU13O2) * xnew(NO3)  &
    + 6.0e-13 * xnew(RO2POOL) * xnew(RU13O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RU13OOH)  &
    + rcphot(IDMEK,k) * xnew(RU13OOH)  &
    + rcphot(IDMEK,k) * xnew(RU13NO3)
  
  L = 7.15e-11 * xnew(OH)  &
    + rct(76,k) * xnew(NO3)  &
    + rcphot(IDHCHO_H,k) *0.6
  
  xnew(CO2C3CHO) = (xold(CO2C3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> CO2C3CO3H

  P = 0.56* rct(129,k) * xnew(CO2C3CO3) * xnew(HO2)
  
  L = 4.18e-12 * xnew(OH)  &
    + rcphot(IDCH3O2H,k)  &
    + rcphot(IDMEK,k)
  
  xnew(CO2C3CO3H) = (xold(CO2C3CO3H) + dt2 * P) / (1.0 + dt2 * L)


!-> RU13NO3

  P = 0.051* rct(109,k) * xnew(RU13O2) * xnew(NO)
  
  L = 1.6e-12 * xnew(OH)  &
    + rcphot(IDMEK,k)
  
  xnew(RU13NO3) = (xold(RU13NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RU12NO3

  P = 0.18* 3.0e-11 * xnew(OH) * xnew(RU14NO3)  &
    + 0.035* rct(109,k) * xnew(RU12O2) * xnew(NO)
  
  L = 2.5e-12 * xnew(OH)  &
    + rcphot(IDMEK,k) *83.
  
  xnew(RU12NO3) = (xold(RU12NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> RU10NO3

  P = 0.48* 3.0e-11 * xnew(OH) * xnew(RU14NO3)  &
    + 0.013* rct(109,k) * xnew(MACRO2) * xnew(NO)  &
    + 0.035* rct(109,k) * xnew(RU10O2) * xnew(NO)
  
  L = 5.26e-13 * xnew(OH)  &
    + rcphot(IDMEK,k) *13.3
  
  xnew(RU10NO3) = (xold(RU10NO3) + dt2 * P) / (1.0 + dt2 * L)


!-> CO2C3PAN

  P = rct(100,k) * xnew(CO2C3CO3) * xnew(NO2)
  
  L = 5.93e-13 * xnew(OH)  &
    + rct(101,k)
  
  xnew(CO2C3PAN) = (xold(CO2C3PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> HCC7CO

  P = 1.27e-10 * xnew(OH) * xnew(RN19OOH)  &
    + 9.6e-11 * xnew(OH) * xnew(RN19NO3)  &
    + 0.722* rct(109,k) * xnew(RN19O2) * xnew(NO)  &
    + rct(113,k) * xnew(RN19O2) * xnew(NO3)  &
    + 2.5e-13 * xnew(RO2POOL) * xnew(RN19O2)  &
    + rcphot(IDCH3O2H,k) * xnew(RN19OOH)  &
    + rcphot(IDMEK,k) *0.778 * xnew(RN19NO3)
  
  L = 1.19e-10 * xnew(OH)
  
  xnew(HCC7CO) = (xold(HCC7CO) + dt2 * P) / (1.0 + dt2 * L)


!-> shipNOx

  P = rcemis(SHIPNOX,k)
  
  L = rct(26,k) * xnew(OH)  &
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


