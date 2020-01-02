!-> ASOC_ng100

  P = rct(76,k) * xnew(OH) * xnew(ASOC_ug1)
  
  ! L = 0.0
  
  xnew(ASOC_ng100) = xold(ASOC_ng100) + dt2 * P


!-> ASOC_ug1

  P = 0.00206* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.010294* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 0.008413* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 0.315476* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + rct(77,k) * xnew(OH) * xnew(ASOC_ug10)
  
  L = rct(76,k) * xnew(OH)
  
  xnew(ASOC_ug1) = (xold(ASOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug10

  P = 0.108* rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + 0.2132* rct(42,k) * xnew(SECC4H9O2) * xnew(HO2)  &
    + 0.01029* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.018529* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 0.820238* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 1.261905* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + rct(78,k) * xnew(OH) * xnew(ASOC_ug1e2)
  
  L = rct(77,k) * xnew(OH)
  
  xnew(ASOC_ug10) = (xold(ASOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e2

  P = 0.078235* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.123529* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 1.261905* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 1.577381* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + rct(79,k) * xnew(OH) * xnew(ASOC_ug1e3)
  
  L = rct(78,k) * xnew(OH)
  
  xnew(ASOC_ug1e2) = (xold(ASOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> ASOC_ug1e3

  P = 0.30882* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.463235* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 1.829762* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 2.208333* rct(47,k) * xnew(OXYO2) * xnew(HO2)
  
  L = rct(79,k) * xnew(OH)
  
  xnew(ASOC_ug1e3) = (xold(ASOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ng100

  P = 0.9* rct(76,k) * xnew(OH) * xnew(ASOC_ug1)  &
    + 1.075* rct(80,k) * xnew(OH) * xnew(NON_C_ASOA_ug1)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ng100) = xold(NON_C_ASOA_ng100) + dt2 * P


!-> NON_C_ASOA_ug1

  P = 0.0173* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.08647* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 0.11105* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 4.16429* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + 0.9* rct(77,k) * xnew(OH) * xnew(ASOC_ug10)  &
    + 1.075* rct(81,k) * xnew(OH) * xnew(NON_C_ASOA_ug10)
  
  L = rct(80,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1) = (xold(NON_C_ASOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug10

  P = 0.90753* rct(41,k) * xnew(NO) * xnew(SECC4H9O2)  &
    + 1.7912* rct(42,k) * xnew(SECC4H9O2) * xnew(HO2)  &
    + 0.08647* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 0.15565* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 10.8271* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 16.6571* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + 0.9* rct(78,k) * xnew(OH) * xnew(ASOC_ug1e2)  &
    + 1.075* rct(82,k) * xnew(OH) * xnew(NON_C_ASOA_ug1e2)
  
  L = rct(81,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug10) = (xold(NON_C_ASOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug1e2

  P = 0.65718* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 1.03765* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 16.6571* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 20.82143* rct(47,k) * xnew(OXYO2) * xnew(HO2)  &
    + 0.9* rct(79,k) * xnew(OH) * xnew(ASOC_ug1e3)  &
    + 1.075* rct(83,k) * xnew(OH) * xnew(NON_C_ASOA_ug1e3)
  
  L = rct(82,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1e2) = (xold(NON_C_ASOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_ASOA_ug1e3

  P = 2.59412* rct(41,k) * xnew(NO) * xnew(PRRO2)  &
    + 3.89118* rct(46,k) * xnew(PRRO2) * xnew(HO2)  &
    + 24.15286* rct(41,k) * xnew(OXYO2) * xnew(NO)  &
    + 29.15* rct(47,k) * xnew(OXYO2) * xnew(HO2)
  
  L = rct(83,k) * xnew(OH)
  
  xnew(NON_C_ASOA_ug1e3) = (xold(NON_C_ASOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ng100

  P = rct(84,k) * xnew(OH) * xnew(BSOC_ug1)
  
  ! L = 0.0
  
  xnew(BSOC_ng100) = xold(BSOC_ng100) + dt2 * P


!-> BSOC_ug1

  P = 0.002833* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.0255* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 0.08* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 0.715333* rct(75,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + rct(85,k) * xnew(OH) * xnew(BSOC_ug10)
  
  L = rct(84,k) * xnew(OH)
  
  xnew(BSOC_ug1) = (xold(BSOC_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug10

  P = 0.065167* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.085* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 0.813333* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 0.612* rct(75,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + rct(86,k) * xnew(OH) * xnew(BSOC_ug1e2)
  
  L = rct(85,k) * xnew(OH)
  
  xnew(BSOC_ug10) = (xold(BSOC_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e2

  P = 0.0425* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.0425* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 1.34* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 2.391333* rct(75,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + rct(87,k) * xnew(OH) * xnew(BSOC_ug1e3)
  
  L = rct(86,k) * xnew(OH)
  
  xnew(BSOC_ug1e2) = (xold(BSOC_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> BSOC_ug1e3

  P = 3.333333* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 4.05* rct(75,k) * xnew(HO2) * xnew(TERPPeroxy)
  
  L = rct(87,k) * xnew(OH)
  
  xnew(BSOC_ug1e3) = (xold(BSOC_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ng100

  P = 0.9* rct(84,k) * xnew(OH) * xnew(BSOC_ug1)  &
    + 1.075* rct(88,k) * xnew(OH) * xnew(NON_C_BSOA_ug1)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ng100) = xold(NON_C_BSOA_ng100) + dt2 * P


!-> NON_C_BSOA_ug1

  P = 0.034* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.306* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 0.672* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 6.009* rct(75,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + 0.9* rct(85,k) * xnew(OH) * xnew(BSOC_ug10)  &
    + 1.075* rct(89,k) * xnew(OH) * xnew(NON_C_BSOA_ug10)
  
  L = rct(88,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1) = (xold(NON_C_BSOA_ug1) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug10

  P = 0.782* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 1.02* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 6.832* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 5.1408* rct(75,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + 0.9* rct(86,k) * xnew(OH) * xnew(BSOC_ug1e2)  &
    + 1.075* rct(90,k) * xnew(OH) * xnew(NON_C_BSOA_ug1e2)
  
  L = rct(89,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug10) = (xold(NON_C_BSOA_ug10) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug1e2

  P = 0.51* rct(41,k) * xnew(ISRO2) * xnew(NO)  &
    + 0.51* rct(48,k) * xnew(ISRO2) * xnew(HO2)  &
    + 11.256* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 20.0872* rct(75,k) * xnew(HO2) * xnew(TERPPeroxy)  &
    + 0.9* rct(87,k) * xnew(OH) * xnew(BSOC_ug1e3)  &
    + 1.075* rct(91,k) * xnew(OH) * xnew(NON_C_BSOA_ug1e3)
  
  L = rct(90,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1e2) = (xold(NON_C_BSOA_ug1e2) + dt2 * P) / (1.0 + dt2 * L)


!-> NON_C_BSOA_ug1e3

  P = 28.* rct(41,k) * xnew(NO) * xnew(TERPPeroxy)  &
    + 34.02* rct(75,k) * xnew(HO2) * xnew(TERPPeroxy)
  
  L = rct(91,k) * xnew(OH)
  
  xnew(NON_C_BSOA_ug1e3) = (xold(NON_C_BSOA_ug1e3) + dt2 * P) / (1.0 + dt2 * L)


!-> Fffuel_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for Fffuel_ng10! xnew(Fffuel_ng10) = max(0.0, xold(Fffuel_ng10))


!-> woodOA_ng10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for woodOA_ng10! xnew(woodOA_ng10) = max(0.0, xold(woodOA_ng10))


!-> H2

  P = rcphot(IDHCHO_H2,k) * xnew(HCHO)
  
  L = rct(14,k) * xnew(OH)
  
  xnew(H2) = (xold(H2) + dt2 * P) / (1.0 + dt2 * L)


!-> CH4

  P = 0.1* rct(45,k) * xnew(O3) * xnew(C3H6)
  
  L = rct(19,k) * xnew(OH)
  
  xnew(CH4) = (xold(CH4) + dt2 * P) / (1.0 + dt2 * L)


!-> POM_f_wood

  P = rcemis(POM_F_WOOD,k)
  
  ! L = 0.0
  
  xnew(POM_f_wood) = xold(POM_f_wood) + dt2 * P


!-> POM_f_ffuel

  P = rcemis(POM_F_FFUEL,k)
  
  ! L = 0.0
  
  xnew(POM_f_ffuel) = xold(POM_f_ffuel) + dt2 * P


!-> POM_c_ffuel

  P = rcemis(POM_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(POM_c_ffuel) = xold(POM_c_ffuel) + dt2 * P


!-> EC_f_wood_new

  P = rcemis(EC_F_WOOD_NEW,k)
  
  L = rct(92,k)
  
  xnew(EC_f_wood_new) = (xold(EC_f_wood_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_wood_age

  P = rcemis(EC_F_WOOD_AGE,k)  &
    + rct(92,k) * xnew(EC_f_wood_new)
  
  ! L = 0.0
  
  xnew(EC_f_wood_age) = xold(EC_f_wood_age) + dt2 * P


!-> EC_c_wood

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for EC_c_wood! xnew(EC_c_wood) = max(0.0, xold(EC_c_wood))


!-> EC_f_ffuel_new

  P = rcemis(EC_F_FFUEL_NEW,k)
  
  L = rct(92,k)
  
  xnew(EC_f_ffuel_new) = (xold(EC_f_ffuel_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_ffuel_age

  P = rcemis(EC_F_FFUEL_AGE,k)  &
    + rct(92,k) * xnew(EC_f_ffuel_new)
  
  ! L = 0.0
  
  xnew(EC_f_ffuel_age) = xold(EC_f_ffuel_age) + dt2 * P


!-> EC_c_ffuel

  P = rcemis(EC_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(EC_c_ffuel) = xold(EC_c_ffuel) + dt2 * P


!-> REMPPM25

  P = rcemis(REMPPM25,k)
  
  ! L = 0.0
  
  xnew(REMPPM25) = xold(REMPPM25) + dt2 * P


!-> REMPPM_c

  P = rcemis(REMPPM_C,k)
  
  ! L = 0.0
  
  xnew(REMPPM_c) = xold(REMPPM_c) + dt2 * P


!-> OM25_BGND

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_BGND! xnew(OM25_BGND) = max(0.0, xold(OM25_BGND))


!-> OM25_p

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_p! xnew(OM25_p) = max(0.0, xold(OM25_p))


!-> ffire_OM

  P = rcemis(FFIRE_OM,k)
  
  ! L = 0.0
  
  xnew(ffire_OM) = xold(ffire_OM) + dt2 * P


!-> ffire_BC

  P = rcemis(FFIRE_BC,k)
  
  ! L = 0.0
  
  xnew(ffire_BC) = xold(ffire_BC) + dt2 * P


!-> ffire_remPPM25

  P = rcemis(FFIRE_REMPPM25,k)
  
  ! L = 0.0
  
  xnew(ffire_remPPM25) = xold(ffire_remPPM25) + dt2 * P


!-> ffire_c

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for ffire_c! xnew(ffire_c) = max(0.0, xold(ffire_c))


!-> SeaSalt_f

  P = rcemis(SEASALT_F,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_f) = xold(SeaSalt_f) + dt2 * P


!-> SeaSalt_c

  P = rcemis(SEASALT_C,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_c) = xold(SeaSalt_c) + dt2 * P


!-> Ash_f

  P = rcemis(ASH_F,K)
  
  ! L = 0.0
  
  xnew(Ash_f) = xold(Ash_f) + dt2 * P


!-> Ash_c

  P = rcemis(ASH_C,K)
  
  ! L = 0.0
  
  xnew(Ash_c) = xold(Ash_c) + dt2 * P


