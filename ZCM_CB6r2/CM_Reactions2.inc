!-> ASOC_ng10

  P = YCOXY(-2)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YCBENZ(-2)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YCTOL(-2)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(ASOC_ng10) = xold(ASOC_ng10) + dt2 * P


!-> ASOC_ng1e2

  P = YCOXY(-1)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YCBENZ(-1)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YCTOL(-1)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(ASOC_ng1e2) = xold(ASOC_ng1e2) + dt2 * P


!-> ASOC_ug1

  P = YCOXY(0)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YCBENZ(0)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YCTOL(0)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(ASOC_ug1) = xold(ASOC_ug1) + dt2 * P


!-> ASOC_ug10

  P = YCOXY(1)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YCBENZ(1)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YCTOL(1)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(ASOC_ug10) = xold(ASOC_ug10) + dt2 * P


!-> ASOC_ug1e2

  P = YCOXY(2)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YCBENZ(2)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YCTOL(2)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(ASOC_ug1e2) = xold(ASOC_ug1e2) + dt2 * P


!-> ASOC_ug1e3

  P = YCOXY(3)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YCBENZ(3)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YCTOL(3)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(ASOC_ug1e3) = xold(ASOC_ug1e3) + dt2 * P


!-> NON_C_ASOA_ng10

  P = YNOXY(-2)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YNBENZ(-2)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YNTOL(-2)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ng10) = xold(NON_C_ASOA_ng10) + dt2 * P


!-> NON_C_ASOA_ng1e2

  P = YNOXY(-1)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YNBENZ(-1)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YNTOL(-1)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ng1e2) = xold(NON_C_ASOA_ng1e2) + dt2 * P


!-> NON_C_ASOA_ug1

  P = YNOXY(0)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YNBENZ(0)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YNTOL(0)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ug1) = xold(NON_C_ASOA_ug1) + dt2 * P


!-> NON_C_ASOA_ug10

  P = YNOXY(1)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YNBENZ(1)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YNTOL(1)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ug10) = xold(NON_C_ASOA_ug10) + dt2 * P


!-> NON_C_ASOA_ug1e2

  P = YNOXY(2)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YNBENZ(2)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YNTOL(2)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ug1e2) = xold(NON_C_ASOA_ug1e2) + dt2 * P


!-> NON_C_ASOA_ug1e3

  P = YNOXY(3)* rct(124,k) * xnew(XYL) * xnew(OH)  &
    + YNBENZ(3)* rct(121,k) * xnew(BENZ) * xnew(OH)  &
    + YNTOL(3)* rct(123,k) * xnew(TOL) * xnew(OH)
  
  ! L = 0.0
  
  xnew(NON_C_ASOA_ug1e3) = xold(NON_C_ASOA_ug1e3) + dt2 * P


!-> BSOC_ng10

  P = YCISOP(-2)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(-2)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YCTERP(-2)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YCTERP(-2)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(BSOC_ng10) = xold(BSOC_ng10) + dt2 * P


!-> BSOC_ng1e2

  P = YCISOP(-1)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(-1)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YCTERP(-1)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YCTERP(-1)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(BSOC_ng1e2) = xold(BSOC_ng1e2) + dt2 * P


!-> BSOC_ug1

  P = YCISOP(0)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(0)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YCTERP(0)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YCTERP(0)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(BSOC_ug1) = xold(BSOC_ug1) + dt2 * P


!-> BSOC_ug10

  P = YCISOP(1)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(1)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YCTERP(1)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YCTERP(1)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(BSOC_ug10) = xold(BSOC_ug10) + dt2 * P


!-> BSOC_ug1e2

  P = YCISOP(2)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(2)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YCTERP(2)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YCTERP(2)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(BSOC_ug1e2) = xold(BSOC_ug1e2) + dt2 * P


!-> BSOC_ug1e3

  P = YCISOP(3)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YCTERP(3)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YCTERP(3)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YCTERP(3)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(BSOC_ug1e3) = xold(BSOC_ug1e3) + dt2 * P


!-> NON_C_BSOA_ng10

  P = YNISOP(-2)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(-2)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YNTERP(-2)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YNTERP(-2)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ng10) = xold(NON_C_BSOA_ng10) + dt2 * P


!-> NON_C_BSOA_ng1e2

  P = YNISOP(-1)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(-1)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YNTERP(-1)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YNTERP(-1)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ng1e2) = xold(NON_C_BSOA_ng1e2) + dt2 * P


!-> NON_C_BSOA_ug1

  P = YNISOP(0)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(0)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YNTERP(0)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YNTERP(0)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ug1) = xold(NON_C_BSOA_ug1) + dt2 * P


!-> NON_C_BSOA_ug10

  P = YNISOP(1)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(1)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YNTERP(1)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YNTERP(1)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ug10) = xold(NON_C_BSOA_ug10) + dt2 * P


!-> NON_C_BSOA_ug1e2

  P = YNISOP(2)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(2)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YNTERP(2)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YNTERP(2)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ug1e2) = xold(NON_C_BSOA_ug1e2) + dt2 * P


!-> NON_C_BSOA_ug1e3

  P = YNISOP(3)* rct(103,k) * xnew(OH) * xnew(C5H8)  &
    + YNTERP(3)* rct(118,k) * xnew(OH) * xnew(TERP)  &
    + YNTERP(3)* rct(119,k) * xnew(O3) * xnew(TERP)  &
    + YNTERP(3)* rct(120,k) * xnew(NO3) * xnew(TERP)
  
  ! L = 0.0
  
  xnew(NON_C_BSOA_ug1e3) = xold(NON_C_BSOA_ug1e3) + dt2 * P


!-> POM_f_wood

  P = rcemis(POM_F_WOOD,k)
  
  ! L = 0.0
  
  xnew(POM_f_wood) = xold(POM_f_wood) + dt2 * P


!-> POM_c_wood

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for POM_c_wood! xnew(POM_c_wood) = max(0.0, xold(POM_c_wood))


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
  
  L = rct(140,k)
  
  xnew(EC_f_wood_new) = (xold(EC_f_wood_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_wood_age

  P = rcemis(EC_F_WOOD_AGE,k)  &
    + rct(140,k) * xnew(EC_f_wood_new)
  
  ! L = 0.0
  
  xnew(EC_f_wood_age) = xold(EC_f_wood_age) + dt2 * P


!-> EC_c_wood

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for EC_c_wood! xnew(EC_c_wood) = max(0.0, xold(EC_c_wood))


!-> EC_f_ffuel_new

  P = rcemis(EC_F_FFUEL_NEW,k)
  
  L = rct(140,k)
  
  xnew(EC_f_ffuel_new) = (xold(EC_f_ffuel_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_ffuel_age

  P = rcemis(EC_F_FFUEL_AGE,k)  &
    + rct(140,k) * xnew(EC_f_ffuel_new)
  
  ! L = 0.0
  
  xnew(EC_f_ffuel_age) = xold(EC_f_ffuel_age) + dt2 * P


!-> EC_c_ffuel

  P = rcemis(EC_C_FFUEL,k)
  
  ! L = 0.0
  
  xnew(EC_c_ffuel) = xold(EC_c_ffuel) + dt2 * P


!-> pSO4f

  P = rcemis(PSO4F,k)
  
  ! L = 0.0
  
  xnew(pSO4f) = xold(pSO4f) + dt2 * P


!-> pSO4c

  P = rcemis(PSO4C,k)
  
  ! L = 0.0
  
  xnew(pSO4c) = xold(pSO4c) + dt2 * P


!-> remPPM25

  P = rcemis(REMPPM25,k)
  
  ! L = 0.0
  
  xnew(remPPM25) = xold(remPPM25) + dt2 * P


!-> remPPM_c

  P = rcemis(REMPPM_C,k)
  
  ! L = 0.0
  
  xnew(remPPM_c) = xold(remPPM_c) + dt2 * P


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


