!-> SO4

  P = rcemis(SO4,k)  &
    + AQRCK(ICLOHSO2,K)*2.00e-12 * xnew(OH) * xnew(SO2)  &
    + AQRCK(ICLRC1,K) * xnew(SO2) * xnew(H2O2)  &
    + AQRCK(ICLRC2,K) * xnew(SO2) * xnew(O3)  &
    + AQRCK(ICLRC3,K) * xnew(SO2)
  
  ! L = 0.0
  
  xnew(SO4) = xold(SO4) + dt2 * P


!-> NH3

  P = rcemis(NH3,k)
  
  ! L = 0.0
  
  xnew(NH3) = xold(NH3) + dt2 * P


!-> NO3_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NO3_f! xnew(NO3_f) = max(0.0, xold(NO3_f))


!-> NO3_c

  P = rct(96,k) * xnew(HNO3)  &
    + rct(97,k) * xnew(HNO3)
  
  ! L = 0.0
  
  xnew(NO3_c) = xold(NO3_c) + dt2 * P


!-> NH4_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NH4_f! xnew(NH4_f) = max(0.0, xold(NH4_f))


!-> Ash_f

  P = rcemis(ASH_F,K)
  
  ! L = 0.0
  
  xnew(Ash_f) = xold(Ash_f) + dt2 * P


!-> Ash_c

  P = rcemis(ASH_C,K)
  
  ! L = 0.0
  
  xnew(Ash_c) = xold(Ash_c) + dt2 * P


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
  
  L = rct(119,k)
  
  xnew(EC_f_WOOD_new) = (xold(EC_f_WOOD_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_WOOD_age

  P = rcemis(EC_F_WOOD_AGE,k)  &
    + rct(119,k) * xnew(EC_f_WOOD_new)
  
  ! L = 0.0
  
  xnew(EC_f_WOOD_age) = xold(EC_f_WOOD_age) + dt2 * P


!-> EC_c_WOOD

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for EC_c_WOOD! xnew(EC_c_WOOD) = max(0.0, xold(EC_c_WOOD))


!-> EC_f_FFUEL_new

  P = rcemis(EC_F_FFUEL_NEW,k)
  
  L = rct(119,k)
  
  xnew(EC_f_FFUEL_new) = (xold(EC_f_FFUEL_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_FFUEL_age

  P = rcemis(EC_F_FFUEL_AGE,k)  &
    + rct(119,k) * xnew(EC_f_FFUEL_new)
  
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


!-> OM25_BGND

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_BGND! xnew(OM25_BGND) = max(0.0, xold(OM25_BGND))


!-> OM25_p

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OM25_p! xnew(OM25_p) = max(0.0, xold(OM25_p))


!-> FFIRE_CO

  P = rcemis(FFIRE_CO,k)
  
  L = rct(120,k) * xnew(OH)
  
  xnew(FFIRE_CO) = (xold(FFIRE_CO) + dt2 * P) / (1.0 + dt2 * L)


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


!-> SeaSalt_f

  P = rcemis(SEASALT_F,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_f) = xold(SeaSalt_f) + dt2 * P


!-> SeaSalt_c

  P = rcemis(SEASALT_C,K)
  
  ! L = 0.0
  
  xnew(SeaSalt_c) = xold(SeaSalt_c) + dt2 * P


