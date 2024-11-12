!-> OH

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OH! xnew(OH) = max(0.0, xold(OH))


!-> HO2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HO2! xnew(HO2) = max(0.0, xold(HO2))


!-> O3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for O3! xnew(O3) = max(0.0, xold(O3))


!-> NO

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NO! xnew(NO) = max(0.0, xold(NO))


!-> NO2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NO2! xnew(NO2) = max(0.0, xold(NO2))


!-> HNO3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HNO3! xnew(HNO3) = max(0.0, xold(HNO3))


!-> SO2

  P = rcemis(SO2,k)
  
  L = 3.5e-5
  
  xnew(SO2) = (xold(SO2) + dt2 * P) / (1.0 + dt2 * L)


!-> SO4

  P = rcemis(SO4,k)  &
    + 3.5e-5 * xnew(SO2)
  
  ! L = 0.0
  
  xnew(SO4) = xold(SO4) + dt2 * P


!-> NH3

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NH3! xnew(NH3) = max(0.0, xold(NH3))


!-> NO3_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NO3_f! xnew(NO3_f) = max(0.0, xold(NO3_f))


!-> NO3_c

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NO3_c! xnew(NO3_c) = max(0.0, xold(NO3_c))


!-> NH4_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NH4_f! xnew(NH4_f) = max(0.0, xold(NH4_f))


!-> CO

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CO! xnew(CO) = max(0.0, xold(CO))


!-> HCHO

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for HCHO! xnew(HCHO) = max(0.0, xold(HCHO))


!-> CH4

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CH4! xnew(CH4) = max(0.0, xold(CH4))


!-> H2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for H2! xnew(H2) = max(0.0, xold(H2))


!-> C2H6

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for C2H6! xnew(C2H6) = max(0.0, xold(C2H6))


!-> NC4H10

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NC4H10! xnew(NC4H10) = max(0.0, xold(NC4H10))


!-> OXYL

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for OXYL! xnew(OXYL) = max(0.0, xold(OXYL))


!-> C5H8

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for C5H8! xnew(C5H8) = max(0.0, xold(C5H8))


!-> APINENE

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for APINENE! xnew(APINENE) = max(0.0, xold(APINENE))


!-> BIOTERP

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for BIOTERP! xnew(BIOTERP) = max(0.0, xold(BIOTERP))


!-> C3H6

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for C3H6! xnew(C3H6) = max(0.0, xold(C3H6))


!-> C2H4

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for C2H4! xnew(C2H4) = max(0.0, xold(C2H4))


!-> C2H5OH

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for C2H5OH! xnew(C2H5OH) = max(0.0, xold(C2H5OH))


!-> CH3OH

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CH3OH! xnew(CH3OH) = max(0.0, xold(CH3OH))


!-> PPM25_FIRE

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for PPM25_FIRE! xnew(PPM25_FIRE) = max(0.0, xold(PPM25_FIRE))


!-> MEK

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MEK! xnew(MEK) = max(0.0, xold(MEK))


!-> GLYOX

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for GLYOX! xnew(GLYOX) = max(0.0, xold(GLYOX))


!-> MGLYOX

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for MGLYOX! xnew(MGLYOX) = max(0.0, xold(MGLYOX))


!-> UNREAC

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for UNREAC! xnew(UNREAC) = max(0.0, xold(UNREAC))


!-> NPP_Cs137

  P = rcemis(NPP_CS137,K)
  
  ! L = 0.0
  
  xnew(NPP_Cs137) = xold(NPP_Cs137) + dt2 * P


!-> NPP_I131

  P = rcemis(NPP_I131,K)
  
  L = 0.994e-06
  
  xnew(NPP_I131) = (xold(NPP_I131) + dt2 * P) / (1.0 + dt2 * L)


!-> NPP_Sr90

  P = rcemis(NPP_SR90,K)
  
  ! L = 0.0
  
  xnew(NPP_Sr90) = xold(NPP_Sr90) + dt2 * P


!-> NPP_Xe133

  P = rcemis(NPP_XE133,K)
  
  L = 0.152e-05
  
  xnew(NPP_Xe133) = (xold(NPP_Xe133) + dt2 * P) / (1.0 + dt2 * L)


!-> NPP_Kr85

  P = rcemis(NPP_KR85,K)
  
  ! L = 0.0
  
  xnew(NPP_Kr85) = xold(NPP_Kr85) + dt2 * P


!-> ASH_1

  P = rcemis(ASH_1,K)
  
  ! L = 0.0
  
  xnew(ASH_1) = xold(ASH_1) + dt2 * P


!-> ASH_2

  P = rcemis(ASH_2,K)
  
  ! L = 0.0
  
  xnew(ASH_2) = xold(ASH_2) + dt2 * P


!-> ASH_3

  P = rcemis(ASH_3,K)
  
  ! L = 0.0
  
  xnew(ASH_3) = xold(ASH_3) + dt2 * P


!-> ASH_4

  P = rcemis(ASH_4,K)
  
  ! L = 0.0
  
  xnew(ASH_4) = xold(ASH_4) + dt2 * P


!-> ASH_5

  P = rcemis(ASH_5,K)
  
  ! L = 0.0
  
  xnew(ASH_5) = xold(ASH_5) + dt2 * P


!-> ASH_6

  P = rcemis(ASH_6,K)
  
  ! L = 0.0
  
  xnew(ASH_6) = xold(ASH_6) + dt2 * P


!-> ASH_7

  P = rcemis(ASH_7,K)
  
  ! L = 0.0
  
  xnew(ASH_7) = xold(ASH_7) + dt2 * P


