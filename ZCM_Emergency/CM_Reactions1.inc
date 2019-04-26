!-> OD

  ! P = 0.0
  
  L = rct(2,k)
  
  xnew(OD) = xold(OD) / (1.0 + dt2 * L)


!-> OP

  ! P = 0.0
  
  L = rct(1,k)
  
  xnew(OP) = xold(OP) / (1.0 + dt2 * L)


!-> OH

  P = 2.* rct(2,k) * xnew(OD)  &
    + rct(6,k) * xnew(NO) * xnew(HO2)
  
  L = rct(4,k) * xnew(VOC)  &
    + rct(5,k) * xnew(CH3CHO)  &
    + 1.11e-11 * xnew(NO2)
  
  xnew(OH) = (xold(OH) + dt2 * P) / (1.0 + dt2 * L)


!-> HO2

  P = rct(7,k) * xnew(NO) * xnew(RO2)
  
  L = rct(6,k) * xnew(NO)
  
  xnew(HO2) = (xold(HO2) + dt2 * P) / (1.0 + dt2 * L)


!-> RO2

  P = rct(4,k) * xnew(OH) * xnew(VOC)  &
    + rct(8,k) * xnew(CH3COO2) * xnew(NO)
  
  L = rct(7,k) * xnew(NO)
  
  xnew(RO2) = (xold(RO2) + dt2 * P) / (1.0 + dt2 * L)


!-> TRACER1

  P = rcemis(TRACER1,k)
  
  ! L = 0.0
  
  xnew(TRACER1) = xold(TRACER1) + dt2 * P


!-> TRACER2

  P = rcemis(TRACER2,k)
  
  ! L = 0.0
  
  xnew(TRACER2) = xold(TRACER2) + dt2 * P


!-> O3

  P = rct(1,k) * xnew(OP)
  
  L = rct(3,k) * xnew(NO)
  
  xnew(O3) = (xold(O3) + dt2 * P) / (1.0 + dt2 * L)


!-> NO

  P = rcemis(NO,k)
  
  L = rct(3,k) * xnew(O3)  &
    + rct(6,k) * xnew(HO2)  &
    + rct(7,k) * xnew(RO2)  &
    + rct(8,k) * xnew(CH3COO2)
  
  xnew(NO) = (xold(NO) + dt2 * P) / (1.0 + dt2 * L)


!-> NO2

  P = rct(3,k) * xnew(O3) * xnew(NO)  &
    + rct(6,k) * xnew(NO) * xnew(HO2)  &
    + rct(7,k) * xnew(NO) * xnew(RO2)  &
    + rct(8,k) * xnew(CH3COO2) * xnew(NO)  &
    + rct(9,k) * xnew(PAN)  &
    + rcemis(NO2,k)
  
  L = 1.11e-11 * xnew(OH)  &
    + 4.7e-12 * xnew(CH3COO2)
  
  xnew(NO2) = (xold(NO2) + dt2 * P) / (1.0 + dt2 * L)


!-> SHIPNOX

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for SHIPNOX! xnew(SHIPNOX) = max(0.0, xold(SHIPNOX))


!-> HNO3

  P = 1.11e-11 * xnew(NO2) * xnew(OH)
  
  ! L = 0.0
  
  xnew(HNO3) = xold(HNO3) + dt2 * P


!-> VOC

  P = rcemis(VOC,k)
  
  L = rct(4,k) * xnew(OH)
  
  xnew(VOC) = (xold(VOC) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3CHO

  P = rct(7,k) * xnew(NO) * xnew(RO2)  &
    + rcemis(CH3CHO,k)
  
  L = rct(5,k) * xnew(OH)
  
  xnew(CH3CHO) = (xold(CH3CHO) + dt2 * P) / (1.0 + dt2 * L)


!-> PAN

  P = 4.7e-12 * xnew(CH3COO2) * xnew(NO2)
  
  L = rct(9,k)
  
  xnew(PAN) = (xold(PAN) + dt2 * P) / (1.0 + dt2 * L)


!-> CH3COO2

  P = rct(5,k) * xnew(OH) * xnew(CH3CHO)  &
    + rct(9,k) * xnew(PAN)
  
  L = rct(8,k) * xnew(NO)  &
    + 4.7e-12 * xnew(NO2)
  
  xnew(CH3COO2) = (xold(CH3COO2) + dt2 * P) / (1.0 + dt2 * L)


!-> SO2

  P = rcemis(SO2,k)
  
  ! L = 0.0
  
  xnew(SO2) = xold(SO2) + dt2 * P


!-> SO4

  P = rcemis(SO4,k)
  
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

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NO3_c! xnew(NO3_c) = max(0.0, xold(NO3_c))


!-> NH4_f

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for NH4_f! xnew(NH4_f) = max(0.0, xold(NH4_f))


!-> CO

  P = rcemis(CO,k)
  
  ! L = 0.0
  
  xnew(CO) = xold(CO) + dt2 * P


!-> HCHO

  P = rcemis(HCHO,k)
  
  ! L = 0.0
  
  xnew(HCHO) = xold(HCHO) + dt2 * P


!-> CH4

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for CH4! xnew(CH4) = max(0.0, xold(CH4))


!-> H2

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for H2! xnew(H2) = max(0.0, xold(H2))


!-> C2H6

  P = rcemis(C2H6,k)
  
  ! L = 0.0
  
  xnew(C2H6) = xold(C2H6) + dt2 * P


!-> NC4H10

  P = rcemis(NC4H10,k)
  
  ! L = 0.0
  
  xnew(NC4H10) = xold(NC4H10) + dt2 * P


!-> OXYL

  P = rcemis(OXYL,k)
  
  ! L = 0.0
  
  xnew(OXYL) = xold(OXYL) + dt2 * P


!-> C5H8

  P = rcemis(C5H8,k)
  
  ! L = 0.0
  
  xnew(C5H8) = xold(C5H8) + dt2 * P


!-> APINENE

  P = rcemis(APINENE,k)
  
  ! L = 0.0
  
  xnew(APINENE) = xold(APINENE) + dt2 * P


!-> BIOTERP

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for BIOTERP! xnew(BIOTERP) = max(0.0, xold(BIOTERP))


!-> C3H6

  P = rcemis(C3H6,k)
  
  ! L = 0.0
  
  xnew(C3H6) = xold(C3H6) + dt2 * P


!-> C2H4

  P = rcemis(C2H4,k)
  
  ! L = 0.0
  
  xnew(C2H4) = xold(C2H4) + dt2 * P


!-> C2H5OH

  P = rcemis(C2H5OH,k)
  
  ! L = 0.0
  
  xnew(C2H5OH) = xold(C2H5OH) + dt2 * P


!-> CH3OH

  P = rcemis(CH3OH,k)
  
  ! L = 0.0
  
  xnew(CH3OH) = xold(CH3OH) + dt2 * P


!-> PPM25_FIRE

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for PPM25_FIRE! xnew(PPM25_FIRE) = max(0.0, xold(PPM25_FIRE))


!-> MEK

  P = rcemis(MEK,k)
  
  ! L = 0.0
  
  xnew(MEK) = xold(MEK) + dt2 * P


!-> GLYOX

  P = rcemis(GLYOX,k)
  
  ! L = 0.0
  
  xnew(GLYOX) = xold(GLYOX) + dt2 * P


!-> MGLYOX

  P = rcemis(MGLYOX,k)
  
  ! L = 0.0
  
  xnew(MGLYOX) = xold(MGLYOX) + dt2 * P


!-> UNREAC

  ! P = 0.0
  
  ! L = 0.0
  
  ! Nothing to do for UNREAC! xnew(UNREAC) = max(0.0, xold(UNREAC))


!-> V1702A02B_1

  P = rcemis(V1702A02B_1,k)
  
  ! L = 0.0
  
  xnew(V1702A02B_1) = xold(V1702A02B_1) + dt2 * P


!-> V1702A02B_2

  P = rcemis(V1702A02B_2,k)
  
  ! L = 0.0
  
  xnew(V1702A02B_2) = xold(V1702A02B_2) + dt2 * P


!-> V1702A02B_3

  P = rcemis(V1702A02B_3,k)
  
  ! L = 0.0
  
  xnew(V1702A02B_3) = xold(V1702A02B_3) + dt2 * P


!-> V1702A02B_4

  P = rcemis(V1702A02B_4,k)
  
  ! L = 0.0
  
  xnew(V1702A02B_4) = xold(V1702A02B_4) + dt2 * P


!-> V1702A02B_5

  P = rcemis(V1702A02B_5,k)
  
  ! L = 0.0
  
  xnew(V1702A02B_5) = xold(V1702A02B_5) + dt2 * P


!-> V1702A02B_6

  P = rcemis(V1702A02B_6,k)
  
  ! L = 0.0
  
  xnew(V1702A02B_6) = xold(V1702A02B_6) + dt2 * P


!-> V1702A02B_7

  P = rcemis(V1702A02B_7,k)
  
  ! L = 0.0
  
  xnew(V1702A02B_7) = xold(V1702A02B_7) + dt2 * P


