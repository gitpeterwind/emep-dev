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


!-> EC_f_new

  P = rcemis(EC_F_NEW,k)
  
  L = rct(80,k)
  
  xnew(EC_f_new) = (xold(EC_f_new) + dt2 * P) / (1.0 + dt2 * L)


!-> EC_f_age

  P = rcemis(EC_F_AGE,k)  &
    + rct(80,k) * xnew(EC_f_new)
  
  ! L = 0.0
  
  xnew(EC_f_age) = xold(EC_f_age) + dt2 * P


