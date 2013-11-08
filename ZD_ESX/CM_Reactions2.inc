
!-> SO4 

      P =   &
          1.0e-5 * xnew(SO2 )  &
        + rcemis(SO4,k) 
      ! L = 0.0
 

      xnew(SO4)=  xold(SO4) + dt2 * P

!-> NH3 

      P =   &
          rcemis(NH3,k) 
      ! L = 0.0
 

      xnew(NH3)=  xold(NH3) + dt2 * P

!-> NO3_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NO3_F! xnew(NO3_F)=  max(0.0, xold(NO3_F))

!-> NO3_C 

      P =   &
          rct(62,k) * xnew(HNO3 ) 
      ! L = 0.0
 

      xnew(NO3_C)=  xold(NO3_C) + dt2 * P

!-> NH4_F 
      ! P = 0.0
 
      ! L = 0.0
 
!Nothing to do for NH4_F! xnew(NH4_F)=  max(0.0, xold(NH4_F))
