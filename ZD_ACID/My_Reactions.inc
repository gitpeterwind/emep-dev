 ! ds - MOD
    do iter = 1, toiter(k)
 



!-> NO 

      P =   &

      rcphot(IDNO2,K) * xnew(NO2 )  &
      + 0.95*rcemis(QRCNO,K) 

      L =  &

      rct(1,k)* xn_2d_bgn(IXBGN_O3,k )  
!if ((i==2).and.(j==2).and.(k==20))then
!write(*,*)'Ozone',xn_2d_bgn(IXBGN_O3,k)*PPBINV/amk(k)
!write(*,*)'H2O2',xn_2d_bgn(IXBGN_H2O2,k)*PPBINV/amk(k)
!endif
       xnew(NO)=  amax1(0.0, ( xold(NO) + dt2 * P))  &
                     /(1.0 + dt2*L ) 


!-> NO2 

      P =   &

      rct(1,k) * xn_2d_bgn(IXBGN_O3,k ) * xnew(NO )  &
!      + rctroe(6,k) * xnew(PAN )  &
      + rcmisc(16,k) * xnew(PAN )  &
      + 0.05*rcemis(QRCNO,K) 

      L =  &

     rct(2,k)* xn_2d_bgn(IXBGN_O3,k )  &
!      + rctroe(5,k)* xn_2d_bgn(IXBGN_CH3COO2,k )  &
!      + rctroe(4,k)* xn_2d_bgn(IXBGN_OH,k )  &
      + rcmisc(14,k)* xn_2d_bgn(IXBGN_OH,k ) &
      + rcmisc(15,k)* xn_2d_bgn(IXBGN_CH3COO2,k )  &
      + rcphot(IDNO2,K) 


       xnew(NO2)=  amax1(0.0, ( xold(NO2) + dt2 * P))  &
                     /(1.0 + dt2*L ) 


!-> PAN 

      P =   &

!     rctroe(5,k) * xn_2d_bgn(IXBGN_CH3COO2,k ) * xnew(NO2 ) 
     rcmisc(15,k) * xn_2d_bgn(IXBGN_CH3COO2,k ) * xnew(NO2 ) 

      L =  &

!     rctroe(6,k) 
     rcmisc(16,k) 
       xnew(PAN)=  amax1(0.0, ( xold(PAN) + dt2 * P))  &
                     /(1.0 + dt2*L ) 


!-> HNO3 

      P =   &

!      rctroe(4,k) * xnew(NO2 ) * xn_2d_bgn(IXBGN_OH,k ) 
     rcmisc(14,k) * xnew(NO2 ) *  xn_2d_bgn(IXBGN_OH,k ) &
!added
      + rct(2,k) * xn_2d_bgn(IXBGN_O3,k) * xnew(NO2)

      L = rcmisc(10,k)    ! Coarse pNO3 formation
       xnew(HNO3)=  amax1(0.0, ( xold(HNO3) + dt2 * P))  &
                     /(1.0 + dt2*L ) 


!-> SO2 

!hf
     psd_h2o2 = xn_2d_bgn(IXBGN_H2O2,k)*xn_2d_bgn(IXBGN_H2O2,k)/&
               (xn_2d_bgn(IXBGN_H2O2,k)+ xnew(SO2))

      P =   &

     (1-0.05)*rcemis(QRCSO2,K) 

      L =  &

!       rcmisc(9,k)  ! Eliassen+Saltbones
      2e-12*aqrck(ICLOHSO2,K)* xn_2d_bgn(IXBGN_OH,k )  &
!hf       + aqrck(ICLRC1,K)*xn_2d_bgn(IXBGN_H2O2,k)  &
       + aqrck(ICLRC1,K)*psd_h2o2  &
       + aqrck(ICLRC2,K)*xn_2d_bgn(IXBGN_O3,k)  &
       + aqrck(ICLRC3,K)    ! Fe-catalysed reaction


       xnew(SO2)=  amax1(0.0, ( xold(SO2) + dt2 * P))  &
                     /(1.0 + dt2*L ) 


 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   end do !! End iterations
 ! Just before SO4
 !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


!-> SO4 

      P =   &

      2e-12*aqrck(ICLOHSO2,K) * xn_2d_bgn(IXBGN_OH,k ) * xnew(SO2 )  &
!hfOH       + aqrck(ICLRC1,K)*xn_2d_bgn(IXBGN_H2O2,k) * xnew(SO2 )  &
       + aqrck(ICLRC1,K)*psd_h2o2 * xnew(SO2 )  &
       + aqrck(ICLRC2,K)*xn_2d_bgn(IXBGN_O3,k) * xnew(SO2 )  &
       + aqrck(ICLRC3,K) * xnew(SO2 )  &
!       rcmisc(9,k) * xnew(SO2) &    !! Eliassen+Saltbones
      + 0.05*rcemis(QRCSO2,K) 


       xnew(SO4)=  amax1(0.0,  xold(SO4) + dt2 * P)



!-> NH3 

      P =   &

     rcemis(QRCNH3,K) 
      ! L = 0.0
 

       xnew(NH3)=  amax1(0.0,  xold(NH3) + dt2 * P)


!-> AMSU 
      ! P = 0.0
 
      ! L = 0.0
 

!ds       xnew(AMSU)=  amax1(0.0, xold(AMSU))


!-> AMNI 
      ! P = 0.0
 
      ! L = 0.0
 

!ds       xnew(AMNI)=  amax1(0.0, xold(AMNI))
!==================================================

!-> PM25

      P =   &
     rcemis(QRCPM25,K) 
      ! L = 0.0
       xnew(PM25)=  amax1(0.0,  xold(PM25) + dt2 * P)


!-> PMCO

      P =   &
     rcemis(QRCPMCO,K) 
      ! L = 0.0
       xnew(PMCO)=  amax1(0.0,  xold(PMCO) + dt2 * P)

!-> pNO3 

      P = rcmisc(10,k)  * xnew(HNO3)
      xnew(pNO3)=  amax1(0.0,  xold(pNO3) + dt2 * P)

!-> Rn222

      P =  rc_Rn222(k)
      L = 1.0/(5.51614*24.0*3600.0)    ! Dietmar, 11/4/2005
       xnew(Rn222)=  amax1(0.0, ( xold(Rn222) + dt2 * P))  &
                     /(1.0 + dt2*L ) 
!dsbug       xnew(Rn222)=  max(0.0,  xold(Rn222) + dt2 * P)

!-> Pb210
      P = 1.0/(5.51614*24.0*3600.0) * xnew(Rn222)    ! Dietmar, 11/4/2005
      ! L = 0.0
       xnew(Pb210)=  max(0.0,  xold(Pb210) + dt2 * P)

!SeaS
   if (SEASALT) then

!-> Sea salt fine

      P =   &
     rcss(QSSFI,K) 
      ! L = 0.0
       xnew(SSfi)=  amax1(0.0,  xold(SSfi) + dt2 * P)

!-> Sea salt coarse

      P =   &
     rcss(QSSCO,K) 
      ! L = 0.0
       xnew(SSco)=  amax1(0.0,  xold(SSco) + dt2 * P)

   endif














