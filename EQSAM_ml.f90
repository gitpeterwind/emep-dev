module EQSAM_ml

 implicit none
 private


 !/- subroutines:
 public   ::  aero



!_______________________________________________________________________________
!EQSAM_v1.0
!programmed by swen metzger 3/11/99
!included in martijn schaaps lotus version 11/10/2000
!
!purpose
!-------
!simplified gas/aerosol partioning and aerosol associated water mass
!calculation based on aerosol composition parameterizations
!
!interface
!---------
!call  aero(c)
!
!method
!------
!equilibrium / internal mixture assumption
!System: NH3,NH4+/H2SO4+,SO4--/HNO3,NO3-, H2O
!
!external
!--------
!none
!
!      reference
!      ---------
!      Swen Metzger Ph.D Thesis, University Utrecht, 2000
!         http://www.mpch-mainz.mpg.de/~metzger
!
!      Metzger, S. M., F. J. Dentener, J. Lelieveld, and S. N. Pandis,
!         GAS/AEROSOL PARTITIONING I: A COMPUTATIONALLY EFFICIENT MODEL, JGR, in press.
!      Metzger, S. M., F. J. Dentener, A. Jeuken, and M. Krol, J. Lelieveld,
!         GAS/AEROSOL PARTITIONING II: GLOBAL MODELING RESULTS, JGR, in press.
!
!      email: metzger@mpch-mainz.mpg.de
!_______________________________________________________________________________



!hfmodule m_aero



 contains

!hfsubroutine aero(c,ah2o)

!>-------------------------------------------------------------------------------<
subroutine aero(so4in, hno3in,no3in,nh3in,nh4in,rh,temp,   &
                  aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, ERRMARK,deb) 

!>-------------------------------------------------------------------------------<

  use ModelConstants_ml,  only :  KMAX_MID, KCHEMTOP

!use dims, only : nx, ny, nz, nspec ! gevlagd door Paul Smeets, okt, 2002
!use indices                        !
!use meteo, only : temp, rh

implicit none
 real, intent(in):: temp(KCHEMTOP:KMAX_MID),rh(KCHEMTOP:KMAX_MID)

!hf real :: c(nx,ny,nz,nspec), ah2o(nx,ny,nz)
 !.. local
  real,intent(in)    :: so4in(KCHEMTOP:KMAX_MID), &
             no3in(KCHEMTOP:KMAX_MID), &
             nh4in(KCHEMTOP:KMAX_MID), &
             hno3in(KCHEMTOP:KMAX_MID), &
             nh3in(KCHEMTOP:KMAX_MID)


   real,intent(out):: &
             aSO4out(KCHEMTOP:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), & 
             gNH3out(KCHEMTOP:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID)

  integer :: i,j,k,errmark,irh
  logical, intent(in) :: deb

   real ::  ph2o, psat, p, tso4, hno3, nh3, & !hf Changes :
                                           !hno3in->hno3, nh3in->nh3 
                                           !in order not to confuse variables
         ghno3, gnh3, aso4, ano3, anh4, &
         frac, tnh4, tno3, sulf, hfac,pref2,tk

real, parameter :: mwso4=96.0576
real, parameter :: mwhno3=63.01287
real, parameter :: mwno3=62.0049
real, parameter :: mwnh3=17.03061
real, parameter :: mwnh4=18.03061
real, parameter :: mwna=23.0
real, parameter :: mwcl=35.5
real, parameter :: tref = 298.15
real, parameter :: pref = 3169.0
real, parameter :: Href = 43990
real, parameter :: R = 8.31451

real, parameter :: RR=82.0567e-6,RHMAX=0.9999,RHMIN=0.0001
real, parameter :: GF0=0.81,GF1=0.25,GF2=0.50,MCH20=1000.*GF0,ZERO=0.
real, parameter :: MWAN=80.,MWSA=98.,MWAB=115.,MWAS=132.,MWLC=MWAB+MWAS
real, parameter :: ZWAN=1.00,ZWSA=0.50,ZWAB=1.00,ZWAS=0.75,ZWLC=1.00
real, parameter :: NWAN=4.00,NWSA=4.00,NWAB=4.00,NWAS=27.0,NWLC=9.00

real :: T0,RHD,KAN,X0,X1,X2,X3,X4,X5,X6,XK10,PY,RATIONS,ZKAN,rhl
real :: ZIONIC,ZFLAG,COEF,T0T,TT,PH,WH2O,GNO3,GSO4,PM
real :: GAMA,CAS0,CAN0,CSA0,CAB0,CLC0,GG,HPLUS,CONVT
real :: w(8)

real :: RHDX(8),RHDZ(8)    ! RHD / MRHD arrays for different aerosol types
! RHD / MRHD values as of ISORROPIA / SCAPE (T=298.15K)
real, parameter :: RHDA(8) = &
(/0.32840,0.4906,0.6183,0.7997,0.67500,0.5000,0.4000,0.0000/)
real, parameter :: RHDE(8) = &
(/-1860.0,-431.0,852.00,80.000,262.000,3951.0,384.00,0.0000/)


T0=tref

!hf do i=1,nx
!hf do j=1,ny
!hf do k=1,nz
  do k=KCHEMTOP,KMAX_MID
   tk = temp(k) !hf temp(i,j,k)
   rhl= rh(k)   !rh(i,j,k)/100

! calculation aerosolconcentrations:
! 1 ppb= 2.46e10 molec/cm**3 = 0.0409 * (MW) * p/pref microgram/m**3
! c(i_so4) is the sulphuric acid concentration in the gasfase

!hf pref2 = 101325
!hf p     = pref2
!hf CONVT = 0.0409 * p/pref2 !hf MARS calc Temp dependent CONVT


! compute gas-phase in umol/m3
!hf  hno3in = c(i,j,k,i_hno3)  * CONVT
!hf  nh3in  = c(i,j,k,i_nh3)   * CONVT
errmark=1

! compute total aerosol
!hf tso4 = c(i,j,k,i_so4)*  CONVT &
!hf     + c(i,j,k,i_so4a)/mwso4
!hf tno3 = c(i,j,k,i_no3a)/mwno3 + hno3in
!hf tnh4 = c(i,j,k,i_nh4a)/mwnh4 + nh3in

!hf already in umol/m3
    tso4 = so4in(k)
    tno3 = no3in(k) + hno3in(k)
    tnh4 = nh4in(k) + nh3in(k)
! All sulphate is assumed to be in the particulate phase, hense H2So4 (g)
!is assumed to condense

!
!______________________________________________
!  meteorology
TT = tk                        ! T                      [K]
!     RH = rh                  ! RH                     [0-1]
!     py = p/100.              ! p                      [hPa]
!
! gas+aerosol:
!
w(1) = 0.0 !c(i,j,k,i_na)/mwna    ! Na+ (ss  + xsod) (a)   [umol/m^3]
w(2) = tso4                  ! H2SO4    + SO4-- (p)   [umol/m^3]
w(3) = tnh4                  ! NH3  (g) + NH4+  (p)   [umol/m^3]
w(4) = tno3                  ! HNO3 (g) + NO3-  (p)   [umol/m^3]
w(5) = 0.0 !c(i,j,k,i_cl)/mwcl    ! HCl  (g)               [umol/m^3]

! ions (account for particles with 0.01 > r < 1 micrometer):
w(6) = 0.                      ! K+   (p) from Dust     [umol/m^3]
w(7) = 0.                      ! Ca++ (p) from Dust     [umol/m^3]
w(8) = 0.                      ! Mg++ (p) from Dust     [umol/m^3]
!______________________________________________
!
zflag=1.
!
! account for sea salt and mineral dust

!     w(2)=w(2) + 0.25 * w(1)  ! sulfate concentration [umol/m^3]
!     w(5)=w(5) + 1.80 * w(1)  ! HCl concentration     [umol/m^3]
!     w(6)=w(6) + 0.00 * w(1)  ! K+ concentration[umol/m^3]
!     w(7)=w(7) + 0.04 * w(1)  ! Ca++ concentration    [umol/m^3]
!     w(8)=w(8) + 0.12 * w(1)  ! Mg++ concentration    [umol/m^3]
!
w(:) = w(:) * 1.0e-6
!
TSO4  = w(2)             ! total input sulfate  (g+p)
TNH4  = w(3)             ! total input ammonium (g+p)
TNO3  = w(4)             ! total input nitrate  (g+p)
!
! SULFATE RICH
!
if((w(1)+w(3)+w(6)+2.*(w(7)+w(8))).le.(2.*w(2))) then
    zflag=3.
endif
!
! SULFATE VERY RICH CASE if (NH4+Na+K+2(Ca+Mg))/SO4 < 1
!
if((w(1)+w(3)+w(6)+2.*(w(7)+w(8))).lt.w(2)) then
    zflag=4.
endif
!
! SULFATE NEUTRAL CASE
!
if((w(1)+w(3)+w(6)+2.*(w(7)+w(8))).gt.(2.*w(2))) then
    zflag=2.
endif
!
! SULFATE POOR AND CATION POOR CASE
!
if((w(1)+w(6)+2.*(w(7)+w(8))).gt.(2.*w(2))) then
!
    zflag=1.
!
endif
!
IF ( rhl .LT. RHMIN ) rhl=RHMIN
IF ( rhl .GT. RHMAX ) rhl=RHMAX

! CALCULATE TEMPERATURE DEPENDENCY FOR SOME RHDs
do irh = 1,8
 RHDX(irh)=RHDA(irh)*exp(RHDE(irh)*(1./TT-1./T0))
 RHDZ(irh)=RHDX(irh)
enddo
!
! GET WATER ACTIVITIES ACCORDING TO METZGER ET AL. 1999
!
CAN0=((MCH20*NWAN*(1./rhl-1.))/MWAN)**ZWAN  !  NH4NO3
CSA0=((MCH20*NWSA*(1./rhl-1.))/MWSA)**ZWSA  !  2H-SO4
CAB0=((MCH20*NWAB*(1./rhl-1.))/MWAB)**ZWAB  !  NH4HSO4
CAS0=((MCH20*NWAS*(1./rhl-1.))/MWAS)**ZWAS  ! (NH4)2SO4
CLC0=((MCH20*NWLC*(1./rhl-1.))/MWLC)**ZWLC  ! (NH4)3H(SO4)2
!
!  GET DELIQUESENCE RELATIVE HUMIDITY (ISORROPIA/SCAPE)
!
RHD = 0.6183
if(TT.ne.298.0.or.TT.ne.T0) RHD=RHD*exp(852.*(1./TT-1./T0))
!
! GET MEAN MOLAL IONIC ACTIVITY COEFF ACCORDING TO METZGER ET AL. 1999
!
GAMA=0.0
IF(rhl.GE.RHD)  GAMA=(rhl**ZFLAG/(1000./ZFLAG*(1.-rhl)+ZFLAG))
!
!     NH4NO3  :              GAMA=GAMA**GF1  (GF1=1/4)
!    (NH4)2S04:              GAMA=GAMA**GF2  (GF2=1/2)
!    (NH4)3H(S04)2  :        GAMA=GAMA**GF3  (GF3=2/5)
!     CaS04   :              GAMA=GAMA**GF4  (GF4=1/1)
!
! GET EQ-const. NH4NO3(s) <==> NH3(g) + HNO3(g) [mol^2/kg] (ISORROPIA)
!
T0T=T0/TT
COEF=1.0+LOG(T0T)-T0T
!
XK10 = 5.746e-17
XK10= XK10 * EXP(-74.38*(T0T-1.0) + 6.120*COEF)
!
KAN = XK10/(RR*TT)/(RR*TT)
!
! DEFINE AQUEOUSE PHASE (NO SOLID NH4NO3 IF NO3/SO4>1, TEN BRINK, ET AL.,
! 1996, ATMOS ENV, 24, 4251-4261)
!
IF(rhl.LT.RHD.AND.ZFLAG.GT.2.) GAMA=0.  ! NH4NO3 EQ ABOVE DRY AEROSOL
!     IF(rhl.LT.RHD) GAMA=0.
ZKAN=1.-GAMA**GF1
!
! ACCOUNT FOR VARIOUS AMMOMIUM SULFATE SALTS ACCORDING TO MEAN VALUE AS OF
!ISORROPIA
GG=2.0                          ! (NH4)2SO4 IS THE PREFFERED SPECIES
!FOR SULFATE DEFICIENT CASES
IF(ZFLAG.EQ.3.) THEN
   IF(rhl.LE.RHDZ(7)) THEN       ! ACCOUNT FOR MIXTURE OF (NH4)2SO4(s) &NH4HSO4(s) & (NH4)3H(SO4)2(s)
      GG=1.5
   ELSEIF(rhl.GT.RHDZ(7).AND.rhl.LE.RHDZ(5)) THEN ! MAINLY (NH4)2SO4(s) &(NH4)3H(SO4)2(s)
      GG=1.75
!     GG=1.5
   ELSEIF(rhl.GE.RHDZ(5)) THEN   ! (NH4)2SO4(S) & NH4HSO4(S) & SO4-- & HSO4-
      GG=1.5
   ENDIF
ENDIF
IF(ZFLAG.EQ.4.) GG=1.0          ! IF SO4 NEUTRALIZED, THEN ONLY AS NH4HSO4(S) OR  HSO4- / H2SO4
!
! CALCULATE RHD DEPENDENT EQ: NH4NO3(s) <==> NH3(g) + HNO3(g) (ISORROPIA)
!
X0   = MAX(ZERO,MIN(TNH4,GG*TSO4))      ! MAX AMMOMIUM SULFATE
X1   = MAX(ZERO,MIN(TNH4-X0,TNO3))      ! MAX AMMOMIUM NITRATE
X2   = MAX(TNH4-X1-X0,ZERO)             ! RES NH3
X3   = MAX(TNO3-X1,ZERO)                ! RES HNO3
X4   = X2 + X3
!
X5   = SQRT(X4*X4+KAN*ZKAN)
X6   = 0.5*(-X4+X5)
X6   = MIN(X1,X6)
!
GSO4 = TSO4-X0/GG                       ! H2SO4
GNH3 = X2+X6
GNO3 = X3+X6
!
! ELIMINATE SMALL NEGATIVE VALUES IF THEY OCCUR DUE TO NUMERICS
!
IF(GSO4.LT.0.)   GSO4=0.
IF(GNH3.LT.0.)   GNH3=0.
IF(GNO3.LT.0.)   GNO3=0.
IF(GSO4.GT.TSO4) GSO4=TSO4
IF(GNH3.GT.TNH4) GNH3=TNH4
IF(GNO3.GT.TNO3) GNO3=TNO3
!
! ASSUME SUPERSATURATED SOLUTIONS RATHER THAN SOLIDS
!
ASO4 = TSO4-GSO4
ANH4 = TNH4-GNH3
ANO3 = TNO3-GNO3
!
! GET AEROSOL WATER/PM + pH
!
!     MWAN=80.  /     MWSA=98.  /    MWAB=115.   /MWAS=132.   / MWLC=MWAB+MWAS
! CAN0 = NH4NO3 / CSA0 = 2H-SO4 / CAB0 = NH4HSO4 / CAS0 = (NH4)2SO4 / CLC0 = (NH4)3H(SO4)2
!
! note that the water fractions (i.e. the ratios e.g. ASO4/CAS0)
! of various aerosol compositions (domains) were determined to yield
! the best agreement with ISORROPOIA (were it is explizitly calculated).
! for consistence, the absolute mass of the total particulate matter (PM),
! is determined in the same way. SM 17/2/2000.
!
! CALCULATE AEROSOL WATER
!
!     MWAN=80.  /     MWSA=98.  /    MWAB=115.   /MWAS=132.   / MWLC=MWAB+MWAS
! CAN0 = NH4NO3 / CSA0 = 2H-SO4 / CAB0 = NH4HSO4 / CAS0 = (NH4)2SO4 / CLC0= (NH4)3H(SO4)2
!
WH2O = 0.
IF(ZFLAG.EQ.1.) WH2O = ASO4/CAS0 + GSO4/CAB0 + ANO3/CAN0
IF(ZFLAG.EQ.2.) WH2O = ASO4/CLC0 + GSO4/CAB0 + ANO3/CAN0
IF(ZFLAG.EQ.3.) WH2O = ASO4/CAB0 + GSO4/CSA0 + ANO3/CAN0
IF(ZFLAG.EQ.4.) WH2O = ASO4/CAB0 + GSO4/CSA0 + ANO3/CAN0
!
! uncomment if needed
!
!     PM   = 0.
!     PM=ANH4*MWNH4+ASO4*MWSO4+ANO3*MWNO3
!     PM   = PM   * 1.e6                       ! total PM, excl. H20 [ug/m^3]
!
!hf uncommented conversion of water CORRECT units???
     WH2O = WH2O * 1.e9                       ! aerosol water [ug/m^3]
     IF(WH2O.LT.1.e-3) WH2O=0.
!
!     PH=7.
!     HPLUS=0.
!     IF(WH2O.GT.0.)  HPLUS=(2.*TSO4+ANO3-ANH4)*1000./WH2O
!     IF(HPLUS.GT.0.) PH=-ALOG10(HPLUS)
!
! IONIC STRENGTH [nmol/m^3 air]
!
!     ZIONIC=0.
!     IF(WH2O.GT.0.) ZIONIC=0.5*(ANH4+ANO3+ASO4*4.)*1.e3/WH2O! ionic strength [moles/kg]
!     ZIONIC=min(ZIONIC,200.0)                         ! limit for output
!     ZIONIC=max(ZIONIC,0.0)
!
!     RATIONS = 0.
!     IF(TSO4.GT.0.) RATIONS = ANO3/TSO4
!
!     GAMAAN=0.0
!     IF(WH2O.GT.0.) GAMAAN = GAMA**GF1                ! activity coefficient (NH4NO3)
!     GAMAAN=min(GAMAAN,1.0)                           ! focus on 0-1 scale
!     GAMAAN=max(GAMAAN,0.0)
!
! Adjust gasfaseconcentrations(ppb) and aerosolconcentrations(microgram/m**3)
!
!hf c(i,j,k,i_nh4a)  = anh4 * mwnh4 * 1.e6
!hf c(i,j,k,i_nh3)  = gnh3 / CONVT *1.e6
!hf c(i,j,k,i_so4a) = TSO4 * mwso4 * 1.e6
!hf c(i,j,k,i_so4)  = 0.0
!hf c(i,j,k,i_no3a) = ano3 * mwno3 * 1.e6
!hf c(i,j,k,i_hno3) = gno3/ CONVT *1.e6
!hf ah2o(i,j,k) = wh2o

!hf outputted as molec/cm3, converted in My_EQSAM
     aSO4out(k)= aso4*1.e6  !hf *1.e6 since total array multiplied by 1.e-6 in the start
     aNO3out(k)= ano3*1.e6
     aNH4out(k)= anh4*1.e6
     gNH3out(k)= gnh3*1.e6
     gNO3out(k)= gno3*1.e6
     aH2Oout(k)= wh2o
enddo !k
!hf enddo
!hf enddo

 return
 end subroutine aero
!>-------------------------------------------------------------------------------<

 end module EQSAM_ml
