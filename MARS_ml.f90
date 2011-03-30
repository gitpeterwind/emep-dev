! <MARS_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2011 met.no
!* 
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*  
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!* 
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!* 
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!
 
module MARS_ml
 ! -----------------------------------------------------------------------
 ! Calculates gas-aerosol equilibrium for SO4, HNO3-NO3 and NH3-NH4 system
 ! Made available by Frank Binkowski (originally from EPA's RPM model)
 ! Presently not in use, EQSAM is used for inorganic equilibrium stuff
 !------------------------------------------------------------------------

 use Io_ml,              only : ios
 use MARS_Aero_water_ml, only:  Awater
 use ModelConstants_ml,  only : NPROC
 use Par_ml,             only : me
 implicit none
 private

  real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  


 !/- subroutines:
 public   ::  rpmares,   &
              cubic,     &
              actcof

 contains

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine rpmares ( SO4, HNO3, NO3, NH3, NH4, RH, TEMP,   &
                           ASO4, ANO3, AH2O, ANH4, GNH3, GNO3,   &
                           ERRMARK,deb) 

!-----------------------------------------------------------------------
!C
!C Description:
!C
!C   ARES calculates the chemical composition of a sulfate/nitrate/
!C   ammonium/water aerosol based on equilibrium thermodynamics.
!C
!C   This code considers two regimes depending upon the molar ratio 
!C   of ammonium to sulfate. 
!C
!C   For values of this ratio less than 2,the code solves a cubic for 
!C   hydrogen ion molality, HPLUS,  and if enough ammonium and liquid
!C   water are present calculates the dissolved nitric acid. For molal
!C   ionic strengths greater than 50, nitrate is assumed not to be present. 
!C   
!C   For values of the molar ratio of 2 or greater, all sulfate is assumed
!C   to be ammonium sulfate and a calculation is made for the presence of
!C   ammonium nitrate.
!C
!C   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!C   obtain the activity coefficients. Abandoned -7/30/97 FSB 

!c   The Bromley method of calculating the activity coefficients is s used
!c    in this version

!c   The calculation of liquid water
!C   is done in subroutine water. Details for both calculations are given
!C   in the respective subroutines.
!C
!C   Based upon MARS due to 
!C   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld, 
!C   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
!C
!C   and SCAPE due to 
!C   Kim, Seinfeld, and Saxeena, Aerosol Ceience and Technology,
!C   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
!C
!C NOTE: All concentrations supplied to this subroutine are TOTAL
!C       over gas and aerosol phases
!C
!C Parameters:
!C 
!C  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate (IN)
!C  HNO3  : Nitric Acid in MICROGRAMS/M**3 as nitric acid (IN)
!C  NO3   : Total nitrate in MICROGRAMS/M**3 as nitric acid (IN)
!C  NH3   : Total ammonia in MICROGRAMS/M**3 as ammonia (IN)
!C  NH4   : Ammonium in MICROGRAMS/M**3 as ammonium (IN)
!C  RH    : Fractional relative humidity (IN)
!C  TEMP  : Temperature in Kelvin (IN)
!C  GNO3  : Gas phase nitric acid in MICROGRAMS/M**3 (OUT)
!C  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 (OUT)
!C  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 (OUT) 
!C  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 (OUT)
!C  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 (OUT)
!C  AH2O  : Aerosol phase water in MICROGRAMS/M**3 (OUT)
!C  NITR  : Number of iterations for obtaining activity coefficients  (OUT) 
!C  NR    : Number of real roots to the cubic in the low ammonia case (OUT)
!C 
!C Revision History:
!C      Who       When        Detailed description of changes
!C   ---------   --------  -------------------------------------------
!C   S.Roselle   11/10/87  Received the first version of the MARS code
!C   S.Roselle   12/30/87  Restructured code
!C   S.Roselle   2/12/88   Made correction to compute liquid-phase 
!C                         concentration of H2O2.  
!C   S.Roselle   5/26/88   Made correction as advised by SAI, for 
!C                         computing H+ concentration.
!C   S.Roselle   3/1/89    Modified to operate with EM2
!C   S.Roselle   5/19/89   Changed the maximum ionic strength from 
!C                         100 to 20, for numerical stability.
!C   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!C                         using equations for nitrate budget.
!C   F.Binkowski 6/18/91   New ammonia poor case which 
!C                         omits letovicite.
!C   F.Binkowski 7/25/91   Rearranged entire code, restructured
!C                         ammonia poor case.
!C   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!C                         as SO4--
!C   F.Binkowski 12/6/91   Changed the ammonia defficient case so that 
!C                         there is only neutralized sulfate (ammonium
!C                         sulfate) and sulfuric acid.
!C   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement 
!C                          with the Cohen et al. (1987)  maximum molality
!C                          of 36.2 in Table III.( J. Phys Chem (91) page
!C                          4569, and Table IV p 4587.)
!C   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
!C                         possibility for denomenator becoming zero; 
!C                         this involved solving for HPLUS first.
!C                         Note that for a relative humidity
!C                          less than 50%, the model assumes that there is no 
!C                          aerosol nitrate.
!C   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)  
!C                          Redid logic as follows
!C                         1. Water algorithm now follows Spann & Richardson
!C                         2. Pitzer Multicomponent method used
!C                         3. Multicomponent practical osmotic coefficient 
!C                            use to close iterations.
!C                         4. The model now assumes that for a water
!C                            mass fraction WFRAC less than 50% there is
!C                            no aerosol nitrate.
!C   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor 
!C                         case, and changed the WFRAC criterion to 40%.
!C                         For ammonium to sulfate ratio less than 1.0 
!C                         all ammonium is aerosol and no nitrate aerosol
!C                         exists.
!C   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
!C                         allow gas-phase ammonia to exist. 
!C   F.Binkowski 7/26/95   Changed equilibrium constants to values from 
!C                         Kim et al. (1993) 
!C   F.Binkowski 6/27/96   Changed to new water format
!c   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent 
!c                         activity coefficients. The binary activity coefficients 
!c                         are the same as the previous version
!c   F.Binkowski 8/1/97    Chenged minimum sulfate from 0.0 to 1.0e-6 i.e.
!c                         1 picogram per cubic meter
!C   I.Ackermann 2/23/98   modification for solving precision problems
!C        	  on workstations
!C   I.Ackermann 2/25/99   changed logic as follows:
!c                         If iterations fail, initial values of nitrate
!c                          are retained.
!c                         Total mass budgets are changed to account for gas
!c                         phase returned. (incorporated from FZB's models3
!c                         framework)
!C                         eliminated ratio=5 !! for low to zero sulfate
!C   I.Ackermann 3/17/99   modified ratio = 5 treatment see RB3,p.125
!C
!C-----------------------------------------------------------------------


!...........ARGUMENTS and their descriptions

  real, intent(in) ::  SO4   &     ! Total sulfate in micrograms / m**3 
                      ,HNO3  &     ! Total nitric acid in micrograms / m**3
                      ,NO3   &     ! Total nitrate in micrograms / m**3
                      ,NH3   &     ! Total ammonia in micrograms / m**3
                      ,NH4   &     ! Total ammonium in micrograms / m**3
                      ,RH    &     ! Fractional relative humidity 
                      ,TEMP        ! Temperature in Kelvin 

  real, intent(out)::  ASO4  &     ! Aerosol sulfate in micrograms / m**3 
                      ,ANO3  &     ! Aerosol nitrate in micrograms / m**3
                      ,AH2O  &     ! Aerosol liquid water content water in micrograms / m**3
                      ,ANH4  &     ! Aerosol ammonium in micrograms / m**3
                      ,GNO3  &     ! Gas-phase nitric acid in micrograms / m**3
                      ,GNH3        ! Gas-phase ammonia in micrograms / m**3

  logical, intent(in) :: deb

!C...........INCLUDES and their descriptions
!!      INCLUDE SUBST_CONST          ! constants

!...........PARAMETERS and their descriptions:

      real        MWNACL           ! molecular weight for NaCl
      parameter ( MWNACL = 58.44277 )

      real        MWNO3            ! molecular weight for NO3
      parameter ( MWNO3  = 62.0049 ) 

      real        MWHNO3           ! molecular weight for HNO3
      parameter ( MWHNO3 = 63.01287 )       

      real        MWSO4            ! molecular weight for SO4
      parameter ( MWSO4 = 96.0576 )

      real        MWHSO4           ! molecular weight for HSO4
      parameter ( MWHSO4 = MWSO4 + 1.0080 ) 

      real        MH2SO4           ! molecular weight for H2SO4
      parameter ( MH2SO4 = 98.07354 ) 

      real        MWNH3            ! molecular weight for NH3
      parameter ( MWNH3 = 17.03061 ) 

      real        MWNH4            ! molecular weight for NH4
      parameter ( MWNH4 = 18.03858 )

      real        MWORG            ! molecular weight for Organic Species
      parameter ( MWORG = 16.0 )

      real        MWCL             ! molecular weight for Chloride  
      parameter ( MWCL = 35.453 )

      real        MWAIR            ! molecular weight for AIR
      parameter ( MWAIR = 28.964 )

      real        MWLCT            ! molecular weight for Letovicite
      parameter ( MWLCT = 3.0 * MWNH4 + 2.0 * MWSO4 + 1.0080 )

      real        MWAS             ! molecular weight for Ammonium Sulfate
      parameter ( MWAS = 2.0 * MWNH4 + MWSO4 )

      real        MWABS            ! molecular weight for Ammonium Bisulfate 
      parameter ( MWABS = MWNH4 + MWSO4 + 1.0080 )


!...........SCRATCH LOCAL VARIABLES and their descriptions:
       
      REAL        irh              ! Index set to percent relative humidity  
      INTEGER     NITR             ! Number of iterations for activity coefficients
      INTEGER     NNN              ! Loop index for iterations 
      INTEGER     NR               ! Number of roots to cubic equation for HPLUS
      INTEGER     ERRMARK

      real          A0             ! Coefficients and roots of 
      real          A1             ! Coefficients and roots of 
      real          A2             ! Coefficients and roots of 
      real        AA               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        BAL              ! internal variables ( high ammonia case)
      real        BB               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        BHAT             ! Variables used for ammonia solubility 
      real        CC               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        CONVT            ! Factor for conversion of units  
      real        DD               ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        DISC             ! Coefficients and discriminant for quadratic equation for ammonium nitrate
      real        EROR             ! Relative error used for convergence test  
      real        FNH3             ! "Free ammonia concentration", that which exceeds TWOSO4       
      real        GAMAAB           ! Activity Coefficient for (NH4+, HSO4-)GAMS( 2,3 )
      real        GAMAAN           ! Activity coefficient for (NH4+, NO3-) GAMS( 2,2 )
      real        GAMAHAT          ! Variables used for ammonia solubility 
      real        GAMANA           ! Activity coefficient for (H+ ,NO3-)   GAMS( 1,2 )
      real        GAMAS1           ! Activity coefficient for (2H+, SO4--) GAMS( 1,1 )
      real        GAMAS2           ! Activity coefficient for (H+, HSO4-)  GAMS( 1,3 )
      real        GAMOLD           ! used for convergence of iteration
      real        GASQD            ! internal variables ( high ammonia case)
      real        HPLUS            ! Hydrogen ion (low ammonia case) (moles / kg water)
      real        K1A              ! Equilibrium constant for ammoniua to ammonium
      real        K2SA             ! Equilibrium constant for sulfate-bisulfate (aqueous)
      real        K3               ! Dissociation constant for ammonium nitrate 
      real        KAN              ! Equilibrium constant for ammonium nitrate (aqueous)
      real        KHAT             ! Variables used for ammonia solubility 
      real        KNA              ! Equilibrium constant for nitric acid (aqueous)   
      real        KPH              ! Henry's Law Constant for ammonia       
      real        KW               ! Equilibrium constant for water dissociation             
      real        KW2              ! Internal variable using KAN 
      real        MAN              ! Nitrate (high ammonia case) (moles / kg water)
      real        MAS              ! Sulfate (high ammonia case) (moles / kg water)
      real        MHSO4            ! Bisulfate (low ammonia case) (moles / kg water)
      real        MNA              ! Nitrate (low ammonia case) (moles / kg water)
      real        MNH4             ! Ammonium (moles / kg water)
      real        MOLNU            ! Total number of moles of all ions
      real        MSO4             ! Sulfate (low ammonia case) (moles / kg water)
      real        PHIBAR           ! Practical osmotic coefficient      
      real        PHIOLD           ! Previous value of practical osmotic coefficient used for convergence of iteration
      real        RATIO            ! Molar ratio of ammonium to sulfate
      real        RK2SA            ! Internal variable using K2SA
      real        RKNA             ! Internal variables using KNA
      real        RKNWET           ! Internal variables using KNA
      real        RR1
      real        RR2
      real        STION            ! Ionic strength
      real        T1               ! Internal variables for temperature corrections
      real        T2               ! Internal variables for temperature corrections
      real        T21              ! Internal variables of convenience (low ammonia case)
      real        T221             ! Internal variables of convenience (low ammonia case)
      real        T3               ! Internal variables for temperature corrections
      real        T4               ! Internal variables for temperature corrections
      real        T6               ! Internal variables for temperature corrections
      real        TNH4             ! Total ammonia and ammonium in micromoles / meter ** 3
      real        TNO3             ! Total nitrate in micromoles / meter ** 3
      real        TOLER1           ! Tolerances for convergence test 
      real        TOLER2           ! Tolerances for convergence test 
      real        TSO4             ! Total sulfate in micromoles / meter ** 3
      real        TWOSO4           ! 2.0 * TSO4  (high ammonia case) (moles / kg water)
      real        WFRAC            ! Water mass fraction 
      real        WH2O             ! Aerosol liquid water content (internally)
                                   ! micrograms / meter **3 on output
                                   ! internally it is 10 ** (-6) kg (water) / meter ** 3
                                   ! the conversion factor (1000 g = 1 kg) is applied 
                                   ! for AH2O output
      real        WSQD             ! internal variables ( high ammonia case)
      real        XNO3             ! Nitrate aerosol concentration in micromoles / meter ** 3
      real        XXQ              ! Variable used in quadratic solution
      real        YNH4             ! Ammonium aerosol concentration in micromoles / meter** 3
      real        ZH2O             ! Water variable saved in case ionic strength too high.
      real        ZSO4             ! Total sulfate molality - mso4 + mhso4 (low ammonia case) (moles / kg water)

      real        CAT( 2 )         ! Array for cations (1, H+); (2, NH4+) (moles / kg water)
      real        AN ( 3 )         ! Array for anions (1, SO4--); (2, NO3-); (3, HSO4-)  (moles / kg water) 
      real        CRUTES( 3 )      ! Coefficients and roots of 
      real        GAMS( 2, 3 )     ! Array of activity coefficients 
      real        MINSO4           ! Minimum value of sulfate laerosol concentration
       parameter( MINSO4 = 1.0E-6 / MWSO4 ) 
      real        MINNO3
       parameter( MINNO3 = 1.0E-6 / MWNO3 )    !2/25/99 IJA
!st      real        FLOOR
!st       parameter( FLOOR = 1.0E-30) ! minimum concentration       
!2/25/99 IJA
! FSB New variables Total ammonia and nitrate mass concentrations
      real  TMASSNH3  ! Total ammonia (gas and particle)
                      !  as ammonia gas [ug /m**3]
      real  TMASSHNO3 ! Total nitrate (vapor and particle) as
                      !  nitric acid vapor [ug /m**3]

!-----------------------------------------------------------------------
!  begin body of subroutine RPMARES
                                                                         
!...convert into micromoles/m**3
 
!..iamodels3 merge NH3/NH4 , HNO3,NO3 here
      TSO4 = MAX( 0.0, SO4 / MWSO4  )
      TNO3 = MAX( 0.0, (NO3 / MWNO3 + HNO3 / MWHNO3) )
      TNH4 = MAX( 0.0, (NH3 / MWNH3 + NH4 / MWNH4)  )

!2/25/99 IJA
!      TMASSNH3  = MAX(0.0, NH3 + (MWNH3 / MWNH4)  * NH4 )
!      TMASSHNO3 = MAX(0.0, NO3 + (MWHNO3 / MWNO3) * NO3 )

      TMASSNH3  = MAX(0.0, NH3 +  NH4 )
      TMASSHNO3 = MAX(0.0, HNO3 + NO3 )
 
!...now set humidity index IRH as a percent

!st      IRH = NINT( 100.0 * RH )
         irh = RH 
!...Check for valid IRH

       irh = MAX( 0.01, IRH )
       irh = MIN( 0.99, IRH )

!...Specify the equilibrium constants at  correct
!...  temperature.  Also change units from ATM to MICROMOLE/M**3 (for KAN,
!...  KPH, and K3 )
!...  Values from Kim et al. (1993) except as noted.
 
      CONVT = 1.0 / ( 0.082 * TEMP ) 
      T6 = 0.082E-9 *  TEMP
      T1 = 298.0 / TEMP
      T2 = ALOG( T1 )
      T3 = T1 - 1.0
      T4 = 1.0 + T2 - T1
      KNA  = 2.511E+06 *  EXP(  29.17 * T3 + 16.83 * T4 ) * T6
      K1A  = 1.805E-05 *  EXP(  -1.50 * T3 + 26.92 * T4 )
      K2SA = 1.015E-02 *  EXP(   8.85 * T3 + 25.14 * T4 )
      KW   = 1.010E-14 *  EXP( -22.52 * T3 + 26.92 * T4 )
      KPH  = 57.639    *  EXP(  13.79 * T3 - 5.39  * T4 ) * T6
!!!      K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6  
      KHAT =  KPH * K1A / KW  
      KAN  =  KNA * KHAT  

!...Compute temperature dependent equilibrium constant for NH4NO3
!...  ( from Mozurkewich, 1993)

      K3 = EXP( 118.87  - 24084.0 / TEMP -  6.025  * ALOG( TEMP ) )

!...Convert to (micromoles/m**3) **2

      K3 = K3 * CONVT * CONVT

      WH2O   = 0.0
      STION  = 0.0
      AH2O   = 0.0
      MAS    = 0.0
      MAN    = 0.0
      HPLUS  = 0.0
      TOLER1 = 0.00001
      TOLER2 = 0.001
      NITR   = 0
      NR     = 0
      RATIO  = 0.0
      GAMAAN = 1.0
      GAMOLD = 1.0

!...set the ratio according to the amount of sulfate and nitrate

      IF ( TSO4 > MINSO4 ) THEN
        RATIO = TNH4 / TSO4

!...If there is no sulfate and no nitrate, there can be no ammonium
!...  under the current paradigm. Organics are ignored in this version.

      ELSE 
      
        IF ( TNO3 <= MINNO3 ) THEN

! *** If there is very little sulfate and nitrate set concentrations
!      to a very small value and return.    
          ASO4 = MAX(FLOOR, ASO4)
          ANO3 = MAX(FLOOR, ANO3 )          
          WH2O = 0.0
          AH2O = 0.0
          GNH3 = MAX(FLOOR,GNH3)
          GNO3 = MAX(FLOOR,GNO3)
          RETURN
       END IF
       
!...For the case of no sulfate and nonzero nitrate, set ratio to 5
!...  to send the code to the high ammonia case if there is more
!...  ammonia than sulfate, otherwise send to low ammonia case.

       IF (TNH4 > TSO4) THEN
         RATIO = 5.0        !this is a high ammonia case with low sulfur
       ELSE
         RATIO = 1.        !this is a low ammonia case with low sulfur
       END IF
 
      END IF 

!....................................
!......... High Ammonia Case ........
!....................................
 
      IF ( RATIO > 2.0 ) THEN
        GAMAAN = 0.1

!...Set up twice the sulfate for future use.

        TWOSO4 = 2.0 * TSO4
        XNO3 = 0.0            
        YNH4 = TWOSO4

!...Treat different regimes of relative humidity 

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  start with ammomium sulfate solution without nitrate

      CALL awater(IRH,TSO4,YNH4,TNO3,AH2O ) !**** note TNO3
        WH2O = 1.0E-3 * AH2O  
        ASO4 = TSO4 * MWSO4
        ANO3 = 0.0
        ANH4 = YNH4 * MWNH4
        WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )
!!!!       IF ( WFRAC == 0.0 )  RETURN   ! No water       
        IF ( WFRAC < 0.2 ) THEN
 
!..."dry" ammonium sulfate and ammonium nitrate
!...  compute free ammonia

          FNH3 = TNH4 - TWOSO4
          CC = TNO3 * FNH3 - K3

!...check for not enough to support aerosol      

          IF ( CC <= 0.0 ) THEN
            XNO3 = 0.0
          ELSE
            AA = 1.0
            BB = -( TNO3 + FNH3 ) 
            DISC = BB * BB - 4.0 * CC

!...Check for complex roots of the quadratic
!...  set nitrate to zero and RETURN if complex roots are found
!2/25/99 IJA

            IF ( DISC < 0.0 ) THEN
              XNO3 = 0.0
              AH2O = 1000.0 * WH2O
              YNH4 = TWOSO4
              GNO3 = HNO3
              ASO4 = TSO4 * MWSO4
              ANO3 = NO3
              ANH4 = YNH4 * MWNH4
              GNH3 = TMASSNH3 - ANH4
              RETURN
            END IF

!...to get here, BB .lt. 0.0, CC .gt. 0.0 always      

            DD = SQRT( DISC )
            XXQ = -0.5 * ( BB + SIGN ( 1.0, BB ) * DD )

!...Since both roots are positive, select smaller root.      

            XNO3 = MIN( XXQ / AA, CC / XXQ )
          
          END IF
!2/25/99 IJA
          AH2O = 1000.0 * WH2O
          YNH4 = TWOSO4 + XNO3
          ASO4 = TSO4 * MWSO4
          ANO3 = XNO3 * MWNO3
          ANH4 = YNH4 * MWNH4
          GNH3 = TMASSNH3 - ANH4
          GNO3 = TMASSHNO3 - ANO3
          RETURN

        END IF

!...liquid phase containing completely neutralized sulfate and
!...  some nitrate.  Solve for composition and quantity.
 
        MAS = TSO4 / WH2O
        MAN = 0.0
        XNO3 = 0.0
        YNH4 = TWOSO4
        PHIOLD = 1.0

!...Start loop for iteration
 
!...The assumption here is that all sulfate is ammonium sulfate,
!...  and is supersaturated at lower relative humidities.
 
        DO 1501 NNN = 1, 150 
          NITR = NNN
          GASQD = GAMAAN * GAMAAN
          WSQD = WH2O * WH2O
          KW2 = KAN * WSQD / GASQD
          AA = 1.0 - KW2
          BB = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
          CC = -KW2 * TNO3 * ( TNH4 - TWOSO4 )

!...This is a quadratic for XNO3 [MICROMOLES / M**3] of nitrate in solution.

          DISC = BB * BB - 4.0 * AA * CC

!...Check for complex roots, retain inital values and RETURN
!2/25/99 IJA

          IF ( DISC < 0.0 ) THEN
            XNO3 = 0.0
            AH2O = 1000.0 * WH2O
            YNH4 = TWOSO4
            GNO3 = HNO3
            ASO4 = TSO4 * MWSO4
            ANO3 = NO3
            ANH4 = YNH4 * MWNH4
            GNH3 = TMASSNH3 - ANH4

!!!            WRITE( 10, * ) ' COMPLEX ROOTS '
            RETURN
          END IF
! 2/25/99 IJA

! Deal with degenerate case (yoj)

          IF ( AA /= 0.0 ) THEN
            DD = SQRT( DISC )
            XXQ = -0.5 * ( BB + SIGN ( 1.0, BB ) * DD )
            RR1 = XXQ / AA
            RR2 = CC / XXQ

!...choose minimum positve root         

            IF ( ( RR1 * RR2 ) < 0.0 ) THEN
              XNO3 = MAX( RR1, RR2 )
            ELSE 
              XNO3 = MIN( RR1, RR2 )
            END IF
          ELSE
             XNO3 = - CC / BB
          END IF

          XNO3 = MIN( XNO3, TNO3 )

!...This version assumes no solid sulfate forms (supersaturated ) 
!...  Now update water

          CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O)

!...ZSR relationship is used to set water levels. Units are
!...  10**(-6) kg water/ (cubic meter of air)
!...  The conversion from micromoles to moles is done by the units of WH2O.

          WH2O = 1.0E-3 * AH2O

!...Ionic balance determines the ammonium in solution.

          MAN = XNO3 / WH2O
          MAS = TSO4 / WH2O
          MNH4 = 2.0 * MAS + MAN
          YNH4 = MNH4 * WH2O

 !st ...  FIXING
   if(MNH4<0. .or. MAS<0. .or. MAN<0.) then
      MNH4 = 1.e-30
      MAS  = 1.e-30
      MAN  = 1.e-30
   endif

!...MAS, MAN and MNH4 are the aqueous concentrations of sulfate, nitrate,
!...  and ammonium in molal units (moles/(kg water) ).

          STION = 3.0 * MAS + MAN
          CAT( 1 ) = 0.0
          CAT( 2 ) = MNH4 
          AN ( 1 ) = MAS
          AN ( 2 ) = MAN
          AN ( 3 ) = 0.0
          CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR , ERRMARK,1,deb)
          GAMAAN = GAMS( 2, 2 )

!...Use GAMAAN for convergence control

          EROR = ABS( GAMOLD - GAMAAN ) / GAMOLD
          GAMOLD = GAMAAN

!...Check to see if we have a solution

          IF ( EROR <= TOLER1 ) THEN 
!!!            WRITE( 11, * ) RH, STION, GAMS( 1, 1 ),GAMS( 1, 2 ), GAMS( 1, 3 ),
!!!     &      GAMS( 2, 1 ), GAMS( 2, 2 ), GAMS( 2, 3 ), PHIBAR
! 2/25/99 IJA
            ASO4 = TSO4 * MWSO4
            ANO3 = XNO3 * MWNO3
            ANH4 = YNH4 * MWNH4
            GNO3 = TMASSHNO3  - ANO3
            GNH3 = TMASSNH3   - ANH4
            AH2O = 1000.0 * WH2O
            RETURN
          END IF

1501    CONTINUE

!...If after NITR iterations no solution is found, then:
! FSB retain the initial values of nitrate particle and vapor
! 2/25/99 IJA
        ASO4 = TSO4 * MWSO4
        ANO3 = NO3
        XNO3 = NO3 / MWNO3
        YNH4 = TWOSO4
        ANH4 = YNH4 * MWNH4
        CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O)
        GNO3 = HNO3
        GNH3 = TMASSNH3 - ANH4
        RETURN

      ELSE
       
!......................................
!......... Low Ammonia Case ...........
!......................................
      
!...coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)
      ! modified 8/28/98

!...All cases covered by this logic
 
        WH2O = 0.0
        CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
        WH2O = 1.0E-3 * AH2O
        ZH2O = AH2O

!...convert 10**(-6) kg water/(cubic meter of air) to micrograms of water
!...  per cubic meter of air (1000 g = 1 kg)
! 2/25/99 IJA 
        ASO4 = TSO4 * MWSO4
        ANH4 = TNH4 * MWNH4
        ANO3 = NO3
        GNO3 = TMASSHNO3 - ANO3
        GNH3 = FLOOR

!...Check for zero water.      

        IF ( WH2O == 0.0 ) RETURN
        ZSO4 = TSO4 / WH2O 

!...ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4      

!!!         IF ( ZSO4 > 11.0 ) THEN

!...do not solve for aerosol nitrate for total sulfate molality
!...  greater than 11.0 because the model parameters break down
!...  greater than  9.0 because the model parameters break down

        IF ( ZSO4 > 9.0 ) THEN   ! 18 June 97
          RETURN
        END IF

!...First solve with activity coeffs of 1.0, then iterate.

        PHIOLD = 1.0
        GAMANA = 1.0
        GAMAS1 = 1.0
        GAMAS2 = 1.0
        GAMAAB = 1.0
        GAMOLD = 1.0

!...All ammonia is considered to be aerosol ammonium. 

        MNH4 = TNH4 / WH2O

!...MNH4 is the molality of ammonium ion.

        YNH4 = TNH4
      
!...loop for iteration
 
        DO 1601 NNN = 1, 150
          NITR = NNN

!...set up equilibrium constants including activities
!...  solve the system for hplus first then sulfate & nitrate

          RK2SA = K2SA * GAMAS2 * GAMAS2 / ( GAMAS1 * GAMAS1 * GAMAS1 )
          RKNA = KNA / ( GAMANA * GAMANA )
          RKNWET = RKNA * WH2O       
          T21  = ZSO4 - MNH4
          T221 = ZSO4 + T21

!...set up coefficients for cubic       

          A2 = RK2SA + RKNWET - T21
          A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET )    &
             - RK2SA * ZSO4 - RKNA * TNO3
          A0 = - (T21 * RK2SA * RKNWET                      &
             + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )   
         
          CALL CUBIC ( A2, A1, A0, NR, CRUTES )
       
!...Code assumes the smallest positive root is in CRUTES(1)
 
          HPLUS = CRUTES( 1 )
          BAL = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0
          MSO4 = RK2SA * ZSO4 / ( HPLUS + RK2SA )   ! molality of sulfate ion
          MHSO4 = ZSO4 - MSO4                       ! molality of bisulfate ion
          MNA = RKNA * TNO3 / ( HPLUS + RKNWET )    ! molality of nitrate ion
          MNA = MAX( 0.0, MNA )
          MNA = MIN( MNA, TNO3 / WH2O )
          XNO3 = MNA * WH2O
          ANO3 = MNA * WH2O * MWNO3
! 2/25/99 IJA
          GNO3 = TMASSHNO3 - ANO3
        
!...Calculate ionic strength      

          STION = 0.5 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0 * MSO4 )
          
!...Update water

          CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

!...Convert 10**(-6) kg water/(cubic meter of air) to micrograms of water
!...  per cubic meter of air (1000 g = 1 kg)                       

          WH2O = 1.0E-3 * AH2O 
          CAT( 1 ) = HPLUS
          CAT( 2 ) = MNH4
          AN ( 1 ) = MSO4
          AN ( 2 ) = MNA
          AN ( 3 ) = MHSO4

          CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR, ERRMARK,2,deb)

          GAMANA = GAMS( 1, 2 )
          GAMAS1 = GAMS( 1, 1 )
          GAMAS2 = GAMS( 1, 3 )
          GAMAAN = GAMS( 2, 2 )

          GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
          BHAT = KHAT * GAMAHAT 
!!!          EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
!!!          PHIOLD = PHIBAR
          EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD 
          GAMOLD = GAMAHAT   

!...write out molalities and activity coefficient
!...  and return with good solution

          IF ( EROR <= TOLER2 ) THEN
!!!            WRITE(12,*) RH, STION,HPLUS,ZSO4,MSO4,MHSO4,MNH4,MNA
!!!            WRITE(11,*) RH, STION, GAMS(1,1),GAMS(1,2),GAMS(1,3),
!!!     &                  GAMS(2,1),GAMS(2,2),GAMS(2,3), PHIBAR
            RETURN
          END IF

1601    CONTINUE     

!...after NITR iterations, failure to solve the system, no ANO3
! 2/25/99 IJA
        ANH4 = TNH4 * MWNH4
        GNH3 = FLOOR
        GNO3 = HNO3
        ANO3 = NO3
        CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )      
        RETURN
            
      END IF   ! ratio .gt. 2.0

      end subroutine rpmares ! end RPMares

!>-------------------------------------------------------------------------------<
!<------------------------------------------------------------------------------->

      subroutine cubic(a2,a1,a0,nr,crutes)

  !.. subroutine  to find the roots of a cubic equation / 3rd order polynomial
  !.. formulae can be found in numer. recip.  on page 145
  !..  kiran  developed  this version on 25/4/1990
  !..  Dr. Francis Binkowski modified the routine on 6/24/91, 8/7/97
  !--------------------------------------------------------------

      implicit none

      real, intent(in)     :: a2,a1,a0
      integer, intent(out) :: nr      
      real, intent(out)    :: crutes(3)
!.. local
      real ::  qq,rr,a2sq,theta, sqrt3, one3rd
      real ::  dum1,dum2,part1,part2,part3,rrsq,phi,yy1,yy2,yy3
      real ::  costh, sinth

      data sqrt3/1.732050808/, one3rd/0.333333333/
!=======

      a2sq=a2*a2
      qq=(a2sq-3.*a1)/9.
      rr=( a2*(2.*a2sq - 9.*a1) + 27.*a0 )/54.
! CASE 1 THREE real ROOTS or  CASE 2 ONLY ONE real ROOT
      dum1=qq*qq*qq 
      rrsq=rr*rr
      dum2=dum1 - rrsq
      
      if(dum2 >= 0.) then
! NOW WE HAVE THREE real ROOTS
         phi=sqrt(dum1)
         if(abs(phi) <= 1.e-20) then 
!           write(10,*) ' cubic phi small, phi = ',phi
            crutes(1) = 0.0 
            crutes(2) = 0.0
            crutes(3) = 0.0
            nr = 0            
           stop 
         end if
         theta=acos(rr/phi)/3.0
         costh = cos(theta)
         sinth = sin(theta)
! *** use trig identities to simplify the expressions 
! *** binkowski's modification
         part1=sqrt(qq)
         yy1=part1*costh
         yy2=yy1-a2/3.0
         yy3=sqrt3*part1*sinth
         crutes(3) = -2.0*yy1 - a2/3.0
         crutes(2) = yy2 + yy3
         crutes(1) = yy2 - yy3
! *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
         if(crutes(1) <= 0.0) crutes(1) = 1.0e9
         if(crutes(2) <= 0.0) crutes(2) =1.0e9
         if(crutes(3) <= 0.0) crutes(3) = 1.0e9
! *** put smallest positive root in crutes(1)
         crutes(1)=min( crutes(1),crutes(2),crutes(3))
         nr=3
      else  ! dum IS NEGATIVE
!     NOW HERE WE HAVE ONLY ONE real ROOT
         part1=sqrt(rrsq-dum1)
         part2=abs(rr)
         part3=(part1+part2)**one3rd
         crutes(1) = -sign(1.0,rr) * ( part3 + (qq/part3) ) - a2/3. 
         crutes(2)=0.
         crutes(3)=0.
!IAREV02...ADDITIONAL CHECK on NEGATIVE ROOTS
! *** SET NEGATIVE ROOTS TO A LARGE POSITIVE VALUE
!IA ACTIONIA
      if(crutes(1) <= 0.0) THEN
         crutes(1) = 1.0e9
  !!    if(deb ) write(6,*) 'WARNING: NEGATIVE ROOTS IN CUBIC', crutes(1)
  !!st       stop
      end if
      nr=1
      end if

   end subroutine cubic
    
!>-------------------------------------------------------------------------------<
!<------------------------------------------------------------------------------->

      subroutine actcof ( CAT, AN, GAMA, MOLNU, PHIMULT , ERRMARK, IA2, deb)

!C-----------------------------------------------------------------------
!C
!C DESCRIPTION:
!C
!C  This subroutine computes the activity coefficients of (2NH4+,SO4--),
!C  (NH4+,NO3-),(2H+,SO4--),(H+,NO3-),AND (H+,HSO4-) in aqueous
!C  multicomponent solution, using Bromley's model and Pitzer's method.
!C
!C REFERENCES:
!C
!C   Bromley, L.A. (1973) Thermodynamic properties of strong electrolytes
!C     in aqueous solutions.  AIChE J. 19, 313-320.
!C
!C   Chan, C.K. R.C. Flagen, & J.H.  Seinfeld (1992) Water Activities of 
!C     NH4NO3 / (NH4)2SO4 solutions, Atmos. Environ. (26A): 1661-1673.
!C
!C   Clegg, S.L. & P. Brimblecombe (1988) Equilibrium partial pressures 
!C     of strong acids over saline solutions - I HNO3, 
!C     Atmos. Environ. (22): 91-100
!C
!C   Clegg, S.L. & P. Brimblecombe (1990) Equilibrium partial pressures
!C     and mean activity and osmotic coefficients of 0-100% nitric acid 
!C     as a function of temperature,   J. Phys. Chem (94): 5369 - 5380
!C
!C   Pilinis, C. and J.H. Seinfeld (1987) Continued development of a
!C     general equilibrium model for inorganic multicomponent atmospheric
!C     aerosols.  Atmos. Environ. 21(11), 2453-2466.
!C


!
!CC ARGUMENT DESCRIPTION:
!
!C     CAT(1) : conc. of H+    (moles/kg)
!C     CAT(2) : conc. of NH4+  (moles/kg)
!C     AN(1)  : conc. of SO4-- (moles/kg)
!C     AN(2)  : conc. of NO3-  (moles/kg)
!C     AN(3)  : conc. of HSO4- (moles/kg)
!C     GAMA(2,1)    : mean molal ionic activity coeff for (2NH4+,SO4--)
!C     GAMA(2,2)    :  "    "     "       "       "    "  (NH4+,NO3-)
!C     GAMA(2,3)    :  "    "     "       "       "    "  (NH4+. HSO4-)
!C     GAMA(1,1)    :  "    "     "       "       "    "  (2H+,SO4--)
!C     GAMA(1,2)    :  "    "     "       "       "    "  (H+,NO3-)
!C     GAMA(1,3)    :  "    "     "       "       "    "  (H+,HSO4-)
!C     MOLNU   : the total number of moles of all ions.
!C     PHIMULT : the multicomponent paractical osmotic coefficient.
!C
!C REVISION HISTORY:
!C      Who       When        Detailed description of changes
!C   ---------   --------  -------------------------------------------
!C   S.Roselle   7/26/89   Copied parts of routine BROMLY, and began this
!C                         new routine using a method described by Pilinis
!C                         and Seinfeld 1987, Atmos. Envirn. 21 pp2453-2466.
!C   S.Roselle   7/30/97   Modified for use in Models-3
!C   F.Binkowski 8/7/97    Modified coefficients BETA0, BETA1, CGAMA
!C
!-----------------------------------------------------------------------




!...........INCLUDES and their descriptions

!      INCLUDE SUBST_XSTAT     ! M3EXIT status codes

!...........ARGUMENTS and their descriptions

            
     real, intent(in)  :: cat(2)    &  ! cation conc in moles/kg (input)
                         ,an (3)    &   ! anion conc in moles/kg (input)
                         ,molnu     &  ! tot # moles of all ions
                         ,phimult      ! multicomponent paractical osmotic coef
     real, intent(out) :: gama(2,3)   ! mean molal ionic activity coefs
     logical, intent(in) :: deb

!....................................................................

      INTEGER    XSTAT0       ! Normal, successful completion
      PARAMETER (XSTAT0 = 0)
      INTEGER    XSTAT1       ! File I/O error
      PARAMETER (XSTAT1 = 1)
      INTEGER    XSTAT2       ! Execution error
      PARAMETER (XSTAT2 = 2)
      INTEGER    XSTAT3       ! Special  error
      PARAMETER (XSTAT3 = 3)
      INTEGER ERRMARK
      INTEGER IA2
      CHARACTER*120 XMSG

!...........PARAMETERS and their descriptions:

      INTEGER      NCAT                 ! number of cations
      PARAMETER  ( NCAT = 2 )

      INTEGER      NAN                  ! number of anions
      PARAMETER  ( NAN = 3 )


!...........SCRATCH LOCAL VARIABLES and their descriptions:

      CHARACTER*16 PNAME            ! driver program name
      SAVE         PNAME

      INTEGER      IAN                  ! anion indX
      INTEGER      ICAT                 ! cation indX

      REAL         FGAMA                ! 
      REAL         I                    ! ionic strength 
      REAL         R                    ! 
      REAL         S                    ! 
      REAL         TA                   ! 
      REAL         TB                   ! 
      REAL         TC                   ! 
      REAL         TEXPV                ! 
      REAL         TRM                  ! 
      REAL         TWOI                 ! 2*ionic strength
      REAL         TWOSRI               ! 2*sqrt of ionic strength
      REAL         ZBAR                 ! 
      REAL         ZBAR2                ! 
      REAL         ZOT1                 ! 
      REAL         SRI                  ! square root of ionic strength 
      REAL         F2( NCAT )           ! 
      REAL         F1( NAN )            ! 
      REAL         ZP( NCAT )           ! absolute value of charges of cation
      REAL         ZM( NAN )            ! absolute value of charges of anion
      REAL         BGAMA ( NCAT, NAN )  ! 
      REAL         X     ( NCAT, NAN )  ! 
      REAL         M     ( NCAT, NAN )  ! molality of each electrolyte
      REAL         LGAMA0( NCAT, NAN )  ! binary activity coefficients
      REAL         Y     ( NAN, NCAT )  ! 
      REAL         BETA0 ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         BETA1 ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         CGAMA ( NCAT, NAN )  ! binary activity coefficient parameter
      REAL         V1    ( NCAT, NAN )  ! number of cations in electrolyte formula
      REAL         V2    ( NCAT, NAN )  ! number of anions in electrolyte formula

      DATA         ZP / 1.0, 1.0 /
      DATA         ZM / 2.0, 1.0, 1.0 /
      DATA         XMSG / ' ' /
      DATA         PNAME / 'ACTCOF' /

! *** Sources for the coefficients BETA0, BETA1, CGAMA:
 
! *** (1,1);(1,3)  - Clegg & Brimblecombe (1988)
! *** (2,3)        - Pilinis & Seinfeld (1987), cgama different 
! *** (1,2)        - Clegg & Brimblecombe (1990)
! *** (2,1);(2,2)  - Chan, Flagen & Seinfeld (1992)                                 
      
! *** now set the basic constants, BETA0, BETA1, CGAMA 

      DATA BETA0(1,1) /2.98E-2/,    BETA1(1,1) / 0.0/,          &
           CGAMA(1,1) / 4.38E-2/                                ! 2H+SO4-
     
      DATA BETA0(1,2) / 1.2556E-1/,   BETA1(1,2) / 2.8778E-1/,  & 
           CGAMA(1,2) / -5.59E-3/                               ! HNO3
     
      DATA BETA0(1,3) / 2.0651E-1/,   BETA1(1,3) / 5.556E-1/,   &
           CGAMA(1,3) /0.0/                                     ! H+HSO4-
     
      DATA BETA0(2,1) /4.6465E-2/,   BETA1(2,1) /-0.54196/,     &    
           CGAMA(2,1) /-1.2683E-3/                              ! (NH4)2SO4
     
      DATA BETA0(2,2) /-7.26224E-3/, BETA1(2,2) /-1.168858/,    &
           CGAMA(2,2) /3.51217E-5/                              ! NH4NO3
     
      DATA BETA0(2,3) / 4.494E-2/,    BETA1(2,3) / 2.3594E-1/,  &
           CGAMA(2,3) /-2.962E-3/                               ! NH4HSO4

      DATA V1(1,1), V2(1,1) / 2.0, 1.0 /     ! 2H+SO4-
      DATA V1(2,1), V2(2,1) / 2.0, 1.0 /     ! (NH4)2SO4
      DATA V1(1,2), V2(1,2) / 1.0, 1.0 /     ! HNO3 
      DATA V1(2,2), V2(2,2) / 1.0, 1.0 /     ! NH4NO3
      DATA V1(1,3), V2(1,3) / 1.0, 1.0 /     ! H+HSO4-
      DATA V1(2,3), V2(2,3) / 1.0, 1.0 /     ! NH4HSO4

!-----------------------------------------------------------------------
!  begin body of subroutine ACTCOF

!...compute ionic strength
      I = 0.0

      DO ICAT = 1, NCAT
        I = I + CAT( ICAT ) * ZP( ICAT ) * ZP( ICAT )
      END DO

      DO IAN = 1, NAN
        I = I + AN( IAN ) * ZM( IAN ) * ZM( IAN )
      END DO

      I = 0.5 * I

!...check for problems in the ionic strength

      IF ( I .EQ. 0.0 ) THEN

        DO IAN = 1, NAN
          DO ICAT = 1, NCAT
            GAMA( ICAT, IAN ) = 0.0
          END DO
        END DO
        
        XMSG = 'Ionic strength is zero...returning zero activities'
       if(deb ) WRITE(6,*) XMSG 
        RETURN

      ELSE IF ( I .LT. 0.0 ) THEN
        XMSG = 'Ionic strength below zero...negative concentrations'
   if(deb ) then
        WRITE(6,*) XMSG
        WRITE(6,*) 'called over ', IA2
        WRITE(6,*) ' I =', I
        WRITE(6,*) 'CAT=', CAT
        WRITE(6,*) 'AN=', AN
        WRITE(6,*) 'GAMA=', GAMA
        WRITE(6,*) 'MOLNU=',MOLNU
        WRITE(6,*) 'PHIMULT=',PHIMULT
   endif
 !!       CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT2 )
 !emep1.2      call stop_test(.true.,me,NPROC,ios,'##MARS-negat.con')
     END IF

!...compute some essential expressions

      SRI    = SQRT( I )
      TWOSRI = 2.0 * SRI
      TWOI   = 2.0 * I
      TEXPV  = 1.0 - EXP( -TWOSRI ) * ( 1.0 + TWOSRI - TWOI )
      R      = 1.0 + 0.75 * I
      S      = 1.0 + 1.5  * I
      ZOT1   = 0.511 * SRI / ( 1.0 + SRI )

!...Compute binary activity coeffs

      FGAMA = -0.392 * ( ( SRI / ( 1.0 + 1.2 * SRI )         &
            + ( 2.0 / 1.2 ) * ALOG( 1.0 + 1.2 * SRI ) ) )

      DO ICAT = 1, NCAT
        DO IAN = 1, NAN

         BGAMA( ICAT, IAN ) = 2.0 * BETA0( ICAT, IAN )        &       
                            + ( 2.0 * BETA1( ICAT, IAN ) / ( 4.0 * I ) )   &
                            * TEXPV

!...compute the molality of each electrolyte for given ionic strength

          M( ICAT, IAN ) = ( CAT( ICAT )**V1( ICAT, IAN )       &
                         *   AN( IAN )**V2( ICAT, IAN ) )**( 1.0     &
                         / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) ) )   

!...calculate the binary activity coefficients

         LGAMA0( ICAT, IAN ) = ( ZP( ICAT ) * ZM( IAN ) * FGAMA       &
                             + M( ICAT, IAN )                         &
                             * ( 2.0 * V1( ICAT, IAN ) * V2( ICAT, IAN )   &
                             / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )       &
                             * BGAMA( ICAT, IAN ) )                    &
                             + M( ICAT, IAN ) * M( ICAT, IAN )         &
                             * ( 2.0 * ( V1( ICAT, IAN )               &
                             * V2( ICAT, IAN ) )**1.5                  &
                             / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )    &
                             * CGAMA( ICAT, IAN ) ) ) / 2.302585093

        END DO
      END DO

!...prepare variables for computing the multicomponent activity coeffs

      DO IAN = 1, NAN
        DO ICAT = 1, NCAT
          ZBAR = ( ZP( ICAT ) + ZM( IAN ) ) * 0.5
          ZBAR2 = ZBAR * ZBAR
          Y( IAN, ICAT ) = ZBAR2 * AN( IAN ) / I
          X( ICAT, IAN ) = ZBAR2 * CAT( ICAT ) / I
        END DO
      END DO

      DO IAN = 1, NAN
        F1( IAN ) = 0.0
        DO ICAT = 1, NCAT
          F1( IAN ) = F1( IAN ) + X( ICAT, IAN ) * LGAMA0( ICAT, IAN )    &
                    + ZOT1 * ZP( ICAT ) * ZM( IAN ) * X( ICAT, IAN )
        END DO
      END DO

      DO ICAT = 1, NCAT
        F2( ICAT ) = 0.0
        DO IAN = 1, NAN
          F2( ICAT ) = F2( ICAT ) + Y( IAN, ICAT ) * LGAMA0( ICAT, IAN )    &
                     + ZOT1 * ZP( ICAT ) * ZM( IAN ) * Y( IAN, ICAT )
        END DO
      END DO

!...now calculate the multicomponent activity coefficients

      DO IAN = 1, NAN
        DO ICAT = 1, NCAT

          TA = -ZOT1 * ZP( ICAT ) * ZM( IAN )
          TB = ZP( ICAT ) * ZM( IAN ) / ( ZP( ICAT ) + ZM( IAN ) )
          TC = ( F2( ICAT ) / ZP( ICAT ) + F1( IAN ) / ZM( IAN ) )
          TRM = TA + TB * TC

          IF ( TRM > 30.0 ) THEN
            GAMA( ICAT, IAN ) = 1.0E+30
            XMSG = 'Multicomponent activity coefficient is >>'
       !!     if(deb )  WRITE(6,*) XMSG, gama(icat,ian)
            ERRMARK=2
 
          ELSE
            GAMA( ICAT, IAN ) = 10.0**TRM
          END IF

        END DO
      END DO
 
    end subroutine actcof 
!>-------------------------------------------------------------------------------<

 end module MARS_ml

