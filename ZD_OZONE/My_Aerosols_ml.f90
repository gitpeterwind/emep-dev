!File:  My_Aerosols_ml
! Contains modules relating to aerosols.
! These modules may be switched and altered as needed. 
! ulic parameters ORGANIC_AEROSOLS AND INORGANIC_AEROSOLS are
! assigned 

   !--------------- Inorganic and Organic  aerosols
   ! .. 1) Ammonium_ml    - SO4, AMSU, AMNI, HNO3 equilib - EMEP traditional
   ! .. 2) Organic_ml     - SOA formation    (here=dummy)


module Ammonium_ml
 !----------------------------------------------------------------------------
 ! Module to set up and process the NH4-NO3-SO4 reaction system
 !
 ! Usage:
 !   "call ammonium()"  - from Runchem
 !     - on the first call this runs the tabulation routines. On all calls
 !     the equilibrium relationships are calculated and run to establish
 !     new values of ammonium sulphate (AMSU), NH3, HNO3, SO4 and 
 !     ammonium nitrate (AMNI).
 !----------------------------------------------------------------------------
 !
 use ModelConstants_ml   , only : CHEMTMIN, CHEMTMAX   &! Temp. range
                                 , PPB             &! unit factors
                                 , KCHEMTOP            &! k=2 - top of chemistry
                                 , KMAX_MID                ! K=20 at ground
 implicit none
 private


 !/- subroutines:
 public   :: ammonium          ! Sets up most tables

 private  :: tabulate         ! Sets up most tables, and calls tab_rct_rates
 private  :: setup_ammonium   ! setup data for 1d column calculation
 private  :: calc_ammonium    ! Equilibrium distribution of NH4-SO4-NO3

 !/Wanted?
 logical, public, parameter :: INORGANIC_AEROSOLS = .true.

 !/- Outputs: - updated xn_2d concentrations after equilibrium


 !/-- Local:

  real, private, dimension(CHEMTMIN:CHEMTMAX), save  :: &
                   tab_rhdel   &  ! RH of deliquescence for ammonium nitrate
                  ,tab_Kp_amni &  ! Equil. constant, nh3 + hno3  <--> nh4no3
                  ,tab_MozP1   &  ! Mozurkewich P1 value for Kaq
                  ,tab_MozP2   &  ! Mozurkewich P2 value for Kaq
                  ,tab_MozP3   ! &  ! Mozurkewich P3 value for Kaq
!                 ,tab_vav_n2o5   ! avg. molecular speed N2O5 
				! Might move elsewhere

   logical, private, save :: my_first_call = .true.


 contains

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine ammonium()
!    integer :: k ! for print-out only

   real, dimension(KCHEMTOP:KMAX_MID)  ::  rcnh4 ! equilib. value
                                                      !was :  miscrc(ICRCNH3,k)

     if ( my_first_call ) then
        call tabulate()
     	my_first_call = .false.
     endif

     call setup_ammonium(rcnh4)
     call calc_ammonium(rcnh4)



 end subroutine ammonium

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine tabulate()
 !
     integer :: i
   real, dimension(CHEMTMIN:CHEMTMAX) :: t  ! temp.(K) for tabulations


     !/-- current temperature range: from 148 K (-125C) ro 333K (+60C):

      t = (/ (real(i),i=CHEMTMIN,CHEMTMAX) /)


    ! Tabulations tab_rhedl, tab_Kp_amni, tab_MozP.., tab_vav_n2o5
    !-------------------------------------------------------------------
    !   relative humidity of deliquescence for ammonium nitrate
    !   Ref:  Mozurkewich (1993)  - Journal???
    !   Units : fraction 0-1
    !   (MADE/MACHO notes  : was miscrcit(ICRHD,it)

       tab_rhdel(:) = exp( 618.3/t(:) - 2.551 )

    !-------------------------------------------------------------------
    !    Equilibrium constant (Kp):  NH3 + HNO3  <-------> NH4NO3   
    !    Ref: Mozurkewich (1993)
    !    Units : (molecule/cm3)^2 for Kp
    !   (MADE/MACHO notes  : was miscrcit(ICRS,it)
    !
    !      lnKp = 118.87 - 24084.0/T - 6.025* ln(T)
    !
    ! st: documentation has + 24084!
    ! c.f. Seinfeld, eqn 9.91, p.532 - suggests minus, if it is relevant?

       tab_Kp_amni(:) = exp( 118.87 - 24084.0/t(:)-6.025*alog(t(:)) )

    !-------------------------------------------------------------------
    !    temp. dependant constrants for calcolating dissos. rate 
    !    for  the formation of ammonium nitrate  
    !    Ref: Mozurkewich (1993)
    !   (MADE/MACHO notes  : was miscrcit(ICXK1,it)..miscrcit(ICXK_3,it)
    !    n.b. EMEP report 2/98 had 2446 in P3, but 24.46 is correct

       tab_MozP1(:) = exp( -135.94 +  8763.0/t(:) + 19.12*alog( t(:) ) )
       tab_MozP2(:) = exp( -122.65 +  9969.0/t(:) + 16.22*alog( t(:) ) )
       tab_MozP3(:) = exp( -182.61 + 13875.0/t(:) + 24.46*alog( t(:) ) )


  !-------------------------------------------------------------------
  end subroutine tabulate

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine setup_ammonium(rcnh4)
   !
   ! Calculates the equilibrium constant for the ammonium-suphate
   ! ammonium nitrate, Kp and Kaq (here denoted rcKaq). 
   ! Ref: EMEP Report 2/98, pages B:3, Mozurkewich (1993)
   !
   !      Kpaq = [ P1 -P2(1-rh/100) + P3(1-rh/100)^2 ] .(1-rh/100)**1.75. Kp
   !
   ! Units :  Kp, Kaq : (molecules cm-3)^2
   !          rc      ????
   ! MADE/MACHO notes.. ds- replaced xk by Kp, miscrc(ICRCNH3) by rc...??
   !--------------------------------------------------------------------------
 use Setup_1dfields_ml   , only : rh, amk, itemp

   real, dimension(KCHEMTOP:KMAX_MID)  ::  rcnh4 ! equilib. value
   real, dimension(KCHEMTOP:KMAX_MID) ::   rhd, Kp     ! deliq. rh, Kp
    real, dimension(KCHEMTOP:KMAX_MID) ::    &
               roappm                     &  ! density in ppm?
              ,humd,humdsqrt,humdsqrt2       ! humd = 1-rh
   !!          rhd, Kp                    &  ! deliq. rh, Kp

      rhd(:) =  tab_rhdel( itemp(:) )   ! was: miscrcit(ICRHD,itk)

      Kp(:)  =  tab_Kp_amni( itemp(:) ) ! was: miscrcit(ICRS,itk)

!hf Initialize rcnh4 to tab_Kp_amni,need roappm
      roappm(:) = amk(:)*PPB
      rcnh4(:)  =  tab_Kp_amni( itemp(:) )*roappm(:)* roappm(:)

!  The lines below are a CPU-efficient way of calculating the
!  power of 1.75  for  Mozurkewich Kp, suggested by  su.

      where ( rh >= rhd ) !  old: if(rh(k) >= rhd) then

        humd = 1.0001 - rh                  ! ds why not 1.0?
        humdsqrt = sqrt(humd)
        humdsqrt2 = sqrt(humdsqrt)*humdsqrt
        Kp = (   tab_MozP1(itemp) &
               - tab_MozP2(itemp)*humd &
               + tab_MozP3(itemp)*humd*humd  ) *humd*humdsqrt2*Kp

        roappm = amk*PPB
        rcnh4  = Kp * roappm * roappm   ! old misrc(ICRCNH3,k)

!hf BUG      elsewhere

!hf BUG        rcnh4 = 0.0
      end where


 end subroutine setup_ammonium
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine calc_ammonium(rcnh4)
     !** DESCRIPTION
     !   Calculates the distribution of NH3, NH4SO4, NH4NO3
     !   - needs more text...
     !-------------------------------------------------------------------------
 use GenSpec_tot_ml      , only : SO4, AMSU, AMNI, NH3, HNO3
 use Setup_1dfields_ml   , only :  xn => xn_2d

   real, dimension(KCHEMTOP:KMAX_MID)  ::  rcnh4 ! equilib. value
     real, dimension(KCHEMTOP:KMAX_MID) :: eqnh3, delteq   !ds, delt


     where ( 1.5*xn(SO4,:) >  xn(NH3,:) ) ! ammonium sulphate in excess of NH3

            xn(AMSU,:) = xn(AMSU,:) +  xn(NH3,:)*2./3.

            xn(SO4,:) = xn(SO4,:)    - xn(NH3,:)*2./3.

            xn(NH3,:) = 0.

     elsewhere

            xn(AMSU,:) = xn(AMSU,:) + xn(SO4,:)

            xn(NH3,:) = xn(NH3,:)   - xn(SO4,:)*1.5

            xn(SO4,:) = 0.

     end where

     ! The equilibrium concentration of NH3 is:

     eqnh3 = (xn(NH3,:) - xn(HNO3,:))*0.5   & 
                + sqrt( 0.25*(xn(NH3,:) -xn(HNO3,:))**2 + rcnh4 )+1.
                                                      !ds - why +1 here?
    !hf eqnh3 er i størrelsesorden 10^20.

     delteq     = eqnh3 - xn(NH3,:)
     delteq     = min(delteq,xn(AMNI,:))   ! ds - used to have delt here

     xn(AMNI,:) = xn(AMNI,:) - delteq
     xn(NH3,:)  = xn(NH3,:)  + delteq
     xn(HNO3,:) = xn(HNO3,:) + delteq

   end subroutine calc_ammonium
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module Ammonium_ml

!=============================================================================
!=============================================================================
!=====    Module 2: Organics =================================================
!=============================================================================
!=============================================================================

module OrganicAerosol_ml
! From BOX3D on 9/8/01
   ! Calculates the amount of condesible species in the gas and aerosol
   ! phases. 
   ! Methodology from Andersson-Sk\"old and Simpson, 2000, Secondary Organic
   ! Aerosol Formation in Northern Europe: a Model Study, to be published
   ! in JGR.
   !
   ! Usage: call OrganicAerosol from Runchem, after setup of 1d-fields
   ! finished.  The subroutine initialises itself on the first call
   ! and thereafter modifies two external variables:
   !   xn(SOA,k) : the concentrations of SOA 
   !   Fgas(X,k) : The fraction of X which is gas and not aeorosl
   !
   ! Dave Simpson, August 2001 
   !--------------------------------------------------------------------------
   ! nb - we use all of GenSpec_tot_ml since different model versions
   !  will have different species names. This module is intended to
   !  insensitive to the actual names one day, so this should be
   !  revised .. one day   - ds.
   !--------------------------------------------------------------------------
   use ModelConstants_ml,    only : CHEMTMIN, CHEMTMAX, &
                                    K2 => KMAX_MID, K1 => KCHEMTOP
   use PhysicalConstants_ml, only : AVOG
   use Setup_1dfields_ml,    only : itemp, xn => xn_2d
   use GenChemicals_ml,      only : species   ! for molwts
   use GenSpec_tot_ml,  A1 => FIRST_SOA , A2 => LAST_SOA
   implicit none

   !/-- subroutines
    public  :: OrganicAerosol


   !/-- public

    logical, public, parameter :: ORGANIC_AEROSOLS = .false.   ! true if wanted

    real, public, dimension(A1:A2,K1:K2), save :: Fgas  ! Fraction in gas-phase

   contains
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !+ Driver routine for Secondary Organic Aerosol  module

   subroutine OrganicAerosol(debug_flag)
   logical, intent(in) :: debug_flag  ! for debugging purposes only

     ! empty 

   end subroutine OrganicAerosol

end module OrganicAerosol_ml
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


