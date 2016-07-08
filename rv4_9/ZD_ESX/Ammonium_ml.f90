! <Ammonium_ml.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2013 met.no
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
module Ammonium_ml
 !----------------------------------------------------------------------------
 ! Module to set up and process the NH4-NO3-SO4 reaction system
 ! This code represents the earliest EMEP work in the area (Hov et al., 
 ! Calculation of the distribution of NOx compounds in Europe. Troposheric
 ! ozone. Regional and global scale interactions, I. Isaksen (Ed.), 
 ! D. Reidel, 1988, 239-262.)
 !
 ! This method was subsequenty replaced with EQSAM and MARS alternatives. 
 ! Still, this code is the simplest and most robust.
 !
 ! Usage:
 !   "call ammonium()"  - from Runchem
 !     - on the first call this runs the tabulation routines. On all calls
 !     the equilibrium relationships are calculated and run to establish
 !     new values of ammonium sulphate (AMSU), NH3, HNO3, SO4 and 
 !     ammonium nitrate (AMNI).
 !
 !     Dec 2002 Routine change to treat SO4-NH3-HNO3-NO3_f-NH4_f system instead
 !     This makes code flexible with regards to which eq solver you choos: 
 !     Ammonium, MARS or EQSAM.
 !     In principle, this is exactly the same as using the old indices,
 !     however, SO4 which goes into the chemical solver is now the total 
 !     sulphate, whereas with the old indices it was only free sulphate.
 !
 ! ESX modifications, October 2013
 !----------------------------------------------------------------------------

 use ChemSpecs,            only : species
 use ModelConstants,       only : CHEMTMIN, CHEMTMAX   ! Temp. range
 use SmallUtils_ml,           only : find_index
 use ZchemData,            only : xChem
 use Zmet_ml,              only : rh, temp, itemp, amk =>M
 implicit none
 private

 !/- subroutines:
 public   :: ammonium          ! Sets up most tables
 private  :: tabulate          ! Sets up most tables, and calls tab_rct_rates

 !/- Outputs: - updated xn_2d concentrations after equilibrium


 !/-- Local:

  real, private, dimension(CHEMTMIN:CHEMTMAX), save  :: &
                   tab_rhdel   &  ! RH of deliquescence for ammonium nitrate
                  ,tab_Kp_amni &  ! Equil. constant, nh3 + hno3  <--> nh4no3
                  ,tab_MozP1   &  ! Mozurkewich P1 value for Kaq
                  ,tab_MozP2   &  ! Mozurkewich P2 value for Kaq
                  ,tab_MozP3   ! &  ! Mozurkewich P3 value for Kaq

   logical, private, save :: my_first_call = .true.

 contains

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !--------------------------------------------------------------------------
   ! Calculates the equilibrium constant for the ammonium-suphate
   ! ammonium nitrate, Kp and Kaq (here denoted rcKaq). 
   ! Ref: EMEP Report 2/98, pages B:3, Mozurkewich (1993)
   !
   !      Kpaq = [ P1 -P2(1-rh/100) + P3(1-rh/100)^2 ] .(1-rh/100)**1.75. Kp
   !
   ! Units :  Kp, Kaq : (molecules cm-3)^2
   !          rc      ????
   !--------------------------------------------------------------------------
   !   Calculates the distribution of NH3, (NH4)1.5SO4, NH4NO3
   !   - needs more text...
   !   nov 2002 hf Changed from NH3-AMSU-AMNI-HNO3
   !                       to   NH3-aNH4-pNO3_f-HNO3
   !   in order to have same structure as with EQSAM and MARS 
   !-------------------------------------------------------------------------

  subroutine ammonium(debug_flag)
    logical, intent(in) :: debug_flag 


  ! Local

   integer, save :: iSO4, iNO3, iNH3, iNH4, iHNO3
   real, pointer, dimension(:) :: SO4, NO3, NH3, NH4, HNO3
   character(len=15) :: fmt="(a,20(a,es8.2))" 

   real, dimension(size(rh)):: &
               eqnh3, delteq, freeSO4 &
             , rcnh4     &  ! equilib. value
             , rhd, Kp                    &   ! deliq. rh, Kp
             , roappm                     &  ! density in ppm?
             , humd,humdsqrt,humdsqrt2       ! humd = 1-rh

   if ( my_first_call ) then

     call tabulate()

     iSO4 = find_index( "SO4", species(:)%name )
     iNO3 = find_index( "NO3_F", species(:)%name )
     iNH3 = find_index( "NH3", species(:)%name )
     iNH4 = find_index( "NH4_F", species(:)%name )
     iHNO3 =find_index( "HNO3", species(:)%name )
     if(debug_flag) write(*,*) "AMMONIUM setup", iSO4,iNH3, iNO3, iNH4, iHNO3
     
      my_first_call = .false.
   endif

  if( any ( (/ iSO4, iNH3, iNO3, iNH4, iHNO3 /) < 1 ) ) return

  SO4  => xChem(iSO4,:)
  NO3  => xChem(iNO3,:)
  NH3  => xChem(iNH3,:)
  NH4  => xChem(iNH4,:)
  HNO3 => xChem(iHNO3,:)

  if(debug_flag) write(*,fmt) "AMMONIUM IN ", " SO4:",  SO4(1), &
    " NO3:",NO3(1), " NH3:",NH3(1), " NH4:",NH4(1)," HNO3:",HNO3(1)

   !--------------------------------------------------------------------------
   ! Calculate the equilibrium constant for the ammonium-suphate
   !--------------------------------------------------------------------------


      rhd(:) =  tab_rhdel( itemp(:) )

      Kp(:)  =  tab_Kp_amni( itemp(:) )

    ! Initialize rcnh4 to tab_Kp_amni,need roappm

      roappm(:) = amk(:)*1.0e-9   !PPB
      rcnh4(:)  =  tab_Kp_amni( itemp(:) )*roappm(:)* roappm(:)

!  The lines below are a CPU-efficient way of calculating the
!  power of 1.75  for  Mozurkewich Kp

      where ( rh >= rhd )

        humd = 1.0001 - rh
        humdsqrt = sqrt(humd)
        humdsqrt2 = sqrt(humdsqrt)*humdsqrt
        Kp = (   tab_MozP1(itemp) &
               - tab_MozP2(itemp)*humd &
               + tab_MozP3(itemp)*humd*humd  ) *humd*humdsqrt2*Kp

        !DOUBLED ? roappm = amk*1.0e-9   !PPB
        rcnh4  = Kp * roappm * roappm 

      end where

   ! Sulfate not in form of (NH4)1.5SO4 or NH4NO3:

     freeSO4=SO4-((NH4-NO3)*2./3.) 
     freeSO4=max(0.0,freeSO4)


     where ( 1.5*freeSO4(:) >  NH3 ) ! free SO4 (not in Amm.S form) in excess of NH3

        NH4 = NH4 +  NH3        !hf

        NH3 = 0.

     elsewhere !NH3 in excess

        NH4 = NH4 + freeSO4*1.5 !hf

        NH3 = NH3 - freeSO4*1.5

              
      ! The equilibrium concentration of NH3 is:
        eqnh3 = (NH3 - HNO3)*0.5 + sqrt( 0.25*(NH3 -HNO3)**2 + rcnh4 )+1.

      ! eqnh3  of order 10^20.

        delteq     = eqnh3 - NH3
        delteq     = min(delteq,NO3)

        NO3     = NO3 - delteq

        NH3     = NH3  + delteq
        HNO3    = HNO3 + delteq

        delteq  = min(delteq,NH4)!in  theory not necessary, 
                              !but numerics make very small neg value possible
        NH4     = NH4  - delteq !hf amsu

     end where

  if(debug_flag) write(*,*) "AMMONIUM T", itemp(1)
  if(debug_flag) write(*,fmt) "AMMONIUM OUT", " SO4:",  SO4(1), &
    " NO3:",NO3(1), " NH3:",NH3(1), " NH4:",NH4(1)," HNO3:",HNO3(1), &
        " TC ", itemp(1)-273.0, " Keq ", Kp(1)
   end subroutine ammonium
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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

       tab_rhdel(:) = exp( 618.3/t(:) - 2.551 )

    !-------------------------------------------------------------------
    !    Equilibrium constant (Kp):  NH3 + HNO3  <-------> NH4NO3   
    !    Ref: Mozurkewich (1993)
    !    Units : (molecule/cm3)^2 for Kp
    !
    !      lnKp = 118.87 - 24084.0/T - 6.025* ln(T)
    !
    ! nb: older documentation had + 24084!
    ! c.f. Seinfeld, eqn 9.91, p.532 

       tab_Kp_amni(:) = exp( 118.87 - 24084.0/t(:)-6.025*alog(t(:)) )

    !-------------------------------------------------------------------
    !    temp. dependant constrants for calcolating dissos. rate 
    !    for  the formation of ammonium nitrate  
    !    Ref: Mozurkewich (1993)
    !    n.b. EMEP report 2/98 had 2446 in P3, but 24.46 is correct

       tab_MozP1(:) = exp( -135.94 +  8763.0/t(:) + 19.12*alog( t(:) ) )
       tab_MozP2(:) = exp( -122.65 +  9969.0/t(:) + 16.22*alog( t(:) ) )
       tab_MozP3(:) = exp( -182.61 + 13875.0/t(:) + 24.46*alog( t(:) ) )


  !-------------------------------------------------------------------
  end subroutine tabulate

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module Ammonium_ml
