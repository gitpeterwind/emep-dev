!>  AerosolCalls.f90 - A component of the EMEP MSC-W Chemical transport Model
!!****************************************************************************! 
!> Options for aerosol-gas equilibrium partitioning:
!!
!! * EMEP - old EMEP scheme with (NH4)1.5 SO4
!! * MARS - run MARS equilibrium model
!! * EQSAM - run EQSAM equilibrium model
!! * ISORROPIA - run ISORROPIA II equilibrium model (in testing)

module AerosolCalls

 use AeroConstants_mod,     only: AERO
 use Ammonium_mod,          only: Ammonium
 use CheckStop_mod,         only: StopAll, CheckStop
 use ChemDims_mod,          only: NSPEC_SHL
 use ChemSpecs_mod,         only: species
 use Chemfields_mod,        only: PM25_water, PM25_water_rh50, & !H2O_eqsam, & !PMwater 
                                   cfac
 use Config_module,         only: KMAX_MID, KCHEMTOP, MasterProc, USES,&
                                  SO4_ix, HNO3_ix, NO3_f_ix, NH3_ix, NH4_f_ix
 use Debug_module,          only: DEBUG   ! -> DEBUG%EQUIB
 use EQSAM4clim_ml,        only :  EQSAM4clim
! use EQSAM_v03d_mod,        only: eqsam_v03d
 use MARS_mod,              only: rpmares, rpmares_2900, DO_RPMARES_new
 use PhysicalConstants_mod, only: AVOG
 use SmallUtils_mod,        only: find_index
 use ZchemData_mod,         only: xn_2d, temp, rh, pp
 implicit none
 private

 !/-- public           !!  true if wanted

 public :: AerosolEquilib
 public :: emep2MARS, emep2EQSAM, Aero_Water, Aero_Water_rh50, Aero_Water_MARS
 private :: emep2isorropia
                    
!    logical, public, parameter :: AERO_DYNAMICS     = .false.  &  
!                                , EQUILIB_EMEP      = .false.  & !old Ammonium stuff
!                                , EQUILIB_MARS      = .true.  & !MARS
!                                , EQUILIB_EQSAM     = .false.     !EQSAM
                                
!.. Number of aerosol sizes (1-fine, 2-coarse, 3-'giant' for sea salt )
!    integer, public, parameter :: NSIZE = 5
!           !   FINE_PM = 1, COAR_NO3 = 2, COAR_SS = 3, COAR DUST = 4,pollen = 5    


  integer, private, save :: iSeaSalt ! zero if no seasalt compounds

contains

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine AerosolEquilib(debug_flag)
    logical, intent(in) :: debug_flag
!    integer, intent(in)  :: i, j
    logical, save :: my_first_call=.true.
    character(len=*),parameter:: dtxt='AeroEqui:'
    
    if( my_first_call ) then
      iSeaSalt = find_index('SeaSalt_f',species(:)%name )
      call CheckStop(USES%SEASALT.and.iSeaSalt<1,dtxt//"iSeaSalt neg")
      if(  MasterProc ) then
        write(*,*) 'AerosolEquilib: chosen: ',AERO%EQUILIB 
        write(*,*) 'AerosolEquilib water: chosen: ',AERO%EQUILIB_WATER
        write(*,*) 'AerosolEquilib seasalt index: ',iSeaSalt
      end if
    end if
    select case ( AERO%EQUILIB )
      case ( 'EMEP' )
        call ammonium()
      case ( 'MARS' , 'MARS_2900', 'GEOSCHEM')
        call emep2MARS(debug_flag)
      case ( 'EQSAM' )
        call emep2EQSAM(debug_flag)
      case ( 'ISORROPIA' )
        !NOV22 call StopAll('Isorropia problems found. Removed for now')
        call emep2Isorropia(debug_flag)
      case default
        if( my_first_call .and. MasterProc ) then
          write(*,*) 'WARNING! AerosolEquilib, nothing valid chosen: '//AERO%EQUILIB
        end if
    end select
    my_first_call = .false.

  end subroutine AerosolEquilib

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 ! Adapted from List 10, p130, Isoropia manual

 subroutine emep2isorropia(debug_flag)
  logical, intent(in) :: debug_flag

  real, dimension(8) :: wi = 0.0, wt
  real, dimension(3) :: gas
  real, dimension(15) :: aerliq
  real, dimension(19) :: aersld
  real, parameter, dimension(2) :: CNTRL =  (/ 0, 0 /)
  real, dimension(9) :: other
  !real :: rhi, tempi
  character(len=15) :: scase

  !EMEP
  real, parameter :: Ncm3_to_molesm3 = 1.0e6/AVOG    ! #/cm3 to moles/m3
  real, parameter :: molesm3_to_Ncm3 = 1.0/Ncm3_to_molesm3
  integer :: k
!??? BUG - WAS NOT INITIALISED. WILL TEST HERE, guessing ?= wt(4) - total nitrate (moles/m3):
   real :: tmpno3, tmpnh3, tmpnhx, tmphno3


  ! WI(1)  = max(FLOOR2, xn_2d(Na,k))  / species(Na)%molwt  * Ncm3_to_molesm3 - QUERY WHY MW Here????
  !from isocon.f90:
  !Input = WI
  !Concentrations, expressed in moles/m3. Depending on the type of
  !     problem solved (specified in CNTRL(1)), WI contains either
  !     GAS+AEROSOL or AEROSOL only concentratios.
  !     WI(1) - sodium    WI(2) - sulfate   WI(3) - ammonium   WI(4) - nitrate
  !     WI(5) - chloride  WI(6) - calcium   WI(7) - potassium  WI(8) - magnesium
  !Output
  !     Total concentrations (GAS+AEROSOL) of species, expressed in moles/m3.
  !     If the foreward probelm is solved (CNTRL(1)=0), array WT is
  !     identical to array WI.
  !     WT(1) - total sodium     WT(2) - total sulfate   WT(3) - total ammonium    WT(4) - total nitrate
  !     WT(5) - total chloride   WT(6) - total calcium   WT(7) - total potassium   WT(8) - total magnesium
  ! 2. [GAS] !     real array of length [03].
  !     Gaseous species concentrations, expressed in moles/m3.
  !     GAS(1) - NH3 !     GAS(2) - HNO3 !     GAS(3) - HCl
  ! 3. [AERLIQ] !     real array of length [15].
  !     Liquid aerosol species concentrations, expressed in moles/m3.
  !     AERLIQ(01) - H+(aq) !     AERLIQ(02) - Na+(aq) !     AERLIQ(03) - NH4+(aq) !     AERLIQ(04) - Cl-(aq)
  !     AERLIQ(05) - SO4--(aq) !     AERLIQ(06) - HSO4-(aq) !     AERLIQ(07) - NO3-(aq) !     AERLIQ(08) - H2O
  !     AERLIQ(09) - NH3(aq) (undissociated) !     AERLIQ(10) - HNCl(aq) (undissociated) !     AERLIQ(11) - HNO3(aq) (undissociated) !     AERLIQ(12) - OH-(aq)
  !     AERLIQ(13) - Ca2+(aq) !     AERLIQ(14) - K+(aq) !     AERLIQ(15) - Mg2+(aq)

  !  4. [AERSLD] !     real array of length [19].
  !     Solid aerosol species concentrations, expressed in moles/m3.
  !     AERSLD(01) - NaNO3(s)  !        AERSLD(02) - NH4NO3(s)    ! AERSLD(03) - NaCl(s)     ! AERSLD(04) - NH4Cl(s)
  !     AERSLD(05) - Na2SO4(s) !        AERSLD(06) - (NH4)2SO4(s) ! AERSLD(07) - NaHSO4(s)   ! AERSLD(08) - NH4HSO4(s)
  !     AERSLD(09) - (NH4)4H(SO4)2(s) ! AERSLD(10) - CaSO4(s)     ! AERSLD(11) - Ca(NO3)2(s) ! AERSLD(12) - CaCl2(s)
  !     AERSLD(13) - K2SO4(s)  !        AERSLD(14) - KHSO4(s)     ! AERSLD(15) - KNO3(s)     ! AERSLD(16) - KCl(s)
  !     AERSLD(17) - MgSO4(s)  !        AERSLD(18) - Mg(NO3)2(s)  ! AERSLD(19) - MgCl2(s)


  do k = KMAX_MID, KMAX_MID  ! TESTING KCHEMTOP, KMAX_MID

    WI(1)  = 0.0 !FINE sum( xn_2d(SS_GROUP,k) ) * Ncm3_to_molesm3
    WI(2)  = xn_2d(SO4_ix,k)             * Ncm3_to_molesm3
    !NOV22 WI(3)  = sum( xn_2d(RDN_GROUP,k) ) * Ncm3_to_molesm3  !NH3, NH4
    WI(3)  = ( xn_2d(NH3_ix,k) + xn_2d(NH4_f_ix,k) ) * Ncm3_to_molesm3  !NH3, NH4
    !FINE WI(4)  = ( xn_2d(NO3_F,k) + xn_2d(NO3_C,k) + xn_2d(HNO3,k) )&
    WI(4)  = ( xn_2d(NO3_f_ix,k) + xn_2d(HNO3_ix,k) )&
               * Ncm3_to_molesm3
    WI(5)  =0.0 !FINE  WI(1)  ! Cl only from sea-salt. Needs consideration!

    !NOV22 testing:
    tmpnh3 = xn_2d(NH3_ix,k)
    tmpnhx = tmpnh3 + xn_2d(NH4_f_ix,k)
    tmphno3 = xn_2d(HNO3_ix,k)
    tmpno3 = tmphno3 + xn_2d(NO3_f_ix,k)

    call isoropia ( wi, rh(k), temp(k), CNTRL,&
                    wt, gas, aerliq, aersld, scase, other)

   ! gas outputs are in moles/m3(air)

    xn_2d(NH3_ix,k)  = max(0.0, gas(1)) * molesm3_to_Ncm3
    xn_2d(HNO3_ix,k) = max(0.0, gas(2)) * molesm3_to_Ncm3
    !xn_2d(HCl,k) = gas(3) * molesm3_to_Ncm3

   ! aerosol outputs are in moles/m3(air)
   ! 1=H+, 2=Na+, 3=NH4+, 4=Cl-, 5=SO42-, 6=HSO4-, 7=NO3-, 8=Ca2+
   ! 9=K+, 10=Mg2+
    !xn_2d(NH4_F,k) = MOLAL(3) ???? Unlinekly MOLAL?
    !NOV22 - get some v.small neg., so use max below. Test properly later.
    xn_2d(NH4_f_ix,k) =  max(0.0, wt(3) - gas(1)) * molesm3_to_Ncm3

   ! Just use those needed:
   ! QUERY: Is NaNO3 always solid? Ans = No!

     !xn_2d(NO3_c,k ) = aeroHCl * molesm3_to_Ncm3 ! assume all HCl from NaNO3 formation?
     !FINE xn_2d(NO3_f,k ) = tmpno3 - xn_2d(NO3_c,k ) - xn_2d(HNO3,k)

     !tmpno3 = wt(4) * molesm3_to_Ncm3  ! NOV22  wt4=nitrate
     !xn_2d(NO3_f_ix,k ) = tmpno3 - xn_2d(HNO3_ix,k)
    !NOV22 - get some v.small neg., so use max below. Test properly later.
     xn_2d(NO3_f_ix,k ) = max(0.0, wt(4) - gas(2))  * molesm3_to_Ncm3 

     !if ( xn_2d(NO3_f_ix,k )  < 0.0 .or.  xn_2d(NH4_f_ix,k) < 0.0 ) then
     !   if ( xn_2d(NO3_f_ix,k ) 
     !   print "(a,99e12.3)", 'NEGNO3 pre',tmpnh3, tmpnhx, tmphno3, tmpno3
     !   print "(a,99e12.3)", 'NEGNO3 xn', xn_2d(NH3_ix,k ), xn_2d(NH4_f_ix,k)+xn_2d(NH3_ix,k ), xn_2d(HNO3_ix,k ), xn_2d(NO3_f_ix,k ), xn_2d(SO4_ix,k )
     !   print "(a,99e12.3)", 'NEGNO3 wi', wi(1:4) ! 1 - total sodium 2 - total sulfate 3 - total ammonium 4 - total nitrate
     !   print "(a,99e12.3)", 'NEGNO3 wt', wt(1:4) ! 1 - total sodium 2 - total sulfate 3 - total ammonium 4 - total nitrate
     !   print "(a,99e12.3)", 'NEGNO3 gas', gas !     GAS(1) - NH3 !     GAS(2) - HNO3 !     GAS(3) - HCl
     !   print "(a,99e12.3)", 'NEGNO3 diffs',  xn_2d(NO3_f_ix,k ),  xn_2d(NH4_f_ix,k)
     !   call StopAll('NEGNO3')
     !end if

    if( debug_flag ) then 
      write(*, "(a,2f8.3,99g12.3)") "ISORROPIA ", rh(k), temp(k), gas
    end if
    !call StopAll("ISOR")
    
  end do
end subroutine emep2isorropia
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine emep2MARS(debug_flag)

 !..................................................................
 ! Pretty old F. Binkowski code from EPA CMAQ-Models3
 ! JGR, 108, D6, 4183
 !..................................................................

 logical, intent(in) :: debug_flag 
 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  
 real, parameter ::    FLOOR2 = 1.0E-9         ! minimum concentration  

 !.. local
  real    :: so4in, no3in, nh4in, hno3in, nh3in,   &
             aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out,   &
             coef
  integer :: k, errmark
 !-----------------------------------
  call CheckStop( SO4_ix<1, "emep2MARS:SO4 not defined" )
  call CheckStop( NH4_f_ix<1, "emep2MARS:NH4_f not defined" )
  call CheckStop( HNO3_ix<1, "emep2MARS:HNO3 not defined" )
  call CheckStop( NO3_f_ix<1, "emep2MARS:NO3_f not defined" )
  call CheckStop( NH3_ix<1, "emep2MARS: NH3 not defined" )

   coef = 1.e12 / AVOG

   do k = KCHEMTOP, KMAX_MID
  

!//.... molec/cm3 -> ug/m3
! Use FLOOR2 = 1.0e-8 molec/cm3 for input. Too many problems
      so4in  = max(FLOOR2, xn_2d(SO4_ix,k)) * species(SO4_ix)%molwt  *coef
      hno3in = max(FLOOR2, xn_2d(HNO3_ix,k))* species(HNO3_ix)%molwt *coef 
      nh3in  = max(FLOOR2, xn_2d(NH3_ix,k)) * species(NH3_ix)%molwt  *coef
      no3in  = max(FLOOR2, xn_2d(NO3_f_ix,k)) * species(NO3_f_ix)%molwt  *coef
      nh4in  = max(FLOOR2, xn_2d(NH4_f_ix,k)) * species(NH4_f_ix)%molwt  *coef

 !--------------------------------------------------------------------------                
      if(AERO%EQUILIB=='MARS')then 
         call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='MARS_2900')then!svn version 2908, for testing if there are significant differences. will be deleted
         call rpmares_2900 (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='GEOSCHEM')then
         call DO_RPMARES_new (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag)
      end if

 !--------------------------------------------------------------------------

      if( DEBUG%EQUIB) then
        call CheckStop(gNO3out< 0.0, "XMARS: gNO3out")
        call CheckStop(gNH3out< 0.0, "XMARS: gNH3out")
        call CheckStop(aNO3out< 0.0, "XMARS: aNO3out")
        call CheckStop(aNH4out< 0.0, "XMARS: aNH4out")
      end if ! DEBUG%EQUIB

      xn_2d(HNO3_ix,k)  = max (FLOOR, gNO3out / (species(HNO3_ix)%molwt *coef) )
      xn_2d(NH3_ix,k)   = max (FLOOR, gNH3out / (species(NH3_ix)%molwt  *coef) )
      xn_2d(NO3_f_ix,k)  = max (FLOOR, aNO3out / (species(NO3_f_ix)%molwt  *coef) )
      xn_2d(NH4_f_ix,k)  = max (FLOOR, aNH4out / (species(NH4_f_ix)%molwt  *coef) )

   end do  ! K-levels

 end subroutine emep2MARS

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine emep2EQSAM(debug_flag)
 logical, intent(in) :: debug_flag 
! integer, intent(in)  :: i, j

    !..................................................................
    ! EQSAM4clim - Equlibrium Simplified Aerosol Model by Swen Metzger
    !             version v10 is implemented here
    ! Metzger, S., B. Steil, M. Abdelkader, K. Klingmüller, L. Xu, 
    ! J.E. Penner, C. Fountoukis, A. Nenes, and J. Lelieveld, 2016; 
    ! Aerosol Water Parameterization: A single parameter framework; 
    ! Atmos. Chem. Phys., 16, 7213–7237, doi:10.5194/acp-16-7213-2016; 
    ! https://www.atmos-chem-phys.net/16/7213/2016/acp-16-7213-2016.html
    !..................................................................


 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  


 !.. local

  real    :: so4in(KCHEMTOP:KMAX_MID),   &
             no3in(KCHEMTOP:KMAX_MID),   &
             nh4in(KCHEMTOP:KMAX_MID),   &
             hno3in(KCHEMTOP:KMAX_MID),  &
             nh3in(KCHEMTOP:KMAX_MID),   &
             aH2Oin(KCHEMTOP:KMAX_MID),  &
! The following input is not in use
             NAin(KCHEMTOP:KMAX_MID),    &
             CLin(KCHEMTOP:KMAX_MID),    &

             aSO4out(KCHEMTOP:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), & 
             gNH3out(KCHEMTOP:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID), &
             aNAout(KCHEMTOP:KMAX_MID),  &
             aCLout(KCHEMTOP:KMAX_MID),  &
             gCLout(KCHEMTOP:KMAX_MID),  &
             gSO4out(KCHEMTOP:KMAX_MID)

 !-----------------------------------


  if ( debug_flag  ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4_ix,20),xn_2d(HNO3_ix,20),&
               xn_2d(NH3_ix,20),xn_2d(NO3_f_ix,20),xn_2d(NH4_f_ix,20)
  end if

!//.... molec/cm3 -> micromoles/m**3
    so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4_ix,KCHEMTOP:KMAX_MID)  *1.e12/AVOG
    hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3_ix,KCHEMTOP:KMAX_MID) *1.e12/AVOG
    nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3_ix,KCHEMTOP:KMAX_MID)  *1.e12/AVOG 
    no3in(KCHEMTOP:KMAX_MID)  = xn_2d(NO3_f_ix,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(NH4_f_ix,KCHEMTOP:KMAX_MID)*1.e12/AVOG
!    aH2Oin(KCHEMTOP:KMAX_MID) = H2O_eqsam(i,j,KCHEMTOP:KMAX_MID)
    if ( iSeaSalt > 0 ) then
      NAin(KCHEMTOP:KMAX_MID)   = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*0.306e12/AVOG
      CLin(KCHEMTOP:KMAX_MID)   = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*0.55e12/AVOG
    end if

 !--------------------------------------------------------------------------                
  
!    call eqsam_v03d (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh,temp,pp, &
!                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,             &
!                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout)
    
    call EQSAM4clim (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh, temp,  & !aH2Oin, &
                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,            &
                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout, KCHEMTOP, KMAX_MID)
  
 !--------------------------------------------------------------------------

!//.... micromoles/m**3  -> molec/cm3 
!      xn_2d(NO3,KCHEMTOP:KMAX_MID)  = FLOOR !different for ACID/OZONE

      xn_2d(HNO3_ix,KCHEMTOP:KMAX_MID)  = max(FLOOR,gNO3out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
      xn_2d(NH3_ix,KCHEMTOP:KMAX_MID)   = max(FLOOR,gNH3out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
      xn_2d(NO3_f_ix,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNO3out(KCHEMTOP:KMAX_MID)*AVOG*1.e-12 ) 
      xn_2d(NH4_f_ix,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNH4out(KCHEMTOP:KMAX_MID)*AVOG*1.e-12 )
      xn_2d(SO4_ix,KCHEMTOP:KMAX_MID)   = max(FLOOR,aSO4out(KCHEMTOP:KMAX_MID) *AVOG*1.e-12 )
!      H2O_eqsam(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )

 if ( debug_flag ) then ! Selected debug cell
    write(*,*)'After EQSAM',xn_2d(SO4_ix,20),xn_2d(HNO3_ix,20),&
               xn_2d(NH3_ix,20),xn_2d(NO3_f_ix,20),xn_2d(NH4_f_ix,20)
  end if

 end subroutine emep2EQSAM


 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


 !water 

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine Aero_water_rh50(i,j, debug_flag)

  !.....................................................................
  ! EQSAM is called before every daily output to calculate aerosol water 
  ! at T=20C and Rh = 50%. This should model the particle water content 
  ! for gravitationally determined PM mass
  ! Tsyro, S. (2005). To what extent can aerosol water explain the 
  ! discrepancy between model calculated and gravimetric PM10 and PM2.5?. 
  ! Atmos. Chem.. Phys., 5, 602, 1-8, 2005.
  !.....................................................................

 implicit none

 integer, intent(in)  :: i, j
 logical, intent(in)  :: debug_flag
 !.. local
! integer, parameter    ::  KCHEMTOP = KMAX_MID

!  AerosolCalls.f90(319): error #6592: This symbol must be a defined parameter, an enumerator, or an argument of 
!  an inquiry function that evaluates to a compile-time constant.   [KMAX_MID]
!  integer, parameter    ::  KCHEMTOP = KMAX_MID
!--------------------------------------^

  real    :: so4in(KMAX_MID:KMAX_MID),   &
             no3in(KMAX_MID:KMAX_MID),   &
             nh4in(KMAX_MID:KMAX_MID),   &
             hno3in(KMAX_MID:KMAX_MID),  &
             nh3in(KMAX_MID:KMAX_MID),   &
             aH2Oin(KMAX_MID:KMAX_MID),  &
             NAin(KMAX_MID:KMAX_MID)  ,  &
             CLin(KMAX_MID:KMAX_MID) ,   &
!-- output
             aSO4out(KMAX_MID:KMAX_MID), &
             aNO3out(KMAX_MID:KMAX_MID), &
             aH2Oout(KMAX_MID:KMAX_MID), &
             aNH4out(KMAX_MID:KMAX_MID), & 
             gNH3out(KMAX_MID:KMAX_MID), &
             gNO3out(KMAX_MID:KMAX_MID), &
             aNAout(KMAX_MID:KMAX_MID),  &
             aCLout(KMAX_MID:KMAX_MID),  &
             gCLout(KMAX_MID:KMAX_MID),  &
             gSO4out(KMAX_MID:KMAX_MID), &

             rlhum(KMAX_MID:KMAX_MID),tmpr(KMAX_MID:KMAX_MID)

 !-----------------------------------


  if ( debug_flag ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4_ix,20),xn_2d(HNO3_ix,20),&
               xn_2d(NH3_ix,20),xn_2d(NO3_f_ix,20),xn_2d(NH4_f_ix,20)
  end if


!//.... molec/cm3 -> micromoles/m**3
      so4in(KMAX_MID:KMAX_MID)  = xn_2d(SO4_ix,KMAX_MID:KMAX_MID)*1.e12/AVOG
      hno3in(KMAX_MID:KMAX_MID) = xn_2d(HNO3_ix,KMAX_MID:KMAX_MID)*1.e12/AVOG
      nh3in(KMAX_MID:KMAX_MID)  = xn_2d(NH3_ix,KMAX_MID:KMAX_MID)*1.e12/AVOG 
      no3in(KMAX_MID:KMAX_MID)  = xn_2d(NO3_f_ix,KMAX_MID:KMAX_MID)*1.e12/AVOG
      nh4in(KMAX_MID:KMAX_MID)  = xn_2d(NH4_f_ix,KMAX_MID:KMAX_MID)*1.e12/AVOG
      NAin(KMAX_MID:KMAX_MID)   = xn_2d(iSeaSalt,KMAX_MID:KMAX_MID)*0.306e12/AVOG
      CLin(KMAX_MID:KMAX_MID)   = xn_2d(iSeaSalt,KMAX_MID:KMAX_MID)*0.55e12/AVOG

!.. Rh = 50% and T=20C
      rlhum(KMAX_MID:KMAX_MID) = 0.5
      tmpr (KMAX_MID:KMAX_MID) = 293.15


 !--------------------------------------------------------------------------                
   
    call EQSAM4clim (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rlhum, tmpr,  & 
                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,               &
                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout, KMAX_MID, KMAX_MID)     

 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 

      PM25_water_rh50 (i,j)             = max(0., aH2Oout(KMAX_MID) )


 end subroutine  Aero_water_rh50
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 subroutine Aero_water(i,j, debug_flag)


    !..................................................................
    ! EQSAM4clim calculates PM water 1. at ambient conditions and 
    !          2. at Rh=50% and T=20C (comparison with gravimentric PM)
    !..................................................................

 integer, intent(in)  :: i, j
 logical, intent(in)  :: debug_flag

 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  

 !===  For ambient conditions =====================

 !.. local

  real    :: so4in(KCHEMTOP:KMAX_MID), so4ins(KMAX_MID:KMAX_MID),  &
             no3in(KCHEMTOP:KMAX_MID), no3ins(KMAX_MID:KMAX_MID),  &
             nh4in(KCHEMTOP:KMAX_MID), nh4ins(KMAX_MID:KMAX_MID),  &
             hno3in(KCHEMTOP:KMAX_MID),hno3ins(KMAX_MID:KMAX_MID), &
             nh3in(KCHEMTOP:KMAX_MID), nh3ins(KMAX_MID:KMAX_MID),  &
             aH2Oin(KCHEMTOP:KMAX_MID),aH2Oins(KMAX_MID:KMAX_MID), &
             NAin(KCHEMTOP:KMAX_MID),  Nains(KMAX_MID:KMAX_MID),   &
             CLin(KCHEMTOP:KMAX_MID),  CLins(KMAX_MID:KMAX_MID),   &
!.. output
             aSO4out(KCHEMTOP:KMAX_MID), aSO4outs(KMAX_MID:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), aNO3outs(KMAX_MID:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID), aH2Oouts(KMAX_MID:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), aNH4outs(KMAX_MID:KMAX_MID), & 
             gNH3out(KCHEMTOP:KMAX_MID), gNH3outs(KMAX_MID:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID), gNO3outs(KMAX_MID:KMAX_MID), &
             aNAout(KCHEMTOP:KMAX_MID),  aNAouts(KMAX_MID:KMAX_MID),  &
             aCLout(KCHEMTOP:KMAX_MID),  aCLouts(KMAX_MID:KMAX_MID),  &
             gCLout(KCHEMTOP:KMAX_MID),  gCLouts(KMAX_MID:KMAX_MID),  &
             gSO4out(KCHEMTOP:KMAX_MID), gSO4outs(KMAX_MID:KMAX_MID), &

             rlhum(KCHEMTOP:KMAX_MID), tmpr(KCHEMTOP:KMAX_MID),  &
             rlhums(KMAX_MID:KMAX_MID),tmprs(KMAX_MID:KMAX_MID)

 !-----------------------------------


  if ( debug_flag  ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4_ix,20),xn_2d(HNO3_ix,20),&
               xn_2d(NH3_ix,20),xn_2d(NO3_f_ix,20),xn_2d(NH4_f_ix,20)
  end if

!//.... molec/cm3 -> micromoles/m**3
    so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4_ix,KCHEMTOP:KMAX_MID)  *1.e12/AVOG
    hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3_ix,KCHEMTOP:KMAX_MID) *1.e12/AVOG
    nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3_ix,KCHEMTOP:KMAX_MID)  *1.e12/AVOG 
    no3in(KCHEMTOP:KMAX_MID)  = xn_2d(NO3_f_ix,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
    nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(NH4_f_ix,KCHEMTOP:KMAX_MID)*1.e12/AVOG
    NAin(KCHEMTOP:KMAX_MID)   = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*0.306e12/AVOG
    CLin(KCHEMTOP:KMAX_MID)   = xn_2d(iSeaSalt,KCHEMTOP:KMAX_MID)*0.55e12/AVOG

    rlhum(KCHEMTOP:KMAX_MID) = rh(:)
    tmpr(KCHEMTOP:KMAX_MID)  = temp(:)

 !--------------------------------------------------------------------------                
 
    call EQSAM4clim (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rlhum, tmpr,  & !aH2Oin, &
                     aSO4out, aNO3out, aNH4out, aNAout, aCLout,            &
                     gSO4out, gNH3out, gNO3out, gClout, aH2Oout, KCHEMTOP, KMAX_MID)     

 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 

      PM25_water(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )



 !===  PM water for Rh=50% and T=20C  conditions ================================

!//.... molec/cm3 -> micromoles/m**3
      so4ins  = so4in(KMAX_MID:KMAX_MID)
      hno3ins = hno3in(KMAX_MID:KMAX_MID)
      nh3ins  = nh3in(KMAX_MID:KMAX_MID)
      no3ins  = no3in(KMAX_MID:KMAX_MID)
      nh4ins  = nh4in(KMAX_MID:KMAX_MID)
      NAins   = NAin(KMAX_MID:KMAX_MID)
      CLins   = CLin(KMAX_MID:KMAX_MID)

!//.... Only for the lowest layer KMAX_MID

      rlhums = 0.5
      tmprs  = 293.15


 !--------------------------------------------------------------------------                
   
    call EQSAM4clim (so4ins, hno3ins,no3ins,nh3ins,nh4ins,NAins,CLins, rlhums, tmprs,  & 
                     aSO4outs, aNO3outs, aNH4outs, aNAouts, aCLouts,               &
                     gSO4outs, gNH3outs, gNO3outs, gClouts, aH2Oouts, KMAX_MID, KMAX_MID)     

 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 

      PM25_water_rh50 (i,j)             = max(0., aH2Oouts(KMAX_MID) )


 end subroutine Aero_water


 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine Aero_water_MARS(i,j, debug_flag)

 !..................................................................
 ! Pretty old F. Binkowski code from EPA CMAQ-Models3
 ! JGR, 108, D6, 4183
 !..................................................................

 integer, intent(in)  :: i, j
 logical, intent(in)  :: debug_flag

 !.. local
  real    :: rlhum(KCHEMTOP:KMAX_MID), tmpr(KCHEMTOP:KMAX_MID)
  real    :: so4in, no3in, nh4in, hno3in, nh3in,   &
             aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out,   &
             coef
  integer :: k, errmark
 !-----------------------------------
  if(AERO%EQUILIB/='MARS' .and. AERO%EQUILIB/='MARS_2900' .and. AERO%EQUILIB/='GEOSCHEM')then
     PM25_water(i,j,:) = 0.0
     return
  endif
 
   coef = 1.e12 / AVOG

 !.. PM2.5 water at ambient conditions (3D)
   rlhum(:) = rh(:) 
   tmpr(:)  = temp(:)

    do k = KCHEMTOP, KMAX_MID
  
!//.... molec/cm3 -> ug/m3
      so4in  = xn_2d(SO4_ix,k) * species(SO4_ix)%molwt  *coef
      hno3in = xn_2d(HNO3_ix,k)* species(HNO3_ix)%molwt *coef 
      nh3in  = xn_2d(NH3_ix,k) * species(NH3_ix)%molwt  *coef
      no3in  = xn_2d(NO3_f_ix,k) * species(NO3_f_ix)%molwt  *coef
      nh4in  = xn_2d(NH4_f_ix,k) * species(NH4_f_ix)%molwt  *coef


 !--------------------------------------------------------------------------                
      if(AERO%EQUILIB=='MARS')then 
         call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='MARS_2900')then!svn version 2908, for testing if there are significant differences. will be deleted
         call rpmares_2900 (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='GEOSCHEM')then
         call DO_RPMARES_new (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      end if
      !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3) 
      PM25_water(i,j,k) = max (0., aH2Oout )

    end do  ! k-loop

!.. PM2.5 water at equilibration conditions for gravimetric PM (Rh=50% and t=20C)
                            
    rlhum(:) = 0.5
    tmpr(:)  = 293.15
    k = KMAX_MID
!//.... molec/cm3 -> ug/m3
      so4in  = xn_2d(SO4_ix,k) * species(SO4_ix)%molwt  *coef *cfac(SO4_ix-NSPEC_SHL,i,j) 
      hno3in = xn_2d(HNO3_ix,k)* species(HNO3_ix)%molwt *coef *cfac(HNO3_ix-NSPEC_SHL,i,j)
      nh3in  = xn_2d(NH3_ix,k) * species(NH3_ix)%molwt  *coef *cfac(NH3_ix-NSPEC_SHL,i,j)
      no3in  = xn_2d(NO3_f_ix,k) * species(NO3_f_ix)%molwt *coef *cfac(NO3_f_ix-NSPEC_SHL,i,j)
      nh4in  = xn_2d(NH4_f_ix,k) * species(NH4_f_ix)%molwt *coef *cfac(NH4_f_ix-NSPEC_SHL,i,j)
!--------------------------------------------------------------------------                
     if(AERO%EQUILIB=='MARS')then 
         call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='MARS_2900')then!svn version 2908, for testing if there are significant differences. will be deleted
         call rpmares_2900 (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      elseif(AERO%EQUILIB=='GEOSCHEM')then
         call DO_RPMARES_new (so4in, hno3in,no3in ,nh3in, nh4in , rlhum(k), tmpr(k),   &
              aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, &
              ERRMARK,debug_flag) 
      end if
  !--------------------------------------------------------------------------

      PM25_water_rh50 (i,j) = max (0., aH2Oout )


 end subroutine  Aero_water_MARS
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end module AerosolCalls
