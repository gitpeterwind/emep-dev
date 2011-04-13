!TOTO- check BGND OC or OM
module OrganicAerosol_ml

  ! Calculates the amount of condensible species in the gas and aerosol phases. 
  ! Simple starting module, with alpha-K methodology for 2 aromatic
  ! and 4 monoterpene species.
  !
  ! References:
  !   AS2001: Andersson-Sk\"old, Y. and Simpson, D., 2001, Secondary Organic
  !     Aerosol Formation in Northern Europe: a Model Study, JGR, 106, D7,
  !     7357-7374
  ! 
  !   CS2002: Chung, S. H., and J. Seinfeld (2002), Global distribution and
  !     climate forcing of carbonaceous aerosols, \textitJ. Geophys. Res.,
  !     107(D19), doi: 101029/2001JD001397.
  !
  !   Pun2003: Pun, B., S. Wu, C. Seigneur, J. Seinfeld, R. Griffin, and 
  !     A. Pandis (2003), Uncertainties in Modeling Secondary Organic 
  !     Aerosols: Three-dimensional Modeling Studies in Nashville/Western
  !     Tennessee, Environ. Sci.  Tech., 37, 3647-3661.
  !
  !   Pankow:  (1994a,b), An absorption model of gas/particle partitioning
  !     of organic compounds in the atmosphere, Atmospheric Environment,
  !     28(2), 185-188   ** AND **  189-193    
  ! 
  !   S2007: Simpson, D. et al., 2007, 
  !
  !
  ! Usage: call OrganicAerosol from Runchem, after setup of column data
  !  (for meteorology, etc.). The subroutine initialises itself on the 
  !  first call and thereafter modifies two external variables:
  !   xn(k,SOA) : the concentrations of SOA   (For a 1-d column array(in EMEP))
  !   Fgas(k,X) : The fraction of X which is gas and not aeorosol
  !
  ! ----------------------------------------------------------------------------
  !
  ! From Gas/Particle theory, A/G = K.COA,
  ! therefore, Fgas = G/(G+A) = 1/(1+K.COA)
  !
  !-----------------------------------------------------------------------------
  ! NB- we exclude use of gamma for now, but leave commented out code
  !-----------------------------------------------------------------------------
  !
  ! Dave Simpson, August 2001 -- Mar. 2009
  !--------------------------------------------------------------------------

  ! Functions + GridValues + PT only for BGNDOC
   use Functions_ml, only: StandardAtmos_kPa_2_km !ds for use in Hz scaling
   use ChemChemicals_ml,      only : species   ! for molwts
   use ChemSpecs_tot_ml,  S1 => FIRST_SEMIVOL , S2 => LAST_SEMIVOL

!dsrb added:
   use ChemGroups_ml, only :    &
       NONVOLOC_GROUP &
      ,NONVOLEC_GROUP &
      ,ASOA => ASOA_GROUP &
      ,BSOA => BSOA_GROUP &
!rb only_for_semi-volatile_POA      ,POA  => POA_GROUP &
!rb only_for_semi-volatile_POA      ,OPOA  => OPOA_GROUP &
      ,FFUELEC  => FFUELEC_GROUP &
      ,FFUELOC  => FFUELOC_GROUP &
      ,WOODOC  => WOODOC_GROUP

!dsrb doesn'rt exist now   use ChemSOA_ml    ! Includes tabulate_VPsoa, tabVPsoa

   use GridValues_ml, only: sigma_mid 
   use ModelConstants_ml,    only :  PT

   use ModelConstants_ml,    only : CHEMTMIN, CHEMTMAX, &
                                    MasterProc, DEBUG => DEBUG_SOA, &
                                    K2 => KMAX_MID, K1 => KCHEMTOP
   use Par_ml,               only : LIDIM => MAXLIMAX, LJDIM => MAXLJMAX
   use PhysicalConstants_ml, only : AVOG, RGAS_J 
   use Setup_1dfields_ml,    only : itemp, xn => xn_2d, Fgas, Fpart
   use TimeDate_ml,       only: current_date
   implicit none
   private


   !/-- subroutines

    private  :: Init_OrganicAerosol
    public   :: OrganicAerosol


   !/-- public

    logical, public, parameter :: ORGANIC_AEROSOLS = .true.

   ! We store some values in 3-D fields, to allow the next G/P partitioning
   ! calculation to  start off with values of COA, mw and Fgas which 
   ! are about right. Ensures that very few iterations are needed.

!  real,public, save, dimension(S1:S2,LIDIM,LJDIM,K1:K2) :: &
!            Grid_SOA_Fgas           !EXC Grid_SOA_gamma

  real,public, save, dimension(LIDIM,LJDIM,K1:K2) :: Grid_COA
!VBS            Grid_avg_mw &

  real, private, dimension(K1:K2), save :: &
        COA           & ! Org. aerosol, ug/m3
!VBS    avg_mw        & ! Avg. mol. wt.
       ,BGND_OC       & ! FAKE FOR NOW, 0.50 ugC/m3 at surface
       ,BGND_OA         ! Assumed OA/OC=2, -> 1 ug/m3


  !VBS real, private, dimension(S1:S2,K1:K2), save :: &
  ! TMP - we assign Fparty for all species for now, since
  ! it makes it easier to code for nonvol and vol
! From Setup_1dfields now
!  real, private, dimension(1:NSPEC_TOT,K1:K2), save :: &
!              Fgas    & ! Fraction in gas-phase
!             ,Fpart!   & ! Fraction in gas-phase
             !VBS ,tabRTpL   ! = 1.0e-6 * R.T(i)/pL(Ti) for all temps i:
               !EXC ,gamma   & ! activity coefficient

  real, parameter, public :: SMALLFN  = 1.0e-20 ! Minimum value of ug allowed

   !/-- private

   ! ug = array for aerosol masses (ug/m3). Includes non-volatile compounds:
   ! TMP??? Excluding NVOL for now?

    real, private, dimension(S1:S2,K1:K2), save :: ug_semivol 
!dsrb - use new NONVOLOC grpup to define:
    integer, private, parameter :: NUM_NONVOLOC = size(NONVOLOC_GROUP)
    integer, private, parameter :: NUM_NONVOLEC = size(NONVOLEC_GROUP)
    integer, private, parameter :: NUM_NONVOL = NUM_NONVOLOC+NUM_NONVOLEC
    integer, private, parameter, dimension(NUM_NONVOL) ::  &
      NONVOL = (/ POC_F_WOOD,POC_F_FFUEL,POC_C_FFUEL,FFIRE_OC, EC_F_WOOD_NEW,EC_F_WOOD_AGE,EC_C_WOOD,EC_F_FFUEL_NEW,EC_F_FFUEL_AGE,EC_C_FFUEL,FFIRE_BC /) ! OC+EC in partitioning OM
    real, private, dimension(NUM_NONVOL,K1:K2), save :: ug_nonvol 

    real,  private, dimension(S1:S2,CHEMTMIN:CHEMTMAX) :: tabCiStar

    integer, private, save :: NITER = 2              ! No. iterations for Ksoa
    real, private, save :: xn2molem3  ! Conversion from molec/cm3 to mole/m3
    real, private, save :: xn2ugC, ugC2xn


! Set some indices for deposition, used in DryDep_ml and Aqueous_ml
!EGU - find indices of SO4-like particles (PMc included for now - tmp!)
! for SOA we need a separate array, since Fgas will be needed

!dsrb - not needed now?
!dsrb  integer, public, parameter, dimension(NUM_NONVOLOC+NUM_NONVOLEC) ::&
!dsrb   SO4LIKE_DEP = (/ NONVOLOC, NONVOLEC /)
!dsrb
!dsrb  integer, public, parameter, dimension(NUM_VOL) ::&
!dsrb   SOALIKE_DEP = (/ VOL /)


   !/-- DEBUG variables
   !    Usually, DEBUG_SOA gives extra outputs. debug_flag is used to allow
   !    some extra outputs for a gven i,j - set in CTM model.

    logical, private, save :: my_first_call = .true.  ! True for 1st call only
    character(len=20), public, save :: soa_errmsg = "ok"

   contains
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine Init_OrganicAerosol(debug_flag)
   logical, intent(in) :: debug_flag
   integer :: i, it, k
   real, parameter :: kJ = 1000.0  
  !ds Use of standard atmosphere 
   real, dimension(K2), save :: p_kPa, h_km


      ! Need to convert aeros to ug/m3 or ugC/m3.  Since xn is in molecules/cm3
      ! we divide by AVOG to get mole/cm3, multiply by 1e6 to get mole/m3,
      ! by  mol. weight to get g/m3 and by 1e6 to get ug/m3

         xn2molem3 = 1.0/AVOG * 1.0e6
         xn2ugC   = xn2molem3 * 12.0 * 1.0e6
         ugC2xn   = 1/xn2ugC

  !ds 27/7/2003 - Use Standard Atmosphere to get average heights of layers
  ! Taken from Unimod.BVKam2X and Warneck

       p_kPa(:) = 0.001*( PT + sigma_mid(:)*(101325.0-PT) ) ! Pressure in kPa
       h_km     = StandardAtmos_kPa_2_km(p_kPa)
       BGND_OC(:)= 0.5 * 1.005 ! ng/m3 !!! will give 0.5 ugC/m3 at z=0 m

       do k = K1, K2
            BGND_OC(k) = BGND_OC(k) * exp( -h_km(k)/9.1 )
            if(DEBUG .and. MasterProc ) write(*,"(a,i4,2f8.3)"), &
               "BGND_OC ", k, h_km(k),  BGND_OC(k)
       end do
       BGND_OA(:) = 2*BGND_OC(:)   ! Assume OA/OC = 2 for bgnd



       ! Tabulate change of Fcond 
       !  Smog-chamber data from 310K
       ! Ci = 1.0e6*P0/RT 
       ! Now, pi(T) = Ai exp(-Hi/RT)
       ! And pi(T) = Pi(Tref) * exp( H/RT * (1/Tref - 1/T) )
       ! ->  Ci(T) = Ci(Tref) * Tref/T * exp(...)

       do i=S1,S2
         do it=CHEMTMIN,CHEMTMAX

           !dsrb tabCiStar(i,it) = VBS(i)%CiStar * 300.0/it * &
           !dsrb        exp( VBS(i)%DeltaH * kJ/RGAS_J * (1.0/300.0 - 1.0/it) )
           tabCiStar(i,it) = species(i)%CiStar * 300.0/it * &
                  exp( species(i)%DeltaH * kJ/RGAS_J * (1.0/300.0 - 1.0/it) )
         end do
       end do


         if ( MasterProc ) then 
            do i = S1, S2
               write(6,"(a27,i4,a16,i5,i3,8es12.3)") &
                " Tab SOA: MW, Carbons, C*:", i, trim(species(i)%name), &
                 species(i)%molwt, species(i)%carbons, & !VBStabVPsoa(i,298)
                 tabCiStar(i,273), tabCiStar(i,303)
            end do
         end if


       ! Initial values only:

         Fpart(:,:)         = 0.0
         !dsrb Fpart(NONVOLOC,:)  = 1.0 ! ds added EC here also. OK?
         Fpart(NONVOL,:)  = 1.0
         !Grid_SOA_Fgas    = 1.0 - 1.0e-10
         do k = K1, K2
            Grid_COA(:,:,k) = BGND_OA(k)  ! Use OA, not OC here
         end do
           !VBS Grid_avg_mw    = 250.0              ! Da
           !EXC Grid_SOA_gamma = 1.0

   end subroutine Init_OrganicAerosol

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine OrganicAerosol(i_pos,j_pos,debug_flag)

   integer, intent(in) :: i_pos, j_pos
   logical, intent(in) :: debug_flag 

   integer :: i, it, k, iter, ispec   ! loop variables 
!   real, save :: molcc2ngm3 = 1.0e9*1.0e6/AVOG  !molecules/cc-> ng/m3
   real, save :: molcc2ugm3 = 1.0e12/AVOG  !molecules/cc-> ng/m3
!   real, save :: ngm32molcc = AVOG/1.0e15  !molecules/cc <- ng/m3
   real, save :: ugm32molcc = AVOG/1.0e12  !molecules/cc <- ng/m3
   real :: Nmoles , Ksoa
   integer :: nmonth, nday, nhour, seconds

  ! Outputs:
   real :: surfOC, surfASOA, surfBSOA, surfFFUELOC, surfWOODOC, surfBGNDOC !for semivolatile POA runs, surfOPOA, surfPOA


   nmonth = current_date%month
   nday   = current_date%day
   nhour  = current_date%hour
   seconds = current_date%seconds

   if( DEBUG .and. debug_flag) write(unit=*,fmt=*) "Into SOA"

   if ( my_first_call ) then

       call Init_OrganicAerosol(debug_flag)
       my_first_call = .false.

   endif



! Note that xn(SOA) is not strictly needed with the method we have, but it
! is a good first guess of the sum of the condensed phases, enables easy output
! and saves the need for iteration. 
!
! Remember also that Fgas is saved, so will initially have been set by the
! preceding call to OrganicAerosol for a different set of i,j.

 !VBS forall(i=S1:S2,k=K1:K2)
 !VBS     tabRTpL( i, k) = 1.0e-6 * RGAS_J * itemp(k) / tabVpsoa(i,itemp(k))
 !VBS end forall


 ! 1st guesses:

!  Fgas(S1:S2,:)   =  Grid_SOA_Fgas(S1:S2,i_pos,j_pos,:)
  COA(:)          =  Grid_COA(i_pos,j_pos,:)
  !VBS avg_mw(:)    =  Grid_avg_mw(i_pos,j_pos,:)
     !EXC gamma(:,:) = SOA_gamma(:,i_pos,j_pos,:)


 ! ============ Non-volatile species first ============================

  do i = 1, NUM_NONVOL  ! OA/OC for POC about 1.333
    !dsrb ispec = NONVOLOC(i)
    ispec = NONVOL(i)

    ug_nonvol(i,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt

!rb: something is strange with the output from sites (at least) for the non-volatile components test setting Fpart to 1 here also (but should not really be necessary since it is already done in the initialisation above)
    Fpart(i,:)=1.0
    Fgas(i,:)=0.0

  end do

  ! ============ SOA species now, iteratrion needed ===================

  do iter = 1, NITER


      ! Fgas = G/(G+A) = 1/(1+K.COA)
      ! K = tabRTpL/(mw*gamma)

       do ispec = S1, S2

          !VBS Fgas(ispec,:) = 1.0/(1.0+tabRTpL(ispec,:)*COA(:)/avg_mw(:) )  
                                         !EXC *gamma(:,ispec)) )

          Fpart(ispec,:) = COA(:)/( COA(:)+tabCiStar(ispec,itemp(:)) )

          ug_semivol(ispec,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt &
                         * Fpart(ispec,:)

       end do ! ispec

       !ng(:,:) = max(SMALLFN, ng(:,:))

     ! New estimate of COA  (in ug/m3) and avg_mw (g/mole):
     ! (nb. xn in molecules/cm3)

       do k = K1,K2

         COA(k) = sum( ug_semivol(:,k) ) + sum( ug_nonvol(:,k) ) + BGND_OA(k)

         !VBS Nmoles = ( sum( Fpart(VOL,k) * xn(VOL,k) ) + sum( xn(NONVOLOC,k) ) &
         !VBS             + BGND_OC_ng(k)*ngm32molcc/250.0 ) &  ! Assumed MW
         !VBS          * xn2molem3
         !VBS avg_mw(k)    = 1.0e-6 * COA(k)/Nmoles 

       end do  !k
     ! ====================================================================

      if( DEBUG  .and. debug_flag ) then

         if( iter == NITER .and. seconds == 0 ) then
           write(unit=6,fmt="(a,i2,a,3i3,f7.2)") "Iteration ", Niter, &
                " ======== ",nmonth, nday, nhour, itemp(K2)
         !NITER end if 
           write(unit=6,fmt="(a3,a15,3a10,a4,4a10)") "SOA","Species", "xn", &
               "Ci* ", "Ki"," ", "Fpart", "ng"

           do i = 1, NUM_NONVOL
              !dsrb ispec = NONVOLOC(i)
              ispec = NONVOL(i)
              write(unit=6,fmt="(a4,i3,a15,es10.2,2f10.3,a4,es10.3,f13.4)")&
                "NVOL", ispec,&
                species(ispec)%name, xn(ispec,K2),-999.999, &
                -999.999, " => ", Fpart(ispec,K2), 1000.0*ug_nonvol(i, K2)
           end do

           do ispec = S1,S2
              ! K = tabRTpL*COA/(mw*gamma) QUERY COA===!!!!
              !Ksoa = tabRTpL(ispec,K2)*COA(K2)/(avg_mw(K2))
              Ksoa = 1.0/tabCiStar(ispec,itemp(K2)) !just for printout
              write(unit=6,fmt="(a4,i3,a15,3es10.2,a4,es10.3,f13.4)") "SOA ",ispec,&
                species(ispec)%name, xn(ispec,K2), &
                 tabCiStar(ispec,itemp(K2)),& !VBStabVpsoa(ispec,298),
                  Ksoa, " => ", Fpart(ispec,K2), 1000.0*ug_semivol(ispec, K2)
          end do ! ispec
         end if 

          write(unit=6,fmt="(a,i2,f12.6)")  "COA: ", iter, COA(K2)

       end if ! DEBUG

   end do ! ITER

 ! The above iteration has now given new values to: 
 !
 ! 1) COA(1:nz)
 ! 3) Fgas(FIRST_SEMIVOL:LAST_SEMIVOL,1:nz)
 ! 4) Fpart(FIRST_SEMIVOL:LAST_SEMIVOL,1:nz)
 ! 5) ng(FIRST_SEMIVOL:LAST_SEMIVOL,1:nz)    
 !
 ! Note:     ng(FIRST_NONVOLOC:LAST_NONVOLOC,1:nz)
 !       and xn(1:LAST_SEMIVOL,1:nz)  are unaffected by the
 !       iteration.
 !=========================================================================

  ! Set Fgas for later chemistry, and eset 3-D fields

   Fgas(S1:S2,:)  = 1.0 - Fpart(S1:S2,:)
   Grid_COA(i_pos,j_pos,:)              = COA(:)

    !Grid_SOA_Fgas(S1:S2, i_pos,j_pos,:)  = Fgas(S1:S2,:)
    !VBS Grid_avg_mw(i_pos,j_pos,:)       = avg_mw(:)
    !SOA_gamma(i_pos,j_pos,:,:)     = gamma(:,:)

  ! Outputs, ugC/m3

   do k = K1, K2
     xn(AER_ASOA,k)  = sum ( Fpart(ASOA,k) *xn(ASOA,k)  *species(ASOA)%carbons )
     xn(GAS_ASOA,k)  = sum ( Fgas(ASOA,k) * xn(ASOA,k)  *species(ASOA)%carbons )
     xn(AER_BSOA,k)  = sum ( Fpart(BSOA,k) *xn(BSOA,k)  *species(BSOA)%carbons )
     xn(GAS_BSOA,k)  = sum ( Fgas( BSOA,k) *xn(BSOA,k)  *species(BSOA)%carbons )
!rb     xn(AER_POA,k)  = sum ( Fpart(POA,k) *xn(POA,k)  *species(POA)%carbons )
!rb     xn(GAS_POA,k)  = sum ( Fgas( POA,k) *xn(POA,k)  *species(POA)%carbons )
!rb     xn(AER_OPOA,k)  = sum ( Fpart(OPOA,k) *xn(OPOA,k)  *species(OPOA)%carbons )
!rb     xn(GAS_OPOA,k)  = sum ( Fgas( OPOA,k) *xn(OPOA,k)  *species(OPOA)%carbons )
     xn(AER_FFUELOC,k) = sum ( Fpart(FFUELOC,k)*xn(FFUELOC,k) *species(FFUELOC)%carbons)
     xn(AER_WOODOC,k)  = sum ( Fpart(WOODOC,k) *xn(WOODOC,k)  *species(WOODOC)%carbons )

!HARDCODE
!   xn(AER_TBSOA,k)  =  xn(AER_BSOA,k)  ! Just in case TBSOA is wanted for kam
!...............................................................................
   xn(AER_TBSOA,k)  = &
     Fpart( TERP_ug1,k) *xn(TERP_ug1,k)  *species(TERP_ug1)%carbons + &
     Fpart( TERP_ug10,k) *xn(TERP_ug10,k)  *species(TERP_ug10)%carbons + &
     Fpart( TERP_ug1e2,k) *xn(TERP_ug1e2,k)  *species(TERP_ug1e2)%carbons + &
     Fpart( TERP_ug1e3,k) *xn(TERP_ug1e3,k)  *species(TERP_ug1e3)%carbons
   xn(AER_IBSOA,k)  = &
     Fpart( ISOP_ug1,k) *xn(ISOP_ug1,k)  *species(ISOP_ug1)%carbons + &
     Fpart( ISOP_ug10,k) *xn(ISOP_ug10,k)  *species(ISOP_ug10)%carbons + &
     Fpart( ISOP_ug1e2,k) *xn(ISOP_ug1e2,k)  *species(ISOP_ug1e2)%carbons + &
     Fpart( ISOP_ug1e3,k) *xn(ISOP_ug1e3,k)  *species(ISOP_ug1e3)%carbons
!   xn(AER_SBSOA,k)  = &
!     Fpart( SESQ_ug1,k) *xn(SESQ_ug1,k)  *species(SESQ_ug1)%carbons + &
!    Fpart( SESQ_ug10,k) *xn(SESQ_ug10,k)  *species(SESQ_ug10)%carbons + &
!    Fpart( SESQ_ug1e2,k) *xn(SESQ_ug1e2,k)  *species(SESQ_ug1e2)%carbons + &
!    Fpart( SESQ_ug1e3,k) *xn(SESQ_ug1e3,k)  *species(SESQ_ug1e3)%carbons 
!...............................................................................
   end do
   xn(AER_BGNDOC,:)    = ugC2xn * BGND_OC(:)      ! FAKE FOR NOW, 0.5 ug/m3 at surface

 ! for convencience:
   xn(AER_POC,:) =  xn(AER_FFUELOC,:) ! + xn(AER_POA,:)
   xn(AER_OC,:)  =  xn(AER_ASOA,:)+xn(AER_BSOA,:)+xn(AER_FFUELOC,:) + &
                    xn(AER_WOODOC,:) + xn(AER_BGNDOC,:) ! +xn(AER_OPOA,:)+xn(AER_POA,:)

   surfASOA  = xn2ugC* xn(AER_ASOA,K2) ! sum ( Fpart(ASOA,K2) * xn(ASOA,K2)*species(ASOA)%carbons )
   surfBSOA  = xn2ugC* xn(AER_BSOA,K2) ! sum ( Fpart(BSOA,K2) * xn(BSOA,K2)*species(BSOA)%carbons )
!rb   surfOPOA  = xn2ugC* xn(AER_OPOA,K2) ! sum ( Fpart(BSOA,K2) * xn(BSOA,K2)*species(BSOA)%carbons )
!rb   surfPOA  = xn2ugC* xn(AER_POA,K2) ! 
   surfFFUELOC = xn2ugC* xn(AER_FFUELOC,K2) ! sum ( Fpart(FFUELOC,K2) * xn(FFUELOC,K2)*species(FFUELOC)%carbons )
   surfWOODOC  = xn2ugC* xn(AER_WOODOC,K2) ! sum ( Fpart(WOODOC,K2) * xn(WOODOC,K2)*species(WOODOC)%carbons )

! FAKEugC*sum ( Fpart(BGND_OC,K2) * &
                  !FAKE       xn(BGND_OC,K2)*species(BGND_OC)%carbons )

   surfBGNDOC  = BGND_OC(K2)
   surfOC    = surfASOA+surfBSOA+surfFFUELOC+surfWOODOC+surfBGNDOC ! +surfOPOA+surfPOA


   !/ Sum of Biogenics -----------------------

   if ( debug_flag .and. seconds == 0 ) then

       k=20
       write(unit=6,fmt="(a,3i3,2f7.2,f5.2,20es9.2)")"xns ug ", &
         nmonth, nday, nhour, &
         COA(20), surfOC,surfBGNDOC, surfASOA, surfBSOA,surfFFUELOC, surfWOODOC!, &
!vbs      xn2ugC*xn(AER_IBSOA,k), xn2ugC*xn(AER_TBSOA,k), xn2ugC*xn(AER_SBSOA,k)

   endif

 end subroutine OrganicAerosol

end module OrganicAerosol_ml
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

