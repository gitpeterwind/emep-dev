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
  ! 
  ! NOTE!! 11 AUg 2011: This version of My_SOA_ml is under construction for use with the EMEP VBS
  !        SOA models. The model is under development and changes are likely without warning.
  !        If you want to use this model it is probably a good idea to contact David Simpson
  !        and/or Robert Bergström first.  
  !--------------------------------------------------------------------------

  ! Functions + GridValues + PT only for BGNDOC
   use Functions_ml, only: StandardAtmos_kPa_2_km !ds for use in Hz scaling
   use ChemChemicals_ml,      only : species   ! for molwts
   use ChemSpecs_tot_ml,  S1 => FIRST_SEMIVOL , S2 => LAST_SEMIVOL

!dsrb added:
   use ChemGroups_ml, only :    &
      NONVOLPCM_GROUP &
      ,NVABSOM_GROUP & ! contains the nonvolatile absorbing species for OA partitioning
      ,ASOA => ASOA_GROUP &
      ,BSOA => BSOA_GROUP &
      ,FFUELOA25 => FFUELOA25_GROUP & ! for summing all semivolatile POA+OPOA from "fossil fuel"
      ,ECFINE => ECFINE_GROUP &                    
      ,VFFUELOA25  => VFFUELOA25_GROUP &   ! for VBS FFUELOA25_C
      ,OFFUELOA25  => OFFUELOA25_GROUP & ! for VBS FFUELOA25_C with aging
      ,OX_OFFUELOA25 => OFFOAEXTRAO25_GROUP &! for VBS FFUELOA25_C with aging
      ,WOODOA25  => WOODOA25_GROUP &    ! for summing all semivolatile WOOD BURNING OA
      ,VWOODOA25  => VWOODOA25_GROUP &    ! for VBS WOOD BURNING OA
      ,OWOODOA25  => OWOODOA25_GROUP &  ! for VBS WOOD BURNING OA with aging
      ,OX_OWOODOA25 => OWOODOAEXTRAO25_GROUP & ! for VBS WOOD BURNING OA with aging
      ,FFIREOA25  => FFIREOA25_GROUP &    ! for summing all semivolatile Wildfire OA
      ,VFFIREOA25  => VFFIREOA25_GROUP &    ! for VBS Wildfire OA
      ,OFFIREOA25  => OFFIREOA25_GROUP &  ! for VBS Wildfire OA with aging
      ,OX_OFFIREOA25 => OFFIREOAEXTRAO25_GROUP & ! for VBS Wildfire OA with aging
      ,FFUELEC  => FFUELEC_GROUP &
      ,NVFFUELOC25  => NVFFUELOC25_GROUP & ! only non-volatile FFUELOC emissions (PM2.5 fraction)
      ,NVFFUELOC10  => NVFFUELOC10_GROUP & ! only non-volatile FFUELOC emissions (PM10 fraction, that is all of it)
      ,NVWOODOC25  => NVWOODOC25_GROUP & ! only non-volatile WOODOC emissions, zero in VBS-PAX type runs
      ,NVFFIREOC25  => NVFFIREOC25_GROUP ! only non-volatile FFIREOC emissions, zero in VBS-PAX type runs

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

  real, private, dimension(K1:K2), save :: &
        COA           & ! Org. aerosol, ug/m3  (this version does not include EC as absorber)
       ,BGND_OC       & ! FAKE FOR NOW, 0.50 ugC/m3 at surface
       ,BGND_OA         ! Assumed OA/OC=2, -> 1 ug/m3


  !VBS real, private, dimension(S1:S2,K1:K2), save :: &
  ! TMP - we assign Fpart for all species for now, since
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
    integer, private, parameter :: NUM_NONVOLPCM = size(NONVOLPCM_GROUP)
    integer, private, parameter :: NUM_NVABSOM = size(NVABSOM_GROUP)
!    integer, private, parameter, dimension(NUM_NONVOL) ::  &
!      NONVOL = (/ NONVOLOC_GROUP, NONVOLEC_GROUP /) ! OC+EC in partitioning OM
    real, private, dimension(NUM_NVABSOM,K1:K2), save :: ug_nonvol 
    !real, private, dimension(K1:K2), save :: ug_ecf ! CityZen added 

    real,  private, dimension(S1:S2,CHEMTMIN:CHEMTMAX) :: tabCiStar

    integer, private, save :: NITER = 2              ! No. iterations for Ksoa
    real, private, save :: xn2molem3  ! Conversion from molec/cm3 to mole/m3
    real, private, save :: xn2ugC, ugC2xn


! Set some indices for deposition, used in DryDep_ml and Aqueous_ml
!EGU - find indices of SO4-like particles (PMc included for now - tmp!)
! for SOA we need a separate array, since Fgas will be needed

   !/-- DEBUG variables
   !    Usually, DEBUG_SOA gives extra outputs. debug_flag is used to allow
   !    some extra outputs for a gven i,j - set in CTM model.

    logical, private, save :: my_first_call = .true.  ! True for 1st call only
    character(len=20), public, save :: soa_errmsg = "ok"
    character(len=*), public, parameter :: SOA_MODULE_FLAG="VBS"

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

!rb: C*-values are given for 298K according to most(?) publications.
           tabCiStar(i,it) = species(i)%CiStar * 298./it * &
                  exp( species(i)%DeltaH * kJ/RGAS_J * (1.0/298. - 1.0/it) )
         end do
       end do


         if ( MasterProc ) then 
            do i = S1, S2
               write(6,"(a,i4,a16,f7.1,i3,8es12.3)") &
                " Tab SOA: MW, Carbons, C*:", i, trim(species(i)%name), &
                 species(i)%molwt, species(i)%carbons, & !VBStabVPsoa(i,298)
                 tabCiStar(i,273), tabCiStar(i,303)
            end do
         end if


       ! Initial values only:

         Fpart(:,:)         = 0.0
         !dsrb Fpart(NONVOLOC,:)  = 1.0 ! ds added EC here also. OK?
         Fpart(NONVOLPCM_GROUP,:)  = 1.0
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
   real :: surfOC25, surfASOA, surfBSOA, surfFFUELOC25, surfWOODOC25, surfBGNDOC, surfOFFUELOA25_C, surfFFUELOA25_C, surfOWOODOA25, surfWOODOA25, surfOFFIREOA25, surfFFIREOA25


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
 ! NVABSOM - Only include fine OM! That is no EC and no coarse OM!

  do i = 1, NUM_NVABSOM  ! OA/OC for POC about 1.333
    !dsrb ispec = NONVOLOC(i)
    ispec = NVABSOM_GROUP(i)

    ug_nonvol(i,:) = molcc2ugm3 * xn(ispec,:)*species(ispec)%molwt

!rb: something is strange with the output from sites (at least) for the non-volatile components test setting Fpart to 1 here also (but should not really be necessary since it is already done in the initialisation above)
    Fpart(i,:)=1.0
    Fgas(i,:)=0.0

  end do
  !do k = K1, K2 
  !  ug_ecf(k) = molcc2ugm3 * sum( xn(ECFINE,k) ) * 12.0 !! *species(ispec)%molwt
  !end do

  ! ============ SOA species now, iteration needed ===================

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

           do i = 1, NUM_NONVOLPCM
              !dsrb ispec = NONVOLOC(i)
              ispec = NONVOLPCM_GROUP(i)
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

!CITYZEN. PCM_F is for output only. Has MW 1 to avoid confusion with OC
!do not use ugC outputs, just ug

   xn(PART_OM_F,:)  =  COA(:) * ugC2xn * 12.0

    !Grid_SOA_Fgas(S1:S2, i_pos,j_pos,:)  = Fgas(S1:S2,:)
    !VBS Grid_avg_mw(i_pos,j_pos,:)       = avg_mw(:)
    !SOA_gamma(i_pos,j_pos,:,:)     = gamma(:,:)

  ! Outputs, ugC/m3

  !RB: Would like to be able to store also total OM (not only OC) 
  !    at least for some components. And for a "total" OM and/or OM2.5 and OM10. 
  !    Total TC and EC (and TC2.5, TC10, EC2.5 and EC10) would also be useful.
  !    Also perhaps the names of these species should reflect 
  !    that they are in units of C?   
   do k = K1, K2
     xn(PART_ASOA_C,k)  = sum ( Fpart(ASOA,k) *xn(ASOA,k)  *species(ASOA)%carbons )
     xn(GAS_ASOA_C,k)  = sum ( Fgas(ASOA,k) * xn(ASOA,k)  *species(ASOA)%carbons )
     xn(PART_BSOA_C,k)  = sum ( Fpart(BSOA,k) *xn(BSOA,k)  *species(BSOA)%carbons )
     xn(GAS_BSOA_C,k)  = sum ( Fgas( BSOA,k) *xn(BSOA,k)  *species(BSOA)%carbons )
     xn(NONVOL_FFUELOC25,k) = sum ( xn(NVFFUELOC25,k) *species(NVFFUELOC25)%carbons)
     xn(NONVOL_FFUELOC10,k) = sum ( xn(NVFFUELOC10,k) *species(NVFFUELOC10)%carbons)
     xn(NONVOL_WOODOC25,k)  = sum ( xn(NVWOODOC25,k)  *species(NVWOODOC25)%carbons )
     xn(NONVOL_FFIREOC25,k)  = sum ( xn(NVFFIREOC25,k)  *species(NVFFIREOC25)%carbons )
!RB want to have PART_FFUELOA25/FFIREOA/WOODOA_C working also with nonvolatile POA emissions, test this hard coded version first
     xn(PART_FFUELOA25_C,k)  = sum ( Fpart(VFFUELOA25,k) *xn(VFFUELOA25,k) *species(VFFUELOA25)%carbons ) + &
	xn(NONVOL_FFUELOC25,k)
     xn(GAS_FFUELOA_C,k)  = sum ( Fgas(VFFUELOA25,k) *xn(VFFUELOA25,k) *species(VFFUELOA25)%carbons )
     xn(PART_OFFUELOA25_C,k)  = sum ( Fpart(OFFUELOA25,k) *xn(OFFUELOA25,k)  *species(OFFUELOA25)%carbons )
     xn(GAS_OFFUELOA_C,k)  = sum ( Fgas(OFFUELOA25,k) *xn(OFFUELOA25,k)  *species(OFFUELOA25)%carbons )
     xn(PART_XO_OFFLOA25_O,k)  = sum ( Fpart(OX_OFFUELOA25,k) *xn(OX_OFFUELOA25,k) )
     xn(GAS_XO_OFFLOA_O,k)  = sum ( Fgas( OX_OFFUELOA25,k) *xn(OX_OFFUELOA25,k) )
     xn(PART_WOODOA25_C,k)  = sum ( Fpart(VWOODOA25,k) *xn(VWOODOA25,k)  *species(VWOODOA25)%carbons ) + &
	xn(NONVOL_WOODOC25,k) 
     xn(PART_OWOODOA25_C,k)  = sum ( Fpart(OWOODOA25,k) *xn(OWOODOA25,k)  *species(OWOODOA25)%carbons )
     xn(PART_XO_OWDOA25_O,k)  = sum ( Fpart(OX_OWOODOA25,k) *xn(OX_OWOODOA25,k) )
     xn(PART_FFIREOA25_C,k)  = sum ( Fpart(VFFIREOA25,k) *xn(VFFIREOA25,k)  *species(VFFIREOA25)%carbons ) + &
	xn(NONVOL_FFIREOC25,k)
     xn(PART_OFFIREOA25_C,k)  = sum ( Fpart(OFFIREOA25,k) *xn(OFFIREOA25,k)  *species(OFFIREOA25)%carbons )
     xn(PART_XO_OFFIOA25_O,k)  = sum ( Fpart(OX_OFFIREOA25,k) *xn(OX_OFFIREOA25,k) )
!Test for storing in ug/m3, Use with caution! 
     xn(PART_ASOA_OM,k)  = sum ( ug_semivol(ASOA,k) ) * ugC2xn * 12.0 
     xn(PART_BSOA_OM,k)  = sum ( ug_semivol(BSOA,k) ) * ugC2xn * 12.0 
!RB want to have PART_FFUELOA25/FFIREOA/WOODOA_OM working also with nonvolatile POA emissions, test this hard coded version first
     xn(PART_FFUELOA25_OM,k)  = sum ( ug_semivol(FFUELOA25,k) ) * ugC2xn * 12.0 + xn(NONVOL_FFUELOC25,k) * 1.25
     xn(PART_WOODOA25_OM,k)  = sum ( ug_semivol(WOODOA25,k) ) * ugC2xn * 12.0 + xn(NONVOL_WOODOC25,k) * 1.7
     xn(PART_FFIREOA25_OM,k)  = sum ( ug_semivol(FFIREOA25,k) ) * ugC2xn * 12.0 + xn(NONVOL_FFIREOC25,k) * 1.7

!HARDCODE
!   xn(AER_TBSOA,k)  =  xn(AER_BSOA,k)  ! Just in case TBSOA is wanted for kam
!...............................................................................
   xn(PART_TBSOA_C,k)  = &
     Fpart( TERP_ng100,k) *xn(TERP_ng100,k)  *species(TERP_ng100)%carbons + &
     Fpart( TERP_ug1,k) *xn(TERP_ug1,k)  *species(TERP_ug1)%carbons + &
     Fpart( TERP_ug10,k) *xn(TERP_ug10,k)  *species(TERP_ug10)%carbons + &
     Fpart( TERP_ug1e2,k) *xn(TERP_ug1e2,k)  *species(TERP_ug1e2)%carbons + &
     Fpart( TERP_ug1e3,k) *xn(TERP_ug1e3,k)  *species(TERP_ug1e3)%carbons
   xn(PART_IBSOA_C,k)  = &
     Fpart( ISOP_ng100,k) *xn(ISOP_ng100,k)  *species(ISOP_ng100)%carbons + &
     Fpart( ISOP_ug1,k) *xn(ISOP_ug1,k)  *species(ISOP_ug1)%carbons + &
     Fpart( ISOP_ug10,k) *xn(ISOP_ug10,k)  *species(ISOP_ug10)%carbons + &
     Fpart( ISOP_ug1e2,k) *xn(ISOP_ug1e2,k)  *species(ISOP_ug1e2)%carbons + &
     Fpart( ISOP_ug1e3,k) *xn(ISOP_ug1e3,k)  *species(ISOP_ug1e3)%carbons
!   xn(PART_SBSOA,k)  = &
!     Fpart( SESQ_ug1,k) *xn(SESQ_ug1,k)  *species(SESQ_ug1)%carbons + &
!    Fpart( SESQ_ug10,k) *xn(SESQ_ug10,k)  *species(SESQ_ug10)%carbons + &
!    Fpart( SESQ_ug1e2,k) *xn(SESQ_ug1e2,k)  *species(SESQ_ug1e2)%carbons + &
!    Fpart( SESQ_ug1e3,k) *xn(SESQ_ug1e3,k)  *species(SESQ_ug1e3)%carbons 
!...............................................................................
   end do
   xn(NONVOL_BGNDOC,:)    = ugC2xn * BGND_OC(:)      ! FAKE FOR NOW, 0.5 ug/m3 at surface

 ! for convencience:
!removed not used?   xn(PART_POC,:) =  xn(PART_FFUELOC,:) + xn(PART_FFUELOA25_C,:)
! Hopefully this works for both VBS-NPNA and VBS-PAx type runs. 
! WARNING! Test changing to WOODOC25, FFIREOC25 and NONVOLOC25 may cause problems here, especially if coarse components are inlcuded later! But this is rather hard coded anyway...
   xn(PART_OC10,:)  =  xn(PART_ASOA_C,:)+xn(PART_BSOA_C,:)+xn(NONVOL_FFUELOC10,:) + &
                    xn(NONVOL_WOODOC25,:)+ xn(NONVOL_BGNDOC,:)+xn(PART_OFFUELOA25_C,:)+ &
                    xn(PART_FFUELOA25_C,:)+xn(PART_OWOODOA25_C,:)+xn(PART_WOODOA25_C,:)+ &
                    xn(PART_OFFIREOA25_C,:)+xn(PART_FFIREOA25_C,:)+ xn(NONVOL_FFIREOC25,:)
!WARNING! The below will NOT work for NPNA (nonvolatile) type runs. These include fine and coarse OC in NONVOL_FFUELOC. So for these runs the FFUELOC-contribution has to be added separately!!!
!Test shidting to NONVOL_FFUELOC25!
   xn(PART_OC25,:)  =  xn(PART_ASOA_C,:)+xn(PART_BSOA_C,:)+ xn(NONVOL_FFUELOC25,:)+&
                    xn(NONVOL_WOODOC25,:)+ xn(NONVOL_BGNDOC,:)+xn(PART_OFFUELOA25_C,:)+ &
                    xn(PART_FFUELOA25_C,:)+xn(PART_OWOODOA25_C,:)+xn(PART_WOODOA25_C,:)+ &
                    xn(PART_OFFIREOA25_C,:)+xn(PART_FFIREOA25_C,:)+ xn(NONVOL_FFIREOC25,:)

   surfASOA  = xn2ugC* xn(PART_ASOA_C,K2) ! sum ( Fpart(ASOA,K2) * xn(ASOA,K2)*species(ASOA)%carbons )
   surfBSOA  = xn2ugC* xn(PART_BSOA_C,K2) ! sum ( Fpart(BSOA,K2) * xn(BSOA,K2)*species(BSOA)%carbons )
! Is this used? Then we need the nonvolatile parts as well!?
   surfOFFUELOA25_C  = xn2ugC* xn(PART_OFFUELOA25_C,K2) !
   surfFFUELOA25_C  = xn2ugC* xn(PART_FFUELOA25_C,K2) ! 
   surfOWOODOA25 = xn2ugC* xn(PART_OWOODOA25_C,K2) !
   surfWOODOA25 = xn2ugC* xn(PART_WOODOA25_C,K2) !
   surfOFFIREOA25 = xn2ugC* xn(PART_OFFIREOA25_C,K2) !
   surfFFIREOA25 = xn2ugC* xn(PART_FFIREOA25_C,K2) !

!CHECK this! Is it used also for partitioning runs? Then the partitioning parts are needed as well!!!
   surfFFUELOC25 = xn2ugC* xn(NONVOL_FFUELOC25,K2) ! sum ( Fpart(FFUELOC,K2) * xn(FFUELOC,K2)*species(FFUELOC)%carbons )
   surfWOODOC25  = xn2ugC* xn(NONVOL_WOODOC25,K2) ! sum ( Fpart(WOODOC,K2) * xn(WOODOC,K2)*species(WOODOC)%carbons )

! FAKEugC*sum ( Fpart(BGND_OC,K2) * &
                  !FAKE       xn(BGND_OC,K2)*species(BGND_OC)%carbons )

   surfBGNDOC  = BGND_OC(K2)
   surfOC25    = surfASOA+surfBSOA+surfFFUELOC25+surfWOODOC25+surfBGNDOC+surfOFFUELOA25_C+surfFFUELOA25_C+surfWOODOA25+surfOWOODOA25+surfFFIREOA25+surfOFFIREOA25 ! 
! DO we need surfOC10 as well? 

   !/ Sum of Biogenics -----------------------

   if ( debug_flag .and. seconds == 0 ) then

       k=20
       write(unit=6,fmt="(a,3i3,2f7.2,f5.2,20es9.2)")"xns ug ", &
         nmonth, nday, nhour, &
         xn(PART_OM_F,20)*xn2ugC, & 
         COA(20), surfOC25,surfBGNDOC, surfASOA, surfBSOA,surfFFUELOA25_C, &
         surfOFFUELOA25_C,surfFFUELOC25, surfWOODOA25, surfOWOODOA25, surfWOODOC25, &
         surfFFIREOA25, surfOFFIREOA25!, &
!vbs      xn2ugC*xn(PART_IBSOA,k), xn2ugC*xn(PART_TBSOA,k), xn2ugC*xn(PART_SBSOA,k)

   endif

 end subroutine OrganicAerosol

end module OrganicAerosol_ml
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

