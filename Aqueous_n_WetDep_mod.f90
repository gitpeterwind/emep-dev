module Aqueous_mod
!-----------------------------------------------------------------------
! Aqueous scavenging and cloud-processing routines.
!
! The scavenging of soluble compounds is based upon the work of Berge
! and Jakobsen (1998) and Eliassen and Saltbones (1983).
! Simple scavenging coefficients are used. A distinction is made
! between in-cloud and sub-cloud scavenging.
!
! Usage:
!
! imports CM_WetDep.inc, which gives the chemistry-dependent species and
! associated rates (generated by GenChem.pl)
!
! Setup_Clouds(i,j)  called from Runchem_mod
! WetDeposition(i,j) called from Runchem_mod if prec. clouds are present
!
! Refs:
! Berge, E., 1993, Coupling of wet scavenging of sulphur to clouds in a
!    numerical weather prediction model, Tellus, 45B, 1-22
! Berge, E. and Jakobsen, H.A., 1998, {A regional scale multi-layer model for
!    the calculation of long-term transport and deposition of air pollution in
!    Europe, Tellus, 50,205-223
! Eliassen, A. and Saltbones, J., 1983, Modelling of long-range transport of
!    sulphur over Europe: a two year model run and some experiments, Atmos.
!    Environ., 17, 1457-1473
! Seland, O. and T. Iversen (1999) A scheme for black carbon and
!    sulphate aerosols tested in a hemispheric scale, Eulerian dispersion
!    model. Atm.  Env. Vol. 33, pp.2853-2879.
!-----------------------------------------------------------------------

use My_Derived_mod,    only: nOutputWdep ! number WDEP used
use CheckStop_mod,     only: CheckStop, StopAll
use ChemDims_mod             
use ChemSpecs_mod,     only: species_adv, FIRST_SEMIVOL, LAST_SEMIVOL
use ChemFunctions_mod,  only: IUPAC_TROE
use ChemGroups_mod,    only: ChemGroups
use Config_module,only: &
    CHEMTMIN, CHEMTMAX      &       ! -> range of temperature
   ,MasterProc              &
   ,KMAX_MID                &       ! -> ground, k=20
   ,KUPPER                  &       ! -> top of cloud-chemistry, k=6
   ,KCHEMTOP                &       ! -> top of chemistry, now k=2
   ,dt => dt_advec          &       ! -> model timestep
   ,IOU_INST                &       ! Index: instantaneous values
   ,USES, WDEP_WANTED       &       ! Which outputs wanted!
   ,SO2_ix, SO4_ix, NH3_ix,NH4_f_ix, NO3_f_ix, HNO3_ix ! indices that may or may not exist in species
use Debug_module,       only: DEBUG  !  => DEBUG%AQUEOUS, DEBUG%MY_WETDEP, DEBUG%pH
use DerivedFields_mod,  only: f_2d, d_2d     ! Contains Wet deposition fields
use GasParticleCoeffs_mod, only: WetCoeffs, WDspec, WDmapping, nwdep
use GridValues_mod,     only: gridwidth_m,xm2,dA,dB,i_fdom,j_fdom,A_mid,B_mid
use Io_mod,             only: IO_DEBUG, datewrite
use LocalFractions_mod, only: lf_wetdep, wetdep_lf
use MassBudget_mod,     only : wdeploss,totwdep
use MetFields_mod,       only: pr, roa, z_bnd, cc3d, cw_met
use MetFields_mod,       only: ps
use OrganicAerosol_mod,  only: ORGANIC_AEROSOLS
use Par_mod,             only: limax,ljmax, me,li0,li1,lj0,lj1
use PhysicalConstants_mod,only: GRAV,AVOG,  &    ! "g" & Avogadro's No.
                               ATWAIR,&         ! Mol. weight of air(Jones,1992)
                               RGAS_ATML,RGAS_J ! Gas-constant
use ZchemData_mod,  only: xn_2d, M, Fpart, Fgas, &
                              temp, itemp        ! temperature (K)
use SmallUtils_mod,      only: find_index
use Units_mod,           only: Group_Scale,group_umap

implicit none
private

! Subroutines:
public :: Init_WetDep       ! Call from emepctm
public :: WetDep_Budget     ! called here
public :: init_aqueous
public :: Setup_Clouds      ! characterises clouds and calls WetDeposition if rain
public :: WetDeposition     ! simplified setup_wetdep
private:: tabulate_aqueous
private:: get_frac
private:: setup_aqurates

! Outputs:
logical, public, save,allocatable, dimension(:) :: &
  incloud              ! True for in-cloud k values
! Variables used in module:
real, private, save,allocatable, dimension(:) :: &
  pr_acc                  ! Accumulated precipitation
! (here for debugging)
real, private, save,allocatable, dimension(:) :: &
  pH,so4_aq,no3_aq,nh4_aq,nh3_aq,hso3_aq,so2_aq,so32_aq,co2_aq,hco3_aq   ! pH in cloud

integer, private, save  :: kcloudtop   ! k-level of highest-cloud
integer, private, save  :: ksubcloud   ! k-level just below cloud

real, private, parameter :: & ! Define limits for "cloud"
  PR_LIMIT = 1.0e-7,  &      ! for accumulated precipitation
  CW_LIMIT = 1.0e-7, &      ! for cloud water, kg(H2O)/kg(air)
  B_LIMIT  = 1.0e-3         ! for cloud cover (fraction)

!hf  real, private, save :: &      ! Set in init below
!hf      INV_Hplus          &      ! = 1.0/Hplus       (1/H+)
!hf     ,INV_Hplus0p4              ! = INV_Hplus**0.4  (1/H+)**0.4

! The Henry's law coefficients, K, given in units of M or M atm-1,
! are calculated as effective. A factor K1fac = 1+K1/H+ is defined
! here also.

integer, public, parameter :: &
!hf pH
  NHENRY  = 5, &  ! No. of species with Henry's law applied
  NK1     = 1, &  ! No. of species needing effective Henry's calc.
  IH_SO2  = 1, &
  IH_H2O2 = 2, &
  IH_O3   = 3, &
  IH_NH3  = 4, &     !hf pH
  IH_CO2  = 5
! Aqueous fractions:
real, save,allocatable, public,  dimension(:,:) :: frac_aq
real, private, dimension(NHENRY,CHEMTMIN:CHEMTMAX), save :: H
!hf NEW
real, private, dimension(CHEMTMIN:CHEMTMAX), save :: &
  K1,           & ! K for SO2->HSO3-
  K2,           & ! HSO3->SO32-
  Knh3,         & ! NH3+H20-> NH4+
  Kw,           & ! K for water
  Kco2
! Aqueous reaction rates for usage in gas-phase chemistry:
integer, private, parameter :: &
  NAQUEOUS = 5, & ! No. aqueous rates
  NAQRCT   = 3, & ! No. Temp. dependent rates
  S_H2O2   = 1,  & ! Numbering for Temp. dependent aq. rates
  S2_O3    = 2,  &
  S3_O3    = 3

real, private, save :: aqrcC_O3
real, private, dimension(NAQRCT,CHEMTMIN:CHEMTMAX), save :: aqrcT
                                                   ! Temp dependent rates 
                                                   ! for aq. so2 oxidn.

real, public, save,allocatable, dimension(:,:) :: aqrck
                                                 ! aq. reaction rates
						 ! transferred to chemistry
						 ! (CM_Reactions1(2).inc)

logical, public,save :: prclouds_present      ! true if precipitating clouds

integer, public, parameter :: &
  ICLOHSO2   = 1, & ! for [oh] + [so2]
  ICLRC1     = 2, & ! for [h2o2] + [so2]
  ICLRC2     = 3, & ! for [o3] + [so2]
  ICLRC3     = 4, & ! for [o3] + [o2] (Fe catalytic)
  ICLHO2H2O2 = 5    ! for HO2g --> 0.5 * [H2O2]

! Incloud scavenging: (only dependant on precipitation (not cloud water)
!-----------------------------------------------------------------------
! The parameterization of the scavenging of soluble chemical
! components is scaled to the precipitation in each layer. For
! incloud scavenging it is based on the parameterization described
! in Berge (1998). The incloud scavenging of a soluble component
! X is given by the expression:
!
!       X * f * pr_acc * W_sca
! Q =  ------------------------
!         Z_sca * rho_water
!
! where pr_acc is the accumulated precipitation in the layer,
! Z_sca is the scavenging depth (scale=1000m) and rho_water is the
! density of water.
! f is an efficiency parameter. The module My_WetDep_mod should
! return the user-defined array WetDep%W_sca, with W_Sca representing
! the value of f.W_Sca/Z_sca/rho_water


! Sub cloud scavenging:
!-----------------------------------------------------------------------
! The sub-cloud scavenging distinguishes between particulate and
! gas-phase components. The scavenging of gases is calculated as
! Q = vw * P, where P is the accumulated precipitation (pr_acc, m)
! from all the layers above. (From setup_1d)
! For particles the scavenging is believed to be much less effective,
! as they follow the air-current around the droplets (Berge, 1993).
! Scavenging for particles is calculated as
!    Q = A.e.P/v
! where A is 5.2 m3 kg-1 s-1, v is the fall-speed of the droplets,
! 5 m s-1, and e is the scavenging efficiency, 0.1.


! Create an array to map from the "calc" to the advected species
! Use zeroth column to store number of species in that row

  integer, private, save :: iwdepPMf 
! arrays for species and groups, e.g. SOX, OXN
integer, private, save :: nwgrpOut = 0, nwspecOut = 0  ! no. groups & specs
integer, private, allocatable, dimension(:), save :: wetGroupOut, wetSpecOut
type(group_umap), private, allocatable, dimension(:), target, save :: wetGroupOutUnits


integer, public, save  :: WDEP_PREC=-1   ! Used in Aqueous_mod
contains

subroutine Init_WetDep()
  integer :: iadv, igrp, icalc, n, nc, f2d, alloc_err
  character(len=30) :: dname


   call WetCoeffs()  ! Sets WDspec(:)% name, W_sca, W_sub

  allocate(incloud(KUPPER:KMAX_MID),pr_acc(KUPPER:KMAX_MID))

!####################### gather indices from My_Derived
! WDEP_WANTED array, and determine needed indices in d_2d

  iwdepPMf = find_index('PMf',WDspec(:)%name )
  nwspecOut=count(WDEP_WANTED(1:nOutputWDEP)%txt2=="SPEC")
  nwgrpOut =count(WDEP_WANTED(1:nOutputWDEP)%txt2=="GROUP")
  allocate(wetSpecOut(nOutputWdep),wetGroupOut(nwgrpOut),wetGroupOutUnits(nwgrpOut),stat=alloc_err)
  call CheckStop(alloc_err, "alloc error wetSpecOut/wetGroupOut")

  nwspecOut=0;nwgrpOut=0  ! resets
  do n = 1, nOutputWdep  ! size(WDEP_WANTED(:)%txt1)
    dname = "WDEP_"//trim(WDEP_WANTED(n)%txt1)
    f2d = find_index(dname,f_2d(:)%name)
    call CheckStop(f2d<1, "AQUEOUS f_2d PROBLEM: "//trim(dname))

    iadv=0;igrp=0
    select case(WDEP_WANTED(n)%txt2)
    case("PREC")
      WDEP_PREC=f2d
      if(WDEP_PREC>0) then
        iadv=-999;igrp=-999 ! just for printout
      elseif(DEBUG%AQUEOUS.and.MasterProc)then
        call CheckStop(WDEP_PREC,find_index(dname,f_2d(:)%name),&
          "Inconsistent WDEP_WANTED/f_2d definition for "//trim(dname))
      end if
    case("SPEC")
      iadv=f_2d(f2d)%index
      if(iadv>0) then
        nwspecOut = nwspecOut  + 1
        wetSpecOut(nwspecOut) = f2d
      elseif(DEBUG%AQUEOUS.and.MasterProc)then
        call CheckStop(iadv,find_index(dname,species_adv(:)%name),&
          "Inconsistent WDEP_WANTED/f_2d definition for "//trim(dname))
      end if
    case("GROUP")
      igrp=f_2d(f2d)%index
      if(igrp>0) then
        nwgrpOut = nwgrpOut + 1
        wetGroupOut(nwgrpOut) = f2d
        wetGroupOutUnits(nwgrpOut) = Group_Scale(igrp,f_2d(f2d)%unit,&
          debug=DEBUG%AQUEOUS.and.MasterProc)
      elseif(DEBUG%AQUEOUS.and.MasterProc)then
        call CheckStop(igrp,find_index(dname,chemgroups(:)%name),&
          "Inconsistent WDEP_WANTED/f_2d definition for "//trim(dname))
      end if
    end select

    if(DEBUG%AQUEOUS.and.MasterProc)  then
      write(*,"(2a,3i5)") "WETPPP ", trim(f_2d(f2d)%name), f2d, iadv, igrp
      if(igrp>0) write(*,*) "WETFGROUP ", nwgrpOut, wetGroupOutUnits(nwgrpOut)%iadv
      if(iadv>0) write(*,*) "WETFSPEC  ", nwspecOut, iadv
    end if
  end do

!####################### END indices here ##########

end subroutine Init_WetDep





!-----------------------------------------------------------------------
subroutine Setup_Clouds(i,j,debug_flag)
!-----------------------------------------------------------------------
! DESCRIPTION
! Called fron runchem
! Define incloud and precipitating clouds.
! The layer must contain at least 1.e-7 kgwater/kg air to
! be considered a cloud.
!
! Also calculates
!  pr_acc - the accumulated precipitation for each k
!  b      - fractional cloud cover for each k
!-----------------------------------------------------------------------
  integer, intent(in) ::  i,j
  logical, intent(in) :: debug_flag

  real, dimension(KUPPER:KMAX_MID) :: &
    b,           &  ! Cloud-area (fraction)
    cloudwater,  &  ! Cloud-water (volume mixing ratio, H2O/air)
                    ! cloudwater = 1.e-6 same as 1.g m^-3
    pres            ! Pressure (Pa)
  integer :: k, SO4


  SO4 = SO4_ix

! Add up the precipitation in the column:
!old defintion:
!  pr_acc(KUPPER) = sum ( pr(i,j,1:KUPPER) ) ! prec. from above
!  do k= KUPPER+1, KMAX_MID
!    pr_acc(k) = pr_acc(k-1) + pr(i,j,k)
!    pr_acc(k) = max( pr_acc(k), 0.0 )
!  end do

!now pr is already defined correctly (>=0)
  do k= KUPPER, KMAX_MID
    pr_acc(k) = pr(i,j,k)
    pres(k)=A_mid(k)+B_mid(k)*ps(i,j,1)
  end do

  prclouds_present=(pr_acc(KMAX_MID)>PR_LIMIT) ! --> precipitation at the surface

! initialise with .false. and 0:
  incloud(:)  = .false.
  cloudwater(:) = 0.
!  aqrck(:,:)=0.0 !set in setup_aqurates
!  aqrck(ICLOHSO2,:) = 1.0 !set in setup_aqurates

  !cc3d is in fraction units (values 0 to 1)
  !cw_met is in kg/kg units

! Loop starting at surface finding the cloud base:
  ksubcloud = KMAX_MID+1       ! k-coordinate of sub-cloud limit
  do k = KMAX_MID, KUPPER, -1
    ! roa(i,j,k,1) * cw_met(i,j,k,1) / b(k) * 1e3 > 0.06 checks that the in-cloud liquid water content is at least that of fog (0.06 g/m3)
    if(cc3d(i,j,k,1) > B_LIMIT .and. cw_met(i,j,k,1) > CW_LIMIT .and. roa(i,j,k,1) * cw_met(i,j,k,1) / b(k) * 1e3 > 0.06) exit
    ksubcloud = k
  end do



  if(ksubcloud /= KUPPER)then 
     !clouds were found under KUPPER
     
     ! Define incloud part of the column requiring that both cloud water
     ! and cloud fractions are above limit values
     kcloudtop = -1               ! k-level of cloud top


     do k = KUPPER, ksubcloud-1
        b(k) = cc3d(i,j,k,1)
        ! Units: kg(w)/kg(air) * kg(air(m^3) / density of water 10^3 kg/m^3
        ! ==> cloudwater (volume mixing ratio of water to air *in cloud*)
        ! (it is divided by cloud fraction b )

        if(b(k) > B_LIMIT .and. cw_met(i,j,k,1) > CW_LIMIT .and. roa(i,j,k,1) * cw_met(i,j,k,1) / b(k) * 1e3 > 0.06) then
           ! value of cloudwater in the cloud fraction of the grid in units
           ! vol(water)/vol(air)
           ! if  FIXED_CLW > 0, this value is used for clouds. Otherwise
           ! calculated from NWP values. (In future NWP will be used by default,
           ! but we are invesigating some pH calculation issues. 
           ! For safety, use FIXED_CLW
           if ( USES%FIXED_CLW > 0.0  ) then
             cloudwater(k) = USES%FIXED_CLW
           else ! calculate from NWP data
             cloudwater(k) = 1.0e-3 * roa(i,j,k,1) * cw_met(i,j,k,1) / b(k)
           end if 
           incloud(k) = .true.
           if(kcloudtop<0) kcloudtop = k
        end if
     end do  ! k inside cloud loop


     if(kcloudtop == -1) then
        if(prclouds_present.and.DEBUG%AQUEOUS) &
             write(*,"(a20,2i5,3es12.4)") "ERROR prclouds sum_cw", &
             i,j, maxval(cw_met(i,j,KUPPER:KMAX_MID,1),1), maxval(pr(i,j,:)), pr_acc(KMAX_MID)
        kcloudtop = KUPPER ! for safety
     end if
  else
     ! No cloud water found below KUPPER
     ! Cloud above KUPPER are likely thin
     ! cirrus clouds, and if included may
     ! need special treatment...
     ! ==> assume no cloud
      kcloudtop = -999 ! used as label for no cloud 
  endif   !  ksubcloud /= KUPPER



  
 ! sets up the aqueous phase reaction rates (SO2 oxidation) and the
 ! fractional solubility

 !need to be called also if no clouds for non-cloud rates
 !DSJ18 Query. Couldn't we just set AQRCK etc tozero if kcloudtop < 1



   call setup_aqurates(b ,cloudwater,incloud,pres)
!  Calculates arrays with temperature dependent aqueous phase
!  equilibrium rates and aqueous phase reaction rates


  if(kcloudtop>KUPPER .and.kcloudtop<KMAX_MID)then
     if(DEBUG%pH .and. debug_flag .and. incloud(kcloudtop)) then
        !   write(*,"(a,l1,2i4,es14.4)") "DEBUG_AQ ",prclouds_present, &
        !            kcloudtop, ksubcloud, pr_acc(KMAX_MID)

        write(*,*) "DEBUG%pH ",prclouds_present, &
             kcloudtop, ksubcloud, (pH(k),k=kcloudtop,ksubcloud-1)
        write(*,*) "CONC (mol/l)",&
             so4_aq(ksubcloud-1),no3_aq(ksubcloud-1),nh4_aq(ksubcloud-1),&
             nh3_aq(ksubcloud-1),hco3_aq(ksubcloud-1),co2_aq(ksubcloud-1)
        write(*,*)"H+(ph_factor) ",&
             hco3_aq(ksubcloud-1)+2.*so4_aq(ksubcloud-1)+hso3_aq(ksubcloud-1)&
             +2.*so32_aq(ksubcloud-1)+no3_aq(ksubcloud-1)-nh4_aq(ksubcloud-1)-nh3_aq(ksubcloud-1)
        write(*,*) "CLW(l_vann/l_luft) ",cloudwater(ksubcloud-1)
        if(SO4>0)write(*,*) "xn_2d(SO4) ugS/m3 ",(xn_2d(SO4,k)*10.e12*32./AVOG,k=kcloudtop,KMAX_MID)
     end if
  end if
     
end subroutine Setup_Clouds
!-----------------------------------------------------------------------





subroutine init_aqueous()
!-----------------------------------------------------------------------
! DESCRIPTION
! Called from main before time itegration 
! Calls initial tabulations, sets frac_aq to zero above cloud level, and
! sets constant rates.
! MTRLIM represents mass transport limitations between the clouds
! and the remainder of the grid-box volume. (so2 will be rapidly
! depleted within the clouds, and must be replenished from the
! surrounding cloudfree volume.
!-----------------------------------------------------------------------
!hf  real, parameter :: &
!hf     Hplus = 5.0e-5                    ! H+ (Hydrogen ion concentration)
!     h_plus = 5.0e-5                    ! H+ (Hydrogen ion concentration)
  real, parameter :: MASSTRLIM = 1.0   ! Mass transport limitation
!hf  INV_Hplus    = 1.0/Hplus             ! 1/H+
!hf  INV_Hplus0p4 = INV_Hplus**0.4        ! (1/H+)**0.4


  allocate(pH(KUPPER:KMAX_MID),so4_aq(KUPPER:KMAX_MID),no3_aq(KUPPER:KMAX_MID))
  allocate(nh4_aq(KUPPER:KMAX_MID),nh3_aq(KUPPER:KMAX_MID),hso3_aq(KUPPER:KMAX_MID))
  allocate(so2_aq(KUPPER:KMAX_MID),so32_aq(KUPPER:KMAX_MID),co2_aq(KUPPER:KMAX_MID))
  allocate(hco3_aq(KUPPER:KMAX_MID))
  allocate(frac_aq(NHENRY,KUPPER:KMAX_MID))
  allocate(aqrck(NAQUEOUS,KCHEMTOP:KMAX_MID))



! tabulations
!========================
  call tabulate_aqueous()
!========================

end subroutine init_aqueous
!-----------------------------------------------------------------------



subroutine tabulate_aqueous()
!-----------------------------------------------------------------------
! DESCRIPTION
! Tabulates Henry's law coefficients over the temperature range
! defined in Tabulations_mod.
! For SO2, the effective Henry's law is given by
!   Heff = H * ( 1 + K1/H+ )
! where k2 is omitted as it is significant only at high pH.
! We tabulate also the factor 1+K1/H+ as K1fac.

! Calculates temperature aqueoua equilibrium and reaction rates in
! for temeratures between CHEMTMIN (148) and CHEMTMAX (333) as defined in
! Config_module
! Constant rates: The rates are in mol-1 l. See code for references.
! These need to be multiplied by 1.0e3/AVOG/Vf,so we perform the
! 1.0e3/AVOG scaling here.
!-----------------------------------------------------------------------


  real, parameter :: MASSTRLIM = 1.0   ! Mass transport limitation
  real, dimension(CHEMTMIN:CHEMTMAX) :: t, tfac ! Temperature, K, factor
  integer :: i
  real    :: to_mol_per_l

  t(:)          = (/ ( real(i), i=CHEMTMIN, CHEMTMAX ) /)
  tfac(:)       = 1.0/t(:) -  1.0/298.0
  to_mol_per_l  = 1.0e3 * MASSTRLIM /AVOG

!  Henry's law constants
  H(IH_SO2 ,:)  = 1.3    * exp(2900.0*tfac(:))  ! Sander(2023)
  H(IH_H2O2,:)  = 8.7e4   * exp(7300.0*tfac(:)) ! Sander(2023)
  H(IH_O3  ,:)  = 1.01e-2 * exp(2800.0*tfac(:)) ! Sander(2023)
  H(IH_NH3 ,:)  = 60.0    * exp(4200.0*tfac(:)) ! Sander(2023)
  H(IH_CO2,:)   = 3.4e-2  * exp(2300.0*tfac(:)) ! Sander(2023)


! Equilibrium ractions. These are further used to calculate
! efficient Henry's law constants. 
  K1(:)   = 1.3e-2  * exp( 1960.0*tfac(:))  ! Seinfeldt&Pandis 2016
  K2(:)   = 6.6e-8  * exp( 1500.0*tfac(:))  ! Seinfeldt&Pandis 2016
  Knh3(:) = 1.7e-5  * exp(-450.0*tfac(:))   ! Seinfeldt&Pandis 2016
  Kw(:)   = 1.0e-14 * exp(-6718.0*tfac(:))  ! Seinfeldt&Pandis 2016
  Kco2(:) = 4.3e-7  * exp(-1000.0*tfac(:))  ! Seinfeldt&Pandis 2016


!! Temperature dependent aqueous phase reaction rates
  aqrcT(S_H2O2,:) = 7.5e7 * exp(-4430.0*tfac(:)) * to_mol_per_l  ! Seinfeldt&Pandis 2016
  aqrcT(S2_O3,:)  = 3.7e5 * exp(-5530.0*tfac(:)) * to_mol_per_l  ! Seinfeldt&Pandis 2016
  aqrcT(S3_O3,:)  = 1.5e9 * exp(-5280.0*tfac(:)) * to_mol_per_l  ! Seinfeldt&Pandis 2016

! SO2aq + O3 not temp. dependent
   aqrcC_O3 = 2.4e4 *  to_mol_per_l   ! Seinfeldt&Pandis 2016


! Old Temperature independent aqueous phase reaction rates
!  aqrcT(S_H2O2,:) = 8.3e5   !!!!7.5e7 * exp(-4430.0*tfac(:))  ! Seinfeldt&Pandis 2016
!  aqrcT(S2_O3,:)  = 3.7e5 * exp(-5530.0*tfac(:))  ! Seinfeldt&Pandis 2016
!  aqrcT(S3_O3,:)  = 1.5e9 * exp(-5280.0*tfac(:))  ! Seinfeldt&Pandis 2016
!!!!!  aqrcT(S1_O3,:) = !!! not Temp. dependent


end subroutine tabulate_aqueous
!-----------------------------------------------------------------------





subroutine setup_aqurates(b ,cloudwater,incloud,pres)
!-----------------------------------------------------------------------
! DESCRIPTION
! called from setup_clouds (witch is called from runchem)
! First pH dependence is calculated in an iterative method.
! When pH the (ph dependent) rate-coefficients for the
! aqueous-phase reactions can be set. 
! 
!-----------------------------------------------------------------------
 !DSJ18 Query. Couldn't we just set AQRCK etc tozero if kcloudtop < 1
  real, parameter :: MASSTRLIM = 1.0   ! Mass transport limitation

  real, dimension(KUPPER:KMAX_MID) :: &
    b,          &  ! cloud-aread (fraction)
    cloudwater, &  ! cloud-water
    pres           ! Pressure(Pa) !hf

  logical, dimension(KUPPER:KMAX_MID) :: &
    incloud  ! True for in-cloud k values

! Outputs  -> aqurates

! local
  real, dimension(KUPPER:KMAX_MID) :: &
    fso2grid, & ! f_aq * b = {f_aq}
    fso2aq,   & ! only so2.h2o part (not hso3- and so32-)
    fhso3,    & ! only hso3.h2o part (not so2aq and so32-)
    fso3        ! only so2.h2o part (not so2aq and hso3-)

! PH
  real, dimension(KUPPER:KMAX_MID) :: &
    phfactor, &
    h_plus

  real, parameter :: CO2conc_ppm = 392 !mix ratio for CO2 in ppm
  real :: CO2conc !Co2 in mol/l
  real :: invhplus, invhplus04, K1K2_fac, Heff,Heff_NH3,pH_old
  real, dimension (KUPPER:KMAX_MID) :: VfRT ! Vf * Rgas * Temp
  real, parameter :: Hplus43=5.011872336272724E-005! 10.0**-4.3
  real, parameter :: Hplus55=3.162277660168379e-06! 10.0**-5.5
  integer k, iter
  integer :: SO2, SO4, NO3_f, NH4_f, NH3, HNO3
  real :: pH1,pH2,pHnew,pHout1,pHout2,pHoutnew,h_plusnew
  
  integer, parameter:: N_ITER=40 !each iteration divides the error in two
  real, parameter :: phThres = 1.0E-6 !just to verify result.
  
  SO2 = SO2_ix
  call CheckStop( SO2_ix<1, "Aqueous: SO2 not defined" )
  SO4 = SO4_ix
  call CheckStop( SO4_ix<1, "Aqueous: SO4 not defined" )
  NH3 = NH3_ix
  !call CheckStop( NH3_ix<1, "Aqueous: NH3 not defined" )
  NH4_f = NH4_f_ix
  !call CheckStop( NH4_f_ix<1, "Aqueous: NH4_f not defined" )
  NO3_f = NO3_f_ix
  call CheckStop( NO3_f_ix<1, "Aqueous: NO3_f not defined" )
  HNO3 = HNO3_ix
  call CheckStop( HNO3_ix<1, "Aqueous: HNO3 not defined" )


! pw note: some frac_aq are recalculated (overwritten) after get_frac
  call get_frac(cloudwater,incloud) ! => frac_aq

! initialize:
  aqrck(:,:) = 0.

! for PH

  pH(:)=4.3!dspw 13082012
  h_plus(:)=Hplus43!dspw 13082012
  pH(:)=5.5!stpw 23082012
  h_plus(:)=Hplus55!stpw 23082012

  ph_old=0.0
! Gas phase ox. of SO2 is "default"
! in cloudy air, only the part remaining in gas phase (not
! dissolved) is oxidized
  aqrck(ICLOHSO2,:) = 1.0
  Fgas(SO2,:) = 1.0

  do k = KUPPER,KMAX_MID
     if(.not.incloud(k)) cycle ! Vf > 1.0e-10)

     !  pwjoj (Jan 2023):
     !  we dissolve only the part which is inside the cloud,
     !  So the concentrations on the lhs are concentrations
     !  within the cloud (0 otherwise) 

!For pH calculations:
!   Assume total uptake of so4,no3,hno3,nh4+
!  For pH below 5, all NH3 will be dissolved, at pH=6 around 50%
!  Effectively all dissolved NH3 will ionize to NH4+ (Seinfeldt)
!  cloudwater volume mix. ratio
!  SO4_aq, NO3_aq, NH4_aq converted from molec cm-3 to mol/l   

   so4_aq(k)= (xn_2d(SO4,k)*1000./AVOG)/cloudwater(k) !xn_2d=molec cm-3
   no3_aq(k)= ( (xn_2d(NO3_F,k)+xn_2d(HNO3,k))*1000./AVOG)/cloudwater(k)
   if (NH4_F>0) then
       nh4_aq(k) =  ( xn_2d(NH4_F,k) *1000./AVOG )/cloudwater(k)  !only nh4+ now
    else
       nh4_aq(k) =  0.0
    end if

!!  old stuff
 !   hso3_aq(k)= 0.0 !initial, before dissolved
 !   so32_aq(k)= 0.0
 !   nh3_aq(k) = 0.0 !nh3 dissolved and ionized to nh4+(aq)
 !   hco3_aq(k) = 0.0 !co2 dissolved and ionized to hco3

    VfRT(k) = cloudwater(k) * RGAS_ATML * temp(k)

    !  dissolve CO2 and SO2 (pH independent)
    !  CO2conc=392 ppm
    CO2conc=CO2conc_ppm * 1e-9 * pres(k)/(RGAS_J *temp(k)) !mol/l


    frac_aq(IH_CO2,k) = 1.0 / ( 1.0+1.0/( H(IH_CO2,itemp(k))*VfRT(k) ) )
    co2_aq(k)=frac_aq(IH_CO2,k)*CO2CONC /cloudwater(k)

    frac_aq(IH_SO2,k) = 1.0 / ( 1.0+1.0/( H(IH_SO2,itemp(k))*VfRT(k) ) )
    so2_aq(k)= frac_aq(IH_SO2,k)*(xn_2d(SO2,k)*1000./AVOG)/cloudwater(k)


!!!  pH iteration starts here  !!!!
!!!! ------------------------  !!!!
    !Iterate between min and max value. St
    !Iterate between min and max value. Start with extremes:
    pH1=1.0
    pH2=7.0

    do iter = 1, N_ITER! binary search
       if (iter==1) then
          pHnew=pH1
       else if (iter==2) then
          pHnew=pH2
       else                  
          pHnew=(pH1+pH2)/2
          !Alternative, faster convergence, but the function may not be stable near 7 (why?).
          !make a straight line between values, 
          !and find where the line cross the diagonal, i.e. pH = pHout
          !pHnew=(pH1*pHout2-pH2*pHout1)&
          !     /(pH1-pH2-pHout1+pHout2)
       end if
       h_plusnew=exp(-pHnew*log(10.))
       
       Heff_NH3= H(IH_NH3,itemp(k))*Knh3(itemp(k))*h_plusnew/Kw(itemp(k))
       frac_aq(IH_NH3,k) = 1.0 / ( 1.0+1.0/( Heff_NH3*VfRT(k) ) )
       if (NH3>0) then     
          nh3_aq(k) = frac_aq(IH_NH3,k)*(xn_2d(NH3,k)*1000./AVOG)/cloudwater(k)
       else
          nh3_aq(k) = 0.0
       end if
       
       hco3_aq(k)= co2_aq(k) * Kco2(itemp(k))/h_plusnew
       hso3_aq(k)= so2_aq(k) * K1(itemp(k))/h_plusnew
       so32_aq(k)= hso3_aq(k) * K2(itemp(k))/h_plusnew
       
       phfactor(k)=hco3_aq(k)+2.*so4_aq(k)+hso3_aq(k)+2.*so32_aq(k)+no3_aq(k)-nh4_aq(k)-nh3_aq(k)
       h_plusnew=0.5*(phfactor(k) + sqrt(phfactor(k)*phfactor(k)+4.*1.e-14) )
       h_plusnew=min(1.e-1,max(h_plusnew,1.e-7))! between 1 and 7
       pHoutnew=-log(h_plusnew)/log(10.)

       !commented out: fixed number of iterations
       !       if(abs(pHnew-pHoutnew)< phThres)then       
       !we have converged enough
       !exit
       !       end if       
 
       if (iter==1) then
          pHout1=pHoutnew
       else if (iter==2) then
          pHout2=pHoutnew
       else
          !we replace either pH1 or pH2. They must be on each side of the diagonal
          if((pH1-pHout1)*(pH2-pHout2)<=0)then
             if((pH1-pHout1)*(pHnew-pHoutnew)<=0)then
                !we keep pH1 and the new pH
                pH2=pHnew
                pHout2=pHoutnew
             else if((pH2-pHout2)*(pHnew-pHoutnew)<=0)then
                !we keep pH2 and the new pH
                pH1=pHnew
                pHout1=pHoutnew
             end if
          else
             !something went wrong
              if(abs(pHnew-pHoutnew)<pHThres)then
                 !we are close enough to convergence
              else
                 pH(k)=5.5
                 write(*,*)iter,'Warning: incompatible pH ',pH1,pH2,pHout1,pHout2,pHnew,pHoutnew
              end if
              exit
           end if
           if(iter==N_ITER)then
             !check that we are close to the solution
             if(abs(pHnew-pHoutnew)>pHThres)then
                write(*,*)'Warning: pH did not converge properly ',pH1,pHout1,pH2,pHout2,pHnew,pHoutnew
             end if
          end if
        endif       
    end do  ! iter   Iteration completed and concentration derived pH
            !        and H+ used below to calculate pH dependent
	    !        effective Henry's law constants (that will also
	    !        affect aqueous oxidation rates.
    pH(k)=pHnew
    h_plus(k)=exp(-pH(k)*log(10.))

!!!!!!  pH iteration fiished   !!!!!!!!!



!after pH determined, final numbers of frac_aq(IH_SO2)
!= effective fraction of S(IV):
!include now also ionization to SO32-
!        K1_fac  =  &
!             1.0 + K1(k)/h_plus(k) !not used
!        H (IH_SO2 ,itemp(k))  = H(IH_SO2,itemp(k)) * K1_fac
    invhplus   = 1.0/h_plus(k)
    invhplus04 = (1.0/h_plus(k))**0.4

    K1K2_fac = 1.0  &   
             + K1(itemp(k)) * invhplus &
	     + K1(itemp(k))*K2(itemp(k)) * (invhplus**2)

    Heff  = H(IH_SO2,itemp(k)) * K1K2_fac

    frac_aq(IH_SO2,k) = 1.0 / ( 1.0+1.0/( Heff*VfRT(k) ) )
    
    fso2grid(k) = b(k) * frac_aq(IH_SO2,k)   !  frac of S(IV) in grid
                                             !  in aqueous phas - Saq/(Saq + Sg)
    fso2aq  (k) = fso2grid(k) / K1K2_fac     ! frac of SO2 in total grid
                                             ! in aqueous phase

    fhso3(k)    = fso2aq(k) * K1(itemp(k)) * invhplus  
                                             ! frac of HSO3- in total grid
                                             ! in aqueous phase

    fso3  (k) = fhso3(k) * K2(itemp(k)) * invhplus 
                                             ! frac of SO3-- in total grid
                                             ! in aqueous phase





!!!!   SET  Aquous phase reaction rates exported to chemistry 
!!-------------------------------------------------------------


  ! oh + so2 gas-phase
    aqrck(ICLOHSO2,k) = ( 1.0-fso2grid(k) ) ! now correction factor!
    Fgas(SO2,k) = 1.0-fso2grid(k)
                                            ! as part of SO2 not in gas phase
  !  aqrck(ICLOHSO2,k) = ( 1.0-fso2grid(k) ) * &
  !    IUPAC_TROE(2.8e-31*EXP(2.6*(LOG(300/temp(k)))),2.0e-12,EXP(-temp(k)/472.0),M(k),0.75-1.27*(-temp(k)/472.0)/LOG(10.0))

  !  HO2g ----> 0.5 H2O2
    aqrck(ICLHO2H2O2,k) = 0.066 * cloudwater(k)* 1.e6 * b(k)


!  Incloud oxidation of Siv to Svi by H2O2
!  When calculating fhso3 devided by Hpluss  ==> divide and multiply by Hpluss cancel out
    aqrck(ICLRC1,k)   =  aqrcT(S_H2O2,itemp(k)) * frac_aq(IH_H2O2,k) * fhso3(k) &
                       * H_plus(k) / cloudwater(k)                   !  only HSO3-
                                                                      ! skip 1/(1+Hplus*K), see S&P 2016
                                                                      ! as only significant pH < 2
!!  Incloud oxidation of Siv to Svi by O3
    aqrck(ICLRC2,k)   = ( aqrcC_O3 * fso2aq(k)    &  ! with SO2aq
                      +   aqrcT(S2_O3,itemp(k)) * fhso3(k)     &  ! with HSO3-
		      +   aqrcT(S3_O3,itemp(k)) * fso3(k) )    &  ! with SO3--
                      *   frac_aq(IH_O3,k)/cloudwater(k) 

!    aqrck(ICLRC2,k)   = 1.8e4 * 1.0e3 * frac_aq(IH_O3,k) &
!                      * invhplus04 * fso2grid(k) / cloudwater(k)/ AVOG * MASSTRLIM
!!!!!! (so2aq + hso3-) + H+ + o3 ---> so4, ref: Martin & Damschen 1981
!!!!!!  aqrc(2) = 1.8e4 * 1.0e3/AVOG * MASSTRLIM


!  Incloud oxidation of Siv to Svi by O2 
    aqrck(ICLRC3,k)   = 3.3e-10 *  fso2grid(k)/ cloudwater(k) * MASSTRLIM

!!!!!! (so2aq + hso3-) + o2 ( + Fe ) --> so4, see documentation below
!!!!!!!  aqrc(3) = 3.3e-10  * MASSTRLIM

! Regarding reaction with O2:
! catalytic oxidation with Fe. The assumption is that 2% of SIV
! is oxidised per hour inside the droplets, corresponding to a
! conversion rate of 5.6^-6 (units s^-1  -- Therfore no conversion
! from mol l^-1)


!!!Old code
!!!!    caqh2o2 (k) = aqrc(1) * frac_aq(IH_H2O2,k) / cloudwater(k)
!    aqrck(ICLRC1,k)   = caqh2o2(k) * fso2aq(k) !only SO2
!    !        aqrck(ICLRC2,k)   = caqo3(k) * INV_Hplus0p4 * fso2grid(k)

!  Incloud oxidation of Siv to Svi by O3
!!!    caqo3   (k) = aqrc(2) * frac_aq(IH_O3,k) / cloudwater(k)
!    aqrck(ICLRC2,k)   = caqo3(k) * invhplus04 * fso2grid(k)

!  Incloud oxidation of Siv to Svi by O2 
!!!    caqsx   (k) = aqrc(3) / cloudwater(k)
!    aqrck(ICLRC3,k)   = caqsx(k) *  fso2grid(k)

  end do

end subroutine setup_aqurates
!-----------------------------------------------------------------------




subroutine get_frac(cloudwater,incloud)
!-----------------------------------------------------------------------
! DESCRIPTION
! Calculating pH dependant solubility fractions: Calculates the fraction
! of each soluble gas in the aqueous phase, frac_aq
!-----------------------------------------------------------------------
! intent in from used modules  : cloudwater and logical incloud
! intent out to rest of module : frac_aq

! local
  real, dimension (KUPPER:KMAX_MID), intent(in) :: &
      cloudwater    ! Volume fraction - see notes above.
  logical, dimension(KUPPER:KMAX_MID), intent(in) :: &
      incloud       ! True for in-cloud k values
  real    :: VfRT   ! Vf * Rgas * Temp
  integer :: ih, k  ! index over species with Henry's law, vertical level k

! Make sure frac_aq is zero outside clouds:
  frac_aq(:,:) = 0.0

  do k = KUPPER, KMAX_MID
    if(.not.incloud(k)) cycle

    VfRT = cloudwater(k) * RGAS_ATML * temp(k)
! Get aqueous fractions:
    do ih = 1, NHENRY
      frac_aq(ih,k) = 1.0 / ( 1.0+1.0/( H(ih,itemp(k))*VfRT ) )
    end do
  end do
end subroutine get_frac
!-----------------------------------------------------------------------




subroutine WetDeposition(i,j,debug_flag)
!-----------------------------------------------------------------------
! DESCRIPTION
! Calculates wet deposition and changes in xn concentrations
! WetDeposition called from RunChem if precipitation reach the surface
!-----------------------------------------------------------------------
! input
  integer, intent(in) ::  i,j
  logical, intent(in) :: debug_flag

! local
  integer :: itot,iadv,is !  index in xn_2d arrays
  integer :: k, icalc

  real    :: invgridarea ! xm2/(h*h)
  real    :: f_rho       ! Factors in rho calculation
  real    :: rho(KUPPER:KMAX_MID)
  real    :: loss    ! conc. loss due to scavenging
  real, dimension(KUPPER:KMAX_MID) :: vw ! Svavenging rates (tmp. array)
  real, dimension(KUPPER:KMAX_MID) :: lossfac ! EGU
  real, dimension(KUPPER:KMAX_MID) :: lossfacPMf ! for particle fraction of semi-volatile (VBS) species

  wdeploss(:) = 0.0
  if(WDEP_PREC>0)d_2d(WDEP_PREC,i,j,IOU_INST) = pr(i,j,KMAX_MID) * dt ! Same for all models
  if ( kcloudtop < 1 ) RETURN !  skip wetdep calcs if no cloud

  invgridarea = xm2(i,j)/( gridwidth_m*gridwidth_m )
  f_rho  = 1.0/(invgridarea*GRAV*ATWAIR)
! Loop starting from above:
  do k=kcloudtop, KMAX_MID           ! No need to go above cloudtop
    rho(k) = f_rho*(dA(k) + dB(k)*ps(i,j,1))/ M(k)
  end do

! calculate concentration after wet deposition and sum up the vertical
! column of the depositions for the fully soluble species.

  if(DEBUG%AQUEOUS.and.debug_flag) write(*,*) "(a15,2i4,es14.4)", &
     "DEBUG_WDEP2", kcloudtop, ksubcloud, pr_acc(KMAX_MID)

! need particle fraction wet deposition for semi-volatile species - here hard
!  coded to use scavenging parameters for PMf

  vw(kcloudtop:ksubcloud-1) = WDspec(iwdepPMf)%W_sca ! Scav. for incloud
  vw(ksubcloud:KMAX_MID  )  = WDspec(iwdepPMf)%W_sub ! Scav. for subcloud
  do k = kcloudtop, KMAX_MID
    lossfacPMf(k)  = exp( -vw(k)*pr_acc(k)*dt )
  enddo 

  do icalc = 1, nwdep  ! Here we loop over "model" species

! Put both in- and sub-cloud scavenging ratios in the array vw:

    vw(kcloudtop:ksubcloud-1) = WDspec(icalc)%W_sca ! Scav. for incloud
    vw(ksubcloud:KMAX_MID  )  = WDspec(icalc)%W_sub ! Scav. for subcloud

    do k = kcloudtop, KMAX_MID
      lossfac(k)  = exp( -vw(k)*pr_acc(k)*dt )

      ! For each "calc" species we have often a number of model
      ! species

      do is = 1, size(WDmapping(icalc)%advspecs)  ! number of species
        iadv = WDmapping(icalc)%advspecs(is)
        itot = iadv+NSPEC_SHL

    ! For semivolatile species only the particle fraction is deposited
    !RB: This assumption needs to be revised. The semi-volatile organics are
    ! likely highly soluble and should wet deposit also in the gas phase
        if(itot>=FIRST_SEMIVOL .and. itot<=LAST_SEMIVOL) then
          loss = xn_2d(itot,k) * ( Fpart(itot,k)*( 1.0 - lossfacPMf(k) ) &
                + (1.0-Fpart(itot,k))*( 1.0 - lossfac(k) ) )
        else
          loss = xn_2d(itot,k) * ( 1.0 - lossfac(k)  )
        endif
        xn_2d(itot,k) = xn_2d(itot,k) - loss
        loss = loss * rho(k) ! concentration -> weight
        wdeploss(iadv) = wdeploss(iadv) + loss
        
        if(USES%LocalFractions .and. wetdep_lf(iadv)) call lf_wetdep(iadv, i, j, k, loss, invgridarea)

        if(DEBUG%AQUEOUS.and.debug_flag.and.pr_acc(KMAX_MID)>1.0e-5) then
          write(*,"(a50,2i4,a,9es12.2)") "DEBUG_WDEP, k, icalc, spec", k, &
           icalc, trim(WDspec(icalc)%name)//':'//trim(species_adv(iadv)%name),&
            vw(k), pr_acc(k), lossfac(k)
        end if ! DEBUG%AQUEOUS
      enddo ! is
    enddo ! k loop

  end do ! icalc loop

  ! add other losses into twetdep and wdep arrays:

  call WetDep_Budget(i,j,invgridarea,debug_flag)

 end subroutine WetDeposition
!-----------------------------------------------------------------------
 subroutine WetDep_Budget(i,j,invgridarea, debug_flag)
  integer, intent(in) ::  i,j
  real,    intent(in) :: invgridarea
  logical, intent(in) :: debug_flag

  integer :: f2d, igrp ,iadv, n, g
  real    :: wdep
  type(group_umap), pointer :: gmap=>null()  ! group unit mapping

  ! Mass Budget: Do not include values on outer frame
  if(.not.(i<li0.or.i>li1.or.j<lj0.or.j>lj1)) &
    totwdep(:) = totwdep(:) + wdeploss(:)

  ! Deriv.Output: individual species (SO4, HNO3, etc.) as needed
  do n = 1, nwspecOut
    f2d  = wetSpecOut(n)
    iadv = f_2d(f2d)%index
    d_2d(f2d,i,j,IOU_INST) = wdeploss(iadv) * invgridarea


    if(DEBUG%MY_WETDEP.and.debug_flag) &
      call datewrite("WET-PPPSPEC: "//species_adv(iadv)%name,&
        iadv,(/wdeploss(iadv)/))
  end do

  ! Deriv.Output: groups of species (SOX, OXN, etc.) as needed
  do n = 1, nwgrpOut
    f2d  =  wetGroupOut(n)
    gmap => wetGroupOutUnits(n)
    igrp = f_2d(f2d)%index

    wdep = dot_product(wdeploss(gmap%iadv),gmap%uconv(:))
    d_2d(f2d,i,j,IOU_INST) = wdep * invgridarea
    if(DEBUG%MY_WETDEP.and.debug_flag)then
      do g=1,size(gmap%iadv)
        iadv=gmap%iadv(g)
        call datewrite("WET-PPPGROUP: "//species_adv(iadv)%name ,&
          iadv,(/wdeploss(iadv)/))
      end do
    end if
  end do
 end subroutine WetDep_Budget
!-----------------------------------------------------------------------
end module Aqueous_mod
