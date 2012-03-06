! <Aqueous_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Aqueous_ml

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
! Setup_Clouds(i,j)  called from Runchem_ml
! WetDeposition(i,j) called from Runchem_ml if prec. clouds are present
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

  use My_Derived_ml,     only : WDEP_WANTED ! Which outputs wanted!

  use CheckStop_ml, only : CheckStop
  use ChemChemicals_ml, only: species
  use ChemSpecs_tot_ml
  use ChemSpecs_adv_ml          ! IXADV_SO2, IXADV_SO4, etc.
  use ChemGroups_ml, only : ChemGroups, INDEX_WDEP_SOX_GROUP, &
         INDEX_WDEP_RDN_GROUP, INDEX_WDEP_OXN_GROUP
  use DerivedFields_ml, only : f_2d, d_2d     ! Contains Wet deposition fields
  use GridValues_ml, only : gridwidth_m,xm2,dA,dB
  use Io_ml,             only : IO_DEBUG, datewrite
  use MassBudget_ml,     only : wdeploss,totwdep
  use ModelConstants_ml, only: &
      CHEMTMIN, CHEMTMAX       &       ! -> range of temperature 
     ,MasterProc               &
     ,DEBUG => DEBUG_AQUEOUS   &       !
     ,atwS, atwN, DEBUG_MY_WETDEP &
     ,DEBUG_pH                 &
     ,KMAX_MID                 &       ! -> ground, k=20
     ,KUPPER                   &       ! -> top of cloud-chemistry, k=6
     ,KCHEMTOP                 &       ! -> top of chemistry, now k=2
     ,dt => dt_advec           &       ! -> model timestep
     ,IOU_INST                 &       ! Index: instantaneous values
     ,ATWAIR                           ! -> atw. air
  use MetFields_ml,              only : pr, roa, z_bnd, cc3d, lwc,cw
  use MetFields_ml,              only : ps
  use OrganicAerosol_ml,    only: ORGANIC_AEROSOLS
  use OwnDataTypes_ml,           only : depmap, typ_i3  ! has adv, calc, vg
  use PhysicalConstants_ml, only: GRAV  &   
                                 ,AVOG  &    ! Avogadro's No.
                                 ,RGAS_ATML,RGAS_J  ! Gas-constant
  use Setup_1dfields_ml,   only : xn_2d, amk, Fpart, Fgas &
                                 ,temp  &    ! temperature (K)
                                 ,itemp      ! temperature (K)
  use SmallUtils_ml,  only : find_index




  implicit none
  private

! Subroutines:

  public :: Init_WetDep         ! Call from Unimod
  public :: WetDep_Budget       ! called here

  public  :: init_aqueous
  public  :: Setup_Clouds       ! characterises clouds and calls
                                ! WetDeposition if rain
  public  :: WetDeposition      ! simplified setup_wetdep
  private :: tabulate_aqueous
  private :: get_frac
  private :: setup_aqurates

! Outputs:

  logical, public, save, dimension(KUPPER:KMAX_MID) :: &
           incloud              ! True for in-cloud k values 


! Variables used in module:

  real, private, save, dimension(KUPPER:KMAX_MID) :: &
        pr_acc                  ! Accumulated precipitation
!hf NEW (here for debugging)
  real, private, save, dimension(KUPPER:KMAX_MID) :: &
        pH,so4_aq,no3_aq,nh4_aq,nh3_aq,hso3_aq,so2_aq,so32_aq,co2_aq,hco3_aq   ! pH in cloud

  integer, private, save  :: kcloudtop   ! k-level of highest-cloud
  integer, private, save  :: ksubcloud   ! k-level just below cloud

  real, private, parameter :: & ! Define limits for "cloud"
      PR_LIMIT = 1.0e-7  &      ! for accumulated precipitation
     ,CW_LIMIT = 1.0e-10 &      ! for cloud water, kg(H2O)/kg(air)
     ,B_LIMIT  = 1.0e-3         ! for cloud cover (fraction)

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

  real, public, dimension ( NHENRY, KUPPER:KMAX_MID ), save  :: frac_aq
  real, private, dimension (NHENRY, CHEMTMIN:CHEMTMAX), save :: H
  real, private, dimension (NK1,CHEMTMIN:CHEMTMAX), save     :: K1fac
!hf NEW
  real, private, dimension (CHEMTMIN:CHEMTMAX), save :: K1 ! K for SO2->HSO3-
  real, private, dimension (CHEMTMIN:CHEMTMAX), save :: K2 ! HSO3->SO32-
  real, private, dimension (CHEMTMIN:CHEMTMAX), save :: Knh3 ! NH3+H20-> NH4+
  real, private, dimension (CHEMTMIN:CHEMTMAX), save :: Kw   ! K for water
  real, private, dimension (CHEMTMIN:CHEMTMAX), save :: Kco2
! Aqueous reaction rates for usage in gas-phase chemistry:

  integer, private, parameter :: NAQUEOUS = 4   ! No. aqueous rates
  integer, private, parameter :: NAQRC    = 3   ! No. constant rates

  real, public, dimension(NAQUEOUS,KCHEMTOP:KMAX_MID), save :: aqrck

  real, private, dimension(NAQRC), save :: aqrc ! constant rates for
                                                ! so2 oxidn.
  real, private, dimension(2), save :: vw       ! constant rates for
  logical, public,save :: prclouds_present      ! true if precipitating
                                                ! clouds

  integer, public, parameter :: &
      ICLOHSO2  = 1   &  ! for [oh] + [so2]
     ,ICLRC1    = 2   &  ! for [h2o2] + [so2]
     ,ICLRC2    = 3   &  ! for [o3] + [so2]
     ,ICLRC3    = 4      ! for [o3] + [o2] (Fe catalytic)


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
! f is an efficiency parameter. The module My_WetDep_ml should
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



  ! ------------ WetDep initialisation (old My_WetDep --------------

  type, public :: WScav
     real  :: W_sca       ! Scavenging ratio/z_Sca/rho = W_sca/1.0e6
     real  :: W_sub       ! same for subcloud
  end type WScav
  

  integer, public, parameter :: NWETDEP_CALC =  14 ! No. of solublity classes

   !  Note - these are for "master" or model species - they do not
   !  need to be present in the chemical scheme. However, the chemical
   !  scheme needs to define wet scavenging after these. If you would
   !  like other characteristics, add them here.
  
  integer, parameter, public :: &
     CWDEP_SO2 = 1,  &
     CWDEP_SO4 = 2,  &
     CWDEP_NH3 = 3,  &
     CWDEP_HNO3 = 4,  &
     CWDEP_H2O2 = 5,  &
     CWDEP_HCHO = 6,  &
     CWDEP_PMf  = 7,  &
     CWDEP_PMc  = 8,  &
     CWDEP_ECfn = 9,  &
     CWDEP_SSf  = 10, &
     CWDEP_SSc  = 11, &
     CWDEP_SSg  = 12, &
     CWDEP_POLLw = 13, &
     CWDEP_ROOH = 14   ! TEST!!



            !===========================================!
            ! Chemistry-dependent mapping:
            ! WdepMap = (/
            !    depmap( HNO3, CWDEP_HNO3, -1) & etc.
            ! .... produced from GenChem, also with e.g. 
            !integer, public, parameter ::  NWETDEP_ADV  = 14
            !===========================================!

                       include 'CM_WetDep.inc'

            !===========================================!
            !===========================================!
            !===========================================!


 ! And create an array to map from the "calc" to the advected species
 ! Use zeroth column to store number of species in that row

   integer, public, dimension(NWETDEP_CALC,0:NWETDEP_ADV) :: Calc2tot

  ! arrays for species and groups, e.g. SOX, OXN
    integer, save, private :: nwgrp = 0, nwspec = 0  ! no. groups & specs
    integer, save, allocatable, dimension(:), private :: wetgroup, wetspec 

    type(typ_i3), save, private, &
          dimension(size(WDEP_WANTED(:)%txt1 )) :: tmpgroup, tmpspec

    type(WScav), public, dimension(NWETDEP_CALC), save  :: WetDep
  
    integer, public, save  :: WDEP_PREC   ! Used in Aqueous_ml


contains

  subroutine Init_WetDep()

    integer :: itot, icalc, n, nc, if2, igr, isp, alloc_err, atw

  !/ SUBCLFAC is A/FALLSPEED where A is 5.2 m3 kg-1 s-1, 
  !  and the fallspeed of the raindroplets is assumed to be 5 m/s. 
    real, parameter ::  FALLSPEED = 5.0               ! m/s 
    real, parameter ::  SUBCLFAC = 5.2 / FALLSPEED

  !/ e is the scavenging efficiency (0.02 for fine particles, 0.4 for course)

    real, parameter ::  EFF25 = 0.02*SUBCLFAC  & 
                      , EFFCO = 0.4*SUBCLFAC   &
                      , EFFGI = 0.7*SUBCLFAC  

   !/.. setup the scavenging ratios for in-cloud and sub-cloud. For
   !    gases, sub-cloud = 0.5 * incloud. For particles, sub-cloud=
   !    efficiency * SUBCLFAC
   !/..                             W_Sca  W_sub
    WetDep(CWDEP_SO2)   = WScav(   0.3,  0.15)  ! Berge+Jakobsen
    WetDep(CWDEP_SO4)   = WScav(   1.0,  EFF25) ! Berge+Jakobsen
    WetDep(CWDEP_NH3)   = WScav(   1.4,  0.5 )  ! subcloud = 1/3 of cloud for gases
    WetDep(CWDEP_HNO3)  = WScav(   1.4,  0.5)   ! 
    WetDep(CWDEP_H2O2)  = WScav(   1.4,  0.5)   ! 
    WetDep(CWDEP_HCHO)  = WScav(   0.1,  0.03)  ! 
    WetDep(CWDEP_ECfn)  = WScav(   0.05,  EFF25)
    WetDep(CWDEP_SSf)   = WScav(   1.6,  EFF25)
    WetDep(CWDEP_SSc)   = WScav(   1.6,  EFFCO)
    WetDep(CWDEP_SSg)   = WScav(   1.6,  EFFGI)
    WetDep(CWDEP_PMf)   = WScav(   1.0,  EFF25) !!
    WetDep(CWDEP_PMc)   = WScav(   1.0,  EFFCO) !!
    WetDep(CWDEP_POLLw)  = WScav(   1.0,  SUBCLFAC) ! pollen
    WetDep(CWDEP_ROOH)   = WScav(  0.05,  0.015) ! assumed half of HCHO

  ! Other PM compounds treated with SO4-LIKE array defined above

   !####################### gather indices from My_Derived
   ! WDEP_WANTED array, and determine needed indices in d_2d

     do n = 1, size( WDEP_WANTED(:)%txt1 )

       if2   = find_index("WDEP_"//WDEP_WANTED(n)%txt1,f_2d(:)%name)
       atw = -999
       if( WDEP_WANTED(n)%txt3 == "mgS" ) atw = atwS
       if( WDEP_WANTED(n)%txt3 == "mgN" ) atw = atwN
       if( WDEP_WANTED(n)%txt3 == "mgSS") atw = 58
       if( WDEP_WANTED(n)%txt3 == "mgP")  atw = 12
       if( WDEP_WANTED(n)%txt3 == "mm")   atw = 999 ! Dummy for precip
       call CheckStop( atw <1 , "AQUEOUS ATW PROBLEM:" // trim(WDEP_WANTED(n)%txt3)  ) 

       if ( WDEP_WANTED(n)%txt2 == "GROUP" ) then

          igr   = find_index("WDEP_"//WDEP_WANTED(n)%txt1,chemgroups(:)%name)

          if(igr>0) then
              nwgrp = nwgrp + 1
              tmpgroup(nwgrp) = typ_i3( igr, if2, atw )  ! link to array of 
                                                         ! species integers
          end if

       else if ( WDEP_WANTED(n)%txt2 == "PREC" ) then

          WDEP_PREC= find_index("WDEP_PREC",f_2d(:)%name)
          igr = -999 ! just for printout
          isp = -999 ! just for printout
          atw = -999

       else ! SPEC

          isp   = find_index(WDEP_WANTED(n)%txt1,species(:)%name)
          if(isp>0) then
              nwspec = nwspec + 1
              tmpspec(nwspec) = typ_i3( isp, if2, atw )
          end if
       end if

       if( DEBUG .and. MasterProc )  then
            write(6,"(2a,4i5)") "WETPPP ", trim(f_2d(if2)%name), if2, igr , atw
            if(igr>0) write(*,*) "WETFGROUP ", nwgrp, chemgroups(igr)%ptr, atw
            if(isp>0) write(*,*) "WETFSPEC  ", nwspec, isp, atw
       end if
     end do
     allocate(wetspec(nwspec),stat=alloc_err)
     call CheckStop( alloc_err /= 0, "alloc error wetspec")
     allocate(wetgroup(nwgrp),stat=alloc_err)
     call CheckStop( alloc_err /= 0, "alloc error wetgroup")

    ! And now we fill these arrays with the right indices:
    ! Simplifies the code a little, but we still use the tmp arrays
    ! to store the f2d and atw info.

     wetspec(:)  = tmpspec(1:nwspec)%int1
     wetgroup(:) = tmpgroup(1:nwgrp)%int1
       
   !####################### END indices here ##########

   ! Now create table to map calc species to actual advected ones:
     Calc2tot = 0
     do n = 1, NWETDEP_ADV
         icalc = WDepMap(n)%calc
         itot  = WDepMap(n)%ind
         Calc2tot(icalc,0) =  Calc2tot(icalc,0)  + 1
         nc = Calc2tot(icalc,0)
     if( MasterProc .and.DEBUG) write(6,"(a,4i5)") "CHECKING WetDep Calc2tot ", &
            n,icalc,itot,nc
         Calc2tot(icalc,nc) = itot
      end do 

     if( MasterProc.and.DEBUG ) then
      write(*,*) "FINAL WetDep Calc2tot "
      do icalc = 1, NWETDEP_CALC
        write(*,"(i3,i4,15(1x,a))") icalc, Calc2tot(icalc,0 ), &
         ( trim( species(Calc2tot(icalc,nc))%name ), &
             nc= 1, Calc2tot(icalc,0 ))
      end do
     end if

  end subroutine Init_WetDep

!-----------------------------------------------------------------------
subroutine Setup_Clouds(i,j,debug_flag)

!-----------------------------------------------------------------------
! DESCRIPTION
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
       b           &  !  Cloud-area (fraction)
      ,cloudwater  &  !  Cloud-water (volume mixing ratio) 
                      !  cloudwater = 1.e-6 same as 1.g m^-3
      ,pres           !Pressure (Pa)
  integer :: k

! Add up the precipitation in the column:
  pr_acc(KUPPER) = sum ( pr(i,j,1:KUPPER) ) ! prec. from above 
  do k= KUPPER+1, KMAX_MID
     pr_acc(k) = pr_acc(k-1) + pr(i,j,k)
     pr_acc(k) = max( pr_acc(k), 0.0 )
  end do

  prclouds_present = .false.  
  if ( pr_acc(KMAX_MID) > PR_LIMIT )  prclouds_present = .true. 
                             ! --> precipitation at the surface

! initialise with .false. and 0:
  incloud(:)  = .false.
  cloudwater(:) = 0.
!hf
  pres(:)=0.0

! Loop starting at surface finding the cloud base:

  ksubcloud = KMAX_MID+1       ! k-coordinate of sub-cloud limit

  do k = KMAX_MID, KUPPER, -1  
     if ( lwc(i,j,k) >  CW_LIMIT ) exit
     ksubcloud = k
  end do
  if ( ksubcloud == 0 ) return ! No cloud water found below level 6
                               ! Cloud above level 6 are likely thin
                               ! cirrus clouds, and if included may
                               ! need special treatment...
                               ! ==> assume no cloud 

! Define incloud part of the column requiring that both cloud water
! and cloud fractions are above limit values 

  kcloudtop = -1               ! k-level of cloud top
  do k = KUPPER, ksubcloud-1

     b(k) = cc3d(i,j,k)

! Units: kg(w)/kg(air) * kg(air(m^3) / density of water 10^3 kg/m^3
! ==> cloudwater (volume mixing ratio of water to air in cloud 
! (when devided by cloud fraction b )
! cloudwater(k) = 1.0e-3 * cw(i,j,k,1) * roa(i,j,k,1) / b(k)   

     if ( lwc(i,j,k) >  CW_LIMIT ) then
        cloudwater(k) = lwc(i,j,k) / b(k) ! value of cloudwater in the
                                          ! cloud fraction of the grid

!hf : alternative if cloudwater exists (and can be used) from met model
!        cloudwater(k) = 1.0e-3 * cw(i,j,k,1) * roa(i,j,k,1) / b(k)
!        cloudwater min 0.03 g/m3 (0.03e-6 mix ratio)   
!        cloudwater(k) = max(0.3e-7, (1.0e-3 * cw(i,j,k,1) * roa(i,j,k,1) ))
!        cloudwater(k) =         cloudwater(k)/ b(k)
        incloud(k) = .true.
!hf
        pres(k)=ps(i,j,1)
        if ( kcloudtop < 0 ) kcloudtop = k
     end if

  end do

  if ( kcloudtop == -1 ) then
     if ( prclouds_present ) then
     if ( DEBUG ) write(6,"(a20,2i5,3es12.4)") &
        "ERROR prclouds sum_cw", &
         i,j, maxval(lwc(i,j,KUPPER:KMAX_MID),1) , &
         maxval(pr(i,j,:)), pr_acc(KMAX_MID)
     end if
     kcloudtop = KUPPER ! for safety
  end if

! sets up the aqueous phase reaction rates (SO2 oxidation) and the 
! fractional solubility 

!hf add pres
  call setup_aqurates(b ,cloudwater,incloud,pres)

  if ( DEBUG_pH .and. debug_flag .and. incloud(kcloudtop) ) then
!     write(6,"(a,l1,2i4,es14.4)") "DEBUG_AQ ",prclouds_present, &
!                kcloudtop, ksubcloud, pr_acc(KMAX_MID)

     write(6,*) "DEBUG_pH ",prclouds_present, &
                kcloudtop, ksubcloud, (pH(k),k=kcloudtop,ksubcloud-1)
     write(6,*) "CONC (mol/l)",so4_aq(ksubcloud-1),no3_aq(ksubcloud-1),nh4_aq(ksubcloud-1),nh3_aq(ksubcloud-1),hco3_aq(ksubcloud-1),co2_aq(ksubcloud-1)
write(6,*)"H+(ph_factor) ",hco3_aq(ksubcloud-1)+2.*so4_aq(ksubcloud-1)+hso3_aq(ksubcloud-1)+2.*so32_aq(ksubcloud-1)+no3_aq(ksubcloud-1)-nh4_aq(ksubcloud-1)-nh3_aq(ksubcloud-1)
     write(6,*) "CLW(l_vann/l_luft) ",cloudwater(ksubcloud-1)
     write(6,*) "xn_2d(SO4) ugS/m3 ",(xn_2d(SO4,k)*10.e12*32./AVOG,k=kcloudtop,KMAX_MID)
  end if

end subroutine Setup_Clouds


!-----------------------------------------------------------------------
subroutine init_aqueous() 

!-----------------------------------------------------------------------
! DESCRIPTION
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

! tabulations
!========================
  call tabulate_aqueous()
!========================

! Constant rates: The rates given in Berge (1993) are in mol-1 l.
! These need to be multiplied by 1.0e3/AVOG/Vf,so we perform the
! 1.0e3/AVOG scaling here.

! so2aq + h2o2   ---> so4, ref: Möller 1980
  aqrc(1) = 8.3e5 * 1.0e3/AVOG * MASSTRLIM

! (so2aq + hso3-) + H+ + o3 ---> so4, ref: Martin & Damschen 1981
  aqrc(2) = 1.8e4 * 1.0e3/AVOG * MASSTRLIM

! (so2aq + hso3-) + o2 ( + Fe ) --> so4, see documentation below
  aqrc(3) = 3.3e-10  * MASSTRLIM  

! Regarding aqrc(3):
! catalytic oxidation with Fe. The assumption is that 2% of SIV
! is oxidised per hour inside the droplets, corresponding to a
! conversion rate of 5.6^-6 (units s^-1  -- Therfore no conversion
! from mol l^-1)

! 5.6e-6 * 0.5e-6 (liquid water fraction) /8.5e-3 (fso2 at 10deg C)

! multiply with the assumed liquid water fraction from Seland and
! Iversen (0.5e-6) and with an assumed fso2 since the reaction is
! scaled by the calculated value for these parameters later.

end subroutine init_aqueous


!-----------------------------------------------------------------------
subroutine tabulate_aqueous()

!-----------------------------------------------------------------------
! DESCRIPTION
! Tabulates Henry's law coefficients over the temperature range
! defined in Tabulations_ml.
! For SO2, the effective Henry's law is given by
!   Heff = H * ( 1 + K1/H+ )  
! where k2 is omitted as it is significant only at high pH.
! We tabulate also the factor 1+K1/H+ as K1fac.
!-----------------------------------------------------------------------

  real, dimension(CHEMTMIN:CHEMTMAX) :: t, tfac ! Temperature, K, factor
  integer :: i

  t(:)           = (/ ( real(i), i=CHEMTMIN, CHEMTMAX ) /)
  tfac(:)        = 1.0/t(:) -  1.0/298.0

  H (IH_SO2 ,:)  = 1.23    * exp(3020.0*tfac(:) ) 
  H (IH_H2O2,:)  = 7.1e4   * exp(6800.0*tfac(:) )
  H (IH_O3  ,:)  = 1.13e-2 * exp(2300.0*tfac(:) )
  H (IH_NH3 ,:)  = 60.0    * exp(4400.0*tfac(:) ) !http://www.ceset.unicamp.br/~mariaacm/ST405/Lei%20de%20Henry.pdf
  H (IH_CO2,:)   = 3.5e-2  * exp(2400.0*tfac(:) ) !http://www.ceset.unicamp.br/~mariaacm/ST405/Lei%20de%20Henry.pdf
!hf NEW:

  K1(:) = 1.23e-2  * exp(2010.0*tfac(:)  )
  K2(:) = 6.6e-8   * exp(1122.0*tfac(:)  )!Seinfeldt&Pandis 1998
  Knh3(:) = 1.7e-5 * exp(-4353.0*tfac(:) )!Seinfeldt&Pandis 1998
  Kw(:) = 1.0e-14  * exp(-6718.0*tfac(:) )!Seinfeldt&Pandis 1998
  Kco2(:) = 4.3e-7 * exp(-921.0*tfac(:)  )!Seinfeldt&Pandis 1998

! Need  effective Henry's coefficient for SO2:
!hf belove moved to setup_clouds because pH is needed 
!hf moved to setup_clouds  K1fac(IH_SO2  ,:)  =  &
!hf moved      ( 1.0 + 1.23e-2 * exp(2010.0*tfac(:) ) * INV_Hplus)
!WAS A BUG HERE. Above=(1+K1)/H+ but should be K1fac=1+K1/H+
!hf moved  H (IH_SO2 ,:)  = H(IH_SO2,:) * K1fac(IH_SO2,:)

end subroutine tabulate_aqueous

!-----------------------------------------------------------------------

subroutine setup_aqurates(b ,cloudwater,incloud,pres)

!-----------------------------------------------------------------------
! DESCRIPTION
! sets the rate-coefficients for thr aqueous-phase reactions
!-----------------------------------------------------------------------


  real, dimension(KUPPER:KMAX_MID) :: &
     b            & ! cloud-aread (fraction)
    ,cloudwater   &  ! cloud-water
    ,pres           !Pressure(Pa) !hf

  logical, dimension(KUPPER:KMAX_MID) :: &
           incloud  ! True for in-cloud k values 

! Outputs  -> aqurates

! local
  real, dimension(KUPPER:KMAX_MID) :: &
         fso2grid & ! f_aq * b = {f_aq}
        ,fso2aq   & ! only so2.h2o part (not hso3- and so32-)
        ,caqh2o2  & ! rate of oxidation of so2 with H2O2
        ,caqo3    & ! rate of oxidation of so2 with H2O2
        ,caqsx      ! rate of oxidation of so2 with o2 ( Fe )

! PH
  real, dimension(KUPPER:KMAX_MID) :: &
   phfactor &
   ,h_plus 

  real, parameter :: CO2conc_ppm = 392 !mix ratio for CO2 in ppm
  real :: CO2conc !Co2 in mol/l
 !real :: invhplus04, K1_fac,K1K2_fac, Heff,Heff_NH3
  real :: invhplus04, K1K2_fac, Heff,Heff_NH3
  integer, parameter :: pH_ITER = 5 ! num iter to calc pH. Could probably be reduced     
  real, dimension (KUPPER:KMAX_MID) :: VfRT ! Vf * Rgas * Temp

  integer k, iter



  call get_frac(cloudwater,incloud) ! => frac_aq

    
! initialize:
  aqrck(:,:)=0.

! for PH
  so4_aq(:)=0.0
  no3_aq(:)=0.0
  nh4_aq(:)=0.0

  phfactor(:)=0.0
  pH(:)=0.0

! Gas phase ox. of SO2 is "default"
! in cloudy air, only the part remaining in gas phase (not
! dissolved) is oxidized 

  aqrck(ICLOHSO2,:) = 1.0 

  do k = KUPPER,KMAX_MID
     if ( incloud(k) ) then ! Vf > 1.0e-10) ! lwc > CW_limit

!For pH calculations:
!Assume total uptake of so4,no3,hno3,nh4+
!For pH below 5, all NH3 will be dissolved, at pH=6 around 50%
!Effectively all dissolved NH3 will ionize to NH4+ (Seinfeldt)

        so4_aq(k)= (xn_2d(SO4,k)*1000./AVOG)/cloudwater(k) !xn_2d=molec cm-3
                                                      !cloudwater volume mix. ratio
                                                      !so4_aq= mol/l
        no3_aq(k)= ( (xn_2d(NO3_F,k)+xn_2d(HNO3,k))*1000./AVOG)/cloudwater(k)
        nh4_aq(k)=  ( xn_2d(NH4_F,k) *1000./AVOG )/cloudwater(k)!only nh4+ now
        hso3_aq(k)= 0.0 !initial, before dissolved
        so32_aq(k)= 0.0
        nh3_aq(k) = 0.0 !nh3 dissolved and ionized to nh4+(aq)
        hco3_aq(k) = 0.0 !co2 dissolved and ionized to hco3

        VfRT(k) = cloudwater(k) * RGAS_ATML * temp(k)

           !dissolve CO2 and SO2 (pH independent)          
           !CO2conc=392 ppm
           CO2conc=CO2conc_ppm * 1e-9 * pres(k)/(RGAS_J *temp(k)) !mol/l

           frac_aq(IH_CO2,k) = 1.0 / ( 1.0+1.0/( H(IH_CO2,itemp(k))*VfRT(k) ) )
           co2_aq(k)=frac_aq(IH_CO2,k)*CO2CONC /cloudwater(k)

           frac_aq(IH_SO2,k) = 1.0 / ( 1.0+1.0/( H(IH_SO2,itemp(k))*VfRT(k) ) )
           so2_aq(k)= frac_aq(IH_SO2,k)*(xn_2d(SO2,k)*1000./AVOG)/cloudwater(k)

        do iter = 1,pH_ITER !iteratively calc pH

           phfactor(k)=hco3_aq(k)+2.*so4_aq(k)+hso3_aq(k)+2.*so32_aq(k)+no3_aq(k)-nh4_aq(k)-nh3_aq(k)
           h_plus(k)=0.5*(phfactor(k) + sqrt(phfactor(k)*phfactor(k)+4.*1.e-14) )
           h_plus(k)=min(1.e-1,max(h_plus(k),1.e-7))! between 1 and 7
           pH(k)=-log(h_plus(k))/log(10.)

           !nh4+, hco3, hso3 and so32 dissolve and ionize
           Heff_NH3= H(IH_NH3,itemp(k))*Knh3(itemp(k))*h_plus(k)/Kw(itemp(k)) 
           frac_aq(IH_NH3,k) = 1.0 / ( 1.0+1.0/( Heff_NH3*VfRT(k) ) )
           nh3_aq(k)= frac_aq(IH_NH3,k)*(xn_2d(NH3,k)*1000./AVOG)/cloudwater(k)

           hco3_aq(k)= co2_aq(k) * Kco2(itemp(k))/h_plus(k)
           hso3_aq(k)= so2_aq(k) * K1(itemp(k))/h_plus(k)
           so32_aq(k)= hso3_aq(k) * K2(itemp(k))/h_plus(k)

        enddo

!after pH determined, final numbers of frac_aq(IH_SO2) 
!= effective fraction of S(IV):
!include now also ionization to SO32-
!        K1_fac  =  &
!             1.0 + K1(k)/h_plus(k) !not used
!        H (IH_SO2 ,itemp(k))  = H(IH_SO2,itemp(k)) * K1_fac
        invhplus04= (1.0/h_plus(k))**0.4
        K1K2_fac=&
             1.0 + K1(itemp(k))/h_plus(k) + K1(itemp(k))*K2(itemp(k))/(h_plus(k)**2)
        Heff  = H(IH_SO2,itemp(k)) * K1K2_fac
        frac_aq(IH_SO2,k) = 1.0 / ( 1.0+1.0/( Heff*VfRT(k) ) )
        fso2grid(k) = b(k) * frac_aq(IH_SO2,k)!frac of S(IV) in grid in 
                                              !aqueous phase
!        fso2aq  (k) = fso2grid(k) / K1_fac 
        fso2aq  (k) = fso2grid(k) / K1K2_fac  !frac of SO2 in total grid 
                                              !in aqueous phase
                                             
        
        caqh2o2 (k) = aqrc(1) * frac_aq(IH_H2O2,k) / cloudwater(k)
        caqo3   (k) = aqrc(2) * frac_aq(IH_O3,k) / cloudwater(k)
        caqsx   (k) = aqrc(3) / cloudwater(k)  
        
        ! oh + so2 gas-phase
        aqrck(ICLOHSO2,k) = ( 1.0-fso2grid(k) ) ! now correction factor!
        aqrck(ICLRC1,k)   = caqh2o2(k) * fso2aq(k) !only SO2
        !        aqrck(ICLRC2,k)   = caqo3(k) * INV_Hplus0p4 * fso2grid(k)
        aqrck(ICLRC2,k)   = caqo3(k) * invhplus04 * fso2grid(k)
        aqrck(ICLRC3,k)   = caqsx(k) *  fso2grid(k) 
        
     end if
  enddo

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
      cloudwater   ! Volume fraction - see notes above.
  logical, dimension(KUPPER:KMAX_MID), intent(in) :: &
      incloud      ! True for in-cloud k values
  real    :: VfRT ! Vf * Rgas * Temp
  integer :: ih, k ! index over species with Henry's law, vertical level k

! Make sure frac_aq is zero outside clouds:
  frac_aq(:,:) = 0.0

  do k = KUPPER, KMAX_MID
     if ( incloud(k) ) then

        VfRT = cloudwater(k) * RGAS_ATML * temp(k)

! Get aqueous fractions:
        do ih = 1, NHENRY
           frac_aq(ih,k) = 1.0 / ( 1.0+1.0/( H(ih,itemp(k))*VfRT ) )
        end do

     end if
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
  integer :: itot,is !  index in xn_2d arrays
  integer :: k, icalc

  real    :: invgridarea ! xm2/(h*h)
  real    :: f_rho       ! Factors in rho calculation
  real    :: rho(KUPPER:KMAX_MID)
  real    :: loss    ! conc. loss due to scavenging
  real, dimension(KUPPER:KMAX_MID) :: vw ! Svavenging rates (tmp. array)
  real, dimension(KUPPER:KMAX_MID) :: lossfac ! EGU

  invgridarea = xm2(i,j)/( gridwidth_m*gridwidth_m )
  f_rho  = 1.0/(invgridarea*GRAV*ATWAIR)
! Loop starting from above:
  do k=kcloudtop, KMAX_MID           ! No need to go above cloudtop
   rho(k)  = f_rho*(dA(k) + dB(k)*ps(i,j,1))/ amk(k)
  end do

  wdeploss(:) = 0.0


! calculate concentration after wet deposition and sum up the vertical
! column of the depositions for the fully soluble species.

  if ( DEBUG .and. debug_flag ) then
     Write(6,*) "(a15,2i4,es14.4)", "DEBUG_WDEP2", &
                   kcloudtop, ksubcloud, pr_acc(KMAX_MID)
  end if ! DEBUG

  do icalc = 1, NWETDEP_CALC  ! Here we loop over "model" species

! Put both in- and sub-cloud scavenging ratios in the array vw:

!TMP xnloss = 0.0
     vw(kcloudtop:ksubcloud-1) = WetDep(icalc)%W_sca ! Scav. for incloud
     vw(ksubcloud:KMAX_MID  )  = WetDep(icalc)%W_sub ! Scav. for subcloud

     do k = kcloudtop, KMAX_MID

        lossfac(k)  = exp( -vw(k)*pr_acc(k)*dt ) 

        ! For each "calc" species we have often a number of model
        ! species
            do is = 1, Calc2tot(icalc,0)  ! number of species
               itot = Calc2tot(icalc,is)

    ! For semivolatile species only the particle fraction is deposited

               if ( itot >= FIRST_SEMIVOL .and. itot <=  LAST_SEMIVOL) then
                  loss = xn_2d(itot,k) * Fpart(itot,k)*( 1.0 - lossfac(k)  )
               else
                  loss = xn_2d(itot,k) * ( 1.0 - lossfac(k)  )
               endif
               xn_2d(itot,k) = xn_2d(itot,k) - loss
               wdeploss(itot) = wdeploss(itot) + loss * rho(k)
            end do ! is


     end do ! k loop
     if ( DEBUG .and. debug_flag .and. pr_acc(20)>1.0e-5 ) then
        do k = kcloudtop, KMAX_MID
           Write(6,"(a,2i4,a,9es12.2)") "DEBUG_WDEP, k, icalc, spec", k, &
            icalc, trim(species(itot)%name), vw(k), pr_acc(k), lossfac(k)
        end do ! k loop
     end if ! DEBUG

  end do ! icalc loop

  d_2d(WDEP_PREC,i,j,IOU_INST) = sum ( pr(i,j,:) ) * dt
                                              ! Same for all models

! add other losses into twetdep and wdep arrays:
  call WetDep_Budget(i,j,invgridarea, debug_flag)

 end subroutine WetDeposition

!-----------------------------------------------------------------------
 subroutine WetDep_Budget(i,j,invgridarea, debug_flag)
     integer,            intent(in) ::  i,j
     real, intent(in)  :: invgridarea
     logical, intent(in)  :: debug_flag

     real :: wdep
     real :: fwt
     integer :: itot, n, n2, igr, f2d

   ! Process groups of species, SOX, OXN, DUST ..... as needed
     do n = 1, nwgrp
       wdep = 0.0
       igr = wetgroup(n)
       f2d = tmpgroup(n)%int2
       fwt = tmpgroup(n)%int3 * invgridarea  ! int3=atw

       do n2 = 1, size( chemgroups(igr)%ptr )
          itot = chemgroups(igr)%ptr(n2)
          wdep = wdep + wdeploss(itot)
          if( DEBUG_MY_WETDEP .and. debug_flag ) &
               call datewrite("WET-PPPGROUP: "//species(itot)%name ,&
               itot, (/ wdeploss(itot) /) )
       end do ! n2

!Hardcoded TEMPORARY!
       if(igr==INDEX_WDEP_SOX_GROUP)totwdep(IXADV_SO4)  = totwdep(IXADV_SO4)+wdep
       if(igr==INDEX_WDEP_OXN_GROUP)totwdep(IXADV_HNO3)  = totwdep(IXADV_HNO3)+wdep
       if(igr==INDEX_WDEP_RDN_GROUP)totwdep(IXADV_NH3)  = totwdep(IXADV_NH3)+wdep
       
       d_2d(f2d,i,j,IOU_INST) = wdep * fwt
     end do ! n

   ! Process individual species, SOX, OXN, DUST ..... as needed

     do n = 1, nwspec

       itot  = wetspec(n)
       f2d   = tmpspec(n)%int2
       fwt   = tmpspec(n)%int3 * invgridarea  ! int3=atw
       d_2d(f2d,i,j,IOU_INST) = wdeploss(itot) * fwt

        if( DEBUG_MY_WETDEP .and. debug_flag ) &
             call datewrite("WET-PPPSPEC: "//species(itot)%name ,&
             itot, (/ wdeploss(itot) /) )

     end do ! n

  end subroutine WetDep_Budget


!-----------------------------------------------------------------------
end module Aqueous_ml
