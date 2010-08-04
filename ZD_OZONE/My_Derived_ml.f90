! <My_Derived_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007 met.no
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

!==============================================================================
module My_Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module specifies the "derived" fields, such as accumulated
  ! precipitation
  ! or sulphate, daily, monthly or yearly averages, depositions. These fields
  ! are all typically output as netCDF fields.
  !
  ! This module provides the user-defined setups which are used in Derived_ml.
  ! Derived fields are identified by a "class", such as "ADV" of "VOC", and
  ! the Derived_ml should perform any integrations for this.
  !
  ! Several often-used routines (e.g. for AOTs, acc. sulphate, are defined
  ! in the Derived_ml.f90, but users can define their own here, since
  ! we do not use "use only" in Derived_ml.
  !
  !   Only text strings used here to define wanted data
  !   All data field characteristics should be defined in Derived_ml, e.g.
  !   in f_2d arrays.
  !   Derived fields such as d_2d only exist in Derived_ml, so are
  !   accessed here through subroutine calls - using just the (i,j) part
  !   of the bigger d_2d arrays
  !---------------------------------------------------------------------------

use CheckStop_ml,  only: CheckStop, StopAll
use Chemfields_ml, only : xn_adv, xn_shl, cfac
use ChemSpecs_adv_ml        ! Use IXADV_ indices...
use ChemSpecs_shl_ml        ! Use IXSHL_ indices...
use ChemSpecs_tot_ml !,  only : SO2, SO4, HCHO, CH3CHO  &   !  For mol. wts.
                   !        ,NO2, pNO3_f, pNO3_c, HNO3, NH3, aNH4, PPM25, PPMCO &
                   !       ,O3, PAN, MPAN, SeaSalt_f, SeaSalt_c  !SS=SeaSalt
use ChemGroups_ml,  only :  OXN_GROUP, DDEP_OXNGROUP, DDEP_SOXGROUP, &
                            DDEP_RDNGROUP, SIA_GROUP, BVOC_GROUP
use ChemChemicals_ml, only : species               !  For mol. wts.
use ChemSpecs_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use EmisDef_ml,     only :  EMIS_NAME
use GridValues_ml, only : debug_li, debug_lj, debug_proc
use LandDefs_ml,  only : LandDefs, LandType, Check_LandCoverPresent ! e.g. "CF"
use MetFields_ml,        only : z_bnd, roa    ! 6c REM: zeta
use ModelConstants_ml, only : atwS, atwN, ATWAIR  &
                        , MasterProc  &
                        , SOURCE_RECEPTOR  &
                        , DEBUG => DEBUG_MY_DERIVED &
                        , KMAX_MID & ! =>  z dimension
                        , PPBINV  &  !   1.0e9
                        , MFAC       ! converts roa (kg/m3 to M, molec/cm3)

use OwnDataTypes_ml, only : Deriv, print_deriv_type, TXTLEN_DERIV
use Par_ml,    only: me, MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     limax, ljmax           ! => used x, y area
use SmallUtils_ml,  only : AddArray, LenArray, NOT_SET_STRING, WriteArray, &
                            find_index
use TimeDate_ml,   only : current_date
implicit none
private

 public  :: Init_My_Deriv
 public  :: My_DerivFunc

 private :: misc_xn             ! Miscelleaneous Sums and fractions of xn_adv


   !/** Depositions are stored in separate arrays for now - to keep size of
   !    derived arrays smaller and to allow possible move to a Deposition
   !    module at a later stage.
   !  Factor 1.0e6 converts from kg/m2/a to mg/m2/a

  !        We normally distinguish source-receptor (SR) stuff from model
  !        evaluation.  The SR runs should use as few as possible outputs
  !        to keep CPU and disc-requirements down. We define first then the
  !        minimum list of outputs for use in SR, then define an extra list
  !        of parameters needed in model evaluation, or even for the base-case
  !        of SR runs.


  !============ parameters for source-receptor modelling: ===================!

    integer, public, parameter :: MAX_NUM_DERIV2D = 200
    integer, public, parameter :: MAX_NUM_DERIV3D =   5
    character(len=TXTLEN_DERIV), public, save, &
         dimension(MAX_NUM_DERIV2D) :: wanted_deriv2d = NOT_SET_STRING
    character(len=TXTLEN_DERIV), public, save, &
         dimension(MAX_NUM_DERIV3D) ::  wanted_deriv3d = NOT_SET_STRING

    integer, private, save :: mynum_deriv2d
    integer, private, save :: mynum_deriv3d



!Mass-outputs of advected species, will be added to Derived
!Modify for SR
   integer, public, parameter, dimension(1) :: SRSURF_UG_S = (/ SO4 /)
   integer, public, parameter, dimension(2) ::   SURF_UG_S = (/ SO2, SO4 /)

   integer, public, parameter, dimension(4) :: SRSURF_UG_N = (/ pNO3_f, pNO3_c, aNH4, NO2 /)
   integer, public, parameter, dimension(3) ::  XSURF_UG_N = (/ NH3, HNO3, NO /)
   integer, public, parameter, dimension(7) ::   SURF_UG_N = (/ SRSURF_UG_N, XSURF_UG_N /)

   integer, public, parameter, dimension(1) ::  SURF_UG_C = (/ HCHO /)

   integer, public, parameter, dimension(3) :: SRSURF_UG = (/ PPM25, PPM25_FIRE, PPMCO /)
   integer, public, parameter, dimension(2) ::  XSURF_UG = (/ SeaSalt_f,SeaSalt_c /)
   integer, public, parameter, dimension(5) ::   SURF_UG = (/ SRSURF_UG, XSURF_UG /)

   integer, public, parameter, dimension(1) :: SRSURF_PPB = (/ O3 /)
   integer, public, parameter, dimension(1) ::  XSURF_PPB = (/ HCHO /)
   integer, public, parameter, dimension(2) ::   SURF_PPB = (/ SRSURF_PPB, XSURF_PPB /)

   integer, public, parameter :: NALL_SURF_UG = &
     size(SURF_UG_S) + size(SURF_UG_N) + size(SURF_UG_C) + size(SURF_UG)

   integer, public, parameter, dimension( NALL_SURF_UG ) :: &
      ALL_SURF_UGX = (/ SURF_UG_S, SURF_UG_N, SURF_UG_C, SURF_UG /)
   character(len=3), public, save, dimension( NALL_SURF_UG ) :: &
      ALL_SURF_UGTXT
   real, public, save, dimension( NALL_SURF_UG ) :: ALL_SURF_ATW

! Tropospheric columns
   integer, public, parameter, dimension(5) :: COLUMN_MOLEC_CM2 = (/ CO, CH4, C2H6, HCHO, NO2 /)
   character(len=3), public, parameter, dimension(1) :: COLUMN_LEVELS = &
      (/  "k20" /) ! , "k16", "k12", "k08" /)

    character(len=TXTLEN_DERIV), public, parameter, dimension(10) :: &
  D2_SR = (/ &
       "SURF_MAXO3  " &
      ,"SURF_ug_SIA " & !ds rv3_5_6 using groups
      ,"SURF_ug_PM25 " & !dsMay2010
      ,"SURF_ug_PMc  " & !dsMay2010
      ,"SURF_ug_PM10 " & !dsMay2010
      ,"SURF_ugN_OXN " & !dsMay2010
      ,"SURF_ugN_RDN " & !dsMay2010
      ,"SURF_ugN_TNO3" & !dsMay2010
      ,"SOMO35      " & !"D2_SOMO0    " &
      ,"PSURF       " &  ! Surface  pressure (for cross section):
  /)
!
!      ,ToDo "D2_PM25_H2O " &

  !============ Extra parameters for model evaluation: ===================!
    !character(len=TXTLEN_DERIV), public, parameter, dimension(19) :: &
    character(len=TXTLEN_DERIV), public, parameter, dimension(18) :: &
  D2_EXTRA = (/ &
       "WDEP_SO2          " &
      ,"WDEP_SO4          " &
      ,"WDEP_HNO3         " &
      ,"WDEP_pNO3_f       " &
      ,"WDEP_pNO3_c       " &
      ,"WDEP_NH3          " &
      ,"WDEP_aNH4         " &
      ,"SURF_ugN_NOX      " &
      ,"SURF_ppbC_VOC     " &
!      ,"D2_REDN           " &
!      ,"D2_SNOW           " &
!      ,"D2_SNratio        " &
      ,"Area_Grid_km2     " &
      ,"Area_Conif_Frac   " &
      ,"Area_Decid_Frac   " &
      ,"Area_Seminat_Frac " &
      ,"Area_Crops_Frac   " &
!      ,"Area_Water_D_Frac " &
      ,"HMIX              " &
!      ,"D2_HMIX00         " &
!      ,"D2_HMIX12         " &
      ,"SoilWater          " &
      ,"SoilWater_deep     " &
      ,"USTAR_NWP         " &
  /)


 ! Emissions
   integer, public, parameter, dimension(1) :: EMIS_OUT   = (/ C5H8 /)

  integer, public, parameter :: &   ! Groups for DDEP and WDEP
    SOX_INDEX = -1, OXN_INDEX = -2, RDN_INDEX = -3

 ! Ecosystem dep output uses receiver land-cover classes (LCs)
 ! which might include several landuse types, e.g. Conif in
!  D2_SO2_m2Conif.


  integer, private, save :: nOutDDep, nOutVg, nOutVEGO3
  integer, private, save :: nOutRG  ! RG = resistances and conductances
  integer, private, save :: nOutMET ! RG = resistances and conductances

 ! Mosaic-specific outputs, e.g. VG_CF_HNO3 or Rns_GR_NH3
  integer, public, save :: nMosaic = 0
  integer, public, parameter :: MAX_MOSAIC_OUTPUTS=150


   ! Specify some species and land-covers we want to output
   ! depositions for in netcdf files. DDEP_ECOS must match one of
   ! the DEP_RECEIVERS  from EcoSystem_ml.
   !
    integer, public, parameter :: NNDRYDEP = 3 ! 7 + 1 !JUST HNO3: size(DDEP_OXNGROUP)
   !integer, public, parameter, dimension(7+size(DDEP_OXNGROUP)) :: &
    integer, public, parameter, dimension(NNDRYDEP) :: &
      DDEP_SPECS = (/ SOX_INDEX, OXN_INDEX, RDN_INDEX /) ! , &
       !    SO2,  SO4, NH3, aNH4, HNO3 /) ! DDEP_OXNGROUP /)

    character(len=TXTLEN_DERIV), public, parameter, dimension(5) :: &
      DDEP_ECOS  = (/ "Grid   " , "Conif  ", "Seminat" &! "Water_D" &
                    , "Decid  ", "Crops  " /)

    integer, public, parameter, dimension(2) :: &
      WDEP_SPECS = (/ SO2,  SO4 /)! , aNH4, NH3, pNO3_f, HNO3, pNO3_c /)
      !WDEP_SPECS = (/ SO2,  SO4, aNH4, NH3, pNO3_f, HNO3, pNO3_c /)


  type(Deriv), public, &
     dimension( MAX_MOSAIC_OUTPUTS ), save :: MosaicOutput

  ! Have many combinations: species x ecosystems
!  type(Deriv), public, &
!     dimension( size(DDEP_SPECS)*size(DDEP_ECOS) ), save :: OutDDep

   !ECO08 - specify some species and land-covers we want to output
   ! dep. velocities for in netcdf files. Set in My_DryDep_ml.

    character(len=TXTLEN_DERIV), public, parameter, dimension(1) :: &
      VG_LABELS = (/ "VG" /)
    integer, public, parameter, dimension(1) :: &
      VG_SPECS = (/ O3 /) ! , NH3, SO2, PPM25,  PPMCO , HNO3/)
    character(len=TXTLEN_DERIV), public, parameter, dimension(1) :: &
      VG_LCS  = (/ "Grid" /) ! , "CF  ", "SNL ", "GR  " /)

!    type(Deriv), public, &
!     dimension( size(VG_LABELS)*size(VG_SPECS)*size(VG_LCS) ),  save :: OutVg

! VEGO3 outputs for AFstY and AOTX
!To avoid many unwanted combinations of land and Y values we just give
! the name here and let the code interpret it later.
! *** to use format f3.1 or i2 for the Y  or X value for fluxes/AOT! ***
!
! AFstY vs Fst - not that the accumulated AFstY is processed here.
! the instantaneous Fst is set as for Canopy O3 in METCONCS
          ! N.B. AOTs have several definitions. We usually want
          ! the ICP-veg Mapping Manual (MM) ones. Other
          ! possibilities are EU (8-20daytime) or UN (May-July for
          ! crops)

    character(len=TXTLEN_DERIV), public, parameter, dimension(15) :: &
     VEGO3_OUTPUTS = (/ "AFST_1.6_IAM_DF", &
                        "AFST_0.0_IAM_DF", &
                        "AFST_0.0_BF    ", &
                        "AFST_1.6_BF    ", &
                        "AFST_0.0_IAM_CR", &
                        "AFST_3.0_IAM_CR", &
!                        "AFST_6.0_IAM_CR", &
                        "AFST_0.0_IAM_MF", &
                        "AFST_1.0_IAM_MF", & ! WGSR POD1 birch
                        "AFST_1.0_IAM_DF", & ! WGSR POD1
                        "AFST_1.6_IAM_MF", &
!                        "MMAOT_30_IAM_DF", &
                        "MMAOT_40_IAM_DF", & ! WGSR beech
                        "MMAOT_40_IAM_MF", & ! WGSR birch
!                        "MMAOT_30_IAM_CR", & ! only iam allowed
                        "MMAOT_40_IAM_CR", &
                        "EUAOT_40_IAM_CR", &
                        "MMAOT_40_IAM_WH" /) !NB -last not found. Could
                                             !just be skipped, but kept
                                             !to show behaviour
!    type(Deriv), public, &
!     dimension( size(VEGO3_OUTPUTS) ),  save :: OutVEGO3

! For resistances and conductances we normally want the same landuse
! outputs, so we use a combined variable:

    character(len=TXTLEN_DERIV), public, parameter, dimension(3) :: &
      RG_LABELS = (/ "Rs ", "Rns", "Gns" /)
    integer, public, parameter, dimension(2) :: &
      RG_SPECS = (/ NH3, SO2/)
    character(len=TXTLEN_DERIV), public, parameter, dimension(1) :: &
      RG_LCS  = (/ "Grid" /) !, "CF  ", "SNL ", "GR  " /)
!    character(len=TXTLEN_DERIV), public, parameter, dimension(6) :: &
!      RG_LCS  = (/ "Grid", "CF", "SNL", "GR" , "TESTRG","TC"/)

!    type(Deriv), public, & !Non-stomatal conductance
!     dimension( size(RG_LABELS)*size( RG_SPECS)*size( RG_LCS) ),  save :: OutRG

! For met-data and canopy concs/fluxes ...

    character(len=TXTLEN_DERIV), public, parameter, dimension(2) :: &
      METCONC_PARAMS = (/ "RH      " &
                         ,"CanopyO3" & !SKIP ,"VPD     ", "FstO3   " &
                        !SKIP  ,"USTAR   ", "INVL    "  &
                       /)
                          ! "g_sto" needs more work - only set as L%g_sto
    integer, public, save :: MMC_USTAR, MMC_INVL, MMC_RH, MMC_CANO3,  &
           MMC_VPD, MMC_FST
    character(len=TXTLEN_DERIV), public, save, dimension(2) :: &
      MET_LCS  = (/ "GR    " , "CF    " /) !, "IAM_DF", "IAM_MF"/)
      !MET_LCS  = (/ "GR    " , "IAM_CR", "IAM_DF", "IAM_MF"/)
    !character(len=TXTLEN_DERIV), public, parameter, dimension(5) :: &
      !MET_LCS  = (/ "CF", "SNL", "TESTCF", "GR" ,"TC"/)
      ! Can also set dim 4:1 to exclude all - gives zero size MET_LCS

   ! We use Dep_type anyway since we can use the LCC elements
!    type(Dep_type), public, & !Non-stomatal conductance
!     dimension( size(METCONC_PARAMS)*size( MET_LCS) ),  save :: OutMET

!----------------------
! Less often needed:
!exc  "D2T_HCHO  ","D2T_CH3CHO","D2_O3TC   ","D2_ACCSU  ",
!"D2_FRNIT  ","D2_MAXOH  "

   !======= MY_DERIVED SYSTEM ======================================

  ! use character arrays to specify which outputs are wanted

   character(len=TXTLEN_DERIV), public, parameter, dimension(4) :: &
       WDEP_WANTED = (/ "WDEP_PREC", "WDEP_SOX ", "WDEP_OXN ", &
                      "WDEP_RDN " /)   ! WDEP_PM not used


    ! For some reason having this as a parameter caused problems for
    ! PC-gfortran runs.

 ! hb new 3D output when ppb
    !integer, public, parameter, dimension(3) ::   D3_PPB = (/ O3, pNO3_f, pNO3_c /)
    integer, public, save, dimension(4:1) ::   D3_PPB ! = (/ O3 /)

! hb other (non-ppb) 3D output
     character(len=TXTLEN_DERIV), public, save, dimension(4:1) :: &
!AMVB 2010-07-21: PM-PPB bug fix
       D3_OTHER != (/ "D3_ug_PM25", "D3_ug_PMc",  "D3_m_TH", "D3_m2s_Kz" /)

    integer, private :: i,j,k,n, ivoc, index    ! Local loop variables

   contains

 !=========================================================================
  subroutine Init_My_Deriv()

    integer :: i, ilab, itot, nDD, nVg, nRG, nMET, nVEGO3, iLC, &
      iadv, ispec, atw, n1, n2
    integer :: isep1, isep2  ! location of seperators in sting
    character(len=TXTLEN_DERIV) :: name ! e.g. DDEP_SO2_m2Conif
    character(len=TXTLEN_DERIV) :: txt, txt2, units, txtnum
    real    :: Y           ! threshold for AFSTY, also used as X in AOT40X
    integer :: Threshold   ! threshold for AFSTY
    character(len=TXTLEN_DERIV), &
    dimension(size(COLUMN_MOLEC_CM2)*size(COLUMN_LEVELS)) :: tmpname ! e.g. DDEP_SO2_m2Conif
    character(len=100) :: errmsg
! hb new 3D output
!    character(len=TXTLEN_DERIV), dimension(NALL_SURF_UG+size(SURF_PPB)) ::&
    character(len=TXTLEN_DERIV), dimension(NALL_SURF_UG+size(SURF_PPB)+size(D3_PPB)) ::&
          tag_name    ! Needed to concatanate some text in AddArray calls
                      ! - older (gcc 4.1?) gfortran's had bug
    logical, parameter :: T=.true., F=.false.
    logical :: dt_scale
    real :: scale


    ! Set indices of mosaic metconc params for later use. Will be zero if 
    ! not found, but that's okay I hope...
      MMC_RH    = find_index("RH",METCONC_PARAMS) 
      MMC_CANO3 = find_index("CanopyO3",METCONC_PARAMS) 
      MMC_VPD   = find_index("VPD",METCONC_PARAMS) 
      MMC_FST   = find_index("FstO3",METCONC_PARAMS) 
      MMC_USTAR = find_index("USTAR",METCONC_PARAMS)
      MMC_INVL  = find_index("INVL",METCONC_PARAMS) 

   ! Build up the array wanted_deriv2d with the required field names

     call AddArray(WDEP_WANTED, wanted_deriv2d, NOT_SET_STRING,errmsg)
     call CheckStop( errmsg, errmsg // "WDEP_WANTED too long" )
     call AddArray( D2_SR,  wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "D2_SR too long" )

  ! Emission sums - we always add these (good policy!)
   do  i = 1, size(EMIS_NAME)
     tag_name(1) = "Emis_mgm2_" // trim(EMIS_NAME(i))
     call AddArray( tag_name(1:1), wanted_deriv2d, NOT_SET_STRING, errmsg)
   end do
   do  i = 1, size(BVOC_GROUP)
     itot = BVOC_GROUP(i)
     tag_name(1) = "Emis_mgm2_" // trim(species(itot)%name)
     call AddArray( tag_name(1:1), wanted_deriv2d, NOT_SET_STRING, errmsg)
   end do

! add surf concs in various units (e.g. ugS/m3 or ppb):
     tag_name(1:size(SURF_UG_S)) = "SURF_ugS_" // species(SURF_UG_S)%name
     call AddArray(  tag_name(1:size(SURF_UG_S)) , wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "SURF_ugS too long" )
     tag_name(1:size(SURF_UG_N)) = "SURF_ugN_"//species(SURF_UG_N)%name
     call AddArray( tag_name(1:size(SURF_UG_N)) ,  wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "SURF_ugN too long" )
     tag_name(1:size(SURF_UG_C)) = "SURF_ugC_"//species(SURF_UG_C)%name
     call AddArray( tag_name(1:size(SURF_UG_C)),  wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "SURF_ugC too long" )
     tag_name(1:size(SURF_UG)) = "SURF_ug_"//species(SURF_UG)%name
     call AddArray( tag_name(1:size(SURF_UG)) ,  wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "SURF_ug too long" )
     tag_name(1:size(SURF_PPB)) = "SURF_ppb_"//species(SURF_PPB)%name
     call AddArray( tag_name(1:size(SURF_PPB)) ,  wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "SURF_ppb too long" )

     if ( .not. SOURCE_RECEPTOR ) then !may want extra?
        call AddArray( D2_EXTRA, wanted_deriv2d, NOT_SET_STRING, errmsg)
        call CheckStop( errmsg, errmsg // "D2_EXTRA too long" )
     end if

     ! Column data:
     n = 0
     do n1 = 1, size(COLUMN_MOLEC_CM2)
     do n2 = 1, size(COLUMN_LEVELS)
       n = n + 1
       tmpname(n) = "COLUMN_" // trim( species(COLUMN_MOLEC_CM2(n1))%name ) // "_" // COLUMN_LEVELS(n2)
       !if(MasterProc) print *, "COL MY_DERIV", trim( species(COLUMN_MOLEC_CM2(n))%name )
     end do
     end do
     call AddArray(tmpname, wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "COLUMN too long" )
     ! Didn't work:
     !call AddArray( "COLUMN_" // trim( species(COLUMN_MOLEC_CM2(:))%name ), &
     !  wanted_deriv2d, NOT_SET_STRING)


     !------------- Dry Depositions for d_2d -------------------------
     ! Add species and ecosystem depositions if wanted:
     ! We find the various combinations of gas-species and ecosystem,
     ! adding them to the derived-type array OutDDep (e.g. => D2_SO4_m2Conif)

      nDD = 0
      do i = 1, size(DDEP_SPECS)
        do n = 1, size(DDEP_ECOS)

          nDD = nDD + 1
          ispec  = DDEP_SPECS(i)  ! Index in ix_tot arrays

          if ( ispec > 0 ) then ! normal species, e.g. DDEP_SO2_m2CF
             name = "DDEP_"  // trim( species(ispec)%name ) // &
                    "_m2" // trim( DDEP_ECOS(n) )

             iadv  = ispec - NSPEC_SHL ! adv for use in Derived

             if ( species(ispec)%sulphurs > 0 ) then
               atw  = species(ispec)%sulphurs * atwS
               units  =  "mgS/m2"

                call CheckStop( species( ispec )%nitrogens > 0 , &
                 " ERROR in DDEP_SPECS: BOTH S AND N!"// &
                   species(DDEP_SPECS(i))%name)

             else if ( species( ispec )%nitrogens > 0 ) then
               atw  = species( ispec )%nitrogens * atwN
               units  =  "mgN/m2"

             else
               call StopAll("ERROR: OutDDep atw failure "// &
                   species( ispec )%name)
             end if
          else ! GROUP
             iadv   = ispec  ! e.g. -1 for SOX
             atw   = atwN   ! first guess
             units = "mgN/m2"
            if ( ispec == SOX_INDEX ) then
                atw   = atwS
                units = "mgS/m2"
                name = "DDEP_SOX_m2"//trim(DDEP_ECOS(n))
             else if ( ispec == OXN_INDEX ) then
                name = "DDEP_OXN_m2"//trim(DDEP_ECOS(n))
             else if ( ispec == RDN_INDEX ) then
                name = "DDEP_RDN_m2"//trim(DDEP_ECOS(n))
             end if
          end if

        ! Set defaults
        ! dep_type( name, LC, index, f2d, class, label, txt, scale, atw, units )
        !            x     d      d    d   a10    a10   a10     f    i    a10
!             OutDDep(nDD) = Dep_type(  &
!              name, -99, iadv, -99,"Mosaic", "DDEP", DDEP_ECOS(n), 1.0, atw, units)

             nMosaic = nMosaic + 1
             call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, DDEP" )
          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC, XYLC, dt_scale, scale, avg? rho Inst Yr Mn Day atw
             MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", "DDEP", DDEP_ECOS(n), units, &
                  iadv,-99,-99, -99.9, F, 1.0e6,  F,   F,   F, T,  T,  F, atw)
!QUERY - why no dt_scale??

          if(DEBUG .and. MasterProc) then
            write(6,*) "DDEP setups"
            call print_deriv_type(MosaicOutput(nMosaic))
          end if
        end do ! DDEP_SPECS
     end do ! DDEP_ECOS
     nOutDDep = nDD

      !------------- Deposition velocities for d_2d -------------------------
      ! Add species and ecosystem depositions if wanted:
      ! We find the various combinations of gas-species and ecosystem,
      ! adding them to the derived-type array LCC_DDep (e.g. => VGO3_CF)

      nVg = 0
      do ilab = 1, size(VG_LABELS)
      do i = 1, size(VG_SPECS)
        VG_LC: do n = 1, size(VG_LCS)

          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "VG_LCS", n, VG_LCS, &
                   write_condition=(i==1 .and. ilab == 1))
          if ( iLC < 0 ) cycle  VG_LC
          !-------------End of Check if LC present in this array ------!
          nVg = nVg + 1

          name = trim( VG_LABELS(ilab) ) // "_"  // &
             trim( species(VG_SPECS(i))%name ) // "_" // trim( VG_LCS(n) )
          iadv = VG_SPECS(i) - NSPEC_SHL
          call CheckStop( iadv < 1 .or. iadv > NSPEC_ADV, &
                 " ERR: DDEP_SPECS: VG_SPECS" )

          nMosaic = nMosaic + 1
          call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, nVg" )

          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC,XYLC, scale, avg? rho Inst Yr Mn Day atw
          MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", VG_LABELS(ilab), VG_LCS(n), "cm/s", &
                iadv, -99,iLC, -99.9,  F, 100.0, T,   F,   F, T, T, T, atw)

        end do VG_LC !n
      end do ! i
      end do ! ilab
      nOutVg = nVg

      !------------- VEGO3 stuff ----------------------------------------------
      ! For fluxes or AOTs we start with a formatted name, eg. AFST_3.0_CF and
      !untangle it to get threshold Y (=3.0) and landcover type

      nVEGO3 = 0
      dt_scale = .false. ! for AFstY and AOT we need to reset this.
      VEGO3_LC: do n = 1, size(VEGO3_OUTPUTS)

         name = VEGO3_OUTPUTS(n)
         isep1 = scan(name,"_")                       ! AFST or AOT
         isep2 = isep1 + scan(name(isep1+1:),"_")   ! 1.6 or 40
         !isep3 = isep2 + scan(name(isep2+1:),"_")   ! IAM_CF or SNL
         txt =name(1:isep1-1)           ! AFST or AOT
         txtnum=name(isep1+1:isep2-1)     ! 1.6 or 40
         txt2=name(isep2+1:)            ! IAM_CF or SNL


         if( txt == "AFST" ) then
            read(txtnum,fmt="(f3.1)") Y
            if(txtnum == ".0") txtnum = txtnum(1:1)  ! 3.0 -> 3
            !WHY??!!!! Threshold = nint( 10*Y)   ! Store Y=1.6 as 16
            Threshold = Y
            units = "mmole/m2"
            scale = 1.0e-6   ! Accumulates nmole/s to mmole (*dt_advec)
            dt_scale = .true.   ! Accumulates nmole/s to mmole (*dt_advec)
         else if( txt(3:5) == "AOT" )  then
            read(txtnum,fmt="(i2)") Threshold ! e.g. 40 ppb
            units = "ppb.h"
            scale = 1.0/3600.0 ! AOT in ppb.hour
            dt_scale = .true.
         end if


          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "VEGO3_LCS", txt2, .true. )
          if(DEBUG .and. MasterProc)  write(*,*) "VEGO3 ", trim(name), &
               "=> Y", trim(txtnum), " iLC, LC ", iLC, trim(txt2)
          if ( iLC < 0 ) cycle  VEGO3_LC
          !DS if ( iLC > 0 .and. .not. LandType(iLC)%is_iam ) cycle  VEGO3_LC
          if ( iLC > 0  ) LandType(iLC)%flux_wanted  = .true. 
          !-------------End of Check if LC present in this array ------!
          nVEGO3 = nVEGO3 + 1

          write(unit=name,fmt="(a)") trim(txt)// trim(txtnum)//"_"//trim(txt2)

           nMosaic = nMosaic + 1
           call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, VEGO3" )
          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC,XYCL, scale dt_scale avg? rho Inst Yr Mn Day atw
           MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", txt, txt2, units, &
                iadv, -99,iLC,Threshold,  T,  scale,  F,   F,   F, T, T, F, atw)

      end do VEGO3_LC !n
      nOutVEGO3 = nVEGO3
          if(DEBUG .and. MasterProc)  print *, "VEGO3 FINAL NUM ", nVEGO3


      !------------- Surface resistance for d_2d -------------------------
      ! We find the various combinations of gas-species and ecosystem,
      ! adding them to the derived-type array LCC_DDep (e.g. => RsO3_CF)

      nRG = 0
      do ilab = 1, size(RG_LABELS)
      do i = 1, size(RG_SPECS)
        RG_LC: do n = 1, size(RG_LCS)

          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "RG_LCS", n, RG_LCS, (i==1 .and. ilab == 1))
          if ( iLC < 0 ) cycle  RG_LC
          !-------------End of Check if LC present in this array ------!

          nRG = nRG + 1
          name = trim ( RG_LABELS(ilab) ) // "_"  // &
           trim( species(RG_SPECS(i))%name ) // "_" // trim( RG_LCS(n) )

          iadv  =   RG_SPECS(i) - NSPEC_SHL

          !OutRG(nRG)%label  = RG_LABELS(ilab)
          !OutRG(nRG)%txt  =   RG_LCS(n)

!          OutRG(nRG) = Dep_type(  &
!             name, iLC, iadv, -99,"Mosaic", RG_LABELS(ilab), RG_LCS(n),&
!                                             1.0, -99,  "-")
           nMosaic = nMosaic + 1
           call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, nRG" )

          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC,XYCL, scale, avg? rho Inst Yr Mn Day atw

          if( RG_LABELS(ilab)(1:1) == "R" )  then
             MosaicOutput(nMosaic) = Deriv( &
               name, "Mosaic", RG_LABELS(ilab), RG_LCS(n), "s/m", &
                iadv, -99,iLC,-99.9, F , 1.0,  T,   F,   F, T, T, T, atw)
          else if( RG_LABELS(ilab)(1:1) == "G" )  then
             MosaicOutput(nMosaic) = Deriv( &
               name, "Mosaic", RG_LABELS(ilab), RG_LCS(n), "cm/s", &
                iadv, -99,iLC,-99.9, F, 100.0,  T,   F,   F, T, T, T, atw)
          end if
!          if( OutRG(nRG)%label(1:1) == "R" )  then
!              OutRG(nRG)%units  =   "s/m"
!          else if( OutRG(nRG)%label(1:1) == "G" )  then
!              OutRG(nRG)%scale  =    100.0
!              OutRG(nRG)%units  =   "cm/s"
!          end if
!          if(DEBUG .and. MasterProc) call print_dep_type(OutRG(nRG))
        end do RG_LC !n
      end do ! i
      end do ! ilab
      nOutRG = nRG

!     call AddArray( OutRG(:)%name, wanted_deriv2d, NOT_SET_STRING, errmsg)
!     call CheckStop( errmsg, errmsg // "OutRG  too long" )

      !------------- Met data for d_2d -------------------------
      ! We find the various combinations of met and ecosystem,
      ! adding them to the derived-type array LCC_Met (e.g. => Met_CF)

      nMET = 0
      do ilab = 1, size(METCONC_PARAMS)
        MET_LC: do n = 1, size(MET_LCS)

          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "MET_LCS", n, MET_LCS, (ilab == 1))
          if ( iLC < 0 ) cycle  MET_LC
          if ( METCONC_PARAMS(ilab)(1:6) == "Canopy" .or. & 
               METCONC_PARAMS(ilab)(1:5) == "FstO3" ) &
                 LandType(iLC)%flux_wanted  = .true.  ! Canopy calc in StoFlux
          !-------------End of Check if LC present in this array ------!

          nMET = nMET + 1
          name = trim ( METCONC_PARAMS(ilab) ) // "_"  // trim( MET_LCS(n) )

          nMosaic = nMosaic + 1
          call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, nMET" )
          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC,XYCL, scale, avg? rho Inst Yr Mn Day atw
           MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", "METCONC", MET_LCS(n), METCONC_PARAMS(ilab), &
                ilab, -99,iLC,-99.9,  F , 1.0,  T,   F,   F, T, T, T, atw)
!MESS          OutMET(nMET) = Dep_type( &
!           name, ilab, -99, iLC, "Mosaic", METCONC_PARAMS(ilab),  &
!                                              MET_LCS(n), 1.0, -99, "-")

          if( METCONC_PARAMS(ilab)(1:5) == "USTAR" )  then
              MosaicOutput(nMosaic)%unit  =   "m/s"
          else if( METCONC_PARAMS(ilab)(1:4) == "INVL" )  then
              MosaicOutput(nMosaic)%unit  =   "m"
          else if( METCONC_PARAMS(ilab)(1:8) == "CanopyO3" )  then
              MosaicOutput(nMosaic)%unit  =   "ppb"
          else if( METCONC_PARAMS(ilab)(1:5) == "FstO3" )  then
              MosaicOutput(nMosaic)%unit  =   "mmole/m2" ! accumulated
          end if

          if(DEBUG .and. MasterProc) call print_deriv_type(MosaicOutput(nMosaic))
        end do MET_LC !n
      end do ! ilab
      nOutMET = nMET

!     call AddArray( OutMET(:)%name, wanted_deriv2d, NOT_SET_STRING, errmsg)
!     call CheckStop( errmsg, errmsg // "OutMET too long" )

      !------------- end LCC data for d_2d -------------------------

     call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics" )
     call AddArray( MosaicOutput(1:nMosaic)%name, &
                        wanted_deriv2d, NOT_SET_STRING, errmsg)
     call CheckStop( errmsg, errmsg // "MosaicOutput too long" )

     mynum_deriv2d  = LenArray( wanted_deriv2d, NOT_SET_STRING )

   ! ditto wanted_deriv3d....
! hb new 3D ouput
     if ( .not. SOURCE_RECEPTOR ) then
       if ( size(D3_PPB) > 0 ) then
         tag_name(1:size(D3_PPB)) = "D3_ppb_" // species(D3_PPB)%name
         call AddArray(  tag_name(1:size(D3_PPB)) , wanted_deriv3d, &
            NOT_SET_STRING, errmsg)
         call CheckStop( errmsg, errmsg // "D3_ppb too long" )
       end if
       if ( size(D3_OTHER) > 0 ) then
         call AddArray( D3_OTHER,  wanted_deriv3d, NOT_SET_STRING, errmsg)
         call CheckStop( errmsg, errmsg // "Wanted D3 too long" )
       end if
     end if
     mynum_deriv3d  = LenArray( wanted_deriv3d, NOT_SET_STRING )


     if ( MasterProc ) then
       write(*,*) "Init_My_Deriv, mynum_deriv2d = ", mynum_deriv2d
       write(*,*) "Init_My_Deriv, mynum_deriv3d = ", mynum_deriv3d
       if(  DEBUG ) then
         do i = 1, mynum_deriv2d
           write(*,*) "DEBUG DERIV2D ", i, mynum_deriv2d, wanted_deriv2d(i)
         end do
       end if
       call WriteArray(wanted_deriv2d,mynum_deriv2d," Wanted 2d array is")
       call WriteArray(wanted_deriv3d,mynum_deriv3d," Wanted 3d array is")

     end if

  end subroutine Init_My_Deriv
 !=========================================================================
  subroutine My_DerivFunc( e_2d, class , density )

    ! We define here here any functions which cannot easily be defined
    ! in the more general Derived_ml. For example, we need the
    ! index for IXADV_O3 for AOTs, and this might not be available in the model
    ! we are running (a PM2.5 model for example), so it is better to define
    ! this function here.

  real, dimension(:,:), intent(inout) :: e_2d  !  (i,j) 2-d extract of d_2d
  character(len=*), intent(in)    :: class       ! Class of data

  real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density
! density = 1 ( or = roa when unit ug)

  select case ( class )

      case ( "OX", "NOX", "NOZ", "TOXN", "TRDN", "FRNIT", "tNO3 ", "SSalt" )

           call misc_xn( e_2d, class, density )

      case  default

            write(*,*) "WARNING - REQUEST FOR UNDEFINED OUTPUT:", n, class
     end select



  end subroutine My_DerivFunc
 !=========================================================================

  subroutine misc_xn( e_2d, class, density)
    real, dimension(:,:), intent(inout) :: e_2d  ! i,j section of d_2d arrays
    character(len=*)    :: class   ! Type of data
    real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density
    integer :: itot, iadv, igrp
! density = 1 ( or = roa when unit ug)


    !/--  adds up sulphate, nitrate, or whatever is defined

    select case ( class )

    case ( "TOXN" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_HNO3,i,j,KMAX_MID) * cfac(IXADV_HNO3,i,j) &
              + xn_adv(IXADV_pNO3_f,i,j,KMAX_MID) * cfac(IXADV_pNO3_f,i,j) &
              + xn_adv(IXADV_pNO3_c,i,j,KMAX_MID) * cfac(IXADV_pNO3_c,i,j)) &
              * density(i,j)
      end forall


! OX for O3 and NO2 trend studies

    case ( "OX" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
                xn_adv(IXADV_O3,i,j,KMAX_MID)  * cfac(IXADV_O3,i,j)   &
              + xn_adv(IXADV_NO2,i,j,KMAX_MID) * cfac(IXADV_NO2,i,j)
      end forall

    case ( "NOX" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_NO,i,j,KMAX_MID) &
              + xn_adv(IXADV_NO2,i,j,KMAX_MID) * cfac(IXADV_NO2,i,j) &
              ) * density(i,j)
      end forall

    case ( "NOZ" )
      e_2d( :,: ) = 0.0
      do i=1,limax
        do j=1,ljmax
          do igrp = 1, size(OXN_GROUP)
            itot = OXN_GROUP(igrp)
            iadv = OXN_GROUP(igrp) - NSPEC_SHL
            e_2d( i,j ) = e_2d( i,j ) + &
              xn_adv(iadv,i,j,KMAX_MID) * &
                cfac(iadv,i,j) * species(itot)%nitrogens
          end do ! n
          e_2d( i,j ) = e_2d( i,j ) * density(i,j)
        end do ! j
      end do ! i

    case ( "TRDN" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
               ( xn_adv(IXADV_NH3,i,j,KMAX_MID) * cfac(IXADV_NH3,i,j)    &
              +  xn_adv(IXADV_aNH4,i,j,KMAX_MID) * cfac(IXADV_aNH4,i,j))  &
               * density(i,j)
      end forall


    case ( "FRNIT" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
             ( xn_adv(IXADV_pNO3_f,i,j,KMAX_MID) * cfac(IXADV_pNO3_f,i,j)  &
            +  xn_adv(IXADV_pNO3_c,i,j,KMAX_MID) * cfac(IXADV_pNO3_c,i,j)) &
       /max(1E-80, (xn_adv(IXADV_HNO3,i,j,KMAX_MID) *  cfac(IXADV_HNO3,i,j))&
            +  xn_adv(IXADV_pNO3_f,i,j,KMAX_MID) * cfac(IXADV_pNO3_f,i,j)    &
            +  xn_adv(IXADV_pNO3_c,i,j,KMAX_MID) * cfac(IXADV_pNO3_c,i,j))
      end forall

    case ( "tNO3" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d(  i,j ) = &
              ( xn_adv(IXADV_pNO3_f,i,j,KMAX_MID) * cfac(IXADV_pNO3_f,i,j) &
              + xn_adv(IXADV_pNO3_c,i,j,KMAX_MID) * cfac(IXADV_pNO3_c,i,j) )&
              * density(i,j)
      end forall

    case ( "SSalt" )
      forall ( i=1:limax, j=1:ljmax )
          e_2d( i,j ) = &
              ( xn_adv(IXADV_SeaSalt_f,i,j,KMAX_MID) * cfac(IXADV_SeaSalt_f,i,j) &
              + xn_adv(IXADV_SeaSalt_c,i,j,KMAX_MID) * cfac(IXADV_SeaSalt_c,i,j) )&
              * density(i,j)
      end forall

    end select
  end subroutine misc_xn
 !=========================================================================
end module My_Derived_ml
