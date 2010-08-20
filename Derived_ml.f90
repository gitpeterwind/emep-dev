! <Derived_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module performs the calculations associated with "derived" 2D and 3D,
  ! such as accumulated precipitation or sulphate, daily, monthly or yearly
  ! averages, depositions. These fields are all typically output as netCDF
  ! fields.
  !
  ! This routine defines many possible derived  outputs.
  ! The names of the derived fields actualy required should have been specified
  !  in the user-defined My_Derived_ml.
  !

  ! User-defined routines and treatments are often needed here. Here there is
  ! added stuff for VOC, AOTs, accsu. In
  ! general such code should be added in such a way that it isn't activated if
  ! not needed. It then doesn't need to be commented out if not used.
  !---------------------------------------------------------------------------

!current definitions:
!SOMO35: daily max is found over 00:00 to 24:00. (not emepday)
!SOMO35: accumulated over one year
!D2_MAXO3 :  daily max is found over an EMEPDAY
!D2_MAXO3 : accumulated into yearly_output from April to September
!D2_O3 is now yearly accumulated

use My_Derived_ml, only : &
            wanted_deriv2d, wanted_deriv3d  &! names of wanted derived fields
           ,Init_My_Deriv, My_DerivFunc
use My_Derived_ml,  only : & !EcoDep
      nMosaic, &
      MosaicOutput, &
      COLUMN_MOLEC_CM2, &  !ds added May 2009
      COLUMN_LEVELS   , &  !ds added Dec 2009
      SURF_UG_S, &  !ds added May 2009
      SURF_UG_N, &  !ds added May 2009
      SURF_UG_C, &  !ds added May 2009
      SURF_UG  , &  !ds added May 2009
      SURF_PPB , &  !ds added May 2009
!AMVB 2010-07-21: PM-PPB bug fix
      D3_PPB,D3_OTHER ! hb new 3D output

use AOTx_ml,           only: Calc_GridAOTx, setaccumulate_2dyear
use Biogenics_ml,       only: EmisNat !dsbvoc
use CheckStop_ml,      only: CheckStop, StopAll
use Chemfields_ml, only : xn_adv, xn_shl, cfac,xn_bgn, PM_water
use ChemGroups_ml     ! SIA_GROUP, PMCO_GROUP -- use tot indices
use ChemSpecs_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use ChemSpecs_shl_ml
use ChemSpecs_tot_ml
use ChemChemicals_ml, only : species
use Chemfields_ml , only : so2nh3_24hr,Grid_snow
use EcoSystem_ml,   only : DepEcoSystem, NDEF_ECOSYSTEMS, &
                           EcoSystemFrac,FULL_GRID
use EmisDef_ml,    only : EMIS_NAME
use Emissions_ml,  only : SumSnapEmis
use GridValues_ml, only : debug_li, debug_lj, debug_proc, xm2, GRIDWIDTH_M&
                         ,GridArea_m2 ! dsbvoc
use MetFields_ml, only :   roa,pzpbl,Kz_m2s,th,zen, ustar_nwp, z_bnd
use MetFields_ml, only :   ps
use MetFields_ml, only :   SoilWater, SoilWater_deep
use ModelConstants_ml, only: &
   KMAX_MID     & ! =>  z dimension
  ,NPROC        & ! No. processors
  ,atwS, atwN, ATWAIR, dt_advec  &
  ,PPBINV       & ! 1.0e9, for conversion of units
  ,PPTINV       & ! 1.0e12, for conversion of units
  ,MFAC         & ! converts roa (kg/m3 to M, molec/cm3)
  ,DEBUG_i, DEBUG_j &
  ,DEBUG => DEBUG_DERIVED, MasterProc &
  ,SOURCE_RECEPTOR &
  ,FORECAST     & ! only dayly (and hourly) output on FORECAST mode
  ,NTDAY        & ! Number of 2D O3 to be saved each day (for SOMO)
  ! output types corresponding to instantaneous,year,month,day
  ,IOU_INST, IOU_YEAR, IOU_MON, IOU_DAY, IOU_HOUR, IOU_HOUR_MEAN

use OwnDataTypes_ml, only: Deriv,TXTLEN_DERIV   ! type & length of names
use Par_ml,    only: MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     me,                &   ! for print outs
                     gi0,gj0,IRUNBEG,JRUNBEG,&! for i_fdom, j_fdom
                     li0,lj0,limax, ljmax    ! => used x, y area
use PhysicalConstants_ml,  only : PI
use SmallUtils_ml, only: find_index, LenArray, NOT_SET_STRING
use TimeDate_ml, only : day_of_year,daynumber,current_date
implicit none
private

 public  :: Init_Derived         !
 public  :: ResetDerived         ! Resets values to zero
 public  :: DerivedProds         ! Calculates any production terms
 public  :: AddDeriv             ! Adds Deriv type to def_2d, def_3d
 public  :: AddNewDeriv          ! Creates & Adds Deriv type to def_2d, def_3d
 private :: Define_Derived       !
 private :: Setups
 private :: write_debugadv
 private :: write_debug

 public :: Derived              ! Calculations of sums, avgs etc.
 private :: voc_2dcalc          ! Calculates sum of VOC for 2d fields
 private :: voc_3dcalc          ! Calculates sum of VOC for 3d fields
 private :: uggroup_calc        ! Calculates sum of groups, e.g. pm25


   INCLUDE 'mpif.h'
   INTEGER STATUS(MPI_STATUS_SIZE),INFO
   logical, private, parameter :: T = .true., F = .false. ! shorthands only

  ! Tip. For unix users, do "grep AddDef | grep -v Is3D | wc" or similar
  ! to help get the number of these:
   integer, private, parameter ::  &
       MAXDEF_DERIV2D =250 & ! Max. No. 2D derived fields to be defined
      ,MAXDEF_DERIV3D = 17   ! Max. No. 3D derived fields to be defined

   integer, public, save :: num_deriv2d, num_deriv3d
   integer, private, save :: Nadded2d = 0, Nadded3d=0 ! No. defined derived

  ! We put definitions of **all** possible variables in def_2d, def_3d
  ! and copy the needed ones into f_xx. The data will go into d_2d, d_3d

    type(Deriv),private, dimension(MAXDEF_DERIV2D), save :: def_2d
    type(Deriv),private, dimension(MAXDEF_DERIV3D), save :: def_3d

    type(Deriv),public, allocatable, dimension(:), save :: f_2d
    type(Deriv),public, allocatable, dimension(:), save :: f_3d


  ! The 2-d and 3-d fields use the above as a time-dimension. We define
  ! LENOUTxD according to how fine resolution we want on output. For 2d
  ! fields we use daily outputs. For the big 3d fields, monthly output
  ! is sufficient.

   integer, public, parameter ::  LENOUT2D = 4  ! Allows INST..DAY for 2d fields
   integer, public, parameter ::  LENOUT3D = 4  ! Allows INST..MON for 3d fields

  !e.g. d_2d( num_deriv2d,MAXLIMAX, MAXLJMAX, LENOUT2D)
  ! &   d_3d( num_deriv3d,MAXLIMAX, MAXLJMAX, KMAX_MID, LENOUT3D )
   real, save, public, allocatable, dimension(:,:,:,:) :: d_2d
   real, save, public, allocatable, dimension(:,:,:,:,:) :: d_3d


   ! save O3 every hour during one day to find running max
    real, save,  public :: &     ! to be used for SOMO35
     D2_O3_DAY( MAXLIMAX, MAXLJMAX, NTDAY) = 0.


  ! Counters to keep track of averaging
  ! Initialise to zero in Init.

    integer, public, allocatable, dimension(:,:), save :: nav_2d
    integer, public, allocatable, dimension(:,:), save :: nav_3d

   ! Note - previous versions did not have the LENOUT2D dimension
   ! for wet and dry deposition. Why not?  Are annual or daily
   ! depositions never printed? Since I prefer to keep all 2d
   ! fields as similar as posisble, I have kept this dimension
   ! for now - ds


   !-- some variables for the VOC sum done for ozone models
   !   (have no effect in non-ozone models - leave in code)

   integer, private, save :: nvoc   ! No. VOCs
   integer, private, dimension(NSPEC_ADV), save :: &
             voc_index, &     ! Index of VOC in xn_adv
             voc_carbon       ! Number of C atoms

   logical, private, save :: debug_flag, Is3D
   character(len=100), private :: errmsg

   integer, private :: i,j,k,n, ivoc, index    ! Local loop variables

   contains

    !=========================================================================
    subroutine Init_Derived()

        integer :: alloc_err
          if(MasterProc .and. DEBUG ) write(*,*) "INIT My DERIVED STUFF"
          call Init_My_Deriv()  !-> wanted_deriv2d, wanted_deriv3d

         ! get lengths of wanted arrays (excludes notset values)
          num_deriv2d = LenArray(wanted_deriv2d,NOT_SET_STRING)
          num_deriv3d = LenArray(wanted_deriv3d,NOT_SET_STRING)

          call CheckStop(num_deriv2d<1,"num_deriv2d<1 !!")

     if ( num_deriv2d > 0 ) then
          if(DEBUG .and. MasterProc ) write(*,*) "Allocate arrays for 2d:",&
                                                  num_deriv2d
          allocate(f_2d(num_deriv2d),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of f_2d")
          allocate(d_2d(num_deriv2d,MAXLIMAX,MAXLJMAX,LENOUT2D),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of d_2d")
          call CheckStop(alloc_err,"Allocation of d_3d")
          allocate(nav_2d(num_deriv2d,LENOUT2D),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of nav_2d")
          nav_2d = 0
     end if
     if ( num_deriv3d > 0 ) then
          if(DEBUG .and. MasterProc ) write(*,*) "Allocate arrays for 3d: ",&
                                                 num_deriv3d
          allocate(f_3d(num_deriv3d),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of f_3d")
          allocate(d_3d(num_deriv3d,MAXLIMAX,MAXLJMAX,KMAX_MID,LENOUT3D),&
                 stat=alloc_err)
          allocate(nav_3d(num_deriv3d,LENOUT3D),stat=alloc_err)
          call CheckStop(alloc_err,"Allocation of nav_3d")
          nav_3d = 0
     end if

          call Define_Derived()
          call Setups()

    end subroutine Init_Derived

   !=========================================================================
    subroutine AddNewDeriv( name,class,subclass,txt,unit,index,f2d,&
           LC,Threshold,dt_scale,scale, avg,rho,inst,year,month,day,atw,Is3D)

       character(len=*), intent(in) :: name ! e.g. DDEP_SO2_m2Conif
       character(len=*), intent(in) :: class ! Type of data, e.g. ADV or VOC
       character(len=*), intent(in) :: subclass !
       character(len=*), intent(in) :: txt ! text where needed, e.g. "Conif"
       character(len=*), intent(in) :: unit ! writen in netCDF output
       integer, intent(in)  :: index    ! index in concentation array, or other
       integer, intent(in) :: f2d       ! index in f_2d arrays
       integer, intent(in) :: LC        ! Index of Receiver land-cover (one
       real, intent(in)    :: Threshold ! Threshold or CL, e.f. AOTx or AFstY
       logical, intent(in) :: dt_scale  !  where scaling by dt_advec needed,
       real, intent(in)    :: scale     !  e.g. use 100.0 to get cm/s
       logical, intent(in)  :: avg      ! True => average data (divide by
                           ! nav at end),  else accumulate over run period
       logical, intent(in)  :: rho      ! True when scale is ug (N or S)
       logical, intent(in)  :: inst     ! True when instantaneous values needed
       logical, intent(in)  :: year     ! True when yearly averages wanted
       logical, intent(in)  :: month    ! True when monthly averages wanted
       logical, intent(in)  :: day      ! True when daily averages wanted
                                        ! of Dep_Receivers)
       integer, intent(in)  :: atw           ! atomic weight where needed

       logical, intent(in), optional :: Is3D
       type(Deriv) :: inderiv

       inderiv = Deriv(trim(name),trim(class),trim(subclass),&
                           trim(txt),trim(unit),index,f2d,LC,Threshold,dt_scale, scale,&
                            avg,rho,inst,year,month,day,atw)

       if ( present(Is3D) ) then
          call AddDeriv(inderiv,Is3D)
       else
          call AddDeriv(inderiv)
       end if

    end subroutine AddNewDeriv
   !=========================================================================
    subroutine AddDeriv(inderiv,Is3D)

       type(Deriv), intent(in) :: inderiv
       logical, intent(in), optional :: Is3D
       logical :: Is3D_local

       Is3D_local = .false.
       if ( present(Is3D) ) Is3D_local = Is3D

       if ( Is3D_local ) then
         Nadded3d = Nadded3d + 1
         N = Nadded3d
         if ( DEBUG .and. MasterProc   ) &
                  write(*,*) "Define 3d deriv ", N, trim(inderiv%name)
         call CheckStop(N>MAXDEF_DERIV3D,"Nadded3d too big!")
         def_3d(N) = inderiv

       else
         Nadded2d = Nadded2d + 1
         N = Nadded2d
         if (  DEBUG .and. MasterProc ) write(*,*) "Define 2d deriv ",&
               N, trim(inderiv%name), " Class:", inderiv%class, " Ind:", inderiv%index
         call CheckStop(N>MAXDEF_DERIV2D,"Nadded2d too big!")
         def_2d(N) = inderiv

       end if

    end subroutine AddDeriv
   !=========================================================================
    subroutine Define_Derived()

   ! Set the parameters for the derived parameters, including the codes
   ! used by DNMI/xfelt and scaling factors. (The scaling factors may
   ! be changed later in Derived_ml.
   ! And, Initialise the fields to zero.

    real, save    :: ugSm3 = atwS*PPBINV/ATWAIR
    real, save    :: ugNm3 = atwN*PPBINV/ATWAIR
    real, save    :: ugCm3 = 12*PPBINV/ATWAIR
    real, save    :: ugXm3 = PPBINV/ATWAIR   ! will be multiplied by molwwt(X)
    real, save    :: ugPM  = PPBINV /ATWAIR  ! No multiplication needed

    character(len=30) :: dname
    character(len=3) :: subclass

  ! - for debug  - now not affecting ModelConstants version
  ! integer, dimension(MAXLIMAX) :: i_fdom
  ! integer, dimension(MAXLJMAX) :: j_fdom
   integer :: ind, iadv, itot, idebug, n, n2, iLC

  ! - And to check if a wanted field has been previously defined.
        integer, dimension(MAXDEF_DERIV2D) :: found_ind2d = 0
        integer, dimension(MAXDEF_DERIV3D) :: found_ind3d = 0


    if(DEBUG .and. MasterProc ) write(6,*) " START DEFINE DERIVED "
    !   same mol.wt assumed for PPM25 and PPMCOARSE

     !ds ugPMad = species(PPM25)%molwt * PPBINV /ATWAIR

!DSGC     ugCH3CHO = species ( CH3CHO )%molwt * PPBINV /ATWAIR

!-- Deposition fields. Define all possible fields and their xfelt codes here:

      !code class  avg? ind scale rho Inst Yr Mn Day   name      unit

Is3D = .false.
       !Deriv(name, class,    subc,  txt,           unit
       !Deriv index, f2d,LC, Threshold, dt_scale, scale, avg? rho Inst Yr Mn Day atw
call AddNewDeriv( "WDEP_PREC","PREC ","-","-", "mm",  &
                -1, -99,-99,  0.0,       F,    1.0,   F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_SOX ","WDEP ","-","-", "mgS/m2", &
                -1, -99,-99,  0.0,     F,    1.0e6,   F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_OXN ","WDEP ","-","-", "mgN/m2",  &
                -1, -99,-99,  0.0,     F,    1.0e6,   F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_RDN ","WDEP ","-","-", "mgN/m2",  &
                -1, -99,-99,  0.0,     F,    1.0e6,   F,   F , F ,T ,T ,T ,-999)
! Hard-coded for ECO08 - will rewrite later as with DDEP
!DONEcall AddDef( "WDEP ", F, -1, 1.0e6, F  , F  ,T ,T ,T ,"WDEP_SO2","mgS/m2")
!ixadv not yet used I think.
call AddNewDeriv( "WDEP_SO2 ","WDEP ","-","-", "mgS/m2", &
         IXADV_SO2, -99,-99,  0.0,     F,    1.0e6,   F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_SO4 ","WDEP ","-","-", "mgS/m2", &
         IXADV_SO4, -99,-99,  0.0,     F,    1.0e6,   F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_HNO3 ","WDEP ","-","-", "mgN/m2", &
         IXADV_HNO3, -99,-99,  0.0,     F,    1.0e6,  F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_pNO3_f","WDEP ","-","-", "mgN/m2", &
         IXADV_pNO3_f, -99,-99,  0.0,     F,    1.0e6,  F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_pNO3_c","WDEP ","-","-", "mgN/m2", &
         IXADV_pNO3_c, -99,-99,  0.0,     F,    1.0e6,  F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_NH3 ","WDEP ","-","-", "mgN/m2", &
         IXADV_NH3, -99,-99,  0.0,     F,    1.0e6,  F,   F , F ,T ,T ,T ,-999)
call AddNewDeriv( "WDEP_aNH4 ","WDEP ","-","-", "mgN/m2", &
         IXADV_aNH4, -99,-99,  0.0,     F,    1.0e6,  F,   F , F ,T ,T ,T ,-999)

      !code class  avg? ind scale rho Inst Yr Mn Day   name      unit

!ECO08:
! Compound-specific depositions:

! We process the various combinations of gas-species and ecosystem:
! stuff from My_Derived
! (e.g. for D2_SO2_m2SN)
!-------------------------------------------------------------------------------
  do n = 1, nMosaic
    call AddDeriv( MosaicOutput(n) )
  end do
!-------------------------------------------------------------------------------
!Much deleted
!  do n = 1, nOutVEGO3 ! NB Adv not used
!                 !code            avg?       ind             scale
!                 !  rho Inst Yr Mn Day  Units
!    call AddDef(OutVEGO3(n)%class, F,OutVEGO3(n)%Index,OutVEGO3(n)%scale,&
!                    F  , F  ,T ,T ,T , OutVEGO3(n)%name,OutVEGO3(n)%units)
!  end do
!-------------------------------------------------------------------------------
! Areas of deposition-related ecosystems. Set externally
  do n = 1, NDEF_ECOSYSTEMS
     print *, "ECODEF ",n, trim( DepEcoSystem(n)%name )
!
!    call AddDef("ECOAREA", F,DepEcoSystem(n)%Index,DepEcosystem(n)%scale,&
!                    F  , F  ,T ,F ,F , DepEcosystem(n)%name,&
!                                        DepEcosystem(n)%units)
     call AddDeriv( DepEcoSystem(n) )
  end do
!!-------------------------------------------------------------------------------
!-- 2-D fields - the complex ones
! (multiplied with roa in layers?? ==>  rho "false" ) !ds - explain!

      !Deriv(name, class,    subc,  txt,           unit
      !Deriv index, f2d,LC, Threshold, dt_scale, scale, avg? rho Inst Yr Mn Day atw
call AddNewDeriv( "AOT40_Grid", "AOT ","subclass","-", "ppb h", &
           IXADV_O3, -99,0, 40.0,   T, 1.0/3600.0, F,  F , F ,T ,T ,T ,-999)
!
!-------------------------------------------------------------------------------
!-- Tropospheric columns
!
  do n = 1, size(COLUMN_MOLEC_CM2)
  do n2 = 1, size(COLUMN_LEVELS)

    itot = COLUMN_MOLEC_CM2(n)
    iadv = itot - NSPEC_SHL
      !Deriv(name, class,    subc,  txt,           unit
      !Deriv index, f2d,LC, Threshold, dt_scale, scale, avg? rho Inst Yr Mn Day atw
     dname = "COLUMN_" //trim(species( itot )%name) //"_"//COLUMN_LEVELS(n2)
     subclass = COLUMN_LEVELS(n2)
     read(unit=subclass(2:3),fmt="(i2)") iLC  ! Faking vertical index with iLC :-(
     call AddNewDeriv( dname, "COLUMN ",subclass,"-", "molec/cm2", &
               iadv, -99,iLC,  0.0,   F,    1.0,     T,  F , F ,T ,T ,T ,-999)
  end do
  end do
!-------------------------------------------------------------------------------
!
!       code class  avg? ind scale rho  Inst  Yr  Mn   Day  name      unit

!REWRITE as groups!
!call AddDef( "NOZ  ", T,   -1  ,ugN ,    T , F,T,T,T,"D2_NOZ","ugN/m3")
!call AddDef( "OX   ", T,   -1  ,PPBINV , F , F,T,T,T,"D2_OX","ppb")

       !Deriv(name, class,    subc,  txt,           unit
      !Deriv index, f2d,LC, Threshold, dt_scale, scale, avg? rho Inst Yr Mn Day atw

! NOT YET: Scale pressure by 0.01 to get hPa
call AddNewDeriv( "PSURF ","PSURF",  "SURF","-",   "hPa", &
               -99,  -99,-99, 0.0,     F,  1.0,  T,  F , T,   T, T, T , -999)

call AddNewDeriv( "HMIX  ","HMIX",  "-","-",   "m", &
               -99,  -99,-99, 0.0,     F,  1.0,  T, T , F, T, T, T ,-999)

!call AddDef( "HMIX00",T,  0 ,         F,  1.0,  T , F, T, T, T ,"D2_HMIX00","m")
!call AddDef( "HMIX12",T,  0 ,         F,  1.0,  T , F, T, T, T ,"D2_HMIX12","m")

call AddNewDeriv( "USTAR_NWP","USTAR_NWP",  "-","-",   "m/s", &
               -99,  -99,-99, 0.0,   F, 1.0,  T, T , F, T, T, T ,-999)
!SoilWater
call AddNewDeriv( "SoilWater_deep","SoilWater_deep",  "-","-",   "m", &
               -99,  -99,-99, 0.0,   F, 1.0,  T, T , F, T, T, T ,-999)
call AddNewDeriv( "SoilWater","SoilWater",  "-","-",   "m", &
               -99,  -99,-99, 0.0,   F, 1.0,  T, T , F, T, T, T ,-999)

!Mass-based outputs
do ind = 1, size(SURF_UG_S)
  itot = SURF_UG_S(ind)
  iadv = itot - NSPEC_SHL
  dname = "SURF_ugS_" //species( itot )%name
       !Deriv(name, class,    subc,  txt,           unit
      !Deriv index, f2d,LC, Threshold, dt_scale, scale, avg? rho Inst Yr Mn Day atw
  call AddNewDeriv( dname, "SURF_UG", "MASS", "-", "ugS/m3", &
              iadv , -99, -99, 0.0,   F,   ugSm3,     T, T , F, T, T, T, -999 ) !?? atw?
end do

do ind = 1, size(SURF_UG_N)
  itot = SURF_UG_N(ind)
  iadv = itot - NSPEC_SHL
  dname = "SURF_ugN_" //species( itot )%name
  call AddNewDeriv( dname, "SURF_UG", "MASS", "-", "ugN/m3", &
         iadv , -99, -99, 0.0, F,   ugNm3,   T, T , F, T, T, T, -999 ) !?? atw?
end do

do ind = 1, size(SURF_UG_C)
  itot = SURF_UG_C(ind)
  iadv = itot - NSPEC_SHL
  dname = "SURF_ugC_" //species( itot )%name
  call AddNewDeriv( dname, "SURF_UG", "MASS", "-", "ugC/m3", &
         iadv , -99, -99, 0.0, F,  ugCm3,   T, T , F, T, T, T, -999 ) !?? atw?
end do

do ind = 1, size(SURF_UG)
  itot = SURF_UG(ind)
  iadv = itot - NSPEC_SHL
  dname = "SURF_ug_" //species( itot )%name
  call AddNewDeriv( dname, "SURF_UG", "MASS", "-", "ug/m3", &
         iadv , -99, -99, 0.0, F, ugXm3*species(itot)%molwt,  T, T , F, T, T, T, -999 ) !?? atw?
end do

do ind = 1, size(SURF_PPB)   ! ppb has rho flag set false
  itot = SURF_PPB(ind)
  iadv = itot - NSPEC_SHL
  dname = "SURF_ppb_" //species( itot )%name
  call AddNewDeriv( dname, "SURF_PPB", "-", "-", "ppb", &
         iadv , -99, -99, 0.0, F, PPBINV,  T, F , F, T, T, T, -999 ) !?? atw?
end do

!call AddDef( "VOC  ", T,  -1    ,PPBINV, F, F, T, T, T,"D2_VOC","ppb")
call AddNewDeriv( "SURF_ppbC_VOC", "VOC", "-", "-", "ppb", &
         -1 , -99, -99, 0.0, F, PPBINV,  T, F , F, T, T, T, -999 ) !?? atw?

!Emissions:
! We use mg/m2 outputs for consistency with depositions
! Would need to multiply by GridArea_m2 later to get ktonne/grid, but not
! done here.
!
! BVOC called every dt_advec, so use dt_scale=1.0e6 to get from kg/m2/s to
!  mg/m2 accumulated (after multiplication by dt_advec)

     !Deriv(name,       class,    subc,  txt,           unit
     !Deriv index, f2d,LC, Threshold, dt_scale, scale, avg? rho Inst Yr Mn Day atw

  do  ind = 1, size(BVOC_GROUP)
     itot = BVOC_GROUP(ind)
     dname = "Emis_mgm2_" // trim(species(itot)%name)
     call AddNewDeriv( dname, "NatEmis", "-", "-", "mg/m2", &
                 ind , -99,-99, 0.0,  T ,    1.0e6,     F, F , F, T, T, T, -999 ) !?? atw?
   end do

! SNAP emissions called every hour, so use scale=3600.0 to get from kg/m2/s to kg/m2,
! and by 1.0e6 to get from kg/m2 to mg/m2 accumulated.
!
! Future option - might make use of Emis_Molwt to get mg(N)/m2
do  ind = 1, size(EMIS_NAME)
  dname = "Emis_mgm2_" // trim(EMIS_NAME(ind))
  call AddNewDeriv( dname, "SnapEmis", "-", "-", "mg/m2", &
                     ind , -99,-99, 0.0, F,  3600.0e6,  F, F , F, T, T, T, -999 ) !?? atw?
end do ! ind


!call AddDef( "SNOW",T,  0 ,       1.0, F , T, T, T, T ,"D2_SNOW","frac")
!surface 0.6*SO2/NH3 ratio where conc are 24 averaged
!call AddDef( "SNratio",T,  0 ,       1.0, F , T, T, T, T ,"D2_SNratio","ratio")
!
! Also, use field for EU definition (8-20 CET) and Mapping Manual/UNECE
!
! --  time-averages - here 8-16
!
!call AddDef( "TADV ", T,IXADV_HCHO  ,ugHCHO,  T, F, T, T, T,"D2T_HCHO","ug/m3")
!DSGCcall AddDef( "TADV ", T,IXADV_CH3CHO,ugCH3CHO,T, F, T, T,T,"D2T_CH3CHO","ug/m3")
!
! -- miscellaneous user-defined functions
!
!      code class   avg? ind scale rho Inst Yr Mn  Day   name      unit
!! ,Deriv( "TSO4 ", T,   -1  ,ugS ,    T , F,T,T,T,"D2_SOX","ugS/m3")
!call AddDef( "FRNIT", T,   -1  ,1.0 ,    F , F,T,T,T,"D2_FRNIT","(1)")
!call AddDef( "MAXSHL", F,IXSHL_OH,1.0e13,F , F,T,F,T,"D2_MAXOH","?")
!call AddDef( "H2O  ", T, -1,   1.0 , T, F, T, T, T,"D2_PM25_H2O ", "ug/m3")
!call AddDef( "SSalt", T, -1, ugSS,   T, F, T, T, T,"D2_SS  ", "ug/m3")

              !Deriv(name, class,    subc,  txt,           unit
              !Deriv index, f2d,LC, Threshold, scale, avg? rho Inst Yr Mn Day atw
call AddNewDeriv("SURF_ug_SIA", "SIAGROUP", "MASS", "-", "ug/m3", &
                      -99 , -99,-99, 0.0, F, ugPM,  T, T , F, T, T, T, -999 ) !?? atw?
call AddNewDeriv("SURF_ugN_OXN", "OXNGROUP", "MASS", "-", "ugN/m3", &
                      -99 , -99,-99, 0.0, F, ugNm3,  T, T , F, T, T, T, -999 ) !?? atw?
call AddNewDeriv("SURF_ugN_NOX", "NOXGROUP", "MASS", "-", "ugN/m3", &
                      -99 , -99,-99, 0.0, F, ugNm3,  T, T , F, T, T, T, -999 ) !?? atw?
call AddNewDeriv("SURF_ugN_RDN", "RDNGROUP", "MASS", "-", "ugN/m3", &
                      -99 , -99,-99, 0.0, F, ugNm3,  T, T , F, T, T, T, -999 ) !?? atw?
call AddNewDeriv("SURF_ugN_TNO3", "TNO3GROUP", "MASS", "-", "ugN/m3", &
                      -99 , -99,-99, 0.0, F, ugNm3,  T, T , F, T, T, T, -999 ) !?? atw?
call AddNewDeriv("SURF_ug_PM25", "PM25GROUP", "MASS", "-", "ug/m3", &
                      -99 , -99,-99, 0.0, F, ugPM,  T, T , F, T, T, T, -999 ) !?? atw?
call AddNewDeriv("SURF_ug_PMc ", "PMcGROUP", "MASS", "-", "ug/m3", &
                      -99 , -99,-99, 0.0, F, ugPM,  T, T , F, T, T, T, -999 ) !?? atw?
call AddNewDeriv("SURF_ug_PM10", "PM10GROUP", "MASS", "-", "ug/m3", &
                      -99 , -99,-99, 0.0, F, ugPM,  T, T , F, T, T, T, -999 ) !?? atw?
call AddNewDeriv( "SOMO35","SOMO",  "SURF","-",   "ppb.day", &
                  IXADV_O3, -99,-99,35.0, F, 1.0,   F,  F , F,   T, T, F , -999)
call AddNewDeriv( "SOMO00","SOMO",  "SURF","-",   "ppb.day", &
                  IXADV_O3, -99,-99, 0.0, F, 1.0,   F,  F , F,   T, T, F , -999)

call AddNewDeriv( "SURF_MAXO3","MAXADV", "O3","-",   "ppb", &
           IXADV_O3, -99,-99, 0.0, F, PPBINV,   F,  F , F,   T, T, T , -999)


!-- 3-D fields

Is3D = .true.
!call AddDef( "TH  ",T,  0 ,       1.0, F , T, T, T, F ,"D3_TH","m",Is3D)
! etc... use D3_PPB sustem
!
! Set Year true to allow debug - change later
!call AddDef( "SHL",   T, IXSHL_OH,  PPTINV, T, F, T, T, F ,"D3_OH","?",Is3D)
!DSGC call AddDef( "ADV",   T, IXADV_CH3COO2, &
!DSGC                                     PPTINV, F, F, T, T, F ,"D3_CH3COO2","?",Is3D)
!call AddDef( "MAX3DSHL", T,IXSHL_OH,PPTINV, T, F, T, T, F ,"D3_MAXOH","?",Is3D)   ! rho true for shl
!DSGC call AddDef( "MAX3DADV", T, IXADV_CH3COO2,&
!DSGC                                     PPTINV, F, F, T, T, F ,"D3_MAXCH3COO2","?",Is3D)
!DSGC call AddDef( "PHNO3   ", T, IXSHL_PHNO3,1.0e8, F, F, T, T, F ,"D3_PHNO3","?",Is3D)

! hb new 3D output
do ind = 1, size(D3_PPB)
  itot = D3_PPB(ind)
  iadv = itot - NSPEC_SHL
  dname = "D3_ppb_" //species( itot )%name
   !Deriv(name, class,    subc,  txt,           unit
   !Deriv index, f2d,LC, Threshold, scale, avg? rho Inst Yr Mn Day atw, Is3D
  call AddNewDeriv( dname, "D3_PPB", "-", "-", "ppb", &
         iadv , -99, -99, 0.0, F,  PPBINV,   T, T , F, T, T, F, -999,Is3D ) !?? atw?
end do
! hb new 3D output
!ds: Should not have PM in ppb! We do not know the Mol. Wt
!ds  others removed for now to save memory
!ds  call AddNewDeriv( "D3_PM25","PM25", "-","-",   "ppb", &
!ds           -1, -99,-99, 0.0, F,  PPBINV,  F,  F , F,   F, F, F , -999,Is3D)
!ds  call AddNewDeriv( "D3_PMco","PMco", "-","-",   "ppb", &
!ds           -1, -99,-99, 0.0, F,  PPBINV,  F,  F , F,   F, F, F , -999,Is3D)
!ds  call AddNewDeriv( "D3_TH","TH", "-","-",   "m", &
!ds           0, -99,-99, 0.0, F,  1.0,      F,  F , F,   F, F, F , -999,Is3D)
!ds  call AddNewDeriv( "D3_Kz","Kz", "-","-",   "-", &
!ds           0, -99,-99, 0.0, F,  1.0,      F,  F , F,   F, F, F , -999,Is3D)

!AMVB 2010-07-21: PM-PPB bug fix
do ind = 1, size(D3_OTHER)
  select case ( trim(D3_OTHER(ind)) )
  case ("D3_ug_PM25")
  call AddNewDeriv("D3_ug_PM25", "PM25GROUP", "MASS", "-", "ug/m3", &
         -99, -99,-99, 0.0, F, ugPM,  T, T , F, T, T, F, -999,Is3D ) !?? atw?
  case ("D3_ug_PMc")
  call AddNewDeriv("D3_ug_PMc ", "PMcGROUP" , "MASS", "-", "ug/m3", &
         -99, -99,-99, 0.0, F, ugPM,  T, T , F, T, T, F, -999,Is3D ) !?? atw?
  case ("D3_m_TH")
  call AddNewDeriv("D3_m_TH","TH", "-","-",   "m", &
         -99, -99,-99, 0.0, F,  1.0,  F, F , F, T, T, F, -999,Is3D )
  case ("D3_m2s_Kz")
  call AddNewDeriv( "D3_Kz","Kz", "-","-",   "-", &
         -99, -99,-99, 0.0, F,  1.0,  F, F , F, T, T, F, -999,Is3D )
  end select
end do

!AMVB 2010-08-02: SOURCE_RECEPTOR in FORECAST mode
     if ( .not. FORECAST .and. &
          SOURCE_RECEPTOR .and. num_deriv2d>0 ) then  ! We assume that no
!    if ( SOURCE_RECEPTOR .and. num_deriv2d>0 ) then  ! We assume that no
        def_2d(:)%day = .false.               ! daily outputs are wanted.
     end if

     if ( FORECAST .and. num_deriv2d>0 ) then ! only dayly (and hourly) output
        def_2d(:)%inst  = .false.             ! on FORECAST mode
        def_2d(:)%year  = .false.
        def_2d(:)%month = .false.
     end if
     if ( FORECAST .and. num_deriv3d>0 ) then
        def_3d(:)%inst  = .false.
        def_3d(:)%year  = .false.
        def_3d(:)%month = .false.
     end if

     ! Get indices of wanted fields in larger def_xx arrays:
      do i = 1, num_deriv2d
          ind = find_index( wanted_deriv2d(i), def_2d(:)%name )
          if ( ind>0) then
               f_2d(i) = def_2d(ind)
               call CheckStop ( found_ind2d(ind) > 0,  &
                  "REQUESTED 2D DERIVED ALREADY DEFINED: "// &
                      def_2d(ind)%name  )
               found_ind2d(ind)  = 1

          else
            print *,   "OOOPS wanted_deriv2d not found: ", wanted_deriv2d(i)
            print *,   "OOOPS N,N :", num_deriv2d, Nadded2d
            if(  MasterProc  ) then
               do idebug = 1, Nadded2d
                write (*,"(a,i4,a)") "Had def_2d: ", idebug, def_2d(idebug)%name
               end do
               call CheckStop("OOPS STOPPED")
            end if
          end if
          if (  DEBUG .and. MasterProc  ) then
               write(*,"(a,i4,a,i4,3a)") "Index f_2d ", i, " = def ", ind, &
                 trim(def_2d(ind)%name),trim(def_2d(ind)%unit), &
                    trim(def_2d(ind)%class)
          end if
      end do

      do i = 1, num_deriv3d
          ind = find_index( wanted_deriv3d(i), def_3d(:)%name )
          call CheckStop ( found_ind3d(ind) > 0,  &
                  "REQUESTED 3D DERIVED ALREADY DEFINED: "// &
                      def_3d(ind)%name  )
          found_ind3d(ind)  = 1
          f_3d(i) = def_3d(ind)
          if ( DEBUG .and. MasterProc ) write(*,*) "Index f_3d ", i,&
                " = def ", ind
      end do


   !Initialise to zero

      if ( num_deriv2d > 0  ) d_2d( :,:,:,:) = 0.0
      if ( num_deriv3d > 0  ) d_3d( :,:,:,:,:) = 0.0

      debug_flag = ( DEBUG  .and. debug_proc )

  end subroutine Define_Derived
 !=========================================================================
  subroutine Setups()
     integer :: n

    !/** flexibility note. By making use of character-based tests such
    !    as for "VOC" below, we achieve code which can stay for
    !    different chemical mechanisms, without having to define non-used indices.

    !/ ** if voc wanted, set up voc_array. Works for all ozone chemistries

      if ( any(  f_2d(:)%class == "VOC" ) ) then !TMP .or. &
          !TMP           any(  f_3d(:)%class == "VOC" )  ) then
          ! was call Setup_VOC(), moved here Mar 2010
          !--------------------------------------------------------
          ! Searches through the advected species and colects the
          ! index and carbon content of nmhc/voc species, as they are
          ! defined in CM_ChemSpecs_ml
          !
          !--------------------------------------------------------

         !====================================================================
          do n = 1, NSPEC_ADV

            if ( species( NSPEC_SHL+n )%carbons > 0 .and. &
                 species( NSPEC_SHL+n )%name   /= "CO"  .and. &
                 species( NSPEC_SHL+n )%name   /= "CH4" ) then

                 nvoc = nvoc + 1
                 voc_index(nvoc) = n
                 voc_carbon(nvoc) = species( NSPEC_SHL+n )%carbons
            end if
          end do
         !====================================================================
          if (DEBUG  .and. MasterProc )then
               write(6,*) "Derived VOC setup returns ", nvoc, "vocs"
               write(6,"(a12,/,(20i3))")  "indices ", voc_index(1:nvoc)
               write(6,"(a12,/,(20i3))")  "carbons ", voc_carbon(1:nvoc)
          endif
      end if


    end subroutine Setups
    !=========================================================================

    subroutine Derived(dt,End_of_Day)

    !/** DESCRIPTION
    !  Integration and averaging of chemical fields. Intended to be
    !  a more flexible version of the old chemint routine.
    !  Includes AOT40, AOT60 if present

      real, intent(in)    :: dt           !  time-step used in intergrations
      logical, intent(in) :: End_of_Day   !  e.g. 6am for EMEP sites

      character(len=len(f_2d%class)) :: typ  !  See defs of f_2d
      real :: thour                          ! Time of day (GMT)
      real :: timefrac                       ! dt as fraction of hour (3600/dt)
      real :: dayfrac              ! fraction of day elapsed (in middle of dt)
      real :: af
      real, save :: km2_grid
      integer :: ntime                       ! 1...NTDAYS
      integer :: klow                        !  lowest extent of column data
      real, dimension(MAXLIMAX,MAXLJMAX) :: density !  roa (kgair m-3 when
                                                    ! scale in ug,  else 1
      real, dimension(MAXLIMAX,MAXLJMAX) :: tmpwork

      real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: inv_air_density3D
                ! Inverse of No. air mols/cm3 = 1/M
                ! where M =  roa (kgair m-3) * MFAC  when ! scale in ug,  else 1
      logical :: accumulate_2dyear !flag to know when to accumulate d_2d (case "EXT")
      logical :: first_call = .true.
      integer :: ipm25, ipmc ! will save some calcs for pm10

      timefrac = dt/3600.0
      thour = current_date%hour+current_date%seconds/3600.0

      daynumber=day_of_year(current_date%year,current_date%month,&
                             current_date%day)


     !/***** 2-D fields **************************

     ipm25 = 0  ! Reset once pm25 calculated
     ipmc  = 0  ! pm-coarse
     do n = 1, num_deriv2d

        accumulate_2dyear=.true.
        typ = f_2d(n)%class
        if( debug_flag .and. first_call ) &
           write(*,"(a,i3,7a)") "Derive2d-typ",&
            n, "T:", typ, "N:", f_2d(n)%name, "C:", f_2d(n)%class,"END"


        if ( f_2d(n)%rho ) then
            forall ( i=1:limax, j=1:ljmax )
                density(i,j) = roa(i,j,KMAX_MID,1)
            end forall
        else
            density(:,:) = 1.0
        end if

        !/** user-defined time-averaging. Here we have defined TADV and TVOC
        !    so that 8-hour daytime averages will be calculated.
        !    Just comment out if not wanted, or (better!) don't define any
        !    f_2d as TADV or TVOC

        if ( typ == "TADV" .or. typ == "TVOC" ) then
             if(thour <= 8.0 .or. thour > 16.0 ) cycle  ! Start next species
        end if

       ! hmix average at 00 and 12:

        if ( typ == "HMIX00" .or. typ == "XKSIG00" ) then
             if(thour /= 0.0 ) cycle  ! Start next species
        end if

        if ( typ == "HMIX12" .or. typ == "XKSIG12" ) then
             if(thour /= 12.0 ) cycle  ! Start next species
        end if

        index = f_2d(n)%index
        if ( DEBUG .and. MasterProc .and. first_call ) then
           write(*,*) "DEBUG Derived 2d", n, f_2d(n)%name, index, trim(typ)
        end if

        select case ( typ )
          case ( "USTAR_NWP" )
            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = ustar_nwp(i,j)
          end forall

          case ( "SoilWater_deep" )
            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = SoilWater_deep(i,j,1)
          end forall
          if ( debug_flag ) call write_debug(n,index, "SoilWater_DEEP")

          case ( "SoilWater" )
            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = SoilWater_deep(i,j,1)
          end forall

          case ( "SNOW" )
            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = Grid_snow(i,j)
            end forall

          case ( "SNratio" )
            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = min(3.0,so2nh3_24hr(i,j))
            end forall

          case ( "PSURF" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = ps(i,j,1)*0.01
              !NOT YET - keep hPa in sites:d_2d( n, i,j,IOU_INST) = ps(i,j,1)
            end forall


          case ( "HMIX", "HMIX00", "HMIX12" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = pzpbl(i,j)
            end forall

            if ( debug_flag ) then
             write(*,fmt="(a12,i4,f12.3)") "HMIX" , n , &
                     d_2d(n,debug_li,debug_lj,IOU_INST)
            end if

         ! Simple advected species:
          case ( "ADV", "TADV" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = xn_adv(index,i,j,KMAX_MID)  &
                                     * cfac(index,i,j) * density(i,j)
            end forall

            if ( debug_flag ) call write_debugadv(n,index, &
                                     density(debug_li,debug_lj), "ADV or TADV")

          case ( "SURF_PPB" )
            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = xn_adv(index,i,j,KMAX_MID) &
                                     * cfac(index,i,j)
            end forall
            if ( debug_flag ) call write_debugadv(n,index, 1.0, "PPB OUTS")

          case ( "SURF_UG" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = xn_adv(index,i,j,KMAX_MID) &
                                     * cfac(index,i,j) * density(i,j)
            end forall
            if ( debug_flag ) call write_debugadv(n,index, &
                                     density(debug_li,debug_lj), "SURF_UG")

          case ( "H2O" )      !water

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = PM_water(i,j,KMAX_MID)
            end forall


          case ( "MAXADV" )


            d_2d( n, 1:limax,1:ljmax,IOU_DAY) = &
                 max( d_2d( n, 1:limax,1:ljmax,IOU_DAY),  &
                      xn_adv(index,1:limax,1:ljmax,KMAX_MID)  &
                     * cfac(index,1:limax,1:ljmax) * density(1:limax,1:ljmax))
            if ( debug_flag ) call write_debugadv(n,index, &
                                     density(debug_li,debug_lj), "MAXADV")


            !Monthly and yearly ARE averaged over days
            if(End_of_Day)then
              d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)
              nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
              if(    current_date%month >= 4 &
                 .or.current_date%month <= 9 )then
              d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY)
              nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
              endif
            endif


          case ( "MAXSHL" )        ! Daily maxima - short-lived

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_DAY) = max( d_2d( n, i,j,IOU_DAY), &
                                xn_shl(index,i,j,KMAX_MID)  &
                                    / (density(i,j)*MFAC) )
            end forall


            if ( debug_flag ) then
               write(*, *) "SHL:MAX.,MFAC ", n, index  , MFAC
               write(*,fmt="(a12,2i4,4es12.3)") "SHL MAX. ", n, index  &
                      , d_2d(n,debug_li,debug_lj,IOU_DAY) &
                      ,  xn_shl(index,debug_li,debug_lj,KMAX_MID)  &
                      ,  density(debug_li,debug_lj), MFAC
            end if

            !Monthly and yearly ARE averaged over days
            if(End_of_Day)then
              d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)
              nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
              if(    current_date%month >= 4 &
                 .or.current_date%month <= 9 )then
              d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY)
              nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
              endif
            endif

          case ( "VOC", "TVOC" )

            call voc_2dcalc()

          case( "AOT" )    !  Hardly used these days. The vegetation-specific
                           !  AOTs are handled in the Mosaic class and as
                           !  part of the dry dep calculations.

            d_2d(n, 1:limax, 1:ljmax, IOU_INST) = Calc_GridAOTx(f_2d(n)%Threshold)

           !if( debug_flag .and. i == debug_li .and. j == debug_lj ) then
           if( debug_flag ) then
              write(*,*) "GROWINDERIV? ", n, f_2d(n)%name
           end if


           case( "SOMO" )


              !dt/7200: half a dt time step in hours
              !dayfrac "points" to the middle of the integration step
              dayfrac= (thour-(dt/7200.))/24. !must be < 1
              ntime=int(dayfrac*NTDAY )+1 !must be >=1 and <= NTDAY
              if(dayfrac<0)ntime=NTDAY !midnight

        !last value  (not averaged):
          D2_O3_DAY( : , : , ntime) =&
           xn_adv(IXADV_O3,:,:,KMAX_MID)*cfac(IXADV_O3,:,:)*PPBINV

              if(dayfrac<0)then !only at midnight: write on d_2d


                 call somo_calc( n, DEBUG .and. debug_proc ) !  accumulate
                 d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_DAY)

                ! if(current_date%month>=4.and.current_date%month<=9)then
                 d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_DAY)
                !NB overwritten anyway D2_O3_DAY = 0.
              endif


          case ( "PREC", "WDEP", "DDEP", "VG" ,"Rs", "Rns", "Gns", "Mosaic" )
!            if ( debug_flag ) write(*,"(2a,i4,a,es12.3)")"PROCESS ",trim(typ),&
!                   n, trim(f_2d(n)%name), d_2d(n,debug_li,debug_lj,IOU_INST)
!            Nothing to do - all set in My_DryDep

          case ( "COLUMN" )   ! Dave May 2009, simplified from Joffen's COLUMN
                            ! MFAC gives #/cm3, 100 is for m -> cm

            klow = f_2d(n)%LC  ! here we have used LC to set vertical limit
            do j = 1, ljmax
            do i = 1, limax
               k = 1
               tmpwork( i,j ) =  &
                    xn_adv(index,i,j,k)  * &
                    roa(i,j,k,1)* (z_bnd(i,j,k)-z_bnd(i,j,k+1))

               do k = 2, klow   !!! KMAX_MID
                  tmpwork( i,j ) = tmpwork( i,j ) + &
                    xn_adv(index,i,j,k)  * &
                    roa(i,j,k,1) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))
                  if( debug_flag .and. i == debug_li .and. j == debug_lj ) then
                     write(*,"(a,3i4,a4,f8.3,f8.1,2es12.3)") &
                        trim(f_2d(n)%name), n, index, k, " => ", &
                        roa(i,j,k,1), z_bnd(i,j,k)-z_bnd(i,j,k+1), &
                        xn_adv(index,i,j,k),tmpwork(i,j)
                  end if
               end do ! k
               d_2d( n, i, j, IOU_INST) = MFAC * 100.0 * tmpwork( i, j ) ! Should be molec/cm2


            end do !i
            end do !j
            if( debug_flag ) &
                write(*,"(a18,es12.3)") "COLUMN d2_2d", d_2d( n, debug_li, debug_lj, IOU_INST)


          !case ( "ECOAREA" )
         case ( "EcoFrac" ) ! ODD TO HAVE FRAC AND AREA BELOW:"ECOAREA" )

            if( .not. first_call ) cycle ! Only need to do once
            if( f_2d(n)%Index == FULL_GRID ) then
              km2_grid = (GRIDWIDTH_M*GRIDWIDTH_M) * 1.0e-6 ! km2
              forall ( i=1:limax, j=1:ljmax )
                  d_2d(n,i,j,IOU_YEAR) =  EcoSystemFrac( f_2d(n)%Index ,i,j)&
                        * KM2_GRID /xm2(i,j)
              end forall
            else
              forall ( i=1:limax, j=1:ljmax )
                  d_2d(n,i,j,IOU_YEAR) =  EcoSystemFrac( f_2d(n)%Index ,i,j)
              end forall
            end if
            if( debug_flag ) &
                write(*,"(a18,a,i4,a,2es12.3)") "ECOD2D ", &
                 " f2d:", f_2d(n)%Index, &
                 " Frac", EcoSystemFrac( f_2d(n)%Index, debug_li,debug_lj), &
                 !!" Index: ", DepEcoSystem(n)%Index, &
                    d_2d( n, debug_li, debug_lj, IOU_YEAR)

          case ( "NatEmis" ) !emissions in kg/m2/s converted??

              forall ( i=1:limax, j=1:ljmax )
                  d_2d(n,i,j,IOU_INST) =  EmisNat( i,j, f_2d(n)%Index)
              end forall
              !Not done, keep mg/m2  * GridArea_m2(i,j)

          case ( "SnapEmis" ) !emissions in kg/m2/s converted??

              forall ( i=1:limax, j=1:ljmax )
                  d_2d(n,i,j,IOU_INST) =  SumSnapEmis( i,j, f_2d(n)%Index)
              end forall
              !not done, to keep mg/m2 * GridArea_m2(i,j)


          case ( "EXT" )

          ! Externally set for IOU_INST (in other routines); so no new work
          ! needed except decision to accumalate to yearly or not.
          ! Used for e.g. AOT40s

             call setaccumulate_2dyear(f_2d(n)%name,accumulate_2dyear)

            !if ( debug_flag ) write(*,"(a18,i4,a12,a4,es12.3)")"EXT d_2d",&
            !       n, f_2d(n)%name, " is ", d_2d(n,debug_li,debug_lj,IOU_INST)

          case ( "SIAGROUP" )
            call uggroup_calc( d_2d(n,:,:,IOU_INST), n, typ, SIA_GROUP, &
                               density, 0)
          case ( "OXNGROUP" )
            call uggroup_calc( d_2d(n,:,:,IOU_INST), n, typ, OXN_GROUP, &
                               density, 0)
          case ( "NOXGROUP" )
            call uggroup_calc( d_2d(n,:,:,IOU_INST), n, typ, NOX_GROUP, &
                               density, 0)
          case ( "RDNGROUP" )
            call uggroup_calc( d_2d(n,:,:,IOU_INST), n, typ, RDN_GROUP, &
                               density, 0)
          case ( "TNO3GROUP" )
            call uggroup_calc( d_2d(n,:,:,IOU_INST), n, typ, TNO3_GROUP, &
                               density, 0)
          case ( "PM25GROUP" )
            ipm25 = n
            call uggroup_calc( d_2d(n,:,:,IOU_INST), n, typ, PM25_GROUP, &
                               density, 0)
          case ( "PMcGROUP" )
            ipmc = n
            call uggroup_calc( d_2d(n,:,:,IOU_INST), n, typ, PMCO_GROUP, &
                               density, 0)

          case ( "PM10GROUP" ) ! Consider doing as sum later
            if ( ipm25 > 0 .and. ipmc > 0 ) then ! We have these already
               d_2d(n,:,:,IOU_INST) = &
                 d_2d(ipm25,:,:,IOU_INST) + d_2d(ipmc,:,:,IOU_INST)
            else
              ! call uggroup_calc( d_2d(n,:,:,IOU_INST), n, typ,PM10_GROUP, density )
               call StopAll("PM10 group is special. Need to define PM25 and PMc first!")
            end if

          case  default

            if ( debug_flag ) then
               if( debug_flag .and. i == debug_li .and. j == debug_lj ) &
                 write(*,"(a,i3,4a)") "My_Deriv Defaults called n=",&
                    n, " Type ",trim(typ), " Name ", trim( f_2d(n)%name )

                 write(*,"(a,i3,i8,i4,a)") &
                    "My_Deriv index?, nav? length?, class? ", index,&
                    nav_2d(n,IOU_INST), len(f_2d%class), trim(f_2d(n)%class)
                 write(*,*) "My_Deriv index?, avg ", f_2d(n)%avg
             end if

             call My_DerivFunc( d_2d(n,:,:,IOU_INST), typ, density )

        end select


        !/** add to daily, monthly and yearly average, and increment counters
        !  Note that the MAXADV and MAXSHL and SOMO needn't be summed here, but
        !  since the INST values are zero it doesn't harm, and the code is
        !  shorter. These d_2d ( MAXADV, MAXSHL, SOMO) are set elsewhere

        af = 1.0 ! accumlation factor
        if( f_2d(n)%dt_scale ) then !need to scale with dt_advec
            af = dt_advec
        end if

        d_2d(n,:,:,IOU_DAY )  = d_2d(n,:,:,IOU_DAY )  + af*d_2d(n,:,:,IOU_INST)
        if ( f_2d(n)%avg ) nav_2d(n,IOU_DAY) = nav_2d(n,IOU_DAY) + 1
        d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + af*d_2d(n,:,:,IOU_INST)
        if ( f_2d(n)%avg ) nav_2d(n,IOU_MON) = nav_2d(n,IOU_MON) + 1
        if(accumulate_2dyear)then
           d_2d(n,:,:,IOU_YEAR ) = &
                 d_2d(n,:,:,IOU_YEAR ) + af*d_2d(n,:,:,IOU_INST)
           if ( f_2d(n)%avg ) nav_2d(n,IOU_YEAR) = nav_2d(n,IOU_YEAR) + 1
        endif

     end do   ! num_deriv2d

     !/***** 3-D fields **************************

       if(debug_flag) then ! RUN through indices etc.
            write(*, "(a12,2i4,f12.3)") "3D3D TIME ",  me, num_deriv3d, &
                     (current_date%hour+current_date%seconds/3600.0)
        end if


     do n = 1, num_deriv3d

        index = f_3d(n)%index

       if ( f_3d(n)%rho ) then
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                inv_air_density3D(i,j,k) = 1.0/( roa(i,j,k,1) * MFAC )
            end forall
        else
            inv_air_density3D(:,:,:) = 1.0
        end if

        select case ( f_3d(n)%class )

         ! Simple advected species:
          case ( "ADV" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
            end forall

         case ( "BGN" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_bgn(index,i,j,k)
            end forall

         case ("XKSIG00", "XKSIG12" ) !hf hmix Kz_m2s

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = Kz_m2s(i,j,k)
            end forall

         case ("TH  " ) !JEJ Pot. temp (needed for cross sections)

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = th(i,j,k,1)
            end forall

!DSGC         case ( "PHNO3" )   !ds-hf  rv1_9_28
!DSGC            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
!DSGC              d_3d( n, i,j,k,IOU_INST) = xn_shl(index,i,j,k)
!DSGC            end forall
!DSGC
!DSGC             if(debug_flag) write(*,"(a12,i4,2es12.3)") "3D3D PHNO3", n, &
!DSGC                  xn_shl(index,debug_li,debug_lj,KMAX_MID), &
!DSGC                  d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST)

         case ( "MAX3DSHL" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )! Daily maxima - short-lived
              d_3d( n, i,j,k,IOU_INST) = max( d_3d( n, i,j,k,IOU_INST),&
                                      xn_shl(index,i,j,k) &
                                     * inv_air_density3D(i,j,k) )
            end forall

            if(debug_flag) write(*,"(a13,i4,f8.3,3es12.3)") "3D3D MAX3DSHL", n, thour, &
              xn_shl(index,debug_li,debug_lj,KMAX_MID), &
              1.0/inv_air_density3D(debug_li,debug_lj,KMAX_MID), &
              d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST)

          case ( "MAX3DADV" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) =  max( d_3d( n, i,j,k,IOU_INST),&
                                               xn_adv(index,i,j,k) )
            end forall

             if(debug_flag) write(*,"(a12,i4,f8.3,4es12.3)") "SET MAX3DADV", n, thour, &
                      xn_adv(index,debug_li,debug_lj,KMAX_MID), &
                      d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST)

          case ( "SHL" )
            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) =   xn_shl(index,i,j,k) * inv_air_density3D(i,j,k)
            end forall


          case ( "VOC" )

            call voc_3dcalc()

! hb new 3D output
        case ( "D3_PPB" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
            end forall
           if ( debug_flag ) call write_debugadv(n,index, 1.0, "D3 PPB OUTS")


! hb new 3D output
! ds Bug - cannot have PM in mixing ratio
!ds         case ( "PM25" )
!ds
!ds            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
!ds              d_3d( n, i,j,k,IOU_INST) = xn_adv(IXADV_SO4,i,j,k) &
!ds                + xn_adv(IXADV_pNO3_f,i,j,k) &
!ds                + xn_adv(IXADV_aNH4,i,j,k) &
!ds                + xn_adv(IXADV_PPM25,i,j,k) &
!ds                + xn_adv(IXADV_SeaSalt_f,i,j,k)
!ds            end forall
!ds
!ds! hb new 3D output
!ds         case ( "PMco" )
!ds
!ds            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
!ds              d_3d( n, i,j,k,IOU_INST) =   &
!ds                + xn_adv(IXADV_pNO3_c,i,j,k) &
!ds                + xn_adv(IXADV_PPMCOARSE,i,j,k) &
!ds                + xn_adv(IXADV_SeaSalt_c,i,j,k)
!ds            end forall

!AMVB 2010-07-21: PM-PPB bug fix
        case ( "PM25GROUP" )
          do k=1,KMAX_MID
            call uggroup_calc( d_3d(n,:,:,k,IOU_INST), n, typ, PM25_GROUP, &
                               roa(i,j,k,1), k)
          enddo
        case ( "PMcGROUP" )
          do k=1,KMAX_MID
            call uggroup_calc( d_3d(n,:,:,k,IOU_INST), n, typ, PMCO_GROUP, &
                               roa(i,j,k,1), k)
          enddo

! hb
         case ( "Kz" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = Kz_m2s(i,j,k)
            end forall


          case  default

           write(unit=errmsg,fmt=*) "Derived 3D class NOT FOUND", n, index, &
                         f_3d(n)%name,f_3d(n)%class
           call CheckStop( errmsg )


        end select


      !/** add to monthly and yearly average, and increment counters
       !    ( no daily averaging done for 3-D fields so far).


       ! For the MAX3D possibilities, we store maximum value of the
       !   current day in the IOU_INST variables.
       !   These are then added into IOU_MON **only** at the end of each day.
       ! (NB there is an error made on 1st day used, since only 1st 6 hours
       !  are checked. Still, not much happens on 1st Jan.... ;-)

        if ( (f_3d(n)%class == "MAX3DSHL")  .or. &
             (f_3d(n)%class == "MAX3DADV") )then
           if (End_of_Day) then
              d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
                                     + d_3d(n,:,:,:,IOU_INST)
              d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
                                     + d_3d(n,:,:,:,IOU_INST)
              if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1 !only collected for end_of_day

              if( debug_flag ) then
                    write(*,fmt="(a20,a9,i4,f8.3,2es12.3)") "END_OF_DAY MAX3D", &
                      f_3d(n)%class, n, thour,  &
                      d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_MON ),&
                      d_3d(n,debug_li,debug_lj,KMAX_MID,IOU_INST )
                    write(*,"(a20,i4,2x,6i6)") "END_OF_DAY NAV ", &
                      n, (nav_3d(n,i), i=1,LENOUT3D)
              end if

              d_3d(n,:,:,:,IOU_INST ) = 0.0  !! Reset d_3d

           endif ! End_of_Day
        else
           d_3d(n,:,:,:,IOU_DAY ) = d_3d(n,:,:,:,IOU_DAY ) &
                + d_3d(n,:,:,:,IOU_INST)
           d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
                + d_3d(n,:,:,:,IOU_INST)
           d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
                + d_3d(n,:,:,:,IOU_INST)
           if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1
        endif



!      !/** add to monthly and yearly average, and increment counters
!      !    ( no daily averaging done for 3-D fields so far).
!
!       d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
!                              + d_3d(n,:,:,:,IOU_INST)
!       d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
!                              + d_3d(n,:,:,:,IOU_INST)
!
!       if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1

      end do

      first_call = .false.

    end subroutine Derived
    !=========================================================================

    subroutine DerivedProds(text,dt)

    !/** DESCRIPTION
    !  Calculates chemical changes by comparing values before and  after
    !  chemistry subroutine. Intended to be a more flexible version of the old
    !  PRODO3  calculation

      character(len=*), intent(in) :: text  ! "Before" or "After"
      real,             intent(in) :: dt    ! timestep (s)

      real :: timefrac                      ! dt as fraction of hour (3600/dt)



!      if ( num_deriv3d < 1 ) print *, "DerivedProds "//text, num_deriv3d
      if ( num_deriv3d < 1 ) return
      if (.not. any( f_3d%class == "PROD" ) ) return

      timefrac = dt/3600.0

     !/***** 3-D fields **************************

     do n = 1, num_deriv3d

        if ( f_3d(n)%class  == "PROD " ) then
           index = f_3d(n)%index

           select case ( text )

               case ( "Before" )   !! Initialise to xn_adv

                 forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                   d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
                 end forall

               case ( "After" )    !! Calculate change

                 forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                   d_3d( n, i,j,k,IOU_INST) = &
                      d_3d( n, i,j,k,IOU_INST) - xn_adv(index,i,j,k)
                 end forall

           end select
        end if
      end do

    end subroutine DerivedProds
    !=========================================================================

    subroutine ResetDerived(period)
      integer, intent(in) :: period   ! Either IOU_DAY or IOU_MON

       if ( period <= LENOUT2D ) then
           nav_2d  (:,period) = 0.0
           d_2d(:,:,:,period) = 0.0
       end if


       if ( num_deriv3d > 0 .and.  period <= LENOUT3D ) then
           nav_3d    (:,period) = 0.0
           d_3d(:,:,:,:,period) = 0.0
       end if

    end subroutine ResetDerived
 !=========================================================================


   subroutine voc_2dcalc()

    !/-- Sums up voc species using the indices defined earlier in Setup_VOCs

     ! We initialise d_2d first, the use a simple loop
     ! over voc. Some CPU could be saved by initialising
     ! with the 1st voc, then looping over 2, nvoc, but who cares...


      d_2d( n, 1:limax,1:ljmax,IOU_INST) =  0.0

      do ivoc = 1, nvoc

         index = voc_index(ivoc)           ! Gives which IXADV_ to use.
         forall ( i=1:limax, j=1:ljmax )
             d_2d( n, i,j,IOU_INST) = d_2d( n, i,j,IOU_INST)      &
                                    + xn_adv(index,i,j,KMAX_MID)  &
                                    * voc_carbon(ivoc) * cfac(index,i,j)
                               ! multiplied by nr. of C and "reduced to surface"
         end forall
      end do ! ivoc
   end subroutine voc_2dcalc

 !=========================================================================
   subroutine voc_3dcalc()

    !/-- as for voc_2dcalc

      d_3d( n, 1:limax,1:ljmax,1:KMAX_MID,IOU_INST) =  0.0

      do ivoc = 1, nvoc

         index = voc_index(ivoc)
         forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
             d_3d( n, i,j,k,IOU_INST) = d_3d( n, i,j,k,IOU_INST) + &
                     xn_adv(index,i,j,k)*voc_carbon(ivoc)
         end forall
      end do ! ivoc

   end subroutine voc_3dcalc
 !=========================================================================
!AMVB 2010-08-16: PM-PPB bug fix
 subroutine uggroup_calc( ug_2d, n, class, group, density, ik)
!subroutine uggroup_calc( ug_2d, n, class, group, density )

  !/--  calulates e.g. SIA = SO4 + pNO3_f + pNO3_c + aNH4
  ! (only SIA converted to new group system so far, rv3_5_6 )
  !/--  calulates also PM10  = SIA + PPM2.5 + PPMCOARSE

  real, dimension(:,:), intent(inout) :: ug_2d  ! i,j section of d_2d arrays
  character(len=*)    :: class   ! Type of data
  integer, intent(in)  :: n   !
  !character(len=*)    :: unit   !
  integer, intent(in), dimension(:)  :: group
  real, intent(in), dimension(MAXLIMAX,MAXLJMAX)  :: density
  integer, intent(in) :: ik   !AMVB 2010-08-16: PM-PPB bug fix
  integer :: ig, iadv, itot,k
  real :: scale
  character(len=10) :: unit=""

  if(DEBUG .and. debug_proc) then
    write(*,"(a,4i4,2es12.3)") "DEBUG GROUP-PM-N", size(group)
  end if

!AMVB 2010-08-16: PM-PPB bug fix
  if     (ik==0 .and. n<=num_deriv2d) then
    k=KMAX_MID
    unit=f_2d(n)%unit
  elseif (ik/=0 .and. n<=num_deriv3d) then
    k=ik
    unit=f_3d(n)%unit
  else
    k=-1
    unit="not_found"
  endif

  ug_2d( :,:) = 0.0
  do ig = 1, size(group)
    itot = group(ig)
    iadv = group(ig) - NSPEC_SHL

    select case (trim(unit))
      case("ug/m3" ); scale = species(itot)%molwt
      case("ugN/m3"); scale = species(itot)%nitrogens!*atwN !AMVB 2010-08-16: PM-PPB bug fix
      case default  ; call StopAll("uggroup called with wrong unit='"//unit//"'!")
    end select

    if(ik==0)then
      forall ( i=1:limax, j=1:ljmax )
        ug_2d( i,j) = ug_2d( i,j) + xn_adv(iadv,i,j,k) *scale * cfac(iadv,i,j)
      end forall
    else
      forall ( i=1:limax, j=1:ljmax )
        ug_2d( i,j) = ug_2d( i,j) + xn_adv(iadv,i,j,k) *scale
      end forall
    endif

   if(DEBUG .and. debug_proc) then
      i=debug_li
      j=debug_lj
      write(*,"(a,4i4,f6.1,2es12.3)") "DEBUG GROUP-PM", ig, &
        itot, iadv, species(itot)%molwt, scale, xn_adv(iadv,i,j,k), ug_2d(i,j)
    end if
  end do !n
  forall ( i=1:limax, j=1:ljmax )
    ug_2d( i,j) = ug_2d( i,j) * density(i,j)
  end forall

 end subroutine uggroup_calc
 !=========================================================================

  subroutine somo_calc( n , debug_flag )


    !/-- Calculates SOMO (8hours) values for input threshold.

    implicit none
    integer, intent(in) :: n           ! index in Derived_ml::d_2d arrays
    logical, intent(in) :: debug_flag

    real    :: threshold               ! Threshold, e.g. 35 (ppb)
    real :: o3                         ! Ozone (ppb) - needed if SOMOs
    real :: sum8h
    integer, parameter :: N8h = (NTDAY*8)/24 !number of periods in 8 hours
    real, parameter :: N8h_inv=1./N8h
    integer :: nh


    !BUG threshold = f_2d(n)%index
    threshold = f_2d(n)%Threshold

      do i=1,limax
        do j=1,ljmax

           !find max running 8h sum O3
           sum8h=0.
           do nh=1,N8h
              sum8h = sum8h + D2_O3_DAY( i , j , nh)
           enddo
           o3=sum8h
           do nh=N8h+1,NTDAY
              sum8h =sum8h-D2_O3_DAY( i , j , nh-N8h)+D2_O3_DAY( i , j , nh)
              o3=max(o3,sum8h)
              if(n<0)write(*,*)o3 !pw fake for compiler!!
           enddo

           !divide by N8h to find 8h mean
           o3=o3*N8h_inv

           if ( debug_flag .and. i==debug_li .and. j==debug_lj ) then
             write(*,"(a,i4,f7.1,f12.3)") "SOMO DEBUG ", n, threshold, o3
           end if
   

           o3 = max( o3 - threshold , 0.0 )   ! Definition of SOMOs

             ! d_2d values will be accumulated in Derived_ml

           d_2d(n, i,j,IOU_DAY ) = o3

        end do
      end do
   end subroutine somo_calc

 !=========================================================================
    subroutine write_debugadv(n,index,rho,txt)
       integer, intent(in) :: n, index
       real, intent(in) :: rho
       character(len=*) :: txt

       write(*,fmt="(2a,2i4,a,4f12.3)") "PROCESS " , txt , n, index  &
                  ,trim(f_2d(n)%name)  &
                  ,d_2d(n,debug_li,debug_lj,IOU_INST)*PPBINV &
                  ,xn_adv(index,debug_li,debug_lj,KMAX_MID)*PPBINV &
                  ,rho, cfac(index,debug_li,debug_lj)
    end subroutine write_debugadv
 !=========================================================================
    subroutine write_debug(n,index,txt)
       integer, intent(in) :: n, index
       character(len=*) :: txt

       write(*,fmt="(2a,2i4,a,4f12.3)") "DERIV: GEN " , txt , n, index  &
                  ,trim(f_2d(n)%name)  &
                  ,d_2d(n,debug_li,debug_lj,IOU_INST)
    end subroutine write_debug

end module Derived_ml
