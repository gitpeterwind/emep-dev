! <My_DryDep_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module My_DryDep_ml    ! DryDep_ml
!+
! Module to define the dry deposition components and rates. We
! define the min (vd_min) and max dep. velocity (Vg) first and then derive the
! daytime addition (vd_day).
!/**************************************************************************
!  Specifies which of the possible species (from Wesely's list)
!  are required in the current air pollution model   
!/**************************************************************************

 use GenChemicals_ml, only : species
 use CheckStop_ml,  only : CheckStop, StopAll
 use Derived_ml,    only : f_2d,   d_2d, IOU_INST
 use My_Derived_ml,  only : &
   nOutDDep, DDEP_LCS , OutDDep, nOutVg, OutVg, nOutRG, OutRG &
  ,SOX_INDEX, OXN_INDEX, RDN_INDEX &  ! Equal -1, -2, -3
  ,DDEP_SOXGROUP, DDEP_OXNGROUP, DDEP_RDNGROUP, DDEP_GROUP

 use GenChemicals_ml,    only : species 
 use GenSpec_adv_ml        !   e.g. NSPEC_ADV,IXADV_O3,IXADV_H2O2,
 use GenSpec_shl_ml,     only : NSPEC_SHL   ! For DDEP_SOXGROUP etc.
 use LandDefs_ml,        only : LandDefs, LandType
 use Landuse_ml,         only : WheatGrowingSeason
 use LocalVariables_ml,  only : Grid  !=> izen  integer of zenith angle
 use ModelConstants_ml , only : atwS, atwN, AOT_HORIZON
use Par_ml, only: me ! Vds Ndep
 use PhysicalConstants_ml, only : AVOG
 use SmallUtils_ml,      only: find_index, NOT_FOUND
 use StoFlux_ml,         only : unit_flux, lai_flux, leaf_flux
 use TimeDate_ml,        only : current_date
 use Wesely_ml
 implicit none
 private

  public :: Init_DepMap
  public :: Add_ddep
  public :: Add_Vg
  public :: Add_RG

  !/** Variables used in deposition calculations
 
  ! DDEP_xx gives the index that will be used in the EMEP model
  ! WES_xx gives the index of the Wesely gas to which this corresponds

  integer, private, save :: &
!    DDEP_SOX,   DDEP_OXN,   DDEP_RDN,  &
!    D2_VddACC, D2_VddCOA, & !ECO08 addition
    D2_AFSTDF0, D2_AFSTDF16, D2_AFSTBF0, D2_AFSTBF16, & 
    D2_AFSTCR0, D2_AFSTCR3, D2_AFSTCR6,&
    iam_medoak, iam_beech, iam_wheat, &   ! For Fluxes
    D2_O3DF,    D2_O3WH, &
    D2_EUAOT30WH, D2_EUAOT40WH, D2_EUAOT30DF, D2_EUAOT40DF, &
    D2_UNAOT30WH, D2_UNAOT40WH, D2_UNAOT30DF, D2_UNAOT40DF, &
    D2_MMAOT40WH, D2_MMAOT30WH


  ! Here we define the minimum set of species which has different
  ! deposition velocities. We calculate Vg for these, and then
  ! can use the rates for other similar species. (e.g. AMSU can use
  ! the Vg for SO4.  Must set NDRYDEP_CALC species

  !/** IMPORTANT: the variables below must match up in the sense that, for 
  ! example, if DDEP_NH3=4 then the 4th element of DRYDEP must be WES_NH3.

  integer, public, parameter :: NDRYDEP_GASES = 10  ! gases
  integer, public, parameter :: NDRYDEP_AER = 2    ! aerosols
  integer, public, parameter :: NDRYDEP_CALC = NDRYDEP_GASES + NDRYDEP_AER


  integer, public, parameter :: &
       CDEP_HNO3 = 1, CDEP_O3  = 2, CDEP_SO2 = 3  &
      ,CDEP_NH3  = 4, CDEP_NO2 = 5, CDEP_PAN  = 6 &
      ,CDEP_H2O2 = 7, CDEP_ALD = 8, CDEP_HCHO = 9, &
       CDEP_OP = 10,  CDEP_FIN = 11, CDEP_COA = 12 

  integer, public, parameter :: CDEP_SET = -99    



 ! WE NEED A FLUX_CDEP, FLUX_ADV FOR OZONE;
 ! (set to one for non-ozone models)

  logical, public, parameter :: STO_FLUXES = .true.
  integer, public, parameter :: FLUX_CDEP  = CDEP_O3
  integer, public, parameter :: FLUX_ADV   = IXADV_O3

 
  integer, public, parameter, dimension(NDRYDEP_GASES) :: &
    DRYDEP_GASES = (/ WES_HNO3, WES_O3,   WES_SO2, &
                     WES_NH3,  WES_NO2 , WES_PAN, &
                     WES_H2O2, WES_ALD, WES_HCHO, WES_OP    /)

  !/** Compensation pount approach from CEH used?:

  logical, public, parameter :: COMPENSATION_PT = .false. 



  ! We define also the number of species which will be deposited in
  ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_GASES
  ! The actual species used and their relation to the CDEP_ indices
  ! above will be defined in Init_DepMap

  integer, public, parameter ::  NDRYDEP_ADV  = 22

  !/-- we define a type to map indices of species to be deposited
  !   to the lesser number of species where Vg is calculated

   type, public :: depmap
      integer :: adv   ! Index of species in IXADV_ arrays
      integer :: calc  ! Index of species in  calculated dep arrays
      real    :: vg    ! if CDEP_SET, give vg in m/s
   end type depmap

   type(depmap), public, dimension(NDRYDEP_ADV):: Dep

   real, public, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost

  ! Maps from adv index to one of calc indices
   real, public, save, dimension(NSPEC_ADV) :: DepAdv2Calc ! ECO08

  ! depositions are calculated to the following landuse classes, where
  ! e.g. conif may include both temperate and Medit.
   character(len=7), private, dimension(0:4), parameter :: &
    DEP_RECEIVERS = (/ "Grid", "Conif", "Decid", "Crops", "Seminat" /)
   integer, private, parameter :: &
    GRID_LC=0, CONIF=1, DECID=2, CROP=3, SEMINAT=4

   logical, private, parameter :: MY_DEBUG = .false.
   logical, private, parameter :: DEBUG_ECO = .false.
   logical, private, save :: first_call = .true.

contains
  subroutine Init_DepMap
   integer, dimension(22) :: check_vals  ! Tmp safety check array
   integer :: icheck, iadv, i, i2, n, ndep, nVg, nRG

 ! .... Define the mapping between the advected species and
 !      the specied for which the calculation needs to be done.

   Dep(1) =  depmap( IXADV_HNO3 , CDEP_HNO3, -1.)
   Dep(2) =  depmap( IXADV_PAN,   CDEP_PAN, -1. ) 
   Dep(3) =  depmap( IXADV_NO2,   CDEP_NO2, -1. )
   Dep(4) =  depmap( IXADV_SO2,   CDEP_SO2, -1. )
   Dep(5) =  depmap( IXADV_SO4,   CDEP_FIN,  -1)
   Dep(6) =  depmap( IXADV_NH3,   CDEP_NH3, -1. )
   Dep(7) =  depmap( IXADV_aNH4,  CDEP_FIN,  -1) 
   Dep(8) =  depmap( IXADV_aNO3,  CDEP_FIN,  -1) 
   Dep(9) =  depmap( IXADV_O3   , CDEP_O3  , -1.)
   Dep(10) =  depmap( IXADV_H2O2 , CDEP_H2O2, -1.)
   Dep(11) =  depmap( IXADV_MPAN , CDEP_PAN , -1.)
   Dep(12) =  depmap( IXADV_HCHO , CDEP_HCHO, -1.)
   Dep(13) =  depmap( IXADV_CH3CHO,CDEP_ALD , -1.)
   Dep(14) =  depmap( IXADV_MAL   ,CDEP_ALD , -1.)
   Dep(15) =  depmap( IXADV_CH3O2H,CDEP_OP  , -1.)
   Dep(16) =  depmap( IXADV_C2H5OOH,CDEP_OP  , -1.)
   Dep(17) =  depmap( IXADV_pNO3,  CDEP_COA, -1.)
   Dep(18) =  depmap( IXADV_PM25,  CDEP_FIN, -1. )
   Dep(19) =  depmap( IXADV_PMco,  CDEP_COA, -1. )
   Dep(20) =  depmap( IXADV_SSfi,  CDEP_FIN, -1. )
   Dep(21) =  depmap( IXADV_SSco,  CDEP_COA, -1. )
   Dep(22) =  depmap( IXADV_Pb210,  CDEP_FIN, -1. )


!#######################  ECO08 new mappng  #######################
  do i = 1, NDRYDEP_ADV  ! 22
      iadv = Dep(i)%adv
      DepAdv2Calc(iadv) = Dep(i)%calc
 end do
  
! We process the various combinations of gas-species and ecosystem:
! starting with DryDep, e.g. DDEP_SO2_m2CF

  do ndep = 1, nOutDDep

    OutDdep(ndep)%f2d  = find_index(OutDdep(ndep)%name ,f_2d(:)%name)
   ! We start LC with Grid at zero, so need -1 offset
    OutDdep(ndep)%LC   = find_index(OutDdep(ndep)%txt ,DEP_RECEIVERS) -1

    if(MY_DEBUG .and. me==0) write(6,*) "OUTDdep ", ndep, &
       OutDdep(ndep)%name,OutDdep(ndep)%txt,"=>"&
          ,OutDdep(ndep)%LC, OutDdep(ndep)%f2d
    call CheckStop( OutDDep(ndep)%f2d < 1, &
        "OutDDep-f2d Error " // OutDDep(ndep)%name)
    call CheckStop( OutDDep(ndep)%LC < 0, & !0=GRID
        "OutDDep-LC Error " // OutDDep(ndep)%name)
  end do
     
  do nVg = 1, nOutVg

     if( OutVg(nVg)%txt == "Grid") then
        OutVg(nVg)%LC = GRID_LC   ! zero
     else 
        OutVg(nVg)%LC = find_index( OutVg(nVg)%txt, LandDefs(:)%code )
     end if
     OutVg(nVg)%f2d = find_index( OutVg(nVg)%name, f_2d(:)%name )

     if(MY_DEBUG .and. me==0) write(6,*) "OUTVG ", nVg, &
         OutVg(nVg)%name,OutVg(nVg)%txt,"=>",&
           OutVg(nVg)%LC, OutVg(nVg)%f2d

     call CheckStop( OutVg(nVg)%LC < GRID_LC , & !ie zero
          "OutVg-LC Error " // OutVg(nVg)%name)
     call CheckStop( OutVg(nVg)%f2d < 1, &
          "OutVg-f2d Error " // OutVg(nVg)%name)
  end do

 do nRG = 1, nOutRG

     if( OutRG(nRG)%txt == "Grid") then
        OutRG(nRG)%LC = GRID_LC    ! zero
     else 
        OutRG(nRG)%LC = find_index( OutRG(nRG)%txt, LandDefs(:)%code )
     end if
     OutRG(nRG)%f2d = find_index( OutRG(nRG)%name, f_2d(:)%name )

     if(MY_DEBUG .and. me==0) write(6,*) "OUTRG ", nRG, &
         OutRG(nRG)%name,OutRG(nRG)%txt,"=>",&
           OutRG(nRG)%LC, OutRG(nRG)%f2d

     call CheckStop( OutRG(nRG)%LC < GRID_LC , & ! <zero
          "OutRG-LC Error " // OutRG(nRG)%name)
     call CheckStop( OutRG(nRG)%f2d < 1, &
          "OutRG-f2d Error " // OutRG(nRG)%name)
  end do



!#######################  NEW define indices here #######################

D2_AFSTDF0 = find_index("D2_AFSTDF0",f_2d(:)%name)
D2_AFSTDF16 = find_index("D2_AFSTDF16",f_2d(:)%name)

D2_AFSTBF0 = find_index("D2_AFSTBF0",f_2d(:)%name)
D2_AFSTBF16 = find_index("D2_AFSTBF16",f_2d(:)%name)

D2_AFSTCR0 = find_index("D2_AFSTCR0",f_2d(:)%name)
D2_AFSTCR3 = find_index("D2_AFSTCR3",f_2d(:)%name)
D2_AFSTCR6 = find_index("D2_AFSTCR6",f_2d(:)%name)

D2_O3DF    = find_index("D2_O3DF   ",f_2d(:)%name)
D2_O3WH    = find_index("D2_O3WH   ",f_2d(:)%name)

D2_EUAOT30WH    = find_index("D2_EUAOT30WH",f_2d(:)%name)
D2_EUAOT40WH    = find_index("D2_EUAOT40WH",f_2d(:)%name)
D2_EUAOT30DF    = find_index("D2_EUAOT30DF",f_2d(:)%name)
D2_EUAOT40DF    = find_index("D2_EUAOT40DF",f_2d(:)%name)

D2_UNAOT30WH    = find_index("D2_UNAOT30WH",f_2d(:)%name)
D2_UNAOT40WH    = find_index("D2_UNAOT40WH",f_2d(:)%name)
D2_UNAOT30DF    = find_index("D2_UNAOT30DF",f_2d(:)%name)
D2_UNAOT40DF    = find_index("D2_UNAOT40DF",f_2d(:)%name)

D2_MMAOT30WH    = find_index("D2_MMAOT30WH",f_2d(:)%name)
D2_MMAOT40WH    = find_index("D2_MMAOT40WH",f_2d(:)%name)

iam_wheat   = find_index("IAM_CR",LandDefs(:)%code)
iam_beech   = find_index("IAM_DF",LandDefs(:)%code)
iam_medoak  = find_index("IAM_MF",LandDefs(:)%code)
!####################### ds END of define indices #######################
  check_vals = (/  &
    D2_AFSTDF0, D2_AFSTDF16, D2_AFSTBF0, D2_AFSTBF16, & !4
    D2_AFSTCR0, D2_AFSTCR3, D2_AFSTCR6,& !3
    iam_medoak, iam_beech, iam_wheat, &   ! For Fluxes !3
    D2_O3DF,    D2_O3WH, & !2
    D2_EUAOT30WH, D2_EUAOT40WH, D2_EUAOT30DF, D2_EUAOT40DF, & !4
    D2_UNAOT30WH, D2_UNAOT40WH, D2_UNAOT30DF, D2_UNAOT40DF, & !4
    D2_MMAOT40WH, D2_MMAOT30WH & !2
  /) ! => 37 - 12 -> 25
    ! Will re-write whole subroutine later to avoid individual indices, but  
    ! for now, do:
     do icheck = 1, 22
        call CheckStop( check_vals(icheck) < 1 , "D2D CHECKVAL! " )
     end do
     call CheckStop( any(check_vals < 1) , "D2D CHECKVAL! " )
 
     if(me==0) write(6,*) "Init_DepMap D2D FINISHED"
     call CheckStop( icheck < 0 , "D2D STOPPER! " )

  end subroutine Init_DepMap

  !<==========================================================================
  subroutine Add_ddep(debug_flag,dt,i,j,convfac,lossfrac,fluxfrac,&
           c_hvegppb,coverage)

  !<==========================================================================
     ! Adds deposition losses to ddep arrays
     logical, intent(in) :: debug_flag
     real,    intent(in) :: dt              ! time-step
     integer, intent(in) :: i,j             ! coordinates
     real,    intent(in) ::  convfac, lossfrac
     real, dimension(:,:), intent(in) ::  fluxfrac   ! dim (NADV, NLANDUSE)
     real, dimension(:), intent(in) ::  c_hvegppb   ! dim (NLANDUSE)
     real, dimension(:), intent(in) ::  coverage    ! dim (NLANDUSE), =fraction
     integer :: n, nadv, nadv2, ihh, idd, imm
     real :: o3WH, o3DF   ! O3 over wheat, decid forest


     real, parameter  :: NMOLE_M3 = 1.0e6*1.0e9/AVOG  ! Converts from 
                                                      ! mol/cm3 to nmole/m3

  !ECO08:  Variables added for new ecosystem dep
     real, dimension(0:size(DEP_RECEIVERS)-1) :: invEcoArea, EcoArea
     integer :: ndep, RLC   ! RLC = receiver land-cover
     real :: Fflux

     real ::  to_nmole, timefrac, fstfrac 

     to_nmole =  NMOLE_M3
     timefrac = dt/3600.0
     fstfrac  = dt*1.0e-6     ! Converts also nmole to mmole


  ! Must match areas given above, e.g. DDEP_CONIF -> Conif

  !ECO08, Ecosystem areas:
     EcoArea(CONIF)   = sum( coverage(:), LandType(:)%is_conif )
     EcoArea(DECID)   = sum( coverage(:), LandType(:)%is_decid )
     EcoArea(CROP)    = sum( coverage(:), LandType(:)%is_crop  )
     EcoArea(SEMINAT) = sum( coverage(:), LandType(:)%is_seminat )
     EcoArea(GRID_LC)    = 1.0

     invEcoArea(:) = 0.0

     do n = 0, size(DEP_RECEIVERS)-1
        if ( EcoArea(n) > 1.0e-39 ) invEcoArea(n) = 1.0/EcoArea(n)
     end do 


   !  Query - crops, outisde g.s. ????
     !if ( DEBUG_ECO .and. first_call .and. debug_flag ) then
     if ( DEBUG_ECO .and. debug_flag ) then
       write(*,*)  "ECOAREAS ", i,j, EcoArea(:)
         do n = 1, 19
           write(*,*)  "ECOCHECK ", n, i,j, LandType(n)%is_conif, coverage(n)
         end do
       write(*,*) "ECOCHECK ========================"
        
     end if
     first_call = .false.


    ! Traditional depositions to whole grid::

!     d_2d(DDEP_SOX,i,j,IOU_INST) = (  &
!          DepLoss(IXADV_SO2) + DepLoss(IXADV_SO4) ) * convfac * atwS
!
!       d_2d(D2_VddCOA,i,j,IOU_INST) =  100.0* Vg3m(CDEP_COA)

    ! Ecosystem depositions, for grouped or individual species:

     do ndep = 1, nOutDDep
        nadv  = OutDDep(ndep)%Adv
        RLC   = OutDDep(ndep)%LC
        if ( nadv > 0 ) then  ! normal advectde species
           nadv2 = 1
           DDEP_GROUP(1) = nadv
        else if ( nadv == SOX_INDEX ) then
           nadv2 = size( DDEP_SOXGROUP )
           DDEP_GROUP(1:nadv2) = DDEP_SOXGROUP - NSPEC_SHL
        else if ( nadv == OXN_INDEX ) then
           nadv2 = size( DDEP_OXNGROUP )
           DDEP_GROUP(1:nadv2) = DDEP_OXNGROUP - NSPEC_SHL
        else  !if ( nadv == RDN_INDEX ) then
           nadv2 = size( DDEP_RDNGROUP )
           DDEP_GROUP(1:nadv2) = DDEP_RDNGROUP - NSPEC_SHL
        end if
   
      Fflux = 0.0
      do n = 1, nadv2
         nadv = DDEP_GROUP(n)
        if ( RLC == GRID_LC ) then 
            Fflux = Fflux + Deploss(nadv) * sum( fluxfrac(nadv,:) )
        else if ( RLC == CONIF ) then 
            Fflux = Fflux + Deploss(nadv) * &
               sum( fluxfrac(nadv,:), LandType(:)%is_conif )
        else if ( RLC == DECID ) then 
            Fflux = Fflux + Deploss(nadv) * &
               sum( fluxfrac(nadv,:), LandType(:)%is_decid )
        else if ( RLC == CROP ) then 
            Fflux = Fflux + Deploss(nadv) * &
               sum( fluxfrac(nadv,:), LandType(:)%is_crop )
        else if ( RLC == SEMINAT ) then 
            Fflux = Fflux + Deploss(nadv) * &
               sum( fluxfrac(nadv,:), LandType(:)%is_seminat )
        else 
            Fflux = -1.0
        end if
      end do ! n
        if ( OutDDep(ndep)%f2d < 1 .or. Fflux < 0.0 ) then
             write(6,*) "CATASTR ", ndep, OutDDep(ndep)%f2d,OutDDep(ndep)%name
             write(6,*) "CATASTR  RLC", RLC
             call CheckStop("CATASTROPHE: "//OutDDep(ndep)%name)
        end if

  !ECO08 
  ! - invEcoAreaCF divides the flux per grid by the landarea of each
  ! ecosystem, to give deposition in units of mg/m2 of ecosystem.

        d_2d( OutDDep(ndep)%f2d,i,j,IOU_INST) =  &
            Fflux * convfac * OutDDep(ndep)%atw * invEcoArea(RLC)

        if ( DEBUG_ECO .and. debug_flag ) then
           write(6,*) "DEBUG_ECO Fflux ", ndep, nadv, RLC, Fflux, &
            d_2d( OutDDep(ndep)%f2d,i,j,IOU_INST), DDEP_SOXGROUP
        end if ! DEBUG_ECO 
     end do ! ndep





!MAPPING_MANUAL CHANGES:
!  Use 1.6 for Beech and 3 for crops

!Beech:
     d_2d(D2_AFSTDF0,i,j,IOU_INST) =  fstfrac*leaf_flux(iam_beech)
     d_2d(D2_AFSTDF16,i,j,IOU_INST) = fstfrac* max(leaf_flux(iam_beech)-1.6,0.0)
!Med. Oak:
     d_2d(D2_AFSTBF0,i,j,IOU_INST) =  fstfrac*leaf_flux(iam_medoak)
     d_2d(D2_AFSTBF16,i,j,IOU_INST) = fstfrac* max(leaf_flux(iam_medoak)-1.6,0.0)
!Crops
     d_2d(D2_AFSTCR0,i,j,IOU_INST) =  fstfrac*leaf_flux(iam_wheat)
     d_2d(D2_AFSTCR3,i,j,IOU_INST) =  fstfrac*max(leaf_flux(iam_wheat)-3.0,0.0)
     d_2d(D2_AFSTCR6,i,j,IOU_INST) =  fstfrac*max(leaf_flux(iam_wheat)-6.0,0.0)

   !--- ecosystem specific concentrations..
   ! - use Conif forest for forests - safer for growing seasons

     imm      =    current_date%month            ! for debugging
     idd      =    current_date%day              ! for debugging
     ihh      =    current_date%hour             ! for debugging


     o3WH = c_hvegppb(iam_wheat)* lossfrac
     o3DF = c_hvegppb(iam_beech)* lossfrac

     d_2d(D2_O3DF,i,j,IOU_INST) =   o3DF
     d_2d(D2_O3WH,i,j,IOU_INST) =   o3WH

     if ( ihh >= 9 .and. ihh <= 21 ) then ! 8-20 CET, assuming summertime

        d_2d(D2_EUAOT30WH,i,j,IOU_INST) =  max(o3WH-30.0,0.0) * timefrac
        d_2d(D2_EUAOT40WH,i,j,IOU_INST) =  max(o3WH-40.0,0.0) * timefrac
        d_2d(D2_EUAOT30DF,i,j,IOU_INST) =  max(o3DF-30.0,0.0) * timefrac
        d_2d(D2_EUAOT40DF,i,j,IOU_INST) =  max(o3DF-40.0,0.0) * timefrac
     else
        d_2d(D2_EUAOT30WH,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT40WH,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT30DF,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT40DF,i,j,IOU_INST) =  0.0
     end if


    !/-- Calcuates AOT values for specific veg. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )


           if ( Grid%izen < AOT_HORIZON ) then

             d_2d(D2_UNAOT30WH,i,j,IOU_INST) =  max(o3WH-30.0,0.0) * timefrac
             d_2d(D2_UNAOT40WH,i,j,IOU_INST) =  max(o3WH-40.0,0.0) * timefrac
             d_2d(D2_UNAOT30DF,i,j,IOU_INST) =  max(o3DF-30.0,0.0) * timefrac
             d_2d(D2_UNAOT40DF,i,j,IOU_INST) =  max(o3DF-40.0,0.0) * timefrac

           else
             d_2d(D2_UNAOT30WH,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT40WH,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT30DF,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT40DF,i,j,IOU_INST) =  0.0
           end if

       ! MM AOT added (same as UNECE, but different growing season)
             d_2d(D2_MMAOT30WH,i,j,IOU_INST) = d_2d(D2_UNAOT30WH,i,j,IOU_INST)&
                     * WheatGrowingSeason(i,j)
             d_2d(D2_MMAOT40WH,i,j,IOU_INST) = d_2d(D2_UNAOT40WH,i,j,IOU_INST)&
                     * WheatGrowingSeason(i,j)

    if ( DEBUG_ECO .and. debug_flag ) then
          write(6,"(a12,3i5,f7.2,5es12.3,i3,es12.3)") "DEBUG_ECO ", &
          imm, idd, ihh, o3WH, &
             leaf_flux(iam_beech), d_2d(D2_AFSTDF0,i,j,IOU_INST), &
             leaf_flux(iam_wheat), d_2d(D2_AFSTCR0,i,j,IOU_INST), &
             d_2d(D2_UNAOT40WH,i,j,IOU_INST), &
             WheatGrowingSeason(i,j), d_2d(D2_MMAOT40WH,i,j,IOU_INST)
    end if ! DEBUG

   !---- end ecosystem specific ----------------------------------------------

  end subroutine  Add_ddep
  !<==========================================================================
  subroutine Add_Vg(debug_flag,i,j,VgLU)
     logical, intent(in) :: debug_flag
     integer, intent(in) :: i,j             ! coordinates
     real, dimension(:,0:), intent(in) :: VgLU  ! dim (NLANDUSE*nlu), 0=Grid
     integer :: n, nVg, cdep

       do nVg = 1, nOutVg
         cdep  =  DepAdv2Calc( OutVg(nVg)%Adv ) ! Convert e.g. IXADV_O3
                                                  ! to CDEP_O3
         d_2d( OutVg(nVg)%f2d,i,j,IOU_INST) = VgLU( cdep, OutVg(nVg)%LC ) 

         if( MY_DEBUG .and. debug_flag ) then
              write(*,"(a,a,i3,i4,f8.3)") "ADD_VG: ",  OutVg(nVg)%name, cdep,  &
                 OutVg(nVg)%LC,  100.0*d_2d( OutVg(nVg)%f2d,i,j,IOU_INST)
         end if
       end do

  end subroutine Add_Vg
  !<==========================================================================
  subroutine Add_RG(debug_flag,i,j, GsLU, GnsLU)
     logical, intent(in) :: debug_flag
     integer, intent(in) :: i,j             ! coordinates
     real, dimension(:,0:), intent(in) :: GsLU, GnsLU  ! dim (0:NLANDUSE*nlu)
     integer :: n, nRG, cdep
     real :: Gs, Gns, Rs, Rns

       do nRG = 1, nOutRG
         cdep  =  DepAdv2Calc( OutRG(nRG)%Adv ) ! Convert e.g. IXADV_O3
                                                  ! to CDEP_O3
         Gs = GsLU( cdep, OutRG(nRG)%LC )
         if( Gs < 1.0e-44 ) then
            Rs = -999.0
         else
            Rs = 1.0/Gs
         end if

         Gns = GnsLU( cdep, OutRG(nRG)%LC )
         if( Gns < 1.0e-44 ) then
            Rns = -999.0
         else
            Rns = 1.0/Gns
         end if

         if( trim( OutRG(nRG)%label ) == "Rs" )  then
           d_2d( OutRG(nRG)%f2d,i,j,IOU_INST) = Rs
         else if( trim( OutRG(nRG)%label ) == "Rns" )  then
           d_2d( OutRG(nRG)%f2d,i,j,IOU_INST) = Rns
         else if( trim( OutRG(nRG)%label ) == "Gns" )  then
           d_2d( OutRG(nRG)%f2d,i,j,IOU_INST) = Gns
         else if( trim( OutRG(nRG)%label ) == "Gs" )  then
           d_2d( OutRG(nRG)%f2d,i,j,IOU_INST) = Gs
         end if

         if( MY_DEBUG .and. debug_flag ) then
              write(*,"(a,a,i3,2i4,es12.3)") "ADD_RG: ", OutRG(nRG)%name,&
                cdep,  OutRG(nRG)%LC,  OutRG(nRG)%f2d, &
                   d_2d( OutRG(nRG)%f2d,i,j,IOU_INST)
         end if
       end do

  end subroutine Add_RG
  !<==========================================================================


end module My_DryDep_ml

