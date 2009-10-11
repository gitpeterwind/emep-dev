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

 use My_Derived_ml,  only : &
   nOutDDep, OutDDep, nOutVg, OutVg, nOutRG, OutRG &
  ,nOutMET, OutMET, nOutVEGO3, OutVEGO3 & !MAR2009
  ,SOX_INDEX, OXN_INDEX, RDN_INDEX &  ! Equal -1, -2, -3
  ,DDEP_SOXGROUP, DDEP_RDNGROUP, DDEP_GROUP

 use CheckStop_ml,  only : CheckStop, StopAll
 use ChemChemicals_ml,    only : species 
 use ChemSpecs_adv_ml        !   e.g. NSPEC_ADV,IXADV_O3,IXADV_H2O2,
 use ChemSpecs_shl_ml,     only : NSPEC_SHL   ! For DDEP_SOXGROUP etc.
 use ChemGroups_ml,       only : DDEP_OXNGROUP  !DSGC
 use Derived_ml,    only : f_2d, d_2d, IOU_INST, AddDef,IOU_YEAR
 use EcoSystem_ml,      only : DEF_ECOSYSTEMS, NDEF_ECOSYSTEMS, &
                               EcoSystemFrac, & 
                                 FULL_GRID, Is_EcoSystem !DSOCT09
 use GridValues_ml,  only: debug_proc, debug_li, debug_lj
 use LandDefs_ml,        only : LandDefs, LandType
 use Landuse_ml,         only : WheatGrowingSeason, LandCover   ! DSOCT09
 use LocalVariables_ml,  only : Grid  !=> izen  integer of zenith angle
 use ModelConstants_ml , only : atwS, atwN, AOT_HORIZON, DEBUG_MY_DRYDEP,MasterProc
 use Par_ml,             only: li0, lj0, li1, lj1
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
  public :: Add_LCC_Met

  !/** Variables used in deposition calculations
 
  ! DDEP_xx gives the index that will be used in the EMEP model
  ! WES_xx gives the index of the Wesely gas to which this corresponds

!  integer, private, save :: &
!    DDEP_SOX,   DDEP_OXN,   DDEP_RDN,  &
!    D2_VddACC, D2_VddCOA, & !ECO08 addition
!    D2_AFSTDF0, D2_AFSTDF16, D2_AFSTBF0, D2_AFSTBF16, & 
!    D2_AFSTCR0, D2_AFSTCR3, D2_AFSTCR6,&
!    iam_medoak, iam_beech, iam_wheat, &   ! For Fluxes
!Oct09    D2_O3DF,    D2_O3WH, &
!Oct09    D2_EUAOT40WH, D2_EUAOT40DF, &
!Oct09    D2_UNAOT40WH, D2_UNAOT40DF, &
!Oct09    D2_MMAOT40WH


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


  !/-- we define a type to map indices of species to be deposited
  !   to the lesser number of species where Vg is calculated

   type, public :: depmap
      integer :: adv   ! Index of species in IXADV_ arrays
      integer :: calc  ! Index of species in  calculated dep arrays
      real    :: vg    ! if CDEP_SET, give vg in m/s
   end type depmap

   real, public, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost

  ! Maps from adv index to one of calc indices
   integer, public, save, dimension(NSPEC_ADV) :: DepAdv2Calc ! ECO08

   logical, private, save :: first_call = .true.

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   ! .... Define the mapping between the advected species and
   !      the specied for which the calculation needs to be done.
   !  We also define the number of species which will be deposited in
   ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_GASES
   ! The actual species used and their relation to the CDEP_ indices
   ! above will be defined in Init_DepMap

       include 'CM_DryDep.inc'

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

contains
  subroutine Init_DepMap
   integer :: icheck, iadv, i, i2, n, ndep, nVg, nRG, nMET

!#######################  ECO08 new mappng  #######################
  do i = 1, NDRYDEP_ADV  ! 22
      iadv = Dep(i)%adv
      if(DEBUG_MY_DRYDEP .and. MasterProc) &
         write(6,*) "DEPMAP   ", Dep(i)%adv, Dep(i)%calc
       if ( iadv < 1 ) then
           print *, "DEPPROB  ",  i, Dep(i)%adv, Dep(i)%calc
           call CheckStop( "ERROR: Negative iadv" )
       end if
      DepAdv2Calc(iadv) = Dep(i)%calc
 end do
  
!##################################################################
! We process the various combinations of gas-species and ecosystem:
! starting with DryDep, e.g. DDEP_SO2_m2CF

  do ndep = 1, nOutDDep

    OutDdep(ndep)%f2d  = find_index(OutDdep(ndep)%name ,f_2d(:)%name)
    OutDdep(ndep)%LC   = find_index(OutDdep(ndep)%txt ,DEF_ECOSYSTEMS)

    if(DEBUG_MY_DRYDEP .and. MasterProc) write(6,*) "OUTDdep ", ndep, &
       trim (OutDdep(ndep)%name),trim(OutDdep(ndep)%txt),"=>"&
          ,OutDdep(ndep)%LC, OutDdep(ndep)%f2d
    call CheckStop( OutDDep(ndep)%f2d < 1, &
        "OutDDep-f2d Error " // OutDDep(ndep)%name)
    call CheckStop( OutDDep(ndep)%LC < 0, & !0=GRID
        "OutDDep-LC Error " // OutDDep(ndep)%name)
  end do

!##################################################################
  do nVg = 1, nOutVg

     icheck = OutVg(nVg)%LC  ! TMP DEBUG
     if( OutVg(nVg)%txt == "Grid") then
        OutVg(nVg)%LC = FULL_GRID   ! zero
     else 
        OutVg(nVg)%LC = find_index( OutVg(nVg)%txt, LandDefs(:)%code )
     end if
     OutVg(nVg)%f2d = find_index( OutVg(nVg)%name, f_2d(:)%name )

     if(DEBUG_MY_DRYDEP .and. MasterProc) write(6,*) "OUTVG ", nVg, &
         trim( OutVg(nVg)%name ),trim( OutVg(nVg)%txt ),"=>",&
           OutVg(nVg)%LC, OutVg(nVg)%f2d

     call CheckStop( OutVg(nVg)%LC < FULL_GRID , & !ie zero
          "OutVg-LC Error " // OutVg(nVg)%name)
     call CheckStop( OutVg(nVg)%f2d < 1, &
          "OutVg-f2d Error " // OutVg(nVg)%name)
  end do

!##################################################################
 do nRG = 1, nOutRG

     if( OutRG(nRG)%txt == "Grid") then
        OutRG(nRG)%LC = FULL_GRID    ! zero
     else 
        OutRG(nRG)%LC = find_index( OutRG(nRG)%txt, LandDefs(:)%code )
     end if
     OutRG(nRG)%f2d = find_index( OutRG(nRG)%name, f_2d(:)%name )

     if(DEBUG_MY_DRYDEP .and. MasterProc) write(6,*) "OUTRG ", nRG, &
         OutRG(nRG)%name,OutRG(nRG)%txt,"=>",&
           OutRG(nRG)%LC, OutRG(nRG)%f2d

!     call CheckStop( OutRG(nRG)%LC < FULL_GRID , & ! <zero
!          "OutRG-LC Error " // OutRG(nRG)%name)
     call CheckStop( OutRG(nRG)%f2d < 1, &
          "OutRG-f2d Error " // OutRG(nRG)%name)
  end do

!##################################################################
 do nMET = 1, nOutMET

     OutMET(nMET)%f2d = find_index( OutMET(nMET)%name, f_2d(:)%name )

     if(DEBUG_MY_DRYDEP .and. MasterProc) write(6,*) "OUTMET ", nMET, &
         trim( OutMET(nMET)%name),trim( OutMET(nMET)%txt),"=>",&
           OutMET(nMET)%LC, OutMET(nMET)%f2d

!     call CheckStop( OutMET(nMET)%LC < FULL_GRID , & ! <zero
!          "OutMET-LC Error " // OutMET(nMET)%name)
     call CheckStop( OutMET(nMET)%f2d < 1, &
          "OutMET-f2d Error " // OutMET(nMET)%name)
  end do

!##################################################################
  do n = 1, nOutVEGO3

    OutVEGO3(n)%f2d  = find_index(OutVEGO3(n)%name ,f_2d(:)%name)

    if(DEBUG_MY_DRYDEP .and. MasterProc) write(6,*) "OUTVEGO3 ", n, &
       trim( OutVEGO3(n)%name)," ",trim( OutVEGO3(n)%txt),"=>"&
          ,OutVEGO3(n)%LC, OutVEGO3(n)%f2d
    call CheckStop( OutVEGO3(n)%f2d < 1, &
        "OutVEGO3-f2d Error " // OutVEGO3(n)%name)
    call CheckStop( OutVEGO3(n)%LC < 0, & !0=GRID
        "OutVEGO3-LC Error " // OutVEGO3(n)%name)
  end do

!#######################  NEW define indices here #######################

!D2_AFSTDF0 = find_index("D2_AFSTDF0",f_2d(:)%name)
!D2_AFSTDF16 = find_index("D2_AFSTDF16",f_2d(:)%name)
!
!D2_AFSTBF0 = find_index("D2_AFSTBF0",f_2d(:)%name)
!D2_AFSTBF16 = find_index("D2_AFSTBF16",f_2d(:)%name)
!
!D2_AFSTCR0 = find_index("D2_AFSTCR0",f_2d(:)%name)
!D2_AFSTCR3 = find_index("D2_AFSTCR3",f_2d(:)%name)
!D2_AFSTCR6 = find_index("D2_AFSTCR6",f_2d(:)%name)

!Oct09D2_O3DF    = find_index("D2_O3DF   ",f_2d(:)%name)
!Oct09D2_O3WH    = find_index("D2_O3WH   ",f_2d(:)%name)

!Oct09D2_EUAOT40WH    = find_index("D2_EUAOT40WH",f_2d(:)%name)
!Oct09D2_EUAOT40DF    = find_index("D2_EUAOT40DF",f_2d(:)%name)

!Oct09D2_UNAOT40WH    = find_index("D2_UNAOT40WH",f_2d(:)%name)
!Oct09D2_UNAOT40DF    = find_index("D2_UNAOT40DF",f_2d(:)%name)

!Oct09D2_MMAOT40WH    = find_index("D2_MMAOT40WH",f_2d(:)%name)

!iam_wheat   = find_index("IAM_CR",LandDefs(:)%code)
!iam_beech   = find_index("IAM_DF",LandDefs(:)%code)
!iam_medoak  = find_index("IAM_MF",LandDefs(:)%code)
!####################### ds END of define indices #######################
!  check_vals = (/  &
!    D2_AFSTDF0, D2_AFSTDF16, D2_AFSTBF0, D2_AFSTBF16, & !4
!    D2_AFSTCR0, D2_AFSTCR3, D2_AFSTCR6,& !3
!    iam_medoak, iam_beech, iam_wheat, &   ! For Fluxes !3
!    D2_O3DF,    D2_O3WH, & !2
!    D2_EUAOT40WH, D2_EUAOT40DF, & !4
!    D2_UNAOT40WH, D2_UNAOT40DF, & !4
!    D2_MMAOT40WH & !2
 
     if(MasterProc) write(6,*) "Init_DepMap D2D FINISHED"

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
     integer :: n, nadv, nadv2, ihh, idd, imm, iLC, iEco
     real :: veg_o3  ! O3 (ppb)  at canopy top, for  e.g. wheat, decid forest
     real :: Y       ! Threshold for flux


     real, parameter  :: NMOLE_M3 = 1.0e6*1.0e9/AVOG  ! Converts from 
                                                      ! mol/cm3 to nmole/m3

  !ECO08:  Variables added for new ecosystem dep
     real, dimension(NDEF_ECOSYSTEMS) :: invEcoFrac, EcoFrac
     integer :: ndep
     real :: Fflux

     real ::  to_nmole, timefrac, fstfrac 

     to_nmole =  NMOLE_M3
     timefrac = dt/3600.0
     fstfrac  = dt*1.0e-6     ! Converts also nmole to mmole

  ! Must match areas given above, e.g. DDEP_CONIF -> Conif

  !ECO08, Ecosystem areas, which were assigned in Init_DryDep:
  !  EcoFrac(CONIF)   = sum( coverage(:), LandType(:)%is_conif )
  !  EcoFrac(FULL_GRID)    = 1.0

     EcoFrac(:)    = EcoSystemFrac(:,i,j)
     invEcoFrac(:) = 0.0

     do n = 1, NDEF_ECOSYSTEMS
        if ( EcoFrac(n) > 1.0e-39 ) invEcoFrac(n) = 1.0/EcoFrac(n)
     end do 


   !  Query - crops, outisde g.s. ????
     !if ( DEBUG_MY_DRYDEP .and. first_call .and. debug_flag ) then
     if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
       write(*,*)  "ECOAREAS ", i,j
         do n = 1,  NDEF_ECOSYSTEMS
           write(*,"(a,i3,a,f14.4,g12.3)")  "ECOCHECK ", n, &
               DEF_ECOSYSTEMS(n), EcoFrac(n), invEcoFrac(n)
         end do
       write(*,*) "Done ECOCHECK ========================"
        
     end if
     first_call = .false.


    ! Traditional depositions to whole grid::

!     d_2d(DDEP_SOX,i,j,IOU_INST) = (  &
!          DepLoss(IXADV_SO2) + DepLoss(IXADV_SO4) ) * convfac * atwS
!
!       d_2d(D2_VddCOA,i,j,IOU_INST) =  100.0* Vg3m(CDEP_COA)

    ! Ecosystem depositions, for grouped or individual species:

     do ndep = 1, nOutDDep
        nadv  = OutDDep(ndep)%Index
        iEco   = OutDDep(ndep)%LC
     if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
         call CheckStop(iEco<1, "IECO ERROR: "//OutDDep(ndep)%name)
     end if
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
         Fflux = Fflux + Deploss(nadv) * &
             sum( fluxfrac(nadv,:), Is_EcoSystem(iEco,:) )
!OCT09        if ( RLC == FULL_GRID ) then 
!OCT09            Fflux = Fflux + Deploss(nadv) * sum( fluxfrac(nadv,:) )
!OCT09        else if ( RLC == CONIF ) then 
!OCT09            Fflux = Fflux + Deploss(nadv) * &
!OCT09               sum( fluxfrac(nadv,:), LandType(:)%is_conif )
!OCT09      .........
      end do ! n
        if ( OutDDep(ndep)%f2d < 1 .or. Fflux < 0.0 ) then
             write(6,*) "CATASTR ", ndep, OutDDep(ndep)%f2d,OutDDep(ndep)%name
             write(6,*) "CATASTR  iEco", iEco
             call CheckStop("CATASTROPHE: "//OutDDep(ndep)%name)
        end if

  !ECO08 
  ! - invEcoFracCF divides the flux per grid by the landarea of each
  ! ecosystem, to give deposition in units of mg/m2 of ecosystem.

        d_2d( OutDDep(ndep)%f2d,i,j,IOU_INST) =  &
            Fflux * convfac * OutDDep(ndep)%atw * invEcoFrac(iEco)

        if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
           write(6,"(a,3i4,3es12.3)") "DEBUG_ECO Fflux ", ndep, nadv, &
              iEco, Fflux, d_2d( OutDDep(ndep)%f2d,i,j,IOU_INST)
        end if ! DEBUG_ECO 
     end do ! ndep

! Fluxes, AFstY 
! o3WH = c_hvegppb(iam_wheat)* lossfrac

do n = 1, nOutVEGO3
     iLC = OutVEGO3(n)%LC
     Y   = OutVEGO3(n)%index    ! Index stores threshold Y x 10, or X as X

     if( OutVEGO3(n)%name(1:4) == "AFST")  then    !############### AFSTs ####
        Y = 0.1 * Y ! Index stores threshold Y x 10
        if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
           write(6,"(a,2i3,f4.1,es12.3)") "DEBUG_FST Fflux ", &
                n, iLC, Y, leaf_flux(iLC)
        end if

       ! Add fluxes:
        d_2d( OutVEGO3(n)%f2d ,i,j,IOU_INST) = fstfrac* max(leaf_flux(iLC)-Y,0.0)

     else if( OutVEGO3(n)%name(1:3) == "AOT")  then !############## AOTs ####


    !/-- Calcuates AOT values for specific veg. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )

        if ( Grid%izen < AOT_HORIZON ) then
          veg_o3 = c_hvegppb(iLC)* lossfrac

          d_2d( OutVEGO3(n)%f2d ,i,j,IOU_INST) = max(veg_o3-Y,0.0) * timefrac

          if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
            write(6,"(a,2i3,2f9.1,i4)") "DEBUG_AOT ", n, iLC, Y, veg_o3 &
               ,WheatGrowingSeason(i,j)  ! NOT USED YET....!
          end if
        else
          d_2d( OutVEGO3(n)%f2d ,i,j,IOU_INST) = 0.0
        end if
     else
        call CheckStop("OUTVEG UNDEF")
     end if
     
end do ! VEGO3

!
   !--- ecosystem specific concentrations..
   ! - use Conif forest for forests - safer for growing seasons

     imm      =    current_date%month            ! for debugging
     idd      =    current_date%day              ! for debugging
     ihh      =    current_date%hour             ! for debugging


! EXCLUDED THE REST PENDING RE_WRITE OF AOTx system - will be similar to
! the other "Out" variables, e.g OutVEGO3
!TMPDS
!TMPDS     if ( ihh >= 9 .and. ihh <= 21 ) then ! 8-20 CET, assuming summertime
!TMPDS        d_2d(D2_EUAOT40DF,i,j,IOU_INST) =  max(o3DF-40.0,0.0) * timefrac
!TMPDS     else
!TMPDS        d_2d(D2_EUAOT40DF,i,j,IOU_INST) =  0.0
!TMPDS     end if
!TMPDS
!????????????????????
!TMPDS       ! MM AOT added (same as UNECE, but different growing season)
!TMPDS             d_2d(D2_MMAOT40WH,i,j,IOU_INST) = d_2d(D2_UNAOT40WH,i,j,IOU_INST)&
!TMPDS                     * WheatGrowingSeason(i,j)
!TMPDS
!TMPDS    if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
!TMPDS          write(6,"(a12,3i5,f7.2,5es12.3,i3,es12.3)") "DEBUG_ECO ", &
!TMPDS          imm, idd, ihh, o3WH, &
!TMPDS             leaf_flux(iam_beech), d_2d(D2_AFSTDF0,i,j,IOU_INST), &
!TMPDS             leaf_flux(iam_wheat), d_2d(D2_AFSTCR0,i,j,IOU_INST), &
!TMPDS             d_2d(D2_UNAOT40WH,i,j,IOU_INST), &
!TMPDS             WheatGrowingSeason(i,j), d_2d(D2_MMAOT40WH,i,j,IOU_INST)
!TMPDS    end if ! DEBUG

   !---- end ecosystem specific ----------------------------------------------

  end subroutine  Add_ddep
  !<==========================================================================
  subroutine Add_Vg(debug_flag,i,j,VgLU)
     logical, intent(in) :: debug_flag
     integer, intent(in) :: i,j             ! coordinates
     real, dimension(:,0:), intent(in) :: VgLU  ! dim (NLANDUSE*nlu), 0=Grid
     integer :: n, nVg, cdep
     logical :: my_first_call = .true.

       if ( my_first_call ) then ! Some safety tests.
         if( MasterProc ) then
            do nVg = 1, nOutVg
              if( OutVg(nVg)%Index < 1 ) then
                call CheckStop( "ERROR: Unallocated OutVg Index "//OutVg(nVg)%name )
              end if
              cdep  =  DepAdv2Calc( OutVg(nVg)%Index ) ! Convert e.g. IXADV_O3
  
              if ( cdep < 1 ) then
                print *, "ERROR: OutVg name", OutVg(nVg)%name
                print *, "ERROR: Negative OutVg cdep", cdep, nVg, OutVg(nVg)%Index
                print *, "ERROR: DEPADV2CALC had size", size(DepAdv2Calc)
                do n = 1, size( DepAdv2Calc)
                   print *, "DEPADVLIST ", n, DepAdv2Calc(n)
                end do
                do n = 1, size( OutVg(:)%Index )
                   print *, "LIST ", n, OutVg(n)%Index
                end do
                call CheckStop( cdep  < 1 , "ERROR: Negative OutVg cdep")
              end if
            end do ! nVg
         end if ! MasterProc
         my_first_call = .false.
       end if

       do nVg = 1, nOutVg
         cdep  =  DepAdv2Calc( OutVg(nVg)%Index ) ! Convert e.g. IXADV_O3
         !call CheckStop( cdep  < 1 , "ERROR: Negative OutVg cdep")
         d_2d( OutVg(nVg)%f2d,i,j,IOU_INST) = VgLU( cdep, OutVg(nVg)%LC ) 

         if( DEBUG_MY_DRYDEP .and. debug_flag ) then
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
         if( OutRG(nRG)%Index < 1 ) then
            print *, "ERROR: OutRG Index ", nRG, OutRG(nRG)%name
            call CheckStop( "ERROR: Unallocated OutRG Index ")
         end if
         cdep  =  DepAdv2Calc( OutRG(nRG)%Index ) ! Convert e.g. IXADV_O3
                                                  ! to CDEP_O3
         call CheckStop( cdep  < 1 , "ERROR: Negative OutRG cdep")
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

         if( DEBUG_MY_DRYDEP .and. debug_flag ) then
              write(*,"(a,a,i3,2i4,es12.3)") "ADD_RG: ", OutRG(nRG)%name,&
                cdep,  OutRG(nRG)%LC,  OutRG(nRG)%f2d, &
                   d_2d( OutRG(nRG)%f2d,i,j,IOU_INST)
         end if
       end do

  end subroutine Add_RG
  !<==========================================================================
  !<==========================================================================
  subroutine Add_LCC_Met(debug_flag,i,j, LCC_Met )
     logical, intent(in) :: debug_flag
     integer, intent(in) :: i,j             ! coordinates
     real, dimension(:,:), intent(in) :: LCC_Met ! dim (nMET, NLANDUSE*nlu)
     integer :: n, nMET
     real :: met

       do nMET = 1, nOutMET
         n    =  OutMET(nMET)%Index  !ind = ustar or ..
         met  =  LCC_Met( n, OutMET(nMET)%LC )

         d_2d( OutMET(nMET)%f2d,i,j,IOU_INST) = met

         if( DEBUG_MY_DRYDEP .and. debug_flag ) then
              write(*,"(a,a,i3,2i4,es12.3)") "ADD_MET: ", OutMET(nMET)%name,&
                n,  OutMET(nMET)%LC,  OutMET(nMET)%f2d, &
                   d_2d( OutMET(nMET)%f2d,i,j,IOU_INST)
         end if
       end do

  end subroutine Add_LCC_Met
  !<==========================================================================


end module My_DryDep_ml

