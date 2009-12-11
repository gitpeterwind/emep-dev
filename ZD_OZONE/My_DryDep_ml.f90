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
   nMosaic, MosaicOutput &
  ,MMC_USTAR, MMC_INVL, MMC_RH, MMC_CANO3, MMC_VPD, MMC_FST &
  ,SOX_INDEX, OXN_INDEX, RDN_INDEX &  ! Equal -1, -2, -3
  ,DDEP_SOXGROUP, DDEP_RDNGROUP, DDEP_GROUP

 use AOTx_ml,             only : Calc_AOTx
 use CheckStop_ml,        only : CheckStop, StopAll
 use ChemChemicals_ml,    only : species 
 use ChemSpecs_adv_ml        !   e.g. NSPEC_ADV,IXADV_O3,IXADV_H2O2,
 use ChemSpecs_shl_ml,    only : NSPEC_SHL   ! For DDEP_SOXGROUP etc.
 use ChemGroups_ml,       only : DDEP_OXNGROUP  !DSGC
 use Derived_ml,          only : f_2d, d_2d
 use EcoSystem_ml,        only : DEF_ECOSYSTEMS, NDEF_ECOSYSTEMS, &
                                 EcoSystemFrac, FULL_GRID, Is_EcoSystem
 use GridValues_ml,       only : debug_proc, debug_li, debug_lj
 use LandDefs_ml,         only : LandDefs, LandType, iLC_grass
 use Landuse_ml,          only : WheatGrowingSeason, LandCover
 use LocalVariables_ml,   only : Grid, Sub  !=> izen  integer of zenith angle
 use ModelConstants_ml ,  only : atwS, atwN, DEBUG_MY_DRYDEP, DEBUG_CLOVER, &
                                    MasterProc , IOU_INST, IOU_YEAR
 use OwnDataTypes_ml,     only : print_Deriv_type
 use Par_ml,              only : li0, lj0, li1, lj1
 use PhysicalConstants_ml,only : AVOG
 use SmallUtils_ml,       only : find_index, NOT_FOUND
 use StoFlux_ml,          only : unit_flux, lai_flux
 use TimeDate_ml,         only : current_date
 use Wesely_ml
 implicit none
 private

  public :: Init_DepMap
  public :: Add_MosaicOutput

  !/** Variables used in deposition calculations
 
  ! DDEP_xx gives the index that will be used in the EMEP model
  ! WES_xx gives the index of the Wesely gas to which this corresponds

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
   integer :: iadv, i, i2, n, imc
   character(len=12) :: subclass

!#######################  ECO08 new mappng  #######################
  do i = 1, NDRYDEP_ADV  ! 22
      iadv = Dep(i)%adv
      if(DEBUG_MY_DRYDEP .and. MasterProc) &
         write(6,*) "DEPMAP   ", Dep(i)%adv, Dep(i)%calc
      call CheckStop( iadv < 1, "ERROR: Negative iadv" )
      DepAdv2Calc(iadv) = Dep(i)%calc
 end do
  
!##################################################################
! We process the various combinations of gas-species and ecosystem:
! starting with DryDep, e.g. DDEP_SO2_m2CF

  do imc = 1, nMosaic
    subclass =  MosaicOutput(imc)%subclass
    MosaicOutput(imc)%f2d  = find_index(MosaicOutput(imc)%name ,f_2d(:)%name)
    if( MasterProc) write(6,*) "Init_MyDry process", trim(MosaicOutput(imc)%name)
    call CheckStop( MosaicOutput(imc)%f2d < 1, &
        "MosaicOutput-f2d Error " // trim(MosaicOutput(imc)%name))

    select case ( subclass ) 
      case ( "DDEP" ) 

        MosaicOutput(imc)%LC   = find_index(MosaicOutput(imc)%txt ,DEF_ECOSYSTEMS)
        !call CheckStop( MosaicOutput(imc)%LC < 0, & !0=GRID (careful,FULL_GRID=1)
        call CheckStop( MosaicOutput(imc)%LC < 1, & !0=GRID (careful,FULL_GRID=1)
        "MosaicOutput-LC Error " // trim(MosaicOutput(imc)%name))

      case ( "VG", "Rs ", "Rns", "Gns", "MMAOT", "EUAOT", "AFST" ) ! could we use RG_LABELS? 

        if( MosaicOutput(imc)%txt == "Grid") then
           MosaicOutput(imc)%LC = FULL_GRID   ! zero
        else 
           MosaicOutput(imc)%LC = &
              find_index( MosaicOutput(imc)%txt, LandDefs(:)%code )
        end if
        call CheckStop( MosaicOutput(imc)%LC < FULL_GRID , & !ie zero
          "OutVg-LC Error " // MosaicOutput(imc)%name)

      case default

        if( MasterProc) write(6,*) "Init_MyDry skips ", trim(MosaicOutput(imc)%name)

    end select ! subclass

    if(DEBUG_MY_DRYDEP .and. MasterProc) then
          write(6,*) "MosaicOuts after DryDep Init: ", imc
          call print_Deriv_type( MosaicOutput(imc) )
    end if

  end do

 
     if(MasterProc) write(6,*) "Init_DepMap D2D FINISHED"

  end subroutine Init_DepMap

  !<==========================================================================
  subroutine Add_MosaicOutput(debug_flag,dt,i,j,convfac,fluxfrac,&
                VgLU, GsLU, GnsLU,LCC_Met)

  !<==========================================================================
     ! Adds deposition losses to ddep arrays
     logical, intent(in) :: debug_flag
     real,    intent(in) :: dt              ! time-step
     integer, intent(in) :: i,j             ! coordinates
     real,    intent(in) ::  convfac  !!???, lossfrac
     real, dimension(:,:), intent(in) :: fluxfrac  ! dim (NADV, NLANDUSE)
     real, dimension(:,0:),intent(in) :: VgLU  ! dim (NLANDUSE*nlu), 0=Grid
     real, dimension(:,0:),intent(in) :: GsLU, GnsLU  ! dim (0:NLANDUSE*nlu)
     real, dimension(:,0:), intent(in) :: LCC_Met   ! dim (nMET, NLANDUSE*nlu)

     integer :: n, nadv, nadv2, ihh, idd, imm, iLC, iEco
     integer :: imc, f2d, cdep
     real :: X, Y       ! Threshold for flux AFstY, AOTX
     real :: output     ! tmp variable
     !character(len=len(MosaicOutput(1)%subclass)) :: subclass
     character(len=12) :: subclass
     logical :: my_first_call = .true.
     logical :: first_vgr_call     ! reset each subroutine call

     real, parameter  :: NMOLE_M3 = 1.0e6*1.0e9/AVOG  ! Converts from 
                                                      ! mol/cm3 to nmole/m3
  ! Variables added for ecosystem dep
     real, dimension(NDEF_ECOSYSTEMS) :: invEcoFrac, EcoFrac
     real :: Fflux, Gs, Gns

     real ::  to_nmole, timefrac, nmole_s_to_mmole 

     to_nmole =  NMOLE_M3
     timefrac = dt/3600.0
     nmole_s_to_mmole  = dt*1.0e-6   ! Accumulates nmole/s to mmole
     cdep = -99                      ! set on first_vgr_call

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

     first_vgr_call = .true.
     do imc = 1, nMosaic
        subclass = MosaicOutput(imc)%subclass
        f2d      = MosaicOutput(imc)%f2d
        nadv     = MosaicOutput(imc)%Index  ! can be negatve for groups
        iLC      = MosaicOutput(imc)%LC

        output = 0.0  ! We only have instantaneous outputs, so can initialise
                      ! here and set d-2d at end

        if ( my_first_call ) then ! Some safety tests.
           ! Land-cover index can be zero, for FULL_GRID
            call CheckStop(iLC<0, "ILC ERROR: "//MosaicOutput(imc)%name)
            call CheckStop(f2d<1, "f2d ERROR:  "//MosaicOutput(imc)%name)
        end if
        if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
           write(6,"(a,a)",advance='no') "Add_Mosaic: "// &
                  trim(MosaicOutput(imc)%name), ", " // trim(subclass)
        end if

        select case ( subclass )
        case ( "DDEP" )

           iEco   = iLC  ! We rename for clarity. Eco landcovers can include several
                         ! land-cover classes, see EcoSystem_ml
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
           end do ! n

           if ( DEBUG_MY_DRYDEP .and. Fflux < 0.0 ) then
             write(6,"(a,3i4,a)") "DDEP Fflux CATASTR ", imc, f2d, iEco, &
                    trim(MosaicOutput(imc)%name)
             call CheckStop("CATASTROPHE: "//MosaicOutput(imc)%name)
           end if

        ! - invEcoFracCF divides the flux per grid by the landarea of each
        ! ecosystem, to give deposition in units of mg/m2 of ecosystem.

          output = Fflux * convfac * MosaicOutput(imc)%atw * invEcoFrac(iEco)

          if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
             write(6,"(3i4,3es12.3)") imc, nadv, iEco, Fflux, output
          end if ! DEBUG_ECO 

        case ( "METCONC" )    ! hard-coded bit n' pieces

          n   =  MosaicOutput(imc)%Index  !ind = ustar or ..
          if( n == MMC_USTAR  ) output = Sub(iLC)%ustar
          if( n == MMC_RH     ) output = Sub(iLC)%rh
          if( n == MMC_INVL   ) output = Sub(iLC)%invL
          if( n == MMC_CANO3  ) output = Sub(iLC)%cano3_ppb
          if( n == MMC_VPD    ) output = Sub(iLC)%vpd
          if( n == MMC_FST    ) output = Sub(iLC)%FstO3

          if ( DEBUG_CLOVER .and. debug_flag .and. &
               n == MMC_FST .and. LandType(iLC)%is_clover ) then
             write(6,"(a,3i4,i5,2i4,i6,es12.3)") "MDCV ", imc, n, iLC, &
                 current_date%month, current_date%day, &
                 current_date%hour, current_date%seconds,   output
          end if ! DEBUG_ECO 

        case ( "AFST" )    ! Fluxes, AFstY 
                           ! o3WH = c_hvegppb(iam_wheat)* lossfrac

          Y   = MosaicOutput(imc)%XYCL   ! threshold Y, nmole/m2/s

          if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
             write(6,"(2i3,f6.1,es12.3)") imc, iLC, Y, Sub(iLC)%FstO3
          end if

         ! Add fluxes:
          !output  = nmole_s_to_mmole * max(Sub(iLC)%FstO3 - Y,0.0)
          output  = max(Sub(iLC)%FstO3 - Y,0.0)

          if ( DEBUG_CLOVER .and. debug_flag .and. &
               LandType(iLC)%is_clover ) then
             write(6,"(a,3i4,i5,2i4,i6,2es12.3)") "MDCVA", imc, n, iLC, &
                 current_date%month, current_date%day, &
                 current_date%hour, current_date%seconds,Y, output
          end if ! DEBUG_ECO 

        case ( "MMAOT", "UNAOT", "EUAOT" )    ! AOTX

          !OLD: veg_o3 = c_hvegppb(iLC)* lossfrac  ! lossfrac????

          X      = MosaicOutput(imc)%XYCL   ! threshold X,  ppb.h

          call Calc_AOTx( subclass,iLC, Sub(iLC)%cano3_ppb, X,output) 

          if ( DEBUG_MY_DRYDEP .and. debug_flag ) then
              write(6,"(i5,2i3,f6.1,f9.3,i4,f12.4)") &
               current_date%hour, n, iLC, X, &
                Sub(iLC)%cano3_ppb , &  ! c_hvegppb(iLC), &
                 WheatGrowingSeason(i,j), output  ! NOT USED YET....!
          end if

        case ( "VG", "Rs ", "Rns", "Gns" ) ! could we use RG_LABELS? 

           if ( first_vgr_call ) then ! We pre-calculate some dep-related stuff
             cdep = DepAdv2Calc( nadv ) ! Convert e.g. IXADV_O3
             Gs   = GsLU( cdep, MosaicOutput(imc)%LC )
             Gns  = GnsLU( cdep, MosaicOutput(imc)%LC )
           end if
           first_vgr_call = .false.

           !It is easy to make mistakes with Vg, so we have som extra checks
           !here

           if ( DEBUG_MY_DRYDEP .and. cdep < 1 ) then
             print *, "ERROR: OutVgR name", MosaicOutput(imc)%name
             print *, "ERROR: Negative cdep", cdep, imc, MosaicOutput(imc)%Index
             print *, "ERROR: DEPADV2CALC had size", size(DepAdv2Calc)
             do n = 1, size( DepAdv2Calc)
                print *, "DEPADVLIST ", n, DepAdv2Calc(n)
             end do
             call CheckStop( cdep  < 1 , "ERROR: Negative cdep")
           end if

           if ( subclass == "VG" ) then

             output = VgLU( cdep, iLC )  ! CHECK iLC

           else if ( subclass == "Gs" ) then

             output = Gs

           else if ( subclass == "Gns" ) then

             output = Gns

           else if ( subclass == "Rs" ) then

             if( Gs < 1.0e-44 ) then
                 output = -999.0 
             else
                 output = 1.0/Gs
             end if

           else if ( subclass == "Rns" ) then

            if( Gns < 1.0e-44 ) then
               output = -999.0
            else
               output = 1.0/Gns
            end if
           end if ! subclass

           if( DEBUG_MY_DRYDEP .and. debug_flag ) then
              write(*,"(2i4,f9.3)") cdep, iLC, output

           end if
        case default
           if ( MasterProc ) then
              print *, "OUTVEG UNDEF name ",   MosaicOutput(imc)%name
              print *, "OUTVEG UNDEF subclass ",   MosaicOutput(imc)%subclass
           call CheckStop("OUTVEG UNDEF" // subclass )
           end if
        end select

        if( DEBUG_MY_DRYDEP .and. debug_flag ) then
              write(*,"(a,es12.3)") "ADDED output: ",  output
        end if
        d_2d( f2d,i,j,IOU_INST) = output
     
     end do ! Mosaic

     my_first_call = .false.

!
   !--- ecosystem specific concentrations..
   ! - use Conif forest for forests - safer for growing seasons

     imm      =    current_date%month            ! for debugging
     idd      =    current_date%day              ! for debugging
     ihh      =    current_date%hour             ! for debugging

   !---- end ecosystem specific ----------------------------------------------

  end subroutine  Add_MosaicOutput
  !<==========================================================================

end module My_DryDep_ml

