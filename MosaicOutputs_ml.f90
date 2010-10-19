! <MosaicOutputs_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2010 met.no
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

module MosaicOutputs_ml
 use AOTx_ml,          only : Calc_AOTx
 use CheckStop_ml,  only: CheckStop, StopAll
 use ChemChemicals_ml, only : species
 use ChemSpecs_shl_ml, only : NSPEC_SHL
 use ChemSpecs_adv_ml, only : NSPEC_ADV
 use ChemGroups_ml,       only : DDEP_OXNGROUP,DDEP_SOXGROUP, &
                                 DDEP_RDNGROUP, NMAX_DDEP
 use DerivedFields_ml, only : f_2d, d_2d
 use EcoSystem_ml,     only : NDEF_ECOSYSTEMS, DEF_ECOSYSTEMS, &
    EcoSystemFrac, FULL_GRID, Is_EcoSystem
 use LandDefs_ml,      only : LandDefs, LandType, &
         Check_LandCoverPresent ! e.g. "CF"
 use LocalVariables_ml,   only : Sub, Grid
 use MetFields_ml
 use ModelConstants_ml, only : MasterProc, DEBUG => DEBUG_MOSAICS,&
   NLANDUSEMAX, IOU_INST, &
   atwS, atwN, SOX_INDEX, OXN_INDEX, RDN_INDEX ! indices for dep groups

 use OwnDataTypes_ml,  only: Deriv,O3cl, print_deriv_type, TXTLEN_DERIV, TXTLEN_SHORT
 use SmallUtils_ml, only: find_index
 use TimeDate_ml, only : current_date
 use Wesely_ml, only : NDRYDEP_CALC
 implicit none
 private

 ! From My_Derived we created:
 public ::  Init_MosaicMMC
 public ::  Add_MosaicMetConcs
 public ::  Add_MosaicRG
 public ::  Add_MosaicVG
 public ::  Add_MosaicVEGO3
 public ::  Add_MosaicDDEP

 ! From My_DryDep_ml we just move:
 public :: Init_MosaicOutputs
 public :: Add_MosaicOutput

 
 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO

 integer, public, save :: MMC_RH, MMC_CANO3, MMC_VPD, MMC_FST, &
     MMC_USTAR, MMC_INVL, MMC_GSTO, MMC_EVAP
 character(len=30),private, save :: errmsg = "ok"


 ! Mosaic-specific outputs, e.g. VG_CF_HNO3 or Rns_GR_NH3
  integer, public, save :: nMosaic = 0
  integer, public, parameter :: MAX_MOSAIC_OUTPUTS=100
  logical, private, parameter :: T=.true., F=.false.

  type(Deriv), public, &
     dimension( MAX_MOSAIC_OUTPUTS ), save :: MosaicOutput


 contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     subroutine Init_MosaicMMC(MOSAIC_METCONCS)
        character(len=*), dimension(:), intent(in) :: MOSAIC_METCONCS

    ! Set indices of mosaic metconc params for later use. Will be zero if 
    ! not found, but that's okay I hope...
      MMC_RH    = find_index("RH",MOSAIC_METCONCS)
      MMC_CANO3 = find_index("CanopyO3",MOSAIC_METCONCS)
      MMC_VPD   = find_index("VPD",MOSAIC_METCONCS)
      MMC_FST   = find_index("FstO3",MOSAIC_METCONCS)
      MMC_USTAR = find_index("USTAR",MOSAIC_METCONCS)
      MMC_INVL  = find_index("INVL",MOSAIC_METCONCS)
      MMC_GSTO  = find_index("GSTO",MOSAIC_METCONCS)
      MMC_EVAP  = find_index("EVAP",MOSAIC_METCONCS)

     end subroutine Init_MosaicMMC

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     subroutine Add_MosaicMetConcs(MOSAIC_METCONCS,MET_LCS,nMET)
        character(len=*), dimension(:), intent(in) :: MOSAIC_METCONCS
        character(len=*), dimension(:), intent(in) :: MET_LCS
        integer, intent(out) :: nMET
        integer :: ilab, n, iLC
        character(len=TXTLEN_DERIV) :: name
        real :: atw 


      !------------- Met data for d_2d -------------------------
      ! We find the various combinations of met and ecosystem,
      ! adding them to the derived-type array Mosaic_Met (e.g. => Met_CF)

      nMET = 0
      atw = -999.0 !fake
      do ilab = 1, size(MOSAIC_METCONCS)
        MET_LC: do n = 1, size(MET_LCS)

          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "MET_LCS", n, MET_LCS, (ilab == 1))
          if ( iLC < 0 ) cycle  MET_LC
          if ( MOSAIC_METCONCS(ilab)(1:6) == "Canopy" .or. & 
               MOSAIC_METCONCS(ilab)(1:5) == "FstO3" ) &
                 LandType(iLC)%flux_wanted  = .true.  ! Canopy calc in StoFlux
          !-------------End of Check if LC present in this array ------!

          nMET = nMET + 1
          name = trim ( MOSAIC_METCONCS(ilab) ) // "_"  // trim( MET_LCS(n) )

          nMosaic = nMosaic + 1
          call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, nMET" )
          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC,Threshold, scale, avg? rho Inst Yr Mn Day atw
           MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", "METCONC", MET_LCS(n), MOSAIC_METCONCS(ilab), &
                ilab, -99,iLC,-99.9,  F , 1.0,  T,   F,   F, T, T, T, atw)
!MESS          OutMET(nMET) = Dep_type( &
!           name, ilab, -99, iLC, "Mosaic", MOSAIC_METCONCS(ilab),  &
!                                              MET_LCS(n), 1.0, -99, "-")

          if( MOSAIC_METCONCS(ilab)(1:5) == "USTAR" )  then
              MosaicOutput(nMosaic)%unit  =   "m/s"
          else if( MOSAIC_METCONCS(ilab)(1:4) == "INVL" )  then
              MosaicOutput(nMosaic)%unit  =   "m"
          else if( MOSAIC_METCONCS(ilab)(1:8) == "CanopyO3" )  then
              MosaicOutput(nMosaic)%unit  =   "ppb"
          else if( MOSAIC_METCONCS(ilab)(1:5) == "FstO3" )  then
              MosaicOutput(nMosaic)%unit  =   "mmole/m2" ! accumulated
          else if( MOSAIC_METCONCS(ilab)(1:4) == "EVAP" )  then
              MosaicOutput(nMosaic)%avg       =  .false. ! accumulate
              MosaicOutput(nMosaic)%unit      =  "mm"
              MosaicOutput(nMosaic)%dt_scale  =  .true.
          end if

          if(DEBUG .and. MasterProc) call print_deriv_type(MosaicOutput(nMosaic))
        end do MET_LC !n
      end do ! ilab

 end subroutine Add_MosaicMetConcs
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine Add_MosaicRG(RG_LABELS,RG_LCS,RG_SPECS,nRG)
        character(len=*), dimension(:), intent(in) :: RG_LABELS ! eg Rs,Gs
        character(len=*), dimension(:), intent(in) :: RG_LCS    ! eg CF, Grid
        integer, dimension(:), intent(in) :: RG_SPECS  ! eg NH3
        integer, intent(out) :: nRG
        integer :: ilab, n, itot, iLC, iadv
        character(len=TXTLEN_DERIV) :: name

      !------------- Surface resistance for d_2d -------------------------
      ! We find the various combinations of gas-species and ecosystem,
      ! adding them to the derived-type array LCC_DDep (e.g. => RsO3_CF)

      nRG = 0
      do ilab = 1, size(RG_LABELS)
      do itot = 1, size(RG_SPECS)
        RG_LC: do n = 1, size(RG_LCS)

          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "RG_LCS", n, RG_LCS, (itot==1 .and. ilab == 1))
          if ( iLC < 0 ) cycle  RG_LC
          !-------------End of Check if LC present in this array ------!

          nRG = nRG + 1
          name = trim ( RG_LABELS(ilab) ) // "_"  // &
           trim( species(RG_SPECS(itot))%name ) // "_" // trim( RG_LCS(n) )

          iadv  =   RG_SPECS(itot) - NSPEC_SHL

          !OutRG(nRG)%label  = RG_LABELS(ilab)
          !OutRG(nRG)%txt  =   RG_LCS(n)

!          OutRG(nRG) = Dep_type(  &
!             name, iLC, iadv, -99,"Mosaic", RG_LABELS(ilab), RG_LCS(n),&
!                                             1.0, -99,  "-")
           nMosaic = nMosaic + 1
           call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, nRG" )

          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC,Threshold, scale, avg? rho Inst Yr Mn Day atw

          if( RG_LABELS(ilab)(1:1) == "R" )  then
             MosaicOutput(nMosaic) = Deriv( &
               name, "Mosaic", RG_LABELS(ilab), RG_LCS(n), "s/m", &
                iadv, -99,iLC,-99.9, F , 1.0,  T,   F,   F, T, T, T, -999.0)
          else if( RG_LABELS(ilab)(1:1) == "G" )  then
             MosaicOutput(nMosaic) = Deriv( &
               name, "Mosaic", RG_LABELS(ilab), RG_LCS(n), "cm/s", &
                iadv, -99,iLC,-99.9, F, 100.0,  T,   F,   F, T, T, T, -999.0)
          end if
        end do RG_LC !n
      end do ! itot
      end do ! ilab
 end subroutine Add_MosaicRG

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine Add_MosaicVG(VG_LABELS,VG_LCS,VG_SPECS,nVG)
        character(len=*), dimension(:), intent(in) :: VG_LABELS ! eg Vg
        character(len=*), dimension(:), intent(in) :: VG_LCS ! eg CF, Grid
        integer, dimension(:), intent(in) :: VG_SPECS  ! eg NH3
        integer, intent(out) :: nVG
        integer :: ilab, n, itot, iLC, iadv
        character(len=TXTLEN_DERIV) :: name


      !------------- Deposition velocities for d_2d -------------------------
      ! Add species and ecosystem depositions if wanted:
      ! We find the various combinations of gas-species and ecosystem,
      ! adding them to the derived-type array LCC_DDep (e.g. => VGO3_CF)

      nVg = 0
      do ilab = 1, size(VG_LABELS)
      do itot = 1, size(VG_SPECS)
        VG_LC: do n = 1, size(VG_LCS)

          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "VG_LCS", n, VG_LCS, &
                   write_condition=(itot==1 .and. ilab == 1))
          if ( iLC < 0 ) cycle  VG_LC
          !-------------End of Check if LC present in this array ------!
          nVg = nVg + 1

          name = trim( VG_LABELS(ilab) ) // "_"  // &
             trim( species(VG_SPECS(itot))%name ) // "_" // trim( VG_LCS(n) )
          iadv = VG_SPECS(itot) - NSPEC_SHL
          call CheckStop( iadv < 1 .or. iadv > NSPEC_ADV, &
                 " ERR: DDEP_SPECS: VG_SPECS" )

          nMosaic = nMosaic + 1
          call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, nVg" )

          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC,XYLC, scale, avg? rho Inst Yr Mn Day atw
          MosaicOutput(nMosaic) = Deriv(  &
              name, "Mosaic", VG_LABELS(ilab), VG_LCS(n), "cm/s", &
                iadv, -99,iLC, -99.9,  F, 100.0, T,   F,   F, T, T, T, -999.0)

        end do VG_LC !n
      end do ! i
      end do ! ilab
 end subroutine Add_MosaicVG


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine Add_MosaicVEGO3(VEGO3_OUTPUTS,nVEGO3)
   ! Diagnoses type, index and vegetation of POD or AOT type
   ! metrics: eg POD_1.0_IAM_MF, AOT_40_CF

        type(O3cl), dimension(:), intent(in) :: VEGO3_OUTPUTS
        !character(len=*), dimension(:), intent(in) :: VEGO3_OUTPUTS
        integer, intent(out) :: nVEGO3
        integer :: n, isep1, isep2, iLC
        character(len=TXTLEN_DERIV) :: name
        character(len=TXTLEN_SHORT) :: txt, txtnum, txt2, units
        real :: Threshold, scale, Y
        type(O3cl) :: veg
        integer :: testIX
        logical :: dt_scale

      !------------- VEGO3 stuff ----------------------------------------------
      ! For fluxes or AOTs we start with a formatted name, eg. POD_3.0_CF and
      !untangle it to get threshold Y (=3.0) and landcover type

      nVEGO3 = 0
      dt_scale = .false. ! for AFstY and AOT we need to reset this.
      VEGO3_LC: do n = 1, size(VEGO3_OUTPUTS%name)

         veg = VEGO3_OUTPUTS(n)
         name = veg%name
         !txt = veg%LC

         if( veg%class == "POD" ) then
            Threshold = veg%Threshold
            units = "mmole/m2"
            scale = 1.0e-6   ! Accumulates nmole/s to mmole (*dt_advec)
            dt_scale = .true.   ! Accumulates nmole/s to mmole (*dt_advec)
         else if( veg%class == "AOT" )  then
            Threshold = veg%Threshold
            units = "ppb.h"
            scale = 1.0/3600.0 ! AOT in ppb.hour
            dt_scale = .true.
            !if( veg%defn == "EU") txt = veg%defn     ! Store MM or EU here
         end if


          !------------------- Check if LC present in this array ------!
          iLC = Check_LandCoverPresent( "VEGO3_LCS", veg%TXTLC, .true. )
          if ( iLC < 0  ) cycle  VEGO3_LC
          if ( iLC > 0  ) LandType(iLC)%flux_wanted  = .true. 
          !-------------End of Check if LC present in this array ------!
          nVEGO3 = nVEGO3 + 1

          name = veg%name

           nMosaic = nMosaic + 1
           call CheckStop( NMosaic >= MAX_MOSAIC_OUTPUTS, &
                       "too many nMosaics, VEGO3" )
          !Deriv(name, class,    subc,  txt,           unit
          !Deriv index, f2d,LC,Threshold, scale dt_scale avg? rho Inst Yr Mn Day atw
          ! Use index for veg array. No need to set iadv for VEGO3. Always O3.
           MosaicOutput(nMosaic) = Deriv(  &
              name, veg%class,  veg%defn, veg%TXTLC, units, &
                n, -99,iLC,Threshold,  T,  scale,  F,   F,   F, T, T, F, -999.0)
              ! name, "Mosaic", veg%class, veg%defn, units, & ! Use defn, not LC
              !TEST name, "Mosaic", veg%class, veg%LC, units, &
              ! Txt has defn
              ! class is AOT or defn here
              !name, veg%class, veg%defn, veg%LC, units, &
                !-99, -99,iLC,Threshold,  T,  scale,  F,   F,   F, T, T, F, -999.0)

      end do VEGO3_LC !n
 end subroutine Add_MosaicVEGO3

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine Add_MosaicDDEP(DDEP_ECOS,DDEP_SPECS,nDD)
   ! Diagnoses type, index and vegetation of POD or AOT type
   ! metrics: eg POD_1.0_IAM_MF, AOT_40_CF

        character(len=*), dimension(:), intent(in) :: DDEP_ECOS
        integer, dimension(:), intent(in) :: DDEP_SPECS  ! eg NH3
        integer, intent(out) :: nDD
        integer :: i, n, ispec, iadv, iLC
        character(len=TXTLEN_DERIV) :: name
        character(len=TXTLEN_SHORT) :: txt, txtnum, txt2, units
        real :: Threshold, scale, Y, atw
        logical :: dt_scale
         


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
 end subroutine Add_MosaicDDEP

!<==========================================================================
 subroutine Init_MosaicOutputs()
   integer :: iadv, i, i2, n, imc
   character(len=TXTLEN_SHORT) :: subclass

!##################################################################
! We process the various combinations of gas-species and ecosystem:
! starting with DryDep, e.g. DDEP_SO2_m2CF

  do imc = 1, nMosaic
    subclass =  MosaicOutput(imc)%subclass
    MosaicOutput(imc)%f2d  = find_index(MosaicOutput(imc)%name ,f_2d(:)%name)
    if(DEBUG .and. MasterProc) &
       write(6,*) "Init_MyDry process", trim(MosaicOutput(imc)%name)
    call CheckStop( MosaicOutput(imc)%f2d < 1, &
        "MosaicOutput-f2d Error " // trim(MosaicOutput(imc)%name))

    select case ( subclass ) 
      case ( "DDEP" ) 

        MosaicOutput(imc)%LC   = find_index(MosaicOutput(imc)%txt ,DEF_ECOSYSTEMS)

        call CheckStop( MosaicOutput(imc)%LC < 1, & !0=GRID (careful,FULL_GRID=1)
        "MosaicOutput-LC Error " // trim(MosaicOutput(imc)%name))

      case ( "EUAOT" ) ! we use txt to specify MM or EU
         if(DEBUG .and. MasterProc) then
             write(6,"(a,i4,a,i3,f8.3)") "MosaicOuts EUAOTDEFN: ", imc, &
              trim(MosaicOutput(imc)%txt), MosaicOutput(imc)%LC, Grid%surf_o3_ppb
         end if

      case ( "VG", "Rs ", "Rns", "Gns", "POD", "AOT" ) ! could we use RG_LABELS? 

        if( MosaicOutput(imc)%txt == "Grid") then
           MosaicOutput(imc)%LC = FULL_GRID   ! zero
        else if( MosaicOutput(imc)%txt == "EU") then
           MosaicOutput(imc)%LC = FULL_GRID   ! zero
        else 
           MosaicOutput(imc)%LC = &
              find_index( MosaicOutput(imc)%txt, LandDefs(:)%code )
        end if
        call CheckStop( MosaicOutput(imc)%LC < FULL_GRID , & !ie zero
          "OutVg-LC Error " // MosaicOutput(imc)%name)

      case default

        if( MasterProc) write(6,*) "Init_MyDry skips ", &
            trim(MosaicOutput(imc)%name), ",subclass:", trim(subclass)

    end select ! subclass

    if(DEBUG .and. MasterProc) then
          write(6,*) "MosaicOuts Init: ", imc
          call print_Deriv_type( MosaicOutput(imc) )
    end if

  end do
 
     if(MasterProc) write(6,*) "Init_MosaicOutputs FINISHED"

  end subroutine Init_MosaicOutputs

!<==========================================================================
 subroutine Add_MosaicOutput(debug_flag,dt,i,j,convfac,DepAdv2Calc,fluxfrac,&
                Deploss)

  !<==========================================================================
     ! Adds deposition losses to ddep arrays
     logical, intent(in) :: debug_flag
     real,    intent(in) :: dt              ! time-step
     integer, intent(in) :: i,j             ! coordinates
     real,    intent(in) ::  convfac  !!???, lossfrac
     integer, dimension(:), intent(in) :: DepAdv2Calc
     real, dimension(:,:), intent(in) :: fluxfrac  ! dim (NADV, NLANDUSE)
     real, dimension(:), intent(in) :: Deploss

     integer, dimension(NMAX_DDEP) :: ddep_code 
     integer :: n, nadv, nadv2, iLC, iEco
     integer :: imc, f2d, cdep
     real :: X, Y       ! Threshold for flux AFstY, AOTX
     real :: output     ! tmp variable
     character(len=TXTLEN_SHORT) :: subclass, class
     logical :: my_first_call = .true.
     logical :: first_call = .true.     ! reset each subroutine call

  ! Variables added for ecosystem dep
     real, dimension(NDEF_ECOSYSTEMS) :: invEcoFrac, EcoFrac
     real :: Fflux, Gs, Gns

     cdep = -99                      ! set on first_vgr_call

  ! Must match areas given above, e.g. DDEP_CONIF -> Conif

  !Ecosystem areas, which were assigned in Init_DryDep:
  !  EcoFrac(CONIF)   = sum( coverage(:), LandType(:)%is_conif )
  !  EcoFrac(FULL_GRID)    = 1.0

     EcoFrac(:)    = EcoSystemFrac(:,i,j)
     invEcoFrac(:) = 0.0

     do n = 1, NDEF_ECOSYSTEMS
        if ( EcoFrac(n) > 1.0e-39 ) invEcoFrac(n) = 1.0/EcoFrac(n)
     end do 

   !  Query - crops, outisde g.s. ????
     if ( DEBUG .and. first_call .and. debug_flag ) then
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
!       d_2d(D2_VddCOA,i,j,IOU_INST) =  100.0* Vg3m(CDDEP_COA)

    ! Ecosystem depositions, for grouped or individual species:

     do imc = 1, nMosaic
        class = MosaicOutput(imc)%class
        subclass = MosaicOutput(imc)%subclass
        f2d      = MosaicOutput(imc)%f2d
        nadv     = MosaicOutput(imc)%Index  ! can be negatve for groups
        iLC      = MosaicOutput(imc)%LC
       if( class == "AOT" ) subclass = class
       if( class == "POD" ) subclass = class

        output = 0.0  ! We only have instantaneous outputs, so can initialise
                      ! here and set d-2d at end

        if ( my_first_call ) then ! Some safety tests.
           ! Land-cover index can be zero, for FULL_GRID
            call CheckStop(iLC<0, "ILC ERROR: "//MosaicOutput(imc)%name)
            call CheckStop(f2d<1, "f2d ERROR:  "//MosaicOutput(imc)%name)
        end if
        if ( DEBUG .and. debug_flag ) then
           write(6,"(a,a)",advance='no') "Add_Mosaic: "// &
                  trim(MosaicOutput(imc)%name), ", " // trim(subclass)
        end if

        select case ( subclass )
        case ( "DDEP" )

           iEco   = iLC  ! We rename for clarity. Eco landcovers can include several
                         ! land-cover classes, see EcoSystem_ml
           if ( nadv > 0 ) then  ! normal advectde species
              nadv2 = 1
              ddep_code(1) = nadv
           else if ( nadv == SOX_INDEX ) then
              nadv2 = size( DDEP_SOXGROUP )
              ddep_code(1:nadv2) = DDEP_SOXGROUP - NSPEC_SHL
           else if ( nadv == OXN_INDEX ) then
              nadv2 = size( DDEP_OXNGROUP )
              ddep_code(1:nadv2) = DDEP_OXNGROUP - NSPEC_SHL
           else  !if ( nadv == RDN_INDEX ) then
              nadv2 = size( DDEP_RDNGROUP )
              ddep_code(1:nadv2) = DDEP_RDNGROUP - NSPEC_SHL
           end if
   
           Fflux = 0.0
           do n = 1, nadv2
              nadv = ddep_code(n)
              Fflux = Fflux + Deploss(nadv) * &
                sum( fluxfrac(nadv,:), Is_EcoSystem(iEco,:) )
           end do ! n

           if ( DEBUG .and. Fflux < 0.0 ) then
             write(6,"(a,3i4,a)") "DDEP Fflux CATASTR ", imc, f2d, iEco, &
                    trim(MosaicOutput(imc)%name)
             call CheckStop("CATASTROPHE: "//MosaicOutput(imc)%name)
           end if

        ! - invEcoFracCF divides the flux per grid by the landarea of each
        ! ecosystem, to give deposition in units of mg/m2 of ecosystem.

          output = Fflux * convfac * MosaicOutput(imc)%atw * invEcoFrac(iEco)

          if ( DEBUG .and. debug_flag ) then
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
          if( n == MMC_GSTO   ) output = Sub(iLC)%g_sto
          if( n == MMC_EVAP   ) output = Sub(iLC)%EvapTransp
          if ( DEBUG .and. debug_flag .and. &
               n==MMC_CANO3 .and. iLC == 2 ) then !DF
             write(6,"(a,4i5,f10.4)") "MYDDEP CANO3 ", &
                 current_date%month, current_date%day, &
                 current_date%hour, current_date%seconds,   output
          end if ! DEBUG_ECO 

          !if ( DEBUG_CLOVER .and. debug_flag .and. &
          !     !n == MMC_FST .and. LandType(iLC)%is_clover ) then
          !     n == MMC_EVAP .and. LandType(iLC)%is_forest ) then
          !   write(6,"(a,3i4,i5,2i4,i6,es12.3)") "MDCV ", imc, n, iLC, &
          !       current_date%month, current_date%day, &
          !       current_date%hour, current_date%seconds,   output
          !end if ! DEBUG_ECO 

        case ( "POD" )    ! Fluxes, AFstY 
                           ! o3WH = c_hvegppb(iam_wheat)* lossfrac

          Y   = MosaicOutput(imc)%Threshold   ! threshold Y, nmole/m2/s

          if ( DEBUG .and. debug_flag ) then
             write(6,"(2i3,f6.1,es12.3)") imc, iLC, Y, Sub(iLC)%FstO3
          end if

         ! Add fluxes:
          output  = max(Sub(iLC)%FstO3 - Y,0.0)

          !if ( DEBUG_CLOVER .and. debug_flag .and. &
          !     LandType(iLC)%is_clover ) then
          !   write(6,"(a,3i4,i5,2i4,i6,2es12.3)") "MDCVA", imc, n, iLC, &
          !       current_date%month, current_date%day, &
          !       current_date%hour, current_date%seconds,Y, output
          !end if ! DEBUG_ECO 

        !!case ( "EUAOT" )    ! AOTX using 3m grid O3

        !  X      = MosaicOutput(imc)%Threshold   ! threshold X,  ppb.h

          !call Calc_AOTx( MosaicOutput(imc)%txt,iLC, Grid%surf_o3_ppb, X,output,&
        !  defn = VEG
        !  call Calc_AOTx( VEGO3_OUTPUTS( MosaicOutput(imc)%index)%defn,iLC, Grid%surf_o3_ppb, X,output,&
        !         debug_flag, "EU-AOT " // trim(subclass) ) 
!
        case ( "AOT" )    ! AOTX

          !OLD: veg_o3 = c_hvegppb(iLC)* lossfrac  ! lossfrac????

          X      = MosaicOutput(imc)%Threshold   ! threshold X,  ppb.h

         ! txt has EU or MM
          call Calc_AOTx( MosaicOutput(imc)%subclass, iLC, X,output,&
                 debug_flag, "AOTLC " // trim(MosaicOutput(imc)%name) ) 
          !call Calc_AOTx( MosaicOutput(imc)%txt,iLC, Sub(iLC)%cano3_ppb, X,output,&
          !TEST call Calc_AOTx( subclass,iLC, Sub(iLC)%cano3_ppb, X,output,&

!          if ( DEBUG_AOT .and. debug_flag ) then
!              write(6,"(a,4i3,i5,f8.3,i3,g12.4)") &
!               "AOTLC " // trim(subclass), iLC, &
!                 current_date%month, current_date%day, &
!                 current_date%hour, current_date%seconds,&
!                Sub(iLC)%cano3_ppb , &  ! c_hvegppb(iLC), &
!                 WheatGrowingSeason(i,j), output  ! NOT USED YET....!
!!               current_date%hour, n, iLC, X, &
!          end if

        case ( "VG", "Rs ", "Rns", "Gns" ) ! could we use RG_LABELS? 

             cdep = DepAdv2Calc( nadv ) ! e.g. IXADV_O3 to calc index
             Gs   = Sub(iLC)%Gsur(cdep)
             Gns  = Sub(iLC)%Gns(cdep)

           !It is easy to make mistakes with Vg, so we have som extra checks
           !here

           if ( DEBUG .and. cdep < 1 ) then
             print *, "ERROR: OutVgR name", MosaicOutput(imc)%name
             print *, "ERROR: Negative cdep", cdep, imc, MosaicOutput(imc)%Index
             print *, "ERROR: DEPADV2CALC had size", size(DepAdv2Calc)
             do n = 1, size( DepAdv2Calc)
                print *, "DEPADVLIST ", n, DepAdv2Calc(n)
             end do
             call CheckStop( cdep  < 1 , "ERROR: Negative cdep")
           end if

           if ( subclass == "VG" ) then

             output = Sub(iLC)%Vg_3m(cdep)  ! CHECK iLC

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

           if( DEBUG .and. debug_flag ) then
              write(*,"(2i4,f9.3)") cdep, iLC, output

           end if
        case default
           if ( MasterProc ) then
              print *, "OUTVEG UNDEF name ",   MosaicOutput(imc)%name
              print *, "OUTVEG UNDEF subclass ",   MosaicOutput(imc)%subclass
           call CheckStop("OUTVEG UNDEF" // subclass )
           end if
        end select

        if( DEBUG .and. debug_flag ) then
              write(*,"(a,es12.3)") "ADDED output: ",  output
        end if
        d_2d( f2d,i,j,IOU_INST) = output
     
     end do ! Mosaic

     my_first_call = .false.


  end subroutine  Add_MosaicOutput
  !<==========================================================================

end module MosaicOutputs_ml
