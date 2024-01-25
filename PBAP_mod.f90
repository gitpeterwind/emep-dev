!> <PBAP_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!  **************************************************************************! 

module PBAP_mod

  !/-- Module to deal with primary biological aerosol pollutants (PBAPs).
  !
  !    Gunnar Felix Lange 2023
  !---------------------------------------------------------------------------

  use CheckStop_mod,      only: CheckStop, StopAll
  use ChemSpecs_mod,         only : species
  use Config_module, only : NPROC, MasterProc, TINY, &
                           NLANDUSEMAX, IOU_INST, & 
                           KT => KCHEMTOP, KG => KMAX_MID, & 
                           EURO_SOILNOX_DEPSCALE, & 
                           MasterProc, &
                           USES, &
                           NATBIO, EmBio, EMEP_EuroBVOCFile
  use Debug_module,       only: DebugCell, DEBUG
  use GridValues_mod,     only: i_fdom,j_fdom, debug_proc,debug_li,debug_lj
  use Io_mod,             only: IO_FORES, open_file, ios, datewrite
  use Io_RunLog_mod,      only: PrintLog
  use KeyValueTypes,      only: KeyVal,KeyValue
  use LandDefs_mod,       only: LandType, LandDefs
  use LandPFT_mod,        only: MapPFT_LAI, pft_lai
  use Landuse_mod,        only: LandCover
  use LocalVariables_mod, only: Grid  ! -> izen, DeltaZ
  use MetFields_mod,      only: t2_nwp, q
  use MetFields_mod,      only: PARdbh, PARdif !WN17, in W/m2
  use NetCDF_mod,         only: ReadField_CDF, printCDF
  use OwnDataTypes_mod,   only: Deriv, TXTLEN_SHORT
!  use Paleo_mod, only : PALEO_modai, PALEO_miso, PALEO_mmon
  use Par_mod,            only: MSG_READ1,me, limax, ljmax
  use PhysicalConstants_mod,  only:  AVOG, GRAV, PI
  use Radiation_mod,      only: PARfrac, Wm2_uE
  use SmallUtils_mod,     only: find_index
  use TimeDate_mod,       only: current_date, daynumber
  use ZchemData_mod,      only: rcemis, rcbio
  use Biogenics_mod,      only: NEMIS_BioNat,EMIS_BioNat,EmisNat
  implicit none
  private

  !/-- subroutines for PBAPs
  public ::  init_PBAPs,setup_PBAPs
  !/-- subroutines for Fungal Spores
  private :: Set_FungalSpores

  real,public, save, allocatable, dimension(:,:) :: FungalSpores 

  integer, public, parameter ::   NPBAP = 1 !Number of PBAPs
                                            !(only fungal spores implemented at the moment)

  integer, private, save :: itot_FUNGALSPORES, inat_FUNGALSPORES !Index of fungal spores in Species and EMIS_BioNat

  real*8, DIMENSION(3), parameter  ::  &
  FUNG_PARAMS = [20.426, 275.82, 39300.0] !From Fungal paramterization, Eq. (2) of
                    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                    !DOI 10.1007/978-3-319-35095-0_121

  real, parameter :: unit_conv = 1e-6
                !Original equation outputs number concentration in m^{-2}s^{-1}, but rcemis
                !wants number molc/cm3/s
                !Conversion: Divide by delta_z (in m), multiply by 1e-6 m^{-3}->cm^{-3}
  real, public, save, allocatable, dimension(:,:,:) :: &
     EmisBPAP       !  will be transferred to d_2d emis sums

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine init_PBAPs()
    logical, save ::  my_first_call = .true.

    if (my_first_call) then
      allocate(FungalSpores(LIMAX,LJMAX)) !Spatial distribution of fungal spores
      FungalSpores = 0.0
      itot_FUNGALSPORES = find_index( "FUNGAL_SPORES", species(:)%name)
      inat_FUNGALSPORES = find_index( "FUNGAL_SPORES", EMIS_BioNat(:))
      my_first_call = .false.
      if( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
        write(*,*)"INIT PBAPs (should only happen once!)"
      end if
    end if
  end subroutine init_PBAPs

  

  subroutine Set_FungalSpores(i,j)
    !!!!!!!
    !Fills FungalSpores(i,j)
    !Based on parameterization from
    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
    !DOI 10.1007/978-3-319-35095-0_121
    integer, intent(in) ::  i,j
    integer :: nlu,iiL,LC,i_d,j_d
    real    :: F_FNG, temp_val

    if( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
      if (i .eq. debug_li .and. j .eq. debug_lj) then
       write(*,*)"PBAP_mod DEBUG FUNGAL_SPORES: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%FUNGAL_SPORES, itot_FUNGALSPORES,inat_FUNGALSPORES
      end if
    end if


          F_FNG = 0.0 !Fungal spores flux

           nlu = LandCover(i,j)%ncodes

           do iiL = 1,nlu
              LC = LandCover(i,j)%codes(iiL)
              if ( LandType(LC)%is_water) then
                cycle
              else if ( LandType(LC)%is_ice) then
                cycle
              else
                temp_val =  LandCover(i,j)%fraction(iiL)*(FUNG_PARAMS(1)*(t2_nwp(i,j,1)-FUNG_PARAMS(2))+FUNG_PARAMS(3)*q(i,j,KG,1)*LandCover(i,j)%LAI(iiL))
                temp_val = (unit_conv*temp_val)/Grid%DeltaZ
                F_FNG = F_FNG + max(0.0, temp_val)
                !Eq.(2) of S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                !DOI 10.1007/978-3-319-35095-0_121, scaled by fraction
                !Leaf-area index (LAI) should be in m2/m2
                !Specific humidity q should be in kg/kg
                !Temperature at 2m (t2_nwp) should be in K

              end if
           end do !LC

           FungalSpores(i,j) = F_FNG

    if ( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
       if (i .eq. debug_li .and. j .eq. debug_lj) then
        write(*,"(a,4i4)") "FUNGAL_SPORES i,j: ",  1, limax, 1, ljmax
        write(*,"(a,2i4)") "FUNGAL_SPORES indices: ",  itot_FUNGALSPORES,inat_FUNGALSPORES
        write(*,"(a,2i4,4f12.4)") "FUNGAL_SPORES_EMISSION: ", &
               daynumber, current_date%hour, t2_nwp(i,j,1), q(i,j,KG,1), FungalSpores(i,j),LandCover(i,j)%LAI(1)
       end if
    end if

  end subroutine Set_FungalSpores

  subroutine setup_PBAPs(i,j)
  !
  !---- Adds PBAPs to rcemis  ------------------------------------------------
  !
  !  So far, only adds Fungal Spores to rcemis
  !
  !  Called from setup_1d_mod, every  advection step.
  !----------------------------------------------------------------------------

  integer, intent(in) ::  i,j

  character(len=*), parameter :: dtxt='PBAPModSetup:' 

  logical :: dbg

  if ( NPBAP == 0  ) return   ! Number of PBAPs

  dbg = ( DEBUG%PBAP .and. debug_proc .and. &
          i==debug_li .and. j==debug_lj .and. current_date%seconds == 0 )


  if ( USES%FUNGAL_SPORES ) then
      call Set_FungalSpores(i,j)
      if (itot_FUNGALSPORES > 0) then
        rcemis(itot_FUNGALSPORES,KG) = rcemis(itot_FUNGALSPORES,KG)+(unit_conv*FungalSpores(i,j))/Grid%DeltaZ
        !Recmis in molec/cm3/s
        if (dbg) then
          write(*,*) "FUNGAL SPORES: rcemis ",rcemis(itot_FUNGALSPORES,KG)
        end if

      end if
      
      if (inat_FUNGALSPORES > 0) then
        EmisNat(inat_FUNGALSPORES,i,j) = FungalSpores(i,j)
        !Emissions in molec/m2/s
        if (dbg) then
          write(*,*) "FUNGAL SPORES: EmisNat ",EmisNat(inat_FUNGALSPORES,i,j)
        end if
      end if
  end if

  end subroutine setup_PBAPs

  !----------------------------------------------------------------------------
end module PBAP_mod
