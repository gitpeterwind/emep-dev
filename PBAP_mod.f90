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
  implicit none
  private

  !/-- subroutines for PBAPs
  public ::  init_PBAPs,setup_PBAPs
  !/-- subroutines for Fungal Spores
  private :: Set_FungalSpores

  real,public, save, allocatable, dimension(:,:) :: FungalSpores 

  integer, public, parameter ::   NPBAP = 1 !Number of PBAPs currently implemented
                                            !(only fungal spores at the moment)

  real*8, DIMENSION(3), parameter  ::  &
  FUNG_PARAMS = [20.426d0, 275.82d0, &
         39300.0d0] !From Fungal paramterization, Eq. (2) of
                    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                    !DOI 10.1007/978-3-319-35095-0_121

  ! - allows rcbio in CM_Reactions, but we access elements with
  ! the natbio indices here. These much match the indices used in rcbio
  ! We only use rcbio for isoprene and terpenes so far,  since
  ! soil NO, NH3 emissions etc are dealt with through rcemis.

  ! We hard-code these indices, but only calculate emissions if needed
  ! Must match order of NATBIO to start with 
  character(len=16), save, dimension(NPBAP), public:: &
      EMIS_PBAP = [character(len=16):: &
             "FUNGAL_SPORES"]


  ! EmisNat is used for BVOC; soil-NO, also in futur for sea-salt etc.
  ! Main criteria is not provided in gridded data-bases, often land-use
  ! dependent.

  real, public, save, allocatable, dimension(:,:,:) :: &
     EmisBPAP       !  will be transferred to d_2d emis sums

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine init_PBAPs()
    allocate(FungalSpores(LIMAX,LJMAX)) !Spatial distribution of fungal spores
    FungalSpores = 0d0
  end subroutine init_PBAPs

  

  subroutine Set_FungalSpores(i,j)
    !!!!!!!
    !Fills FungalSpores(i,j)
    !Based on parameterization from
    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
    !DOI 10.1007/978-3-319-35095-0_121
    integer, intent(in) ::  i,j
    integer :: nLC, iLC, LC
    logical :: my_first_call = .true.
    real    :: f, ft, fn, ftn, F_FNG, LAIFac
    real    :: eFGS !, emissions Fungal spores TODO: Fix units
    real :: beta, bmin, bmax, bx, by ! for beta function
    real :: hfac

    if( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
       write(*,*)"PBAP_mod DEBUG FUNGAL_SPORES: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%FUNGAL_SPORES
    end if


          F_FNG = 0d0 !Fungal spores flux

           nlu = LandCover(i,j)%ncodes

           do iiL = 1,nlu
              LC = LandCover(i,j)%codes(iiL)
              if ( LandType(LC)%is_water) then
                cycle
              else if ( LandType(LC)%is_ice) then
                cycle
              else
                F_FNG = F_FNG + &
                LandCover(i,j)%fraction(iiL)*(FUNG_PARAMS(1)*(t2_nwp(i,j,1)-FUNG_PARAMS(2))+FUNG_PARAMS(3)*q(i,j,1)*LandCover(i,j)%LAI(iiL))
                !Eq.(2) of S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                !DOI 10.1007/978-3-319-35095-0_121, scaled by fraction
              end if
           end do !LC

           FungalSpores(i,j) = F_FNG

    if ( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
       i = debug_li
       j = debug_lj
       write(*,"(a,4i4)") "FUNGAL_SPORES_EMISSION: ",  1, limax, 1, ljmax
       write(*,"(a,2i4,2f12.4,es12.4)") "FUNGAL_SPORES_EMISSION: ", &
               daynumber, current_date%hour, t2_nwp(i,j,1), q(i,j,1)(i,j), Fungal_emiss(i,j)
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

  real, parameter :: unit_conv = 1.0d0 !TODO: Add conversion
  logical :: dbg

  if ( size(NPBAP) == 0  ) return   ! Number of PBAPs

  dbg = ( DEBUG%PBAB .and. debug_proc .and. &
          i==debug_li .and. j==debug_lj .and. current_date%seconds == 0 )


  
  rcbio(NATBIO%TERP,KG)    = (E_MTL+E_MTP) * biofac_TERP/Grid%DeltaZ
  EmisNat(NATBIO%TERP,i,j) = (E_MTL+E_MTP) * 1.0e-9/3600.0

  if ( USES%FUNGAL_SPORES ) then
      call Set_FungalSpores(i,j)
      rcemis(itot_FUNGALSPORES,KG) = FungalSpores(i,j)

  end if

  end subroutine setup_PBAPs

  !----------------------------------------------------------------------------
end module PBAP_mod
