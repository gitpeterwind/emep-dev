!----------------------------------------------------------------------------!
!> MODULE esx_Variables
!! Stores and configs various run settings.
!----------------------------------------------------------------------------!
module esx_Variables
  use CheckStop_ml, only : CheckStop
  use ChemSpecs, only : species
  use KeyValueTypes, only: KeyValInt, KeyValReal
  use LocalVariables, only : LocDat   !! Data for Sub
  use ModelConstants, only : UNDEF_R
  use PhysicalConstants_ml, only : T0    !! 273.15
  use SmallUtils_ml, only: LenArray, find_indices
  use TimeDate_ml,       only : date, current_date
  implicit none
  private

  public :: Config_esx  !> Sets run params

  !> Some maximum limits:
  !! (Set these limits lowish to stop the output namelist
  !!  arrays becoming excessive. No harm done, but harder to read.)
  integer, parameter, public :: ESX_MAXNZ =   200 ! No. vertical layers
  integer, parameter, public :: ESX_MAXNDIFF = 50 ! No. diffusive species
  integer, parameter, public :: ESX_MAXNOUT  = 50 ! No. output species

  !> For diffusing species we will allow Vd, Ve, Fb, Ft, so define a type
  !! For units explanation, see config file.

  type DiffSpec_t
    character(len=30) :: name = '-'
    integer           :: ind  = -999
    real              :: Vd   = 0.0   ! Deposition to ground, m/s
    real              :: Ve   = 0.0   ! Escape Vel. at top, m/s
    real              :: Fb   = 0.0   ! mass flux at ground, units as conc
    real              :: Ft   = 0.0   ! mass flux at top, units as conc
  end type DiffSpec_t


  !> ESX general  information. Can be reset in Config

  type esx_t


   !> Some variables are calculated. No need to set directly in config file

    integer :: nz = 0          !> Number of layers. Calculated from esx%zbnd
    integer :: nhveg = 0       !> Number of canopy layers (<= nz).
    integer :: nDiffSpecs = 0  !> Number of diffusing species
    integer :: nOutSpecs = 0   !> Number of output species
  
   !> Use settings. Can skip various modules if not required. May use
   !! character strings one day, but for now simple logicals

    logical :: uses_chem = .false.
    logical :: uses_diff = .false.
    logical :: uses_veg  = .false.
    logical :: uses_plotting = .false.

   !> The z-grid variables

    character(len=10) :: zmethod = 'explicit'   ! explicit, or linear, or ...
    real :: zfirst = -999., zlast =  -999.      ! used if zbnd calculated ('linear' or)

    real, dimension(ESX_MAXNZ) :: & ! zn
       z     &!< grid centre point heights: z_1, z_2,..., z_n
      ,dz    &!< layer heights between boundaries: dz_1, dz_2,..., dz_n (dz_n=zbnd_n+1-zbnd_n)
  !  real :: & ! only zn-1 used, but we dimension zn anyway
     ,zbnd = UNDEF_R  &!< grid boundary heights: z_1½, z_2½,..., z_n-½
     ,dzmid  !< layer height between grid mid-points: dz_1½, dz_2½,..., dz_n-½ (dzmid_n=z_n+1-z_n)

   !> Remaining variables can be set in config_esx.nml. See that file for more info

    character(len=30) :: exp_name    = "-"  !> Run label
    character(len=30) :: units       = "-"  !> ppb or - so far
    character(len=30) :: DataSource  =  "-"
    character(len=30) :: odir        =  "." !> output directory

   !> Start, finish, step, etc.  - basic run times
   !! Usually reset in config_esx.nml

    type(date) :: startdate, indate
    real ::&
      startTime = 0, endTime=0, stepTime=1.0e10, Time=0   ! dummy values for safety
    real :: dt_extern =    0.0  ! Master time step, for e.g. met data change
    real :: dt_phychem=    0.0  ! time step around chem/diffusion loops 
    real :: dt_Zdiff = 9.99e10 ! time step for diffusion. Set in config or calculated

   !>  Top boundary

    logical :: fixedBC = .false.    ! Set true for fixed BC at top, false for free

    integer :: &
       nsteps = 0 &!! Total number of data records
      ,year   = 0 &!!
      ,daynumber = 0 &!! Ordinal day number
      ,ns     = 0  !! counter in 1.. nsteps
   !
   !> Some times for output (allow 100 for now)
   !! and debug settings: 0 for none, 1 for light, 2 for extra
    real, dimension(100) :: printTime

  !> Debugging flags. Use 0 for none, 1 or 2 for 
    integer :: debug_driver= 0   !> Extra output
    integer :: debug_Zdiff = 0   !> Extra output from KdiffSolver
    integer :: debug_Zchem = 0   !> Extra output from ChemSolver
    integer :: debug_Zmet  = 0   !> Extra output
    integer :: debug_Zveg  = 0   !> Extra output from 1-D , e.g. PAR
    integer :: debug_gleaf = 0   !> Extra output, e.g. fphen

   !> Diffusing species have name and index:

    type(DiffSpec_t), dimension(ESX_MAXNDIFF) :: &
       DiffSpecs = DiffSpec_t() !> Diffusing species

   !> Output species:

    type(KeyValInt), dimension(ESX_MAXNOUT) :: &
       OutSpecs = KeyValInt("-",-999) !> Diffusing species


   !> Kz settings. We use key-word arguments to pass numerical
   !! values, since different profile methods have very different
   !! parameters. Allow here for max 5 kwargs:
    character(len=30) :: Kz_method = "-"
    type(KeyValReal), dimension(5) :: Kz_kwargs = KeyValReal("-",0)

   !> Stomatal conductance settings. Has to be set in config
    character(len=30) :: gsto_method = '-'
    character(len=30) :: gleaf_method = '-'

   !> Plot settings, is esx%use_plotting
    character(len=100) :: plot_cmds = ""

  end type esx_t

  type(esx_t), public, target, save :: esx

 !> Location-dependent scalar (z-less) data, e.g. lat, long, ustar, etc.
 !!  are kept in L. 
 !! Also easy to change via namelist.

  type(LocDat), public, save :: Loc


 !> ============ end of config_esx ==========================================


 !> z-dependent meteo data kept here
 !! Will use surface values unless reset.

  type, public :: Zmet_t
    real :: &
       rh  = UNDEF_R     &!< RH (fraction)
      ,Hum = UNDEF_R     &! Specific humidity (UNIT=
      ,VPD = UNDEF_R     &!< VPD (kPa)
      ,Pa  = UNDEF_R     &!< pressure (Pa)
      ,tzK = UNDEF_R     &!< T (degrees K)
      ,tzC =  UNDEF_R    &!< T (degrees C)
      ,tleafC = UNDEF_R  &!< Leaf temperature (degrees C)
      ,Kz  = 0.0         &!< Eddy diffusivity at grid boundaries
      ,Kz2  = 0.0        &!< Eddy diffusivity - TESTING above hSL
      ,Kz3  = 0.0        &!< Eddy diffusivity - TESTING above hSL
      ,Ra   = UNDEF_R    &!< Resistance between layers,Ra1 between z1 to z2
      ,uz  = 0.1         &!< Wind speed (m/s)
      ,CO2 = 392.0       &!< CO2 concentration (ppm)
      ,M  =  UNDEF_R      ! Air concentration, 3rd body for chem reactions
  end type

  type(Zmet_t), target, dimension(ESX_MAXNZ), public, save :: Zmet = Zmet_t()

  !> Generic per-layer vegetation data.

  type, public :: Zveg_t
    real :: &
       dLAI   = 0.0    &!< leaf area index per layer
      ,cumLAI = 0.0  &!< top-down cumulative dLAI
      ,PARz   = 0.0    &!< PAR
      ,gleaf  = 0.0    &!< total leaf conductance
      ,gsto   = 0.0    &!< stomatal conductance
      ,gns    = 0.0      !< non-stomatal conductance
  end type

  type(ZVeg_t), target, dimension(ESX_MAXNZ), public, save :: Zveg = Zveg_t()


contains
  subroutine Config_esx(io, writelog)
    integer, intent(in) :: io
    logical, intent(in) :: writelog
    integer :: ispec, ilog, nspecs, iz
    real :: zinc
    namelist /esxDriver_config/ esx, Loc

    read (io, nml=esxDriver_config)

   ! We convert the date to daynumber
     
    current_date = esx%startdate

   ! Now we can set the z-grid

    if ( esx%zmethod == 'linear' ) then
      print *, "Zbnd from linear range: ",  esx%nz
      call CheckStop( esx%zlast<0.0 .or. esx%nz<=0,&
       "zmethod linear range needs zfirst, zlast and zn")

      esx%zbnd(:) =  0.0   ! clears any other namelist settings
      zinc = ( esx%zlast - esx%zfirst ) / esx%nz

      do iz = 1, esx%nz
        esx%zbnd(iz) = esx%zfirst + zinc *  iz
        if ( esx%debug_driver > 1 ) then
          print *, "Zbnd from linear range: ", iz, esx%zbnd(iz), &
            maxval( esx%zbnd), esx%nz, maxloc(esx%zbnd,dim=1)
        end if
      end do
      
    else  ! derived from explicit config data:

print *, "ZBND 1 ", esx%zbnd(1), UNDEF_R
      call CheckStop( esx%zbnd(1) == UNDEF_R, "zbnds not set. Check config file")

      ! Check for mal-formed input of zbnds (e.g. double defn)
      do iz = 2, ESX_MAXNZ
         if ( esx%zbnd(iz) > 0.0 .and. (esx%zbnd(iz) < esx%zbnd(iz-1)) ) then
              print *, "zbnds ERRROR", iz, esx%zbnd(iz),  esx%zbnd(iz-1)
              call CheckStop( "zbnds not increasing linearly. Check config file")
         end if
      end do
      esx%nz     = maxloc(esx%zbnd,dim=1)


    end if 
    !! find the indices in the chemical mechanism corresponding
    !! to esx DiffSpecs

    nspecs = LenArray( esx%DiffSpecs%name, "-" )
    esx%nDiffSpecs = nspecs

    esx%DiffSpecs(1:nspecs)%ind =  &
       find_indices( esx%DiffSpecs(1:nspecs)%name, species(:)%name )
    
    nspecs = LenArray( esx%OutSpecs%key, "-" )
    esx%nOutSpecs = nspecs
    esx%OutSpecs(1:nspecs)%int =  &
       find_indices( esx%OutSpecs(1:nspecs)%key, species(:)%name )

    do ispec = 1, nspecs
       if( esx%DiffSpecs(ispec)%ind < 1 ) & 
         call CheckStop("DiffSpec not in CM:"//esx%DiffSpecs(ispec)%name )
       if( esx%OutSpecs(ispec)%int < 1 ) &
         call CheckStop("OutSpec not in CM:"//esx%OutSpecs(ispec)%key )
    end do

    if( writelog ) then
      call system("mkdir " // esx%odir)
      open(newunit=ilog,file=trim(esx%odir)//"/LogConfig.esx")
      write(ilog,"(a)") "CONFIG ESX===================================="
      write(ilog, nml=esxDriver_config)
      close(ilog)
    end if

    if ( esx%debug_driver > 0 ) then
      write(*,"(a)"   ) "CONFIG ESX===================================="
      write(*,    nml=esxDriver_config)
    end if

  end subroutine Config_esx

end module esx_Variables

!program test_config
!  use iso_fortran_env, only: INPUT_UNIT
!  use config
!
!  ! Read namelist with the derived type instance from config
!  namelist /esx_config/ conf
!  read (unit=INPUT_UNIT, nml=esx_config) ! Read from stdin
!  print *, "exp_name=", conf%exp_name
!  print *, "tstart=", conf%tstart
!end program
