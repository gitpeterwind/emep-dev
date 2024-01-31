!> TEMPORARY MODULE for consistency with ecosx
!  Will move all DEBUG from Config_module here one day
module Debug_module

  ! DebugCell is set in Solver_mod when DEBUG%RUNCHEM is true, i=debug_li,
  ! j=debug_lj,  and k==20. Allows extra debug for one cell

 logical, public, save ::  DebugCell  = .false.

 type, public :: emep_debug
  logical :: &
     AOT             = .false. &
    ,ADV             = .false. & ! 
    ,AEROSOL         = .false. & ! ...needed for intended debugs are to work
    ,AQUEOUS         = .false. &
    ,BCS             = .false. & ! BoundaryConditions
    ,BIO             = .false. & ! Biogenic emissions
    ,BIDIR           = .false. & ! FUTURE Bi-directional exchange
    ,BLM             = .false. & ! Produces matrix of differnt Kz and Hmix
    ,COLUMN          = .false. & ! Used in Derived_mod for column integration
    ,COLSRC          = .false. & ! Volcanic emissions and Emergency scenarios
    ,DERIVED         = .false. & !
    ,DRYDEP          = .false. & ! 
    ,DRYRUN          = .false. & ! Skips fast chemistry to save some CPU
    ,DUST            = .false. & ! Skips fast chemistry to save some CPU
    ,ECOSYSTEMS      = .false. &
    ,EMISSIONS       = .false. & ! 
    ,EMISSTACKS      = .false. & ! 
    ,EMISTIMEFACS    = .false. &
    ,EQUIB           = .false. &   !MARS, EQSAM etc.
    ,FORESTFIRE      = .false. &
    ,GETEMIS         = .false. &
    ,GLOBBC          = .false. &
    ,GRIDVALUES      = .false. &
    ,HOURLY_OUTPUTS  = .false. & !
    ,IOPROG          = .false. &
    ,Kz              = .false. &
    ,LANDDEFS        = .false. &
    ,LANDIFY         = .false. &
    ,MAINCODE        = .false. & !< debugs main code (emepctm) driver
    ,MASS            = .false. &
    ,MET             = .false. &
    ,MOSAICS         = .false. &
    ,MY_DERIVED      = .false. &
    ,NEST            = .false. &
    ,NEST_ICBC       = .false. & ! IFS-MOZART/C-IFS BC
    ,NETCDF          = .false. &
    ,NETCDF_RF       = .false. & ! ReadField_CDF in NetCDF_mod
    ,OUTPUTCHEM      = .false. & ! Output of netcdf results
    ,pH              = .false. &
    ,PHYCHEM         = .false. &
    ,POLLEN          = .false. &
    ,ROADDUST        = .false. &
    ,RSUR            = .false. & ! Surface resistance
    ,RUNCHEM         = .false. & ! DEBUG%RUNCHEM is SPECIAL, need for some other debugs
    ,MY_WETDEP    = .false. &
    ,SEASALT         = .false. &
    ,SETUP_1DCHEM    = .false. &
    ,SETUP_1DBIO     = .false. &
    ,SITES           = .false. & ! set also DEBUG%SITE below
    ,SOILNOX         = .false. &
    ,SOILWATER       = .false. &
    ,SOLVER          = .false. &
    ,STOFLUX         = .false. &
    ,VERT_DIFF       = .false. &
    ,VDS             = .false. &
    ,PBAP            = .false. &
    ,FUNGAL_SPORES   = .false. &
    ,BACTERIA        = .false.
  ! integer debug options allow different levels of verbosity
   integer               :: &
      PFT_MAPS  = 0         & !< Future option
     ,LANDUSE   = 0         & !
     ,DO3SE     = 0         & !
     ,SOA       = 0         &
     ,SUBMET    = 2         & ! Use 999 for all land-cover, otherwise LC index
     ,STOP_HH   = -1          ! If positive, code will quite when hh==STOP_HH
  !----------------------------------------------------------
   integer, dimension(2) :: IJ = [-999,-999]  ! index for debugging print out
   character(len=20)     :: SPEC = 'O3'       ! default.
   character(len=20)     :: datetxt = '-'     ! default.
   character(len=20)     :: SITE = 'NOT_SET'  !  e.g. Birkenes. (Full name not essential)
   integer               :: ISPEC = -999      ! Will be set after NML
end type emep_debug
type(emep_debug), public, save :: DEBUG


end module Debug_module
