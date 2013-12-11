module MetFields_ml

  implicit none
  private


!----------------- basic met fields ----------------------------------!
!  Here we declare the meteorological fields used in the model        !
!
! Horizonal alignments....
! Placement of q(i,j), u(i,j), v(i,j) 
!
!    --------------- --------------- --------------- ---------------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u03   q13      u13    q23      u23     q33     u33    q43      u43 ... u(LIMAX,3)
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v12------- -----v22------- ------v32------ -----v42-------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u02   q12      u12    q22      u22     q32     u32    q42      u42
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v11------- -----v21------- ------v31------ -----v41-------
!    |              |               |               |               | 
!    |              |               |               |               | 
!   u01   q11      u11    q21      u21     q31     u31    q41      u41
!    |              |               |               |               | 
!    |              |               |               |               | 
!    -----v10------- -----v20------- ------v30------ -----v40-------
!
  !---------------------------------------------------------------------!
  !
  !
  ! Vertical levels: z_mid,  z_bnd, sigma_mid, sigma_bnd
  !=============================================================================
  !*   "mid" and "bnd" are used as suffixes on z and sigma as shown in
  !*   the sketch below. "bnd" is the boundary between two layers and
  !*   "mid" the midddle of the layer. The numbering of layers starts
  !*   from 1 at the surface.
  !*
  !*
  !*
  !* ---------------------------
  !*
  !*
  !* - - - - - - - - - - - -     KMAX_MID -1
  !*
  !*
  !* --------------------------  KMAX_BDN-1       (z_bnd)   (sigma_bnd)
  !*
  !*
  !* - - - - - - - - -           KMAX_MID(old kmax2) = 20    (z_mid)   (sigma_mid)   (old z2)
  !*
  !* ------------------------    KMAX_BND = 21    (z_bnd)                 (old z1)
  !* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  surface \\\\\\\\\\\\\\\\
  !*
  ! RESULT FROM ONE TEST
  !Tested zm3d (=zm) which is sometimes used against
  !z_mid. Seems almost identical. Diff in exner functions maybe?
  !                   z_bnd        z_mid       zm3d 
  !     DEBUG_Z  2  1  16276.2359  15244.9023  15244.7131
  !     DEBUG_Z  2  2  14319.8682  13536.3354  13531.2981
  !     DEBUG_Z  2  3  12815.7707  12185.9090  12177.9740
  !     DEBUG_Z  2  4  11598.2202  11011.4342  11004.0286
  !     DEBUG_Z  2  5  10461.1309   9792.2354   9788.6193
  !     DEBUG_Z  2  6   9168.4701   8431.0622   8430.8546
  !     DEBUG_Z  2  7   7747.1534   7010.9878   7013.5860
  !     DEBUG_Z  2  8   6324.9278   5696.2857   5700.0226
  !     DEBUG_Z  2  9   5103.2253   4595.6436   4598.2893
  !     DEBUG_Z  2 10   4110.6562   3690.7686   3692.3391
  !     DEBUG_Z  2 11   3286.0770   2933.6493   2934.3594
  !     DEBUG_Z  2 12   2591.8972   2296.1632   2296.2730
  !     DEBUG_Z  2 13   2007.7924   1761.4943   1761.6936
  !     DEBUG_Z  2 14   1520.4434   1316.9528   1317.4565
  !     DEBUG_Z  2 15   1116.9477    950.9995    951.4678
  !     DEBUG_Z  2 16    787.3321    655.5900    655.8842
  !     DEBUG_Z  2 17    525.1633    424.5627    424.7443
  !     DEBUG_Z  2 18    324.8232    253.9725    254.0498
  !     DEBUG_Z  2 19    183.6019    137.4135    137.4284
  !     DEBUG_Z  2 20     91.3017     45.6234     45.6146


  !   Vertical level geopotential heights:

  real,public, save, allocatable,&
       dimension(:,:,:) :: z_bnd ! height of full layers
  real,public, save,allocatable, &
       dimension(:,:,:) :: z_mid ! height of half layers

  !   Two sets of Met. fields are read in, and a linear interpolation is made
  !   between these two points in time. NMET  == 2 (two points in time)
  !   note u_xmj, v_xmi are not "real" m/s wind speeds
  !   - they are actually divided by the mapping factor in the perpendicular direction).
  !
  real,target,public, save,allocatable, dimension(:,:,:,:) :: u_xmj 
  real,target,public, save,allocatable, dimension(:,:,:,:) :: v_xmi

  real,target,allocatable,public, dimension(:,:,:,:) :: q

  real,target,public, save,allocatable, dimension(:,:,:,:) :: &
       th      &  ! Potential teperature  ( deg. k )
!       ,q      &  ! Specific humidity
       ,roa    &  ! kg/m3
       ,cw        ! cloudwater
  real,target,public, save,allocatable, dimension(:,:,:,:) :: &
        EtaKz    &! vertical diffusivity in Eta coords
       ,SigmaKz  &! vertical diffusivity in sigma coords
       ,Etadot     &! vertical velocity, Eta coords, Pa/s
       ,sdot     &! vertical velocity, sigma coords, 1/s
       ,Kz_met    ! vertical diffusivity in sigma coordinates from meteorology


  ! since pr,cc3d,cc3dmax,cnvuf,cnvdf used only for 1 time layer - define without NMET
  real,target,public, save,allocatable, dimension(:,:,:) :: &
        pr      & ! Precipitation
       ,cc3d    & ! 3-d cloud cover (cc3d),
       ,cc3dmax & ! and maximum for layers above a given layer
       ,lwc     & !liquid water content
  ! QUERY - should xksig be MID, not BND? Is it needed at all?
       ,Kz_m2s     ! estimated Kz, in intermediate sigma levels, m2/s

  real,target,public, save,allocatable, dimension(:,:,:) :: &
        cnvuf   & ! convective_updraft_flux (kg/s/m2)
       ,cnvdf    ! convective_downdraft_flux (kg/s/m2)


 ! We don't need to calculate u,v for RiB, Kz for all layer in future maybe
 ! Still, for safety  we let this extent to K=1 for now

  real,target,public, save,allocatable, dimension(:,:,:) :: &
        u_mid   & ! wind u-compnent, m/s (real, not projected)
       ,v_mid     ! wind v-compnent, m/s
  

 real,target,public,save,allocatable, dimension(:,:,:) :: &
       tau        ! surf. stress  N/m^2

! Surface fields, interpolated:
 real,target,public, save,allocatable, dimension(:,:,:) :: &
        ps        &! Surface pressure Pa
       ,t2_nwp    & ! Temp 2 m   deg. K
       ,fh        & ! surf.flux.sens.heat W/m^2
       ,fl        & ! latent heat flux W/m^2
!       ,tau       & ! surf. stress  N/m^2
  ! These fields only available for EMEP/PARLAM from 2002 on
       ,rh2m            & !  RH at 2m
       ,SoilWater_uppr  & !  Shallow  (Upper 7.2cm in PARLAM)
       ,SoilWater_deep  & !  Deep (Next 6x7cm in PARLAM), converted to relative value 
       ,sdepth          & !  Snowdepth, m
       ,ice_nwp             & ! QUERY why real?
       ,sst     &  ! SST Sea Surface Temprature- ONLY from 2002 in PARLAM
       ,ws_10m    ! wind speed 10m
 

 real,target,public, save,allocatable, dimension(:,:) :: &
     u_ref             & ! wind speed m/s at 45m (real, not projected)
    ,rho_surf          & ! Surface density
    ,surface_precip    & ! Surface precip mm/hr
    ,Tpot2m            & ! Potential temp at 2m
    ,ustar_nwp         & ! friction velocity m/s ustar^2 = tau/roa
    ,invL_nwp          & ! friction velocity m/s ustar^2 = tau/roa
    ,pzpbl             & ! stores H(ABL) for averaging and plotting purposes, m
    ,pwp               & ! Permanent Wilting Point
    ,fc                  ! Field Capacity

!  temporary placement of solar radiation variations QUERY?
 
  real,target, public,allocatable, dimension(:,:), save:: &
       zen          &  ! Zenith angle (degrees)
      ,coszen       &  ! cos of zenith angle
      ,Idiffuse     &  ! diffuse solar radiation (W/m^2)
      ,Idirect         ! total direct solar radiation (W/m^2)



  real,target,public, save,allocatable, dimension(:,:) :: &   !st-dust
       clay_frac  &  ! clay fraction (%) in the soil
      ,sand_frac     ! sand fraction (%) in the soil

  ! Different NWP outputs for soil water are possible. We can currently
  ! cope with two:
  character(len=10), public, save  :: SoilWaterSource  ! IFS or PARLAM

  real,target,public, save, allocatable,dimension(:,:) :: &
    fSW     ! fSW= f(relative extractable water) =  (sw-swmin)/(swFC-swmin)

  real,target, public, dimension(:,:), save,allocatable  ::&
         xwf  ! extension of water fraction, save after 1st call

  integer, parameter, public :: NEXTEND = 2 ! no. box to side of (i,j) 

  integer, public, save   :: Nhh &         ! number of field stored per 24 hours
       ,nhour_first  ! time of the first meteo stored
! Logical flags, used to determine if some met fields are present in the
! input or not:
  logical, public, save :: &
     foundustar     & ! Used for MM5-type, where u_xmj* but not tau
    ,foundsdot      & ! If not found: compute using divergence=0
    ,sdot_at_mid    & ! set false if sdot is defined
    ,foundSST       & ! false if no SeaSurfaceT in metdata
    ,foundSoilWater_uppr  & ! false if no SW-shallow
    ,foundSoilWater_deep  & ! false if no SW-deep
    ,foundsdepth    & ! false if no snow_flag depth in metdata
    ,foundice       & ! false if no ice_nwp coverage (%) in metdata
    ,foundKz_met    & ! false if no Kz from meteorology
    ,foundconv      & ! false if convection not found or not used
  ! Introduced for FUTURE NH3, but also sea-salt
    ,foundws10_met   & ! false if no u10 from meteorology
    ,foundu10_met   & ! false if no u10 from meteorology
    ,foundv10_met   & ! false if no v10 from meteorology
    ,foundprecip    & ! false if no precipitationfrom meteorology
    ,foundcloudwater& !false if no cloudwater found
    ,foundSMI1& ! false if no Soil Moisture Index level 1 (shallow)
    ,foundSMI3 ! false if no Soil Moisture Index level 3 (deep)

  type,  public :: metfield
     character(len = 100) :: name = 'empty' !name as defined in external meteo file
     integer :: dim = 3 !number of dimension (2 for 2D, 3 for 3D)
     integer :: frequency =3  ! How many hours between two fields
     logical :: time_interpolate = .true. ! Interpolate in time  
     logical :: read_meteo = .false. ! The values will be looked for in the external meteo file
     logical :: found= .false. ! The values have been found in the external meteo file
!note that it is not allowed in fortran to define a target in a derived type
     real, pointer :: field(:,:,:,:) => null() !actual values for the fields; must be pointed to
     integer :: zsize = 1 ! field, size of third index
     integer :: msize = 1 ! field, size of fourth index
  endtype metfield

  integer, public, parameter   :: NmetfieldsMax=100 !maxnumber of metfields
  type(metfield),  public :: met(NmetfieldsMax)!size can be larger
  integer, public, save   :: Nmetfields! number of fields defined in met
  integer, public, save   :: N3Dmetfields! number of 3D fields defined in met
  real,target, public,save,allocatable, dimension(:,:,:) :: uw,ue
  real,target, public,save,allocatable, dimension(:,:,:) :: vs,vn

  public :: Alloc_MetFields !allocate arrays

contains

subroutine Alloc_MetFields(MAXLIMAX,MAXLJMAX,KMAX_MID,KMAX_BND,NMET)
!allocate MetFields arrays arrays
  implicit none
  
  integer, intent(in) ::MAXLIMAX,MAXLJMAX,KMAX_MID,KMAX_BND,NMET
  integer ::ix

  ix=1
  met(ix)%name             = 'u_wind'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(u_xmj(0:MAXLIMAX,1:MAXLJMAX,KMAX_MID,NMET))
  met(ix)%field(0:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:NMET)  => u_xmj
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'v_wind'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(v_xmi(1:MAXLIMAX,0:MAXLJMAX,KMAX_MID,NMET))
  met(ix)%field(1:MAXLIMAX,0:MAXLJMAX,1:KMAX_MID,1:NMET)  => v_xmi
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'specific_humidity'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(q(MAXLIMAX,MAXLJMAX,KMAX_MID,NMET))
  met(ix)%field => q 
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'potential_temperature'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(th(MAXLIMAX,MAXLJMAX,KMAX_MID,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:NMET)  => th
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'sigma_dot'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(sdot(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_BND,1:NMET)  => sdot
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = '3D_cloudcover'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(cc3d(MAXLIMAX,MAXLJMAX,KMAX_MID))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:1)  => cc3d
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'precipitation'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(pr(MAXLIMAX,MAXLJMAX,KMAX_MID))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:1)  => pr
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'cloudwater'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(cw(MAXLIMAX,MAXLJMAX,KMAX_MID,NMET))!should not be NMET?
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:NMET)  => cw
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'convective_updraft_flux'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(cnvuf(MAXLIMAX,MAXLJMAX,KMAX_BND))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_BND,1:1)  => cnvuf
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'convective_downdraft_flux'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(cnvdf(MAXLIMAX,MAXLJMAX,KMAX_BND))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_BND,1:1)  => cnvdf
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'eddy_diffusion_coefficient'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(Kz_met(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_BND,1:NMET)  => Kz_met
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'air_density'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(roa(MAXLIMAX,MAXLJMAX,KMAX_MID,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:NMET)  => roa
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'Kz_sigmacoordinates'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(SigmaKz(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_BND,1:NMET)  => SigmaKz
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'Kz_Etacoordinates'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(EtaKz(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_BND,1:NMET)  => EtaKz
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET


  ix=ix+1
  met(ix)%name             = 'etadot'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(Etadot(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_BND,1:NMET)  => Etadot
  met(ix)%zsize = KMAX_BND
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'max_cloudcover'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(cc3dmax(MAXLIMAX,MAXLJMAX,KMAX_MID))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:1)  => cc3dmax
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'cloud_liquid_water'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(lwc(MAXLIMAX,MAXLJMAX,KMAX_MID))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:1)  => lwc
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'Kz'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(Kz_m2s(MAXLIMAX,MAXLJMAX,KMAX_MID))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:1)  => Kz_m2s
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'u_wind_3D'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(u_mid(MAXLIMAX,MAXLJMAX,KMAX_MID))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:1)  => u_mid
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'v_wind_3D'
  met(ix)%dim              = 3
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(v_mid(MAXLIMAX,MAXLJMAX,KMAX_MID))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:KMAX_MID,1:1)  => v_mid
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = 1

  N3Dmetfields=ix


  ix=ix+1
  met(ix)%name             = 'surface_pressure'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(ps(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => ps
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'temperature_2m'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(t2_nwp(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => t2_nwp
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'relative_humidity_2m'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(rh2m(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => rh2m
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'surface_flux_sensible_heat'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(fh(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => fh
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'surface_flux_latent_heat'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(fl(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => fl
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'surface_stress'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(tau(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => tau
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'ustar_nwp'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(ustar_nwp(MAXLIMAX,MAXLJMAX))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:1)  => ustar_nwp
  met(ix)%zsize = 1
  met(ix)%msize = 1


  ix=ix+1
  met(ix)%name             = 'sea_surface_temperature'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(sst(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => sst
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'SMI1'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(SoilWater_uppr(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => SoilWater_uppr
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'SMI3'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(SoilWater_deep(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => SoilWater_deep
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'snow_depth'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(sdepth(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => sdepth
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'fraction_of_ice'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(ice_nwp(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => ice_nwp
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'u10'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(ws_10m(MAXLIMAX,MAXLJMAX,NMET))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:NMET)  => ws_10m
  met(ix)%zsize = 1
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'large_scale_precipitations'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .false.
  met(ix)%read_meteo       = .true.
  met(ix)%found            = .false.
  allocate(surface_precip(MAXLIMAX,MAXLJMAX))
  met(ix)%field(1:MAXLIMAX,1:MAXLJMAX,1:1,1:1)  => surface_precip
  met(ix)%zsize = 1
  met(ix)%msize = 1

  ix=ix+1
  met(ix)%name             = 'neigbors_wind-uw'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(uw(MAXLJMAX,KMAX_MID,NMET))
  met(ix)%field(1:1,1:MAXLJMAX,1:KMAX_MID,1:NMET)  => uw
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'neigbors_wind-ue'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(ue(MAXLJMAX,KMAX_MID,NMET))
  met(ix)%field(1:1,1:MAXLJMAX,1:KMAX_MID,1:NMET)  => ue
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'neigbors_wind-vs'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(vs(MAXLIMAX,KMAX_MID,NMET))
  met(ix)%field(1:MAXLIMAX,1:1,1:KMAX_MID,1:NMET)  => vs
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  ix=ix+1
  met(ix)%name             = 'neigbors_wind-vn'
  met(ix)%dim              = 2
  met(ix)%frequency        = 3
  met(ix)%time_interpolate = .true.
  met(ix)%read_meteo       = .false.
  met(ix)%found            = .false.
  allocate(vn(MAXLIMAX,KMAX_MID,NMET))
  met(ix)%field(1:MAXLIMAX,1:1,1:KMAX_MID,1:NMET)  => vn
  met(ix)%zsize = KMAX_MID
  met(ix)%msize = NMET

  Nmetfields=ix
  if(Nmetfields>NmetfieldsMax)then
     write(*,*)"Increase NmetfieldsMax! "
     stop
  endif

    allocate(u_ref(MAXLIMAX,MAXLJMAX))
    allocate(rho_surf(MAXLIMAX,MAXLJMAX))
    allocate(Tpot2m(MAXLIMAX,MAXLJMAX))
    allocate(invL_nwp(MAXLIMAX,MAXLJMAX))
    allocate(pzpbl(MAXLIMAX,MAXLJMAX))
    allocate(pwp(MAXLIMAX,MAXLJMAX))
    allocate(fc(MAXLIMAX,MAXLJMAX))
    allocate(xwf(MAXLIMAX+2*NEXTEND,MAXLJMAX+2*NEXTEND)) 
    allocate(fSW(MAXLIMAX,MAXLJMAX))
    fSW = 1.0
    allocate(zen(MAXLIMAX, MAXLJMAX))
    allocate(coszen(MAXLIMAX, MAXLJMAX))
    coszen=0.0
    allocate(Idiffuse(MAXLIMAX, MAXLJMAX))
    allocate(Idirect(MAXLIMAX, MAXLJMAX))
    allocate(clay_frac(MAXLIMAX, MAXLJMAX))
    allocate(sand_frac(MAXLIMAX, MAXLJMAX))


  end subroutine Alloc_MetFields

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


end module MetFields_ml
