!*****************************************************************************!
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#define IF_MPI_NOT_OK_RETURN(action) if (status/=MPI_SUCCESS) then; call MPI_Error_String(status,gol,ngol,status); call goErr; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
!*****************************************************************************!

module DA_Obs_ml

  use DA_Util_ml,       only: ngol, gol, goPr, goErr
  use MPI   ,           only: MPI_SUCCESS, MPI_Error_String
  use NetCDF,           only: NF90_NOERR, NF90_StrError
  use TimeDate_ml,      only: current_date
  use CheckStop_ml,     only: CheckStop
  use Functions_ml,     only: great_circle_distance
  use ChemSpecs_adv_ml, only: NSPEC_ADV
  use ChemChemicals_ml, only: species_adv
  use ChemGroups_ml,    only: chemgroups
  use Io_Progs_ml,      only: PrintLog
  use ModelConstants_ml,only: MasterProc,runlabel1
  use SmallUtils_ml,    only: find_index
  use TimeDate_ExtraUtil_ml, only: date2string
  use DA_ml           , only : debug => DEBUG_DA
  use DA_ml           , only : dafmt => da_fmt_msg
  use DA_ml           , only : damsg => da_msg
  use DA_ml           , only : nlev

  implicit none

  
  ! --- const -----------------------------------------

  character(len=*), parameter  ::  mname = 'DA_Obs_ml'
  
  ! max length:
  integer, parameter  ::  LEN_LABEL = 64
  integer, parameter  ::  LEN_SCODE = 8
  
  ! level types:                               extract from ...
  integer, parameter  ::  LEVTYPE_3D_ML      = 1  ! 3D model layers
  integer, parameter  ::  LEVTYPE_3D_ML_SFC  = 2  ! 3D model layers, surface simulation (3m)
  integer, parameter  ::  LEVTYPE_3D_ML_TC   = 3  ! 3D model layers, sum to total column
  !integer, parameter  ::  LEVTYPE_3D_ML_K1   = 4  ! 3D model layers, apply kernel for column
  integer, parameter  ::  LEVTYPE_3D_ML_K2   = 4  ! 3D model layers, apply kernel for column, trop
  integer, parameter  ::  LEVTYPE_2D_ML_SFC  = 5  ! 2D lowest model layer, surface simulation (3m)
  integer, parameter  ::  LEVTYPE_2D_OBS_SFC = 6  ! 2D simulated observations, surface simulation (3m)
  ! counter:
  integer, parameter  ::  NLEVTYPE = 6
  ! names:
  character(len=10), parameter   ::  LEVTYPE_NAME(NLEVTYPE) = (/ &
                                        '3D_ML     ', &
                                        '3D_ML_SFC ', &
                                        '3D_ML_TC  ', &
                                        '3D_ML_K2  ', &
                                        '2D_ML_SFC ', &
                                        '2D_OBS_SFC' /)
  
  
  ! --- in/out ---------------------------------------
  
  private

  public    ::  DA_Obs_Init, DA_Obs_Done
  
  public    ::  Read_Obs
  
  public    ::  T_ObsOper
  public    ::  T_ObsOpers
  
  public    ::  nObsData, obsData
  public    ::  nObsComp, ObsCompInfo
  
  public    ::  LEVTYPE_3D_ML
  public    ::  LEVTYPE_3D_ML_SFC
  public    ::  LEVTYPE_3D_ML_TC
  !public    ::  LEVTYPE_3D_ML_K1
  public    ::  LEVTYPE_3D_ML_K2
  public    ::  LEVTYPE_2D_ML_SFC
  public    ::  LEVTYPE_2D_OBS_SFC
  
  !public    ::  dbg_cell, dbg_i, dbg_j


  ! --- types ----------------------------------------
  
  ! collection of observation errors per station code
  type T_ObsError_List
    ! number of stations:
    integer                  ::  n
    ! station codes:
    character(len=LEN_SCODE), allocatable   ::  scode(:)   ! (n)
    ! contributions to total error:
    integer                                 ::  nval
    real, allocatable                       ::  err(:,:)   ! (n,nval)
  contains
    procedure   ::  Init         => ObsError_List_Init
    procedure   ::  Done         => ObsError_List_Done
    procedure   ::  ReadFiles    => ObsError_List_ReadFiles
    procedure   ::  GetError     => ObsError_List_GetError
  end type T_ObsError_List
  
  ! *
  
  ! potential storage for kernels:
  type T_ObsKernels
    ! dims:
    integer                   ::  npixel
    integer                   ::  nlayer
    ! kernels and pressures, only allocated if necessary:
    real, allocatable         ::  kernel(:,:)  ! (layer,pixel)
    character(len=32)         ::  kernel_units
    real, allocatable         ::  phlev(:,:)   ! (layer_interface,pixel)
    character(len=32)         ::  phlev_units
  contains
    procedure   ::  Init         => ObsKernels_Init
    procedure   ::  Done         => ObsKernels_Done
    procedure   ::  Clear        => ObsKernels_Clear
    procedure   ::  Fill         => ObsKernels_Fill
  end type T_ObsKernels
  
  ! *
  
  ! observation data collection:
  type obs_data
    !
    ! * elements defined in configuration namelist
    !
    ! file types:
    !   'txt'  : text files used for operational processing
    !   'omi'  : OMI netcdf files (converted from HDFEOS5)
    character(len=32)         ::  ftype      = 'txt'
    ! observation file (analyzed):
    character(len=256)        ::  file       = ''
    ! observation files for validation:
    character(len=256)        ::  file_val   = ''
    !
    ! observed component:
    character(len=LEN_LABEL)  ::  name       = ''    ! 'O3', 'NO2', 'PM10', ...
    ! keyword for how to compute from 3D model fields:
    character(len=LEN_LABEL)  ::  deriv      = ''
    !
    ! observation data units:
    character(len=LEN_LABEL)  ::  unit       = ''
    ! allowed range:
    real                      ::  min        = 0.0
    real                      ::  max        = 1e3
    !
    ! obs.repr. type: 'fraction' | 'estim'
    character(len=LEN_LABEL)  ::  error_type = ''
    !~ fraction and maximum:
    real                      ::  error_rel  = -1.0e0
    real                      ::  error_rep  = -1.0e0
    !~ estimated: files with obs.error and spatial error per station code:
    character(len=1024)       ::  error_obs  = ''
    character(len=1024)       ::  error_spat = ''
    !
    ! background covar:
    character(len=1024)       ::  Bfile      = ''
    character(len=LEN_LABEL)  ::  Bunit      = ''
    ! scaling factor for sigma:
    real                      ::  Bsigfac    = 1.0
    ! feedback analysis to model?
    logical                   ::  feedback      = .true.
    !~ how to feedback?
    character(len=LEN_LABEL)  ::  feedback_type = ''
    !
    ! * auxilary
    !
    ! tag filled for output:
    character(len=256)        ::  tag        = ''
    ! used for mapping from model tracers:
    integer                   ::  iObsComp   = -999
    ! some index range in txt input:
    integer                   ::  iobs(0:1)  = -1
  end type

  !-----------------------------------------------------------------------
  ! H_jac is a matrix with components H_{n;i,j,l,k}, where
  !   n=1,...,nobs     (observation index)
  !   i=1,...,nx       (longitude index)
  !   j=1,...,ny       (latitude index)
  !   l=1,...,nlev     (altitude index)
  !   k=1,...,nchemobs (observed chemical component)
  ! Due to the 4-point interpolation in the horizontal direction
  ! of the model results to the observation point, the Jacobian is
  ! highly sparse wrt i,j,l. For this reason, the matrix is stored using
  ! the following compact storage scheme:
  ! H_{n;i,j,l,k}=H_jac(n,p,l,k), where n=1:nobs,p=1:4,l=l(0):l(1)
  !   pidx(n)%i(1)=i_n,   pidx(n)%j(1)=j_n,   pidx(n)%l(0)=l0_n,
  !   pidx(n)%i(2)=i_n+1, pidx(n)%j(2)=j_n,   pidx(n)%l(1)=l1_n,
  !   pidx(n)%i(3)=i_n,   pidx(n)%j(3)=j_n+1,
  !   pidx(n)%i(4)=i_n+1, pidx(n)%j(4)=j_n+1,
  !-----------------------------------------------------------------------

  ! observation operator 
  ! NOTE: when changing this, also adapt:
  !  - ObsOper_Copy
  !  - ObsOpers_Swap
  type T_ObsOper
    ! observation and std.dev.:
    real                        ::  obs
    real                        ::  obsstddev
    ! model simulations:
    real                        ::  xf       ! forecast
    real                        ::  xa       ! analysis
    real                        ::  sb       ! background sigma
    real                        ::  sf       ! forecast sigma
    ! innovation:
    real                        ::  innov    ! obs - forecast
    ! Jacobian dH/dx :
    real, allocatable           ::  H_jac(:)  ! (nlev)
    ! TESTING: scale factor for profile
    real                        ::  H_jac1
    ! grid cell:
    integer                     ::  i
    integer                     ::  j
    ! vertical:
    integer                     ::  levtype
    integer                     ::  l(0:1)
    ! id for testing:
    integer                     ::  id
    ! observed component:
    integer                     ::  iObsComp
    ! index in 'obsData' array:
    integer                     ::  iObsData
    ! meta data:
    integer                     ::  stnid
    real                        ::  lon, lat
    real                        ::  alt
    character(len=LEN_SCODE)    ::  stncode
    logical                     ::  analyse
  contains
    procedure   ::  Init         => ObsOper_Init
    procedure   ::  Done         => ObsOper_Done
    procedure   ::  Copy         => ObsOper_Copy
    procedure   ::  Fill         => ObsOper_Fill
    procedure   ::  Evaluate     => ObsOper_Evaluate
  end type T_ObsOper

  ! observation operators:
  type T_ObsOpers
    ! number of observations:
    integer                         ::  nobs
    ! info per observations:
    type(T_ObsOper), allocatable    ::  obs(:)
    ! scale factors from online tuning:
    real, allocatable               ::  BScale(:)  ! (nObsComp)
  contains
    procedure   ::  Init          => ObsOpers_Init
    procedure   ::  Done          => ObsOpers_Done
    procedure   ::  Alloc         => ObsOpers_Alloc
    procedure   ::  DeAlloc       => ObsOpers_DeAlloc
    procedure   ::  Set_IDs       => ObsOpers_Set_IDs
    procedure   ::  Swap          => ObsOpers_Swap
    procedure   ::  SelectTracers => ObsOpers_SelectTracers
    procedure   ::  Evaluate      => ObsOpers_Evaluate
    procedure   ::  SetBScale     => ObsOpers_SetBScale
    procedure   ::  WriteToFile   => ObsOpers_WriteToFile
  end type T_ObsOpers

  ! *

  ! observed component is linear combination of advected species:
  type T_ObsCompInfo
    ! component name and units:
    character(len=32)                ::  name
    character(len=32)                ::  units
    ! vertical: 
    !   'mod-lev'   :  extract from (range of) model level(s)
    !   'surface'   :  from bottom model layer to surface using 'cfac'
    character(len=32)                ::  deriv
    integer                          ::  levtype
    ! covariance file:
    character(len=1024)              ::  Bfile
    ! adhoc scale factor for sigma:
    real                             ::  Bsigfac
    ! number of model species:
    integer                          ::  nspec
    ! indices of model advected species:
    integer, allocatable             ::  ispec(:)     ! (nspec)
    ! weights in sum:
    real, allocatable                ::  w(:)         ! (nspec)
    ! unit conversion:
    real, allocatable                ::  unitconv(:)  ! (nspec)
    logical, allocatable             ::  unitroa(:)
    ! add coarse NO3 aerosol fraction (assigned to PM25) ?
    logical                          ::  with_no3c_frac
    ! add aerosol water?
    logical                          ::  with_pmwater
    ! feed back into model?
    logical                          ::  feedback
    !~ how to feed back?
    character(len=32)                ::  feedback_type = ''
  contains
    procedure :: Init                         => ObsCompInfo_Init
    procedure :: Done                         => ObsCompInfo_Done
    procedure :: Show                         => ObsCompInfo_Show
    procedure :: FillFields                   => ObsCompInfo_FillFields
    procedure :: DistributeIncrement          => ObsCompInfo_DistributeIncrement_3d
    procedure :: DistributeIncrement1         => ObsCompInfo_DistributeIncrement1_3d
    procedure :: ChangeProfile                => ObsCompInfo_ChangeProfile
    !procedure :: ChangeProfileFineCoarseRatio => ObsCompInfo_ChangeProfile_FineCoarseRatio
    procedure :: AdjustFineCoarseRatio        => ObsCompInfo_AdjustFineCoarseRatio
    procedure :: ChangeML                     => ObsCompInfo_ChangeML
    procedure :: ChangeVCD                    => ObsCompInfo_ChangeVCD
    !
  end type T_ObsCompInfo


  !--- local -------------------------------------------------------------

  ! maximum number of data sets:
  integer, parameter :: nObsDataMax = 10
  ! storage for dataset definitions:
  type(obs_data)     :: obsData(nObsDataMax)
  ! actual number:
  integer            :: nObsData=0

  ! storage for error lists:
  type(T_ObsError_List), allocatable    ::  ObsError_Lists(:)   ! (nObsData)
  
  ! storage for kernels:
  type(T_ObsKernels), allocatable       ::  ObsKernels(:)  ! (nObsData)

  
  ! number of observed components:
  integer                            ::  nObsComp
  ! info:
  type(T_ObsCompInfo), allocatable   ::  ObsCompInfo(:)  ! (nObsComp)
  
  !! testing ...
  !logical            ::  dbg_cell = .false.
  !integer            ::  dbg_i    = -999
  !integer            ::  dbg_j    = -999
  

contains


  !=========================================================================
  !===
  !=== observed components
  !===
  !=========================================================================

  
  subroutine ObsCompInfo_Init( self, name, units, deriv, Bfile, indices, status, &
                                 with_no3c_frac, with_pmwater, feedback, feedback_type, Bsigfac )
  
    use Units_ml        , only : Units_Scale
    use SmallUtils_ml   , only : find_index
    use ChemChemicals_ml, only : species_adv

    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(out)   ::  self
    character(len=*), intent(in)        ::  name
    character(len=*), intent(in)        ::  units
    character(len=*), intent(in)        ::  deriv  ! 'mod-lev', 'surface'
    character(len=*), intent(in)        ::  Bfile
    integer, intent(in)                 ::  indices(:)
    integer, intent(out)                ::  status
    
    logical, intent(in), optional           ::  with_no3c_frac
    logical, intent(in), optional           ::  with_pmwater
    logical, intent(in), optional           ::  feedback
    character(len=*), intent(in), optional  ::  feedback_type
    real, intent(in), optional              ::  Bsigfac
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_Init'
    
    ! --- local -----------------------------
    
    integer       ::  i
    
    ! --- begin -----------------------------
    
    ! store name:
    self%name  = trim(name)
    self%units = trim(units)
    self%deriv = trim(deriv)
    self%Bfile = trim(Bfile)

    ! scale factor:
    self%Bsigfac = 1.0
    if ( present(Bsigfac) ) self%Bsigfac = Bsigfac

    ! add aerosol water?
    self%with_pmwater = .false.
    if ( present(with_pmwater) ) self%with_pmwater = with_pmwater
    
    ! feed back into model?
    self%feedback = .true.
    if ( present(feedback) ) self%feedback = feedback
    ! how to feed back:
    self%feedback_type = ''
    if ( present(feedback_type) ) self%feedback_type = feedback_type
    
    ! include fraction of coarse nitrate aerosol ? used for PM25:
    self%with_no3c_frac = .false.
    if ( present(with_no3c_frac) ) self%with_no3c_frac = with_no3c_frac

    ! number:
    self%nspec = size(indices)
    ! extra?
    if ( self%with_no3c_frac ) self%nspec = self%nspec + 1
    
    ! storage:
    allocate( self%ispec(self%nspec), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%w(self%nspec), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%unitconv(self%nspec), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%unitroa(self%nspec), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! init with equal weights:
    self%w = 1.0
    
    ! copy:
    self%ispec(1:size(indices)) = indices

    ! extra?
    if ( self%with_no3c_frac ) then
      ! target element:
      i = size(indices) + 1
      ! find index of coarse nitrate:
      self%ispec(i) = find_index( 'NO3_C', species_adv(:)%name )
      ! add 27% of weight:
      self%w(i) = 0.27
    end if
    
    ! loop over input species:
    do i = 1, self%nspec
      ! obtain unit conv factor to convert from ispec units to "units":
#ifdef with_ajs
      call Units_Scale( units, self%ispec(i), self%unitconv(i), &
                           needroa=self%unitroa(i), status=status )
      IF_NOT_OK_RETURN(status=1)
#else
      call Units_Scale( units, self%ispec(i), self%unitconv(i), &
                           needroa=self%unitroa(i) )
#endif
    end do
    
    ! ok
    status = 0
    
  end subroutine ObsCompInfo_Init


  ! ***
  
  
  subroutine ObsCompInfo_Done( self, status )
  
    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(inout)   ::  self
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_Done'
    
    ! --- local -----------------------------
    
    ! --- begin -----------------------------
    
    ! clear:
    deallocate( self%ispec, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%w, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%unitconv, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%unitroa, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine ObsCompInfo_Done


  ! ***
  
  
  subroutine ObsCompInfo_Show( self, status )

    use ChemChemicals_ml, only : species_adv
  
    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(in)      ::  self
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_Show'
    
    ! --- local -----------------------------
    
    integer               ::  i
    integer               ::  ispec
    character(len=6)      ::  rholab
    
    ! --- begin -----------------------------
    
    ! intro ...
    write (gol,'(a,": observed component : ",a)') rname, trim(self%name); call goPr
    write (gol,'(a,":   units            : ",a)') rname, trim(self%units); call goPr
    write (gol,'(a,":   derivation       : ",a)') rname, trim(self%deriv); call goPr
    write (gol,'(a,":   source species   : ",i0)') rname, self%nspec; call goPr
    ! loop:
    do i = 1, self%nspec
      ispec = self%ispec(i)
      ! multiply with density?
      rholab = '      '
      if ( self%unitroa(i) ) rholab = ' * rho'
      ! show:
      write (gol,'(a,":   ",i3," ",a20," * ",es12.4,a," * ",f8.2)') rname, &
               ispec, trim(species_adv(ispec)%name), self%unitconv(i), rholab, self%w(i); call goPr
    end do
    ! aerosol water?
    if ( self%with_pmwater ) then
      ! dummy:
      ispec = -99
      rholab = '      '
      ! extra line:
      write (gol,'(a,":   ",i3," ",a20,"               ",a," * ",f8.2)') rname, &
               ispec, 'PM25_water', rholab, 1.0; call goPr
    end if
    
    ! ok
    status = 0
    
  end subroutine ObsCompInfo_Show


  ! ***
  
  !
  ! Simulate observed component from model species in "xn_adv"
  ! and return this in arguments:
  !   xn_obs_sfc   : sfc field (or 2D field of columns)
  !   xn_obs_ml    : model level field
  !   xn_obs_units
  !
  
  subroutine ObsCompInfo_FillFields( self, xn_adv, xn_adv_units, &
                                       xn_obs_sfc, xn_obs_ml, xn_obs_units, status )
  
    use MetFields_ml , only : roa
    use ChemFields_ml, only : cfac
    use Chemfields_ml, only : PM25_water, PM25_water_rh50

    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(in)      ::  self
    real, intent(in)                      ::  xn_adv(:,:,:,:)    ! (nspec_adv,lnx,lny,nlev)
    character(len=*), intent(in)          ::  xn_adv_units(:)    ! (nspec_adv)
    real, intent(out)                     ::  xn_obs_sfc(:,:)    ! (lnx,lny)
    real, intent(out)                     ::  xn_obs_ml(:,:,:)   ! (lnx,lny,nlev)
    character(len=*), intent(out)         ::  xn_obs_units
    integer, intent(out)                  ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_FillFields'
    
    ! --- local -----------------------------
    
    integer             ::  i
    integer             ::  ispec
    integer             ::  k
    real, allocatable   ::  xn_obs_ml__in(:,:,:)  ! (lnx,lny,nlev)
    real, allocatable   ::  xn_adv1(:,:,:)    ! (lnx,lny,nlev)
    
    ! --- begin -----------------------------
    
    ! storage for single tracer field:
    allocate( xn_adv1(size(xn_obs_ml,1),size(xn_obs_ml,2),size(xn_obs_ml,3)), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! switch:
    select case ( trim(self%deriv) )
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! accumulated fields:
      case ( '2D-OBS-SFC' )
      
        ! units is the same as observations:
        xn_obs_units = trim(self%units)
        
        ! init sum for surface field:
        xn_obs_sfc = 0.0
        xn_obs_ml  = 0.0
        ! loop:
        do i = 1, self%nspec
          ! current:
          ispec = self%ispec(i)
          ! copy original field:
          xn_adv1 = xn_adv(ispec,:,:,:)
          ! include density conversion for this tracer?
          if ( self%unitroa(i) ) then
            ! surface from bottom layer:
            k = size(xn_adv,4)
            ! add weighted contribution:
            xn_obs_sfc = xn_obs_sfc + self%w(i) * xn_adv1(:,:,k) &
                                                * self%unitconv(i)    &  ! unit conversion
                                                * roa(:,:,k,1)        &  ! density
                                                * cfac(ispec,:,:)        ! 50 m -> 3 m
            ! idem for 3D field:
            xn_obs_ml  = xn_obs_ml  + self%w(i) * xn_adv1(:,:,:) &
                                                * self%unitconv(i)    &  ! unit conversion
                                                * roa(:,:,:,1)           ! density
          else
            ! surface from bottom layer:
            k = size(xn_adv,4)
            ! add weighted contribution:
            xn_obs_sfc = xn_obs_sfc + self%w(i) * xn_adv1(:,:,k) &
                                                * self%unitconv(i)    &  ! unit conversion
                                                * cfac(ispec,:,:)        ! 50 m -> 3 m
            ! idem for 3D field:
            xn_obs_ml  = xn_obs_ml  + self%w(i) * xn_adv1(:,:,:) &
                                                * self%unitconv(i)       ! unit conversion
          end if  ! density conversion
        end do  ! specs
        
        ! add water?
        if ( self%with_pmwater ) then
          ! check ..
          if ( self%units /= 'ug/m3' ) then
            write (gol,'("unexpected target units `",a,"` for adding pm water")') trim(self%units); call goErr
            TRACEBACK; status=1; return
          end if
          ! add:
          xn_obs_sfc = xn_obs_sfc + PM25_water_rh50
          xn_obs_ml  = xn_obs_ml  + PM25_water
        endif
        
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! For following deriv's, H will include factors etc to convert 
      ! from model units (mole/mole) to the observation unit and type:
      !  - 'unitconv' factor (for example incl molemass of air, kg to ug, etc)
      !  - eventually 'roa' (air density)
      !  - for surface: 'cfac' (due to deposition profile)
      !
      !case ( '3D-ML-SFC', '2D-ML-SFC', '3D-ML-TC', '3D-ML-K1', '3D-ML-K2' )
      case ( '3D-ML-SFC', '2D-ML-SFC', '3D-ML-TC', '3D-ML-K2' )
      
        ! check ...
        if ( self%nspec /= 1 ) then
          write (gol,'("observation operator on original model concentrations")'); call goErr
          write (gol,'("not implemented for accumulated tracer (",a,") yet ..")') &
                   trim(self%name); call goErr
          TRACEBACK; status=1; return
        end if
        ! current (and only) species:
        ispec = self%ispec(1)

        ! native model units:
        xn_obs_units = 'mol mol^-1'
        
        ! copy original field:
        xn_adv1 = xn_adv(ispec,:,:,:)

        ! surface from bottom layer:
        k = size(xn_adv,4)
        ! copy bottom layer:
        xn_obs_sfc = xn_adv1(:,:,k)
        
        ! copy model levels:
        xn_obs_ml = xn_adv1
        
        ! add water?
        if ( self%with_pmwater ) then
          write (gol,'("observation operator on original model concentrations")'); call goErr
          write (gol,'("not implemented for addition of pmwater to tracer (",a,") yet ..")') &
                   trim(self%name); call goErr
          TRACEBACK; status=1; return
        endif
        
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! not yet
      case default
        write (gol,'("unsupported observation deriv key : ",a)') trim(self%deriv); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! clear:
    deallocate( xn_adv1, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine ObsCompInfo_FillFields


  ! *


  !
  ! Change xn_adv array (all tracers) given observed field (total pm?) and increment.
  !
  
  subroutine ObsCompInfo_DistributeIncrement_3d( self, xn_adv, xn_adv_units, &
                                                   xn_obs, dx_obs, obs_units, status, &
                                                   maxratio, verbose )

    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(in)      ::  self
    real, intent(inout)                   ::  xn_adv(:,:,:,:)  ! (nspec_adv,lnx,lny,nz)
    character(len=*), intent(in)          ::  xn_adv_units(:)  ! (nspec_adv)
    real, intent(in)                      ::  xn_obs(:,:,:)    ! (lnx,lny,nz)
    real, intent(in)                      ::  dx_obs(:,:,:)    ! (lnx,lny,nz)
    character(len=*), intent(in)          ::  obs_units
    integer, intent(out)                  ::  status
    
    real, intent(in), optional            ::  maxratio
    logical, intent(in), optional         ::  verbose
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_DistributeIncrement_3d'
    
    ! --- local -----------------------------
    
    integer               ::  ispec_adv
    logical               ::  verb

    real                  ::  tpm
    integer               ::  i
    integer               ::  ispec
    
    ! --- begin -----------------------------
    
    ! verbose?
    verb = .false.
    if ( present(verbose) ) verb = verbose

    !! testing ...
    !if ( dbg_cell ) then
    !  !nlev = size(xn_adv,4)
    !  tpm = 0.0
    !  do i = 1, self%nspec
    !    ispec = self%ispec(i)
    !    tpm = tpm + xn_adv(ispec,dbg_i,dbg_j,nlev)
    !    write (gol,*) 'yyy before xn_adv = ', xn_adv(ispec,dbg_i,dbg_j,nlev); call goPr
    !  end do
    !  write (gol,*) 'yyy before tpm    = ', tpm; call goPr
    !end if

    ! apply per advected tracer:
    do ispec_adv = 1, size(xn_adv,1)
    
      ! apply per tracer:
      call self%DistributeIncrement1( ispec_adv, xn_adv(ispec_adv,:,:,:), xn_adv_units(ispec_adv), &
                                           xn_obs, dx_obs, obs_units, status, &
                                           maxratio=maxratio, verbose=verbose )
      IF_NOT_OK_RETURN(status=1)
      
    end do  ! ispec_adv

!    ! testing ...
!    if ( dbg_cell ) then
!      tpm = 0.0
!      do i = 1, self%nspec
!        ispec = self%ispec(i)
!        tpm = tpm + xn_adv(ispec,dbg_i,dbg_j,nlev)
!        write (gol,*) 'yyy after  xn_adv = ', xn_adv(ispec,dbg_i,dbg_j,nlev); call goPr
!      end do
!      write (gol,*) 'yyy after  tpm    = ', tpm; call goPr
!    end if

    ! ok
    status = 0
    
  end subroutine ObsCompInfo_DistributeIncrement_3d


  ! *

  
  !
  ! Change single tracer in xn_adv array given observed field (total pm?) and increment.
  ! This routine was made to avoid the need for a copy of the full tracer array in case
  ! assimilation result is requested as post-processing only (not fed back).
  !
  
  subroutine ObsCompInfo_DistributeIncrement1_3d( self, ispec_adv, xn_adv1, xn_adv1_units, &
                                                   xn_obs, dx_obs, obs_units, status, &
                                                   maxratio, verbose )
  
    use MetFields_ml , only : roa
    use Chemfields_ml, only : PM25_water

    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(in)      ::  self
    integer, intent(in)                   ::  ispec_adv
    real, intent(inout)                   ::  xn_adv1(:,:,:)   ! (lnx,lny,nz) valid for ispec_adv
    character(len=*), intent(in)          ::  xn_adv1_units
    real, intent(in)                      ::  xn_obs(:,:,:)    ! (lnx,lny,nz)
    real, intent(in)                      ::  dx_obs(:,:,:)    ! (lnx,lny,nz)
    character(len=*), intent(in)          ::  obs_units
    integer, intent(out)                  ::  status
    
    real, intent(in), optional            ::  maxratio
    logical, intent(in), optional         ::  verbose
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_DistributeIncrement1_3d'
    
    ! --- local -----------------------------
    
    real, allocatable     ::  ratio(:,:,:)
    integer               ::  i
    integer               ::  ispec
    logical               ::  verb
    
    ! --- begin -----------------------------
    
    ! verbose?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! single spec? then just add the increment:
    if ( self%nspec == 1 ) then
    
      ! single spec:
      i = 1
      ! current:
      ispec = self%ispec(i)
      ! target?
      if ( ispec == ispec_adv ) then
        !! testing ..
        !if ( verb ) then
        !  write (gol,*) 'xxx reset tracer ', ispec, ' from ', minval(xn_adv1), ' - ', maxval(xn_adv1); call goPr
        !end if
        ! check units ..
        if ( trim(xn_adv1_units) /= trim(obs_units) ) then
          write (gol,'("xn_adv units `",a,"` while obs units `",a,"`")') trim(xn_adv1_units), trim(obs_units); call goErr
          TRACEBACK; status=1; return
        end if
        ! add increment:
        xn_adv1 = max( 0.0, xn_adv1 + dx_obs )
        !! testing ..
        !if ( verb ) then
        !  write (gol,*) 'xxx                             to ', minval(xn_adv1), ' - ', maxval(xn_adv1); call goPr
        !end if
      end if
      
    else
    
      ! composed of sum, distribute according to ratio;
      ! note that if sum is zero, no distribution can be made ...
      ! also no need to check units, ratio will have units 1 anyway
    
      ! loop over original species:
      do i = 1, self%nspec
        ! current:
        ispec = self%ispec(i)
        ! target?
        if ( ispec == ispec_adv ) then

          ! storage for relative change:
          allocate( ratio(size(xn_obs,1),size(xn_obs,2),size(xn_obs,3)), stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! fill ratio, or set to unity if undefined (zero tracer sum):
          where( xn_obs > 0.0 )
            ratio = max( 0.0, ( xn_obs + dx_obs ) / xn_obs )
          elsewhere
            ratio = 1.0
          end where

          ! maximum?
          if ( present(maxratio) ) ratio = min( ratio, maxratio )
          
          !! testing ..
          !if ( verb ) then
          !  write (gol,*) 'xxx scale tracer ', ispec, ' with ratios ', minval(ratio), ' - ', maxval(ratio); call goPr
          !end if

!          ! testing ...
!          if ( dbg_cell ) then
!            write (gol,*) 'yyy update: ', xn_adv1(dbg_i,dbg_j,nlev), ratio(dbg_i,dbg_j,nlev), xn_adv1(dbg_i,dbg_j,nlev) * ratio(dbg_i,dbg_j,nlev); call goPr
!          end if

          ! scale
          xn_adv1 = xn_adv1 * ratio
          
          ! also apply to aerosol water?
          ! select first species to perform this:
          if ( self%with_pmwater .and. (i == 1) ) then
            !! testing ..
            !if ( verb ) then
            !  write (gol,*) 'xxx scale pm water ...'; call goPr
            !end if
            ! apply same ratio:
            PM25_water = PM25_water * ratio
          end if
          
          ! clear:
          deallocate( ratio, stat=status )
          IF_NOT_OK_RETURN(status=1)

        end if ! target sepc
      end do ! specs

    end if  ! single spec or sum

    ! ok
    status = 0
    
  end subroutine ObsCompInfo_DistributeIncrement1_3d


  ! *


  !
  ! Input:
  !   xn_adv    :  model traces (all)
  !   xn_obs    :  2D simulated observation field
  !   dx_obs    :  3D analysis increment, layer nz has the increment towards xn_obs,
  !                other layers filled using vertical correlations in B
  ! Change tracers in 3D fields "xn_adv" following 2D ratio:
  !    (xn_obs+dx_obs(:,:,nz))/xn_obs
  ! This is the ratio for the observation field (surface?),
  ! with decay to zero for higher levels (1,..,nz-1)
  !
  
  subroutine ObsCompInfo_ChangeProfile( self, xn_adv, xn_adv_units, &
                                       xn_obs, dx_obs, obs_units, status, &
                                       maxratio, verbose, bl )

    use MetFields_ml , only : z_bnd    ! (nx,ny,1:nz+1) gph bounds relative to surface
    use MetFields_ml , only : pzpbl    ! (nx,ny) boundary layer height
    use Chemfields_ml, only : PM25_water         ! (i,j,k) model levels
    use Chemfields_ml, only : PM25_water_rh50    ! (i,j) surface

    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(in)      ::  self
    real, intent(inout)                   ::  xn_adv(:,:,:,:)  ! (nspec_adv,lnx,lny,nz)
    character(len=*), intent(in)          ::  xn_adv_units(:)  ! (nspec_adv)
    real, intent(in)                      ::  xn_obs(:,:)      ! (lnx,lny)
    real, intent(in)                      ::  dx_obs(:,:,:)    ! (lnx,lny,nz)
    character(len=*), intent(in)          ::  obs_units
    integer, intent(out)                  ::  status
    
    real, intent(in), optional            ::  maxratio
    logical, intent(in), optional         ::  verbose
    logical, intent(in), optional         ::  bl
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_ChangeProfile'
    
    ! --- local -----------------------------
    
    logical               ::  verb
    integer               ::  nspec, lnx, lny, nz
    integer               ::  ispec
    integer               ::  i, j, k
    real, allocatable     ::  ratio(:,:,:)  ! (lnx,lny,nz)
    logical               ::  only_bl

    ! --- begin -----------------------------
    
    ! verbose?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! boundary layer only?
    only_bl = .false.
    if ( present(bl) ) only_bl = bl
    
    ! shape:
    nspec = size(xn_adv,1)
    lnx   = size(xn_adv,2)
    lny   = size(xn_adv,3)
    nz    = size(xn_adv,4)
    
    ! storage for relative change:
    allocate( ratio(lnx,lny,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! only bl, or decay with height?
    if ( only_bl ) then
    
      ! init as no change:
      ratio = 1.0
    
      ! loop over cells:
      do i = 1, lnx
        do j = 1, lny
          ! loop over levels
          do k = 1, nz
            ! (partly) in boundary layer?
            if ( z_bnd(i,j,k+1) < pzpbl(i,j) ) then
              ! relative change at lowest layer:
              if ( xn_obs(i,j) > 0.0 ) then
                !              (   "analyzed" wrt obs        ) /  fg
                ratio(i,j,k) = ( xn_obs(i,j) + dx_obs(i,j,nz ) ) / xn_obs(i,j)
              end if
            end if ! in bl
          end do ! k
        end do ! j
      end do ! i
    
    else
    
      ! loop over levels (layer nz is equivalent to obs)
      do k = 1, nz
        ! fill ratio, or set to unity if undefined (zero tracer sum):
        where( xn_obs > 0.0 )
          !              (   "analyzed" wrt obs   ) /  fg
          ratio(:,:,k) = ( xn_obs + dx_obs(:,:,k ) ) / xn_obs
          !! TESTING: same ratio at all levels:
          !ratio(:,:,k) = ( xn_obs + dx_obs(:,:,nz) ) / xn_obs
        elsewhere
          ratio(:,:,k) = 1.0
        end where
      end do  ! levels 

      !! testing ...
      !if ( dbg_cell ) then
      !  write (gol,*) 'xxx ratio = ', ratio(dbg_i,dbg_j,:); call goPr
      !end if

    end if

    ! use minimum ratio to avoid that all aersols are removed,
    ! this would give errors from aerosol equilibrium routines:
    ratio = max( 0.01, ratio )

    ! maximum?
    if ( present(maxratio) ) ratio = min( ratio, maxratio )
    
    ! loop over contributing species:
    do i = 1, self%nspec
      ! global index:
      ispec = self%ispec(i)
      ! check units ..
      !if ( trim(xn_adv_units(ispec)) /= trim(obs_units) ) then
      !  write (gol,'("xn_adv units `",a,"` while obs units `",a,"`")') trim(xn_adv_units(ispec)), trim(obs_units); call goErr
      ! TRACEBACK; status=1; return
      !end if
      ! scale 3D:
      xn_adv(ispec,:,:,:) = xn_adv(ispec,:,:,:) * ratio
      !! testing ...
      !if ( dbg_cell ) then
      !  write (gol,*) 'xxx xn_adv(',ispec,')=',xn_adv(ispec,dbg_i,dbg_j,:); call goPr
      !end if
    end do ! specs

    ! also apply to aerosol water?
    !   Aero_water from AerosolCalls
    !   PM25_water at model levels (ambient)
    !   PM25_water_rh50 at surface (50% rh)
    if ( self%with_pmwater ) then
      ! apply same ratios:
      PM25_water      = PM25_water      * ratio          ! model levels
      PM25_water_rh50 = PM25_water_rh50 * ratio(:,:,nz)  ! surface
      !! testing ...
      !if ( dbg_cell ) then
      !  write (gol,*) 'xxx PM25_water      = ', PM25_water(dbg_i,dbg_j,:); call goPr
      !  write (gol,*) 'xxx PM25_water_rh50 = ', PM25_water_rh50(dbg_i,dbg_j); call goPr
      !end if
    end if
    
    ! clear:
    deallocate( ratio, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine ObsCompInfo_ChangeProfile


  ! *


!  !
!  ! Input:
!  !   xn_adv    :  model traces (all)
!  !   xn_obs    :  2D simulated observation field
!  !   dx_obs    :  3D analysis increment, layer nz has the increment towards xn_obs,
!  !                other layers filled using vertical correlations in B
!  !
!  ! Change ratio of fine/coarse aerosol in 3D fields "xn_adv"
!  ! to reach at the surface (layer nz) the analysed total fine part:
!  !    xn_obs+dx_obs(:,:,nz)
!  ! The increment "dx_obs" decays to zero for higher levels (1,..,nz-1)
!  !
!  ! If bl=.true. this decay is not used but fine/coarse ratio in boundary layer
!  ! are scaled with the same factor as the surface field.
!  !
!  
!  subroutine ObsCompInfo_ChangeProfile_FineCoarseRatio( self, xn_adv, xn_adv_units, &
!                                                     xn_obs, dx_obs, obs_units, status, &
!                                                     verbose, bl )
!
!    use ChemGroups_ml    , only : PMFINE_GROUP, PM10_GROUP   ! requires offset NSPEC_SHL !
!    use ChemSpecs_shl_ml , only : NSPEC_SHL  ! number of short-lived tracers
!    use ChemChemicals_ml , only : species_adv
!    use Units_ml         , only : Units_Scale
!    use MetFields_ml     , only : roa
!    use MetFields_ml     , only : z_bnd    ! (nx,ny,1:nz+1) gph bounds relative to surface
!    use MetFields_ml     , only : pzpbl    ! (nx,ny) boundary layer height
!    use Chemfields_ml    , only : cfac               ! (i,j) 50m -> 3m factor
!    use Chemfields_ml    , only : PM25_water         ! (i,j,k) model levels
!    use Chemfields_ml    , only : PM25_water_rh50    ! (i,j) surface
!    use ModelConstants_ml, only : lev3d, num_lev3d
!
!    ! --- in/out ----------------------------
!    
!    class(T_ObsCompInfo), intent(in)      ::  self
!    real, intent(inout)                   ::  xn_adv(:,:,:,:)  ! (nspec_adv,lnx,lny,nz)
!    character(len=*), intent(in)          ::  xn_adv_units(:)  ! (nspec_adv)
!    real, intent(in)                      ::  xn_obs(:,:)      ! (lnx,lny)
!    real, intent(in)                      ::  dx_obs(:,:,:)    ! (lnx,lny,nz)
!    character(len=*), intent(in)          ::  obs_units
!    integer, intent(out)                  ::  status
!    
!    logical, intent(in), optional         ::  verbose
!    logical, intent(in), optional         ::  bl
!    
!    ! --- const ----------------------------
!    
!    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_ChangeProfile_FineCoarseRatio'
!    
!    ! fraction of NO3c assigned to PM2.5:
!    real, parameter    ::  alfa = 0.27
!    
!    ! pm collections:
!    integer, parameter        ::  npm = 3
!    integer, parameter        ::  ipm_fine   = 1
!    integer, parameter        ::  ipm_coarse = 2
!    integer, parameter        ::  ipm_no3c   = 3
!    
!    ! --- local -----------------------------
!    
!    logical               ::  verb
!    logical               ::  only_bl
!    integer               ::  nspec, lnx, lny, nz
!    real, allocatable     ::  tpm_sfc(:,:,  :)  ! (lnx,lny   ,npm)
!    real, allocatable     ::  tpm_ml (:,:,:,:)  ! (lnx,lny,nz,npm)
!    real, allocatable     ::  xn_25(:,:)        ! (lnx,lny)
!    real, allocatable     ::  xn_10(:,:)        ! (lnx,lny)
!    real, allocatable     ::  factor_sfc(:,:  ,:)  ! (lnx,lny   ,npm)
!    real, allocatable     ::  factor_ml (:,:,:,:)  ! (lnx,lny,nz,npm)
!    logical, allocatable  ::  ispm(:,:)            ! (nspec_adv,npm)
!    integer               ::  ipm
!    integer, pointer      ::  out_group(:)
!    integer               ::  out_offset
!    integer               ::  ispec
!    integer               ::  ispec_no3c
!    real                  ::  fscale
!    logical               ::  needroa
!    integer               ::  k
!
!      integer   :: i,j
!
!    ! --- begin -----------------------------
!    
!    ! verbose?
!    verb = .false.
!    if ( present(verbose) ) verb = verbose
!    
!    ! boundary layer only?
!    only_bl = .false.
!    if ( present(bl) ) only_bl = bl
!    
!    ! shape:
!    nspec = size(xn_adv,1)
!    lnx   = size(xn_adv,2)
!    lny   = size(xn_adv,3)
!    nz    = size(xn_adv,4)
!    
!    ! target obs:
!    allocate( xn_25(lnx,lny), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( xn_10(lnx,lny), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    ! storage for current total fine/coarse aerosol:
!    allocate( tpm_ml (lnx,lny,nz,npm), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( tpm_sfc(lnx,lny   ,npm), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    ! storage for factors:
!    allocate( factor_sfc(lnx,lny   ,npm), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    allocate( factor_ml (lnx,lny,nz,npm), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    
!    ! Coarse nitrate should be treated special ...
!    ! Fraction alfa is assigned to "fine" and (1-alfa) to "coarse",
!    ! thus fine/coarse ratio "within" no3c is fixed.
!    ! This fixed ratio should not be changed by the assimilation;
!    ! actually it is imposible to maintain the ratio while changing
!    ! the total fine/coarse fractions.
!    ! Therefore do not apply scale factors to coarse nitrate,
!    ! only assimilation of pm10 will change its amount.
!    ! Find index of coarse nitrate:
!    ispec_no3c = find_index( 'NO3_C', species_adv(:)%name )
!
!    ! weights of tracers in fine and corase:
!    allocate( ispm(nspec,npm), stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    ! init with no contribution:
!    ispm = .false.    
!    ! loop over pm collections:
!    do ipm = 1, npm
!      ! use group defined in 'ChemGroups_ml':
!      if ( ipm == ipm_fine ) then
!        out_group => PMFINE_GROUP
!      else
!        out_group => PM10_GROUP
!      end if
!      ! offset:
!      out_offset = NSPEC_SHL
!      ! loop over group members:
!      do k = 1, size(out_group)
!        ! extract index of specie:
!        ispec = out_group(k) - out_offset
!        ! store relative weight in fine and coarse parts:
!        if ( ipm == ipm_fine ) then
!          ! contributes to fine part:
!          ispm(ispec,ipm) = .true.
!        !  
!        else if ( ipm == ipm_coarse ) then
!          ! only if exclusevely coarse and not coarse nitrate:
!          ispm(ispec,ipm) = .not. ( ispm(ispec,ipm_fine) .or. (ispec == ispec_no3c) )
!        !
!        else if ( ipm == ipm_no3c ) then
!          ! only if coarse nitrate:
!          ispm(ispec,ipm) = ispec == ispec_no3c
!        end if
!        !
!      end do ! group elements
!    end do  ! total fine and coarse
!    ! special:
!    ispm(ispec_no3c,ipm_no3c) = .true.
!    
!    !! info ...
!    !write (gol,*) 'rrr weights: '; call goPr
!    !do ispec = 1, nspec
!    !  if ( any( ispm(ispec,:) ) ) then
!    !    write (gol,'("rrr ",3l2," ",a)') ispm(ispec,:), trim(species_adv(ispec)%name); call goPr
!    !  end if
!    !end do
!
!    ! init sums:
!    tpm_ml  = 0.0
!    tpm_sfc = 0.0
!    ! loop over fine (=1) and coarse-without-no3c (=2):
!    do ipm = 1, 3
!      !! testing ...
!      !if ( dbg_cell ) then
!      !  write (gol,*) 'yyy ipm = ', ipm; call goPr
!      !end if
!      ! loop over species:
!      do ispec = 1, nspec
!        ! skip if no contribution:
!        if ( .not. ispm(ispec,ipm) ) cycle
!        ! info on conversion to target units:
!#ifdef with_ajs
!        call Units_Scale( obs_units, ispec, fscale, &
!                                 needroa=needroa, status=status )
!        IF_NOT_OK_RETURN(status=1)
!#else
!        call Units_Scale( obs_units, ispec, fscale, &
!                                 needroa=needroa )
!#endif
!        !! testing ...
!        !write (gol,*) 'qqq ', ipm, ' ', obs_units, ispec, ' ', trim(species_adv(ispec)%name), fscale, needroa; call goPr
!        ! convert using air density?
!        if ( needroa ) then
!          ! add contribution, multiply with air density:
!          tpm_ml(:,:,:,ipm) = tpm_ml(:,:,:,ipm) + &
!                               xn_adv(ispec,:,:,lev3d(:num_lev3d)) &
!                               * fscale &
!                               * roa(:,:,lev3d(:num_lev3d),1)    ! density
!          ! surface fraction:
!          tpm_sfc(:,:,ipm)  = tpm_sfc(:,:,ipm) + &
!                               xn_adv(ispec,:,:,nlev) &
!                               * fscale &
!                               * roa(:,:,nlev,1)    &  ! density
!                               * cfac(ispec,:,:)       ! 50 m -> 3 m
!          !! testing ...
!          !if ( dbg_cell ) then
!          !  write (gol,*) 'yyy   ispec ', ispec, xn_adv(ispec,dbg_i,dbg_j,nlev) * w(ispec,ipm) &
!          !                                       * fscale * roa(dbg_i,dbg_j,nlev,1) * cfac(ispec,dbg_i,dbg_j); call goPr
!          !end if
!                               
!        else
!          ! add contribution:
!          tpm_ml(:,:,:,ipm) = tpm_ml(:,:,:,ipm) + &
!                               xn_adv(ispec,:,:,lev3d(:num_lev3d)) &
!                               * fscale
!          ! surface fraction:
!          tpm_sfc(:,:,ipm)  = tpm_sfc(:,:,ipm) + &
!                               xn_adv(ispec,:,:,nlev) &
!                               * fscale &
!                               * cfac(ispec,:,:)       ! 50 m -> 3 m
!          !! testing ...
!          !if ( dbg_cell ) then
!          !  write (gol,*) 'yyy   ispec ', ispec, xn_adv(ispec,dbg_i,dbg_j,nlev) * w(ispec,ipm) &
!          !                                       * fscale * cfac(ispec,dbg_i,dbg_j); call goPr
!          !end if
!        end if
!      end do  ! spec
!    end do  ! pm collections
!    ! add 3D aersosol water to fine fraction:
!    tpm_ml (:,:,:,ipm_fine) = tpm_ml (:,:,:,ipm_fine) + PM25_water
!    tpm_sfc(:,:  ,ipm_fine) = tpm_sfc(:,:  ,ipm_fine) + PM25_water_rh50
!    !! testing ...
!    !if ( dbg_cell ) then
!    !  write (gol,*) 'yyy   rh50           = ', PM25_water_rh50(dbg_i,dbg_j); call goPr
!    !end if
!    
!    !! testing ...
!    !if ( dbg_cell ) then
!    !  write (gol,*) 'yyy before tpm fine,course   = ', tpm_sfc(dbg_i,dbg_j,:); call goPr
!    !  write (gol,*) 'yyy before tpm fine+coarse   = ', sum(tpm_sfc(dbg_i,dbg_j,:)); call goPr
!    !  write (gol,*) 'yyy before tpm f + a*no3c    = ', tpm_sfc(dbg_i,dbg_j,1)+alfa*tpm_sfc(:,:,3)(dbg_i,dbg_j); call goPr
!    !  write (gol,*) 'yyy before tpm f + c + no3c  = ', sum(tpm_sfc(dbg_i,dbg_j,:)) + tpm_sfc(:,:,3)(dbg_i,dbg_j); call goPr
!    !end if
!
!    !! testing ..
!    !write (gol,*) 'ppp0 tpm = ', tpm_sfc(1,1,:), sum(tpm_sfc(1,1,:)); call goPr
!    
!    ! analyzed total pm25 at the surcface:  xn_obs + dx_obs(:,:,nz)
!    ! in theory this could not exceed total pm10,
!    ! but total pm10 was analyzed independently,
!    ! and sometimes observations are found with pm25 > pm10 ;
!    ! limit the analyzed field to the current total 
!    ! (which should remain unchanged in this update):
!    xn_10 = sum( tpm_sfc, dim=3 )
!    xn_25 = min( max( 0.0, xn_obs + dx_obs(:,:,nz) ), xn_10 )
!    !! without the water part:
!    !xn_ana = min( xn_obs + dx_obs - PM25_water, tpm(:,:,:,1)+tpm(:,:,:,ipm_coarse) )
!    
!    !! testing ..
!    !write (gol,*) 'ppp0 xn_ana = ',xn_obs(1,1), dx_obs(1,1), xn_ana(1,1); call goPr
!
!    !! testing ...
!    !if ( dbg_cell ) then
!    !  write (gol,*) 'yyy xn_obs, dx_obs  = ', xn_obs(dbg_i,dbg_j), dx_obs(dbg_i,dbg_j,nz); call goPr
!    !  write (gol,*) 'yyy xn_obs + dx_obs = ', xn_obs(dbg_i,dbg_j) + dx_obs(dbg_i,dbg_j,nz); call goPr
!    !  write (gol,*) 'yyy xn_ana          = ', xn_ana(dbg_i,dbg_j); call goPr
!    !end if
!    
!    !
!    ! Determine scale factors for fine and coarse part such that:
!    !
!    !   total_fine * fine_factor +  alfa*no3c                 * coarse_factor = xn_25
!    !   total_fine * fine_factor + (     no3c + total_coarse) * coarse_factor = xn_10
!    !
!    ! First solve coarse factor:
!    !
!    !   xn_25 - alfa*no3c * coarse_factor + ( no3c + total_coarse) * coarse_factor = xn_10
!    !   ( (1-alfa)*no3c + total_coarse ) * coarse_factor = xn_10 - xn_25
!    !     coarse_factor = (xn_10 - xn_25)/( (1-alfa)*no3c + total_coarse )
!    !
!    ! Thus, scale the coarse fraction to match with analyzed coarse fraction.
!    !
!    where ( xn_10 > xn_25 )
!      ! scale factor for coarse fraction:
!      factor_sfc(:,:,ipm_coarse) = (xn_10 - xn_25) / ( (1.0-alfa)*tpm_sfc(:,:,ipm_no3c) + tpm_sfc(:,:,ipm_coarse) )
!    elsewhere
!      ! no coarse fraction in analysis present, dummy factor:
!      factor_sfc(:,:,ipm_coarse) = 1.0
!    end where
!
!        ! testing ...
!        if ( dbg_cell ) then
!          write (gol,*) 'xxx xn_10             = ', xn_10(dbg_i,dbg_j); call goPr
!          write (gol,*) 'xxx xn_25             = ', xn_25(dbg_i,dbg_j); call goPr
!          write (gol,*) 'xxx tpm_sfc coarse    = ', tpm_sfc(dbg_i,dbg_j,ipm_coarse); call goPr
!          write (gol,*) 'xxx tpm_sfc no3c      = ', tpm_sfc(dbg_i,dbg_j,ipm_no3c); call goPr
!          write (gol,*) 'xxx factor_sfc coarse = ', factor_sfc(dbg_i,dbg_j,ipm_coarse); call goPr
!        end if
!    
!    ! scale factor to make fine fraction (incl part of no3c) equal to observations:
!    !   total_fine * fine_factor + alfa*no3c * coarse_factor =  xn_25
!    ! thus:
!    !   fine_factor  =  (xn_25 - alfa*no3c * coarse_factor) / total_fine
!    ! any fine pm?
!    where ( tpm_sfc(:,:,ipm_fine) > 0.0 )
!      ! fine pm present, obtain scale factor;
!      ! might become negative which means that all fine pm shoulde be removed
!      ! except the no3c fration that is also considered as part of the fine mode;
!      ! keep minimum fraction to avoid that all fine aerosols are removed:
!      factor_sfc(:,:,ipm_fine) =  max( 0.01, ( xn_25 - alfa*tpm_sfc(:,:,ipm_no3c) * factor_sfc(:,:,ipm_coarse) ) / tpm_sfc(:,:,ipm_fine) )
!    elsewhere
!      ! no fine pm present, dummy factor:
!      factor_sfc(:,:,ipm_fine) = 1.0
!    end where
!    !! testing ...
!    !if ( dbg_cell ) then
!    !  write (gol,*) 'yyy fine_factor 1 = ', factor_sfc(:,:,ipm_fine)(dbg_i,dbg_j); call goPr
!    !end if
!    ! use minimum factor to avoid that all aersols are removed,
!    ! this would give errors from aerosol equilibrium routines:
!    !factor_sfc(:,:,ipm_fine) = max( 0.01, factor_sfc(:,:,ipm_fine) )
!    !! testing ...
!    !if ( dbg_cell ) then
!    !  write (gol,*) 'yyy fine_factor 2 = ', factor_sfc(:,:,ipm_fine)(dbg_i,dbg_j); call goPr
!    !end if
!
!        ! testing ...
!        if ( dbg_cell ) then
!          write (gol,*) 'xxx tpm_sfc fine   = ', tpm_sfc(dbg_i,dbg_j,ipm_fine); call goPr
!          write (gol,*) 'xxx tpm_sfc coarse = ', tpm_sfc(dbg_i,dbg_j,ipm_coarse); call goPr
!          write (gol,*) 'xxx tpm_sfc no3c   = ', tpm_sfc(dbg_i,dbg_j,ipm_no3c); call goPr
!          write (gol,*) 'xxx factor_sfc fine = ', factor_sfc(dbg_i,dbg_j,ipm_fine); call goPr
!        end if
!    
!    ! if for the fine factor the minimum value is used,
!    ! then both the modelled pm25 and pm10 will exceed the analyzed pm25 and pm10 
!    ! since pm25 has too much of the coarse nitrate fraction;
!    ! remove the surplus of coarse nitrate:
!    where ( ( (1-alfa)*tpm_sfc(:,:,ipm_no3c) + tpm_sfc(:,:,ipm_coarse) )*factor_sfc(:,:,ipm_coarse) > xn_10 )
!      ! special scale factor for coarse nitrate:
!      !   total_fine * fine_factor + no3c * no3c_factor + total_coarse * coarse_factor = xn_10
!      !   no3c_factor = ( xn_10 - total_fine * fine_factor - total_coarse * coarse_factor ) / no3c
!      ! keep minimum ammount:
!      factor_sfc(:,:,ipm_no3c) = ( xn_10 - tpm_sfc(:,:,ipm_fine)*factor_sfc(:,:,ipm_fine) - tpm_sfc(:,:,ipm_coarse)*factor_sfc(:,:,ipm_coarse) ) / ( (1-alfa)*tpm_sfc(:,:,ipm_no3c) )
!    elsewhere
!      ! keep already computed coarse factor:
!      factor_sfc(:,:,ipm_no3c) = factor_sfc(:,:,ipm_coarse)
!    endwhere
!
!        ! testing ...
!        if ( dbg_cell ) then
!          write (gol,*) 'xxx factor_sfc no3c = ', factor_sfc(dbg_i,dbg_j,ipm_no3c); call goPr
!        end if
!
!    ! expand to higher layers; same for all collections:
!    do ipm = 1, npm
!
!      ! copy surface:
!      factor_ml(:,:,nz,ipm) = factor_sfc(:,:,ipm)
!
!      ! only boundary layer, or decay with dx ?
!      if ( only_bl ) then
!    
!        ! loop over cells:
!        do i = 1, lnx
!          do j = 1, lny
!            ! loop over levels except surface:
!            do k = 1, nz-1
!              ! (partly) in boundary layer?
!              if ( z_bnd(i,j,k+1) < pzpbl(i,j) ) then
!                ! use surface ratio:
!                factor_ml(i,j,k,ipm) = factor_ml(i,j,nz,ipm)
!              else
!                ! no change:
!                factor_ml(i,j,k,ipm) = 1.0
!              end if ! in bl
!            end do ! k
!          end do ! j
!        end do ! i
!
!      else
!
!        ! other layers:
!        do k = 1, nz-1
!          ! use profile following ratio:
!          !   dx(:,:,k)/dx(:,:,nz)
!          ! which is 1 at the bottom layer (nz) and decays to zero with height;
!          ! note that also values >1 might be present in the first layers above the surface,
!          ! and tiny negative values at the top, therefore use abs values;
!          ! also use minimum factor:
!          where ( abs(dx_obs(:,:,nz)) > 0.0 )
!            factor_ml(:,:,k,ipm) = max( 0.01, 1.0 + ( factor_ml(:,:,nz,ipm) - 1.0 ) * abs(dx_obs(:,:,k)/dx_obs(:,:,nz)) )
!          elsewhere
!            factor_ml(:,:,k,ipm) = 1.0
!          end where
!        end do ! levels
!
!        ! testing ...
!        if ( dbg_cell ) then
!          write (gol,*) 'xxx factor_ml(',ipm,')'; call goPr
!          write (gol,*) 'xxx factor_ml at nz ', factor_ml(dbg_i,dbg_j,nz,ipm); call goPr
!          do k = 1, nz-1
!            write (gol,*) 'xxx  layer ',k, dx_obs(dbg_i,dbg_j,k), dx_obs(dbg_i,dbg_j,nz), &
!                              abs(dx_obs(dbg_i,dbg_j,k)/dx_obs(dbg_i,dbg_j,nz)), &
!                               factor_ml(dbg_i,dbg_j,k,ipm); call goPr
!          end do
!        end if
!
!      end if ! bl only or dx
!
!    end do ! ipm
!    
!    !! info ..
!    !write (gol,*) 'ppp2 fine: ', xn_ana(1,ipm_fine), tpm_sfc(1,1,ipm_fine), xn_ana(1,ipm_fine) / tpm_sfc(1,1,ipm_fine); call goPr
!    !write (gol,*) 'ppp2 fine   factor = ', fine_factor  (1,ipm_fine); call goPr
!    !write (gol,*) 'ppp2 coarse factor = ', coarse_factor(1,ipm_fine); call goPr
!
!    ! loop over fine/coarse:
!    do ipm = 1, npm
!      ! loop over all species:
!      do ispec = 1, nspec
!        ! skip if no contribution to pm:
!        if ( .not. ispm(ispec,ipm) ) cycle
!        !! testing ...
!        !if ( dbg_cell ) then
!        !  write (gol,*) 'yyy ispec ', ispec; call goPr
!        !  !write (gol,*) 'yyy   wf ml  ', w(ispec,ipm_fine) * factor_ml(:,:,ipm_fine)(dbg_i,dbg_j,nz), w(ispec,ipm_coarse) * factor_ml(:,:,ipm_coarse)(dbg_i,dbg_j,nz); call goPr
!        !  write (gol,*) 'yyy   wf ml  ', w(ispec,ipm_fine) * factor_ml(:,:,ipm_fine)(dbg_i,dbg_j,nz) + w(ispec,ipm_coarse) * factor_ml(:,:,ipm_coarse)(dbg_i,dbg_j,nz); call goPr
!        !  write (gol,*) 'yyy   xn_adv ', xn_adv(ispec,dbg_i,dbg_j,nz); call goPr
!        !end if
!        ! apply factor fields:
!        xn_adv(ispec,:,:,:) = xn_adv(ispec,:,:,:) * factor_ml(:,:,:,ipm)
!        ! testing ...
!        if ( dbg_cell ) then
!          !write (gol,*) 'yyy   xn_adv ', xn_adv(ispec,dbg_i,dbg_j,nz); call goPr
!          write (gol,*) 'xxx xn_adv(',ispec,')=',xn_adv(ispec,dbg_i,dbg_j,:); call goPr
!        end if
!      end do ! ispec
!    end do ! ipm
!    ! also apply (fine) factor to aersol water;
!    ! better would be to recompute water content after changing the aerosol load ...
!    PM25_water      = PM25_water      * factor_ml (:,:,:,ipm_fine)   ! model levels
!    PM25_water_rh50 = PM25_water_rh50 * factor_sfc(:,:  ,ipm_fine)   ! sfc
!    
!      ! testing ...
!      if ( dbg_cell ) then
!        write (gol,*) 'xxx PM25_water      = ', PM25_water(dbg_i,dbg_j,:); call goPr
!        write (gol,*) 'xxx PM25_water_rh50 = ', PM25_water_rh50(dbg_i,dbg_j); call goPr
!      end if
!
!!    ! .... testing ........................
!!
!!    ! init sums:
!!    tpm_ml  = 0.0
!!    tpm_sfc = 0.0
!!    ! loop over fine (=1) and coarse (=2):
!!    do ipm = 1, 2
!!      !! testing ...
!!      !if ( dbg_cell ) then
!!      !  write (gol,*) 'yyy ipm = ', ipm; call goPr
!!      !end if
!!      ! loop over species:
!!      do ispec = 1, nspec
!!        ! skip if no contribution:
!!        if ( w(ispec,ipm) <= 0.0 ) cycle
!!        ! info on conversion to target units:
!!#ifdef with_ajs
!!        call Units_Scale( obs_units, ispec, fscale, &
!!                                 needroa=needroa, status=status )
!!        IF_NOT_OK_RETURN(status=1)
!!#else
!!        call Units_Scale( obs_units, ispec, fscale, &
!!                                 needroa=needroa )
!!#endif
!!        !! testing ...
!!        !write (gol,*) 'qqq ', ipm, ' ', obs_units, ispec, ' ', trim(species_adv(ispec)%name), fscale, needroa; call goPr
!!        ! convert using air density?
!!        if ( needroa ) then
!!          ! add contribution, multiply with air density:
!!          tpm_ml(:,:,:,ipm) = tpm_ml(:,:,:,ipm) + &
!!                               xn_adv(ispec,:,:,lev3d(:num_lev3d)) * w(ispec,ipm) &
!!                               * fscale &
!!                               * roa(:,:,lev3d(:num_lev3d),ipm_fine)    ! density
!!          ! surface fraction:
!!          tpm_sfc(:,:,ipm)  = tpm_sfc(:,:,ipm) + &
!!                               xn_adv(ispec,:,:,nlev) * w(ispec,ipm) &
!!                               * fscale &
!!                               * roa(:,:,nlev,ipm_fine)    &  ! density
!!                               * cfac(ispec,:,:)       ! 50 m -> 3 m
!!          !! testing ...
!!          !if ( dbg_cell ) then
!!          !  write (gol,*) 'yyy   r ispec ', ispec, xn_adv(ispec,dbg_i,dbg_j,nlev) * w(ispec,ipm) &
!!          !                                       * fscale * roa(dbg_i,dbg_j,nlev,ipm_fine) * cfac(ispec,dbg_i,dbg_j); call goPr
!!          !end if
!!                               
!!        else
!!          ! add contribution:
!!          tpm_ml(:,:,:,ipm) = tpm_ml(:,:,:,ipm) + &
!!                               xn_adv(ispec,:,:,lev3d(:num_lev3d)) * w(ispec,ipm) &
!!                               * fscale
!!          ! surface fraction:
!!          tpm_sfc(:,:,ipm)  = tpm_sfc(:,:,ipm) + &
!!                               xn_adv(ispec,:,:,nlev) * w(ispec,ipm) &
!!                               * fscale &
!!                               * cfac(ispec,:,:)       ! 50 m -> 3 m
!!          !! testing ...
!!          !if ( dbg_cell ) then
!!          !  write (gol,*) 'yyy   p ispec ', ispec, xn_adv(ispec,dbg_i,dbg_j,nlev) * w(ispec,ipm) &
!!          !                                       * fscale * cfac(ispec,dbg_i,dbg_j); call goPr
!!          !end if
!!        end if
!!      end do  ! spec
!!      ! add 3D aersosol water to fine fraction:
!!      if ( ipm == 1 ) then
!!        tpm_ml (:,:,:,ipm) = tpm_ml (:,:,:,ipm) + PM25_water
!!        tpm_sfc(:,:  ,ipm) = tpm_sfc(:,:  ,ipm) + PM25_water_rh50
!!        !! testing ...
!!        !if ( dbg_cell ) then
!!        !  write (gol,*) 'yyy   rh50           = ', PM25_water_rh50(dbg_i,dbg_j); call goPr
!!        !end if
!!      end if
!!    end do  ! total fine and coarse
!!    
!!    ! coarse nitrate:
!!    ispec = ispec_no3c
!!    ! info on conversion to target units:
!!#ifdef with_ajs
!!    call Units_Scale( obs_units, ispec, fscale, &
!!                             needroa=needroa, status=status )
!!    IF_NOT_OK_RETURN(status=1)
!!#else
!!    call Units_Scale( obs_units, ispec, fscale, &
!!                             needroa=needroa )
!!#endif
!!    !! testing ...
!!    !write (gol,*) 'qqq ', ipm, ' ', obs_units, ispec, ' ', trim(species_adv(ispec)%name), fscale, needroa; call goPr
!!    ! convert using air density?
!!    if ( needroa ) then
!!      ! surface fraction:
!!      tpm_sfc(:,:,ipm_no3c) = xn_adv(ispec,:,:,nlev) * fscale &
!!                 * roa(:,:,nlev,ipm_fine)    &  ! density
!!                 * cfac(ispec,:,:)       ! 50 m -> 3 m
!!      !! testing ...
!!      !if ( dbg_cell ) then
!!      !  write (gol,*) 'yyy   r ispec ', ispec, tpm_sfc(:,:,ipm_no3c)(dbg_i,dbg_j); call goPr
!!      !end if
!!
!!    else
!!      ! surface fraction:
!!      tpm_sfc(:,:,ipm_no3c) = xn_adv(ispec,:,:,nlev) * fscale &
!!                 * cfac(ispec,:,:)       ! 50 m -> 3 m
!!      !! testing ...
!!      !if ( dbg_cell ) then
!!      !  write (gol,*) 'yyy   p ispec ', ispec, tpm_sfc(:,:,ipm_no3c)(dbg_i,dbg_j); call goPr
!!      !end if
!!    end if
!!
!!    ! testing ...
!!    if ( dbg_cell ) then
!!      write (gol,*) 'yyy after tpm_sfc      = ', tpm_sfc(dbg_i,dbg_j,:); call goPr
!!      write (gol,*) 'yyy after tpm_sfc(:,:,ipm_no3c)     = ', tpm_sfc(:,:,ipm_no3c)(dbg_i,dbg_j); call goPr
!!      write (gol,*) 'yyy after f + a*no3c   = ', tpm_sfc(dbg_i,dbg_j,ipm_fine)+alfa*tpm_sfc(:,:,ipm_no3c)(dbg_i,dbg_j); call goPr
!!      write (gol,*) 'yyy after f + c + no3c = ', sum(tpm_sfc(dbg_i,dbg_j,:))+tpm_sfc(:,:,ipm_no3c)(dbg_i,dbg_j); call goPr
!!      do ispec = 1, nspec
!!        !if ( any(w(ispec,:) > 0.0) .or. (ispec == ispec_no3c) ) then
!!        !  write (gol,*) 'yyy ispec ', ispec, ' ', trim(species_adv(ispec)%name); call goPr
!!        !  do k = 1, nz
!!        !    write (gol,*) k, xn_adv(ispec,dbg_i,dbg_j,k); call goPr
!!        !  end do
!!        !end if
!!        write (gol,*) 'yyy ispec ', ispec, ' ', trim(species_adv(ispec)%name), xn_adv(ispec,dbg_i,dbg_j,nlev); call goPr
!!      end do
!!    end if
!!    
!!    ! ...............................
!    
!    !
!    ! Recompute water content:
!    !   Aero_water from AerosolCalls
!    !   PM25_water at model levels (ambient)
!    !   PM25_water_rh50 at surface (50% rh)
!    !
!    
!    ! clear:
!    deallocate( xn_25, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( xn_10, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( ispm, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( tpm_ml, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( tpm_sfc, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( factor_sfc, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!    deallocate( factor_ml, stat=status )
!    IF_NOT_OK_RETURN(status=1)
!
!    ! ok
!    status = 0
!    
!  end subroutine ObsCompInfo_ChangeProfile_FineCoarseRatio



  ! *


  
  subroutine ObsCompInfo_AdjustFineCoarseRatio( self, xn_adv, xn_adv_units, &
                                                     xn_obs, dx_obs, obs_units, status, &
                                                     verbose, bl )

    use ChemGroups_ml    , only : PMFINE_GROUP, PM10_GROUP   ! requires offset NSPEC_SHL !
    use ChemSpecs_shl_ml , only : NSPEC_SHL  ! number of short-lived tracers
    use ChemChemicals_ml , only : species_adv
    use Units_ml         , only : Units_Scale
    use MetFields_ml     , only : roa
    use MetFields_ml     , only : z_bnd    ! (nx,ny,1:nz+1) gph bounds relative to surface
    use MetFields_ml     , only : pzpbl    ! (nx,ny) boundary layer height
    use Chemfields_ml    , only : cfac               ! (i,j) 50m -> 3m factor
    use Chemfields_ml    , only : PM25_water         ! (i,j,k) model levels
    use Chemfields_ml    , only : PM25_water_rh50    ! (i,j) surface

    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(in)      ::  self
    real, intent(inout)                   ::  xn_adv(:,:,:,:)  ! (nspec_adv,lnx,lny,nz)
    character(len=*), intent(in)          ::  xn_adv_units(:)  ! (nspec_adv)
    real, intent(in)                      ::  xn_obs(:,:)      ! (lnx,lny)
    real, intent(in)                      ::  dx_obs(:,:,:)    ! (lnx,lny,nz)
    character(len=*), intent(in)          ::  obs_units
    integer, intent(out)                  ::  status
    
    logical, intent(in), optional         ::  verbose
    logical, intent(in), optional         ::  bl
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_AdjustFineCoarseRatio'
    
    ! fraction of NO3c assigned to PM2.5:
    real, parameter    ::  alfa = 0.27
    
    ! pm collections:
    integer, parameter        ::  npm = 2
    integer, parameter        ::  ipm_fine   = 1
    integer, parameter        ::  ipm_coarse = 2
    
    ! --- local -----------------------------
    
    logical               ::  verb
    logical               ::  only_bl
    integer               ::  nspec, lnx, lny, nz
    real, allocatable     ::  no3(:,:,:)  ! (lnx,lny,nz)
    real, allocatable     ::  tpm_sfc(:,:,  :)  ! (lnx,lny   ,npm)
    real, allocatable     ::  tpm_ml (:,:,:,:)  ! (lnx,lny,nz,npm)
    real, allocatable     ::  xn_25(:,:)        ! (lnx,lny)
    real, allocatable     ::  xn_10(:,:)        ! (lnx,lny)
    real, allocatable     ::  q(:,:,:)  ! (lnx,lny,npm)
    real, allocatable     ::  factor_sfc(:,:,  :)  ! (lnx,lny   ,npm)
    real, allocatable     ::  factor_ml (:,:,:,:)  ! (lnx,lny,nz,npm)
    logical, allocatable  ::  ispm(:,:)            ! (nspec_adv,npm)
    integer               ::  ipm
    integer, pointer      ::  out_group(:)
    integer               ::  out_offset
    integer               ::  ispec
    integer               ::  ispec_no3f
    integer               ::  ispec_no3c
    real                  ::  fscale
    logical               ::  needroa
    integer               ::  k
!
!      integer   :: i,j

    ! --- begin -----------------------------
    
    ! verbose?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! boundary layer only?
    only_bl = .false.
    if ( present(bl) ) only_bl = bl
    
    ! shape:
    nspec = size(xn_adv,1)
    lnx   = size(xn_adv,2)
    lny   = size(xn_adv,3)
    nz    = size(xn_adv,4)
    
    ! target obs:
    allocate( xn_25(lnx,lny), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( xn_10(lnx,lny), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for no3 fraction:
    allocate( no3(lnx,lny,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for current total fine/coarse aerosol:
    allocate( tpm_ml (lnx,lny,nz,npm), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( tpm_sfc(lnx,lny   ,npm), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for factors:
    allocate( q(lnx,lny,npm), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( factor_sfc(lnx,lny,npm), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( factor_ml (lnx,lny,nz,npm), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! Coarse nitrate should be treated special ...
    ! Fraction alfa is assigned to "fine" and (1-alfa) to "coarse",
    ! thus fine/coarse ratio "within" no3c is fixed.
    ! Here add fraction alfa from course to fine mode,
    ! and at the end assign part of fine mode to coarse mode again.
    ! Find indices of nitrate aerosols:
    ispec_no3f = find_index( 'NO3_F', species_adv(:)%name )
    ispec_no3c = find_index( 'NO3_C', species_adv(:)%name )
    ! fraction of coarse part assigned to fine:
    no3 = alfa * xn_adv(ispec_no3c,:,:,:)
    ! switch:
    xn_adv(ispec_no3f,:,:,:) = xn_adv(ispec_no3f,:,:,:) + no3
    xn_adv(ispec_no3c,:,:,:) = xn_adv(ispec_no3c,:,:,:) - no3

    ! storage for flags to know if spec is fine or coarse:
    allocate( ispm(nspec,npm), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! init with no contribution:
    ispm = .false.    
    ! loop over pm collections:
    do ipm = 1, npm
      ! use group defined in 'ChemGroups_ml':
      if ( ipm == ipm_fine ) then
        out_group => PMFINE_GROUP
      else
        out_group => PM10_GROUP
      end if
      ! offset:
      out_offset = NSPEC_SHL
      ! loop over group members:
      do k = 1, size(out_group)
        ! extract index of specie:
        ispec = out_group(k) - out_offset
        ! store relative weight in fine and coarse parts:
        if ( ipm == ipm_fine ) then
          ! contributes to fine part:
          ispm(ispec,ipm) = .true.
        !  
        else if ( ipm == ipm_coarse ) then
          ! only if exclusevely coarse and not coarse nitrate:
          ispm(ispec,ipm) = .not. ispm(ispec,ipm_fine)
        end if
        !
      end do ! group elements
    end do  ! total fine and coarse
    
    !! info ...
    !write (gol,*) 'rrr weights: '; call goPr
    !do ispec = 1, nspec
    !  if ( any( ispm(ispec,:) ) ) then
    !    write (gol,'("rrr ",2l2," ",a)') ispm(ispec,:), trim(species_adv(ispec)%name); call goPr
    !  end if
    !end do

    ! init sums:
    tpm_ml  = 0.0
    tpm_sfc = 0.0
    ! loop over fine (=1), coarse-without-no3c (=2), and no3c (=3):
    do ipm = 1, npm
      !! testing ...
      !if ( dbg_cell ) then
      !  write (gol,*) 'yyy ipm = ', ipm; call goPr
      !end if
      ! loop over species:
      do ispec = 1, nspec
        ! skip if no contribution:
        if ( .not. ispm(ispec,ipm) ) cycle
        ! info on conversion to target units:
#ifdef with_ajs
        call Units_Scale( obs_units, ispec, fscale, &
                                 needroa=needroa, status=status )
        IF_NOT_OK_RETURN(status=1)
#else
        call Units_Scale( obs_units, ispec, fscale, &
                                 needroa=needroa )
#endif
        !! testing ...
        !write (gol,*) 'qqq ', ipm, ' ', obs_units, ispec, ' ', trim(species_adv(ispec)%name), fscale, needroa; call goPr
        ! convert using air density?
        if ( needroa ) then
          ! add contribution, multiply with air density:
          tpm_ml(:,:,:,ipm) = tpm_ml(:,:,:,ipm) + &
                               xn_adv(ispec,:,:,:) &
                               * fscale &
                               * roa(:,:,:,1)    ! density
          ! surface fraction:
          tpm_sfc(:,:,ipm)  = tpm_sfc(:,:,ipm) + &
                               xn_adv(ispec,:,:,nz) &
                               * fscale &
                               * roa(:,:,nz,1)    &  ! density
                               * cfac(ispec,:,:)       ! 50 m -> 3 m
          !! testing ...
          !if ( dbg_cell ) then
          !  write (gol,*) 'yyy   ispec ', ispec, xn_adv(ispec,dbg_i,dbg_j,nz) * w(ispec,ipm) &
          !                                       * fscale * roa(dbg_i,dbg_j,nz,1) * cfac(ispec,dbg_i,dbg_j); call goPr
          !end if
                               
        else
          ! add contribution:
          tpm_ml(:,:,:,ipm) = tpm_ml(:,:,:,ipm) + &
                               xn_adv(ispec,:,:,:) &
                               * fscale
          ! surface fraction:
          tpm_sfc(:,:,ipm)  = tpm_sfc(:,:,ipm) + &
                               xn_adv(ispec,:,:,nz) &
                               * fscale &
                               * cfac(ispec,:,:)       ! 50 m -> 3 m
          !! testing ...
          !if ( dbg_cell ) then
          !  write (gol,*) 'yyy   ispec ', ispec, xn_adv(ispec,dbg_i,dbg_j,nz) * w(ispec,ipm) &
          !                                       * fscale * cfac(ispec,dbg_i,dbg_j); call goPr
          !end if
        end if
      end do  ! spec
    end do  ! pm collections
    ! add 3D aersosol water to fine fraction:
    tpm_ml (:,:,:,ipm_fine) = tpm_ml (:,:,:,ipm_fine) + PM25_water
    tpm_sfc(:,:  ,ipm_fine) = tpm_sfc(:,:  ,ipm_fine) + PM25_water_rh50
    !! testing ...
    !if ( dbg_cell ) then
    !  write (gol,*) 'yyy   rh50           = ', PM25_water_rh50(dbg_i,dbg_j); call goPr
    !end if
    
!    !! testing ...
!    !if ( dbg_cell ) then
!    !  write (gol,*) 'yyy before tpm fine,course   = ', tpm_sfc(dbg_i,dbg_j,:); call goPr
!    !  write (gol,*) 'yyy before tpm fine+coarse   = ', sum(tpm_sfc(dbg_i,dbg_j,:)); call goPr
!    !  write (gol,*) 'yyy before tpm f + a*no3c    = ', tpm_sfc(dbg_i,dbg_j,1)+alfa*tpm_sfc(:,:,3)(dbg_i,dbg_j); call goPr
!    !  write (gol,*) 'yyy before tpm f + c + no3c  = ', sum(tpm_sfc(dbg_i,dbg_j,:)) + tpm_sfc(:,:,3)(dbg_i,dbg_j); call goPr
!    !end if
!
!    !! testing ..
!    !write (gol,*) 'ppp0 tpm = ', tpm_sfc(1,1,:), sum(tpm_sfc(1,1,:)); call goPr
    
    ! analyzed total pm25 at the surcface:  xn_obs + dx_obs(:,:,nz)
    ! in theory this could not exceed total pm10,
    ! but total pm10 was analyzed independently,
    ! and sometimes observations are found with pm25 > pm10 ;
    ! limit the analyzed field to the current total 
    ! (which should remain unchanged in this update):
    xn_10 = sum( tpm_sfc, dim=3 )
    xn_25 = min( max( 0.0, xn_obs + dx_obs(:,:,nz) ), xn_10 )
    !! without the water part:
    !xn_ana = min( xn_obs + dx_obs - PM25_water, tpm(:,:,:,1)+tpm(:,:,:,ipm_coarse) )
    
!    !! testing ...
!    !if ( dbg_cell ) then
!    !  write (gol,*) 'yyy xn_obs, dx_obs  = ', xn_obs(dbg_i,dbg_j), dx_obs(dbg_i,dbg_j,nz); call goPr
!    !  write (gol,*) 'yyy xn_obs + dx_obs = ', xn_obs(dbg_i,dbg_j) + dx_obs(dbg_i,dbg_j,nz); call goPr
!    !  write (gol,*) 'yyy xn_ana          = ', xn_ana(dbg_i,dbg_j); call goPr
!    !end if

    ! pmf transformed into pmc or reverse?
    where ( xn_25 < tpm_sfc(:,:,ipm_fine) )
    
      !
      !   |
      !   |  pmf | pmc |
      !   |      |     |
      !   +---*--o-----o---------
      !      xn  
      !
      ! fraction of fine mode remaining fine:
      !    qf = xn / pmf
      !
      ! new amounts:
      !   pmf  :=     qf  pmf
      !   pmc  :=  (1-qf) pmf + pmc
      !
      
      ! fraction of fine mode remaining fine:
      where ( tpm_sfc(:,:,ipm_fine) > 0.0 )
        q(:,:,ipm_fine  ) = max( 0.01, xn_25 / tpm_sfc(:,:,ipm_fine) )
      elsewhere
        q(:,:,ipm_fine  ) = 0.01
      endwhere
      ! all coarse mode remains coarse:
      q(:,:,ipm_coarse) = 1.0
      
    elsewhere

      !
      !   |
      !   |  pmf | pmc |
      !   |      |     |
      !   +------o--*--o---------
      !            xn  
      !
      ! fraction of coarse mode remaining coarse:
      !    qf = (pmf+pmc-xn) / pmc
      !
      ! new amounts:
      !   pmf  :=  pmf + (1-qc) pmc
      !   pmc  :=           qc  pmc
      !
      
      ! all fine aerosol remains fine:
      q(:,:,ipm_fine  ) = 1.0
      ! fraction of coarse mode remaining coarse:
      where ( tpm_sfc(:,:,ipm_coarse) > 0.0 )
        q(:,:,ipm_coarse) = max( 0.01, ( tpm_sfc(:,:,ipm_fine) + tpm_sfc(:,:,ipm_coarse) - xn_25 ) / tpm_sfc(:,:,ipm_coarse) )
      elsewhere
        q(:,:,ipm_coarse) = 0.01
      endwhere
      
    endwhere

    ! scale factors for fine mode at surface level:
    where ( tpm_sfc(:,:,ipm_fine) > 0.0  )
      factor_sfc(:,:,ipm_fine  ) = (    q(:,:,ipm_fine) *tpm_sfc(:,:,ipm_fine) + (1.0-q(:,:,ipm_coarse))*tpm_sfc(:,:,ipm_coarse) ) / tpm_sfc(:,:,ipm_fine  )
    elsewhere
      factor_sfc(:,:,ipm_fine  ) = 1.0
    endwhere
    ! scale factors for coarse mode at surface level:
    where ( tpm_sfc(:,:,ipm_coarse) > 0.0 )
      factor_sfc(:,:,ipm_coarse) = ( (1-q(:,:,ipm_fine))*tpm_sfc(:,:,ipm_fine) +      q(:,:,ipm_coarse) *tpm_sfc(:,:,ipm_coarse) ) / tpm_sfc(:,:,ipm_coarse)
    elsewhere
      factor_sfc(:,:,ipm_coarse) = 1.0
    endwhere

    ! loop over layers:
    do k = 1, nz

      ! scale factors for fine mode in this layer:
      where ( tpm_ml(:,:,k,ipm_fine) > 0.0  )
        factor_ml(:,:,k,ipm_fine  ) = (    q(:,:,ipm_fine) *tpm_ml(:,:,k,ipm_fine) + (1.0-q(:,:,ipm_coarse))*tpm_ml(:,:,k,ipm_coarse) ) / tpm_ml(:,:,k,ipm_fine  )
      elsewhere
        factor_ml(:,:,k,ipm_fine  ) = 1.0
      endwhere
      ! scale factors for coarse mode in this layer:
      where ( tpm_ml(:,:,k,ipm_coarse) > 0.0 )
        factor_ml(:,:,k,ipm_coarse) = ( (1-q(:,:,ipm_fine))*tpm_ml(:,:,k,ipm_fine) +      q(:,:,ipm_coarse) *tpm_ml(:,:,k,ipm_coarse) ) / tpm_ml(:,:,k,ipm_coarse)
      elsewhere
        factor_ml(:,:,k,ipm_coarse) = 1.0
      endwhere
      
      ! boundary layer only?
      if ( only_bl ) then
        ! reset if bottom above blh:
        where ( z_bnd(:,:,k+1) >= pzpbl )
          factor_ml(:,:,k,ipm_fine  ) = 1.0
          factor_ml(:,:,k,ipm_coarse) = 1.0
        endwhere
      end if
      
    end do ! levels
   
    ! loop over fine/coarse:
    do ipm = 1, npm
      ! loop over all species:
      do ispec = 1, nspec
        ! skip if no contribution to pm:
        if ( .not. ispm(ispec,ipm) ) cycle
        !! testing ...
        !if ( dbg_cell ) then
        !  write (gol,*) 'yyy ispec ', ispec; call goPr
        !  !write (gol,*) 'yyy   wf ml  ', w(ispec,ipm_fine) * factor_ml(:,:,ipm_fine)(dbg_i,dbg_j,nz), w(ispec,ipm_coarse) * factor_ml(:,:,ipm_coarse)(dbg_i,dbg_j,nz); call goPr
        !  write (gol,*) 'yyy   wf ml  ', w(ispec,ipm_fine) * factor_ml(:,:,ipm_fine)(dbg_i,dbg_j,nz) + w(ispec,ipm_coarse) * factor_ml(:,:,ipm_coarse)(dbg_i,dbg_j,nz); call goPr
        !  write (gol,*) 'yyy   xn_adv ', xn_adv(ispec,dbg_i,dbg_j,nz); call goPr
        !end if
        ! apply factor fields:
        xn_adv(ispec,:,:,:) = xn_adv(ispec,:,:,:) * factor_ml(:,:,:,ipm)
        !! testing ...
        !if ( dbg_cell ) then
        !  !write (gol,*) 'yyy   xn_adv ', xn_adv(ispec,dbg_i,dbg_j,nz); call goPr
        !  write (gol,*) 'xxx xn_adv(',ispec,')=',xn_adv(ispec,dbg_i,dbg_j,:); call goPr
        !end if
      end do ! ispec
    end do ! ipm

    ! also apply (fine) factor to aersol water;
    ! better would be to recompute water content after changing the aerosol load ...
    PM25_water      = PM25_water      * factor_ml (:,:,:,ipm_fine)   ! model levels
    PM25_water_rh50 = PM25_water_rh50 * factor_sfc(:,:  ,ipm_fine)   ! sfc
    
!      ! testing ...
!      if ( dbg_cell ) then
!        write (gol,*) 'xxx PM25_water      = ', PM25_water(dbg_i,dbg_j,:); call goPr
!        write (gol,*) 'xxx PM25_water_rh50 = ', PM25_water_rh50(dbg_i,dbg_j); call goPr
!      end if

    ! part of fine mode that could be assigned to coarse mode,
    ! remove max 99% of the fine mode:
    no3 = min( alfa * xn_adv(ispec_no3c,:,:,:), 0.99 * xn_adv(ispec_no3f,:,:,:) )
    ! switch:
    xn_adv(ispec_no3f,:,:,:) = xn_adv(ispec_no3f,:,:,:) - no3
    xn_adv(ispec_no3c,:,:,:) = xn_adv(ispec_no3c,:,:,:) + no3
    
    ! clear:
    deallocate( xn_25, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( xn_10, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( no3, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( ispm, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( tpm_ml, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( tpm_sfc, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( q, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( factor_sfc, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( factor_ml, stat=status )
    IF_NOT_OK_RETURN(status=1)

    !! testing ...
    !write (gol,'("break")'); call goErr
    !TRACEBACK; status=1; return

    ! ok
    status = 0
    
  end subroutine ObsCompInfo_AdjustFineCoarseRatio



  ! *


  !
  ! Input:
  !   xn_adv    :  model traces (all)
  !   xn_obs    :  simulated observation field
  !   dx_obs    :  analysis increments
  !
  ! Change tracers in "xn_adv" following ratio:
  !    (xn_obs+dx_obs)/xn_obs
  !
  ! The 3D increment is used to ensure that factors applied
  ! to the levels decrease to 1 at the top.
  !
  ! If bl=.true. layers are changed from bottom (nz) up to layer
  ! including boundary layer top ("pzpbl" from "MetFields_ml")
  !
  
  subroutine ObsCompInfo_ChangeML( self, xn_adv, xn_adv_units, &
                                       xn_obs, dx_obs, obs_units, status, &
                                       maxratio, verbose, bl )

    use MetFields_ml, only : z_bnd    ! (nx,ny,1:nz+1) gph bounds relative to surface
    use MetFields_ml, only : pzpbl    ! (nx,ny) boundary layer height

    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(in)      ::  self
    real, intent(inout)                   ::  xn_adv(:,:,:,:)  ! (nspec_adv,lnx,lny,nz)
    character(len=*), intent(in)          ::  xn_adv_units(:)  ! (nspec_adv)
    real, intent(in)                      ::  xn_obs(:,:,:)    ! (lnx,lny,nz)
    real, intent(in)                      ::  dx_obs(:,:,:)    ! (lnx,lny,nz)
    character(len=*), intent(in)          ::  obs_units
    integer, intent(out)                  ::  status
    
    real, intent(in), optional            ::  maxratio
    logical, intent(in), optional         ::  verbose
    logical, intent(in), optional         ::  bl
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_ChangeML'
    
    ! --- local -----------------------------
    
    logical               ::  verb
    logical               ::  only_bl
    integer               ::  nspec, lnx, lny, nz
    real, allocatable     ::  ratio(:,:)  ! (lnx,lny)
    integer               ::  k
    integer               ::  is, ispec

    ! --- begin -----------------------------
    
    ! verbose?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! boundary layer only?
    only_bl = .false.
    if ( present(bl) ) only_bl = bl
    
    ! shape:
    nspec = size(xn_adv,1)
    lnx   = size(xn_adv,2)
    lny   = size(xn_adv,3)
    nz    = size(xn_adv,4)
    
    !! testing ...
    !i = 1
    !j = 1
    !write (gol,*) 'xxx testing pbl'; call goPr
    !write (gol,*) 'xxx z_bnd ', lbound(z_bnd), ubound(z_bnd); call goPr
    !do k = 1, nz
    !  write (gol,*) 'xxx   layer top,bottom : ', z_bnd(i,j,k:k+1); call goPr
    !end do
    !write (gol,*) 'xxx pbl', pzpbl(i,j); call goPr
    !TRACEBACK; status=1; return
    
    ! storage for relative change:
    allocate( ratio(lnx,lny), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over levels, skip top layer (boundary conditions):
    do k = 2, nz

      ! fill ratio, or set to unity if undefined (zero tracer sum):
      where( xn_obs(:,:,k) > 0.0 )
        ratio = ( xn_obs(:,:,k) + dx_obs(:,:,k) ) / xn_obs(:,:,k)
      elsewhere
        ratio = 1.0
      end where
      
      ! boundary layer only?
      if ( only_bl ) then
        ! reset if bottom is above boundary layer height:
        where ( z_bnd(:,:,k+1) >= pzpbl )
          ratio = 1.0
        endwhere
      end if

      ! use minimum ratio to avoid that all tracers are removed,
      ! this would give errors from aerosol equilibrium routines:
      ratio = max( 0.01, ratio )
      ! maximum?
      if ( present(maxratio) ) ratio = min( ratio, maxratio )

      ! loop over contributing species:
      do is = 1, self%nspec
        ! global index:
        ispec = self%ispec(is)
        ! check units ..
        if ( trim(xn_adv_units(ispec)) /= trim(obs_units) ) then
          write (gol,'("xn_adv units `",a,"` while obs units `",a,"`")') trim(xn_adv_units(ispec)), trim(obs_units); call goErr
        TRACEBACK; status=1; return
        end if

        ! apply:
        xn_adv(ispec,:,:,k) = xn_adv(ispec,:,:,k) * ratio
    
      end do ! specs

    end do ! levels

    ! clear:
    deallocate( ratio, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine ObsCompInfo_ChangeML



  ! *


  !
  ! Input:
  !   xn_adv    :  model traces (all)
  !   xn_obs    :  3D simulated observation field
  !   dx_obs    :  2D analysis increments
  !
  ! Change tracers in "xn_adv" following vcd profile:
  !    xn_obs(:) + vcd(:)/sum(vcd) dx_obs
  !
  ! No changes in layers above 200 hPa.
  !
  
  subroutine ObsCompInfo_ChangeVCD( self, xn_adv, xn_adv_units, &
                                       xn_obs, dx_obs, obs_units, status, &
                                       verbose )


    use PhysicalConstants_ml, only : Avog   ! 6.02e23 mlc/mole
    use PhysicalConstants_ml, only : AtwAir ! 28.964 (g air)/mole
    use GridValues_ml       , only : A_bnd, B_bnd
    use MetFields_ml        , only : z_bnd    ! (nx,ny,1:nz+1) gph bounds relative to surface
    use MetFields_ml        , only : roa
    use MetFields_ml        , only : ps

    ! --- in/out ----------------------------
    
    class(T_ObsCompInfo), intent(in)      ::  self
    real, intent(inout)                   ::  xn_adv(:,:,:,:)  ! (nspec_adv,lnx,lny,nz)
    character(len=*), intent(in)          ::  xn_adv_units(:)  ! (nspec_adv)
    real, intent(in)                      ::  xn_obs(:,:,:)    ! (lnx,lny,nz)
    real, intent(in)                      ::  dx_obs(:,:)      ! (lnx,lny)
    character(len=*), intent(in)          ::  obs_units
    integer, intent(out)                  ::  status
    
    logical, intent(in), optional         ::  verbose
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/ObsCompInfo_ChangeVCD'
    
    ! --- local -----------------------------
    
    logical               ::  verb
    integer               ::  nspec, lnx, lny, nz
    real, allocatable     ::  vcd(:,:,:)        ! (lnx,lny,nz)
    real, allocatable     ::  vc(:,:)           ! (lnx,lny)
    real, allocatable     ::  fraction(:,:)     ! (lnx,lny)
    real, allocatable     ::  ph(:,:,:)         ! (lnx,lny,nz)
    integer               ::  ilev
    integer               ::  i, j
    integer               ::  is, ispec

    ! --- begin -----------------------------
    
    ! verbose?
    verb = .false.
    if ( present(verbose) ) verb = verbose
    
    ! shape:
    nspec = size(xn_adv,1)
    lnx   = size(xn_adv,2)
    lny   = size(xn_adv,3)
    nz    = size(xn_adv,4)
    
    !! testing ...
    !i = 1
    !j = 1
    !write (gol,*) 'xxx testing pbl'; call goPr
    !write (gol,*) 'xxx z_bnd ', lbound(z_bnd), ubound(z_bnd); call goPr
    !do k = 1, nz
    !  write (gol,*) 'xxx   layer top,bottom : ', z_bnd(i,j,k:k+1); call goPr
    !end do
    !write (gol,*) 'xxx pbl', pzpbl(i,j); call goPr
    !TRACEBACK; status=1; return
    
    ! storage for column densities:
    allocate( vcd(lnx,lny,nz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! total:
    allocate( vc(lnx,lny), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for fractions
    allocate( fraction(lnx,lny), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! half leel pressure:
    allocate( ph(lnx,lny,0:nz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill using hybride coeff;
    ! first record of surface pressure is said to be the correct one ..
    do j = 1, lny
      do i = 1, lnx
        ph(i,j,:) = A_bnd + B_bnd * ps(i,j,1)
      end do
    end do

    ! [1e15 mlc/cm2]/m = 
    !   (mole tracer)/(mole air)   ! input units (ppv)
    !   * (kg air)/m3              ! applied later
    !   * 1e-15   mlc/(mole tracer)  (mole air)/(kg air)   m2/cm2  
    vcd = xn_obs * & ! ppv
          roa(:,:,:,1) * & ! kg/m3
          1e-15 *      Avog        /    (AtwAir*1e-3)    *  1e-4  ! <units>/m

    ! loop over levels:
    do ilev = 1, nz
      ! multiply with layer thickness:
      !                 <units>/m    *                   m
      vcd(:,:,ilev) = vcd(:,:,ilev)  * ( z_bnd(:,:,ilev) - z_bnd(:,:,ilev+1) )  ! <units>
    end do
    
    ! skip above tropo:
    where( ph(:,:,1:nz) < 200.0e2 )
      vcd = 0.0  ! <units>
    end where
    
    ! sum:
    vc = sum(vcd,dim=3)  ! <units>
    
    !write (gol,*) 'vvv1 dx_obs = ', dx_obs(60,1); call goPr
    !write (gol,*) 'vvv1 vc = ', vc(60,1); call goPr
    
    ! loop over levels, skip top layer (boundary conditions):
    do ilev = 2, nz
    
      ! fractions:
      where ( vc > 0.0 )
        fraction = vcd(:,:,ilev) / vc
      elsewhere
        fraction = 1.0
      end where

      ! loop over contributing species:
      do is = 1, self%nspec
        ! global index:
        ispec = self%ispec(is)
        ! check units ..
        if ( trim(xn_adv_units(ispec)) /= trim(obs_units) ) then
          write (gol,'("xn_adv units `",a,"` while obs units `",a,"`")') trim(xn_adv_units(ispec)), trim(obs_units); call goErr
          TRACEBACK; status=1; return
        end if
        
        !write (gol,*) 'vvv1 lev ', ilev, fraction(60,1), xn_adv(ispec,60,1,ilev), &
        !          xn_adv(ispec,60,1,ilev) + fraction(60,1) * dx_obs(60,1); call goPr

        ! apply:
        xn_adv(ispec,:,:,ilev) = max( 0.0, xn_adv(ispec,:,:,ilev) + fraction * dx_obs )
    
      end do ! specs

    end do ! levels
    
    !print *, 'break'
    !TRACEBACK; status=1; return

    ! clear:
    deallocate( vcd, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( vc, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( fraction, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( ph, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine ObsCompInfo_ChangeVCD



  !=========================================================================
  !===
  !=== data sets
  !===
  !=========================================================================

  
  subroutine DA_Obs_Init( status )
  
    use Io_ml           , only : IO_NML
    use SmallUtils_ml   , only : find_index
    use ChemChemicals_ml, only : species_adv
    use ChemGroups_ml   , only : PPM25_GROUP, PPM10_GROUP, PMFINE_GROUP, PM10_GROUP
    use ChemSpecs_shl_ml, only : NSPEC_SHL  ! number of short-lived tracers, offseet to PM groups
    
    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/DA_Obs_Init'

    ! --- namelists -------------------------
    
    ! name list that defines the actual values:
    namelist /OBSERVATIONS/ nObsData,obsData

    ! --- local -----------------------------
    
    integer             ::  iObsData
    integer             ::  ispec
    character(len=16)   ::  ObsComp(nObsDataMax)
    character(len=16)   ::  ObsBUnit(nObsDataMax)
    character(len=16)   ::  ObsType(nObsDataMax)
    character(len=1024) ::  ObsBfile(nObsDataMax)
    real                ::  Bsigfac(nObsDataMax)
    logical             ::  ObsFeedback(nObsDataMax)
    character(len=32)   ::  ObsFeedbackType(nObsDataMax)
    integer             ::  iObsComp
    integer             ::  i

    ! --- begin -----------------------------
    
    ! * configuration by namelist
    
    !
    ! Read definition of data sets to be used,
    ! variables and namelist defined as module variables:
    !
    !   nObsData                ! number of datasets
    !   obsData(1:nObsData)%..  ! types with dataset properties
    !
    ! back to start:
    rewind(IO_NML)
    ! read using namelist:
    read( unit=IO_NML, nml=OBSERVATIONS, iostat=status )
    if ( status /= 0 ) then
      write (gol,'("reading namelist `OBSERVATIONS` from: config_emep.nml")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! info ...
    write (gol,'(a,": number of observation data sets : ",i0)') rname, nObsData; call goPr
    
    ! * kernels
    
    ! potential storage for kernels:
    allocate( ObsKernels(nObsData), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do iObsData = 1, nObsData
      ! init:
      call ObsKernels(iObsData)%Init( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! * flags, indices, ..
    
    ! loop over datasets:
    do iObsData = 1, nObsData
    
      ! from model level or at surface?
      if ( obsData(iObsData)%deriv == '' ) then
        ! default is surface ..
        obsData(iObsData)%deriv = 'surface'
      end if
      ! add tag:
      write (obsData(iObsData)%tag,"(2(A,1X),'(',A,')')") &
               trim(obsData(iObsData)%name), trim(obsData(iObsData)%deriv), trim(obsData(iObsData)%unit)
      
      ! * obs.repr.error
                
      ! how is observation representation error set?
      select case ( trim(obsData(iObsData)%error_type) )

        !~ original:
        case ( 'fraction', 'retrieval' )
          ! no special tasks ...

        !~ values per site:
        case ( 'estim' )
          ! not initialized yet?
          if ( .not. allocated(ObsError_Lists) ) then
            ! storage:
            allocate( ObsError_Lists(nObsData), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! loop over elements:
            do i = 1, nObsData
              ! basic setup:
              call ObsError_Lists(i)%Init( status )
              IF_NOT_OK_RETURN(status=1)
            end do
          end if
          ! read data files:
          call ObsError_Lists(iObsData)%ReadFiles( (/ obsData(iObsData)%error_obs , &
                                                      obsData(iObsData)%error_spat /), status )
          IF_NOT_OK_RETURN(status=1)

        !~
        case default
          write (gol,'("unsupported obs.repr.error type `",a,"`")') trim(obsData(iObsData)%error_type); call goErr
          TRACEBACK; status=1; return
      end select

    end do ! obs data sets

    !
    ! Tracer info defined in "CM_ChemSpecs_ml.f90", module "ChemChemicals_ml",
    ! array "species(NSPEC_TOT)" of type "Chemical" with fields:
    !    character   :  name
    !
    ! Species indices 's' in  'xn_adv(s,x,y,z)' are in species_adv:
    !    xn_adv( 1,x,y,z)  :  species_adv( 1)%name = 'O3'
    !    xn_adv( 2,x,y,z)  :  species_adv( 2)%name = 'NO2'
    !    xn_adv(53,x,y,z)  :  species_adv(53)%name = 'SO2'
    !
    ! info ...
    write (gol,'(a,": components observed by dataset :")') rname; call goPr
    ! init counter:
    nObsComp = 0
    ! loop over datasets:
    do iObsData = 1, nObsData
      ! info ..
      write (gol,'(a,":   data set ",i0)') rname, iObsData; call goPr
      write (gol,'(a,":     component : ",a)') rname, trim(obsData(iObsData)%name); call goPr
      write (gol,'(a,":     unit      : ",a)') rname, trim(obsData(iObsData)%unit); call goPr
      write (gol,'(a,":     type      : ",a)') rname, trim(obsData(iObsData)%deriv); call goPr
      ! compare with current:
      iObsComp = -999
!TESTING: all individually
!      do i = 1, nObsComp
!        ! match?
!        if ( trim(obsData(iObsData)%name) == trim(ObsComp(i)) ) then
!          ! store index:
!          iObsComp = i
!          ! check ...
!          if ( trim(ObsBUnit(iObsComp)) /= trim(obsData(iObsData)%Bunit) ) then
!            write (gol,'("could not handle observations of the same component with different Bunits yet")'); call goErr
!            write (gol,'("obsdata Bunit   : ",a)') trim(obsData(iObsData)%Bunit); call goErr
!            write (gol,'("collected Bunit : ",a)') trim(ObsBUnit(iObsComp)); call goErr
!            TRACEBACK; status=1; return
!          end if
!          ! check ...
!          if ( trim(ObsBfile(iObsComp)) /= trim(obsData(iObsData)%Bfile) ) then
!            write (gol,'("could not handle observations of the same component for different Bfiles yet")'); call goErr
!            TRACEBACK; status=1; return
!          end if
!          ! leave:
!          exit
!        end if
!      end do
      ! not defined yet?
      if ( iObsComp < 0 ) then
        ! increase counter:
        nObsComp = nObsComp + 1
        ! store name etc:
        ObsComp        (nObsComp) = trim(obsData(iObsData)%name)
        ObsType        (nObsComp) = trim(obsData(iObsData)%deriv)
        ObsBfile       (nObsComp) = trim(obsData(iObsData)%Bfile)
        ObsBUnit       (nObsComp) = trim(obsData(iObsData)%Bunit)
        Bsigfac        (nObsComp) = obsData(iObsData)%Bsigfac
        ObsFeedback    (nObsComp) = obsData(iObsData)%feedback
        ObsFeedbackType(nObsComp) = obsData(iObsData)%feedback_type
        ! set index:
        iObsComp = nObsComp
      end if
      ! store index of observed component with data set:
      obsData(iObsData)%iObsComp = iObsComp
      ! info ...
      write (gol,'(a,":   assigned to ",a," (",i0,")")') &
               rname, trim(ObsComp(iObsComp)), iObsComp; call goPr

    end do ! data sets
    
    ! check ...
    if ( nObsComp == 0 ) then
      write (gol,'("no observed components found, something wrong in configuration?")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! info ...
    write (gol,'(a,": observed components:")') rname; call goPr
    do iObsComp = 1, nObsComp
      write (gol,'(a,":   ",i0," ",a)') rname, iObsComp, trim(ObsComp(iObsComp)); call goPr
    end do
    
    ! storage:
    allocate( ObsCompInfo(nObsComp), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! info ...
    write (gol,'(a,": observed component originating species:")') rname; call goPr
    ! loop:
    do iObsComp = 1, nObsComp
      ! info ...
      write (gol,'(a,":   species name      : ",a)') rname, trim(ObsComp(iObsComp)); call goPr
      write (gol,'(a,":   target unit (B)   : ",a)') rname, trim(ObsBUnit(iObsComp)); call goPr
      write (gol,'(a,":   observation type  : ",a)') rname, trim(ObsType(iObsComp)); call goPr
      ! switch:
      select case ( trim(ObsComp(iObsComp)) )
        !
        !!~ total primary fine pm:
        !case ( 'PPM25' )
        !  ! use group defined in 'CM_ChemGroups_ml',
        !  ! for index in xn_adv remove the offset for short-lived species:
        !  call ObsCompInfo(iObsComp)%Init( trim(ObsComp(iObsComp)), trim(ObsBUnit(iObsComp)), &
        !                                   trim(ObsType(iObsComp)), trim(ObsBfile(iObsComp)), &
        !                                   PPM25_GROUP-NSPEC_SHL, status, &
        !                                   Bsigfac=Bsigfac(iObsComp) )
        !  IF_NOT_OK_RETURN(status=1)
        !  !
        !!~ total primary coarse pm:
        !case ( 'PPM10' )
        !  ! use group defined in 'CM_ChemGroups_ml',
        !  ! for index in xn_adv remove the offset for short-lived species:
        !  call ObsCompInfo(iObsComp)%Init( trim(ObsComp(iObsComp)), trim(ObsBUnit(iObsComp)), &
        !                                   trim(ObsType(iObsComp)), trim(ObsBfile(iObsComp)), &
        !                                   PPM10_GROUP-NSPEC_SHL, status, &
        !                                   Bsigfac=Bsigfac(iObsComp) )
        !  IF_NOT_OK_RETURN(status=1)
        !
        !~ total fine pm:
        case ( 'PM25' )
          ! use group defined in 'CM_ChemGroups_ml',
          ! for index in xn_adv remove the offset for short-lived species;
          ! also add part of coarse nitrate and aerosol water:
          call ObsCompInfo(iObsComp)%Init( trim(ObsComp(iObsComp)), trim(ObsBUnit(iObsComp)), &
                                           trim(ObsType(iObsComp)), trim(ObsBfile(iObsComp)), &
                                           PMFINE_GROUP-NSPEC_SHL, status, &
                                           Bsigfac=Bsigfac(iObsComp), &
                                           with_no3c_frac=.true., with_pmwater=.true., &
                                           feedback=ObsFeedback(iObsComp), &
                                           feedback_type=ObsFeedbackType(iObsComp) )
          IF_NOT_OK_RETURN(status=1)
          !
        !~ total coarse pm:
        case ( 'PM10' )
          ! use group defined in 'CM_ChemGroups_ml',
          ! for index in xn_adv remove the offset for short-lived species;
          ! also add aerosol water:
          call ObsCompInfo(iObsComp)%Init( trim(ObsComp(iObsComp)), trim(ObsBUnit(iObsComp)), &
                                           trim(ObsType(iObsComp)), trim(ObsBfile(iObsComp)), &
                                           PM10_GROUP-NSPEC_SHL, status, &
                                           Bsigfac=Bsigfac(iObsComp), &
                                           with_pmwater=.true., &
                                           feedback=ObsFeedback(iObsComp), &
                                           feedback_type=ObsFeedbackType(iObsComp) )
          IF_NOT_OK_RETURN(status=1)
          !
        !~ expected a single model species ...
        case default
          ! find index:
          ispec = find_index( trim(ObsComp(iObsComp)), species_adv(:)%name )
          ! check ...
          if ( ispec < 1 ) then
            write (gol,'(a,": observed species not found:",a)') rname, trim(obsData(iObsData)%name); call goErr
            TRACEBACK; status=1; return
          end if
          ! init mapping:
          call ObsCompInfo(iObsComp)%Init( trim(ObsComp(iObsComp)), trim(ObsBUnit(iObsComp)), &
                                           trim(ObsType(iObsComp)), trim(ObsBfile(iObsComp)), &
                                           (/ispec/), status, &
                                           Bsigfac=Bsigfac(iObsComp), &
                                           feedback=ObsFeedback(iObsComp), &
                                           feedback_type=ObsFeedbackType(iObsComp) )
          IF_NOT_OK_RETURN(status=1)
      end select
      ! info ...
      call ObsCompInfo(iObsComp)%Show( status )
      IF_NOT_OK_RETURN(status=1)
    end do ! observed components
    
    ! *

    ! ok
    status = 0
    
  end subroutine DA_Obs_Init
  
  
  ! ***
  
  
  subroutine DA_Obs_Done( status )

    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/DA_Obs_Done'
    
    ! --- local ----------------------------
    
    integer     ::  i
    
    ! --- begin -----------------------------
    
    ! error lists allocated?
    if ( allocated(ObsError_Lists) ) then
      ! loop over datasets:
      do i = 1, size(ObsError_Lists)
        ! done:
        call ObsError_Lists(i)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do
      ! clear:
      deallocate( ObsError_Lists, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! kernels allocated?
    if ( allocated(ObsKernels) ) then
      ! loop:
      do i = 1, size(ObsKernels)
        ! clear:
        call ObsKernels(i)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do
      ! clear:
      deallocate( ObsKernels, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine DA_Obs_Done


  !=========================================================================
  !===
  !=== observation errors
  !===
  !=========================================================================


  subroutine ObsError_List_Init( self, status )
  
    ! --- in/out -----------------------------
    
    class(T_ObsError_List), intent(out)   ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsError_List_Init'

    ! --- local -----------------------------
    
    ! --- begin -----------------------------
    
    ! no data yet:
    self%n = 0
    
    ! ok
    status = 0
    
  end subroutine ObsError_List_Init


  ! ***


  subroutine ObsError_List_Done( self, status )

    ! --- in/out -----------------------------
    
    class(T_ObsError_List), intent(inout)   ::  self
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsError_List_Done'

    ! --- begin -----------------------------
    
    ! any records?
    if ( self%n > 0 ) then
      ! clear:
      deallocate( self%err, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine ObsError_List_Done
  
  
  ! ***


  subroutine ObsError_List_ReadFiles( self, descriptions, status )
  
    use GO, only : goVarValue, goReadFromLine
    use GO, only : goGetFU
  
    ! --- in/out -----------------------------
    
    class(T_ObsError_List), intent(inout)   ::  self
    character(len=*), intent(in)            ::  descriptions(:)
    integer, intent(out)                    ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsError_List_ReadFiles'

    ! --- local -----------------------------
    
    integer                 ::  ifile
    character(len=1024)     ::  filename
    character(len=64)       ::  header_code
    character(len=64)       ::  header_err
    logical                 ::  exist
    integer                 ::  fu
    integer                 ::  nloop, iloop
    integer                 ::  nline
    character(len=1)        ::  sep
    character(len=1024)     ::  line
    integer                 ::  ncol, icol
    integer                 ::  icol_code, icol_err
    character(len=64)       ::  header
    character(len=64)       ::  scode
    character(len=64)       ::  value
    real                    ::  error
    integer                 ::  irec
    integer                 ::  i

    ! --- begin -----------------------------
    
    ! info ..
    write (gol,'("read observation representation errors ...")'); call goPr
    
    ! empty:
    self%n = 0
    ! number of values:
    self%nval = size(descriptions)
    
    ! number of loops needed, first 2 (number of records and values), then 1
    nloop = 2
    ! loop over input files:
    do ifile = 1, self%nval
    
      ! info ..
      write (gol,'("  description ",i0,"/",i0," : ",a)') &
               ifile, self%nval, trim(descriptions(ifile)); call goPr
    
      ! example description:
      !   file=/data/errors.csv,code=code,err=ozone_err
      ! extract elements:
      filename = 'None'
      call goVarValue( descriptions(ifile), ',', 'file', '=', filename, status )
      IF_NOT_OK_RETURN(status=1)
      header_code = 'code'
      call goVarValue( descriptions(ifile), ',', 'code', '=', header_code, status )
      IF_NOT_OK_RETURN(status=1)
      header_err = 'err'
      call goVarValue( descriptions(ifile), ',', 'err', '=', header_err, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! seperation character in file:
      sep = ';'
    
      ! info ..
      write (gol,'("    file             : ",a)') trim(filename); call goPr
      write (gol,'("    header for code  : ",a)') trim(header_code); call goPr
      write (gol,'("    header for error : ",a)') trim(header_err); call goPr
    
      ! check ...
      inquire( file=filename, exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found: ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if
      
      ! new file unit:
      call goGetFU( fu, status )
      IF_NOT_OK_RETURN(status=1)
      ! loops:
      do iloop = 1, nloop
      
        !! info ...
        !write (gol,'("    loop ",i0,"/",i0)') iloop, nloop; call goPr
      
        ! init counter:
        nline = 0
        
        ! open for reading:
        open( fu, file=filename, iostat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over lines:
        do
          ! read line:
          read (fu,'(a)',iostat=status) line
          if ( status < 0 ) exit ! eof
          IF_NOT_OK_RETURN(status=1)

          ! increase counter:
          nline = nline + 1
          ! only line counting?
          if ( (ifile == 1) .and. (iloop == 1) ) cycle
          
          ! header line?
          if ( nline == 1 ) then

            ! init column numbers:
            icol_code = -999
            icol_err  = -999
            ! loop over columns:
            ncol = 0
            do
              ! empty?
              if ( len_trim(line) == 0 ) exit
              ! extract next header:
              call goReadFromLine( line, header, status, sep=sep )
              IF_NOT_OK_RETURN(status=1)
              ! increase counter:
              ncol = ncol + 1
              ! match?
              if ( trim(header) == trim(header_code) ) icol_code = ncol
              if ( trim(header) == trim(header_err ) ) icol_err  = ncol
            end do
            ! check ..
            if ( icol_code < 0 ) then
              write (gol,'("no code  column `",a,"` found in file: ",a)') trim(header_code), trim(filename); call goErr
              TRACEBACK; status=1; return
            end if
            ! check ..
            if ( icol_err < 0 ) then
              write (gol,'("no error column `",a,"` found in file: ",a)') trim(header_err), trim(filename); call goErr
              TRACEBACK; status=1; return
            end if
            
          else
          
            ! loop over columns:
            do icol = 1, ncol
              ! special ?
              if ( icol == icol_code ) then
                ! read value:
                call goReadFromLine( line, scode, status, sep=sep )
                IF_NOT_OK_RETURN(status=1)
              else if ( icol == icol_err ) then
                ! read value:
                call goReadFromLine( line, error, status, sep=sep )
                IF_NOT_OK_RETURN(status=1)
              else
                ! dummy
                call goReadFromLine( line, value, status, sep=sep )
                IF_NOT_OK_RETURN(status=1)
              end if
            end do ! columns

            ! store if first file, otherwise search record:
            if ( ifile == 1 ) then
              ! record number:
              irec = nline - 1
              ! store code:
              self%scode(irec) = trim(scode)
            else
              ! search in existing records:
              irec = -999
              do i = 1, self%n
                ! match?
                if ( trim(self%scode(i)) == trim(scode) ) then
                  ! store record:
                  irec = i
                  ! leave:
                  exit
                end if
              end do ! existing records
              ! check:
              if ( irec < 0 ) then
                write (gol,'("station code `",a,"` not found in records")') trim(scode); call goErr
                TRACEBACK; status=1; return
              end if
            end if

            !! info ...
            !write (gol,'("    record ",i4," ",a," value ",i0,f8.2)') irec, trim(self%scode(irec)), ifile, error; call goPr
            ! store value:
            self%err(irec,ifile) = error

          end if  ! header or record
          
        end do ! lines
        ! close:
        close( fu, iostat=status )
        IF_NOT_OK_RETURN(status=1)
        
        ! storage?
        if ( (ifile == 1) .and. (iloop == 1) ) then
          ! number of records:
          self%n = nline - 1
          ! info ..
          write (gol,'("    allocate storage for ",i0," records")') self%n; call goPr
          ! storage:
          allocate( self%scode(self%n), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( self%err(self%n,self%nval), stat=status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
      end do ! counting loop or storage
      
      ! only one loop needed from now on:
      nloop = 1
      
    end do ! files
    
    ! info ..
    write (gol,'("  ok.")'); call goPr

    ! ok
    status = 0

  end subroutine ObsError_List_ReadFiles

  
  
  ! ***

  
  !
  ! Return total error for provided station code:
  !
  !    sqrt(sum err^2)
  !
  ! Error if some values are not defined.
  !

  subroutine ObsError_List_GetError( self, scode, error, status )
  
    ! --- in/out -----------------------------
    
    class(T_ObsError_List), intent(in)    ::  self
    character(len=*), intent(in)          ::  scode
    real, intent(out)                     ::  error
    integer, intent(out)                  ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsError_List_ReadFiles'
    
    ! default station code (median value for suburban type):
    character(len=*), parameter  ::  scode0 = 'XX9999S'

    ! --- local -----------------------------
    
    integer                 ::  i
    integer                 ::  irec
    integer                 ::  ival

    ! --- begin -----------------------------

    ! search in existing records;
    ! init target index:
    irec = -999
    ! loop over records:
    do i = 1, self%n
      ! match?
      if ( trim(self%scode(i)) == trim(scode) ) then
        ! store record:
        irec = i
        ! leave:
        exit
      end if ! match
    end do ! records

    ! not found? try default for station type "suburban" ...
    if ( irec < 0 ) then
      ! init target index:
      irec = -999
      ! loop over records:
      do i = 1, self%n
        ! match?
        if ( trim(self%scode(i)) == trim(scode0) ) then
          !! info ...
          !write (gol,'("WARNING - station code `",a,"` not found in obserr records, use default `",a,"`")') &
          !              trim(scode), trim(scode0); call goPr
          ! store record:
          irec = i
          ! leave:
          exit
        end if
      end do ! records
    end if  ! station code not found

    ! check:
    if ( irec < 0 ) then
      write (gol,'("station code `",a,"` nor default `",a,"` found in records")') &
                        trim(scode), trim(scode0); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check values ..
    if ( any( self%err(irec,:) < 0.0 ) ) then
      write (gol,'("some error values undefined for station `",a,"` :")') trim(scode); call goErr
      do ival = 1, self%nval
        write (gol,'("  element ",i0," value : ",f12.4)') ival, self%err(irec,ival); call goErr
      end do
      TRACEBACK; status=1; return
    end if
    
    ! total error:
    error = sqrt( sum( self%err(irec,:)**2 ) )
    
    ! ok
    status = 0
    
  end subroutine ObsError_List_GetError


  !=========================================================================
  !===
  !=== observation errors
  !===
  !=========================================================================


  subroutine ObsKernels_Init( self, status )
  
    ! --- in/out -----------------------------
    
    class(T_ObsKernels), intent(out)   ::  self
    integer, intent(out)               ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsKernels_Init'

    ! --- local -----------------------------
    
    ! --- begin -----------------------------
    
    ! dummy:
    self%npixel = 0
    self%nlayer = 0
    
    ! ok
    status = 0
    
  end subroutine ObsKernels_Init


  ! ***


  subroutine ObsKernels_Done( self, status )

    ! --- in/out -----------------------------
    
    class(T_ObsKernels), intent(inout)   ::  self
    integer, intent(out)                 ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsKernels_Done'

    ! --- begin -----------------------------

    ! clear storage:
    call self%Clear( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ObsKernels_Done


  ! ***


  subroutine ObsKernels_Clear( self, status )

    ! --- in/out -----------------------------
    
    class(T_ObsKernels), intent(inout)   ::  self
    integer, intent(out)                 ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsKernels_Clear'

    ! --- begin -----------------------------

    ! storage allocated?
    if ( allocated(self%kernel) ) then
      ! clear:
      deallocate( self%kernel, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%phlev, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! dummy:
    self%npixel = 0
    self%nlayer = 0

    ! ok
    status = 0

  end subroutine ObsKernels_Clear


  ! ***


  subroutine ObsKernels_Fill( self, kernel, kernel_units, phlev, phlev_units, status )

    ! --- in/out -----------------------------
    
    class(T_ObsKernels), intent(inout)   ::  self
    real, intent(in)                     ::  kernel(:,:)    ! (layer,pixels)
    character(len=*), intent(in)         ::  kernel_units
    real, intent(in)                     ::  phlev (:,:)    ! (0:layer,pixel)
    character(len=*), intent(in)         ::  phlev_units
    integer, intent(out)                 ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsKernels_Fill'

    ! --- begin -----------------------------
    
    ! clear current storage if necessary:
    call self%Clear( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! new dimensions:
    self%nlayer = size(kernel,1)
    self%npixel = size(kernel,2)
    
    ! check ..
    if ( any( shape(phlev) /= (/self%nlayer+1,self%npixel/) ) ) then
      write (gol,'("kernels shaped (",i0,",",i0,") but phlev shaped (",i0,",",i0,")")') &
               shape(kernel), shape(phlev); call goErr
      TRACEBACK; status=1; return
    end if

    ! new storage:
    allocate( self%kernel(self%nlayer,self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%phlev(0:self%nlayer,self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    self%kernel       = kernel
    self%kernel_units = kernel_units
    self%phlev        = phlev
    self%phlev_units  = phlev_units

    ! ok
    status = 0

  end subroutine ObsKernels_Fill


  !=========================================================================
  !===
  !=== observation operators
  !===
  !=========================================================================


  subroutine ObsOper_Init( self, status )
  
    use DA_ml, only : nlev

    ! --- in/out -----------------------------
    
    class(T_ObsOper), intent(out)   ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOper_Init'

    ! --- begin -----------------------------
    
    ! dummy values:
    self%obs       = -999.9
    self%xf        = -999.9
    self%xa        = -999.9
    self%sb        = -999.9
    self%sf        = -999.9
    self%innov     = -999.9
    self%obsstddev = -999.9
    self%i         = -999
    self%j         = -999
   
    ! level type:
    self%levtype   = -999
    ! storage for vertical profile:
    allocate( self%H_jac(nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! dummy values:
    self%H_jac     = -999.9
    self%H_jac1    = -999.9
    self%l         = -999
    
    ! dummy values:
    self%id        = -999
    self%iObsData  = -999
    self%iObsComp  = -999
    self%lon       = -999.9
    self%lat       = -999.9
    self%alt       = -999.9
    self%stncode   = '-'
    self%analyse   = .true.

    ! ok
    status = 0

  end subroutine ObsOper_Init


  ! ***


  subroutine ObsOper_Done( self, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOper), intent(inout)   ::  self
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOper_Done'

    ! --- begin -----------------------------
    
    ! clear:
    deallocate( self%H_jac, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ObsOper_Done


  ! ***

  ! copy values from src into own

  subroutine ObsOper_Copy( self, src, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOper), intent(inout)   ::  self
    class(T_ObsOper), intent(in)      ::  src
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOper_Copy'

    ! --- local -----------------------------

    ! --- begin -----------------------------
    
    ! copy scalar values:
    self%obs       = src%obs
    self%obsstddev = src%obsstddev
    self%xf        = src%xf
    self%xa        = src%xa
    self%sb        = src%sb
    self%sf        = src%sf
    self%innov     = src%innov
    self%H_jac     = src%H_jac
    self%H_jac1    = src%H_jac1
    self%i         = src%i
    self%j         = src%j
    self%levtype   = src%levtype
    self%l         = src%l
    self%id        = src%id
    self%iObsComp  = src%iObsComp
    self%iObsData  = src%iObsData
    self%stnid     = src%stnid
    self%lon       = src%lon
    self%lat       = src%lat
    self%alt       = src%alt
    self%stncode   = src%stncode
    self%analyse   = src%analyse

    ! ok
    status = 0

  end subroutine ObsOper_Copy


  ! ***


  subroutine ObsOper_Fill( self, &
                            !xn_adv_b, xn_adv_units, &
                            sf_obs, xn_obs, xn_obs_units, &
                            ipar, stnid, lat,lon, alt, &
                            yn, status )
  !-----------------------------------------------------------------------
  ! @description
  ! H operator mapping from model to observation space;
  ! Further: interpolate air density to observation point
  !
  ! The following 4D concentration fields are currently used;
  ! seems that content of 'xn_b' is also included in 'xn_adv',
  ! so maybe these can be combined in future ?
  !
  !   xn_adv_b :  Background concentration on local grid ("xn_adv");
  !               all tracers, used for indirect observations only.
  !               Shape: (nspec_adv,lnx,lny,nlev)
  !               Used from module 'Chemfields_ml'.
  !               Include scale factor 'FGSCALE' when using elements
  !               of this array.
  !
  !   xn_b     :  Background concentration on local grid,
  !               only tracers involved in assimilation.
  !               Shape: (lnx,lny,nlev,nchem)
  !               Passed as argument.
  !               Scale factor 'FGSCALE' is already included.
  !
  ! Output:
  !   base class:
  !     self      : innovation, H_jac, etc
  !   arguments:
  !     yn        : simulated value, ~  H_jac * xn_b
  !     rho       : air density (for unit conversions)
  !
  ! Note that the interpolation is on the local domain only!
  ! In future, halo of 1 cell should be passed to have same
  ! results when docomposition changes.
  !
  ! @author AMVB
  !-----------------------------------------------------------------------

    use PhysicalConstants_ml, only : Avog   ! 6.02e23 mlc/mole
    use PhysicalConstants_ml, only : AtwAir ! 28.964 (g air)/mole
    use GridValues_ml       , only : glon,glat,coord_in_domain
    use GridValues_ml       , only : A_bnd, B_bnd
    use MetFields_ml        , only : z_bnd    ! (nx,ny,1:nz+1) gph bounds relative to surface
    use MetFields_ml        , only : pzpbl    ! (nx,ny) boundary layer height
    use MetFields_ml        , only : roa
    use MetFields_ml        , only : ps
    use ChemFields_ml       , only : cfac
    use ChemFields_ml       , only : PM25_water
    use ChemFields_ml       , only : pm25_water_rh50

    ! --- in/out ----------------------------
    
    class(T_ObsOper), intent(inout)    ::  self
!    real, intent(in)              :: xn_adv_b(:,:,:,:)  ! (nspec_adv,lnx,lny,nlev)
!    character(len=*), intent(in)  :: xn_adv_units(:)    ! (nspec_adv)
    real, intent(in)              :: sf_obs(:,:  ,:)      ! (lnx,lny     ,nObsComp)
    real, intent(in)              :: xn_obs(:,:,:,:)      ! (lnx,lny,nlev,nObsComp)
    character(len=*), intent(in)  :: xn_obs_units(:)      ! (nObsComp)
    integer, intent(in)           :: ipar         ! index in 'obsData' array
    integer, intent(in)           :: stnid        ! station id (1-based record number in file)
    real, intent(inout)           :: lat,lon      ! location
    real, intent(in)              :: alt          ! altitude
    real, intent(out)             :: yn           ! simulation
    integer, intent(out)          :: status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOper_Fill'

    ! --- local -----------------------------

    logical              ::  local_model_grid
    integer              ::  i, j
    integer              ::  k
    integer              ::  ilev
    integer              ::  itr
    integer              ::  ispec
    real                 ::  xn, sf
    real                 ::  uconv
    logical              ::  needroa
    real, allocatable    ::  ph(:)  ! (0:nlev)
    real, allocatable    ::  H0(:)  ! (1:nlev)
    real, allocatable    ::  vcd(:)  ! (1:nlev)
    real                 ::  Hi
    real                 ::  xn1
    real                 ::  yn_blh

    ! --- begin -----------------------------
    
    ! store meta data:
    self%stnid = stnid
    self%lon   = lon
    self%lat   = lat
    self%alt   = alt
    ! default values, will be filled outside this routine:
    self%stncode = ''
    self%analyse = .true.
    ! store index in 'obsData' (and 'ObsKernels') array:
    self%iObsData = ipar

    ! get local grid indices
    local_model_grid = coord_in_domain("processor",lon,lat,iloc=i,jloc=j)

    ! check, expecting currently local indices only!
    if(.not.local_model_grid) then
      write (gol,'("found interpolation indices outside local domain:")'); call goErr
      write (gol,'("  local grid shape : ",2i4)') size(sf_obs,1), size(sf_obs,2); call goErr
      write (gol,'("  i index : ",i4)') i; call goErr
      write (gol,'("  j index : ",i4)') j; call goErr
      TRACEBACK; status=1; return
    end if

    !-----------------------------------------------------------------------
    ! setup obsData element
    !-----------------------------------------------------------------------

    ! store cell indices:
    self%i = i
    self%j = j

    ! copy value from obsData element into data type:
    self%iObsComp = obsData(self%iObsData)%iObsComp
    
    ! switch:
    select case ( trim(obsData(self%iObsData)%deriv) )
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! derive from 3D model levels (or 3D lowest layer) via surface (cfac)
      case ( '3D-ML-SFC', '2D-ML-SFC' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! cfac is for single species ...
        if ( ObsCompInfo(self%iObsComp)%nspec /= 1 ) then
          write (gol,'(a,": surface observation not implemented for accumulated tracer (",a,") yet ..")') &
                   rname, trim(ObsCompInfo(self%iObsComp)%name); call goErr
          TRACEBACK; status=1; return
        end if
        ! exctract from first (and only) index:
        itr = 1
        ! spec index in xn_adv array:
        ispec = ObsCompInfo(self%iObsComp)%ispec(itr)

        ! defined on model levels:
        select case ( trim(obsData(self%iObsData)%deriv) )
          case ( '3D-ML-SFC' )
            self%levtype = LEVTYPE_3D_ML_SFC
          case ( '2D-ML-SFC' )
            self%levtype = LEVTYPE_2D_ML_SFC
          case default
            write (gol,'("deriv `",a,"` not supported")') trim(obsData(self%iObsData)%deriv); call goErr
            TRACEBACK; status=1; return
        end select  ! deriv
        
        ! selected level:
        ilev = nlev    ! bottom layer (top-down order!)
        ! range of single layer only:
        self%l = (/ilev,ilev/)

        ! init with unit conversion:
        self%H_jac(ilev) = ObsCompInfo(self%iObsComp)%unitconv(itr)
        ! use density conversion?
        if ( ObsCompInfo(self%iObsComp)%unitroa(itr) ) then
          self%H_jac(ilev) = self%H_jac(ilev) * roa(i,j,ilev,1)
        end if
        ! extra factor for 50m->3m conversion:
        self%H_jac(ilev) = self%H_jac(ilev) * cfac(ispec,i,j)   ! 50m->3m conversion
        
        ! not used ...
        self%H_jac1 = -999.9

        ! check ...
        if ( trim(xn_obs_units(self%iObsComp)) /= 'mol mol^-1' ) then
          write (gol,'("xn_obs_units is `",a,"` while expected `",a,"`")') &
                  trim(xn_obs_units(self%iObsComp)), 'mol mol^-1'; call goErr
          TRACEBACK; status=1; return
        end if

        ! extract concentration:
        xn = xn_obs(i,j,ilev,self%iObsComp)
        ! simulate observation:
        yn = self%H_jac(ilev) * xn
    
!      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      ! derive from 2D surface field
!      case ( '2D-ML-SFC' )
!      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!        ! cfac is for single species ...
!        if ( ObsCompInfo(self%iObsComp)%nspec /= 1 ) then
!          write (gol,'(a,": surface observation not implemented for accumulated tracer (",a,") yet ..")') &
!                   rname, trim(ObsCompInfo(self%iObsComp)%name); call goErr
!          TRACEBACK; status=1; return
!        end if
!        ! exctract from first (and only) index:
!        itr = 1
!        ! spec index in xn_adv array:
!        ispec = ObsCompInfo(self%iObsComp)%ispec(itr)
!
!        ! defined on 2D field:
!        self%levtype = LEVTYPE_2D_ML_SFC
!        ! no level range involved:
!        self%l = (/-999,-999/)
!        ! selected level:
!        ilev = nlev    ! bottom layer (top-down order!)
!
!        ! store in first layer of H_jac:
!        k = 1
!        ! init with unit conversion:
!        self%H_jac(k) = ObsCompInfo(self%iObsComp)%unitconv(itr)
!        ! use density conversion?
!        if ( ObsCompInfo(self%iObsComp)%unitroa(itr) ) then
!          self%H_jac(k) = self%H_jac(k) * roa(i,j,ilev,1)
!        end if
!        ! extra factor for 50m->3m conversion:
!        self%H_jac(k) = self%H_jac(k) * cfac(ispec,i,j)   ! 50m->3m conversion
!
!        ! check ...
!        if ( trim(xn_obs_units(self%iObsComp)) /= 'mol mol^-1' ) then
!          write (gol,'("xn_obs_units is `",a,"` while expected `",a,"`")') &
!                  trim(xn_obs_units(self%iObsComp)), 'mol mol^-1'; call goErr
!          TRACEBACK; status=1; return
!        end if
!
!        ! extract concentration:
!        xn = sf_obs(i,j,self%iObsComp)   <-- should be xn_obs? cfac already applied ...
!        ! simulate observation:
!        yn = self%H_jac(k) * xn
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! derive from 3D model levels as total column
      case ( '3D-ML-TC' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! only for single species yet ..
        if ( ObsCompInfo(self%iObsComp)%nspec /= 1 ) then
          write (gol,'(a,": total column observation not implemented for accumulated tracer (",a,") yet ..")') &
                   rname, trim(ObsCompInfo(self%iObsComp)%name); call goErr
          TRACEBACK; status=1; return
        end if
        ! exctract from first (and only) index:
        itr = 1
        ! spec index in xn_adv array:
        ispec = ObsCompInfo(self%iObsComp)%ispec(itr)

        ! defined on 3D field:
        self%levtype = LEVTYPE_3D_ML_TC
        ! all levels are involved:
        self%l = (/1,nlev/)

        ! check ...
        if ( trim(xn_obs_units(self%iObsComp)) /= 'mol mol^-1' ) then
          write (gol,'("xn_obs_units is `",a,"` while expected `",a,"`")') &
                  trim(xn_obs_units(self%iObsComp)), 'mol mol^-1'; call goErr
          TRACEBACK; status=1; return
        end if

        ! conversion factor to <units>/m/(kg/m3),
        ! integral over layer height is included per level:
        select case ( trim(obsData(self%iObsData)%unit) )
          !~
          case ( '1e15molec/cm2' )
            ! [1e15 mlc/cm2]/m = 
            !     (mol tracer)/(mole air)   ! input units
            !     * (kg air)/m3             ! applied later
            !     * 1e-15   mlc/(mole tracer)  (mole air)/(kg air)   m2/cm2  
            uconv = 1e-15 *      Avog        /    (AtwAir*1e-3)    *  1e-4  ! <units>/m/(kg/m3)
            ! multiply with air density later:
            needroa = .true.
          !~
          case default
            write (gol,'("unsupported obsdata units `",a,"`")') &
                             trim(obsData(self%iObsData)%unit); call goErr
            TRACEBACK; status=1; return
        end select
        
        ! loop over layers:
        do ilev = self%l(0), self%l(1)
          ! init with factor:
          self%H_jac(ilev) = uconv   ! <units>/m/(kg/m3)
          ! use density conversion?
          if ( needroa ) then
            !                 <units>/m/(kg/m3)      (kg/m3)
            self%H_jac(ilev) = self%H_jac(ilev) * roa(i,j,ilev,1)   ! <units>/m
          end if
          ! layer thickness:
          !                      <units>/m    *                   m
          self%H_jac(ilev) = self%H_jac(ilev) * ( z_bnd(i,j,ilev) - z_bnd(i,j,ilev+1) )  ! <units>
        end do ! ilev
        
        ! not used ...
        self%H_jac1 = -999.9

        ! init simulation:
        yn = 0.0
        ! loop:
        do ilev = self%l(0), self%l(1)
          ! extract concentration:
          xn = xn_obs(i,j,ilev,self%iObsComp)
          ! add contribution to simulated observation:
          yn = yn + self%H_jac(ilev) * xn
        end do  ! ilev
    
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! derive from 3D model levels using kernels
      !case ( '3D-ML-K1', '3D-ML-K2' )
      case ( '3D-ML-K2' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! only for single species yet ..
        if ( ObsCompInfo(self%iObsComp)%nspec /= 1 ) then
          write (gol,'(a,": total column observation not implemented for accumulated tracer (",a,") yet ..")') &
                   rname, trim(ObsCompInfo(self%iObsComp)%name); call goErr
          TRACEBACK; status=1; return
        end if
        ! exctract from first (and only) index:
        itr = 1
        ! spec index in xn_adv array:
        ispec = ObsCompInfo(self%iObsComp)%ispec(itr)

        ! all levels are involved:
        self%l = (/1,nlev/)

        ! check ...
        if ( trim(xn_obs_units(self%iObsComp)) /= 'mol mol^-1' ) then
          write (gol,'("xn_obs_units is `",a,"` while expected `",a,"`")') &
                  trim(xn_obs_units(self%iObsComp)), 'mol mol^-1'; call goErr
          TRACEBACK; status=1; return
        end if

        ! conversion factor to <units>/m/(kg/m3),
        ! integral over layer height is included per level:
        select case ( trim(obsData(self%iObsData)%unit) )
          !~
          case ( '1e15molec/cm2' )
            ! [1e15 mlc/cm2]/m = 
            !     (mole tracer)/(mole air)   ! input units (ppv)
            !     * (kg air)/m3              ! applied later
            !     * 1e-15   mlc/(mole tracer)  (mole air)/(kg air)   m2/cm2  
            uconv = 1e-15 *      Avog        /    (AtwAir*1e-3)    *  1e-4  ! <units>/m/(kg/m3)/((mole tracer)/(mole air))
            ! multiply with air density later:
            needroa = .true.
          !~
          case default
            write (gol,'("unsupported obsdata units `",a,"`")') &
                             trim(obsData(self%iObsData)%unit); call goErr
            TRACEBACK; status=1; return
        end select
        
        !! testing ..
        !write (gol,*) 'xxx ps    ', lbound(ps), ';', ubound(ps); call goPr
        !write (gol,*) 'xxx A_bnd ', lbound(A_bnd), ';', ubound(A_bnd); call goPr
        !write (gol,*) 'xxx B_bnd ', lbound(B_bnd), ';', ubound(B_bnd); call goPr
        !write (gol,*) 'xxx cell  ', self%i, self%j; call goPr
        
        ! Limit change of model profile to OMI trop.column top at 200 hPa?
        ! First 17 layers are below this top ..
        
        ! storage for model half level pressures:
        allocate( ph(0:nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill using hybride coeff;
        ! first record of surface pressure is said to be the correct one ..
        ph = A_bnd + B_bnd * ps(self%i,self%j,1)
        
        ! check ...
        if ( trim(ObsKernels(self%iObsData)%phlev_units) /= 'Pa' ) then
          write (gol,'("unsupported retrieval layer units `",a,"`")') trim(ObsKernels(self%iObsData)%phlev_units); call goErr
          TRACEBACK; status=1; return
        end if
        
        !! testing ...
        !write (gol,'("ppp0 surface pressures: ",2f12.2)') ps(self%i,self%j,:); call goPr
        !write (gol,'("ppp0 model layers: ",i0)') nlev; call goPr
        !do ilev = 0, nlev
        !  write (gol,'("ppp0 model half level ",i2,3f12.2)') ilev, A_bnd(ilev+1), B_bnd(ilev+1), ph(ilev); call goPr
        !end do
        !write (gol,'("ppp0 ")'); call goPr
        
        ! jacobian without kernel:
        allocate( H0(1:nlev), stat=status, source=0.0 )
        IF_NOT_OK_RETURN(status=1)
        ! partial columns:
        allocate( vcd(1:nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)
        
        ! loop over layers:
        do ilev = self%l(0), self%l(1)
          ! init with factor:
          H0(ilev) = uconv   ! <units> / m / (kg/m3) / ppv
          ! use density conversion?
          if ( needroa ) then
            !          <units>/m/(kg/m3)/ppv      (kg/m3)
            H0(ilev) =        H0(ilev)        * roa(i,j,ilev,1)   ! <units> / m / ppv
          end if
          ! multiply with layer thickness:
          !          <units>/m/ppv   *                   m
          H0(ilev) =    H0(ilev)     * ( z_bnd(i,j,ilev) - z_bnd(i,j,ilev+1) )  ! <units> / ppv
          
          ! contribution of this layer in kernel-weighted column:
          call PartialColumnContribution( ph(ilev-1), ph(ilev), ObsKernels(self%iObsData)%nlayer, &
                                           ObsKernels(self%iObsData)%phlev (:,self%stnid), &
                                           ObsKernels(self%iObsData)%kernel(:,self%stnid), &
                                           Hi, status )
          IF_NOT_OK_RETURN(status=1)
          ! apply:
          !                   <units>/m/ppv  * 1
          self%H_jac(ilev) =      H0(ilev)   * Hi  ! <units> / ppv
          
          ! reset above 200 hPa:
          if ( ph(ilev) < 200.0e2 ) H0(ilev) = 0.0
        
          !! testing ..
          !if ( modulo(self%stnid-1,10) == 0 ) then
          !  write (gol,'("ppp0 layer weight ",i2,f8.2)') ilev, Hi; call goPr
          !end if

        end do ! ilev
        
        !! testing ...
        !write (gol,'("ppp1 model layers: ",i0)') nlev; call goPr
        !do ilev = 1, nlev
        !  write (gol,'("ppp1 layer ",i2,2f12.2,f16.2)') ilev, ph(ilev-1), ph(ilev), H0(ilev); call goPr
        !end do
        !write (gol,'("ppp1 ")'); call goPr
        
        !write (gol,*) 'xxx break'; call goPr
        !TRACEBACK; status=1; return

        ! init simulation:
        yn = 0.0
        yn_blh = 0.0
        ! loop:
        do ilev = self%l(0), self%l(1)
          ! extract concentration:
          xn = xn_obs(i,j,ilev,self%iObsComp)  ! ppv
          ! add contribution to simulated observation:
          yn = yn + self%H_jac(ilev) * xn      ! units
          ! (partly) in boundary layer?
          if ( z_bnd(i,j,ilev+1) < pzpbl(i,j) ) then
            yn_blh = yn_blh + self%H_jac(ilev) * xn      ! units
          end if
        end do  ! ilev
        
        !! testing ..
        !write (gol,*) 'www1f cell = ', i, j, self%l(0), self%l(1); call goPr
        !write (gol,*) '  w1f xn   = ', xn_obs(i,j,self%l(0):self%l(1),self%iObsComp); call goPr
        !write (gol,*) '  w1f Hj   = ', self%H_jac(self%l(0):self%l(1)); call goPr
        !write (gol,*) '  w1f yn   = ', yn, 'xf'; call goPr

        !! check ...
        !if ( yn <= 1.0e-2 ) then
        !  write (gol,'("unexpected small simulated OMI column ",e12.4)') yn; call goErr
        !  do ilev = self%l(0), self%l(1)
        !    write (gol,'("  layer ",i2," xn ",e12.4," H_jac ",e12.4)') &
        !             ilev, xn_obs(i,j,ilev,self%iObsComp), self%H_jac(ilev); call goErr
        !  end do
        !  TRACEBACK; status=1; return
        !end if
        
!        if ( trim(obsData(self%iObsData)%deriv) == '3D-ML-K1' ) then
!
!          ! analysis should scale concentration in boundary layer
!          ! relative to change at lowest model layer
!          ! (which is also used to simulate surface (at 3m) concentrations)
!          !     vcd       [units]
!          !     --- dx    -------  [ppv]
!          !     xn1        [ppv]
!          ! where both xn1 and dx are in model concentration units
!          ! (observation operator H_jac converts to observations units, 
!          ! for example ug/m3 or 1e15 mlc/cm2)   
!          ! level from which surface is computed:
!          ilev = nlev    ! bottom layer (top-down order!)
!          ! conventration value for lowest layer (in model units):
!          xn1 = xn_obs(i,j,ilev,self%iObsComp)  ! ppv
!          !
!          !! the "ObsCompInfo(self%iObsComp)%units" is the units of B (probably "ug/m3")
!          !! and therefore the units of the analysis increment ;
!          !! the "unitconv(itr)" factor and "unitroa(itr)" flag convert from model units (ppb?)
!          !! to these "increment" units:
!          !xn1 = xn1 * ObsCompInfo(self%iObsComp)%unitconv(itr)
!          !if ( ObsCompInfo(self%iObsComp)%unitroa(itr) ) then
!          !  xn1 = xn1 * roa(i,j,ilev,1)
!          !end if
!
!          ! defined on 3D field:
!          self%levtype = LEVTYPE_3D_ML_K1
!          ! ratio between boundary layer column and surface value;
!          ! skip values with very small concentrations in lowest layer (in ppv):
!          if ( xn1 < 1.0e-12 ) then
!            self%H_jac1 = 0.0
!          else
!            !! TESTING: blame boundary layer only:
!            !self%H_jac1 = yn_blh / xn1   ! units/ppv
!            ! blame entire profile:
!            self%H_jac1 = yn / xn1   ! units/ppv
!          end if
!        
!          !! testing ...
!          !write (gol,*) 'xxx1 H_jac1 : ', yn_blh, xn1, self%H_jac1; call goPr
!          
!        else 
        if ( trim(obsData(self%iObsData)%deriv) == '3D-ML-K2' ) then

          ! defined on 3D field:
          self%levtype = LEVTYPE_3D_ML_K2
          
          ! assume that increment dx is distrbuted relative to current vcd:
          !   x(:) +  vcd(:)/sum(vcd(:)) dx
          ! jacobian of observation operator w.r.t. dx:
          !   H_jac^T vcd(:)/sum(vcd)
          
          ! current columns:
          !     units/ppv             ppv
          vcd =     H0    * xn_obs(i,j,:,self%iObsComp)  ! units
          
          ! jacobian for application of dx:
          if ( sum(vcd) > 0.0 ) then
            !                   units/ppv   units/units
            self%H_jac1 = sum( self%H_jac * vcd )/sum(vcd)   ! units/ppv
          else
            ! no change:
            self%H_jac1 = 0.0   ! units/ppv
          end if
        
          !! testing ...
          !write (gol,*) 'xxx1 H_jac1 : ', i,j, sum(vcd), self%H_jac1; call goPr
          !TRACEBACK; status=1; return
          
        else
          write (gol,'("unsupported deriv `",a,"`")') trim(obsData(self%iObsData)%deriv); call goErr
          TRACEBACK; status=1; return
        end if

        ! clear:
        deallocate( ph, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( H0, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( vcd, stat=status )
        IF_NOT_OK_RETURN(status=1)
                         
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! 2D field of observation simulations
      case ( '2D-OBS-SFC' )
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! defined on 2D field:
        self%levtype = LEVTYPE_2D_OBS_SFC
        ! no level range involved:
        self%l = (/-999,-999/)

        ! check ...
        if ( trim(obsData(self%iObsData)%unit) /= trim(ObsCompInfo(self%iObsComp)%units) ) then
          write (gol,'("obsdata units are `",a,"` while observed component units are `",a,"`")') &
                  trim(obsData(self%iObsData)%unit), trim(ObsCompInfo(self%iObsComp)%units); call goErr
          TRACEBACK; status=1; return
        end if   
        ! check ...
        if ( trim(xn_obs_units(self%iObsComp)) /= trim(ObsCompInfo(self%iObsComp)%units) ) then
          write (gol,'("xn_obs_units is `",a,"` while expected `",a,"`")') &
                  trim(xn_obs_units(self%iObsComp)), trim(ObsCompInfo(self%iObsComp)%units); call goErr
          TRACEBACK; status=1; return
        end if
        
        ! level from which surface is computed:
        ilev = nlev    ! bottom layer (top-down order!)
        ! extract simulation first model level:
        xn = xn_obs(i,j,ilev,self%iObsComp)
        ! extract simulation from surface field:
        sf = sf_obs(i,j,self%iObsComp)

!        ! store weight in first layer of H_jac:
!        k = 1
!        ! set only first value of profile,
!        ! here use ratio between surface and model layer:
!        self%H_jac = 0.0
!        self%H_jac(k) = sf / xn
!
!        ! combine:
!        yn = self%H_jac(1) * xn

        ! store weight in first layer of H_jac:
        k = 1
        ! set only first value of profile only;
        ! not that it will extract from the surface field:
        self%H_jac    = 0.0
        self%H_jac(k) = 1.0
        
        ! not used ...
        self%H_jac1 = -999.9

        ! combine:
        yn = self%H_jac(k) * sf
    

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      case default
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        write (gol,'("deriv `",a,"` not supported")') trim(obsData(self%iObsData)%deriv); call goErr
        TRACEBACK; status=1; return

    end select  ! deriv

    ! ok
    status = 0

  end subroutine ObsOper_Fill
  
  
  ! ***
  
  ! 
  ! Contribution of pressure interval [p1,p2] to kernel weighted column.
  !
  ! Total column is:
  !        nr
  !    y = sum K(r) vcd(r)
  !        r=1
  ! with:
  !    nr         : number of retrieval layers
  !    K(1:nr)    : kernel
  !    vcd(1:nr)  : partial columns, for example in units "1e15 mlc/cm2"
  !                 valid for retrieval layers defined by boundaries phlev(0:nr)
  !
  ! The interval [p1,p2] is distributed over the retrieval layers:
  !
  !               4
  !    phlev(3) --------
  !               3   ------ p2
  !    phlev(2) --------  |      25% to retrieval layer 3
  !               2       |      50% to retrieval layer 2
  !    phlev(1) --------  |      25% to retrieval layer 1
  !               1   ------ p1
  !    phlev(0) --------
  !
  ! The total column is then computed from:
  !
  !     y = K(1) * ( 0.25 * vcd([p1,p2]) + ... fractions from other vcd's ... )
  !       + K(2) * ( 0.50 * vcd([p1,p2]) )
  !       + K(3) * ( 0.25 * vcd([p1,p2]) + ... fractions from other vcd's ... )
  !       + K(4) ...
  !
  ! The contribution of [p1,p2] to the total column is therefore:
  !
  !     y =  K(1) * 0.25  +  K(2) * 0.50  +  K(3) * 0.25
  
  subroutine PartialColumnContribution( p1, p2, nr, phlev, kernel, H, status )

    ! --- in/out -----------------------------
    
    real, intent(in)                  ::  p1, p2
    integer, intent(in)               ::  nr
    real, intent(in)                  ::  phlev(0:nr)
    real, intent(in)                  ::  kernel(1:nr)
    real, intent(out)                 ::  H
    integer, intent(out)              ::  status
    
    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/PartialColumnContribution'

    ! --- local -----------------------------
    
    integer     ::  ir
    real        ::  pmin, pmax
    real        ::  rmin, rmax
    real        ::  op1, op2
    real        ::  fraction
    
    ! --- begin -----------------------------
    
    !! testing ...
    !write (gol,'("ppp1 model interval: [",f9.2,",",f9.2,"]")') p1, p2; call goPr
    !write (gol,'("ppp1 retrieval layers: ",i0)') nr; call goPr
    
    ! pressure range:
    pmin = min( p1, p2 )
    pmax = max( p1, p2 )
    
    ! init output:
    H = 0
    ! loop over retrieval layers:
    do ir = 1, nr
      ! retrieval layer pressure range:
      rmin = min( phlev(ir-1),phlev(ir) )
      rmax = max( phlev(ir-1),phlev(ir) )
      ! no overlap?
      if ( (pmax <= rmin) .or. (pmin >= rmax) ) then
        ! no coverage:
        fraction = 0.0
      else
        ! overlapping interval:
        op1 = max( pmin, rmin )
        op2 = min( pmax, rmax )
        ! fraction of [p1,p2] that overlaps with retrieval layer:
        fraction = (op2-op1)/(pmax-pmin)
        ! add weighted 
        H = H + kernel(ir) * fraction
      end if
      !! testing:
      !write (gol,'("ppp1 retrieval layer ",i2," [",f9.2,",",f9.2,"] ",f12.6,f8.2)') &
      !                ir, phlev(ir-1), phlev(ir), kernel(ir), fraction; call goPr
    end do
    !! testing ...
    !write (gol,'("ppp1 ")'); call goPr
  
    ! ok
    status = 0
    
  end subroutine PartialColumnContribution


  ! ***


  subroutine ObsOper_Evaluate( self, key, sf_b, xn_b, status, debug )

    ! --- in/out -----------------------------
    
    class(T_ObsOper), intent(inout)   ::  self
    character(len=*), intent(in)      ::  key
    real, intent(in)                  ::  sf_b(:,:  ,:)      ! (lnx,lny     ,nObsComp)  2D field (sfc obs?)
    real, intent(in)                  ::  xn_b(:,:,:,:)      ! (lnx,lny,nlev,nObsComp)  3D field
    integer, intent(out)              ::  status
    
    logical, intent(in), optional     ::  debug

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOper_Evaluate'

    ! --- local -----------------------------
    
    !logical         ::  dbg
    
    integer         ::  i, j
    integer         ::  l0, l1
    integer         ::  iObsComp
    real            ::  yn
    
    ! --- begin -----------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug

    ! involved grid cell:
    i = self%i
    j = self%j

    ! observed tracer index:
    iObsComp = self%iObsComp

    ! switch:
    select case ( self%levtype )
      !~ surface:
      case ( LEVTYPE_2D_OBS_SFC )
        ! extract simulation for this cell:
        yn = self%H_jac(1) * sf_b(i,j,iObsComp)
        !! testing ..
        !if (dbg) then
        !  write (gol,*) 'yyy1 station ', trim(self%stncode); call goErr
        !  write (gol,*) '  y1 comp ', self%iObsComp, trim(ObsCompInfo(self%iObsComp)%name); call goErr
        !  write (gol,*) '  y1 cell ', i, j; call goErr
        !  write (gol,*) '  y1 sf_b ', self%H_jac(1), sf_b(i,j,iObsComp); call goErr
        !  write (gol,*) '  y1 yn ', yn; call goErr
        !  TRACEBACK; status=1; return
        !end if
      !~ model levels:
      case ( LEVTYPE_2D_ML_SFC, LEVTYPE_3D_ML_SFC, LEVTYPE_3D_ML_TC, LEVTYPE_3D_ML_K2 )
        ! level range:
        l0 = self%l(0)
        l1 = self%l(1)
        ! extract simulation for this cell, sum over layers:
        yn = sum( self%H_jac(l0:l1) * xn_b(i,j,l0:l1,iObsComp) )
        !write (gol,*) 'www1w cell = ', i, j, l0, l1; call goPr
        !write (gol,*) '  w1w xn   = ', xn_b(i,j,l0:l1,iObsComp); call goPr
        !write (gol,*) '  w1w Hj   = ', self%H_jac(l0:l1); call goPr
        !write (gol,*) '  w1w yn   = ', yn, key; call goPr
      !~ unknown
      case default
        write (gol,'("unsupported levtype ",i0)') self%levtype; call goErr
        TRACEBACK; status=1; return
    end select

    ! store:
    select case ( key )
      case ( 'xf' )
        self%xf = yn
      case ( 'xa' )
        self%xa = yn
      case ( 'sb' )
        self%sb = yn
      case ( 'sf' )
        self%sf = yn
      case default
        write (gol,'("unsupported key `",a,"`")') trim(key); call goErr
        TRACEBACK; status=1; return
    end select

    ! ok
    status = 0
    
  end subroutine ObsOper_Evaluate


  !=========================================================================
  !===
  !=== observation operators
  !===
  !=========================================================================


  subroutine ObsOpers_Init( self, status )

    !-----------------------------------------------------------------------
    ! @description
    ! Allocate the arrays innov and obsstddev. This step is
    ! required, since the number of observations can vary from
    ! one time step to the next
    ! @author AMVB
    !-----------------------------------------------------------------------

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(out)  ::  self
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Init'

    ! --- begin -----------------------------
    
    ! no values yet:
    self%nobs = -999

    ! ok
    status = 0

  end subroutine ObsOpers_Init


  ! ***


  subroutine ObsOpers_Done( self, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(inout)  ::  self
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Done'

    ! --- begin -----------------------------
    
    ! clear storage for observations:
    call self%DeAlloc( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear storage of over-all scale factors:
    if ( allocated(self%BScale) ) then
      deallocate( self%BScale, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0

  end subroutine ObsOpers_Done


  ! ***


  subroutine ObsOpers_Alloc( self, nobs, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(out)  ::  self
    integer, intent(in)             ::  nobs
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Alloc'

    ! --- local -----------------------------
    
    integer   ::  iobs

    ! --- begin -----------------------------
    
    ! different from current ?
    if ( self%nobs /= nobs ) then
      ! clear (if necessary):
      call self%DeAlloc( status )
      IF_NOT_OK_RETURN(status=1)
      ! new; trap empty:
      if ( nobs == 0 ) then
        ! dummy with single element to avoid problems:
        allocate( self%obs(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! storage:
        allocate( self%obs(nobs), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! loop over new values:
      do iobs = 1, size(self%obs)
        ! init value:
        call ObsOper_Init( self%obs(iobs), status )
        IF_NOT_OK_RETURN(status=1)
      end do
    end if
    
    ! store:
    self%nobs = nobs
    
    ! ok
    status = 0

  end subroutine ObsOpers_Alloc


  ! ***


  subroutine ObsOpers_DeAlloc( self, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(inout)  ::  self
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_DeAlloc'

    ! --- local -----------------------------
    
    integer   ::  iobs

    ! --- begin -----------------------------
    
    ! initialized ?
    if ( allocated(self%obs) ) then
      ! loop over new values:
      do iobs = 1, self%nobs
        ! done with value:
        call ObsOper_Done( self%obs(iobs), status )
        IF_NOT_OK_RETURN(status=1)
      end do
      ! clear:
      deallocate( self%obs, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! empty:
      self%nobs = -999
    end if

    ! ok
    status = 0

  end subroutine ObsOpers_DeAlloc


  ! ***


  subroutine ObsOpers_Set_IDs( self, status )

  !-----------------------------------------------------------------------
  ! @description
  ! Set unique for each id equal to global order in all domains.
  !-----------------------------------------------------------------------
  
    use MPI   , only : MPI_INTEGER
    use MPIF90, only : MPIF90_Displacements
    
    use MPI_Groups_ml    , only : MPI_COMM_CALC
    use ModelConstants_ml, only : nproc
    use Par_ml           , only : me

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(inout)  ::  self
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Set_IDs'

    ! --- local -----------------------------
    
    integer                         ::  iobs
    integer, allocatable            ::  recvcounts(:)  ! (0:nproc-1)
    integer, allocatable            ::  rdispls(:)     ! (0:nproc-1)
    
    ! --- begin -----------------------------

    ! storage:
    allocate( recvcounts(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! exchange number of observations:
    call MPI_AllGather( self%nobs , 1, MPI_INTEGER, &
                        recvcounts, 1, MPI_INTEGER, &
                        MPI_COMM_CALC, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! create displacements from these:
    call MPIF90_Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    ! fill id's:
    do iobs = 1, self%nobs
      self%obs(iobs)%id = rdispls(me) + iobs
    end do
    
    ! info ...
    if ( self%nobs > 0 ) then
      write (gol,'(a,": proc observations : ",i0," - ",i0)') rname, &
                    self%obs(1)%id, self%obs(self%nobs)%id; call goPr
    else
      write (gol,'(a,": proc observations : none")') rname; call goPr
    end if
    
    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine ObsOpers_Set_IDs


  ! ***


  subroutine ObsOpers_Swap( self, doms, Hops_new, doms_new, status )

  !-----------------------------------------------------------------------
  ! @description
  ! Fill initialized obseration operator decomposed over domains 'doms_new'
  ! from similar structure 'self' decomposed over domains 'doms'.
  !-----------------------------------------------------------------------
  
    ! NOTE: This routine was originaly implemented by defining
    ! a new MPI type for an instance of the 'T_ObsOper' class.
    ! However, this requires that the compiler aligns the variables
    ! contigeously but that is dificult to ensure for all compilers around ...
    ! Therefore here we simply copy the data into arrays and swap
    ! these to the other domains.
  
    use MPIF90       , only : MPIF90_AllToAll
    use MPIF90       , only : MPIF90_AllToAllV
    use MPIF90       , only : MPIF90_Displacements
    
    use MPI_Groups_ml    , only : MPI_COMM_CALC
    use ModelConstants_ml, only : nproc
    use Par_ml           , only : me
    use DA_Util_ml       , only : T_Domains
    
    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(in)     ::  self
    type(T_Domains), intent(in)       ::  doms
    type(T_ObsOpers), intent(inout)   ::  Hops_new
    type(T_Domains), intent(in)       ::  doms_new
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Swap'

    ! --- local -----------------------------
    
    integer                         ::  iobs, iobs_in
    integer                         ::  iproc, iproc_in
    integer                         ::  ind(4)
    !integer                         ::  nrr, sHj(2), nii
    integer                         ::  nrr, nii
    real, allocatable               ::  srt_rr(:,:)
    !real, allocatable               ::  srt_Hj(:,:,:)  ! (nlev,nchemobs,nobs)
    real, allocatable               ::  srt_Hj(:,:)     ! (nlev,nobs)
    real, allocatable               ::  srt_ii(:,:)
    real, allocatable               ::  new_rr(:,:)
    !real, allocatable               ::  new_Hj(:,:,:)  ! (nlev,nchemobs,nobs)
    real, allocatable               ::  new_Hj(:,:)     ! (nlev,nobs)
    real, allocatable               ::  new_ii(:,:)
    integer, allocatable            ::  sendcounts(:)  ! (0:nproc-1)
    integer, allocatable            ::  sdispls(:)     ! (0:nproc-1)
    integer, allocatable            ::  recvcounts(:)  ! (0:nproc-1)
    integer, allocatable            ::  rdispls(:)     ! (0:nproc-1)
    integer                         ::  nobs_new
    integer                         ::  nval
    
    ! --- begin -----------------------------
    
    ! storage:
    allocate( sendcounts(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( sdispls(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( recvcounts(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for sorted version that can be passed through MPI;
    ! decomposition over observation, so this should be last index:
    nrr = 3
    allocate( srt_rr(nrr,self%nobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    nii = 8
    allocate( srt_ii(nii,self%nobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! Hj
    allocate( srt_Hj(nlev,self%nobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init counter for sorted obs on this proc:
    iobs = 0
    ! loop over target procs:
    do iproc = 0, nproc-1
      ! set offset, init counter:
      sdispls   (iproc) = iobs  ! zero based offset
      sendcounts(iproc) = 0
      ! loop over input observations:
      do iobs_in = 1, self%nobs
        ! global index of first grid cell involved:
        ind = [doms%off(1,me) + self%obs(iobs_in)%i, &
               doms%off(2,me) + self%obs(iobs_in)%j, &
               1, 1 ]
        ! target proc:
        call doms_new%Find( ind, iproc_in, status )
        IF_NOT_OK_RETURN(status=1)
        ! current ?
        if ( iproc_in == iproc ) then
          ! increase counters:
          iobs = iobs + 1
          sendcounts(iproc) = sendcounts(iproc) + 1
          ! copy values, fill with global horizontal indices:
          srt_rr(1,iobs)   = self%obs(iobs_in)%innov
          srt_rr(2,iobs)   = self%obs(iobs_in)%obsstddev
          srt_rr(3,iobs)   = self%obs(iobs_in)%H_jac1
          srt_Hj(:  ,iobs) = self%obs(iobs_in)%H_jac
          srt_ii(1  ,iobs) = doms%off(1,me) + self%obs(iobs_in)%i
          srt_ii(2  ,iobs) = doms%off(2,me) + self%obs(iobs_in)%j
          srt_ii(3  ,iobs) = self%obs(iobs_in)%levtype
          srt_ii(4:5,iobs) = self%obs(iobs_in)%l
          srt_ii(6  ,iobs) = self%obs(iobs_in)%id
          srt_ii(7  ,iobs) = self%obs(iobs_in)%iObsComp
          if ( self%obs(iobs_in)%analyse ) then
            srt_ii(8  ,iobs) = 1
          else
            srt_ii(8  ,iobs) = 0
          end if
        end if
      end do  ! input obs
    end do  ! target procs
    ! check ...
    if ( iobs /= self%nobs ) then
      write (gol,'("number of obs in sorted operator ",i0," differs from original ",i0)') iobs, self%nobs; call goErr
      TRACEBACK; status=1; return
    end if

    ! collect number of observations to be received:
    call MPIF90_AllToAll( sendcounts, 1, &
                          recvcounts, 1, &
                          MPI_COMM_CALC, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! compute displacements:
    call MPIF90_Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! total number to be received:
    nobs_new = sum(recvcounts)

    ! storage for arrays to be received; 
    ! allocate at least something to avoid bound check errors:
    nval = max( 1, nobs_new )
    allocate( new_rr(nrr,nval), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( new_Hj(nlev,nval), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( new_ii(nii,nval), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! swap real values:
    nval = nrr
    call MPIF90_AllToAllV( srt_rr, nval*sendcounts, nval*sdispls, &
                           new_rr, nval*recvcounts, nval*rdispls, &
                           MPI_COMM_CALC, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! swap Hj values:
    nval = nlev
    call MPIF90_AllToAllV( srt_Hj, nval*sendcounts, nval*sdispls, &
                           new_Hj, nval*recvcounts, nval*rdispls, &
                           MPI_COMM_CALC, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    ! swap integer values:
    nval = nii
    call MPIF90_AllToAllV( srt_ii, nval*sendcounts, nval*sdispls, &
                           new_ii, nval*recvcounts, nval*rdispls, &
                           MPI_COMM_CALC, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! storage for structures:
    call Hops_new%Alloc( nobs_new, status )
    IF_NOT_OK_RETURN(status=1)
    ! unpack:
    do iobs = 1, nobs_new
      ! copy values, convert from global to local horizontal indices:
      Hops_new%obs(iobs)%innov     = new_rr(1,iobs)
      Hops_new%obs(iobs)%obsstddev = new_rr(2,iobs)
      Hops_new%obs(iobs)%H_jac1    = new_rr(3,iobs)
      Hops_new%obs(iobs)%H_jac     = new_Hj(:  ,iobs)
      Hops_new%obs(iobs)%i         = new_ii(  1,iobs) - doms_new%off(1,me)
      Hops_new%obs(iobs)%j         = new_ii(  2,iobs) - doms_new%off(2,me)
      Hops_new%obs(iobs)%levtype   = new_ii(  3,iobs)
      Hops_new%obs(iobs)%l         = new_ii(4:5,iobs)
      Hops_new%obs(iobs)%id        = new_ii(  6,iobs)
      Hops_new%obs(iobs)%iObsComp  = new_ii(  7,iobs)
      Hops_new%obs(iobs)%analyse   = new_ii(  8,iobs) == 1
    end do
    
    ! clear:
    deallocate( new_rr, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( new_Hj, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( new_ii, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( srt_rr, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( srt_Hj, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( srt_ii, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( sdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ObsOpers_Swap


  ! ***


  !-----------------------------------------------------------------------
  ! @description
  ! Fill self (already initialized) with 
  ! subset of H with only observed components
  !-----------------------------------------------------------------------
  
  subroutine ObsOpers_SelectTracers( self, H, iObsComps, status, nana )

    use MPIF90           , only : MPIF90_AllReduce, MPI_SUM
    use MPI_Groups_ml    , only : MPI_COMM_CALC
    use ModelConstants_ml, only : MasterProc

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(inout)  ::  self
    type(T_ObsOpers), intent(in)      ::  H
    integer, intent(in)               ::  iObsComps(:)  ! (nselected)
    integer, intent(out)              ::  status
    integer, intent(out), optional    ::  nana

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_SelectTracers'

    ! --- local -----------------------------
    
    integer                         ::  nobs
    integer                         ::  iobs
    integer                         ::  k
    
    ! counters ...
    integer, allocatable    ::  tab_loc(:,:,:)  ! (nObsComp,NLEVTYPE,ana/val)
    integer, allocatable    ::  tab_tot(:,:,:)  ! (nObsComp,NLEVTYPE,ana/val)
    integer                 ::  iobscomp
    integer                 ::  ilevtype
    integer                 ::  iav
    
    ! --- begin -----------------------------
    
    ! init counter:
    if ( present(nana) ) nana = 0
    
    ! storage:
    allocate( tab_loc(nObsComp,NLEVTYPE,2), stat=status )    
    IF_NOT_OK_RETURN(status=1)
    allocate( tab_tot(nObsComp,NLEVTYPE,2), stat=status )    
    IF_NOT_OK_RETURN(status=1)
    ! init:
    tab_loc = 0

    ! init new number of observations:
    nobs = 0
    ! loop over existing observations:
    do k = 1, H%nobs
      ! match?
      if ( any( H%obs(k)%iObsComp == iObsComps ) ) then
        ! increase counter:
        nobs = nobs + 1
        ! increase specific counter?
        if ( present(nana) ) then
          if ( H%obs(k)%analyse ) nana = nana + 1
        end if ! nana
        ! increase stat:
        iobscomp = H%obs(k)%iObsComp
        ilevtype = H%obs(k)%levtype
        iav = 2
        if ( H%obs(k)%analyse ) iav = 1
        tab_loc(iobscomp,ilevtype,iav) = tab_loc(iobscomp,ilevtype,iav) + 1
      end if ! match comp
    end do ! obs
    
    ! storage for structures:
    call self%Alloc( nobs, status )
    IF_NOT_OK_RETURN(status=1)
    ! init new index:
    iobs = 0
    ! loop over existing observations:
    do k = 1, H%nobs
      ! match?
      if ( any( H%obs(k)%iObsComp == iObsComps ) ) then
        ! increase index:
        iobs = iobs + 1
        ! copy values:
        call self%obs(iobs)%Copy( H%obs(k), status )
        IF_NOT_OK_RETURN(status=1)
      end if ! selected?
    end do ! original obs
    ! check ...
    if ( iobs /= nobs ) then
      write (gol,'("something wrong, filled ",i0," of ",i0," selected observations")') iobs, nobs; call goErr
      TRACEBACK; status=1; return
    end if

    ! collect from other processes:
    call MPIF90_AllReduce( tab_loc, tab_tot, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do iobscomp = 1, nObsComp
      do ilevtype = 1, NLEVTYPE
        ! overview:
        if ( any( tab_tot(iobscomp,ilevtype,:) > 0 ) ) then
          write (gol,'(a,":     number of ",a6," ",a10," observations: ",i6," (",i6," analyzed, ",i6," validated)")') &
                              rname, trim(ObsCompInfo(iobscomp)%name), &
                              trim(LEVTYPE_NAME(ilevtype)), sum(tab_tot(iobscomp,ilevtype,:)), &
                              tab_tot(iobscomp,ilevtype,1), tab_tot(iobscomp,ilevtype,2); call goPr
        end if  ! any obs
        ! for logfile:
        if ( MasterProc ) then
          if ( tab_tot(iobscomp,ilevtype,1) > 0 ) then
            write (gol,'(a,": number of analyzed ",a6," ",a10," observations: ",i6)') &
                              rname, trim(ObsCompInfo(iobscomp)%name), &
                              trim(LEVTYPE_NAME(ilevtype)), tab_tot(iobscomp,ilevtype,1); call goPr
            !call PrintLog(gol)
          end if  ! any obs
        end if  ! master
      end do ! levtype
    end do ! comp

    ! clear:
    deallocate( tab_loc, stat=status )    
    IF_NOT_OK_RETURN(status=1)
    deallocate( tab_tot, stat=status )    
    IF_NOT_OK_RETURN(status=1)


    ! ok
    status = 0

  end subroutine ObsOpers_SelectTracers


  ! ***


  subroutine ObsOpers_Evaluate( self, key, sf_b, xn_b, status, debug )

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(inout)  ::  self
    character(len=*), intent(in)      ::  key
    real, intent(in)                  ::  sf_b(:,:  ,:)      ! (lnx,lny     ,nObsComp)
    real, intent(in)                  ::  xn_b(:,:,:,:)      ! (lnx,lny,nlev,nObsComp)
    integer, intent(out)              ::  status
    
    logical, intent(in), optional     ::  debug

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Evaluate'

    ! --- local -----------------------------
    
    !logical         ::  dbg    
    integer         ::  n
    
    ! --- begin -----------------------------
    
    !! debug?
    !dbg = .false.
    !if ( present(debug) ) dbg = debug

    ! any local observations?
    if ( self%nobs > 0 ) then
      ! loop over observations:
      do n = 1, self%nobs
        ! evaluate:
        call self%obs(n)%Evaluate( key, sf_b, xn_b, status )
                                    !debug=dbg .and. (trim(self%obs(n)%stncode)=='PT07001') &
                                    !      .and. (self%obs(n)%iObsComp==2) )
        IF_NOT_OK_RETURN(status=1)
      end do ! n
    end if  ! nobs > 0
    
    ! ok
    status = 0
    
  end subroutine ObsOpers_Evaluate


  ! ***


  subroutine ObsOpers_SetBScale( self, status )

    use MPIF90           , only : MPIF90_AllGather
    use MPIF90           , only : MPIF90_Displacements
    use MPIF90           , only : MPIF90_GatherV
    use MPIF90           , only : MPIF90_BCast
    use MPI_Groups_ml    , only : MPI_COMM_CALC
    use MPI_Groups_ml    , only : MasterPE
    use ModelConstants_ml, only : nproc
    use Par_ml           , only : me

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(inout)  ::  self
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_SetBScale'
    
    ! minimum number of chi values for stats:
    integer, parameter    ::  nchi_min = 10

    ! scale factors in [1,3]:
    real, parameter       ::  fscale0 = 1.0
    real, parameter       ::  dscale  = 0.1
    integer, parameter    ::  nscale  = 20+1

    ! --- local -----------------------------
    
    integer, allocatable      ::  recvcounts(:)  ! (0:nproc-1)
    integer, allocatable      ::  rdispls(:)     ! (0:nproc-1)

    integer                   ::  nobs_all
    integer, allocatable      ::  ivalues(:)
    real, allocatable         ::  rvalues(:)

    integer, allocatable      ::  iObsData_all(:)
    integer, allocatable      ::  analyse_all(:)
    real, allocatable         ::  yy_all(:)
    real, allocatable         ::  rr_all(:)
    real, allocatable         ::  xf_all(:)
    real, allocatable         ::  sb_all(:)

    integer                   ::  iObsComp
    integer                   ::  iobs
    integer, allocatable      ::  chi(:)
    integer                   ::  nchi
    real                      ::  chi_mean, chi_stdv

    integer                   ::  iscale
    real                      ::  fscale
    real                      ::  fscale_opt
    real                      ::  chi_stdv_opt
    character(len=1)          ::  mark

    ! --- begin -----------------------------
      
    ! info ..
    write (gol,'("Compute B sigma scaling from Chi2 test")'); call goPr
    
    ! setup storage if not done yet:
    if ( .not. allocated(self%BScale) ) then
      allocate( self%BScale(nObsComp), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! init result to no scaling ..
    self%Bscale = 1.0
    
    ! storage:
    allocate( recvcounts(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! exchange number of observations:
    call MPIF90_AllGather( (/self%nobs/), recvcounts, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! create displacements from these:
    call MPIF90_Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! total:
    nobs_all = sum(recvcounts)

    ! check ...
    if ( nobs_all == 0 ) then
      
      ! info ..
      write (gol,'("WARNING - no observations at all, default scaling ...")'); call goPr

    else

      ! reset for procs that do not need to receive:
      if ( .not. MasterProc )  nobs_all = 1  

      ! target arrays:
      allocate( iObsData_all(nobs_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( analyse_all(nobs_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( yy_all(nobs_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( rr_all(nobs_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( xf_all(nobs_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( sb_all(nobs_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! ~ gather on root
      
      ! storage for local data:
      allocate( ivalues(max(1,self%nobs)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( rvalues(max(1,self%nobs)), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! fill local values:
      do iobs = 1, self%nobs
        ivalues(iobs) = self%obs(iobs)%iObsData
      end do
      ! gather in 'values_all' on root:
      call MPIF90_GatherV( ivalues, self%nobs, &
                           iObsData_all, recvcounts, rdispls, &
                           MasterPE, MPI_COMM_CALC, status )
      IF_MPI_NOT_OK_RETURN(status=1)

      ! fill local values:
      do iobs = 1, self%nobs
        if ( self%obs(iobs)%analyse ) then
          ivalues(iobs) = 1
        else
          ivalues(iobs) = 0
        end if
      end do
      ! gather in 'values_all' on root:
      call MPIF90_GatherV( ivalues, self%nobs, &
                           analyse_all, recvcounts, rdispls, &
                           MasterPE, MPI_COMM_CALC, status )
      IF_MPI_NOT_OK_RETURN(status=1)

      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%obs
      end do
      ! gather in 'values_all' on root:
      call MPIF90_GatherV( rvalues, self%nobs, &
                           yy_all, recvcounts, rdispls, &
                           MasterPE, MPI_COMM_CALC, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%obsstddev
      end do
      ! gather in 'values_all' on root:
      call MPIF90_GatherV( rvalues, self%nobs, &
                           rr_all, recvcounts, rdispls, &
                           MasterPE, MPI_COMM_CALC, status )
      IF_MPI_NOT_OK_RETURN(status=1)

      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%xf
      end do
      ! gather in 'values_all' on root:
      call MPIF90_GatherV( rvalues, self%nobs, &
                           xf_all, recvcounts, rdispls, &
                           MasterPE, MPI_COMM_CALC, status )
      IF_MPI_NOT_OK_RETURN(status=1)

      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%sb
      end do
      ! gather in 'values_all' on root:
      call MPIF90_GatherV( rvalues, self%nobs, &
                           sb_all, recvcounts, rdispls, &
                           MasterPE, MPI_COMM_CALC, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( ivalues, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( rvalues, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! ~ chi2 test
      
      ! all data collected on root ...
      if ( MasterProc ) then
      
        ! info ...
        write (gol,'("determine sigma scale factors ...")'); call goPr
      
        ! storage:
        allocate( chi(nobs_all), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over observed components:
        do iObsComp = 1, nObsComp
      
          ! info ...
          write (gol,'("  observed component ",i0," `",a,"`")') &
                             iObsComp, trim(ObsCompInfo(iObsComp)%name); call goPr
                             
          ! loop over scale factors:
          do iscale = 1, nscale
          
            ! current value:
            fscale = fscale0 + (iscale-1)*dscale

            ! counter for collected chi values:
            nchi = 0
            ! loop over observations:
            do iobs = 1, nobs_all
              ! skip if not analysed:
              if ( analyse_all(iobs) /= 1 ) cycle
              ! this component?
              if ( obsdata(iObsData_all(iobs))%iObsComp == iObsComp ) then
                ! increase counter:
                nchi = nchi + 1
                ! compute and store chi value:
                chi(nchi) = ( yy_all(iobs) - xf_all(iobs) ) / sqrt( rr_all(iobs)**2 + (fscale*sb_all(iobs))**2 )
              end if ! matched component
            end do ! obs

            ! not enough observations for this component? then next:
            if ( nchi < nchi_min ) then
              ! info ...
              write (gol,'("    found only ",i0," observations, at least ",i0," needed; keep default scaling")') nchi, nchi_min; call goPr
              ! default:
              fscale_opt = 1.0
              ! next component:
              exit
            end if
            
            ! mean value:
            chi_mean = sum( chi(1:nchi) ) / nchi
            ! std.dev:
            chi_stdv = sqrt( sum( (chi(1:nchi) - chi_mean)**2 )/(nchi-1.0) )
            ! closer to 1.0?
            if ( iscale == 1 ) then
              ! store first:
              fscale_opt = fscale
              chi_stdv_opt = chi_stdv
              ! current optimum:
              mark = ' '
            else
              ! closer to 1 ?
              if ( abs(chi_stdv-1.0) < abs(chi_stdv_opt-1.0) ) then
                ! reset optimal values:
                fscale_opt = fscale
                chi_stdv_opt = chi_stdv
                ! current optimum:
                mark = '*'
              else
                ! keep previous:
                mark = ' '
              end if
            end if
            ! info ...
            write (gol,'("    chi stdv : ",f6.1,f12.4," ",a)') fscale, chi_stdv, mark; call goPr
            
          end do
          
          ! store result:
          self%Bscale(iObsComp) = fscale_opt

        end do ! observed component
            
      end if ! root
    
      ! clear:
      deallocate( iObsData_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( analyse_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( yy_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( rr_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xf_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( sb_all, stat=status )
      IF_NOT_OK_RETURN(status=1)

    end if  ! any obs
    
    ! send scale factors from root to other processes:
    call MPIF90_BCast( self%Bscale, MasterPE, MPI_COMM_CALC, status )
    IF_MPI_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine ObsOpers_SetBScale


  ! ***


  subroutine ObsOpers_WriteToFile( self, cdate, status )

    use MPIF90           , only : MPIF90_AllGather
    use MPIF90           , only : MPIF90_Displacements
    use MPIF90           , only : MPIF90_GatherV
    use MPI_Groups_ml    , only : MPI_COMM_CALC
    use MPI_Groups_ml    , only : MasterPE
    use ModelConstants_ml, only : nproc
    use Par_ml           , only : me
    use TimeDate_ml      , only : date
    use C3PO             , only : Datafile
    use C3PO             , only : Dimension
    use C3PO             , only : IntegerCoordinate
    use NetCDF           , only : NF90_FLOAT, NF90_INT, NF90_CHAR
    use NetCDF           , only : NF90_Def_Var, NF90_Put_Var
    use NetCDF           , only : NF90_Put_Att

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(in)     ::  self
    type(date), intent(in)            ::  cdate     ! current date
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_WriteToFile'

    ! --- local -----------------------------
    
    integer, allocatable      ::  recvcounts(:)  ! (0:nproc-1)
    integer, allocatable      ::  rdispls(:)     ! (0:nproc-1)
    
    character(len=1024)       ::  fname
    type(Datafile)            ::  F
    integer                   ::  varid
    integer                   ::  varid_name, varid_units, varid_obstype
    integer                   ::  varid_iobsdata
    integer                   ::  varid_stnid
    integer                   ::  varid_scode
    integer                   ::  varid_analyse
    integer                   ::  varid_lon, varid_lat, varid_alt
    integer                   ::  varid_y, varid_r
    integer                   ::  varid_xf, varid_xa, varid_sb, varid_sf

    integer                   ::  iObsComp
    integer                   ::  varid_sigma_scale

    integer, allocatable      ::  ivalues(:), ivalues_all(:)
    real, allocatable         ::  rvalues(:)
    character(len=LEN_SCODE), allocatable         ::  svalues(:)
    integer                   ::  nobs_all
    integer                   ::  iobs
    type(IntegerCoordinate)   ::  obs_coor
    
    integer                   ::  iobsdata
    type(Dimension)           ::  obsdata_dim
    type(Dimension)           ::  obsdata_len_dim
    type(Dimension)           ::  scode_len_dim
    character(len=LEN_LABEL)  ::  label
    character(len=LEN_SCODE)  ::  scode
    integer                   ::  vlen
    
    ! --- begin -----------------------------
    
    ! storage:
    allocate( recvcounts(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( rdispls(0:nproc-1), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! exchange number of observations:
    call MPIF90_AllGather( (/self%nobs/), recvcounts, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! create displacements from these:
    call MPIF90_Displacements( recvcounts, rdispls, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! total:
    nobs_all = sum(recvcounts)

    ! check ...
    if ( nobs_all == 0 ) then
      
      ! info ..
      write (gol,'("WARNING - could not write, no observations at all ...")'); call goPr

    else

      ! storage for local data:
      allocate( ivalues(max(1,self%nobs)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( rvalues(max(1,self%nobs)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( svalues(max(1,self%nobs)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! write on root:
      if (MasterProc) then

        ! new coordinate:
        call obs_coor%Init( 'obs', status )
        IF_NOT_OK_RETURN(status=1)
        ! set total size:
        call obs_coor%Set_Dim( status, n=nobs_all )
        IF_NOT_OK_RETURN(status=1)
        ! set attributes:
        call obs_coor%Set_Attrs( status, units='1', long_name='observation id' )
        IF_NOT_OK_RETURN(status=1)

        ! new dim:
        call obsdata_dim%Init( 'obsdata', status )
        IF_NOT_OK_RETURN(status=1)
        ! set total size:
        call obsdata_dim%Set_Dim( status, n=nObsData )
        IF_NOT_OK_RETURN(status=1)

        ! new dim:
        call obsdata_len_dim%Init( 'obsdata_len', status )
        IF_NOT_OK_RETURN(status=1)
        ! set total size:
        call obsdata_len_dim%Set_Dim( status, n=LEN_LABEL )
        IF_NOT_OK_RETURN(status=1)

        ! new dim:
        call scode_len_dim%Init( 'stncode_len', status )
        IF_NOT_OK_RETURN(status=1)
        ! set total size:
        call scode_len_dim%Set_Dim( status, n=LEN_SCODE )
        IF_NOT_OK_RETURN(status=1)
        
        ! target file:
        write (fname,'("obs_",i4.4,2i2.2,"_",2i2.2,".nc")') &
                 cdate%year, cdate%month, cdate%day, cdate%hour, nint(cdate%seconds/60.0)
        ! create file:
        call F%Create( trim(fname), status )
        IF_NOT_OK_RETURN(status=1)

        ! define coordinate in output file:
        call obs_coor%Def( F, status )
        IF_NOT_OK_RETURN(status=1)
        ! define dimension in output file:
        call obsdata_dim%Def( F, status )
        IF_NOT_OK_RETURN(status=1)
        ! define dimension in output file:
        call obsdata_len_dim%Def( F, status )
        IF_NOT_OK_RETURN(status=1)
        ! define dimension in output file:
        call scode_len_dim%Def( F, status )
        IF_NOT_OK_RETURN(status=1)

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'tracer', NF90_CHAR, &
                                 (/obsdata_len_dim%dimid,obsdata_dim%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'tracer name' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! save:
        varid_name = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'units', NF90_CHAR, &
                                 (/obsdata_len_dim%dimid,obsdata_dim%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'tracer units' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! save:
        varid_units = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'obstype', NF90_CHAR, &
                                 (/obsdata_len_dim%dimid,obsdata_dim%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'observation type' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! save:
        varid_obstype = varid
        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'iobsdata', NF90_INT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'observation dataset index' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_iobsdata = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'sigma_scale', NF90_FLOAT, (/obsdata_dim%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'calibration factor for sigma' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_sigma_scale = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'stnid', NF90_INT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'station id (1-based record number)' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_stnid = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'stncode', NF90_CHAR, &
                                 (/scode_len_dim%dimid,obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'station code' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_scode = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'analyse', NF90_INT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'analysis flag' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'description', '1 : include in analysis, 0 : for validation only' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_analyse = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'longitude', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'standard_name', 'longitude' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', 'degrees_east' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_lon = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'latitude', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'standard_name', 'latitude' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', 'degrees_north' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_lat = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'altitude', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'standard_name', 'altitude' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', 'm' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_alt = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'y', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'observed value' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_y = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'r', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'observation representation error std.dev.' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_r = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'xf', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'simulated value forecast' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_xf = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'sb', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'sigma background' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_sb = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'sf', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'sigma forecast' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_sf = varid

        ! define new variable:
        status = NF90_Def_Var( F%ncid, 'xa', NF90_FLOAT, (/obs_coor%dimid/), varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! extra attributes:
        status = NF90_Put_Att( F%ncid, varid, 'long_name', 'simulated value analysis' )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Put_Att( F%ncid, varid, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! store:
        varid_xa = varid

        ! end of defintions:
        call F%EndDef( status )
        IF_NOT_OK_RETURN(status=1)

        ! storage for collected data:
        allocate( ivalues_all(nobs_all), stat=status )
        IF_NOT_OK_RETURN(status=1)

      else

        ! dummy ...
        allocate( ivalues_all(1), stat=status )
        IF_NOT_OK_RETURN(status=1)

      end if  ! master
      
      ! fill local observation ids:
      do iobs = 1, self%nobs
        ivalues(iobs) = self%obs(iobs)%id
      end do
      ! gather:
      call MPIF90_GatherV( ivalues, self%nobs, &
                           ivalues_all, recvcounts, rdispls, &
                           MasterPE, MPI_COMM_CALC, status )
      IF_MPI_NOT_OK_RETURN(status=1)
      ! write on root:
      if (MasterProc) then
        ! set values:
        call obs_coor%Set_Values( status, values=ivalues_all )
        IF_NOT_OK_RETURN(status=1)
        ! write coordinate:
        call obs_coor%Write( F, status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! write obsdata info from root, same everywhere:
      if (MasterProc) then

        ! loop:
        do iobsdata = 1, nObsData
          ! whitespace:
          label = repeat( ' ', LEN_LABEL )
          ! insert value:
          vlen      = len_trim(obsData(iobsdata)%name)
          label(1:vlen) = trim(obsData(iobsdata)%name)
          ! write padded value:
          status = NF90_Put_Var( F%ncid, varid_name, label, &
                                  start=(/1,iobsdata/), count=(/LEN_LABEL,1/) )
          IF_NF90_NOT_OK_RETURN(status=1)
        end do ! obsdata
        
        ! loop:
        do iobsdata = 1, nObsData
          ! whitespace:
          label = repeat( ' ', LEN_LABEL )
          ! insert value:
          vlen      = len_trim(obsData(iobsdata)%unit)
          label(1:vlen) = trim(obsData(iobsdata)%unit)
          ! write padded value:
          status = NF90_Put_Var( F%ncid, varid_units, label, &
                                  start=(/1,iobsdata/), count=(/LEN_LABEL,1/) )
          IF_NF90_NOT_OK_RETURN(status=1)
        end do ! obsdata

        ! loop:
        do iobsdata = 1, nObsData
          ! whitespace:
          label = repeat( ' ', LEN_LABEL )
          ! insert value:
          vlen      = len_trim(obsData(iobsdata)%deriv)
          label(1:vlen) = trim(obsData(iobsdata)%deriv)
          ! write padded value:
          status = NF90_Put_Var( F%ncid, varid_obstype, label, &
                                  start=(/1,iobsdata/), count=(/LEN_LABEL,1/) )
          IF_NF90_NOT_OK_RETURN(status=1)
        end do ! obsdata
        
        ! write scale factors from root since same on all processes ;
        ! available per observed component in self%BScale(1:nObsComp) ;
        ! loop over target elemtns:
        do iobsdata = 1, nObsData
          ! curent:
          iObsComp = obsData(iobsdata)%iObsComp
          ! put element:
          status = NF90_Put_Var( F%ncid, varid_sigma_scale, (/self%BScale(iObsComp)/), &
                                    start=(/iobsdata/), count=(/1/) )
          IF_NF90_NOT_OK_RETURN(status=1)
        end do ! obsdata

      end if  ! root
      
      ! fill local values:
      do iobs = 1, self%nobs
        ivalues(iobs) = self%obs(iobs)%iObsData
      end do
      ! write:
      call ObsOpers_Write_i( F%ncid, varid_iobsdata, ivalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        ivalues(iobs) = self%obs(iobs)%stnid
      end do
      ! write:
      call ObsOpers_Write_i( F%ncid, varid_stnid, ivalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        ! whitespace:
        scode = repeat( ' ', LEN_SCODE )
        ! insert value:
        vlen          = len_trim(self%obs(iobs)%stncode)
        scode(1:vlen) = trim(self%obs(iobs)%stncode)
        ! store:
        svalues(iobs) = scode
      end do
      ! write:
      call ObsOpers_Write_s( F%ncid, varid_scode, svalues, self%nobs, LEN_SCODE, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        if ( self%obs(iobs)%analyse ) then
          ivalues(iobs) = 1
        else
          ivalues(iobs) = 0
        end if
      end do
      ! write:
      call ObsOpers_Write_i( F%ncid, varid_analyse, ivalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%lon
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_lon, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%lat
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_lat, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%alt
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_alt, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%obs
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_y, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%obsstddev
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_r, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%xf
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_xf, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%sb
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_sb, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%sf
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_sf, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%xa
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_xa, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      IF_NOT_OK_RETURN(status=1)

      ! only root has written:
      if (MasterProc) then
        ! close:
        call F%Close( status )
        IF_NOT_OK_RETURN(status=1)
        ! done:
        call obs_coor%Done( status )
        IF_NOT_OK_RETURN(status=1)
        call obsdata_dim%Done( status )
        IF_NOT_OK_RETURN(status=1)
        call obsdata_len_dim%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end if

      ! clear:
      deallocate( ivalues, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( rvalues, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( svalues, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! clear:
      deallocate( ivalues_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! any obs
    
    ! clear:
    deallocate( recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( rdispls, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine ObsOpers_WriteToFile
  
  ! *
  
  subroutine ObsOpers_Write_i( ncid, varid, values, n, &
                                recvcounts, rdispls, &
                                me, root, comm, status )

    use MPIF90, only : MPIF90_GatherV
    use NetCDF, only : NF90_Put_Var

    ! --- in/out -----------------------------
    
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  varid
    integer, intent(in)               ::  values(:)      ! (n) or dummy if n==0
    integer, intent(in)               ::  n              ! number of valid values
    integer, intent(in)               ::  recvcounts(:)  ! (nproc)
    integer, intent(in)               ::  rdispls   (:)  ! (nproc)
    integer, intent(in)               ::  me
    integer, intent(in)               ::  root
    integer, intent(in)               ::  comm
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Write_i'

    ! --- local -----------------------------
    
    integer, allocatable      ::  values_all(:)
    integer                   ::  n_all
    
    ! --- begin -----------------------------
      
    ! root?
    if ( me == root ) then                          
      ! total number:
      n_all = sum(recvcounts)
      ! storage:
      allocate( values_all(n_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! dummy:
      allocate( values_all(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! gather in 'values_all' on root:
    call MPIF90_GatherV( values, n, &
                         values_all, recvcounts, rdispls, &
                         root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! write on root:
    if ( me == root ) then
      ! write:
      status = NF90_Put_Var( ncid, varid, values_all )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! clear:
    deallocate( values_all, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ObsOpers_Write_i
  
  ! *
  
  subroutine ObsOpers_Write_r( ncid, varid, values, n, &
                                recvcounts, rdispls, &
                                me, root, comm, status )

    use MPIF90, only : MPIF90_GatherV
    use NetCDF, only : NF90_Put_Var

    ! --- in/out -----------------------------
    
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  varid
    real, intent(in)                  ::  values(:)      ! (n) or dummy if n==0
    integer, intent(in)               ::  n              ! number of valid values
    integer, intent(in)               ::  recvcounts(:)  ! (nproc)
    integer, intent(in)               ::  rdispls   (:)  ! (nproc)
    integer, intent(in)               ::  me
    integer, intent(in)               ::  root
    integer, intent(in)               ::  comm
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Write_r'

    ! --- local -----------------------------
    
    real, allocatable         ::  values_all(:)
    integer                   ::  n_all
    
    ! --- begin -----------------------------
      
    ! root?
    if ( me == root ) then                          
      ! total number:
      n_all = sum(recvcounts)
      ! storage:
      allocate( values_all(n_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! dummy:
      allocate( values_all(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! gather in 'values_all' on root:
    call MPIF90_GatherV( values, n, &
                         values_all, recvcounts, rdispls, &
                         root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! write on root:
    if ( me == root ) then
      ! write:
      status = NF90_Put_Var( ncid, varid, values_all )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! clear:
    deallocate( values_all, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ObsOpers_Write_r
  
  ! *
  
  subroutine ObsOpers_Write_s( ncid, varid, values, n, maxlen, &
                                recvcounts, rdispls, &
                                me, root, comm, status )

    use MPIF90, only : MPIF90_GatherV
    use NetCDF, only : NF90_Put_Var

    ! --- in/out -----------------------------
    
    integer, intent(in)               ::  ncid
    integer, intent(in)               ::  varid
    character(len=*), intent(in)      ::  values(:)      ! (n) or dummy if n==0
    integer, intent(in)               ::  n              ! number of valid values
    integer, intent(in)               ::  maxlen         ! max length of character strings
    integer, intent(in)               ::  recvcounts(:)  ! (nproc)
    integer, intent(in)               ::  rdispls   (:)  ! (nproc)
    integer, intent(in)               ::  me
    integer, intent(in)               ::  root
    integer, intent(in)               ::  comm
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Write_s'

    ! --- local -----------------------------
    
    character(len=maxlen), allocatable      ::  values_all(:)
    integer                                 ::  n_all
    
    ! --- begin -----------------------------

    ! root?
    if ( me == root ) then                          
      ! total number:
      n_all = sum(recvcounts)
      ! storage:
      allocate( values_all(n_all), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! dummy:
      allocate( values_all(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! gather in 'values_all' on root:
    call MPIF90_GatherV( values, n, maxlen, &
                         values_all, recvcounts, rdispls, &
                         root, comm, status )
    IF_MPI_NOT_OK_RETURN(status=1)

    ! write on root:
    if ( me == root ) then
      ! write:
      status = NF90_Put_Var( ncid, varid, values_all )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! clear:
    deallocate( values_all, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ObsOpers_Write_s
  
  
  !=========================================================================
  !===
  !=== raw data
  !===
  !=========================================================================


  !
  ! Collected observations (location, values, errors) in 1D arrays.
  !
  
  subroutine read_obs( domain, maxobs, &
                        stnid, flat,flon,falt, y,stddev, scode, analyse, &
                        ipar, nobs, status )
    use GridValues_ml, only : coord_in_domain
    use Io_ml        , only : IO_TMP
    use AOD_PM_ml    , only : AOD_init,wavelength
    use MPI          , only : MPI_IN_PLACE,MPI_INTEGER,MPI_SUM
    use MPI_Groups_ml, only : MPI_COMM_CALC
  
    ! --- in/out ---------------------------------

    character(len=*), intent(in)    ::  domain ! domain/scope for observations
!    "full":      all obs in full model domain
!    "processor": all obs in this processor sub-domain
    integer, intent(in)             ::  maxobs
    integer, intent(out)            ::  stnid(maxobs)  ! 1-based station number
    real, intent(out)               ::  flat(maxobs)
    real, intent(out)               ::  flon(maxobs)
    real, intent(out)               ::  falt(maxobs)
    real, intent(out)               ::  y(maxobs)
    real, intent(out)               ::  stddev(maxobs)
    character(len=*), intent(out)   ::  scode(maxobs)
    logical, intent(out)            ::  analyse(maxobs)
    integer, intent(out)            ::  ipar(maxobs)
    integer, intent(out)            ::  nobs
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/read_obs'
    
    ! --- local ----------------------------------

    integer             ::  nd
    integer             ::  iav
    logical             ::  avflag
    character(len=256)  ::  file
    logical             ::  exist
    integer             ::  no
!    character(len=1024) ::  line
!    real                ::  flon0, flat0, falt0, y0, stddev0
!    character(len=16)   ::  scode0
!    integer             ::  slen
                                            
    ! --- begin ----------------------------------

    ! init counter:
    nobs = 0

    !-----------------------------------------------------------------------
    ! open observational data file in standard format, if applicable
    !-----------------------------------------------------------------------
    ! loop over data sets:
    do nd = 1, nobsData
    
      ! init index range:
      obsData(nd)%iobs(0)=0
      obsData(nd)%iobs(1)=0

      ! loop:
      do iav = 1, 2
      
        ! switch:
        !~ analysed
        if ( iav == 1 ) then
    
          ! undefined ? 
          if ( len_trim(obsData(nd)%file) == 0 ) then
            ! info ...
            write (gol,'("WARNING - no analysis data file specified for Obsdata ",i6)') nd; call goPr
            ! next:
            cycle
          end if
          ! replace time values in template:
          file = date2string( obsData(nd)%file, current_date )
          ! set flag:
          avflag = .true.
          
        !~ validation
        else if ( iav == 2 ) then
    
          ! undefined ? 
          if ( len_trim(obsData(nd)%file_val) == 0 ) then
            ! info ...
            write (gol,'("WARNING - no validation data file specified for Obsdata ",i6)') nd; call goPr
            ! next:
            cycle
          end if
          ! replace time values in template:
          file = date2string( obsData(nd)%file_val, current_date )
          ! set flag:
          avflag = .false.
          
        !~
        else
          write (gol,'("unsupported iav ",i0)') iav; call goErr
          TRACEBACK; status=1; return
        end if        
      
        ! might not be present ...
        inquire( file=trim(file), exist=exist )
        if ( .not. exist ) then
          write (gol,'("WARNING - observation file `",a,"` not found; continue ...")') trim(file); call goPr
          cycle
        end if
        
        ! switch:
        select case ( trim(obsData(nd)%ftype) )
        
          ! good old text files:
          case ( 'txt' )
         
            ! read:
            call read_obs_txt( obsData(nd), nd, avflag, trim(file), &
                                domain, maxobs, &
                                stnid, flat,flon,falt, y,stddev, scode, analyse, ipar, &
                                nobs, status )
            IF_NOT_OK_RETURN(status=1)
        
          ! OMI extract:
          case ( 'omi' )
         
            ! read:
            call read_obs_omi( obsData(nd), nd, avflag, trim(file), &
                                domain, maxobs, &
                                stnid, flat,flon,falt, y,stddev, scode, analyse, ipar, &
                                nobs, status )
            IF_NOT_OK_RETURN(status=1)
            
          ! unknown ...
          case default
            write (gol,'("unsupported file type `",a,"`")') trim(obsData(nd)%ftype); call goErr
            TRACEBACK; status=1; return
        end select ! file type

      end do  ! ana, val

    end do  ! data sets

    !! Write #obs to log file
    !! NOT USEFUL! Only observations at local domain are present ...
    !if ( MasterProc ) then
    !  do nd = 1, nobsData
    !    file = date2string(obsData(nd)%file,current_date)
    !    no = max(obsData(nd)%iobs(1)-obsData(nd)%iobs(0)+1,0)
    !    write (damsg,"('obsData(',I0,') contains ',I0,3(1X,A))") &
    !        nd, no, trim(obsData(nd)%name), &
    !        trim(obsData(nd)%deriv), "observations"
    !    write (damsg,dafmt) trim(damsg)
    !    call PrintLog(damsg)
    !    call PrintLog(file)
    !  end do
    !end if
    
    ! ok
    status = 0

  end subroutine read_obs
  
  
  ! *
  
  
  subroutine read_obs_txt( obsd, nd, avflag, filename, &
                            domain, maxobs, &
                            stnid, flat,flon,falt, y,stddev, scode, analyse, ipar, &
                            nobs, status )

    use GridValues_ml, only : coord_in_domain
    use Io_ml        , only : IO_TMP
  
    ! --- in/out ---------------------------------

    type(obs_data), intent(inout)     ::  obsd  ! data set info
    integer, intent(in)               ::  nd    ! id of 'obsd' argument, stored in 'ipar' array
    logical, intent(in)               ::  avflag  ! analysis or validation set? copied into 'analyse' array
    character(len=*), intent(in)      ::  filename
    character(len=*), intent(in)      ::  domain ! domain/scope for observations
                                                 !    "full":      all obs in full model domain
                                                 !    "processor": all obs in this processor sub-domain
    integer, intent(in)               ::  maxobs
    integer, intent(inout)            ::  stnid(maxobs)  ! 1-based station number
    real, intent(inout)               ::  flat(maxobs)
    real, intent(inout)               ::  flon(maxobs)
    real, intent(inout)               ::  falt(maxobs)
    real, intent(inout)               ::  y(maxobs)
    real, intent(inout)               ::  stddev(maxobs)
    character(len=*), intent(inout)   ::  scode(maxobs)
    logical, intent(inout)            ::  analyse(maxobs)
    integer, intent(inout)            ::  ipar(maxobs)
    integer, intent(inout)            ::  nobs       ! number values filled in (maxobs) arrays on input and output
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/read_obs_text'
    
    ! --- local ----------------------------------

    character(len=1024) ::  line
    integer             ::  no
    real                ::  flon0, flat0, falt0, y0, stddev0
    character(len=16)   ::  scode0
    integer             ::  slen
                                            
    ! --- begin ----------------------------------
  
    ! open:
    open( IO_TMP, file=trim(filename), form='formatted', status='old', iostat=status )
    if ( status /= 0 ) then
      write (gol,'("could not open obsdata file : ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! info ...
    write (gol,'("Obsdata opened ",a)') trim(filename); call goPr

    !-----------------------------------------------------------------------
    ! read data
    !-----------------------------------------------------------------------

    ! init index to collected observation array:
    obsd%iobs(0) = nobs+1

    ! loop over lines:
    do

      ! expected format:
      !     0    48.391670   13.671114      0.00       5.19      0.519 #AT0ENK1

      ! read line:
      read (IO_TMP,'(a)',iostat=status) line
      if ( status < 0 ) then

        ! eof reached
        exit

      elseif(status==0) then

        ! extract parts:
        read (line,*,iostat=status) no,flat0,flon0,falt0,y0,stddev0!,scode0
        IF_NOT_OK_RETURN(status=1)         

        ! read station code, starts with '#' .. (comment on line);
        ! first search for hash:
        no = index(line,'#')
        ! present?
        if ( no > 0 ) then
          ! total lenth:
          slen = len_trim(line)
          ! station code after # to end:
          scode0 = line(no+1:slen)
        else
          ! dummy:
          write (scode0,'("xx",i5.5)') nobs+1
        end if

        ! if ouside domain/scope, next record:
        if ( .not. coord_in_domain(domain,flon0,flat0) ) cycle

        ! check if in range:
        if ( (y0 < obsd%min) .or. (y0 > obsd%max) ) then
          ! next record:
          cycle
        end if
        ! accepted!
        ! increse counter:
        nobs = nobs + 1
        ! check if storage is sufficient:
        if ( nobs > maxobs ) then
          write (gol,'("maxobs too small, cannot store all observations")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! store:
        stnid (nobs) = no + 1  ! covert from 0-based to 1-based
        flat  (nobs) = flat0
        flon  (nobs) = flon0
        falt  (nobs) = falt0
        y     (nobs) = y0
        ! copy station code, skip leading '#'
        scode (nobs) = trim(scode0)
        ! store data set number:
        ipar(nobs) = nd
        ! store analyse/validation flag:
        analyse(nobs) = avflag

        ! how is observation representation error set?
        select case ( trim(obsd%error_type) )
          !~ original:
          case ( 'fraction' )
            ! error estimate read from file (fraction of value?)
            stddev(nobs) = stddev0
            ! truncate std.dev. to relative or absolute maximum:
            stddev(nobs) = max( stddev(nobs), &
                                obsd%error_rel * y(nobs), &
                                obsd%error_rep            )
            ! check .. 
            if ( stddev(nobs) <= 0e0 ) then
              print dafmt,'WARNING obs stddev <= 0'
              stddev(nobs) = 1e-9
            end if
          !~ values per site:
          case ( 'estim' )
            ! extract for current site:
            call ObsError_Lists(nd)%GetError( trim(scode(nobs)), stddev(nobs), status )
            IF_NOT_OK_RETURN(status=1)
          !~
          case default
            write (gol,'("unsupported obs.repr.error type `",a,"`")') trim(obsd%error_type); call goErr
            TRACEBACK; status=1; return
        end select

        !! testing ...
        !write (gol,'(a,": observation ",2i6,2f8.2)') rname, nd, nobs, flon(nobs), flat(nobs); call goPr

      else
        write (gol,'("while reading observations")'); call goErr
        TRACEBACK; status=1; return
      end if

    end do  ! lines

    ! store end number:
    obsd%iobs(1) = nobs

    ! close:
    close( IO_TMP, iostat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
        
  end subroutine read_obs_txt
  
  
  ! *
  
  
  subroutine read_obs_omi( obsd, nd, avflag, filename, &
                            domain, maxobs, &
                            stnid, flat,flon,falt, y,stddev, scode, analyse, ipar, &
                            nobs, status )

    use GridValues_ml, only : coord_in_domain
    use EMIP_OMI     , only : T_EMIP_OMI
  
    ! --- in/out ---------------------------------

    type(obs_data), intent(inout)     ::  obsd  ! data set info
    integer, intent(in)               ::  nd    ! id of 'obsd' argument, stored in 'ipar' array
    logical, intent(in)               ::  avflag  ! analysis or validation set? copied into 'analyse' array
    character(len=*), intent(in)      ::  filename
    character(len=*), intent(in)      ::  domain ! domain/scope for observations
                                                 !    "full":      all obs in full model domain
                                                 !    "processor": all obs in this processor sub-domain
    integer, intent(in)               ::  maxobs
    integer, intent(inout)            ::  stnid(maxobs)  ! 1-based station number
    real, intent(inout)               ::  flat(maxobs)
    real, intent(inout)               ::  flon(maxobs)
    real, intent(inout)               ::  falt(maxobs)
    real, intent(inout)               ::  y(maxobs)
    real, intent(inout)               ::  stddev(maxobs)
    character(len=*), intent(inout)   ::  scode(maxobs)
    logical, intent(inout)            ::  analyse(maxobs)
    integer, intent(inout)            ::  ipar(maxobs)
    integer, intent(inout)            ::  nobs       ! number values filled in (maxobs) arrays on input and output
    integer, intent(out)              ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/read_obs_omi'
    
    ! minimum allowed retrieval:
    real, parameter   ::  y0_minval = -0.2

    ! --- local ----------------------------------

    type(T_EMIP_OMI)    ::  df
    integer             ::  ipix
    real                ::  flon0, flat0
    real                ::  y0, stddev0
                                            
    ! --- begin ----------------------------------
  
    ! info ...
    write (gol,'(a,": open ",a)') rname, trim(filename); call goPr
    
    ! init index to collected observation array:
    obsd%iobs(0) = nobs+1

    ! read file content:
    call df%Init( filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy kernels:
    call ObsKernels(nd)%Fill( df%kernel, df%kernel_units, &
                              df%phlev, df%phlev_units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ..
    write (gol,'(a,":   number of pixels: ",i0)') rname, df%npixel; call goPr
    write (gol,'(a,":   skip retrievals below ",f8.2)') rname, y0_minval; call goPr

    ! loop over pixels:
    do ipix = 1, df%npixel
    
      ! extract:
      flon0   = df%longitude(ipix)
      flat0   = df%latitude (ipix)
      y0      = df%vcd_trop (ipix)
      stddev0 = df%sigma_vcd_trop(ipix)

      ! if ouside domain/scope, next record:
      if ( .not. coord_in_domain(domain,flon0,flat0) ) cycle

      ! OMI NO2 could be slightly negative, but some values are really strange ...
      ! skip the values that are very negative:
      if ( y0 < y0_minval ) then
        ! next record:
        cycle
      end if

      ! accepted!
      ! increse counter:
      nobs = nobs + 1
      ! check if storage is sufficient:
      if ( nobs > maxobs ) then
        write (gol,'("maxobs too small, cannot store all observations")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! store:
      stnid (nobs) = ipix  ! use pixel number as 'station' id
      flat  (nobs) = flat0
      flon  (nobs) = flon0
      falt  (nobs) = 0.0   ! dummy
      y     (nobs) = y0
      ! use 'station' code with pixel number:
      write (scode(nobs),'("p",i5.5)') ipix
      ! store data set number:
      ipar(nobs) = nd
      ! store analyse/validation flag:
      analyse(nobs) = avflag

      ! how is observation representation error set?
      select case ( trim(obsd%error_type) )
        !~ as in product:
        case ( 'retrieval' )
          ! copy:
          stddev(nobs) = stddev0
        !~ relative to observed value:
        case ( 'fraction' )
          ! error estimate read from file (fraction of value?)
          stddev(nobs) = stddev0
          ! truncate std.dev. to relative or absolute maximum:
          stddev(nobs) = max( stddev(nobs), &
                              obsd%error_rel * y(nobs), &
                              obsd%error_rep            )
          ! check .. 
          if ( stddev(nobs) <= 0e0 ) then
            print dafmt,'WARNING obs stddev <= 0'
            stddev(nobs) = 1e-9
          end if
        !~
        case default
          write (gol,'("unsupported obs.repr.error type `",a,"`")') trim(obsd%error_type); call goErr
          TRACEBACK; status=1; return
      end select
      
    end do  ! pixels

    ! store end number:
    obsd%iobs(1) = nobs

    ! clear:
    call df%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
        
  end subroutine read_obs_omi

end module DA_Obs_ml
