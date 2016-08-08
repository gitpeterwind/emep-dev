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
  use ChemChemicals_ml, only: species
  use ChemGroups_ml,    only: chemgroups
  use Io_ml,            only: IO_TMP
  use Io_Progs_ml,      only: PrintLog
  use ModelConstants_ml,only: KMAX_MID,MasterProc,runlabel1
  use SmallUtils_ml,    only: find_index
  use TimeDate_ExtraUtil_ml, only: date2string
  use Units_ml        , only : Units_Scale, Group_Units
  use DA_ml           , only : debug => DEBUG_DA
  use DA_ml           , only : dafmt => da_fmt_msg
  use DA_ml           , only : damsg => da_msg
  use DA_ml           , only : FGSCALE
  use DA_ml           , only : nlev
  use DA_ml           , only : nChem, iChemInv
  use DA_ml           , only : nChemObs, iChemObs
  use ModelConstants_ml, only: nproc
  use Par_ml,            only: me

  implicit none

  
  ! --- const -----------------------------------------

  character(len=*), parameter  ::  mname = 'DA_Obs_ml'
  
  ! max length:
  integer, parameter  ::  LEN_LABEL = 64
  integer, parameter  ::  LEN_SCODE = 8
  
  
  ! --- in/out ---------------------------------------
  
  private

  public    ::  Read_Obs
  
  public    ::  T_ObsOper
  public    ::  T_ObsOpers
  
  public    ::  varName, varSpec, varSpecInv
  public    ::  obsVarName, observedVar
  public    ::  OBSERVATIONS, nobsData, obsData
  public    ::  nchemobs, ichemObs
  

  ! --- types ----------------------------------------
  
  type obs_data
    logical            :: set       = .false.
    logical            :: found     = .false.
    logical            :: unitroa   = .false.
    integer            :: ispec     = -1
    integer            :: ichem     = -1
    integer            :: ichemObs  = -1
    integer            :: iobs(0:1) = -1
    real               :: min       = 0.0
    real               :: max       = 1e3
    real               :: error_rel = -1e0
    real               :: error_rep = -1e0
    real               :: unitconv  = -1e0
    character(len=LEN_LABEL)  :: name      = ''    ! tracer name
    character(len=LEN_LABEL)  :: unit      = ''
    character(len=LEN_LABEL)  :: deriv     = ''
    character(len=256)        :: file      = ''
    character(len=256)        :: tag       = ''
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

  ! observation operator :
  ! NOTE: when changing this, adapt the 'ObsOpers_Swap' routine too!
  type T_ObsOper
    ! observation and std.dev.:
    real                        ::  obs
    real                        ::  obsstddev
    ! model simulations:
    real                        ::  xf       ! forecast
    real                        ::  xa       ! analysis
    ! innovation:
    real                        ::  innov    ! obs - forecast
    ! Jacobian dH/dx :
    real, allocatable           ::  H_jac(:,:)  ! (nlev,nChemObs)
    ! involved cells:
    integer                     ::  i
    integer                     ::  j
    integer                     ::  l(0:1)
    ! id for testing:
    integer                     ::  id
    ! index in 'obsData' array:
    integer                     ::  iObsData
    ! meta data:
    integer                     ::  stnid
    real                        ::  lon, lat
    real                        ::  alt
    character(len=LEN_SCODE)    ::  stncode
    ! copy of values from original observation info:
    integer                     ::  ichemObs
    ! ... or new index:
    integer                     ::  itracer
  contains
    procedure   ::  Init         => ObsOper_Init
    procedure   ::  Done         => ObsOper_Done
    procedure   ::  FillTracers  => ObsOper_FillTracers
    procedure   ::  Fill         => ObsOper_Fill
    procedure   ::  Evaluate     => ObsOper_Evaluate
  end type T_ObsOper

  ! observation operators:
  type T_ObsOpers
    integer                         ::  nobs
    type(T_ObsOper), allocatable    ::  obs(:)
  contains
    procedure   ::  Init          => ObsOpers_Init
    procedure   ::  Done          => ObsOpers_Done
    procedure   ::  Alloc         => ObsOpers_Alloc
    procedure   ::  DeAlloc       => ObsOpers_DeAlloc
    procedure   ::  Set_IDs       => ObsOpers_Set_IDs
    procedure   ::  Swap          => ObsOpers_Swap
    procedure   ::  SelectTracers => ObsOpers_SelectTracers
    procedure   ::  Evaluate      => ObsOpers_Evaluate
    procedure   ::  WriteToFile   => ObsOpers_WriteToFile
  end type T_ObsOpers


  !--- local -------------------------------------------------------------

  integer            :: nobsData=0
  integer, parameter :: nobsDataMax=50
  type(obs_data)     :: obsData(nobsDataMax)

  namelist /OBSERVATIONS/nobsData,obsData

  !-----------------------------------------------------------------------

  integer               :: varSpec(NSPEC_ADV),varSpecInv(NSPEC_ADV)
  character(len=016)    :: varName(NSPEC_ADV)='',obsVarName(NSPEC_ADV)=''
  logical               :: observedVar(NSPEC_ADV)=.false.


contains


  !=========================================================================
  !===
  !=== observation operators
  !===
  !=========================================================================


  subroutine ObsOper_Init( self, nchm, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOper), intent(out)   ::  self
    integer, intent(in)             ::  nchm  ! (nChemObs) or (ntracer)
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOper_Init'

    ! --- begin -----------------------------
    
    ! storage:
    allocate( self%H_jac(nlev,nchm), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! dummy values:
    self%H_jac     = -999.9
    
    ! dummy values:
    self%obs       = -999.9
    self%xf        = -999.9
    self%xa        = -999.9
    self%innov     = -999.9
    self%obsstddev = -999.9
    self%i         = -999
    self%j         = -999
    self%l         = -999
    self%id        = -999
    self%iObsData  = -999
    self%lon       = -999.9
    self%lat       = -999.9
    self%alt       = -999.9
    self%iChemObs  = -999
    self%itracer   = -999

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


  ! copy tracer selection from 'nChemObs' arrays to 'ntracer' arrays
  
  subroutine ObsOper_FillTracers( self, src, itracers, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOper), intent(inout)   ::  self
    class(T_ObsOper), intent(in)      ::  src      ! with iChemObs and H_jac(:,:,nChemObs)
    integer, intent(in)               ::  itracers(:)   ! (nChemObs)
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOper_FillTracers'

    ! --- local -----------------------------
    
    integer       ::  ichemobs
    integer       ::  itracer

    ! --- begin -----------------------------
    
    ! copy scalar values:
    self%obs       = src%obs
    self%xf        = src%xf
    self%xa        = src%xa
    self%innov     = src%innov
    self%obsstddev = src%obsstddev
    self%i         = src%i
    self%j         = src%j
    self%l         = src%l
    self%id        = src%id
    self%iObsData  = src%iObsData
    self%lon       = src%lon
    self%lat       = src%lat
    self%alt       = src%alt
    ! no target obs anymore:
    self%iChemObs   = -999
    ! instead, new index:
    self%itracer    = itracers(src%iChemObs)
    
    ! loop over original chemobs:
    do ichemobs = 1, nChemObs
      ! target tracer:
      itracer = itracers(ichemobs)
      ! defined?
      if ( itracer > 0 ) then
        ! copy:
        self%H_jac(:,itracer) = src%H_jac(:,ichemobs)
      end if
    end do ! original tracers

    ! ok
    status = 0

  end subroutine ObsOper_FillTracers


  ! ***


subroutine ObsOper_Fill( self, xn_adv_b, &
                          ipar, stnid, lat,lon, alt, &
                          yn, status )
!-----------------------------------------------------------------------
! @description
! H operator mapping from model to observation space;
! Further: interpolate air density to observation point
!
! The following 4D concentration fields are currently used;
!
!   xn_adv_b :  Background concentration on local grid ("xn_adv");
!               all tracers, used for indirect observations only.
!               Shape: (nspec_adv,lnx,lny,kmax_mid)
!               Used from module 'Chemfields_ml'.
!               Include scale factor 'FGSCALE' when using elements
!               of this array.
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
!-----------------------------------------------------------------------

use GridValues_ml, only: glon,glat,coord_in_domain
use MetFields_ml , only: z_bnd,roa
use AOD_PM_ml,     only: aod_grp,SpecExtCross
use ChemFields_ml, only: cfac,PM25_water,pm25_water_rh50

! --- in/out ----------------------------

class(T_ObsOper), intent(inout)    ::  self
real, intent(in)              :: xn_adv_b(:,:,:,:)  ! (nspec_adv,lnx,lny,kmax_mid)
integer, intent(in)           :: ipar         ! index in 'obsData' array
integer, intent(in)           :: stnid        ! station id (1-based record number in file)
real, intent(inout)           :: lat,lon      ! location
real, intent(in)              :: alt          ! location
real, intent(out)             :: yn           ! simulation
integer, intent(out)          :: status

! --- const ----------------------------

character(len=*), parameter  ::  rname = mname//'/ObsOper_Fill'

! --- local -----------------------------

integer :: i,j,l,ll,k,kk,igrp,g,ispec,nn
real    :: xn,Hj,T1,QML1,PS1
real    :: Hchem(nchem)
real    :: unitconv
integer, pointer, dimension(:) :: gspec=>null()      ! group array of indexes
real,    pointer, dimension(:) :: gunit_conv=>null() ! group array of unit conv. factors
real                           :: z_middle
logical :: local_model_grid

! --- begin -----------------------------
    
  ! store meta data:
  self%stnid = stnid
  self%lon   = lon
  self%lat   = lat
  self%alt   = alt
  self%stncode = ''
  ! store index in 'obsData' array:
  self%iObsData = ipar
  ! get local grid indices
  local_model_grid=coord_in_domain("processor",lon,lat,iloc=i,jloc=j)
  self%i=i
  self%j=j
  !! testing ...
  !!write (gol,'(a,": H_op ",2i6,2f8.2,2i6)') rname, ipar, self%id, lon, lat, i, j; call goPr

  ! check, expecting currently local indices only!
  if(.not.local_model_grid) then
    write (gol,'("found interpolation indices outside local domain:")'); call goErr
    write (gol,'("  local grid shape : ",2i4)') size(xn_adv_b,2), size(xn_adv_b,3); call goErr
    write (gol,'("  i index : ",i4)') i; call goErr
    write (gol,'("  j index : ",i4)') j; call goErr
    TRACEBACK; status=1; return
  end if
  
  !-----------------------------------------------------------------------
  ! z grid coordinate of observations:
  !-----------------------------------------------------------------------
  ll=KMAX_MID
  do while(alt>z_bnd(i,j,ll).and.ll>=1)
    ll=ll-1
  enddo
  l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
  select case(obsData(ipar)%deriv)
    case('mod-lev'  );self%l(:)=l
    case('surface'  );self%l(:)=l
    case('Trop.Col.');self%l(:)=[l,nlev]
    case default
      write (gol,'("obsData not yet supported")'); call goErr
      TRACEBACK; status=1; return
  endselect

  !! testing ...
  !write (gol,'(a,":   level ll=",i0," l=",i0)') rname, ll, l; call goPr


  !-----------------------------------------------------------------------
  ! info
  !-----------------------------------------------------------------------
  if(debug)then
    ! cell with observation location:
    ! adhoc, original 'z_mid' is not swapped ...
    z_middle = 0.5*( z_bnd(i,j,ll) + z_bnd(i,j,ll+1) )
    ! info ...
    write (gol,dafmt) 'Observation/Model Location:Weight[%],lon[deg],lat[deg],alt[km]'; call goPr
    write (gol,"(5(1X,A3,':',F0.2,',',F0.2,',',F0.2))") &
      'Observation',lon,lat,alt,&
      'Model',glon(i,j),glat(i,j),z_middle*1e-3; call goPr
  endif

  ! obsData element should be setup in ObsOper_Fill
  call CheckStop(.not.obsData(ipar)%set,&
    "obsData not set: "//trim(obsData(ipar)%name))

  ! copy value from obsData element into data type:
  self%iChemObs = obsData(ipar)%iChemObs
  
  !! info ...
  !write (gol,'(a,":   shape(xn_b)=",4i4)') rname, shape(xn_b); call goPr
  !write (gol,'(a,":   deriv=",a)') rname, trim(obsData(ipar)%deriv); call goPr

  !if ( n == obsData(ipar)%iobs(0) ) then
  !  print dafmt,"Analysis of "//trim(obsData(ipar)%tag)
  !end if

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  yn=0e0
!-----------------------------------------------------------------------
! analysis of direct observations: NO2, O3, SO2, etc
!-----------------------------------------------------------------------
  if(obsData(ipar)%found)then
    unitconv=obsData(ipar)%unitconv
    ispec=obsData(ipar)%ispec
    k =obsData(ipar)%ichem
    select case(obsData(ipar)%deriv)
!-----------------------------------------------------------------------
! direct observations: model mid-levels
!-----------------------------------------------------------------------
    case('mod-lev')
      if(obsData(ipar)%unitroa)then   ! unit conversion        
        Hj=unitconv*roa(i,j,ll,1)
      else
        Hj=unitconv
      endif
      xn=xn_adv_b(ispec,i,j,ll)*FGSCALE
      yn=yn+Hj*xn
      self%H_jac(l,k)=Hj
      !! testing ...
      !write (gol,'(a,":   lev H_jac(",3i4,") * xn(",4i4,") = ",2es12.4," = ",es12.4)') &
      !        rname, l,k, i,j,l,kk, self%H_jac(l,k), xn, self%H_jac(l,k)*xn; call goPr
      !write (gol,'(a,":     yn = ",es12.4)') rname, yn; call goPr
!-----------------------------------------------------------------------
! direct observations: surface level
!-----------------------------------------------------------------------
    case('surface')
      if(obsData(ipar)%unitroa)then  ! unit conversion,50m->3m
        Hj=unitconv*cfac(ispec,i,j)*roa(i,j,ll,1)
      else
        Hj=unitconv*cfac(ispec,i,j)
      endif
      xn=xn_adv_b(ispec,i,j,ll)*FGSCALE
      yn=yn+Hj*xn
      self%H_jac(l,k)=Hj
      !! testing ...
      !write (gol,'(a,":   sfc H_jac(",3i4,") * xn(",4i4,") = ",2es12.4," = ",es12.4)') &
      !        rname, p,l,k, i,j,l,kk, self%H_jac(l,k), xn, self%H_jac(l,k)*xn; call goPr
      !write (gol,'(a,":     yn = ",es12.4)') rname, yn; call goPr
!-----------------------------------------------------------------------
! direct observations: COLUMN (eg NO2tc)
!-----------------------------------------------------------------------
    case('Trop.Col.')
      call CheckStop(.not.obsData(ipar)%unitroa,&
        "Unsupported obsData unit for: "//trim(obsData(ipar)%tag))     
      do ll=1,KMAX_MID
        Hj=unitconv                       & ! unit conversion
          * roa(i,j,ll,1)                 & ! density.
          * (z_bnd(i,j,ll)-z_bnd(i,j,ll+1)) ! level thickness
        xn=xn_adv_b(ispec,i,j,ll)*FGSCALE
        yn=yn+Hj*xn
        l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
        if(l>0)self%H_jac(l,k)=Hj
      enddo
    endselect
!-----------------------------------------------------------------------
! analysis of indirect observations: PM* (eg PM10)
!-----------------------------------------------------------------------
  elseif(index(obsData(ipar)%name,'PM')>0)then
    igrp=obsData(ipar)%ispec
!!  print *,igrp,chemgroups(igrp)%name
    call Group_Units(igrp,obsData(ipar)%unit,gspec,gunit_conv,debug)
    select case(obsData(ipar)%deriv)
!-----------------------------------------------------------------------
! indirect observations: model mid-levels
!-----------------------------------------------------------------------
    case('mod-lev')
      if(obsData(ipar)%unitroa)&    ! unit conversion
        gunit_conv(:)=gunit_conv(:)*roa(i,j,ll,1)
      do g=1,size(gspec)
        Hj=gunit_conv(g)
        xn=xn_adv_b(gspec(g),i,j,ll)*FGSCALE
        yn=yn+Hj*xn
        kk=varSpecInv(gspec(g))
        k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
        if(k>0)self%H_jac(l,k)=Hj
      enddo
      ! PM water content in ug/m3 at model mid-levels
      if(index(obsData(ipar)%name,'PM')>0)&
        yn=yn+PM25_water(i,j,ll)*FGSCALE
!-----------------------------------------------------------------------
! indirect observations: surface level
!-----------------------------------------------------------------------
    case('surface')
      if(obsData(ipar)%unitroa)&    ! unit conversion
        gunit_conv(:)=gunit_conv(:)*roa(i,j,ll,1)
      do g=1,size(gspec)
        Hj=gunit_conv(g)          & ! unit conversion
          *cfac(gspec(g),i,j)       ! 50m->3m conversion
        xn=xn_adv_b(gspec(g),i,j,ll)*FGSCALE
        yn=yn+Hj*xn
        kk=varSpecInv(gspec(g))
        k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
        if(k>0)self%H_jac(l,k)=Hj
      enddo
      ! PM water content in ug/m3 at surface level
      if(index(obsData(ipar)%name,'PM')>0)&
        yn=yn+PM25_water_rh50(i,j)*FGSCALE
!-----------------------------------------------------------------------
! indirect observations: COLUMN
!-----------------------------------------------------------------------
    case('Trop.Col.')
      do ll=1,KMAX_MID
        Hj=roa(i,j,ll,1)                 & ! interpolation,density
          *(z_bnd(i,j,ll)-z_bnd(i,j,ll+1)) ! level thickness
        l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
        do g=1,size(gspec)
          xn=xn_adv_b(gspec(g),i,j,ll)*gunit_conv(g)*FGSCALE
          yn=yn+Hj*xn
          kk=varSpecInv(gspec(g))
          k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
          if(k>0.and.l>0)self%H_jac(l,k)=Hj*gunit_conv(g)
        enddo
      enddo
      ! PM water content in ug/m3 at model mid-levels
      if(index(obsData(ipar)%name,'PM')>0)then
        xn=0.0
        do ll=1,KMAX_MID
          xn=xn+PM25_water(i,j,ll)&
            *(z_bnd(i,j,ll)-z_bnd(i,j,ll+1)) ! level thickness
        enddo
        yn=yn+xn*FGSCALE
      endif
    endselect
!-----------------------------------------------------------------------
! analysis of indirect observations:
!   BSC and EXT GROUPs (wavelength in %unit and wlen index in %ichem)
!-----------------------------------------------------------------------
! elseif(any(obsData(ipar)%name==['BSC','EXT','AOD'])then
  elseif(any(obsData(ipar)%name==['EXT','AOD']))then
    gspec=>aod_grp          ! AOD group
    nn=obsData(ipar)%ichem  ! wavelength index
    select case(obsData(ipar)%deriv)
!-----------------------------------------------------------------------
! indirect observations: model mid-levels
!-----------------------------------------------------------------------
    case('mod-lev')
      gunit_conv=>SpecExtCross(:,ll,i,j,nn)
      do g=1,size(gspec)
        Hj=gunit_conv(g)
        xn=xn_adv_b(gspec(g),i,j,ll)*FGSCALE
        yn=yn+Hj*xn
        kk=varSpecInv(gspec(g))
        k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
        if(k>0)self%H_jac(l,k)=Hj
      enddo
!-----------------------------------------------------------------------
! indirect observations: surface level
!-----------------------------------------------------------------------
    case('surface')
      gunit_conv=>SpecExtCross(:,ll,i,j,nn)
      do g=1,size(gspec)
        Hj=gunit_conv(g)          & ! unit conversion
          *cfac(gspec(g),i,j)       ! 50m->3m conversion
        xn=xn_adv_b(gspec(g),i,j,ll)*FGSCALE
        yn=yn+Hj*xn
        kk=varSpecInv(gspec(g))
        k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
        if(k>0)self%H_jac(l,k)=Hj
      enddo
!-----------------------------------------------------------------------
! indirect observations: COLUMN (eg AOD)
!-----------------------------------------------------------------------
    case('Trop.Col.')
      do ll=1,KMAX_MID
        Hj=(z_bnd(i,j,ll)-z_bnd(i,j,ll+1)) ! level thickness
        l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
        gunit_conv=>SpecExtCross(:,ll,i,j,nn)
        do g=1,size(gspec)
          xn=xn_adv_b(gspec(g),i,j,ll)*gunit_conv(g)*FGSCALE
          yn=yn+Hj*xn
          kk=varSpecInv(gspec(g))
          k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
          if(k>0.and.l>0)self%H_jac(l,k)=Hj*gunit_conv(g)
        enddo
      enddo
    endselect
!-----------------------------------------------------------------------
  else
    write (gol,'("obsData not yet supported")'); call goErr
    TRACEBACK; status=1; return
  endif
!-----------------------------------------------------------------------
  status = 0
endsubroutine ObsOper_Fill


  ! ***


  subroutine ObsOper_Evaluate( self, key, xn_b, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOper), intent(inout)   ::  self
    character(len=*), intent(in)      ::  key
    real, intent(in)                  ::  xn_b(:,:,:,:)      ! (lnx,lny,nlev,nchem)
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOper_Evaluate'

    ! --- local -----------------------------

    integer         ::  i, j
    integer         ::  l0, l1
    integer         ::  ichem
    real            ::  yn
    
    ! --- begin -----------------------------

    ! init sum:
    yn = 0.0e0
    ! involved grid cell:
    i = self%i
    j = self%j
    ! level range:
    l0 = self%l(0)
    l1 = self%l(1)
    ! observed tracer index:
    ichem = self%ichemObs
    ! add contribution for this cell, sum over layers:
    yn = yn + sum(self%H_jac(l0:l1,ichem) * xn_b(i,j,l0:l1,ichem))

    ! store:
    select case ( key )
      case ( 'xf' )
        self%xf = yn
      case ( 'xa' )
        self%xa = yn
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
    
    ! clear:
    call self%DeAlloc( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine ObsOpers_Done


  ! ***


  subroutine ObsOpers_Alloc( self, nobs, nchm, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(out)  ::  self
    integer, intent(in)             ::  nobs
    integer, intent(in)             ::  nchm
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
        call ObsOper_Init( self%obs(iobs), nchm, status )
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
  
    use MPI, only : MPI_INTEGER
    use MPI, only : MPI_AllGather
    use MPIF90, only : MPIF90_Displacements
    
    use MPI_Groups_ml, only : MPI_COMM_CALC

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
      write (gol,'(a,": proc observations : ",i0," - ",i0)') rname, self%obs(1)%id, self%obs(self%nobs)%id; call goPr
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
    use DA_Util_ml   , only : T_Domains
    use MPI_Groups_ml, only : MPI_COMM_CALC

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
    integer                         ::  nrr, sHj(2), nii
    real, allocatable               ::  srt_rr(:,:),srt_Hj(:,:,:),srt_ii(:,:)
    real, allocatable               ::  new_rr(:,:),new_Hj(:,:,:),new_ii(:,:)
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
    
    !! info ...
    !do iobs = 1, self%nobs
    !  write (gol,*) rname//': obs ', self%obs(iobs)%id, &
    !    self%obs(iobs)%innov, self%obs(iobs)%obsstddev, ';', &
    !    sum(self%obs(iobs)%H_jac), ';', &
    !    self%obs(iobs)%i(1), self%obs(iobs)%j(1), self%obs(iobs)%l; call goPr
    !end do
    
    ! storage for sorted version that can be passed through MPI;
    ! decomposition over observation, so this should be last index:
    nrr = 2
    allocate( srt_rr(nrr,self%nobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    sHj = [ nlev, nchemobs ]
    allocate( srt_Hj(sHj(1),sHj(2),self%nobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    nii = 6
    allocate( srt_ii(nii,self%nobs), stat=status )
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
          srt_Hj(:,:,iobs) = self%obs(iobs_in)%H_jac
          srt_ii(1  ,iobs) = doms%off(1,me) + self%obs(iobs_in)%i
          srt_ii(2  ,iobs) = doms%off(2,me) + self%obs(iobs_in)%j
          srt_ii(3:4,iobs) = self%obs(iobs_in)%l
          srt_ii(5  ,iobs) = self%obs(iobs_in)%id
          srt_ii(6  ,iobs) = self%obs(iobs_in)%iChemObs
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
    
    !! info ...
    !write (gol,'(a,": exchange H ...")') rname; call goPr
    !write (gol,*) rname//':    sendcounts   : ', sendcounts; call goPr
    !write (gol,*) rname//':    sdispls      : ', sdispls; call goPr
    !write (gol,*) rname//':    recvcounts   : ', recvcounts; call goPr
    !write (gol,*) rname//':    rdispls      : ', rdispls; call goPr

    ! total number to be received:
    nobs_new = sum(recvcounts)

    ! storage for arrays to be received; 
    ! allocate at least something to avoid bound check errors:
    nval = max( 1, nobs_new )
    allocate( new_rr(nrr,nval), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( new_Hj(sHj(1),sHj(2),nval), stat=status )
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
    nval = product( sHj )
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
    call Hops_new%Alloc( nobs_new, nChemObs, status )
    IF_NOT_OK_RETURN(status=1)
    ! unpack:
    do iobs = 1, nobs_new
      ! copy values, convert from global to local horizontal indices:
      Hops_new%obs(iobs)%innov     = new_rr(1,iobs)    
      Hops_new%obs(iobs)%obsstddev = new_rr(2,iobs)    
      Hops_new%obs(iobs)%H_jac     = new_Hj(:,:,iobs)
      Hops_new%obs(iobs)%i         = new_ii(1  ,iobs) - doms_new%off(1,me)
      Hops_new%obs(iobs)%j         = new_ii(2  ,iobs) - doms_new%off(2,me)
      Hops_new%obs(iobs)%l         = new_ii(3:4,iobs)
      Hops_new%obs(iobs)%id        = new_ii(5  ,iobs)
      Hops_new%obs(iobs)%iChemObs  = new_ii(6  ,iobs)
    end do
    
    !! info ...
    !do iobs = 1, Hops_new%nobs
    !  write (gol,*) rname//': obs ', Hops_new%obs(iobs)%id, &
    !    Hops_new%obs(iobs)%innov, Hops_new%obs(iobs)%obsstddev, ';', &
    !    sum(Hops_new%obs(iobs)%H_jac), ';', &
    !    Hops_new%obs(iobs)%i(1), Hops_new%obs(iobs)%j(1), Hops_new%obs(iobs)%l; call goPr
    !end do
    
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


  subroutine ObsOpers_SelectTracers( self, H, itracers, status )

  !-----------------------------------------------------------------------
  ! @description
  ! Fill (already initialized) observation structure
  ! with subset of self with only selected tracers
  !-----------------------------------------------------------------------
  
    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(inout)  ::  self
    type(T_ObsOpers), intent(in)      ::  H
    integer, intent(in)               ::  itracers(:)   ! (nChemObs)
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Swap'

    ! --- local -----------------------------
    
    integer                         ::  ntracer
    integer                         ::  nobs
    integer                         ::  iobs
    integer                         ::  k
    
    ! --- begin -----------------------------
    
    ! check ...
    if ( size(itracers) /= nChemObs ) then
      write (gol,'("size of itracer is ",i0," while nChemObs is ",i0)') size(itracers), nChemObs; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! count:
    ntracer = count( itracers > 0 )
    
    ! init new number of observations:
    nobs = 0
    ! loop over existing observations:
    do iobs = 1, H%nobs
      ! check index in 'nChemObs' arrays:
      if ( itracers(H%obs(iobs)%iChemObs) > 0 ) nobs = nobs + 1
    end do
    
    ! storage for structures:
    call self%Alloc( nobs, ntracer, status )
    IF_NOT_OK_RETURN(status=1)
    ! init new index:
    iobs = 0
    ! loop over existing observations:
    do k = 1, H%nobs
      ! check index in 'nChemObs' arrays:
      if ( itracers(H%obs(k)%iChemObs) > 0 ) then
        ! increase index:
        iobs = iobs + 1
        ! copy values:
        call self%obs(iobs)%FillTracers( H%obs(k), itracers, status )
        IF_NOT_OK_RETURN(status=1)
      end if ! selected?
    end do ! original obs

    ! ok
    status = 0

  end subroutine ObsOpers_SelectTracers


  ! ***


  subroutine ObsOpers_Evaluate( self, key, xn_b, status )

    ! --- in/out -----------------------------
    
    class(T_ObsOpers), intent(inout)  ::  self
    character(len=*), intent(in)      ::  key
    real, intent(in)                  ::  xn_b(:,:,:,:)      ! (lnx,lny,nlev,nchem)
    integer, intent(out)              ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Evaluate'

    ! --- local -----------------------------
    
    integer         ::  n
    
    ! --- begin -----------------------------

    ! any local observations?
    if ( self%nobs > 0 ) then
      ! loop over observations:
      do n = 1, self%nobs
        ! evaluate:
        call self%obs(n)%Evaluate( key, xn_b, status )
        IF_NOT_OK_RETURN(status=1)
      end do ! n
    end if  ! nobs > 0
    
    ! ok
    status = 0
    
  end subroutine ObsOpers_Evaluate


  ! ***


  subroutine ObsOpers_WriteToFile( self, cdate, status )

    use MPIF90           , only : MPIF90_AllGather
    use MPIF90           , only : MPIF90_Displacements
    use MPIF90           , only : MPIF90_GatherV
    use MPI_Groups_ml    , only : MPI_COMM_CALC
    use MPI_Groups_ml    , only : MasterPE
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
    integer                   ::  varid_lon, varid_lat, varid_alt
    integer                   ::  varid_y, varid_r
    integer                   ::  varid_xf, varid_xa

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

      ! any local obs?
      if ( self%nobs > 0 ) then
        ! storage for local data:
        allocate( ivalues(self%nobs), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( rvalues(self%nobs), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( svalues(self%nobs), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! dummy ...
        allocate( ivalues(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( rvalues(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( svalues(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
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

      end if  ! root
      
      ! fill local values:
      do iobs = 1, self%nobs
        ivalues(iobs) = self%obs(iobs)%iObsData
      end do
      ! write:
      call ObsOpers_Write_i( F%ncid, varid_iobsdata, ivalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      
      ! fill local values:
      do iobs = 1, self%nobs
        ivalues(iobs) = self%obs(iobs)%stnid
      end do
      ! write:
      call ObsOpers_Write_i( F%ncid, varid_stnid, ivalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      
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
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%lon
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_lon, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%lat
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_lat, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%alt
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_alt, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%obs
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_y, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%obsstddev
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_r, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%xf
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_xf, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )
      
      ! fill local values:
      do iobs = 1, self%nobs
        rvalues(iobs) = self%obs(iobs)%xa
      end do
      ! write:
      call ObsOpers_Write_r( F%ncid, varid_xa, rvalues, self%nobs, &
                                recvcounts, rdispls, &
                                me, MasterPE, MPI_COMM_CALC, status )

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
      
    end if
    
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

    character(len=*), parameter  ::  rname = mname//'/ObsOpers_Write_r'

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


  subroutine read_obs( domain, maxobs, &
                        stnid, flat,flon,falt, y,stddev, scode, &
                        ipar, nobs, status )
  !-----------------------------------------------------------------------
  ! @description
  ! Read in observations y in standard format
  ! (Station-No flat flon falt obs-value stddev).
  ! @author M.Kahnert
  !-----------------------------------------------------------------------
  
    use GridValues_ml, only : coord_in_domain
    use Io_ml        , only : IO_TMP
    use AOD_PM_ml    , only : AOD_init,wavelength
    use MPI          , only : MPI_IN_PLACE,MPI_INTEGER,MPI_SUM
    use MPI_Groups_ml, only : MasterPE,MPI_COMM_CALC
  
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
    integer, intent(out)            ::  ipar(maxobs)
    integer, intent(out)            ::  nobs
    integer, intent(out)            ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/read_obs'
    
    ! --- local ----------------------------------

    integer             ::  nd
    character(len=256)  ::  file
    logical             ::  exist
    integer             ::  no
    character(len=1024) ::  line
    real                ::  flon0, flat0, falt0, y0, stddev0
    character(len=16)   ::  scode0
    integer             ::  slen
                                            
    ! --- begin ----------------------------------

    ! init counter:
    nobs = 0
    ! init index range:
    obsData%iobs(0)=0
    obsData%iobs(1)=0

    !-----------------------------------------------------------------------
    ! open observational data file in standard format, if applicable
    !-----------------------------------------------------------------------
    ! loop over data sets:
    do nd = 1, nobsData
    
      ! undefined ? probably a mistake ..
      if ( len_trim(obsData(nd)%file) == 0 ) then
        ! info ...
        write (gol,'("no data file specified for Obsdata ",i6)') nd; call goErr
        TRACEBACK; status=1; return
      end if

      ! replace time values in templage:
      file = date2string( obsData(nd)%file, current_date )
      
      ! might not be present ...
      inquire( file=trim(file), exist=exist )
      if ( .not. exist ) then
        write (gol,'("WARNING - observation file `",a,"` not found; continue ...")') trim(file); call goPr
        cycle
      end if
      
      ! open:
      open( IO_TMP, file=trim(file), form='formatted', status='old', iostat=status )
      if ( status /= 0 ) then
        write (gol,'("could not open obsdata file : ",a)') trim(file); call goErr
        TRACEBACK; status=1; return
      end if

      ! info ...
      write (gol,'("Obsdata opened ",a)') trim(file); call goPr

      !-----------------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------------

      ! init index to collected observation array:
      obsData(nd)%iobs(0) = nobs+1

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
          read(line,*,iostat=status)no,flat0,flon0,falt0,y0,stddev0!,scode0
          IF_NOT_OK_RETURN(status=1)         

          ! read station code, code starts with '#', skip leading hash:
          no=index(line,'#')
          slen=len_trim(line)
          if(no>0) scode0=line(no+1:slen)

          ! if ouside domain/scope, next record:
          if(.not.coord_in_domain(domain,flon0,flat0))cycle

          ! check if in range:
          if ( (y0 < obsData(nd)%min) .or. &
               (y0 > obsData(nd)%max)      ) then
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
          stddev(nobs) = stddev0
          ! copy station code, skip leading '#'
          scode (nobs) = trim(scode0)
          ! store data set number:
          ipar(nobs) = nd
          ! truncate std.dev. to relative or absolute maximum:
          stddev(nobs) = max( stddev(nobs), &
                              obsData(nd)%error_rel * y(nobs), &
                              obsData(nd)%error_rep            )
          ! check .. 
          if ( stddev(nobs) <= 0e0 ) then
            print dafmt,'WARNING obs stddev <= 0'
            stddev(nobs) = 1e-9
          end if
          
          !! testing ...
          !write (gol,'(a,": observation ",2i6,2f8.2)') rname, nd, nobs, flon(nobs), flat(nobs); call goPr

        else
          write (gol,'("while reading observations")'); call goErr
          TRACEBACK; status=1; return
        end if

      end do  ! lines

      ! store end number:
      obsData(nd)%iobs(1) = nobs

      ! close:
      close( IO_TMP, iostat=status )
      IF_NOT_OK_RETURN(status=1)

    end do  ! data sets
!-----------------------------------------------------------------------
! Setup obsData &  Write #obs to log file
!-----------------------------------------------------------------------
    do nd=1,nobsData
! fill obsData(nd) structure
      if(.not.obsData(nd)%set)then
        no=find_index(obsData(nd)%name,obsVarName(:))
        obsData(nd)%found=(no>0)
        if(obsData(nd)%found)then
          obsData(nd)%ichem=no                    ! index observed
          obsData(nd)%ichemObs=ichemobs(no)       ! index observed/unobserved
          obsData(nd)%ispec=varSpec(ichemobs(no)) ! index species
          call Units_Scale(obsData(nd)%unit,obsData(nd)%ispec,obsData(nd)%unitconv,&
                           needroa=obsData(nd)%unitroa)
        else
          select case(obsData(nd)%name)
          case("BSC","EXT","AOD")
            call AOD_init("DA_Obs:"//trim(obsData(nd)%name),wlen=obsData(nd)%unit)
            obsData(nd)%ichem=find_index(obsData(nd)%unit,wavelength)
          case default
            obsData(nd)%ichem=find_index(obsData(nd)%name,chemgroups(:)%name)
            call Units_Scale(obsData(nd)%unit,-1,obsData(nd)%unitconv,&
                             needroa=obsData(nd)%unitroa)
          endselect
        endif
        call CheckStop(obsData(nd)%deriv=='',"Inconsistent obsData%deriv")
        write(obsData(nd)%tag,"(2(A,1X),'(',A,')')")&
          trim(obsData(nd)%name),trim(obsData(nd)%deriv),trim(obsData(nd)%unit)
        call CheckStop(obsData(nd)%ichem<1,&
          "Unsupported obsData "//trim(obsData(nd)%tag))
        obsData(nd)%set=.true.
      endif
! check #obs
      no=0
      if(any(obsData(nd)%iobs/=0))&
        no=DIM(obsData(nd)%iobs(1)+1,obsData(nd)%iobs(0))
      call CheckStop(no,count(ipar(:nobs)==nd),"Inconsistent obsData%nobs")
! total #obs
      if(domain=='processor')then
        CALL MPI_REDUCE(MPI_IN_PLACE,no,1,MPI_INTEGER,MPI_SUM,&
          MasterPE,MPI_COMM_CALC,status)
        IF_MPI_NOT_OK_RETURN(status=1)
      endif
! log #obs
      if(MasterProc)then
        file=date2string(obsData(nd)%file,current_date)
        write(damsg,"('obsData(',I0,') contains ',I0,2(1X,A),' observations')")&
          nd,no,trim(obsData(nd)%name),trim(obsData(nd)%deriv)
        write(damsg,dafmt)trim(damsg)
        call PrintLog(damsg)
        call PrintLog(file)
      endif
    enddo
    
    ! ok
    status = 0

  end subroutine read_obs
end module DA_Obs_ml
