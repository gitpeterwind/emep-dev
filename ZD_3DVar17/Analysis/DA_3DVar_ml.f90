!*****************************************************************************!
!
! EMEP 3D-var data assimilation module
!
! HISTORY
!   <2014 Origianal code by M. Kahnert, adapated for EMEP model by Alvaro Valdebenito
!   2014-2016 Re-implementation by Arjo Segers
!
!*****************************************************************************!
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!*****************************************************************************!

module DA_3DVar_ml

#ifdef with_ajs
  use GO, only : gol, goPr, goErr
  use GO, only : GO_Print_Set
  !use GO, only : goMem
#else
  use DA_Util_ml     , only : gol, goPr, goErr
#endif
  use TimeDate_ml    , only : date
  use EMEP_BCovarSqrt, only : T_BCovarSqrt
  use DA_Util_ml     , only : T_Domains


  implicit none

  ! --- in/out -----------------------------------------

  private

  public  ::  NTIMING_3DVAR
  public  ::  T_3DVAR

  public  ::  DA_3DVar_Init, DA_3DVar_Done
  public  ::  Main_3DVar


  ! --- const -----------------------------------------

  character(len=*), parameter  ::  mname = 'DA_3DVar_ml'

  ! timing parameters:
  integer, parameter  ::  NTIMING_3DVAR = 0
  integer, parameter  ::  T_3DVAR       = 0
  

  ! --- local ----------------------------------------

  ! Analysis date patterns:
  integer, parameter :: ANALYSIS_NDATE_MAX = 24   ! MAX patterns to check
  integer, save      :: ANALYSIS_NDATE = 1        ! Patterns to check
  type(date), dimension(ANALYSIS_NDATE_MAX) :: &  ! When to perform the analysis,
    analysis_date=date(-1,-1,-1,-1,00)            ! every hour by default

  ! link from observation element to original observation dataset:
  integer, allocatable  ::  iObsData(:)

  ! Limit dx,du to 500%
  real, parameter       ::  ANALYSIS_RELINC_MAX = 5.0

  ! Solver (m1qn3) settings
  integer, save :: &
    solver_scaling =  1,& ! 0=DIS, 1=SIS
    solver_nupdates= 20,& ! number of updates stored in work array
    solver_maxiter =500,& ! maximum number of solver iterations
    solver_maxsim  =  2   ! maximum number of simulations per solver iteration
  logical, save ::  &
    solver_logall = .false. ! save log from all processors
  ! summary table:
  integer               ::  solver_table
  ! counter:
  integer               ::  solver_number

#ifdef with_ajs
  ! timers:
  integer               ::  itim_read_obs, itim_innov, itim_swap
  integer               ::  itim_loop, itim_optimizer, itim_costfunc
  integer               ::  itim_chi2x, itim_innov_adj, itim_gradient, itim_x2chi
  integer               ::  itim_post
#endif

  ! lengths:
  integer, parameter    ::  TRACER_NAME_LEN = 32
  integer, parameter   ::  FILE_NAME_LEN   = 1024

  ! info on assim with single covar:
  type T_BInfo
    ! short description
    character(len=FILE_NAME_LEN)      ::  name
    ! level selection:
    integer                           ::  nlev_all
    integer                           ::  nlev
    integer, allocatable              ::  levels(:)  ! (nlev)
    !integer                           ::  nlev_sfc
    !integer, allocatable              ::  levels_sfc(:)  ! (nlev_sfc)
    ! number of tracers involved:
    integer                           ::  ntracer
    ! tracer names:
    character(len=TRACER_NAME_LEN), allocatable ::  tracers(:)  ! (ntracer)
    ! mapping from tracer to iObsComp:
    integer, allocatable              ::  iObsComp(:)  ! (ntracer)
    ! mapping from iObsComp to itracer, could be undefined:
    integer, allocatable              ::  itracer(:)  ! (nObsComp)
    ! scale factors for sigma:
    logical                           ::  sigma_scale_online
    real, allocatable                 ::  sigma_scale_factors(:)   ! (ntracer)
    ! input file
    character(len=FILE_NAME_LEN)      ::  filename
    ! covariance square root:
    type(T_BCovarSqrt)                ::  BCovarSqrt
  end type T_BInfo

  ! collection of B matrices and associated tracers:
  integer                             ::  nBmat
  type(T_BInfo), allocatable          ::  Bmats(:)  ! (nBmat)

  ! local domain definitions:
  type(T_Domains)       ::  doms_adv_m   ! (s,x,y,z) all tracers, model xy decomposition
  type(T_Domains)       ::  doms_an_m    ! (x,y,z,s) analysis tracers, model xy decomposition
  type(T_Domains)       ::  doms_an_fg   ! (x,y,z,s) analysis tracers, fft y decomposition
  
  ! testing ..
  integer ::  nmemtest = 0


contains


  !-----------------------------------------------------------------------
  ! module init/done
  !-----------------------------------------------------------------------

  subroutine DA_3DVar_Init( status )

    use DA_Util_ml       , only : DA_Util_Init
#ifdef with_ajs
    use DA_Util_ml       , only : GO_Timer_Def
    use DA_Util_ml       , only : GO_Print_Set
#endif
    use DA_Util_ml       , only : goGetFU
    use ModelConstants_ml, only : masterProc
    use ModelConstants_ml, only : KMAX_MID
    use MPI_Groups_ml    , only : MPI_COMM_CALC
    use Io_Nums_ml       , only : IO_NML
    use DA_ml            , only : nlev
    use DA_Obs_ml        , only : DA_Obs_Init

    ! --- in/out ---------------------------------

    integer, intent(out)           ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/DA_3DVar_Init'

    ! --- namelists ------------------------------

    namelist /DA_CONFIG/ &
      analysis_ndate, analysis_date, &
      solver_scaling,solver_nupdates,solver_maxiter,solver_maxsim,solver_logall

    ! --- local ----------------------------------

    character(len=1024)     ::  fname

    ! --- begin ----------------------------------

    ! check ...
#ifndef _MPI
    write (gol,'("MPI code should be enabled, define _MPI macro!")'); call goErr
    TRACEBACK; status=1; return
#endif

    ! initialize utilities:
    call DA_Util_Init( status )
    IF_NOT_OK_RETURN(status=1)

#ifdef with_ajs
    ! info as error message to have it in the log file ..
    write (gol,'(a,": enable GO logging from this routine if necessary ...")') rname; call goPr
    !! log from all:
    !call GO_Print_Set( status, apply=.true. )
    !IF_NOT_OK_RETURN(status=1)
    !! info ..
    !write (gol,'(a,": enabled GO logging on all processes ...")') rname; call goPr

    ! define timers:
    call GO_Timer_Def( itim_read_obs , 'read observations', status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_innov    , 'innovations'      , status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_swap     , 'swap'             , status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_loop     , 'optimization loop', status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_optimizer, 'optimizer'        , status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_costfunc , 'costfunction'     , status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_chi2x    , 'chi to x'         , status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_x2chi    , 'x to chi'         , status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_innov_adj, 'innovation adj'   , status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_gradient , 'gradient'         , status )
    IF_NOT_OK_RETURN(status=1)
    call GO_Timer_Def( itim_post     , 'post'             , status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! Read settings for:
    ! - analysis date
    ! - solver
    ! back to start:
    rewind(IO_NML)
    ! read using namelist:
    read (unit=IO_NML,nml=DA_CONFIG,iostat=status)
    if ( status /= 0 ) then
      write (gol,'(a,": reading namelist DA_CONFIG from: config_emep.nml")') rname; call goErr
      TRACEBACK; status=1; return
    end if

        !! testing ...
        !solver_maxiter = 5
        !write (gol,*) 'qqq reset maxiter to ', solver_maxiter; call goPr

    ! table with norms and costs will be written by root:
    if ( masterProc ) then
      ! target file:
      fname = 'm1qn3.csv'
      ! select free file unit:
      call goGetFU( solver_table, status )
      IF_NOT_OK_RETURN(status=1)
      ! open:
      open( unit=solver_table, file=trim(fname), form='formatted', iostat=status )
      if ( status /= 0 ) then
        write (gol,'("could not open solver table file (directory not present?) : ",a)') trim(fname); call goErr
        TRACEBACK; status=1; return
      end if
      ! header:
      write (solver_table,'("number,step,J,Jb,Jo,gn")')
    end if

    ! number of levels:
    nlev = KMAX_MID  ! 20

    ! init observations:
    call DA_Obs_Init( status )
    IF_NOT_OK_RETURN(status=1)


    ! ok
    status = 0

  end subroutine DA_3DVar_Init


  ! ***


  subroutine DA_3DVar_Done( status )

    use DA_Util_ml       , only : DA_Util_Done
    use ModelConstants_ml, only : MasterProc
    use DA_Obs_ml        , only : DA_Obs_Done

    ! --- in/out ----------------------------

    integer, intent(out)           ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/DA_3DVar_Done'

    ! --- local ----------------------------

    integer     ::  iB

    ! --- begin -----------------------------

    ! done with observations:
    call DA_Obs_Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over covar matrices:
    do iB = 1, nBmat
      ! clear:
      deallocate( Bmats(iB)%levels, stat=status )
      IF_NOT_OK_RETURN(status=1)
      !deallocate( Bmats(iB)%levels_sfc, stat=status )
      !IF_NOT_OK_RETURN(status=1)
      deallocate( Bmats(iB)%tracers, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( Bmats(iB)%iObsComp, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( Bmats(iB)%itracer, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( Bmats(iB)%sigma_scale_factors, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! done with covariance structure:
      call Bmats(iB)%BCovarSqrt%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! clear:
    deallocate( Bmats, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! only on root ...
    if ( masterProc ) then
      ! close data file:
      close( unit=solver_table, iostat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! done with domains on model decomposition:
    call doms_adv_m%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! done with domains for increment swap:
    call doms_an_m%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call doms_an_fg%Done( status )
    IF_NOT_OK_RETURN(status=1)
    !call doms_an_fs%Done( status )
    !IF_NOT_OK_RETURN(status=1)

    ! done with utilities:
    call DA_Util_Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine DA_3DVar_Done


  ! ***


  !-----------------------------------------------------------------------
  subroutine Main_3DVar( status )
  !-----------------------------------------------------------------------
  ! @description
  ! Main module for starting 3dvar program.
  ! @author M.Kahnert & AMVB
  !-----------------------------------------------------------------------

    use ModelConstants_ml    , only : ANALYSIS
    use ModelConstants_ml    , only : MasterProc
    use TimeDate_ml          , only : current_date
    use TimeDate_ExtraUtil_ml, only : date2string,compare_date
    use DA_ml                , only : dafmt     => da_fmt_msg
    use DA_ml                , only : DAFMT_DEF => DA_FMT_DEF
    use DA_ml                , only : debug => DEBUG_DA

    ! --- in/out ----------------------------------

    integer, intent(out)    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter  ::  rname = mname//'/main_3dvar'

    ! --- saved --------------------------------------

    logical, save :: first_call=.true.

    ! --- local --------------------------------------

    ! --- begin --------------------------------------

    ! not enabled ? then leave:
    if ( .not. ANALYSIS ) return
    
    ! reset format for messages:
    dafmt = date2string(DAFMT_DEF,current_date)

    ! need to initialize ?
    if ( first_call ) then
      ! info ...
      if ( debug .and. MasterProc) print dafmt,'Initialisation'
      ! initialize:
      call Init_3DVar( status )
      IF_NOT_OK_RETURN(status=1)
      ! reset flag:
      first_call=.false.
    end if

    ! leave if no analysis time:
    if(.not.compare_date(ANALYSIS_NDATE,current_date,analysis_date,wildcard=-1))then
      status=0; return
    endif

    ! info ...
    if(debug.and.MasterProc) print dafmt,'Start analysis'
    ! run:
    call generic3dvar( current_date, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine main_3dvar


  ! ***

  !-----------------------------------------------------------------------
  subroutine Init_3DVar( status )
  !-----------------------------------------------------------------------
  ! @description
  ! Initialise 3dvar variables.
  ! @author AMVB
  !-----------------------------------------------------------------------

    use MPIF90               , only : MPIF90_BCast
    use MPI_Groups_ml        , only : MPI_COMM_CALC
    use ModelConstants_ml    , only : MasterProc, nproc
    use ModelConstants_ml    , only : RUNDOMAIN
    use ModelConstants_ml    , only : KMAX_MID
    use ModelConstants_ml    , only : KMAX_BND, KCHEMTOP
    use MPI_Groups_ml        , only : MasterPE
    use Par_ml               , only : me
    use Par_ml               , only : tlimax, tgi0, tgi1, tljmax, tgj0, tgj1
    use ChemGroups_ml        , only : chemgroups       ! group  names
    use ChemSpecs_adv_ml     , only : NSPEC_ADV
    use ChemSpecs_shl_ml     , only : NSPEC_SHL        ! Maps indices
    use SmallUtils_ml        , only : find_index
    use Io_Nums_ml           , only : IO_NML
    use DA_ml                , only : dafmt => da_fmt_msg
    use DA_ml                , only : debug => DEBUG_DA
    use DA_Obs_ml            , only : nObsComp, ObsCompInfo
    use DA_ml                , only : nlev

    ! --- in/out ----------------------------------

    integer, intent(out)    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Init_3DVar'

    ! --- local ----------------------------------

    integer               :: nvar
    integer               :: k
    integer               :: ny_local, iy_offset

    integer               ::  nx, ny
    integer               ::  itracer

    integer                                     ::  iObsComp
    character(len=TRACER_NAME_LEN), allocatable ::  Bname(:)       ! (nChemObs)
    logical, allocatable                        ::  Bml(:)         ! (nChemObs)
    integer, allocatable                        ::  tracer2B(:)    ! (nChemObs)
    real, allocatable                           ::  Bsigfac(:)     ! (nChemObs)
    integer                                     ::  iB

    ! --- begin ------------------------------

    ! info ...
    write (gol,'(a,": initalize 3D var variables ...")') rname; call goPr

    ! assumed grid shape:
    nx   = RUNDOMAIN(2) - RUNDOMAIN(1) + 1
    ny   = RUNDOMAIN(4) - RUNDOMAIN(3) + 1
    nlev = kmax_mid

    !-----------------------------------------------------------------------
    ! Read covariance matrix
    !-----------------------------------------------------------------------

    ! * assign B to tracer(s)

    ! info ...
    write (gol,'(a,": setup B matrices ..")') rname; call goPr

    ! init counter:
    nBmat = 0
    ! storage for keywords defining B for individual or groups of tracers:
    allocate( Bname(nObsComp), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bml(nObsComp), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( tracer2B(nObsComp), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bsigfac(nObsComp), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! defaults:
    Bname(:)    = ''
    Bml(:)      = .true.
    tracer2B(:) = ''
    Bsigfac(:)  = 1.0

    ! loop over tracers to be analysed:
    do iObsComp = 1, nObsComp
      ! info ...
      write (gol,'(a,":   observed tracer ",i0," : ",a)') rname, iObsComp, trim(ObsCompInfo(iObsComp)%name); call goPr
      write (gol,'(a,":     associated B : ",a)') rname, trim(ObsCompInfo(iObsComp)%Bfile); call goPr

      ! search in list of already defined keys:
      iB = find_index( trim(ObsCompInfo(iObsComp)%Bfile), Bname(:) )
      ! new?
      if ( iB < 0 ) then
        ! increase counter
        nBmat = nBmat + 1
        ! store name:
        Bname(nBmat) = trim(ObsCompInfo(iObsComp)%Bfile)
        ! copy factor:
        Bsigfac(nBmat) = ObsCompInfo(iObsComp)%Bsigfac
        ! use B with model levels?
        select case ( trim(ObsCompInfo(iObsComp)%deriv) )
          !~ 2D operators:
          case ( '2D-ML-SFC', '2D-OBS-SFC' )
            Bml(nBmat) = .false.
          !~ 3D operators:
          case ( '3D-ML', '3D-ML-SFC', '3D-ML-TC' )
            Bml(nBmat) = .true.
          !~
          case default
            write (gol,'("unsupported deriv: ",a)') trim(ObsCompInfo(iObsComp)%deriv); call goErr
            TRACEBACK; status=1; return
        end select
        ! assign index:
        iB = nBmat
        ! info ...
        write (gol,'(a,":     assign new B matrix nr. ",i0)') rname, iB; call goPr
      end if
      ! store index in mapping from tracer to covariance:
      tracer2B(iObsComp) = iB
    end do ! tracers

    ! info ..
    write (gol,'(a,":   number of B matrices: ",i0)') rname, nBmat; call goPr

    ! storage for B matrices and associated info:
    allocate( Bmats(nBMat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over covar matrices:
    do iB = 1, nBmat

      ! store name:
      Bmats(iB)%name = trim(Bname(iB))

      ! all model levels (3D) ? otherwise 2D field
      if ( Bml(iB) ) then
        ! all model levels:
        Bmats(iB)%nlev_all = nlev
        Bmats(iB)%nlev     = nlev
        ! storage:
        allocate( Bmats(iB)%levels(nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill:
        do k = 1, nlev
          Bmats(iB)%levels(k) = k
        end do
      else
        ! only bottom layer ...
        Bmats(iB)%nlev_all = nlev
        Bmats(iB)%nlev     = 1
        ! storage:
        allocate( Bmats(iB)%levels(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! fill:
        Bmats(iB)%levels = (/nlev/)  ! bottom layer
      end if
      
      ! number of tracers:
      Bmats(iB)%ntracer = count( tracer2B == iB )
      ! storage for tracer names in B:
      allocate( Bmats(iB)%tracers(Bmats(iB)%ntracer), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! storage for online scale factors:
      allocate( Bmats(iB)%sigma_scale_factors(Bmats(iB)%ntracer), stat=status )
      IF_NOT_OK_RETURN(status=1)
      !~ init with no scaling:
      Bmats(iB)%sigma_scale_factors = 1.0
      ! storage for mapping from B tracers to observed components:
      allocate( Bmats(iB)%iObsComp(Bmats(iB)%ntracer), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! mapping from nChemObs to tracer:
      allocate( Bmats(iB)%itracer(nObsComp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! init as undefined:
      Bmats(iB)%itracer = -999
      ! init counter:
      itracer = 0
      ! loop over tracers:
      do iObsComp = 1, nObsComp
        ! current?
        if ( tracer2B(iObsComp) == iB ) then
          ! increase counter:
          itracer = itracer + 1
          ! store name:
          Bmats(iB)%tracers(itracer) = trim(ObsCompInfo(iObsComp)%name)
          ! store index:
          Bmats(iB)%iObsComp(itracer) = iObsComp
          ! store reverse index:
          Bmats(iB)%itracer(iObsComp) = itracer
          ! copy covariance filename:
          Bmats(iB)%filename = trim(ObsCompInfo(iObsComp)%Bfile)
        end if ! current
      end do  ! observed variables

      ! info ...
      write (gol,'(a,":     B matrix ",i0," (",a,")")') rname, iB, trim(Bname(iB)); call goPr
      do itracer = 1, Bmats(iB)%ntracer
        write (gol,'(a,":       tracer ",i2," ",a)') rname, itracer, trim(Bmats(iB)%tracers(itracer)); call goPr
      end do ! tracers
      write (gol,'(a,":       file: ",a)') rname, trim(Bmats(iB)%filename); call goPr

      ! init covariance structure,
      ! setup mapping to the subset of tracers
      ! and subset of observed levels:
      call Bmats(iB)%BCovarSqrt%Init( trim(Bmats(iB)%filename), status, &
                                       tracers=Bmats(iB)%tracers, &
                                       obs_levels=Bmats(iB)%levels, &
                                       comm=MPI_COMM_CALC )
      IF_NOT_OK_RETURN(status=1)

      ! info ...
      write (gol,'(a,":   adhoc scaling factor for sigma : ",f6.1)') &
                        rname, Bsigfac(iB); call goPr
      ! if negative, online scaling based on obs-sim statistics:
      Bmats(ib)%sigma_scale_online = Bsigfac(iB) < 0.0
      ! switch ...
      if ( Bmats(ib)%sigma_scale_online ) then
        ! info ...
        write (gol,'(a,":     apply online scaling per hour ...")') rname; call goPr
      else if ( Bsigfac(iB) == 1.0 ) then
        ! info ...
        write (gol,'(a,":     no scaling specified ...")') rname; call goPr
      else
        ! check ...
        if ( Bmats(iB)%ntracer > 1 ) then
          write (gol,'("no sigma scale factor yet for covariance with multiple tracers")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! write ...
        write (gol,'(a,":     apply factor ...")') rname; call goPr
        ! apply scale factor:
        Bmats(iB)%BCovarSqrt%S = Bmats(iB)%BCovarSqrt%S * Bsigfac(iB)
      end if

    end do ! covars

    ! clear:
    deallocate( Bname, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( Bml, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( tracer2B, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! *

    ! loop over B matrices:
    do iB = 1, nBmat

      ! check ...
      call Bmats(iB)%BCovarSqrt%Check( status, nlon=nx, nlat=ny, nlev=nlev )
      IF_NOT_OK_RETURN(status=1)

      ! extract info on grid and decomposition used by covariance:
      call Bmats(iB)%BCovarSqrt%Get( status, nlat_ex_local=ny_local, ilat_ex_offset=iy_offset )
      IF_NOT_OK_RETURN(status=1)

    end do ! B matrices

    !-----------------------------------------------------------------------
    ! model domain decomposition
    !-----------------------------------------------------------------------

    ! info ...
    write (gol,'(a,": model decomposition")') rname; call goPr
    do k = 0, nproc-1
      write (gol,'(a,":   domain ",i2," cells  [",i3,",",i3,"] x [",i3,",",i3,"]")') &
               rname, k, tgi0(k), tgi1(k), tgj0(k), tgj1(k); call goPr
    end do

    ! info ...
    if ( ny_local > 0 ) then
      write (gol,'(a,":   ufft local y-slab : ",i0,":",i0)') rname, iy_offset+1, iy_offset+ny_local; call goPr
    else
      write (gol,'(a,":   ufft without local y-slab; ny_local = ",i0)') rname, ny_local; call goPr
    end if

    ! original model decomposition on standard domain:
    !                          var        lon       lat       lev
    call doms_adv_m%Init( (/         1, tgi0(me), tgj0(me),     1    /), &
                          (/ nspec_adv, tgi1(me), tgj1(me), kmax_mid /), status )
    IF_NOT_OK_RETURN(status=1)

    ! analyis tracers and order;
    ! original model decomposition on standard domain:
    !                         lon       lat      lev
    call doms_an_m%Init( (/ tgi0(me), tgj0(me),    1,        1 /), &
                         (/ tgi1(me), tgj1(me), nlev, nObsComp /), status )
    IF_NOT_OK_RETURN(status=1)

    ! analyis tracers and order;
    ! define new local decomposition on standard domain ;
    ! some domains are completely in the extension and therefore undefined:
    if ( (iy_offset >= ny) .or. (ny_local == 0) ) then
      ! undefined, use dummy non-existing range with correct zero size:
      !                       var lon lat lev
      call doms_an_fg%Init( (/  0,  0,  0,  0 /), &
                            (/ -1, -1, -1, -1 /), status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! define slab, truncate at y bound:
      !                       lon           lat                lev      var
      call doms_an_fg%Init( (/  1,     iy_offset+       1    ,  1  ,        1 /), &
                            (/ nx, min(iy_offset+ny_local,ny), nlev, nObsComp /), status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! info ...
    write (gol,'(a,": end")') rname; call goPr

    ! ok
    status = 0

  end subroutine init_3dvar


  ! ***


  !-----------------------------------------------------------------------
  ! @description
  ! Generic version of 3d variational analysis.
  ! @author M.Kahnert
  !-----------------------------------------------------------------------

  subroutine generic3dvar( cdate, status )

    use MPIF90               , only : MPIF90_AllReduce, MPI_SUM
#ifdef with_ajs
    use DA_Util_ml           , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
#endif
    use MPI_Groups_ml        , only : MPI_COMM_CALC
    use ModelConstants_ml    , only : MasterProc, NPROC
    use ModelConstants_ml    , only : RUNDOMAIN
    use ModelConstants_ml    , only : KMAX_MID
    use MPI_Groups_ml        , only : MasterPE
    use My_Timing_ml         , only : Code_timer, Add_2timing
    use TimeDate_ml          , only : date  ! date/time structure
    use Par_ml               , only : me
    use Par_ml               , only : MAXLIMAX, MAXLJMAX   ! local x, y dimensions
    use Par_ml               , only : tlimax, tgi0, tgi1, tljmax, tgj0, tgj1
    use Par_ml               , only : limax, ljmax
    use GridValues_ml        , only : glon, glat
    use ChemSpecs_adv_ml     , only : NSPEC_ADV
    use MetFields_ml         , only : z_bnd
    use ChemFields_ml        , only : xn_adv
    use Chemfields_ml        , only : PM25_water
    use Chemfields_ml        , only : PM25_water_rh50
    use DA_ml                , only : debug => DEBUG_DA
    use DA_ml                , only : dafmt => da_fmt_msg
    use DA_ml                , only : damsg => da_msg
    use DA_ml                , only : tim_before => datim_before, tim_after => datim_after
    use DA_Obs_ml            , only : T_ObsOpers
    use DA_Obs_ml            , only : Read_Obs
    use DA_Obs_ml            , only : nObsComp, ObsCompInfo
    use DA_ml                , only : nlev

    !-----------------------------------------------------------------------
    ! Formal parameters
    !-----------------------------------------------------------------------

    type(date), intent(in)    ::  cdate   ! current date
    integer, intent(out)      ::  status

    !-----------------------------------------------------------------------
    ! constants
    !-----------------------------------------------------------------------

    character(len=*), parameter  ::  rname = mname//'/generic3dvar'

    !-----------------------------------------------------------------------
    ! Local parameters
    !-----------------------------------------------------------------------

    integer               ::  iout
    integer, pointer      ::  out_group(:)
    integer, target       ::  out_group1(1)
    integer               ::  ispec
    logical               ::  needroa

    integer               ::  i,j,ilev,ilev0,k,kk,ii, maxobs !, nvar
    integer               ::  nobs
    integer               ::  nobs_tot
    integer               ::  iObsComp
    
    ! collected observations:
    integer, allocatable            ::  stnid(:)        ! (nobs)
    real, allocatable               ::  flat(:)         ! (nobs)
    real, allocatable               ::  flon(:)         ! (nobs)
    real, allocatable               ::  falt(:)         ! (nobs)
    real, allocatable               ::  obs(:)          ! (nobs)
    real, allocatable               ::  obsstddev1(:)   ! (nobs)
    character(len=32), allocatable  ::  stncodes(:)     ! (nobs)
    logical, allocatable            ::  obsanalyse(:)   ! (nobs)

    integer               ::  lnx, lny
    logical               ::  has_local_domain

    ! idem operator on model and assimilation decomposition:
    type(T_ObsOpers)      ::  Hops_m
    type(T_ObsOpers)      ::  Hops_f
    type(T_ObsOpers)      ::  Hops_f_B
    integer               ::  nobs_B
    integer               ::  nana, nana_tot

    real, allocatable     ::  sf_an (:,:  ,:)     ! (limax,ljmax     ,nObsComp)
    real, allocatable     ::  xn_an (:,:,:,:)     ! (limax,ljmax,nlev,nObsComp)
    real, allocatable     ::  dx_an (:,:,:,:)     ! (limax,ljmax,nlev,nObsComp)
    real, allocatable     ::  dx_loc(:,:,:,:)     ! (lnx  ,lny  ,nlev,nObsComp)
    character(len=16), allocatable    ::  xn_obs_units(:)   ! (nObsComp)
    character(len=16), allocatable    ::  xn_adv_units(:)   ! (nspec_adv)
    character(len=16)                 ::  tracer_name
    character(len=16)                 ::  tracer_units

    real, allocatable     ::  sigma_m(:,:,:,:)     ! (limax,ljmax,nlev,nObsComp)

    integer               ::  iB
    real, allocatable     ::  dx_loc_B    (:,:,:,:)    ! (lnx  ,lny  ,B_nlev,ntracer)
    real, allocatable     ::  du_loc_B    (:,:,:,:)    ! (lnx  ,lny  ,nlev  ,ntracer)
    integer               ::  itracer
    integer               ::  iobs
    integer               ::  itime

    !-----------------------------------------------------------------------
    ! begin
    !-----------------------------------------------------------------------

#ifdef with_ajs
    !! testing ..
    !call goMem( rname//' - MEMORY begin', status )
    !IF_NOT_OK_RETURN(status=1)
#endif

    !-----------------------------------------------------------------------
    ! first-guess output
    !-----------------------------------------------------------------------
    
    ! fill arrays for requested subclass using concentrations:
    call Fill_Output_xn_adv( '3DVAR_FG', xn_adv, status )
    IF_NOT_OK_RETURN(status=1)
    
#ifdef with_ajs
    !! testing ..
    !call goMem( rname//' - MEMORY 1    ', status )
    !IF_NOT_OK_RETURN(status=1)
#endif

    !-----------------------------------------------------------------------
    ! local grid
    !-----------------------------------------------------------------------

    ! extract local grid size for convenience:
    lnx = doms_an_fg%shp(1,me)
    lny = doms_an_fg%shp(2,me)
    ! flag ...
    has_local_domain = all( (/lnx,lny/) > 0 )
    ! part of model domain in this part of spectral grid ?
    if ( has_local_domain ) then
      ! info ..
      write (gol,'(a,": local 3D-var domain, size : ",2i6)') rname, lnx, lny; call goPr
    else
      ! use dummy values for propper allocation:
      lnx = 1
      lny = 1
      ! info ...
      write (gol,'(a,": no local 3D-var domain")') rname; call goPr
    end if

#ifdef with_ajs
    !! testing ..
    !call goMem( rname//' - MEMORY 2a   ', status )
    !IF_NOT_OK_RETURN(status=1)
#endif

    !-----------------------------------------------------------------------
    ! read observations
    !-----------------------------------------------------------------------

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_read_obs, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! obs regridded to model resolution ?
    !maxobs = nx*ny
    maxobs = (RUNDOMAIN(2)-RUNDOMAIN(1)+1) * (RUNDOMAIN(4)-RUNDOMAIN(3)+1)

    ! info ..
    if (debug) print dafmt,'Read observations'

    ! storage per observation:
    !~ index in 'obsData' array:
    allocate( iObsData(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    !~ other:
    allocate( stnid(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( stncodes(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( flat(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( flon(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( falt(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( obs(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( obsstddev1(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( obsanalyse(maxobs), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! selection on local grid point domain used by model:
    call read_obs( "processor", maxobs, &
                     stnid, flat,flon,falt, obs, obsstddev1, stncodes, obsanalyse, &
                     iObsData, nobs, status )
    IF_NOT_OK_RETURN(status=1)

    ! number of analysed observations:
    if ( nobs > 0 ) then
      ! count:
      nana = count(obsanalyse(1:nobs))
    else
      ! dummy:
      nana = 0
    end if
    
    ! info ...
    write (gol,'(a,": local number of observations read : ",i6," (analysed ",i0,", validation ",i0,")")') &
                  rname, nobs, nana, nobs-nana; call goPr
    
    ! total number of observations over all domains:
    call MPIF90_AllReduce( nobs, nobs_tot, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    call MPIF90_AllReduce( nana, nana_tot, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'(a,": total number of observations read : ",i6," (analysed ",i0,", validation ",i0,")")') &
                   rname, nobs_tot, nana_tot, nobs_tot-nana_tot; call goPr

#ifdef with_ajs
    ! end timing:
    call GO_Timer_End( itim_read_obs, status )
    IF_NOT_OK_RETURN(status=1)

    !! testing ..
    !call goMem( rname//' - MEMORY 3a   ', status )
    !IF_NOT_OK_RETURN(status=1)
    !write (gol,'("WARNING - reset n_obs_tot=0 ...")'); call goErr
    !nobs_tot = 0
#endif

    ! no observations at all ?
    if ( nobs_tot == 0 ) then

      ! info ...
      if(MasterProc) write (*,'("WARNING: No obserations found")')

    else

      ! info ..
      if(debug) write (*,'(1X,I0,1X,A)') nobs_tot, 'observations read.'

      !-----------------------------------------------------------------------
      ! extract background fields
      !-----------------------------------------------------------------------

      ! info ..
      write (gol,'(a,": collect tracer fields involved in analysis ...")') rname; call goPr
      write (gol,'(a,": limax,ljmax,nlev,nObsComp = ",6i4)') rname, limax,ljmax,nlev,nObsComp; call goPr

      ! NOTE: This is only needed for not changing 'H_op' too much,
      !       which extracts concentrations from an array with shape:
      !           (lon,lat,lev,ichemobs)
      !       In future, extract everthing directly from 'xn_adv',
      !       which is already used in 'H_op' for the indirect observations

      ! storage for extract of concentrations involved in analysis:
      allocate( sf_an(limax,ljmax     ,nObsComp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( xn_an(limax,ljmax,nlev,nObsComp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! storage for increments:
      allocate( dx_an(limax,ljmax,nlev,nObsComp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! storage for units:
      allocate( xn_obs_units(nObsComp), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! storage for units:
      allocate( xn_adv_units(size(xn_adv,1)), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! all the same ...
      do itracer = 1, size(xn_adv_units)
        xn_adv_units(itracer) = 'mol mol^-1'
      end do  ! tracers

      ! loop over variables involved in analysis:
      do iObsComp = 1, nObsComp
        ! info ...
        write (gol,'(a,": copy xn_adv species into xn_an var ",i0)') rname, iObsComp; call goPr
        !! extract slab:
        !!     x,y,z,s                     s     ,    x  ,   y   ,          z
        !xn_an(:,:,:,iObsComp) = xn_adv(varSpec(iObsComp),1:limax,1:ljmax,KMAX_MID-nlev+1:KMAX_MID)
        ! weighted sum:
        call ObsCompInfo(iObsComp)%FillFields( xn_adv, &
                    sf_an(:,:,iObsComp), xn_an(:,:,:,iObsComp), xn_obs_units(iObsComp), &
                    status )
        IF_NOT_OK_RETURN(status=1)
        ! store if needed:
        call Fill_Output_sf_xn( '3DVAR_OBS_FG', trim(ObsCompInfo(iObsComp)%name), &
                                   sf_an(:,:,iObsComp), xn_an(:,:,:,iObsComp), &
                                   trim(ObsCompInfo(iObsComp)%units), status )
        IF_NOT_OK_RETURN(status=1)
      end do

      !! info ...
      !write (gol,'(a,": xn_an range : ",2e16.6)') rname, minval(xn_an), maxval(xn_an); call goPr


      !-----------------------------------------------------------------------
      ! allocate work arrrays
      !-----------------------------------------------------------------------

      ! local y-slabs ?
      if ( lny > 0 ) then
        ! storage for local concentration increments, new decomposition:
        !allocate( xn_loc(lnx,lny,nlev,nObsComp), stat=status )
        !IF_NOT_OK_RETURN(status=1)
        allocate( dx_loc(lnx,lny,nlev,nObsComp), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! dummy:
        !allocate( xn_loc(1,1,1,1), stat=status )
        !IF_NOT_OK_RETURN(status=1)
        allocate( dx_loc(1,1,1,1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if

#ifdef with_ajs
      !! testing ..
      !call goMem( rname//' - MEMORY 4a   ', status )
      !IF_NOT_OK_RETURN(status=1)
#endif


      !-----------------------------------------------------------------------
      ! compute innovations H(xb)-y
      !-----------------------------------------------------------------------

      ! info ..
      write (gol,'(a,": compute innovations ...")') rname; call goPr

#ifdef with_ajs
      ! start timing:
      call GO_Timer_Start( itim_innov, status )
      IF_NOT_OK_RETURN(status=1)
#endif

      ! start timing:
      call Code_timer(tim_before)

      ! init operator:
      call Hops_m%Init( status )
      IF_NOT_OK_RETURN(status=1)

      ! fill innovations and Jacobian of observation operator ;
      ! use full concentration arrays on model decomposition:
      call get_innovations( Hops_m, nobs, &
                              xn_adv, xn_adv_units, &
                              sf_an, xn_an, xn_obs_units, &
                              maxobs, stnid, flat,flon,falt, obs, obsstddev1, stncodes, obsanalyse, &
                              'xf', status )
      IF_NOT_OK_RETURN(status=1)

      !!... TESTING ...
      !! save all:
      !call Hops_m%WriteToFile( cdate, status )
      !IF_NOT_OK_RETURN(status=1)
      !! break
      !write (gol,'("break after write obs")'); call goErr
      !TRACEBACK; status=1; return
      !!... TESTING ...

      ! assign unique id's to observations:
      call Hops_m%Set_IDs( status )
      IF_NOT_OK_RETURN(status=1)

      ! ...................................

      ! info ..
      write (gol,'(a,": evaluate sigma at observation locations ...")') rname; call goPr

      ! loop over B matrices:
      do iB = 1, nBmat

        ! get index of record holding current hour:
        call Bmats(iB)%BcovarSqrt%FindTime( cdate%hour, itime, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! loop over tracers in covar:
        do itracer = 1, Bmats(iB)%ntracer
          ! observed component index:
          iObsComp = Bmats(iB)%iObsComp(itracer)
          ! copy sigma, fill in dx_loc which is not used yet:
          if ( lny > 0 ) then
            dx_loc(:,:,:,iObsComp) = Bmats(iB)%Bcovarsqrt%S(:,:,:,itracer,itime)
          end if
        end do  ! tarcers in covar

      end do  ! Bmat
      
      ! storage for sigma at model grid:
      allocate( sigma_m(limax,ljmax,nlev,nObsComp), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! swap sigma field from analysis decomposition to model decomposition:
      call doms_an_fg%Swap( dx_loc, doms_an_m, sigma_m, status )
      IF_NOT_OK_RETURN(status=1)

      ! evaluate observation operator on sigma fields
      ! (lowest level for 2D surface field);
      ! store in 'sb' element as background-sigma:
      call Hops_m%Evaluate( 'sb', sigma_m(:,:,nlev,:), sigma_m, status )!, debug=.true. )
      IF_NOT_OK_RETURN(status=1)

      ! info ..
      write (gol,'(a,": obtain scale factors for sigma (might not be used) ...")') rname; call goPr
      
      ! determine scale factors for B.sigma based on Chi2 test,
      ! stored in Hops_m%BScale(1:nObsComp)
      call Hops_m%SetBScale( status )
      IF_NOT_OK_RETURN(status=1)
      
      ! info ..
      write (gol,'(a,":   optimal scale factors for B sigma:")') rname; call goPr
      do iObsComp = 1, nObsComp
        write (gol,'(a,":     comp ",i0," `",a,"` : ",f6.1)') &
            rname, iObsComp, trim(ObsCompInfo(iObsComp)%name), Hops_m%Bscale(iObsComp); call goPr
      end do

      ! store factors in B matrices:
      do iB = 1, nBmat
        ! info ...
        write (gol,'(a,":   B matrix ",i0)') rname, iB; call goPr
        ! loop over tracers in covar:
        do itracer = 1, Bmats(iB)%ntracer
          ! observed component index:
          iObsComp = Bmats(iB)%iObsComp(itracer)
          ! use online computed scale factors?
          if ( Bmats(iB)%sigma_scale_online ) then
            ! info ...
            write (gol,'(a,":     comp `",a,"` online scaled with factor ",f6.1)') &
                           rname, trim(ObsCompInfo(iObsComp)%name), Hops_m%Bscale(iObsComp); call goPr
            ! store:
            Bmats(iB)%sigma_scale_factors(itracer) = Hops_m%Bscale(iObsComp)
            ! apply to put out 'sf':
            sigma_m(:,:,:,iObsComp) = sigma_m(:,:,:,iObsComp) * Hops_m%Bscale(iObsComp)
          else
            ! info ...
            write (gol,'(a,":     comp `",a,"` not scaled")') &
                           rname, trim(ObsCompInfo(iObsComp)%name); call goPr
            ! default:
            Bmats(iB)%sigma_scale_factors(itracer) = 1.0
          end if
        end do ! tracers in covar
      end do  ! Bmat

      ! evaluate observation operator on new sigma fields,
      ! store as 'sigma forecast':
      call Hops_m%Evaluate( 'sf', sigma_m(:,:,nlev,:), sigma_m, status )!, debug=.true. )
      IF_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( sigma_m, stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! ...................................

#ifdef with_ajs
      !! testing ..
      !call goMem( rname//' - MEMORY 5a   ', status )
      !IF_NOT_OK_RETURN(status=1)
#endif

      ! initialize observation operator on analysis decomposition:
      call Hops_f%Init( status )
      IF_NOT_OK_RETURN(status=1)
      ! swap observation operator data (Hops_m) on model decomposition (doms_an_m)
      ! onto new observation operator (Hops_f) on analysis decomposition (doms_an_fg):
      call Hops_m%Swap( doms_an_m, Hops_f, doms_an_fg, status )
      IF_NOT_OK_RETURN(status=1)

      ! update timer:
      call Add_2timing(41,tim_after,tim_before,'3DVar: Get innovations from observations.')

#ifdef with_ajs
      ! end timing:
      call GO_Timer_End( itim_innov, status )
      IF_NOT_OK_RETURN(status=1)
#endif

      !-----------------------------------------------------------------------
      ! perform variational analysis
      !-----------------------------------------------------------------------

      ! info ...
      write (gol,'(a,": call var3d ...")') rname; call goPr

#ifdef with_ajs
      !! testing ..
      !call goMem( rname//' - MEMORY 6a   ', status )
      !IF_NOT_OK_RETURN(status=1)
#endif

      ! init array with all increments to zero:
      dx_loc = 0.0

      ! loop over B matrices:
      do iB = 1, nBmat

        ! info ...
        write (gol,'(a,":   B matrix ",i0," (",a,")")') rname, iB, trim(Bmats(iB)%name); call goPr
      
        ! init selected observations:
        call Hops_f_B%Init( status )
        IF_NOT_OK_RETURN(status=1)
        ! fill with copy of selected observed components only:
        call Hops_f_B%SelectTracers( Hops_f, Bmats(iB)%iObsComp, status, nana=nana )
        IF_NOT_OK_RETURN(status=1)

        ! total number of observations over all domains:
        call MPIF90_AllReduce( Hops_f_B%nobs, nobs_B, MPI_SUM, MPI_COMM_CALC, status )
        IF_NOT_OK_RETURN(status=1)
        call MPIF90_AllReduce( nana, nana_tot, MPI_SUM, MPI_COMM_CALC, status )
        IF_NOT_OK_RETURN(status=1)
        ! info ...
        write (gol,'(a,":     number of observations : ",i6," (analysed ",i0,", validation ",i0,")")') &
                            rname, nobs_B, nana_tot, nobs_B-nana_tot; call goPr

        ! any observations for this B ?
        if ( nobs_B > 0 ) then

          ! only tracers associated with this covar:
          if ( lny > 0 ) then
            ! storage for local concentration increments, new decomposition:
            allocate( dx_loc_B(lnx,lny,Bmats(iB)%nlev    ,Bmats(iB)%ntracer), stat=status )
            IF_NOT_OK_RETURN(status=1)
            allocate( du_loc_B(lnx,lny,Bmats(iB)%nlev_all,Bmats(iB)%ntracer), stat=status )
            IF_NOT_OK_RETURN(status=1)
          else
            ! dummy:
            allocate( dx_loc_B(1,1,1,1), stat=status )
            IF_NOT_OK_RETURN(status=1)
            allocate( du_loc_B(1,1,1,1), stat=status )
            IF_NOT_OK_RETURN(status=1)
          end if
          ! set to zero:
          dx_loc_B = 0.0
          du_loc_B = 0.0

          ! check units in B; loop over observed components:
          do iObsComp = 1, size(Bmats(iB)%itracer)
            ! mapping to tracer in this covariance:
            itracer = Bmats(iB)%itracer(iObsComp)
            ! undefined? then not in this covar:
            if ( itracer < 0 ) cycle
            ! get properties:
            call Bmats(iB)%BCovarSqrt%TracerGet( itracer, status, name=tracer_name, units=tracer_units )
            IF_NOT_OK_RETURN(status=1)
            ! compare:
            if ( trim(tracer_units) /= trim(xn_obs_units(iObsComp)) ) then
              write (gol,'("covariance tracer `",a,"` has units `",a,"` while xn_obs units are `",a,"`")') &
                  trim(tracer_name), trim(tracer_units), trim(xn_obs_units(iObsComp)); call goErr
              TRACEBACK; status=1; return
            end if
          end do  ! observed components

          ! perform analysis, all processes involved for operations with H;
          ! only the tracer values assigned to this B matrix are changed;
          ! dx might have only 1 level (sfc, or total column),
          ! but du is the update on all levels according to vertical correlations:
          call var3d( cdate, Hops_f_B, Bmats(iB), dx_loc_B, status, &
                         du_loc=du_loc_B )
          IF_NOT_OK_RETURN(status=1)

          ! info ...
          if ( Bmats(iB)%nlev == 1 ) then
            write (gol,'(a,":     copy component update (2D to 3D) ...")') rname; call goPr
          else
            write (gol,'(a,":     copy component update (3D) ...")') rname; call goPr
          end if
          ! restore:
          if ( lny > 0 ) then
            ! copy:
            do itracer = 1, Bmats(iB)%ntracer
              ! index of observed component for this tracer:
              iobs = Bmats(iB)%iObsComp(itracer)
              ! 2D to 3D ?
              if ( Bmats(iB)%nlev == 1 ) then
                ! copy from 'unobserved' levels:
                dx_loc(:,:,:,iobs) = du_loc_B(:,:,:,itracer)
              else
                ! copy, dx on all levels already:
                dx_loc(:,:,:,iobs) = dx_loc_B(:,:,:,itracer)
              end if
            end do ! tracers
          end if ! local slab?

          ! clear:
          deallocate( dx_loc_B, stat=status )
          IF_NOT_OK_RETURN(status=1)
          deallocate( du_loc_B, stat=status )
          IF_NOT_OK_RETURN(status=1)

        end if  ! nobs_B > 0

        ! done:
        call Hops_f_B%Done( status )
        IF_NOT_OK_RETURN(status=1)

      end do ! B matrices

#ifdef with_ajs
      !! testing ..
      !call goMem( rname//' - MEMORY 6b   ', status )
      !IF_NOT_OK_RETURN(status=1)
#endif


      !-----------------------------------------------------------------------
      ! read result for \delta x and \delta u and add to background field;
      !-----------------------------------------------------------------------
      ! Only update whithin the non-extended zone (since there are no obs. in
      ! that zone => zero increment). Leave BCs (lateral & top) unchanged.
      !-----------------------------------------------------------------------

      ! info ...
      write (gol,'(a,": swap increments to model decomposition ...")') rname; call goPr

      !! rescale if necessary:
      !if ( FGSCALE_INV /= 1e0 ) dx_loc = dx_loc * FGSCALE_INV

      ! swap to original model decomposition:
      call doms_an_fg%Swap( dx_loc, doms_an_m, dx_an, status )
      IF_NOT_OK_RETURN(status=1)

      ! info ...
      write (gol,'(a,": reset concentrations in model array ...")') rname; call goPr

      ! loop over variables involved in analysis:
      do iObsComp = 1, nObsComp
      
        ! NOTE: need routine that distributes dx over sf and xn
        !  given the obsdata type (2D-OBS, 2D-ML, or 3D-ML)

        ! store if needed ;
        ! corrected fields are original plus delta,
        ! for surface just use the lowest level
        call Fill_Output_sf_xn( '3DVAR_OBS_AN', trim(ObsCompInfo(iObsComp)%name), &
                                   sf_an(:,:,iObsComp)+dx_an(:,:,size(dx_an,3),iObsComp), &
                                   xn_an(:,:,:,iObsComp)+dx_an(:,:,:,iObsComp), &
                                   trim(ObsCompInfo(iObsComp)%units), status )
        IF_NOT_OK_RETURN(status=1)

        ! simulate observed components (surface and 3D)
        ! from xn_adv (first guess) plus increments (and also xn_an if relative change is needed):
        call ObsCompInfo(iObsComp)%FillFields( xn_adv, &
                    sf_an(:,:,iObsComp), xn_an(:,:,:,iObsComp), xn_obs_units(iObsComp), status, &
                    dx_obs=dx_an(:,:,:,iObsComp) )
        IF_NOT_OK_RETURN(status=1)

        ! apply feedback? 
        if ( ObsCompInfo(iObsComp)%feedback ) then
          ! info ...
          write (gol,'(a,":   feedback ",a," ...")') rname, trim(ObsCompInfo(iObsComp)%name); call goPr
          ! distribute increment over original species:
          call ObsCompInfo(iObsComp)%DistributeIncrement( xn_adv, &
                        xn_an(:,:,:,iObsComp), dx_an(:,:,:,iObsComp), status, &
                        maxratio=ANALYSIS_RELINC_MAX )
          IF_NOT_OK_RETURN(status=1)
        else
          ! info ...
          write (gol,'(a,":   do not feedback ",a," (yet) ...")') rname, trim(ObsCompInfo(iObsComp)%name); call goPr
        end if

      end do


      ! analysed output
      ! - surface and 3D fields
      ! - might not have been feeded back to the model)

      ! fill arrays for requested subclass using analysed concentrations:
      call Fill_Output_xn_adv( '3DVAR_AN', xn_adv, status )
      IF_NOT_OK_RETURN(status=1)

#ifdef with_ajs
      !! testing ..
      !call goMem( rname//' - MEMORY 6c   ', status )
      !IF_NOT_OK_RETURN(status=1)
#endif


      !-----------------------------------------------------------------------
      ! innovations after analysis
      !-----------------------------------------------------------------------

      ! info ...
      write (gol,'(a,": evaluate observations after analysis ...")') rname; call goPr

      ! evaluate observation operator, store in it's 'xa' field:
      call Hops_m%Evaluate( 'xa', sf_an, xn_an, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! info ...
      write (gol,'(a,": write innovations ...")') rname; call goPr

      ! save all:
      call Hops_m%WriteToFile( cdate, status )
      IF_NOT_OK_RETURN(status=1)

#ifdef with_ajs
      !! testing ..
      !write (gol,'("break after Hops write")'); call goErr
      !TRACEBACK; status=1; return

      !! testing ..
      !call goMem( rname//' - MEMORY 6d   ', status )
      !IF_NOT_OK_RETURN(status=1)
#endif


      !-----------------------------------------------------------------------
      ! deallocate analysis arrays
      !-----------------------------------------------------------------------

      ! info ...
      write (gol,'(a,": clear analysis arrays ...")') rname; call goPr

      ! done with obs operators:
      call Hops_f%Done( status )
      IF_NOT_OK_RETURN(status=1)

#ifdef with_ajs
      !! testing ..
      !call goMem( rname//' - MEMORY 5b   ', status )
      !IF_NOT_OK_RETURN(status=1)
#endif

      ! done with obs operators:
      call Hops_m%Done( status )
      IF_NOT_OK_RETURN(status=1)

#ifdef with_ajs
      !! testing ..
      !call goMem( rname//' - MEMORY 4b   ', status )
      !IF_NOT_OK_RETURN(status=1)
#endif

      ! clear:
      deallocate( dx_loc, stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! clear:
      deallocate( sf_an, stat=status  )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xn_an, stat=status  )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dx_an, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xn_obs_units, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xn_adv_units, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! observations present

#ifdef with_ajs
    !! testing ..
    !call goMem( rname//' - MEMORY 3b   ', status )
    !IF_NOT_OK_RETURN(status=1)
#endif


    !-----------------------------------------------------------------------
    ! deallocate observation arrays
    !-----------------------------------------------------------------------

    ! info ...
    write (gol,'(a,": clear observation arrays ...")') rname; call goPr

    ! clear:
    deallocate( stnid, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( stncodes, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( flat, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( flon, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( falt, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( obs, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( obsstddev1, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( obsanalyse, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( iObsData, stat=status )
    IF_NOT_OK_RETURN(status=1)

#ifdef with_ajs
    !! testing ..
    !call goMem( rname//' - MEMORY 2b   ', status )
    !IF_NOT_OK_RETURN(status=1)

    !! testing ..
    !call goMem( rname//' - MEMORY end  ', status )
    !IF_NOT_OK_RETURN(status=1)
#endif

    !! testing ...
    !nmemtest = nmemtest + 1
    !if ( nmemtest >= 5 ) then
    !  write (gol,'(a," - MEMORY break")') rname; call goPr
    !  TRACEBACK; status=1; return
    !end if

    ! ok
    status = 0

  end subroutine generic3dvar
  

  ! ***



  subroutine Fill_Output_xn_adv( subclass, xn_adv, status )

    use DerivedFields_ml     , only : d_2d, f_2d, d_3d, f_3d
    use ModelConstants_ml    , only : IOU_INST, num_lev3d, lev3d
    use ModelConstants_ml    , only : KMAX_MID
    use ChemGroups_ml        , only : PMFINE_GROUP, PM10_GROUP   ! requires offset NSPEC_SHL !
    use ChemSpecs_shl_ml     , only : NSPEC_SHL  ! number of short-lived tracers
    use ChemChemicals_ml     , only : species_adv
    use Chemfields_ml        , only : cfac
    use Chemfields_ml        , only : PM25_water_rh50  ! (i,j) surface
    use Chemfields_ml        , only : PM25_water       ! (i,j,k) model levels
    use MetFields_ml         , only : roa
    use Units_ml             , only : Units_Scale
    use SmallUtils_ml        , only : find_index

    ! --- in/out ----------------------------
    
    character(len=*), intent(in)  :: subclass
    real, intent(in)              :: xn_adv(:,:,:,:)
    integer, intent(out)          :: status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Fill_Output_xn_adv'

    ! --- local -----------------------------

    integer               ::  iout
    integer, pointer      ::  out_group(:)
    integer, target       ::  out_group1(1)
    integer               ::  out_offset
    integer               ::  ispec
    real                  ::  fscale
    logical               ::  needroa
    logical               ::  addwater
    integer               ::  k
    integer               ::  ilev

    ! --- begin -----------------------------

    ! loop over 2D output fields:
    do iout = 1, size(f_2d)
      ! filter:
      if ( trim(f_2d(iout)%class   ) /= 'USET'         ) cycle
      if ( trim(f_2d(iout)%subclass) /= trim(subclass) ) cycle
      ! set list of contributing species:
      select case ( trim(f_2d(iout)%txt) )
        !~ total fine pm:
        case ( 'PM25' )
          ! use group defined in 'ChemGroups_ml':
          out_group => PMFINE_GROUP
          ! offset:
          out_offset = NSPEC_SHL
          ! include aerosol water:
          addwater = .true.
        !~ total coarse pm:
        case ( 'PM10' )
          ! use group defined in 'ChemGroups_ml':
          out_group => PM10_GROUP
          ! offset:
          out_offset = NSPEC_SHL
          ! include aerosol water:
          addwater = .true.
        !~ expected a single model species ...
        case default
          ! find index:
          out_group1(1) = find_index( trim(f_2d(iout)%txt), species_adv(:)%name )
          ! assign:
          out_group => out_group1
          ! no offset:
          out_offset = 0
          ! no aerosol water:
          addwater = .false.
      end select
      ! no bulk unit conversion, apply scalings below:
      f_2d(iout)%scale = 1.0
      ! init sum:
      d_2d(iout,:,:,IOU_INST) = 0.0
      ! loop over group members:
      do k = 1, size(out_group)
        ! extract index of specie:
        ispec = out_group(k) - out_offset
        ! info on conversion to target units:
#ifdef with_ajs
        call Units_Scale( f_2d(iout)%unit, ispec, fscale, &
                                 needroa=needroa, status=status )
        IF_NOT_OK_RETURN(status=1)
#else
        call Units_Scale( f_2d(iout)%unit, ispec, fscale, &
                                 needroa=needroa )
#endif
        ! extract from bottom level:
        ilev = KMAX_MID
        ! convert using air density?
        if ( needroa ) then
          ! add contribution, multiply with air density:
          d_2d(iout,:,:,IOU_INST) = d_2d(iout,:,:,IOU_INST) + &
                                      xn_adv(ispec,:,:,ilev) &
                                      * fscale &
                                      * roa(:,:,ilev,1) &
                                      * cfac(ispec,:,:)
        else
          ! add contribution:
          d_2d(iout,:,:,IOU_INST) = d_2d(iout,:,:,IOU_INST) + &
                                      xn_adv(ispec,:,:,ilev) &
                                      * fscale &
                                      * cfac(ispec,:,:)
        end if
      end do  ! spec
      ! add aerosol water?
      if ( addwater ) then
        ! add surface aersosol water:
        d_2d(iout,:,:,IOU_INST) = d_2d(iout,:,:,IOU_INST) + PM25_water_rh50
      end if ! add water
    end do  ! 2D output fields

    ! loop over 3D output fields:
    do iout = 1, size(f_3d)
      ! filter:
      if ( trim(f_3d(iout)%class   ) /= 'USET'         ) cycle
      if ( trim(f_3d(iout)%subclass) /= trim(subclass) ) cycle
      ! set list of contributing species:
      select case ( trim(f_2d(iout)%txt) )
        !~ total fine pm:
        case ( 'PM25' )
          ! use group defined in 'ChemGroups_ml':
          out_group => PMFINE_GROUP
          ! offset:
          out_offset = NSPEC_SHL
          ! include aerosol water:
          addwater = .true.
        !~ total coarse pm:
        case ( 'PM10' )
          ! use group defined in 'ChemGroups_ml':
          out_group => PM10_GROUP
          ! offset:
          out_offset = NSPEC_SHL
          ! include aerosol water:
          addwater = .true.
        !~ expected a single model species ...
        case default
          ! find index:
          out_group1(1) = find_index( trim(f_2d(iout)%txt), species_adv(:)%name )
          ! assign:
          out_group => out_group1
          ! no offset:
          out_offset = 0
          ! no aerosol water:
          addwater = .false.
      end select
      ! no bulk unit conversion, apply scalings below:
      f_3d(iout)%scale = 1.0
      ! init sum:
      d_3d(iout,:,:,:,IOU_INST) = 0.0
      ! loop over group members:
      do k = 1, size(out_group)
        ! extract index of specie:
        ispec = out_group(k) - out_offset
        ! info on conversion to target units:
#ifdef with_ajs
        call Units_Scale( f_3d(iout)%unit, ispec, fscale, &
                                needroa=needroa, status=status )
        IF_NOT_OK_RETURN(status=1)
#else
        call Units_Scale( f_3d(iout)%unit, ispec, fscale, &
                                needroa=needroa )
#endif
        ! convert using air density?
        if ( needroa ) then
          ! add contribution, multiply with air density:
          d_3d(iout,:,:,:,IOU_INST) = d_3d(iout,:,:,:,IOU_INST) + &
                                        xn_adv(ispec,:,:,lev3d(:num_lev3d)) &
                                        * fscale &
                                        * roa(:,:,lev3d(:num_lev3d),1)
        else
          ! add contribution:
          d_3d(iout,:,:,:,IOU_INST) = d_3d(iout,:,:,:,IOU_INST) + &
                                        xn_adv(ispec,:,:,lev3d(:num_lev3d)) &
                                        * fscale
        end if
      end do  ! spec
      ! add aerosol water?
      if ( addwater ) then
        ! add surface aersosol water:
        d_3d(iout,:,:,:,IOU_INST) = d_3d(iout,:,:,:,IOU_INST) + PM25_water
      end if ! add water
    end do  ! 3D output fields
  
    ! ok
    status = 0

  end subroutine Fill_Output_xn_adv
  

  ! ***



  subroutine Fill_Output_sf_xn( subclass, txt, sf, xn, units, status )

    use DerivedFields_ml     , only : d_2d, f_2d, d_3d, f_3d
    use ModelConstants_ml    , only : IOU_INST, num_lev3d, lev3d

    ! --- in/out ----------------------------
    
    character(len=*), intent(in)  :: subclass
    character(len=*), intent(in)  :: txt
    character(len=*), intent(in)  :: units
    real, intent(in)              :: sf(:,:)
    real, intent(in)              :: xn(:,:,:)
    integer, intent(out)          :: status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/Fill_Output_sf_xn'

    ! --- local -----------------------------

    integer               ::  iout

    ! --- begin -----------------------------

    ! loop over 2D output fields:
    do iout = 1, size(f_2d)
      ! filter:
      if ( trim(f_2d(iout)%class   ) /= 'USET'         ) cycle
      if ( trim(f_2d(iout)%subclass) /= trim(subclass) ) cycle
      if ( trim(f_2d(iout)%txt     ) /= trim(txt     ) ) cycle
      ! check ...
      if ( trim(units) /= trim(f_2d(iout)%unit) ) then
        write (gol,'("units `",a,"1 of txt `",a,"` from subclass `",a,"` does not match with output units `",a,"`")') &
                  trim(units), trim(txt), trim(subclass), trim(f_2d(iout)%unit); call goErr
        TRACEBACK; status=1; return
      end if
      ! no bulk unit conversion, units already checked:
      f_2d(iout)%scale = 1.0
      ! store field:
      d_2d(iout,:,:,IOU_INST) = sf
      ! found
      exit
    end do ! 2D output fields

    ! loop over 3D output fields:
    do iout = 1, size(f_3d)
      ! filter:
      if ( trim(f_3d(iout)%class   ) /= 'USET'         ) cycle
      if ( trim(f_3d(iout)%subclass) /= trim(subclass) ) cycle
      if ( trim(f_3d(iout)%txt     ) /= trim(txt     ) ) cycle
      ! check ...
      if ( trim(units) /= trim(f_3d(iout)%unit) ) then
        write (gol,'("units `",a,"1 of txt `",a,"` from subclass `",a,"` does not match with output units `",a,"`")') &
                  trim(units), trim(txt), trim(subclass), trim(f_3d(iout)%unit); call goErr
        TRACEBACK; status=1; return
      end if
      ! no bulk unit conversion, units already checked:
      f_3d(iout)%scale = 1.0
      ! store field:
      d_3d(iout,:,:,:,IOU_INST) = xn(:,:,lev3d(:num_lev3d))
      ! found
      exit
    end do  ! 3D output fields
  
    ! ok
    status = 0

  end subroutine Fill_Output_sf_xn
  

  ! ***


  !-----------------------------------------------------------------------
  subroutine var3d( cdate, Hops, Bmat, dx_loc, status, &
                       du_loc )
  !-----------------------------------------------------------------------
  ! @description
  ! 3-D variational analysis
  ! @author M.Kahnert
  !-----------------------------------------------------------------------

    use MPIF90               , only : MPIF90_BCast
    use MPIF90               , only : MPIF90_AllReduce, MPI_SUM
#ifdef with_ajs
    use DA_Util_ml           , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
#endif
    use ModelConstants_ml    , only : MasterProc, NPROC
    use MPI_Groups_ml        , only : MasterPE
    use My_Timing_ml         , only : Code_timer, Add_2timing
    use DA_Util_ml           , only : norm
    use Io_ml                , only : m1qn3_io=>IO_TMP  ! write unit for m1qn3 printouts
    use TimeDate_ml          , only : date,current_date ! date/time structure
    use TimeDate_ExtraUtil_ml, only : date2string
    use Par_ml               , only : me
    use DA_ml                , only : debug => DEBUG_DA
    use DA_ml                , only : dafmt => da_fmt_msg, damsg => da_msg
    use DA_ml                , only : tim_before => datim_before, tim_after => datim_after
    use DA_Obs_ml            , only : T_ObsOpers
    use DA_ml                , only : nlev
    !use DA_ml                , only : nChemObs

    !-----------------------------------------------------------------------
    ! Formal parameters
    !-----------------------------------------------------------------------

    type(date), intent(in)            ::  cdate
    type(T_ObsOpers), intent(inout)   ::  Hops
    type(T_BInfo), intent(inout)      ::  Bmat
    real, intent(out)                 ::  dx_loc(:,:,:,:)  ! (limax,ljmax,nlev,ntracer)
    integer, intent(out)              ::  status

    real, intent(out), optional       ::  du_loc(:,:,:,:)  ! (limax,ljmax,nlev_all,ntracer)

    !-----------------------------------------------------------------------
    ! external routines
    !-----------------------------------------------------------------------

    ! functions defined in m1qn3 :
    external simul_rc,euclid,ctonbe,ctcabe

    !-----------------------------------------------------------------------
    ! parameters
    !-----------------------------------------------------------------------

    character(len=*), parameter   ::  rname = mname//'/var3d'

    character(len=*), parameter   ::  normtype='dfn'

    !! configuration:
    !integer, parameter ::   imode(3)=(/ 1,&   ! run in SIS mode
    !                                    0,&   ! cold start
    !                                    0/)   ! calculate Jcost,gradJcost every interation

    !
    !integer, parameter            ::  nupdates = 4

    ! norm resolution
    ! ~ originally copied from 'DAPREC' parameter in spectralcov.f90,
    !   where it was set to a value based on real number represenation:
    !     ! 2**-22=2.38e-07 1/real4 wont't overflow
    !     real, parameter :: DAPREC = 1e0/2**22
    !real(kind=8), parameter       ::  dxmin = DAPREC
    ! ~ same value, but parameterized here:
    real(kind=8), parameter       ::  dxmin = 1e0/2**22

    character(len=*),parameter :: LOGFILE="YYYYMMDDhh-PPP_m1qn3.log"

    !-----------------------------------------------------------------------
    ! Local parameters
    !-----------------------------------------------------------------------

    integer                   ::  itime
    integer                   ::  nw, nw_local
    integer                   ::  nv_hcr
    real, allocatable         ::  chi_hcr(:)
    real, allocatable         ::  gradJcost_hcr(:)
    integer, allocatable      ::  l2w_hcr(:)

    ! global arrays on each pe needed for calls to m1qn3:
    integer                   ::  nhcr
    integer                   ::  nv_hcr_all

    real                      ::  Jcost0, Jcost, Jcost_b, Jcost_obs
    real                      ::  gradnorm_loc, gradnorm

    ! m1qn3 variables:
    integer                   ::  iz(5)
    integer                   ::  ndz, ndz_all
    integer                   ::  indic
    integer                   ::  imode(3)
    integer                   ::  omode
    integer                   ::  niter
    integer                   ::  nsim
    integer                   ::  reverse
    integer                   ::  impres
    integer                   ::  istep
    real                      ::  epsg,df1
    real(kind=8), allocatable ::  dz(:)
    integer                   ::  izs(1) = 0
    real(kind=4)              ::  rzs(1) = 0.0
    real(kind=8)              ::  dzs(1) = 0.0

    integer                   ::  n,p,l,k

    !-----------------------------------------------------------------------
    ! storage
    !-----------------------------------------------------------------------

    ! get index of record holding current hour:
    call Bmat%BcovarSqrt%FindTime( cdate%hour, itime, status )
    IF_NOT_OK_RETURN(status=1)

    ! current size of half-complex state,
    ! local as well as global:
    call Bmat%BcovarSqrt%TimeGet( itime, status, &
                                obs_nw_local=nw_local, obs_nw=nw )
    IF_NOT_OK_RETURN(status=1)

    ! size of half-complex-real state is double:
    nv_hcr     = 2 * nw_local
    nv_hcr_all = 2 * nw

    !! info ...
    !write (gol,'(a,": horizontal size of half-complex-to-real (local) : ",i0," (",i0,")")') rname, nhcr, lnhcr; call goPr
    !write (gol,'(a,": total      size of half-complex-to-real (local) : ",i0," (",i0,")")') rname, nv_hcr_all, nv_hcr; call goPr

    ! trap zero length:
    if ( nv_hcr > 0 ) then
      ! storage for state provided to optimizer:
      allocate( chi_hcr(nv_hcr), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! weight for l2 norm:
      allocate( l2w_hcr(nv_hcr), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! gradient vector:
      allocate( gradJcost_hcr(nv_hcr), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! storage for state provided to optimizer:
      allocate( chi_hcr(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! weight for l2 norm:
      allocate( l2w_hcr(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! gradient vector:
      allocate( gradJcost_hcr(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if


    !-----------------------------------------------------------------------
    ! Make one first call to the objective function
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! optimizer loop
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Call minimization routine (m1qn3-3.3;2009)
    ! subroutine m1qn3 (simul, prosca, ctonb, ctcab, n, x, f, g,
    !                   dxmin, df1, epsg, normtype, impres, io,
    !                   imode, omode, niter, nsim, iz, dz, ndz,
    !                   reverse, indic, izs, rzs, dzs)
    ! Very old call minimization routine (m1qn3-2.0d;1996)
    ! subroutine m1qn3 (simul, prosca, ctonb, ctcab, n, x, f, g,
    !                   dxmin, df1, epsg,           impres, io,
    !                   mode,         niter, nsim, iz, rz, nrz)
    !-----------------------------------------------------------------------

    ! initial state:
    chi_hcr = 0.0
    l2w_hcr = 0

    ! thresholds:
    epsg = dxmin**0.9

    ! reverse communication, call cost function (in parallel) outside optimizer:
    reverse = 1

    ! max counters:
    niter = solver_maxiter
    nsim  = solver_maxiter*solver_maxsim

    ! size of work array:
    ndz     = 4*nv_hcr     + solver_nupdates*(2*nv_hcr    +1)
    ndz_all = 4*nv_hcr_all + solver_nupdates*(2*nv_hcr_all+1)
    ! storage:
    allocate( dz(ndz), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! settings:
    imode(1) = solver_scaling ! 0=DIS, 1=SIS
    imode(2) = 0              ! cold start
    imode(3) = 0              ! calculate Jcost,gradJcost every interation

    if ( masterProc .or. solver_logall ) then
      impres = 3  ! verbosity level: 0=No print .. 5=Full verbosity
    else
      impres = 0  ! No print
    endif

    ! info ...
    if(debug.and.MasterProc.and.impres>0)then
      write(damsg,"(A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
        'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
      print dafmt,'Calling m1qn3 '//trim(damsg)
    endif
    if(impres>0)then
      open(m1qn3_io,file=date2string(LOGFILE,current_date),iostat=status)
      IF_NOT_OK_RETURN(status=1)
    endif

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_loop, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! start timing:
    call Code_timer(tim_before)

    ! increase counter:
    solver_number = solver_number + 1

    ! reverse mode, 'costFunction' is called externally;
    ! perform loop over optimer steps:
    istep = 0
    do

      ! increase counter:
      istep = istep + 1

      !! info ...
      !write (gol,'(a,": solver step ",i4," ...")') rname, istep; call goPr
      !! testing ...
      !if ( istep > 3 ) then
      !  write (gol,'(a,": WARNING - break iteration loop")') rname; call goPr
      !  exit
      !end if

#ifdef with_ajs
      ! start timing:
      call GO_Timer_Start( itim_costfunc, status )
      IF_NOT_OK_RETURN(status=1)
#endif

      !! info ...
      !write (gol,'(a,":   evaluate cost function and gradient ...")') rname; call goPr
      ! next evaluation of costfunction and gradient;
      ! for operations with H this involves all processes :
      call costFunction( itime, nv_hcr, chi_hcr, Hops, Bmat, &
                           Jcost, Jcost_b, Jcost_obs, gradJcost_hcr, l2w_hcr, &
                           dx_loc, .false., status )
      IF_NOT_OK_RETURN(status=1)

      !! info ...
      !write (gol,'(a,":     Jcost    : ",e16.6)') rname, Jcost; call goPr
      !write (gol,'(a,":     dx range : ",2e16.6)') rname, minval(dx_loc), maxval(dx_loc); call goPr

      ! only on root ...
      if ( masterProc ) then
        ! add cost function values to summary:
        write (solver_table,'(i,",",i,4(",",e))',iostat=status) &
                 solver_number, istep, Jcost, Jcost_b, Jcost_obs, norm(gradJcost_hcr)
        IF_NOT_OK_RETURN(status=1)
      end if

      !! info ...
      !write (gol,'(a,":     call m1qn3 ...")') rname; call goPr
      !write (gol,'(a,":       state size : ",i0," (local ",i0,")")') rname, nv_hcr_all, nv_hcr; call goPr
      !write (gol,'(a,":       work  size : ",i0," (local ",i0,")")') rname, ndz_all, ndz; call goPr

#ifdef with_ajs
      ! switch timing:
      call GO_Timer_Switch( itim_costfunc, itim_optimizer, status )
      IF_NOT_OK_RETURN(status=1)
#endif

      ! first step ?
      if ( istep == 1 ) then
        ! store initial evaluation of cost function:
        Jcost0 = Jcost
        ! set threshold:
        df1 = Jcost * 0.5
      end if

      ! Set indicator to inform m1qn3 about latest evaluation:
      !  < 0  : cost function could not be evaluated at given chi ;
      !  = 0  : m1qn3 has to stop ;
      !  > 1  : cost funcation and gradient have been evaluated correctly.
      ! Here we assume that evaluation was always possible,
      ! or that otherwise errors have been trapped already:
      indic = 1

      !
      ! NOTE on the half-complex storage!
      !
      !  The spectral fields have the size of the extended grid.
      !  For feeding them to m1qn3 routine these have to be stored
      !  as a single 1D array of real numbers, thus having size:
      !    2 * nxex * nyex * (nlev*nchemobs)
      !  where the factor 2 is need to store the real and imag parts.
      !  However, since the complex values are Fourrier transforms of real values,
      !  about half of them are conjugates of the other half.
      !  It is therefor not necessary to store all of them,
      !  half is enough to describe the complete field ("half-complex storage").
      !  This was used in the original code too when reformatting the
      !  full complex array into a 1D real vector ; only half of the values
      !  in the "y" direction (actually a spectral direction) were stored.
      !  Including the complex-to-2-reals packing the size became:
      !     2 * nxex * nyex/2 * (nlev*nchemobs)
      !
      !  A problem with the half-complex-to-two-reals packing is the computation
      !  of an equivalent for the Hermitian product:
      !       x^H y = conj(x)^T y
      !  This is the base for the l2-norm over a vector with complex numbers:
      !      ||x||^2 = x^H x
      !  The norm is used by the minimization algorithm to decide on convergence.
      !  For a full-complex field with the half-is-conjugate feature it can be shown that:
      !      x^H y  = (Re[x],Im[x])^T (Re[y],Im[y])
      !  thus the Hermitian-product can be computed as the regular dot-product
      !  over vectors containing the real and imaginary parts of the vectors.
      !  Thus, a complex-to-real packing maintains the Hermitian product.
      !  To remain this property for a half-complex-to-real packing,
      !  an element wise weight should be used in the dot product:
      !      x^H y = (Re[x_hc],Im[x_hc])^T D_hc (Re[y_hc],Im[y_hc])
      !  where D_hc is a diagonal matrix with values :
      !      1  if the element originates from a real    value in the Fourrier spectrum
      !      2  "  "   "       "          "    " complex "     "  "   "        "
      !
      !  In the original code this was ignored ;
      !  now introduced, but results will therefore be slightly different.
      !

      !
      ! Optimizer step;
      !  - provided 'simul_rc' is dummy for reverse communication;
      !  - special l2norm for half-complex-to-real states ;
      !    provide integer weights to 'izs' user array which is used
      !    by the 'euclid_hcr' norm routine
      !  - write messages to seperate log file;
      !  - provide only local part of the states 'chi' and 'gradJcost',
      !    and use for the dot-product a 'euclid_hcr' function which
      !    computes a global sum;
      !  - size of the problem should at least 1, even if on this
      !    processor no elements are present ;
      !    therefore pass the actual size of the state vector
      !
      call m1qn3( simul_rc, euclid_hcr, ctonbe, ctcabe, &
                  nv_hcr_all, nv_hcr, chi_hcr, Jcost, gradJcost_hcr, &
                  dxmin, df1, epsg, normtype, impres, m1qn3_io, &
                  imode, omode, niter, nsim, iz, &
                  dz, ndz_all, ndz, &
                  reverse, indic, &
                  l2w_hcr, rzs, dzs )

#ifdef with_ajs
      ! end timing:
      call GO_Timer_End( itim_optimizer, status )
      IF_NOT_OK_RETURN(status=1)
#endif

      !! info ...
      !write (gol,'(a,":     returned reverse : ",i4)') rname, reverse; call goPr
      !
      ! check value:
      if ( reverse < 0 ) then
        ! interupt call loop:
        exit
      else if ( reverse == 1 ) then
        ! continue
      else
        write (gol,'("received unsupported reverse value : ",i0)') reverse; call goErr
        TRACEBACK; status=1; return
      end if

      !! info ...
      !write (gol,'(a,":     returned indic   : ",i4)') rname, indic; call goPr
      ! check value:
      select case ( indic )
        !~ do not call simulator; this should have been trapped already:
        case ( 1 )
          write (gol,'("unexpected indic ",i0," should not occure")') indic; call goErr
          TRACEBACK; status=1; return
        !~ request new J(chi) and nablaJ(chi):
        case ( 4 )
          ! continue
        !~ unknown:
        case default
          write (gol,'("unexpected indic value ",i0)') indic; call goErr
          TRACEBACK; status=1; return
      end select

    end do  ! optimizer steps

    ! info ....
    if(debug.and.MasterProc.and.impres>0)then
      write(damsg,"(A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
        'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
      print dafmt,'Finish  m1qn3 '//trim(damsg)
    endif
    if(impres>0)then
      close(m1qn3_io,iostat=status)
      IF_NOT_OK_RETURN(status=1)
    endif

    !! info ...
    !write (gol,'(a,":     returned omode   : ",i4)') rname, omode; call goPr
    ! check ...
    select case ( omode )
      case ( 0 )
        write (gol,'("omode==0: The simulator asks to stop by returning the value (indic=0)")'); call goErr
        TRACEBACK; status=1; return
      case ( 1 )
        write (gol,'("omode==1: Normal m1qn3 exit: successfull gradient test")'); call goPr
      case ( 2 )
        write (gol,'("omode==2: One of the input arguments is not well initialized")'); call goErr
        TRACEBACK; status=1; return
      case ( 3 )
        write (gol,'("omode==3: Line-search blocked on tmax = 10**20")'); call goErr
        TRACEBACK; status=1; return
      case ( 4 )
        write (gol,'("omode==4: Reached maximal number of iterations (maxiter)")'); call goPr
      case ( 5 )
        write (gol,'("omode==5: Reached maximal number of simulations (maxsim)")'); call goPr
      case ( 6 )
        write (gol,'("omode==6: Stop on dxmin during the line-search")'); call goPr
      case ( 7 )
        write (gol,'("omode==7: Either <g,d> is nonnegative or <y,s> is nonpositive")'); call goErr
        TRACEBACK; status=1; return
      case default
        write (gol,'("unsupported omode==",i0)') omode; call goErr
        TRACEBACK; status=1; return
    end select

    ! finish timing:
    call Add_2timing(42,tim_after,tim_before,'3DVar: Optimization.')

#ifdef with_ajs
    ! end timing:
    call GO_Timer_End( itim_loop, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    !-----------------------------------------------------------------------
    ! Make a final call to costFunction for converting \chi to \delta x.
    !-----------------------------------------------------------------------

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_post, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! postprocessing evaluation,
    ! computes final dx and evaluates du :
    call costFunction( itime, nv_hcr, chi_hcr, Hops, Bmat, &
                        Jcost, Jcost_b, Jcost_obs, gradJcost_hcr, l2w_hcr, &
                        dx_loc, .true., status, &
                        du_loc=du_loc )
    IF_NOT_OK_RETURN(status=1)

    ! info ...
    if ( MasterProc ) then
      write (damsg,*) Jcost0, '-->', Jcost, '=', (1.0-Jcost/Jcost0)*100
      print dafmt,'Cost function '//trim(ADJUSTL(damsg))//'% Reduction'
    end if

#ifdef with_ajs
    ! end timing:
    call GO_Timer_End( itim_post, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    !-----------------------------------------------------------------------
    ! done
    !-----------------------------------------------------------------------

    ! clear:
    deallocate( chi_hcr, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( l2w_hcr, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( gradJcost_hcr, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( dz, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine var3d


  ! ***


  !
  ! Dot-product function to be used by M1QN3.
  !
  ! Original template:
  !    m1qn3.f/euclid
  !
  ! Here result weighted sum over the element products,
  ! where the weights are either 1 or 2 to compensate
  ! for the contribution of the complex conjugates
  ! that are orignally present.
  ! The weights should be passed to the 'izs' argument
  ! with integer user data:
  !   izs(1:n)   : weights
  !
  ! Implemented as collective call over all processes,
  ! each processor computes a local sum first which is
  ! then summed over all.
  !

  subroutine euclid_hcr( n, x, y, ps, izs, rzs, dzs )

    use MPIF90       , only : MPIF90_AllReduce, MPI_SUM
    use MPI_Groups_ml, only : MPI_COMM_CALC

    ! --- in/out ---------------------------------

    integer, intent(in)     ::  n
    real(8), intent(in)     ::  x(n)
    real(8), intent(in)     ::  y(n)
    real(8), intent(out)    ::  ps
    integer, intent(in)     ::  izs(*)
    real(4), intent(in)     ::  rzs(*)
    real(8), intent(in)     ::  dzs(*)

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/euclid_hcr'

    ! --- local ----------------------------------

    integer               ::  i
    real(8)               ::  ps_loc
    integer               ::  status

    ! --- begin ----------------------------------

    ! init local sum:
    ps_loc = 0.d0
    ! add elements:
    do i = 1, n
      ps_loc = ps_loc + x(i) * y(i) * izs(i)
    end do

    ! collect global sum:
    call MPIF90_AllReduce( ps_loc, ps, MPI_SUM, MPI_COMM_CALC, status )
    if ( status /=0 ) then
      write (*,'("ERROR in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__
      stop
    end if

  end subroutine euclid_hcr


  ! ***


  subroutine get_innovations( Hops, nobs, &
                                xn_adv_b, xn_adv_units, &
                                sf_b, xn_b, xn_obs_units, &
                                maxobs, stnid, flat,flon,falt, &
                                obs, obs_stddev, stncodes, obs_analyse, &
                                state, status )
  !-----------------------------------------------------------------------
  ! @description
  ! Compute innovations H(xb)-y, where xb is
  ! the background field, y denotes the observations, and
  ! H is an operator mapping from model to observation space.
  ! On input, innov contains the observations y.
  ! On output, innov contains the innovations.
  ! @author M.Kahnert
  !-----------------------------------------------------------------------
    use Par_ml               , only : me
    use ModelConstants_ml    , only : MasterProc
    use GridValues_ml        , only : coord_in_processor
    use DA_ml                , only : debug => DEBUG_DA
    use DA_Obs_ml            , only : T_ObsOpers
    use DA_Obs_ml            , only : obsData
    !use DA_ml                , only : FGSCALE, FGSCALE_INV
    !use DA_ml                , only : nx, ny
    use DA_ml                , only : nlev
    !use DA_ml                , only : nChemObs

    ! --- in/out ----------------------------

    type(T_ObsOpers), intent(inout)     ::  Hops
    integer, intent(in)                 ::  nobs
    real, intent(in)                    ::  xn_adv_b(:,:,:,:)  ! (nspec_adv,lnx,lny,kmax_mid)
    character(len=*), intent(in)        ::  xn_adv_units(:)    ! (nspec_adv)
    real, intent(in)                    ::  sf_b(:,:  ,:)      ! (lnx,lny     ,nObsComp)
    real, intent(in)                    ::  xn_b(:,:,:,:)      ! (lnx,lny,nlev,nObsComp)
    character(len=*), intent(in)        ::  xn_obs_units(:)    ! (nObsComp)
    integer, intent(in)                 ::  maxobs
    integer, intent(in)                 ::  stnid(maxobs)
    real, intent(inout)                 ::  flon(maxobs)  ! normalized to [-180,180]
    real, intent(inout)                 ::  flat(maxobs)
    real, intent(in)                    ::  falt(maxobs)
    real, intent(in)                    ::  obs(maxobs)
    real, intent(in)                    ::  obs_stddev(maxobs)
    character(len=*), intent(in)        ::  stncodes(maxobs)
    logical, intent(in)                 ::  obs_analyse(maxobs)
    character(len=*), intent(in)        ::  state     ! xf (forecast), xa (analysis)
    integer, intent(out)                ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/H_op'

    ! --- local -----------------------------

    integer      ::  i,j,l,k,n,ilev,err,ierr
    integer      ::  ipar
    real         ::  yn
    !real         ::  alt(nx,ny,nlev)
    !real         ::  rho0

    ! --- begin -----------------------------

    ! info ...
    write (gol,'(a,": nobs = ",i0)') rname, nobs; call goPr

    ! allocate storage, dummy in case nobs is zero:
    call Hops%Alloc( nobs, status )
    IF_NOT_OK_RETURN(status=1)

    ! any obs on this domain ?
    if ( nobs > 0 ) then

      !-----------------------------------------------------------------------
      ! x,y grid coordinates of observations:
      !-----------------------------------------------------------------------
      do n = 1, nobs
        ! index in obsdata array:
        ipar = iObsData(n)
        ! check range ...
        if ( (obs(n) < obsData(ipar)%min) .or. &
             (obs(n) > obsData(ipar)%max)     ) then
          write (gol,'("Observation ",i6," from ",i6," has value ",e16.6," outside accepted value range ",2e16.6)') &
                   n, ipar, obs(n), obsData(ipar)%min, obsData(ipar)%max; call goErr
          TRACEBACK; status=1; return
        end if

        !-----------------------------------------------------------------------
        ! mapping from model to obs-space:
        !-----------------------------------------------------------------------

        ! testing ...
        write (gol,'(a,":   fill ",i0)') rname, n; call goPr

        ! only local field:
        call Hops%obs(n)%Fill( xn_adv_b, xn_adv_units, &
                               sf_b, xn_b, xn_obs_units, &
                               iObsData(n), stnid(n), flat(n),flon(n),falt(n), &
                               yn, status )
        IF_NOT_OK_RETURN(status=1)

        ! store observation:
        Hops%obs(n)%obs = obs(n)

        ! simulation, scale towards observation units:
        select case ( trim(state) )
          case ( 'xf' )
            Hops%obs(n)%xf = yn
          case ( 'xa' )
            Hops%obs(n)%xa = yn
          case default
            write (gol,'("unsupported state `",a,"`")') trim(state); call goErr
            TRACEBACK; status=1; return
        end select

        ! compute innovation:
        Hops%obs(n)%innov = yn - obs(n)

        ! scale obs.error std.dev.:
        Hops%obs(n)%obsstddev = obs_stddev(n)

        ! store station code:
        Hops%obs(n)%stncode = trim(stncodes(n))

        ! testing ...
        write (gol,*) rname//':     flag1 ', size(obs_analyse), n
        write (gol,*) rname//':     flag2 ', obs_analyse(n)
        
        ! store analysis/validation flag:
        Hops%obs(n)%analyse = obs_analyse(n)
        
        ! testing ...
        write (gol,*) rname//':     flag3 '

        ! info ..
        if ( debug .and. MasterProc ) then
          write (gol, '("#",I0,2(1X,A3,":",E12.3))') &
                   n, 'Observation', obs(n), 'Model', yn; call goPr
        end if

      end do  ! observations

    end if  ! nobs > 0

    ! info ...
    write (gol,'(a,": ok")') rname; call goPr

    ! ok
    status = 0

  end subroutine get_innovations


  ! ***


  !-----------------------------------------------------------------------
  ! @description
  ! Compute the costfunction and its gradient
  !
  ! @author M.Kahnert
  !-----------------------------------------------------------------------

  subroutine costFunction( itime, nv_hcr, chi_hcr, &
                            Hops, Bmat, &
                            Jcost, Jb, Jobs, gradJcost_hcr, l2w_hcr, &
                            dx_loc, post, status, &
                            du_loc )

    use MPIF90               , only : MPIF90_AllReduce, MPI_SUM
    use MPIF90               , only : MPIF90_BCast
#ifdef with_ajs
    use DA_Util_ml           , only : GO_Timer_Start, GO_Timer_End
#endif
    use MPI_Groups_ml        , only : MPI_COMM_CALC
    use MPI_Groups_ml        , only : MasterPE
    use ModelConstants_ml    , only : MasterProc, NPROC
    use Par_ml               , only : me
    use Par_ml               , only : limax, ljmax
    use Par_ml               , only : tgi0, tgi1, tgj0, tgj1
    use My_Timing_ml         , only : Add_2timing
    use DA_ml                , only : debug => DEBUG_DA
    use DA_ml                , only : dafmt => da_fmt_msg, damsg => da_msg
    use DA_ml                , only : tim_before => datim_before, tim_after => datim_after
    use DA_Obs_ml            , only : T_ObsOpers
    use DA_Obs_ml            , only : LEVTYPE_3D_ML_SFC, LEVTYPE_3D_ML_TC
    use DA_Obs_ml            , only : LEVTYPE_2D_ML_SFC, LEVTYPE_2D_OBS_SFC

    !-----------------------------------------------------------------------
    ! Formal parameters
    !-----------------------------------------------------------------------

    integer, intent(in)               ::  itime            ! record for hour-of-the-day
    integer, intent(in)               ::  nv_hcr           ! number of elements in half-complex-to-real state
    real, intent(in)                  ::  chi_hcr(nv_hcr)  ! input state
    type(T_ObsOpers), intent(in)      ::  Hops
    type(T_BInfo), intent(inout)      ::  Bmat
    real, intent(out)                 ::  Jcost, Jb, Jobs
    real, intent(out)                 ::  gradJcost_hcr(nv_hcr)
    integer, intent(out)              ::  l2w_hcr(nv_hcr)
    real, intent(out)                 ::  dx_loc(:,:,:,:)  ! (lnx,lny,nlev,ntracer)
    logical, intent(in)               ::  post             ! post processing call ?
    integer, intent(out)              ::  status

    real, intent(out), optional       ::  du_loc(:,:,:,:)  ! (lnx,lny,nlev_all,ntracer)

    !-----------------------------------------------------------------------
    ! parameters
    !-----------------------------------------------------------------------

    character(len=*), parameter  ::  rname = mname//'/costFunction'

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    integer                         ::  nxg, lnyg, iyg0
    integer                         ::  nw_local
    complex(kind=8), allocatable    ::  chi_hc(:)      ! (nw_local)
    complex(kind=8), allocatable    ::  chi2_hc(:)      ! (nw_local)
    integer, allocatable            ::  hcn_local(:)   ! (nw_local)

    integer                         ::  i, j, k
    integer                         ::  n
    integer                         ::  l0, l1
    real, allocatable               ::  yn(:), dep(:), OinvDep(:)
    real, allocatable               ::  HT_OinvDep_loc(:,:,:,:)  ! (limax,ljmax,nlev,ntracer)
    integer                         ::  itracer

    real                            ::  Jb_loc
    real                            ::  Jobs_loc
    real                            ::  gradJb_hcr(nv_hcr)
    real                            ::  gradJobs_hcr(nv_hcr)
    real                            ::  gradnorm, gradnorm_loc

    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------

    !! info ...
    !write (gol,'(a,":   evaluate costFunction ...")') rname; call goPr

    ! init to avoid errors on undefined output:
    dx_loc = 0.0

    ! add contribution to timing:
    call Add_2timing(42,tim_after,tim_before,'3DVar: Optimization.')

    !-----------------------------------------------------------------------
    ! Copy chi (in compact storage format) into chi_arr (in matrix format),
    ! taking into account the relation chi_arr(i,m,n)=conjg(chi_arr(i,-m,-n))
    !-----------------------------------------------------------------------

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_chi2x, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! complex state array
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! local grid size (standard model domain, Bsqrt handles extension):
    call Bmat%BcovarSqrt%Get( status, nlon=nxg, nlat_local=lnyg, ilat_offset=iyg0 )
    IF_NOT_OK_RETURN(status=1)

    ! local y-slab defined ?
    if ( lnyg > 0 ) then
      ! check ...
      if ( any( shape(dx_loc) /= (/nxg,lnyg,Bmat%nlev,Bmat%ntracer/) ) ) then
        write (gol,'("unexpected input shape:")'); call goErr
        write (gol,'("  nxg,lnyg,nlev,ntracer   : ",4i6)') nxg,lnyg,Bmat%nlev,Bmat%ntracer; call goErr
        write (gol,'("  dx_loc                  : ",4i6)') shape(dx_loc); call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! size of current complexe state:
    call Bmat%BcovarSqrt%TimeGet( itime, status, obs_nw_local=nw_local )
    IF_NOT_OK_RETURN(status=1)

    ! check ..
    if ( nv_hcr /= (2*nw_local) ) then
      write (gol,'("number of real values in state (",i0,")")') nv_hcr; call goErr
      write (gol,'("is expected to be twice the number of complex values (",i0,")")') nw_local; call goErr
      TRACEBACK; status=1; return
    end if

    ! local slab defined?
    if ( nw_local > 0 ) then

      ! storage:
      allocate( chi_hc(nw_local), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( chi2_hc(nw_local), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! copy from half-complex-real into half-complex,
      ! loop over complex elements:
      do k = 1, nw_local
        ! copy pair as real and imag part:
        chi_hc(k) = cmplx( chi_hcr(2*k-1), chi_hcr(2*k) )
        chi2_hc(k) = cmplx( chi_hcr(2*k-1), chi_hcr(2*k) )
      end do ! complex elements

      ! storage:
      allocate( hcn_local(nw_local), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! number of values represented in half-complex storage:
      call Bmat%BcovarSqrt%TimeGet( itime, status, obs_hcn_local=hcn_local )
      IF_NOT_OK_RETURN(status=1)

      ! loop over complex elements:
      do k = 1, nw_local
        l2w_hcr(2*k-1:2*k) = hcn_local(k)
      end do
      ! clear:
      deallocate( hcn_local, stat=status )
      IF_NOT_OK_RETURN(status=1)

    else

      ! dummy ...
      allocate( chi_hc(1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( chi2_hc(1), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! set to zero for safety:
      chi_hc = cmplx(0e0,0e0)
      chi2_hc = cmplx(0e0,0e0)

    end if


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! conversion from spectral to model space
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !! info ...
    !write (gol,'(a,":   evaluate dx = B^{1/2} w ...")') rname; call goPr
    
    ! transform:
    call Bmat%BcovarSqrt%Reverse_Obs( itime, chi_hc, dx_loc, status, &
                                       sigma_scale_factors=Bmat%sigma_scale_factors )
    IF_NOT_OK_RETURN(status=1)
    
    !! testing ...
    !write (gol,*) rname, ': xxx dx_loc shape ', shape(dx_loc); call goPr
    !do k = 1, size(dx_loc,3)
    !  write (gol,*) rname, ': xxx dx_loc lev ', k, sum(abs(dx_loc(:,:,k,:)))/size(dx_loc(:,:,k,:)); call goPr
    !end do

    ! also unobserved levels ?
    if ( present(du_loc) ) then
      ! local y-slab defined ?
      if ( lnyg > 0 ) then
        ! check size ...
        if ( any( shape(du_loc) /= (/nxg,lnyg,Bmat%nlev_all,Bmat%ntracer/) ) ) then
          write (gol,'("unexpected input shape:")'); call goErr
          write (gol,'("  nxg,lnyg,nlev,ntracer   : ",4i6)') nxg,lnyg,Bmat%nlev_all,Bmat%ntracer; call goErr
          write (gol,'("  du_loc                  : ",4i6)') shape(du_loc); call goErr
          TRACEBACK; status=1; return
        end if
      end if ! local slab present
      ! other levels are unobserved:
      call Bmat%BcovarSqrt%Reverse_UnObs( itime, chi_hc, du_loc, status, &
                                            sigma_scale_factors=Bmat%sigma_scale_factors )
      IF_NOT_OK_RETURN(status=1)

    end if


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! end
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef with_ajs
    ! end timing:
    call GO_Timer_End( itim_chi2x, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    !-----------------------------------------------------------------------
    ! update unobserved species:
    !-----------------------------------------------------------------------

    ! postprocessing run ?
    if ( post ) then

      !! info ..
      !write (gol,'(a,": posteriori evaluation")') rname; call goPr

      ! clear:
      call my_deallocate( MasterProc, "Postprocessing call to costFunction; no Chi2 need, so leave now." )
      ! ok
      status=0; return

    endif ! post

    !-----------------------------------------------------------------------
    ! projection from model to observation space
    !-----------------------------------------------------------------------

    ! any local observations?
    if ( Hops%nobs > 0 ) then

      !! check ..
      !if ( size(Hops%obs(1)%H_jac,2) /= Bsqrt%ntracer ) then
      !  write (gol,'("covariance for ",i0," tracers, but H_jac for ",i0)') &
      !           Bsqrt%ntracer, size(Hops%obs(1)%H_jac,2); call goErr
      !  TRACEBACK; status=1; return
      !end if

      ! storage:
      allocate( yn(Hops%nobs), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dep(Hops%nobs), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( OinvDep(Hops%nobs), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over observations:
      do n = 1, Hops%nobs

        ! mapping from observed component to covariance tracer index:
        itracer = Bmat%itracer(Hops%obs(n)%iObsComp)

        ! involved grid cell:
        i = Hops%obs(n)%i
        j = Hops%obs(n)%j

        ! switch:
        select case ( Hops%obs(n)%levtype )

          ! simulate from levels:
          case ( LEVTYPE_3D_ML_SFC, LEVTYPE_3D_ML_TC )
            ! level range:
            l0 = Hops%obs(n)%l(0)
            l1 = Hops%obs(n)%l(1)
            ! check ...
            if ( (l0 < 1) .or. (l1 > Bmat%nlev) ) then
              write (gol,'("level range ",i0,":",i0," out of array bounds 1:",i0)') &
                             l0, l1, Bmat%nlev; call goErr
              TRACEBACK; status=1; return
            end if
            ! simulation for this cell, sum over layers:
            yn(n) = sum( Hops%obs(n)%H_jac(l0:l1) * dx_loc(i,j,l0:l1,itracer) )

          ! surface
          case ( LEVTYPE_2D_ML_SFC, LEVTYPE_2D_OBS_SFC )
            ! check ...
            if ( Bmat%nlev /= 1 ) then
              write (gol,'("expected 2D fields for levtype sfc")'); call goErr
              TRACEBACK; status=1; return
            end if
            ! simulation for this cell:
            yn(n) = Hops%obs(n)%H_jac(1) * dx_loc(i,j,1,itracer)

          ! unknown ...
          case default
            write (gol,'("unsupported levtype ",i0)') Hops%obs(n)%levtype; call goErr
            TRACEBACK; status=1; return
        end select
        
        !! testing ...
        !write (gol,*) rname//':    yyy1 H_jac*dx ', n, yn(n); call goPr

      end do  ! n (obs index)

    end if  ! nobs > 0

    !-----------------------------------------------------------------------
    ! adjoint forcing:  H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
    !-----------------------------------------------------------------------

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! departures: O^{-1} * [ H(xb)+H_jac*dx-y ]
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !! info ..
    !write (gol,'(a,":   compute departures ...")') rname; call goPr

    ! any local observations ?
    if ( Hops%nobs > 0 ) then

      ! loop over observations:
      do n = 1, Hops%nobs
        ! to be analysed?
        if ( Hops%obs(n)%analyse ) then
          ! departures:
          !  dep =  [h(xb) + H_jac*dx] - obs
          !      =        h(xb) - obs      +   H_jac*dx
          dep(n) =     Hops%obs(n)%innov   +    yn(n)
        else
          ! set to zero, then no impact in J_obs and grad_J_obs:
          dep(n) = 0.0
        end if
        ! O^{-1} * [ H(xb)+H_jac*dx - obs ]
        !          = [ H(xb)+H_jac*dx - obs ] /       sigma_obs**2
        OinvDep(n) =         dep(n)           / (Hops%obs(n)%obsstddev**2)
        
        !! testing ...
        !write (gol,*) rname//':    yyy1 dep ', n, dep(n), OinvDep(n); call goPr

      end do  ! n (obs)

    end if


    !-----------------------------------------------------------------------
    ! adjoint forcing:  H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
    !-----------------------------------------------------------------------

    !! info ..
    !write (gol,'(a,":   compute forcing ...")') rname; call goPr

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_innov_adj, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! local y-slab defined ?
    if ( lnyg > 0 ) then
      ! storage on extended domain, local y-slab:
      allocate( HT_OinvDep_loc(nxg,lnyg,Bmat%nlev,Bmat%ntracer), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! dummy:
      allocate( HT_OinvDep_loc(1,1,1,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! init to zero:
    HT_OinvDep_loc = 0.0

    ! loop over observations ;
    ! result remains zero if locally no observations present:
    do n = 1, Hops%nobs
      ! grid cell indices:
      i = Hops%obs(n)%i
      j = Hops%obs(n)%j

      ! mapping from observed component to covariance tracer index:
      itracer = Bmat%itracer(Hops%obs(n)%iObsComp)

      ! switch:
      select case ( Hops%obs(n)%levtype )

        ! model levels:
        case ( LEVTYPE_3D_ML_SFC, LEVTYPE_3D_ML_TC )
          ! level range:
          l0 = Hops%obs(n)%l(0)
          l1 = Hops%obs(n)%l(1)
          ! check ...
          if ( (l0 < 1) .or. (l1 > Bmat%nlev) ) then
            write (gol,'("level range ",i0,":",i0," out of array bounds 1:",i0)') &
                           l0, l1, Bmat%nlev; call goErr
            TRACEBACK; status=1; return
          end if
          ! fill departure in level range:
          HT_OinvDep_loc(i,j,l0:l1,itracer) = HT_OinvDep_loc(i,j,l0:l1,itracer) + Hops%obs(n)%H_jac(l0:l1) * OinvDep(n)
          !! testing ...
          !write (gol,'("HT_OinvDep_loc n=",i0," p=",i0," l0:l1=",i0,":",i0," ichem=",i0," H_jac=",es12.4," OinvDep=",es12.4)') &
          !         n, p, l0, l1, ichem, Hops%obs(n)%H_jac(p,l0:l1,ichem), OinvDep(n); call goPr

        ! surface
        case ( LEVTYPE_2D_ML_SFC, LEVTYPE_2D_OBS_SFC )
          ! check ...
          if ( Bmat%nlev /= 1 ) then
            write (gol,'("expected 2D fields for levtype sfc")'); call goErr
            TRACEBACK; status=1; return
          end if
          ! fill departure:
          HT_OinvDep_loc(i,j,1,itracer) = HT_OinvDep_loc(i,j,1,itracer) + Hops%obs(n)%H_jac(1) * OinvDep(n)

        ! unknown ...
        case default
          write (gol,'("unsupported levtype ",i0)') Hops%obs(n)%levtype; call goErr
          TRACEBACK; status=1; return
      end select

    end do ! n obs

#ifdef with_ajs
    ! end timing:
    call GO_Timer_End( itim_innov_adj, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    !-----------------------------------------------------------------------
    ! grad J = grad Jb + grad Jobs
    !   grad Jb   = Chi
    !   grad Jobs = U^{-+}*H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ],
    !-----------------------------------------------------------------------

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_gradient, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! grad Jb
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !! info ..
    !write (gol,'(a,":   compute gradJb ...")') rname; call goPr

    ! gradient to elements of chi:
    gradJb_hcr = chi_hcr

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! grad Jobs
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !! info ..
    !write (gol,'(a,":   HT_OinvDep_loc : ",2es12.4)') rname, minval(HT_OinvDep_loc), maxval(HT_OinvDep_loc); call goPr

    !! info ...
    !write (gol,'(a,":   evaluate w = B^{H/2} (H^T R^{-1} dy) ...")') rname; call goPr

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_x2chi, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! tranform for observed levels/tracers:
    call Bmat%BcovarSqrt%Forward_Obs( itime, HT_OinvDep_loc, chi_hc, status, &
                                         sigma_scale_factors=Bmat%sigma_scale_factors )
    IF_NOT_OK_RETURN(status=1)

    !! info ..
    !write (gol,'(a,":   chi_hc         : ",2es12.4)') rname, minval(abs(chi_hc)), maxval(abs(chi_hc)); call goPr

#ifdef with_ajs
    ! end timing:
    call GO_Timer_End( itim_x2chi, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! any local elements?
    if ( nw_local > 0 ) then
      ! loop over complex elements:
      do k = 1, nw_local
        ! add gradients:
        gradJobs_hcr(2*k-1) =  real(chi_hc(k))
        gradJobs_hcr(2*k  ) = aimag(chi_hc(k))
      end do
    end if  ! nw_local > 0

    !! info ..
    !write (gol,'(a,":   gradJobs_hcr   : ",2es12.4)') rname, minval(gradJobs_hcr), maxval(gradJobs_hcr); call goPr

    ! combine:
    gradJcost_hcr = gradJb_hcr + gradJobs_hcr

    !! info ..
    !write (gol,'(a,":     sum l2w : ",i0)') rname, sum(l2w_hcr); call goPr

    ! info: local l2 sum:
    gradnorm_loc = sum( l2w_hcr * gradJb_hcr**2 )
    ! collect global sum:
    call MPIF90_AllReduce( gradnorm_loc, gradnorm, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! convert to norm:
    gradnorm = sqrt(gradnorm)
    !! info ..
    !write (gol,'(a,":     norm(gradJb  ) : ",e16.6)') rname, gradnorm; call goPr

    ! info: local l2 sum:
    gradnorm_loc = sum( l2w_hcr * gradJobs_hcr**2 )
    ! collect global sum:
    call MPIF90_AllReduce( gradnorm_loc, gradnorm, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! convert to norm:
    gradnorm = sqrt(gradnorm)
    !! info ..
    !write (gol,'(a,":     norm(gradJobs) : ",e16.6)') rname, gradnorm; call goPr

    ! info: local l2 sum:
    gradnorm_loc = sum( l2w_hcr * gradJcost_hcr**2 )
    ! collect global sum:
    call MPIF90_AllReduce( gradnorm_loc, gradnorm, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! convert to norm:
    gradnorm = sqrt(gradnorm)
    !! info ..
    !write (gol,'(a,":     norm(gradJ   ) : ",e16.6)') rname, gradnorm; call goPr

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! end
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef with_ajs
    ! end timing:
    call GO_Timer_End( itim_gradient, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    !-----------------------------------------------------------------------
    ! J = Jb + Jobs
    !  Jb = 0.5 * \chi^{\dagger}*\chi
    !  Jobs = ... (see above)
    !-----------------------------------------------------------------------

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Jb
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !! info ..
    !write (gol,'(a,":   compute Jb ...")')  rname; call goPr

    ! norm over local elements of chi (complex!)
    ! which is the same as the norm over the half-complex-to-real elements
    ! taking into account weights for the redundant values:
    Jb_loc = 0.5 * sum( l2w_hcr * chi_hcr**2 )
    ! total over all domains:
    call MPIF90_AllReduce( Jb_loc, Jb, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)

    !! info ...
    !write (gol,'(a,":     Jb       = ",e16.6)') rname, Jb; call goPr


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Jobs = 0.5 * [ H(xb)+H_jac*dx-y ]^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
    !      = 0.5 *          dep             *         OinvDep
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !! info ..
    !write (gol,'(a,":   compute Jobs ...")') rname; call goPr

    ! any local observations ?
    if ( Hops%nobs > 0 ) then
      ! cost of local observations:
      Jobs_loc = 0.5 * dot_product(dep,OinvDep)
    else
      ! no local costs ...
      Jobs_loc = 0.0
    end if  ! nobs > 0
    ! total over all domains:
    call MPIF90_AllReduce( Jobs_loc, Jobs, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)

    !! info ...
    !write (gol,'(a,":     Jobs     = ",e16.6)') rname, Jobs; call goPr


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! J = Jb + Jobs
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! combine:
    Jcost = Jb + Jobs

    !! info ...
    !write (gol,'(a,":     Jcost    = ",e16.6)') rname, Jcost; call goPr


    !-----------------------------------------------------------------------
    ! end
    !-----------------------------------------------------------------------

    ! stop timing:
    call Add_2timing(43,tim_after,tim_before,'3DVar: costFunction.')

    ! done:
    call my_deallocate(.false.,"NoMessage")

    ! ok
    status = 0

  !-----------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------

    subroutine my_deallocate(verb,msg)
      logical, intent(in)           ::  verb
      character(len=*), intent(in)  ::  msg
      if (verb) print dafmt,msg
      if ( allocated(chi_hc        ) ) deallocate( chi_hc         )
      if ( allocated(hcn_local     ) ) deallocate( hcn_local      )
      if ( allocated(yn            ) ) deallocate( yn             )
      if ( allocated(dep           ) ) deallocate( dep            )
      if ( allocated(OinvDep       ) ) deallocate( OinvDep        )
      if ( allocated(HT_OinvDep_loc) ) deallocate( HT_OinvDep_loc )
    end subroutine my_deallocate

  end subroutine costFunction

end module DA_3DVar_ml
