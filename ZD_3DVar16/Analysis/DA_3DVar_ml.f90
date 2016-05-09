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

  use GO             , only : gol, goPr, goErr
  use TimeDate_ml    , only : date
  use EMEP_BCovarSqrt, only : T_BCovarSqrt
  use GO             , only : T_Domains

  implicit none
  
  ! --- in/out -----------------------------------------
  
  private
  
  public  ::  NTIMING_3DVAR, T_3DVAR
  
  public  ::  DA_3DVar_Init, DA_3DVar_Done
  public  ::  Main_3DVar

  
  ! --- const -----------------------------------------

  character(len=*), parameter  ::  mname = 'DA_3DVar_ml'

  ! settings used by 3D-var:
  character(len=*), parameter   ::  rcfile = 'emo.rc'

  ! timing parameters:
  integer, parameter  ::  NTIMING_3DVAR = 0
  integer, parameter  ::  T_3DVAR = 0
  
  ! debug flags:
  integer, parameter    :: debug_n = 1
  integer, parameter    :: debug_p = 1
  integer, parameter    :: debug_k = 1
  

  ! --- local ----------------------------------------
  
  !! Number of analysis to perform
  !integer, parameter     ::  ANALYSIS_NDATE = 4
  !! When to perform the analysis
  !! By default an analysis will be done every 00,06,12,18 UTC :
  !type(date) :: analysis_date(ANALYSIS_NDATE)= (/ &
  !                 date(-1,-1,-1,00,0), &
  !                 date(-1,-1,-1,06,0), &
  !                 date(-1,-1,-1,12,0), &
  !                 date(-1,-1,-1,18,0)    /)
  ! ... Fill from rcfile settings:
  ! Number of analysis to perform
  integer                   ::  ANALYSIS_NDATE
  ! When to perform the analysis:
  !   (/ date(-1,-1,-1,00,00), date(-1,-1,-1,06,00), ... /)
  type(date), allocatable   ::  analysis_date(:)

  ! link from observation element to original observation dataset:
  integer, allocatable  ::  iObsData(:)
  
  ! Limit dx,du to 500%
  real, parameter       ::  ANALYSIS_RELINC_MAX = 5.0

  !! interpolation type, see 'DA_Obs_ml/H_op' :
  !character(len=16)     ::  H_op_interp
  
  ! m1qn3 verbose ?
  logical               ::  m1qn3_log_all
  ! loging unit:
  integer               ::  m1qn3_io
  ! summary table:
  integer               ::  m1qn3_table
  ! counter:
  integer               ::  m1qn3_number
  ! settings:
  integer               ::  m1qn3_scaling   ! 0=DIS, 1=SIS
  ! number of updates stored in work array:
  integer               ::  m1qn3_nupdates
  
  ! counters:
  integer               ::  m1qn3_maxiter
  integer               ::  m1qn3_maxsim
  
#ifdef with_ajs
  ! timers:
  integer               ::  itim_read_obs, itim_innov, itim_swap
  integer               ::  itim_loop, itim_optimizer, itim_costfunc
  integer               ::  itim_chi2x, itim_innov_adj, itim_gradient, itim_x2chi
#endif
  
  ! info on assim with single covar:
  type T_BInfo
    ! short description
    character(len=32)                 ::  name
    ! number of tracers involved:
    integer                           ::  ntracer
    ! tracer names:
    character(len=32), allocatable    ::  tracers(:)  ! (ntracer)
    ! mapping from tracer to iChemObs:
    integer, allocatable              ::  iChemObs(:)  ! (ntracer)
    ! mapping from iChemObs to itracer,
    ! could be undefined:
    integer, allocatable              ::  itracer(:)  ! (nChemObs)
    ! input file
    character(len=1024)               ::  filename
    ! covariance square root:
    type(T_BCovarSqrt)                ::  BCovarSqrt
  end type T_BInfo
  
  ! collection of B matrices and associated tracers:
  integer                             ::  nBmat
  type(T_BInfo), allocatable          ::  Bmat(:)  ! (nBmat)
  
  !! adhoc factor:
  !real                             ::  BCovarSqrt_sigma_factor

  ! local domain definitions:
  type(T_Domains)       ::  doms_adv_m   ! (s,x,y,z) all tracers, model xy decomposition
  type(T_Domains)       ::  doms_an_m    ! (x,y,z,s) analysis tracers, model xy decomposition
  type(T_Domains)       ::  doms_an_fg   ! (x,y,z,s) analysis tracers, fft y decomposition
  !type(T_Domains)       ::  doms_an_fs   ! (m,nloc,nv1) complex coeff, local slab


contains


  !-----------------------------------------------------------------------
  ! module init/done
  !-----------------------------------------------------------------------

  
  subroutine DA_3DVar_Init( status )

#ifdef with_ajs
#else
    use GO               , only : GO_Init
    use GO               , only : GO_Par_Setup
    use GO               , only : GO_Print_Set
    use GO               , only : GO_Timer_Def
#endif
    use GO               , only : TrcFile, Init, Done, ReadRc
    use GO               , only : goGetFU, goStdErr
    use GO               , only : me
    use ModelConstants_ml, only : masterProc
    use MPI_Groups_ml    , only : MPI_COMM_CALC

    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/DA_3DVar_Init'
    
    ! --- local -----------------------------
    
    character(len=1024)     ::  fname
    character(len=32)       ::  key
    type(TRcFile)           ::  rcF
    integer                 ::  analysis_dhour
    integer                 ::  idate
    integer                 ::  hour
    
    ! --- begin -----------------------------
    
    ! check ...
#ifndef _MPI
    write (gol,'("MPI code should be enabled, define _MPI macro!")'); call goErr
    TRACEBACK; status=1; return
#endif
  
#ifdef with_ajs
#else
    ! initialize GO tools:
    call GO_Init( status )
    if ( status /= 0 ) then
      write (*,'("in ",a," (line",i5,")")') __FILE__, __LINE__
      stop
    end if
    ! from now on, the gol/goPr/goPr logging could be used ...

    ! log from root only:
    call GO_Print_Set( status, apply=MasterProc )
    IF_NOT_OK_RETURN(status=1)
  
    ! setup parallel tools in GO modules:
    call GO_Par_Setup( MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
#endif
  
    ! extra settings:
    call Init( rcF, rcfile, status )
    IF_NOT_OK_RETURN(status=1)
     
    ! init counter for number of times that optimzation loop is performed:
    m1qn3_number = 0
    
    ! should all m1qn3 process shout info ?
    call ReadRc( rcF, 'emo.3dvar.m1qn3.log.all', m1qn3_log_all, status )
    IF_NOT_OK_RETURN(status=1)
    ! m1qn3 messages: only on root or on all ?
    if ( masterProc .or. m1qn3_log_all ) then
      ! log file for m1qn3 internal messages:
      call ReadRc( rcF, 'emo.3dvar.m1qn3.log.file', fname, status )
      IF_NOT_OK_RETURN(status=1)
      ! extend with processor number:
      if ( me > 0 ) write (fname,'(a,".",i0)') trim(fname), me
      ! select free file unit:
      call goGetFU( m1qn3_io, status )
      IF_NOT_OK_RETURN(status=1)
      ! open:
      open( unit=m1qn3_io, file=trim(fname), form='formatted', iostat=status )
      if ( status /= 0 ) then
        write (gol,'("could not open m1qn3 log file (directory not present?) : ",a)') trim(fname); call goErr
        TRACEBACK; status=1; return
      end if      
    else    
      ! m1qn3 is not supposed to print messages on other processes
      ! since below 'impres' is set to '0' for these ;
      ! if it tries to print anyway, write to std.error :
      m1qn3_io = goStdErr
    end if
    
    ! table with norms and costs will be written by root:
    if ( masterProc ) then
      ! csv file for m1qn3 results:
      call ReadRc( rcF, 'emo.3dvar.m1qn3.table', fname, status )
      IF_NOT_OK_RETURN(status=1)
      ! select free file unit:
      call goGetFU( m1qn3_table, status )
      IF_NOT_OK_RETURN(status=1)
      ! open:
      open( unit=m1qn3_table, file=trim(fname), form='formatted', iostat=status )
      if ( status /= 0 ) then
        write (gol,'("could not open m1qn3 table file (directory not present?) : ",a)') trim(fname); call goErr
        TRACEBACK; status=1; return
      end if
      ! header:
      write (m1qn3_table,'("number,step,J,Jb,Jo,gn")')
    end if
  
    ! scaling:
    call ReadRc( rcF, 'emo.3dvar.m1qn3.scaling', key, status )
    IF_NOT_OK_RETURN(status=1)
    ! set integer flag:
    select case ( trim(key) )
      case ( 'DIS' ) ; m1qn3_scaling = 0
      case ( 'SIS' ) ; m1qn3_scaling = 1
      case default
        write (gol,'("unsupported scaling `",a,"`")') trim(key); call goErr
        TRACEBACK; status=1; return
    end select
  
    ! maximum number of iterations:
    call ReadRc( rcF, 'emo.3dvar.m1qn3.maxiter', m1qn3_maxiter, status )
    IF_NOT_OK_RETURN(status=1)
    ! same for maximum number of simulations:
    m1qn3_maxsim = m1qn3_maxiter
    
    ! number of updates stored in memory:
    call ReadRc( rcF, 'emo.3dvar.m1qn3.nupdates', m1qn3_nupdates, status )
    IF_NOT_OK_RETURN(status=1)

#ifdef with_ajs
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
#endif
    
    ! time step between analyses:
    call ReadRc( rcF, 'emo.3dvar.analysis.dhour', analysis_dhour, status )
    IF_NOT_OK_RETURN(status=1)
    ! check ..
    if ( modulo(24,analysis_dhour) /= 0 ) then
      write (gol,'("number of analysis times per day should be integer number;")'); call goErr
      write (gol,'("requested time step in hours between analysis: ",i0)') analysis_dhour; call goErr
      TRACEBACK; status=1; return
    end if
    ! number of times per day:
    ANALYSIS_NDATE = 24 / analysis_dhour
    ! storage for dates:
    allocate( analysis_date(ANALYSIS_NDATE), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over analysis times:
    do idate = 1, ANALYSIS_NDATE
      ! hour from 00:00 with step dhour:
      hour = (idate-1)*analysis_dhour
      ! fill:
      analysis_date(idate) = date(-1,-1,-1,hour,00)
    end do
  
    ! done with settings:
    call Done( rcF, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine DA_3DVar_Init
  
  
  ! ***
  
  
  subroutine DA_3DVar_Done( status )

#ifdef with_ajs
#else
    use GO               , only : GO_Done
#endif
    use ModelConstants_ml, only : MasterProc

    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/DA_3DVar_Done'
    
    ! --- local ----------------------------
    
    integer     ::  iB
    
    ! --- begin -----------------------------

    ! loop over covar matrices:
    do iB = 1, nBmat
      ! clear:
      deallocate( Bmat(iB)%tracers, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( Bmat(iB)%iChemObs, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( Bmat(iB)%itracer, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! done with covariance structure:
      call Bmat(iB)%BCovarSqrt%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! clear:
    deallocate( Bmat, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! only on root ...
    if ( masterProc ) then

      ! close data file:
      close( unit=m1qn3_table, iostat=status )
      IF_NOT_OK_RETURN(status=1)

      ! close logfile:
      close( unit=m1qn3_io, iostat=status )
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

    ! clear:
    deallocate( analysis_date, stat=status )
    IF_NOT_OK_RETURN(status=1)
  
#ifdef with_ajs
#else
    ! done with GO modules:
    call GO_Done( status )
    if ( status /= 0 ) then
      write (*,'("in ",a," (line",i5,")")') __FILE__, __LINE__
      stop
    end if
#endif

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
    use TimeDate_ExtraUtil_ml, only : date2string
    use Util_ml              , only : compare_date
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
    if(first_call) then
      ! info ...
      if ( debug .and. MasterProc) print dafmt,'Initialisation'
      ! initialize:
      call Init_3DVar( status )
      IF_NOT_OK_RETURN(status=1)
      ! reset flag:
      first_call=.false.
    endif
    
    ! info ..
    !if(debug.and.MasterProc)print dafmt,'Test analysis'
    
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
  
    use GO                   , only : TrcFile, Init, Done, ReadRc
    use MPIF90               , only : MPIF90_BCast
    
    use MPI_Groups_ml        , only : MPI_COMM_CALC
    use ModelConstants_ml    , only : MasterProc, nproc
    use ModelConstants_ml    , only : RUNDOMAIN
    use ModelConstants_ml    , only : KMAX_MID, KMAX_BND, KCHEMTOP
    use MPI_Groups_ml        , only : MasterPE
    use Par_ml               , only : me
    use Par_ml               , only : tlimax, tgi0, tgi1, tljmax, tgj0, tgj1
    use CheckStop_ml         , only : CheckStop
    use ChemChemicals_ml     , only : species          ! Gives names
    use ChemGroups_ml        , only : chemgroups       ! group  names
    use ChemSpecs_adv_ml     , only : NSPEC_ADV
    use ChemSpecs_shl_ml     , only : NSPEC_SHL        ! Maps indices
    use SmallUtils_ml        , only : find_index
    use DA_ml                , only : dafmt     => da_fmt_msg
    use DA_ml                , only : debug => DEBUG_DA
    use DA_Obs_ml            , only : varName, varSpec, varSpecInv
    use DA_Obs_ml            , only : obsVarName, observedVar
    use DA_Obs_ml            , only : OBSERVATIONS
    use DA_ml                , only : nlev
    use DA_ml                , only : nchem, nchemObs
    use DA_ml                , only : iChemObs, iChemInv
    use DA_ml                , only : FGSCALE, FGSCALE_INV

    ! --- in/out ----------------------------------
    
    integer, intent(out)    ::  status

    ! --- const ------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Init_3DVar'

    ! --- local ----------------------------------
    
    type(TRcFile)         ::  rcF
    character(len=1024)   :: namelistfile
    integer(4)            :: inNml = 172
    integer               :: nvar, k
    integer               :: ny_local, iy_offset
    
    ! dummy for namelist, actual variables are
    ! 'm1qn3_maxiter' and 'm1qn3_maxsim' read from rcfile:
    integer               :: maxiter, maxsim
    
    integer               ::  nx, ny
    integer               ::  iobs, inoobs
    integer               ::  itracer
    character(len=16)     ::  tracer_name
    character(len=16)     ::  tracer_units
    character(len=16)     ::  species_units
    integer               ::  ivar

    integer                           ::  nB
    character(len=32), allocatable    ::  Bname(:)       ! (nChemObs)
    integer, allocatable              ::  tracer2B(:)    ! (nChemObs)
    integer                           ::  itracer2B
    integer                           ::  iB
    character(len=32)                 ::  key

    ! --- namelists ------------------------------

    namelist /DA_CONFIG/ analysis_date, nChem, nChemObs,&
                         varName, obsVarName, observedVar,&
                         maxiter, maxsim

    ! --- begin ------------------------------
    
    ! info ...
    write (gol,'(a,": initalize 3D var variables ...")') rname; call goPr

    !-----------------------------------------------------------------------
    ! Read config: variable names (observed & unobserved)
    !-----------------------------------------------------------------------

    ! input file:
    namelistfile = 'namelist.nml'
    ! open:
    open( unit=inNml, file=trim(namelistfile), status='OLD', action='READ', &
            form='FORMATTED',iostat=status)
    if ( status /= 0 ) then
      write (gol,'("openening namelist file: ",a)') trim(namelistfile ); call goErr
      TRACEBACK; status=1; return
    end if
    
    !+------------------------------------------------------------------
    ! observed species
    !+------------------------------------------------------------------

    ! find index of 'DAOBS' in 'chemgroups' array:
    k = find_index("DAOBS",chemgroups(:)%name)
    call CheckStop(k<1,'DA group not found: "DAOBS".')
    ! count:
    nChemObs = size(chemgroups(k)%ptr)
    ! copy names from this chemgroup:
    obsVarName(1:nChemObs) = species(chemgroups(k)%ptr(:))%name
    ! copy:
    varName(1:nChemObs) = obsVarName(1:nChemObs)

    ! find index of 'DAUNOBS' in 'chemgroups' array:
    k=find_index("DAUNOBS",chemgroups(:)%name)
    if ( k > 0 ) then
      nChem = nChemObs+size(chemgroups(k)%ptr(:))
      varName(nChemObs+1:nChem) = species(chemgroups(k)%ptr(:))%name
    else
      nChem = nChemObs
      if ( MasterProc ) print dafmt,'WARNING: DA group not found: "DAUNOBS".'
    end if

    ! info ...
    write (gol,'(a,":   settings from DA_CONFIG namelist:")') rname; call goPr
    write (gol,'(a,":     total obs/unobs species (nChem)    : ",i6)') rname, nChem; call goPr
    write (gol,'(a,":     observed species        (nChemObs) : ",i6)') rname, nChemObs; call goPr
    do k = 1, nChemObs
      write (gol,'(a,":       ",i2," ",a)') rname, k, trim(varName(k)); call goPr
    end do
    write (gol,'(a,":     unobserved species                 : ",i6)') rname, nChem-nChemObs; call goPr
    if ( nChemObs < nChem ) then
      do k = 1, nChemObs+1, nChem
        write (gol,'(a,":       ",i2," ",a)') rname, k, trim(varName(k)); call goPr
      end do
    end if
    
    ! set flags:
    observedVar(:)=.false.

    !+------------------------------------------------------------------
    !
    !+------------------------------------------------------------------

    ! read:
    read (unit=inNml,nml=DA_CONFIG,iostat=status)
    if ( status /= 0 ) then
      write (gol,'("reading namelist `DA_CONFIG` from: ",a)') trim(namelistfile); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    call CheckStop(nChemObs>0.eqv.any(observedVar),&
      'Incomplete/Redundant definition of nChemObs and observedVar on DA_CONFIG namelist.')

    !+------------------------------------------------------------------
    !
    !+------------------------------------------------------------------
    
    ! loop over observed variables:
    do nvar=1,nChemObs
      ! check if is in variable list 
      ! (always true, since 'varName' is copied from 'obsVarName' ...)
      k=find_index( obsVarName(nvar), varName(:nChem) )
      call CheckStop(k<1,'Unknown observed variable: '//trim(obsVarName(nvar)))
      ! set flag:
      observedVar(k) = .true.
    enddo
    
    ! at least something should be ovserved ...
    call CheckStop(.not.any(observedVar),'No observed variables found (observedVar).')
    
    ! sort obsVarName following varName order
    nChemObs=0 ! count(observedVar)
    do nvar=1,nChem
      if(observedVar(nvar))then
        nChemObs=nChemObs+1
        obsVarName(nChemObs)=varName(nvar)
      endif
    enddo
    call CheckStop(nChemObs<1,'No observed variables found (nChemObs).')

    ! info ...
    if(debug.and.MasterProc) then
      print dafmt,'B matrix description'
      print "(2(A,:,'(',I0,')'))",&
        'Variable: Observed',nChemObs,'/Unobserved',nChem-nChemObs
      do nvar=1,nChem
        if(observedVar(nvar))then
          print "(I4,': ',A10,'=O',I3.3,$)",nvar,trim(varName(nvar)),count(observedVar(:nvar))
        else
          print "(I4,': ',A10,'=U',I3.3,$)",nvar,trim(varName(nvar)),count(.not.observedVar(:nvar))
        endif
        if(mod(nvar,5)==0)print *,''
      enddo
      if(mod(nChem,5)/=0)print *,''
    endif
    
    ! info ...
#ifdef gFortran
    if(debug) write(*,nml=DA_CONFIG)
#else
    if(debug) write(*,nml=DA_CONFIG,delim='QUOTE')
#endif
    
    !-----------------------------------------------------------------------
    ! Read observation parameters
    !-----------------------------------------------------------------------

    !
    ! read settings, namelist defined in 'DA_Obs_ml' :
    !
    !   nobsData                ! number of datasets
    !   obsData(1:nobsData)%..  ! types with dataset properties
    !   interpolation_method    ! 'nearest-neighbor'
    !
    read (unit=inNml,nml=OBSERVATIONS,iostat=status)
    if ( status /= 0 ) then
      write (gol,'("reading namelist `OBSERVATIONS` from: ",a)') trim(namelistfile); call goErr
      TRACEBACK; status=1; return
    end if


    !-----------------------------------------------------------------------
    ! Read covariance matrix
    !-----------------------------------------------------------------------

    ! info ...
    write (gol,'(a,": setup observed/unobserved variables")') rname; call goPr

    ! Set numbers of observed/unobserved variables and
    ! mapping from observed/unobserved index to analysis variables:
    !   nchemObs  , ichemObs  (1:nchemObs  ) -> ivar
    !   nchemNoObs, ichemNoObs(1:nchemNoObs) -> ivar
    !               ichemInv  (1:nvar)       -> iChem (in iChemObs or iChemNoObs)
    ! (Original implementation in "covmat_ml/set_chemobs_idx")
    ! Count:
    nChemObs   = count(observedVar)
    ! storage:
    allocate( ichemObs(nchem), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( ichemInv(nchem), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! init indices:
    iobs   = 0
    ! loop over DA variables:
    do ivar = 1, nChem
      ! observed ?
      if ( observedVar(ivar) ) then
        ! increase counter:
        iobs = iobs + 1
        ! store mapping:
        iChemObs(iobs) = ivar
        iChemInv(ivar) = iobs
      else
        call CheckStop(rname//", unobserved variables are no longer supported")
      end if
    end do  ! DA variables
    
    ! * assign B to tracer(s)

    ! info ...
    write (gol,'(a,": setup B matrices ..")') rname; call goPr
      
    ! extra settings:
    call Init( rcF, rcfile, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init counter:
    nBmat = 0
    ! storage for keywords defining B for individual or groups of tracers:
    allocate( Bname(nChemObs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( tracer2B(nChemObs), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over tracers to be analysed:
    do ivar = 1, nChemObs
      ! read name that identifies B:
      call ReadRc( rcF, 'emo.3dvar.tracer2B.'//trim(varName(ivar)), key, status )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'(a,":   observed tracer ",i0," : ",a)') rname, ivar, trim(varName(ivar)); call goPr
      write (gol,'(a,":     associated B : ",a)') rname, key; call goPr
      ! search in list of already defined keys:
      itracer2B = -999
      do k = 1, nBmat
        ! tracer should use this B?
        if ( trim(key) == trim(Bname(k)) ) then
          ! store index:
          itracer2B = k
          ! info ...
          write (gol,'(a,":     use B matrix nr. ",i0)') rname, itracer2B; call goPr
          ! leave:
          exit
        end if  ! tracer uses this B
      end do  ! B matrices
      ! new?
      if ( itracer2B < 0 ) then
        ! increase counter:
        nBmat = nBmat + 1
        ! store name:
        Bname(nBmat) = trim(key)
        ! assign index:
        itracer2B = nBmat
        ! info ...
        write (gol,'(a,":     assign new B matrix nr. ",i0)') rname, itracer2B; call goPr
      end if
      ! store index in mapping from tracer to covariance:
      tracer2B(ivar) = itracer2B
    end do ! tracers
    
    ! info ..
    write (gol,'(a,":   number of B matrices: ",i0)') rname, nBmat; call goPr

    ! storage for B matrices and associated info:
    allocate( Bmat(nBMat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over covar matrices:
    do iB = 1, nBmat

      ! store name:
      Bmat(iB)%name = trim(Bname(iB))
      ! number of tracers:
      Bmat(iB)%ntracer = count( tracer2B == iB )
      ! storage:
      allocate( Bmat(iB)%tracers(Bmat(iB)%ntracer), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( Bmat(iB)%iChemObs(Bmat(iB)%ntracer), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! mapping from nChemObs to tracer:
      allocate( Bmat(iB)%itracer(nChemObs), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! init as undefined:
      Bmat(iB)%itracer = -999
      ! init counter:
      itracer = 0
      ! loop over tracers:
      do ivar = 1, nChemObs
        ! current?
        if ( tracer2B(ivar) == iB ) then
          ! increase counter:
          itracer = itracer + 1
          ! store name:
          Bmat(iB)%tracers(itracer) = trim(varName(ivar))
          ! store index:
          Bmat(iB)%iChemObs(itracer) = ivar
          ! store reverse index:
          Bmat(iB)%itracer(ivar) = itracer
        end if ! current
      end do  ! observed variables

      ! read filename:
      call ReadRc( rcF, 'emo.3dvar.B.'//trim(Bname(iB))//'.filename', Bmat(iB)%filename, status )
      IF_NOT_OK_RETURN(status=1)

      ! info ...
      write (gol,'(a,":     B matrix ",i0," (",a,")")') rname, iB, trim(Bname(iB)); call goPr
      do itracer = 1, Bmat(iB)%ntracer
        write (gol,'(a,":       tracer ",i2," ",a)') rname, itracer, trim(Bmat(iB)%tracers(itracer)); call goPr
      end do ! tracers
      write (gol,'(a,":       file: ",a)') rname, trim(Bmat(iB)%filename); call goPr

      ! init covariance structure,
      ! setup mapping to the subset of observed variables:
      call Bmat(iB)%BCovarSqrt%Init( trim(Bmat(iB)%filename), status, &
                                       tracers=Bmat(iB)%tracers, &
                                       comm=MPI_COMM_CALC )
      IF_NOT_OK_RETURN(status=1)

      !! info ...
      !write (gol,'(a,":   adhoc scaling factor for sigma : ",f6.1)') &
      !                  rname, BCovarSqrt_sigma_factor; call goPr
      !! apply scale factor:
      !BCovarSqrt%S = BCovarSqrt%S * BCovarSqrt_sigma_factor

    end do ! covars

    ! clear:
    deallocate( Bname, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( tracer2B, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! done with settings:
    call Done( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! *
    
    ! assumed grid shape:
    nx   = RUNDOMAIN(2) - RUNDOMAIN(1) + 1
    ny   = RUNDOMAIN(4) - RUNDOMAIN(3) + 1
    ! number of levels:
    nlev = KMAX_MID  ! 20
    
    ! loop over B matrices:
    do iB = 1, nBmat
    
      ! check ...
      call Bmat(iB)%BCovarSqrt%Check( status, nlon=nx, nlat=ny, nlev=nlev )
      IF_NOT_OK_RETURN(status=1)

      !! get extended grid size:
      !call BCovarSqrt%Get( status, nlon_ex=nxex, nlat_ex=nyex )
      !IF_NOT_OK_RETURN(status=1)

      ! extract info on grid and decomposition used by covariance:
      call Bmat(iB)%BCovarSqrt%Get( status, nlat_ex_local=ny_local, ilat_ex_offset=iy_offset )
      IF_NOT_OK_RETURN(status=1)

      !! check of tracers in B:
      !call BCovarSqrt%Check( status, ntracer=nChem )
      !IF_NOT_OK_RETURN(status=1)

      ! check units; for the moment assume that all model tracers
      ! are in 'mix_ratio', and the covariance should therefore be in 'mol mol^-1' :
      species_units = 'mol mol^-1'
      ! loop over tracers in B:
      do ivar = 1, nChem
        ! mapping to tracer in this covariance:
        itracer = Bmat(iB)%itracer(ivar)
        ! undefined? then not in this covar:
        if ( itracer < 0 ) cycle
        ! get properties:
        call Bmat(iB)%BCovarSqrt%TracerGet( itracer, status, name=tracer_name, units=tracer_units )
        IF_NOT_OK_RETURN(status=1)
        ! compare:
        if ( trim(tracer_units) /= trim(species_units) ) then
          write (gol,'("covariance tracer `",a,"` has units `",a,"` while model has `",a,"`")') &
              trim(tracer_name), trim(tracer_units), trim(species_units); call goErr
          TRACEBACK; status=1; return
        end if
      end do  ! analysis tracers
      
    end do ! B matrices

    ! scaling factor to convert from model units to BCovarSqrt units,
    ! might become tracer depenend in future; 
    ! with the above assumed units it is simply unity:
    FGSCALE     = 1.0
    FGSCALE_INV = 1.0 / FGSCALE

    ! info ...
    write (gol,'(a,": DA variables: ",i0)') rname, nChem; call goPr
    do ivar = 1, nChem
      write (gol,'(a,":   ",i2," ; observed : ",l1," ; inv : ",i0)') &
              rname, ivar, observedVar(ivar), iChemInv(ivar); call goPr
    end do
    write (gol,'(a,": DA Obs variables: ",i0)') rname, nChemObs; call goPr
    do ivar = 1, nChemObs
      write (gol,'(a,":   ",i2," ; ivar : ",i0)') &
              rname, ivar, iChemObs(ivar); call goPr
    end do
    
    !-----------------------------------------------------------------------
    ! mapping from 3D-var variables to model species
    !-----------------------------------------------------------------------

    ! info ...
    write (gol,'(a,":   model species involved in analysis ...")') rname; call goPr

    ! init mapping arrays:
    varSpec    = -1  ! (1:nvar) ->  global tracer index k
    varSpecInv = -1  ! (k)      ->  ivar
    ! loop over specs involved in analysis:
    do nvar = 1, nChem
      ! info ...
      write (gol,'(a,":     analysis var ",i2," `",a,"`")') rname, nvar, trim(varName(nvar)); call goPr
      ! check ...
      if ( count( varName(1:nvar) == varName(nvar) ) > 1 ) then
        write (gol,'(a,": multiple definitions of variable name: ",a)') rname, trim(varName(nvar)); call goErr
        TRACEBACK; status=1; return
      end if
      ! search ...
      k = find_index(varName(nvar),species(NSPEC_SHL+1:)%name)
      ! check ...
      if ( k < 1 ) then
        write (gol,'(a,": Wrong variable name: ",a)') rname, trim(varName(nvar)); call goErr
        TRACEBACK; status=1; return
      end if
      ! info ...
      write (gol,'(a,":       model tracer index: ",i4)') rname, k; call goPr
      ! for DA variable index 'nvar' provide global tracer index 'k' :
      varSpec(nvar) = k
      ! for global tracer index 'k' provide DA variable index 'nvar' :
      varSpecInv(k) = nvar
    end do
    
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
    call doms_an_m%Init( (/ tgi0(me), tgj0(me),    1,     1 /), &
                         (/ tgi1(me), tgj1(me), nlev, nChem /), status )
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
      !                       lon           lat                lev     var
      call doms_an_fg%Init( (/  1,     iy_offset+       1    ,  1  ,     1 /), &
                            (/ nx, min(iy_offset+ny_local,ny), nlev, nChem /), status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! info ...
    write (gol,'(a,": end")') rname; call goPr
    
    ! ok
    status = 0
    
  end subroutine init_3dvar


  ! ***


  !-----------------------------------------------------------------------
  subroutine generic3dvar( cdate, status )
  !-----------------------------------------------------------------------
  ! @description
  ! Generic version of 3d variational analysis.
  ! @author M.Kahnert
  !-----------------------------------------------------------------------
    use MPIF90               , only : MPIF90_AllReduce, MPI_SUM
    use GO                   , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use MPI_Groups_ml        , only : MPI_COMM_CALC
    use ModelConstants_ml    , only : MasterProc, NPROC
    use ModelConstants_ml    , only : RUNDOMAIN
    use ModelConstants_ml    , only : KMAX_MID, KMAX_BND, KCHEMTOP
    use MPI_Groups_ml        , only : MasterPE
    use CheckStop_ml         , only : CheckStop
    use My_Timing_ml         , only : Code_timer, Add_2timing
    use TimeDate_ml          , only : date  ! date/time structure
    use Par_ml               , only : me
    use Par_ml               , only : MAXLIMAX, MAXLJMAX   ! local x, y dimensions
    use Par_ml               , only : tlimax, tgi0, tgi1, tljmax, tgj0, tgj1
    use Par_ml               , only : limax, ljmax
    use GridValues_ml        , only : glon, glat
    use ChemSpecs_adv_ml     , only : NSPEC_ADV
    use MetFields_ml         , only : z_bnd
    use MetFields_ml         , only : roa
    use ChemFields_ml        , only : xn_adv
    use Chemfields_ml        , only : cfac
    use Chemfields_ml        , only : PM25_water
    use Chemfields_ml        , only : PM25_water_rh50
    !use exd_domain_ml        , only : EXT_DOMAIN, EXT_DOMAIN_INV
    use DA_ml                , only : debug => DEBUG_DA
    use DA_ml                , only : dafmt => da_fmt_msg
    use DA_ml                , only : damsg => da_msg
    use DA_ml                , only : tim_before => datim_before, tim_after => datim_after
    use DA_Obs_ml            , only : T_ObsOpers
    use DA_Obs_ml            , only : varName, varSpec
    use DA_Obs_ml            , only : Read_Obs
    use DA_ml                , only : nlev
    use DA_ml                , only : nChem
    use DA_ml                , only : nchemObs  , ichemObs
    use DA_ml                , only : FGSCALE, FGSCALE_INV

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

    integer               ::  i,j,ilev,ilev0,k,kk,ii, maxobs, nvar
    integer               ::  nobs
    integer               ::  nobs_tot
    integer               ::  ichem
    integer, allocatable  ::  stnid(:)        ! (nobs)
    real, allocatable     ::  flat(:)         ! (nobs)
    real, allocatable     ::  flon(:)         ! (nobs)
    real, allocatable     ::  falt(:)         ! (nobs)
    real, allocatable     ::  obs(:)          ! (nobs)
    real, allocatable     ::  obsstddev1(:)   ! (nobs)
    character(len=32), allocatable  ::  stncodes(:)     ! (nobs)
    
    integer               ::  lnx, lny
    logical               ::  has_local_domain

    ! idem operator on model and assimilation decomposition:
    type(T_ObsOpers)      ::  Hops_m
    type(T_ObsOpers)      ::  Hops_f
    type(T_ObsOpers)      ::  Hops_f_B
    integer               ::  nobs_B

    real, allocatable     ::  xn_an (:,:,:,:)     ! (limax,ljmax,nlev,nChem         )
    real, allocatable     ::  xn_loc(:,:,:,:)     ! (lnx  ,lny  ,nlev,nChem         )
    real, allocatable     ::  dx_loc(:,:,:,:)     ! (lnx  ,lny  ,nlev,nChemObs      )
    
    integer               ::  iB
    real, allocatable     ::  dx_loc_B(:,:,:,:)     ! (lnx  ,lny  ,nlev,ntracer)
    integer               ::  itracer
    integer               ::  ichm
    
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

    ! storage per observed grid cell:
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

    ! selection on local grid point domain used by model:
    call read_obs( doms_adv_m, maxobs, &
                     stnid, flat,flon,falt, obs, obsstddev1, stncodes, &
                     iObsData, nobs, status )
    IF_NOT_OK_RETURN(status=1)

    ! info ...
    write (gol,'(a,": local number of observations read : ",i6)') rname, nobs; call goPr

    ! total number of observations over all domains:
    call MPIF90_AllReduce( nobs, nobs_tot, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! info ...
    write (gol,'(a,": total number of observations read : ",i6)') rname, nobs_tot; call goPr

#ifdef with_ajs
    ! end timing:
    call GO_Timer_End( itim_read_obs, status )
    IF_NOT_OK_RETURN(status=1)
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

      ! NOTE: This is only needed for not changing 'H_op' too much,
      !       which extracts concentrations from an array with shape:
      !           (lon,lat,lev,ichemobs)
      !       In future, extract everthing directly from 'xn_adv',
      !       which is already used in 'H_op' for the indirect observations

      ! storage for extract of concentrations involved in analysis:
      allocate( xn_an(limax,ljmax,nlev,nChem), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! loop over variables involved in analysis:
      do nvar = 1, nChem
        ! info ...
        write (gol,'(a,": copy xn_adv spec ",i0," into xn_an var ",i0)') rname, varSpec(nvar), nvar; call goPr
        ! extract slab:
        !     x,y,z,s                     s     ,    x  ,   y   ,          z
        xn_an(:,:,:,nvar) = xn_adv(varSpec(nvar),1:limax,1:ljmax,KMAX_MID-nlev+1:KMAX_MID)
      end do

      ! info ...
      write (gol,'(a,": xn_an range : ",2e16.6)') rname, minval(xn_an), maxval(xn_an); call goPr

      ! scale if necessary:
      if ( FGSCALE /= 1.0e0 ) xn_an = xn_an * FGSCALE

      ! info ...
      write (gol,'(a,": xn_an range : ",2e16.6)') rname, minval(xn_an), maxval(xn_an); call goPr


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
      call get_innovations( Hops_m, nobs, xn_adv, xn_an, &
                              maxobs, stnid, flat,flon,falt, obs, obsstddev1, stncodes, &
                              'xf', status )
      IF_NOT_OK_RETURN(status=1)
      
      !! save all:
      !call Hops_m%WriteToFile( cdate, status )
      !IF_NOT_OK_RETURN(status=1)
      !! break 
      !write (gol,'("break after write obs")'); call goErr
      !TRACEBACK; status=1; return
      
      ! assign unique id's to observations:
      call Hops_m%Set_IDs( status )
      IF_NOT_OK_RETURN(status=1)

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
      ! allocate work arrrays
      !-----------------------------------------------------------------------

      ! local y-slabs ?
      if ( lny > 0 ) then
        ! storage for local concentration increments, new decomposition:
        allocate( xn_loc(lnx,lny,nlev,nChem     ), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( dx_loc(lnx,lny,nlev,nChemObs  ), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! dummy:
        allocate( xn_loc(1,1,1,1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( dx_loc(1,1,1,1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if

      ! dumy values:
      xn_loc = 0.0
      dx_loc = 0.0

      !-----------------------------------------------------------------------
      ! perform variational analysis
      !-----------------------------------------------------------------------

      ! info ...
      write (gol,'(a,": call var3d ...")') rname; call goPr
      
      ! loop over B matrices:
      do iB = 1, nBmat

        ! info ...
        write (gol,'(a,":   B matrix ",i0," (",a,")")') rname, iB, trim(Bmat(iB)%name); call goPr
        
        ! init selected observations:
        call Hops_f_B%Init( status )
        IF_NOT_OK_RETURN(status=1)
        ! fill with copy for tracer selection:
        call Hops_f_B%SelectTracers( Hops_f, Bmat(iB)%itracer, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! total number of observations over all domains:
        call MPIF90_AllReduce( Hops_f_B%nobs, nobs_B, MPI_SUM, MPI_COMM_CALC, status )
        IF_NOT_OK_RETURN(status=1)
        ! info ...
        write (gol,'(a,":     number of observations : ",i6)') rname, nobs_B; call goPr

        ! any?
        if ( nobs_B > 0 ) then
        
          ! only tracers associated with this covar:
          if ( lny > 0 ) then
            ! storage for local concentration increments, new decomposition:
            allocate( dx_loc_B(lnx,lny,nlev,Bmat(iB)%ntracer), stat=status )
            IF_NOT_OK_RETURN(status=1)
            ! copy:
            do itracer = 1, Bmat(iB)%ntracer
              ichm = Bmat(iB)%iChemObs(itracer)
              dx_loc_B(:,:,:,itracer) = dx_loc(:,:,:,ichm)
            end do
          else
            ! dummy:
            allocate( dx_loc_B(1,1,1,1), stat=status )
            IF_NOT_OK_RETURN(status=1)
          end if

          ! perform analysis, all processes involved for operations with H;
          ! only the tracer values assigned to this B matrix are changed:
          call var3d( cdate, Hops_f_B, Bmat(iB)%Bcovarsqrt, &
                        dx_loc_B, status )
          IF_NOT_OK_RETURN(status=1)

          ! restore:
          if ( lny > 0 ) then
            ! copy:
            do itracer = 1, Bmat(iB)%ntracer
              ichm = Bmat(iB)%iChemObs(itracer)
              dx_loc(:,:,:,ichm) = dx_loc_B(:,:,:,itracer)
            end do
          end if

          ! clear:
          deallocate( dx_loc_B, stat=status )
          IF_NOT_OK_RETURN(status=1)
          
        end if  ! nobs_B > 0

        ! done:
        call Hops_f_B%Done( status )
        IF_NOT_OK_RETURN(status=1)

      end do ! B matrices
        
      ! info ...
      write (gol,'(a,": back from var3d ...")') rname; call goPr


      !-----------------------------------------------------------------------
      ! read result for \delta x and \delta u and add to background field;
      !-----------------------------------------------------------------------
      ! Only update whithin the non-extended zone (since there are no obs. in
      ! that zone => zero increment). Leave BCs (lateral & top) unchanged.
      !-----------------------------------------------------------------------

      ! info ...
      write (gol,'(a,": swap xn_an to analysis domains ...")') rname; call goPr

      ! swap analysed concentrations to new decomposition:
      call doms_an_m%Swap( xn_an, doms_an_fg, xn_loc, status )
      IF_NOT_OK_RETURN(status=1)

      ! info ...
      write (gol,'(a,": add analysis increments ...")') rname; call goPr

      ! y-slab defined ?
      if ( lny > 0 ) then

        ! any observed ?
        if ( nChemObs > 0 ) then
          ! loop over tracers:
          do nvar = 1, nChemObs
            ! index:
            ichem = iChemObs(nvar)
            ! add, truncate between factors [0,fac] :
            xn_loc(:,:,:,ichem) = min( max( 0.0, xn_loc(:,:,:,ichem) + dx_loc(:,:,:,nvar) ), &
                                         xn_loc(:,:,:,ichem) * ANALYSIS_RELINC_MAX )
          end do  ! tracers
        end if

        ! rescale if necessary:
        if ( FGSCALE_INV /= 1e0 ) xn_loc = xn_loc * FGSCALE_INV

      end if  ! y-slab

      ! info ...
      write (gol,'(a,": swap back to model decomposition ...")') rname; call goPr

      ! swap to original model decomposition:
      call doms_an_fg%Swap( xn_loc, doms_an_m, xn_an, status )
      IF_NOT_OK_RETURN(status=1)

      ! info ...
      write (gol,'(a,": reset concentrations in model array ...")') rname; call goPr

      ! loop over variables involved in analysis:
      do nvar = 1, nChem
        ! restore in original location:
        !xn_adv(varSpec(nvar),1:limax,1:ljmax,KMAX_MID-nlev+1:KMAX_MID) = xn_an_bm(:,:,:,nvar)
        xn_adv(varSpec(nvar),1:limax,1:ljmax,KMAX_MID-nlev+1:KMAX_MID) = xn_an(:,:,:,nvar)
      end do

      !-----------------------------------------------------------------------
      ! innovations after analysis
      !-----------------------------------------------------------------------

      ! info ...
      write (gol,'(a,": evaluate observations after analysis ...")') rname; call goPr

      ! evaluate observation operator, store in it's 'xa' field:
      call Hops_m%Evaluate( 'xa', xn_an, status )
      IF_NOT_OK_RETURN(status=1)

      ! info ...
      write (gol,'(a,": write innovations ...")') rname; call goPr
      
      ! save all:
      call Hops_m%WriteToFile( cdate, status )
      IF_NOT_OK_RETURN(status=1)


      !-----------------------------------------------------------------------
      ! deallocate analysis arrays
      !-----------------------------------------------------------------------

      ! info ...
      write (gol,'(a,": clear analysis arrays ...")') rname; call goPr

      ! done with obs operators:
      call Hops_m%Done( status )
      IF_NOT_OK_RETURN(status=1)
      call Hops_f%Done( status )
      IF_NOT_OK_RETURN(status=1)

      ! clear:
      deallocate( xn_an, stat=status  )
      IF_NOT_OK_RETURN(status=1)
      deallocate( xn_loc, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dx_loc, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! observations present

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

    ! clear:
    deallocate( iObsData, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0

  end subroutine generic3dvar

  
  ! ***

  
  !-----------------------------------------------------------------------
  subroutine var3d( cdate, Hops, Bsqrt, dx_loc, status )
  !-----------------------------------------------------------------------
  ! @description
  ! 3-D variational analysis
  ! @author M.Kahnert
  !-----------------------------------------------------------------------
  
    use MPIF90               , only : MPIF90_BCast
    use MPIF90               , only : MPIF90_AllReduce, MPI_SUM
    use GO                   , only : GO_Timer_Start, GO_Timer_End, GO_Timer_Switch
    use ModelConstants_ml    , only : MasterProc, NPROC
    use MPI_Groups_ml        , only : MasterPE
    use CheckStop_ml         , only : CheckStop
    use My_Timing_ml         , only : Code_timer, Add_2timing
    use Util_ml              , only : norm
    use TimeDate_ml          , only : date  ! date/time structure
    use Par_ml               , only : me
    use DA_ml                , only : debug => DEBUG_DA
    use DA_ml                , only : dafmt => da_fmt_msg, damsg => da_msg
    use DA_ml                , only : tim_before => datim_before, tim_after => datim_after
    use DA_Obs_ml            , only : T_ObsOpers
    use DA_ml                , only : nlev
    use DA_ml                , only : nChemObs
    
    !-----------------------------------------------------------------------
    ! Formal parameters
    !-----------------------------------------------------------------------

    type(date), intent(in)            ::  cdate
    type(T_ObsOpers), intent(inout)   ::  Hops
    type(T_BCovarSqrt), intent(inout) ::  Bsqrt
    real, intent(out)                 ::  dx_loc(:,:,:,:)  ! (limax,ljmax,nlev,nChemObs)
    integer, intent(out)              ::  status

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
    call Bsqrt%FindTime( cdate%hour, itime, status )
    IF_NOT_OK_RETURN(status=1)

    ! current size of half-complex state,
    ! local as well as global:
    call Bsqrt%TimeGet( itime, status, &
                                nw_local=nw_local, nw=nw )
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
    niter = m1qn3_maxiter
    nsim  = m1qn3_maxsim

    ! size of work array:
    ndz     = 4*nv_hcr     + m1qn3_nupdates*(2*nv_hcr    +1)
    ndz_all = 4*nv_hcr_all + m1qn3_nupdates*(2*nv_hcr_all+1)
    ! storage:
    allocate( dz(ndz), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! settings:
    imode(1) = m1qn3_scaling   ! 0=DIS, 1=SIS
    imode(2) = 0               ! cold start
    imode(3) = 0               ! calculate Jcost,gradJcost every interation

    ! verbosity level:
    if ( masterProc .or. m1qn3_log_all ) then
      ! shout:
      impres = 5    ! Full verbosity
    else
      ! shut up ...
      impres = 0     ! No print
    end if
    
    ! info ...
    if(debug.and.MasterProc.and.impres>0)then
      write(damsg,"(A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
        'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
      print dafmt,'Calling m1qn3 '//trim(damsg)
    endif

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_loop, status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! start timing:
    call Code_timer(tim_before)
    
    ! increase counter:
    m1qn3_number = m1qn3_number + 1

    ! reverse mode, 'costFunction' is called externally;
    ! perform loop over optimer steps:
    istep = 0
    do

      ! increase counter:
      istep = istep + 1
      !! info ...
      !write (gol,'(a,": solver step ",i4," ...")') rname, istep; call goPr
    
#ifdef with_ajs
      ! start timing:
      call GO_Timer_Start( itim_costfunc, status )
      IF_NOT_OK_RETURN(status=1)
#endif

      !! info ...
      !write (gol,'(a,":   evaluate cost function and gradient ...")') rname; call goPr
      ! next evaluation of costfunction and gradient;
      ! for operations with H this involves all processes :
      call costFunction( itime, nv_hcr, chi_hcr, Hops, Bsqrt, &
                           Jcost, Jcost_b, Jcost_obs, gradJcost_hcr, l2w_hcr, &
                           dx_loc, .false., status )
      IF_NOT_OK_RETURN(status=1)

      !! info ...
      !write (gol,'(a,":     Jcost    : ",e16.6)') rname, Jcost; call goPr
      !write (gol,'(a,":     dx range : ",2e16.6)') rname, minval(dx_loc), maxval(dx_loc); call goPr
      
      ! only on root ...
      if ( masterProc ) then
        ! add cost function values to summary:
        write (m1qn3_table,'(i,",",i,4(",",e))',iostat=status) &
                 m1qn3_number, istep, Jcost, Jcost_b, Jcost_obs, norm(gradJcost_hcr)
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

    ! info ...
    write (gol,'(a,":     returned omode   : ",i4)') rname, omode; call goPr
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
    call costFunction( itime, nv_hcr, chi_hcr, Hops, Bsqrt, &
                        Jcost, Jcost_b, Jcost_obs, gradJcost_hcr, l2w_hcr, &
                        dx_loc, .true., status )
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


  subroutine get_innovations( Hops, nobs, xn_adv_b, xn_b, &
                                maxobs, stnid, flat,flon,falt, &
                                obs, obs_stddev, stncodes, state, status )
  !-----------------------------------------------------------------------
  ! @description
  ! Compute innovations H(xb)-y, where xb is
  ! the background field, y denotes the observations, and
  ! H is an operator mapping from model to observation space.
  ! On input, innov contains the observations y. 
  ! On output, innov contains the innovations.
  ! @author M.Kahnert
  !-----------------------------------------------------------------------
    use ModelConstants_ml    , only : MasterProc
    use CheckStop_ml         , only : CheckStop
    use GridValues_ml        , only : coord_in_processor
    use DA_ml                , only : debug => DEBUG_DA
    use DA_Obs_ml            , only : T_ObsOpers
    use DA_Obs_ml            , only : obsData
    use DA_ml                , only : FGSCALE, FGSCALE_INV
    !use DA_ml                , only : nx, ny
    use DA_ml                , only : nlev
    use DA_ml                , only : nChemObs
    
    ! --- in/out ----------------------------

    type(T_ObsOpers), intent(inout)     ::  Hops
    integer, intent(in)                 ::  nobs
    real, intent(in)                    ::  xn_adv_b(:,:,:,:)  ! (nspec_adv,lnx,lny,kmax_mid)
    real, intent(in)                    ::  xn_b(:,:,:,:)      ! (lnx,lny,nlev,nChem)
    integer, intent(in)                 ::  maxobs
    integer, intent(in)                 ::  stnid(maxobs)
    real, intent(inout)                 ::  flon(maxobs)  ! normalized to [-180,180]
    real, intent(inout)                 ::  flat(maxobs)
    real, intent(in)                    ::  falt(maxobs)
    real, intent(in)                    ::  obs(maxobs)
    real, intent(in)                    ::  obs_stddev(maxobs)
    character(len=*), intent(in)        ::  stncodes(maxobs)
    character(len=*), intent(in)        ::  state     ! xf (forecast), xa (analysis)
    integer, intent(out)                ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/H_op'

    ! --- local -----------------------------

    integer      ::  i,j,l,k,n,ilev,err,ierr
    integer      ::  ipar
    real         ::  yn
    !real         ::  alt(nx,ny,nlev)
    real         ::  rho0

    ! --- begin -----------------------------
    
    ! info ...
    write (gol,'(a,": nobs = ",i0)') rname, nobs; call goPr
    
    ! allocate storage, dummy in case nobs is zero:
    call Hops%Alloc( nobs, nChemObs, status )
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

        ! only local field ;
        !  'xn_b' contains scale factor 'FGSCALE' already,
        !  simulation 'yn' is therefore including this factor too:
        call Hops%obs(n)%Fill( xn_adv_b, xn_b, &
                               iObsData(n), stnid(n), flat(n),flon(n),falt(n), &
                               yn, rho0, status )
        IF_NOT_OK_RETURN(status=1)

        !if(obsData(ipar(n))%error_rel>0)&
        !  obs_stddev(n)=max(obs_stddev(n),obs(n)*obsData(ipar(n))%error_rel)
        !if(obsData(ipar(n))%error_rep>0)&
        !  obs_stddev(n)=max(obs_stddev(n),obsData(ipar(n))%error_rep)

        !-----------------------------------------------------------------------
        !**************** check  Unimod units **********************************
        ! Innovation = yn-y
        !  for chemical observations: on input, obs(n)=y=observation [g/m3];
        !    on output innov(n)=innovation [g/kg].
        !    Convert standard deviation from g/m3 to g/kg
        !  for NO2 Trop.Col. observation are converted to [g/kg] uppon reading
        !    with the use of sfalling factors
        !  for lidar observations: on input, obs(n)=y [1/m/sr]
        !    likewise on output
        !-----------------------------------------------------------------------
        !if(trim(obspararr(ipar(n)))=='BSCA')then
        !  innov(n)=yn-obs(n)
        !  obsstddev(n)=obs_stddev(n)
        !else
        !  innov(n)=yn-obs(n)/rho0
        !  obsstddev(n)=obs_stddev(n)/rho0
        !endif
        
        ! store observation, scale towards model units:
        Hops%obs(n)%obs = obs(n)*FGSCALE
        
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
        Hops%obs(n)%innov = yn - obs(n)*FGSCALE

        ! scale obs.error std.dev.:
        Hops%obs(n)%obsstddev = obs_stddev(n)*FGSCALE
        
        ! store station code:
        Hops%obs(n)%stncode = trim(stncodes(n))

        ! info ..
        if ( debug .and. MasterProc ) then
          write (gol, '("#",I0,2(1X,A3,":",E12.3))') &
                   n, 'Observation', obs(n), 'Model', yn*FGSCALE_INV; call goPr
        end if

      end do  ! observations
      
    end if  ! nobs > 0
    
    ! ok
    status = 0

  end subroutine get_innovations


  ! ***


  subroutine costFunction( itime, nv_hcr, chi_hcr, &
                            Hops, Bsqrt, &
                            Jcost, Jb, Jobs, gradJcost_hcr, l2w_hcr, &
                            dx_loc, post, status )
  !-----------------------------------------------------------------------
  ! @description
  ! Compute the costfunction and its gradient
  !
  ! @author M.Kahnert
  !-----------------------------------------------------------------------
    use MPIF90               , only : MPIF90_AllReduce, MPI_SUM
    use MPIF90               , only : MPIF90_BCast
    use GO                   , only : GO_Timer_Start, GO_Timer_End
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
    use DA_Obs_ml            , only : nchemobs, ichemObs
    use DA_Obs_ml            , only : varName, varSpec
    use DA_Obs_ml            , only : obsVarName
    use DA_Obs_ml            , only : obsData
    use DA_ml                , only : nlev

    !-----------------------------------------------------------------------
    ! Formal parameters
    !-----------------------------------------------------------------------

    integer, intent(in)               ::  itime            ! record for hour-of-the-day
    integer, intent(in)               ::  nv_hcr           ! number of elements in half-complex-to-real state
    real, intent(in)                  ::  chi_hcr(nv_hcr)  ! input state
    type(T_ObsOpers), intent(in)      ::  Hops
    type(T_Bcovarsqrt), intent(inout) ::  Bsqrt
    real, intent(out)                 ::  Jcost, Jb, Jobs
    real, intent(out)                 ::  gradJcost_hcr(nv_hcr)
    integer, intent(out)              ::  l2w_hcr(nv_hcr)
    real, intent(out)                 ::  dx_loc(:,:,:,:)  ! (lnx,lny,nlev,ntracer)
    logical, intent(in)               ::  post             ! post processing call ?
    integer, intent(out)              ::  status
    
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
    integer, allocatable            ::  hcn_local(:)   ! (nw_local)

    integer                         ::  i, j, k
    integer                         ::  n
    integer                         ::  p
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

    ! info ...
    write (gol,'(a,":   evaluate costFunction ...")') rname; call goPr
    
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
    call Bsqrt%Get( status, nlon=nxg, nlat_local=lnyg, ilat_offset=iyg0 )
    IF_NOT_OK_RETURN(status=1)

    ! size of current complexe state:
    call Bsqrt%TimeGet( itime, status, nw_local=nw_local )
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

      ! copy from half-complex-real into half-complex,
      ! loop over complex elements:
      do k = 1, nw_local
        ! copy pair as real and imag part:
        chi_hc(k) = cmplx( chi_hcr(2*k-1), chi_hcr(2*k) )
      end do ! complex elements

      ! storage:
      allocate( hcn_local(nw_local), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! number of values represented in half-complex storage:
      call Bsqrt%TimeGet( itime, status, hcn_local=hcn_local )
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

      ! set to zero for safety:
      chi_hc = cmplx(0e0,0e0)

    end if


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! conversion from spectral to model space
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! info ...
    write (gol,'(a,":   evaluate dx = B^{1/2} w ...")') rname; call goPr
    
    ! transform ...
    call Bsqrt%Reverse( itime, chi_hc, dx_loc, status )
    IF_NOT_OK_RETURN(status=1)


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

      ! info ..
      write (gol,'(a,": posteriori evaluation")') rname; call goPr

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

      ! check ..
      if ( size(Hops%obs(1)%H_jac,3) /= Bsqrt%ntracer ) then
        write (gol,'("covariance for ",i0," tracers, but H_jac for ",i0)') Bsqrt%ntracer, size(Hops%obs(1)%H_jac,3); call goErr
        TRACEBACK; status=1; return
      end if
    
      ! storage:    
      allocate( yn(Hops%nobs), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dep(Hops%nobs), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( OinvDep(Hops%nobs), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! loop over observations:
      do n = 1, Hops%nobs
        ! init sum:
        yn(n) = 0.0e0
        ! observed covariance tracer index:
        itracer = Hops%obs(n)%itracer
        ! level range:
        l0 = Hops%obs(n)%l(0)
        l1 = Hops%obs(n)%l(1)
        ! loop over points involved in horizontal interpolation:
        do p = 1, 4
          ! involved grid cell:
          i = Hops%obs(n)%i(p)
          j = Hops%obs(n)%j(p)
          ! add contribution for this cell, sum over layers:
          yn(n) = yn(n) + sum( Hops%obs(n)%H_jac(p,l0:l1,itracer) * dx_loc(i,j,l0:l1,itracer) )
        end do
      end do

    end if  ! nobs > 0

    !-----------------------------------------------------------------------
    ! adjoint forcing:  H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
    !-----------------------------------------------------------------------

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! departures: O^{-1} * [ H(xb)+H_jac*dx-y ]
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! info ..
    write (gol,'(a,":   compute departures ...")') rname; call goPr
    
    ! any local observations ?
    if ( Hops%nobs > 0 ) then
      ! loop:
      do n = 1, Hops%nobs
        ! departures: 
        !  dep =  [h(xb) + H_jac*dx] - obs 
        !      =        h(xb) - obs      +   H_jac*dx 
        dep(n) =     Hops%obs(n)%innov   +    yn(n)
        ! O^{-1} * [ H(xb)+H_jac*dx - obs ]
        !          = [ H(xb)+H_jac*dx - obs ] /       sigma_obs**2
        OinvDep(n) =         dep(n)           / (Hops%obs(n)%obsstddev**2)
      end do
    end if


    !-----------------------------------------------------------------------
    ! adjoint forcing:  H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
    !-----------------------------------------------------------------------

    ! info ..
    write (gol,'(a,":   compute forcing ...")') rname; call goPr

#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_innov_adj, status )
    IF_NOT_OK_RETURN(status=1)
#endif
    
    ! local y-slab defined ?
    if ( lnyg > 0 ) then
      ! storage on extended domain, local y-slab:
      allocate( HT_OinvDep_loc(nxg,lnyg,nlev,Bsqrt%ntracer), stat=status )
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
      ! loop over involved grid cells:
      do p = 1, 4
        ! grid cell indices:
        i = Hops%obs(n)%i(p)
        j = Hops%obs(n)%j(p)
        ! level range:
        l0 = Hops%obs(n)%l(0)
        l1 = Hops%obs(n)%l(1)
        ! covariance tracer index:
        itracer = Hops%obs(n)%itracer
        ! first on local array only:
        HT_OinvDep_loc(i,j,l0:l1,itracer) = HT_OinvDep_loc(i,j,l0:l1,itracer) + Hops%obs(n)%H_jac(p,l0:l1,itracer) * OinvDep(n)
        !! testing ...
        !write (gol,'("HT_OinvDep_loc n=",i0," p=",i0," l0:l1=",i0,":",i0," ichem=",i0," H_jac=",es12.4," OinvDep=",es12.4)') &
        !         n, p, l0, l1, ichem, Hops%obs(n)%H_jac(p,l0:l1,ichem), OinvDep(n); call goPr
      end do
    end do
    
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
    
    ! info ..
    write (gol,'(a,":   compute gradJb ...")') rname; call goPr
    
    ! gradient to elements of chi:
    gradJb_hcr = chi_hcr

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! grad Jobs
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! info ..
    write (gol,'(a,":   HT_OinvDep_loc : ",2es12.4)') rname, minval(HT_OinvDep_loc), maxval(HT_OinvDep_loc); call goPr

    ! info ...
    write (gol,'(a,":   evaluate w = B^{H/2} (H^T R^{-1} dy) ...")') rname; call goPr
    
#ifdef with_ajs
    ! start timing:
    call GO_Timer_Start( itim_x2chi, status )
    IF_NOT_OK_RETURN(status=1)
#endif
    
    ! transform ...
    call Bsqrt%Forward( itime, HT_OinvDep_loc, chi_hc, status )
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
    ! info ..
    write (gol,'(a,":     norm(gradJb  ) : ",e16.6)') rname, gradnorm; call goPr
    
    ! info: local l2 sum:
    gradnorm_loc = sum( l2w_hcr * gradJobs_hcr**2 )
    ! collect global sum:
    call MPIF90_AllReduce( gradnorm_loc, gradnorm, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! convert to norm:
    gradnorm = sqrt(gradnorm)
    ! info ..
    write (gol,'(a,":     norm(gradJobs) : ",e16.6)') rname, gradnorm; call goPr
    
    ! info: local l2 sum:
    gradnorm_loc = sum( l2w_hcr * gradJcost_hcr**2 )
    ! collect global sum:
    call MPIF90_AllReduce( gradnorm_loc, gradnorm, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
    ! convert to norm:
    gradnorm = sqrt(gradnorm)
    ! info ..
    write (gol,'(a,":     norm(gradJ   ) : ",e16.6)') rname, gradnorm; call goPr

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

    ! info ..
    write (gol,'(a,":   compute Jb ...")')  rname; call goPr
    
    ! norm over local elements of chi (complex!)
    ! which is the same as the norm over the half-complex-to-real elements
    ! taking into account weights for the redundant values:
    Jb_loc = 0.5 * sum( l2w_hcr * chi_hcr**2 )
    ! total over all domains:
    call MPIF90_AllReduce( Jb_loc, Jb, MPI_SUM, MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)

    ! info ...
    write (gol,'(a,":     Jb       = ",e16.6)') rname, Jb; call goPr


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Jobs = 0.5 * [ H(xb)+H_jac*dx-y ]^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
    !      = 0.5 *          dep             *         OinvDep
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! info ..
    write (gol,'(a,":   compute Jobs ...")') rname; call goPr
    
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
    
    ! info ...
    write (gol,'(a,":     Jobs     = ",e16.6)') rname, Jobs; call goPr


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! J = Jb + Jobs
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! combine:
    Jcost = Jb + Jobs
    
    ! info ...
    write (gol,'(a,":     Jcost    = ",e16.6)') rname, Jcost; call goPr


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
