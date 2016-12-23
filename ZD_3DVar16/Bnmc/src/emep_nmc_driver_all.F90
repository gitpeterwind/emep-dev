!######################################################################
!
! EMEP DA NMC - driver routines
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!######################################################################

module EMEP_NMC_Driver_ALL

  use GO                  , only : gol, goPr, goErr
  use GO                  , only : TDate
  use EMEP_NMC_ModelOutput, only : T_ModelOutputSeries

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  ParseArguments
  public  ::  T_ModelRuns
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver'
  
  ! number of model runs, probably 2 for NMC test:
  integer, parameter  ::  nrun = 2

  ! assumed maximum number of variables:
  integer, parameter  ::  maxvar = 10

  ! assumed maximum number of times within days:
  integer, parameter  ::  maxhours = 24

  ! assumed maximum number of point:
  integer, parameter  ::  maxpoint = 10


  ! --- types ----------------------------------------
  
  type T_ModelField
    character(len=32)           ::  name
    character(len=32)           ::  units
  end type T_ModelField
  
  ! *
  
  type T_Point
    character(len=32)           ::  name
    real                        ::  lon, lat   ! deg
    integer                     ::  ilon, ilat, ilev
  end type T_Point
  
  ! *
  
  type T_ModelRuns
    ! model run info:
    type(T_ModelOutputSeries)         ::  run(nrun)
    ! time range:
    type(TDate)                       ::  tday1, tday2
    ! hours within day:
    integer                           ::  nhour
    integer                           ::  hours(maxhours)
    ! variables:
    integer                           ::  nvar
    type(T_ModelField), allocatable   ::  var(:)
    ! points:
    integer                           ::  npoint
    type(T_Point), allocatable        ::  point(:)
    ! statistics files created or read:
    character(len=1024)               ::  eps_samples_filename
    character(len=1024)               ::  eps_stats_filename
    character(len=1024)               ::  eta_samples_filename
    character(len=1024)               ::  eta_f_samples_filename
    character(len=1024)               ::  C_f_filename
    character(len=1024)               ::  D_f_filename
    character(len=1024)               ::  gamma_filename
    character(len=1024)               ::  B_f_filename
    character(len=1024)               ::  XLX_filename
    character(len=1024)               ::  BB_filename
  contains
    procedure   ::  Init            => ModelRuns_Init
    procedure   ::  Done            => ModelRuns_Done
    procedure   ::  Eps_Samples     => ModelRuns_Eps_Samples
    procedure   ::  Evaluate_Eps    => ModelRuns_Evaluate_Eps
    procedure   ::  Eta_Samples     => ModelRuns_Eta_Samples
    procedure   ::  Eta_f_Samples   => ModelRuns_Eta_f_Samples
    procedure   ::  Compute_C_f     => ModelRuns_Compute_C_f
    procedure   ::  Compute_D_f     => ModelRuns_Compute_D_f
    procedure   ::  Compute_B_f     => ModelRuns_Compute_B_f
    procedure   ::  Compute_XLX     => ModelRuns_Compute_XLX
    procedure   ::  Collect_BB      => ModelRuns_Collect_BB
    procedure   ::  Evaluate_BB     => ModelRuns_Evaluate_BB
  end type T_ModelRuns
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** arguments
  ! ***
  ! ********************************************************************


  !
  ! Return status:
  !     -1 : usage displayed
  !      0 : ok
  !  other : error
  !
  
  subroutine ParseArguments( rcfile, status )
  
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(out)    ::  rcfile
    integer, intent(out)             ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ParseArguments'
    
    ! assumed target name of executable fromed:
    character(len=*), parameter  ::  xname = 'emep_nmc.x'
    
    ! --- local ----------------------------------
    
    integer                         ::  narg, iarg
    integer                         ::  arglength
    character(len=1024)             ::  argvalue
    integer                         ::  k
    
    ! --- begin ----------------------------------
    
    ! dummy values at start:
    rcfile = ''
    
    ! count:
    narg = Command_Argument_Count()

    ! loop:
    do iarg = 1, narg
      ! get length:
      call Get_Command_Argument( iarg, length=arglength, status=status )
      if ( status /= 0 ) then
        write (gol,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength; call goErr
        TRACEBACK; stop 1
      end if
      ! get value: 
      call Get_Command_Argument( iarg, value=argvalue, status=status )
      if ( status /= 0 ) then
        write (gol,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength; call goErr
        TRACEBACK; stop 1
      end if
      ! switch:
      select case ( trim(argvalue) )
        !~ help:
        case ( '-h', '--help' )
          write (gol,'("")'); call goPr
          write (gol,'("Usage:")'); call goPr
          write (gol,'("  ",a," rcfile")') xname; call goPr
          write (gol,'("  ",a," -h|--help")') xname; call goPr
          write (gol,'("")'); call goPr
          status = -1; return
        !~ values:
        case default
          ! switch:
          if ( len_trim(rcfile) == 0 ) then
            rcfile = trim(argvalue)
          else
            write (*,'("unsupported argument `",a,"`")') trim(argvalue)
            TRACEBACK; stop 1
          end if
      end select
    end do   ! arguments

    ! check ...
    if ( len_trim(rcfile) == 0 ) then
      write (*,'("no input rcfile argument specified")')
      TRACEBACK; stop 1
    end if
    
    ! ok
    status = 0
  
  end subroutine ParseArguments


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine ModelRuns_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine, goSplitString
  
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(out)           ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Init'
    
    ! --- local ----------------------------------
    
    integer               ::  irun
    character(len=32)     ::  rckey
    character(len=32)     ::  id
    character(len=32)     ::  tvalue
    character(len=1024)   ::  line
    character(len=32)     ::  vnames(maxvar)
    integer               ::  ivar
    character(len=32)     ::  pnames(maxpoint)
    integer               ::  ipoint
    
    ! --- begin ----------------------------------
    
    ! start time:
    call ReadRc( rcF, 'nmc.timerange.t1', tvalue, status )
    IF_NOT_OK_RETURN(status=1)
    ! extract:
    call goReadFromLine( tvalue, self%tday1, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! end time:
    call ReadRc( rcF, 'nmc.timerange.t2', tvalue, status )
    IF_NOT_OK_RETURN(status=1)
    ! extract:
    call goReadFromLine( tvalue, self%tday2, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! steps within day as list:
    call ReadRc( rcF, 'nmc.hours', line, status )
    IF_NOT_OK_RETURN(status=1)
    ! extract:
    self%nhour = 0
    do
      ! empty ?
      if ( len_trim(line) == 0 ) exit
      ! next:
      self%nhour = self%nhour + 1
      ! check ...
      if ( self%nhour > size(self%hours) ) then
        write (gol,'("number of time steps exceeds maximum storage ",i6)') size(self%hours); call goErr
        TRACEBACK; status=1; return
      end if
      ! extract:
      call goReadFromLine( line, self%hours(self%nhour), status, sep=' ' )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! check ...
    if ( self%nhour == 0 ) then
      write (gol,'("no time steps specified")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! variable names:
    call ReadRc( rcF, 'nmc.variables', line, status )
    IF_NOT_OK_RETURN(status=1)
    ! split:
    call goSplitString( trim(line), self%nvar, vnames, status )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%var(self%nvar) )
    ! fill names:
    do ivar = 1, self%nvar
      self%var(ivar)%name = trim(vnames(ivar))
    end do
    
    ! loop over runs:
    do irun = 1, nrun
      ! info ...
      write (gol,'("run ",i2," ...")') irun; call goPr
      ! write key for id:
      write (rckey,'("nmc.run.id",i1)') irun
      ! read id:
      call ReadRc( rcF, trim(rckey), id, status )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'("run id : ",a)') trim(id); call goPr
      ! init runs:
      call self%run(irun)%Init( rcF, 'nmc.run', id, vnames(1:self%nvar), status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! point names:
    call ReadRc( rcF, 'nmc.points', line, status )
    IF_NOT_OK_RETURN(status=1)
    ! split:
    call goSplitString( trim(line), self%npoint, pnames, status )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%point(self%npoint) )
    ! loop:
    do ipoint = 1, self%npoint
      ! store:
      self%point(ipoint)%name = trim(pnames(ipoint))
      ! location:
      call ReadRc( rcF, 'nmc.point.'//trim(pnames(ipoint))//'.lon', self%point(ipoint)%lon, status )
      IF_NOT_OK_RETURN(status=1)
      call ReadRc( rcF, 'nmc.point.'//trim(pnames(ipoint))//'.lat', self%point(ipoint)%lat, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! input/output:
    call ReadRc( rcF, 'nmc.eps.samples.filename', self%eps_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.eps.stats.filename', self%eps_stats_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.eta.samples.filename', self%eta_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.eta_f.samples.filename', self%eta_f_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.C_f.filename', self%C_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.D_f.filename', self%D_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.gamma.filename', self%gamma_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.B_f.filename', self%B_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.XLX.filename', self%XLX_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.BB.filename', self%BB_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine ModelRuns_Init


  ! ***
  

  subroutine ModelRuns_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! loop over runs:
    do irun = 1, nrun
      ! done with runs:
      call self%run(irun)%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! clear variable info:
    deallocate( self%var )
    
    ! clear point info:
    deallocate( self%point )
    
    ! ok
    status = 0
  
  end subroutine ModelRuns_Done


  ! ***
  
  
  !
  ! Collect samples of NMC differences:
  !
  !   eps_comp(sample,hour,lev,lat,lon)  =  x1_comp(sample,hour,lev,lat,lon) - x2_comp(sample,hour,lev,lat,lon)
  !
  ! Meanwhile also compute statistics of eps:
  !
  !   eps_comp_mean(hour,lev,lat,lon) =  sample_mean( eps_comp )
  !   eps_comp_stdv(hour,lev,lat,lon) =  sample_stdv( eps_comp )
  !

  subroutine ModelRuns_Eps_Samples( self, rcF, status )
  
    use GO             , only : TDate, IncrDate
    use GO             , only : operator(+), operator(>), wrtgol
    use GO             , only : T_SampStat
    use GO             , only : TrcFile, ReadRc
    use GO             , only : goReplace
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : RealCoordinate, HybrideLevelCoordinate, TimeCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Eps_Samples'
    
    ! --- local ----------------------------------
    
    integer                          ::  isample
    type(TDate)                      ::  tday, t
    integer                          ::  ihour, ihr
    integer                          ::  irun
    logical                          ::  setup
    logical                          ::  isetup(nrun)
    integer                          ::  ivar
    integer                          ::  nlon, nlat, nlev
    integer                          ::  ilon, ilat, ilev
    real, allocatable                ::  lons(:)         ! (nlon)
    real, allocatable                ::  lats(:)         ! (nlat)
    real, allocatable                ::  hyam(:), hybm(:)     ! (nlev)
    real, allocatable                ::  hyai(:), hybi(:)     ! (nlev+1)
    real                             ::  p0
    character(len=32)                ::  vname_ps
    character(len=32)                ::  units_p
    character(len=32)                ::  units
    real, allocatable                ::  samp_ps  (:,:    ,:)   ! (nlon,nlat     ,    ,nrun)
    real, allocatable                ::  samp_data(:,:,:,:,:)   ! (nlon,nlat,nlev,nvar,nrun)
    real, allocatable                ::  aver_ps  (:,:  )       ! (nlon,nlat     )
    real, allocatable                ::  eps_g     (:,:,:)      ! (nlon,nlat,nlev)

    type(T_SampStat), allocatable    ::  sampstat_eps_g(:,:)    ! (nhour,nvar)
    type(T_SampStat), allocatable    ::  sampstat_ps(:)         ! (nhour)
    
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(TimeCoordinate)             ::  time_coor, ctime_coor, sample_coor

    character(len=1024)              ::  filename
    character(len=1024)              ::  filename_curr
    type(T_NMC_Output)               ::  eps_file

    real, allocatable                ::  ps(:,:)        ! (nlon,nlat,ntime)
    real, allocatable                ::  dmean(:,:,:)   ! (nlon,nlat,nlev)
    real, allocatable                ::  dstdv(:,:,:)   ! (nlon,nlat,nlev)
    integer, allocatable             ::  fixed(:)   ! (nlev)
    integer                          ::  ps__varid
    integer, allocatable             ::  eps__varids(:)     ! (nvar)
    integer, allocatable             ::  dmean__varid(:)   ! (nvar)
    integer, allocatable             ::  dstdv__varid(:)   ! (nvar)
    integer, allocatable             ::  fixed__varids(:)   ! (nvar)
    type(TDate)                      ::  climat(2)

    character(len=128)               ::  run_labels(nrun)
    character(len=1024)              ::  msg
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Collect NMC samples **")'); call goPr
    write (gol,'("")'); call goPr

    ! start time:
    tday = self%tday1
    ! set flag:
    setup = .false.
    isetup = .false.
    ! no filename defined yet:
    filename_curr = ''
    ! time loop:
    do
      ! info ...
      call wrtgol( 'day ', tday ); call goPr
      
      ! loop over hours:
      do ihour = 1, self%nhour

        ! info ...
        write (gol,'("  hour ",i4," ...")') ihour; call goPr
      
        ! current time:
        t = tday + IncrDate(hour=self%hours(ihour))
      
        ! loop over runs:
        do irun = 1, nrun

          ! info ...
          write (gol,'("    run ",i4," ...")') irun; call goPr

          ! read into buffer if necessary:
          call self%run(irun)%ReadBuffer( t, status )
          IF_NOT_OK_RETURN(status=1)

          ! store ?
          if ( .not. isetup(irun) ) then
            ! get label:
            call self%run(irun)%Get( status, label=run_labels(irun) )
            IF_NOT_OK_RETURN(status=1)
            ! loop over variables:
            do ivar = 1, self%nvar
              ! set or check units:
              if ( irun == 1 ) then
                ! copy units:
                self%var(ivar)%units = trim(self%run(irun)%var(ivar)%units)
              else
                ! check ...
                if ( trim(self%var(ivar)%units) /= trim(self%run(irun)%var(ivar)%units) ) then
                  write (gol,'("units do not match for ",i0," ",a)') ivar, trim(self%var(ivar)%name); call goErr
                  write (gol,'("  first run   : ",a)')  trim(self%var(ivar)%units); call goErr
                  write (gol,'("  current run : ",a)')  trim(self%run(irun)%var(ivar)%units); call goErr
                  TRACEBACK; status=1; return
                end if
              end if
            end do  ! var
            ! reset flag:
            isetup(irun) = .false.
          end if
          
        end do  ! nmc runs

        ! setup data ?
        if ( .not. setup ) then

          ! info ...
          write (gol,'("    initial setup ...")'); call goPr

          ! use first run:
          irun = 1
          ! get dimensions etc:
          call self%run(irun)%Get( status, nlon=nlon, nlat=nlat, nlev=nlev )
          IF_NOT_OK_RETURN(status=1)

          ! storage:
          allocate( lons (nlon), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( lats (nlat), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( hyam(nlev), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( hybm(nlev), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( hyai(nlev+1), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( hybi(nlev+1), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( samp_ps   (nlon,nlat               ,nrun), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( samp_data (nlon,nlat,nlev,self%nvar,nrun), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( aver_ps   (nlon,nlat     ), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( eps_g     (nlon,nlat,nlev), stat=status )
          IF_NOT_OK_RETURN(status=1)

          ! coordinates:
          call self%run(irun)%Get( status, lons=lons, lats=lats, &
                                      hyam=hyam, hybm=hybm, hyai=hyai, hybi=hybi, &
                                      p0=p0, vname_ps=vname_ps, units_p=units_p )
          IF_NOT_OK_RETURN(status=1)

          ! create longitude coordinate:
          call lon_coor%Init( 'longitude', status )
          IF_NOT_OK_RETURN(status=1)
          ! set size:
          call lon_coor%Set_Dim( status, n=nlon )
          IF_NOT_OK_RETURN(status=1)
          ! set attributes:
          call lon_coor%Set_Attrs( status, units='degrees_east', standard_name='longitude' )
          IF_NOT_OK_RETURN(status=1)
          ! fill values:
          call lon_coor%Set_Values( status, values=lons )
          IF_NOT_OK_RETURN(status=1)

          ! create latitude coordinate:
          call lat_coor%Init( 'latitude', status )
          IF_NOT_OK_RETURN(status=1)
          ! set size:
          call lat_coor%Set_Dim( status, n=nlat )
          IF_NOT_OK_RETURN(status=1)
          ! set attributes:
          call lat_coor%Set_Attrs( status, units='degrees_north', standard_name='latitude' )
          IF_NOT_OK_RETURN(status=1)
          ! fill values:
          call lat_coor%Set_Values( status, values=lats )
          IF_NOT_OK_RETURN(status=1)
          
          ! create level coordinate:
          call lev_coor%Init( 'lev', status )
          IF_NOT_OK_RETURN(status=1)
          ! set size:
          call lev_coor%Set_Dim( status, n=nlev )
          IF_NOT_OK_RETURN(status=1)
          ! fill attributes and values:
          call lev_coor%Set_Values( status, ap=hyam, b=hybm, api=hyai, bi=hybi, &
                                       p0=p0, units_p=units_p, &
                                       vname_ps=trim(vname_ps) )
          IF_NOT_OK_RETURN(status=1)
          
          ! create hour coordinate:
          call time_coor%Init( 'time', status )
          IF_NOT_OK_RETURN(status=1)
          ! set size:
          call time_coor%Set_Dim( status, n=self%nhour )
          IF_NOT_OK_RETURN(status=1)
          ! fill attributes:
          call time_coor%Set_Values( status, units_step='hours', units_t0=self%tday1 )
          IF_NOT_OK_RETURN(status=1)
          
          ! create sample coordinate:
          call sample_coor%Init( 'sample', status )
          IF_NOT_OK_RETURN(status=1)
          ! set unlimited size:
          call sample_coor%Set_Dim( status, unlimited=.true. )
          IF_NOT_OK_RETURN(status=1)
          ! fill attributes:
          call sample_coor%Set_Values( status, units_step='days', units_t0=self%tday1 )
          IF_NOT_OK_RETURN(status=1)
          
          ! for average pressure:
          allocate( sampstat_ps(self%nhour), stat=status )
          IF_NOT_OK_RETURN(status=1)
          do ihr = 1, self%nhour
            call sampstat_ps(ihr)%Init( (/nlon,nlat/), status )
            IF_NOT_OK_RETURN(status=1)
          end do  ! hour

          ! access:
          allocate( eps__varids(self%nvar) )

          ! stats:
          allocate( sampstat_eps_g(self%nhour,self%nvar), stat=status )
          IF_NOT_OK_RETURN(status=1)
          do ihr = 1, self%nhour
            do ivar = 1, self%nvar
              call sampstat_eps_g(ihr,ivar)%Init( (/nlon,nlat,nlev/), status )
              IF_NOT_OK_RETURN(status=1)
            end do  ! var
          end do  ! hour

          ! reset flag:
          setup = .true.
        end if  ! not setup yet
        
        ! ~ init file for eps

        ! target file:
        filename = trim(self%eps_samples_filename)
        ! replace keys:
        call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
        IF_NOT_OK_RETURN(status=1)
        call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
        IF_NOT_OK_RETURN(status=1)
        call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
        IF_NOT_OK_RETURN(status=1)

        ! new ?
        if ( trim(filename) /= trim(filename_curr) ) then

          ! info ...
          write (gol,'("    create output file: ",a)') trim(filename); call goPr
          
          ! close existing?
          if ( len_trim(filename_curr) > 0 ) then
            ! close output file:
            call eps_file%Close( status )
            IF_NOT_OK_RETURN(status=1)
          end if

          ! copy:
          filename_curr = trim(filename)

          ! create new result file:
          call eps_file%Create( trim(filename), status )
          IF_NOT_OK_RETURN(status=1)
      
          ! description:
          write (msg,'("run labels `",a,"` and `",a,"`")') &
                   trim(run_labels(1)), trim(run_labels(2))
          ! add:
          call eps_file%Extend_History( trim(msg), status )
          IF_NOT_OK_RETURN(status=1)

          ! define coordinates in output file:
          call lon_coor%Def( eps_file, status )
          IF_NOT_OK_RETURN(status=1)
          call lat_coor%Def( eps_file, status )
          IF_NOT_OK_RETURN(status=1)

          ! define level coordinate in output file:
          call lev_coor%Def( eps_file, status )
          IF_NOT_OK_RETURN(status=1)

          ! define hour coordinate in output file:
          call time_coor%Def( eps_file, status )
          IF_NOT_OK_RETURN(status=1)

          ! define sample coordinate in output file:
          call sample_coor%Def( eps_file, status )
          IF_NOT_OK_RETURN(status=1)
          
          ! init sample counter:
          isample = 1

          ! define:
          call eps_file%Def_Sample2D( trim(vname_ps), trim(units_p), &
                                       lon_coor, lat_coor, time_coor, sample_coor, &
                                       ps__varid, status )
          IF_NOT_OK_RETURN(status=1)
          call eps_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
          IF_NOT_OK_RETURN(status=1)

          ! loop over variables:
          do ivar = 1, self%nvar
            ! define:
            call eps_file%Def_Sample3D( trim(self%var(ivar)%name), &
                                    trim(self%var(ivar)%units), &
                                    lon_coor, lat_coor, lev_coor, time_coor, sample_coor, &
                                    eps__varids(ivar), status )
            IF_NOT_OK_RETURN(status=1)
            call eps_file%Put_Att( eps__varids(ivar), 'long_name', trim(self%var(ivar)%name), status )
            IF_NOT_OK_RETURN(status=1)
          end do  ! var

          ! end:
          call eps_file%EndDef( status )
          IF_NOT_OK_RETURN(status=1)

          ! write grid coordinate variables:
          call lon_coor%Write( eps_file, status )
          IF_NOT_OK_RETURN(status=1)
          call lat_coor%Write( eps_file, status )
          IF_NOT_OK_RETURN(status=1)
          
          ! write vertical coordinate variables:
          call lev_coor%Write( eps_file, status )
          IF_NOT_OK_RETURN(status=1)
          
        end if  ! new output file
        
        ! ~

        ! info ...
        write (gol,'("    get records ...")'); call goPr
      
        ! loop over runs:
        do irun = 1, nrun

          ! loop over variables:
          do ivar = 1, self%nvar
            ! read data:
            call self%run(irun)%GetRecord( t, samp_ps(:,:,irun), ivar, &
                                            samp_data(:,:,:,ivar,irun), &
                                            self%var(ivar)%units, status )
            IF_NOT_OK_RETURN(status=1)
          end do ! var

        end do  ! nmc runs

        ! average surface pressure:
        aver_ps = 0.5*( samp_ps(:,:,1) + samp_ps(:,:,2) )
        ! update statistics:
        call sampstat_ps(ihour)%AddSample( aver_ps, status )
        IF_NOT_OK_RETURN(status=1)

        ! first sample? then write hour info:
        if ( isample == 1 ) then
          ! time for this hour on first day:
          t = self%tday1 + IncrDate(hour=self%hours(ihour))
          ! write time value:
          call time_coor%Write_Value( eps_file, ihour, status, t=t )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! first hour ? then write sample info:
        if ( ihour == 1 ) then
          ! write sample time:
          call sample_coor%Write_Value( eps_file, isample, status, t=tday )
          IF_NOT_OK_RETURN(status=1)
        end if

        ! write surface pressure:
        call eps_file%Put_Sample2D( ps__varid, ihour, isample, aver_ps, status )
        IF_NOT_OK_RETURN(status=1)

        ! vars:
        do ivar = 1, self%nvar
          ! difference field:
          eps_g = samp_data(:,:,:,ivar,2) - samp_data(:,:,:,ivar,1)
          ! write sample:
          call eps_file%Put_Sample3D( eps__varids(ivar), ihour, isample, eps_g, status )
          IF_NOT_OK_RETURN(status=1)
          ! update statistics:
          call sampstat_eps_g(ihour,ivar)%AddSample( eps_g, status )
          IF_NOT_OK_RETURN(status=1)
        end do ! vars
        
      end do  ! hours

      ! next:
      tday = tday + IncrDate(day=1)
      ! end ?
      if ( tday > self%tday2 ) exit
      
      ! increase counter:
      isample = isample + 1
      
    end do  ! days

    ! close output file:
    call eps_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! info ...
    write (gol,'("post processing ...")'); call goPr

    ! data present ?
    if ( setup ) then

      ! clear:
      deallocate( samp_ps )
      deallocate( samp_data )
      deallocate( aver_ps )
      deallocate( eps_g )

      ! ~ output

      ! info ...
      write (gol,'("  write statistics ...")'); call goPr
      
      ! storage:
      allocate( ps(nlon,nlat), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dmean(nlon,nlat,nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dstdv(nlon,nlat,nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( fixed(nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! access:
      allocate( dmean__varid(self%nvar) )
      allocate( dstdv__varid(self%nvar) )
      allocate( fixed__varids(self%nvar) )
      
      ! create new result file:
      call eps_file%Create( trim(self%eps_stats_filename), status )
      IF_NOT_OK_RETURN(status=1)
      
      ! description:
      write (msg,'("run labels `",a,"` and `",a,"`")') &
               trim(run_labels(1)), trim(run_labels(2))
      ! add:
      call eps_file%Extend_History( trim(msg), status )
      IF_NOT_OK_RETURN(status=1)
      
      ! define coordinates in output file:
      call lon_coor%Def( eps_file, status )
      IF_NOT_OK_RETURN(status=1)
      call lat_coor%Def( eps_file, status )
      IF_NOT_OK_RETURN(status=1)

      ! define hour coordinate, now climatology:
      call ctime_coor%Init( 'time', status )
      IF_NOT_OK_RETURN(status=1)
      ! set size:
      call ctime_coor%Set_Dim( status, n=self%nhour )
      IF_NOT_OK_RETURN(status=1)
      ! fill attributes:
      call ctime_coor%Set_Values( status, climatology=.true., &
                                   units_step='days', units_t0=self%tday1 )
      IF_NOT_OK_RETURN(status=1)

      ! define hour coordinate in output file:
      call ctime_coor%Def( eps_file, status )
      IF_NOT_OK_RETURN(status=1)

      ! define level coordinate in output file:
      call lev_coor%Def( eps_file, status )
      IF_NOT_OK_RETURN(status=1)

      ! define:
      call eps_file%Def_Field2D_Series( trim(vname_ps), trim(units_p), &
                                   lon_coor, lat_coor, time_coor, &
                                   ps__varid, status )
      IF_NOT_OK_RETURN(status=1)
      call eps_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over variables:
      do ivar = 1, self%nvar
        ! define:
        call eps_file%Def_Field3D_Series( trim(self%var(ivar)%name)//'_mean', &
                               trim(self%var(ivar)%units), &
                               lon_coor, lat_coor, lev_coor, ctime_coor, &
                               dmean__varid(ivar), status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( dmean__varid(ivar), 'long_name', trim(self%var(ivar)%name)//' mean', status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( dmean__varid(ivar), 'cell_methods', 'time: point within days time: mean over days', status )
        IF_NOT_OK_RETURN(status=1)
        ! define:
        call eps_file%Def_Field3D_Series( trim(self%var(ivar)%name)//'_stdv', &
                               trim(self%var(ivar)%units), &
                               lon_coor, lat_coor, lev_coor, ctime_coor, &
                               dstdv__varid(ivar), status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( dstdv__varid(ivar), 'long_name', trim(self%var(ivar)%name)//' standard deviation', status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( dstdv__varid(ivar), 'cell_methods', 'time: point within days time: standard_deviation over days', status )
        IF_NOT_OK_RETURN(status=1)
        ! define:
        call eps_file%Def_IField1D( trim(self%var(ivar)%name)//'_fixed', '1', &
                               lev_coor, fixed__varids(ivar), status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( fixed__varids(ivar), 'long_name', trim(self%var(ivar)%name)//' fixed layer (1=fixed, 0=dynamic)', status )
        IF_NOT_OK_RETURN(status=1)
      end do  ! var

      ! end:
      call eps_file%EndDef( status )
      IF_NOT_OK_RETURN(status=1)

      ! write grid coordinate variables:
      call lon_coor%Write( eps_file, status )
      IF_NOT_OK_RETURN(status=1)
      call lat_coor%Write( eps_file, status )
      IF_NOT_OK_RETURN(status=1)

      ! write vertical coordinate variables:
      call lev_coor%Write( eps_file, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! loop over hours:
      do ihour = 1, self%nhour

        ! time for this hour on first day:
        t = self%tday1+IncrDate(hour=self%hours(ihour))
        ! climatological boundaries:
        climat(1) = self%tday1
        climat(2) = self%tday2 + IncrDate(hour=24)
        ! fill time and climatological bounds:
        call ctime_coor%Write_Value( eps_file, ihour, status, t=t, climat=climat )
        IF_NOT_OK_RETURN(status=1)

        ! extract:
        call sampstat_ps(ihour)%Get2d( status, mean=ps )
        IF_NOT_OK_RETURN(status=1)
        ! write:
        call eps_file%Put_Field2D_Series( ps__varid, ihour, ps, status )
        IF_NOT_OK_RETURN(status=1)
      
        ! loop over variables:
        do ivar = 1, self%nvar

          ! extract:
          call sampstat_eps_g(ihour,ivar)%Get3d( status, mean=dmean, stdv=dstdv )
          IF_NOT_OK_RETURN(status=1)
          
          ! write mean:
          call eps_file%Put_Field3D_Series( dmean__varid(ivar), ihour, dmean, status )
          IF_NOT_OK_RETURN(status=1)
          ! write std.dev.:
          call eps_file%Put_Field3D_Series( dstdv__varid(ivar), ihour, dstdv, status )
          IF_NOT_OK_RETURN(status=1)
          
          ! first hour ?
          if ( ihour == 1 ) then
            ! fixed layers ?
            do ilev = 1, nlev
              if ( any( dstdv(:,:,ilev) > 0.0 ) ) then
                fixed(ilev) = 0
              else
                fixed(ilev) = 1
              end if
            end do
            ! write record:
            call eps_file%Put_IField1D( fixed__varids(ivar), fixed, status )
            IF_NOT_OK_RETURN(status=1)
          end if

        end do  ! var

      end do  ! hours

      ! close:
      call eps_file%Close( status )
      IF_NOT_OK_RETURN(status=1)
      
      ! clear:
      deallocate( ps )
      deallocate( dmean )
      deallocate( dstdv )
      deallocate( fixed )

      ! ~

      ! info ...
      write (gol,'("done.")'); call goPr

      ! clear:
      deallocate( lons, lats )
      deallocate( hyam, hybm, hyai, hybi )
      
      ! clear:
      call lon_coor%Done( status )
      IF_NOT_OK_RETURN(status=1)
      call lat_coor%Done( status )
      IF_NOT_OK_RETURN(status=1)
      call lev_coor%Done( status )
      IF_NOT_OK_RETURN(status=1)
      call time_coor%Done( status )
      IF_NOT_OK_RETURN(status=1)
      call ctime_coor%Done( status )
      IF_NOT_OK_RETURN(status=1)
      call sample_coor%Done( status )
      IF_NOT_OK_RETURN(status=1)

      ! clear access:
      deallocate( eps__varids )
      deallocate( dmean__varid )
      deallocate( dstdv__varid )
      deallocate( fixed__varids )

      ! done with stats:
      do ihr = 1, self%nhour
        call sampstat_ps(ihr)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do
      deallocate( sampstat_ps )
      ! done with stats:
      do ihr = 1, self%nhour
        do ivar = 1, self%nvar
          call sampstat_eps_g(ihr,ivar)%Done( status )
          IF_NOT_OK_RETURN(status=1)
        end do
      end do
      deallocate( sampstat_eps_g )
      
    else
    
      ! info ...
      write (gol,'("WARNING - no records in time range ...")'); call goPr

    end if
    
    ! ok
    status = 0
  
  end subroutine ModelRuns_Eps_Samples


  ! ***
  
  
  !
  ! Evaluate sample covariance matrix at selected locations.
  !

  subroutine ModelRuns_Evaluate_Eps( self, rcF, status )
  
    use GO             , only : TDate, IncrDate, operator(+), operator(>), wrtgol
    use GO             , only : TrcFile, ReadRc
    use GO             , only : goReplace
    use GO             , only : T_SampStat, T_SampCovr_3D_1D, T_SampCovr_2D_2D
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Evaluate_Eps'
  
    ! --- local ----------------------------------

    type(TDate)                      ::  tday, t
    character(len=1024)              ::  filename
    character(len=1024)              ::  eps_filename
    integer                          ::  eps_isample
    logical                          ::  setup
    type(T_NMC_Output)               ::  eps_file
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(TimeCoordinate)             ::  ctime_coor
    type(TimeCoordinate)             ::  eps_sample_coor
    integer                          ::  nlon, nlat, nlev
    integer                          ::  ntime
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    real, allocatable                ::  ps(:,:)
    character(len=32)                ::  ps__units
    integer, allocatable             ::  eps__varids(:)   ! (nvar)
    character(len=32), allocatable   ::  eps__units(:)  ! (nvar)
    real, allocatable                ::  eps(:,:,:,:)   ! (nlon,nlat,nlev,nvar)
    character(len=32)                ::  vname
    integer                          ::  ilon, ilat, ilev
    integer                          ::  ivar
    integer                          ::  itime
    
    logical                          ::  correlation
    character(len=32)                ::  cv_label
    character(len=1)                 ::  cv_id
    character(len=32)                ::  cv_units

    type(LabelCoordinate)            ::  tracer_coor
    character(len=32)                ::  tracer_name
    character(len=32)                ::  tracer_units
    integer                          ::  ipoint
    real, allocatable                ::  eps_profiles(:,:,:)   ! (npoint,nlev,nvar)

    type(T_SampStat), allocatable         ::  sampstat_ps(:)        ! (nhour)
    type(T_SampStat), allocatable         ::  sampstat_eps(:,:)     ! (nhour,nvar)
    type(T_SampCovr_3D_1D), allocatable   ::  sampcovr_eps(:,:,:)   ! (nhour,nvar,nvar)
    integer                               ::  ivar2
    
    character(len=1024)              ::  B_point_filename
    type(T_NMC_Output)               ::  out_file
    character(len=1024)              ::  msg
    integer                          ::  ps__varid
    integer, allocatable             ::  S__varids(:)  ! (ntracer)
    real, allocatable                ::  S(:,:,:)         ! (nlon,nlat,nlev)
    integer, allocatable             ::  Bx__varids(:,:)  ! (npoint,ntracer)
    real, allocatable                ::  Bx(:,:,:)        ! (nlon,nlat,nlev)
    real                             ::  stdv

    logical                               ::  levtr
    type(T_SampCovr_2D_2D), allocatable   ::  sampcovr_Bzs(:,:)     ! (nhour,npoint)
    real, allocatable                     ::  Vzs(:,:,:,:)  ! (nlev,ntracer,nlev,ntracer)
    real, allocatable                     ::  Bzs(:,:,:,:)  ! (nlev,nlev,ntracer,ntracer)
    integer, allocatable                  ::  Bzs__varids(:)  ! (npoint)
    integer                               ::  k
    
    ! --- begin ----------------------------------

    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Evaluate 1/N sum (eps_j-<eps_j>)(eps_j-<eps_j>)^T e_i **")'); call goPr
    write (gol,'("")'); call goPr
    
    ! flag:
    call ReadRc( rcF, 'nmc.sample-B.correlation', correlation, status )
    IF_NOT_OK_RETURN(status=1)
    ! set label:
    if ( correlation ) then
      cv_label = 'correlation'
      cv_id    = 'C'
      cv_units = '1'
    else
      cv_label = 'covariance'
      cv_id    = 'B'
      cv_units = 'tracer*tracer'
    end if
    
    ! flag:
    call ReadRc( rcF, 'nmc.sample-B.levtr', levtr, status )
    IF_NOT_OK_RETURN(status=1)

    ! start time:
    tday = self%tday1
    ! set flag:
    setup = .false.
    ! init current filenames:
    eps_filename = ''
    ! time loop:
    do
      ! info ...
      call wrtgol( 'day ', tday ); call goPr
      
      ! ~ input samples

      ! input file:
      filename = trim(self%eps_samples_filename)
      ! replace keys:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
      IF_NOT_OK_RETURN(status=1)
      
      ! new ?
      if ( trim(filename) /= trim(eps_filename) ) then
      
        ! open?
        if ( len_trim(eps_filename) > 0 ) then
          ! close:
          call eps_file%Close( status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! copy:
        eps_filename = trim(filename)

        ! open file with eps samples:
        call eps_file%Open( trim(filename), status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init counter:
        eps_isample = 1
        
      else
      
        ! next sample for already open file:
        eps_isample = eps_isample + 1
        
      end if
      
      ! ~ setup

      ! initial setup?
      if ( .not. setup ) then

        ! create longitude coordinate:
        call lon_coor%Init( 'longitude', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lon_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lon_coor%Get_Dim( status, n=nlon )
        IF_NOT_OK_RETURN(status=1)

        ! create latitude coordinate:
        call lat_coor%Init( 'latitude', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lat_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lat_coor%Get_Dim( status, n=nlat )
        IF_NOT_OK_RETURN(status=1)

        ! create level coordinate:
        call lev_coor%Init( 'lev', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lev_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lev_coor%Get_Dim( status, n=nlev )
        IF_NOT_OK_RETURN(status=1)

        ! create hour coordinate:
        call ctime_coor%Init( 'time', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file, here includes climatology bounds:
        call ctime_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call ctime_coor%Get_Dim( status, n=ntime )
        IF_NOT_OK_RETURN(status=1)

        ! create sample coordinate:
        call eps_sample_coor%Init( 'sample', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call eps_sample_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! get name of pressure variable:
        call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( eps__varids(self%nvar) )
        ! loop over variables:
        do ivar = 1, self%nvar
          ! get id:
          call eps_file%Get_VarID( trim(self%var(ivar)%name), eps__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do

        ! storage for 3D fields:
        allocate( ps(nlon,nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( eps(nlon,nlat,nlev,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        
        ! storage for meta data:
        allocate( eps__units(self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ points

        ! loop over points:
        do ipoint = 1, self%npoint
          ! nearby index:
          call lon_coor%Get_Index( self%point(ipoint)%lon, self%point(ipoint)%ilon, status )
          IF_NOT_OK_RETURN(status=1)
          ! nearby index:
          call lat_coor%Get_Index( self%point(ipoint)%lat, self%point(ipoint)%ilat, status )
          IF_NOT_OK_RETURN(status=1)
          ! surface (top-down order!)
          self%point(ipoint)%ilev = nlev
        end do

        ! storage for eps profiles:
        allocate( eps_profiles(self%npoint,nlev,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ statistics and covariance

        ! for average pressure:
        allocate( sampstat_ps(ntime), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over hours:
        do itime = 1, ntime
            ! init mean/sigma computation over 2D field:
          call sampstat_ps(itime)%Init( (/nlon,nlat/), status )
          IF_NOT_OK_RETURN(status=1)
        end do  ! hour

        ! for sigma:
        allocate( sampstat_eps(ntime,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over hours:
        do itime = 1, ntime
          ! loop over variables:
          do ivar = 1, self%nvar
            ! init mean/sigma computation over 3D field:
            call sampstat_eps(itime,ivar)%Init( (/nlon,nlat,nlev/), status )
            IF_NOT_OK_RETURN(status=1)
          end do  ! var
        end do  ! hour

        ! covariance of 3D fields with point values:
        allocate( sampcovr_eps(ntime,self%nvar,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over hours:
        do itime = 1, ntime
          ! loop over variables:
          do ivar = 1, self%nvar
            do ivar2 = 1, self%nvar
              ! init covariance between 3D field and point values:
              call sampcovr_eps(itime,ivar,ivar2)%Init( (/nlon,nlat,nlev/), self%npoint, status )
              IF_NOT_OK_RETURN(status=1)
            end do
          end do  ! var
        end do  ! hour

        ! apply?
        if ( levtr ) then
          ! covariance between tracers/levels :
          allocate( sampcovr_Bzs(ntime,self%npoint), stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop over hours:
          do itime = 1, ntime
            ! loop over points:
            do ipoint = 1, self%npoint
              ! init covariance between pair of 2D fields
              call sampcovr_Bzs(itime,ipoint)%Init( (/nlev,self%nvar/), (/nlev,self%nvar/), status )
              IF_NOT_OK_RETURN(status=1)
            end do ! points
          end do  ! hour
        end if  ! levtr

        ! ~

        ! reset flag:
        setup = .true.

      end if  ! not setup yet

      ! ~
      
      ! get sample info:
      call eps_sample_coor%Get_Value( eps_isample, status, t=t )
      IF_NOT_OK_RETURN(status=1)

      ! loop over hours:
      do itime = 1, ntime

        ! info ...
        write (gol,'("  time ",i4," ...")') itime; call goPr
        
        !! testing ...
        !if ( itime /= 2 ) then
        !  write (gol,'("    WARNING - skip this time")'); call goErr
        !  cycle
        !end if
    
        ! read sample:
        call eps_file%Get_Sample2D( vname_ps, itime, eps_isample, ps, ps__units, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! update statistics:
        call sampstat_ps(itime)%AddSample( ps, status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%nvar
    
          ! info ...
          write (gol,'("    variable ",i0," ",a)') itime, trim(self%var(ivar)%name); call goPr
        
          ! read sample:
          call eps_file%Get_Sample3D( trim(self%var(ivar)%name), itime, eps_isample, &
                                        eps(:,:,:,ivar), eps__units(ivar), status )
          IF_NOT_OK_RETURN(status=1)
          
          ! update statistics:
          call sampstat_eps(itime,ivar)%AddSample( eps(:,:,:,ivar), status )
          IF_NOT_OK_RETURN(status=1)
          
        end do ! var

        ! collect sample differences at point locations:
        do ipoint = 1, self%npoint

          ! grid cell indices of point:
          ilon = self%point(ipoint)%ilon
          ilat = self%point(ipoint)%ilat
          ilev = self%point(ipoint)%ilev

          ! profiles:
          eps_profiles(ipoint,:,:) = eps(ilon,ilat,:,:)
          
        end do ! points

        ! loop over variables:
        do ivar = 1, self%nvar
          do ivar2 = 1, self%nvar

            ! update sample covariances:
            call sampcovr_eps(itime,ivar,ivar2)%AddSample( eps(:,:,:,ivar), eps_profiles(:,ilev,ivar2), status )
            IF_NOT_OK_RETURN(status=1)

          end do ! var2
        end do ! var

        ! level/tracer covar?
        if ( levtr ) then
          ! loop over points:
          do ipoint = 1, self%npoint
            ! update sample covariances:
            call sampcovr_Bzs(itime,ipoint)%AddSample( eps_profiles(ipoint,:,:), eps_profiles(ipoint,:,:), status )
            IF_NOT_OK_RETURN(status=1)
          end do ! points
        end if ! levtr
        
      end do  ! hour

      ! next:
      tday = tday + IncrDate(day=1)
      ! end ?
      if ( tday > self%tday2 ) exit
      
    end do  ! samples
    
    ! info ...
    write (gol,'("ok")'); call goPr
    
    ! ~ output
    
    ! info ...
    write (gol,'("create B_point file ...")'); call goPr
    
    ! target file:
    call ReadRc( rcF, 'nmc.sample-B.filename', B_point_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! create new result file:
    call out_file%Create( trim(B_point_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("evaluate ",a," from `",a,"`")') &
            trim(cv_label), trim(self%BB_filename)
    ! add:
    call out_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define coordinates in output file:
    call lon_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define level coordinate in output file:
    call lev_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! set size:
    call tracer_coor%Set_Dim( status, n=self%nvar )
    IF_NOT_OK_RETURN(status=1)
    ! set attributes:
    call tracer_coor%Set_Attrs( status, long_name='tracer' )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do ivar = 1, self%nvar
      ! fill tracer name:
      call tracer_coor%Set_Value( ivar, status, value=trim(self%var(ivar)%name) )
      IF_NOT_OK_RETURN(status=1)
    end do  ! var

    ! define tracer coordinate in output file:
    call tracer_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call out_file%Def_Field2D_Series( trim(vname_ps), trim(units_p), &
                                  lon_coor, lat_coor, ctime_coor, &
                                  ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call out_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( S__varids(self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bx__varids(self%npoint,self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over tracers:
    do ivar = 1, self%nvar

      ! copy:
      tracer_name  = trim(self%var(ivar)%name)
      tracer_units = trim(eps__units(ivar))
      
      ! variable:
      write (vname,'("S_",a)') trim(tracer_name)
      ! define:
      call out_file%Def_Field3D_Series( trim(vname), trim(tracer_units), &
                              lon_coor, lat_coor, lev_coor, ctime_coor, &
                              S__varids(ivar), status )
      IF_NOT_OK_RETURN(status=1)
      call out_file%Put_Att( S__varids(ivar), 'long_name', &
                     trim(tracer_name)//' grid point standard deviation', status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over points:
      do ipoint = 1, self%npoint

        ! variable:
        write (vname,'(a,"_",a,"_",a)') trim(cv_id), trim(tracer_name), trim(self%point(ipoint)%name)
        ! define:
        call out_file%Def_Sample3D( trim(vname), trim(cv_units), &
                                lon_coor, lat_coor, lev_coor, ctime_coor, tracer_coor, &
                                Bx__varids(ipoint,ivar), status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'long_name', &
                       trim(cv_label)//' with '//trim(tracer_name)//' in '//trim(self%point(ipoint)%name), status )
        IF_NOT_OK_RETURN(status=1)
        ! write as attributes:
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'lon', self%point(ipoint)%lon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'lat', self%point(ipoint)%lat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'ilon', self%point(ipoint)%ilon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'ilat', self%point(ipoint)%ilat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'ilev', self%point(ipoint)%ilev, status )
        IF_NOT_OK_RETURN(status=1)

      end do ! points

    end do ! tracers
    
    ! evaluate level/tracer covariance?
    if ( levtr ) then

      ! storage for result:
      allocate( Bzs(nlev,nlev,self%nvar,self%nvar), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! idem in different order:
      allocate( Vzs(nlev,self%nvar,nlev,self%nvar), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! output:
      allocate( Bzs__varids(self%npoint), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over points:
      do ipoint = 1, self%npoint

        ! variable:
        write (vname,'(a,"zs_",a)') trim(cv_id), trim(self%point(ipoint)%name)
        ! define:
        call out_file%Def_ACovar1( trim(vname), trim(cv_units), &
                                    lev_coor, tracer_coor, ctime_coor, &
                                    Bzs__varids(ipoint), status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'long_name', &
                       'vertical '//trim(cv_label)//' in '//trim(self%point(ipoint)%name), status )
        IF_NOT_OK_RETURN(status=1)
        ! write as attributes:
        call out_file%Put_Att( Bzs__varids(ipoint), 'lon', self%point(ipoint)%lon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'lat', self%point(ipoint)%lat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'ilon', self%point(ipoint)%ilon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'ilat', self%point(ipoint)%ilat, status )
        IF_NOT_OK_RETURN(status=1)
        
      end do  ! points
      
    end if  ! level/tracer covar

    ! end:
    call out_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lon_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lev_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call ctime_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call tracer_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! temporary storage:
    allocate( S(nlon,nlat,nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bx(nlon,nlat,nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over hours:
    do itime = 1, ntime

      ! extract:
      call sampstat_ps(itime)%Get2d( status, mean=ps )
      IF_NOT_OK_RETURN(status=1)
      ! write:
      call out_file%Put_Field2D_Series( ps__varid, itime, ps, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! ~ std.dev fields

      ! loop over variables:
      do ivar = 1, self%nvar
      
        ! extract 3D std.dev. field:
        call sampstat_eps(itime,ivar)%Get3D( status, stdv=S )
        IF_NOT_OK_RETURN(status=1)
        ! write:
        call out_file%Put_Field3d_Series( S__varids(ivar), itime, S, status )
        IF_NOT_OK_RETURN(status=1)  
          
      end do  ! var

      ! ~ 3D covariance with point

      ! loop over variables:
      do ivar = 1, self%nvar

        ! second tracer dim:
        do ivar2 = 1, self%nvar

          ! loop over points:
          do ipoint = 1, self%npoint
        
            ! correlation or covariance?
            if ( correlation ) then
              ! extract 3D correlation field with point:
              call sampcovr_eps(itime,ivar,ivar2)%Get3D( ipoint, status, corr=Bx )
              IF_NOT_OK_RETURN(status=1)
            else
              ! extract 3D covariance field with point:
              call sampcovr_eps(itime,ivar,ivar2)%Get3D( ipoint, status, covr=Bx )
              IF_NOT_OK_RETURN(status=1)
            end if  ! correlation or covariance
            ! write:
            call out_file%Put_Sample3D( Bx__varids(ipoint,ivar), itime, ivar2, Bx, status )
            IF_NOT_OK_RETURN(status=1)
          
          end do ! points

        end do  ! var2
          
      end do  ! var
 
      ! ~ level/tracer covariance

      ! evaluate full level/tracer covar?
      if ( levtr ) then 

        ! loop over points:
        do ipoint = 1, self%npoint
        
          ! correlation or covariance?
          if ( correlation ) then
            ! extract sample result:
            call sampcovr_Bzs(itime,ipoint)%GetXY( status, corr=Vzs )
            IF_NOT_OK_RETURN(status=1)
          else
            ! extract sample result:
            call sampcovr_Bzs(itime,ipoint)%GetXY( status, covr=Vzs )
            IF_NOT_OK_RETURN(status=1)
          end if
          ! re-order:
          do ivar = 1, self%nvar
            do ivar2 = 1, self%nvar
              Bzs(:,:,ivar,ivar2) = Vzs(:,ivar,:,ivar2)
            end do
          end do

          ! write sample:
          call out_file%Put_ACovar1( Bzs__varids(ipoint), itime, Bzs, status  )
          IF_NOT_OK_RETURN(status=1)

        end do ! points

      end if ! level/tracer covariance

      ! ~

    end do ! times
    
    ! clear:
    deallocate( S, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( Bx, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( Bx__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! evaluate level/tracer covariance?
    if ( levtr ) then
      ! clear:
      deallocate( Bzs, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! clear:
      deallocate( Vzs, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! clear:
      deallocate( Bzs__varids, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! level/tracer covar
      
    ! clear:
    call tracer_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close output:
    call out_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ clear stats
    
    ! loop over time:
    do itime = 1, ntime
      ! done:
      call sampstat_ps(itime)%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! clear:
    deallocate( sampstat_ps, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over time and var:
    do itime = 1, ntime
      do ivar = 1, self%nvar
        ! done:
        call sampstat_eps(itime,ivar)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do
    end do
    ! clear:
    deallocate( sampstat_eps, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over time and var:
    do itime = 1, ntime
      do ivar = 1, self%nvar
        do ivar2 = 1, self%nvar
          ! done:
          call sampcovr_eps(itime,ivar,ivar2)%Done( status )
          IF_NOT_OK_RETURN(status=1)
        end do ! ivar2
      end do ! ivar
    end do ! itime
    ! clear:
    deallocate( sampcovr_eps, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( eps_profiles, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ close input
    
    ! clear:
    deallocate( eps__units, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( eps, stat=status )
    IF_NOT_OK_RETURN(status=1)
        
    ! clear:
    deallocate( eps__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    call lon_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call eps_sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! close input:
    call eps_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ 

    ! ok
    status = 0
  
  end subroutine ModelRuns_Evaluate_Eps


  ! ***
  
  
  !
  ! Compute eta samples:
  ! 
  !    eta_comp(sample,hour,lev,lat,lon)  =  eps_comp(sample,hour,lev,lat,lon) / eps_comp_stdv(hour,lev,lat,lon)
  !

  subroutine ModelRuns_Eta_Samples( self, rcF, status )
  
    use GO             , only : TDate, IncrDate, operator(+), operator(>), wrtgol
    use GO             , only : TrcFile, ReadRc
    use GO             , only : goReplace
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : RealCoordinate, HybrideLevelCoordinate, TimeCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Eta_Samples'
    
    ! --- local ----------------------------------

    type(TDate)                      ::  tday, t
    character(len=1024)              ::  filename
    character(len=1024)              ::  eps_filename
    character(len=1024)              ::  eta_filename
    integer                          ::  eps_isample
    integer                          ::  eta_isample
    logical                          ::  setup
    type(T_NMC_Output)               ::  eps_stats_file
    type(T_NMC_Output)               ::  eps_file
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(TimeCoordinate)             ::  ctime_coor
    type(TimeCoordinate)             ::  eps_sample_coor
    type(TimeCoordinate)             ::  eta_sample_coor
    integer                          ::  nlon, nlat, nlev
    integer                          ::  ntime
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    integer, allocatable             ::  eps__varids(:)   ! (nvar)
    real, allocatable                ::  ps(:,:)
    real, allocatable                ::  eps(:,:,:)
    real, allocatable                ::  eta(:,:,:)
    character(len=32)                ::  vname, units
    real, allocatable                ::  eps_mean(:,:,:,:,:)
    real, allocatable                ::  eps_stdv(:,:,:,:,:)
    integer, allocatable             ::  fixed(:,:)  ! (nlev,nvar)
    character(len=32)                ::  fixed__units
    integer, allocatable             ::  fixed__varids(:)
    integer                          ::  ivar
    integer                          ::  itime
    
    type(T_NMC_Output)               ::  eta_file
    character(len=1024)              ::  msg
    integer                          ::  ps__varid
    integer, allocatable             ::  eta__varids(:)   ! (nvar)
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Compute eta samples **")'); call goPr
    write (gol,'("")'); call goPr

    ! start time:
    tday = self%tday1
    ! set flag:
    setup = .false.
    ! init current filenames:
    eps_filename = ''
    eta_filename = ''
    ! time loop:
    do
      ! info ...
      call wrtgol( 'day ', tday ); call goPr
      
      ! ~ input samples

      ! input file:
      filename = trim(self%eps_samples_filename)
      ! replace keys:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
      IF_NOT_OK_RETURN(status=1)
      
      ! new ?
      if ( trim(filename) /= trim(eps_filename) ) then
      
        ! open?
        if ( len_trim(eps_filename) > 0 ) then
          ! close:
          call eps_file%Close( status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! copy:
        eps_filename = trim(filename)

        ! open file with eps samples:
        call eps_file%Open( trim(filename), status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init counter:
        eps_isample = 1
        
      else
      
        ! next sample for already open file:
        eps_isample = eps_isample + 1
        
      end if
      
      ! ~ setup

      ! initial setup?
      if ( .not. setup ) then

        ! create longitude coordinate:
        call lon_coor%Init( 'longitude', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lon_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lon_coor%Get_Dim( status, n=nlon )
        IF_NOT_OK_RETURN(status=1)

        ! create latitude coordinate:
        call lat_coor%Init( 'latitude', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lat_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lat_coor%Get_Dim( status, n=nlat )
        IF_NOT_OK_RETURN(status=1)

        ! create level coordinate:
        call lev_coor%Init( 'lev', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lev_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lev_coor%Get_Dim( status, n=nlev )
        IF_NOT_OK_RETURN(status=1)

        ! create hour coordinate:
        call ctime_coor%Init( 'time', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file, here includes climatology bounds:
        call ctime_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call ctime_coor%Get_Dim( status, n=ntime )
        IF_NOT_OK_RETURN(status=1)

        ! create sample coordinate:
        call eps_sample_coor%Init( 'sample', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call eps_sample_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! get name of pressure variable:
        call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( eps__varids(self%nvar) )
        ! loop over variables:
        do ivar = 1, self%nvar
          ! get id:
          call eps_file%Get_VarID( trim(self%var(ivar)%name), eps__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do

        ! storage:
        allocate( ps(nlon,nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( eps(nlon,nlat,nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ eps statistics:

        ! open file with eps samples:
        call eps_stats_file%Open( trim(self%eps_stats_filename), status )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( eps_mean(nlon,nlat,nlev,ntime,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( eps_stdv(nlon,nlat,nlev,ntime,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( fixed(nlev,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%nvar
          !
          ! variable:
          write (vname,'(a,"_mean")') trim(self%var(ivar)%name)
          ! loop over hours:
          do itime = 1, ntime
            ! read:
            call eps_stats_file%Get_Field3D_Series( trim(vname), itime, eps_mean(:,:,:,itime,ivar), units, status )
            IF_NOT_OK_RETURN(status=1)
          end do ! time
          !
          ! variable:
          write (vname,'(a,"_stdv")') trim(self%var(ivar)%name)
          ! loop over times:
          do itime = 1, ntime
            ! read:
            call eps_stats_file%Get_Field3D_Series( trim(vname), itime, eps_stdv(:,:,:,itime,ivar), units, status )
            IF_NOT_OK_RETURN(status=1)
          end do ! time
          !
          ! variable:
          write (vname,'(a,"_fixed")') trim(self%var(ivar)%name)
          call eps_stats_file%Get_IField1D( trim(vname), fixed(:,ivar), fixed__units, status )
          IF_NOT_OK_RETURN(status=1)
          !
        end do ! var

        ! close:
        call eps_stats_file%Close( status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ output storage

        ! access:
        allocate( eta__varids(self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( fixed__varids(self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! ~

        ! reset flag:
        setup = .true.

      end if  ! not setup yet

      ! ~ output
        
      ! input file:
      filename = trim(self%eta_samples_filename)
      ! replace keys:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
      IF_NOT_OK_RETURN(status=1)
      
      ! new ?
      if ( trim(filename) /= trim(eta_filename) ) then
      
        ! open?
        if ( len_trim(eta_filename) > 0 ) then
          ! write times of samples:
          call eta_sample_coor%Write( eta_file, status )
          IF_NOT_OK_RETURN(status=1)
          ! close:
          call eta_file%Close( status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! copy:
        eta_filename = trim(filename)

        ! create new result file:
        call eta_file%Create( trim(filename), status )
        IF_NOT_OK_RETURN(status=1)

        ! description:
        write (msg,'("eta samples from `",a,"` and `",a,"`")') &
                 trim(self%eps_samples_filename), trim(self%eps_stats_filename)
        ! add:
        call eta_file%Extend_History( trim(msg), status )
        IF_NOT_OK_RETURN(status=1)

        ! define coordinates in output file:
        call lon_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lev_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call ctime_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! output sample coor:
        call eta_sample_coor%Init( 'sample', status )
        IF_NOT_OK_RETURN(status=1)
        ! set unlimited size:
        call eta_sample_coor%Set_Dim( status, unlimited=.true. )
        IF_NOT_OK_RETURN(status=1)
        ! fill attributes:
        call eta_sample_coor%Set_Values( status, units_step='days', units_t0=self%tday1 )
        IF_NOT_OK_RETURN(status=1)
        ! define sample coordinate in output file:
        call eta_sample_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init counter:
        eta_isample = 1

        ! define:
        call eta_file%Def_Sample2D( trim(vname_ps), trim(units_p), &
                                      lon_coor, lat_coor, ctime_coor, eta_sample_coor, &
                                      ps__varid, status )
        IF_NOT_OK_RETURN(status=1)
        call eta_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%nvar
          ! define, result will be unit-less:
          call eta_file%Def_Sample3D( trim(self%var(ivar)%name), '1', &
                                        lon_coor, lat_coor, lev_coor, ctime_coor, eta_sample_coor, &
                                        eta__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
          call eta_file%Put_Att( eta__varids(ivar), 'long_name', trim(self%var(ivar)%name)//' (eps-eps_mean)/eps_stdv', status )
          IF_NOT_OK_RETURN(status=1)

          ! define:
          call eta_file%Def_IField1D( trim(self%var(ivar)%name)//'_fixed', trim(fixed__units), &
                                        lev_coor, fixed__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
          call eta_file%Put_Att( fixed__varids(ivar), 'long_name', trim(self%var(ivar)%name)//' fixed layer (1=fixed, 0=dynamic)', status )
          IF_NOT_OK_RETURN(status=1)

        end do  ! var

        ! end:
        call eta_file%EndDef( status )
        IF_NOT_OK_RETURN(status=1)

        ! write coordinate variables:
        call lon_coor%Write( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_coor%Write( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lev_coor%Write( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call ctime_coor%Write( eta_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%nvar
          ! write:
          call eta_file%Put_IField1D( fixed__varids(ivar), fixed(:,ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do
        
      else
      
        ! next record in already opened file:
        eta_isample = eta_isample + 1
        
      end if  ! new output file

      ! ~
      
      ! get sample info:
      call eps_sample_coor%Get_Value( eps_isample, status, t=t )
      IF_NOT_OK_RETURN(status=1)
      ! copy to output coordinate:
      call eta_sample_coor%Set_Value( eta_isample, status, t=t )
      IF_NOT_OK_RETURN(status=1)

      ! loop over hours:
      do itime = 1, self%nhour

        ! info ...
        write (gol,'("  time ",i4," ...")') itime; call goPr
    
        ! read sample:
        call eps_file%Get_Sample2D( vname_ps, itime, eps_isample, ps, units, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! write sample:
        call eta_file%Put_Sample2D( ps__varid, itime, eta_isample, ps, status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%nvar
    
          ! info ...
          write (gol,'("    variable ",i0," ",a)') itime, trim(self%var(ivar)%name); call goPr
        
          ! read sample:
          call eps_file%Get_Sample3D( trim(self%var(ivar)%name), itime, eps_isample, eps, units, status )
          IF_NOT_OK_RETURN(status=1)
          
          ! convert:
          where ( eps_stdv(:,:,:,itime,ivar) > 0.0 )
            eta = ( eps - eps_mean(:,:,:,itime,ivar) ) / eps_stdv(:,:,:,itime,ivar)
          elsewhere
            eta = 0.0
          end where
          
          ! write sample:
          call eta_file%Put_Sample3D( eta__varids(ivar), itime, eta_isample, eta, status )
          IF_NOT_OK_RETURN(status=1)
          
        end do ! var
        
      end do  ! hour

      ! next:
      tday = tday + IncrDate(day=1)
      ! end ?
      if ( tday > self%tday2 ) exit
      
    end do  ! samples
    
    ! info ...
    write (gol,'("ok")'); call goPr

    ! ~ close output
    
    ! write times of samples:
    call eta_sample_coor%Write( eta_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! close ouptut:
    call eta_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! done:
    call eta_sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( eta__varids )
    deallocate( fixed__varids )
    
    ! ~ clear stats

    ! clear stats:
    deallocate( eps_mean, eps_stdv )
    deallocate( fixed )

    ! ~ close input
    
    ! clear:
    deallocate( ps )
    deallocate( eps )

    ! clear:
    deallocate( eps__varids )
      
    ! clear:
    call lon_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call eps_sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! close input:
    call eps_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~

    ! ok
    status = 0
  
  end subroutine ModelRuns_Eta_Samples


  ! ***
  
  
  !
  ! Compute Fourrier transforms of eta samples:
  ! 
  !    eta_f = F eta
  !

  subroutine ModelRuns_Eta_f_Samples( self, rcF, status )
  
    use GO             , only : TDate, IncrDate, operator(+), operator(>), wrtgol
    use GO             , only : TrcFile, ReadRc
    use GO             , only : goReplace
    use GO             , only : T_SampStat
    use EMEP_NMC_Output, only : T_NMC_Output
    use UFFTW3         , only : T_UFFTW3_2d
    use C3PO           , only : Dimension
    use C3PO           , only : RealCoordinate, HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : IntegerCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Eta_f_Samples'
    
    ! --- local ----------------------------------

    type(TDate)                      ::  tday, t
    character(len=1024)              ::  filename
    logical                          ::  setup

    character(len=1024)              ::  eta_filename
    integer                          ::  eta_isample
    type(T_NMC_Output)               ::  eta_file
    type(TimeCoordinate)             ::  eta_sample_coor
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(TimeCoordinate)             ::  ctime_coor
    integer                          ::  nlon, nlat, nlev
    integer                          ::  ntime
    real                             ::  p0
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    integer                          ::  ivar
    integer, allocatable             ::  eta__varids(:)   ! (nvar)
    real, allocatable                ::  eta(:,:,:)
    character(len=32)                ::  units
    integer, allocatable             ::  fixed(:,:)  ! (nlev,nvar)
    character(len=32)                ::  fixed__units
    integer, allocatable             ::  fixed__varids(:)  ! (nvar)
    
    character(len=64)                ::  rckey
    integer                          ::  nlon_ex, nlat_ex
    real, allocatable                ::  eta_ex(:,:,:)

    type(T_UFFTW3_2d)                ::  ufft
    integer                          ::  nlon_f, nlat_f
    integer, allocatable             ::  lons_f(:), lats_f(:)
    real, allocatable                ::  kstar(:,:)
    real                             ::  kstar__fill_value
    integer                          ::  kstar__varid
    character(len=256)               ::  kstar__formula
    real, allocatable                ::  theta(:,:)
    real                             ::  theta__fill_value
    integer                          ::  theta__varid
    integer, allocatable             ::  hcn(:,:)
    integer                          ::  hcn__varid
    type(IntegerCoordinate)          ::  lon_f_coor, lat_f_coor
    type(Dimension)                  ::  cplx
    type(Dimension)                  ::  lon_ex_dim, lat_ex_dim

    character(len=1024)              ::  eta_f_filename
    integer                          ::  eta_f_isample
    type(T_NMC_Output)               ::  eta_f_file
    type(TimeCoordinate)             ::  eta_f_sample_coor
    character(len=1024)              ::  msg
    complex, allocatable             ::  eta_f(:,:,:)   ! (nlon_f,nlat_f,nlev)
    real                             ::  eta_f__fill_value
    integer, allocatable             ::  eta_f__varids(:)   ! (nvar)
    real, allocatable                ::  ps(:,:)
    integer                          ::  ps__varid
    integer                          ::  itime
    integer                          ::  ilev
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Compute eta spectral transforms **")'); call goPr
    write (gol,'("")'); call goPr

    ! start time:
    tday = self%tday1
    ! set flag:
    setup = .false.
    ! init current filenames:
    eta_filename = ''
    eta_f_filename = ''
    ! time loop:
    do
      ! info ...
      call wrtgol( 'day ', tday ); call goPr
    
      ! ~ input samples

      ! input file:
      filename = trim(self%eta_samples_filename)
      ! replace keys:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
      IF_NOT_OK_RETURN(status=1)
      
      ! new ?
      if ( trim(filename) /= trim(eta_filename) ) then
      
        ! open?
        if ( len_trim(eta_filename) > 0 ) then
          ! close:
          call eta_file%Close( status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! copy:
        eta_filename = trim(filename)
          
        ! open file with eta samples:
        call eta_file%Open( trim(filename), status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init counter:
        eta_isample = 1
        
      else
      
        ! next record in already opened file:
        eta_isample = eta_isample + 1
        
      end if
      
      ! ~ setup

      ! initial setup?
      if ( .not. setup ) then

        ! create longitude coordinate:
        call lon_coor%Init( 'longitude', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lon_coor%Read( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lon_coor%Get_Dim( status, n=nlon )
        IF_NOT_OK_RETURN(status=1)

        ! create latitude coordinate:
        call lat_coor%Init( 'latitude', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lat_coor%Read( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lat_coor%Get_Dim( status, n=nlat )
        IF_NOT_OK_RETURN(status=1)

        ! create level coordinate:
        call lev_coor%Init( 'lev', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lev_coor%Read( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lev_coor%Get_Dim( status, n=nlev )
        IF_NOT_OK_RETURN(status=1)

        ! create (climatological) time coordinate:
        call ctime_coor%Init( 'time', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call ctime_coor%Read( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call ctime_coor%Get_Dim( status, n=ntime )
        IF_NOT_OK_RETURN(status=1)

        ! create sample coordinate:
        call eta_sample_coor%Init( 'sample', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call eta_sample_coor%Read( eta_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! get name of pressure variable:
        call lev_coor%Get_Values( status, vname_ps=vname_ps, &
                                   units_p=units_p, p0=p0 )
        IF_NOT_OK_RETURN(status=1)

        ! storage for variable id's, needed to copy units:
        allocate( eta__varids(self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over variables:
        do ivar = 1, self%nvar
          ! get id:
          call eta_file%Get_VarID( trim(self%var(ivar)%name), &
                                     eta__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do

        ! storage:
        allocate( eta(nlon,nlat,nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( fixed(nlev,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop:
        do ivar = 1, self%nvar
          ! read:
          call eta_file%Get_IField1D( trim(self%var(ivar)%name)//'_fixed', &
                                        fixed(:,ivar), fixed__units, status )
          IF_NOT_OK_RETURN(status=1)
        end do

        ! ~ fft

        ! padded dimensions, read from rcfile:
        !~ lon
        write (rckey,'("nmc.spectral.padding.",i3.3)') nlon
        call ReadRc( rcf, trim(rckey), nlon_ex, status )
        IF_NOT_OK_RETURN(status=1)
        !~ lat
        write (rckey,'("nmc.spectral.padding.",i3.3)') nlat
        call ReadRc( rcf, trim(rckey), nlat_ex, status )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( eta_ex(nlon_ex,nlat_ex,nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! init transform for extended grid:
        call ufft%Init( nlon_ex, nlat_ex, status, unitairy=.true. )
        IF_NOT_OK_RETURN(status=1)

        ! spectral size in half-complex representation:
        call ufft%Get( status, nxf=nlon_f, nyf=nlat_f )
        IF_NOT_OK_RETURN(status=1)

        ! storage for spectral indices:
        allocate( lons_f(nlon_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( lats_f(nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! storage for wave number:
        allocate( kstar(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! storage for spectral angle:
        allocate( theta(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! storage for complex counter:
        allocate( hcn(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! get spectral indices etc:
        call ufft%Get( status, xf=lons_f, yf=lats_f,&
                          kstar=kstar, kstar_formula=kstar__formula, &
                          theta=theta, hcn=hcn )
        IF_NOT_OK_RETURN(status=1)

        ! create longitude coordinate:
        call lon_f_coor%Init( 'lon_f', status )
        IF_NOT_OK_RETURN(status=1)
        ! set size:
        call lon_f_coor%Set_Dim( status, n=nlon_f )
        IF_NOT_OK_RETURN(status=1)
        ! set attributes:
        call lon_f_coor%Set_Attrs( status, units='1', long_name='longitude frequency m' )
        IF_NOT_OK_RETURN(status=1)
        ! fill values:
        call lon_f_coor%Set_Values( status, values=lons_f )
        IF_NOT_OK_RETURN(status=1)

        ! create latitude coordinate:
        call lat_f_coor%Init( 'lat_f', status )
        IF_NOT_OK_RETURN(status=1)
        ! set size:
        call lat_f_coor%Set_Dim( status, n=nlat_f )
        IF_NOT_OK_RETURN(status=1)
        ! set attributes:
        call lat_f_coor%Set_Attrs( status, units='1', long_name='latitude frequency n' )
        IF_NOT_OK_RETURN(status=1)
        ! fill values:
        call lat_f_coor%Set_Values( status, values=lats_f )
        IF_NOT_OK_RETURN(status=1)

        ! create extended longitude dimension:
        call lon_ex_dim%Init( 'lon_ex', status )
        IF_NOT_OK_RETURN(status=1)
        ! set size:
        call lon_ex_dim%Set_Dim( status, n=nlon_ex )
        IF_NOT_OK_RETURN(status=1)

        ! create extended latitude dimension:
        call lat_ex_dim%Init( 'lat_ex', status )
        IF_NOT_OK_RETURN(status=1)
        ! set size:
        call lat_ex_dim%Set_Dim( status, n=nlat_ex )
        IF_NOT_OK_RETURN(status=1)

        ! complex dimension:
        call cplx%Init( 'cmplx', status )
        IF_NOT_OK_RETURN(status=1)
        ! set size for (real,imag):
        call cplx%Set_Dim( status, n=2 )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( eta_f(nlon_f,nlat_f,nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( ps(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! access:
        allocate( eta_f__varids(self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( fixed__varids(self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! reset flag:
        setup = .true.

      end if  ! initial setup


      ! ~ output
        
      ! input file:
      filename = trim(self%eta_f_samples_filename)
      ! replace keys:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
      IF_NOT_OK_RETURN(status=1)
      
      ! new ?
      if ( trim(filename) /= trim(eta_f_filename) ) then
      
        ! open?
        if ( len_trim(eta_f_filename) > 0 ) then
          ! write times of samples:
          call eta_f_sample_coor%Write( eta_f_file, status )
          IF_NOT_OK_RETURN(status=1)
          ! close:
          call eta_f_file%Close( status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! copy:
        eta_f_filename = trim(filename)
    
        ! create new result file:
        call eta_f_file%Create( trim(filename), status )
        IF_NOT_OK_RETURN(status=1)

        ! description:
        write (msg,'("spectral transforms from `",a,"`")') trim(self%eta_samples_filename)
        ! add:
        call eta_f_file%Extend_History( trim(msg), status )
        IF_NOT_OK_RETURN(status=1)

        ! define complex dimension in output file:
        call cplx%Def( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! define coordinates in output file:
        call lon_f_coor%Def( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_f_coor%Def( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! define dimensions in output file:
        call lon_ex_dim%Def( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_ex_dim%Def( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! define level coordinate in output file:
        call lev_coor%Def( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! define hour coordinate in output file:
        call ctime_coor%Def( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! output sample coor:
        call eta_f_sample_coor%Init( 'sample', status )
        IF_NOT_OK_RETURN(status=1)
        ! set unlimited size:
        call eta_f_sample_coor%Set_Dim( status, unlimited=.true. )
        IF_NOT_OK_RETURN(status=1)
        ! fill attributes:
        call eta_f_sample_coor%Set_Values( status, units_step='days', units_t0=self%tday1 )
        IF_NOT_OK_RETURN(status=1)
        ! define sample coordinate in output file:
        call eta_f_sample_coor%Def( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init output record counter:
        eta_f_isample = 1

        ! define:
        call eta_f_file%Def_Field2D( 'kstar', '1', lon_f_coor, lat_f_coor, &
                                       kstar__varid, status, &
                                       fill_value=kstar__fill_value )
        IF_NOT_OK_RETURN(status=1)
        call eta_f_file%Put_Att( kstar__varid, 'long_name', 'wave number', status )
        IF_NOT_OK_RETURN(status=1)
        call eta_f_file%Put_Att( kstar__varid, 'formula', trim(kstar__formula), status )
        IF_NOT_OK_RETURN(status=1)

        ! define:
        call eta_f_file%Def_Field2D( 'theta', 'rad', lon_f_coor, lat_f_coor, &
                                       theta__varid, status, &
                                       fill_value=theta__fill_value )
        IF_NOT_OK_RETURN(status=1)
        call eta_f_file%Put_Att( theta__varid, 'long_name', 'spectral angle', status )
        IF_NOT_OK_RETURN(status=1)

        ! define:
        call eta_f_file%Def_IField2D( 'hcn', '1', lon_f_coor, lat_f_coor, &
                                       hcn__varid, status )
        IF_NOT_OK_RETURN(status=1)
        call eta_f_file%Put_Att( hcn__varid, 'long_name', 'number of represented complex coefficients', status )
        IF_NOT_OK_RETURN(status=1)

        ! define:
        call eta_f_file%Def_Sample2D( trim(vname_ps), trim(units_p), &
                                        lon_f_coor, lat_f_coor, ctime_coor, eta_f_sample_coor, &
                                        ps__varid, status )
        IF_NOT_OK_RETURN(status=1)
        call eta_f_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
        IF_NOT_OK_RETURN(status=1)

        ! loop:
        do ivar = 1, self%nvar

          ! get units:
          call eta_file%Get_Att( eta__varids(ivar), 'units', units, status )
          IF_NOT_OK_RETURN(status=1)
          ! define:
          call eta_f_file%Def_CSample3D( trim(self%var(ivar)%name), trim(units), &
                       cplx, lon_f_coor, lat_f_coor, lev_coor, ctime_coor, eta_f_sample_coor, &
                       eta_f__varids(ivar), status, &
                       fill_value=eta_f__fill_value )
          IF_NOT_OK_RETURN(status=1)
          ! annote:
          call eta_f_file%Put_Att( eta_f__varids(ivar), 'long_name', trim(self%var(ivar)%name)//' spectral transform', status )
          IF_NOT_OK_RETURN(status=1)

          ! define:
          call eta_f_file%Def_IField1D( trim(self%var(ivar)%name)//'_fixed', trim(fixed__units), &
                                        lev_coor, &
                                        fixed__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
          call eta_f_file%Put_Att( fixed__varids(ivar), 'long_name', trim(self%var(ivar)%name)//' fixed layer (1=fixed, 0=dynamic)', status )
          IF_NOT_OK_RETURN(status=1)

        end do ! var

        ! end:
        call eta_f_file%EndDef( status )
        IF_NOT_OK_RETURN(status=1)

        ! write coordinate variables:
        call lon_f_coor%Write( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_f_coor%Write( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lev_coor%Write( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call ctime_coor%Write( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! write extended dimensions:
        call lon_ex_dim%Write( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_ex_dim%Write( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%nvar
          ! write:
          call eta_f_file%Put_IField1D( fixed__varids(ivar), fixed(:,ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do
        
      else
      
        ! next record in already open file:
        eta_f_isample = eta_f_isample + 1
        
      end if  ! new output file
        

      ! ~ convert
      
      ! get sample info:
      call eta_sample_coor%Get_Value( eta_isample, status, t=t )
      IF_NOT_OK_RETURN(status=1)
      ! copy to output coordinate:
      call eta_f_sample_coor%Set_Value( eta_f_isample, status, t=t )
      IF_NOT_OK_RETURN(status=1)
    
      ! loop over hours:
      do itime = 1, ntime
    
        ! info ...
        write (gol,'("  time ",i0)') itime; call goPr
        
        ! dumy value defined at spectral grid:
        ps = p0
        ! write sample:
        call eta_f_file%Put_Sample2D( ps__varid, itime, eta_f_isample, ps, status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%nvar
    
          ! info ...
          write (gol,'("    variable ",i0," ",a)') itime, trim(self%var(ivar)%name); call goPr
        
          ! read sample:
          call eta_file%Get_Sample3D( self%var(ivar)%name, itime, eta_isample, eta, units, status )
          IF_NOT_OK_RETURN(status=1)
          
          ! padding:
          eta_ex = 0.0
          eta_ex(1:nlon,1:nlat,1:nlev) = eta
          
          ! loop over levels:
          do ilev = 1, nlev
          
            ! forward real-to-(half-)complex fft:
            call ufft%Forward( eta_ex(:,:,ilev), eta_f(:,:,ilev), status )
            IF_NOT_OK_RETURN(status=1)
            
            ! eliptic truncation (only needed to have nicely masked output):
            call ufft%Truncate( eta_f(:,:,ilev), status, &
                                  fill_value=eta_f__fill_value )
            IF_NOT_OK_RETURN(status=1)
            
            ! truncate spectral info arrays:
            where( eta_f(:,:,ilev) == eta_f__fill_value )
              kstar = kstar__fill_value
              theta = theta__fill_value
            endwhere
            
          end do  ! levels
          
          ! write sample:
          call eta_f_file%Put_CSample3D( eta_f__varids(ivar), itime, eta_f_isample, &
                                           eta_f, status )
          IF_NOT_OK_RETURN(status=1)

        end do ! var
        
      end do  ! time

      ! write static output?
      if ( eta_f_isample == 1 ) then
      
        ! write wavenumbers (truncated):
        call eta_f_file%Put_Field2D( kstar__varid, kstar, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! write spectral angles (truncated):
        call eta_f_file%Put_Field2D( theta__varid, theta, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! write half-complex weights:
        call eta_f_file%Put_IField2D( hcn__varid, hcn, status )
        IF_NOT_OK_RETURN(status=1)
      
      end if  ! static output

      ! next:
      tday = tday + IncrDate(day=1)
      ! end ?
      if ( tday > self%tday2 ) exit
    
    end do  ! samples
    
    ! ~ done with output
    
    ! write times of samples:
    call eta_f_sample_coor%Write( eta_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! close output:
    call eta_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
        
    ! done:
    call eta_f_sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( eta_f__varids )
    deallocate( fixed__varids )

    ! ~ done with fft
    
    ! clear:
    deallocate( eta_ex )
    deallocate( eta_f )
    deallocate( lons_f, lats_f )
    deallocate( kstar )
    deallocate( theta )
    deallocate( hcn )
    
    ! init transform for extended grid:
    call ufft%Done( status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    call lon_f_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_f_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ done with input
    
    ! clear:
    deallocate( ps )
    deallocate( eta )
    deallocate( eta__varids )
    deallocate( fixed )
      
    ! clear:
    call lon_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call eta_sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close input:
    call eta_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~

    ! ok
    status = 0
  
  end subroutine ModelRuns_Eta_f_Samples


  ! ***
  
  
  !
  ! Compute spectral covariances
  !

  subroutine ModelRuns_Compute_C_f( self, rcF, status )
  
    use GO             , only : TDate, IncrDate, operator(+), operator(>), wrtgol
    use GO             , only : TrcFile, ReadRc
    use GO             , only : goReplace
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : IntegerCoordinate, LabelCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Compute_C_f'
    
    ! --- local ----------------------------------

    type(TDate)                      ::  tday, t
    character(len=1024)              ::  filename
    logical                          ::  setup

    character(len=1024)              ::  eta_f_filename
    type(T_NMC_Output)               ::  eta_f_file
    type(Dimension)                  ::  cplx
    type(IntegerCoordinate)          ::  lon_f_coor, lat_f_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(TimeCoordinate)             ::  ctime_coor, sample_coor
    integer                          ::  ncplx
    integer                          ::  nlon_f, nlat_f
    integer                          ::  nlev
    integer                          ::  ntime
    integer                          ::  nsample
    integer                          ::  isample
    real                             ::  p0
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    character(len=32)                ::  units
    real, allocatable                ::  kstar(:,:)
    character(len=32)                ::  kstar__units
    character(len=256)               ::  kstar__formula
    real                             ::  kstar__fill_value
    integer                          ::  kstar__varid
    real, allocatable                ::  theta(:,:)
    character(len=32)                ::  theta__units
    real                             ::  theta__fill_value
    integer                          ::  theta__varid
    integer, allocatable             ::  hcn(:,:)
    character(len=32)                ::  hcn__units
    integer                          ::  hcn__varid
    complex, allocatable             ::  eta_f_sample(:,:,:,:)  ! (lon_f,lat_f,lev,tracer)
    integer                          ::  ilev
    integer                          ::  ivar
    integer                          ::  itime
    real                             ::  eta_f_sample__fill_value
    integer, allocatable             ::  fixed(:,:)  ! (nlev,ntime,nvar)
    character(len=32)                ::  fixed__units
    integer                          ::  fixed__varid
    
    type(T_NMC_Output)               ::  C_f_file
    character(len=1024)              ::  msg
    type(LabelCoordinate)            ::  tracer_coor
    real, allocatable                ::  C_f(:,:,:,:,:,:,:)  ! (nlon_f,nlat_f,nlev,nlev,ntr,nr,nhour)
    real                             ::  C_f__fill_value
    integer                          ::  C_f__varid
    real, allocatable                ::  ps(:,:)   ! (nlon_f,nlat_f,nhour)
    integer                          ::  ps__varid
    integer                          ::  ilev1, ilev2, ilev2_first
    integer                          ::  ivar1, ivar2
    real, allocatable                ::  sigma(:,:,:,:)   ! (nlon_f,nlat_f,nlev,nvar)
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Compute level/comp covariances **")'); call goPr
    write (gol,'("")'); call goPr

    ! start time:
    tday = self%tday1
    ! set flag:
    setup = .false.
    ! init current filenames:
    eta_f_filename = ''
    ! init counter:
    nsample = 0
    ! time loop:
    do
      ! info ...
      call wrtgol( 'day ', tday ); call goPr
    
      ! ~ input samples

      ! input file:
      filename = trim(self%eta_f_samples_filename)
      ! replace keys:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
      IF_NOT_OK_RETURN(status=1)
      
      ! new ?
      if ( trim(filename) /= trim(eta_f_filename) ) then
      
        ! open?
        if ( len_trim(eta_f_filename) > 0 ) then
          ! close:
          call eta_f_file%Close( status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! copy:
        eta_f_filename = trim(filename)
   
        ! open file with eta spectral samples:
        call eta_f_file%Open( trim(filename), status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init counter:
        isample = 1
        
      else
      
        ! next record in already opened file:
        isample = isample + 1
        
      end if
      
      ! ~ setup

      ! initial setup?
      if ( .not. setup ) then
    
        ! init real/imag dimension:
        call cplx%Init( 'cmplx', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call cplx%Read( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call cplx%Get_Dim( status, n=ncplx )
        IF_NOT_OK_RETURN(status=1)
    
        ! create longitude coordinate:
        call lon_f_coor%Init( 'lon_f', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lon_f_coor%Read( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lon_f_coor%Get_Dim( status, n=nlon_f )
        IF_NOT_OK_RETURN(status=1)
    
        ! create latitude coordinate:
        call lat_f_coor%Init( 'lat_f', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lat_f_coor%Read( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lat_f_coor%Get_Dim( status, n=nlat_f )
        IF_NOT_OK_RETURN(status=1)
        
        ! storage for wave numbers:
        allocate( kstar(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! read:
        call eta_f_file%Get_Field2D( 'kstar', kstar, kstar__units, status )
        IF_NOT_OK_RETURN(status=1)
        ! obtain id:
        call eta_f_file%Get_VarID( 'kstar', kstar__varid, status )
        IF_NOT_OK_RETURN(status=1)
        ! info:
        call eta_f_file%Get_Att( kstar__varid, 'formula', kstar__formula, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! storage for spectral angles:
        allocate( theta(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! read:
        call eta_f_file%Get_Field2D( 'theta', theta, theta__units, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! storage for half-complex weights:
        allocate( hcn(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! read:
        call eta_f_file%Get_IField2D( 'hcn', hcn, hcn__units, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! create level coordinate:
        call lev_coor%Init( 'lev', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lev_coor%Read( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lev_coor%Get_Dim( status, n=nlev )
        IF_NOT_OK_RETURN(status=1)
    
        ! create (climatological) time coordinate:
        call ctime_coor%Init( 'time', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call ctime_coor%Read( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call ctime_coor%Get_Dim( status, n=ntime )
        IF_NOT_OK_RETURN(status=1)
    
        ! create sample coordinate:
        call sample_coor%Init( 'sample', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call sample_coor%Read( eta_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! get name of pressure variable:
        call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p, p0=p0 )
        IF_NOT_OK_RETURN(status=1)
    
        ! storage:
        allocate( eta_f_sample(nlon_f,nlat_f,nlev,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
    
        ! storage:
        allocate( fixed(nlev,self%nvar), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over variables:
        do ivar = 1, self%nvar
          ! read:
          call eta_f_file%Get_IField1D( trim(self%var(ivar)%name)//'_fixed', fixed(:,ivar), fixed__units, status )
          IF_NOT_OK_RETURN(status=1)
        end do
        
        ! ~ sample averages
        
        ! storage:
        allocate( ps(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( C_f(nlon_f,nlat_f,nlev,nlev,self%nvar,self%nvar,self%nhour), stat=status )
        IF_NOT_OK_RETURN(status=1)
    
        ! init sums:
        ps = 0.0
        C_f = 0.0
        
        ! ~ output file
        !   (init now to obtain fill_values)
        
        ! create new result file:
        call C_f_file%Create( trim(self%C_f_filename), status )
        IF_NOT_OK_RETURN(status=1)
    
        ! description:
        write (msg,'("spectral level/tracer covariances from `",a,"`")') trim(self%eta_f_samples_filename)
        ! add:
        call C_f_file%Extend_History( trim(msg), status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define coordinates in output file:
        call lon_f_coor%Def( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_f_coor%Def( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define, prepare for no-data values where truncated:
        call C_f_file%Def_Field2D( 'kstar', trim(kstar__units), &
                                      lon_f_coor, lat_f_coor, kstar__varid, status, &
                                      fill_value=kstar__fill_value )
        IF_NOT_OK_RETURN(status=1)
        call C_f_file%Put_Att( kstar__varid, 'long_name', 'wave number', status )
        IF_NOT_OK_RETURN(status=1)
        call C_f_file%Put_Att( kstar__varid, 'formula', trim(kstar__formula), status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define, prepare for no-data values where truncated:
        call C_f_file%Def_Field2D( 'theta', trim(theta__units), &
                                      lon_f_coor, lat_f_coor, theta__varid, status, &
                                      fill_value=theta__fill_value )
        IF_NOT_OK_RETURN(status=1)
        call C_f_file%Put_Att( theta__varid, 'long_name', 'spectral angle', status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define:
        call C_f_file%Def_IField2D( 'hcn', trim(hcn__units), &
                                      lon_f_coor, lat_f_coor, hcn__varid, status )
        IF_NOT_OK_RETURN(status=1)
        call C_f_file%Put_Att( hcn__varid, 'long_name', 'number of represented complex coefficients', status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define level coordinate in output file:
        call lev_coor%Def( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! init tracer coordinate:
        call tracer_coor%Init( 'tracer', status )
        IF_NOT_OK_RETURN(status=1)
        ! set size:
        call tracer_coor%Set_Dim( status, n=self%nvar )
        IF_NOT_OK_RETURN(status=1)
        ! set attributes:
        call tracer_coor%Set_Attrs( status, long_name='tracer' )
        IF_NOT_OK_RETURN(status=1)
        ! loop:
        do ivar = 1, self%nvar
          ! fill tracer name:
          call tracer_coor%Set_Value( ivar, status, value=trim(self%var(ivar)%name) )
          IF_NOT_OK_RETURN(status=1)
        end do  ! var
        ! define tracer coordinate in output file:
        call tracer_coor%Def( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define hour coordinate in output file:
        call ctime_coor%Def( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define:
        call C_f_file%Def_Field2D_Series( trim(vname_ps), trim(units_p), &
                                            lon_f_coor, lat_f_coor, ctime_coor, &
                                            ps__varid, status )
        IF_NOT_OK_RETURN(status=1)
        call C_f_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define:
        call C_f_file%Def_IField2D( 'fixed', trim(fixed__units), &
                                            lev_coor, tracer_coor, &
                                            fixed__varid, status )
        IF_NOT_OK_RETURN(status=1)
        call C_f_file%Put_Att( fixed__varid, 'long_name', 'fixed layer (1=fixed, 0=dynamic)', status )
        IF_NOT_OK_RETURN(status=1)
    
        ! define, obtain fill value for no-data:
        call C_f_file%Def_Covar( 'C_f', '1', &
                                  lon_f_coor, lat_f_coor, lev_coor, tracer_coor, ctime_coor, &
                                  C_f__varid, status, &
                                  fill_value=C_f__fill_value )
        IF_NOT_OK_RETURN(status=1)
        ! annote:
        call C_f_file%Put_Att( C_f__varid, 'long_name', 'covariance of spectral values', status )
        IF_NOT_OK_RETURN(status=1)
    
        ! end:
        call C_f_file%EndDef( status )
        IF_NOT_OK_RETURN(status=1)
    
        ! write coordinate variables:
        call lon_f_coor%Write( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_f_coor%Write( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lev_coor%Write( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call tracer_coor%Write( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
        call ctime_coor%Write( C_f_file, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! write record
        call C_f_file%Put_IField2D( fixed__varid, fixed, status )
        IF_NOT_OK_RETURN(status=1)
    
        ! ~

        ! reset flag:
        setup = .true.
        
      end if  ! initial setup
    
    
      ! loop over hours:
      do itime = 1, ntime
    
        ! info ...
        write (gol,'("  time ",i0)') itime; call goPr
    
        ! info ...
        write (gol,'("    read variables ...")'); call goPr
        
        ! loop over variables:
        do ivar = 1, self%nvar
          ! read sample:
          call eta_f_file%Get_CSample3D( self%var(ivar)%name, itime, isample, &
                                            eta_f_sample(:,:,:,ivar), units, status, &
                                          fill_value=eta_f_sample__fill_value )
          IF_NOT_OK_RETURN(status=1)
        end do
    
        ! info ...
        write (gol,'("    add contribution to covariance ...")'); call goPr
        
        ! loops over upper triangle of tracers:
        do ivar1 = 1, self%nvar
          do ivar2 = ivar1, self%nvar
          
            ! loops over upper triangle of levels:
            do ilev1 = 1, nlev
              if ( ivar2 == ivar1 ) then
                ilev2_first = ilev1
              else
                ilev2_first = 1
              end if
              do ilev2 = ilev2_first, nlev
              
                ! only valid data:
                where ( eta_f_sample(:,:,ilev1,ivar1) == eta_f_sample__fill_value )
                
                  ! no data:
                  C_f(:,:,ilev1,ilev2,ivar1,ivar2,itime) = C_f__fill_value
                  
                elsewhere
                
                  ! update sum:
                  C_f(:,:,ilev1,ilev2,ivar1,ivar2,itime) = C_f(:,:,ilev1,ilev2,ivar1,ivar2,itime) &
                       +  real(eta_f_sample(:,:,ilev1,ivar1)) *  real(eta_f_sample(:,:,ilev2,ivar2)) &
                       + aimag(eta_f_sample(:,:,ilev1,ivar1)) * aimag(eta_f_sample(:,:,ilev2,ivar2))

                endwhere
                
                ! copy to lower triangle:
                C_f(:,:,ilev2,ilev1,ivar2,ivar1,itime) = C_f(:,:,ilev1,ilev2,ivar1,ivar2,itime)
                
              end do ! ilev2
            end do ! ilev1
            
          end do ! ivar2
        end do ! ivar1
        
        !! testing ...
        !write (gol,'("WARNING - only first hour")'); call goPr
        !exit
          
      end do  ! hour
        
      ! increase counter:
      nsample = nsample + 1
      
      !! testing ...
      !write (gol,'("WARNING - only first sample")'); call goPr
      !exit

      ! next:
      tday = tday + IncrDate(day=1)
      ! end ?
      if ( tday > self%tday2 ) exit
    
    end do  ! samples
    
    ! info ...
    write (gol,'("sample averages ...")'); call goPr
    
    ! check ...
    if ( nsample < 2 ) then
      write (gol,'("at least 2 samples required to compute correlations, found ",i0)') nsample; call goErr
      TRACEBACK; status=1; return
    end if

    ! sample covariance requires division by N-1 :
    where ( C_f /= C_f__fill_value )
      C_f = C_f  / (nsample-1.0)
    endwhere
    
    ! info ...
    write (gol,'("write averages ...")'); call goPr
    
    ! loop over hours:
    do itime = 1, ntime

      ! fill with reference value:
      ps = p0
      ! write average:
      call C_f_file%Put_Field2D_Series( ps__varid, itime, ps, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! write average:
      call C_f_file%Put_Covar( C_f__varid, itime, &
                                  C_f(:,:,:,:,:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
          
    end do ! hours

    ! write wave numbers:
    call C_f_file%Put_Field2D( kstar__varid, kstar, status )
    IF_NOT_OK_RETURN(status=1)

    ! write spectral angles:
    call C_f_file%Put_Field2D( theta__varid, theta, status )
    IF_NOT_OK_RETURN(status=1)

    ! write half-complex weights:
    call C_f_file%Put_IField2D( hcn__varid, hcn, status )
    IF_NOT_OK_RETURN(status=1)

    ! done with coordinates:
    call tracer_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close output:
    call C_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ done with output
    
    ! clear:
    deallocate( C_f )
    deallocate( ps )
    
    ! ~ done with input
    
    ! clear:
    deallocate( kstar )
    deallocate( theta )
    deallocate( hcn )
    deallocate( eta_f_sample )
    deallocate( fixed )
      
    ! clear:
    call lon_f_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_f_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close input:
    call eta_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ 

    ! ok
    status = 0
  
  end subroutine ModelRuns_Compute_C_f


  ! ***
  
  
  !
  ! Compute angular averages.
  ! 

  subroutine ModelRuns_Compute_D_f( self, rcF, status )
  
    use GO             , only : TrcFile, ReadRc
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : IntegerCoordinate, LabelCoordinate
    use C3PO           , only : RealCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Compute_D_f'
  
    ! gonio:
    real, parameter     ::  pi = 4.0 * atan(1.0)
    
    ! --- local ----------------------------------

    type(T_NMC_Output)               ::  C_f_file
    type(IntegerCoordinate)          ::  lon_f_coor, lat_f_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(LabelCoordinate)            ::  tracer_coor
    type(TimeCoordinate)             ::  ctime_coor
    integer                          ::  nlon_f, nlat_f
    integer                          ::  nlev
    integer                          ::  ntracer
    integer                          ::  ntime
    real                             ::  p0
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    character(len=32)                ::  units
    real, allocatable                ::  kstar(:,:)
    real                             ::  kstar__fill_value
    character(len=32)                ::  kstar__units
    character(len=256)               ::  kstar__formula
    integer                          ::  kstar__varid
    real, allocatable                ::  theta(:,:)
    real                             ::  theta__fill_value
    character(len=32)                ::  theta__units
    integer                          ::  theta__varid
    integer, allocatable             ::  hcn(:,:)
    character(len=32)                ::  hcn__units
    integer, allocatable             ::  fixed(:,:)  ! (nlev,nvar)
    character(len=32)                ::  fixed__units
    integer                          ::  fixed__varid
    real, allocatable                ::  C_f(:,:,:,:,:,:)  ! (nlon_f,nlat_f,nlev,nlev,ntr,nr)
    integer                          ::  itime

    real                             ::  ks(2)
    real                             ::  kstar_max
    real                             ::  akstar_step
    integer                          ::  nakstar
    real, allocatable                ::  akstar_bnds(:,:)  ! (nv,nakstar)
    integer                          ::  ik
    integer, allocatable             ::  iakstar(:,:)      ! (nlon_f,nlat_f)
    integer                          ::  iakstar__fill_value
    integer                          ::  iakstar__varid
    integer                          ::  i, j
    type(RealCoordinate)             ::  akstar_coor
    
    integer                          ::  ntheta90, ntheta
    real                             ::  DeltaTheta0
    real                             ::  theta0
    integer                          ::  ith
    integer, allocatable             ::  itheta(:,:)      ! (nlon_f,nlat_f)
    integer, allocatable             ::  ntheta_r(:,:)    ! (nakstar,ntheta)
    integer, allocatable             ::  ntheta_v(:,:)    ! (nakstar,ntheta)
    real, allocatable                ::  DeltaTheta(:,:)  ! (nakstar,ntheta)
    integer                          ::  it, ds1, ds2
    real                             ::  dth_r, dth_v
    real, allocatable                ::  dtheta_sum(:)  ! (nakstar)
    real, allocatable                ::  dtheta(:,:)    ! (nlon_f,nlat_f)
    real                             ::  dtheta__fill_value
    integer                          ::  dtheta__varid
    
    type(T_NMC_Output)               ::  D_f_file
    character(len=1024)              ::  msg
    real, allocatable                ::  D_f(:,:,:,:,:)  ! (nakstar,nlev,nlev,ntr,ntr)
    real                             ::  D_f__fill_value
    integer                          ::  D_f__varid
    real, allocatable                ::  ps(:,:)   ! (nakstar,ntime)
    integer                          ::  ps__varid
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Compute level/comp covariances angular average **")'); call goPr
    write (gol,'("")'); call goPr
    
    ! ~ input samples
          
    ! open file with C_f:
    call C_f_file%Open( trim(self%C_f_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create longitude coordinate:
    call lon_f_coor%Init( 'lon_f', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lon_f_coor%Read( C_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lon_f_coor%Get_Dim( status, n=nlon_f )
    IF_NOT_OK_RETURN(status=1)

    ! create latitude coordinate:
    call lat_f_coor%Init( 'lat_f', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lat_f_coor%Read( C_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lat_f_coor%Get_Dim( status, n=nlat_f )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for wave numbers:
    allocate( kstar(nlon_f,nlat_f), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read, obtain value used for no-data:
    call C_f_file%Get_Field2D( 'kstar', kstar, kstar__units, status, &
                                 fill_value=kstar__fill_value )
    IF_NOT_OK_RETURN(status=1)
    ! obtain id:
    call C_f_file%Get_VarID( 'kstar', kstar__varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! info:
    call C_f_file%Get_Att( kstar__varid, 'formula', kstar__formula, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for spectral angles:
    allocate( theta(nlon_f,nlat_f), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read, obtain value used for no-data:
    call C_f_file%Get_Field2D( 'theta', theta, theta__units, status, &
                                 fill_value=theta__fill_value )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for half-complex weights:
    allocate( hcn(nlon_f,nlat_f), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call C_f_file%Get_IField2D( 'hcn', hcn, hcn__units, status )
    IF_NOT_OK_RETURN(status=1)

    ! create level coordinate:
    call lev_coor%Init( 'lev', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lev_coor%Read( C_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lev_coor%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( C_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=ntracer )
    IF_NOT_OK_RETURN(status=1)

    ! create (climatological) time coordinate:
    call ctime_coor%Init( 'time', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call ctime_coor%Read( C_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call ctime_coor%Get_Dim( status, n=ntime )
    IF_NOT_OK_RETURN(status=1)
    
    ! get name of pressure variable:
    call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p, p0=p0 )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage:
    allocate( C_f(nlon_f,nlat_f,nlev,nlev,ntracer,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( fixed(nlev,self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call C_f_file%Get_IField2D( 'fixed', fixed, fixed__units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ average kstar
    
    ! maximum kstar along dimensions:
    ks(1) = maxval( kstar(:,1), mask=(kstar(:,1) /= kstar__fill_value) )
    ks(2) = maxval( kstar(1,:), mask=(kstar(1,:) /= kstar__fill_value) )
    ! elliptic values expected, check this:
    if ( ks(1) /= ks(2) ) then
      write (gol,'("expected similar k* values along dimensions, found: ",2f16.8)') ks; call goErr
      TRACEBACK; status=1; return
    end if
    ! for safety, copy minimum value to ensure circles,
    ! but should be the same as tested above:
    kstar_max = maxval( ks )
    
    ! read stepping:
    call ReadRc( rcF, 'nmc.D_f.akstar_step', akstar_step, status )
    IF_NOT_OK_RETURN(status=1)

    ! number of steps:
    nakstar = int(ceiling( kstar_max / akstar_step ))

    ! storage for boundaries:
    allocate( akstar_bnds(2,nakstar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do ik = 1, nakstar
      ! lower and upper bound:
      akstar_bnds(1,ik) = (ik-1) * akstar_step
      akstar_bnds(2,ik) =  ik    * akstar_step
    end do
    
    ! storage for (lon_f,lat_f) to akstar index:
    allocate( iakstar(nlon_f,nlat_f), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! undefined, later on replaced by fill value:
    iakstar = -999
    ! loop over spectral coeff:
    do j = 1, nlat_f
      ! loop over spectral coeff:
      do i = 1, nlon_f
        ! skip no-data values:
        if ( kstar(i,j) == kstar__fill_value ) cycle
        ! index:
        ik = int(ceiling( kstar(i,j) / akstar_step ))
        ! first is zero, add to first bin:
        if ( kstar(i,j) == 0.0 ) ik = 1
        ! check ..
        if ( (ik < 1) .or. (ik > nakstar) ) then
          write (gol,'("unexpected index ",i0," in ak* array with step ",f8.3)') ik,akstar_step; call goErr
          write (gol,'("  spectral index : (",i0,",",i0,")")') i,j; call goErr
          write (gol,'("  k* value       : ",f8.3)') kstar(i,j); call goErr
          TRACEBACK; status=1; return
        end if
        ! store:
        iakstar(i,j) = ik
      end do  ! i
    end do  ! j
    
    ! define ak* coordinate:
    call akstar_coor%Init( 'akstar', status )
    IF_NOT_OK_RETURN(status=1)
    call akstar_coor%Set_Dim( status, n=nakstar )
    IF_NOT_OK_RETURN(status=1)
    call akstar_coor%Set_Attrs( status, units=kstar__units, &
                                  long_name='average wave number', &
                                  formula=trim(kstar__formula) )
    IF_NOT_OK_RETURN(status=1)
    call akstar_coor%Set_Values( status, value_bnds=akstar_bnds )
    IF_NOT_OK_RETURN(status=1)

    
    ! ~ angular weights
    
    ! number of theta segments per 90 deg:
    call ReadRc( rcF, 'nmc.D_f.ntheta', ntheta90, status )
    IF_NOT_OK_RETURN(status=1)
    ! bins are centered around -pi/2,..,0,..,pi/2,
    ! thus segments [-pi/2-DeltaTheta/2,-pi/2+DeltaTheta/2] etc
    ! number of segments:
    ntheta = 2*ntheta90 + 1
    ! initial spacing:
    DeltaTheta0 = (0.5*pi) / real(ntheta90)   ! rad
    ! offset:
    theta0 = (-pi/2.0) - (0.5*DeltaTheta0)
    
    ! storage segment index:
    allocate( itheta(nlon_f,nlat_f), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for weights:
    allocate( ntheta_r(nakstar,ntheta), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( ntheta_v(nakstar,ntheta), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! init counters:
    ntheta_r = 0
    ntheta_v = 0
    ! loop over spectral coeff:
    do j = 1, nlat_f
      do i = 1, nlon_f
        ! skip undefined:
        if ( kstar(i,j) == kstar__fill_value ) cycle
        ! average k* bin:
        ik = iakstar(i,j)
        ! theta bin:
        ith = (theta(i,j) - theta0)/DeltaTheta0 + 1
        ! store:
        itheta(i,j) = ith
        ! increase counters:
        ntheta_r(ik,ith) = ntheta_r(ik,ith) + 1
        if ( hcn(i,j) >= 2 ) ntheta_v(ik,ith) = ntheta_v(ik,ith) + 1
      end do ! i
    end do ! j
    
    ! storage for segment size, might be larger then initial size
    ! to account for segments without data:
    allocate( DeltaTheta(nakstar,ntheta), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! initial values:
    DeltaTheta = DeltaTheta0
    ! loop over segments:
    do ik = 1, nakstar
      do ith = 1, ntheta
        ! no values ?
        if ( ntheta_r(ik,ith) == 0 ) then
          ! check ...
          if ( (ith == 1) .or. (ith == ntheta) ) then
            write (gol,'("found zero ntheta_r at k* band ",i0," segment ",i0)') ik, ith; call goErr
            TRACEBACK; status=1; return
          end if
          ! find distances to  neigbours at each side:
          ds1 = -999
          do it = ith-1, 1, -1
            if ( ntheta_r(ik,it) > 0 ) then
              ds1 = ith - it
              exit
            end if
          end do
          ds2 = -999
          do it = ith+1, ntheta
            if ( ntheta_r(ik,it) > 0 ) then
              ds2 = it - ith
              exit
            end if
          end do
          ! check ...
          if ( (ds1 < 0) .or. (ds2 < 0) ) then
            write (gol,'("could not find neighbour with non-zero number of elements")'); call goPr
            write (gol,'("  k* band ",i0," segment ",i0)') ik, ith; call goErr
            write (gol,*) '  ntheta_r = ', ntheta_r(ik,:); call goErr
            TRACEBACK; status=1; return
          end if
          ! distribute over neighbours, 
          ! smaller fraction to neighbours further away;
          ! segments around -pi/2 and pi/2 receive double:
          if ( ith-ds1 == 1 ) then
            DeltaTheta(ik,ith-ds1) = DeltaTheta(ik,ith-ds1) + 2 * real(ds2)/real(ds1+ds2) * DeltaTheta(ik,ith)
          else
            DeltaTheta(ik,ith-ds1) = DeltaTheta(ik,ith-ds1) +     real(ds2)/real(ds1+ds2) * DeltaTheta(ik,ith)
          end if
          if ( ith+ds2 == ntheta ) then
            DeltaTheta(ik,ith+ds2) = DeltaTheta(ik,ith+ds2) + 2 * real(ds1)/real(ds1+ds2) * DeltaTheta(ik,ith)
          else
            DeltaTheta(ik,ith+ds2) = DeltaTheta(ik,ith+ds2) +     real(ds1)/real(ds1+ds2) * DeltaTheta(ik,ith)
          end if
          ! reset for safety ...
          DeltaTheta(ik,ith) = 0
        end if  ! reset
      end do
    end do
    
    ! storage for end result:
    allocate( dtheta(nlon_f,nlat_f), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! for check:
    allocate( dtheta_sum(nakstar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! init sum:
    dtheta_sum = 0.0
    
    ! init result, negative values will be replaced by fill_value:
    dtheta = -999.9
    ! counted, fill dth
    do j = 1, nlat_f
      do i = 1, nlon_f
        ! skip undefined:
        if ( kstar(i,j) == kstar__fill_value ) cycle
        ! bin indices:
        ik  = iakstar(i,j)
        ith = itheta(i,j)
        ! check ...
        if ( DeltaTheta(ik,ith) <= 0.0 ) then
          write (gol,'("found zero DeltaTheta(",i0,",",i0,") = ",f0.4)') DeltaTheta(ik,ith); call goErr
          TRACEBACK; status=1; return
        end if
        ! contributions, different for -pi and pi segments:
        if ( (ith == 1) .or. (ith == ntheta) ) then
          dth_r = DeltaTheta(ik,ith) / ( ntheta_r(ik,ith) + ntheta_v(ik,ith) )
          dth_v = dth_r
        else
          dth_r = DeltaTheta(ik,ith) / ntheta_r(ik,ith)
          dth_v = DeltaTheta(ik,ith) / ntheta_v(ik,ith)
        end if
        ! assign portion of circle to coeff:
        if ( hcn(i,j) == 1 ) then
          dtheta(i,j) = dth_r
        else if ( hcn(i,j) == 2 ) then
          dtheta(i,j) = dth_r + dth_v
        else
          write (gol,'("unsupported hcn value ",i0)') hcn(i,j); call goErr
          TRACEBACK; status=1; return
        end if
        ! add contribution:
        dtheta_sum(ik) = dtheta_sum(ik) + dtheta(i,j)
      end do ! i
    end do ! j
    
    ! check ...
    do ik = 1, nakstar
      if ( abs( dtheta_sum(ik) - 2*pi ) > 0.01 ) then
        write (gol,'("dtheta sum for k* band ",i0," equal to ",f0.4," while 2pi is ",f0.4)') &
                   dtheta_sum(ik), 2*pi; call goErr
      end if
    end do
    ! clear:
    deallocate( dtheta_sum, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ output

    ! create new result file:
    call D_f_file%Create( trim(self%D_f_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("angular averages of spectral level/tracer covariances from `",a,"`")') trim(self%C_f_filename)
    ! add:
    call D_f_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define coordinates in output file:
    call lon_f_coor%Def( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_f_coor%Def( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define:
    call D_f_file%Def_Field2D( 'kstar', trim(kstar__units), &
                                    lon_f_coor, lat_f_coor, kstar__varid, status, &
                                    fill_value=kstar__fill_value )
    IF_NOT_OK_RETURN(status=1)
    call D_f_file%Put_Att( kstar__varid, 'long_name', 'wave number', status )
    IF_NOT_OK_RETURN(status=1)
    call D_f_file%Put_Att( kstar__varid, 'formula', trim(kstar__formula), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define:
    call D_f_file%Def_Field2D( 'theta', trim(theta__units), &
                                    lon_f_coor, lat_f_coor, theta__varid, status, &
                                    fill_value=theta__fill_value )
    IF_NOT_OK_RETURN(status=1)
    call D_f_file%Put_Att( theta__varid, 'long_name', 'wave number', status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define:
    call D_f_file%Def_Field2D( 'dtheta', trim(theta__units), &
                                    lon_f_coor, lat_f_coor, dtheta__varid, status, &
                                    fill_value=dtheta__fill_value )
    IF_NOT_OK_RETURN(status=1)
    call D_f_file%Put_Att( dtheta__varid, 'long_name', 'weight in angular average', status )
    IF_NOT_OK_RETURN(status=1)
    ! reset undefined values:
    where ( dtheta < 0.0 )
      dtheta = dtheta__fill_value
    end where
    
    ! define ak* coordinate in output file:
    call akstar_coor%Def( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define mapping field
    call D_f_file%Def_IField2D( 'iakstar', '1', &
                                    lon_f_coor, lat_f_coor, iakstar__varid, status, &
                                    fill_value=iakstar__fill_value )
    IF_NOT_OK_RETURN(status=1)
    call D_f_file%Put_Att( iakstar__varid, 'long_name', 'index in '//trim(akstar_coor%name), status )
    IF_NOT_OK_RETURN(status=1)

    ! define level coordinate in output file:
    call lev_coor%Def( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define tracer coordinate in output file:
    call tracer_coor%Def( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call D_f_file%Def_Field2D( trim(vname_ps), trim(units_p), &
                                 akstar_coor, ctime_coor, &
                                 ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call D_f_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call D_f_file%Def_IField2D( 'fixed', trim(fixed__units), &
                                  lev_coor, tracer_coor, &
                                  fixed__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call D_f_file%Put_Att( fixed__varid, 'long_name', 'fixed layer (1=fixed, 0=dynamic)', status )
    IF_NOT_OK_RETURN(status=1)

    ! define, obtain fill value for no-data:
    call D_f_file%Def_ACovar( 'D_f', '1', &
                              akstar_coor, lev_coor, tracer_coor, ctime_coor, &
                              D_f__varid, status, &
                              fill_value=D_f__fill_value )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call D_f_file%Put_Att( D_f__varid, 'long_name', 'angular average of covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)
    
    ! end:
    call D_f_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lon_f_coor%Write( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_f_coor%Write( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call akstar_coor%Write( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Write( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Write( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Write( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write wave numbers:
    call D_f_file%Put_Field2D( kstar__varid, kstar, status )
    IF_NOT_OK_RETURN(status=1)

    ! write angles:
    call D_f_file%Put_Field2D( theta__varid, theta, status )
    IF_NOT_OK_RETURN(status=1)

    ! write weights in agular average:
    call D_f_file%Put_Field2D( dtheta__varid, dtheta, status )
    IF_NOT_OK_RETURN(status=1)

    ! reset undefined indices to fill value:
    where ( iakstar < 0 )
      iakstar = iakstar__fill_value
    endwhere
    ! write:
    call D_f_file%Put_IField2D( iakstar__varid, iakstar, status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( ps(nakstar,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill with reference pressure:
    ps = p0
    ! write:
    call D_f_file%Put_Field2D( ps__varid, ps, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write flags:
    call D_f_file%Put_IField2D( fixed__varid, fixed, status )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( D_f(nakstar,nlev,nlev,ntracer,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ time loop

    ! info ...
    write (gol,'("dimensions:")'); call goPr
    write (gol,'("  number of time values : ",i0)') ntime; call goPr
    write (gol,'("  number of tracers     : ",i0)') ntracer; call goPr
    
    ! loop over hours:
    do itime = 1, ntime
    
      ! info ...
      write (gol,'("  itime ",i0)') itime; call goPr

      ! read spectral covariances:
      call C_f_file%Get_Covar( 'C_f', itime, C_f, units, status )
      IF_NOT_OK_RETURN(status=1)

      ! init sum:
      D_f = 0.0
      ! loop over spectral coeff:
      do j = 1, nlat_f
        do i = 1, nlon_f
          ! skip undefined:
          if ( kstar(i,j) == kstar__fill_value ) cycle
          ! target bin:
          ik = iakstar(i,j)
          ! add contribution:
          D_f(ik,:,:,:,:) = D_f(ik,:,:,:,:) + C_f(i,j,:,:,:,:) * dtheta(i,j)
        end do
      end do
      
      ! loop over sums:
      do ik = 1, nakstar
        ! average:
        D_f(ik,:,:,:,:) = D_f(ik,:,:,:,:) / (2*pi)
      end do  ! sums
      
      ! write average:
      call D_f_file%Put_ACovar( D_f__varid, itime, D_f, status )
      IF_NOT_OK_RETURN(status=1)
          
    end do  ! times

    ! ~ done with output
    
    ! clear:
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( D_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
   
    ! close output:
    call D_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ done with avers
    
    ! clear:
    deallocate( akstar_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( iakstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( itheta, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( ntheta_r, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( ntheta_v, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( DeltaTheta, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( dtheta, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ done with input
    
    ! clear:
    deallocate( kstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( theta, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( C_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( fixed, stat=status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    call lon_f_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_f_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call akstar_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close input:
    call C_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ 

    ! ok
    status = 0
  
  end subroutine ModelRuns_Compute_D_f


  ! ***
  
  
  !
  ! Compute B_f = D_f/sqrt(gamma*gamma)
  ! 

  subroutine ModelRuns_Compute_B_f( self, rcF, status )
  
    use GO             , only : TrcFile, ReadRc
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Compute_B_f'
  
    ! gonio:
    real, parameter     ::  pi = 4.0 * atan(1.0)
    
    ! --- local ----------------------------------

    type(T_NMC_Output)               ::  D_f_file
    type(RealCoordinate)             ::  akstar_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(LabelCoordinate)            ::  tracer_coor
    type(TimeCoordinate)             ::  ctime_coor
    integer                          ::  nakstar
    integer                          ::  nlev
    integer                          ::  ntracer
    integer                          ::  ntime
    real                             ::  p0
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    integer, allocatable             ::  fixed(:,:)  ! (nlev,nvar)
    character(len=32)                ::  fixed__units
    integer                          ::  fixed__varid
    real, allocatable                ::  akstar(:)  ! (nakstar)
    real, allocatable                ::  D_f(:,:,:,:,:)  ! (nakstar,nlev,nlev,ntr,ntr)
    character(len=32)                ::  D_f__units
    integer                          ::  itime

    character(len=1024)              ::  msg
    real, allocatable                ::  ps(:,:)   ! (nakstar,ntime)
    integer                          ::  ps__varid

    type(T_NMC_Output)               ::  gamma_file
    real, allocatable                ::  gamma(:,:,:)  ! (nakstar,nlev,ntr)
    integer                          ::  gamma__varid
    
    type(T_NMC_Output)               ::  B_f_file
    real, allocatable                ::  B_f(:,:,:,:,:)  ! (nakstar,nlev,nlev,ntr,ntr)
    integer                          ::  B_f__varid

    integer                          ::  ilev, ilev1, ilev2
    integer                          ::  itracer, itracer1, itracer2
    integer                          ::  ik
    real                             ::  D_f_diagelem_sum
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Compute gamma and B_f **")'); call goPr
    write (gol,'("")'); call goPr
    
    ! ~ input samples
          
    ! open file with D_f:
    call D_f_file%Open( trim(self%D_f_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create average k* coordinate:
    call akstar_coor%Init( 'akstar', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call akstar_coor%Read( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call akstar_coor%Get_Dim( status, n=nakstar )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( akstar(nakstar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! extract:
    call akstar_coor%Get_Values( status, values=akstar )
    IF_NOT_OK_RETURN(status=1)

    ! create level coordinate:
    call lev_coor%Init( 'lev', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lev_coor%Read( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lev_coor%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=ntracer )
    IF_NOT_OK_RETURN(status=1)

    ! create (climatological) time coordinate:
    call ctime_coor%Init( 'time', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call ctime_coor%Read( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call ctime_coor%Get_Dim( status, n=ntime )
    IF_NOT_OK_RETURN(status=1)
    
    ! get name of pressure variable:
    call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p, p0=p0 )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage:
    allocate( D_f(nakstar,nlev,nlev,ntracer,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( fixed(nlev,self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call D_f_file%Get_IField2D( 'fixed', fixed, fixed__units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ output general
    
    ! storage:
    allocate( ps(nakstar,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill with reference pressure:
    ps = p0

    ! ~ output: gamma

    ! create new result file:
    call gamma_file%Create( trim(self%gamma_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("normalization of angular averages of spectral level/tracer covariances from `",a,"`")') trim(self%D_f_filename)
    ! add:
    call gamma_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define ak* coordinate in output file:
    call akstar_coor%Def( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define level coordinate in output file:
    call lev_coor%Def( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define tracer coordinate in output file:
    call tracer_coor%Def( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call gamma_file%Def_Field2D( trim(vname_ps), trim(units_p), &
                                 akstar_coor, ctime_coor, &
                                 ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call gamma_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call gamma_file%Def_IField2D( 'fixed', trim(fixed__units), &
                                        lev_coor, tracer_coor, &
                                        fixed__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call gamma_file%Put_Att( fixed__varid, 'long_name', 'fixed layer (1=fixed, 0=dynamic)', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call gamma_file%Def_AVar( 'gamma', '1', &
                              akstar_coor, lev_coor, tracer_coor, ctime_coor, &
                              gamma__varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call gamma_file%Put_Att( gamma__varid, 'long_name', 'scaling for angular average of covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)
    
    ! end:
    call gamma_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call akstar_coor%Write( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Write( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Write( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Write( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write:
    call gamma_file%Put_Field2D( ps__varid, ps, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write flags:
    call gamma_file%Put_IField2D( fixed__varid, fixed, status )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( gamma(nakstar,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ output: B_f

    ! create new result file:
    call B_f_file%Create( trim(self%B_f_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("angular averages of spectral level/tracer covariances from `",a,"`")') trim(self%C_f_filename)
    write (msg,'(a," scalled with `",a,"`")') trim(msg), trim(self%gamma_filename)
    ! add:
    call B_f_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define ak* coordinate in output file:
    call akstar_coor%Def( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define level coordinate in output file:
    call lev_coor%Def( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define tracer coordinate in output file:
    call tracer_coor%Def( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call B_f_file%Def_Field2D( trim(vname_ps), trim(units_p), &
                                 akstar_coor, ctime_coor, &
                                 ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call B_f_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call B_f_file%Def_IField2D( 'fixed', trim(fixed__units), &
                                   lev_coor, tracer_coor, &
                                   fixed__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call B_f_file%Put_Att( fixed__varid, 'long_name', 'fixed layer (1=fixed, 0=dynamic)', status )
    IF_NOT_OK_RETURN(status=1)

    ! define, obtain fill value for no-data:
    call B_f_file%Def_ACovar( 'B_f', '1', &
                              akstar_coor, lev_coor, tracer_coor, ctime_coor, &
                              B_f__varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call B_f_file%Put_Att( B_f__varid, 'long_name', 'scaled angular average of covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)
    
    ! end:
    call B_f_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call akstar_coor%Write( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Write( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Write( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Write( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write:
    call B_f_file%Put_Field2D( ps__varid, ps, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write flags:
    call B_f_file%Put_IField2D( fixed__varid, fixed, status )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( B_f(nakstar,nlev,nlev,ntracer,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ time loop

    ! info ...
    write (gol,'("dimensions:")'); call goPr
    write (gol,'("  number of k* values   : ",i0)') nakstar; call goPr
    write (gol,'("  number of time values : ",i0)') ntime; call goPr
    write (gol,'("  number of tracers     : ",i0)') ntracer; call goPr
    
    ! loop over hours:
    do itime = 1, ntime
    
      ! info ...
      write (gol,'("  itime ",i0)') itime; call goPr

      ! read angular averaged spectral covariances:
      call D_f_file%Get_ACovar( 'D_f', itime, D_f, D_f__units, status )
      IF_NOT_OK_RETURN(status=1)
      
      !
      !                        D_f(k,l;k,l;k*)
      ! gamma(k,l;k*) = ---------------------------------
      !                 sum_k*' 2 pi k*' D_f(k,l;k,l;k*')
      !
      ! loop over diagonal elements:
      do itracer = 1, ntracer
        do ilev = 1, nlev
        
          ! sum:
          D_f_diagelem_sum = sum( 2*pi*akstar * D_f(:,ilev,ilev,itracer,itracer) )
        
          ! non-zero variance?
          if ( D_f_diagelem_sum > 0.0 ) then
            ! fill scale factor:
            gamma(:,ilev,itracer) = D_f(:,ilev,ilev,itracer,itracer) / D_f_diagelem_sum
          else
            ! no scaling:
            gamma(:,ilev,itracer) = 0.0
          end if
          
        end do ! levels
      end do  ! tracers
      
      !
      !                                D_f(k,l;k',l';k*)
      !  B_f(k,l;k',l';k*) = -------------------------------------
      !                      sqrt( gamma(k,l;k*) gamma(k',l';k*) )
      !
      ! loop over covariance elements:
      do itracer1 = 1, ntracer
        do itracer2 = 1, ntracer
          do ilev1 = 1, nlev
            do ilev2 = 1, nlev
            
              ! loop over wave numbers:
              do ik = 1, nakstar
              
                ! trap zero variances:
                if ( (gamma(ik,ilev1,itracer1) > 0.0) .and. &
                     (gamma(ik,ilev2,itracer2) > 0.0)       ) then
                  ! fill scalled version:
                  B_f(ik,ilev1,ilev2,itracer1,itracer2) = &
                    D_f(ik,ilev1,ilev2,itracer1,itracer2) / &
                    sqrt( gamma(ik,ilev1,itracer1) * gamma(ik,ilev2,itracer2) )
                else
                  ! no variance (top level?), so remain zero:
                  B_f(ik,ilev1,ilev2,itracer1,itracer2) = 0.0
                end if  ! zero variance
                
              end do  ! k*
              
            end do  ! ilev2
          end do ! ilev1
        end do ! itracer2
      end do ! itracer1

      ! write scaling:
      call gamma_file%Put_AVar( gamma__varid, itime, gamma, status )
      IF_NOT_OK_RETURN(status=1)
      ! write covar:
      call B_f_file%Put_ACovar( B_f__varid, itime, B_f, status )
      IF_NOT_OK_RETURN(status=1)

    end do  ! times

    ! ~ done with output
    
    ! clear:
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( gamma, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! close output:
    call gamma_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( B_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! close output:
    call B_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ done with input
    
    ! clear:
    deallocate( akstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( fixed, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( D_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    call akstar_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close input:
    call D_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ 

    ! ok
    status = 0
  
  end subroutine ModelRuns_Compute_B_f


  ! ***
  
  
  !
  ! Compute eigenvalue decompo:
  !   B_f = X Lambda X'
  ! 
  ! Intel MKL on Lapack eigenvalue solvers:
  !
  !   https://software.intel.com/en-us/node/469188#DE6F1082-B885-4C07-BAC2-4E844ACC14E3
  !

  subroutine ModelRuns_Compute_XLX( self, rcF, status )
  
    use Lapack95, only : SyEvR
  
    use GO             , only : TrcFile, ReadRc
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    use C3PO           , only : IntegerCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Compute_XLX'
  
    ! --- local ----------------------------------

    type(T_NMC_Output)               ::  B_f_file
    type(RealCoordinate)             ::  akstar_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(LabelCoordinate)            ::  tracer_coor
    type(TimeCoordinate)             ::  ctime_coor
    integer                          ::  nakstar
    integer                          ::  nlev
    integer                          ::  ntracer
    integer                          ::  ntime
    real                             ::  p0
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    integer, allocatable             ::  fixed(:,:)  ! (nlev,nvar)
    character(len=32)                ::  fixed__units
    integer                          ::  fixed__varid
    real, allocatable                ::  B_f(:,:,:,:,:)  ! (nakstar,nlev,nlev,ntr,ntr)
    character(len=32)                ::  B_f__units
    integer                          ::  itime
    
    character(len=1024)              ::  msg
    real, allocatable                ::  ps(:,:)   ! (nakstar,ntime)
    integer                          ::  ps__varid

    type(T_NMC_Output)               ::  XLX_file
    type(IntegerCoordinate)          ::  eval_coor
    integer, allocatable             ::  evals(:)     ! (nv)
    integer, allocatable             ::  nev(:)       ! (nakstar)
    integer                          ::  nev__varid
    real, allocatable                ::  Lambda(:,:)  ! (nakstar,nv)
    real                             ::  Lambda__fill_value
    integer                          ::  Lambda__varid
    real, allocatable                ::  X(:,:,:,:)   ! (nakstar,nlev,ntr,nv)
    real                             ::  X__fill_value
    integer                          ::  X__varid
    
    integer                          ::  ik
    integer                          ::  nv, iv, iv1, iv2
    integer                          ::  ilev, ilev1, ilev2
    integer                          ::  itracer, itracer1, itracer2
    real, allocatable                ::  A(:,:)  ! (nv,nv)
    real, allocatable                ::  Z(:,:)  ! (nv,nv)
    real, allocatable                ::  w(:)  ! (nv)
    integer                          ::  m, m1, m2
    integer                          ::  nzero
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Compute eigenvalue decompositions **")'); call goPr
    write (gol,'("")'); call goPr
    
    ! ~ input samples
          
    ! open file with B_f:
    call B_f_file%Open( trim(self%B_f_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create average k* coordinate:
    call akstar_coor%Init( 'akstar', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call akstar_coor%Read( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call akstar_coor%Get_Dim( status, n=nakstar )
    IF_NOT_OK_RETURN(status=1)

    ! create level coordinate:
    call lev_coor%Init( 'lev', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lev_coor%Read( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lev_coor%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=ntracer )
    IF_NOT_OK_RETURN(status=1)

    ! create (climatological) time coordinate:
    call ctime_coor%Init( 'time', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call ctime_coor%Read( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call ctime_coor%Get_Dim( status, n=ntime )
    IF_NOT_OK_RETURN(status=1)
    
    ! get name of pressure variable:
    call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p, p0=p0 )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage:
    allocate( B_f(nakstar,nlev,nlev,ntracer,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( fixed(nlev,self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call B_f_file%Get_IField2D( 'fixed', fixed, fixed__units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ eigenvalue dimension

    ! maximum number of eigenvalues:
    nv = nlev * ntracer

    ! create longitude coordinate:
    call eval_coor%Init( 'eigenvalue', status )
    IF_NOT_OK_RETURN(status=1)
    ! set size:
    call eval_coor%Set_Dim( status, n=nv )
    IF_NOT_OK_RETURN(status=1)
    ! set attributes:
    call eval_coor%Set_Attrs( status, units='1', long_name='eigen value' )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( evals(nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do iv = 1, nv
      evals(iv) = iv
    end do
    ! fill values:
    call eval_coor%Set_Values( status, values=evals )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( evals, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ output
    
    ! storage:
    allocate( ps(nakstar,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill with reference pressure:
    ps = p0

    ! create new result file:
    call XLX_file%Create( trim(self%XLX_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("eigenvalue decomposition of normalized angular averages of spectral level/tracer covariances from `",a,"`")') trim(self%B_f_filename)
    ! add:
    call XLX_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define ak* coordinate in output file:
    call akstar_coor%Def( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define level coordinate in output file:
    call lev_coor%Def( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define tracer coordinate in output file:
    call tracer_coor%Def( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define eigenvalue coordinate in output file:
    call eval_coor%Def( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call XLX_file%Def_Field2D( trim(vname_ps), trim(units_p), &
                                 akstar_coor, ctime_coor, &
                                 ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call XLX_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call XLX_file%Def_IField2D( 'fixed', trim(fixed__units), &
                                 lev_coor, tracer_coor, &
                                 fixed__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call XLX_file%Put_Att( fixed__varid, 'long_name', 'fixed layer (1=fixed, 0=dynamic)', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call XLX_file%Def_Nev( 'nev', '1', &
                              akstar_coor, ctime_coor, &
                              nev__varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call XLX_file%Put_Att( nev__varid, 'long_name', 'number of significant eigenvalues', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call XLX_file%Def_EVal( 'Lambda', '1', &
                              akstar_coor, eval_coor, ctime_coor, &
                              Lambda__varid, status, &
                              fill_value=Lambda__fill_value )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call XLX_file%Put_Att( Lambda__varid, 'long_name', 'eigenvalues of normalized angular averaged covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call XLX_file%Def_EVec( 'X', '1', &
                              akstar_coor, lev_coor, tracer_coor, eval_coor, ctime_coor, &
                              X__varid, status, &
                              fill_value=X__fill_value )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call XLX_file%Put_Att( X__varid, 'long_name', 'eigenvectors of normalized angular averaged covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)

    ! end:
    call XLX_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call akstar_coor%Write( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Write( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Write( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)
    call eval_coor%Write( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Write( XLX_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write:
    call XLX_file%Put_Field2D( ps__varid, ps, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write record
    call XLX_file%Put_IField2D_Series( fixed__varid, itime, fixed, status )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( X(nakstar,nlev,ntracer,nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Lambda(nakstar,nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( nev(nakstar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ time loop

    ! info ...
    write (gol,'("dimensions:")'); call goPr
    write (gol,'("  number of k* values   : ",i0)') nakstar; call goPr
    write (gol,'("  number of time values : ",i0)') ntime; call goPr
    write (gol,'("  number of tracers     : ",i0)') ntracer; call goPr
    
    ! storage for eigenvalue decompo:
    allocate( A(nv,nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( w(nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Z(nv,nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over hours:
    do itime = 1, ntime
    
      ! info ...
      write (gol,'("  itime ",i0)') itime; call goPr

      ! read angular averaged spectral covariances:
      call B_f_file%Get_ACovar( 'B_f', itime, B_f, B_f__units, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! init output:
      X      = X__fill_value
      Lambda = Lambda__fill_value

      ! loop over wavenumbers:
      do ik = 1, nakstar
      
        ! pack covar in matrix array;
        A = -999.9
        ! first dim:
        iv1 = 0
        do itracer1 = 1, ntracer
          do ilev1 = 1, nlev
            iv1 = iv1 + 1
            ! second dim:
            iv2 = 0
            do itracer2 = 1, ntracer
              do ilev2 = 1, nlev
                iv2 = iv2 + 1
                ! copy:
                A(iv1,iv2) = B_f(ik,ilev1,ilev2,itracer1,itracer2)
              end do ! ilev2
            end do ! itracer2
          end do ! ilev1
        end do ! itracer1
        
        ! count number of zero rows:
        nzero = 0
        do iv = 1, nv
          if ( all( A(iv,:) == 0.0 ) ) nzero = nzero + 1
        end do
        
        ! eigenvalue decomposition,
        ! only positive eigenvalues:
        call SyEvr( A, w, Z=Z, m=m, vl=0.0, abstol=1.0e-4, info=status )
        if ( status < 0 ) then
          write (gol,'("SyEvr: illegal value for parameter ",i0)') abs(status); call goErr
          TRACEBACK; status=1; return
        else if ( status > 0 ) then
          write (gol,'("SyEvr: internal error, info=",i0)') status; call goErr
          TRACEBACK; status=1; return
        end if
        
        ! skip nearly zero values:
        if ( m > nv-nzero ) then
          m1 = m - (nv-nzero) + 1
          m2 = m
        else
          m1 = 1
          m2 = m
        end if

        !! testing ...
        !print *, ''
        !print *, 'ik = ', ik, nzero, m, m1, m2
        !print *, '  m    = ', m
        !print *, '  cond = ', w(m)/w(1)
        !print *, '  w    = ', w(1:m)
        !write (*,'("vvv ",i4,2f20.6)') m, w(m)/w(1), w(m)/w(m-37)
        
        ! reset counter:
        m = m2 - m1 + 1
        
        ! testing ...
        !exit
        
        ! number of valid eigenvalues:
        nev(ik) = m
        
        ! loop over valid eigenvalues
        do iv = 1, m
          ! lev/tracer dim:
          iv2 = 0
          do itracer = 1, ntracer
            do ilev = 1, nlev
              iv2 = iv2 + 1
              ! copy eigenvector, reverse order:
              X(ik,ilev,itracer,iv) = Z(iv2,m2+1-iv)
            end do
          end do
          ! copy eigenvalue, reverse order:
          Lambda(ik,iv) = w(m2+1-iv)
        end do  ! eigenvalues
        
      end do  ! ik

      ! write number of eigenvalues:
      call XLX_file%Put_Nev( nev__varid, itime, nev, status )
      IF_NOT_OK_RETURN(status=1)
      ! write eigenvalues:
      call XLX_file%Put_EVal( Lambda__varid, itime, Lambda, status )
      IF_NOT_OK_RETURN(status=1)
      ! write eigenvectors:
      call XLX_file%Put_EVec( X__varid, itime, X, status )
      IF_NOT_OK_RETURN(status=1)
      
      !! testing ...
      !exit
      
    end do  ! times
    
    ! clear:
    deallocate( A, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( w, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( Z, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ done with output
    
    ! clear:
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( X, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( Lambda, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( nev, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! close output:
    call XLX_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ done with input
    
    ! clear:
    deallocate( fixed, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( B_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    call akstar_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close input:
    call B_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ 

    ! ok
    status = 0
  
  end subroutine ModelRuns_Compute_XLX


  ! ***
  
  
  !
  ! Collect all entitities needed for background error factors:
  !    S
  !    E   (dimension of extended grid)
  !    gamma
  !    Lambda
  !    X
  ! and dimensions, coordinates,etc.
  !

  subroutine ModelRuns_Collect_BB( self, rcF, status )
  
    use GO             , only : TrcFile, ReadRc
    use GO             , only : goReplace
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    use C3PO           , only : IntegerCoordinate
    use EMEP_BCovarSqrt, only : T_BCovarSqrt
    
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Collect_BB'
  
    ! --- local ----------------------------------

    character(len=1024)              ::  filename
    type(T_NMC_Output)               ::  inp_file
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(LabelCoordinate)            ::  tracer_coor
    type(TimeCoordinate)             ::  ctime_coor
    type(RealCoordinate)             ::  akstar_coor
    integer                          ::  nlon, nlat, nlev
    integer                          ::  ntime
    integer                          ::  nakstar
    integer                          ::  ntracer
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    real, allocatable                ::  ps(:,:,:)   ! (nlon,nlat,ntime)
    integer                          ::  nlon_ex, nlat_ex
    type(Dimension)                  ::  lon_ex_dim, lat_ex_dim
    integer, allocatable             ::  fixed(:,:)  ! (nlev,nvar)
    character(len=32)                ::  fixed__units
    integer                          ::  fixed__varid
    character(len=32)                ::  vname, units
    real, allocatable                ::  S(:,:,:,:,:)  ! (nlon,nlat,nlev,ntime,nvar)
    integer, allocatable             ::  S__varids(:)
    character(len=32), allocatable   ::  S__units(:)
    integer                          ::  itime
    integer                          ::  ivar
    integer                          ::  itracer
    integer                          ::  ilev
    
    type(T_NMC_Output)               ::  BB_file
    character(len=1024)              ::  msg
    integer                          ::  ps__varid
    real, allocatable                ::  gamma(:,:,:,:)  ! (nakstar,nlev,ntr,ntime)
    integer                          ::  gamma__varid
    character(len=32)                ::  gamma__units
    type(IntegerCoordinate)          ::  eval_coor
    integer                          ::  nv
    integer, allocatable             ::  nev(:,:)       ! (nakstar)
    character(len=32)                ::  nev__units
    integer                          ::  nev__varid
    real, allocatable                ::  Lambda(:,:,:)  ! (nakstar,nv)
    character(len=32)                ::  Lambda__units
    real                             ::  Lambda__fill_value
    integer                          ::  Lambda__varid
    real, allocatable                ::  X(:,:,:,:,:)   ! (nakstar,nlev,ntr,nv)
    character(len=32)                ::  X__units
    real                             ::  X__fill_value
    integer                          ::  X__varid
    
    real, allocatable                ::  phi(:,:,:)   ! (nlev,ntracer,ntime)
    character(len=32)                ::  phi__units
    integer                          ::  phi__varid
    type(T_BCovarSqrt)               ::  BCovarSqrt
    integer                          ::  ilon, ilat
    real, allocatable                ::  x1(:,:,:,:)  ! (nlon,nlat,nlev,ntracer)
    real, allocatable                ::  Cx1(:,:,:,:)  ! (nlon,nlat,nlev,ntracer)
    real                             ::  c1
    integer                          ::  nw
    complex, allocatable             ::  w(:)  ! (nw)
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Collect B^{1/2} B^{H/2} entities **")'); call goPr
    write (gol,'("")'); call goPr

    ! ~ eps statistics:
    
    ! info ...
    write (gol,'("read eps stats ...")'); call goPr
          
    ! open file with eps samples:
    call inp_file%Open( trim(self%eps_stats_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create longitude coordinate:
    call lon_coor%Init( 'longitude', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lon_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lon_coor%Get_Dim( status, n=nlon )
    IF_NOT_OK_RETURN(status=1)

    ! create latitude coordinate:
    call lat_coor%Init( 'latitude', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lat_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lat_coor%Get_Dim( status, n=nlat )
    IF_NOT_OK_RETURN(status=1)

    ! create level coordinate:
    call lev_coor%Init( 'lev', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lev_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lev_coor%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)
    ! get name of pressure variable:
    call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p )
    IF_NOT_OK_RETURN(status=1)
    
    ! create hour coordinate:
    call ctime_coor%Init( 'time', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file, here includes climatology bounds:
    call ctime_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call ctime_coor%Get_Dim( status, n=ntime )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( ps(nlon,nlat,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( fixed(nlev,self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( S(nlon,nlat,nlev,ntime,self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( S__units(self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over times:
    do itime = 1, ntime
      ! read surface pressure:
      call inp_file%Get_Field2D_Series( trim(vname_ps), itime, ps(:,:,itime), units, status )
      IF_NOT_OK_RETURN(status=1)
    end do

    ! loop over variables:
    do ivar = 1, self%nvar
      !
      ! variable:
      write (vname,'(a,"_stdv")') trim(self%var(ivar)%name)
      ! loop over times:
      do itime = 1, ntime
        ! read:
        call inp_file%Get_Field3D_Series( trim(vname), itime, &
                             S(:,:,:,itime,ivar), S__units(ivar), status )
        IF_NOT_OK_RETURN(status=1)
      end do ! time
      !
      ! variable:
      write (vname,'(a,"_fixed")') trim(self%var(ivar)%name)
      call inp_file%Get_IField1D( trim(vname), fixed(:,ivar), fixed__units, status )
      IF_NOT_OK_RETURN(status=1)
      !
    end do ! var
    
    ! close:
    call inp_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ extended grid
    
    ! info ...
    write (gol,'("read extended grid dimenions ...")'); call goPr
        
    ! input file:
    filename = trim(self%eta_f_samples_filename)
    ! replace keys:
    call goReplace( filename, '%{yyyy}', '(i4.4)', self%tday1%year , status )
    IF_NOT_OK_RETURN(status=1)
    call goReplace( filename, '%{mm}'  , '(i2.2)', self%tday1%month, status )
    IF_NOT_OK_RETURN(status=1)
    call goReplace( filename, '%{dd}'  , '(i2.2)', self%tday1%day  , status )
    IF_NOT_OK_RETURN(status=1)
          
    ! open file with scaling:
    call inp_file%Open( trim(filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! create extended dimension:
    call lon_ex_dim%Init( 'lon_ex', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lon_ex_dim%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! create extended dimension:
    call lat_ex_dim%Init( 'lat_ex', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lat_ex_dim%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! close:
    call inp_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ gamma
    
    ! info ...
    write (gol,'("read gamma ...")'); call goPr
          
    ! open file with scaling:
    call inp_file%Open( trim(self%gamma_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create average k* coordinate:
    call akstar_coor%Init( 'akstar', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call akstar_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call akstar_coor%Get_Dim( status, n=nakstar )
    IF_NOT_OK_RETURN(status=1)

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=ntracer )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( gamma(nakstar,nlev,ntracer,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over times:
    do itime = 1, ntime
      ! read scaling:
      call inp_file%Get_AVar( 'gamma', itime, gamma(:,:,:,itime), gamma__units, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! close:
    call inp_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ XLX
    
    ! info ...
    write (gol,'("read XLX ...")'); call goPr
          
    ! open file with scaling:
    call inp_file%Open( trim(self%XLX_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create eigenvalue coordinate:
    call eval_coor%Init( 'eigenvalue', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call eval_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call eval_coor%Get_Dim( status, n=nv )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( X(nakstar,nlev,ntracer,nv,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Lambda(nakstar,nv,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( nev(nakstar,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over times:
    do itime = 1, ntime
      ! read number of significant eigenvalues:
      call inp_file%Get_Nev( 'nev', itime, nev(:,itime), nev__units, status )
      IF_NOT_OK_RETURN(status=1)
      ! read eigenvalues:
      call inp_file%Get_EVal( 'Lambda', itime, Lambda(:,:,itime), &
                                Lambda__units, status, &
                                fill_value=Lambda__fill_value )
      IF_NOT_OK_RETURN(status=1)
      ! read eigenvectors:
      call inp_file%Get_EVec( 'X', itime, X(:,:,:,:,itime), &
                                X__units, status, &
                                fill_value=X__fill_value )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! close:
    call inp_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ phi
    
    ! storage for compenstation factors:
    allocate( phi(nlev,ntracer,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ output
    
    ! info ...
    write (gol,'("write BB file ...")'); call goPr
    
    ! create new result file:
    call BB_file%Create( trim(self%BB_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("covariance entities from `",a,"`, `",a,"`, `",a,"`, and `",a,"`")') &
            trim(self%eps_stats_filename), &
            trim(self%eta_f_samples_filename), &
            trim(self%gamma_filename), &
            trim(self%XLX_filename)
    ! add:
    call BB_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define coordinates in output file:
    call lon_coor%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define dimensions in output file:
    call lon_ex_dim%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_ex_dim%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define level coordinate in output file:
    call lev_coor%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define ak* coordinate in output file:
    call akstar_coor%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define tracer coordinate in output file:
    call tracer_coor%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define eigenvalue coordinate in output file:
    call eval_coor%Def( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call BB_file%Def_Field2D_Series( trim(vname_ps), trim(units_p), &
                                  lon_coor, lat_coor, ctime_coor, &
                                  ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call BB_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( S__varids(self%nvar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over variables:
    do ivar = 1, self%nvar
      ! variable:
      write (vname,'("S_",a)') trim(self%var(ivar)%name)
      ! define:
      call BB_file%Def_Field3D_Series( trim(vname), trim(S__units(ivar)), &
                              lon_coor, lat_coor, lev_coor, ctime_coor, &
                              S__varids(ivar), status )
      IF_NOT_OK_RETURN(status=1)
      call BB_file%Put_Att( S__varids(ivar), 'long_name', trim(self%var(ivar)%name)//' grid point standard deviation', status )
      IF_NOT_OK_RETURN(status=1)
    end do ! var

    ! define:
    call BB_file%Def_IField2D( 'fixed', trim(fixed__units), &
                                  lev_coor, tracer_coor, &
                                  fixed__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call BB_file%Put_Att( fixed__varid, 'long_name', 'fixed layer (1=fixed, 0=dynamic)', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call BB_file%Def_AVar( 'gamma', trim(gamma__units), &
                              akstar_coor, lev_coor, tracer_coor, ctime_coor, &
                              gamma__varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call BB_file%Put_Att( gamma__varid, 'long_name', 'scaling for angular average of covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call BB_file%Def_Nev( 'nev', trim(nev__units), &
                              akstar_coor, ctime_coor, &
                              nev__varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call BB_file%Put_Att( nev__varid, 'long_name', 'number of significant eigenvalues', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call BB_file%Def_EVal( 'Lambda', trim(Lambda__units), &
                              akstar_coor, eval_coor, ctime_coor, &
                              Lambda__varid, status, &
                              fill_value=Lambda__fill_value )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call BB_file%Put_Att( Lambda__varid, 'long_name', 'eigenvalues of normalized angular averaged covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call BB_file%Def_EVec( 'X', trim(X__units), &
                              akstar_coor, lev_coor, tracer_coor, eval_coor, ctime_coor, &
                              X__varid, status, &
                              fill_value=X__fill_value )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call BB_file%Put_Att( X__varid, 'long_name', 'eigenvectors of normalized angular averaged covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call BB_file%Def_Field2D_Series( 'phi', '1', &
                                  lev_coor, tracer_coor, ctime_coor, &
                                  phi__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call BB_file%Put_Att( phi__varid, 'long_name', 'compensation factor for truncations', status )
    IF_NOT_OK_RETURN(status=1)

    ! end:
    call BB_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lon_coor%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write extended dimensions:
    call lon_ex_dim%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_ex_dim%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lev_coor%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call ctime_coor%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call akstar_coor%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call tracer_coor%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call eval_coor%Write( BB_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over hours:
    do itime = 1, ntime
      ! write:
      call BB_file%Put_Field2D_Series( ps__varid, itime, ps(:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! times

    ! loop over variables:
    do ivar = 1, self%nvar
      ! loop over hours:
      do itime = 1, ntime
        ! write; note different writing order!
        call BB_file%Put_Field3D_Series( S__varids(ivar), itime, S(:,:,:,itime,ivar), status )
        IF_NOT_OK_RETURN(status=1)
      end do  ! times
    end do ! vars
    
    ! write:
    call BB_file%Put_IField2D( fixed__varid, fixed, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over times:
    do itime = 1, ntime
      ! write scaling:
      call BB_file%Put_AVar( gamma__varid, itime, gamma(:,:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! times

    ! loop over times:
    do itime = 1, ntime
      ! write number of eigenvalues:
      call BB_file%Put_Nev( nev__varid, itime, nev(:,itime), status )
      IF_NOT_OK_RETURN(status=1)
      ! write eigenvalues:
      call BB_file%Put_EVal( Lambda__varid, itime, Lambda(:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
      ! write eigenvectors:
      call BB_file%Put_EVec( X__varid, itime, X(:,:,:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! times
    
    ! loop over hours:
    do itime = 1, ntime
      ! init as unity factor, final value will be evaluated afterwards:
      phi(:,:,itime) = 1.0
      ! write:
      call BB_file%Put_Field2D_Series( phi__varid, itime, phi(:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! times

    ! close output:
    call BB_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ done with input
      
    ! clear:
    call lon_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    call lon_ex_dim%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_ex_dim%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( gamma, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( S, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( S__units, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( S__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    call akstar_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    !
    ! ~ compensation factors
    !
    
    ! info ...
    write (gol,'("  compute compensation factors phi ...")'); call goPr

    ! init covariance structure:
    call BCovarSqrt%Init( trim(self%BB_filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! sample cell:
    ilon = nlon/2
    ilat = nlat/2

    ! storage for input state:
    allocate( x1(nlon,nlat,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for evaluated correlations:
    allocate( Cx1(nlon,nlat,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over hours:
    do itime = 1, ntime
    
      ! info ...
      write (gol,'("    time record ",i0)') itime; call goPr

      ! size:
      call BCovarSqrt%TimeGet( itime, status, nw=nw )
      IF_NOT_OK_RETURN(status=1)
      ! storage for complex 1D state:
      allocate( w(nw), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! loop over tracers:
      do itracer = 1, ntracer
    
        ! info ...
        write (gol,'("      tracer ",i0," / ",i0)') itracer, ntracer; call goPr
      
        ! loop over levels:
        do ilev = 1, nlev
        
          ! fixed layer?
          if ( fixed(ilev,itracer) == 1 ) then
          
            ! no scaling needed, set to default value:
            phi(ilev,itracer,itime) = 1.0
            
          else
    
            ! fill as unity vector:
            x1 = 0.0
            x1(ilon,ilat,ilev,itracer) = 1.0
          
            ! evaluate square roots:
            call BCovarSqrt%Forward( itime, x1, w, status, correlation=.true. )
            IF_NOT_OK_RETURN(status=1)
            call BCovarSqrt%Reverse( itime, w, Cx1, status, correlation=.true. )
            IF_NOT_OK_RETURN(status=1)
            
            ! top:
            c1 = Cx1(ilon,ilat,ilev,itracer)
    
            ! info ...
            write (gol,'("          level ",i3," : ",f8.2," max ",f8.2)') &
                             ilev, c1, maxval(Cx1); call goPr

            ! check ...
            if ( c1 <= 0.0 ) then
              write (gol,'("evaluated correlation should be positive, found ",es16.6)') c1; call goErr
              TRACEBACK; status=1; return
            end if
            
            ! top of correlation should become 1.0,
            ! define scale factor to compensate for deviation:
            phi(ilev,itracer,itime) = 1.0 / c1

          end if  ! fixed layer
      
        end do ! levels
        
      end do  ! tracers
      
      ! clear:
      deallocate( w, stat=status )
      IF_NOT_OK_RETURN(status=1)
    
    end do ! hours
    
    ! clear:
    deallocate( x1, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( Cx1, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    call BCovarSqrt%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! data file contains initial version filled with 1.0,
    ! now replace by final version ...
    
    ! open existing file:
    call BB_file%Open( trim(self%BB_filename), status, writable=.true. )
    IF_NOT_OK_RETURN(status=1)
    
    ! access to variable:
    call BB_file%Get_VarID( 'phi', phi__varid, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over hours:
    do itime = 1, ntime
      ! write:
      call BB_file%Put_Field2D_Series( phi__varid, itime, phi(:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! times

    ! close output:
    call BB_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( phi, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ 
    
    ! clear:
    deallocate( fixed, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ 

    ! ok
    status = 0
  
  end subroutine ModelRuns_Collect_BB


  ! ***
  
  
  !
  ! Evaluate parameterized covariance matrix at selected locations.
  !

  subroutine ModelRuns_Evaluate_BB( self, rcF, status )
  
    use GO             , only : TrcFile, ReadRc
    use EMEP_NMC_Output, only : T_NMC_Output
    use EMEP_BCovarSqrt, only : T_BCovarSqrt
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_ModelRuns), intent(inout)           ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelRuns_Evaluate_BB'
  
    ! --- local ----------------------------------
    
    type(T_BCovarSqrt)               ::  BCovarSqrt
    logical                          ::  correlation
    character(len=32)                ::  cv_label
    character(len=1)                 ::  cv_id
    character(len=32)                ::  cv_units

    type(T_NMC_Output)               ::  inp_file
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(LabelCoordinate)            ::  tracer_coor
    type(TimeCoordinate)             ::  ctime_coor
    integer                          ::  nlon, nlat, nlev
    integer                          ::  ntime
    integer                          ::  ntracer
    character(len=32)                ::  tracer_name
    character(len=32)                ::  tracer_units
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    real, allocatable                ::  ps(:,:,:)   ! (nlon,nlat,ntime)
    character(len=32)                ::  vname, units
    real, allocatable                ::  S(:,:,:,:)  ! (nlon,nlat,nlev,ntracer)
    integer, allocatable             ::  S__varids(:)
    integer                          ::  itime
    integer                          ::  ipoint
    integer                          ::  ilon, ilat, ilev
    integer                          ::  itracer, itracer2
    
    character(len=1024)              ::  B_point_filename
    type(T_NMC_Output)               ::  out_file
    character(len=1024)              ::  msg
    integer                          ::  ps__varid
    real, allocatable                ::  x(:,:,:,:)  ! (nlon,nlat,nlev,ntracer)
    real, allocatable                ::  Bx(:,:,:,:)  ! (nlon,nlat,nlev,ntracer)
    integer, allocatable             ::  Bx__varids(:,:)  ! (npoint,ntracer)

    integer                          ::  nw
    complex, allocatable             ::  w(:)
    
    logical                          ::  levtr
    real, allocatable                ::  Bzs(:,:,:,:)  ! (nlev,nlev,ntracer,ntracer)
    integer, allocatable             ::  Bzs__varids(:)  ! (npoint)
    integer                          ::  k
    
    ! --- begin ----------------------------------

    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Evaluate B e_i **")'); call goPr
    write (gol,'("")'); call goPr
    
    ! ~ covar sqrt
    
    ! info ...
    write (gol,'("read B sqrt ...")'); call goPr

    ! init covariance structure:
    call BCovarSqrt%Init( trim(self%BB_filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! flag:
    call ReadRc( rcF, 'nmc.approx-B.correlation', correlation, status )
    IF_NOT_OK_RETURN(status=1)
    ! set label:
    if ( correlation ) then
      cv_label = 'correlation'
      cv_id    = 'C'
      cv_units = '1'
    else
      cv_label = 'covariance'
      cv_id    = 'B'
      cv_units = 'tracer*tracer'
    end if
    
    ! flag:
    call ReadRc( rcF, 'nmc.approx-B.levtr', levtr, status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ coordiantes etc:
    
    ! info ...
    write (gol,'("read coordinates ...")'); call goPr
          
    ! open file with covar sqrt:
    call inp_file%Open( trim(self%BB_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create longitude coordinate:
    call lon_coor%Init( 'longitude', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lon_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lon_coor%Get_Dim( status, n=nlon )
    IF_NOT_OK_RETURN(status=1)

    ! create latitude coordinate:
    call lat_coor%Init( 'latitude', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lat_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lat_coor%Get_Dim( status, n=nlat )
    IF_NOT_OK_RETURN(status=1)

    ! create level coordinate:
    call lev_coor%Init( 'lev', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lev_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lev_coor%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)
    ! get name of pressure variable:
    call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p )
    IF_NOT_OK_RETURN(status=1)
    
    ! create hour coordinate:
    call ctime_coor%Init( 'time', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file, here includes climatology bounds:
    call ctime_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call ctime_coor%Get_Dim( status, n=ntime )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( ps(nlon,nlat,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over times:
    do itime = 1, ntime
      ! read surface pressure:
      call inp_file%Get_Field2D_Series( trim(vname_ps), itime, ps(:,:,itime), units, status )
      IF_NOT_OK_RETURN(status=1)
    end do

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=ntracer )
    IF_NOT_OK_RETURN(status=1)
      
    ! close:
    call inp_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ points
    
    ! loop over points:
    do ipoint = 1, self%npoint
      ! nearby index:
      call lon_coor%Get_Index( self%point(ipoint)%lon, self%point(ipoint)%ilon, status )
      IF_NOT_OK_RETURN(status=1)
      ! nearby index:
      call lat_coor%Get_Index( self%point(ipoint)%lat, self%point(ipoint)%ilat, status )
      IF_NOT_OK_RETURN(status=1)
      ! surface (top-down order!)
      self%point(ipoint)%ilev = nlev
    end do
    
    ! ~ output
    
    ! info ...
    write (gol,'("create B_point file ...")'); call goPr
    
    ! target file:
    call ReadRc( rcF, 'nmc.approx-B.filename', B_point_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! create new result file:
    call out_file%Create( trim(B_point_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("evaluate ",a," from `",a,"`")') &
            trim(cv_label), trim(self%BB_filename)
    ! add:
    call out_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define coordinates in output file:
    call lon_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define level coordinate in output file:
    call lev_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define tracer coordinate in output file:
    call tracer_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call out_file%Def_Field2D_Series( trim(vname_ps), trim(units_p), &
                                  lon_coor, lat_coor, ctime_coor, &
                                  ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call out_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( S__varids(ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bx__varids(self%npoint,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over tracers:
    do itracer = 1, ntracer

      ! name and units:
      call BCovarSqrt%TracerGet( itracer, status, &
                                 name=tracer_name, units=tracer_units )
      IF_NOT_OK_RETURN(status=1)
      
      ! variable:
      write (vname,'("S_",a)') trim(tracer_name)
      ! define:
      call out_file%Def_Field3D_Series( trim(vname), trim(tracer_units), &
                              lon_coor, lat_coor, lev_coor, ctime_coor, &
                              S__varids(itracer), status )
      IF_NOT_OK_RETURN(status=1)
      call out_file%Put_Att( S__varids(itracer), 'long_name', &
                     trim(tracer_name)//' grid point standard deviation', status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over points:
      do ipoint = 1, self%npoint

        ! variable:
        write (vname,'(a,"_",a,"_",a)') trim(cv_id), trim(tracer_name), trim(self%point(ipoint)%name)
        ! define:
        call out_file%Def_Sample3D( trim(vname), trim(cv_units), &
                                lon_coor, lat_coor, lev_coor, ctime_coor, tracer_coor, &
                                Bx__varids(ipoint,itracer), status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'long_name', &
                       trim(cv_label)//' with '//trim(tracer_name)//' in '//trim(self%point(ipoint)%name), status )
        IF_NOT_OK_RETURN(status=1)
        ! write as attributes:
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'lon', self%point(ipoint)%lon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'lat', self%point(ipoint)%lat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'ilon', self%point(ipoint)%ilon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'ilat', self%point(ipoint)%ilat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'ilev', self%point(ipoint)%ilev, status )
        IF_NOT_OK_RETURN(status=1)

      end do ! points

    end do ! tracers
    
    ! evaluate level/tracer covariance?
    if ( levtr ) then

      ! storage for result:
      allocate( Bzs(nlev,nlev,ntracer,ntracer), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! output:
      allocate( Bzs__varids(self%npoint), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over points:
      do ipoint = 1, self%npoint

        ! variable:
        write (vname,'(a,"zs_",a)') trim(cv_id), trim(self%point(ipoint)%name)
        ! define:
        call out_file%Def_ACovar1( trim(vname), trim(cv_units), &
                                    lev_coor, tracer_coor, ctime_coor, &
                                    Bzs__varids(ipoint), status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'long_name', &
                       'vertical '//trim(cv_label)//' in '//trim(self%point(ipoint)%name), status )
        IF_NOT_OK_RETURN(status=1)
        ! write as attributes:
        call out_file%Put_Att( Bzs__varids(ipoint), 'lon', self%point(ipoint)%lon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'lat', self%point(ipoint)%lat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'ilon', self%point(ipoint)%ilon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'ilat', self%point(ipoint)%ilat, status )
        IF_NOT_OK_RETURN(status=1)
        
      end do  ! points
      
    end if  ! level/tracer covar

    ! end:
    call out_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lon_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lev_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call ctime_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call tracer_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over hours:
    do itime = 1, ntime
      ! write:
      call out_file%Put_Field2D_Series( ps__varid, itime, ps(:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! times

    ! storage for std.dev.
    allocate( S(nlon,nlat,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for real state:
    allocate( x(nlon,nlat,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bx(nlon,nlat,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over times:
    do itime = 1, ntime

      ! info ...
      write (gol,'("  time ",i0," / ",i0," ...")') itime, ntime; call goPr
      
      ! extract:
      call BCovarSqrt%TimeGet( itime, status, S=S )
      IF_NOT_OK_RETURN(status=1)
      ! loop over tracers:
      do itracer = 1, ntracer
        ! write sample:
        call out_file%Put_Field3D_Series( S__varids(itracer), itime, &
                                           S(:,:,:,itracer), status )
        IF_NOT_OK_RETURN(status=1)
      end do

      ! count:
      call BCovarSqrt%TimeGet( itime, status, nw=nw )
      IF_NOT_OK_RETURN(status=1)

      ! storage:
      if ( allocated(w) .and. (size(w) /= nw) ) then
        deallocate( w, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      if ( .not. allocated(w) ) then
        allocate( w(nw), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! loop over selected points:
      do ipoint = 1, self%npoint

        ! info ...
        write (gol,'("    point ",i0," / ",i0," ...")') ipoint, self%npoint; call goPr
      
        ! grid cell index:
        ilon = self%point(ipoint)%ilon
        ilat = self%point(ipoint)%ilat
        ilev = self%point(ipoint)%ilev
        
        ! loop over tracers (observed in point):
        do itracer = 1, ntracer

          ! info ...
          write (gol,'("      tracer ",i0," / ",i0," ...")') itracer, ntracer; call goPr

          ! unity vector:
          x = 0.0
          x(ilon,ilat,ilev,itracer) = 1.0

          ! evaluate:  B x =  B^{1/2} B^{H/2} x
          call BCovarSqrt%Forward( itime, x, w, status, correlation=correlation )
          IF_NOT_OK_RETURN(status=1)
          call BCovarSqrt%Reverse( itime, w, Bx, status, correlation=correlation )
          IF_NOT_OK_RETURN(status=1)
          
          ! info ...
          write (gol,'("        value range : ",f8.2," - ",f8.2)') minval(Bx), maxval(Bx); call goPr
          
          ! loop over correlated tracers:
          do itracer2 = 1, ntracer
            ! write sample:
            call out_file%Put_Sample3D( Bx__varids(ipoint,itracer), itime, itracer2, &
                                               Bx(:,:,:,itracer2), status )
            IF_NOT_OK_RETURN(status=1)
          end do  ! correlated tracers
          
        end do ! tracers
        
        ! ~
        
        ! evaluate full level/tracer covar?
        if ( levtr ) then

          ! info ...
          write (gol,'("      evaluate level/tracer covariance ...")'); call goPr

          ! loop over tracers (observed in point):
          do itracer = 1, ntracer
            ! loop:
            do k = 1, nlev
              ! info ...
              write (gol,'("        level ",i0," / ",i0)') k, nlev; call goPr
              ! unity vector:
              x = 0.0
              x(ilon,ilat,k,itracer) = 1.0
              ! evaluate:  B x =  B^{1/2} B^{H/2} x
              call BCovarSqrt%Forward( itime, x, w, status, correlation=correlation )
              IF_NOT_OK_RETURN(status=1)
              call BCovarSqrt%Reverse( itime, w, Bx, status, correlation=correlation )
              IF_NOT_OK_RETURN(status=1)
              ! info ...
              write (gol,'("          value range : ",f8.2," - ",f8.2)') minval(Bx), maxval(Bx); call goPr
              ! copy:
              Bzs(k,:,itracer,:) = Bx(ilon,ilat,:,:)
              Bzs(:,k,:,itracer) = Bx(ilon,ilat,:,:)
            end do ! levels
          end do ! tracers

          ! write sample:
          call out_file%Put_ACovar1( Bzs__varids(ipoint), itime, Bzs, status  )
          IF_NOT_OK_RETURN(status=1)

        end if ! level/tracer covariance
        
        ! ~
          
      end do  ! points

    end do ! times
    
    ! clear:
    deallocate( w, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( Bx, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( S, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! evaluate level/tracer covariance?
    if ( levtr ) then
      ! clear:
      deallocate( Bzs, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! clear:
      deallocate( Bzs__varids, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! level/tracer covar

    ! ~ done with input
      
    ! clear:
    call lon_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ done with output

    ! close output:
    call out_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    deallocate( Bx__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( S__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ done with covar
    
    ! clear:
    call BCovarSqrt%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ 

    ! ok
    status = 0
  
  end subroutine ModelRuns_Evaluate_BB



end module EMEP_NMC_Driver_ALL

