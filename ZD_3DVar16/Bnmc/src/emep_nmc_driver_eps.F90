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

module EMEP_NMC_Driver_Eps

  use GO                  , only : gol, goPr, goErr
  use GO                  , only : TDate
  use EMEP_NMC_Common     , only : nrun
  use EMEP_NMC_Common     , only : T_ModelVars
  use EMEP_NMC_ModelOutput, only : T_ModelOutputSeries

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_Eps
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_Eps'
  
  ! assumed maximum number of variables:
  integer, parameter  ::  maxvar = 10

  ! assumed maximum number of times within days:
  integer, parameter  ::  maxhours = 24


  ! --- types ----------------------------------------
  
  ! *
  
  type T_Driver_Eps
    ! model run info:
    type(T_ModelOutputSeries)         ::  run(nrun)
    ! time range:
    type(TDate)                       ::  tday1, tday2
    ! hours within day:
    integer                           ::  nhour
    integer                           ::  hours(maxhours)
    ! variables:
    type(T_ModelVars)                 ::  mvars
    ! statistics files created or read:
    character(len=1024)               ::  eps_samples_filename
    character(len=1024)               ::  eps_stats_filename
    !! testing ...
    !logical                           ::  test_shift
  contains
    procedure   ::  Init            => Driver_Eps_Init
    procedure   ::  Done            => Driver_Eps_Done
    procedure   ::  Compute         => Driver_Eps_Compute
  end type T_Driver_Eps
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_Eps_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine, goSplitString
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eps), intent(out)           ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eps_Init'
    
    ! --- local ----------------------------------
    
    integer               ::  irun
    character(len=32)     ::  rckey
    character(len=32)     ::  id
    character(len=32)     ::  tvalue
    character(len=1024)   ::  line
    
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
    call self%mvars%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    
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
      call self%run(irun)%Init( rcF, 'nmc.run', id, self%mvars, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! input/output:
    call ReadRc( rcF, 'nmc.eps.samples.filename', self%eps_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.eps.stats.filename', self%eps_stats_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    !! test ...
    !call ReadRc( rcF, 'nmc.eps.test.shift', self%test_shift, status )
    !IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine Driver_Eps_Init


  ! ***
  

  subroutine Driver_Eps_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eps), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eps_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! loop over runs:
    do irun = 1, nrun
      ! done with runs:
      call self%run(irun)%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! clear:
    call self%mvars%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_Eps_Done


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

  subroutine Driver_Eps_Compute( self, status )
  
    use GO             , only : TDate, IncrDate
    use GO             , only : operator(+), operator(>), wrtgol
    use GO             , only : T_SampStat
    use GO             , only : goReplace
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : RealCoordinate, HybrideLevelCoordinate, TimeCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eps), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eps_Compute'
    
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
    real, allocatable                ::  eps_g2    (:,:,:)      ! (nlon,nlat,nlev)

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
    
        ! >> tst >>>>>>>>>
        !integer       ::  shift
        !integer       ::  ip
    
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
            do ivar = 1, self%mvars%n
              ! set or check units:
              if ( irun == 1 ) then
                ! copy units:
                self%mvars%value(ivar)%units = trim(self%run(irun)%var(ivar)%units)
              else
                ! check ...
                if ( trim(self%mvars%value(ivar)%units) /= trim(self%run(irun)%var(ivar)%units) ) then
                  write (gol,'("units do not match for ",i0," ",a)') ivar, trim(self%mvars%value(ivar)%name); call goErr
                  write (gol,'("  first run   : ",a)')  trim(self%mvars%value(ivar)%units); call goErr
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
          allocate( samp_ps   (nlon,nlat                         ,nrun), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( samp_data (nlon,nlat,nlev,self%mvars%n,nrun), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( aver_ps   (nlon,nlat     ), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( eps_g     (nlon,nlat,nlev), stat=status )
          IF_NOT_OK_RETURN(status=1)
          allocate( eps_g2    (nlon,nlat,nlev), stat=status )
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
          allocate( eps__varids(self%mvars%n) )

          ! stats:
          allocate( sampstat_eps_g(self%nhour,self%mvars%n), stat=status )
          IF_NOT_OK_RETURN(status=1)
          do ihr = 1, self%nhour
            do ivar = 1, self%mvars%n
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
        call goReplace( filename, '%{hh}'  , '(i2.2)', tday%hour , status )
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
          do ivar = 1, self%mvars%n
            ! define:
            call eps_file%Def_Sample3D( trim(self%mvars%value(ivar)%name), &
                                    trim(self%mvars%value(ivar)%units), &
                                    lon_coor, lat_coor, lev_coor, time_coor, sample_coor, &
                                    eps__varids(ivar), status )
            IF_NOT_OK_RETURN(status=1)
            call eps_file%Put_Att( eps__varids(ivar), 'long_name', trim(self%mvars%value(ivar)%name), status )
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
          do ivar = 1, self%mvars%n
            ! read data:
            call self%run(irun)%GetRecord( t, samp_ps(:,:,irun), ivar, &
                                            samp_data(:,:,:,ivar,irun), &
                                            self%mvars%value(ivar)%units, status )
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
        do ivar = 1, self%mvars%n
          ! difference field:
          eps_g = samp_data(:,:,:,ivar,2) - samp_data(:,:,:,ivar,1)

            !! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            !if ( self%test_shift ) then
            !  ! save original:
            !  eps_g2 = eps_g
            !  ! test different shifts per hour:
            !  if ( ihour == 1 ) then
            !    ! correlated in x direction:
            !    do ilev = 1, nlev
            !      do ilat = 1, nlat
            !        eps_g(:,ilat,ilev) = sum(eps_g2(:,ilat,ilev))/size(eps_g2(:,ilat,ilev))
            !      end do ! ilat
            !    end do ! ilev
            !  else if ( ihour == 2 ) then
            !    ! correlated in y direction:
            !    do ilev = 1, nlev
            !      do ilon = 1, nlon
            !        eps_g(ilon,:,ilev) = sum(eps_g2(ilon,:,ilev))/size(eps_g2(ilon,:,ilev))
            !      end do ! ilon
            !    end do ! ilev
            !  else if ( ihour == 3 ) then
            !    ! correlated in xy direction:
            !    do ilev = 1, nlev
            !      do ilat = 1, nlat
            !        do ilon = 1, nlon
            !          ! diagonal location:
            !          ip = min( min( max( 1, (ilon+ilat)/2 ), nlon ), nlat )
            !          eps_g(ilon,ilat,ilev) = eps_g2(ip,ip,ilev)
            !        end do ! ilon
            !      end do ! ilat
            !    end do ! ilev
            !  else if ( ihour == 4 ) then
            !    ! correlated in xy direction:
            !    do ilev = 1, nlev
            !      do ilat = 1, nlat
            !        do ilon = 1, nlon
            !          ! diagonal location:
            !          ip = min( min( max( 1, max(ilon,ilat) - min(ilon,ilat) ), nlon ), nlat )
            !          eps_g(ilon,ilat,ilev) = eps_g2(ip,ip,ilev)
            !        end do ! ilon
            !      end do ! ilat
            !    end do ! ilev
            !  end if ! ihour
            !end if ! test shift
            !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
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
      deallocate( samp_ps, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( samp_data, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( aver_ps, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( eps_g, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( eps_g2, stat=status )
      IF_NOT_OK_RETURN(status=1)

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
      allocate( dmean__varid(self%mvars%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( dstdv__varid(self%mvars%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( fixed__varids(self%mvars%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      
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
      do ivar = 1, self%mvars%n
        ! define:
        call eps_file%Def_Field3D_Series( trim(self%mvars%value(ivar)%name)//'_mean', &
                               trim(self%mvars%value(ivar)%units), &
                               lon_coor, lat_coor, lev_coor, ctime_coor, &
                               dmean__varid(ivar), status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( dmean__varid(ivar), 'long_name', trim(self%mvars%value(ivar)%name)//' mean', status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( dmean__varid(ivar), 'cell_methods', 'time: point within days time: mean over days', status )
        IF_NOT_OK_RETURN(status=1)
        ! define:
        call eps_file%Def_Field3D_Series( trim(self%mvars%value(ivar)%name)//'_stdv', &
                               trim(self%mvars%value(ivar)%units), &
                               lon_coor, lat_coor, lev_coor, ctime_coor, &
                               dstdv__varid(ivar), status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( dstdv__varid(ivar), 'long_name', trim(self%mvars%value(ivar)%name)//' standard deviation', status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( dstdv__varid(ivar), 'cell_methods', 'time: point within days time: standard_deviation over days', status )
        IF_NOT_OK_RETURN(status=1)
        ! define:
        call eps_file%Def_IField1D( trim(self%mvars%value(ivar)%name)//'_fixed', '1', &
                               lev_coor, fixed__varids(ivar), status )
        IF_NOT_OK_RETURN(status=1)
        call eps_file%Put_Att( fixed__varids(ivar), 'long_name', trim(self%mvars%value(ivar)%name)//' fixed layer (1=fixed, 0=dynamic)', status )
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
        do ivar = 1, self%mvars%n

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
      deallocate( ps, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dmean, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dstdv, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( fixed, stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! ~

      ! info ...
      write (gol,'("done.")'); call goPr

      ! clear:
      deallocate( lons, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( lats, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( hyam, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( hybm, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( hyai, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( hybi, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
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
      deallocate( eps__varids, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dmean__varid, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( dstdv__varid, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( fixed__varids, stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! done with stats:
      do ihr = 1, self%nhour
        call sampstat_ps(ihr)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do
      deallocate( sampstat_ps, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! done with stats:
      do ihr = 1, self%nhour
        do ivar = 1, self%mvars%n
          call sampstat_eps_g(ihr,ivar)%Done( status )
          IF_NOT_OK_RETURN(status=1)
        end do
      end do
      deallocate( sampstat_eps_g, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    else
    
      ! info ...
      write (gol,'("WARNING - no records in time range ...")'); call goPr

    end if
    
    ! ok
    status = 0
  
  end subroutine Driver_Eps_Compute


end module EMEP_NMC_Driver_Eps
