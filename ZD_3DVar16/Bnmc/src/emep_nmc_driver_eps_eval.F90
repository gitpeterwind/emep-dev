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

module EMEP_NMC_Driver_Eps_Eval

  use GO                  , only : gol, goPr, goErr
  use EMEP_NMC_Common     , only : T_ModelVars
  use EMEP_NMC_Common     , only : T_Points
  use GO                  , only : TDate

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_Eps_Eval
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_Eps_Eval'
  
  ! --- types ----------------------------------------
  
  type T_Driver_Eps_Eval
    ! time range:
    type(TDate)                       ::  tday1, tday2
    ! variables:
    type(T_ModelVars)                 ::  mvars
    ! points:
    type(T_Points)                    ::  points
    ! flags:
    logical                           ::  correlation
    logical                           ::  levtr
    ! statistics files created or read:
    character(len=1024)               ::  eps_samples_filename
    character(len=1024)               ::  B_point_filename
  contains
    procedure   ::  Init            => Driver_Eps_Eval_Init
    procedure   ::  Done            => Driver_Eps_Eval_Done
    procedure   ::  Evaluate        => Driver_Eps_Eval_Evaluate
  end type T_Driver_Eps_Eval
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_Eps_Eval_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine, goSplitString
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eps_Eval), intent(out)           ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eps_Eval_Init'
    
    ! --- local ----------------------------------
    
    character(len=32)     ::  tvalue
    
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
    
    ! variable names:
    call self%mvars%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! point locations:
    call self%points%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! input/output:
    call ReadRc( rcF, 'nmc.eps.samples.filename', self%eps_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! flag:
    call ReadRc( rcF, 'nmc.sample-B.correlation', self%correlation, status )
    IF_NOT_OK_RETURN(status=1)
    ! flag:
    call ReadRc( rcF, 'nmc.sample-B.levtr', self%levtr, status )
    IF_NOT_OK_RETURN(status=1)

    ! target file:
    call ReadRc( rcF, 'nmc.sample-B.filename', self%B_point_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_Eps_Eval_Init


  ! ***
  

  subroutine Driver_Eps_Eval_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eps_Eval), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eps_Eval_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%mvars%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    call self%points%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_Eps_Eval_Done


  ! ***
  
  
  !
  ! Evaluate sample covariance matrix at selected locations.
  !

  subroutine Driver_Eps_Eval_Evaluate( self, status )
  
    use GO             , only : TDate, IncrDate, operator(+), operator(>), wrtgol
    use GO             , only : TrcFile, ReadRc
    use GO             , only : goReplace
    use GO             , only : T_SampStat, T_SampCovr_3D_1D, T_SampCovr_2D_2D
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eps_Eval), intent(inout)            ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eps_Eval_Evaluate'
  
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
    
    type(T_NMC_Output)               ::  out_file
    character(len=1024)              ::  msg
    integer                          ::  ps__varid
    integer, allocatable             ::  S__varids(:)  ! (ntracer)
    real, allocatable                ::  S(:,:,:)         ! (nlon,nlat,nlev)
    integer, allocatable             ::  Bx__varids(:,:)  ! (npoint,ntracer)
    real, allocatable                ::  Bx(:,:,:)        ! (nlon,nlat,nlev)
    real                             ::  stdv

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
    
    ! set label:
    if ( self%correlation ) then
      cv_label = 'correlation'
      cv_id    = 'C'
      cv_units = '1'
    else
      cv_label = 'covariance'
      cv_id    = 'B'
      cv_units = 'tracer*tracer'
    end if
    
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
      call goReplace( filename, '%{hh}'  , '(i2.2)', tday%hour , status )
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
        allocate( eps__varids(self%mvars%n) )
        ! loop over variables:
        do ivar = 1, self%mvars%n
          ! get id:
          call eps_file%Get_VarID( trim(self%mvars%value(ivar)%name), eps__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do

        ! storage for 3D fields:
        allocate( ps(nlon,nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( eps(nlon,nlat,nlev,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        
        ! storage for meta data:
        allocate( eps__units(self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ points

        ! loop over points:
        do ipoint = 1, self%points%n
          ! nearby index:
          call lon_coor%Get_Index( self%points%value(ipoint)%lon, self%points%value(ipoint)%ilon, status )
          IF_NOT_OK_RETURN(status=1)
          ! nearby index:
          call lat_coor%Get_Index( self%points%value(ipoint)%lat, self%points%value(ipoint)%ilat, status )
          IF_NOT_OK_RETURN(status=1)
          ! surface (top-down order!)
          self%points%value(ipoint)%ilev = nlev
        end do

        ! storage for eps profiles:
        allocate( eps_profiles(self%points%n,nlev,self%mvars%n), stat=status )
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
        allocate( sampstat_eps(ntime,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over hours:
        do itime = 1, ntime
          ! loop over variables:
          do ivar = 1, self%mvars%n
            ! init mean/sigma computation over 3D field:
            call sampstat_eps(itime,ivar)%Init( (/nlon,nlat,nlev/), status )
            IF_NOT_OK_RETURN(status=1)
          end do  ! var
        end do  ! hour

        ! covariance of 3D fields with point values:
        allocate( sampcovr_eps(ntime,self%mvars%n,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over hours:
        do itime = 1, ntime
          ! loop over variables:
          do ivar = 1, self%mvars%n
            do ivar2 = 1, self%mvars%n
              ! init covariance between 3D field and point values:
              call sampcovr_eps(itime,ivar,ivar2)%Init( (/nlon,nlat,nlev/), self%points%n, status )
              IF_NOT_OK_RETURN(status=1)
            end do
          end do  ! var
        end do  ! hour

        ! apply?
        if ( self%levtr ) then
          ! covariance between tracers/levels :
          allocate( sampcovr_Bzs(ntime,self%points%n), stat=status )
          IF_NOT_OK_RETURN(status=1)
          ! loop over hours:
          do itime = 1, ntime
            ! loop over points:
            do ipoint = 1, self%points%n
              ! init covariance between pair of 2D fields
              call sampcovr_Bzs(itime,ipoint)%Init( (/nlev,self%mvars%n/), (/nlev,self%mvars%n/), status )
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
        do ivar = 1, self%mvars%n
    
          ! info ...
          write (gol,'("    variable ",i0," ",a)') itime, trim(self%mvars%value(ivar)%name); call goPr
        
          ! read sample:
          call eps_file%Get_Sample3D( trim(self%mvars%value(ivar)%name), itime, eps_isample, &
                                        eps(:,:,:,ivar), eps__units(ivar), status )
          IF_NOT_OK_RETURN(status=1)
          
          ! update statistics:
          call sampstat_eps(itime,ivar)%AddSample( eps(:,:,:,ivar), status )
          IF_NOT_OK_RETURN(status=1)
          
        end do ! var

        ! collect sample differences at point locations:
        do ipoint = 1, self%points%n

          ! grid cell indices of point:
          ilon = self%points%value(ipoint)%ilon
          ilat = self%points%value(ipoint)%ilat
          ilev = self%points%value(ipoint)%ilev

          ! profiles:
          eps_profiles(ipoint,:,:) = eps(ilon,ilat,:,:)
          
        end do ! points

        ! loop over variables:
        do ivar = 1, self%mvars%n
          do ivar2 = 1, self%mvars%n

            ! update sample covariances:
            call sampcovr_eps(itime,ivar,ivar2)%AddSample( eps(:,:,:,ivar), eps_profiles(:,ilev,ivar2), status )
            IF_NOT_OK_RETURN(status=1)

          end do ! var2
        end do ! var

        ! level/tracer covar?
        if ( self%levtr ) then
          ! loop over points:
          do ipoint = 1, self%points%n
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
    
    ! create new result file:
    call out_file%Create( trim(self%B_point_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("evaluate ",a," from `",a,"`")') &
            trim(cv_label), trim(self%eps_samples_filename)
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
    call tracer_coor%Set_Dim( status, n=self%mvars%n )
    IF_NOT_OK_RETURN(status=1)
    ! set attributes:
    call tracer_coor%Set_Attrs( status, long_name='tracer' )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do ivar = 1, self%mvars%n
      ! fill tracer name:
      call tracer_coor%Set_Value( ivar, status, value=trim(self%mvars%value(ivar)%name) )
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
    allocate( S__varids(self%mvars%n), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bx__varids(self%points%n,self%mvars%n), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over tracers:
    do ivar = 1, self%mvars%n

      ! copy:
      tracer_name  = trim(self%mvars%value(ivar)%name)
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
      do ipoint = 1, self%points%n

        ! variable:
        write (vname,'(a,"_",a,"_",a)') trim(cv_id), trim(tracer_name), trim(self%points%value(ipoint)%name)
        ! define:
        call out_file%Def_Sample3D( trim(vname), trim(cv_units), &
                                lon_coor, lat_coor, lev_coor, ctime_coor, tracer_coor, &
                                Bx__varids(ipoint,ivar), status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'long_name', &
                       trim(cv_label)//' with '//trim(tracer_name)//' in '//trim(self%points%value(ipoint)%name), status )
        IF_NOT_OK_RETURN(status=1)
        ! write as attributes:
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'lon', self%points%value(ipoint)%lon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'lat', self%points%value(ipoint)%lat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'ilon', self%points%value(ipoint)%ilon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'ilat', self%points%value(ipoint)%ilat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,ivar), 'ilev', self%points%value(ipoint)%ilev, status )
        IF_NOT_OK_RETURN(status=1)

      end do ! points

    end do ! tracers
    
    ! evaluate level/tracer covariance?
    if ( self%levtr ) then

      ! storage for result:
      allocate( Bzs(nlev,nlev,self%mvars%n,self%mvars%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! idem in different order:
      allocate( Vzs(nlev,self%mvars%n,nlev,self%mvars%n), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! output:
      allocate( Bzs__varids(self%points%n), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over points:
      do ipoint = 1, self%points%n

        ! variable:
        write (vname,'(a,"zs_",a)') trim(cv_id), trim(self%points%value(ipoint)%name)
        ! define:
        call out_file%Def_ACovar1( trim(vname), trim(cv_units), &
                                    lev_coor, tracer_coor, ctime_coor, &
                                    Bzs__varids(ipoint), status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'long_name', &
                       'vertical '//trim(cv_label)//' in '//trim(self%points%value(ipoint)%name), status )
        IF_NOT_OK_RETURN(status=1)
        ! write as attributes:
        call out_file%Put_Att( Bzs__varids(ipoint), 'lon', self%points%value(ipoint)%lon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'lat', self%points%value(ipoint)%lat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'ilon', self%points%value(ipoint)%ilon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'ilat', self%points%value(ipoint)%ilat, status )
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
      do ivar = 1, self%mvars%n
      
        ! extract 3D std.dev. field:
        call sampstat_eps(itime,ivar)%Get3D( status, stdv=S )
        IF_NOT_OK_RETURN(status=1)
        ! write:
        call out_file%Put_Field3d_Series( S__varids(ivar), itime, S, status )
        IF_NOT_OK_RETURN(status=1)  
          
      end do  ! var

      ! ~ 3D covariance with point

      ! loop over variables:
      do ivar = 1, self%mvars%n

        ! second tracer dim:
        do ivar2 = 1, self%mvars%n

          ! loop over points:
          do ipoint = 1, self%points%n
        
            ! correlation or covariance?
            if ( self%correlation ) then
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
      if ( self%levtr ) then 

        ! loop over points:
        do ipoint = 1, self%points%n
        
          ! correlation or covariance?
          if ( self%correlation ) then
            ! extract sample result:
            call sampcovr_Bzs(itime,ipoint)%GetXY( status, corr=Vzs )
            IF_NOT_OK_RETURN(status=1)
          else
            ! extract sample result:
            call sampcovr_Bzs(itime,ipoint)%GetXY( status, covr=Vzs )
            IF_NOT_OK_RETURN(status=1)
          end if
          ! re-order:
          do ivar = 1, self%mvars%n
            do ivar2 = 1, self%mvars%n
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
    if ( self%levtr ) then
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
      do ivar = 1, self%mvars%n
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
      do ivar = 1, self%mvars%n
        do ivar2 = 1, self%mvars%n
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
  
  end subroutine Driver_Eps_Eval_Evaluate


end module EMEP_NMC_Driver_Eps_Eval

