!######################################################################
!
! Read S5p superobservation file provided by Henk Eskes (KNMI)
!
! EXAMPLE CONTENT
!
!    netcdf s5p_no2_offl_superobs_0.5deg_emep_global_20190101_test_d {
!    dimensions:
!	    lon = 720 ;
!	    lat = 360 ;
!	    layer = 34 ;
!	    vertices = 2 ;
!    variables:
!	    float lon(lon) ;
!		    lon:standard_name = "longitude" ;
!		    lon:long_name = "longitude" ;
!		    lon:units = "degrees_east" ;
!	    float lat(lat) ;
!		    lat:standard_name = "latitude" ;
!		    lat:long_name = "latitude" ;
!		    lat:units = "degrees_north" ;
!	    float layer(layer) ;
!		    layer:standard_name = "layer" ;
!		    layer:long_name = "layer" ;
!		    layer:units = "1" ;
!	    float no2_superobs(lat, lon) ;
!		    no2_superobs:standard_name = "nitrogen_dioxide_tropospheric_column_superobservation" ;
!		    no2_superobs:long_name = "NO2 tropospheric column superobservation, TROPOMI sensor" ;
!		    no2_superobs:units = "1e-6 mol m^-2" ;
!		    no2_superobs:_FillValue = 9.96921e+36f ;
!	    float no2_superobs_uncertainty(lat, lon) ;
!		    no2_superobs_uncertainty:standard_name = "nitrogen_dioxide_tropospheric_column_superobservation_uncertainty" ;
!		    no2_superobs_uncertainty:long_name = "NO2 tropospheric column superobservation uncertainty, TROPOMI sensor, including representation error" ;
!		    no2_superobs_uncertainty:units = "1e-6 mol m^-2" ;
!		    no2_superobs_uncertainty:_FillValue = 9.96921e+36f ;
!	    float no2_superobs_uncertainty_representation(lat, lon) ;
!		    no2_superobs_uncertainty_representation:standard_name = "nitrogen_dioxide_tropospheric_column_superobservation_representativity_uncertainty" ;
!		    no2_superobs_uncertainty_representation:long_name = "NO2 tropospheric column superobservation uncertainty, TROPOMI sensor, only representation error" ;
!		    no2_superobs_uncertainty_representation:units = "1e-6 mol m^-2" ;
!		    no2_superobs_uncertainty_representation:_FillValue = 9.96921e+36f ;
!	    float average_uncertainty(lat, lon) ;
!		    average_uncertainty:standard_name = "nitrogen_dioxide_tropospheric_column_average_uncertainty" ;
!		    average_uncertainty:long_name = "NO2 tropospheric column superobservation uncertainty, TROPOMI sensor, average over individual observation" ;
!		    average_uncertainty:units = "1e-6 mol m^-2" ;
!		    average_uncertainty:_FillValue = 9.96921e+36f ;
!	    float surface_pressure(lat, lon) ;
!		    surface_pressure:standard_name = "surface_pressure" ;
!		    surface_pressure:long_name = "superobservation surface pressure" ;
!		    surface_pressure:units = "Pa" ;
!		    surface_pressure:_FillValue = 9.96921e+36f ;
!	    float cloud_radiance_fraction(lat, lon) ;
!		    cloud_radiance_fraction:standard_name = "cloud_radiance_fraction" ;
!		    cloud_radiance_fraction:long_name = "superobservation cloud radiance fraction" ;
!		    cloud_radiance_fraction:units = "1" ;
!		    cloud_radiance_fraction:_FillValue = 9.96921e+36f ;
!	    float kernel_troposphere(layer, lat, lon) ;
!		    kernel_troposphere:standard_name = "averaging_kernel_troposphere" ;
!		    kernel_troposphere:long_name = "NO2 tropospheric column averaging kernel of the superobservation, TROPOMI sensor" ;
!		    kernel_troposphere:units = "1" ;
!		    kernel_troposphere:_FillValue = 9.96921e+36f ;
!	    float covered_area_fraction(lat, lon) ;
!		    covered_area_fraction:standard_name = "covered_area_fraction" ;
!		    covered_area_fraction:long_name = "Fraction of the grid box covered by satellite observations" ;
!		    covered_area_fraction:units = "1" ;
!		    covered_area_fraction:_FillValue = 9.96921e+36f ;
!	    float no2_scd_geo_superobs(lat, lon) ;
!		    no2_scd_geo_superobs:standard_name = "nitrogen_dioxide_slant_amfgeo_superobservation" ;
!		    no2_scd_geo_superobs:long_name = "NO2 slant column divided by AMFgeo superobservation, TROPOMI sensor" ;
!		    no2_scd_geo_superobs:units = "1e-6 mol m^-2" ;
!		    no2_scd_geo_superobs:_FillValue = 9.96921e+36f ;
!	    float no2_strat_superobs(lat, lon) ;
!		    no2_strat_superobs:standard_name = "nitrogen_dioxide_stratospheric_column_superobservation" ;
!		    no2_strat_superobs:long_name = "NO2 stratospheric column superobservation, TROPOMI sensor" ;
!		    no2_strat_superobs:units = "1e-6 mol m^-2" ;
!		    no2_strat_superobs:_FillValue = 9.96921e+36f ;
!	    float amf_trop_superobs(lat, lon) ;
!		    amf_trop_superobs:standard_name = "tropospheric_air_mass_factor_superobservation" ;
!		    amf_trop_superobs:long_name = "NO2 tropospheric air mass factor, TROPOMI sensor" ;
!		    amf_trop_superobs:units = "1" ;
!		    amf_trop_superobs:_FillValue = 9.96921e+36f ;
!	    float amf_strat_superobs(lat, lon) ;
!		    amf_strat_superobs:standard_name = "stratospheric_air_mass_factor_superobservation" ;
!		    amf_strat_superobs:long_name = "NO2 stratospheric air mass factor, TROPOMI sensor" ;
!		    amf_strat_superobs:units = "1" ;
!		    amf_strat_superobs:_FillValue = 9.96921e+36f ;
!	    float tm5_constant_a(layer, vertices) ;
!		    tm5_constant_a:standard_name = "tm5_constant_a" ;
!		    tm5_constant_a:long_name = "TM5 hybrid A coefficient at upper and lower interface levels" ;
!		    tm5_constant_a:units = "Pa" ;
!		    tm5_constant_a:_FillValue = 9.96921e+36f ;
!	    float tm5_constant_b(layer, vertices) ;
!		    tm5_constant_b:standard_name = "tm5_constant_b" ;
!		    tm5_constant_b:long_name = "TM5 hybrid B coefficient at upper and lower interface levels" ;
!		    tm5_constant_b:units = "1" ;
!		    tm5_constant_b:_FillValue = 9.96921e+36f ;
!
!    // global attributes:
!		    :authors = "Henk Eskes, eskes@knmi.nl" ;
!		    :institution = "Royal Netherlands Meteorological Institute (KNMI)" ;
!		    :product_version = "1.2.2, reprocessing, RPRO" ;
!		    :project = "S5P, phase E2" ;
!		    :reference = "http://www.tropomi.eu/" ;
!		    :source = "TROPOMI / Sentinel-5P" ;
!		    :title = "Nitrogen dioxide (NO2) column data, gridded fields, 1 Jan 2019, India overpass" ;
!		    :vertical_column_processor_version = "TM5-MP-Domino3 v3.5.3" ;
!		    :date_created = "2020-02-19T09:12:47Z" ;
!		    :filter_min_QA_value = 0.75f ;
!    }
!
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
!######################################################################

module HESK_SuperObs

  use GO, only : gol, goPr, goErr
  use GO, only : TDate

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_HESK_Listing
  public  ::  T_HESK_SuperObs
  

  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'HESK_SuperObs'


  ! --- types ----------------------------------------
  
  type T_HESK_Listing
    ! originating file:
    character(len=1024)                 ::  dirname
    character(len=1024)                 ::  filename
    ! number of records:
    integer                             ::  nrec
    ! data values:
    character(len=128), allocatable     ::  orbitfile(:)  ! (nrec)
    character(len=8), allocatable       ::  orbit(:)      ! (nrec)
    type(TDate), allocatable            ::  tt(:,:)       ! (2,nrec)
    !
  contains
    procedure   ::  Init          => HESK_Listing_Init
    procedure   ::  Done          => HESK_Listing_Done
    procedure   ::  Search        => HESK_Listing_Search
  end type T_HESK_Listing
  
  ! *
  
  type T_HESK_SuperObs
    ! originating file:
    character(len=1024)   ::  filename
    ! dimensions:
    integer               ::  npixel
    integer               ::  nlon, nlat
    integer               ::  nlev
    integer               ::  nhlev
    ! data values:
    integer, allocatable  ::  ii(:)
    integer, allocatable  ::  jj(:)
    real, allocatable     ::  longitude(:)
    character(len=32)     ::  longitude_units
    real, allocatable     ::  latitude(:)
    character(len=32)     ::  latitude_units
    real, allocatable     ::  vcd_trop(:)
    character(len=32)     ::  vcd_trop_units
    real, allocatable     ::  sigma_vcd_trop(:)
    character(len=32)     ::  sigma_vcd_trop_units
    real, allocatable     ::  kernel(:,:)
    character(len=32)     ::  kernel_units
    real, allocatable     ::  phlev(:,:)
    character(len=32)     ::  phlev_units
  contains
    procedure   ::  Init          => HESK_SuperObs_Init
    procedure   ::  Done          => HESK_SuperObs_Done
  end type T_HESK_SuperObs


contains



  ! ====================================================================
  ! ===
  ! === Listing file
  ! ===
  ! ====================================================================
  
  !
  ! Listing file is a csv table with filenames and timerange of pixels included.
  ! Sample content:
  !
  !     orbit;filename;time_coverage_start;time_coverage_end
  !     2833;201805_0.5deg/compressed/s5p_no2_rpro_v1.2.2_superobs_ll_0.5deg_02833.nc;2018-05-01T02:03:58Z;2018-05-01T02:56:14Z
  !     2834;201805_0.5deg/compressed/s5p_no2_rpro_v1.2.2_superobs_ll_0.5deg_02834.nc;2018-05-01T03:45:28Z;2018-05-01T04:43:52Z
  !        :
  !

  subroutine HESK_Listing_Init( self, filename, status )
  
    use GO, only : goGetFU
    use GO, only : goDirname
    use GO, only : goSplitString
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_HESK_Listing), intent(out)        ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/HESK_Listing_Init'
    
    ! --- local ----------------------------------
    
    logical                 ::  exist
    integer                 ::  fu
    integer                 ::  nloop, iloop
    integer                 ::  iline
    character(len=1)        ::  sep
    character(len=1024)     ::  line
    integer                 ::  irec
    integer                 ::  nheader, iheader
    character(len=64)       ::  headers(10)

    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("read HESK listing file: ",a)') trim(filename); call goPr
    
    ! check ...
    inquire( file=filename, exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found: ",a)') trim(filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! store:
    self%filename = filename
    
    ! directory part:
    call goDirname( self%filename, self%dirname, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! seperation character in file:
    sep = ';'
    
    ! new file unit:
    call goGetFU( fu, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of loops needed, first 2 (number of records and values), then 1
    nloop = 2
    ! loops:
    do iloop = 1, nloop

      ! open for reading:
      open( fu, file=trim(self%filename), iostat=status )
      IF_NOT_OK_RETURN(status=1)

      ! init counters:
      iline = 0
      irec = 0
      ! loop over lines:
      do
        ! read line:
        read (fu,'(a)',iostat=status) line
        if ( status < 0 ) exit ! eof
        if ( status > 0 ) then
          write (gol,'("reading line ",i0," from ",a)') iline, trim(self%filename); call goErr
          TRACEBACK; status=1; return
        end if

        ! increase counter:
        iline = iline + 1
        ! header line?
        if ( iline == 1 ) then
          ! split:
          call goSplitString( line, nheader, headers, status, sep=sep )
          IF_NOT_OK_RETURN(status=1)
          ! next:
          cycle
        end if
        
        ! increase counter:
        irec = irec + 1
        ! counting only?
        if ( iloop == 1 ) cycle

        ! read values:
        do iheader = 1, nheader
          ! switch:
          select case ( trim(headers(iheader)) )
            !~
            case ( 'orbit' )
              call goReadFromLine( line, self%orbit(irec), status, sep=sep )
              IF_NOT_OK_RETURN(status=1)
            !~
            case ( 'filename' )
              call goReadFromLine( line, self%orbitfile(irec), status, sep=sep )
              IF_NOT_OK_RETURN(status=1)
            !~
            case ( 'time_coverage_start' )
              call goReadFromLine( line, self%tt(1,irec), status, sep=sep )
              IF_NOT_OK_RETURN(status=1)
            !~
            case ( 'time_coverage_end' )
              call goReadFromLine( line, self%tt(2,irec), status, sep=sep )
              IF_NOT_OK_RETURN(status=1)
            !~
            case default
              write (gol,'("unsupported header `",a,"` in file: ",a)') &
                       trim(headers(iheader)), trim(self%filename); call goErr
              TRACEBACK; status=1; return
          end select
        end do ! headers
          
      end do ! lines
      ! close:
      close( fu, iostat=status )
      IF_NOT_OK_RETURN(status=1)
      
      ! storage?
      if ( iloop == 1 ) then
        ! number of records:
        self%nrec = irec
        ! storage:
        allocate( self%orbitfile(self%nrec), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%orbit(self%nrec), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%tt(2,self%nrec), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
    end do ! iloop
    
    ! ok
    status = 0
  
  end subroutine HESK_Listing_Init


  ! ***
  

  subroutine HESK_Listing_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_HESK_Listing), intent(inout)      ::  self
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/HESK_Listing_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%orbitfile, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%orbit, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%tt, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine HESK_Listing_Done


  ! ***
  

  subroutine HESK_Listing_Search( self, cdate, orbitfile, orbit, status )

    use GO            , only : TDate, NewDate, IncrDate, operator(-)
    use GO            , only : operator(<), operator(<=), wrtgol
    use GO            , only : pathsep
    use TimeDate_mod  , only : date
  
    ! --- in/out ---------------------------------
    
    class(T_HESK_Listing), intent(inout)      ::  self
    type(date), intent(in)                    ::  cdate
    character(len=*), intent(out)             ::  orbitfile
    character(len=*), intent(out)             ::  orbit
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/HESK_Listing_Search'
    
    ! --- local ----------------------------------
    
    type(TDate)     ::  t, t0
    integer         ::  irec
    integer         ::  nfound
    
    ! --- begin ----------------------------------
    
    ! convert ...
    t = NewDate( cdate%year, cdate%month, cdate%day, cdate%hour, int(floor(cdate%seconds/60.0)), modulo(cdate%seconds,60) )
    ! info ...
    call wrtgol( 'search oribit files including (or entirely in hour prior to) ', t ); call goPr
    
    ! previous hour:
    t0 = t - IncrDate(hour=1)
    
    ! init result:
    orbitfile = ''
    ! no matches yet ...
    nfound = 0
    ! loop over records:
    do irec = 1, self%nrec
      ! time matches?
      if ( ( (self%tt(1,irec) < t ) .and. (t <= self%tt(2,irec)) ) .or. &
           ( (t0 < self%tt(1,irec)) .and. (self%tt(2,irec) <= t) )     ) then
        ! increase counter:
        nfound = nfound + 1
        ! copy, include path:
        orbitfile = trim(self%dirname)//pathsep//trim(self%orbitfile(irec))
        ! copy:
        orbit = trim(self%orbit(irec))
        ! info ...
        write (gol,'("  found orbitfile : ",a)') trim(orbitfile); call goPr
        call wrtgol( '  valid for       : ', self%tt(:,irec) ); call goPr
      end if
    end do

    ! check ..
    if ( nfound > 1 ) then
      write (gol,'("found ",i0," matching orbit file, something wrong in selection ...")') nfound; call goPr
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
  
  end subroutine HESK_Listing_Search



  ! ====================================================================
  ! ===
  ! === Access nc file
  ! ===
  ! ====================================================================
  
  subroutine HESK_SuperObs_Init( self, filename, status )
  
    use C3PO      , only : Datafile

    ! --- in/out ---------------------------------
    
    class(T_HESK_SuperObs), intent(out)       ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/HESK_SuperObs_Init'
    
    ! --- local ----------------------------------
    
    ! data file:
    type(Datafile)            ::  df
    ! temporary field2d:
    integer                   ::  nlon, nlat
    real, allocatable         ::  field1d(:)
    real, allocatable         ::  field2d(:,:)    ! (nlon,nlat)
    real, allocatable         ::  field3d(:,:,:)  ! (nlon,nlat,nlev)
    character(len=32)         ::  units
    real                      ::  fill_value
    integer                   ::  ipix
    integer                   ::  i, j
    real, allocatable         ::  a(:,:), b(:,:)   ! (2,nlev)

    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("read HESK SuperObs file ...")'); call goPr
    write (gol,'("  input file: ",a)') trim(filename); call goPr
    
    ! store:
    self%filename = filename
    
    ! open file:
    call df%Open( trim(self%filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! get dimensions:
    call df%Get_Dim( 'lon', status, n=nlon )
    IF_NOT_OK_RETURN(status=1)
    call df%Get_Dim( 'lat', status, n=nlat )
    IF_NOT_OK_RETURN(status=1)
    ! number of levels:
    call df%Get_Dim( 'layer', status, n=self%nlev )
    IF_NOT_OK_RETURN(status=1)
    ! number of interfaces:
    self%nhlev = self%nlev + 1
    
    ! storage:
    allocate( field2d(nlon,nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( field3d(nlon,nlat,self%nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! read observations:
    call df%Get_Var( 'no2_superobs', field2d, units, status, fill_value=fill_value )
    IF_NOT_OK_RETURN(status=1)
    ! number of valid observations:
    self%npixel = count( field2d /= fill_value )
    ! check ...
    if ( self%npixel == 0 ) then
      write (gol,'("no valid values in ",a)') trim(self%filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage:
    allocate( self%ii(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%longitude(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%vcd_trop(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%sigma_vcd_trop(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%kernel(self%nlev,self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%phlev(self%nhlev,self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init counter:
    ipix = 0
    ! loop:
    do j = 1, nlat
      do i = 1, nlon
        ! valid?
        if ( field2d(i,j) /= fill_value ) then
          ! increase counter:
          ipix = ipix + 1
          ! store:
          self%ii(ipix) = i
          self%jj(ipix) = j
        end if ! valid
      end do ! i
    end do ! j

    ! storage:
    allocate( field1d(nlon), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call df%Get_Var( 'lon', field1d, self%longitude_units, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    do ipix = 1, self%npixel
      self%longitude(ipix) = field1d(self%ii(ipix))
    end do
    ! clear:
    deallocate( field1d, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( field1d(nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call df%Get_Var( 'lat', field1d, self%latitude_units, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    do ipix = 1, self%npixel
      self%latitude(ipix) = field1d(self%jj(ipix))
    end do
    ! clear:
    deallocate( field1d, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! read:
    call df%Get_Var( 'no2_superobs', field2d, self%vcd_trop_units, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    do ipix = 1, self%npixel
      self%vcd_trop(ipix) = field2d(self%ii(ipix),self%jj(ipix))
    end do

    ! read:
    call df%Get_Var( 'no2_superobs_uncertainty', field2d, self%sigma_vcd_trop_units, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    do ipix = 1, self%npixel
      self%sigma_vcd_trop(ipix) = field2d(self%ii(ipix),self%jj(ipix))
    end do

    ! read:
    call df%Get_Var( 'kernel_troposphere', field3d, self%kernel_units, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    do ipix = 1, self%npixel
      self%kernel(:,ipix) = field3d(self%ii(ipix),self%jj(ipix),:)
    end do

    ! read:
    call df%Get_Var( 'surface_pressure', field2d, self%phlev_units, status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for hybride coeff:
    allocate( a(2,self%nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( b(2,self%nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call df%Get_Var( 'tm5_constant_a', a, units, status )
    IF_NOT_OK_RETURN(status=1)
    call df%Get_Var( 'tm5_constant_b', b, units, status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    do ipix = 1, self%npixel
      self%phlev(1           ,ipix) = a(1,1) + b(1,1) * field2d(self%ii(ipix),self%jj(ipix))
      self%phlev(2:self%nhlev,ipix) = a(2,:) + b(2,:) * field2d(self%ii(ipix),self%jj(ipix))
    end do
    ! clear:
    deallocate( a, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( b, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( field2d, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( field3d, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! close:
    call df%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine HESK_SuperObs_Init


  ! ***
  

  subroutine HESK_SuperObs_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_HESK_SuperObs), intent(inout)     ::  self
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/HESK_SuperObs_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%ii, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%jj, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%longitude, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%latitude, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%vcd_trop, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%sigma_vcd_trop, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%kernel, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%phlev, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine HESK_SuperObs_Done



end module HESK_SuperObs


!!=======================================================================
!!===
!!=== test
!!===
!!=======================================================================
!
!program test
!
!  use EMIP_OMI, only : T_HESK_SuperObs
!  
!
!  type(T_HESK_SuperObs)      ::  omi
!  character(len=1024)   ::  filename
!  integer               ::  status
!  
!  ! sample:
!  filename = '/home/metno/alvarov/CWF/obs/2017_OMI/OMI_201701.nc'
!  
!end program test

