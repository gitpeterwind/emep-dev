!######################################################################
!
! Read OMI/TROPOMI file created by EMIP (EMEP Input Processor).
!
! EXAMPLE CONTENT
!
!   netcdf OMI-Aura_NO2_20120701_0100 {
!     dimensions:
!	      track_image = 737 ;
!	      track_pixel = 601 ;
!	      corner = 4 ;
!	      pixel = 4435 ;
!	      layer = 34 ;
!	      layer_interface = 35 ;
!     variables:
!	      float track_longitude(track_image, track_pixel) ;
!	        track_longitude:units = "degrees_east" ;
!	        track_longitude:long_name = "longitudes" ;
!	      float track_latitude(track_image, track_pixel) ;
!	        track_latitude:units = "degrees_north" ;
!	        track_latitude:long_name = "latitudes" ;
!	      float track_corner_longitudes(track_image, track_pixel, corner) ;
!	        track_corner_longitudes:units = "degrees_east" ;
!	        track_corner_longitudes:long_name = "corner longitudes" ;
!	      float track_corner_latitudes(track_image, track_pixel, corner) ;
!	        track_corner_latitudes:units = "degrees_north" ;
!	        track_corner_latitudes:long_name = "corner latitudes" ;
!	      float longitude(pixel) ;
!	        longitude:units = "degrees_east" ;
!	        longitude:long_name = "longitudes" ;
!	      float latitude(pixel) ;
!	        latitude:units = "degrees_north" ;
!	        latitude:long_name = "latitudes" ;
!	      float corner_longitudes(pixel, corner) ;
!	        corner_longitudes:units = "degrees_east" ;
!	        corner_longitudes:long_name = "corner longitudes" ;
!	      float corner_latitudes(pixel, corner) ;
!	        corner_latitudes:units = "degrees_north" ;
!	        corner_latitudes:long_name = "corner latitudes" ;
!	      int pixel(pixel) ;
!	        pixel:units = "1" ;
!	        pixel:long_name = "compressed coordinate" ;
!	        pixel:compress = "track_image track_pixel" ;
!	        pixel:description = "original zero-based indices: `track_image` = `pixel`/len(`track_pixel`), `track_pixel` = `pixel` mod len(`track_pixel`)" ;
!	      int layer(layer) ;
!	        layer:units = "1" ;
!	        layer:standard_name = "atmosphere_hybrid_sigma_pressure_coordinate" ;
!	        layer:formula = "p(n,k,i) = ap(k) + b(k)*ps(n,i)" ;
!	        layer:formula_terms = "ap: hyam b: hybm ps: surface_pressure" ;
!	      double time(pixel) ;
!	        time:units = "seconds since 1993-01-01 00:00:00" ;
!	        time:long_name = "Time at Start of Scan (s, TAI93)" ;
!	      float vcd_trop(pixel) ;
!	        vcd_trop:_FillValue = 1.e+20f ;
!	        vcd_trop:units = "1e15 cm**-2" ;
!	        vcd_trop:long_name = "NO2 troposheric vertical column density" ;
!	      float sigma_vcd_trop(pixel) ;
!	        sigma_vcd_trop:_FillValue = 1.e+20f ;
!	        sigma_vcd_trop:units = "1e15 cm**-2" ;
!	        sigma_vcd_trop:long_name = "Error in the NO2 tropospheric vertical column density" ;
!	      float sigma_vcd_trop_ak(pixel) ;
!	        sigma_vcd_trop_ak:_FillValue = 1.e+20f ;
!	        sigma_vcd_trop_ak:units = "1e15 cm**-2" ;
!	        sigma_vcd_trop_ak:long_name = "Error in the NO2 tropospheric vertical column density using averaging kernel information (w/o profile error contribution)" ;
!	      float kernel(pixel, layer) ;
!	        kernel:units = "1" ;
!	        kernel:long_name = "averaging kernel" ;
!	      float pressure_levels(pixel, layer_interface) ;
!	        pressure_levels:units = "Pa" ;
!	        pressure_levels:long_name = "half level pressures" ;
!	      float cloud_top_pressure(pixel) ;
!	        cloud_top_pressure:_FillValue = 1.e+20f ;
!	        cloud_top_pressure:units = "hPa" ;
!	        cloud_top_pressure:long_name = "Effective cloud pressure" ;
!	      float cloud_fraction(pixel) ;
!	        cloud_fraction:_FillValue = 1.e+20f ;
!	        cloud_fraction:units = "1" ;
!	        cloud_fraction:long_name = "Effective cloud fraction" ;
!	      float cloud_radiance_fraction(pixel) ;
!	        cloud_radiance_fraction:_FillValue = 1.e+20f ;
!	        cloud_radiance_fraction:units = "%" ;
!	        cloud_radiance_fraction:long_name = "Fraction of the radiance coming from the cloudy part of the pixel (%)" ;
!	      float image_number(pixel) ;
!	        image_number:_FillValue = -999.9f ;
!	        image_number:units = "1" ;
!	        image_number:long_name = "scan number" ;
!	      float pixel_number(pixel) ;
!	        pixel_number:_FillValue = -999.9f ;
!	        pixel_number:units = "1" ;
!	        pixel_number:long_name = "pixel number" ;
!	      float tropopause_layer(pixel) ;
!	        tropopause_layer:_FillValue = 1.e+20f ;
!	        tropopause_layer:units = "1" ;
!	        tropopause_layer:long_name = "TM4 pressure level where the tropopause occurs" ;
!	      float surface_pressure(pixel) ;
!	        surface_pressure:_FillValue = 1.e+20f ;
!	        surface_pressure:units = "Pa" ;
!	        surface_pressure:long_name = "Model surface pressure of the center of the ground pixel, as used in the AMF calculation (Zhou et al., AMT, 2, 401-416, 2009)" ;
!	      float hyam(layer) ;
!	        hyam:_FillValue = 1.e+20f ;
!	        hyam:units = "Pa" ;
!	        hyam:long_name = "TM4 pressure level p = a + p_surf*b" ;
!	      float hybm(layer) ;
!	        hybm:_FillValue = 1.e+20f ;
!	        hybm:units = "1" ;
!	        hybm:long_name = "TM4 pressure level p = a + p_surf*b" ;
!	      float hyai(layer_interface) ;
!	        hyai:units = "Pa" ;
!	        hyai:long_name = "half level hybride coefficients" ;
!	      float hybi(layer_interface) ;
!	        hybi:units = "1" ;
!	        hybi:long_name = "half level hybride coefficients" ;
!	      float orbit_number(pixel) ;
!	        orbit_number:_FillValue = -999.9f ;
!	        orbit_number:units = "1" ;
!	        orbit_number:long_name = "orbit number" ;
!
!      // global attributes:
!	        :format = "1.0" ;
!	        :Conventions = "CF-1.6" ;
!	        :author = "Arjo Segers" ;
!	        :institution = "MetNorway, Oslo, Norway" ;
!	        :email = "Arjo.Segers@met.no" ;
!	        :history = "Tue Nov 13 15:15:18 2018, averaged OMI-Aura_NO2_20120701_o42344.nc to regular grid\n",
!	    	    "Mon Nov 05 17:53:48 2018, selected pixels with \'Longitude\' in [-30.0625,45.0625], selected pixels with \'Latitude\' in [29.96875,76.03125], selected pixels with valid \'TroposphericVerticalColumn\' values, selected pixels with \'TroposphericColumnFlag\' == 0, selected pixels with \'SurfaceAlbedo\' <= 0.3, added 253 pixels from OMI-Aura_L2-OMDOMINO_2012m0701t0017-o42344_v003-2012m0705t000410.he5" ;
!      }
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

module EMIP

  use GO, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_EMIP
  

  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'EMIP'


  ! --- types ----------------------------------------
  
  type T_EMIP
    ! originating file:
    character(len=1024)   ::  filename
    ! dimensions:
    integer               ::  npixel
    integer               ::  nlev
    integer               ::  nhlev
    ! data values:
    real, allocatable     ::  longitude(:)
    character(len=32)     ::  longitude_units
    real, allocatable     ::  latitude(:)
    character(len=32)     ::  latitude_units
    real, allocatable     ::  vcd(:)
    character(len=32)     ::  vcd_units
    real, allocatable     ::  sigma_vcd(:)
    character(len=32)     ::  sigma_vcd_units
    real, allocatable     ::  kernel(:,:)
    character(len=32)     ::  kernel_units
    real, allocatable     ::  phlev(:,:)
    character(len=32)     ::  phlev_units
  contains
    procedure   ::  Init          => EMIP_Init
    procedure   ::  Done          => EMIP_Done
  end type T_EMIP


contains



  !
  ! Access EMEP OMI extract
  !
  
  subroutine EMIP_Init( self, filename, status )
  
    use C3PO      , only : Datafile

    ! --- in/out ---------------------------------
    
    class(T_EMIP), intent(out)                ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/EMIP_Init'
    
    ! --- local ----------------------------------
    
    ! data file:
    type(Datafile)            ::  df

    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("read EMEP OMI file ...")'); call goPr
    write (gol,'("  input file: ",a)') trim(filename); call goPr
    
    ! store:
    self%filename = filename
    
    ! open file:
    call df%Open( trim(self%filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! get dimensions:
    call df%Get_Dim( 'pixel', status, n=self%npixel )
    IF_NOT_OK_RETURN(status=1)
    call df%Get_Dim( 'layer', status, n=self%nlev )
    IF_NOT_OK_RETURN(status=1)
    call df%Get_Dim( 'layer_interface', status, n=self%nhlev )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage:
    allocate( self%longitude(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%latitude(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%vcd(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%sigma_vcd(self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%kernel(self%nlev,self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%phlev(self%nhlev,self%npixel), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read locations:
    call df%Get_Var( 'longitude', self%longitude, self%longitude_units, status )
    IF_NOT_OK_RETURN(status=1)
    call df%Get_Var( 'latitude', self%latitude, self%latitude_units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! different names are used ...
    call df%Inquire( status, varname='vcd_trop' )
    if ( status == 0 ) then
      call df%Get_Var( 'vcd_trop', self%vcd, self%vcd_units, status )
      IF_NOT_OK_RETURN(status=1)
      call df%Get_Var( 'sigma_vcd_trop', self%sigma_vcd, self%sigma_vcd_units, status )
      IF_NOT_OK_RETURN(status=1)
    else if ( status < 0 ) then
      call df%Get_Var( 'vcd', self%vcd, self%vcd_units, status )
      IF_NOT_OK_RETURN(status=1)
      call df%Get_Var( 'sigma_vcd', self%sigma_vcd, self%sigma_vcd_units, status )
      IF_NOT_OK_RETURN(status=1)
    else
      write (gol,'("could not read (find?) vcd variables")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! kernel and levels:
    call df%Get_Var( 'kernel', self%kernel, self%kernel_units, status )
    IF_NOT_OK_RETURN(status=1)
    call df%Get_Var( 'pressure_levels', self%phlev, self%phlev_units, status )
    IF_NOT_OK_RETURN(status=1)

    ! close:
    call df%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine EMIP_Init


  ! ***
  

  subroutine EMIP_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_EMIP), intent(inout)          ::  self
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/EMIP_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%longitude, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%latitude, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%vcd, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%sigma_vcd, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%kernel, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%phlev, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine EMIP_Done



end module EMIP


!!=======================================================================
!!===
!!=== test
!!===
!!=======================================================================
!
!program test
!
!  use EMIP, only : T_EMIP
!  
!
!  type(T_EMIP)      ::  omi
!  character(len=1024)   ::  filename
!  integer               ::  status
!  
!  ! sample:
!  filename = '/home/metno/alvarov/CWF/obs/2017_OMI/OMI_201701.nc'
!  
!end program test

