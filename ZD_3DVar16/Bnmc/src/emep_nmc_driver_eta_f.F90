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

module EMEP_NMC_Driver_Eta_f

  use GO                  , only : gol, goPr, goErr
  use EMEP_NMC_Common     , only : T_ModelVars
  use GO                  , only : TDate

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_Eta_f
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_Eta_f'
  
  ! --- types ----------------------------------------
  
  type T_Driver_Eta_f
    ! time range:
    type(TDate)                       ::  tday1, tday2
    ! variables:
    type(T_ModelVars)                 ::  mvars
    ! statistics files created or read:
    character(len=1024)               ::  eta_samples_filename
    character(len=1024)               ::  eta_f_samples_filename
  contains
    procedure   ::  Init            => Driver_Eta_f_Init
    procedure   ::  Done            => Driver_Eta_f_Done
    procedure   ::  Compute         => Driver_Eta_f_Compute
  end type T_Driver_Eta_f
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_Eta_f_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine, goSplitString
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eta_f), intent(out)          ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eta_f_Init'
    
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
    
    ! input/output:
    call ReadRc( rcF, 'nmc.eta.samples.filename', self%eta_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.eta_f.samples.filename', self%eta_f_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_Eta_f_Init


  ! ***
  

  subroutine Driver_Eta_f_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eta_f), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eta_f_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%mvars%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_Eta_f_Done


  ! ***
  
  
  !
  ! Compute Fourrier transforms of eta samples:
  ! 
  !    eta_f = F eta
  !

  subroutine Driver_Eta_f_Compute( self, rcF, status )
  
    use GO             , only : TDate, IncrDate, operator(+), operator(>), wrtgol
    use GO             , only : TrcFile, ReadRc
    use GO             , only : goReplace
    use GO             , only : T_SampStat
    use EMEP_NMC_Output, only : T_NMC_Output
    use UFFTW3         , only : T_UFFTW3_2d
    use C3PO           , only : Dimension
    use C3PO           , only : RealCoordinate, HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : IntegerCoordinate
    
      ! testing ...
      use file_nc
    
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eta_f), intent(inout)        ::  self
    type(TrcFile), intent(in)                   ::  rcF
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eta_f_Compute'
    
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
    
      ! testing ...
      character(len=256)  ::  fname
      integer             ::  ilon, ilat
    
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
        allocate( eta__varids(self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over variables:
        do ivar = 1, self%mvars%n
          ! get id:
          call eta_file%Get_VarID( trim(self%mvars%value(ivar)%name), &
                                     eta__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do

        ! storage:
        allocate( eta(nlon,nlat,nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( fixed(nlev,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop:
        do ivar = 1, self%mvars%n
          ! read:
          call eta_file%Get_IField1D( trim(self%mvars%value(ivar)%name)//'_fixed', &
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
        allocate( eta_f__varids(self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( fixed__varids(self%mvars%n), stat=status )
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
        do ivar = 1, self%mvars%n

          ! get units:
          call eta_file%Get_Att( eta__varids(ivar), 'units', units, status )
          IF_NOT_OK_RETURN(status=1)
          ! define:
          call eta_f_file%Def_CSample3D( trim(self%mvars%value(ivar)%name), trim(units), &
                       cplx, lon_f_coor, lat_f_coor, lev_coor, ctime_coor, eta_f_sample_coor, &
                       eta_f__varids(ivar), status, &
                       fill_value=eta_f__fill_value )
          IF_NOT_OK_RETURN(status=1)
          ! annote:
          call eta_f_file%Put_Att( eta_f__varids(ivar), 'long_name', trim(self%mvars%value(ivar)%name)//' spectral transform', status )
          IF_NOT_OK_RETURN(status=1)

          ! define:
          call eta_f_file%Def_IField1D( trim(self%mvars%value(ivar)%name)//'_fixed', trim(fixed__units), &
                                        lev_coor, &
                                        fixed__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
          call eta_f_file%Put_Att( fixed__varids(ivar), 'long_name', trim(self%mvars%value(ivar)%name)//' fixed layer (1=fixed, 0=dynamic)', status )
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
        do ivar = 1, self%mvars%n
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
        do ivar = 1, self%mvars%n
    
          ! info ...
          write (gol,'("    variable ",i0," ",a)') itime, trim(self%mvars%value(ivar)%name); call goPr
        
          ! read sample:
          call eta_file%Get_Sample3D( self%mvars%value(ivar)%name, itime, eta_isample, eta, units, status )
          IF_NOT_OK_RETURN(status=1)
          
          ! padding:
          eta_ex = 0.0
          eta_ex(1:nlon,1:nlat,1:nlev) = eta
          
          ! loop over levels:
          do ilev = 1, nlev
          
            !! testing ..
            !if ( itime == 2 ) then
            !  do ilon = 1, nlon
            !    eta_ex(ilon,:,ilev) = eta_ex(ilon     ,1,ilev)
            !  end do
            !  do ilon = nlon+1, nlon_ex
            !    eta_ex(ilon,:,ilev) = eta_ex(ilon-nlon,1,ilev)
            !  end do
            !end if
          
            ! forward real-to-(half-)complex fft:
            call ufft%Forward( eta_ex(:,:,ilev), eta_f(:,:,ilev), status )
            IF_NOT_OK_RETURN(status=1)
            
            !! testing ...
            !if ( (ivar == 1) .and. (ilev == nlev) ) then
            !  write (fname,'("dump-",i2.2,"-eta.nc")') itime
            !  call nc_dump( trim(fname), eta_ex(:,:,ilev), 'eta', (/'x','y'/), status )
            !  IF_NOT_OK_RETURN(status=1)
            !  write (fname,'("dump-",i2.2,"-eta-f-r.nc")') itime
            !  call nc_dump( trim(fname), real(eta_f(:,:,ilev)), 'eta_r', (/'x','y'/), status )
            !  IF_NOT_OK_RETURN(status=1)
            !  write (fname,'("dump-",i2.2,"-eta-f-i.nc")') itime
            !  call nc_dump( trim(fname), aimag(eta_f(:,:,ilev)), 'eta_i', (/'x','y'/), status )
            !  IF_NOT_OK_RETURN(status=1)
            !end if
            
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
      
      !! testing ...
      !write (gol,'("break after first day")'); call goErr
      !exit

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
    deallocate( eta_f__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( fixed__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ done with fft
    
    ! clear:
    deallocate( eta_ex, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( eta_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( lons_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( lats_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( kstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( theta, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( hcn, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
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
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( eta, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( eta__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( fixed, stat=status )
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
    call eta_sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close input:
    call eta_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    !write (gol,'("break after first eta_f")'); call goErr
    !RACEBACK; status=1; return
    
    ! ~

    ! ok
    status = 0
  
  end subroutine Driver_Eta_f_Compute
  
  
end module EMEP_NMC_Driver_Eta_f
