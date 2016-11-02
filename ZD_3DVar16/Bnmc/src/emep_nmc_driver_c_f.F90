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

module EMEP_NMC_Driver_C_f

  use GO                  , only : gol, goPr, goErr
  use EMEP_NMC_Common     , only : T_ModelVars
  use GO                  , only : TDate

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_C_f
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_C_f'
  
  ! --- types ----------------------------------------
  
  type T_Driver_C_f
    ! time range:
    type(TDate)                       ::  tday1, tday2
    ! variables:
    type(T_ModelVars)                 ::  mvars
    ! statistics files created or read:
    character(len=1024)               ::  eta_f_samples_filename
    character(len=1024)               ::  C_f_filename
  contains
    procedure   ::  Init            => Driver_C_f_Init
    procedure   ::  Done            => Driver_C_f_Done
    procedure   ::  Compute         => Driver_C_f_Compute
  end type T_Driver_C_f
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_C_f_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_C_f), intent(out)          ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_C_f_Init'
    
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
    call ReadRc( rcF, 'nmc.eta_f.samples.filename', self%eta_f_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.C_f.filename', self%C_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_C_f_Init


  ! ***
  

  subroutine Driver_C_f_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_C_f), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_C_f_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%mvars%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_C_f_Done


  ! ***
  
  
  !
  ! Compute spectral covariances
  !

  subroutine Driver_C_f_Compute( self, status )
  
    use GO             , only : TDate, IncrDate, operator(+), operator(>), wrtgol
    use GO             , only : goReplace
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : IntegerCoordinate, LabelCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_Driver_C_f), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_C_f_Compute'
    
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
    character(len=64)                ::  kstar__formula
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
      call goReplace( filename, '%{hh}'  , '(i2.2)', tday%hour , status )
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
        allocate( eta_f_sample(nlon_f,nlat_f,nlev,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
    
        ! storage:
        allocate( fixed(nlev,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over variables:
        do ivar = 1, self%mvars%n
          ! read:
          call eta_f_file%Get_IField1D( trim(self%mvars%value(ivar)%name)//'_fixed', fixed(:,ivar), fixed__units, status )
          IF_NOT_OK_RETURN(status=1)
        end do
        
        ! ~ sample averages
        
        ! storage:
        allocate( ps(nlon_f,nlat_f), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( C_f(nlon_f,nlat_f,nlev,nlev,self%mvars%n,self%mvars%n,ntime), stat=status )
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
        do ivar = 1, self%mvars%n
          ! read sample:
          call eta_f_file%Get_CSample3D( self%mvars%value(ivar)%name, itime, isample, &
                                            eta_f_sample(:,:,:,ivar), units, status, &
                                          fill_value=eta_f_sample__fill_value )
          IF_NOT_OK_RETURN(status=1)
        end do
    
        ! info ...
        write (gol,'("    add contribution to covariance ...")'); call goPr
        
        ! loops over upper triangle of tracers:
        do ivar1 = 1, self%mvars%n
          do ivar2 = ivar1, self%mvars%n
          
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
  
  end subroutine Driver_C_f_Compute

  
end module EMEP_NMC_Driver_C_f
  
