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

module EMEP_NMC_Driver_BB

  use GO                  , only : gol, goPr, goErr
  use GO                  , only : TDate

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_BB
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_BB'
  
  ! --- types ----------------------------------------
  
  type T_Driver_BB
    ! start time, used to open eta-f file:
    type(TDate)                       ::  tday1
    ! statistics files created or read:
    character(len=1024)               ::  eps_stats_filename
    character(len=1024)               ::  eta_f_samples_filename
    character(len=1024)               ::  gamma_filename
    character(len=1024)               ::  XLX_filename
    character(len=1024)               ::  BB_filename
  contains
    procedure   ::  Init            => Driver_BB_Init
    procedure   ::  Done            => Driver_BB_Done
    procedure   ::  Collect         => Driver_BB_Collect
  end type T_Driver_BB
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_BB_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_BB), intent(out)          ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_BB_Init'
    
    ! --- local ----------------------------------
    
    character(len=32)     ::  tvalue
    
    ! --- begin ----------------------------------
    
    ! start time:
    call ReadRc( rcF, 'nmc.timerange.t1', tvalue, status )
    IF_NOT_OK_RETURN(status=1)
    ! extract:
    call goReadFromLine( tvalue, self%tday1, status )
    IF_NOT_OK_RETURN(status=1)

    ! input/output:
    call ReadRc( rcF, 'nmc.eps.stats.filename', self%eps_stats_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.eta_f.samples.filename', self%eta_f_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.gamma.filename', self%gamma_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.XLX.filename', self%XLX_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.BB.filename', self%BB_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_BB_Init


  ! ***
  

  subroutine Driver_BB_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_BB), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_BB_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
  
  end subroutine Driver_BB_Done


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

  subroutine Driver_BB_Collect( self, status )
  
    use GO             , only : goReplace
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    use C3PO           , only : IntegerCoordinate
    use EMEP_BCovarSqrt, only : T_BCovarSqrt
    
    ! --- in/out ---------------------------------
    
    class(T_Driver_BB), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_BB_Collect'
  
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
    character(len=32)                ::  tracer
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

    ! ~ gamma
    !   (read this first to setup some coordiantes)
    
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

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=ntracer )
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
    allocate( fixed(nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call inp_file%Get_IField2D( 'fixed', fixed, fixed__units, status )
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

    ! storage:
    allocate( ps(nlon,nlat,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( S(nlon,nlat,nlev,ntime,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( S__units(ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over times:
    do itime = 1, ntime
      ! read surface pressure:
      call inp_file%Get_Field2D_Series( trim(vname_ps), itime, ps(:,:,itime), units, status )
      IF_NOT_OK_RETURN(status=1)
    end do

    ! loop over variables:
    do itracer = 1, ntracer
      !
      ! current:
      call tracer_coor%Get_Value( itracer, status, value=tracer )
      IF_NOT_OK_RETURN(status=1)
      ! variable:
      write (vname,'(a,"_stdv")') trim(tracer)
      ! loop over times:
      do itime = 1, ntime
        ! read:
        call inp_file%Get_Field3D_Series( trim(vname), itime, &
                             S(:,:,:,itime,itracer), S__units(itracer), status )
        IF_NOT_OK_RETURN(status=1)
      end do ! time
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
    
    ! create new result file, prepare for parallel i/o:
    call BB_file%Create( trim(self%BB_filename), status, fmt='netcdf4' )
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
    allocate( S__varids(ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over variables:
    do itracer = 1, ntracer
      ! current:
      call tracer_coor%Get_Value( itracer, status, value=tracer )
      IF_NOT_OK_RETURN(status=1)
      ! variable:
      write (vname,'("S_",a)') trim(tracer)
      ! define:
      call BB_file%Def_Field3D_Series( trim(vname), trim(S__units(itracer)), &
                              lon_coor, lat_coor, lev_coor, ctime_coor, &
                              S__varids(itracer), status )
      IF_NOT_OK_RETURN(status=1)
      call BB_file%Put_Att( S__varids(itracer), 'long_name', trim(tracer)//' grid point standard deviation', status )
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
    do itracer = 1, ntracer
      ! loop over hours:
      do itime = 1, ntime
        ! write; note different writing order!
        call BB_file%Put_Field3D_Series( S__varids(itracer), itime, S(:,:,:,itime,itracer), status )
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
            write (gol,'("          level ",i3," : ",f8.2,"    (max ",f8.2,")")') &
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
  
  end subroutine Driver_BB_Collect
  
  

end module EMEP_NMC_Driver_BB
