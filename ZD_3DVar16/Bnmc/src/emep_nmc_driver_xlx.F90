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

module EMEP_NMC_Driver_XLX

  use GO                  , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_XLX
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_XLX'
  
  ! --- types ----------------------------------------
  
  type T_Driver_XLX
    ! statistics files created or read:
    character(len=1024)               ::  B_f_filename
    character(len=1024)               ::  XLX_filename
  contains
    procedure   ::  Init            => Driver_XLX_Init
    procedure   ::  Done            => Driver_XLX_Done
    procedure   ::  Compute         => Driver_XLX_Compute
  end type T_Driver_XLX
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_XLX_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_XLX), intent(out)          ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_XLX_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! input/output:
    call ReadRc( rcF, 'nmc.B_f.filename', self%B_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.XLX.filename', self%XLX_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_XLX_Init


  ! ***
  

  subroutine Driver_XLX_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_XLX), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_XLX_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
  
  end subroutine Driver_XLX_Done


  ! ***
  
  
  !
  ! Compute eigenvalue decompo:
  !   B_f = X Lambda X'
  ! 

  subroutine Driver_XLX_Compute( self, status )
  
    use EMEP_NMC_LinAlg, only : Decompose_Bhat
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    use C3PO           , only : IntegerCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_Driver_XLX), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_XLX_Compute'
  
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
    integer                          ::  nv, iv
    integer                          ::  ik
    
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
    allocate( fixed(nlev,ntracer), stat=status )
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
      
        ! eigenvalue decomposition:
        call Decompose_Bhat( B_f(ik,:,:,:,:), &
                             X(ik,:,:,:), Lambda(ik,:), nev(ik), &
                             status )
        IF_NOT_OK_RETURN(status=1)
        
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
      
    end do  ! times

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
  
  end subroutine Driver_XLX_Compute


end module EMEP_NMC_Driver_XLX
