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

module EMEP_NMC_Driver_D_f

  use GO                  , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_D_f
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_D_f'
  
  ! --- types ----------------------------------------
  
  type T_Driver_D_f
    ! step size:
    real                              ::  akstar_step
    integer                           ::  ntheta90
    ! statistics files created or read:
    character(len=1024)               ::  C_f_filename
    character(len=1024)               ::  D_f_filename
  contains
    procedure   ::  Init            => Driver_D_f_Init
    procedure   ::  Done            => Driver_D_f_Done
    procedure   ::  Compute         => Driver_D_f_Compute
  end type T_Driver_D_f
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_D_f_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_D_f), intent(out)          ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_D_f_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! read stepping:
    call ReadRc( rcF, 'nmc.D_f.akstar_step', self%akstar_step, status )
    IF_NOT_OK_RETURN(status=1)

    ! number of theta segments per 90 deg:
    call ReadRc( rcF, 'nmc.D_f.ntheta', self%ntheta90, status )
    IF_NOT_OK_RETURN(status=1)

    ! input/output:
    call ReadRc( rcF, 'nmc.C_f.filename', self%C_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.D_f.filename', self%D_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_D_f_Init


  ! ***
  

  subroutine Driver_D_f_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_D_f), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_D_f_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
  
  end subroutine Driver_D_f_Done


  ! ***
  
  
  !
  ! Compute angular averages.
  ! 

  subroutine Driver_D_f_Compute( self, status )
  
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : IntegerCoordinate, LabelCoordinate
    use C3PO           , only : RealCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_Driver_D_f), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_D_f_Compute'
  
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
    character(len=64)                ::  kstar__formula
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
    integer                          ::  nakstar
    real, allocatable                ::  akstar_bnds(:,:)  ! (nv,nakstar)
    integer                          ::  ik
    integer, allocatable             ::  iakstar(:,:)      ! (nlon_f,nlat_f)
    integer                          ::  iakstar__fill_value
    integer                          ::  iakstar__varid
    integer                          ::  i, j
    type(RealCoordinate)             ::  akstar_coor
    
    integer                          ::  ntheta
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
    allocate( fixed(nlev,ntracer), stat=status )
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
    
    ! number of steps:
    nakstar = int(ceiling( kstar_max / self%akstar_step ))

    ! storage for boundaries:
    allocate( akstar_bnds(2,nakstar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do ik = 1, nakstar
      ! lower and upper bound:
      akstar_bnds(1,ik) = (ik-1) * self%akstar_step
      akstar_bnds(2,ik) =  ik    * self%akstar_step
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
        ik = int(ceiling( kstar(i,j) / self%akstar_step ))
        ! first is zero, add to first bin:
        if ( kstar(i,j) == 0.0 ) ik = 1
        ! check ..
        if ( (ik < 1) .or. (ik > nakstar) ) then
          write (gol,'("unexpected index ",i0," in ak* array with step ",f8.3)') ik, self%akstar_step; call goErr
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
    
    ! bins are centered around -pi/2,..,0,..,pi/2,
    ! thus segments [-pi/2-DeltaTheta/2,-pi/2+DeltaTheta/2] etc
    ! number of segments:
    ntheta = 2*self%ntheta90 + 1
    ! initial spacing:
    DeltaTheta0 = (0.5*pi) / real(self%ntheta90)   ! rad
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
  
  end subroutine Driver_D_f_Compute

end module EMEP_NMC_Driver_D_f
