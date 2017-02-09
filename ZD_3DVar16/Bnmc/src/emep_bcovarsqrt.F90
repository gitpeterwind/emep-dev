!######################################################################
!
! EMEP DA Background Covriance SquareRoot
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!######################################################################

module EMEP_BCovarSqrt

  use GO    , only : gol, goPr, goErr
  use UFFTW3, only : T_UFFTW3_2d

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_BCovarSqrt
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_BCovarSqrt'
  

  ! --- types ----------------------------------------
  
  type T_BCovarSqrt
    ! time values:
    integer                           ::  ntime
    integer, allocatable              ::  hour(:)    ! (ntime)
    ! tracers:
    integer                           ::  ntracer
    character(len=16), allocatable    ::  tracer_name(:)  ! (ntracer)
    ! subset of tracers, might be in other order:
    integer                           ::  ntr
    integer, allocatable              ::  itr(:)      ! (ntracer)
    integer, allocatable              ::  itracer(:)  ! (ntr)
    ! grid dimensions:
    integer                           ::  nlon, nlat
    integer                           ::  nlat_local, ilat_offset
    ! extended grid:
    integer                           ::  nlon_ex, nlat_ex
    integer                           ::  nlat_ex_local, ilat_ex_offset
    ! levels:
    integer                           ::  nlev
    integer, allocatable              ::  fixed(:,:)  ! (nlev,ntracer)
    ! std.dev.
    real, allocatable                 ::  S(:,:,:,:,:) ! (nlon,nlat,nlev,ntracer,ntime)
    character(len=32), allocatable    ::  S__units(:)  ! (ntracer)
    ! wavenumbers:
    integer                           ::  nakstar
    real, allocatable                 ::  akstar(:)         ! (nakstar)
    real, allocatable                 ::  akstar_bnds(:,:)  ! (2,nakstar)
    ! fft:
    integer                           ::  nlon_f, nlat_f
    integer                           ::  nlat_f_local, ilat_f_offset
    type(T_UFFTW3_2d)                 ::  ufft
    real, allocatable                 ::  kstar(:,:)    ! (nlon_f,nlat_f)
    integer, allocatable              ::  iakstar(:,:)  ! (nlon_f,nlat_f)
    ! eigenvalue decompo:
    integer                           ::  nv
    integer, allocatable              ::  nev(:,:)    ! (nakstar,ntime)
    real, allocatable                 ::  GXL(:,:,:,:,:)  ! (nakstar,nlev,ntracer,nv,ntime)
    ! size of 1D transformed state:
    integer, allocatable              ::  nw(:)         ! (ntime)
    integer, allocatable              ::  nw_local(:)   ! (ntime)
    ! number of complex numbers represented in half-complex storage:
    integer, allocatable              ::  hcn_local(:,:)      ! (max(nw_local),ntime)
    ! compensation factor for truncations:
    real, allocatable                 ::  sqrt_phi(:,:,:)  ! (nlev,ntracer,ntime)
    !
  contains
    procedure   ::  Init            => BCovarSqrt_Init
    procedure   ::  Done            => BCovarSqrt_Done
    procedure   ::  Check           => BCovarSqrt_Check
    procedure   ::  Get             => BCovarSqrt_Get
    procedure   ::  FindTime        => BCovarSqrt_FindTime
    procedure   ::  LevelTracerGet  => BCovarSqrt_LevelTracerGet
    procedure   ::  TracerGet       => BCovarSqrt_TracerGet
    procedure   ::  TimeGet         => BCovarSqrt_TimeGet
    procedure   ::  Forward         => BCovarSqrt_Forward
    procedure   ::  Reverse         => BCovarSqrt_Reverse
  end type T_BCovarSqrt
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  !
  ! Read covariance entities from file.
  !
  ! Optional arguments:
  !
  !   comm     : MPI comunicator to allow parallel i/o
  !
  !   tracers=(/'O3 ','NO2'/)  : 
  !      Subset of tracers, or same tracers in other order.
  !      If not provided, the original tracer (order) is used.
  !      Names are matched with input tracers to check and set mapping.
  !      The grid array passed to 'Forward' should have tracers in
  !      the order of the subset (or the original tracers).
  !      Similar, the grid array returned from Reverse has this order too.
  !

  subroutine BCovarSqrt_Init( self, filename, status, &
                                comm, tracers )

    use GO             , only : TDate, Get
    use GO             , only : goMatchValue
    use C3PO           , only : Datafile
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    use C3PO           , only : IntegerCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(out)          ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  comm
    character(len=*), intent(in), optional    ::  tracers(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_Init'
    
    ! --- local ----------------------------------
    
    type(Datafile)                   ::  inp_file
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(Dimension)                  ::  lon_ex_dim, lat_ex_dim
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(LabelCoordinate)            ::  tracer_coor
    type(TimeCoordinate)             ::  ctime_coor
    type(RealCoordinate)             ::  akstar_coor
    integer                          ::  itime
    type(TDate)                      ::  t
    integer                          ::  itracer
    integer                          ::  itr
    integer                          ::  iakstar
    integer                          ::  ilev
    character(len=32)                ::  vname, units
    integer                          ::  i, j, k
    integer                          ::  j_local
    type(IntegerCoordinate)          ::  eval_coor
    real, allocatable                ::  gamma(:,:,:,:)  ! (nakstar,nlev,ntracer,ntime)
    real, allocatable                ::  Lambda(:,:,:)  ! (nakstar,nv,ntime)
    real, allocatable                ::  X(:,:,:,:,:)   ! (nakstar,nlev,ntracer,nv,ntime)
    integer                          ::  iev
    integer, allocatable             ::  hcn_local(:,:)  ! (nxf,nyf_local)
    integer                          ::  nwmax
    integer                          ::  iw
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("read covariance square root entities ...")'); call goPr
    
    ! ~ open file
    
    !! open file with BB entities:
    !if ( present(comm) ) then
    !  ! parallel input:
    !  call inp_file%Open( trim(filename), status, comm=comm )
    !  IF_NOT_OK_RETURN(status=1)
    !else
      ! sequential input:
      call inp_file%Open( trim(filename), status )
      IF_NOT_OK_RETURN(status=1)
    !end if
    
    ! ~ coordinates
    
    ! create hour coordinate:
    call ctime_coor%Init( 'time', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call ctime_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call ctime_coor%Get_Dim( status, n=self%ntime )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%hour(self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over time values:
    do itime = 1, self%ntime
      ! get time value:
      call ctime_coor%Get_Value( itime, status, t=t )
      IF_NOT_OK_RETURN(status=1)
      ! extract hour:
      call Get( t, hour=self%hour(itime) )
    end do

    ! create longitude coordinate:
    call lon_coor%Init( 'longitude', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lon_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lon_coor%Get_Dim( status, n=self%nlon )
    IF_NOT_OK_RETURN(status=1)

    ! create latitude coordinate:
    call lat_coor%Init( 'latitude', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lat_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lat_coor%Get_Dim( status, n=self%nlat )
    IF_NOT_OK_RETURN(status=1)

    ! create level coordinate:
    call lev_coor%Init( 'lev', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lev_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lev_coor%Get_Dim( status, n=self%nlev )
    IF_NOT_OK_RETURN(status=1)

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=self%ntracer )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%tracer_name(self%ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over tracers:
    do itracer = 1, self%ntracer
      ! get tracer name:
      call tracer_coor%Get_Value( itracer, status, value=self%tracer_name(itracer) )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! ~ fixed layers
      
    ! storage:
    allocate( self%fixed(self%nlev,self%ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    !call inp_file%Get_IField2D( 'fixed', self%fixed, units, status )
    call inp_file%Get_Var( 'fixed', self%fixed, units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ subset of tracers
    
    ! number of sub-tracer selection:
    if ( present(tracers) ) then
      self%ntr = size(tracers)
    else
      self%ntr = self%ntracer
    end if
    ! storage for mapping:
    allocate( self%itr(self%ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%itracer(self%ntr), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill mapping:
    do itr = 1, self%ntr
      ! obtain index of original tracer:
      if ( present(tracers) ) then
        ! match with original names:
        call goMatchValue( tracers(itr), self%tracer_name, itracer, status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! 1-1 copy:
        itracer = itr
      end if
      ! store:
      self%itracer(itr) = itracer
      self%itr(itracer) = itr
    end do
    
    ! ~ extended grid coordinates
    
    ! create extended dimensin:
    call lon_ex_dim%Init( 'lon_ex', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lon_ex_dim%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lon_ex_dim%Get_Dim( status, n=self%nlon_ex )
    IF_NOT_OK_RETURN(status=1)
    
    ! create extended dimension:
    call lat_ex_dim%Init( 'lat_ex', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lat_ex_dim%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lat_ex_dim%Get_Dim( status, n=self%nlat_ex )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ average wavenumber coordinate
    
    ! create average k* coordinate:
    call akstar_coor%Init( 'akstar', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call akstar_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call akstar_coor%Get_Dim( status, n=self%nakstar )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%akstar(self%nakstar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%akstar_bnds(2,self%nakstar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! get:
    call akstar_coor%Get_Values( status, values=self%akstar, value_bnds=self%akstar_bnds )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ fft
    
    ! init transform for extended grid:
    call self%ufft%Init( self%nlon_ex, self%nlat_ex, status, &
                           comm=comm, unitairy=.true. )
    IF_NOT_OK_RETURN(status=1)

    ! spectral size in half-complex representation:
    call self%ufft%Get( status, nxf=self%nlon_f, nyf=self%nlat_f )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ domain decomposition
    
    ! decomposition in 2nd dimension:
    call self%ufft%Get( status, nyf_local=self%nlat_f_local, iyf_offset=self%ilat_f_offset )
    IF_NOT_OK_RETURN(status=1)  

    ! copy, same decompositions used for grid fields:
    self%nlat_ex_local  = self%nlat_f_local
    self%ilat_ex_offset = self%ilat_f_offset
    
    ! decomposition of grid fields;
    ! local size on normal domain, some slabs might be completely in extension:
    self%nlat_local = max( 0, min(self%ilat_ex_offset+self%nlat_ex_local,self%nlat) - self%ilat_ex_offset )
    ! use same offset as in extended grid, or to zero if in extension:
    if ( self%nlat_local > 0 ) then
      self%ilat_offset = self%ilat_ex_offset
    else
      self%ilat_offset = 0
    end if
    
    ! ~ S

    ! storage:
    allocate( self%S(self%nlon,max(self%nlat_local,1),self%nlev,self%ntracer,self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%S__units(self%ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over tracers:
    do itracer = 1, self%ntracer
      ! variable:
      write (vname,'("S_",a)') trim(self%tracer_name(itracer))
      ! read values and units (collective call, each slab has at least 1 y):
      call inp_file%Get_Var( trim(vname), self%S(:,:,:,itracer,:), self%S__units(itracer), status, &
                               start=(/        1,    self%ilat_offset+1,        1,           1/), &
                               count=(/self%nlon,max(self%nlat_local,1),self%nlev,self%ntime/) )
      IF_NOT_OK_RETURN(status=1)
    end do ! tracers
    
    ! ~ kstar

    ! storage:
    allocate( self%kstar(self%nlon_f,self%nlat_f), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    call self%ufft%Get( status, kstar=self%kstar )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ iakstar
    
    ! storage:
    allocate( self%iakstar(self%nlon_f,self%nlat_f), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over spectral elements:
    do j = 1, self%nlat_f
      do i = 1, self%nlon_f
        ! in range?
        if ( self%kstar(i,j) <= self%akstar_bnds(2,self%nakstar) ) then
          ! for safety check ...
          self%iakstar(i,j) = -999
          ! loop over wavenumber bands:
          do k = 1, self%nakstar
            ! match ?
            if ( (self%kstar(i,j) >= self%akstar_bnds(1,k)) .and. &
                 (self%kstar(i,j) <= self%akstar_bnds(2,k))       ) then
              ! store:
              self%iakstar(i,j) = k
              ! leave:
              exit
            end if  ! match
          end do  ! wavenumber bands
          ! check ...
          if ( self%iakstar(i,j) < 1 ) then
            write (gol,'("no wavenumber band found for coeff (",i0,",",i0,") with k* ",f0.2)') &
                     i, j, self%iakstar(i,j); call goErr
            TRACEBACK; status=1; return
          end if
        else
          ! outside eliptic truncation:
          self%iakstar(i,j) = -999
        end if
      end do  ! i
    end do ! j
    
    ! ~ gamma

    ! temporary storage:
    allocate( gamma(self%nakstar,self%nlev,self%ntracer,self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read scaling:
    call inp_file%Get_Var( 'gamma', gamma, units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ eigenvalue decomposition
    
    ! create eigenvalue coordinate:
    call eval_coor%Init( 'eigenvalue', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call eval_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call eval_coor%Get_Dim( status, n=self%nv )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( self%nev(self%nakstar,self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
      
    ! temporary storage:
    allocate( Lambda(self%nakstar,self%nv,self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( X(self%nakstar,self%nlev,self%ntracer,self%nv,self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! read number of significant eigenvalues:
    call inp_file%Get_Var( 'nev', self%nev, units, status )
    IF_NOT_OK_RETURN(status=1)
    ! read eigenvalues:
    call inp_file%Get_Var( 'Lambda', Lambda, units, status )
    IF_NOT_OK_RETURN(status=1)
    ! read eigenvectors:
    call inp_file%Get_Var( 'X', X, units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! precompue  GXL = sqrt(Gamma) X sqrt(Lambda) ;
    ! storage:
    allocate( self%GXL(self%nakstar,self%nlev,self%ntracer,self%nv,self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! init to zero:
    self%GXL = 0.0
    ! loop over time:
    do itime = 1, self%ntime
      ! loop over wavenumbers:
      do iakstar = 1, self%nakstar
        ! loop over eigenvalues:
        do iev = 1, self%nev(iakstar,itime)
          ! loop over tracer:
          do itracer = 1, self%ntracer
            ! loop over levels:
            do ilev = 1, self%nlev
              ! only non-zero for non-fixed layers:
              if ( self%fixed(ilev,itracer) == 0 ) then
                ! fill element of matrix product:
                self%GXL(iakstar,ilev,itracer,iev,itime) = &
                       sqrt( gamma(iakstar,ilev,itracer    ,itime)) * &
                                 X(iakstar,ilev,itracer,iev,itime)  * &
                       sqrt(Lambda(iakstar             ,iev,itime))
              end if ! zero levels
            end do ! level
          end do ! tracer
        end do ! eigenvalue
      end do ! ak*
    end do ! time
    
    ! clear temporary storage:
    deallocate( Lambda, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( X, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ clear gamma
    
    ! clear temporary storage:
    deallocate( gamma, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ phi
    
    ! storage for compensation factors, might differe per hour:
    allocate( self%sqrt_phi(self%nlev,self%ntracer,self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! read compenstation factor, this is the 'squared' version!
    call inp_file%Get_Var( 'phi', self%sqrt_phi, units, status )
    IF_NOT_OK_RETURN(status=1)
    ! convert:
    self%sqrt_phi = sqrt( self%sqrt_phi )
    
    ! ~ close
    
    ! close:
    call inp_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ 1D states
    
    ! info ...
    write (gol,'("  count number of elements in 1D complex state ...")'); call goPr
    write (gol,'("    local spectral slab : ",i0,"+1:",i0)') self%ilat_f_offset, self%nlat_f_local; call goPr
    ! storage:
    allocate( self%nw(self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%nw_local(self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! init to zero:
    self%nw       = 0
    self%nw_local = 0
    ! loop over time
    do itime = 1, self%ntime
      ! loop over spectral coeff:
      do j = 1, self%nlat_f
        do i = 1, self%nlon_f
          ! current:
          iakstar = self%iakstar(i,j)
          ! skip if outside eliptic truncation:
          if ( iakstar < 1 ) cycle
          ! loop over eigenvalues:
          do iev = 1, self%nev(iakstar,itime)
            ! increase global counter:
            self%nw(itime) = self%nw(itime) + 1
            ! local slab?
            if ( self%nlat_f_local > 0 ) then
              ! increase local counter if in slab:
              if ( (self%ilat_f_offset+1 <= j) .and. (j <= self%ilat_f_offset+self%nlat_f_local) ) then
                self%nw_local(itime) = self%nw_local(itime) + 1
              end if
            end if ! local slab
          end do ! eigenvalues
        end do ! i
      end do ! j
      ! info ...
      write (gol,'("    time record ",i0)') itime; call goPr
      write (gol,'("      number of grid point elements : ",i0)') &
              self%nlon * self%nlat * self%nlev * self%ntracer; call goPr
      write (gol,'("      number of spectral   elements : ",i0," (local ",i0,")")') &
              self%nw(itime), self%nw_local(itime); call goPr
    end do ! time
    
    ! storage:
    allocate( hcn_local(self%nlon_f,self%nlat_f_local), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! 2D field with weights:
    call self%ufft%Get( status, hcn_local=hcn_local )
    IF_NOT_OK_RETURN(status=1)
    ! maximum number of local hcr elements, minimum of 1:
    nwmax = max( 1, maxval( self%nw_local ) )
    ! storage:
    allocate( self%hcn_local(nwmax,self%ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! init for safety:
    self%hcn_local = 0
    ! loop over time
    do itime = 1, self%ntime
      ! init counter in 1D complex state:
      iw = 0
      ! loop over local 2nd dimension of spectral coeff:
      do j_local = 1, self%nlat_f_local
        ! global index:
        j = self%ilat_f_offset + j_local
        ! loop over 1st dimension:
        do i = 1, self%nlon_f
          ! index of wavenumber band:
          iakstar = self%iakstar(i,j)
          ! skip if outside eliptic truncation:
          if ( iakstar < 1 ) cycle
          ! loop over eigenvalues:
          do iev = 1, self%nev(iakstar,itime)
            ! increase counter:
            iw = iw + 1
            ! copy weight:
            self%hcn_local(iw,itime) = hcn_local(i,j_local)
          end do ! eigenvalues
        end do ! i
      end do ! j
      ! check ...
      if ( iw /= self%nw_local(itime) ) then
        write (gol,'("only ",i0," of ",i0," elements of hcn_local filled")') iw, self%nw_local(itime); call goErr
        TRACEBACK; status=1; return
      end if
    end do ! time    
    ! clear:
    deallocate( hcn_local, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~
    
    ! ok
    status = 0
  
  end subroutine BCovarSqrt_Init


  ! ***
  

  subroutine BCovarSqrt_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%hour, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%tracer_name, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%itr, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%itracer, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%fixed, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%S, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%S__units, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%akstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%akstar_bnds, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%kstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%iakstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%nev, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%GXL, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%nw, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%nw_local, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( self%hcn_local, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%sqrt_phi, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine BCovarSqrt_Done


  ! ***
  
  
  !
  ! Return time record 'itime' in 1,..,ntime that corresponds to hour
  !

  subroutine BCovarSqrt_FindTime( self, hour, itime, status )
  
    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(in)                         ::  hour
    integer, intent(out)                        ::  itime
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_FindTime'
    
    ! --- local ----------------------------------
    
    integer    ::  k
    
    ! --- begin ----------------------------------
    
    ! init result:
    itime = -999
    ! loop:
    do k = 1, self%ntime
      ! match?
      if ( self%hour(k) == hour ) then
        ! copy:
        itime = k
        ! leave:
        exit
      end if
    end do ! time records
    ! check ...
    if ( itime < 1 ) then
      write (gol,'("could not find hour ",i2.2," in time records:")') hour; call goErr
      do k = 1, self%ntime
        write (gol,'("  record ",i2," hour ",i2.2)') k, self%hour(k); call goErr
      end do
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
  
  end subroutine BCovarSqrt_FindTime


  ! ***
  
  
  subroutine BCovarSqrt_Check( self, status, &
                                nlon, nlat, nlon_ex, nlat_ex, &
                                nlev, ntracer )
  
    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(out)                        ::  status
    
    integer, intent(in), optional               ::  nlon
    integer, intent(in), optional               ::  nlat
    integer, intent(in), optional               ::  nlon_ex
    integer, intent(in), optional               ::  nlat_ex
    integer, intent(in), optional               ::  nlev
    integer, intent(in), optional               ::  ntracer

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_Check'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check?
    if ( present(nlon) ) then
      ! compare:
      if ( nlon /= self%nlon ) then
        write (gol,'("nlon in covariance is ",i0," while expected ",i0)') &
                self%nlon, nlon; call goErr
        status=1
      end if
    end if
    
    ! check?
    if ( present(nlat) ) then
      ! compare:
      if ( nlat /= self%nlat ) then
        write (gol,'("nlat in covariance is ",i0," while expected ",i0)') &
                self%nlat, nlat; call goErr
        status=1
      end if
    end if
    
    ! check?
    if ( present(nlon_ex) ) then
      ! compare:
      if ( nlon_ex /= self%nlon_ex ) then
        write (gol,'("nlon_ex in covariance is ",i0," while expected ",i0)') &
                self%nlon_ex, nlon_ex; call goErr
        status=1
      end if
    end if
    
    ! check?
    if ( present(nlat_ex) ) then
      ! compare:
      if ( nlat_ex /= self%nlat_ex ) then
        write (gol,'("nlat_ex in covariance is ",i0," while expected ",i0)') &
                self%nlat_ex, nlat_ex; call goErr
        status=1
      end if
    end if
    
    ! check?
    if ( present(nlev) ) then
      ! compare:
      if ( nlev /= self%nlev ) then
        write (gol,'("nlev in covariance is ",i0," while expected ",i0)') &
                self%nlev, nlev; call goErr
        status=1
      end if
    end if
    
    ! check?
    if ( present(ntracer) ) then
      ! compare:
      if ( ntracer /= self%ntracer ) then
        write (gol,'("ntracer in covariance is ",i0," while expected ",i0)') &
                self%ntracer, ntracer; call goErr
        status=1
      end if
    end if
    
    ! any errors ?
    if ( status /= 0 ) then
      TRACEBACK; status=1; return
    end if
    
    ! ok
    status = 0
  
  end subroutine BCovarSqrt_Check


  ! ***
  

  subroutine BCovarSqrt_Get( self, status, &
                                nlon, nlat, nlat_local, ilat_offset, &
                                nlon_ex, nlat_ex, nlat_ex_local, ilat_ex_offset, &
                                ntracer )
  
    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(out)                        ::  status
    
    integer, intent(out), optional              ::  nlon
    integer, intent(out), optional              ::  nlat
    integer, intent(out), optional              ::  nlat_local
    integer, intent(out), optional              ::  ilat_offset
    integer, intent(out), optional              ::  nlon_ex
    integer, intent(out), optional              ::  nlat_ex
    integer, intent(out), optional              ::  nlat_ex_local
    integer, intent(out), optional              ::  ilat_ex_offset
    integer, intent(out), optional              ::  ntracer

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! grid:
    if ( present(nlon       ) ) nlon        = self%nlon
    if ( present(nlat       ) ) nlat        = self%nlat
    if ( present(nlat_local ) ) nlat_local  = self%nlat_local
    if ( present(ilat_offset) ) ilat_offset = self%ilat_offset
    
    ! extended grid:
    if ( present(nlon_ex       ) ) nlon_ex        = self%nlon_ex
    if ( present(nlat_ex       ) ) nlat_ex        = self%nlat_ex
    if ( present(nlat_ex_local ) ) nlat_ex_local  = self%nlat_ex_local
    if ( present(ilat_ex_offset) ) ilat_ex_offset = self%ilat_ex_offset
    
    ! number of tracer elements:
    if ( present(ntracer) ) ntracer = self%ntracer
    
    ! ok
    status = 0
  
  end subroutine BCovarSqrt_Get


  ! ***
  

  subroutine BCovarSqrt_LevelTracerGet( self, ilev, itracer, status, fixed )
  
    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(in)                         ::  ilev
    integer, intent(in)                         ::  itracer
    integer, intent(out)                        ::  status
    logical, intent(out), optional              ::  fixed

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_LevelTracerGet'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (itracer < 1) .or. (itracer > self%ntracer) ) then
      write (gol,'("tracer index ",i0," not in range 1 .. ",i0)') itracer, self%ntracer; call goErr
      TRACEBACK; status=1; return
    end if    
    
    ! check ...
    if ( (ilev < 1) .or. (ilev > self%nlev) ) then
      write (gol,'("level index ",i0," not in range 1 .. ",i0)') ilev, self%nlev; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! return flag?
    if ( present(fixed) ) fixed = self%fixed(ilev,itracer) == 1
    
    ! ok
    status = 0
  
  end subroutine BCovarSqrt_LevelTracerGet


  ! ***
  

  subroutine BCovarSqrt_TimeGet( self, itime, status, &
                                    nw, nw_local, hcn_local, &
                                    S )
  
    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(in)                         ::  itime
    integer, intent(out)                        ::  status
    
    integer, intent(out), optional              ::  nw, nw_local
    integer, intent(out), optional              ::  hcn_local(:)  ! (nw_local)
    real, intent(out), optional                 ::  S(:,:,:,:) ! (nlon,nlat,nlev,ntracer)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_TimeGet'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (itime < 1) .or. (itime > self%ntime) ) then
      write (gol,'("time index ",i0," not in range 1 .. ",i0)') itime, self%ntime; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! return 1D sizes?
    if ( present(nw        ) ) nw         = self%nw        (itime)
    if ( present(nw_local  ) ) nw_local   = self%nw_local  (itime)
    
    ! return weights?
    if ( present(hcn_local) ) then
      ! only if local slab is present:
      if ( self%nw_local(itime) > 0 ) then
        ! check ...
        if ( size(hcn_local) /= self%nw_local(itime) ) then
          write (gol,'("size of hcn argument is ",i0," while expected ",i0)') &
                  size(hcn_local), self%nw_local(itime); call goErr
          TRACEBACK; status=1; return
        end if
        ! copy:
        hcn_local = self%hcn_local(1:self%nw_local(itime),itime)
      end if
    end if
    
    ! return std.dev.?
    if ( present(S) ) S = self%S(:,:,:,:,itime)
    
    ! ok
    status = 0
  
  end subroutine BCovarSqrt_TimeGet


  ! ***
  

  subroutine BCovarSqrt_TracerGet( self, itracer, status, name, units )
  
    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(in)                         ::  itracer
    integer, intent(out)                        ::  status
    character(len=*), intent(out), optional     ::  name
    character(len=*), intent(out), optional     ::  units

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_TracerGet'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (itracer < 1) .or. (itracer > self%ntracer) ) then
      write (gol,'("tracer index ",i0," not in range 1 .. ",i0)') itracer, self%ntracer; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! return tracer name?
    if ( present(name) ) name = trim(self%tracer_name(itracer))
    
    ! return tracer units?
    if ( present(units) ) units = trim(self%S__units(itracer))
    
    ! ok
    status = 0
  
  end subroutine BCovarSqrt_TracerGet


  ! ***
  

  !
  ! w =                 B^{H/2}                  x
  !
  !   = [ (Gamma^{1/2} X^T Lambda^{1/2}) F E S ] x
  !
  ! With 'correlation=.true.' the std.dev. matrix "S" is ommitted.
  !
  
  subroutine BCovarSqrt_Forward( self, itime, x, w, status, correlation )

    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(in)                         ::  itime
    real, intent(in)                            ::  x(:,:,:,:) ! (nlon,nlat,nlev,ntr)
    complex, intent(out)                        ::  w(:)  ! (nw)
    integer, intent(out)                        ::  status
    logical, intent(in), optional               ::  correlation

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_Forward'
    
    ! --- local ----------------------------------
    
    logical               ::  corr
    real, allocatable     ::  x_ex(:,:,:,:)  ! (nlon_ex,nlat_ex,nlev,ntr)
    complex, allocatable  ::  x_f(:,:,:,:)  ! (nlon_f,nlat_f,nlev,ntr)
    integer               ::  i, j
    integer               ::  j_local
    integer               ::  iakstar
    integer               ::  itracer, itr
    integer               ::  ilev
    integer               ::  iev
    integer               ::  iw
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (itime < 1) .or. (itime > self%ntime) ) then
      write (gol,'("time index ",i0," not in range 1 .. ",i0)') itime, self%ntime; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( self%nlat_local > 0 ) then
      if ( any( shape(x) /= (/self%nlon,self%nlat_local,self%nlev,self%ntr/) ) ) then
        write (gol,'("shape of input state is (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
                 shape(x), (/self%nlon,self%nlat_local,self%nlev,self%ntr/); call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! check ...
    if ( self%nw_local(itime) > 0 ) then
      if ( size(w) /= self%nw_local(itime) ) then
        write (gol,'("size of 1D complex output (local) is ",i0," while expected ",i0)') &
                 size(w), self%nw_local(itime); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! ~ E S x
    
    ! flag:
    corr = .false.
    if ( present(correlation) ) corr = correlation
    
    ! storage for extended grid:
    allocate( x_ex(self%nlon_ex,self%nlat_ex_local,self%nlev,self%ntr), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! init everything with zero for padding:
    x_ex = 0.0
    ! local slab ?
    if ( self%nlat_local > 0 ) then
      ! fill first part with x or S x :
      if ( corr ) then
        x_ex(1:self%nlon,1:self%nlat_local,:,:) =                         x
      else
        x_ex(1:self%nlon,1:self%nlat_local,:,:) = self%S(:,:,:,:,itime) * x
      end if
    end if ! local slab
    
    ! ~ F (E S x)
    
    ! storage for spectral transform:
    allocate( x_f(self%nlon_f,self%nlat_f_local,self%nlev,self%ntr), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! spectral transforms;
    ! loop over trs:
    do itr = 1, self%ntr
      ! original:
      itracer = self%itracer(itr)
      ! loop over layers:
      do ilev = 1, self%nlev
        ! skip fixed layers:
        if ( self%fixed(ilev,itracer) == 1 ) then
          ! padd with zeros:
          x_f(:,:,ilev,itr) = (0.0,0.0)
        else
          ! forward real-to-(half-)complex fft:
          call self%ufft%Forward( x_ex(:,:,ilev,itr), x_f(:,:,ilev,itr), status )
          IF_NOT_OK_RETURN(status=1)
          ! NOTE: elliptic truncation is not explicitly needed since:
          !  - following operations skip truncated elements
          !  - we do not need to save this variable, so no no nice mask needed
          ! While testing ...
          call self%ufft%Truncate( x_f(:,:,ilev,itr), status, fill_value=0.0 )
          IF_NOT_OK_RETURN(status=1)
        end if
      end do  ! levels
    end do  ! tracers
    
    ! compensate for truncations; 
    ! loop over subtracers:
    do itr = 1, self%ntr
      ! original:
      itracer = self%itracer(itr)
      ! loop over levels:
      do ilev = 1, self%nlev
        ! skip if fixed:
        if ( self%fixed(ilev,itr)==1 ) cycle
        ! scale:
        x_f(:,:,ilev,itr) = x_f(:,:,ilev,itr) * self%sqrt_phi(ilev,itracer,itime)
      end do ! levels
    end do ! tracers
    
    ! ~  (Lambda^{1/2} X^T Gamma^{1/2})  #        (F E S x)
    !     (nev x nlev*ntracer)(k*(i,j))  #  (nlev*ntracer x 1)(i,j)

    ! init 1D index:
    iw = 0
    ! loop over 2nd dimension of local spectral coeff:
    do j_local = 1,self%nlat_f_local
      ! global 2nd index:
      j = self%ilat_f_offset + j_local
      ! loop over 1st dimension:
      do i = 1, self%nlon_f
        ! wavenumber band:
        iakstar = self%iakstar(i,j)
        ! skip if outside eliptic truncation:
        if ( iakstar < 1 ) cycle
        ! loop over eigenvalues:
        do iev = 1, self%nev(iakstar,itime)
          ! increase index:
          iw = iw + 1
          ! init sum:
          w(iw) = (0.0,0.0)
          ! loop over subtracers:
          do itr = 1, self%ntr
            ! original:
            itracer = self%itracer(itr)
            ! loop over levels:
            do ilev = 1, self%nlev
              ! skip fixed layers:
              if ( self%fixed(ilev,itracer) == 1 ) cycle
              ! add contribution for this level/tracer:
              w(iw) = w(iw) + self%GXL(iakstar,ilev,itracer,iev,itime) * x_f(i,j_local,ilev,itr)
            end do ! levels
          end do ! tracers
        end do ! eigenvalues
      end do ! i
    end do ! j    
    ! check ...
    if ( iw /= self%nw_local(itime) ) then
      write (gol,'("only ",i0," of ",i0," elements of w filled")') iw, self%nw_local(itime); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! ~
    
    ! clear:
    deallocate( x_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( x_ex, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine BCovarSqrt_Forward


  ! ***
  

  !
  ! x =                   B^{1/2}                     w
  !
  !   = [ S E^T F^{-1} (Gamma^{1/2} X Lambda^{1/2}) ] w
  !
  ! With 'correlation=.true.' the std.dev. matrix "S" is ommitted.
  !
  
  subroutine BCovarSqrt_Reverse( self, itime, w, x, status, correlation )

    ! --- in/out ---------------------------------
    
    class(T_BCovarSqrt), intent(inout)          ::  self
    integer, intent(in)                         ::  itime
    complex, intent(in)                         ::  w(:)  ! (nw)
    real, intent(out)                           ::  x(:,:,:,:) ! (nlon,nlat,nlev,ntr)
    integer, intent(out)                        ::  status
    logical, intent(in), optional               ::  correlation

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/BCovarSqrt_Reverse'
    
    ! --- local ----------------------------------
    
    complex, allocatable  ::  x_f(:,:,:,:)  ! (nlon_f,nlat_f,nlev,ntr)
    real, allocatable     ::  x_ex(:,:,:,:)  ! (nlon_ex,nlat_ex,nlev,ntr)
    integer               ::  i, j
    integer               ::  j_local
    integer               ::  iakstar
    integer               ::  itracer, itr
    integer               ::  ilev
    integer               ::  iev
    integer               ::  iw
    logical               ::  corr
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (itime < 1) .or. (itime > self%ntime) ) then
      write (gol,'("time index ",i0," not in range 1 .. ",i0)') itime, self%ntime; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! check ...
    if ( self%nw_local(itime) > 0 ) then
      if ( size(w) /= self%nw_local(itime) ) then
        write (gol,'("size of 1D complex input (local) is ",i0," while expected ",i0)') &
                 size(w), self%nw_local(itime); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! check ...
    if ( self%nlat_local > 0 ) then
      if ( any( shape(x) /= (/self%nlon,self%nlat_local,self%nlev,self%ntr/) ) ) then
        write (gol,'("shape of output state is (",i0,3(",",i0),") while expected (",i0,3(",",i0),")")') &
                 shape(x), (/self%nlon,self%nlat_local,self%nlev,self%ntr/); call goErr
        TRACEBACK; status=1; return
      end if
    end if

    ! ~   (Gamma^{1/2} X Lambda^{1/2})   #   w
    !     (nlev*ntracer x nev)(k*(i,j))  #  (nev)(i,j)
    
    ! storage for spectral transform:
    allocate( x_f(self%nlon_f,self%nlat_f_local,self%nlev,self%ntr), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! init 1D index:
    iw = 0
    ! loop over spectral coeff:
    do j_local = 1, self%nlat_f_local
      j = self%ilat_f_offset + j_local
      do i = 1, self%nlon_f
        ! init sum:
        x_f(i,j_local,:,:) = (0.0,0.0)
        ! wavenumber band:
        iakstar = self%iakstar(i,j)
        ! skip if outside eliptic truncation:
        if ( iakstar < 1 ) cycle
        ! loop over eigenvalues:
        do iev = 1, self%nev(iakstar,itime)
          ! increase index:
          iw = iw + 1
          ! loop over subtracers:
          do itr = 1, self%ntr
            ! original:
            itracer = self%itracer(itr)
            ! loop over levels:
            do ilev = 1, self%nlev
              ! skip fixed layers:
              if ( self%fixed(ilev,itracer) == 1 ) cycle
              ! add contribution of eigenvalue to this level/tracer:
              x_f(i,j_local,ilev,itr) = x_f(i,j_local,ilev,itr) + &
                    self%GXL(iakstar,ilev,itracer,iev,itime) * w(iw)
            end do ! levels
          end do ! tracers
        end do ! eigenvalues
      end do ! i
    end do ! j    
    ! check ...
    if ( iw /= self%nw_local(itime) ) then
      write (gol,'("only ",i0," of ",i0," elements of w used")') iw, self%nw_local(itime); call goErr
      TRACEBACK; status=1; return
    end if

    ! ~ F^{-1} (Gamma^{1/2} X Lambda^{1/2} w)
    
    ! storage for extended grid:
    allocate( x_ex(self%nlon_ex,max(self%nlat_ex_local,1),self%nlev,self%ntr), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! inverse spectral transform;
    ! loop over trs:
    do itr = 1, self%ntr
      ! original:
      itracer = self%itracer(itr)
      ! loop over levels:
      do ilev = 1, self%nlev
        ! skip fixed layers:
        if ( self%fixed(ilev,itracer) == 1 ) then
          ! pad with zeros:
          x_ex(:,:,ilev,itr) = 0.0
        else
          ! inverse (half-complex)-to-real transform:
          call self%ufft%Inverse( x_f(:,:,ilev,itr), x_ex(:,:,ilev,itr), status )
          IF_NOT_OK_RETURN(status=1)
        end if
      end do  ! levels
    end do  ! tracers

    ! compensate for truncations;
    ! loop over subtracers:
    do itr = 1, self%ntr
      ! original:
      itracer = self%itracer(itr)
      ! loop over levels:
      do ilev = 1, self%nlev
        ! skip if fixed:
        if ( self%fixed(ilev,itracer) == 1 ) cycle
        ! scale:
        x_ex(:,:,ilev,itr) = x_ex(:,:,ilev,itr) * self%sqrt_phi(ilev,itracer,itime)
      end do ! levels
    end do ! tracers

    ! ~ S E^T (F^{-1} Gamma^{1/2} X Lambda^{1/2} w)
    
    ! flag:
    corr = .false.
    if ( present(correlation) ) corr = correlation
    
    ! local slab ?
    if ( self%nlat_local > 0 ) then
      ! copy first part of extended grid, and multiply with std.dev. if necessary:
      if ( corr ) then
        x =                         x_ex(1:self%nlon,1:self%nlat_local,:,:)
      else
        x = self%S(:,:,:,:,itime) * x_ex(1:self%nlon,1:self%nlat_local,:,:)
      endif
    end if
    
    ! ~
    
    ! clear:
    deallocate( x_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( x_ex, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine BCovarSqrt_Reverse


end module EMEP_BCovarSqrt

