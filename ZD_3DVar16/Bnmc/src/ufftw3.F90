!#######################################################################
!
! Tools arround FFTW3 library.
!
! Performs real-to-complex FFT or complex-to-real inverse FFT.
!
! In parallel, the FFTW3 library suggest an optimal decomposition 
! of the last dimension, thus y in 2D and z in 3D.
!
! Usage:
!
!   ! storage for dimensions, transformation plan, etc:
!   type(T_UFFTW3_2d)       ::  ufft
!   ! gridded fields:
!   real, allocatable       ::  xx(:,:)
!   real, allocatable       ::  yy(:,:)
!   ! spectral fields:
!   complex, allocatable    ::  ww(:,:)
!
!   ! initialize for specified grid size ;
!   ! eventually provide MPI communicator and logfile unit:
!   ufft%Init( 100, 160, status [,comm=MPI_COMM_WORLD] [,uout=6] )
!   if (status/=0) stop
!
!   ! grid size: nx=100, ny=160
!   call ufft%Get( status, nx=.., ny=.. )
!   if (status/=0) stop
!
!   ! spectral size in half-complex representation: nx=51, ny=160
!   call ufft%Get( status, nxf=.., nyf=.. )
!   if (status/=0) stop
!
!   ! spectral decomposition of last dimension;
!   ! on 4 pe: nyf_local=40/40/40/40, iyf_offset=0/40/80/120
!   call ufft%Get( status, nyf_local=.., iyf_offset=.. )
!   if (status/=0) stop
!
!   ! required grid decomposition of last dimension
!   ! (same as required spectral decomposition ):
!   ! on 4 pe: ny_local=40/40/40/40, iy_offset=0/40/80/120
!   call ufft%Get( status, ny_local=.., iy_offset=.. )
!   if (status/=0) stop
!
!   ! number of (local) real values in half-complex-to-real packing,
!   ! on 4 pe: nhcr=32000, nhcr_local=8000
!   call ufft%Get( status, nhcr=.., nhcr_local=.. )
!   if (status/=0) stop
!
!   ! storage:
!   allocate( xx(nx ,ny_local ) )
!   allocate( yy(nx ,ny_local ) )
!   allocate( ww(nxf,nyf_local) )
!
!   ! forward real-to-(half-)complex fft:
!   call ufft%Forward( xx, ww, status )
!   if (status/=0) stop
!
!   ! inverse (half-)complex-to-real fft:
!   call ufft%Inverse( ww, yy, status )
!   if (status/=0) stop
!
!   ! pack half-complex field into 1D array of real values,
!   ! return weights of dot product (see below):
!   call uuft%HC2HCR( ww, rr, status, l2w=.. )
!
!   ! clear:
!   calll ufft%Done( status )
!   if (status/=0) stop
!
!
! To pad or not to pad ...
! ------------------------
!
! Following the documentation, for real-to-complex tranforms
! with MPI it is necessary to use real input arrays that are
! about twice the logical size, with the remainder padded with
! zeros ("6.5 Multi-dimensional MPI DFTs of Real Data").
! For a serial run this does not seem to be true however;
! tests with an input field filled with diagonals showed that
! the returned complex field is not alligned as expected.
! Just allocating an input array with the logical shape
! solved this. The implementation has therefore different 
! storage for input and output depending on MPI enabled or not,
! until the documentation is better understood and a
! uniform code can be made.
!
!
! Half-complex representation and Hermitian producut
! --------------------------------------------------
!
! An FFT of a real field results in a "full" complex field
! where half of the elements are complex conjugates of the 
! other half. This leads to the property that the 
! Hermitian-product is equal to the dot product of
! vectors holding the real and imaginary parts:
!             _
!    x^H y =  x^T y 
!          = ..
!          = (Re[x],Im[x])^T (Re[y],Im[y])
!
! In half-complex storage only half of the complex number,
! since the other half consists of complex conjugates and
! are therefore implied by the first half.
! To be able to compute the Hermitian product as a dot product
! of real and imaginary parts, an element-wise weight of
! '1' and '2' values should be used:
!
!    x^H y =  l2w  o  (Re[x_hc],Im[x_hc])^T (Re[y_hc],Im[y_hc])
!
!             nhcr
!          =  sum  l2w(i) * x_hcr(i) * y_hcr(i)
!             i=1
!
!#######################################################################
!
#define TRACEBACK write (self%uout,'("ERROR in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__
!
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!#######################################################################

module UFFTW3

  use, intrinsic :: iso_c_binding
  
  implicit none
  

  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_UFFTW3_2d
!  public  ::  T_UFFTW3_3d
  

  ! --- const ------------------------------------
  
  character(len=*), parameter   ::  mname = 'UFFTW3_ml'
  
  ! gonio:
  real, parameter     ::  pi = 4.0 * atan(1.0)
  
  
  ! --- types ------------------------------------
  
  type T_UFFTW3_2d
    ! logfile unit:
    integer                         ::  uout
#ifdef _MPI
    ! physical problem and total size:
    integer(C_SIZE_T)               ::  nx, ny
    integer(C_SIZE_T)               ::  n
#else
    ! physical problem and total size:
    integer(C_INT)                  ::  nx, ny
    integer(C_INT)                  ::  n
#endif
    ! unitairy scale factor:
    real(C_DOUBLE)                  ::  ufactor
    ! spectral size:
    integer(C_INTPTR_T)             ::  nxf, nyf
    ! spectral indices:
    integer                         ::  Kx, Ky
    integer, allocatable            ::  xf(:)
    integer, allocatable            ::  yf(:), yf_local(:)
    ! wave numbers:
    integer                         ::  Ks
    real, allocatable               ::  kstar(:,:)
    real, allocatable               ::  kstar_local(:,:)
    character(len=256)              ::  kstar_formula
    ! spectral angle:
    real, allocatable               ::  theta(:,:)
    real, allocatable               ::  theta_local(:,:)
    ! number of represented complex numbers
    ! in half-complex storage:
    integer, allocatable            ::  hcn(:,:)
    integer, allocatable            ::  hcn_local(:,:)
    ! decomposition:
    integer(C_INTPTR_T)             ::  ny_local, iy_offset
    integer(C_INTPTR_T)             ::  nyf_local, iyf_offset
    ! half-complex-to-real storage:
    integer(C_INTPTR_T)             ::  nhcr, nhcr_local
    ! input and output arrays:
#ifdef _MPI
    integer(C_INTPTR_T)             ::  alloc_local
    real(C_DOUBLE), pointer         ::  input (:,:)
    real(C_DOUBLE), pointer         ::  input2(:,:)
    complex(C_DOUBLE), pointer      ::  output(:,:)
    type(C_PTR)                     ::  input_p
    type(C_PTR)                     ::  input2_p
    type(C_PTR)                     ::  output_p
#else
    real(C_DOUBLE), allocatable     ::  input (:,:)
    real(C_DOUBLE), allocatable     ::  input2(:,:)
    complex(C_DOUBLE), allocatable  ::  output(:,:)
#endif
    ! transformation plans:
    type(C_PTR)                     ::  plan_fwd
    type(C_PTR)                     ::  plan_inv
    !
  contains
    procedure   ::  Init         => UFFTW3_2d_Init
    procedure   ::  Done         => UFFTW3_2d_Done
    procedure   ::  Get          => UFFTW3_2d_Get
    procedure   ::  Forward      => UFFTW3_2d_Forward
    procedure   ::  Inverse      => UFFTW3_2d_Inverse
    procedure   ::  Truncate     => UFFTW3_2d_Truncate
    procedure   ::  TruncateCompensationFactor  => UFFTW3_2d_TruncateCompensationFactor
    procedure   ::  HC2HCR       => UFFTW3_2d_HC2HCR
    procedure   ::  HCR2HC       => UFFTW3_2d_HCR2HC
  end type T_UFFTW3_2d
  


contains


  ! ********************************************************************
  ! ***
  ! *** FFTW3 2D
  ! ***
  ! ********************************************************************


  subroutine UFFTW3_2d_Init( self, nx, ny, status, comm, uout, unitairy )

    use FFTW3, only : FFTW_Alloc_Real, FFTW_Alloc_Complex
    use FFTW3, only : FFTW_ESTIMATE, FFTW_MEASURE
#ifdef _MPI
    use FFTW3, only : FFTW_MPI_Local_Size_2d
    use FFTW3, only : FFTW_MPI_Plan_DFT_r2c_2d
    use FFTW3, only : FFTW_MPI_Plan_DFT_c2r_2d
#else
    use FFTW3, only : FFTW_Plan_DFT_r2c_2d
    use FFTW3, only : FFTW_Plan_DFT_c2r_2d
#endif
  
    ! --- in/out ---------------------------------
    
    class(T_UFFTW3_2d), intent(out)   ::  self
    integer, intent(in)               ::  nx, ny
    integer, intent(out)              ::  status
    integer, intent(in), optional     ::  comm
    integer, intent(in), optional     ::  uout
    logical, intent(in), optional     ::  unitairy

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_Init'
    
    ! --- local ----------------------------------

    integer(C_INT)              ::  fftw_flags
    integer                     ::  i, j, j_local
    logical                     ::  nxf_is_even
    integer                     ::  hcn
    
    ! --- begin ----------------------------------
    
    ! logfile unit:
    if ( present(uout) ) then 
      self%uout = uout
    else
      self%uout = 6  ! std.output
    end if
    
    ! save:
    self%nx = nx
    self%ny = ny
    
    ! total:
    self%n = self%nx * self%ny
    
    ! default scale factor:
    self%ufactor = 1.0
    ! reset ?
    if ( present(unitairy) .and. unitairy ) then
      ! scale factor for unitairy transformations:
      self%ufactor = 1.0 / sqrt( real(self%n) )
    end if
    
    ! spectral dimensions ; for real-to-complex in Fortran, 
    ! first dimension has roughly half of the complex numbers,
    ! since other half is implicitly equal to the conjugate:
    self%nxf = nx/2 + 1  ! ok for both odd and even
    self%nyf = ny
    
    ! how to optimize:
    !fftw_flags = FFTW_ESTIMATE
    fftw_flags = FFTW_MEASURE

#ifdef _MPI
    ! check communicator presence:
    if ( .not. present(comm) ) then
      write (self%uout,'("UFFTW3 compiled with _MPI macro defined, but no comm argument passed")')
      TRACEBACK; status=1; return
    end if

    ! data distribution is determined by shape of complex data:
    ! - provide spectral size in reversed order as 'integer(C_INTPTR_T)'
    ! - local total size is returned as 'integer(C_INTPTR_T)' too;
    ! - returns span and offset of last dimension as 'integer(C_INTPTR_T)'
    ! note reverse order of dimensions:
    self%alloc_local = FFTW_MPI_Local_Size_2d( self%nyf, self%nxf, comm, &
                                            self%nyf_local, self%iyf_offset )

    ! grid point decomposition is the same as spectral ...
    self%ny_local  = self%nyf_local
    self%iy_offset = self%iyf_offset
    
    ! real input requires array twice as large as the complex output (in elements);
    ! number of real values is therefore slightly larger than original field
    ! (FFTW3 doc, section 6.5 "Multi-dimensional MPI DFTs of Real Data")
    self%input_p = FFTW_Alloc_Real( 2*self%alloc_local )
    ! allow reference as Fortran array:
    call C_F_Pointer( self%input_p, self%input, (/2*self%nxf,self%nyf_local/) )
    
    ! other array of same size for inversion result:
    self%input2_p = FFTW_Alloc_Real( 2*self%alloc_local )
    ! allow reference as Fortran array:
    call C_F_Pointer( self%input2_p, self%input2, (/2*self%nxf,self%nyf_local/) )
    
    ! complex output field:
    self%output_p = FFTW_Alloc_Complex( self%alloc_local )
    ! allow reference as Fortran array:
    call C_F_Pointer( self%output_p, self%output, (/self%nxf,self%nyf_local/) )

    ! setup real-to-complex (half-complex) transphorm;
    ! dimension orders need to be reversed !
    self%plan_fwd = FFTW_MPI_Plan_DFT_r2c_2d( self%ny, self%nx, &
                                self%input, self%output, &
                                comm, fftw_flags )

    ! setup inverse:
    self%plan_inv = FFTW_MPI_Plan_DFT_c2r_2d( self%ny, self%nx, &
                                self%output, self%input, &
                                comm, fftw_flags )
#else

    ! no decomposition, single slab:
    self%nyf_local  = self%nyf
    self%iyf_offset = 0

    ! grid point decomposition is the same as spectral,
    ! here no decomposition actually:
    self%ny_local  = self%nyf_local
    self%iy_offset = self%iyf_offset
    
    ! storage for input (gridded):
    allocate( self%input(self%nx,self%ny_local), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for inversion result (gridded):
    allocate( self%input2(self%nx,self%ny_local), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage for output (spectral):
    allocate( self%output(self%nxf,self%nyf_local), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! setup real-to-complex (half-complex) transphorm;
    ! dimension orders need to be reversed !
    self%plan_fwd = FFTW_Plan_DFT_r2c_2d( self%ny, self%nx, &
                                self%input, self%output, &
                                fftw_flags )

    ! setup inverse:
    self%plan_inv = FFTW_Plan_DFT_c2r_2d( self%ny, self%nx, &
                                self%output, self%input2, &
                                fftw_flags )
#endif

    ! size of half-complex to real storage ;
    ! twice the number of complex values to store 
    ! individual real and imaginary values:
    self%nhcr       = 2 * self%nxf * self%nyf
    self%nhcr_local = 2 * self%nxf * self%nyf_local

    ! spectral indices in x-direction ; in full storage:
    !    -nx/2-1, ..., -2, -1, 0, 1, 2, ..., nx/2-1 [,nx/2]
    ! in half complex storage only the positive half:
    !                 0, 1, 2, ...
    ! storage:
    allocate( self%xf(self%nxf), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do i = 1, self%nxf
      self%xf(i) = i-1
    end do
    ! maximum:
    self%Kx = maxval( self%xf )
    
    ! spectral indices in y-direction in original storage:
    !    -ny/2-1, ..., -2, -1, 0, 1, 2, ..., ny/2-1 [,ny/2]
    ! here order is changed:
    !     0, 1, 2, ..., ny/2-1 [,ny/2], -ny/2-1, ..., -2, -1
    ! storage:
    allocate( self%yf(self%nyf), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do j = 1, self%nyf
      if ( j <= self%nyf/2+1 ) then
        self%yf(j) = j-1
      else
        self%yf(j) = -self%nyf-1 + j
      end if
    end do
    ! maximum:
    self%Ky = maxval( self%yf )
    ! storage for local slab:
    allocate( self%yf_local(self%nyf_local), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do j_local = 1, self%nyf_local
      j = self%iyf_offset + j_local
      if ( j <= self%nyf/2+1 ) then
        self%yf_local(j_local) = j-1
      else
        self%yf_local(j_local) = -self%nyf-1 + j
      end if
    end do
    
    ! Wave vector:
    !   _
    !   k* = Ns ( m/M, n/N )
    !
    ! where Ns is arbitrary scale factor, here:
    !   Ns = max( M, N )
    !
    ! Wave number is norm:
    !   k* = Ns sqrt( (m/M)**2 + (n/N)**2 )
    !
    ! scale:
    self%Ks = max( self%Kx, self%Ky )
    ! storage:
    allocate( self%kstar(self%nxf,self%nyf), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%kstar_local(self%nxf,self%nyf_local), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do j = 1, self%nyf
      do i = 1, self%nxf
        ! norm:
        self%kstar(i,j) = self%Ks * sqrt( (real(self%xf(i))/self%Kx)**2 + &
                                          (real(self%yf(j))/self%Ky)**2     )
        ! copy:
        j_local = j - self%iyf_offset
        if ( (0 < j_local) .and. (j_local <= self%nyf_local) ) then
          self%kstar_local(i,j_local) = self%kstar(i,j)
        end if
      end do ! i
    end do ! j
    ! store info:
    write (self%kstar_formula,'("kstar = Ks sqrt( (m/Kx)^2 + (n/Ky)^2 ) ; Kx = ",i0,", Ky = ",i0," ; Ks = max(Kx,Ky) = ",i0)') &
                         self%Kx, self%Ky, self%Ks
    
    
    ! spectral angle:
    !   
    !   tan(theta) = (n/N) / (m/M)
    !
    ! storage:
    allocate( self%theta(self%nxf,self%nyf), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%theta_local(self%nxf,self%nyf_local), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do j = 1, self%nyf
      do i = 1, self%nxf
        ! trap divison by zero:
        if ( self%xf(i) == 0 ) then
          ! on m==0, set to +/- 90 degrees:
          if ( self%yf(j) < 0 ) then
            self%theta(i,j) = - 0.5 * pi
          else if ( self%yf(j) > 0 ) then
            self%theta(i,j) = 0.5 * pi
          else
            self%theta(i,j) = 0.0
          end if
        else
          ! angle:
          self%theta(i,j) = atan( (real(self%yf(j))/self%Ky) / (real(self%xf(i))/self%Kx) )
        end if
        ! copy:
        j_local = j - self%iyf_offset
        if ( (0 < j_local) .and. (j_local <= self%nyf_local) ) then
          self%theta_local(i,j_local) = self%theta(i,j)
        end if
      end do ! i
    end do ! j
    !
    ! numbe of represented values in half-complex storage
    !
    ! storage:
    allocate( self%hcn(self%nxf,self%nyf), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%hcn_local(self%nxf,self%nyf_local), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! half-complex definition is used in x-dimension ;
    ! check if this is odd or even:
    nxf_is_even = mod(self%nxf,2) == 0
    ! fill:
    do j = 1, self%nyf
      do i = 1, self%nxf
        ! weight is 1 for original real values:
        if ( (i == 1) .or. ((i==self%nxf) .and. nxf_is_even) ) then
          self%hcn(i,j) = 1
        else
          self%hcn(i,j) = 2
        end if
        ! copy to local array:
        j_local = j - self%iyf_offset
        if ( (0 < j_local) .and. (j_local <= self%nyf_local) ) then
          self%hcn_local(i,j_local) = self%hcn(i,j)
        end if
      end do ! i
    end do ! j

    ! ok
    status = 0

  end subroutine UFFTW3_2d_Init


  ! ***
  
  
  subroutine UFFTW3_2d_Done( self, status )

    use FFTW3, only : FFTW_Free
    use FFTW3, only : FFTW_Destroy_Plan
    
    ! --- in/out ---------------------------------
    
    class(T_UFFTW3_2d), intent(inout)       ::  self
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! done:
    call FFTW_Destroy_Plan( self%plan_fwd )
    call FFTW_Destroy_Plan( self%plan_inv )

#ifdef _MPI
    ! release:
    call FFTW_Free( self%input_p )
    call FFTW_Free( self%input2_p )
    call FFTW_Free( self%output_p )
    ! reset:
    nullify( self%input )
    nullify( self%input2 )
    nullify( self%output )
#else
    ! clear:
    deallocate( self%input, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%input2, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%output, stat=status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! clear:
    deallocate( self%xf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%yf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%yf_local, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%kstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%kstar_local, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%theta, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%theta_local, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%hcn, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%hcn_local, stat=status )
    IF_NOT_OK_RETURN(status=1)
  
    ! dummy:
    self%nx = -999
    self%ny = -999
    
    ! ok
    status = 0

  end subroutine UFFTW3_2d_Done


  ! ***


  subroutine UFFTW3_2d_Get( self, status, &
                              nx, ny, ny_local, iy_offset, &
                              nxf, nyf, nyf_local, iyf_offset, &
                              nhcr, nhcr_local, &
                              xf, yf, kstar, kstar_formula, theta, &
                              hcn, hcn_local )

    ! --- in/out ---------------------------------

    class(T_UFFTW3_2d), intent(inout)       ::  self
    integer, intent(out)                    ::  status
    integer, intent(out), optional          ::  nx, ny
    integer, intent(out), optional          ::  ny_local, iy_offset
    integer, intent(out), optional          ::  nxf, nyf
    integer, intent(out), optional          ::  nyf_local, iyf_offset
    integer, intent(out), optional          ::  nhcr, nhcr_local
    integer, intent(out), optional          ::  xf(:), yf(:)
    real, intent(out), optional             ::  kstar(:,:)
    character(len=*), intent(out), optional ::  kstar_formula
    real, intent(out), optional             ::  theta(:,:)
    integer, intent(out), optional          ::  hcn(:,:)
    integer, intent(out), optional          ::  hcn_local(:,:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_Get'

    ! --- local ----------------------------------
    
    integer         ::  i

    ! --- begin ----------------------------------

    ! return values:
    if ( present(nx        ) ) nx         = self%nx
    if ( present(ny        ) ) ny         = self%ny
    if ( present(ny_local  ) ) ny_local   = self%ny_local
    if ( present(iy_offset ) ) iy_offset  = self%iy_offset
    if ( present(nxf       ) ) nxf        = self%nxf
    if ( present(nyf       ) ) nyf        = self%nyf
    if ( present(nyf_local ) ) nyf_local  = self%nyf_local
    if ( present(iyf_offset) ) iyf_offset = self%iyf_offset
    if ( present(nhcr      ) ) nhcr       = self%nhcr
    if ( present(nhcr_local) ) nhcr_local = self%nhcr_local
    
    ! spectral indices:
    if ( present(xf) ) then
      ! check ..
      if ( size(xf) /= self%nxf ) then
        write (self%uout,'("ERROR - size of xf does not match:")')
        write (self%uout,'("ERROR -   argument   : ",i6)') size(xf)
        write (self%uout,'("ERROR -   expected   : ",i6)') size(self%xf)
        TRACEBACK; status=1; return
      end if
      ! copy:
      xf = self%xf
    end if
    
    ! spectral indices:
    if ( present(yf) ) then
      ! check ..
      if ( size(yf) /= self%nyf ) then
        write (self%uout,'("ERROR - size of yf does not match:")')
        write (self%uout,'("ERROR -   argument   : ",i6)') size(yf)
        write (self%uout,'("ERROR -   expected   : ",i6)') size(self%yf)
        TRACEBACK; status=1; return
      end if
      ! copy:
      yf = self%yf
    end if
    
    ! wave number:
    if ( present(kstar) ) then
      ! check ..
      if ( (size(kstar,1) /= self%nxf) .or. (size(kstar,2) /= self%nyf) ) then
        write (self%uout,'("ERROR - size of kstar does not match:")')
        write (self%uout,'("ERROR -   argument   : ",2i6)') shape(yf)
        write (self%uout,'("ERROR -   expected   : ",2i6)') self%nxf, self%nyf
        TRACEBACK; status=1; return
      end if
      ! copy:
      kstar = self%kstar
    end if
    ! info:
    if ( present(kstar_formula) ) kstar_formula = trim(self%kstar_formula)
    
    ! spectral angle:
    if ( present(theta) ) then
      ! check ..
      if ( (size(theta,1) /= self%nxf) .or. (size(theta,2) /= self%nyf) ) then
        write (self%uout,'("ERROR - size of theta does not match:")')
        write (self%uout,'("ERROR -   argument   : ",2i6)') shape(yf)
        write (self%uout,'("ERROR -   expected   : ",2i6)') self%nxf, self%nyf
        TRACEBACK; status=1; return
      end if
      ! copy:
      theta = self%theta
    end if
    
    ! number of represented complex numbers in half-complex storage
    if ( present(hcn) ) then
      ! check ..
      if ( any( shape(hcn) /= (/self%nxf,self%nyf/) ) ) then
        write (self%uout,'("ERROR - size of hcn does not match:")')
        write (self%uout,'("ERROR -   argument   : ",2i6)') shape(hcn)
        write (self%uout,'("ERROR -   expected   : ",2i6)') self%nxf, self%nyf
        TRACEBACK; status=1; return
      end if
      ! copy:
      hcn = self%hcn
    end if
    
    ! number of represented complex numbers in half-complex storage
    if ( present(hcn_local) ) then
      ! check ..
      if ( any( shape(hcn_local) /= (/self%nxf,self%nyf_local/) ) ) then
        write (self%uout,'("ERROR - size of hcn_local does not match:")')
        write (self%uout,'("ERROR -   argument   : ",2i6)') shape(hcn_local)
        write (self%uout,'("ERROR -   expected   : ",2i6)') self%nxf, self%nyf_local
        TRACEBACK; status=1; return
      end if
      ! copy:
      hcn_local = self%hcn_local
    end if
    
    ! ok
    status = 0

  end subroutine UFFTW3_2d_Get


  ! ***
  
  
  subroutine UFFTW3_2d_Forward( self, input, output, status )

#ifdef _MPI
    use FFTW3, only : FFTW_MPI_Execute_DFT_r2c
#else
    use FFTW3, only : FFTW_Execute_DFT_r2c
#endif
    
    ! --- in/out ---------------------------------
    
    class(T_UFFTW3_2d), intent(inout)       ::  self
    real(C_DOUBLE), intent(in)              ::  input(:,:)
    complex(C_DOUBLE), intent(out)          ::  output(:,:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_Forward'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! local gridpoint slab defined ?
    if ( self%ny_local > 0 ) then
      ! check size ..
      if ( any( shape(input) /= (/int(self%nx),int(self%ny_local)/) ) ) then
        write (self%uout,'("ERROR - real input sizes do not match:")')
        write (self%uout,'("ERROR -   argument : ",2i6)') shape(input)
        write (self%uout,'("ERROR -   intern   : ",2i6)') self%nx, self%ny_local
        TRACEBACK; status=1; return
      end if
      ! fill input, pad extra elements with zeros:
      self%input = 0.0
      self%input(1:self%nx,1:self%ny_local) = input
    end if

    ! transform:
#ifdef _MPI
    call FFTW_MPI_Execute_DFT_r2c( self%plan_fwd, self%input, self%output )
#else
    call FFTW_Execute_DFT_r2c( self%plan_fwd, self%input, self%output )
#endif

    ! local spectral slab defined ?
    if ( self%nyf_local > 0 ) then
      ! check size ..
      if ( any( shape(output) /= shape(self%output) ) ) then
        write (self%uout,'("ERROR - complex output sizes do not match:")')
        write (self%uout,'("ERROR -   argument : ",2i6)') shape(output)
        write (self%uout,'("ERROR -   intern   : ",2i6)') shape(self%output)
        TRACEBACK; status=1; return
      end if
      ! copy result, scale:
      output = self%output * self%ufactor
    end if

    ! ok
    status = 0

  end subroutine UFFTW3_2d_Forward


  ! ***
  
  
  subroutine UFFTW3_2d_Inverse( self, output, input, status )

#ifdef _MPI
    use FFTW3, only : FFTW_MPI_Execute_DFT_c2r
#else
    use FFTW3, only : FFTW_Execute_DFT_c2r
#endif
    
    ! --- in/out ---------------------------------
    
    class(T_UFFTW3_2d), intent(inout)       ::  self
    complex(C_DOUBLE), intent(in)           ::  output(:,:)
    real(C_DOUBLE), intent(out)             ::  input(:,:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_Inverse'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! local spectral slab defined ?
    if ( self%nyf_local > 0 ) then
      ! check size ..
      if ( any( shape(output) /= shape(self%output) ) ) then
        write (self%uout,'("ERROR - complex input sizes do not match:")')
        write (self%uout,'("ERROR -   argument : ",2i6)') shape(output)
        write (self%uout,'("ERROR -   intern   : ",2i6)') shape(self%output)
        TRACEBACK; status=1; return
      end if
      ! copy into alligned memory:
      self%output = output
    end if
    
    ! inverse transform:
#ifdef _MPI
    call FFTW_MPI_Execute_DFT_c2r( self%plan_inv, self%output, self%input2 )
#else
    call FFTW_Execute_DFT_c2r( self%plan_inv, self%output, self%input2 )
#endif
    
    ! local gridpoint slab defined ?
    if ( self%ny_local > 0 ) then
      ! check size ..
      if ( any( shape(input) /= (/int(self%nx),int(self%ny_local)/) ) ) then
        write (self%uout,'("ERROR - real output do not match:")')
        write (self%uout,'("ERROR -   argument : ",2i6)') shape(input)
        write (self%uout,'("ERROR -   intern   : ",2i6)') self%nx, self%ny_local
        TRACEBACK; status=1; return
      end if
      ! extract from padded array, scale result:
      input = self%input2(1:self%nx,1:self%ny_local) * self%ufactor
    end if

    ! ok
    status = 0

  end subroutine UFFTW3_2d_Inverse


  ! ***
  
  
  !
  ! Eliptic truncation of spectral transform:
  !   (m,n) := 0.0  where  (m/M)**2 + (n/N)**2 <= 1
  ! where:
  !    m  in  -M, .., -2, -1, 0, 1, 2, .., M
  !    n  in  -N, .., -2, -1, 0, 1, 2, .., N
  !
  ! Instead of 0.0, truncated values are set to the fill_value if present.
  !
  
  subroutine UFFTW3_2d_Truncate( self, output, status, fill_value )

    ! --- in/out ---------------------------------
    
    class(T_UFFTW3_2d), intent(in)          ::  self
    complex(C_DOUBLE), intent(inout)        ::  output(:,:)
    integer, intent(out)                    ::  status
    
    real, intent(in), optional              ::  fill_value

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_Truncate'
    
    ! --- local ----------------------------------
    
    integer   ::  ixf, iyf_local
    real      ::  nodata_value
    
    ! --- begin ----------------------------------
    
    ! local spectral slab defined ?
    if ( self%nyf_local > 0 ) then
      ! check size ..
      if ( any( shape(output) /= shape(self%output) ) ) then
        write (self%uout,'("ERROR - complex output sizes do not match:")')
        write (self%uout,'("ERROR -   argument : ",2i6)') shape(output)
        write (self%uout,'("ERROR -   intern   : ",2i6)') shape(self%output)
        TRACEBACK; status=1; return
      end if
      ! set the no-data value:
      nodata_value = 0.0
      if ( present(fill_value) ) nodata_value = fill_value
      ! loop:
      do ixf = 1, self%nxf
        do iyf_local = 1, self%nyf_local
          if ( self%kstar_local(ixf,iyf_local) > self%Ks ) then
            output(ixf,iyf_local) = nodata_value
          end if
        end do  ! i
      end do ! j_local
    end if

    ! ok
    status = 0

  end subroutine UFFTW3_2d_Truncate


  ! *
  
  !
  ! Compensation factor (>= 1) for loss due to eliptic truncation
  !
  
  subroutine UFFTW3_2d_TruncateCompensationFactor( self, factor, status )

    ! --- in/out ---------------------------------
    
    class(T_UFFTW3_2d), intent(inout)       ::  self
    real, intent(out)                       ::  factor
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_TruncateCompensationFactor'
    
    ! --- local ----------------------------------
    
    real(C_DOUBLE), allocatable         ::  x(:,:)
    complex(C_DOUBLE), allocatable      ::  w(:,:)
    
    ! --- begin ----------------------------------
    
    ! local gridpoint slab defined ?
    if ( self%ny_local > 0 ) then
      ! storage for input:
      allocate( x(self%nx,self%ny_local), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! dummy:
      allocate( x(self%nx,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! local spectral slab defined ?
    if ( self%nyf_local > 0 ) then
      ! storage for output:
      allocate( w(self%nxf,self%nyf_local), stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      ! dummy:
      allocate( w(self%nxf,1), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! fill input with unity vector:
    x = 0.0
    if ( self%iy_offset+1 == 1 ) x(1,self%iy_offset+1) = 1.0
    
    ! forward transform:
    call self%Forward( x, w, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! eliptic truncation:
    call self%Truncate( w, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! backward transform:
    call self%Inverse( w, x, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! factor to compensate for peak smaller than 1.0 :
    factor = 1.0 / maxval(x)
    
    ! clear:
    deallocate( x, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( w, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine UFFTW3_2d_TruncateCompensationFactor
  


  ! ***
  
  
  ! Pack half-complex field into 1D array of real values.
  !
  ! Eventually return weights in dot product to maintain
  ! the property that for this type of complex fields
  ! (with half the values the complex conjugate of the other)
  ! the hermitian-product is equal to the dot product of
  ! vectors holding the real and imaginary parts:
  !
  !    x^H y =  l2w  o  (Re[x_hc],Im[x_hc])^T (Re[y_hc],Im[y_hc])
  !
  !             nhcr
  !          =  sum  l2w(i) * x_hcr(i) * y_hcr(i)
  !             i=1
  !
  
  subroutine UFFTW3_2d_HC2HCR( self, hc, hcr, status, l2w )

    ! --- in/out ---------------------------------
    
    class(T_UFFTW3_2d), intent(in)          ::  self
    complex(C_DOUBLE), intent(in)           ::  hc(:,:)
    real(C_DOUBLE), intent(out)             ::  hcr(:)
    integer, intent(out)                    ::  status
    integer, intent(out), optional          ::  l2w(:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_HC2HCR'
    
    ! --- local ----------------------------------
    
    integer     ::  i, j
    integer     ::  i0
    
    ! --- begin ----------------------------------
    
    ! local spectral slab defined ?
    if ( self%nyf_local > 0 ) then

      ! check size ..
      if ( any( shape(hc) /= shape(self%output) ) ) then
        write (self%uout,'("ERROR - complex input sizes do not match:")')
        write (self%uout,'("ERROR -   argument : ",2i6)') shape(hc)
        write (self%uout,'("ERROR -   intern   : ",2i6)') shape(self%output)
        TRACEBACK; status=1; return
      end if

      ! check size ..
      if ( size(hcr) /= self%nhcr_local ) then
        write (self%uout,'("ERROR - real output size does not match:")')
        write (self%uout,'("ERROR -   argument : ",i12)') size(hcr)
        write (self%uout,'("ERROR -   expected : ",i12)') self%nhcr
        TRACEBACK; status=1; return
      end if

      ! check size ..
      if ( present(l2w) ) then
        if ( size(l2w) /= self%nhcr_local ) then
          write (self%uout,'("ERROR - l2w output size does not match:")')
          write (self%uout,'("ERROR -   argument : ",i12)') size(l2w)
          write (self%uout,'("ERROR -   expected : ",i12)') self%nhcr
          TRACEBACK; status=1; return
        end if
      end if

      ! init offset in 1D array:
      i0 = 0
      ! loop over local cells:
      do j = 1, self%nyf_local
        do i = 1, self%nxf
          ! copy from 2D complex field into real array:
          hcr(i0+1) =  real( hc(i,j) )
          hcr(i0+2) = aimag( hc(i,j) )

          ! set weight in l2 norm ?
          if ( present(l2w) ) then
            ! copy from number of represented values:
            l2w(i0+1:i0+2) = self%hcn_local(i,j)
          end if

          ! next pair of elements:
          i0 = i0 + 2

        end do  ! i
      end do  ! j
      
    end if  ! nyf_local > 0

    ! ok
    status = 0

  end subroutine UFFTW3_2d_HC2HCR
  
  ! *
  
  subroutine UFFTW3_2d_HCR2HC( self, hcr, hc, status )

    ! --- in/out ---------------------------------
    
    class(T_UFFTW3_2d), intent(in)          ::  self
    real(C_DOUBLE), intent(in)              ::  hcr(:)
    complex(C_DOUBLE), intent(out)          ::  hc(:,:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/UFFTW3_2d_HC2HCR'
    
    ! --- local ----------------------------------
    
    integer     ::  i, j
    integer     ::  i0
    
    ! --- begin ----------------------------------
    
    ! local spectral slab defined ?
    if ( self%nyf_local > 0 ) then

      ! check size ..
      if ( size(hcr) /= self%nhcr_local ) then
        write (self%uout,'("ERROR - real input size does not match:")')
        write (self%uout,'("ERROR -   argument : ",i12)') size(hcr)
        write (self%uout,'("ERROR -   expected : ",i12)') self%nhcr_local
        TRACEBACK; status=1; return
      end if

      ! check size ..
      if ( any( shape(hc) /= shape(self%output) ) ) then
        write (self%uout,'("ERROR - complex output sizes do not match:")')
        write (self%uout,'("ERROR -   argument : ",2i6)') shape(hc)
        write (self%uout,'("ERROR -   intern   : ",2i6)') shape(self%output)
        TRACEBACK; status=1; return
      end if

      ! init offset in 1D array:
      i0 = 0
      ! loop over local cells:
      do j = 1, self%nyf_local
        do i = 1, self%nxf
          ! copy values from 1D real state to 3D complex array:
          hc(i,j) = cmplx( hcr(i0+1), hcr(i0+2) )
          ! next pair of elements:
          i0 = i0 + 2
        end do  ! i
      end do  ! j
      
    end if  ! nyf_local > 0

    ! ok
    status = 0

  end subroutine UFFTW3_2d_HCR2HC


end module UFFTW3

