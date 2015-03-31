!#######################################################################
! Fortran module for FFTW,
! following instructions in documentation:
!
!   FFTW User Manual (v3.3.4)
!   7 Calling FFTW from Modern Fortran
!   7.7 Defining an FFTW module
!
! http://www.fftw.org/fftw3_doc/Defining-an-FFTW-module.html
!#######################################################################
!
! Bindings and Tools arround FFTW3 library.
! @author A Segers
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
!   call ufft%Forward( ww, yy, status )
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
! Half-complex representation and Hermitian producut
! --------------------------------------------------
!
! An FFT of a real field results in a "full" complex field
! where half of the elements are complex conjugates of the 
! other half. This leads to the property that the 
! Hermitian-product is equal to the dot product of
! vectors holding the real and imaginary parts:
!
!    x^H y =  (Re[x],Im[x])^T (Re[y],Im[y])
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
#define STRING2(x) #x
#define STRING(x) STRING2(x)
#define HERE(MSG) MSG//" ("//__FILE__//":"//STRING(__LINE__)//")."
module FFTW3_MPI
use, intrinsic :: iso_c_binding
use mpi, only : MPI_COMM_WORLD, MPI_SUM, MPI_DOUBLE_PRECISION !,MPI_ALLREDUCE
implicit none
private
public :: &
  FFTW_Alloc_Real, FFTW_Alloc_Complex, &
  FFTW_Free, FFTW_Destroy_Plan, FFTW_FLAGS
public :: &
  FFTW_MPI_Local_Size_2d,   FFTW_MPI_Local_Size_3d, &
  FFTW_MPI_Plan_DFT_r2c_2d, FFTW_MPI_Plan_DFT_c2r_2d, &
  FFTW_MPI_Plan_DFT_r2c_3d, FFTW_MPI_Plan_DFT_c2r_3d, &
  FFTW_MPI_Execute_DFT_r2c, FFTW_MPI_Execute_DFT_c2r
include 'fftw3-mpi.f03'
integer, parameter :: FFTW_FLAGS=FFTW_ESTIMATE  ! how to optimize
!nteger, parameter :: FFTW_FLAGS=FFTW_MEASURE   ! how to optimize
endmodule FFTW3_MPI
! ********************************************************************
! ***
! *** FFTW3 2D
! ***
! ********************************************************************
module UFFTW3_2d
use, intrinsic :: iso_c_binding
use CheckStop_ml, only: CheckStop
use FFTW3_MPI, only: &
  FFTW_Alloc_Real, FFTW_Alloc_Complex, &
  FFTW_Free, FFTW_Destroy_Plan, FFTW_FLAGS, &
  FFTW_Local_Size   => FFTW_MPI_Local_Size_2d,  &
  FFTW_Plan_DFT_r2c => FFTW_MPI_Plan_DFT_r2c_2d,&
  FFTW_Exec_DFT_r2c => FFTW_MPI_Execute_DFT_r2c,&
  FFTW_Plan_DFT_c2r => FFTW_MPI_Plan_DFT_c2r_2d,&
  FFTW_Exec_DFT_c2r => FFTW_MPI_Execute_DFT_c2r
implicit none
! --- in/out -----------------------------------
private
public :: UFFTW3_T
! --- const ------------------------------------ 
character(len=*), parameter   ::  mname = 'UFFTW3_2d'
! --- types ------------------------------------ 
type UFFTW3_T
  integer(C_SIZE_T)           ::  nx, ny, n ! physical problem and total size
  real(C_DOUBLE)              ::  ufactor   ! unitairy scale factor
  integer(C_INTPTR_T)         ::  nxf, nyf  ! spectral size
  ! decomposition:
  integer(C_INTPTR_T)         ::  alloc_local, &
                                  ny_local,  iy_offset, &
                                  nyf_local, iyf_offset
  ! half-complex-to-real storage:
  integer(C_INTPTR_T)         ::  nhcr, nhcr_local
  ! input and output arrays:
  real(C_DOUBLE), pointer     ::  input(:,:),input2(:,:)
  complex(C_DOUBLE), pointer  ::  output(:,:)
  type(C_PTR)                 ::  input_p,input2_p,output_p
  ! transformation plans:
  type(C_PTR)                 ::  plan_fwd, plan_inv
contains
  procedure :: Init    => Init
  procedure :: Done    => Done
  procedure :: Get     => Get
  procedure :: Forward => Forward
  procedure :: Inverse => Inverse
  procedure :: HC2HCR  => HC2HCR
  procedure :: HCR2HC  => HCR2HC
endtype UFFTW3_T
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
subroutine Init(self, nx, ny, comm)
! --- in/out --------------------------------- 
  class(UFFTW3_T), intent(out)  ::  self
  integer, intent(in)           ::  nx, ny
  integer, intent(in)           ::  comm
! --- const --------------------------------------
  character(len=*), parameter ::  rname = mname//'%Init'
! --- begin ----------------------------------
    
  ! save:
  self%nx = nx
  self%ny = ny
    
  ! total:
  self%n = self%nx * self%ny
    
  !! scale factor for unitairy transformations:
  !self%ufactor = 1.0 / sqrt( real(self%n) )
  ! not yet ...
  self%ufactor = 1.0 
    
  ! spectral dimensions ; for real-to-complex in Fortran, 
  ! first dimension has roughly half of the complex numbers,
  ! since other half is implicitly equal to the conjugate:
  self%nxf = nx/2 + 1  ! ok for both odd and even
  self%nyf = ny
    
  ! data distribution is determined by shape of complex data:
  ! - provide spectral size in reversed order as 'integer(C_INTPTR_T)'
  ! - local total size is returned as 'integer(C_INTPTR_T)' too;
  ! - returns span and offset of last dimension as 'integer(C_INTPTR_T)'
  ! note reverse order of dimensions:
  self%alloc_local = FFTW_Local_Size(self%nyf, self%nxf, comm, &
                                     self%nyf_local, self%iyf_offset)

  ! grid point decomposition is the same as spectral ...
  self%ny_local  = self%nyf_local
  self%iy_offset = self%iyf_offset
    
  ! real input requires array twice as large as the complex output (in elements);
  ! number of real values is therefore slightly larger than original field
  ! (FFTW3 doc, section 6.5 "Multi-dimensional MPI DFTs of Real Data")
  self%input_p = FFTW_Alloc_Real(2*self%alloc_local)
  ! allow reference as Fortran array:
  call C_F_Pointer(self%input_p, self%input, [2*self%nxf,self%nyf_local])
    
  ! extra:
  self%input2_p = FFTW_Alloc_Real(2*self%alloc_local)
  ! allow reference as Fortran array:
  call C_F_Pointer(self%input2_p,self%input2,[2*self%nxf,self%nyf_local])

  ! complex output field:
  self%output_p = FFTW_Alloc_Complex(self%alloc_local)
  ! allow reference as Fortran array:
  call C_F_Pointer(self%output_p, self%output, [self%nxf,self%nyf_local])

  ! setup real-to-complex (half-complex) transphorm;
  ! dimension orders need to be reversed !
  self%plan_fwd = FFTW_Plan_DFT_r2c(self%ny, self%nx, &
                              self%input, self%output, comm, FFTW_FLAGS)

  ! setup inverse:
  self%plan_inv = FFTW_Plan_DFT_c2r(self%ny, self%nx, &
                              self%output, self%input2, comm, FFTW_FLAGS)

  ! size of half-complex to real storage ;
  ! twice the number of complex values to store 
  ! individual real and imaginary values:
  self%nhcr       = 2 * self%nxf * self%nyf
  self%nhcr_local = 2 * self%nxf * self%nyf_local
endsubroutine Init
!-----------------------------------------------------------------------
subroutine Done(self)
  ! --- in/out ---------------------------------
  class(UFFTW3_T), intent(inout)  :: self
  ! --- const --------------------------------------
  character(len=*), parameter     ::  rname = mname//'%Done'
  ! --- begin ----------------------------------
    
  ! done:
  call FFTW_Destroy_Plan(self%plan_fwd)
  call FFTW_Destroy_Plan(self%plan_inv)

  ! clear:
  call FFTW_Free(self%input_p )
  call FFTW_Free(self%input2_p)
  call FFTW_Free(self%output_p)

  ! dummy:
  self%nx = -999
  self%ny = -999
endsubroutine Done
!-----------------------------------------------------------------------
subroutine Get(self, &
               nx, ny, ny_local, iy_offset, &
               nxf, nyf, nyf_local, iyf_offset, &
               nhcr, nhcr_local, l2w_r, l2w_c)
  ! --- in/out ---------------------------------
  class(UFFTW3_T), intent(inout)  :: self
  integer, intent(out), optional  :: nx, ny, &
                                     ny_local, iy_offset, &
                                     nxf, nyf, &
                                     nyf_local, iyf_offset, &
                                     nhcr, nhcr_local
  integer, intent(out), pointer, optional  :: l2w_r(:),l2w_c(:,:)
  ! --- const --------------------------------------
  character(len=*), parameter     ::  rname = mname//'%Get'
  ! --- begin ----------------------------------

  ! return values:
  if(present(nx        )) nx         = self%nx
  if(present(ny        )) ny         = self%ny
  if(present(ny_local  )) ny_local   = self%ny_local
  if(present(iy_offset )) iy_offset  = self%iy_offset
  if(present(nxf       )) nxf        = self%nxf
  if(present(nyf       )) nyf        = self%nyf
  if(present(nyf_local )) nyf_local  = self%nyf_local
  if(present(iyf_offset)) iyf_offset = self%iyf_offset
  if(present(nhcr      )) nhcr       = self%nhcr
  if(present(nhcr_local)) nhcr_local = self%nhcr_local

  ! l2 norm weights for real array
  if(present(l2w_r))then
    if(associated(l2w_r))deallocate(l2w_r)
    if(self%nyf_local > 0)then
      allocate(l2w_r(self%nhcr_local))
      ! half-complex definition is used in x-dimension
      l2w_r(:)=2
      l2w_r(1::self%nxf*2)=1
      l2w_r(2::self%nxf*2)=1
      if(mod(self%nxf,2)==0)then  ! nxf is even
        l2w_r(self%nxf*2-1::self%nxf*2)=1
        l2w_r(self%nxf*2-0::self%nxf*2)=1
      endif
    else
      allocate(l2w_r(1))          ! dummy
      l2w_r(:)=0
    endif                         ! nyf_local > 0
  endif

  ! l2 norm weights for complex array
  if(present(l2w_c))then 
    if(associated(l2w_c))deallocate(l2w_c)
    if(self%nyf_local > 0)then
      allocate(l2w_c(self%nxf,self%nyf_local))
      ! half-complex definition is used in x-dimension
      l2w_c(:,:)=2
      l2w_c(1,:)=1
      if(mod(self%nxf,2)==0)&     ! nxf is even
        l2w_c(self%nxf,:)=1
    else
      allocate(l2w_c(1,1))        ! dummy
      l2w_c(:,:)=0
    endif                         ! nyf_local > 0
  endif
endsubroutine Get
!-----------------------------------------------------------------------
subroutine Forward(self, input, output)
  ! --- in/out --------------------------------- 
  class(UFFTW3_T), intent(inout)  :: self
  real(C_DOUBLE), intent(in)      :: input(:,:)
  complex(C_DOUBLE), intent(out)  :: output(:,:)
  ! --- const --------------------------------------
  character(len=*), parameter     :: rname = mname//'%Forward'    
  ! --- begin ----------------------------------
    
  ! local gridpoint slab defined ?
  if(self%ny_local>0 ) then
    ! check size ..
    if(any(shape(input)/=[self%nx,self%ny_local])) then
      write(*,'("ERROR - real input sizes do not match:")')
      write(*,'("ERROR -   argument : ",2i6)') shape(input)
      write(*,'("ERROR -   intern   : ",2i6)') self%nx,self%ny_local
      call CheckStop(HERE(rname))
    end if
    ! fill input, pad extra elements with zeros:
    self%input = 0.0
    self%input(1:self%nx,1:self%ny_local) = input
  endif

  ! transform:
  call FFTW_Exec_DFT_r2c(self%plan_fwd, self%input, self%output)

  ! local spectral slab defined ?
  if(self%nyf_local>0) then
    ! check size ..
    if(any(shape(output)/=shape(self%output))) then
      write(*,'("ERROR - complex output sizes do not match:")')
      write(*,'("ERROR -   argument : ",2i6)') shape(output)
      write(*,'("ERROR -   intern   : ",2i6)') shape(self%output)
      call CheckStop(HERE(rname))
    endif
    ! copy result, scale:
    output = self%output * self%ufactor
  endif
endsubroutine Forward
!-----------------------------------------------------------------------
subroutine Inverse(self, output, input)
  ! --- in/out ---------------------------------    
  class(UFFTW3_T), intent(inout)  :: self
  complex(C_DOUBLE), intent(in)   :: output(:,:)
  real(C_DOUBLE), intent(out)     :: input(:,:)
  ! --- const --------------------------------------
  character(len=*), parameter     :: rname = mname//'%Inverse'
  ! --- begin ----------------------------------
    
  ! local spectral slab defined ?
  if(self%nyf_local>0) then
    ! check size ..
    if(any(shape(output)/=shape(self%output))) then
      write(*,'("ERROR - complex input sizes do not match:")')
      write(*,'("ERROR -   argument : ",2i6)') shape(output)
      write(*,'("ERROR -   intern   : ",2i6)') shape(self%output)
      call CheckStop(HERE(rname))
    endif
    ! copy into alligned memory:
    self%output = output
  endif
    
  ! inverse transform:
   call FFTW_Exec_DFT_c2r(self%plan_inv,self%output,self%input2)
    
  ! local gridpoint slab defined ?
  if(self%ny_local>0) then
    ! check size ..
    if(any(shape(input)/=[self%nx,self%ny_local])) then
      write(*,'("ERROR - real output do not match:")')
      write(*,'("ERROR -   argument : ",2i6)') shape(input)
      write(*,'("ERROR -   intern   : ",2i6)') self%nx,self%ny_local
      call CheckStop(HERE(rname))
    endif
    ! extract from padded array, scale result:
    input = self%input2(1:self%nx,1:self%ny_local) * self%ufactor
  endif
endsubroutine Inverse
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
subroutine HC2HCR(self, hc, hcr)
  ! --- in/out ---------------------------------   
  class(UFFTW3_T), intent(in)    :: self
  complex(C_DOUBLE), intent(in)  :: hc(:,:)
  real(C_DOUBLE), intent(out)    :: hcr(:)
  ! --- const --------------------------------------
  character(len=*), parameter :: rname = mname//'%HC2HCR'   
  ! --- local ----------------------------------   
  integer ::  i, j, i0
  ! --- begin ----------------------------------
  
  ! local spectral slab defined ?
  if ( self%nyf_local > 0 ) then
    ! check size ..
    if(any(shape(hc)/=shape(self%output))) then
      write(*,'("ERROR - complex input sizes do not match:")')
      write(*,'("ERROR -   argument : ",2i6)') shape(hc)
      write(*,'("ERROR -   intern   : ",2i6)') shape(self%output)
      call CheckStop(HERE(rname))
    endif

    ! check size ..
    if(size(hcr)/=self%nhcr_local) then
      write(*,'("ERROR - real output size does not match:")')
      write(*,'("ERROR -   argument : ",i12)') size(hcr)
      write(*,'("ERROR -   expected : ",i12)') self%nhcr
      call CheckStop(HERE(rname))
    endif

    i0 = 0                    ! init offset in 1D array
    do j = 1, self%nyf_local  ! loop over local cells
      do i = 1, self%nxf
        ! copy from 2D complex field into real array:
        hcr(i0+1) =  real(hc(i,j))
        hcr(i0+2) = aimag(hc(i,j))
        i0 = i0+2             ! next pair of elements
      enddo
    enddo
  endif  ! nyf_local > 0
endsubroutine HC2HCR
!-----------------------------------------------------------------------
subroutine HCR2HC(self, hcr, hc)
  ! --- in/out ---------------------------------
  class(UFFTW3_T), intent(in)    :: self
  real(C_DOUBLE), intent(in)     :: hcr(:)
  complex(C_DOUBLE), intent(out) :: hc(:,:)
  ! --- const --------------------------------------
  character(len=*), parameter :: rname = mname//'%HC2HCR'
  ! --- local ----------------------------------
  integer ::  i, j, i0   
  ! --- begin ----------------------------------
    
  ! local spectral slab defined ?
  if(self%nyf_local>0) then
    ! check size ..
    if(size(hcr)/=self%nhcr_local) then
      write(*,'("ERROR - real input size does not match:")')
      write(*,'("ERROR -   argument : ",i12)') size(hcr)
      write(*,'("ERROR -   expected : ",i12)') self%nhcr_local
      call CheckStop(HERE(rname))
    endif

    ! check size ..
    if(any(shape(hc)/=shape(self%output))) then
      write(*,'("ERROR - complex output sizes do not match:")')
      write(*,'("ERROR -   argument : ",2i6)') shape(hc)
      write(*,'("ERROR -   intern   : ",2i6)') shape(self%output)
      call CheckStop(HERE(rname))
    endif

    i0 = 0                    ! init offset in 1D array
    do j = 1, self%nyf_local  ! loop over local cells
      do i = 1, self%nxf
        ! copy values from 1D real state to 3D complex array:
        hc(i,j) = cmplx(hcr(i0+1),hcr(i0+2))
        i0 = i0+2             ! next pair of elements
      enddo
    enddo
  endif ! nyf_local > 0
endsubroutine HCR2HC
!-----------------------------------------------------------------------
endmodule UFFTW3_2d
! ********************************************************************
! ***
! *** FFTW3 3D
! ***
! ********************************************************************
module UFFTW3_3d
use, intrinsic :: iso_c_binding
use CheckStop_ml, only: CheckStop
use FFTW3_MPI, only: &
  FFTW_Alloc_Real, FFTW_Alloc_Complex, &
  FFTW_Free, FFTW_Destroy_Plan, FFTW_FLAGS, &
  FFTW_Local_Size   => FFTW_MPI_Local_Size_3d,  &
  FFTW_Plan_DFT_r2c => FFTW_MPI_Plan_DFT_r2c_3d,&
  FFTW_Exec_DFT_r2c => FFTW_MPI_Execute_DFT_r2c,&
  FFTW_Plan_DFT_c2r => FFTW_MPI_Plan_DFT_c2r_3d,&
  FFTW_Exec_DFT_c2r => FFTW_MPI_Execute_DFT_c2r
implicit none
! --- in/out -----------------------------------
private
public :: UFFTW3_T
! --- const ------------------------------------ 
character(len=*), parameter   ::  mname = 'UFFTW3_3d'
! --- types ------------------------------------ 
type UFFTW3_T
  integer                     ::  uout          ! logfile unit
  integer(C_SIZE_T)           ::  nx, ny, nz, n ! physical problem and total size
  real(C_DOUBLE)              ::  ufactor       ! unitairy scale factor
  integer(C_INTPTR_T)         ::  nxf, nyf, nzf ! spectral size
  ! decomposition:
  integer(C_INTPTR_T)         ::  alloc_local, &
                                  nz_local,  iz_offset, &
                                  nzf_local, izf_offset
  ! input and output arrays:
  real(C_DOUBLE), pointer     ::  input (:,:,:),input2(:,:,:)
  complex(C_DOUBLE), pointer  ::  output(:,:,:)
  type(C_PTR)                 ::  input_p,input2_p,output_p
  ! transformation plans:
  type(C_PTR)                 ::  plan_fwd,plan_inv
contains
  procedure :: Init    => Init
  procedure :: Done    => Done
  procedure :: Get     => Get
  procedure :: Forward => Forward
  procedure :: Inverse => Inverse
endtype UFFTW3_T
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
subroutine Init(self, nx, ny, nz, comm)
  ! --- in/out ---------------------------------
  class(UFFTW3_T), intent(out)  :: self
  integer, intent(in)           :: nx, ny, nz
  integer, intent(in)           :: comm
  ! --- const --------------------------------------
  character(len=*), parameter   ::  rname = mname//'%Init'
    
  ! save:
  self%nx = nx
  self%ny = ny
  self%nz = nz
    
  ! total:
  self%n = self%nx * self%ny * self%nz
  
  !! scale factor for unitairy transformations:
  !self%ufactor = 1.0 / sqrt( real(self%n) )
  ! not yet ...
  self%ufactor = 1.0
    
  ! spectral dimensions ; for real-to-complex in Fortran, 
  ! first dimension has roughly half of the complex numbers,
  ! since other half is implicitly equal to the conjugate:
  self%nxf = nx/2 + 1  ! ok for both odd and even
  self%nyf = ny
  self%nzf = nz
    
  ! data distribution is determined by shape of complex data:
  ! - provide spectral size in reversed order as 'integer(C_INTPTR_T)'
  ! - local total size is returned as 'integer(C_INTPTR_T)' too;
  ! - returns span and offset of last dimension as 'integer(C_INTPTR_T)'
  ! note reverse order of dimensions:
  self%alloc_local = FFTW_Local_Size(self%nzf, self%nyf, self%nxf, comm, &
                                     self%nzf_local, self%izf_offset )
  ! adhoc ...
  if(self%nzf_local==0)then
    write(*,'("ufftw3 3D not prepared for zero sized slabs yet ...")')
    call CheckStop(HERE(rname))
  endif

  ! grid point decomposition is the same as spectral ...
  self%nz_local  = self%nzf_local
  self%iz_offset = self%izf_offset
    
  ! real input requires array twice as large as the complex output (in elements);
  ! number of real values is therefore slightly larger than original field
  ! (FFTW3 doc, section 6.5 "Multi-dimensional MPI DFTs of Real Data")
  self%input_p = FFTW_Alloc_Real( 2*self%alloc_local )
  ! allow reference as Fortran array:
  call C_F_Pointer( self%input_p, self%input, (/2*self%nxf,self%nyf,self%nzf_local/) )
    
  ! extra:
  self%input2_p = FFTW_Alloc_Real( 2*self%alloc_local )
  ! allow reference as Fortran array:
  call C_F_Pointer( self%input2_p, self%input2, (/2*self%nxf,self%nyf,self%nzf_local/) )
    
  ! complex output field:
  self%output_p = FFTW_Alloc_Complex( self%alloc_local )
  ! allow reference as Fortran array:
  call C_F_Pointer( self%output_p, self%output, (/self%nxf,self%nyf,self%nzf_local/) )

  ! setup real-to-complex (half-complex) transphorm;
  ! dimension orders need to be reversed !
  self%plan_fwd = FFTW_Plan_DFT_r2c(self%nz, self%ny, self%nx, &
                              self%input, self%output, comm, FFTW_FLAGS)

  ! setup inverse:
  self%plan_inv = FFTW_Plan_DFT_c2r(self%nz, self%ny, self%nx, &
                              self%output, self%input2, comm, FFTW_FLAGS)
endsubroutine Init
!-----------------------------------------------------------------------
subroutine Done(self)
  ! --- in/out ---------------------------------
  class(UFFTW3_T), intent(inout)  ::  self
  ! --- const --------------------------------------
  character(len=*), parameter     ::  rname = mname//'%Done'
  ! --- begin ----------------------------------
    
  ! done:
  call FFTW_Destroy_Plan(self%plan_fwd)
  call FFTW_Destroy_Plan(self%plan_inv)

  ! clear:
  call FFTW_Free(self%input_p)
  call FFTW_Free(self%input2_p)
  call FFTW_Free(self%output_p)

  ! dummy:
  self%nx = -999
  self%ny = -999
  self%nz = -999
endsubroutine Done
!-----------------------------------------------------------------------
subroutine Get(self, &
               nx, ny, nz, nz_local, iz_offset, &
               nxf, nyf, nzf, nzf_local, izf_offset)  
  ! --- in/out ---------------------------------  
  class(UFFTW3_T), intent(inout)  :: self
  integer, intent(out), optional  :: nx, ny, nz, nz_local, iz_offset, &
                                     nxf, nyf, nzf, nzf_local, izf_offset
  ! --- const --------------------------------------
  character(len=*), parameter     ::  rname = mname//'%Get'
  ! --- begin ----------------------------------
    
  ! return values:
  if(present(nx        )) nx         = self%nx
  if(present(ny        )) ny         = self%ny
  if(present(nz        )) nz         = self%nz
  if(present(nz_local  )) nz_local   = self%nz_local
  if(present(iz_offset )) iz_offset  = self%iz_offset
  if(present(nxf       )) nxf        = self%nxf
  if(present(nyf       )) nyf        = self%nyf
  if(present(nzf       )) nzf        = self%nzf
  if(present(nzf_local )) nzf_local  = self%nzf_local
  if(present(izf_offset)) izf_offset = self%izf_offset
endsubroutine Get
!-----------------------------------------------------------------------
subroutine Forward(self, input, output)
  ! --- in/out ---------------------------------
  class(UFFTW3_T), intent(inout)  :: self
  real(C_DOUBLE), intent(in)      :: input(:,:,:)
  complex(C_DOUBLE), intent(out)  :: output(:,:,:)
  ! --- const --------------------------------------
  character(len=*), parameter     :: rname = mname//'%Forward'
  ! --- local ----------------------------------
  ! check size ..
  if(any(shape(input)/=[self%nx,self%ny,self%nzf_local])) then
    write(*,'("ERROR - real input sizes do not match:")')
    write(*,'("ERROR -   argument : ",3i6)') shape(input)
    write(*,'("ERROR -   intern   : ",3i6)') self%nx,self%ny,self%nz
    call CheckStop(HERE(rname))
  endif
  ! fill input, pad extra elements with zeros:
  self%input = 0.0
  self%input(1:self%nx,1:self%ny,1:self%nzf_local) = input

  ! transform:
  call FFTW_Exec_DFT_r2c(self%plan_fwd, self%input, self%output)

  ! check size ..
  if(any(shape(output)/=shape(self%output))) then
    write(*,'("ERROR - complex output sizes do not match:")')
    write(*,'("ERROR -   argument : ",3i6)') shape(output)
    write(*,'("ERROR -   intern   : ",3i6)') shape(self%output)
    call CheckStop(HERE(rname))
  endif
  ! copy result, scale:
  output = self%output * self%ufactor
endsubroutine Forward
!-----------------------------------------------------------------------
subroutine Inverse(self, output, input)
  ! --- in/out ---------------------------------
  class(UFFTW3_T), intent(inout)  :: self
  complex(C_DOUBLE), intent(in)   :: output(:,:,:)
  real(C_DOUBLE), intent(out)     :: input(:,:,:)
  ! --- const --------------------------------------
  character(len=*), parameter     ::  rname = mname//'%Inverse'
  ! --- begin ----------------------------------
    
  ! check size ..
  if(any(shape(output)/=shape(self%output))) then
    write(*,'("ERROR - complex input sizes do not match:")')
    write(*,'("ERROR -   argument : ",3i6)') shape(output)
    write(*,'("ERROR -   intern   : ",3i6)') shape(self%output)
    call CheckStop(HERE(rname))
  endif
  ! copy into alligned memory:
  self%output = output

  ! inverse transform:
  call FFTW_Exec_DFT_c2r(self%plan_inv, self%output, self%input2)

  ! check size ..
  if(any(shape(input)/=[self%nx,self%ny,self%nzf_local])) then
    write(*,'("ERROR - real output do not match:")')
    write(*,'("ERROR -   argument : ",3i6)') shape(input)
    write(*,'("ERROR -   intern   : ",3i6)') self%nx,self%ny,self%nz
    call CheckStop(HERE(rname))
  endif
  ! extract from padded array, scale result:
  input = self%input2(1:self%nx,1:self%ny,1:self%nzf_local) * self%ufactor
endsubroutine Inverse
!-----------------------------------------------------------------------
endmodule UFFTW3_3d
!-----------------------------------------------------------------------
module UFFTW3_ml
use UFFTW3_2d, only: UFFTW3_T2d=>UFFTW3_T
use UFFTW3_3d, only: UFFTW3_T3d=>UFFTW3_T
endmodule UFFTW3_ml
