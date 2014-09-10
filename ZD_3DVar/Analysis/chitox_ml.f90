module ChiToX_ml
use CheckStop_ml,only: CheckStop
use spectralcov, only: nx,ny,stddev,nxex,nyex,nlev,nchemobs,ichemobs,&
                       nv1,sqrt_lambda,vt,sqrt_gamma,nkstar,ikstar,&
                       wsave, lensav, nmin, nmax
use UFFTW3_ml,   only: UFFTW3_T2d
implicit none
private 
public :: ufft, ChiToX, ChiToX_adj
type(UFFTW3_T2d) ::  ufft
contains
!-----------------------------------------------------------------------
subroutine ChiToX(chi_hc,x)
!-----------------------------------------------------------------------
! @description
! Transform control variable chi to x
! Inverse preconditioning of the control variable
! @author M Kahnert, A Segers
!-----------------------------------------------------------------------
! --- in/out ----------------------------------
  complex, intent(in):: chi_hc(:,:,:) ! (nxf,lnyf,nv1)
  real, intent(out)  :: x(:,:,:,:)    ! (nxex,nyex,nlev,nchemobs)
! --- local -----------------------------------
  integer :: k, kk, l, ibox, m, n, ln, ik, ndim
  complex, pointer ::  temp1(:), temp2(:), temp3(:,:,:)
  integer ::  nxg, nyg, lnyg, iyg0, nxf, nyf, lnyf, iyf0
  integer ::  iy1, iy2, nys
  complex, allocatable:: ww(:,:)      ! (nxf,nyf)
  real, allocatable   :: c2(:,:)      ! (nxex,nyex)
!-----------------------------------------------------------------------
! Compute:
!           c = S * Finv * Linv * X * \Lambda^{1/2} * \chi
! where:
!   S       = diag.matrix with grid point std.dev.
!   Finv    = inverse 2D Fourrier transform
!   Linv    = diag(sqrt(gamma))
!   \Lambda = diagonal matrix containing the reduced set 
!              of eigenvalues of the spectral covariance matrix, 
!   X       = matrix containing the corresponding eigenvectors,
!               here stored in 'vt'
!-----------------------------------------------------------------------
  call ufft%Get(nx=nxg ,ny=nyg ,ny_local=lnyg ,iy_offset=iyg0 ,&
                nxf=nxf,nyf=nyf,nyf_local=lnyf,iyf_offset=iyf0) ! sizes

  if(lnyf>0)then              ! defined spectral y-slab 
    if(any(shape(chi_hc)/=[nxf,lnyf,nv1]))then    ! check input
      print "(A,2(1X,A,3(I0,',')))","ChiToX: shape(chi_hc)",&
        "chi_hc:",shape(chi_hc),"[nxf,lnyf,nv1]:",nxf,lnyf,nv1
      call CheckStop("ChiToX: shape(chi_hc)")
    endif
    ndim = nchemobs * nlev                        ! size of diagonal blocks
    allocate(temp1(nv1), temp2(ndim))             ! temporary storage
    allocate(temp3(ndim,nxf,lnyf), ww(nxf,lnyf))  ! storage
    temp3(:,:,:)=0.0
    do ln = 1, lnyf           ! loop over spectral coeffs
      n = iyf0 + ln           ! full number
      do m = 1, nxf           ! loop over (half of) spectral coeffs
        ik = ikstar(m,n)      ! index in 1D spectral coefs:
        if(ik>nkstar) cycle   ! elliptic truncation
    !-- temp3 =  Gamma^1/2 X Lambda^1/2 Chi
        temp1(:)      = chi_hc(m,ln,:) * sqrt_lambda(:,ik)
        temp2(:)      = matmul(vt(:,:,ik), temp1)
        temp3(:,m,ln) = temp2(:) * sqrt_gamma(:,ik)
      enddo
    enddo
    deallocate(temp1,temp2)   ! cleanup
  else
    allocate(ww(1,1))         ! dummy
  endif

  if(lnyg>0)then              ! defined gridpoint y-slab
    if(iyg0<ny)then           ! output overlap ?
      iy1 =     iyg0+1        ! global index space
      iy2 = min(iyg0+lnyg,ny)
      nys = iy2-iy1+1         ! expected size
      if(any(shape(x)/=[nx,nys,nlev,nChemObs]))then ! check output
        print "(A,2(1X,A,4(I0,',')))","ChiToX: shape(x)",&
          "x:",shape(x),"[nx,nys,nlev,nChemObs]:",nx,nys,nlev,nChemObs
        call CheckStop("ChiToX: shape(x)")
      endif
    endif
    allocate(c2(nxg,lnyg))    ! temporary storage
  else
    allocate(c2(1,1))         ! dummy
  endif

  x = 0.0                     ! init result
  do k = 1, nchemobs          ! loop over tracer
    kk = ichemobs(k)          ! local index
    do l = 1, nlev            ! loop over levels
      ibox = (k-1)*nlev+l     ! index in block matrix
!-----------------------------------------------------------------------
! Inverse 2d-FFT:
!-----------------------------------------------------------------------
      if(lnyf>0)&             ! local y-slab defined ?         
        ww = temp3(ibox,:,:)  ! copy non-redundant part of complex coeffs:

      call ufft%Inverse(ww,c2) ! tranform

      ! copy elements in model domain, extension zone remains zero
      if((lnyg>0).and.(iyg0<ny)) then
        iy1 =     iyg0+1      ! global index space of local slab
        iy2 = min(iyg0+lnyg,ny)
        nys = iy2-iy1+1       ! overlap size
                              ! store multiplication with gridpoint std.dev.
        x(1:nx,1:nys,l,k) = c2 * stddev(1:nx,iy1:iy2,l,kk)
      endif
    enddo   ! level l
  enddo     ! component k

  if(associated(temp3)) deallocate(temp3)  ! cleanup
  if(allocated(ww))     deallocate(ww)
  if(allocated(c2))     deallocate(c2)
endsubroutine ChiToX
!-----------------------------------------------------------------------
subroutine ChiToX_adj(chi_hc,x)
!-----------------------------------------------------------------------
! @description
! Adjoint of the transformation ChiToX
! @author M Kahnert, A Segers
!-----------------------------------------------------------------------
! --- in/out ----------------------------------
  real, intent(in)                :: x(:,:,:,:)     ! (lnx ,lny ,nlev,nchemobs)
  complex, intent(out)            :: chi_hc(:,:,:)  ! (nxf,lnyf,nv1)
! --- local -----------------------------------
  integer :: i,j,k,kk,l,ik,m,n,ir,iostat,ibox, nxy
  integer :: nxg, lnyg, iyg0, nxf, lnyf, iyf0
  integer :: iy1, iy2, ln
  real, allocatable    ::  eta(:,:)  ! (nxex,nyex)
  complex, allocatable ::  ww(:,:)   ! (nxf,nyf)

  ! sizes
  call ufft%Get(nx=nxg ,ny_local=lnyg ,iy_offset=iyg0, &
                nxf=nxf,nyf_local=lnyf,iyf_offset=iyf0)

  if(lnyg>0)then                    ! defined local gridpoint y-slab
    if(any(shape(x)/=[nxg,lnyg,nlev,nChemObs]))then ! check input
      print "(A,2(1X,A,4(I0,',')))","ChiToX_adj: shape(x)",&
        "x:",shape(x),"[nxg,lnyg,nlev,nChemObs]:",nxg,lnyg,nlev,nChemObs
      call CheckStop("ChiToX_adj: shape(x)")       
    endif
    allocate(eta(nxg,lnyg))         ! storage for testing
  else
    allocate(eta(1,1))              ! dummy
  endif

  if(lnyf>0)then                    ! defined local spectral y-slab
    if(any(shape(chi_hc)/=[nxf,lnyf,nv1]))then  ! check output
      print "(A,2(1X,A,3(I0,',')))","ChiToX_adj: shape(chi_hc)",&
        "chi_hc:",shape(chi_hc),"[nxf,lnyf,nv1]:",nxf,lnyf,nv1
      call CheckStop("ChiToX_adj: shape(chi_hc)")  
    endif
    allocate(ww(nxf,lnyf))          ! storage for testing
  else
    allocate(ww(1,1))               ! dummy
  endif


  nxy = real( nxex * nyex ) ! factor for unitary transforms
  chi_hc = cmplx(0e0,0e0)   ! init result
    
  do k = 1, nchemobs        ! loop over tracers
    kk = ichemobs(k)        ! local index
    do l = 1, nlev          ! loop over levels
      ibox = (k-1)*nlev + l ! index in block matrix        
      eta = 0.0       ! extension zone is padded with zero's, fill this first
      if(lnyg>0)then  ! defined local gridpoint y-slab
        iy1 = iyg0+1    ! local slab
        iy2 = iyg0+lnyg
        ! Multiply from the left with S, where S is a diagonal matrix
        ! containing the standard deviations sqrt(sigma):
        eta(1:nxg,1:lnyg) = x(1:nxg,1:lnyg,l,k) * stddev(1:nxg,iy1:iy2,l,kk)
      endif

      call ufft%Forward(eta,ww)  ! transform
      ! To have the same results as using fftpack, 
      ! it seems necessary to devide by 'nxex*nyex' ... Strange!
      ! The factor applied under "b)" below originates from the theory that:
      !        F^{-H} = (nxex*nyex) * F
      ! From the comment in 'cmplx_fft_2d.F' it is understood that
      ! that fftpack's "cfft2f" routine has a scaling 1/N already ;
      ! therefor apply this here too now:
      ww = ww / nxy

      if(lnyf>0)then ! defined local spectral y-slab
        do m = 1, nxf           ! loop over spectral coefs
          do ln = 1, lnyf       ! loop over local spectral coefs
            n = iyf0 + ln       ! global index
            ik = ikstar(m,n)    ! index in 1D spectral coeff array            
            if(ik>nkstar) cycle ! elliptic truncation
            !-----------------------------------------------------------------------
            !     b) multiply with (nxex*nyex);
            !        note that this resets the above division by 'nxy' !
            !-----------------------------------------------------------------------
            ww(m,ln) = ww(m,ln) * nxy
            !-----------------------------------------------------------------------
            !     Multiply from the left with L^{-1}, where L is a diagonal matrix
            !     with elements 1/sqrt(\gamma_ibox) along the diagonal, where
            !     \gamma_ibox denotes the bi-Fourier spectral density of the variance
            !     for level/component ibox (ibox is a level/component super-index):
            !-----------------------------------------------------------------------
            ww(m,ln) = ww(m,ln) * sqrt_gamma(ibox,ik)
            !-----------------------------------------------------------------------
            !     Compute chi = \Lambda^{1/2} * X^{T} * \chi, where \Lambda is
            !     the diagonal matrix containing the reduced set of eigenvalues of
            !     the spectral covariance matrix, and X is the matrix containing the
            !     corresponding eigenvectors:
            !-----------------------------------------------------------------------
            do ir = 1, nv1
              chi_hc(m,ln,ir) = chi_hc(m,ln,ir) + &
                    vt(ibox,ir,ik) * sqrt_lambda(ir,ik) * ww(m,ln)
            enddo
          enddo ! spectral component n
        enddo   ! spectral component m
      endif
    enddo ! altitude l
  enddo   ! component k
  
  if(allocated(eta))  deallocate(eta) ! cleanup
  if(allocated(ww))   deallocate(ww)
endsubroutine ChiToX_adj
!-----------------------------------------------------------------------
endmodule ChiToX_ml


