#define STRING2(x) #x
#define STRING(x) STRING2(x)
#define HERE(MSG) MSG//" ("//__FILE__//":"//STRING(__LINE__)//")."
module ChiToX_ml
use ModelConstants_ml,only: MasterProc,NPROC,NPROCY,IIFULLDOM,JJFULLDOM
use Par_ml,           only: me,limax,ljmax,gi0,gi1,gj0,gj1
use CheckStop_ml,only: CheckStop
use spectralcov, only: nx,ny,stddev,nxex,nyex,nlev,nchemobs,ichemobs,&
                       nv1,sqrt_lambda,vt,sqrt_gamma,nkstar,ikstar,&
                       wsave, lensav, nmin, nmax
use spectralcov, only: kx,ky,kxmin,kymin,mm,nn
use DA_ml,            only: debug=>DA_DEBUG,dafmt=>da_fmt_msg
use UFFTW3_ml,   only: UFFTW3_T2d
use mpi,         only: MPI_COMM_WORLD,MPI_BARRIER,MPI_BCAST,MPI_ALLREDUCE,MPI_LAND,MPI_SUM,&
                       MPI_IN_PLACE,MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_COMPLEX
implicit none
private 
public :: ChiToX, ChiToX_adj, chi_Init, chi_Done, &
          chi_vec2fc, chi_fc2vec, chi_fc2hc, chi_hc2fc 
type(UFFTW3_T2d) :: ufft
logical, public  :: matched_domain
integer, allocatable, public :: iyf(:,:)
!-----------------------------------------------------------------------
INTEGER :: INFO
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
subroutine chi_Init()
! --- local -----------------------------------
  integer :: i, lnx,lny,iy0,iy1, lnxf,lnyf,iyf0,iyf1
  character(len=*), parameter :: &
    DOMFMT="(A,I3,3(1X,A5,'(',I3,'x',I3,')',I3,'--',I3,'x',I3,'--',I3))"
!-----------------------------------------------------------------------
  call ufft%Init(nxex,nyex,comm=MPI_COMM_WORLD)
  ! extract info on y-decomposition used by fftw ;
  ! note that this is defined on the extended domain
  call ufft%Get(nx=lnx,  ny_local=lny,  iy_offset=iy0,&   ! real
               nxf=lnxf,nyf_local=lnyf,iyf_offset=iyf0) ! complex
  iy0 =iy0 +1;iy1 =iy0 +lny -1
  iyf0=iyf0+1;iyf1=iyf0+lnyf-1
  if(debug)then
    if(MasterProc)then
      print dafmt,"chi_Init Domain"
      print domfmt,'@',me,&
        'model',IIFULLDOM,JJFULLDOM,1,IIFULLDOM,1,JJFULLDOM,&
        'trim',nx,ny,1,nx,1,ny,'ext.',nxex,nyex,1,nxex,1,nyex
    endif
    do i=0,NPROC-1
      if(i==me) print DOMFMT,'@',me,&
        'local',limax,ljmax,gi0,gi1 ,gj0 ,gj1 ,&
        'dx|du',lnx  ,lny  ,1  ,lnx ,iy0 ,iy1 ,&
        'fftw3',lnxf ,lnyf ,1  ,lnxf,iyf0,iyf1
      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    enddo
  endif
  matched_domain=(lnyf==ljmax)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,matched_domain,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,INFO)
  allocate(iyf(0:1,0:NPROC-1))
  iyf(:,:)=0
  iyf(:,me)=[iyf0,iyf1]
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,iyf,NPROC*2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,INFO)
endsubroutine chi_Init
!-----------------------------------------------------------------------
subroutine chi_Done()
!-----------------------------------------------------------------------
  call ufft%Done()
  if(allocated(iyf))deallocate(iyf)
endsubroutine chi_Done
!-----------------------------------------------------------------------
subroutine chi_vec2fc(nv2,chi,chi_arr)
!-----------------------------------------------------------------------
! Copy chi (in compact storage format) into chi_arr (in matrix format),
! taking into account the relation chi_arr(i,m,n)=conjg(chi_arr(i,-m,-n))
!----------------------------------------------------------------------
! --- in/out ----------------------------------
  integer, intent(in)           :: nv2
  real, intent(in)              :: chi(nv2)
  complex, intent(out), pointer :: chi_arr(:,:,:)
! --- local -----------------------------------
  integer :: kx1,ky1,kxm1,kym1,kx2
  integer :: i,j,k,m,n,m1,n1,mneg,nneg
!-----------------------------------------------------------------------
  kx1=kx+1
  ky1=ky+1
  kxm1=kxmin+1
  kym1=kymin+1
  kx2=nxex-kxmin+1
  if(associated(chi_arr))deallocate(chi_arr)
  allocate(chi_arr(nv1,nxex,nyex))  
  chi_arr=cmplx(0e0,0e0)
!-----------------------------------------------------------------------
! independent elements:
!-----------------------------------------------------------------------
  k=0
  do i=1,nv1
    do m=1,kx1
      k=k+1
      chi_arr(i,m,1)=cmplx(chi(k),0e0)
    enddo
    do m=1,kx1
      do n=2,ky1
        k=k+1
        chi_arr(i,m,n)=cmplx(chi(k),0e0)
      enddo
    enddo
    do m=kx2,nxex
      do n=2,kym1
        k=k+1
        chi_arr(i,m,n)=cmplx(chi(k),0e0)
      enddo
    enddo
    do m=2,kxm1
      k=k+1
      chi_arr(i,m,1)=chi_arr(i,m,1)+cmplx(0e0,chi(k))
    enddo
    do n=2,kym1
      k=k+1
      chi_arr(i,1,n)=chi_arr(i,1,n)+cmplx(0e0,chi(k))
    enddo
    do m=2,kx1
      do n=2,ky1
        if(mod(nxex,2)==0.and.mod(nyex,2)==0.and.m==kx1.and.n==ky1)then
          continue
        else
          k=k+1
          chi_arr(i,m,n)=chi_arr(i,m,n)+cmplx(0e0,chi(k))
        endif
      enddo
    enddo
    do m=kx2,nxex
      do n=2,kym1
        k=k+1
        chi_arr(i,m,n)=chi_arr(i,m,n)+cmplx(0e0,chi(k))
      enddo
    enddo
  enddo
!-----------------------------------------------------------------------
! dependent elements:
!-----------------------------------------------------------------------
  do i=1,nv1
    do m1=-kxmin,kx
      m=mm(m1)
      mneg=mm(-m1)
      do n1=-kymin,-1
        n=nn(n1)
        nneg=nn(-n1)
        chi_arr(i,m,n)=conjg(chi_arr(i,mneg,nneg))
      enddo
    enddo
    do m1=-kxmin,-1
      m=mm(m1)
      mneg=mm(-m1)
      n=nn(0)
      chi_arr(i,m,n)=conjg(chi_arr(i,mneg,n))
    enddo
    if(mod(nyex,2)==0)then
      do m1=-kxmin,-1
        m=mm(m1)
        mneg=mm(-m1)
        n=nn(ky)
        chi_arr(i,m,n)=conjg(chi_arr(i,mneg,n))
      enddo
    endif
! zero elements:
    m=mm(kx);n=nn(ky)
    chi_arr(i,1,1)=cmplx(real(chi_arr(i,1,1)),0e0)
    if(mod(nxex,2)==0)&
      chi_arr(i,m,1)=cmplx(real(chi_arr(i,m,1)),0e0)
    if(mod(nyex,2)==0)&
      chi_arr(i,1,n)=cmplx(real(chi_arr(i,1,n)),0e0)
    if(mod(nxex,2)==0.and.mod(nyex,2)==0)&
      chi_arr(i,m,n)=cmplx(real(chi_arr(i,m,n)),0e0)
  enddo
endsubroutine chi_vec2fc
!-----------------------------------------------------------------------
subroutine chi_fc2vec(nv2,chi,chi_arr)
! only add independent elements of chi_arr to array chi
! --- in/out ----------------------------------
  integer, intent(in) :: nv2
  real, intent(inout) :: chi(nv2)
  complex, intent(in) :: chi_arr(nv1,nxex,nyex)
! --- local -----------------------------------
  integer :: kx1,ky1,kxm1,kym1,kx2  
  integer :: i,j,k,m,n
!-----------------------------------------------------------------------
  kx1=kx+1
  ky1=ky+1
  kxm1=kxmin+1
  kym1=kymin+1
  kx2=nxex-kxmin+1
!-----------------------------------------------------------------------
  k=0
  do i=1,nv1
    do m=1,kx1
      k=k+1
      chi(k)=chi(k)+real(chi_arr(i,m,1))
    enddo
    do m=1,kx1
      do n=2,ky1
        k=k+1
        chi(k)=chi(k)+real(chi_arr(i,m,n))
      enddo
    enddo
    do m=kx2,nxex
      do n=2,kym1
        k=k+1
        chi(k)=chi(k)+real(chi_arr(i,m,n))
      enddo
    enddo
    do m=2,kxm1
      k=k+1
      chi(k)=chi(k)+aimag(chi_arr(i,m,1))
    enddo
    do n=2,kym1
      k=k+1
      chi(k)=chi(k)+aimag(chi_arr(i,1,n))
    enddo
    do m=2,kx1
      do n=2,ky1
        if(mod(nxex,2)==0.and.mod(nyex,2)==0.and.m==kx1.and.n==ky1)then
          continue
        else
          k=k+1
          chi(k)=chi(k)+aimag(chi_arr(i,m,n))
        endif
      enddo
    enddo
    do m=kx2,nxex
      do n=2,kym1
        k=k+1
        chi(k)=chi(k)+aimag(chi_arr(i,m,n))
      enddo
    enddo
  enddo
!-----------------------------------------------------------------------
endsubroutine chi_fc2vec
!-----------------------------------------------------------------------
subroutine chi_fc2hc(chi_arr,chi_hc)
!----------------------------------------------------------------------
! Fill non-redundant part, or above called "independent" part.
! Use same real-to-complex decoding to transform 'chi' into 'chi_hc' ;
! experiments showd that in a serial run the first part of the
! full complex field 'chi_arr' is a proper input to the fftw routines,
! that is, the first half of the first dimension is sufficient:
!     chi_arr(iv,1:nxf,1:nyf)
! --- in/out ----------------------------------
  complex, intent(in) :: chi_arr(nv1,nxex,nyex)
  complex, intent(out), pointer :: chi_hc(:,:,:)
! --- local -----------------------------------
  integer :: i,nxf,nyf,lnyf,iyf0,iyf1
!-----------------------------------------------------------------------
  call ufft%Get(nxf=nxf,nyf=nyf,nyf_local=lnyf,iyf_offset=iyf0) ! sizes
  if(associated(chi_hc))deallocate(chi_hc)
  if(lnyf>0)then            ! defined local y-slab
    allocate(chi_hc(nxf,lnyf,nv1))
    iyf0=iyf0+1
    iyf1=iyf0+lnyf-1
    forall(i=1:nv1) &                          ! copy non-redundant slab, 
      chi_hc(:,:,i) = chi_arr(i,:nxf,iyf0:iyf1)! 'chi_arr' is the same on every proc
  else
    allocate(chi_hc(1,1,1)) ! dummy
  endif
!-----------------------------------------------------------------------
endsubroutine chi_fc2hc
!-----------------------------------------------------------------------
subroutine chi_hc2fc(chi_arr,chi_hc)
!----------------------------------------------------------------------
! full complex form from FFTW3's half-complex
! --- in/out ----------------------------------
  complex, intent(in)  :: chi_hc(:,:,:)
  complex, intent(out) :: chi_arr(nv1,nxex,nyex)
! --- local -----------------------------------
  integer :: i,m,n,nxf,nyf,lnyf,iyf0,iyf1
  complex, pointer :: chi_hc_glb(:,:,:)=>null()
!-----------------------------------------------------------------------
  call ufft%Get(nxf=nxf,nyf=nyf,nyf_local=lnyf,iyf_offset=iyf0) ! sizes
 !!collect chi_hc on every process:
  allocate(chi_hc_glb(nxf,nyf,nv1))
  if(iyf(0,me)<iyf(1,me))&
    chi_hc_glb(:,iyf(0,me):iyf(1,me),:)=chi_hc
  do n=0,NPROC-1
    if(iyf(0,n)==iyf(1,n))cycle
    CALL MPI_BCAST(chi_hc_glb(:,iyf(0,n):iyf(1,n),:),&
              size(chi_hc_glb(:,iyf(0,n):iyf(1,n),:)),&
      MPI_DOUBLE_COMPLEX,n,MPI_COMM_WORLD,INFO)
  enddo
  if(lnyf>0)then            ! defined local y-slab
    if(any(shape(chi_hc)/=[nxf,lnyf,nv1]))then    ! check input
      print "(A,2(1X,A,3(I0,',')))","chi_fc2hc: shape(chi_hc)",&
        "chi_hc:",shape(chi_hc),"[nxf,lnyf,nv1]:",nxf,lnyf,nv1
      call CheckStop(HERE("chi_fc2hc: shape(chi_hc)"))
    endif
    do i=1,nv1
      ! half complex X-domain
      do m=1,nxf
        chi_arr(i,m,:)=chi_hc_glb(m,:,i)
      enddo
      ! full complex X-domain
      do m=nxex,nxf+1,-1
        do n=nyex,2,-1
          chi_arr(i,m,n)=conjg(chi_hc_glb(nxex-m+2,nyex-n+2,i))
        enddo
      enddo
    enddo
  endif
  deallocate(chi_hc_glb)
!-----------------------------------------------------------------------
endsubroutine chi_hc2fc
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
  integer :: k, kk, l, ibox, m, n, ln, ik, ndim, ierr=0
  integer :: nxg, nyg, lnyg, iyg0, nxf, nyf, lnyf, iyf0
  integer :: iy1, iy2, nys
  complex, allocatable:: ww(:,:)  ! (nxf,nyf)
  real, allocatable   :: c2(:,:)  ! (nxex,nyex)
  complex, pointer :: t1(:)=>null(), t2(:)=>null(), t3(:,:,:)=>null()
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
      call CheckStop(HERE("shape(chi_hc)"))
    endif
    ndim = nchemobs * nlev                        ! size of diagonal blocks
    allocate(t1(nv1), t2(ndim))             ! temporary storage
    allocate(t3(ndim,nxf,lnyf), ww(nxf,lnyf))  ! storage
    t3(:,:,:)=0.0
    do ln = 1, lnyf           ! loop over spectral coeffs
      n = iyf0 + ln           ! full number
      do m = 1, nxf           ! loop over (half of) spectral coeffs
        ik = ikstar(m,n)      ! index in 1D spectral coefs:
        if(ik>nkstar) cycle   ! elliptic truncation
    !-- t3 =  Gamma^1/2 X Lambda^1/2 Chi
        t1(:)      = chi_hc(m,ln,:) * sqrt_lambda(:,ik)
        t2(:)      = matmul(vt(:,:,ik), t1)
        t3(:,m,ln) = t2(:) * sqrt_gamma(:,ik)
      enddo
    enddo
    deallocate(t1,t2,stat=ierr);call CheckStop(ierr,HERE(''))   ! cleanup
  else
    allocate(ww(1,1))         ! dummy
  endif

  if(lnyg>0)then              ! defined gridpoint y-slab
    if(iyg0<ny)then           ! output overlap ?
      iy1 =     iyg0+1        ! global index space
      iy2 = min(iyg0+lnyg,ny)
      nys = iy2-iy1+1         ! expected size
      if(any(shape(x)<[nx,nys,nlev,nChemObs]))then ! check output
        print "(A,2(1X,A,4(I0,',')))","ChiToX: shape(x)",&
          "x:",shape(x),"[nx,nys,nlev,nChemObs]:",nx,nys,nlev,nChemObs
        call CheckStop(HERE("shape(x)"))
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
        ww = t3(ibox,:,:)     ! copy non-redundant part of complex coeffs:

      call ufft%Inverse(ww,c2)! transform

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

  if(associated(t3)) deallocate(t3)  ! cleanup
  if(allocated(ww))  deallocate(ww)
  if(allocated(c2))  deallocate(c2)
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
      call CheckStop(HERE("shape(x)"))       
    endif
    allocate(eta(nxg,lnyg))         ! storage for testing
  else
    allocate(eta(1,1))              ! dummy
  endif

  if(lnyf>0)then                    ! defined local spectral y-slab
    if(any(shape(chi_hc)/=[nxf,lnyf,nv1]))then  ! check output
      print "(A,2(1X,A,3(I0,',')))","ChiToX_adj: shape(chi_hc)",&
        "chi_hc:",shape(chi_hc),"[nxf,lnyf,nv1]:",nxf,lnyf,nv1
      call CheckStop(HERE("shape(chi_hc)"))  
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


