#define STRING2(x) #x
#define STRING(x) STRING2(x)
#define HERE(MSG) MSG//" ("//__FILE__//":"//STRING(__LINE__)//")."
module covmat_ml
use ModelConstants_ml,only: MasterProc
use Par_ml,           only: me
use DA_ml,            only: debug=>DEBUG_DA,dafmt=>da_fmt_msg,damsg=>da_msg
use CheckStop_ml,     only: CheckStop
use spectralcov,      only: covmat,ucovmat,ncorr,&
                            nx,ny,nlev,nchem,nex,nxex,nyex,&
                            ichemInv,nchemobs,ichemobs,nchemnoobs,ichemnoobs,&
                            nbg1,nv1,nstarmax,&
                            nkstar,ikstar,kstar,nmin,nmax,mm,nn,&
                            kx,kxmin,ky,kymin,lensav,wsave,&
                            stddev,sqrt_gamma,vt,sqrt_lambda,spec_allocate,&
                            nstar,w=>intweight,mt,nt,DAPREC
#ifdef MKL
! use MKL95_LAPACK,     only: lamch=>dlamch, spevx
use My_LAPACK77_ml,   only: lamch=>dlamch
use MKL95_LAPACK,     only: spevx
#else
use My_LAPACK77_ml,   only: lamch=>dlamch, spevx=>dspevx
#endif
implicit none
real,public,save :: dxdim=0.25, dydim=0.125
integer :: nv=500, nbg=500
contains
subroutine set_chemobs_idx(nchem,observedVar)
implicit none
integer, intent(in) :: nchem
logical, intent(in) :: observedVar(nchem)
integer :: ierr,nvar
  nchemobs=count(observedVar)
  nchemnoobs=nchem-nchemobs
  ncorr=((nlev*nchemobs)*(nlev*nchemobs+1))/2
! call allocvar(ichemInv,nchem,(/1,nchem/),'chemInv',dafmt)
! call allocvar(ichemObs,nchemObs,(/1,nchemObs/),'ichemObs',dafmt)
! call allocvar(ichemNoObs,nchemNoObs,(/0,nchemNoObs/),'ichemNoObs',dafmt)
  if(allocated(ichemInv))then
    call CheckStop(size(ichemInv)/=nchem,HERE('Wrong size ichemInv'))
  else
    allocate(ichemInv(nchem),stat=ierr)
    call CheckStop(ierr,HERE('Allocate ICHEMINV'))
    ichemInv=0
  endif
  call CheckStop(nchemobs<=0.or.nchemobs>nchem,HERE('Wrong nchemObs'))
  if(allocated(ichemObs))then
    call CheckStop(size(ichemObs)/=nchemObs,HERE('Wrong size ichemObs'))
  else
    allocate(ichemObs(nchemObs),stat=ierr)
    call CheckStop(ierr,HERE('Allocate ICHEMOBS'))
    ichemObs=0
  endif
  call CheckStop(nchemNoObs<0.or.nchemNoObs>nchem,HERE('Wrong nchemNoObs'))
  if(allocated(ichemNoObs))then
    call CheckStop(size(ichemNoObs)/=nchemNoObs,HERE('Wrong size ichemNoObs'))
  elseif(nchemNoObs==0) then
    if(MasterProc)print dafmt,'WARNING nchemNoObs==0'
  else
    allocate(ichemNoObs(nchemNoObs),stat=ierr)
    call CheckStop(ierr,HERE('Allocate ICHEMNOOBS'))
    ichemNoObs=0
  endif
  nchemObs=0
  nchemNoObs=0
  do nvar=1,nchem
    if(observedVar(nvar))then
      nchemObs=nchemObs+1
      ichemObs(nchemObs)=nvar
      ichemInv(nvar)=nchemObs
    else
      nchemNoObs=nchemNoObs+1
      ichemNoObs(nchemNoObs)=nvar
      ichemInv(nvar)=nchemNoObs
    endif
  enddo
end subroutine set_chemobs_idx
subroutine allocate_covmat(nex,nlev,nkstar,selective)
integer, intent(in) :: nex,nlev,nkstar
logical, optional :: selective
integer :: ierr
logical :: allocate_all
  allocate_all=.true.;if(present(selective))allocate_all=.not.selective
  call CheckStop(.not.(allocated(ichemobs).and.(allocated(ichemnoobs).or.nchemNoObs==0)),&
     HERE('Allocate COVMAT/UCOVMAT. Frist call set_chemobs_idx'))
  if(.not.allocated(covmat).and.allocate_all)then
   !allocate(covmat(ncorr,nex),stat=ierr)
    allocate(covmat(ncorr,nkstar),stat=ierr)
    call CheckStop(ierr,HERE('Allocate COVMAT'))
    covmat=0e0
  endif
  if(.not.allocated(ucovmat).and.nchemNoObs>0)then
    allocate(ucovmat(nlev*nchemnoobs,nlev*nchemobs,nkstar),stat=ierr)
    call CheckStop(ierr,HERE('Allocate UCOVMAT'))
    ucovmat=0e0
  endif
end subroutine allocate_covmat
!subroutine update_covmat(nxex,nyex,nex,nlev,nconcmlp,nchem,nttot,&
!  concmlp,nmin,nmax,ikstar,nkstar,nstar,w,mt,nt)
subroutine update_covmat(nxex,nyex,nlev,nconcmlp,nchem,nttot,&
  concmlp)
!-----------------------------------------------------------------------
! @description
! Update the spectral error covariance matrix, assuming
! horizontal isotropy
! @author M.Kahnert
!-----------------------------------------------------------------------
  integer nxex,nyex,nlev,nconcmlp,nchem,nttot
  real concmlp(nxex,nyex,nlev,nconcmlp,nchem)

  integer m,n,l1,l2,k1,k2,kk1,kk2,ir,ii,ik,ibox,jbox,i
  integer bindex,bidx
  real dc!,denom,num

  bindex(m,n)=m+((n-1)*(2*nlev*nchemobs-n))/2

  ir=nconcmlp-1
  ii=nconcmlp

  nttot=nttot+1

  do l1=1,nlev
    do k1=1,nchemobs
      kk1=ichemobs(k1)
      ibox=(k1-1)*nlev+l1
      do l2=1,nlev
        do k2=1,nchemobs
          kk2=ichemobs(k2)
          jbox=(k2-1)*nlev+l2
          bidx=bindex(jbox,ibox) ! super-index
!-----------------------------------------------------------------------
! The covariance matrix is symmetric
! => only need to compute the lower triangle of the matrix:
!-----------------------------------------------------------------------
          if(jbox>=ibox) then
!-----------------------------------------------------------------------
! integration over circle in k-space with integration weights w:
!-----------------------------------------------------------------------
            do ik=1,nkstar
              do i=1,nstar(ik)
                m=mt(i,ik)
                n=nt(i,ik)
                dc=concmlp(m,n,l1,ir,kk1)*concmlp(m,n,l2,ir,kk2)&
                  +concmlp(m,n,l1,ii,kk1)*concmlp(m,n,l2,ii,kk2)
                covmat(bidx,ik)=covmat(bidx,ik)+w(i,ik)*dc
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
  enddo
end subroutine update_covmat
!subroutine update_unobs_covmat(nxex,nyex,nex,nlev,nconcmlp,nchem,concmlp,&
!     nmin,nmax,nkstar,nstar,w,mt,nt)
subroutine update_unobs_covmat(nxex,nyex,nlev,nconcmlp,nchem,concmlp)
!-----------------------------------------------------------------------
! @description
! Update that part of the spectral error covariance matrix, that contains
! the cross-covariances between observed and unobserved species
! @author M Kahnert
!-----------------------------------------------------------------------
  integer nxex,nyex,nlev,nconcmlp,nchem
  real concmlp(nxex,nyex,nlev,nconcmlp,nchem)

  integer m,n,l1,l2,k1,k2,kk1,kk2,ir,ii,ik,ibox,jbox,i
 !integer bindex,bidx
  real dc!,denom,num

  if(nchemnoobs<1)then
    if(MasterProc)print dafmt,'WARNING no unobserved to update'
    return
  endif

  ir=nconcmlp-1
  ii=nconcmlp

  do l1=1,nlev
    do k1=1,nchemnoobs
      kk1=ichemnoobs(k1)
      ibox=(k1-1)*nlev+l1
      do l2=1,nlev
        do k2=1,nchemobs
          kk2=ichemobs(k2)
          jbox=(k2-1)*nlev+l2
!-----------------------------------------------------------------------
! integration over circle in k-space with integration weights w:
!-----------------------------------------------------------------------
          do ik=1,nkstar
            do i=1,nstar(ik)
              m=mt(i,ik)
              n=nt(i,ik)
              dc=concmlp(m,n,l1,ir,kk1)*concmlp(m,n,l2,ir,kk2)&
                +concmlp(m,n,l1,ii,kk1)*concmlp(m,n,l2,ii,kk2)
              ucovmat(ibox,jbox,ik)=ucovmat(ibox,jbox,ik)+w(i,ik)*dc
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
end subroutine update_unobs_covmat
!ubroutine write_stat(nex,nxex,nyex,nlev,nkstar,kstar,ikstar,gamma,avgstddev)
subroutine write_stat(nex,nxex,nyex,nlev,gamma,avgstddev)

!-----------------------------------------------------------------------
! @description
! Write various statistical functions and fields to file
! @author M Kahnert
!-----------------------------------------------------------------------
  character(len=*), parameter :: &
    FMT1001="(17e15.4)",&   !"(18e15.4)"
    FMT1002="(22e15.4)",&
    FMT1003="(23e15.4)",&
    FMT1004="(18e15.4)",&   !"(19e15.4)"
    FMT1005="(i4,17e15.4)"  !"(i4,18e15.4)"

  integer :: nex,nxex,nyex,nlev!,nkstar
  real :: gamma(nlev*nchemobs,nex),avgstddev(nlev,nchemobs)!,kstar(nex)

  integer :: ik,l,k,kk,m,n,ibox,jbox,lref,ib,jb
  real :: pi,vcorr(nlev),ccorr(nchemobs),num,denom,lscale(nchemobs)
  real :: sps,PS(nlev,nchemobs,nex),rdim
  character(len=19) :: file=""

  integer bindex,bidx
  bindex(m,n)=m+((n-1)*(2*nlev*nchemobs-n))/2
  pi=acos(-1e0)
!-----------------------------------------------------------------------
! power spectra of spectral densities of horizontal
! correlation functions, normalised to unity:
!-----------------------------------------------------------------------
  do l=1,nlev
    do k=1,nchemobs
      ibox=(k-1)*nlev+l
      sps=0e0
      do ik=1,nkstar
        PS(l,k,ik)=kstar(ik)*gamma(ibox,ik)/float(nex**2)
        sps=sps+PS(l,k,ik)
      enddo
      do ik=1,nkstar
        PS(l,k,ik)=PS(l,k,ik)/sps
      enddo
    enddo
  enddo
!-----------------------------------------------------------------------
! write horizontal auto-correlation spectra
!-----------------------------------------------------------------------
  open(10,file='spectr_horiz_corr_by_lev.dat')
  write(10,*)'spectral horizontal correlations for different levels and chemical species'
  write(10,*)'(one table for each species, one column per table for each level)'
  write(10,*)
  do k=1,nchemobs
    write(10,*)'chemical species ',k,':'
    do ik=1,nkstar
     !write(10,*)kstar(ik),(PS(l,k,ik),l=1,nlev)
      write(10,FMT1003)kstar(ik),(PS(l,k,ik),l=1,nlev)
    enddo
  enddo
  close(10)
  open(10,file='spectr_horiz_corr_by_species.dat')
  write(10,*)'spectral horizontal correlations for different levels and chemical species'
  write(10,*)'(one table for each level, one column per table for each species)'
  write(10,*)
  do l=1,nlev
    write(10,*)'level ',l,':'
    do ik=1,nkstar
      !write(10,*)kstar(ik),(PS(l,k,ik),k=1,nchemobs)
      write(10,FMT1004)kstar(ik),(PS(l,k,ik),k=1,nchemobs)
    enddo
  enddo
  close(10)
!-----------------------------------------------------------------------
! write horizontally averaged standard deviations
!-----------------------------------------------------------------------
  open(10,file='stddev.dat')
  write(10,*)'horizontally averaged standard deviations as a'
  write(10,*)'function of level'
  write(10,*)'(one column for each species)'
  write(10,*)
  do l=1,nlev
   !write(10,*)l,(avgstddev(l,k),k=1,nchemobs)
    write(10,FMT1005)l,(avgstddev(l,k),k=1,nchemobs)
  enddo
  close(10)
!-----------------------------------------------------------------------
! write vertical correlations as a function of horizontal scale:
!-----------------------------------------------------------------------
  lref=min(18,nlev) !reference level for correlations:
!-------

  open(10,file='vertical_corr.dat')
  write(10,*)'vertical covariance wrt reference level l=',lref
  write(10,*)'as a function of wavenumber'
  write(10,*)'(one column for each level, one row for each kstar,'
  write(10,*)'one table for each species)'

  do k=1,nchemobs
    write(10,*)
    write(10,*)'chemical species ',k,':'
    do ik=1,nkstar
      do l=1,nlev
        ibox=(k-1)*nlev+l
        jbox=(k-1)*nlev+lref
        !covmat only contains the lower triangle => use symmetry of matrix:
        if(jbox>=ibox)then
          bidx=bindex(jbox,ibox)
        else
          bidx=bindex(ibox,jbox)
        endif
        vcorr(l)=covmat(bidx,ik)
      enddo
      write(10,FMT1002)(vcorr(l),l=1,nlev)
      !write(10,*)(vcorr(l),l=1,nlev)
    enddo
  enddo
  close(10)
!-----------------------------------------------------------------------
! write chemical correlations for different horizontal scales:
!-----------------------------------------------------------------------
! write correlations at different levels (one file for each wavenumber):
!-------
  do ik=1,nkstar
    if(ik.gt.99999)then
      print*,'ERROR: ik too large'
      return
    endif

    write (file,'(a,"_",i5.5,".",a)')"chem_corr",ik,"dat"
    open(10,file=file)
    write(10,*)'chemical covariance at k=',kstar(ik)
    write(10,*)'(one matrix for each level)'

    do l=1,nlev
      write(10,*)
      write(10,*)'level=',l,':'
      do k=1,nchemobs
        ibox=(k-1)*nlev+l
        do kk=1,nchemobs
          jbox=(kk-1)*nlev+l
          !covmat only contains the lower triangle => use symmetry of matrix:
          if(jbox>=ibox)then
            bidx=bindex(jbox,ibox)
          else
            bidx=bindex(ibox,jbox)
          endif
          ccorr(kk)=covmat(bidx,ik)
        enddo
       !write(10,*)(ccorr(kk),kk=1,nchemobs)
        write(10,FMT1001)(ccorr(kk),kk=1,nchemobs)
      enddo
    enddo
    close(10)
  enddo
!-----------------------------------------------------------------------
! write covmat matrix
!-----------------------------------------------------------------------
  open(10,file='chem_corr_matrix.dat')
  write(10,*)'chemical correlations (nchemobs*nlev)x(nchemobs*nlev)'
  do ik=1,nkstar
    write(10,*)'k=',kstar(ik)
    do ibox=1,nchemobs*nlev
      do k=1,nchemobs
        do l=1,nlev
          jbox=(k-1)*nlev+l
          !covmat only contains the lower triangle => use symmetry of matrix:
          if(jbox>=ibox)then
            bidx=bindex(jbox,ibox)
          else
            bidx=bindex(ibox,jbox)
          endif
          ib=bindex(ibox,ibox)
          jb=bindex(jbox,jbox)
          vcorr(l)=covmat(bidx,ik)/sqrt(covmat(ib,ik)*covmat(jb,ik))
        enddo
        write(10,FMT1002)(vcorr(l),l=1,nlev)
      enddo
    enddo
  enddo
  close(10)
!-----------------------------------------------------------------------
! write horizontal length scales as a function of altitude:
!-----------------------------------------------------------------------
  if(nxex.gt.nyex) then
    rdim=dxdim/1e3
  else
    rdim=dydim/1e3
  endif
  open(10,file='horiz_lengthscale.dat')
  write(10,*)'horizontal length scale as a function of level height'
  write(10,*)'(one column for each species)'
  do l=1,nlev
    do k=1,nchemobs
      ibox=(k-1)*nlev+l
      bidx=bindex(ibox,ibox)
      num=0e0
      denom=0e0
      do m=1,nxex
        do n=1,nyex
! AMVB 2014-11-10
if(n<=nmin(m).or.n>=nmax(m))then
          ik=ikstar(m,n)
          num=num+covmat(bidx,ik)
          denom=denom+covmat(bidx,ik)*kstar(ik)**2
endif
        enddo
      enddo
      lscale(k)=rdim*nex * sqrt(2e0*num/denom) /(2e0*pi)
    enddo
   !write(10,*)l,(lscale(kk),kk=1,nchemobs)
    write(10,FMT1005)l,(lscale(kk),kk=1,nchemobs)
  enddo

end subroutine write_stat
!subroutine normalise_covmat(nex,nxex,nyex,nlev,nttot,nkstar,nstar,kstar,ikstar)
subroutine normalise_covmat(nxex,nyex,nlev,nttot)
!-----------------------------------------------------------------------
! @description
! Normalisation of the spectral background error covariance
! matrix by use of the spectral densities of horizontal error correlations
! @author M Kahnert
!-----------------------------------------------------------------------
  integer nxex,nyex,nlev,nttot!,nex,nkstar
  integer ik,l,k,m,n,ibox,jbox
  real g,twopi!,gamma(nlev*nchemobs,nex),avgstddev(nlev,nchemobs)
  real, allocatable :: gamma(:,:),avgstddev(:,:)

  integer bindex,bidx

  bindex(m,n)=m+((n-1)*(2*nlev*nchemobs-n))/2
  twopi=2e0*acos(-1e0)

  !divide covariances by number of time steps and number of discrete
  !wave vectors in spectral space:
! do ik=1,nkstar
!   do n=1,ncorr
!     covmat(n,ik)=covmat(n,ik)/float(nttot)/twopi
!   enddo
!   do n=1,nchemobs*nlev
!     do m=1,nchemnoobs*nlev
!       ucovmat(m,n,ik)=ucovmat(m,n,ik)/float(nttot)/twopi
!     enddo
!   enddo
! enddo
  covmat (:,:)  =covmat (:,:)  /float(nttot)/twopi
  if(nchemNoObs>0) &
  ucovmat(:,:,:)=ucovmat(:,:,:)/float(nttot)/twopi
  allocate(gamma(nlev*nchemobs,nex),avgstddev(nlev,nchemobs))

  !compute spectral densities gamma of horizontal error correlations:
  do l=1,nlev
    do k=1,nchemobs
      ibox=(k-1)*nlev+l
      bidx=bindex(ibox,ibox) ! super-index for diagonal elements
      g=0e0
      do ik=1,nkstar
        g=g+covmat(bidx,ik)*twopi*kstar(ik)
      enddo
      g=float(nex**2)/g
      avgstddev(l,k)=sqrt(1e0/g)
      do ik=1,nkstar
        gamma(ibox,ik)=g*covmat(bidx,ik)
      enddo
    enddo
  enddo

  !write statistics files:
 !call write_stat(nex,nxex,nyex,nlev,nkstar,kstar,ikstar,gamma,avgstddev)
  call write_stat(nex,nxex,nyex,nlev,                    gamma,avgstddev)

  !convert gamma to sqrt(gamma):
  forall(ik=1:nkstar,ibox=1:nlev*nchemobs) &
    gamma(ibox,ik)=sqrt(gamma(ibox,ik))

  !use gamma to normalise covariance matrix:
  do ik=1,nkstar
    do ibox=1,nlev*nchemobs
      do jbox=ibox,nlev*nchemobs
        bidx=bindex(jbox,ibox)
        if(gamma(ibox,ik)==0e0.or.gamma(jbox,ik)==0e0)then
          write(*,*)'WARNING: gamma=0 in normalise_covmat',ibox,jbox,ik
        else
          covmat(bidx,ik)=covmat(bidx,ik)/(gamma(ibox,ik)*gamma(jbox,ik))
        endif
      enddo
    enddo
  enddo

  !save gamma in file:
  open(10,file='gamma.cov')
  write(10,*) gamma
  close(10)
  !save ucovmat in file:
  open(10,file='ucovmat.cov')
  if(nchemNoObs>0)write(10,*) ucovmat
  close(10)

  deallocate(gamma,avgstddev)
end subroutine normalise_covmat

!subroutine write_covmat(nex,nxex,nyex,nbg,nv,ik,kstar,ikstar,nkstar,vt,lambda)
subroutine write_covmat(nbg,nv,ik,vt,lambda)
!-----------------------------------------------------------------------
! @description
! Write diagonalised covariance matrix and other output
! parameters to output files
! @author M Kahnert
!-----------------------------------------------------------------------
  integer ik,nbg,nv!,nex,nxex,nyex
  real vt(nbg,nv),lambda(nv)
  integer :: k,nvpos,nvneg,nv01,nv02,nv05,nv10!,i,j
  character(len=18) :: file
!-----------------------------------------------------------------------
! Diagnostics of the diagonalisation of the covariance matrix:
!-----------------------------------------------------------------------
  nvpos=nv
  nvneg=0
  nv01=nv
  nv02=nv
  nv05=nv
  nv10=nv
  do k=1,nv
    if(lambda(k)<0.01*lambda(nv))nv01=nv01-1
    if(lambda(k)<0.02*lambda(nv))nv02=nv02-1
    if(lambda(k)<0.05*lambda(nv))nv05=nv05-1
    if(lambda(k)<0.10*lambda(nv))nv10=nv10-1
    if(lambda(k)<0.00)then
      nvneg=nvneg+1
      nvpos=nvpos-1
    endif
  enddo

  write(6,*)'ik=',ik,'  kstar=',kstar(ik)
  write(6,*)

  if(nvneg>0)then
    write(6,*)
    write(6,*)'WARNING: Negative egenvalues found '
    write(6,*)'Number of negative eigenvalues ',nvneg
    write(6,*)'Number of positive eigenvalues ',nvpos
    write(6,*)'Total of eigenvalues ',nv
  endif
  if(float(nvneg)/nv>0.1)then
    write(6,*)'Number of negative egenvalues larger than 10%.'
    write(6,*)'No analysis is made'
    write(6,*)'lambda'
    write(6,*)(lambda(k),k=1,nv)
    return
  endif
  write(6,*)'Number eigenvalues ',nv
  write(6,*)'Number of positive eigenvalues ',nvpos
  write(6,*)'Number of eigenvalues larger than 10 percent of max',nv10
  write(6,*)'Number of eigenvalues larger than 5 percent of max ',nv05
  write(6,*)'Number of eigenvalues larger than 2 percent of max ',nv02
  write(6,*)'Number of eigenvalues larger than 1 percent of max ',nv01

  write(6,*)'Conditioning number: ',&
    maxval(abs(lambda))/max(minval(abs(lambda)),1e-30)

  write(6,*)'lambda'
  write(6,*)(lambda(k),k=1,nv)
!-----------------------------------------------------------------------
! Write \Lambda^{1/2} and X to files, where \Lambda is the diagonal
! matrix of eigenvalues, and X is the orthogonal matrix of eigenvectors:
!-----------------------------------------------------------------------
  lambda(:)=sqrt(max(0e0,lambda(:)))

  call CheckStop(ik>99999,HERE('ik too large'))

  write (file,'(a,"_",i5.5,".",a)')"eigenvec",ik,"cov"
  open(10,file=file)
  write(10,*)vt
  close(10)

  write (file,'(a,"_",i5.5,".",a)')"eigenval",ik,"cov"
  open(11,file=file)
  write(11,*)lambda
  close(11)
end subroutine write_covmat

!subroutine diagonalise_covmat(nex,nxex,nyex,nx,ny,nlev,nchem,&
!           nkstar,kstar,ikstar,nstarmax)
subroutine diagonalise_covmat(nex,nxex,nyex,nx,ny,nlev,nchem,nkstar,nstarmax)
!-----------------------------------------------------------------------
! @description
! Determine eigenvalues and eigenvectors of the
! background error covariance matrix
! @author M Kahnert
!-----------------------------------------------------------------------
  integer :: nex,nxex,nyex,nx,ny,nlev,nchem,nkstar,nstarmax
  logical :: ierrN=.false.,ierrP=.false.
  integer, allocatable :: ifail(:)
  real, allocatable :: vt(:,:),lambda(:)
#ifdef MKL
  real ::accu
  integer(kind=4) :: il,iu,iostat ! for compatibility with mkl_lapack95_lp64
  integer :: nfound,ik
#else
  real :: accu,work(8*nlev*nchemobs)
  integer :: iw(5*nlev*nchemobs),il,iu,nfound,ik,iostat
#endif

  nbg=nlev*nchemobs
  nv=min(nv,nbg)
  il=nbg-nv+1
  iu=nbg
 !accu=2*lamch('S')
#ifdef MKL
  accu=1e0/float(2**22)
#else
 !accu=1e0/float(2**22)
  accu=2*lamch('S')
#endif
 !accu=max(accu,1e-3*DAPREC)<-- TRY THIS IN FUTURE

  allocate(vt(nbg,nv),lambda(nbg),ifail(nbg),stat=iostat)
  call CheckStop(iostat,HERE('Allocate vt,lambda,ifail'))

  do ik=1,nkstar
    vt=0e0
!##call sspev('V','L',nbg,covmat(1,ik),lambda,vt,nbg,work,iostat)
!##call CheckStop(iostat,HERE("ssyev"))
   !write(6,*)'calling sspevx with accuracy ',accu
   !write(6,*)'il,iu=',il,iu
!       [s/d]spevx
!  call dspevx('V','I','L',nbg,covmat(1,ik),0,0,il,iu,&
!        accu,nfound,lambda,vt,nbg,work,iw,ifail,iostat)
!  call CheckStop(iostat,HERE("sspevx"))
   !accu=1e-6
    write(damsg,"(A,'=',E12.3,2(:,',',1X,A,'=',I0))"),&
        'accuracy',accu,'il',il,'iu',iu
    print dafmt,'Valling spevx '//trim(damsg)
#ifdef MKL
    call spevx(covmat(:,ik),lambda,uplo='L',z=vt,il=il,iu=iu,&
               abstol=accu,m=nfound,ifail=ifail,info=iostat)
#else
    call spevx('V','I','L',nbg,covmat(:,ik),0.0,0.0,il,iu,&
               accu,nfound,lambda,vt,nbg,work,iw,ifail,iostat)
#endif
   !where(abs(lambda)<1e-15)lambda=0.0!accu/2
    lambda=ANINT(lambda/(accu*0.5))*accu*0.5  ! <-- NEEDED? -- 2011-11-23
   !lambda=ANINT(lambda*1d15)*1d-15
    if((nfound<nv).or.(iostat>0))then
      if(nfound<nv)&
      print *,"WARNING: only",nfound,"eigenvalues converged"
      if(iostat>0)&
      print *,"WARNING:",iostat,"eigenvectors failed to converge:",ifail(1:iostat)
      print *,"TRY TO REDUCEING 'nv' to",minval((/nv,nfound,nv-iostat/))
    endif
    ierrN=ierrN.or.(iostat<0)
    ierrP=ierrP.or.(iostat>0)
!-----------------------------------------------
!     Write output files and diagnostics
!-----------------------------------------------
   !call write_covmat(nex,nxex,nyex,nbg,nv,ik,kstar,ikstar,nkstar,vt,lambda)
    call write_covmat(nbg,nv,ik,vt,lambda(:nv))
  enddo
  deallocate(vt,lambda,ifail,stat=iostat)
  call CheckStop(iostat,HERE('Deallocate vt,lambda,ifail'))

  open(12,file='dim.tmp')
  write(12,*)nex,nxex,nyex,nbg,nv,nx,ny,nlev,nchem,nchemobs,nstarmax
  close(12)
  call CheckStop(ierrN,HERE("spevx, i-th argument had an illegal value"))
  call CheckStop(ierrP,HERE("spevx, i eigenvectors failed to converge"))
endsubroutine diagonalise_covmat
subroutine read_speccov
!-----------------------------------------------------------------------
! @description
! Read several arrays used for book-keeping of indices
!   in spectral space, assuming isotropic covariance matrix (kstar.tmp)
! Read FFT paramenters (fft.tmp)
! Read normalised spectral background error covariance
!   ucovmat:
!   sqrt_gamma:
!   sqrt_lambda:
!   vt:
! @author AMVB
!-----------------------------------------------------------------------
  integer :: i,j,k,ik,ierr
  real, allocatable :: sq0(:),vt0(:,:)
  character(len=18) file

  open(10,file='dim.tmp',status='OLD',iostat=ierr)
  call CheckStop(ierr,HERE('open dim.tmp'))
  if(debug.and.me==0)print dafmt,'Reading extended domain size from file'
  read(10,*)nex,nxex,nyex,nbg1,nv1,nx,ny,nlev,nchem,nchemobs,nstarmax
  close(10)
  if(debug.and.me==0)print "(11(1X,A,'=',I0))",&
    'nex',nex,'nxex',nxex,'nyex',nyex,'nbg1',nbg1,'nv1',nv1,'nx',nx,'ny',ny,&
    'nlev',nlev,'nchem',nchem,'nchemobs',nchemobs,'nstarmax',nstarmax
  call spec_allocate(selective=.true.)
  open(10,file='kstar.tmp',status='OLD',iostat=ierr)
  call CheckStop(ierr,HERE('open kstar.tmp'))
  read(10,*)nkstar
  read(10,*)ikstar
  read(10,*)kstar
  read(10,*)nmin
  read(10,*)nmax
  read(10,*)mm
  read(10,*)nn
  read(10,*)kx,kxmin,ky,kymin
  close(10)

  call allocate_covmat(nex,nlev,nkstar,selective=.true.)
  open(10,file='ucovmat.cov',status='OLD',iostat=ierr)
  call CheckStop(ierr,HERE('open ucovmat.cov'))
  if(nchemNoObs>0)read(10,*)ucovmat
  close(10)

  allocate(sqrt_gamma(nbg1,nex),stat=ierr)
  call CheckStop(ierr,HERE('Allocation error: SQRT_GAMMA.'))
  open(10,file='gamma.cov',status='OLD',iostat=ierr)
  call CheckStop(ierr,HERE('open gamma.cov'))
  read(10,*)sqrt_gamma
  close(10)

  allocate(stddev(nxex,nyex,nlev,nchem),stat=ierr)
  call CheckStop(ierr,HERE('Allocate STDDEV'))
  open(10,file='stddev.tmp',status='OLD',iostat=ierr)
  call CheckStop(ierr,HERE('open stddev.tmp'))
  read(10,*)stddev
  close(10)

  open(10,file='fft.tmp',status='OLD',iostat=ierr)
  read(10,*)lensav
  allocate(wsave(lensav),stat=ierr)
  call CheckStop(ierr,HERE('Allocate WSAVE'))
  read(10,*)wsave
  close(10)

  allocate(vt(nbg1,nv1,nex),sqrt_lambda(nv1,nex),&
           vt0(nbg1,nv1),sq0(nv1),stat=ierr)
  call CheckStop(ierr,HERE('Allocate VT,..'))

  call CheckStop(nkstar>99999,HERE('nkstar too large'))
  do ik=1,nkstar ! start loop over wavenumbers
    write(file,'(A,"_",I5.5,".",A)')"eigenvec",ik,"cov"
    open(10,file=file,status='OLD',iostat=ierr)
    call CheckStop(ierr,HERE('open '//trim(file)))
    read(10,*)vt0
    close(10)
    vt(:,:,ik)=vt0(:,:)

    write(file,'(A,"_",I5.5,".",A)')"eigenval",ik,"cov"
    open(10,file=file,status='OLD',iostat=ierr)
    call CheckStop(ierr,HERE('open '//trim(file)))
    read(10,*)sq0
    close(10)
    sqrt_lambda(:,ik)=sq0(:)
  enddo
  deallocate(vt0,sq0,stat=ierr)
  call CheckStop(ierr,HERE('Deallocate VT0,SQ0'))

end subroutine read_speccov
! function bindex(k1,l1,k2,l2,upper_only) result(bidx)
!   integer, intent(in) :: k1,l1,k2,l2
!   logical, intent(in) :: upper_only
!   integer :: bidx,b1,b2,m,n
!   logical :: upp
!   b1=(k1-1)*nlev+l1
!   b2=(k2-1)*nlev+l2
!   upp=.true.;if(present(upper_only))upp=upper_only.and.b1>=b2
!   if(upp)then
!     m=max(b1,b2)
!     n=min(b1,b2)
!     bidx=m+((n-1)*(2*nlev*nchemobs-n))/2
!   else
!     bidx=-1
!   endif
! end function bindex
end module covmat_ml
