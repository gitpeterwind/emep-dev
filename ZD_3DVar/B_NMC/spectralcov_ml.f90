module spectralcov
use ModelConstants_ml,only: PPT,PPTINV
use CheckStop_ml,     only: CheckStop
implicit none
!real, parameter :: FGSCALE=PPBINV,FGSCALE_INV=PPB
!real, parameter :: FGSCALE=PPB,FGSCALE_INV=PPBINV
 real, parameter :: FGSCALE=1e9,FGSCALE_INV=1e-9
 real, parameter :: DAPREC=PPT*FGSCALE_INV !1e0/2**22; 5.96e-08=2**-24 (single precission)

integer, save :: nchem,nchemobs,nchemnoobs,nx,ny,nlev
integer, save :: nex,nxex,nyex,nbg1,nv1
integer, save :: nkstar,kx,kxmin,ky,kymin,lensav
integer, dimension(:,:), allocatable, save :: ikstar
integer, dimension(:),   allocatable, save :: nmin,nmax
integer, dimension(:),   allocatable, save :: mm,nn
integer, dimension(:),   allocatable, save :: ichemInv,ichemObs,ichemNoObs
real, dimension(:),      allocatable, save :: kstar,wsave
real, dimension(:,:),    allocatable, save :: sqrt_gamma
real, dimension(:,:,:,:),allocatable, save :: stddev
real, dimension(:,:,:),  allocatable, save :: vt
real, dimension(:,:),    allocatable, save :: sqrt_lambda
real, dimension(:,:,:),  allocatable, save :: ucovmat
real, dimension(:,:),    allocatable, save :: covmat

integer :: nstarmax,nttot,ncorr
integer, dimension(:),   allocatable :: nstar
integer, dimension(:,:), allocatable :: mt, nt
real, dimension(:,:),    allocatable :: intweight
!#define DEBUG_B

contains
subroutine spec_allocate(selective)
logical, optional :: selective
integer :: ierr
logical :: allocate_all
  allocate_all=.true.;if(present(selective))allocate_all=.not.selective
  if(.not.allocated(intweight).and.allocate_all)then
    allocate(intweight(nxex*nyex,nex),stat=ierr)
    call CheckStop(ierr,'Allocation error: INTWEIGHT.')
    intweight=0e0
  endif
  if(.not.allocated(nstar).and.allocate_all)then
    allocate(nstar(nex),stat=ierr)
    call CheckStop(ierr,'Allocation error: NSTAR.')
    nstar=0
  endif
  if(.not.allocated(mt).and.allocate_all)then
    allocate(mt(nex,nex),stat=ierr)
    call CheckStop(ierr,'Allocation error: MT.')
    mt=0
  endif
  if(.not.allocated(nt).and.allocate_all)then
    allocate(nt(nex,nex),stat=ierr)
    call CheckStop(ierr,'Allocation error: NT.')
    nt=0
  endif
  if(.not.allocated(kstar))then
    allocate(kstar(nex),stat=ierr)
    call CheckStop(ierr,'Allocation error: KSTAR.')
    kstar=0e0
  endif
  if(.not.allocated(nmin))then
    allocate(nmin(nxex),stat=ierr)
    call CheckStop(ierr,'Allocation error: NMIN.')
    nmin=0
  endif
  if(.not.allocated(nmax))then
    allocate(nmax(nxex),stat=ierr)
    call CheckStop(ierr,'Allocation error: NMAX.')
    nmax=0
  endif
  if(.not.allocated(ikstar))then
    allocate(ikstar(nxex,nyex),stat=ierr)
    call CheckStop(ierr,'Allocation error: IKSTAR.')
    ikstar=0
  endif
  if(.not.allocated(mm))then
    allocate(mm(-nxex:nxex),stat=ierr)
    call CheckStop(ierr,'Allocation error: MM.')
    mm=0
  endif
  if(.not.allocated(nn))then
    allocate(nn(-nyex:nyex),stat=ierr)
    call CheckStop(ierr,'Allocation error: NN.')
    nn=0
  endif
end subroutine spec_allocate
subroutine write_kstar(nstarmax)
!-----------------------------------------------------------------------
! from MATCH-3DVar
! @description
! Initialise several arrays used for book-keeping of indices
! in spectral space, assuming isotropic covariance matrix
! @author M.Kahnert
!-----------------------------------------------------------------------
implicit none
integer, intent(in) :: nstarmax
integer :: m,n,m1,n1,ik,i,j,k,l,l1,j1,j0
real :: pi,pi2,twopi,dy,dx!,sum
real :: theta(nstarmax+1)
#ifdef DEBUG_B
  character(len=18) :: file
#endif
  pi=acos(-1e0)
  pi2=pi/2e0
  twopi=2e0*pi

!-----------------------------------------------------------------------
! For each wavenumber kstar, find the set of wavevector components
! (m,n) with corresponding angles theta=arctan[(n/Nyex)/(m/Nxex)]
! in ascending order of theta
!-----------------------------------------------------------------------
  theta=0e0
  intweight=0e0
  do ik=1,Nkstar !loop over wavenumbers
    i=0
    !upper right quadrant of wavevector-circle
    do n=0,Ky
      n1=nn(n)
      dy=float(n)/float(Nyex)
      do m=Kx,0,-1
        m1=mm(m)
        if(n==0.and.m==0)then
          if(ik==1)then
            i=i+1
            mt(i,ik)=m1
            nt(i,ik)=n1
            theta(i)=0e0
          endif
        else
          dx=float(m)/float(Nxex)
          if(ikstar(m1,n1)==ik.and.(n1<=nmin(m1).or.n1>=nmax(m1)))then
            i=i+1
            mt(i,ik)=m1
            nt(i,ik)=n1
!           if(m==0)then
!             theta(i)=pi2
!           else
!             theta(i)=atan(dy/dx)
!           endif
            theta(i)=atan2(dy,dx)
          endif
        endif
      enddo
    enddo
    !upper left quadrant
    do n=Ky,0,-1
      n1=nn(n)
      dy=float(n)/float(Nyex)
      do m=-1,-Kxmin,-1
        m1=mm(m)
        dx=float(m)/float(Nxex)
        if(ikstar(m1,n1)==ik.and.(n1<=nmin(m1).or.n1>=nmax(m1)))then
          i=i+1
          mt(i,ik)=m1
          nt(i,ik)=n1
!         theta(i)=pi+atan(dy/dx)
          theta(i)=atan2(dy,dx)
        endif
      enddo
    enddo
    !lower left quadrant
    do n=-1,-Kymin,-1
      n1=nn(n)
      dy=float(n)/float(Nyex)
      do m=-Kxmin,-1
        m1=mm(m)
        dx=float(m)/float(Nxex)
        if(ikstar(m1,n1)==ik.and.(n1<=nmin(m1).or.n1>=nmax(m1)))then
          i=i+1
          mt(i,ik)=m1
          nt(i,ik)=n1
!         if(m==0)then
!           theta(i)=3*pi2
!         else
!           theta(i)=pi+atan(dy/dx)
!         endif
          theta(i)=atan2(dy,dx)
       endif
      enddo
    enddo
    !lower right quadrant
    do n=-Kymin,-1
      n1=nn(n)
      dy=float(n)/float(Nyex)
      do m=0,Kx
        if(n==0.and.m==0)then
        else
          m1=mm(m)
          dx=float(m)/float(Nxex)
          if(ikstar(m1,n1)==ik.and.(n1<=nmin(m1).or.n1>=nmax(m1)))then
            i=i+1
            mt(i,ik)=m1
            nt(i,ik)=n1
!           theta(i)=twopi+atan(dy/dx)
            theta(i)=atan2(dy,dx)
          endif
        endif
      enddo
    enddo
    !last point = first point (to close the circle)
    i=i+1
    theta(i)=twopi
!-----------------------------------------------------------------------
! sort angles in ascending order
!-----------------------------------------------------------------------
   !call quicksort(theta,1,i,mt(:,ik),nt(:,ik))
    call Qsort(theta(:i),mt(:i,ik),nt(:i,ik))
!-----------------------------------------------------------------------
! determine trapezoid weights for isotropic integration
!-----------------------------------------------------------------------
    !exeptional point (m,n)=(0,0):
    if(ik==1)then
      j0=1
    else
      j0=0
    endif

    k=1
    if(theta(2+j0)/=theta(1+j0))then
      intweight(1+j0,ik)=0.5e0*(theta(2+j0)-theta(1+j0)) ! OBS: Th is in radians
    else
      do while(theta(1+j0)==theta(1+j0+k))
        k=k+1
      enddo
      do j=1+j0,k+j0
        intweight(j,ik)=0.5e0*(theta(1+j0+k)-theta(1+j0))/float(k)
      enddo
    endif

    l=1
    if(theta(i)/=theta(i-1))then
      intweight(i,ik)=0.5e0*(theta(i)-theta(i-1))
    else
      do while(theta(i)==theta(i-l))
        l=l+1
      enddo
      do j=i-l+1,i
        intweight(j,ik)=0.5e0*(theta(i)-theta(i-l))/float(l)
      enddo
    endif

    j1=1
    do j=k+1,i-l
      if(j>=j1)then
        if(theta(j+1)/=theta(j))then
          intweight(j,ik)=0.5e0*(theta(j+1)-theta(j-1))
        else
          l1=1
          do while(theta(j)==theta(j+l1))
            l1=l1+1
          enddo
          intweight(j,ik)=0.5e0*(theta(j+l1)-theta(j-1))/float(l1)
          j1=j+l1
        endif
      else !hop over equal theta-values
        intweight(j,ik)=intweight(j-1,ik)
      endif
    enddo
!-------------------------------------------
! close the integration circle by adding the first and last
! integration weight, i.e.
! intweight1=0.5*(th_2-th_1)+0.5(th_n-th_n-1)=0.5*(th2-th_n1), since th_n=th_1
!-------------------------------------------
    do j=1+j0,k+j0
      intweight(j,ik)=intweight(j,ik)+intweight(i,ik)/float(k)
    enddo
!-------------------------------------------
!   The special point (m,n)=(0,0) contributes
!   isotropically for ik=1
!-------------------------------------------
    if(ik==1)then
      do j=2,nstar(1)
        intweight(j,1)=intweight(j,1)*float(nstar(1)-1)/float(nstar(1))
      enddo
      intweight(1,1)=twopi/float(nstar(1))
    endif
#ifdef DEBUG_B
  write (file,'(a,"_",i5.5,".",a)')"w",ik,"aux"
  open(11,file=file)
  write(11,*)nxex,nyex,nex,ik,nstar(ik),nstarmax
  if(nstar(ik)>0) write(11,*)intweight(:nstar(ik),ik)
  close(11)
  write (file,'(a,"_",i5.5,".",a)')"t",ik,"aux"
  open(11,file=file)
  write(11,*)ik,i,nstarmax+1
  if(nstar(i)>0) write(11,*)theta(:i)
  close(11)
#endif
  enddo ! end of loop over wavenumbers
!-----------------------------------------------------------------------
!     Store results in file:
!-----------------------------------------------------------------------
  open(13,file='kstar.tmp')
  write(13,*)nkstar
  write(13,*)ikstar
  write(13,*)kstar
  write(13,*)nmin
  write(13,*)nmax
  write(13,*)mm
  write(13,*)nn
  write(13,*)kx,kxmin,ky,kymin
  close(13)
  contains
! http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
  RECURSIVE SUBROUTINE Qsort(a,m,n)
    REAL,    INTENT(INOUT) :: a(:)
    INTEGER, INTENT(INOUT) :: m(:),n(:)
    INTEGER :: split
    IF(size(a) > 1) THEN
      CALL Partition(a,m,n,split)
      CALL Qsort(a(:split-1),m(:split-1),n(:split-1))
      CALL Qsort(a(split:)  ,m(split:)  ,n(split:))
    ENDIF
  END SUBROUTINE Qsort
  SUBROUTINE Partition(a,m,n,marker)
    REAL,    INTENT(INOUT) :: a(:)
    INTEGER, INTENT(INOUT) :: m(:),n(:)
    INTEGER, INTENT(OUT)   :: marker
    INTEGER :: left, right, itemp
    REAL :: pivot, rtemp

    pivot = (a(1)+a(size(a)))/2  ! Average of first and last elements to prevent quadratic
    left  = 0                    ! behavior with sorted or reverse sorted data
    right = size(a)+1
    DO WHILE(left < right)
      right = right-1
      DO WHILE(a(right) > pivot)
        right = right-1
      ENDDO
      left = left+1
      DO WHILE(a(left) < pivot)
        left = left+1
      ENDDO
      IF(left < right) THEN
        rtemp=a(left);a(left)=a(right);a(right)=rtemp
        itemp=m(left);m(left)=m(right);m(right)=itemp
        itemp=n(left);n(left)=n(right);n(right)=itemp
      ENDIF
    ENDDO

    IF(left == right) THEN
      marker = left+1
    ELSE
      marker = left
    ENDIF
  END SUBROUTINE Partition
end subroutine write_kstar
subroutine initfft()
!-----------------------------------------------------------------------
! @description
! Initialise FFT
! @author AMVB
!-----------------------------------------------------------------------
implicit none
integer :: ierr
  lensav=3*NXEX+3*NYEX+8
  allocate(wsave(lensav),stat=ierr)
  call CheckStop(ierr,'Allocation error: WSAVE.')
  call cfft2i(nxex,nyex,wsave,lensav,ierr)
  call CheckStop(ierr,'ERROR IN CFFT2I')
  open(10,file='fft.tmp');write(10,*)lensav;write(10,*)wsave;close(10)
end subroutine initfft
subroutine initspec()
!-----------------------------------------------------------------------
! from MATCH-3DVar
! @description
! Initialise several arrays used for book-keeping of indices
! in spectral space, assuming isotropic covariance matrix
! @author M Kahnert
!-----------------------------------------------------------------------
implicit none

integer m,n,ntrunc!,Kx,Ky,Kxmin,Kymin
!real theta(nxex*nyex+1)
integer m1,n1,ik!,i,j,k,l,l1,j1,j0
real rkstar, dkstar,kstar_min!,pi,pi2,twopi,dy,dx,sum

! pi=acos(-1e0)
! pi2=pi/2e0
! twopi=2e0*pi
  dkstar=3e0
  call spec_allocate()

!-----------------------------------------------------------------------
! convert Fourier indices mm=1,...,Kx' to m=-Kxmin,...,Kx
!-----------------------------------------------------------------------
  mm(:)=0
  nn(:)=0
  if(mod(nxex,2)==0)then
    Kx=nxex/2
    Kxmin=Kx-1
  else
    Kx=(nxex-1)/2
    Kxmin=Kx
  endif
  if(mod(nyex,2)==0)then
    Ky=nyex/2
    Kymin=Ky-1
  else
    Ky=(nyex-1)/2
    Kymin=Ky
  endif

  do m=-Kxmin,Kx
    if(m<0)then
      mm(m)=m+nxex+1
    else
      mm(m)=m+1
    endif
  enddo
  do n=-Kymin,Ky
    if(n<0)then
      nn(n)=n+nyex+1
    else
      nn(n)=n+1
    endif
  enddo

!-----------------------------------------------------------------------
! elliptical truncation,
! such that sqrt(m**2/Kx**2 + n**2/Ky**2) <= 1 :
! nn(n)=1,...,nmin,nmax,nmax+1,...,Nyex
! where nmin=Ntrunc+1, nmax=Nyex-Ntrunc+1,
! Ntrunc(m) = Ky*sqrt(1-m**2/Kx**2)
!-----------------------------------------------------------------------
  do m=-Kxmin,Kx
    ntrunc=int(float(Ky)*sqrt(1e0-float(m**2)/float(Kx**2)))
    nmin(mm(m))=ntrunc+1
    nmax(mm(m))=max(1,-ntrunc+nyex+1)
  enddo

!-----------------------------------------------------------------------
! Let rkstar=nex * sqrt(m**2/Kx**2 + n**2/Ky**2). Further, subdivide
! the interval from kstar_min to kstar_max into subintervals
! I(1),...,I(nmax) of width dkstar. Then ikstar(mm,nn) denotes the
! number of the subinterval to which the pair (mm,nn) belongs, i.e.
! nex * sqrt(m**2/Kx**2 + n**2/Ky**2) \in I(ikstar).
!-----------------------------------------------------------------------
  kstar_min=0e0
! kstar_min=nex*sqrt(1e0/float(Kx)**2+1e0/float(Ky)**2)

  Nkstar=-9999

! do ik=1,nex
!   Nstar(ik)=0
! enddo
  Nstar(:)=0

!-----------------------------------------------------------------------
!     rkstar:       wavenumber (saved in kstar)
!     ikstar(m,n): integer index of wavenumber belonging to wavevector
!                  components (m,n), i.e. rkstar=kstar(ikstar(m,n))
!     Nkstar:      max number of wavenumbers
!     Nstar(i):    max number of wavevectors belonging to wavenumber
!                  kstar(i)
!     Nstarmax:    maximum of Nstar array, i.e. max number of
!                  wavevectors per wavenumber
!-----------------------------------------------------------------------
  do m=-Kxmin,Kx
    m1=mm(m)
    do n=-Kymin,Ky
      n1=nn(n)
      rkstar=nex*sqrt((float(m)/float(kx))**2+(float(n)/float(ky))**2)
      ikstar(m1,n1)=int((rkstar-kstar_min)/dkstar+1e0)
      if(n1<=nmin(m1).or.n1>=nmax(m1))then
        Nstar(ikstar(m1,n1))=Nstar(ikstar(m1,n1))+1
        if(ikstar(m1,n1)>Nkstar)Nkstar=ikstar(m1,n1)
      endif
    enddo
  enddo

  do m=1,nex
    kstar(m)=kstar_min+(float(m)-0.5e0)*dkstar
  enddo

! Nstarmax=-9999
! do m=1,Nkstar
!   if(Nstar(m).gt.Nstarmax) Nstarmax=Nstar(m)
! enddo
  Nstarmax=maxval(Nstar(:Nkstar))

  call write_kstar(nstarmax)

end subroutine initspec
end module spectralcov
