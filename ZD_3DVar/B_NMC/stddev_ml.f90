MODULE stddev_ml
use CheckStop_ml, only : CheckStop
use spectralcov,  only : stddev
implicit none
integer, parameter :: nHH=6, HH0=0, HH1=24/nHH-1 ! HH1=23
real, parameter :: KTRUNC=0e0
real, dimension(:,:,:,:,:), allocatable :: bias,variance
real, dimension(:,:),       allocatable :: nt0
logical, parameter, private :: debug_ij=.true.
!nteger, parameter, private :: debug_i=10, debug_j=14
integer, parameter, private :: debug_i=21, debug_j=32, debug_k=4, debug_c=1
CONTAINS
subroutine allocate_stddev(nxex,nyex,nlev,nchem)
implicit none
integer, intent(in) :: nxex,nyex,nlev,nchem
integer :: ierr
  if(.not.allocated(bias))then
    allocate(bias    (nxex,nyex,nlev,nchem,HH0:HH1),stat=ierr)
    call CheckStop(ierr,'Allocation error: BIAS.')
    bias=0e0
  endif
  if(.not.allocated(variance))then
    allocate(variance(nxex,nyex,nlev,nchem,HH0:HH1),stat=ierr)
    call CheckStop(ierr,'Allocation error: VARIANCE.')
    variance=0e0
  endif
  if(.not.allocated(stddev))then
    allocate(stddev  (nxex,nyex,nlev,nchem),        stat=ierr)
    call CheckStop(ierr,'Allocation error: STDDEV.')
    stddev=0e0
  endif
  if(.not.allocated(nt0))then
    allocate(nt0     (nchem,HH0:HH1),               stat=ierr)
    call CheckStop(ierr,'Allocation error: NT0.')
    nt0=0e0
  endif
end subroutine allocate_stddev
subroutine update_stddev(nxex,nyex,nlev,nconcmlp,nchem,k,HOUR,concmlp)
!-----------------------------------------------------------------------
!     @description
!     Add up field differences at times t and t+dt_rochon at
!     every time step; count the number of time steps.
!     @author M Kahnert
!-----------------------------------------------------------------------
implicit none

integer, intent(in) :: nxex,nyex,nlev,nconcmlp,nchem,k,hour
real,    intent(in) :: concmlp(nxex,nyex,nlev,nconcmlp)

real epsilon
integer i,j,ilev,it

  if(.not.allocated(bias)    .or.&
     .not.allocated(variance).or.&
     .not.allocated(nt0)) call allocate_stddev(nxex,nyex,nlev,nchem)

  it=hour/nHH
  nt0(k,it)=nt0(k,it)+1e0
!-----------------------------------------------------------------------
!     Update bias and variance:
!     This subroutine performs the adding of the fields at the current
!     time step. The division by nt0 and the conversion from variance
!     to standard deviation will later be performed by the
!     subroutine normalise_stddev.
!
!     NOTE: To remove diurnal and semi-diurnal signals, fields are
!     grouped according to the time of the day HOUR, and the mean for
!     each HOUR is subtracted. Thus
!     variance = 1/Ntbound sum_HOUR {
!     1/nt0(HOUR) sum_j epsilon_{HOUR,j}**2 -
!     [ 1/nt0(HOUR) sum_j epsilon_{HOUR,j} ]**2 }.
!-----------------------------------------------------------------------
  do ilev=1,nlev
    do j=1,nyex
      do i=1,nxex
        epsilon=CONCMLP(i,j,ilev,2)-CONCMLP(i,j,ilev,1)
        bias    (i,j,ilev,k,it)=bias    (i,j,ilev,k,it)+epsilon
        variance(i,j,ilev,k,it)=variance(i,j,ilev,k,it)+epsilon**2
      enddo
    enddo
  enddo
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  if(debug_ij)then
    write(*,*)'xxx1 b,s:',bias(debug_i,debug_j,max(nlev-debug_k,1),min(debug_c,nchem),it),&
                      variance(debug_i,debug_j,max(nlev-debug_k,1),min(debug_c,nchem),it)
  endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
end subroutine update_stddev
subroutine normalise_stddev(nxex,nyex,nlev,nchem)
!-----------------------------------------------------------------------
!     @description
!     Normalise variance and bias fields by dividing by the
!     number of time steps. Subtract biases from variances
!     and convert to standard deviation.
!     @author M Kahnert
!-----------------------------------------------------------------------
implicit none

integer, intent(in) :: nxex,nyex,nlev,nchem

integer i,j,k,ilev,it
real sigmamin,nt
character(len=18) file

  call CheckStop(.not.allocated(bias),    'Variable not allocated: BIAS.')
  call CheckStop(.not.allocated(variance),'Variable not allocated: VARIANCE.')
  call CheckStop(.not.allocated(nt0),     'Variable not allocated: NT0.')
  if(.not.allocated(stddev))call allocate_stddev(nxex,nyex,nlev,nchem)

  stddev=0e0
  sigmamin=9e20

  do k=1,nchem
    do ilev=1,nlev
      do j=1,nyex
        do i=1,nxex
          nt=0e0
          do it=HH0,HH1
            if(nt0(k,it).gt.0.0)then
              nt=nt+1e0
  !     if(abs(nt0(k,it)-1e0).lt.1e-8.or.abs(nt0(k,it)).lt.1e-8)then
  !     write(*,*)'WARNING: too few data points for determining'
  !     write(*,*)'standard deviation and bias for k=',k
  !     endif
              bias(i,j,ilev,k,it)=bias(i,j,ilev,k,it)/nt0(k,it)
              variance(i,j,ilev,k,it)=variance(i,j,ilev,k,it)/nt0(k,it)
              stddev(i,j,ilev,k)=stddev(i,j,ilev,k)+abs(variance(i,j,ilev,k,it)-bias(i,j,ilev,k,it)**2)
            end if
          enddo
          stddev(i,j,ilev,k)=sqrt(stddev(i,j,ilev,k)/nt)
          if(stddev(i,j,ilev,k).gt.0e0.and.stddev(i,j,ilev,k).lt.sigmamin)sigmamin=stddev(i,j,ilev,k)
        enddo
      enddo
    enddo
  enddo
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  if(debug_ij)then
    write(*,*)'xxx2 b,s:',  bias(debug_i,debug_j,max(nlev-debug_k,1),min(debug_c,nchem),0),&
                        variance(debug_i,debug_j,max(nlev-debug_k,1),min(debug_c,nchem),0)
    write(*,*)'xxx2a sd:',stddev(debug_i,debug_j,max(nlev-debug_k,1),min(debug_c,nchem))
  endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  do k=1,nchem
    do ilev=1,nlev
      call spectral_filter(ilev,k,KTRUNC)
      do j=1,nyex
        do i=1,nxex
          if(stddev(i,j,ilev,k).le.0e0)then
  !    write(*,*)'WARNING: standard deviation <= 0'
            stddev(i,j,ilev,k)=sigmamin
          endif
        enddo
      enddo
    enddo
  enddo
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  if(debug_ij)then
    write(*,*)'xxx2b sd:',stddev(debug_i,debug_j,max(nlev-debug_k,1),min(debug_c,nchem))
  endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!     write standard deviations to temporary file
  open(10,file='stddev.tmp')
  write(10,*)stddev
  close(10)
end subroutine normalise_stddev
subroutine get_moderr(nxex,nyex,nlev,nconcmlp,nchem,k,HOUR,concmlp)
!-----------------------------------------------------------------------
! @description
! Compute the model errors in physical space and normalise
! them by dividing by the standard deviation
! @author M.Kahnert
!-----------------------------------------------------------------------
implicit none
integer, intent(in) :: nxex,nyex,nlev,nconcmlp,nchem,k,hour
real,  intent(inout):: concmlp(nxex,nyex,nlev,nconcmlp)
real epsilon
integer i,j,ilev,it
!-----------------------------------------------------------------------
! The error epsilon is defined as epsilon= ( x-x' ), where x' is some perturbed field
! The normalised errors are stored in CONCMLP(i,j,ilev,nconcmlp-1):
!-----------------------------------------------------------------------
  call CheckStop(k<0.or.k>nchem,'get_moderr: index k out of range')
  it=hour/nHH
  do ilev=1,nlev
    do j=1,nyex
      do i=1,nxex
        epsilon=CONCMLP(i,j,ilev,2)-CONCMLP(i,j,ilev,1)-bias(i,j,ilev,k,it)
        CONCMLP(i,j,ilev,1)=epsilon/stddev(i,j,ilev,k)
      enddo
    enddo
  enddo
end subroutine get_moderr
END MODULE stddev_ml
