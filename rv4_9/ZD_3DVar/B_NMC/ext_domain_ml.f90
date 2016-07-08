module exd_domain_ml
use CheckStop_ml, only : CheckStop
#ifndef OLD_SPL
use splines, only : surf1=>fitp_surf1,surf2=>fitp_surf2
#endif
implicit none
contains
SUBROUTINE EXT_DOMAIN(NX,NY,NLEV,NXEX,NYEX,FIELD,EXT_FIELD)
!-----------------------------------------------------------------------
! from MATCH-3DVar
! @description
! Extrapolate the state vector to an extended domain in 2d with periodic
! boundary conditions as a preparatory step to a 2d-FFT from physical
! to spectral space
! @author M.Kahnert
!-----------------------------------------------------------------------
implicit none
integer, intent(in) :: nx, ny, nlev, nxex, nyex
real,    intent(in) :: field(nx,ny,nlev)
real,    intent(out):: ext_field(nxex,nyex,nlev)

real x(nx+1),y(ny+1),zx1(ny+1),zxm(ny+1),zy1(nx+1),zyn(nx+1),&
     zxy11,zxym1,zxy1n,zxymn,zp(nx+1,ny+1,3),temp(2*nxex+nyex),&
     sigma,xx(nx+1:nxex),yy(ny+1:nyex)
real z(nx+1,ny+1)
integer i,j,l,ierr

!real amax,amin,scale
!data scale/1e0/
data sigma/0.1e0/ ! tension factor of spline interpolation

#ifdef OLD_SPL
real surf2
#endif

  do i=1,nx
    x(i)=float(i)
  enddo
  do j=1,ny
    y(j)=float(j)
  enddo
  x(nx+1)=nxex
  y(ny+1)=nyex

  do i=nx+1,nxex
    xx(i)=i
  enddo
  do j=ny+1,nyex
    yy(j)=j
  enddo


!-----------------------------------------------------------------------
!     Loop over layers
!-----------------------------------------------------------------------
  do l=1,nlev
    do i=1,nx
      do j=1,ny
!       field(i,j,l)=scale*field(i,j,l)
        ext_field(i,j,l)=field(i,j,l)
      enddo
    enddo
    do i=1,nx
      ext_field(i,nyex,l)=field(i,1,l) ! northern=southern boundary
    enddo
    do j=1,ny
      ext_field(nxex,j,l)=field(1,j,l) ! eastern=western boundary
    enddo
    ext_field(nxex,nyex,l)=field(1,1,l) ! NE=SW corner
!-----------------------------------------------------------------------
!     Step 1: determine partial derivatives at domain boundary
!-----------------------------------------------------------------------

    call surf1(nx,ny,x,y,field(:,:,l),nx,zx1,zxm,zy1,zyn,zxy11,&
               zxym1,zxy1n,zxymn,255,zp,temp,sigma,ierr)
    call CheckStop(ierr,'ERROR1 in ext_domain returned from surf1')

!     Copy field values from model domain to z:

    do i=1,nx
      do j=1,ny
        z(i,j)=field(i,j,l)
      enddo
    enddo
!     periodic continuation:
    do i=1,nx
      z(i,ny+1)=z(i,1)
      zyn(i)=zy1(i)
    enddo
    do j=1,ny
      z(nx+1,j)=z(1,j)
      zxm(j)=zx1(j)
    enddo
    z(nx+1,ny+1)=z(1,1)
    zyn(nx+1)=zyn(1)
    zxm(ny+1)=zxm(1)

    zxym1=zxy11
    zxy1n=zxy11
    zxymn=zxy11

!-----------------------------------------------------------------------
!     Step 2: interpolation of extended domain:
!-----------------------------------------------------------------------
    call surf1(nx+1,ny+1,x,y,z,nx+1,zx1,zxm,zy1,zyn,&
               zxy11,zxym1,zxy1n,zxymn,0,zp,temp,sigma,ierr)
    call CheckStop(ierr,'ERROR2 in ext_domain returned from surf1')
    do i=nx+1,nxex-1
      do j=1,ny
        ext_field(i,j,l)=surf2(xx(i),y(j),nx+1,ny+1,x,y,z,nx+1,zp,sigma)
      enddo
    enddo

    do i=nx+1,nxex
      do j=ny+1,nyex
        ext_field(i,j,l)=surf2(xx(i),yy(j),nx+1,ny+1,x,y,z,nx+1,zp,sigma)
      enddo
    enddo
    do i=1,nx
      do j=ny+1,nyex-1
        ext_field(i,j,l)=surf2(x(i),yy(j),nx+1,ny+1,x,y,z,nx+1,zp,sigma)
      enddo
    enddo

  enddo ! end of loop over layers
END SUBROUTINE EXT_DOMAIN

SUBROUTINE EXT_DOMAIN_INV(NX,NY,NLEV,NXEX,NYEX,CONCMLP,MLPFIELD)
!-----------------------------------------------------------------------
! from MATCH-3DVar
! @description
! cutout model domain from extended domain
! @author M.Kahnert
!-----------------------------------------------------------------------
implicit none
integer, intent(in) :: nx,ny,nlev,nxex,nyex
real,    intent(out):: mlpfield(nx,ny,nlev)
real,    intent(in) :: concmlp(nxex,nyex,nlev)
integer i,j,ilev

  do i=1,nx
    do j=1,ny
      do ilev=1,nlev
        mlpfield(i,j,ilev)=concmlp(i,j,ilev)
      enddo
    enddo
  enddo

END SUBROUTINE EXT_DOMAIN_INV
end module exd_domain_ml
