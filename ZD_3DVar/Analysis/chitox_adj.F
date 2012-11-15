C---- Comment Start --------------------------------
C
C     xtochi
C
C     @shortdescr Adjoint to transformation from x to chi
C
C     @keywords
C
C     @description
C     Adjoint of the transformation chitox
C
C     @author M Kahnert
C
C
C
C---- Comment End ---------------------------------
      subroutine chitox_adj(chi,x)
      use spectralcov
      implicit none

      real x(nxex,nyex,nlev,nchemobs)
      complex chi(nv1,nxex,nyex)

      integer i,j,k,kk,l,ik,m,n,ir,iostat,ibox
      real reta(nxex,nyex),ieta(nxex,nyex),nxy

      nxy=float(nxex*nyex)

      chi=cmplx(0e0,0e0)
      do k=1,nchemobs
        kk=ichemobs(k)
        do l=1,nlev
          ibox=(k-1)*nlev+l
          do i=1,nxex
            do j=1,nyex
c-----------------------------------------------------------------------
c     Multiply from the left with S, where S is a diagonal matrix
c     containing the standard deviations sqrt(sigma):
c-----------------------------------------------------------------------
              reta(i,j)=x(i,j,l,k) *stddev(i,j,l,kk)
              ieta(i,j)=0e0
            enddo
          enddo
c-----------------------------------------------------------------------
c     Apply adjoint inverse FFT, and use that
c     F^{-\dagger}=(nxex*nyex)*F
c
c     a) Perform forward 2d-FFT:
c-----------------------------------------------------------------------
          call cplx_fft_2d(nxex,nyex,nmin,nmax,reta,ieta,wsave,lensav)
          do m=1,nxex
            do n=1,nyex
              ik=ikstar(m,n)
              if(ik.le.nkstar)then ! elliptic truncation
c-----------------------------------------------------------------------
c     b) multiply with (nxex*nyex);
c-----------------------------------------------------------------------
                reta(m,n)=reta(m,n)*nxy
                ieta(m,n)=ieta(m,n)*nxy
c-----------------------------------------------------------------------
c     Multiply from the left with L^{-1}, where L is a diagonal matrix
c     with elements 1/sqrt(\gamma_ibox) along the diagonal, where
c     \gamma_ibox denotes the bi-Fourier spectral density of the variance
c     for level/component ibox (ibox is a level/component super-index):
c-----------------------------------------------------------------------
                reta(m,n)=reta(m,n)*sqrt_gamma(ibox,ik)
                ieta(m,n)=ieta(m,n)*sqrt_gamma(ibox,ik)
c-----------------------------------------------------------------------
c     Compute chi = \Lambda^{1/2} * X^{T} * \chi, where \Lambda is
c     the diagonal matrix containing the reduced set of eigenvalues of
c     the spectral covariance matrix, and X is the matrix containing the
c     corresponding eigenvectors:
c-----------------------------------------------------------------------
                do ir=1,nv1
                  chi(ir,m,n)=chi(ir,m,n)+vt(ibox,ir,ik)*
     &                     sqrt_lambda(ir,ik)*cmplx(reta(m,n),ieta(m,n))
                enddo
              endif
            enddo            ! spectral component n
          enddo               ! spectral component m
        enddo                  ! altitude l
      enddo                     ! component k

      end


