C---- Comment Start --------------------------------
C
C     chitox
C
C     @shortdescr Transform control variable chi to x
C
C     @keywords
C
C     @description
C     Inverse preconditioning of the control variable
C
C     @author M Kahnert
C
C
C
C---- Comment End ---------------------------------
      subroutine chitox(chi,x)
      use spectralcov
      implicit none

      real x(nxex,nyex,nlev,nchemobs)
      complex chi(nv1,nxex,nyex)

      integer i,j,k,kk,l,ik,m,n,lenwork,iostat,ibox
      real work(2*nxex*nyex)

      integer :: ndim
      complex, pointer :: temp1(:), temp2(:), temp3(:,:,:), c(:,:)

      lenwork=2*nxex*nyex
      ndim=nchemobs*nlev
      allocate(temp1(nv1),temp2(ndim),temp3(ndim,nxex,nyex))
      do n=1,nyex
        do m=1,nxex
          ik=ikstar(m,n)
          if(ik<=nkstar)then
            temp1(:)=chi(:,m,n)*sqrt_lambda(:,ik)
            temp2(:)=matmul(vt(:,:,ik),temp1)
            temp3(:,m,n)=temp2(:)*sqrt_gamma(:,ik)
          else
            temp3(:,m,n)=0.0
          endif
        enddo
      enddo
      deallocate(temp1,temp2)

c-----------------------------------------------------------------------
c     Compute c = X * \Lambda^{1/2} * \chi, where \Lambda is the
c     diagonal matrix containing the reduced set of eigenvalues of the
c     spectral covariance matrix, and X is the matrix containing the
c     corresponding eigenvectors:
c-----------------------------------------------------------------------
      x=0e0
      do k=1,nchemobs
        kk=ichemobs(k)
        do l=1,nlev
          ibox=(k-1)*nlev+l
          c=>temp3(ibox,:,:)
c-----------------------------------------------------------------------
c     Perform inverse 2d-FFT:
c-----------------------------------------------------------------------
          call cfft2b(nxex,nxex,nyex,c,wsave,lensav,work,lenwork,iostat)
c     nested double loop over horizontal physical space coordinates:
c-----------------------------------------------------------------------
c     Only update whithin the non-extended zone (since there are no
c     observations in that zone => zero increment):
c-----------------------------------------------------------------------
c     Multiply from the left with S, where S is a diagonal matrix
c     containing the standard deviations sqrt(sigma).
c-----------------------------------------------------------------------
          forall(i=1:nx,j=1:ny) x(i,j,l,k)=real(c(i,j))*stddev(i,j,l,kk)
        enddo ! level l
      enddo    ! component k

      deallocate(temp3)

      end


