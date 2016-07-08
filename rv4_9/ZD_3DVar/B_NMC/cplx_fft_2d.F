C---- Comment Start --------------------------------
C
C     cplx_fft_2d
C
C     @shortdescr driver for 2dim fast Fourier transform
C
C     @keywords
C
C     @description
C     Transform the normalised model errors from physical to
C     spectral space by performing a 2D fast Fourier transform
C
C     on input, RE_ETA contains the normalised model errors
C        in physical space (real numbers!)
C
C     on output, RE_ETA and IM_ETA contain the real and imaginary
C        parts of the model error field in spectral space
C
C     @author M Kahnert
C
C
C
C---- Comment End ---------------------------------
      subroutine cplx_fft_2d(NXEX,NYEX,NMIN,NMAX,RE_ETA,IM_ETA,
     &     WSAVE,LENSAV)
      implicit none
      integer NXEX,NYEX,LENSAV
      integer nmin(nxex),nmax(nxex)
      real RE_ETA(NXEX,NYEX), IM_ETA(NXEX,NYEX), WSAVE(LENSAV)

      complex C(NXEX,NYEX)
      integer i,j,m,n,lenwork,iostat,nx
      real work(2*nxex*nyex)

      lenwork=2*nxex*nyex

c-----------------------------------------------------------------------
c     copy physical input field RE_ETA to C
c-----------------------------------------------------------------------
      do i=1,nxex
        do j=1,nyex
          c(i,j)=cmplx(re_eta(i,j),0e0)
        enddo
      enddo

c-----------------------------------------------------------------------
c     perform 2d-FFT
c-----------------------------------------------------------------------
      call cfft2f(nxex,nxex,nyex,c,wsave,lensav,work,lenwork,iostat)
      if(iostat.ne.0) print *,'ERROR in cfft2f, iostat=',iostat
c-----------------------------------------------------------------------
c     OBS: fftpack uses the following indexing convention for the FFT g
c     of a function f:
c       g(m) = (1/N) * \sum_{j=1}^N f(j) * \exp[-2\pi i * (j-1)(m-1)/N ]
c       f(j) = \sum_{m=1}^N g(m) * \exp[+2\pi i * (j-1)(m-1)/N ].
c
c     Consequently, if the input field f is real, the output field g
c     has the following symmetry: g(nxex-m+2,nyex-n+2)=conjg(g(m,n)).
c     Owing to limited numerical precision, this symmetry may not be
c     perfectly fulfilled in the output
c     => enforce symmetry:
c-----------------------------------------------------------------------
      if(mod(nyex,2).eq.0)then
        nx=(nyex+2)/2
      else
        nx=(nyex+1)/2
      endif
      do m=1,nxex
        i=nxex-m+2
        if(i.gt.nxex)i=i-nxex
        do n=1,nx
          j=nyex-n+2
          if(j.gt.nyex)j=j-nyex
!         if((n.le.nmin(m).or.n.ge.nmax(m)).and.
!    &       (i.le.nmin(j).or.i.ge.nmax(j)))c(i,j)=conjg(c(m,n))
          if((n.le.nmin(m).or.n.ge.nmax(m)).and.
     &       (j.le.nmin(i).or.j.ge.nmax(i)))c(i,j)=conjg(c(m,n))
        enddo
      enddo
      c(1,1)=cmplx(real(c(1,1)),0e0)
      if(mod(nxex,2).eq.0.and.mod(nyex,2).eq.0)then
        m=(nxex+2)/2
        n=(nyex+2)/2
        if((n.le.nmin(m).or.n.ge.nmax(m)))c(m,n)=cmplx(real(c(m,n)),0e0)
      endif
      if(mod(nxex,2).eq.0)then
        m=(nxex+2)/2
        c(m,1)=cmplx(real(c(m,1)),0e0)
      endif
      if(mod(nyex,2).eq.0)then
        n=(nyex+2)/2
        if((n.le.nmin(1).or.n.ge.nmax(1)))c(1,n)=cmplx(real(c(1,n)),0e0)
      endif

c-----------------------------------------------------------------------
c     save spectral output results in re_eta and im_eta:
c-----------------------------------------------------------------------
      do m=1,nxex
        do n=1,nyex
          re_eta(m,n)=real(c(m,n))
          im_eta(m,n)=aimag(c(m,n))
        enddo
      enddo

      return
      end
