C---- Comment Start --------------------------------
C
C     spectral_filter
C
C     @shortdescr 2D spectral filter
C
C     @keywords
C
C     @description
C     Perform spectral filtering of a 2d field
C
C
C     @author M Kahnert
C
C
C
C---- Comment End ---------------------------------
      subroutine spectral_filter(ilev,k)!,ktrunc)
      use spectralcov
      implicit none

c     input:
      integer ilev,k
!     real ktrunc

c     local variables:
      complex c(nxex,nyex)
      integer m,n,ik,iostat,lenwork,n0,i,j!,m1,n1
      real work(2*nxex*nyex),avg

      if(ilev.eq.1.and.k.eq.1) print *,'ktrunc:',ktrunc

      lenwork=2*nxex*nyex

      if(ktrunc.eq.0e0)then
c     compute horizontally isotropic standard deviations:
        if(ilev.eq.1.and.k.eq.1)then
          write(*,*)
          write(*,*)
     &      'computing horizontally homogeneous, isotropic stddev'
          write(*,*)
        endif
        avg=0e0
        do j=1,ny ! nyex
          do i=1,nx ! nxex
            avg=avg+stddev(i,j,ilev,k)**2
          enddo
        enddo
        avg=sqrt(avg/float(nx*ny)) ! sqrt(avg/float(nxex*nyex))
        do j=1,nyex
          do i=1,nxex
            stddev(i,j,ilev,k)=avg
          enddo
        enddo
      elseif(ktrunc.gt.kstar(nkstar))then
c     do not perform spectral filtering:
        if(ilev.eq.1.and.k.eq.1)then
          write(*,*)
          write(*,*)'no spectral filtering of stddev'
          write(*,*)
        endif
        return

      else
c     perform spectral filtering:
        if(ilev.eq.1.and.k.eq.1)then
          write(*,*)
          write(*,*)'spectral filtering of stddev with k_trunc=',ktrunc
          write(*,*)
        endif

c     stddev array:
        do i=1,nxex
          do j=1,nyex
            c(i,j)=cmplx(stddev(i,j,ilev,k),0e0)
          enddo
        enddo

c     2d FFT:
        call cfft2f(nxex,nxex,nyex,c,wsave,lensav,work,lenwork,iostat)
        if(iostat.ne.0) print *,'ERROR in cfft2f, iostat=',iostat
c     truncation in spectral space:
        do m=1,nxex
          do n=1,nyex
            ik=ikstar(m,n)
            if(kstar(ik).gt.ktrunc)then
              c(m,n)=cmplx(0e0,0e0)
            endif
          enddo
        enddo
c
c         do m=-Kxmin,Kx
c            m1=mm(m)
c            do n=-Kymin,Ky
c               n1=nn(n)
c               ik=ikstar(m1,n1)
c               if(kstar(ik).gt.ktrunc)then
c                  c(m1,n1)=cmplx(0e0,0e0)
c               endif
c            enddo
c         enddo

c-----------------------------------------------------------------------
c     OBS: fftpack uses the following indexing convention for the FFT g
c     of a function f:
c     g(m) = (1/N) * \sum_{j=1}^N f(j) * \exp[-2\pi i * (j-1)(m-1)/N ]
c     f(j) = \sum_{m=1}^N g(m) * \exp[+2\pi i * (j-1)(m-1)/N ].
c
c     Consequently, if the input field f is real, the output field g
c     has the following symmetry: g(nxex-m+2,nyex-n+2)=conjg(g(m,n)).
c     Owing to limited numerical precision, this symmetry may not be
c     perfectly fulfilled in the output
c     => enforce symmetry:
c-----------------------------------------------------------------------
        if(mod(nyex,2).eq.0)then
          n0=(nyex+2)/2
        else
          n0=(nyex+1)/2
        endif
        do m=1,nxex
          i=nxex-m+2
          if(i.gt.nxex)i=i-nxex
          do n=1,n0
            j=nyex-n+2
            if(j.gt.nyex)j=j-nyex
            if((n.le.nmin(m).or.n.ge.nmax(m)).and.
     &         (i.le.nmin(j).or.i.ge.nmax(j))) c(i,j)=conjg(c(m,n))
          enddo
        enddo
        c(1,1)=cmplx(real(c(1,1)),0e0)
        if(mod(nxex,2).eq.0.and.mod(nyex,2).eq.0)then
          m=(nxex+2)/2
          n=(nyex+2)/2
          if((n.le.nmin(m).or.n.ge.nmax(m)))
     &      c(m,n)=cmplx(real(c(m,n)),0e0)
        endif
        if(mod(nxex,2).eq.0)then
          m=(nxex+2)/2
          c(m,1)=cmplx(real(c(m,1)),0e0)
        endif
        if(mod(nyex,2).eq.0)then
          n=(nyex+2)/2
          if((n.le.nmin(1).or.n.ge.nmax(1)))
     &      c(1,n)=cmplx(real(c(1,n)),0e0)
        endif


c     inverse 2d FFT:
        call cfft2b(nxex,nxex,nyex,c,wsave,lensav,work,
     &        lenwork,iostat)
c     copy spectrally filtered field to stddev:
        do i=1,nxex
          do j=1,nyex
            stddev(i,j,ilev,k)=real(c(i,j))
          enddo
        enddo
      endif

      end
