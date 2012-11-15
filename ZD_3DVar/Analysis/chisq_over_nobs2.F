C---- Comment Start --------------------------------
C
C     chisq_over_nobs2
C
C     @shortdescr Compute CHI^2/nobs for final analysis
C
C     @description Compute CHI^2/nobs on an observation by observation
C     basis by taken a subset of the total number of observations. This
C     should be of order 0.5.
C
C     <p>The code is not all too readable as observations are spread
C     over the processors, and some variables are hidden in include
C     files.
C
C     <p>Should CHI^2 be determined for each observation type separately?
C
C     @author M Kahnert (SMHI)
C     @author L Robertson (SMHI)
C
C---- Comment End ---------------------------------
      subroutine chisq_over_nobs2(nobs,dep,h)
C     use spectralcov
      implicit none

      integer nobs                    ! number of observations
      real dep(nobs)                  ! Analysis departure from background
      real h(nxex,nyex,nlev,nchemobs) ! Works array for observation operator

      complex uh(nv1,nxex,nyex)

      real :: chisq,chisq_,hbh
      integer :: i,j,k,n,nm,no4,no8,gi,gj,gn
      integer :: p,kk,ilev,frac
      character(len=128) file

      integer inside

      no4=4*nobs2
      no8=8*nobs2
      frac=nobs/max(nobs/5,min(200,nobs))

      ON_MONITOR  write(*,*)'CHI^2 calculations'
c-----------------------------------------------------------------------
c     Determine (H^T)_n; compute U^{-\dagger}*(H^T)_n, n: nth component
c     of H^T, n=1,...,nobs
c-----------------------------------------------------------------------
      fchisq(n)=-1
      chisq=0e0
      n=min(nobs,1)
      nm=0
      do gn=1,nobs
        h=0e0
c     We only try for every frac'th observation
        if(mod(gn,frac).eq.0)then
          nm=nm+1
          if(n.gt.0.and.obs_order(n).eq.gn)then
c     print*,"using",my_proc,gn,n
            do p=1,no4
              gi=pidx(p)
              gj=pidx(p+no4)
              i = gi
              j = gj + g2l_y(gj)
              ilev=pidx(p+no8)
              do k=1,nchemobs
                kk=ichemobs(k)
C Set up observation operator for this observations
                if(inside(float(i),float(j),nxex,nyex))
     +            h(i,j,ilev,k)=H_jac(n,p,kk)
              enddo
            enddo
          endif
C Transform to spectral space
          call chitox_adj(uh,h)
          chisq_=0
          hbh=0
C Determine HBH
          do j=1,nyex
            do i=1,nxex
              do ilev=1,nv1
                hbh=hbh+conjg(uh(ilev,i,j))*uh(ilev,i,j)
              enddo
            enddo
          enddo
C Determine HGH+R
          if(n.gt.0 .and. obs_order(n).eq.gn)then
            if(obs_owner(n).eq.my_proc)then
              hbh=hbh+obsstddev(n)**2
            endif
          endif
C Determine chi-square for this observation
          if(n.gt.0 .and. obs_order(n).eq.gn)then
            if(obs_owner(n).eq.my_proc)then
              chisq_=dep(n)*dep(n)/hbh
              fchisq(n)=chisq_
            endif
            n=min(nobs,n+1)
          endif
C Sum up chi-square
          chisq=chisq+chisq_
        elseif(n.gt.0.and.obs_order(n).eq.gn)then
          n=min(nobs,n+1)
        endif
      enddo

      BEGIN_MONITOR
C Normalize by number of used observations
      chisq=chisq/nm
      write(*,*)
      write(*,*)'CHI^2/nobs=',chisq
      write(*,*)"nobs=",nobs,nm
      write(*,*)
      write(file,'("chisqr_",i4,3i2.2,".dat")')year,month,day,hour
      open(10,file=file,form="formatted")
      write(10,'(2i3.2,i5.4,3i3.2,g12.4,i7)')month,day,year,hour,0,0,
     +     chisq,nm
      close(10)
      END_MONITOR

      end

