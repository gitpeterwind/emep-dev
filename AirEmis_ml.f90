module AirEmis_ml
  ! Emissions from aircraft and lightning. Replaces eulcon_mach
  !
   use Par_ml,            only : MAXLIMAX, MAXLJMAX, NPROC, me
   use ModelConstants_ml, only : KCHEMTOP, KMAX_MID, KMAX_BND
!hf u2   use My_Runmode_ml,     only : stop_test
   implicit none
   private

   real, public, dimension(KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX), save :: &
                  airn                 & ! aircraft NOx emissions
                 ,airlig                 ! lightning NOx emissions

   public :: aircraft_nox
   public :: lightning

   private :: air_inter

   integer,private ,parameter :: ILEV=18

 contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine aircraft_nox(newseason)

	use Io_ml  ,           only : IO_AIRN, ios, open_file

!	input
      integer,intent(in):: newseason

!	local
	integer ILON, IGL
	parameter (ILON=128, IGL=32)
	integer i, j, k, nlon, ngl, nlev, level,nlevb
	integer intnox(ILON,64)
	REAL  ygrdum(IGL), ygrida(IGL), flux(ILON,IGL,-1:ILEV)		&
		,rlon(ILON+1), area(IGL)
	real DLON,RLON0	
	parameter (DLON = 4.21875-1.40625,RLON0 = -1.40625)
	real zrmin, zfak, secmonth

        character*20 fname

	data ygrdum / 87.86379884, 85.09652699, 82.31291295, 79.52560657 &
	,76.73689968, 73.94751515, 71.15775201, 68.36775611, 65.5776070	 &
	,62.78735180, 59.99702011, 57.20663153, 54.41619953, 51.6257336	 &
	,48.83524097, 46.04472663, 43.25419467, 40.46364818, 37.6730896	 &
	,34.88252099, 32.09194388, 29.30135962, 26.51076933, 23.7201739	 &
	,20.92957425, 18.13897099, 15.34836476, 12.55775612, 9.76714556	 &
	,6.97653355, 4.18592053, 1.39530/

	data area /3.4123374E+09,  8.2558925E+09,  1.2956985E+10 	&
		,1.7624316E+10,  2.2249359E+10,  2.6821497E+10 		&
		,3.1329966E+10,  3.5764105E+10,  4.0113402E+10 		&
		,4.4367553E+10,  4.8516461E+10,  5.2550296E+10 		&
		,5.6459481E+10,  6.0234756E+10,  6.3867154E+10		&
		,6.7348070E+10,  7.0669238E+10,  7.3822790E+10 		&
		,7.6801245E+10,  7.9597535E+10,  8.2205024E+10 		&
		,8.4617535E+10,  8.6829343E+10,  8.8835195E+10 		&
		,9.0630349E+10,  9.2210528E+10,  9.3572006E+10 		&
		,9.4711529E+10,  9.5626412E+10,  9.6314483E+10		&
		,9.6774103E+10,  9.7004184E+10/

! ---- Defines the ANCAT grid ----------------------------------------------
!   WARNING!!!  This is not the correct t42 grid, but the aircraft 
!               emissions are defined in this grid

	secmonth = 3600.*24.*31.
	flux = 0.

	if(me == 0)then

! --- LESER UTSLIPPS DATA FRA FILER HENTET FRA DLR

          write(fname,fmt='(''ancat'',i2.2,''.dat'')') newseason

          call open_file(IO_AIRN,"r",fname,needed=.true.,skip=1)
          !odin if (ios /= 0 ) call stop_all("ios error: ancat")
          if (ios /= 0) call gc_abort(me,NPROC,"ios error: ancat")
        end if ! me == 0


!hf u2	call stop_test(.true.,me,NPROC,ios,"ios error: ancat")

	if(me == 0)then
	  read(IO_AIRN,'(3i4,2e22.13)') nlon,ngl,nlev,zrmin,zfak
	  write(6,*) nlon,ngl,nlev,zrmin,zfak
	  do 10 k=0,NLEV-1
	    read(IO_AIRN,'(i2)') level
	    read(IO_AIRN,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
	    do 10 i=1,IGL
	      do 10 j=1,nlon
	        flux(j,i,k)=(float(intnox(j,i))*zfak)+zrmin
10	  continue

	  close(IO_AIRN)

          call open_file(IO_AIRN,"r","ancatmil.dat",needed=.true.,skip=1)
          !odin if (ios /= 0 ) call stop_all("ios error: ancatmil")
          if (ios /= 0) call gc_abort(me,NPROC,"ios error: ancatmil")
        end if ! me == 0

!hf u2	call stop_test(.true.,me,NPROC,ios,"ios error: ancatmil")

	if(me == 0)then
	  read(IO_AIRN,'(3i4,2e22.13)') nlon,ngl,nlevb,zrmin,zfak
	  write(6,*) nlon,ngl,nlevb,zrmin,zfak
	  do 11 k=1,nlevb
	    read(IO_AIRN,'(i2)') level
	    read(IO_AIRN,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
	    do 11 i=1,IGL
	      do 11 j=1,nlon
	        flux(j,i,k)=flux(j,i,k)+(float(intnox(j,i))*zfak)+zrmin
11	  continue

	  close(IO_AIRN)

!su	sumnox = 0.
!su	do 20 k=0,nlev-1
!su	  do 20 i=1,ngl
!su	    do 20 j=1,nlon
!su	      flux(j,i,k)=(float(intnox(j,i,k))*zfak)+zrmin
!su	      sumnox = sumnox + flux(j,i,k)
!su20	continue
!su	  write(6,*) 'SUMNOX, ANCAT:',sumnox

	endif

	call air_inter(ILON,IGL,0			&
		,flux,airn				&
		,ygrdum, ygrida,DLON,RLON0		&
		,rlon,area,secmonth)

    end subroutine aircraft_nox

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine lightning()

	use Io_ml,         only : IO_LIGHT, ios, open_file
        use ModelConstants_ml ,    only : current_date
	implicit none

	integer ILON, IGL
	parameter (ILON=64, IGL=16)

	integer i, j, k, nlon, ngl, nlev ,level

	integer intnox(ILON,32)
	real flux(ILON,IGL,-1:ILEV)
	real ygrdum(IGL), ygrida(IGL)		&
		,rlon(ILON+1), area(IGL)
	real  zrmin, zfak, secmonth, sumnox
	real DLON,RLON0	
	parameter (DLON = 8.4375 - 2.8125,RLON0 = -2.8125)

        character*20 fname

	data ygrdum / 85.76058712, 80.26877907, 74.74454037	&
		,69.21297617, 63.67863556, 58.14295405		&
		,52.60652603, 47.06964206, 41.53246125		&
		,35.99507841, 30.45755396, 24.91992863		&
		,19.38223135, 13.84448373,  8.30670286		&
		,2.76890301/

	data area / 26851631001.37,  64778953547.94, 101133318591.76	&
		,136519495343.02, 170627187541.37, 203140714921.94	&
		,233757180440.66, 262190945894.60, 288176619318.18	&
		,311471618334.64, 331858463070.53, 349146817322.73	&
		,363175270173.37, 373812845114.06, 380960223957.41	&
		,384550674664.23/

! ---- Defines the ANCAT grid ----------------------------------------------
!   WARNING!!!  This is not the correct t21 grid, 
!     Talk to Jostein about obtaining the original grid

	secmonth = 1.

        flux = 0.

! --- LESER UTSLIPPS DATA FRA FILER HENTET FRA DLR

	if(me == 0)then

	  sumnox = 0.

	  write(fname,fmt='(''lightn'',i2.2,''.dat'')') 	&
		current_date%month

          ! ds - open and read 1 line of header
          call open_file(IO_LIGHT,"r",fname,needed=.true.,skip=1)
          !odin if (ios /= 0 ) call stop_all("ios error: lightning")
          if (ios /= 0 ) call gc_abort(me,NPROC,"ios error: lightning")
         end if
!hf u2         call stop_test(.true.,me,NPROC,ios,"ios error: lightning")

	if(me == 0)then
	  read(IO_LIGHT,'(3i4,2e22.13)') nlon,ngl,nlev,zrmin,zfak
	  write(6,*) nlon,ngl,nlev,zrmin,zfak
	  do 10 k=1,nlev
	    read(IO_LIGHT,'(i2)') level
	    read(IO_LIGHT,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
	    do 10 i=1,IGL
	      do 10 j=1,nlon
		flux(j,i,k)=(float(intnox(j,i))*zfak)+zrmin
		sumnox = sumnox + flux(j,i,k)
10	  continue

	  close(IO_LIGHT)

	  write(6,*) 'lightningumnox',sumnox
	endif


	call air_inter(ILON,IGL,1			&
		,flux,airlig				&
		,ygrdum, ygrida,DLON,RLON0		&
		,rlon,area,secmonth)


    end subroutine lightning

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine air_inter(ilon,igl,iktop		&
		,flux,airem				&
		,ygrdum, ygrida,dlon,rlon0		&
		,rlon,area,secmonth)

	use Par_ml   , only: NPROC,me,limax,ljmax
	use GridValues_ml         , only : gl,gb
        use PhysicalConstants_ml , only : AVOG

	implicit none

	integer, parameter :: KMAX_BND_AIR = 21 !pw u3
	integer, intent(in) :: ilon,igl,iktop
	real, intent(out) ::   ygrida(igl)
	real, intent(inout) ::   flux(ilon,igl,-1:ILEV)
      real, dimension(KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX), intent(out) :: airem
	real, intent(out) :: rlon(ilon+1)
	real, intent(in) :: area(igl), ygrdum(igl),dlon,rlon0,secmonth

!	local
	integer info
	real height,atwno2
	real sum,sumnox,volcm,sum2
	integer ixn(MAXLIMAX,MAXLJMAX)			&
		,jxn(MAXLIMAX,MAXLJMAX)			&
            ,ilevel(KMAX_MID)
	integer lon,lat,i,j,ig,jg,kg,k
      real frac, fraca(KMAX_MID), fracb(KMAX_MID), above,below, vol(KMAX_MID)
	real glij

!pw u3      real, dimension(KMAX_BND):: zz1
      real, dimension(KMAX_BND_AIR):: zz1


   data zz1 /16175.41, 14138.79, 12562.57, 11274.08, 10058.51		&
		,8710.764, 7294.476, 5944.395, 4806.321, 3880.146	&
		,3106.688, 2456.548, 1909.316, 1450.109, 1068.357  	&
		,755.2472, 505.0488, 313.8083, 178.2393, 88.75819	&
		,0.0000000E+00/

   if(KMAX_BND_AIR.ne.KMAX_BND)then
      write(*,*)'AirEmis : KMAX_BND hardcoded!',KMAX_BND,KMAX_BND_AIR
      call gc_abort(me,NPROC,'please, modify subroutine air_inter!')
   endif

	height = 1.e5
	atwno2 = 46.
	if(me.eq.0)then

	  sum = 0.
	  sumnox = 0.
	  do k = iktop,ILEV
	    do lat = 7,13
	      do lon = 1,5
	        sum = sum + flux(lon,lat,k)
	      end do
	    end do

	    do lat=1,igl
	      volcm = area(lat)*1.e4*height
	      do lon=1,ilon
	        sumnox = sumnox + flux(lon,lat,k)
	        flux(lon,lat,k)=flux(lon,lat,k)*1.e3*AVOG	&
			/volcm/secmonth/atwno2
	      end do
	    end do
	  end do
	  write(6,*) 'SUMNOX, ANCAT:',sumnox

	  call gc_rbcast(317, igl*ilon*(ILEV+1-iktop)	&
		,0, NPROC,info,flux(1,1,iktop))

	endif		!me=0
! -- N/S
	ygrida(1) = 90.
	do i=2,igl
	  ygrida(i) = (ygrdum(i-1)+ygrdum(i))*0.5
	enddo

! -  E/V
	rlon(1) = rlon0
	do i=2,ilon
	  rlon(i) = rlon(i-1) + dlon
	end do
	rlon(ilon+1) = rlon(1)+360.

!    Assign Ines t42 gridpoints to the EMEP grid

	jg = igl-1

	do j = 1,ljmax
	  do i = 1,limax
	    if(gb(i,j).eq.90. .and. gl(i,j).eq.0.0) then
	      ixn(i,j) = 1
	      jxn(i,j) = 1
	    else
	      do while(gb(i,j).lt.ygrida(jg+1))
		jg = jg+1
	      enddo
	      do while(gb(i,j).ge.ygrida(jg))
		jg = jg-1
	      enddo
	      jxn(i,j) = jg
	      glij = gl(i,j)
	      if(glij.le.rlon(1))glij = glij+360.
	      ig = int((glij-rlon(1))/dlon)+1
	      if(glij.eq.rlon(ig))ig=ig+1
	      ixn(i,j) = ig
	    end if
	  end do
	end do

	if(me.ne.0)then
	  call gc_rbcast(317, igl*ilon*(ILEV+1-iktop)		&
			,0, NPROC,info,flux(1,1,iktop))
	endif

	kg = 1
	above = 1.e3
	below = 0.
      do k = KMAX_MID,KCHEMTOP,-1
	  do while(zz1(k+1).gt.below+1.e3) 
	    below = below+1.e3
	  end do
	  do while (zz1(k).gt.above) 
	    kg = kg+1
	    above = above+1.e3
	  end do
	  ilevel(k) = kg
	  fraca(k) = 1.
	  if(above-below.gt.1.1e3) 				&
		fraca(k) = (zz1(k)-(above-1.e3))/(zz1(k) - zz1(k+1))
	  fracb(k) = 0.
	  if(above-below.gt.2.1e3)				&
		fracb(k) = (below+1.e3 - zz1(k+1))/(zz1(k) - zz1(k+1))
	end do

	do j = 1,ljmax
	  do i = 1,limax
	    lon = ixn(i,j)
	    lat = jxn(i,j)
          do k = KCHEMTOP,KMAX_MID
	      frac = 1. - fraca(k) - fracb(k)
	      kg = ilevel(k)
	      airem(k,i,j) = flux(lon,lat,kg)*fraca(k) 		&
			+ flux(lon,lat,kg-1)*frac		&
			+ flux(lon,lat,kg-2)*fracb(k)
	    end do

!  surface emissions

	    if(iktop.eq.0)		&
            airem(KMAX_MID,i,j) = airem(KMAX_MID,i,j)+flux(lon,lat,0)
	  end do
	end do

	sum2 = 0.
      do k=KCHEMTOP,KMAX_MID
	  vol(k) = 5.e6*5.e6*(zz1(k)-zz1(k+1))*1.e2
	end do
	do j = 1,ljmax
	  do i = 1,limax
	    if(gl(i,j).gt.rlon(1) .and. gl(i,j).lt.rlon(6) .and.	&
		gb(i,j).lt.ygrida(7) .and. gb(i,j).gt.ygrida(14)) then
            do k=KCHEMTOP,KMAX_MID
		sum2 = sum2 + airem(k,i,j)*atwno2/AVOG*		&
			vol(k)*secmonth*1.e-3
	      end do
	    end if
	  end do
	end do

	call gc_rsum(1,NPROC,info,sum2)
	if(me.eq.0) write(6,*) 'ancat on limited area:',sum,sum2

    end subroutine air_inter
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


end module AirEmis_ml
