module AirEmis_ml
  ! Emissions from aircraft and lightning.
  ! Ancat emissions converted to flux as molecules cm-3 s-1 in ancat grid. 
  ! Emissions on (finer) model grid then assigned from the acat grid that 
  ! model grid falls within. Note that aircraft emissions is given on a t42 
  ! and lightning on t21.
  ! 
   use Par_ml,            only : MAXLIMAX, MAXLJMAX, limax,ljmax, NPROC, me
   use ModelConstants_ml, only : KCHEMTOP, KMAX_MID, KMAX_BND, current_date
   use Io_ml  ,           only : IO_AIRN, IO_LIGHT, ios, open_file
   use GridValues_ml        , only : gl,gb, GRIDWIDTH_M
   use PhysicalConstants_ml , only : AVOG
   use Met_ml,                only : z_bnd  

   implicit none
   private

   real, public, dimension(KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX), save :: &
                  airn                 & ! aircraft NOx emissions
                 ,airlig                 ! lightning NOx emissions

   public :: aircraft_nox !reads in the raw data
   public :: lightning

   private :: air_inter !interpolate the data into required grid

   integer,private ,parameter :: ILEV=18

 contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      subroutine aircraft_nox(newseason)


!	input
      integer,intent(in):: newseason

!	local
   integer, parameter ::  ILON = 128, IGL = 32, GGL = 64
   real,    parameter ::  DLON = 4.21875-1.40625,RLON0 = -1.40625

	integer i, j, k, nlon, ngl, nlev, level,nlevb
	real zrmin, zfak, secmonth

! Definition of the ancat grid   ~ t42 :
! Data read in from N --> S  and from longitude 0  ( not from +-180 )  
   real, dimension(ggl) ::     ygrida       ! grid mid. pt. N-S
   real, dimension(igl) ::     ygrdum   &   ! grid mid. pt. N-S
                               ,area        ! grid area N-S
   real, dimension(ilon+1) ::   rlon        !  

   integer, dimension(ILON,GGL) ::          intnox ! global emission flux  
   real,    dimension(ILON,GGL,-1:ILEV) :: flux   ! emission flux converted to

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
	flux(:,:,:) = 0.

	if(me == 0)then

! --- Open and read ancat data (originally from DLR, EU project POLINAT)
! --- Commercial aircraft emissions every season
! --- Military aircraft emission read in as annual data


          write(fname,fmt='(''ancat'',i2.2,''.dat'')') newseason

          call open_file(IO_AIRN,"r",fname,needed=.true.,skip=1)
          if (ios /= 0) call gc_abort(me,NPROC,"ios error: ancat")
        end if ! me == 0



	if(me == 0)then
	  read(IO_AIRN,'(3i4,2e22.13)') nlon,ngl,nlev,zrmin,zfak
	  write(6,*) nlon,ngl,nlev,zrmin,zfak
	  do k = 0,NLEV-1
	    read(IO_AIRN,'(i2)') level
	    read(IO_AIRN,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
	    do  i = 1,NGL
	       do j = 1,nlon
       	          flux(j,i,k)=(float(intnox(j,i))*zfak)+zrmin
               end do
            end do
          end do

	  close(IO_AIRN)

          call open_file(IO_AIRN,"r","ancatmil.dat",needed=.true.,skip=1)
          if (ios /= 0) call gc_abort(me,NPROC,"ios error: ancatmil")
        end if ! me == 0


	if(me == 0)then
	  read(IO_AIRN,'(3i4,2e22.13)') nlon,ngl,nlevb,zrmin,zfak
	  write(6,*) nlon,ngl,nlevb,zrmin,zfak
	  do  k = 1,nlevb
	    read(IO_AIRN,'(i2)') level
	    read(IO_AIRN,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
	    do  i = 1,NGL
	      do  j = 1,nlon
	        flux(j,i,k)=flux(j,i,k)+(float(intnox(j,i))*zfak)+zrmin
              end do
            end do
          end do

	  close(IO_AIRN)


	endif

	call air_inter(ILON,IGL,GGL,0			&
		,flux,airn				&
		,ygrdum, ygrida,DLON,RLON0		&
		,rlon,area,secmonth)

    end subroutine aircraft_nox

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine lightning()

	implicit none

   integer, parameter ::  ILON = 64, IGL = 16, GGL = 32
   real,    parameter ::  DLON =  8.4375 - 2.8125, RLON0 = -2.8125


	integer i, j, k, nlon, ngl, nlev ,level

	real  zrmin, zfak, secmonth, sumnox



! Definition of the ancat grid   ~ t21 :
! Data read in from N --> S  and from longitude 0  ( not from +-180 )  
! NB!!  note the difference between lightning and aircraft emission grid


   real, dimension(ggl) ::     ygrida       ! grid mid. pt. N-S
   real, dimension(igl) ::     ygrdum   &   ! grid mid. pt. N-S
                               ,area        ! grid area N-S
   real, dimension(ilon+1) ::   rlon        !  

   integer, dimension(ILON,GGL) ::          intnox ! global emission flux  
   real,    dimension(ILON,GGL,-1:ILEV) :: flux   ! emission flux converted to



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

	secmonth = 1.

        flux(:,:,:) = 0.

! --- LESER UTSLIPPS DATA FRA FILER HENTET FRA DLR

	if(me == 0)then

	  sumnox = 0.

	  write(fname,fmt='(''lightn'',i2.2,''.dat'')') 	&
		current_date%month

          ! ds - open and read 1 line of header
          call open_file(IO_LIGHT,"r",fname,needed=.true.,skip=1)
          if (ios /= 0 ) call gc_abort(me,NPROC,"ios error: lightning")
         end if

	if(me == 0)then
	  read(IO_LIGHT,'(3i4,2e22.13)') nlon,ngl,nlev,zrmin,zfak
	  write(6,*) nlon,ngl,nlev,zrmin,zfak
	  do  k = 1,nlev
	    read(IO_LIGHT,'(i2)') level
	    read(IO_LIGHT,'(12i6)') ((intnox(j,i),j=1,nlon),i=1,ngl)
	    do  i = 1,NGL
	      do  j = 1,nlon
		flux(j,i,k)=(float(intnox(j,i))*zfak)+zrmin
      		sumnox = sumnox + flux(j,i,k)
              end do
            end do
          end do


	  close(IO_LIGHT)

	  write(6,*) 'lightningumnox',sumnox
	endif


	call air_inter(ILON,IGL,GGL,1			&
		,flux,airlig				&
		,ygrdum, ygrida,DLON,RLON0		&
		,rlon,area,secmonth)


    end subroutine lightning

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine air_inter(ilon,igl,ggl,iktop		&
         ,flux,airem				&
         ,ygrdum, ygrida,dlon,rlon0		&
         ,rlon,area,secmonth)


      implicit none

      integer, parameter :: KMAX_BND_AIR = 21 !pw u3

      integer, intent(in) :: ilon,igl,ggl,iktop
      real, intent(in) :: area(igl), ygrdum(igl),dlon,rlon0,secmonth


      real, intent(inout) ::   flux(ilon,ggl,-1:ILEV)

      real, dimension(KCHEMTOP:KMAX_MID,MAXLIMAX,MAXLJMAX), intent(out) :: airem
      real, intent(out) ::   ygrida(ggl)
      real, intent(out) :: rlon(ilon+1)

      !	local
      integer info
      integer lon,lat,i,j,ig,jg,kg,k, i_sh
      integer la_tst1, la_tst2, lo_tst1, lo_tst2   !  test area for sums
      real height,     &  !  height of the emission levels
           atwno2,     &  !  atomic weight of NO2
           vol            !  volume of model grid boxes
      real frac, above, below, glij
      real sum,sumnox,volcm,sum2
      integer, dimension(MAXLIMAX,MAXLJMAX) :: ixn  & !  mapping of emission 
           ,jxn    !  grid to model grid
      integer, dimension(KMAX_MID)          :: ilevel
      real  fraca(KMAX_MID), fracb(KMAX_MID)


      height = 1.e5
      atwno2 = 46.

      !  print out values on a sub-domain for comparison with dirrect model innput
      !  NB!  due to different resolution the subdomain will be different for 
      !       aircraft and lightning emissions

      la_tst1 = 7
      la_tst2 = 13
      lo_tst1 = 1
      lo_tst2 = 5

      if(me.eq.0)then
         sum = 0.
         sumnox = 0.
         do k = iktop,ILEV
	    do lat = la_tst1,la_tst2
               do lon = lo_tst1,lo_tst2
                  sum = sum + flux(lon,lat,k)!test
               end do
	    end do

	    do lat=1,ggl
               if(lat<=igl)then
                  volcm = area(lat)*1.e4*height
               else
                  !area not defined for Southern Hemisphere
                  volcm = area(ggl-lat+1)*1.e4*height
               endif
               do lon=1,ilon
                  sumnox = sumnox + flux(lon,lat,k)
                  flux(lon,lat,k)=flux(lon,lat,k)*1.e3*AVOG	&
                       /volcm/secmonth/atwno2
               end do
	    end do
         end do
         write(6,*) 'SUMNOX, ANCAT:',sumnox

      endif		!me=0


      call gc_rbcast(317, ggl*ilon*(ILEV+1-iktop)	&
           ,0, NPROC,info,flux(1,1,iktop))

      ! -- N/S
      ygrida(1) = 90.
      ygrida(GGL) = -90.
      do i=2,igl
         i_sh = GGL + 1 - i
         ygrida(i) = (ygrdum(i-1)+ygrdum(i))*0.5
         ygrida(i_sh) = - ygrida(i)
      enddo

      ! -  E/V
      rlon(1) = rlon0
      do i=2,ilon
         rlon(i) = rlon(i-1) + dlon
      end do
      rlon(ilon+1) = rlon(1)+360.

      !    Assign gridpoints to the EMEP grid

      jg = ggl-1

      do j = 1,ljmax
         do i = 1,limax
!pw	    if(gb(i,j).eq.90. .and. gl(i,j).eq.0.0) then
	    if(abs(gb(i,j)-90.0)<0.001) then
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
!pw            if(glij.eq.rlon(ig))ig=ig+1
               if(ig>ilon)ig=ig-ilon !pw
               ixn(i,j) = ig

	    end if
         end do
      end do

      !	if(me.eq.0)then
      !           write(*,*) 'ygrida  ', ggl
      !           do i = 1,ggl
      !           write(*,*) 'ygrida  ', i,ygrida(i)
      !           end do
      !	endif

      do j = 1,ljmax
         do i = 1,limax

 	    kg = 1
	    above = 1.e3
	    below = 0.

            do k = KMAX_MID,KCHEMTOP,-1
               do while(z_bnd(i,j,k+1).gt.below+1.e3) 
                  below = below+1.e3
               end do
               do while (z_bnd(i,j,k).gt.above) 
                  kg = kg+1
                  above = above+1.e3
               end do
               ilevel(k) = kg
               fraca(k) = 1.
               if(above-below.gt.1.1e3) 				&
                    fraca(k) = (z_bnd(i,j,k)-(above-1.e3))  &
                    /(z_bnd(i,j,k) - z_bnd(i,j,k+1))
               fracb(k) = 0.
               if(above-below.gt.2.1e3)				&
                    fracb(k) = (below+1.e3 - z_bnd(i,j,k+1))           &
                    /(z_bnd(i,j,k) - z_bnd(i,j,k+1))
            end do

            lon = ixn(i,j)
            lat = jxn(i,j)
if(lon>128)write(*,*)'longg',i,j,ixn(i,j)
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


      !! Print out on a limited part of the dommain both raw data ( flux ) and 
      !! the re-gridded emissions ( airem ).  Expect only an opproximate match! 
      !! Summation of aircraft emission and lightning are on different grids 
      !! and hence for different domains. 
      sum2 = 0.
      do j = 1,ljmax
         do i = 1,limax
            if(gl(i,j).gt.rlon(lo_tst1) .and. gl(i,j).lt.rlon(lo_tst2+1) .and. &
                 gb(i,j).lt.ygrida(la_tst1) .and. gb(i,j).gt.ygrida(la_tst2+1)) then
               do k=KCHEMTOP,KMAX_MID
                  vol = GRIDWIDTH_M*GRIDWIDTH_M   &
                       *(z_bnd(i,j,k)-z_bnd(i,j,k+1))*1.e6
                  sum2 = sum2 + airem(k,i,j)*atwno2/AVOG*		&
                       vol*secmonth*1.e-3
               end do
            end if
         end do
      end do

      call gc_rsum(1,NPROC,info,sum2)
      if(me.eq.0) write(6,*) 'ancat on limited area:',sum,sum2

    end subroutine air_inter
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


end module AirEmis_ml
