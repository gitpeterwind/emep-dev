!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module Met_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! ds rv1.2 MetModel_LandUse added here for snow and iclass
!  - combined from hf and pw Met_ml
! October 2001 hf added call to ReadField
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
use GridValues_ml,     only : xm,xmd, sigma_bnd,sigma_mid
use ModelConstants_ml, only : PASCAL, PT, CLOUDTHRES, METSTEP, &
                              KMAX_BND,KMAX_MID,NMET
use Par_ml           , only : MAXLIMAX,MAXLJMAX,NPROC,me  &
                               ,limax,ljmax,li0,li1,lj0,lj1  &
                               ,neighbor,WEST,EAST,SOUTH,NORTH,NOPROC  &
                               ,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2
use PhysicalConstants_ml, only : KARMAN, XKAP, R, CP, GRAV, ROWATER
use Tabulations_ml , only : TPI,PBAS,PINC
implicit none
private

    !----------------- basic met fields ----------------------------------!
    !  Here we declare the metoeorlogical fields used in the model        !
    !  From old eulmc.inc=eulmet.inc
    !---------------------------------------------------------------------!
!
! Vertical levels: z_mid,  z_bnd, sigma_mid, sigma_bnd
!=============================================================================
!  In the original MADE/MACHO setup the notation k1, k2, scor1, scor2 was
!*  used for heights and sigma coordinates.
!*  In this sytstem the "1" index is for the half layers, the "2" index is for 
!*  the full layers. Half layers (as they used to be defined in NORLAM) would
!*  then also represent the fixed upper and lower boundaries. The same index
!*  is used for height z. 
!*  According to this   - - - - -    is the full levels and   -------- the 1/2
!*  levels in the sketch below.
!* 
!*   In order to avoid confusion, this notation has now been changed, so 
!*   that "mid" and "bnd" are used as suffixes on z and sigma. 
!* 
!* 
!* 
!* ---------------------------
!* 
!* 
!* - - - - - - - - - - - -     KMAX_MID -1
!* 
!* 
!* --------------------------  KMAX_BDN-1       (z_bnd)   (sigma_bnd)
!* 
!* 
!* - - - - - - - - -           KMAX_MID(old kmax2) = 20    (z_mid)   (sigma_mid)   (old z2)
!* 
!* ------------------------    KMAX_BND = 21    (z_bnd)                 (old z1)
!* \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  surface \\\\\\\\\\\\\\\\
!* 


!
  real,public, save, dimension(0:MAXLIMAX,MAXLJMAX,KMAX_MID,NMET) :: u  ! m/s
  real,public, save, dimension(MAXLIMAX,0:MAXLJMAX,KMAX_MID,NMET) :: v  ! m/s
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET) :: sdot ! dp/dt
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID,NMET) :: &
                           th      &  ! Potential teperature  ( deg. k )
                           ,q      &  ! Specific humidity
                          ,cw         ! Cloud water           (kg water)/(kg air)


!    since pr,cc3d,cc3dmax used only for 1 time layer - define without NMET
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: &
                     pr      & ! Precipitation
                    ,cc3d    & ! 3-d cloud cover (cc3d),
                    ,cc3dmax & ! and maximum for layers above a given layer
                    ,trw    & ! total rainwater kg/kg 
                               ! NB: trw: accumulated during metstep. ( Direct output 
                               ! from MM5 gives accumulated during prognosis length )
                    ,lwc       !liquid water content


! surface fields
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,NMET) :: &
                     ps      & ! Surface pressure hPa (or Pa- CHECK!)
                    ,th2m    & ! Temp 2 m   deg. K
                    ,fh      & ! surf.flux.sens.heat W/m^2   ! ds u7.4vg added
                    ,fl      & ! latent heat flux W/m^2
                    ,fm      & ! surf. stress  N/m^2
                    ,ustar     ! pw u3 friction velocity m/s ustar^2 = fm/roa



! Derived met
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET) :: skh
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID,NMET) :: roa ! kg/m^3
  real,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &
                  psurf & !u7.4lu psa  Surface pressure hPa
                 ,surface_precip    & ! Surface precip mm/hr   ! ds rv1.6.2
                 ,t2      !u7.4vg temp2m  Temp 2 m   deg. K

  real,public, save, &
      dimension(MAXLIMAX,MAXLJMAX,KMAX_BND) :: z_bnd ! height of full layers
  real,public, save, &
      dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: z_mid ! height of half layers

!ds rv1.2  keep HIRLAM/xx met model landuse stuff in same routine
  integer,public, save, dimension(MAXLIMAX,MAXLJMAX) :: & 
       snow,    &  ! monthly snow (1=true), read in MetModel_LandUse
       iclass      ! roughness class ,         "       "       "


  logical, public, save :: foundclouds,foundustar,foundpreta,mm5 !pw u3


!hf tiphys
!check dimension
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: &
                xksig ! estimated exchange coefficient, Kz, in intermediate 
                      ! sigma levels, m2/s

  real,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &
       pzpbl,   &  !stores H(ABL) for averaging and plotting purposes, m
       Kz_min      ! Min Kz below hmix  !hf Hilde&Anton

!ds  real, dimension(MAXLIMAX,MAXLJMAX) :: ven !ventilation coefficient, m3
!hf end tiphys

 public :: infield
 public :: MetModel_LandUse    ! rv1.2 combines old in_isnowc and inpar
 !ds rv1.2 public :: in_isnowc
 public :: metvar
 public :: metint
 public :: tiphys           !hf NEW

 contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      subroutine infield(numt)

!	the subroutine reads meteorological fields and parameters every
!	six-hour from fields-files, divide the fields into domains
!	and sends subfields to the processors using nx calls

	use ModelConstants_ml , only : identi	&
                                       ,current_date & ! date-type
                                       , METSTEP
	use Par_ml ,  only :  limax,ljmax
        use Dates_ml, only : nmdays,nydays,date,Init_nmdays	&
                            ,add_dates       ! No. days per year
        use GridValues_ml , only : sigma_bnd,sigma_mid, xp, yp, &
                            fi, GRIDWIDTH_M
	use Io_ml ,only : IO_INFIELD, ios

	implicit none

	integer, intent(in):: numt

!	local

	integer ierr,fid,nr, i, j, k, ident(20)
	integer*8 itmp(6+(MAXLIMAX*MAXLJMAX+3)/4)
	character*20 fname
	integer nyear,nmonth,nday,nhour,nprognosis
	type(date) addhours_to_input
	type(date) next_inptime, buf_date

	real dumhel(MAXLIMAX,MAXLJMAX)

!	definition of the parameter number of the meteorological variables
!	read from field-files:

!ds u7.4vg fl added
!	u(2),v(3),q(9),sdot(11),th(18),cw(22),pr(23),cc3d(39),
!	ps(8),th2m(31),fh(36),fl(37),fm(38),ustar(53),trw(845)

!	scaling factor to be read from field file.

!	data  par/2, 3, 9, 11, 18, 22, 23, 39, 8, 31, 36, 38/

	nr=2
	if(numt == 1)then
	  nr = 1     
	  u  = 0.0
	  v  = 0.0
	endif

!***********************
	if(me.eq.0) then
!***********************

	  fid=IO_INFIELD
	  fname='filxxxx'

	  write(fname,fmt='(''fil'',i4.4)') numt

	  open (fid,file=fname				&
		,form='unformatted',access='sequential'	&
		,status='old',iostat=ios)
          if( ios /= 0 ) then
             write(6,*) "Error opening", fid, fname
             call gc_abort(me,NPROC,"ERROR infield")
          endif !ios
	  ierr=0

	endif ! me == 0

        foundclouds = .false.
        foundustar = .false.
        foundpreta = .false.
        mm5 = .false.

	do while(.true.)

	  call getflti2(fid,ident,itmp,ierr)
	  if(ierr == 2)goto 998

	  k = ident(7)

	  if (ident(6) .eq. 2.and.k.eq.1)then

	    nyear=ident(12)
	    nmonth=ident(13)/100
	    nday=ident(13)-(ident(13)/100)*100
	    nhour=ident(14)/100
	    xp = ident(15)/100.
	    yp = ident(16)/100.
            if(ident(2).eq.1841)then
               GRIDWIDTH_M = 50000.0 ! =~ 1000.*abs(ident(17))/10.
               fi = ident(18)
            elseif(ident(2).eq.1600)then
               mm5 = .true.
               xp = 41.006530761718750 !=~ ident(15)
               yp = 3234.5815429687500 !=~ ident(16)
               !AN = 11888.44824218750
               GRIDWIDTH_M = 1000.*abs(ident(17))/10.
               fi = 10.50000 ! =~ ident(18)
            endif
            nprognosis = ident(4)

	    identi(:)=ident(:)

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!..hj..ko
!.. the name of the fltqhh.yymmdd-file input contains the "ifelt"
!.. time parameters of the analysis, thus the 3 hour prognosis 
!.. data are valid 9-12 hours later!!!!
!
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

!ko... go 9 or 12 hours forward!
!ko    this shift is necessary as the date in the heading is inconsistent 
!ko    with the file name
!ds      Also, the same date appears in two files in the input NWP data,
!ds    hence we need to add 9 hours for nr=even, 12 hours for nr=odd to
!ds    the right dates... :-(

!g11 - introduce new date type
!	initialise first nmonth etc, then add 12(9) hours

!pw use nprognosis=ident(4) for determining the date of the prognosis


          if(numt == 1) then   !!! initialise 
              current_date = date(nyear, nmonth, nday, nhour, 0 )
    	      call Init_nmdays( current_date )
	      addhours_to_input = date(0, 0, 0, nprognosis, 0 )
!pw	      addhours_to_input = date(0, 0, 0, 12, 0 )
	      current_date = add_dates(current_date,addhours_to_input)

              !  if we start 1. of January, then nyear is the year before
              !  so we have to rerun Init_nmdays

	      if(current_date%year.ne.nyear)	&
			call Init_nmdays( current_date )

              ! for printout assign current_date to next_inptime

	      next_inptime = add_dates(current_date,0)

	  else

!   find time for which meteorology is valid:

	      next_inptime = date(nyear, nmonth, nday, nhour, 0 )
              addhours_to_input = date(0, 0, 0, nprognosis, 0 )
	      next_inptime = add_dates(next_inptime,addhours_to_input) 

!	compare the input time with current_date, it should be METSTEP hours later
!	check if current_date+METSTEP = next_inptime

               addhours_to_input = date(0, 0, 0, METSTEP, 0 )
               buf_date = add_dates(current_date,addhours_to_input)

!	now check:

              if(buf_date%year    .ne. next_inptime%year   .or.        &
                 buf_date%month   .ne. next_inptime%month  .or.        &
                 buf_date%day     .ne. next_inptime%day    .or.        &
                 buf_date%hour    .ne. next_inptime%hour   .or.        &
                 buf_date%seconds .ne. next_inptime%seconds     ) then
                     print *,'error in infield, wrong next input time'
                     print *,'current time',current_date
                     print *,'next input time',next_inptime
                     print *,'meteo step',METSTEP

                     call gc_abort(me,NPROC,"infield - time")
               endif

!pw u3	      next_inptime = date(nyear, nmonth, nday, nhour, 0 )
!	now make the correction for the correct input time, but with
!	3 hours less to come to current_date
!
!	      if(mod(numt,2) == 1 )then
!	        addhours_to_input = date(0, 0, 0, 9, 0 )
!	      else
!	        addhours_to_input = date(0, 0, 0, 6, 0 )
!	      endif
!	      next_inptime = add_dates(next_inptime,addhours_to_input) 
!
!	now check:
!
!              if(current_date%year    .ne. next_inptime%year   .or.        &
!                 current_date%month   .ne. next_inptime%month  .or.        &
!                 current_date%day     .ne. next_inptime%day    .or.        &
!                 current_date%hour    .ne. next_inptime%hour   .or.        &
!                 current_date%seconds .ne. next_inptime%seconds     ) then
!                     print *,'error in infield, wrong next input time'
!                     print *,'current time',current_date
!                     print *,'next input time - 3 hours',next_inptime
!                     call gc_abort(me,NPROC,"infield - time")
!               endif
!
!	  add 3 hours to get the real time for printout
!
!               addhours_to_input = date(0, 0, 0, 3, 0 )
!               next_inptime = add_dates(next_inptime,addhours_to_input)
!pw u3


!	now the last check, if we have reached a new year, i.e. current date
!	is 1.1. midnight, then we have to rerun Init_nmdays

               if(  current_date%month == 1 .and.         &
                    current_date%day   == 1 .and.         &
                    current_date%hour  == 0 )	          &
         		call Init_nmdays( current_date )

          end if  ! numt

          if(me == 0) write(6,*) 					&
                '*** nyear,nmonth,nday,nhour,numt,nmdays2,nydays'	&
                ,next_inptime%year,next_inptime%month,next_inptime%day	&
      	        ,next_inptime%hour,numt,nmdays(2),nydays

          endif

	  select case (ident(6))

	  case (2)

          sigma_mid(k)=ident(19)/1.0e+4

	    call getmetfield(ident(20),itmp,dumhel)

	    do j = 1,ljmax
	      do i = 1,limax
		u(i,j,k,nr) = dumhel(i,j)
	      enddo
	    enddo

	  case (3)

	    call getmetfield(ident(20),itmp,dumhel)

	    do j = 1,ljmax
	      do i = 1,limax
		v(i,j,k,nr) = dumhel(i,j)
	      enddo
	    enddo

	  case (9)

	      call getmetfield(ident(20),itmp,q(1,1,k,nr))

	  case (11)

	      call getmetfield(ident(20),itmp,sdot(1,1,k,nr))

	  case (810) !pw u3 MM5 SIGMADOT

	      call getmetfield(ident(20),itmp,sdot(1,1,k,nr))

	  case (18)

	      call getmetfield(ident(20),itmp,th(1,1,k,nr))

	  case (22)

	      call getmetfield(ident(20),itmp,cw(1,1,k,nr))

	  case (23)

              foundpreta = .true. !pw u3
              call getmetfield(ident(20),itmp,pr(1,1,k))

	  case (845) ! pw u3 MM5 TOTALRW

              call getmetfield(ident(20),itmp,trw(1,1,k))

	  case (39)

              foundclouds = .true.
	      call getmetfield(ident(20),itmp,cc3d(1,1,k))

!..2D fields!

	  case (8)

	      call getmetfield(ident(20),itmp,ps(1,1,nr))

	  case (31)

	      call getmetfield(ident(20),itmp,th2m(1,1,nr))

	  case (36)

	      call getmetfield(ident(20),itmp,fh(1,1,nr))

	  case (37)       !ds u7.4vg fl added

              call getmetfield(ident(20),itmp,fl(1,1,nr))

	  case (38)

	      call getmetfield(ident(20),itmp,fm(1,1,nr))

	  case (53) !pw u3

              foundustar = .true.
	      call getmetfield(ident(20),itmp,ustar(1,1,nr))

	  end select

	enddo

998	continue

!     definition of the half-sigma levels from the full levels.

      sigma_bnd(KMAX_BND) = 1.

      do k = KMAX_MID,2,-1
        sigma_bnd(k) = 2.*sigma_mid(k) - sigma_bnd(k+1)
	enddo

      sigma_bnd(1) = 0.

!              set sdot equal to zero at the boundaries.

!ko    sdot is interpolated to get values at half levels; see metvar.f
!ko	    sdot(:,:,1,nr)=0.
!ko
      sdot(:,:,KMAX_BND,nr)=0.

   end subroutine infield

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   subroutine metvar(numt)

	use Par_ml , only : limax,ljmax,li0,li1,lj0,lj1,me		&
			,neighbor,WEST,EAST,SOUTH,NORTH,NOPROC		&
			,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2
        use GridValues_ml , only : xm,xmd, sigma_bnd,sigma_mid
	use ModelConstants_ml, only : PASCAL, PT, CLOUDTHRES, METSTEP,&
                         V_RAIN  !rv1.2
	use PhysicalConstants_ml, only : KARMAN, XKAP, R, CP, GRAV      &
			,ROWATER
	use Tabulations_ml , only : TPI,PBAS,PINC
	implicit none

	integer, intent(in):: numt

!	local

      real prhelp(KMAX_MID)
      real exf1(KMAX_BND), exf2(KMAX_MID) 
 	real bm, cm, dm, divt, x1,x2, xkmin, p1, p2, uvh2, uvhs
	real ri, z00, a2, cdh, fac, fac2, ro, xkh, dvdz, dvdzm, xlmix
	real ric, arg, sl2,dz2k,dex12
	integer i, j, k, lx1,lx2, nr,info
!su	send in x-direction
      real usnd(MAXLJMAX,KMAX_MID)      
!su	send in y-direction
      real vsnd(MAXLIMAX,KMAX_MID)      
!ko
	real prhelp_sum         !u1 - jej ,pr_factor
        real :: inv_METSTEP     ! ds

	nr = 2
	if (numt.eq.1) nr = 1


        divt = 1./(3600.0*METSTEP)
!su
!	define u(i=0,v(j=0)

	if (neighbor(EAST) .ne. NOPROC) then
        do k = 1,KMAX_MID
	    do j = 1,ljmax
	      usnd(j,k) = u(limax,j,k,nr)
	    enddo
	  enddo
        call gc_rsend(MSG_WEST2, MAXLJMAX*KMAX_MID, neighbor(EAST),      &
			info, usnd, usnd)
	endif

	if (neighbor(NORTH) .ne. NOPROC) then
        do k = 1,KMAX_MID
	    do i = 1,limax
	      vsnd(i,k) = v(i,ljmax,k,nr)
	    enddo
	  enddo
        call gc_rsend(MSG_SOUTH2, MAXLIMAX*KMAX_MID, neighbor(NORTH),      &
			info, vsnd, vsnd)
	endif

!     receive from WEST neighbor if any

	if (neighbor(WEST) .ne. NOPROC) then

        call gc_rrecv(MSG_WEST2, MAXLJMAX*KMAX_MID, neighbor(WEST),      &
		info, usnd, usnd)
        do k = 1,KMAX_MID
	    do j = 1,ljmax
	      u(0,j,k,nr) = usnd(j,k)
	    enddo
	  enddo

	else

        do k = 1,KMAX_MID
	    do j = 1,ljmax
	      u(0,j,k,nr) = u(1,j,k,nr)
	    enddo
	  enddo

	endif

!     receive from SOUTH neighbor if any

	if (neighbor(SOUTH) .ne. NOPROC) then

        call gc_rrecv(MSG_SOUTH2, MAXLIMAX*KMAX_MID, neighbor(SOUTH),      &
			info, vsnd, vsnd)
        do k = 1,KMAX_MID
	    do i = 1,limax
	      v(i,0,k,nr) = vsnd(i,k)
	    enddo
	  enddo

	else

        do k = 1,KMAX_MID
	    do i = 1,limax
	      v(i,0,k,nr) = v(i,1,k,nr)
	    enddo
	  enddo

	endif

        inv_METSTEP = 1.0/METSTEP   !ds

	do j = 1,ljmax
	  do i = 1,limax

!     conversion of pressure from mb to Pascal.

	    ps(i,j,nr) = ps(i,j,nr)*PASCAL
	    psurf(i,j) = ps(i,j,nr) !u7.4vg - was psa

!ds rv1.6.2
! surface precipitation, mm/hr

	    surface_precip(i,j) = pr(i,j,KMAX_MID) * inv_METSTEP

!     surface temperature to potential temperature 

	    t2(i,j) = th2m(i,j,nr)  !u7.4vg not -273.15
!su	    th2m(i,j,nr) = th2m(i,j,nr)*(1.e+5/ps(i,j,nr))**(XKAP)
	    th2m(i,j,nr) = th2m(i,j,nr)*exp(-XKAP*log(ps(i,j,nr)*1.e-5))
	    prhelp_sum = 0.0
	    prhelp(1) = max(pr(i,j,1),0.)

	    prhelp_sum = prhelp_sum + prhelp(1)
!ko	    if(numt.ge.2)deriv_2d(i,j,IPRDEP) = deriv_2d(i,j,IPRDEP)
!ko     &			 + prhelp(1)
!ko

!     pr is 3 hours accumulated precipitation in mm in each
!     layer summed from above. This is first converted to precipitation
!     release in each layer.

          do k = 2,KMAX_MID

	      !u4 prhelp(k) = max(pr(i,j,k) - pr(i,j,k-1),0.)
              !u4 - finally fix error reported by hf:
	      prhelp(k) = pr(i,j,k) - pr(i,j,k-1)
!               prhelp(i,j,k) = max(prhelp(i,j,k), 0.)
!ko
	      prhelp_sum = prhelp_sum + prhelp(k)
!ko                if(numt.ge.2)deriv_2d(i,j,IPRDEP) = 
!ko     &			deriv_2d(i,j,IPRDEP) + prhelp(k)
!ko

	    enddo

!ko   accumulated deposition over 3 hour interval 
!ko   k=KMAX_MID now includes accumulated precipitation over all layers
!ko   evaporation has been set to zero as it is not accounted for in the
!ko   wet deposition

    !u1 - error reported by jej, Feb 2002
    !u1 - jej err    if(prhelp_sum.gt.0.0001)then
    !u1 - jej err       pr_factor = pr(i,j,KMAX_MID)/prhelp_sum
    !u1 - jej err    else
    !u1 - jej err      pr_factor = 1.0
    !u1 - jej err    endif

    !u1 - jej err    prhelp(:) = prhelp(:) * pr_factor

!hf I add up in WetDeposition, to have the prec used in the model

	    pr(i,j,:) = prhelp(:)*divt

!ko  interpolation of sigma dot for half layers

          do k = KMAX_MID,2,-1

!ko	  if(me.eq.0.and.numt.eq.1)write(6,*)' k i j sdot(k) sdot(k-1'
!ko     &       ,i,j,k,sdot(i,j,k,1),sdot(i,j,k-1,1)

	      sdot(i,j,k,nr) = sdot(i,j,k-1,nr) 		&
		+ (sdot(i,j,k,nr)-sdot(i,j,k-1,nr))		&
            *(sigma_bnd(k)-sigma_mid(k-1)) / (sigma_mid(k)-sigma_mid(k-1))

!ko	  if(me.eq.0.and.numt.eq.1)write(6,*)' k i j sdot(k) sdot(k-1'
!ko     &       ,i,j,k,sdot(i,j,k,1),sdot(i,j,k-1,1)

  	    enddo
	    sdot(i,j,1,nr)=0.0

            if(foundclouds)then

!ko	conversion from % to fractions (<0,1>) for cloud cover
!ko	calculation of cc3dmax (see eulmc.inc) - 
!ko	maximum of cloud fractions for layers above a given layer

	    cc3d(i,j,1) = 0.01 * cc3d(i,j,1)
	    cc3dmax(i,j,1) = cc3d(i,j,1)
!hf
            lwc(i,j,:)=0.
          do k=2,KMAX_MID
	      cc3d(i,j,k) = 0.01 * cc3d(i,j,k)
	      cc3dmax(i,j,k) = amax1(cc3dmax(i,j,k-1),cc3d(i,j,k-1))
              lwc(i,j,k)=0.6e-6*cc3d(i,j,k)
	    enddo
            else

!pw u3     make cloud cover from cloud water!
               where(cw(i,j,:,nr).gt.CLOUDTHRES)
                  cc3d(i,j,:) = 1.0
               else where
                  cc3d(i,j,:) = 0.0
               end where

               cc3dmax(i,j,1) = cc3d(i,j,1)               
               do k=2,KMAX_MID
                  cc3dmax(i,j,k) = amax1(cc3dmax(i,j,k-1),cc3d(i,j,k-1))
               enddo

            endif


	  enddo
	enddo

!ko 	lines associated with computation of surface diffusion 
!ko	coefficient are commented
 
	xkmin = 1.e-3

!     derive the meteorological parameters from the basic parameters
!     read from field-files.

	do j = 1,ljmax
	  do i = 1,limax
          p1 = sigma_bnd(KMAX_BND)*(ps(i,j,nr) - PT) + PT
	    x1 = (p1 - PBAS)/PINC
	    lx1 = x1
          exf1(KMAX_BND) = tpi(lx1)                  &
			+ (x1-lx1)*(tpi(lx1+1) - tpi(lx1))
          z_bnd(i,j,KMAX_BND) = 0.
          do k = KMAX_MID,1,-1

!     eddy diffusivity in the surface-layer follows the formulation used 
!     in the nwp-model which is based on Louis (1979), (see mc7e.f).

!     the shorter loop is the inner loop to save memory. the order 
!     of the do loops will be changed on a vector machine.

!     exner-function of the half-layers

            p1 = sigma_bnd(k)*(ps(i,j,nr) - PT) + PT
	      x1 = (p1 - PBAS)/PINC
	      lx1 = x1
	      exf1(k) = tpi(lx1) + (x1-lx1)*(tpi(lx1+1) - tpi(lx1))
            p2 = sigma_mid(k)*(ps(i,j,nr) - PT) + PT
	      x2 = (p2 - PBAS)/PINC
	      lx2 = x2

!     exner-function of the full-levels

	      exf2(k) = tpi(lx2) + (x2-lx2)*(tpi(lx2+1) - tpi(lx2))

!     height of the half-layers ! full layers ???

            z_bnd(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
			(exf1(k+1) - exf1(k)))/GRAV

!     height of the full levels. !half layers??

            z_mid(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
			(exf1(k+1) - exf2(k)))/GRAV

            roa(i,j,k,nr) = CP*((ps(i,j,nr) - PT)*sigma_mid(k) + PT)/      &
			(R*th(i,j,k,nr)*exf2(k))

	    enddo

!-----------------------------------------------------------------------

!     local k above the surface layer

	    fac = GRAV/(ps(i,j,nr) - PT)
	    fac2 = fac*fac

          do k = 2,KMAX_MID

            dz2k = z_mid(i,j,k-1)-z_mid(i,j,k)
	      dex12 = th(i,j,k-1,nr)*(exf2(k) - exf1(k)) + 	&
			th(i,j,k,nr)*(exf1(k) - exf2(k-1))

	      dvdzm = 0.00001*dz2k*dz2k
!	      dvdz = (u(i,j,k-1,nr)-u(i,j,k,nr))**2
!     &		 + (v(i,j,k-1,nr) - v(i,j,k,nr))**2
!su	use centred velocities, not yet divided by map factor

	      dvdz = 0.25*((u(i,j,k-1,nr)+u(i-1,j,k-1,nr)	&
			-u(i,j,k,nr)-u(i-1,j,k,nr))**2		&
		+ (v(i,j,k-1,nr)+v(i,j-1,k-1,nr) 		&
			- v(i,j,k,nr) - v(i,j-1,k,nr))**2)	
!pw Error			*xmd(i,j)
          
	      dvdz=amax1(dvdz,dvdzm)
     
!     xlmix, the local mixing height is estimated.
  
            xlmix=KARMAN*z_bnd(i,j,k)
     
	      if(xlmix.gt.70.) xlmix=70.

!     ri is local richardson number

	      ri = GRAV*(th(i,j,k-1,nr)-th(i,j,k,nr))*(exf2(k) -	&
			exf2(k-1))*dz2k/(dvdz*dex12)

!    dvdz is local, vertical windshear

	      dvdz=sqrt(dvdz)/dz2k
  
!     ric is critical richardson number (=1/4 in continuum)

!su            ric=0.115*((z_mid(i,j,k-1) - z_mid(i,j,k))*100.)**0.175
	      ric=0.115*exp(0.175*log(dz2k*100.))

!     xkh is exchange-coeff. for air pollution

	      arg = ri/ric
	      xkh = xkmin
	      sl2 = xlmix*xlmix*dvdz
!hf add
	      if(arg.le.0.)then
!hf ested this .or.( (th(i,j,k-1,nr)-th(i,j,k,nr))<0.05 .and. ri<(1.1/87.) )) then
		xkh=sl2*sqrt(1.1-87.*ri) + xkmin
	      elseif(arg.le.0.5) then
		xkh = sl2*(1.1-1.2*arg) + xkmin
	      elseif(arg.le.1.0) then
		xkh = sl2*(1.-arg) + xkmin
!                write(*,*)'arg le 1',arg,xkh
              endif

!              if(( (th(i,j,k-1,nr)-th(i,j,k,nr))<=0.3).and.(arg.gt.0.0))then
!                xkh=sl2*0.5 + xkmin
!                write(*,*)'New xkh',i,j,k,me,xkh,th(i,j,k-1,nr)-th(i,j,k,nr)
!              endif
!Super adiabat: Theta ved bakken høyere enn Theta i lag k, der 
!l agene k-1,..har lavere pot. T 
!pass på lag 20. ok
!    if ( th(i,j,k,nr)<th(i,j,KMAX_MID,nr).and. &
!         th(i,j,k,nr)>th(i,j,k+1,nr))then
!                xkh = sl2*5.
!vil at alle lagene under også skal ha rask utveksling
!Ri=-2 is large=>(Ri)=13
!               endif

!     exhange coefficient in sigma-coordinates.

!su            ro = ((ps(i,j,1) - PT)*sigma_bnd(k) + PT)*CP*(exf2(k) -
            ro = ((ps(i,j,nr) - PT)*sigma_bnd(k) + PT)*CP*(exf2(k) -      &
			exf2(k-1))/(R*exf1(k)*dex12)
!      if(me == 0)then
!      write(*,*)'Richardson',i,j,k,ri,th(i,j,k-1,nr)-th(i,j,k,nr),ric,xkh
!      endif
!su	      xkz(i,j,k-1) = xkh
!hf TEST	      skh(i,j,k,nr) = xkh*ro*ro*fac2
!Increase diffusivity by 50 percent
!hf NEW SKH!!!	      skh(i,j,k,nr) = xkh*ro*ro*fac2

!     stores the vertical diffusivities

!	      skh(i,j,k,nr) = xkh

!            if (i.eq.it.and.j.eq.jt) write(6,*)'metvar-free troposphere 
!     1           layer:i,j,k,kh,dvdz,ro,ri',i,j,k,xkhw,dvdz,ro,ri
!            if (i.eq.it.and.j.eq.jt) write(6,*)'th,u,v,dvdzm',
!     1           th(i,j,k,nr),u(i,j,k,nr),v(i,j,k,nr),dvdzm

!pw u3 (corrected emep1.2beta)
!    derive precipitations (in mm/s) from rainwater (in kg(water)/kg(air))
!                  pr = V_RAIN * trw * roa  /rowater 
!    trw*roa = kgwater/m3
!    trw * roa  /rowater = vol(water)/vol
!    in one second a "volume" will change height by RAINV. 
!    The volume which comes through the level boundary contains   
!     RAINV*trw*roa/rowater meters of water.
!
! TO BE IMPROVED!

            if(.not.foundpreta)then

              pr(i,j,k) = 1000.*V_RAIN*max(0.0, trw(i,j,k) * roa(i,j,k,nr)  &
                            /ROWATER  )

            endif



	    enddo ! k

            if(.not.foundpreta)then !pw u3
               k=1
               pr(i,j,k) = 1000.*divt*max(0.0, trw(i,j,k) * roa(i,j,k,nr) * &
                            (z_bnd(i,j,k)-z_bnd(i,j,k+1))/ROWATER  )
            endif
               
! ------------------------------------------------------------------

!ko computation of skh in the surface layer is commented
      
!ko          uvh2 = u(i,j,KMAX_MID,nr)*u(i,j,KMAX_MID,nr) +
!ko     1           v(i,j,KMAX_MID,nr)*v(i,j,KMAX_MID,nr)
!su	centred velocities, not yet divided by map factor
!ko          uvh2 = 0.25*( (u(i,j,KMAX_MID,nr)+u(i-1,j,KMAX_MID,nr))**2 +
!ko     1           (v(i,j,KMAX_MID,nr)+v(i,j-1,KMAX_MID,nr))**2 )
!ko     &		*xmd(i,j)

!ko	    if(uvh2.lt.0.0001) uvh2=0.0001

!ko	    uvhs=sqrt(uvh2)

!     local richardson number in the surface layer

!ko          ri = 2.*G*z_mid(i,j,KMAX_MID)*(th(i,j,KMAX_MID,nr)-th2m(i,j,nr))/
!ko     1           ((th(i,j,KMAX_MID,nr) + th2m(i,j,nr))*uvh2)

!     roughness parameter

!ko          z00 = z_mid(i,j,KMAX_MID)/z0(i,j)
!ko	    a2 = alog(z00)
!ko	    a2 = XKAR*XKAR/(a2*a2)

!ko	    if(ri.lt.0.) then
!ko	      cdh = a2*(1.-2.*bm*ri
!ko     1              /(1.+3.*bm*cm*a2*sqrt(-ri*(1.+z00))))
!ko	    else
!ko	      cdh = bm*ri/sqrt(1.+dm*ri)
!ko	      cdh = a2/(1.+3.*cdh)
!ko	    endif

!     convert to sigma-coordinates.

!ko          ro = CP*((ps(i,j,nr) - PT)*sigma_mid(KMAX_MID) + PT)/
!ko     1           (R*th(i,j,KMAX_MID,nr)*exf1(KMAX_BND))
!ko          xkh = cdh*uvhs*z_mid(i,j,KMAX_MID)
!ko          xkz(i,j,KMAX_MID) = xkh
!eb            xkh = (fac*ro*(1.-sigma_mid(KMAX_MID)))
!ko          skh(i,j,KMAX_BND,nr) = xkh*fac2*ro*ro
!            xkh = cdh*uvhs*(fac*ro*(1.-sigma_mid(KMAX_MID)))

!     stores the vertical diffusivity in the surface layer.

!          skh(i,j,KMAX_BND,nr) = xkh

!            if (i.eq.it.and.j.eq.jt) write(6,*)'metvar-surface 
!     1           layer:i,j,cdh,uvhs,ps,z_bnd,z_mid,ri,th2m',i,j,cdh,uvhs,
!     2           ps(i,j,nr),z_bnd(KMAX_MID),z_mid(KMAX_MID),ri,th2m(i,j,nr)

	  enddo
	enddo

!su	do j=lj0,lj1
!su	  do i=li0,li1
!su          roa(i,j,1,nr) = CP*((ps(i,j,nr) - PT)*sigma_mid(1) + PT)/
!su     1           (R*th(i,j,1,nr)*exf2(1))
!su	  enddo
!su	enddo

! 1000 continue

!     Horizontal velocity divided by map-factor.
     
      do k = 1,KMAX_MID

	  do j = 1,ljmax
	    do i = 0,limax

	      u(i,j,k,nr) = 2.*u(i,j,k,nr)/(xm(i,j)+xm(i+1,j))

	    enddo
	  enddo

	  do j = 0,ljmax
	    do i = 1,limax

	      v(i,j,k,nr) = 2.*v(i,j,k,nr)/(xm(i,j)+xm(i,j+1))

	    enddo
	  enddo

	enddo

	return

	end subroutine metvar

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	subroutine metint

!     this routine does the forward linear stepping of the meteorological
!     fields read or derived every 3 hours.

	use ModelConstants_ml , only : nmax,nstep

	implicit none

	real div
     
	if (nstep.lt.nmax) then
 
	  div = 1./real(nmax-(nstep-1))

	  u(:,:,:,1) = u(:,:,:,1) 				&
			+ (u(:,:,:,2) - u(:,:,:,1))*div
	  v(:,:,:,1) = v(:,:,:,1) 				&
			+ (v(:,:,:,2) - v(:,:,:,1))*div
	  sdot(:,:,:,1) = sdot(:,:,:,1) 			&
			+ (sdot(:,:,:,2) - sdot(:,:,:,1))*div
	  th(:,:,:,1) = th(:,:,:,1) 				&
			+ (th(:,:,:,2) - th(:,:,:,1))*div
	  q(:,:,:,1) = q(:,:,:,1) 				&
			+ (q(:,:,:,2) - q(:,:,:,1))*div
	  cw(:,:,:,1) = cw(:,:,:,1) 				&
			+ (cw(:,:,:,2) - cw(:,:,:,1))*div
	  skh(:,:,:,1) = skh(:,:,:,1) 				&
			+ (skh(:,:,:,2) - skh(:,:,:,1))*div
	  roa(:,:,:,1) = roa(:,:,:,1) 				&
			+ (roa(:,:,:,2) - roa(:,:,:,1))*div

	  psurf(:,:) = ps(:,:,1)  !u7.4vg was psa

	  ps(:,:,1) = ps(:,:,1) 				&
			+ (ps(:,:,2) - ps(:,:,1))*div
	  th2m(:,:,1) = th2m(:,:,1) 				&
			+ (th2m(:,:,2) - th2m(:,:,1))*div
!u7.4vg - note we need pressure first
          t2(:,:)    =   th2m(:,:,1) * exp(XKAP*log(psurf(:,:)*1.e-5))

	  fh(:,:,1) = fh(:,:,1) 				&
			+ (fh(:,:,2) - fh(:,:,1))*div
!ds u7.4vg fl added
	  fl(:,:,1) = fl(:,:,1) 				&
			+ (fl(:,:,2) - fl(:,:,1))*div

	  fm(:,:,1) = fm(:,:,1) 				&
			+ (fm(:,:,2) - fm(:,:,1))*div
          if(foundustar)then
	  ustar(:,:,1) = ustar(:,:,1) 				&
			+ (ustar(:,:,2) - ustar(:,:,1))*div
          endif

!ko  precipitation and cloud cover are no longer interpolated

!gv	  if(MADE)then				!MADE
!gv	    vd1komp(:,:,:,1) = vd1komp(:,:,:,1) 
!gv     &			+ (vd1komp(:,:,:,2) - vd1komp(:,:,:,1))*div
!gv	    vd50komp(:,:,:,1) = vd50komp(:,:,:,1) 
!gv     &			+ (vd50komp(:,:,:,2) - vd50komp(:,:,:,1))*div
!gv	  endif					!MADE

	else

!     assign the the meteorological data at time-level 2 to level 1 for
!     the next 6 hours integration period before leaving the inner loop.

	  u(:,:,:,1) = u(:,:,:,2)
	  v(:,:,:,1) = v(:,:,:,2)
	  sdot(:,:,:,1) = sdot(:,:,:,2)
	  th(:,:,:,1) = th(:,:,:,2)
	  q(:,:,:,1) = q(:,:,:,2)
	  cw(:,:,:,1) = cw(:,:,:,2)
	  skh(:,:,:,1) = skh(:,:,:,2)
	  roa(:,:,:,1) = roa(:,:,:,2)

!su	don't forget psa !!!!
	  !7.4vg, but put after ps update psa(:,:) = ps(:,:,1)

	  ps(:,:,1) = ps(:,:,2)
	  psurf(:,:) = ps(:,:,1)   ! u7.4vg
	  th2m(:,:,1) = th2m(:,:,2)

!u7.4vg - note we need pressure first
          t2(:,:)    =   th2m(:,:,1) * exp(XKAP*log(psurf(:,:)*1.e-5))

	  fh(:,:,1) = fh(:,:,2)
	  fm(:,:,1) = fm(:,:,2)
!ds u7.4vg fl added
	  fl(:,:,1) = fl(:,:,2)
	  if(foundustar) ustar(:,:,1) = ustar(:,:,2)

!gv	  if(MADE)then				!MADE
!gv	    vd1komp(:,:,:,1) = vd1komp(:,:,:,2)
!gv	    vd50komp(:,:,:,1) = vd50komp(:,:,:,2)
!gv	  endif					!MADE

	endif

	end subroutine metint

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   subroutine MetModel_LandUse(callnum)

  !ds rv1.2 combines old subroutines in_isnowc and inpar
  ! (and commented out SetZ0)
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   !     This subroutine reads parameterfields from file
   !     reading surface roughness classes from file: rough.170
   !     reading snow                      from file: rough.170
   !
   !     ... fields as used in meteorological model
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    use ModelConstants_ml, only : current_date
    use Io_ml, only :   IO_SNOW, IO_ROUGH,  ios, open_file  
    use ReadField_ml, only : ReadField ! reads ascii fields
    implicit none

    integer, intent(in) :: callnum
    integer ::  i,j, err
    real, allocatable, dimension(:,:) :: r_class  ! Roughness (real) 
    character*20 fname

    ios = 0

    if ( callnum == 1  ) then 

        if ( me == 0  ) then
           write(fname,fmt='(''rough.170'')') 
           write(6,*) 'filename for landuse ',fname
        end if

        allocate(r_class(MAXLIMAX,MAXLJMAX),stat=err)
        if ( err /= 0 ) call gc_abort(me,NPROC,"alloc err:rough")

        call ReadField(IO_ROUGH,fname,r_class)
   
       ! And convert from real to integer field
      
        do j=1,ljmax
           do i=1,limax
              iclass(i,j)=nint(r_class(i,j)) 
           enddo
        enddo
        deallocate(r_class,stat=err)
        if ( err /= 0 ) call gc_abort(me,NPROC,"dealloc err:rough")
    else ! callnum == 2
        if (me == 0) then
           write(fname,fmt='(''snowc'',i2.2,''.dat'')') current_date%month
           write(6,*) 'filename for snow ',fname

        endif !me==0

        call ReadField(IO_SNOW,fname,snow)

    end if ! callnum == 1 

  !ds rv1.2 commented out code put here for possible future use.
  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !  subroutine SetZ0()
   !
   !    ! This routine is not used? (emep1.2beta)
   !    !
   !    !  Set the surface roughness equal to the lam50e values.
   !    !  (ds comment - all this rougness stuff has to be reviewed/replaced 
   !    !   with new  deposition modules, at least that from RIVM)
   !
   !    integer              :: i, j, icl
   !    real, dimension(0:6) ::  class   ! Some surface roughnesses ??
   !
   !    class =  (/ 1.0e-4,1.0e-3,3.0e-1,3.0e-1,3.0e-1,3.0e-1,1.0e-3 /)
   !
   !    do j = 1,ljmax
   !      do i = 1,limax
   !         icl = iclass(i,j)
   !         z0(i,j) = class(icl)
   !      end do
   !    end do
   !
   !  end subroutine SetZ0
 end subroutine MetModel_LandUse
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!Originally tiphys.f in hmix.f
!imax=>limax(sums) ,MAXLIMAX(dimensions)
!ljmax=>lljmax,MAXLJMAX
!z2    =>z_mid
!z1    =>z_bnd
!kmax2 =>KMAX_MID
!kmax1 =>KMAX_BND
!scor1 =>sigma_bnd
!scor2 =>sigma_mid
!pt=>PT
!g=>GRAV
!xkar=>KARMAN
!sm(1)=>sm
!cp=>CP
!pi=>PI
!xkar=>KARMAN

      subroutine tiphys(numt)
        use ModelConstants_ml, only : KMAX_MID ,KMAX_BND,PT
        use GridValues_ml , only : sigma_bnd,sigma_mid
        use Par_ml, only :limax,ljmax,MAXLIMAX,MAXLJMAX
        use PhysicalConstants_ml, only :CP,PI,KARMAN,GRAV,XKAP,R
!z_mid, z_bnd Met_ml
!c
!c	file: met.eulmet-mnd.f
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c	written by Trond Iversen,  modified by Hugo Jakobsen, 060994
!c
!c	Called from: eulmain.f
!c
!c-----------------------------------------------------------------
!c
!!        This routine calculates the exner function, 
!!        the  geopotential height, and the vertical exchange coefficient 
!!        in sigma surfaces.
!!        The height zi of the "well mixed layer" or ABL-height 
!!        is also calculated.
!c
!c
!c	if nroa = 1 also roa is calculated.
!c	if nfirst=1 for the initial timelevel
!c
!c
!c-----------------------------------------------------------------
!c	files included:
!c
!c	metpar.inc	: parameters for LAM50E-meteorology
!c	metvar.inc	: dependent meteorological variables
!c	metcon.inc	: model constants 
!c
!c-----------------------------------------------------------------
!c	routines called:
!c
!c		smoosp
!c
!c
!c-----------------------------------------------------------------
!c
!c    DescriPTion of the parameters/variables defined in this file:  
!c
!c
!c	absfac	: |xfac|
!c	abshd	: |fm|
!c	amax1	: fortran function, choosing largest value
!c	amin1	: fortran function, choosing smallest value
!c	CP	: heat capaciyt of air at constant pressure, J/(kg K)
!c	delq	: available heat flux for developing the unstable ABL, J/m2
!!              : heat-input per m2 from the ground during unstable BL
!c	deltaz	: zm(i,k) - zm(i,k+1), m
!c	dpidth	: heat increasement in accordance with temp. increasement, J/m2
!c	dth	: iterative increament in potential temperature
!c	dth0	: accumulated increament in iterative temperature
!c	dtz	: time interwall for integration of surface heat fluxes
!c		  in the ABL-height calculations, s
!c	dvdz	: Wind shear, 1/s
!c	eps	: small number avoiding ri to become infinitely large
!c	exfrco	: parameter in the Kz model
!c	exnm	: exner function in the full sigma-levels, J/(kg K)
!c	exns	: exner function in the half sigma-levels, J/(kg K)
!c	fh	: surface flux of sensible heat, W/m2
!c	fl	: surface flux of sensible heat, W/m2 ! ds u7.4vg
!c	fm	: surface stress (flux of momentum), N/m2
!c	g	: gravitational acceleration, m/s2
!c	hs	: height of surface layer (i.e. prandtl-layer), m
!c	hsl	: (= hs/l, where l is the monin-obhukov length)
!c	i	: grid index in x-direction
!c	iip	: limax + 1
!c	limax	: max number of grid points in x-direction
!c	iznew	: index for new value of the ABL-height, m
!c	izold	: index for previous value of the ABL-height, m
!c	j	: grid index in y-direction
!c	jjp	: ljmax + 1
!c	ljmax	: max number of grid points in y-direction
!c	k	: grid index in vertical-direction
!c	kkk	: helping index for the cycling of ABL-height
!c	kkm	: number of full s-levels, *** not used ***
!c	KMAX_BND	: max number of vertical half levels
!c		  in sigma coordinates
!c	KMAX_MID	: max number of vertical full levels
!c		  in sigma coordinates
!c       kzmax   : maximum value of xksig, m2/s
!c       kzmin   : minimum value of xksig, m2/s
!c	ndth	: do variable for convective ABL-height iteration loop
!c	nh1	: counts number of layers below zlimax
!c	nh2	: counts number of layers with Kz > ( Kz )limit
!c	nr	: number of met.fields stored in arrays (= 1 or 2)
!c	nt	: time counting variable of the outer time-loop
!c	p	: local pressure, hPa (mb)
!c       pi      : pi = 4.*atan(1.) = 3.14 ...
!c	pidth	: heat used to adjust air temperature, J/m2
!c	pref	: refference pressure (at ground level), 1.e+5 Pa
!c	ps	: surface pressure, hPa
!c	PT	: pressure at the top of the model atmosphere, hPa (mb)
!c	pz	: local pressure in half sigma levels, hPa (mb),
!c		  helping array (j - slices) for pressure
!c	pzpbl	: stores H(ABL) for averaging and plotting purposes, m
!c	ri	: richardson`s number
!c	ri0	: critical richardson`s number
!c	risig	: richardson's number in sigmas-levels
!c	roas	: air density at surface, kg/m3
!c	sigma_bnd	: height of the half-sigma layers
!c	sigma_mid	: height of the full-sigma layers
!c	sm	: height of the surface layer in s-coordinates (4% of H(ABL), m
!c	th	: potensial temperature (theta), K
!c	th2m	: potensial temperature at 2m height, K
!c	thadj	: adjustable surface temperature, K
!c	thsrf	: potensial temperature at the surface, K
!c	trc	: helping variable telling whether or not unstable ABL exists
!!              :       0 => no need for further calc. of ziu
!!              :       1 => ziu not found yet. 
!c	u	: wind speed in the x-direction, m/s
!c       umax    : maximum value of u and v, m/s
!c	ustar	: friction velocity, m/s
!c	v	: wind speed in the y-direction, m/s
!c       ven     : ventilation coefficient, m3
!c       venav   : time averaged ventilation coefficient, m3
!c       venmax  : maximum value of ven, m3
!c       venmin  : minimum value of ven, m3
!c       ven00   : averaged ventilation coefficient at 00 UTC, m3
!c       ven06   : averaged ventilation coefficient at 06 UTC, m3
!c       ven12   : averaged ventilation coefficient at 12 UTC, m3
!c       ven18   : averaged ventilation coefficient at 18 UTC, m3
!c	vdfac	: factor for reduction of vD(1m) to vD(hs)
!!              : i.e. factor for aerodynamic resistance towards dry deposition
!!              : vd(50m) = vd(1m)/(1 + vd(1m)*vdfac)
!c	x12	: mixing length squared, m2
!c	xfac	: helping variable for reducing concentrations to 1m values
!c	xfrco	: parameter in the Kz model
!c	XKAP	: r/CP (-)
!c	KARMAN	: von Karmans constant
!c	xkdz	: the vertical derivative of xkhs at hs, m/s
!!              : i.e. vertical gradient of xkhs
!c	xkhs	: diffusivity at hs (in surface layer), m2/s
!!              : i.e. vertical exchange coeff. on top of prandtl-layer 
!c	xksig	: estimated exchange coefficient, Kz,  in intermediate 
!c		  sigma levels, m2/s
!c	xksm	: spacially smoothed Kz in z direction, m2/s.
!!              : xksig smoothed over three adjacent layers
!c	xkzi	: local helping array for the vertical diffusivity, m2/s
!!              : i.e. vertical exchange coeff. on top of ABL for unstable BL
!c       xtime   : 6.*3600. (seconds in one term, six hours)
!c	zi	: Height of ABL (final value), m
!c	zlimax	: maximum value of ABL-height, zi, (2000), m
!c	zimhs	: ziu - hs
!c	zimin	: minimum value of ABL-height, zi, (200), m
!c	zimz	: ziu - zs1
!c	zis	: height of the stable ABL, m
!c	ziu	: height of the unstable ABL, m
!c	zixx	: Height og ABL (intermediate value), m
!c	zm	: geopotential height of full sigma levels above topography, m
!c	zmhs	: zs1 - hs
!c	zs1	: geopotential height of  half sigma levels above topography, m
!c	ztop	: height of the uppermost layer in s-coordinates
!c
!c-------------------------------------------------------------------
!c
!c	underlying models for (including litterature references):
!c	----------------------------------------------------------
!c
!c	local air density (based on the ideal gas law):
!c
!c	roa = ro = p/(r*T) = (CP*p)/(r*th*exf)
!c
!c
!c	local pressure:
!c
!c	p = PT + scor_*(ps - PT)
!c
!c
!c	local values of the exner function:
!c
!c	exnm = CP*(p/pref)**XKAP
!c
!c
!c	in the unstable (convective) surface layer:
!c	___________________________________________
!c
!c
!c			Iversen and Nordeng(1987), pp. 16-17, eq.(3.18)
!c
!c
!c	xfrco = 0.5*(sqrt(6859.0)-1), 	eq.(3.22), p.17
!c
!c	exfrco = 1.0/3.0
!c
!c	ustar = (fm/roas)**(1/2),		uxo = amax1(ustar, 0.00001)
!c
!c
!c	Monin-Obhukov length:
!c
!c	defined
!c 
!c
!c 	L = theta*ustar**3/(KARMAN*GRAV*<w*theta>),
!c	
!c		where 
!c			<w*theta> = ustar*theta0
!c
!c	  = theta*ustar**2/(KARMAN*GRAV*theta0),
!c
!c	    	where 
!c			fh = CP*ro*<w*theta> = CP*ro*ustar*theta0 
!c			   = - CP*ro*Kz*d(theta)/dz
!c
!c	  = ps*ustar**3*CP/(xhar*GRAV*fh*r)
!c
!c		where theta is set equal to T at the surface; ps=ro*r*T
!c
!c
!c	
!c	the local eddy diffusivity, Iversen and Nordeng(1987), 
!c				    eq.(3.21)
!c
!c	z/L > -2 :
!c	-----------
!c
!c	xkhs = ustar*KARMAN*z/phi(z/L)
!c		
!c		where	phi(z/L) = 0.74/(1-9*((z/L))**0.5
!c
!c	xkdz  = xkhs*(1.0-4.5*(z/L)/(1.0-9.0*(z/L)))/z
!c 
!c
!c	z/L <= -2 :
!c	-----------
!c
!c	xkhs = ustar*KARMAN*z/phi(z/L)
!c		
!c		where	phi(z/L) = 0.74/(1-xfrco*((z/L))**exfrco
!c
!c	xkdz  = xkhs*(1.0-xfrco*(z/L)/(3.0*(1.0-xfrco*(z/L)))/z
!c
!c
!c
!c-->	in the unstable boundary layer:
!c	_______________________________
!c
!c
!c	the profile formulae for the exchange coefficients used
!c	when the surface layer heat flux is directed upwards
!c	(fh<0, upwards, away from the surface), is the 
!c	O'Brien (1970)-relation.
!c
!c
!c			Iversen and Nordeng(1987), pp. 16, eq.(3.17)
!c
!c
!c	Kz	=	z*Kz(hs)/hs	if z < hs
!c 	
!c		=	zs1*xkhs/hs	if zs1 <= hs
!c
!c	and
!c
!c	Kz	=	((zi-z)/(zi-hs))**2*{dK/dz*(hs)-Kz(zi)
!c			+ (z-hs)[Kz(hs)+2*(Kz(hs)-Kz(zi))/(zi-hss)]}
!c
!c					if hs <= z < zi
!c
!c		=	xkzi*((zimz/zimhs)**2*(xkhs-xkzi
!c			+ zmhs*(xkdz + 2.0*(xkhs-xkzi/zimhs))
!c
!c					if hs <= zs1 < ziu
!c
!c
!c
!c
!c-->	in the stable surface layer (Stull (1988), p.361, eq.(9.4.1d)):
!c	_______________________________________________________________
!c
!c	the atmospheric boundary layer is defined as stable as a whole 
!c	when the surface layer heat flux is directed downwards.
!c	However, local unstable layers are possible.
!c
!c	xkhs = ustar*KARMAN*z/phi(z/L)
!c
!c		where	phi(z/L) = 0.74 + 4.7*(z/L)
!c
!c
!c
!c-->	above the boundary layer (i.e. in the free atmosphere),
!c	and in the stable boundary layer:
!c	________________________________________________________
!c
!c
!c	local mixing height, Iversen and Nordeng(1987), 
!c				p. 15, eq.(3.13)
!c
!c	xlmix = KARMAN*z	, if z <= zm
!c	xlmix = KARMAN*zm	, if z >  zm,	where zm = 200 m
!c
!c
!c	local richardson number (Louis (1979), eq.(23), p.193):
!c
!c	ri = GRAV*d(z)*d(th)/(th*d(v)**2)
!c
!c
!c	critical richardson number, Iversen and Nordeng(1987), 
!c				    p. 16, eq.(3.16); modified !!!
!c
!c	ric0 = A*(d(z)/z0)**B, 	where	A = 0.115, B = 0.175,
!c					d(z) is the model layer thickness
!c					d(z0) = (1/100) m
!c
!c	ric = amax1(0.25, ric0), 	i.e. ric,min = 0.25
!c
!c
!c	local eddy diffusivity, Iversen and Nordeng(1987), p. 15
!c				and Louis(1979), eq.(21), p.193
!c
!c
!c	xksig = xl2*|u|*F(ri)/|z|
!c
!c	where
!c		xl2   = (KARMAN*amin1(zs1, zimin))**2
!c
!c		F(ri) = sqrt(1.1-87*ri),   if ri <= 0
!c		F(ri) = 1.1-ri/ric,	   if 0 < ri <= 0.5*ric
!c		F(ri) = 1.0-ri/ric,	   if 0.5*ric < ri < ric
!c		F(ri) = 0.001		   if ri >= ric
!c
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c.. the ventilation coefficient:
!c
!c
!c     ven = |U|*A*dt
!c
!c     where
!c
!!          |U| = sqrt(u**2 + v**2)
!!           A  = 2.*h/sqrt(pi)*zi
!!           dt = 6.*60.*60.
!c
!c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c............................................
!c..the following sketches the sigma-surfaces:
!c
!c
!!                ///////////////////
!c    sigma_bnd(1) = 0 - -sigmas - - - - - sdot(1) = 0, xksig(1)=xksm(1)=0, 
!!                                           pr(1)=0,PT,exns(1), zs1(1)
!c
!c        sigma_mid(1) ---sigmam---------- u, v, th, q, cw, exnm (1)
!c
!c
!c        sigma_bnd(2) - - - - s - - - - - sdot(2), xksig(2), exns(2), pr(2) 
!!                                                 zs1(2), xksm(2)
!c
!c        sigma_mid(2) --------m---------- u, v, th, q, cw, exnm (2)
!c
!c
!c        sigma_bnd(3) - - - - s - - - - - sdot(3), xksig(3), exns(3), pr(3)
!!                                                 zs1(3), xksm(3)
!c
!c        sigma_mid(3) --------m---------- u, v, th, q, cw, exnm (3)
!c
!c
!c        sigma_bnd(4) - - - - s - - - - - sdot(4), xksig(4), exns(4), pr(4)
!!                                                  zs1(4), xksm(4)
!c
!c        sigma_mid(4) --------m---------- u, v, th, q, cw, exnm (4)
!c
!c
!c        sigma_bnd(5) - - - - s - - - - - sdot(5), xksig(5), exns(5), pr(5)
!!                                                  zs1(5), xksm(5)
!c
!!                        :
!!                        :
!c
!c  sigma_bnd(KMAX_BND-1) - - - - s - - - - - sdot(KMAX_BND-1), xksig(KMAX_MID), 
!!                                    exns(KMAX_BND-1),zs1(KMAX_BND-1), 
!!                                    pr(KMAX_BND-1),xksm(KMAX_MID)
!c
!c    sigma_mid(KMAX_MID) --------m---------- u, v, th, q, cw, exnm (KMAX_MID); 
!!                                    this level is assumed to be
!!                                    the top of Prandtl-layer (LAM50E)
!c
!c sigma_bnd(KMAX_BND) = 1- - - - s - - - - - sdot(KMAX_BND) = 0, ps, th2m, fh, 
!!                ///////////////////        fm, mslp, xksig(KMAX_MID)=0, 
!!                                           exns(KMAX_BND), zs1(KMAX_BND), 
!!                                           pr(KMAX_BND),xksm(KMAX_BND)=0.
!c
!c
!c..alternativ names:   kkin = KMAX_MID=20 (number of sigma-layers)
!!                     kkinp = kkin+1 = KMAX_BND=21 (number of 
!!                                        level-bounds for layers) 
!!                     kkinm = kkin-1 = KMAX_MID-1                
!c
!c**********************************************************************
!ds rv1_6_x
 logical, parameter :: DEBUG_KZ = .false.
 logical, parameter :: PIELKE_KZ = .true.
 logical :: debug_flag   ! set true when i,j match DEBUG_i, DEBUG_j



!definer alle dimensjoner med MAXLIMAX,MAXLJMAX
      real exnm(MAXLIMAX,MAXLJMAX,KMAX_MID),exns(MAXLIMAX,MAXLJMAX,KMAX_BND),zm(MAXLIMAX,KMAX_MID),&
          zs1(MAXLIMAX,MAXLJMAX,KMAX_BND),risig(MAXLIMAX,KMAX_BND),xksm(MAXLIMAX,KMAX_BND),&
          zis(MAXLIMAX),ziu(MAXLIMAX,MAXLJMAX),delq(MAXLIMAX),thsrf(MAXLIMAX),trc(MAXLIMAX),&
          pidth(MAXLIMAX),dpidth(MAXLIMAX),lim,help(MAXLIMAX,MAXLJMAX),a(MAXLIMAX,MAXLJMAX),&
          zixx(MAXLIMAX,MAXLJMAX),xdthdz,dthdz(MAXLIMAX,KMAX_MID),zmmin,&
          deltaz(MAXLIMAX,KMAX_MID),pz(MAXLIMAX,KMAX_BND),thc(MAXLIMAX,KMAX_MID),&
          roas(MAXLIMAX,MAXLJMAX),zimin,zlimax,kzmin,kzmax,uabs(MAXLIMAX,MAXLJMAX),&
          vdfac(MAXLIMAX,MAXLJMAX),xkhs(MAXLIMAX),xkdz(MAXLIMAX),xkzi(MAXLIMAX),hs(MAXLIMAX),&
!hf new
          sm,pref,xtime,umax,eps,ric,ric0,dthdzm,dthc,xdth,xfrco,exfrco,hsl,dtz,p,&
          dvdz,xl2,uvhs,zimhs,zimz,zmhs,ux0,fac,fac2,dex12,ro,xkh100(MAXLIMAX)
!hf Hilde&ANton
      real hsurfl
!hf new
      integer i,j,k,km,km1,km2,kabl,iip,jjp,numt,kp
         


      integer nh1(MAXLIMAX),nh2(MAXLIMAX),nr

!hf new
      iip = limax+1
      jjp = ljmax+1

!
!     Preliminary definitions
!
      nr = 2
      if (numt.eq.1) nr = 1
      pref = 1.e+5
!from ModelC      pi = 4.*atan(1.)

      xtime = 6.*3600.
      umax=+70.
      zlimax = 3000.
      zimin = 100.
      zmmin = 200.
      kzmin = 1.e-3
      kzmax = 1.e3
!hf      venmin=0.
!hf      venmax= sqrt(umax**2 + umax**2)*xtime*2.*h/sqrt(PI)*zlimax
      eps = 0.01
      dtz = 3600.
      sm = 0.04
!c
!c
!c..preset=zero:
      do 5 k=1,KMAX_MID
      do 5 i=1,MAXLIMAX
            xksm(i,k)=0.
            risig(i,k)=0.
      do 5 j=1,MAXLJMAX
            xksig(i,j,k)=0.
 5    continue
!c..................................
!c..exner-function in the full sigma-levels..
!c
      do 10 k=1,KMAX_MID
      do 10 j=1,ljmax
      do 10 i=1,limax
!c
!c..pressure (pa)
        p = PT + sigma_mid(k)*(ps(i,j,nr) - PT)
!c..exner (j/k kg)
        exnm(i,j,k)=CP*(p/pref)**XKAP
!c..temp (K):
!hf never used        temp(i,j,k)=th(i,j,k,nr)*exnm(i,j,k)/CP
!c
!c..thadj (K):
!c
!hf never used         thadj(i,j,k)=th(i,j,k,nr)
!hf      if (k.eq.KMAX_MID) thadj(i,j,KMAX_BND) = th2m(i,j,nr)
 10   continue
!c

!c.........................................
!c..procedure to arrive at mixing height..:
!c
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c
!c     Start j-slice here.
!c
      lim = 0.1
!c..hj..test lim=1.
      do 40 j=1,ljmax
!c
!c..exner in half-sigma levels:
         do 15 k=1,KMAX_BND
         do 15 i=1,limax
            p = PT + sigma_bnd(k)*(ps(i,j,nr) - PT)
            pz(i,k) = p
            exns(i,j,k)=CP*(p/pref)**XKAP  
!c
 15      continue
!c
!c.. exns(KMAX_BND), th(KMAX_BND) and height of sigmas:
         do 16 i=1,limax
            zs1(i,j,KMAX_BND)=0.
 16      continue
!c
!c     Height of the half levels
!c
         do 17 k=KMAX_BND-1,1,-1 
         do 17 i=1,limax
            zs1(i,j,k)=zs1(i,j,k+1)+th(i,j,k,nr)*&
                 (exns(i,j,k+1)-exns(i,j,k))/GRAV

 17      continue
!c
!c..height of sigma:
         do 18 k=1,KMAX_MID
         do 18 i=1,limax
            zm(i,k)=((exnm(i,j,k)-exns(i,j,k))*zs1(i,j,k+1)&
                   +(exns(i,j,k+1)-exnm(i,j,k))*zs1(i,j,k))&
                  /(exns(i,j,k+1)-exns(i,j,k))    
!c
 18      continue
!c
!c
!c
!c
!c----------------------------------------------------------------------
!c...........................................
!c..the following variables in sigmas-levels:
!c
         do 19 k=2,KMAX_MID
            km=k-1
         do 19 i=1,limax
!c
!c.........................
!c..wind sheare
!c
!HF Slightly different formulation of dvdz than in metvar
        dvdz = ( (u(i,j,km,nr)-u(i,j,k,nr))**2 &
              + (v(i,j,km,nr)-v(i,j,k,nr))**2 + eps)
!c
!c.........................
!c..the richardsons number:
!c
!!           risig(i,k)=(2.*GRAV/(th(i,j,km,nr)+th(i,j,k,nr)))*
!c     +          (th(i,j,km,nr)-th(i,j,k,nr))*(zm(i,km)-zm(i,k))
!c     +       / ( (u(i,j,km,nr)-u(i,j,k,nr))**2 + 
!c     +         (v(i,j,km,nr)-v(i,j,k,nr))**2 + eps )
!c
!c..modified
!c
            risig(i,k)=(2.*GRAV/(th(i,j,km,nr)+th(i,j,k,nr)))*&
               (th(i,j,km,nr)-th(i,j,k,nr))*(zm(i,km)-zm(i,k))&
               /dvdz
!c........................
!c..mixing length squared:
!c
            xl2=(KARMAN*amin1(zs1(i,j,k),zmmin))**2

!c
!c..............................
!c..critical richardsons number:
!c
            ric0=0.115*((zm(i,km)-zm(i,k))*100.)**0.175
            ric=amax1(0.25,ric0)
!c

!c.............
!c..wind shear:
!c
!!           dvdz = ( (u(i,j,km,nr)-u(i,j,k,nr))**2 
!c     +              + (v(i,j,km,nr)-v(i,j,k,nr))**2 )**0.5
!c     +             /(zm(i,km)-zm(i,k))
!c
!c..modified
!c
            dvdz = sqrt(dvdz)/(zm(i,km)-zm(i,k))

!c..................................................................
!c..exchange coefficient (Pielke,...)
 !ds alternative
     if ( PIELKE_KZ  ) then
          if (risig(i,k) > ric ) then
              !xkds(k) = 0.1
              xksig(i,j,k) = 0.1
          else
              !xkds(k) = 1.1 * (ric-risig(i,k)) * xl2 * dvdz /ric
              xksig(i,j,k) = 1.1 * (ric-risig(i,k)) * xl2 * dvdz /ric
          end if
     else

!c..exchange coefficient (blackadar, 1979; iversen & nordeng, 1987):
!c
            if(risig(i,k).le.0.) then
               xksig(i,j,k)=xl2*dvdz*sqrt(1.1-87.*risig(i,k))
            elseif(risig(i,k).le.0.5*ric) then
               xksig(i,j,k)=xl2*dvdz*(1.1-1.2*risig(i,k)/ric)
            elseif(risig(i,k).le.ric) then
               xksig(i,j,k)=xl2*dvdz*(1.-risig(i,k)/ric)
            else
               xksig(i,j,k)=0.001
            endif
    end if ! Pielke or Blackadar
!c
 19      continue
!c
!ctttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
!c
!c---------------------------------------------------------------------
!c..................................
!c..height of stable boundary layer:
!c
!c.........................................................
!c..vertical smoothing of xksig over three adjacent layers:
!c
            k=2 
            km=1
            kp=3
            do i=1,limax
               xksm(i,k)=( (zm(i,km)-zm(i,k))*xksig(i,j,k)&
                         + (zm(i,k)-zm(i,kp))*xksig(i,j,kp) )&
                       / ( zm(i,km) - zm(i,kp) )
            enddo
!c
            k=KMAX_MID 
            km2=k-2
            km1=k-1
            do i=1,limax
               xksm(i,k)=( (zm(i,km2)-zm(i,km1))*xksig(i,j,km1)&
                         + (zm(i,km1)-zm(i,k))*xksig(i,j,k) )&
                       / ( zm(i,km2) - zm(i,k) )
            enddo   
!c
            do k = 3,KMAX_MID-1
               km1=k-1
               km2=k-2
               kp=k+1
               do i=1,limax
               xksm(i,k)=(  (zm(i,km2)-zm(i,km1))*xksig(i,j,km1)&
                         + (zm(i,km1)-zm(i,k))*xksig(i,j,k)&
                         + (zm(i,k)-zm(i,kp))*xksig(i,j,kp) )& 
                       / ( zm(i,km2) - zm(i,kp) )
               enddo
            enddo


!c
!c............................................................
!c..The height of the stable BL is the lowest level for which:
!c..xksm .le. 1 m2/s (this limit may be changed):
!c
	do i = 1,limax
	   zis(i)=zimin
           nh1(i) = KMAX_MID
           nh2(i) = 1        	
        enddo
!c
	do 25 k=KMAX_MID,2,-1
	do 25 i=1,limax

	if(xksm(i,k).ge.lim .and. nh2(i).eq.1) then
           nh1(i)=k
        else
           nh2(i)=0
        endif
!c
  25	continue
!c
        do 26 i=1,limax
!c
        k=nh1(i)
!c
	if(zs1(i,j,nh1(i)).ge.zimin) then

           if( abs(xksm(i,k)-xksm(i,k-1)) .gt. eps) then 

              zis(i)=((xksm(i,k)-lim)*zs1(i,j,k-1) &
                  + (lim-xksm(i,k-1))*zs1(i,j,k))&
                  /(xksm(i,k)-xksm(i,k-1))
           else

              zis(i)=zimin
           endif

	endif
!c
  26    continue
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c---------------------------------------------------------------------
!c....................................
!c..height of unstable boundary layer:
!c
!c
!c..assuring that th is increasing with height.
!c..adjusted th-sounding is assigned to thc-array.
!c..This adjusted th is not meant to be used in
!c..other parts of the model program
!c
	dthdzm = 1.e-4 
	do i =1,limax
	   thc(i,KMAX_MID)=th(i,j,KMAX_MID,nr)
		do k=KMAX_MID-1,1,-1
!c
		   dthc=(th(i,j,k,nr)-th(i,j,k+1,nr))&
     			/(zm(i,k)-zm(i,k+1))
!c
		   dthdz(i,k)=amax1(dthc,dthdzm)
!c
	   thc(i,k)=thc(i,k+1)+dthdz(i,k)*(zm(i,k)-zm(i,k+1))
!c
		enddo
	enddo

!c
!c
!c..estimated as the height to which an hour's input
!c..of heat from the ground is vertically distributed,
!c..assuming dry adiabatic adjustment.
!c
!c
         do 260 i=1,limax
!c
            delq(i)=-amin1((fh(i,j,nr)),0.)*dtz
!c
            thsrf(i)=0.
            ziu(i,j)=0.
!c
!c.................................
!c..trc=1 for unstable BL (delq>0):
!c..   =0 for stable BL (delq=0):
!c
            if(delq(i).gt.0.00001) then
               trc(i)=1.
            else
               trc(i)=0.
            endif
!c
!c------------------------------------------------------------
!c calculating the height of unstable ABL
!c
	if(trc(i).eq.1.) then
!c
	kabl = KMAX_MID
!c	
 28		if(trc(i).eq.1.) then
!c
		kabl = kabl-1
		pidth(i)=0.
!c
!c	endres til 29 etter at algoritmene er sammenlignet?
!c
		do 281 k=KMAX_MID,kabl,-1
		xdth = thc(i,kabl)-thc(i,k)
		dpidth(i) = exnm(i,j,k)*xdth*(pz(i,k+1)-pz(i,k))/GRAV
		pidth(i) = pidth(i) + dpidth(i)
 281		continue
!c
!c
	  		if(pidth(i).ge.delq(i).and.trc(i).eq.1.) then

!c	at level kabl or below level kabl and above level kabl+1


			  	thsrf(i) = thc(i,kabl)-&
     				         (thc(i,kabl)-thc(i,KMAX_MID))*&
     				         (pidth(i)-delq(i))/pidth(i)

                    	  	xdthdz=(thc(i,kabl)-thc(i,kabl+1))/&
                                      (zm(i,kabl)-zm(i,kabl+1))

                          	ziu(i,j) = zm(i,kabl+1) +  &
                                    (thsrf(i)-thc(i,kabl+1))/xdthdz

			trc(i)=0.

			endif


			if(kabl.le.4 .and. trc(i).eq.1.) then

				write(6,*)'ziu calculations failed!'

				ziu(i,j)=zlimax

				trc(i)=0.
			endif


		go to 28			  

		endif

	endif

 260      continue

!c..iteration finished
!c.....................................................................


         do 35 i=1,limax

            zixx(i,j)=amax1(ziu(i,j),zis(i))
            zixx(i,j)=amin1(zlimax,zixx(i,j))

 35      continue


!     End j-slice


40   continue
            
!..spatial smoothing of new zi:

44 format(I2,30F5.0)
      call smoosp(zixx,zimin,zlimax)

      do j=1,ljmax
         do i=1,limax
            pzpbl(i,j) = zixx(i,j)
         enddo
      enddo


!cttttttttttttttttttttttttttttttttttttttttttttttttttttttt
!c..height of ABL finished..............................
!c------------------------------------------------------
!c......................................................
!c..exchange coefficients for convective boundary layer:
!c..o'brien's profile formula:
!c..and the air density at ground level:
!c
!c..constants for free-convection limit:
!c
      xfrco=0.5*(sqrt(6859.)-1)
      exfrco=1./3.

      do 70 j=1,ljmax

!c..exchange parameter and its vertical derivative at z = hs

      do 60 i=1,limax

         xkh100(i)=0.  !Hilde&Anton
         xkhs(i)=0.                                            
         xkdz(i)=0.
         xkzi(i)=0.
!c
!c
!c...................................................................
!c..air density at ground level is always calculated diagnostically:
!c
         roas(i,j)=ps(i,j,nr)/(XKAP*th2m(i,j,nr)*exns(i,j,KMAX_BND))


         if(ziu(i,j).ge.zimin) then
!c
!c..........................
!c..unstable surface-layer.:
!co
!c..height of surface layer
            hs(i)=sm*ziu(i,j)
!c..u*
!c
         if(foundustar)then                      
            ux0 = ustar(i,j,1)
         else
            ux0 = sqrt(fm(i,j,nr)/roas(i,j))
         endif

            ux0=amax1(ux0,0.00001)

!c..hsl=hs/l where l is the monin-obhukov length
            hsl=KARMAN*GRAV*hs(i)*fh(i,j,nr)*XKAP &
             /(ps(i,j,nr)*ux0*ux0*ux0)


        !ds rv1_7_2 changes: use simple Garratt \Phi function
        !   instead of "older" Businge and Iversen/Nordeng stuff:

	! old code distinguished free convection from less unstable
	! we don't bother now.
            !dsif(hsl >= 0  ) then           ! Unstable
               xkhs(i)=ux0*KARMAN*hs(i)*sqrt(1.0-16.0*hsl)/1.00   
               xkdz(i)=xkhs(i)*(1.-5.0*hsl/(1.0-16.0*hsl))/hs(i)        
            !dselse
            !dsif(hsl.ge.-2.) then           
            !ds   xkhs(i)=ux0*KARMAN*hs(i)*sqrt(1.-9.*hsl)/0.74   
            !ds   xkdz(i)=xkhs(i)*(1.-4.5*hsl/(1.-9.*hsl))/hs(i)        
            !dselse

!                                              
!c..free convection :                                            
            !ds   xkhs(i)=ux0*KARMAN*hs(i)*(1.-xfrco*hsl)**exfrco/0.74  
            !ds   xkdz(i)=xkhs(i)*(1.-xfrco*hsl/(3.*(1.-xfrco*hsl)))/hs(i)
            !dsendif
!Hilde&Anton
!pw & hf            hsurfl=KARMAN*GRAV*100.*amax1(0.001,fh(i,j,nr))*XKAP&
!pw & hf                 &             /(ps(i,j,nr)*ux0*ux0*ux0)
            hsurfl=KARMAN*GRAV*100.*fh(i,j,nr)*XKAP&
                 &             /(ps(i,j,nr)*ux0*ux0*ux0)

            !ds rv1_7_2
            !if(hsurfl >= -2.) then
               xkh100(i)=ux0*KARMAN*100.*sqrt(1.-16.*hsurfl)/1.00
            !else
            !   xkh100(i)=ux0*KARMAN*100.*(1.-xfrco*hsurfl)**exfrco/0.74
            !endif

            Kz_min(i,j)=xkh100(i)
            xksig(i,j,KMAX_MID)=xkhs(i)

         else
!c
!c..........................
!c..stable surface-layer...:
!c
!c..height of surface layer
            hs(i)=sm*zixx(i,j)
!c..u*
            ux0=sqrt(fm(i,j,nr)/roas(i,j))
            ux0=amax1(ux0,0.00001)
!c
!c..hsl=hs/l where l is the monin-obhukov length
            hsl=KARMAN*GRAV*hs(i)*amax1(0.001,fh(i,j,nr))*XKAP&
             /(ps(i,j,nr)*ux0*ux0*ux0)


            !xksig(i,j,KMAX_MID)=ux0*KARMAN*hs(i)/(0.74+4.7*hsl)   
            xksig(i,j,KMAX_MID)=ux0*KARMAN*hs(i)/(1.00+5.0*hsl)   

         endif
!hf Hilde&Anton
            hsurfl=KARMAN*GRAV*100.*amax1(0.001,fh(i,j,nr))*XKAP&
                 &             /(ps(i,j,nr)*ux0*ux0*ux0)
            !Kz_min(i,j)=1.35*ux0*KARMAN*100./(0.74+4.7*hsurfl)
            Kz_min(i,j)=ux0*KARMAN*100./(1.00+5.0*hsurfl)
!c
!c...............................................................
!c..factor for reduction of dry-deposition speed from 1m to hs..:
!c
!NOT NEEDED
!         abshd=abs(fh(i,j,nr))
!         if(abshd.lt.0.1) then
!            uvhs=(u(i,j,KMAX_MID,nr)*u(i,j,KMAX_MID,nr)&
!               + v(i,j,KMAX_MID,nr)*v(i,j,KMAX_MID,nr))**0.5
!            vdfac(i,j)=0.74*uvhs/(ux0*ux0)
!         else
!            xfac=(th(i,j,KMAX_MID,nr)/th2m(i,j,nr) - 1.)*ps(i,j,nr)&
!                /(XKAP*abshd)
!            vdfac(i,j)=abs(xfac)
!         endif

 60   continue
!c
!c
!c..exchange parameter at z = ziu
!c
      do 65 k=1,KMAX_MID
      do 65 i=1,limax

         if(ziu(i,j).gt.zimin .and. zs1(i,j,k).ge.ziu(i,j)) then
            xkzi(i)=xksig(i,j,k)
         elseif (ziu(i,j).gt.zimin) then
!c
!c.....................................................   
!c..the obrien-profile for z<ziu                      . 
!c.....................................................                  
!c
            if(zs1(i,j,k).le.hs(i)) then   
               xksig(i,j,k)=zs1(i,j,k)*xkhs(i)/hs(i)          
            else                                                      
               zimhs=ziu(i,j)-hs(i)   
               zimz=ziu(i,j)-zs1(i,j,k)                     
               zmhs=zs1(i,j,k)-hs(i)                  
               xksig(i,j,k)=xkzi(i)+(zimz/zimhs)*(zimz/zimhs)  &  
                 *(xkhs(i)-xkzi(i)+zmhs*(xkdz(i)+&
                 2.*(xkhs(i)-xkzi(i))/zimhs))
            endif  
         endif

 65   continue

 70   continue
!c
!c..spatial smoothing of xksig:
!c
!c

!pw emep1.2      do 80 k=1,KMAX_MID
      do 80 k=2,KMAX_MID

         do i=1,limax
            do j=1,ljmax
!hf Anton&Hilde
!pw emep1.2               if ( (pzpbl(i,j)>z_mid(i,j,k+1)) .and. k>1 )then
!hf new               if ( (pzpbl(i,j)>z_mid(i,j,k-1)) )then
               if ( (pzpbl(i,j)>z_mid(i,j,k)) )then
                  xksig(i,j,k)=max(xksig(i,j,k),Kz_min(i,j))
               endif 
               help(i,j) = xksig(i,j,k)
            enddo
         enddo

       call smoosp(help,kzmin,kzmax)

         do i=1,limax
            do j=1,ljmax
               xksig(i,j,k) = help(i,j)
!hf I need to convert to sigma coordinates
! Kz(sigma)=Kz*ro^2*(GRAV/p*)
! I have exns=exf1,exnm=exf2,
! and use the formulation from original Unified model
!Ikke saerlig effektivt, men enkelt....

	    fac = GRAV/(ps(i,j,nr) - PT)
	    fac2 = fac*fac
            dex12 = th(i,j,k-1,nr)*(exnm(i,j,k) - exns(i,j,k)) + 	&
			th(i,j,k,nr)*(exns(i,j,k) - exnm(i,j,k-1))
            ro = ((ps(i,j,nr) - PT)*sigma_bnd(k) + PT)*CP*(exnm(i,j,k) -      &
			exnm(i,j,k-1))/(R*exns(i,j,k)*dex12)
!hf NEW
               skh(i,j,k,nr) = xksig(i,j,k)*ro*ro*fac2
            enddo
         enddo

 80   continue


!c
!c...............................................................
!c..mixing-layer parameterization finished.......................
!c...............................................................
!c
!c
!c..calculating the ventilation coefficient, ven(i,j)
!hf NOT NEEDED
!      do i=1,limax
!         do j=1,ljmax
!
!            a(i,j) = 2.0*h/sqrt(PI)
!
!            uabs(i,j)=0.
!
!            do k=KMAX_MID,2,-1
!               if(zixx(i,j).ge.zs1(i,j,k)) then
!
!                  dz = zs1(i,j,k)-zs1(i,j,k+1)
!
!                  u2 = ( u(i,j,k,nr)**2 + v(i,j,k,nr)**2 )
!
!                  uabs(i,j)=uabs(i,j)+sqrt(u2)*dz
!
!               elseif (zs1(i,j,k).gt.zixx(i,j) &
!                     .and. zs1(i,j,k+1).lt.zixx(i,j)) then
!
!                    dz = zixx(i,j)-zs1(i,j,k+1)
!
!                    u2 = ( u(i,j,k,nr)**2 + v(i,j,k,nr)**2 )
!
!                    uabs(i,j)=uabs(i,j)+sqrt(u2)*dz
!
!                endif
!
!            enddo
!
!            ven(i,j) = uabs(i,j)*xtime*a(i,j)
!
!         enddo
!      enddo
!c
!c..spatial smoothing of ven:
!c
!hf not needed      call smoosp(ven,limax,ljmax,iip,jjp,venmin,venmax)

      end subroutine tiphys

!c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c
!c
      subroutine smoosp(f,rmin,rmax)      

!c	file: eulmet-mnd.f
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c	written by Trond Iversen,  modified by Hugo Jakobsen, 080994
!       parallellized and modified by Peter February 2003
!
!c
!c	Called from: tiphys.f
!c
!c----------------------------------------------------------------------
!c
!c	This routine applies the shapiro filter with s=0.5 and s=-0.5         
!c	to the field f usinh h as a work space also the boundaries 
!c	are smoothed. f contains the smoothed field upon return.
!c

!c    Definition of the variables:  
!c
!c
!c	f	: data to be smoothed
!c	iif	: =limax
!c	jjf	: =ljmax
!c	h1,h2	: = help variable
!c	rmin	: min allowed
!c	rmax	: max allowed
!c
implicit none
!       real f(MAXLIMAX,MAXLJMAX),h1(0:MAXLIMAX+1,0:MAXLJMAX+1),rmin,rmax
       real, intent(inout):: f(MAXLIMAX,MAXLJMAX)
       real, intent(in):: rmin,rmax
       real, dimension(MAXLIMAX+4,MAXLJMAX+4):: h1, h2

    integer iif,jjf,is,i,j,ii,jj,iifl,jjfl
    real s
    real, dimension(MAXLIMAX,2) ::f_south,f_north
    real, dimension(MAXLJMAX+2*2,2) ::f_west,f_east
    integer  thick

    iif=limax
    jjf=ljmax

    thick=2  !we fetch 2 neighbors at once, so that we don't need to call readneighbours twice
    iifl=iif+2*thick
    jjfl=jjf+2*thick

    call readneighbors(f,f_south,f_north,f_west,f_east,thick)

    do j=1,jjf
       jj=j+thick
       do i=1,iif
          ii=i+thick
          h1(ii,jj) = f(i,j)
       enddo
    enddo
    do j=1,thick
       do i=1,iif
          ii=i+thick
          h1(ii,j) = f_south(i,j)
       enddo
    enddo
    do j=1,thick
       jj=j+jjf+thick
       do i=1,iif
          ii=i+thick
          h1(ii,jj) = f_north(i,j)
       enddo
    enddo
    do j=1,jjfl
       do i=1,thick
          h1(i,j) = f_west(j,i)
       enddo
    enddo
    do j=1,jjfl
       do i=1,thick
          ii=i+iif+thick
          h1(ii,j) = f_east(j,i)
       enddo
    enddo

    do j=1,jjfl
       h2(1,j) = 0.
       h2(iifl,j) = 0.
    enddo
    do i=1,iifl
       h2(i,1) = 0.
       h2(i,jjfl) = 0.
    enddo
 44 format(I2,30F5.0)

     do 30 is=2,1,-1                                          

         s=is-1.5  !s=0,5 s=-0.5
         if(is /= 2)h1=h2
                                           
!..the smoothing

         do 20 j=2,jjfl-1
            do 20 i=2,iifl-1
               h2(i,j)=(1.-2.*s+s*s)*h1(i,j)&                              
                    +0.5*s*(1.-s)*(h1(i+1,j)+h1(i-1,j)+h1(i,j+1)+h1(i,j-1))   &  
                    +s*s*(h1(i+1,j+1)+h1(i-1,j-1)+h1(i+1,j-1)+h1(i-1,j+1))/4. 
               h2(i,j) = amax1(h2(i,j),rmin)
               h2(i,j) = amin1(h2(i,j),rmax)
         20 continue                                        
               
30  continue                                      

    do j=1,jjf
       jj=j+thick
       do i=1,iif
          ii=i+thick
          f(i,j)=h2(ii,jj) 
       enddo
    enddo
                          
      return  
      end subroutine smoosp                                                             
!  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   subroutine readneighbors(data,data_south,data_north,data_west,data_east,thick)
     
     
! Read data at the other side of the boundaries 
!
! thick is the number of gridcells in each direction to be transferred
! Note that we also fetch data from processors in the "diagonal"
! directions
!
! Written by Peter February 2003
!
!Comment from Peter after ds bug-fix:
!The data_west(jj,:)=data(1,j) is not a bug: when there is no west neighbour, 
!the data is simply copied from the nearest points: data_west(jj,:) should 
!be =data(-thick+1:0,j), but since this data does not exist, we 
!put it =data(1,j).

!     use Par_ml , only : me,NPROC,limax,ljmax,MAXLIMAX,MAXLJMAX &
!          ,neighbor,SOUTH,NORTH,WEST,EAST,NOPROC

     implicit none

     integer, intent(in) :: thick
     real,intent(in), dimension(MAXLIMAX,MAXLJMAX) ::data
     real,intent(out), dimension(MAXLIMAX,thick) ::data_south,data_north
     real,intent(out), dimension(MAXLJMAX+2*thick,thick) ::data_west,data_east

     integer :: msgnr,info
     integer :: i,j,tj,jj,jt
     
 !check that limax and ljmax are large enough
     if(limax < thick .or. ljmax < thick)then
        print *,'me, limax, ljmax, thick ', me, limax, ljmax, thick
        call gc_abort(me,NPROC,"ERROR readneighbors")
     endif
     
     msgnr=1

     data_south(:,:)=data(:,1:thick)
     data_north(:,:)=data(:,ljmax-thick+1:ljmax)
     if(neighbor(SOUTH) >= 0 )then
        call gc_rsend(msgnr,MAXLIMAX*thick,		&
             neighbor(SOUTH), info, data_south, data_south)
     endif
     if(neighbor(NORTH) >= 0 )then
        call gc_rsend(msgnr+9,MAXLIMAX*thick,		&
             neighbor(NORTH), info, data_north, data_north)
     endif
     
     if(neighbor(SOUTH) >= 0 )then
        call gc_rrecv(msgnr+9,MAXLIMAX*thick,		&
             neighbor(SOUTH), info, data_south, data_south)
     else
        do tj=1,thick
           data_south(:,tj)=data(:,1)
        enddo
     endif
 44 format(I2,30F5.0)
     if(neighbor(NORTH) >= 0 )then
        call gc_rrecv(msgnr,MAXLIMAX*thick,		&
             neighbor(NORTH), info, data_north, data_north)
     else
        do tj=1,thick
           data_north(:,tj)=data(:,ljmax)
        enddo
     endif
     
     jj=0
     do jt=1,thick
        jj=jj+1
        !ds bug ? data_west(jj,:)=data_south(:,jt)
        ! may be wrong! Check also assignments below.
        data_west(jj,:)=data_south(1:thick,jt)
        data_east(jj,:)=data_south(limax-thick+1:limax,jt)
     enddo
     do j=1,ljmax
        jj=jj+1
        !ds bug ? data_west(jj,:)=data(:,j)
        data_west(jj,:)=data(1:thick,j)
        data_east(jj,:)=data(limax-thick+1:limax,j)
     enddo
     do jt=1,thick
        jj=jj+1
        !ds bug data_west(jj,:)=data_north(:,jt)
        data_west(jj,:)=data_north(1:thick,jt)
        data_east(jj,:)=data_north(limax-thick+1:limax,jt)
     enddo
     
     if(neighbor(WEST) >= 0 )then
        call gc_rsend(msgnr+3,(MAXLJMAX+2*thick)*thick,		&
             neighbor(WEST), info, data_west, data_west)
     endif
     if(neighbor(EAST) >= 0 )then
        call gc_rsend(msgnr+7,(MAXLJMAX+2*thick)*thick,		&
             neighbor(EAST), info, data_east, data_east)
     endif
     
     if(neighbor(WEST) >= 0 )then
       call gc_rrecv(msgnr+7,(MAXLJMAX+2*thick)*thick,		&
             neighbor(WEST), info, data_west, data_west)
     else
        jj=0
        do jt=1,thick
           jj=jj+1
           data_west(jj,:)=data_south(1,jt)
        enddo
        do j=1,ljmax
           jj=jj+1
           data_west(jj,:)=data(1,j)
        enddo
        do jt=1,thick
           jj=jj+1
           data_west(jj,:)=data_north(1,jt)
        enddo
     endif
     if(neighbor(EAST) >= 0 )then
        call gc_rrecv(msgnr+3,(MAXLJMAX+2*thick)*thick,		&
             neighbor(EAST), info, data_east, data_east)
     else
        jj=0
        do jt=1,thick
           jj=jj+1
           data_east(jj,:)=data_south(limax,jt)
        enddo
        do j=1,ljmax
           jj=jj+1
           data_east(jj,:)=data(limax,j)
        enddo
        do jt=1,thick
           jj=jj+1
           data_east(jj,:)=data_north(limax,jt)
        enddo
     endif
     
   end subroutine readneighbors
end module met_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
!  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



