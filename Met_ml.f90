!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module Met_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! ds rv1.2 MetModel_LandUse added here for snow and iclass
!  - combined from hf and pw Met_ml
! October 2001 hf added call to ReadField:::: WITH TKE scheme (Octobar 2004)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
use Functions_ml,      only : Exner_tab, Exner_nd    ! ds apr2005
use GridValues_ml,     only : xm,xmd, sigma_bnd,sigma_mid, &
                               i_glob, j_glob
use ModelConstants_ml, only : PASCAL, PT, CLOUDTHRES, METSTEP, &
                              KMAX_BND,KMAX_MID,NMET, &
                              DEBUG_i, DEBUG_j
use Par_ml           , only : MAXLIMAX,MAXLJMAX,NPROC,me  &
                               ,limax,ljmax,li0,li1,lj0,lj1  &
                               ,neighbor,WEST,EAST,SOUTH,NORTH,NOPROC  &
                               ,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2
use PhysicalConstants_ml, only : KARMAN, KAPPA, R, RGAS_KG, CP, GRAV, ROWATER
!ds apr2005 use Tabulations_ml , only : TPI,PBAS,PINC
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
                           ,q         ! Specific humidity


!    since pr,cc3d,cc3dmax used only for 1 time layer - define without NMET
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: &
                     pr      & ! Precipitation
                    ,cc3d    & ! 3-d cloud cover (cc3d),
                    ,cc3dmax & ! and maximum for layers above a given layer
                    ,lwc     & !liquid water content
                    ,sst       ! SST Sea Surface Temprature- ONLY from 2002



! surface fields
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,NMET) :: &
                     ps      & ! Surface pressure hPa (or Pa- CHECK!)
                    !ds apr2005 ,th2m    & ! Temp 2 m   deg. K
                    ,t2_nwp    & ! Temp 2 m   deg. K
                    ,fh      & ! surf.flux.sens.heat W/m^2   ! ds u7.4vg added
                    ,fl      & ! latent heat flux W/m^2
                    ,tau       ! surf. stress  N/m^2
                               ! (ds - was fm - renamed with apr2005 system)





! Derived met
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND,NMET) :: skh
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID,NMET) :: roa ! kg/m^3
  real,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &
                  psurf & !u7.4lu psa  Surface pressure hPa
                 ,surface_precip    & ! Surface precip mm/hr   ! ds rv1.6.2
                 !ds apr2005 ,t2&      !u7.4vg temp2m  Temp 2 m   deg. K
                 ,u_ref !wind speed

!ds apr2005 - now calculated in met_derived:

  real,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &
                  rho_surf    & ! Surface density
                 ,invL_nwp    & ! inverse Obukhov length
                 ,Tpot2m      & ! Potential temp at 2m
                 ,ustar_nwp     ! friction velocity m/s ustar^2 = tau/roa

  real,public, save, &
      dimension(MAXLIMAX,MAXLJMAX,KMAX_BND) :: z_bnd ! height of full layers
  real,public, save, &
      dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: z_mid ! height of half layers

!ds rv1.2  keep HIRLAM/xx met model landuse stuff in same routine
  integer,public, save, dimension(MAXLIMAX,MAXLJMAX) :: & 
       snow,    &  ! monthly snow (1=true), read in MetModel_LandUse
       iclass      ! roughness class ,         "       "       "

  logical, private, parameter ::  MY_DEBUG = .false.
  logical, private, save      ::  debug_proc = .false.
  integer, private, save      :: debug_iloc, debug_jloc  ! local coords

  logical, public, save :: foundustar !pw Used for MM5-type, where u* but not tau
  logical, public, save :: sdot_at_mid !pw rv1_9_24 . set false if sdot 
               !is defined (when read) at level  boundaries and therefore 
               !do not need to be interpolated.

!hf tiphys
!check dimension
  real,public, save, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: &
                xksig ! estimated exchange coefficient, Kz, in intermediate 
                      ! sigma levels, m2/s

  real,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &
       pzpbl,   &  !stores H(ABL) for averaging and plotting purposes, m
       Kz_min      ! Min Kz below hmix  !hf Hilde&Anton


 !ds apr2005 - temnporary placement of solar radiation variations

    real, public, dimension(MAXLIMAX, MAXLJMAX), save:: &
        zen          &  ! Zenith angle (degrees)
       ,coszen=0.0   &  ! cos of zenith angle
       ,Idiffuse     &  ! diffuse solar radiation (W/m^2)
       ,Idirect         ! total direct solar radiation (W/m^2)

 public :: infield,Meteoread_CDF
 public :: MetModel_LandUse    ! rv1.2 combines old in_isnowc and inpar
 public :: metvar
 public :: metint
 public :: tiphys           !hf NEW

 contains

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      subroutine Meteoread_CDF(numt)

!	the subroutine reads meteorological fields and parameters (every
!	METSTEP-hours) from NetCDF fields-files, divide the fields into 
!       domains	and sends subfields to the processors

	use ModelConstants_ml , only : current_date & ! date-type
                                       , METSTEP      ! Meteo read timestep
        use GridValues_ml , only : sigma_bnd,sigma_mid, sigma_mid, xp, yp, &
                            fi, GRIDWIDTH_M,xm,ref_latitude
        use Dates_ml, only : date,Init_nmdays,add_dates,nmdays
 
	implicit none

	integer, intent(in):: numt

        character (len = 100),save ::meteoname !name of the meteofile
        character (len = 100) :: namefield & ! name of the requested field
                   ,period_of_validity !field is either instaneous or averaged
        integer ::ndim,nyear,nmonth,nday,nhour,k
        integer,save :: Nhh &    !number of field stored per 24 hours
                   ,nhour_first& !time of the first meteo stored
                   ,nrec     !nrec=record in meteofile, for example 
                             !(Nhh=8): 1=00:00 2=03:00 ... 8=21:00
                             !if nhour_first=3 then 1=03:00 2=06:00...8=24:00

        integer :: nr!Fields are interpolate in time: between nr=1 and nr=2
   	type(date) next_inptime,addhours_to_input

	nr=2 !set to one only when the first time meteo is read
	if(numt == 1)then !first time meteo is read
	  nr = 1     
          
          sdot_at_mid = .true.
          foundustar = .false.
          
          nyear=2001
          nmonth=1
          nday=1
          nhour=0
          current_date = date(nyear, nmonth, nday, nhour, 0 )
          call Init_nmdays( current_date )
          
         !*********initialize grid parameters*********
56        FORMAT(a5,i4.4,i2.2,i2.2,a3)
          write(meteoname,56)'meteo',nyear,nmonth,nday,'.nc'
          write(*,*)'reading meteo from ',trim(meteoname)

          call Getgridparams(meteoname,GRIDWIDTH_M,xp,yp,fi,xm,&
               ref_latitude,sigma_mid,Nhh,nyear,nmonth,nday,nhour,nhour_first)
          if(me==0)then
             write(*,*)'sigma_mid:',(sigma_mid(k),k=1,20)
             write(*,*)'grid resolution:',GRIDWIDTH_M
             write(*,*)'xcoordinate of North Pole, xp:',xp
             write(*,*)'ycoordinate of North Pole, yp:',yp
             write(*,*)'longitude rotation of grid, fi:',fi
             write(*,*)'true distances latitude, ref_latitude:',ref_latitude
          endif

! If origin of meteodomain does not coincide with origin of large domain,
! xp and yp should be shifted here, and coordinates must be shifted when 
! meteofields are read (not yet implemented)

	endif ! numt==1

        if(numt==1)then
           next_inptime = current_date
        else
           addhours_to_input = date(0, 0, 0, METSTEP, 0 )
           next_inptime = add_dates(current_date,addhours_to_input)
        endif
        nyear=next_inptime%year
        nmonth=next_inptime%month
        nday=next_inptime%day 
        nhour=next_inptime%hour
        if(  current_date%month == 1 .and.         &
             current_date%day   == 1 .and.         &
             current_date%hour  == 0 )	          &
             call Init_nmdays( current_date )
        
        if(me == 0) write(6,*) 					&
             '*** nyear,nmonth,nday,nhour,numt,nmdays2'	&
             ,next_inptime%year,next_inptime%month,next_inptime%day	&
             ,next_inptime%hour,numt,nmdays(2)
        
!Read rec=1 in case 00:00 from 1st January is missing
        if((numt-1)*METSTEP<=nhour_first)nrec=0 
        nrec=nrec+1
        if(nrec>Nhh.or.nrec==1)then! define a new meteo input file
           write(meteoname,56)'meteo',nyear,nmonth,nday,'.nc'
           if(me==0)write(*,*)'reading ',trim(meteoname)
           nrec = 1
           !could open and close file here instead of in Getmeteofield
        endif
        if(me==0)write(*,*)'nrec,nhour=',nrec,nhour

! 3D fields (i,j,k)
        ndim=3
        namefield='u_wind'
!note that u and v have dimensions 0:MAXLIJMAX instead of 1:MAXLIJMAX  
!u(i=0) and v(j=0) are set in metvar
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity,u(1:MAXLIMAX,1:MAXLJMAX,:,nr))
        namefield='v_wind'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity,v(1:MAXLIMAX,1:MAXLJMAX,:,nr))
        namefield='specific_humidity'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, q(:,:,:,nr))
        namefield='sigma_dot'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, sdot(:,:,:,nr))
        namefield='potential_temperature'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, th(:,:,:,nr))
        namefield='precipitation'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, pr(:,:,:))
        namefield='3D_cloudcover'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, cc3d(:,:,:))
        if(trim(period_of_validity)/='averaged')then
           if(me==0)write(*,*)'WARNING: 3D cloud cover is not averaged'
        endif
! 2D fields (surface) (i,j)
        ndim=2
        namefield='surface_pressure'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, ps(:,:,nr))
        namefield='temperature_2m'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, t2_nwp(:,:,nr))
        namefield='surface_flux_sensible_heat'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, fh(:,:,nr))
        if(period_of_validity=='averaged')fh(:,:,1)=fh(:,:,nr)
        namefield='surface_flux_latent_heat'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, fl(:,:,nr))
        if(period_of_validity=='averaged')fl(:,:,1)=fl(:,:,nr)
        namefield='surface_stress'
        call Getmeteofield(meteoname,namefield,nrec,ndim,&
             period_of_validity, tau(:,:,nr))
         if(period_of_validity=='averaged')tau(:,:,1)=tau(:,:,nr)
      

      end subroutine Meteoread_CDF

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
                            fi, GRIDWIDTH_M,ref_latitude
	use Io_ml ,only : IO_INFIELD, ios

	implicit none

	integer, intent(in):: numt

!	local

	integer ierr,fid,nr, i, j, k, ij, ident(20)
	integer*8 itmp(6+(MAXLIMAX*MAXLJMAX+3)/4)
	character*20 fname
	integer nyear,nmonth,nday,nhour,nprognosis
	type(date) addhours_to_input
	type(date) next_inptime, buf_date

	real dumhel(MAXLIMAX,MAXLJMAX)
        real ::xrand(20)

!	definition of the parameter number of the meteorological variables
!	read from field-files:

!ds u7.4vg fl added
!	u(2),v(3),q(9),sdot(11),th(18),cw(22),pr(23),cc3d(39),
!	ps(8),t2_nwp(31),fh(36),fl(37),tau(38),ustar(53),trw(845), sst(103)

!	scaling factor to be read from field file.

!	data  par/2, 3, 9, 11, 18, 22, 23, 39, 8, 31, 36, 38, 103/

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

        sdot_at_mid = .true.
        foundustar = .false.
        
	do while(.true.)

	  call getflti2(fid,ident,itmp,ierr)
	  if(ierr == 2)goto 998

	  k = ident(7)

	  if (ident(6) .eq. 2.and.k.eq.1)then

	    nyear=ident(12)
	    nmonth=ident(13)/100
	    nday=ident(13)-(ident(13)/100)*100
	    nhour=ident(14)/100
            if(ident(17)>0.0)then
               xp = ident(15)/100.
               yp = ident(16)/100.
            else
               xp = ident(15)/1.
               yp = ident(16)/1.
            endif
            fi = ident(18)
            if(ident(2).eq.1841)then
               GRIDWIDTH_M = 50000.0 ! =~ 1000.*abs(ident(17))/10.
            else
               GRIDWIDTH_M = 1000.*abs(ident(17))/10.
               if(me==0)write(*,*)'GRIDWIDTH_M=' ,GRIDWIDTH_M ,&
                    'AN= ',6.370e6*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M
            endif
            if(ident(2).eq.1600)then
               xp = 41.006530761718750 !=~ ident(15)
               yp = 3234.5815429687500 !=~ ident(16)
               fi = 10.50000 ! =~ ident(18)
               if(me==0)write(*,*)ident(15),ident(16),ident(18),xp,yp,fi
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

!pw use nprognosis=ident(4) for determining the date of the prognosis


          if(numt == 1) then   !!! initialise 
              current_date = date(nyear, nmonth, nday, nhour, 0 )
    	      call Init_nmdays( current_date )
	      addhours_to_input = date(0, 0, 0, nprognosis, 0 )
!pw	      addhours_to_input = date(0, 0, 0, 12, 0 )
	      current_date = add_dates(current_date,addhours_to_input)!

              !  if we start 1. of January, then nyear is the year before
              !  so we have to rerun Init_nmdays!

	      if(current_date%year.ne.nyear)	&
			call Init_nmdays( current_date )

              ! for printout assign current_date to next_inptime!

	      next_inptime = add_dates(current_date,0)!

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
		u(i,j,k,nr) = dumhel(i,j)-1.E-9
!pw the "-1.E-9" is put in order to avoid possible different 
!roundings on different machines. 
	      enddo
	    enddo

	  case (3)

	    call getmetfield(ident(20),itmp,dumhel)

	    do j = 1,ljmax
	      do i = 1,limax
		v(i,j,k,nr) = dumhel(i,j)-1.E-9
!pw the "-1.E-9" is put in order to avoid possible different 
!roundings on different machines. 
	      enddo
	    enddo

	  case (9)

	      call getmetfield(ident(20),itmp,q(1,1,k,nr))

	  case (11)

	      call getmetfield(ident(20),itmp,sdot(1,1,k,nr))

	  case (-11) !pw (not standard convention)

	      call getmetfield(ident(20),itmp,sdot(1,1,k,nr))

              sdot_at_mid = .false.

 	  case (810) !pw u3 MM5 SIGMADOT

	      call getmetfield(ident(20),itmp,sdot(1,1,k,nr))

	  case (18)

	      call getmetfield(ident(20),itmp,th(1,1,k,nr))

!	  case (22)
!
!	      call getmetfield(ident(20),itmp,cw(1,1,k,nr))

!          case (26)                                          !ASSYCON
!              call getmetfield(ident(20),itmp,ccc(1,1,k,nr)) !ASSYCON
               
	  case (23)

              call getmetfield(ident(20),itmp,pr(1,1,k))

!	  case (845) ! pw u3 MM5 TOTALRW
!
!              call getmetfield(ident(20),itmp,trw(1,1,k))

	  case (39)

	      call getmetfield(ident(20),itmp,cc3d(1,1,k))

!..2D fields!

	  case (8)

	      call getmetfield(ident(20),itmp,ps(1,1,nr))

	  case (31)

	      call getmetfield(ident(20),itmp,t2_nwp(1,1,nr))

	  case (36)

	      call getmetfield(ident(20),itmp,fh(1,1,nr))

	  case (37)       !ds u7.4vg fl added

              call getmetfield(ident(20),itmp,fl(1,1,nr))

	  case (38)

	      call getmetfield(ident(20),itmp,tau(1,1,nr)) ! ds was fm

	  case (53) !pw u3

              foundustar = .true.
	      !ds apr2005 call getmetfield(ident(20),itmp,ustar_nwp(1,1,nr))
	      call getmetfield(ident(20),itmp,ustar_nwp(1,1))

	  case (103) ! SST

	      call getmetfield(ident(20),itmp,sst(1,1,nr))


	  end select

	enddo

998	continue

!     definition of the half-sigma levels from the full levels.

!pw rv2_1_10: moved into DegGrid
!      sigma_bnd(KMAX_BND) = 1.!
!
!      do k = KMAX_MID,2,-1
!        sigma_bnd(k) = 2.*sigma_mid(k) - sigma_bnd(k+1)
!	enddo
!
!      sigma_bnd(1) = 0.

!              set sdot equal to zero at the boundaries.

!ko    sdot is interpolated to get values at half levels; see metvar.f
!ko	    sdot(:,:,1,nr)=0.
!ko

!pw rv2_1_9: moved into metvar
!      sdot(:,:,KMAX_BND,nr)=0.

   end subroutine infield

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   subroutine metvar(numt)

! This routines postprocess the meteo fields:
! Unit changes, special definitions etc...


	use Par_ml , only : limax,ljmax,li0,li1,lj0,lj1,me		&
			,neighbor,WEST,EAST,SOUTH,NORTH,NOPROC		&
			,MSG_NORTH2,MSG_EAST2,MSG_SOUTH2,MSG_WEST2
        use GridValues_ml , only : xm,xm2,xmd, sigma_bnd,sigma_mid,GRIDWIDTH_M
	use ModelConstants_ml, only : PASCAL, PT, CLOUDTHRES, METSTEP,&
                         V_RAIN  !rv1.2
	use PhysicalConstants_ml, only : KARMAN, KAPPA, R, CP, GRAV      &
			,ROWATER
	!ds apr2005 use Tabulations_ml , only : TPI,PBAS,PINC
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
	if (numt.eq.1) then

          nr = 1

          !-------------------------------------------------------------------
          !ds apr2005  Initialisations:

          call Exner_tab()

          ! Look for processor containing debug coordinates
          debug_iloc    = -999
          debug_jloc    = -999

          do i = 1, limax
          do j = 1, ljmax
               if (MY_DEBUG .and. &
                   i_glob(i) == DEBUG_I .and. j_glob(j) == DEBUG_J ) then
                       debug_proc = .true.
                       debug_iloc    = i
                       debug_jloc    = j
               end if
          end do
          end do
           if( debug_proc ) write(*,*) "DEBUG EXNER me", me, Exner_nd(99500.0)
          !-------------------------------------------------------------------
          
        end if


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

	    !ds apr2005 t2(i,j) = th2m(i,j,nr)  !u7.4vg not -273.15
!su	    th2m(i,j,nr) = th2m(i,j,nr)*(1.e+5/ps(i,j,nr))**(KAPPA)
	   !ds apr2005  th2m(i,j,nr) = th2m(i,j,nr)*exp(-KAPPA*log(ps(i,j,nr)*1.e-5))

            rho_surf(i,j)  = psurf(i,j)/(RGAS_KG * t2_nwp(i,j,nr) ) 

!ds apr2005: For MM5 we get u*, not tau. Since it seems better to
!             interpolate tau than u*  between time-steps we convert

            if ( foundustar) then
               tau(i,j,nr)    = ustar_nwp(i,j)*ustar_nwp(i,j)* rho_surf(i,j) 
            end if


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


!hf I add up in WetDeposition, to have the prec used in the model

	    pr(i,j,:) = prhelp(:)*divt

         ! interpolation of sigma dot for half layers

          if(sdot_at_mid)then !pw rv1_9_24
             do k = KMAX_MID,2,-1

	      sdot(i,j,k,nr) = sdot(i,j,k-1,nr) 		&
		+ (sdot(i,j,k,nr)-sdot(i,j,k-1,nr))		&
            *(sigma_bnd(k)-sigma_mid(k-1)) / (sigma_mid(k)-sigma_mid(k-1))

  	      enddo
            endif
	  
           ! set sdot equal to zero at the top and bottom of atmosphere. 
           !ds pw correction here

            sdot(i,j,KMAX_BND,nr)=0.0
            sdot(i,j,1,nr)=0.0

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

          !ds apr2005 remove:
	  ! x1 = (p1 - PBAS)/PINC
	  ! lx1 = x1
          ! exf1(KMAX_BND) = tpi(lx1) + (x1-lx1)*(tpi(lx1+1) - tpi(lx1))

          exf1(KMAX_BND) = CP * Exner_nd(p1)   !ds apr2005tpi(lx1)                  &

          z_bnd(i,j,KMAX_BND) = 0.0

          do k = KMAX_MID,1,-1

!     eddy diffusivity in the surface-layer follows the formulation used 
!     in the nwp-model which is based on Louis (1979), (see mc7e.f).

!     the shorter loop is the inner loop to save memory. the order 
!     of the do loops will be changed on a vector machine.

!     exner-function of the half-layers

            p1 = sigma_bnd(k)*(ps(i,j,nr) - PT) + PT

             !ds apr2005 remove:
	     !x1 = (p1 - PBAS)/PINC
	     !lx1 = x1
	     !exf1(k) = tpi(lx1) + (x1-lx1)*(tpi(lx1+1) - tpi(lx1))

	      exf1(k) = CP * Exner_nd( p1 ) !ds apr2005

            p2 = sigma_mid(k)*(ps(i,j,nr) - PT) + PT

             !ds apr2005 remove:
	     !x2 = (p2 - PBAS)/PINC
	     !lx2 = x2
	     !exf2(k) = tpi(lx2) + (x2-lx2)*(tpi(lx2+1) - tpi(lx2))

!     exner-function of the full-levels

	      exf2(k) = CP * Exner_nd(p2)  !ds apr2005

!     height of the half-layers ! full layers ???

            z_bnd(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
			(exf1(k+1) - exf1(k)))/GRAV


!     height of the full levels. !half layers??

            z_mid(i,j,k) = z_bnd(i,j,k+1) + (th(i,j,k,nr)*            &
			(exf1(k+1) - exf2(k)))/GRAV

            roa(i,j,k,nr) = CP*((ps(i,j,nr) - PT)*sigma_mid(k) + PT)/      &
			(R*th(i,j,k,nr)*exf2(k))

          enddo  ! k

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
!commented out: pw 6/3-2005: ro not used!
!            ro = ((ps(i,j,nr) - PT)*sigma_bnd(k) + PT)*CP*(exf2(k) -      &
!			exf2(k-1))/(R*exf1(k)*dex12)
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


	    enddo ! k


	  enddo
	enddo

        if( MY_DEBUG .and. debug_proc ) then
           write(*,*) "DEBUG meIJ" , me, limax, ljmax
           do k = 1, KMAX_MID
              write(6,"(a12,2i3,2f12.4)") "DEBUG Z",me, k, &
                  z_bnd(debug_iloc,debug_jloc,k), z_mid(debug_iloc,debug_jloc,k)
           end do
        end if

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
       
       call met_derived !compute derived meteo fields

       call tiphys(numt) 

     end subroutine metvar

  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

     subroutine metint

!     this routine does the forward linear stepping of the meteorological
!     fields read or derived every 3 hours.

	use ModelConstants_ml , only : nmax,nstep

	implicit none

        integer :: i,j
	real :: div,ii
     
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
!          ccc(:,:,:,1) = ccc(:,:,:,1)                           & !ASSYCON
!                        + (ccc(:,:,:,2) - ccc(:,:,:,1))*div       !ASSYCON
	  skh(:,:,:,1) = skh(:,:,:,1) 				&
			+ (skh(:,:,:,2) - skh(:,:,:,1))*div
	  roa(:,:,:,1) = roa(:,:,:,1) 				&
			+ (roa(:,:,:,2) - roa(:,:,:,1))*div

	  psurf(:,:) = ps(:,:,1)  !u7.4vg was psa

	  ps(:,:,1) = ps(:,:,1) 				&
			+ (ps(:,:,2) - ps(:,:,1))*div
	  t2_nwp(:,:,1) = t2_nwp(:,:,1) 				&
			+ (t2_nwp(:,:,2) - t2_nwp(:,:,1))*div
!u7.4vg - note we need pressure first
!ds apr2005 t2(:,:)    =   th2m(:,:,1) * exp(KAPPA*log(psurf(:,:)*1.e-5))

	  fh(:,:,1) = fh(:,:,1) 				&
			+ (fh(:,:,2) - fh(:,:,1))*div
!ds u7.4vg fl added
	  fl(:,:,1) = fl(:,:,1) 				&
			+ (fl(:,:,2) - fl(:,:,1))*div

	  tau(:,:,1) = tau(:,:,1) 				&
			+ (tau(:,:,2) - tau(:,:,1))*div

!SST
	  sst(:,:,1)   = sst(:,:,1) 				&
			+ (sst(:,:,2)   - sst(:,:,1))*div


!ko  precipitation and cloud cover are no longer interpolated

	else

!     assign the the meteorological data at time-level 2 to level 1 for
!     the next 6 hours integration period before leaving the inner loop.

	  u(:,:,:,1) = u(:,:,:,2)
	  v(:,:,:,1) = v(:,:,:,2)
	  sdot(:,:,:,1) = sdot(:,:,:,2)
	  th(:,:,:,1) = th(:,:,:,2)
	  q(:,:,:,1) = q(:,:,:,2)
!          ccc(:,:,:,1) = ccc(:,:,:,2)   !ASSYCON
	  skh(:,:,:,1) = skh(:,:,:,2)
	  roa(:,:,:,1) = roa(:,:,:,2)

!su	don't forget psa !!!!
	  !7.4vg, but put after ps update psa(:,:) = ps(:,:,1)

	  ps(:,:,1) = ps(:,:,2)
	  psurf(:,:) = ps(:,:,1)   ! u7.4vg
	  t2_nwp(:,:,1) = t2_nwp(:,:,2)

!u7.4vg - note we need pressure first
!ds apr2005 t2(:,:)    =   th2m(:,:,1) * exp(KAPPA*log(psurf(:,:)*1.e-5))

	  fh(:,:,1) = fh(:,:,2)
	  tau(:,:,1) = tau(:,:,2)
!ds u7.4vg fl added
	  fl(:,:,1) = fl(:,:,2)

	  sst(:,:,1) = sst(:,:,2)          !SST

	endif

        call met_derived !update derived meteo fields

      end subroutine metint

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   subroutine met_derived

! This routine calculates fields derived from meteofields.
! The interpolation in time is done for the meteofields and the
! fields here are derived from the interpolated fields after 
! each interpolation (i.e. every dt_advec). 
! CPU costly fields (those with special functions like log )
! can be computed in metvar only once every METSTEP and interpolated 
! in metint.

!horizontal wind speed (averaged over the four edges)       
!Note that u and v are wind velocities divided by xm
!At present u_ref is defined at KMAX_MID

     implicit none
     integer ::i,j
     logical :: DEBUG_DERIV = .false.

     do j = 1,ljmax
        do i = 1,limax
           u_ref(i,j)=0.125*(&
                sqrt(0.25*( u(i,j,KMAX_MID,1)*(xm(i,j)+xm(i+1,j))&
                +u(i-1,j,KMAX_MID,1)*(xm(i-1,j)+xm(i,j)) )**2&
                +( v(i,j,KMAX_MID,1)*(xm(i,j)+xm(i,j+1) ))**2)&
                +sqrt(0.25*( u(i,j,KMAX_MID,1)*(xm(i,j)+xm(i+1,j))&
                +u(i-1,j,KMAX_MID,1)*(xm(i-1,j)+xm(i,j)) )**2&
                +( v(i,j-1,KMAX_MID,1)*(xm(i,j-1)+xm(i,j)) )**2)&
                +sqrt(( u(i,j,KMAX_MID,1)*(xm(i,j)+xm(i+1,j)) )**2&
                +0.25*( v(i,j,KMAX_MID,1)*(xm(i,j)+xm(i,j+1))&
                +v(i,j-1,KMAX_MID,1)*(xm(i,j-1)+xm(i,j)) )**2)&
                +sqrt((u(i-1,j,KMAX_MID,1)*(xm(i-1,j)+xm(i,j)))**2&
                +0.25*( v(i,j,KMAX_MID,1)*(xm(i,j)+xm(i,j+1))&
                +v(i,j-1,KMAX_MID,1)*(xm(i,j-1)+xm(i,j)) )**2) )
           
        enddo  
     enddo

   !ds apr2005
   ! Tmp ustar solution. May need re-consideration for MM5 etc., but
   ! basic principal should be that fm is interpolated with time, and
   ! ustar derived from this.

     !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

      forall( i=1:limax, j=1:ljmax ) 
           rho_surf(i,j)  = psurf(i,j)/(RGAS_KG * t2_nwp(i,j,1) ) 
      end forall

      if(.not. foundustar)then
         forall( i=1:limax, j=1:ljmax ) 
            ustar_nwp(i,j)   = sqrt( tau(i,j,1)/rho_surf(i,j) )
         end forall
      endif


      forall( i=1:limax, j=1:ljmax ) 


          ustar_nwp(i,j) = max( ustar_nwp(i,j), 1.0e-5 )

          invL_nwp(i,j)  = KARMAN * GRAV * fh(i,j,1)/ &
             (CP*rho_surf(i,j)* ustar_nwp(i,j)*ustar_nwp(i,j)*ustar_nwp(i,j) * t2_nwp(i,j,1))

      !.. we limit the range of 1/L to prevent numerical and printout problems
      !.. and because I don't trust HIRLAM enough.
      !   This range is very wide anyway.

         invL_nwp(i,j)  = max( -1.0, invL_nwp(i,j) ) !! limit very unstable
         invL_nwp(i,j)  = min(  1.0, invL_nwp(i,j) ) !! limit very stable

      end forall
      if ( DEBUG_DERIV .and. debug_proc ) then
          i = debug_iloc
          j = debug_jloc
          write(*,*) "MET_DERIV DONE ", me, ustar_nwp(i,j), invl_nwp(i,j), rho_surf(i,j)
      end if
     !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa


   end subroutine met_derived

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

!Variable names now more intuitive. This is a list with some new and ols 
!variable names
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
        use PhysicalConstants_ml, only :CP,PI,KARMAN,GRAV,KAPPA,R
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
!c	t2_nwp	: potensial temperature at 2m height, K
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
!c	KAPPA	: r/CP (-)
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
!c	zimz	: ziu - zs_bnd
!c	zis	: height of the stable ABL, m
!c	ziu	: height of the unstable ABL, m
!c	zixx	: Height og ABL (intermediate value), m
!c	zm	: geopotential height of full sigma levels above topography, m
!c	zmhs	: zs_bnd - hs
!c	zs_bnd	: geopotential height of  half sigma levels above topography, m
!c	ztop	: height of the uppermost layer in s-coordinates
!c
!c-------------------------------------------------------------------
!c..the following sketches the sigma-surfaces:
!c
!c
!!                ///////////////////
!c    sigma_bnd(1) = 0 - -sigmas - - - - - sdot(1) = 0, xksig(1)=xksm(1)=0, 
!!                                           pr(1)=0,PT,exns(1), zs_bnd(1)
!c
!c        sigma_mid(1) ---sigmam---------- u, v, th, q, cw, exnm (1)
!c
!c
!c        sigma_bnd(2) - - - - s - - - - - sdot(2), xksig(2), exns(2), pr(2) 
!!                                                 zs_bnd(2), xksm(2)
!c
!c        sigma_mid(2) --------m---------- u, v, th, q, cw, exnm (2)
!c
!c
!c        sigma_bnd(3) - - - - s - - - - - sdot(3), xksig(3), exns(3), pr(3)
!!                                                 zs_bnd(3), xksm(3)
!c
!c        sigma_mid(3) --------m---------- u, v, th, q, cw, exnm (3)
!c
!c
!c        sigma_bnd(4) - - - - s - - - - - sdot(4), xksig(4), exns(4), pr(4)
!!                                                  zs_bnd(4), xksm(4)
!c
!c        sigma_mid(4) --------m---------- u, v, th, q, cw, exnm (4)
!c
!c
!c        sigma_bnd(5) - - - - s - - - - - sdot(5), xksig(5), exns(5), pr(5)
!!                                                  zs_bnd(5), xksm(5)
!c
!!                        :
!!                        :
!c
!c  sigma_bnd(KMAX_BND-1) - - - - s - - - - - sdot(KMAX_BND-1), xksig(KMAX_MID), 
!!                                    exns(KMAX_BND-1),zs_bnd(KMAX_BND-1), 
!!                                    pr(KMAX_BND-1),xksm(KMAX_MID)
!c
!c    sigma_mid(KMAX_MID) --------m---------- u, v, th, q, cw, exnm (KMAX_MID); 
!!                                    this level is assumed to be
!!                                    the top of Prandtl-layer (LAM50E)
!c
!c sigma_bnd(KMAX_BND) = 1- - - - s - - - - - sdot(KMAX_BND) = 0, ps, t2_nwp, fh, 
!!                ///////////////////        fm, mslp, xksig(KMAX_MID)=0, 
!!                                           exns(KMAX_BND), zs_bnd(KMAX_BND), 
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
 logical, parameter :: PIELKE_KZ = .true.    ! Default
 logical, parameter :: TKE_DIFF = .false.  !!! CODE NEEDS TESTING/TIDY UP STILL!!!!


!definer alle dimensjoner med MAXLIMAX,MAXLJMAX
 real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID)::exnm
 real, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND)::exns,zs_bnd
 real, dimension(MAXLIMAX,KMAX_MID)::zm,dthdz,deltaz,thc
 real, dimension(MAXLIMAX,KMAX_BND)::risig,xksm,pz
 real, dimension(MAXLIMAX)::zis,delq,thsrf,trc,pidth,dpidth,xkhs,xkdz,xkzi,&
                            hs,xkh100
 !ds apr2005 real, dimension(MAXLIMAX,MAXLJMAX)::ziu,help,a,zixx,roas,uabs,vdfac
 real, dimension(MAXLIMAX,MAXLJMAX)::ziu,help,a,zixx,uabs,vdfac
 real ::lim,xdthdz,zmmin,zimin,zlimax,kzmin,kzmax,sm,pref,xtime,umax,eps,ric,&
        ric0,dthdzm,dthc,xdth,xfrco,exfrco,hsl,dtz,p,dvdz,xl2,uvhs,zimhs,&
        zimz,zmhs,ux0,fac,fac2,dex12,ro


!hf Hilde&ANton
      real hsurfl
!hf new
      integer i,j,k,km,km1,km2,kabl,iip,jjp,numt,kp
!ds Kz-tests
      real, parameter :: KZ_MINIMUM = 0.001   ! m2/s
      real, parameter :: KZ_MAXIMUM = 1.0e3 ! m2/s - as old kzmax
      real, parameter :: KZ_SBL_LIMIT = 0.1 ! m2/s - Defines stable BL height

      integer nh1(MAXLIMAX),nh2(MAXLIMAX),nr
      real :: h100 ! Top of lowest layer - replaces 100.0 

!ds Check:
      if(  KZ_SBL_LIMIT < 1.01*KZ_MINIMUM ) &
             call gc_abort(me,NPROC,"SBL limit too low!!")
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

      eps = 0.01
      dtz = 3600.
      sm = 0.04
!c
!c
!c..preset=zero:
   xksm(:,:)   = 0
   risig(:,:)  = 0.
   xksig(:,:,:)= 0.


!c..................................
!c..exner-function in the full sigma-levels..
!c
  do  k=1,KMAX_MID
    do j=1,ljmax
     do i=1,limax
!c
!c..pressure (pa)
        p = PT + sigma_mid(k)*(ps(i,j,nr) - PT)
!c..exner (j/k kg)
        !ds apr2005   exnm(i,j,k)=CP*(p/pref)**KAPPA
        exnm(i,j,k)= CP * Exner_nd(p)
      end do
    end do
  end do

!c.........................................
!c..procedure to arrive at mixing height..:
!c
!c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!c
!c     Start j-slice here.
!c
!ds-Kz      lim = 0.1
!c..hj..test lim=1.


   do 40 j=1,ljmax

!c..exner in half-sigma levels:
     do  k=1,KMAX_BND
       do i=1,limax
          p = PT + sigma_bnd(k)*(ps(i,j,nr) - PT)
          pz(i,k) = p
          !ds apr2005 exns(i,j,k)=CP*(p/pref)**KAPPA  
          exns(i,j,k)= CP * Exner_nd(p)
       end do
     end do
!c
!c
!c.. exns(KMAX_BND), th(KMAX_BND) and height of sigmas:
     do  i=1,limax
       zs_bnd(i,j,KMAX_BND)=0.
     end do
!c
!c     Height of the half levels
!c
     do  k=KMAX_BND-1,1,-1 
       do i=1,limax
            zs_bnd(i,j,k)=zs_bnd(i,j,k+1)+th(i,j,k,nr)*&
                 (exns(i,j,k+1)-exns(i,j,k))/GRAV

       end do
     end do
!c
!c..height of sigma:
     do  k=1,KMAX_MID
       do i=1,limax
          zm(i,k) = ((exnm(i,j,k)-exns(i,j,k))*zs_bnd(i,j,k+1)&
                  + (exns(i,j,k+1)-exnm(i,j,k))*zs_bnd(i,j,k))&
                  / (exns(i,j,k+1)-exns(i,j,k))    
       end do
     end do


!c----------------------------------------------------------------------
!c...........................................
!c..the following variables in sigmas-levels:
!c
     do  k=2,KMAX_MID
       km=k-1
       do i=1,limax
!c
!c.........................
!c..wind sheare
!c
!HF Slightly different formulation of dvdz than in metvar
        dvdz = ( (u(i,j,km,nr)-u(i,j,k,nr))**2 &
              + (v(i,j,km,nr)-v(i,j,k,nr))**2 + eps)
!c
            risig(i,k)=(2.*GRAV/(th(i,j,km,nr)+th(i,j,k,nr)))*&
               (th(i,j,km,nr)-th(i,j,k,nr))*(zm(i,km)-zm(i,k))&
               /dvdz
!c........................
!c..mixing length squared:
!c
            xl2=(KARMAN*amin1(zs_bnd(i,j,k),zmmin))**2

!c
!c..............................
!c..critical richardsons number:
!c
            ric0=0.115*((zm(i,km)-zm(i,k))*100.)**0.175
            ric=amax1(0.25,ric0)
!c

!c
            dvdz = sqrt(dvdz)/(zm(i,km)-zm(i,k))

!c..................................................................
!c..exchange coefficient (Pielke,...)
 !ds alternative
       if ( PIELKE_KZ  ) then
          if (risig(i,k) > ric ) then
              xksig(i,j,k) = KZ_MINIMUM
          else
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
      end do
    end do
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
  do k=KMAX_MID,2,-1
     do i=1,limax

	!ds-Kz if(xksm(i,k).ge.lim .and. nh2(i).eq.1) then
	if(xksm(i,k) >= KZ_SBL_LIMIT .and. nh2(i) == 1) then
           nh1(i)=k   ! Still unstable
        else
           nh2(i)=0   ! Now stable
        endif
!c
     end do
  end do
!c
   do i=1,limax
!c
        k=nh1(i)
!c
	if(zs_bnd(i,j,nh1(i)).ge.zimin) then

           if( abs(xksm(i,k)-xksm(i,k-1)) .gt. eps) then 

              zis(i)=((xksm(i,k)-KZ_SBL_LIMIT )*zs_bnd(i,j,k-1) &
                  + (KZ_SBL_LIMIT -xksm(i,k-1))*zs_bnd(i,j,k))&
                  /(xksm(i,k)-xksm(i,k-1))
           else

              zis(i)=zimin
           endif

	endif
!c
    end do
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
      do k=KMAX_MID,kabl,-1
	 xdth = thc(i,kabl)-thc(i,k)
	 dpidth(i) = exnm(i,j,k)*xdth*(pz(i,k+1)-pz(i,k))/GRAV
	 pidth(i) = pidth(i) + dpidth(i)
      end do
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


    do i=1,limax

      zixx(i,j)=amax1(ziu(i,j),zis(i))
      zixx(i,j)=amin1(zlimax,zixx(i,j))
    end do

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




!----------------------------------------------------------!
       if( TKE_DIFF ) then
          call tkediff (nr)                            ! guta
       else 
          call O_Brian(nr, KZ_MINIMUM, KZ_MAXIMUM, zimin, zs_bnd, ziu  & 
                     , exns, exnm, zixx )
       end if
!----------------------------------------------------------!       

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

    thick=2  !we fetch 2 neighbors at once, so that we don't need to call
             ! readneighbours twice
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
                    +0.5*s*(1.-s)*(h1(i+1,j)+h1(i-1,j)+h1(i,j+1)+h1(i,j-1))  &
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

!************************************************************************!
     subroutine tkediff (nr)                                             !
!************************************************************************!
!                                                                        !   
!    This routine computes vertical eddy diffusivities as a function     !
!    altitude, height of PBL, and a velocity scale, square root of       !
!    turbulent kinetic energy (TKE). This is a non-local scheme.         !
!    The TKE at the surface is diagnosed using scales for horizontaland  !
!    vertical velocities (ustar and wstar) in the surface layer          !
!    (Alapaty 2004; Holstag et al. 1990 and Mihailovic et al. 2004)      !
!    PBL ht is calculated using the EMEP formulation                     !
!                                                                        !
!    Written by DT Mihailovic (October 2004)                             !
!    EMEP polishing and comments: JE Jonson and P Wind                   !
!************************************************************************!

      implicit none
      integer i, j, k, l, kcbl

!     Local constants
      real   , parameter :: SZKM=1600.     &    ! Constant (Blackadar, 1976)
             ,              CKZ=0.001      &    ! Constant (Zhang and Athens, 1982)
             ,              REFPR=1.0E+05  &    ! Referent pressure
             ,              KZ0LT=1.0E-04  &    ! Constant (Alapaty et al., 1997)
             ,              RIC=0.10       &    ! Critical Richardson number (Holstlag et al., 1993)
             ,              ROVG=R/GRAV         ! Used in Calculation of R-number
      integer, parameter :: KLM =KMAX_MID-1   

!     INPUT      
      integer                              nr              ! Number of meteorological stored
                                                           ! in arrays (1 or 2)
     
!     OUTPUT
!     skh(i,j,k,nr) array    
!     Values of the Kz coefficients (eddyz (i,j,k)) are transformed nto sigma system and
!     then they stored in this array which is later used in ADVECTION module     


!     Local arrays 
   
      integer, dimension(MAXLIMAX,MAXLJMAX)      ::iblht      ! Level of the PBL top
      real, dimension(MAXLIMAX,MAXLJMAX,KMAX_BND):: eddyz     ! Eddy coefficients (m2/s)
      real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID):: &
                    t_virt    &! Potential temperature (K)
                   ,e         &! Kinetic energy with respect to height (m2/s2)
                   ,dzq       &! Thickness of sigma interface layers (m)
                   ,u_mid     &! Wind speed in x-direction (m/s)  
                   ,v_mid      ! Wind speed in y-direction (m/s)       
      real, dimension(MAXLIMAX,MAXLJMAX,KLM):: &
                   dza       ! Thickness of half sigma layers (m)

      real, dimension(MAXLIMAX,MAXLJMAX):: &
                   pblht ,            &! PBL (Holstag, 1990) (m)    
                   h_flux,            &! Sensible heat flux  (W/m2)
                   ust_r ,            &! Friction velocity (m/s) 
                   !ds apr2005 ro_sur,            &! Air density (kg/m3)
                   mol   ,            &! Monin-obukhov length (m)
                   wstar               ! Convective velocity (m/s)
                                                       
      real, dimension(KMAX_BND) :: rib         ! Bulk Richardson number 
      real, dimension(KMAX_MID) :: &
                   rich,           &! Richardson number 
                   psi_zi           ! Used in the vertical integration
      real, dimension (10) :: psi_z, zovh  ! Used for calculating the TKE

!     Local variables
      real dtmp, tog, wssq1, wssq2, wssq, tconv, wss, wst, PSI_TKE,    &
           dusq, dvsq, ri, ss, dthdz, busfc, zvh,                      &
           part1, part2, fract1, fract2, apbl, fac0, fac02, kz0,       &
           cell, dum1, rpsb, press, teta_h, u_s, goth, pressure

!     Functions for averaging the vertical turbulent kinetic energy 
!      (Alapaty, 2003)
      data psi_z /0.00,2.00,1.85,1.51,1.48,1.52,1.43,1.10,1.20,0.25/
      data zovh  /0.00,0.05,0.10,0.20,0.40,0.60,0.80,1.00,1.10,1.20/

!     Store the NMW meteorology and variables derived from its
   
!     Change the sign
      h_flux(1:limax,1:ljmax)=-fh(1:limax,1:ljmax,nr)

!     Avoid devision by zero later in the code
     
      where (ABS(h_flux(1:limax,1:ljmax))<0.0001) h_flux(1:limax,1:ljmax)=0.0001 

!     Check PBL height   !ds - strange tests! Negative pzpbl check? From 1 to 100m
                         !   - odd!
      do i=1,limax
       do j=1,ljmax
            if(ABS(pzpbl(i,j)) < 1.) then
              pzpbl(i,j)=100. 
            endif
        enddo 
       enddo

!     Calculate velocity components in the (h) poits (Arakawa notation)
      do k=1,KMAX_MID
       do i=1,limax
        do j=1,ljmax
!          u_mid(i,j,k)=0.5*(u(i-1,j  ,k,nr)+u(i,j,k,nr))
!          v_mid(i,j,k)=0.5*(v(i  ,j-1,k,nr)+v(i,j,k,nr))
 
          u_mid(i,j,k)=u(i,j  ,k,nr)
          v_mid(i,j,k)=v(i  ,j,k,nr)
               
       enddo
       enddo
      enddo

!     Avoid small values
      where (ABS(u_mid(1:limax,1:ljmax,1:KMAX_MID))<0.001) &
                                u_mid(1:limax,1:ljmax,1:KMAX_MID)=0.001
      where (ABS(v_mid(1:limax,1:ljmax,1:KMAX_MID))<0.001) &
                                v_mid(1:limax,1:ljmax,1:KMAX_MID)=0.001

!     Initialize eddy difussivity arrays
      eddyz(1:limax,1:ljmax,1:KMAX_MID)=0.

!     Calculate tickness of the layers
      dzq(1:limax,1:ljmax,1:KMAX_MID)= &
            z_bnd(1:limax,1:ljmax,1:KMAX_MID)-  &
            z_bnd(1:limax,1:ljmax,2:KMAX_BND)     ! Full levels

      dza(1:limax,1:ljmax,1:KLM)= &
            z_mid(1:limax,1:ljmax,1:KLM)-          &
            z_mid(1:limax,1:ljmax,2:KMAX_MID)          ! Half-sigma lavels

!     Calculate virtual temperature

      t_virt(1:limax,1:ljmax,1:KMAX_MID)=&
             th(1:limax,1:ljmax,1:KMAX_MID,nr)* &
                    (1.0+0.622*q(1:limax,1:ljmax,1:KMAX_MID,nr))


!     Calculate Monin-Obuhkov length   (Garratt, 1994)

      do i=1,limax
       do j=1,ljmax
         !ds apr2005 BUG anyway? for ro_sur ?
         !ds apr2005 ro_sur(i,j)=ps(i,j,nr)/(KAPPA*t2_nwp(i,j,nr)*                     &
                     !ds apr2005 CP*(ps(i,j,nr)/REFPR)**KAPPA)      ! Surface density
           !ds apr2005 if(foundustar) then
           !ds apr2005   ust_r(i,j)=ustar(i,j,1)
           !ds apr2005 else
           !ds apr2005   ust_r(i,j)=SQRT(fm(i,j,nr)/ro_sur(i,j))
           !ds apr2005   ust_r(i,j)=AMAX1(ust_r(i,j),0.00001)
           !ds apr2005 endif
         !ds apr2005 u_s=ust_r(i,j)
         u_s=ustar_nwp(i,j)
         mol(i,j)=-(ps(i,j,nr)*u_s*u_s*u_s)/                        &
                   (KARMAN*GRAV*h_flux(i,j)*KAPPA)
       enddo
      enddo

!     Calculate the convective velocity (wstar)
      do i=1,limax
       do j=1,ljmax
        !ds apr2005 wstar(i,j)=GRAV*h_flux(i,j)*pzpbl(i,j)/ro_sur(i,j)/CP/th(i,j,KMAX_MID,nr)
        wstar(i,j)=GRAV*h_flux(i,j)*pzpbl(i,j)/rho_surf(i,j)/CP/th(i,j,KMAX_MID,nr)
           if(wstar(i,j) < 0.) then
            wstar(i,j)=-ABS(wstar(i,j))**(0.3333)
           else
            wstar(i,j)=(wstar(i,j))**(0.3333)
           endif
       enddo
      enddo

!                            ------------------------------------------>
!     Start with a long loop ------------------------------------------>
!                            ------------------------------------------>
      DO  i=1,limax
      DO  j=1,ljmax

       rib(1:KMAX_MID) = 0.0                        ! Initialize bulk Richardson number

       part1=ust_r(i,j)*ust_r(i,j)*ust_r(i,j)
       wst=AMAX1(wstar(i,j),1.0E-20)
       part2=0.6*wst*wst*wst
       wss=AMAX1(1.0E-4,(part1+part2))
       wss=EXP(0.333333*ALOG(wss))

       if (h_flux(i,j) < 0.0) then
        tconv=0.0                                   ! Holstlag et al. (1990)
       else
        !ds apr2005 tconv=8.5*h_flux(i,j)/ro_sur(i,j)/CP/wss    !Conversion to kinematic flux
        tconv=8.5*h_flux(i,j)/rho_surf(i,j)/CP/wss    !Conversion to kinematic flux
      endif

      do k=KMAX_MID,1,-1
        dtmp=t_virt(i,j,k)-t_virt(i,j,KMAX_MID)-tconv
        tog=0.5*(t_virt(i,j,k)+t_virt(i,j,KMAX_MID))/GRAV
        wssq1=u_mid(i,j,k)*u_mid(i,j,k)
        wssq2=v_mid(i,j,k)*v_mid(i,j,k)
        wssq=wssq1+wssq2
        wssq=AMAX1(wssq,1.0E-4)
        rib(k)=z_mid(i,j,k)*dtmp/(tog*wssq)
        if(rib(k).ge.RIC) go to 9001
      enddo
9001  continue

!     Calculate PBL height according to Holtslag et al. (1993)
      pblht(i,j)=0.
          if(k.ne.KMAX_MID) then
            fract1=(RIC-rib(k+1))/(rib(k)-rib(k+1))
            fract2=1.-fract1
            apbl=z_mid(i,j,k)*fract1
            pblht(i,j)=apbl+z_mid(i,j,k+1)*fract2
                if(pblht(i,j) > z_bnd(i,j,k+1)) then
                   kcbl=k
                else
                   kcbl=k+1
                endif
           endif
      iblht(i,j)=kcbl

      if(pblht(i,j)<z_bnd(i,j,KMAX_MID)) then
         pblht(i,j)=z_bnd(i,j,KMAX_MID)
         iblht(i,j)=KMAX_MID
      endif


      if(pblht(i,j).le.100.) then              !Minimum of PBL height
         pblht(i,j)=100.
      endif
            
!     Find the critical Richardson number (Shir and Borestein, 1976)
      do k=2,iblht(i,j)-1
       rich(k)=0.257*dza(i,j,k)**0.175
      enddo

!     Free troposphere and cloudy case Kz values estimation
      do k=2,iblht(i,j)-1
        dusq=(u_mid(i,j,k-1)-u_mid(i,j,k))*(u_mid(i,j,k-1)-u_mid(i,j,k))
        dvsq=(v_mid(i,j,k-1)-v_mid(i,j,k))*(v_mid(i,j,k-1)-v_mid(i,j,k))
        ss=(dusq+dvsq)/(dza(i,j,k-1)*dza(i,j,k-1))+1.E-9
        goth=2.*GRAV/(t_virt(i,j,k-1)+t_virt(i,j,k))
        dthdz=(t_virt(i,j,k-1)-t_virt(i,j,k))/dza(i,j,k-1)
        ri=goth*dthdz/ss

!     (Duran and Clemp, 1982)

      kz0=CKZ*dzq(i,j,k)
       if (ri-rich(k) > 0.) then
         eddyz(i,j,k)=kz0
       else
         eddyz(i,j,k)=kz0+SZKM*SQRT(ss)*(rich(k)-ri)/rich(k)
       endif
          eddyz(i,j,k)=AMIN1(eddyz(i,j,k),100.)
      enddo

!     Eddy diffusivity coefficients for all regimes in the mixed layer

      do  k=iblht(i,j),KMAX_MID
          if (mol(i,j) < 0.0) then                 !Unstable conditions
             ri=(1.0-15.*z_mid(i,j,k)/mol(i,j))**(-0.25)
             ri=ri/KARMAN/z_mid(i,j,k)
             ri=ri*AMAX1(0.0,pblht(i,j)-z_mid(i,j,k))
             dthdz=ri*ust_r(i,j)**3.
             goth=AMAX1(wstar(i,j),0.0)
             dusq=0.4*goth**3.
             ri=(dthdz+dusq)**(2./3.)
             e(i,j,k)=0.5*ri*(2.6)**(2./3.)        !Moeng and Sullivan (1994)
          else
             ri=z_bnd(i,j,k)/pblht(i,j)               !Stable
             ri=z_mid(i,j,k)/pblht(i,j)               !New
             ri=(1.0-ri)
             ri=AMAX1(0.0,ri)
             ri=(1.0-ri)**1.75
             e(i,j,k)=6.*ust_r(i,j)*ust_r(i,j)*ri  !Lenshow(1988)
          endif

!     Calculate Ksi function using interpolation in the vertical
!     Alapaty (2001, 2003)

      zvh=z_mid(i,j,k)/pblht(i,j)
          do l=1,9
             if (zvh > zovh(l).and.                                    &
                 zvh < zovh(l+1)) then
               psi_zi(k)=(psi_z(l+1)-psi_z(l))/(zovh(l+1)-zovh(l))
               psi_zi(k)=psi_zi(k)*(zvh-zovh(l))
               psi_zi(k)=psi_zi(k)+psi_z(l)
               psi_zi(k)=psi_zi(k)/2.0               !Normalized the value
             endif
          enddo
       enddo

!      Calculate integral for Ksi
       psi_tke=0.
        do k=KMAX_MID,iblht(i,j),-1
          psi_tke=psi_tke+psi_zi(k)*dzq(i,j,k)*sqrt(e(i,j,k))
        enddo
       psi_tke=psi_tke/pblht(i,j)

       do k=iblht(i,j),KMAX_MID                           !Calculate coefficients
        goth=psi_tke
        goth=goth*KARMAN*z_mid(i,j,k)
        dthdz=z_mid(i,j,k)/pblht(i,j)
        dthdz=1.0-dthdz
        dthdz=AMAX1(1.0E-2,dthdz)
          if(mol(i,j) > 0.0) then                         !Stable
             goth=sqrt(e(i,j,iblht(i,j)))                 ! Mihailovic (2004)
             goth=goth*KARMAN*z_mid(i,j,k)                 ! -----------------
             dthdz=z_mid(i,j,k)/pzpbl(i,j)                 ! -----------------
             dthdz=1.0-dthdz
             dthdz=AMAX1(1.0E-2,dthdz)
             busfc=0.74+4.7*z_mid(i,j,KMAX_MID)/mol(i,j)     
             busfc=AMAX1(busfc,1.0)
!             dthdz=dthdz*dthdz
              dthdz=dthdz**1.50                                  !test (2004)
              eddyz(i,j,k)=goth*dthdz/busfc
          else
             dthdz=dthdz*dthdz
             busfc=1.0
             eddyz(i,j,k)=goth*dthdz/busfc
          endif
       enddo

!      Checking procedure
       do k=2,iblht(i,j)-1
         if(eddyz(i,j,k).le.0.0) THEN
          eddyz(i,j,k)= KZ0LT
         endif
       enddo

!      Avoid phisically unrealistic values
       do k=2,KMAX_MID
        IF(eddyz(i,j,k).le.0.1) then
           eddyz(i,j,k)=0.1
        endif
       enddo

!     To avoid loss of mass/energy through top of the model
!     put eddyz (I,J,K) to zero at the last  level from top
      eddyz(i,j,KMAX_BND)=0.0

!     Calculate eddy coefficients at the interfaces
      do k=2,KMAX_MID
       eddyz(i,j,k)=0.5*(eddyz(i,j,k-1)+eddyz(i,j,k)) !!

!              if(i.eq.10.and.j.eq.10.) then
!             if (abs(u(i,j  ,k,nr)-u_mid(i,j,k)).gt.5.) then
!
!       print *,"NEW ",i,j,u(i,j  ,KMAX_MID,nr),u_mid(i,j,KMAX_MID)
!              endif
      enddo
 
!     Transform values of the eddy coeficients into the the sigma coordinate
    
      do k=2,KMAX_MID
       eddyz(i,j,k)=eddyz(i,j,k)*((sigma_mid(k)-sigma_mid(    k-1))/   &
                                 (    z_mid(i,j,k)-z_mid(i,j,k-1)))**2.

      enddo

      ENDDO                  !---------------------------------------->
      ENDDO                  !---------------------------------------->
                             !---------------------------------------->

!     Store diffusivity coefficients into skh(i,j,k,nr) array
      do k=2,KMAX_MID
       do i=1,limax
        do j=1,ljmax
            skh(i,j,k,nr)=eddyz(i,j,k)
         enddo
        enddo
       enddo

   ! For plotting set pblht  =  pzpbl
     
   pzpbl(:,:) = pblht(:,:)
                     
       RETURN 
       end subroutine tkediff
!---------------------------------------------------------------


!************************************************************************!
     subroutine O_Brian(nr, KZ_MINIMUM, KZ_MAXIMUM, zimin, zs_bnd, ziu  &
                      , exns, exnm, zixx )                            !
!************************************************************************!

!c......................................................
!c..exchange coefficients for convective boundary layer:
!c..o'brien's profile formula:
!c..and the air density at ground level:
!c
!c..constants for free-convection limit:
!c

 integer, intent(in) :: nr

 real, intent(in) :: zimin        &
                    ,KZ_MINIMUM   &
                    ,KZ_MAXIMUM

 real,intent(in), dimension(MAXLIMAX,MAXLJMAX,KMAX_BND) :: zs_bnd &
                                                          ,exns   
 real,intent(in), dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: exnm 

 real,intent(in), dimension(MAXLIMAX,MAXLJMAX) :: ziu    &
                                                 ,zixx


  real :: h100   & ! Top of lowest layer - replaces 100.0 
         ,xfrco  &
         ,exfrco &
         ,sm     &
         ,ux0    &   ! local ustar
         ,ux3    &   ! ustar**3, ds apr2005
         ,hsl    &
         ,hsurfl &
         ,zimhs  &
         ,zimz   &
         ,zmhs   &
         ,fac    &
         ,fac2   &
         ,dex12  &
         ,ro


  integer :: i,j,k

! local arrays:
  real, dimension(MAXLIMAX)::   xkh100  &
                                ,xkhs    &
                                ,xkdz    &
                                ,xkzi    &
                                ,hs      


  real, dimension(MAXLIMAX,MAXLJMAX) :: & !ds apr2005 roas   &
                                       help

      sm = 0.04


      xfrco=0.5*(sqrt(6859.)-1)
      exfrco=1./3.

      do 70 j=1,ljmax

!c..exchange parameter and its vertical derivative at z = hs

      do 60 i=1,limax

         xkh100(i)=0.  !Hilde&Anton
         xkhs(i)=0.                                            
         xkdz(i)=0.
         xkzi(i)=0.
         h100 = zs_bnd(i,j,KMAX_MID)
!c
!c
!c...................................................................
!c..air density at ground level is always calculated diagnostically:
!c
         !ds apr2005 roas(i,j)=ps(i,j,nr)/(KAPPA*th2m(i,j,nr)*exns(i,j,KMAX_BND))

         ux0 = ustar_nwp(i,j)   !ds apr2005
         ux3 = ux0*ux0*ux0      !ds apr2005


         if(ziu(i,j) >= zimin) then
!c
!c..........................
!c..unstable surface-layer.:
!co
!c..height of surface layer
            hs(i)=sm*ziu(i,j)
!c..u*
!c
            !ds apr2005 if(foundustar)then                      
            !ds apr2005    ux0 = ustar(i,j,1)
            !ds apr2005 else
            !ds apr2005    ux0 = sqrt(fm(i,j,nr)/roas(i,j))
            !ds apr2005 endif
            !ds apr2005 ux0=amax1(ux0,0.00001)

!c..hsl=hs/l where l is the monin-obhukov length
            hsl=KARMAN*GRAV*hs(i)*fh(i,j,nr)*KAPPA &
             /(ps(i,j,nr)*ux3)


           !ds rv1_7_2 changes: use simple Garratt \Phi function
           !   instead of "older" Businge and Iversen/Nordeng stuff:

               xkhs(i)=ux0*KARMAN*hs(i)*sqrt(1.0-16.0*hsl)  ! /Pr=1.00   
               xkdz(i)=xkhs(i)*(1.-0.5*16.0*hsl/(1.0-16.0*hsl))/hs(i)        

!Hilde&Anton
!pw & hf            hsurfl=KARMAN*GRAV*100.*amax1(0.001,fh(i,j,nr))*KAPPA&
!pw & hf                 &             /(ps(i,j,nr)*ux0*ux0*ux0)
            !ds hsurfl=KARMAN*GRAV*100.*fh(i,j,nr)*KAPPA&
            hsurfl=KARMAN*GRAV*h100*fh(i,j,nr)*KAPPA&
                 &             /(ps(i,j,nr)*ux3)
                 !ds apr2005 &             /(ps(i,j,nr)*ux0*ux0*ux0)

               !ds xkh100(i)=ux0*KARMAN*100.*sqrt(1.-16.*hsurfl) !/Pr=1.00
               xkh100(i)=ux0*KARMAN*h100*sqrt(1.-16.*hsurfl)

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
            !ds apr2005 ux0=sqrt(fm(i,j,nr)/roas(i,j))
            !ds apr2005 ux0=amax1(ux0,0.00001)
!c
!c..hsl=hs/l where l is the monin-obhukov length
            hsl=KARMAN*GRAV*hs(i)*amax1(0.001,fh(i,j,nr))*KAPPA&
              /(ps(i,j,nr)*ux3)
             !ds apr2005 /(ps(i,j,nr)*ux0*ux0*ux0)


            !xksig(i,j,KMAX_MID)=ux0*KARMAN*hs(i)/(0.74+4.7*hsl)   
            xksig(i,j,KMAX_MID)=ux0*KARMAN*hs(i)/(1.00+5.0*hsl)   

   !ds TEST:
    ! Should apply PhiM for all layers below:
    !       if (hs(i) > zs_bnd(i,j,KMAX_MID) ) then
    !         do k = KMAX_MID, nh1(i), -1
    !            hsl=KARMAN*GRAV*zs_bnd(i,j,k)*amax1(0.001,fh(i,j,nr))*KAPPA&
    !             /(ps(i,j,nr)*ux0*ux0*ux0)
    !            xksig(i,j,KMAX_MID)=ux0*KARMAN*hs(i)/(1.00+5.0*hsl)   
    !         end do
    !       end if
           

         endif
!hf Hilde&Anton
            hsurfl=KARMAN*GRAV*100.*amax1(0.001,fh(i,j,nr))*KAPPA&
                 &             /(ps(i,j,nr)*ux3)   !ds apr2005 - ux3
            !Kz_min(i,j)=1.35*ux0*KARMAN*100./(0.74+4.7*hsurfl)
            !ds Kz_min(i,j)=ux0*KARMAN*100./(1.00+5.0*hsurfl)
            Kz_min(i,j)=ux0*KARMAN*h100/(1.00+5.0*hsurfl)
!c
!c...............................................................

 60   continue
!c
!c
!c..exchange parameter at z = ziu
!c
      do 65 k=1,KMAX_MID
      do 65 i=1,limax

         if(ziu(i,j).gt.zimin .and. zs_bnd(i,j,k).ge.ziu(i,j)) then
            xkzi(i)=xksig(i,j,k)
         elseif (ziu(i,j).gt.zimin) then
!c
!c.....................................................   
!c..the obrien-profile for z<ziu                      . 
!c.....................................................                  
!c
            if(zs_bnd(i,j,k).le.hs(i)) then   
               xksig(i,j,k)=zs_bnd(i,j,k)*xkhs(i)/hs(i)          
            else                                                      
               zimhs=ziu(i,j)-hs(i)   
               zimz=ziu(i,j)-zs_bnd(i,j,k)                     
               zmhs=zs_bnd(i,j,k)-hs(i)                  
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
!ds DaveTest         if ( (pzpbl(i,j)>z_mid(i,j,k)) )then
!ds Restored:
               if ( (pzpbl(i,j)>z_mid(i,j,k)) )then
                  xksig(i,j,k)=max(xksig(i,j,k),Kz_min(i,j))
               endif 
               help(i,j) = xksig(i,j,k)
            enddo
         enddo

       !ds-Kz call smoosp(help,kzmin,kzmax)
       call smoosp(help,KZ_MINIMUM ,KZ_MAXIMUM )

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
!               if(zixx(i,j).ge.zs_bnd(i,j,k)) then
!
!                  dz = zs_bnd(i,j,k)-zs_bnd(i,j,k+1)
!
!                  u2 = ( u(i,j,k,nr)**2 + v(i,j,k,nr)**2 )
!
!                  uabs(i,j)=uabs(i,j)+sqrt(u2)*dz
!
!               elseif (zs_bnd(i,j,k).gt.zixx(i,j) &
!                     .and. zs_bnd(i,j,k+1).lt.zixx(i,j)) then
!
!                    dz = zixx(i,j)-zs_bnd(i,j,k+1)
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


!c++++++++++++++++
                     
       RETURN 
       end subroutine O_Brian


     subroutine Getmeteofield(meteoname,namefield,nrec,&
                              ndim,period_of_validity,field)
!
! Read the meteofields and distribute to nodes
!

     use netcdf
      use Par_ml,                only : me,NPROC,MSG_READ4,ISMBEG,JSMBEG&
                                       ,GIMAX,GJMAX,MAXLIMAX,MAXLJMAX
     use ModelConstants_ml,     only : KMAX_MID

     implicit none

     real, dimension(*),intent(out) :: field !dimensions: (MAXLIMAX,MAXLJMAX)
                                             ! or (MAXLIMAX,MAXLJMAX,KMAX)
     character (len = *),intent(in) ::meteoname,namefield
     character (len = *),intent(out) ::period_of_validity
     integer,intent(in) :: nrec,ndim

     integer*2 :: var_local(MAXLIMAX,MAXLJMAX,KMAX_MID)
     integer*2, allocatable ::var_global(:,:,:) !faster if defined with fixed dimensions for all nodes?
     real :: scalefactors(2)
     integer :: gc_info,KMAX,ijk,i,k,j,nfetch

     period_of_validity=''

     if(ndim==3)KMAX=KMAX_MID
     if(ndim==2)KMAX=1
     if(me==0)then
        allocate(var_global(GIMAX,GJMAX,KMAX))
        nfetch=1
        call GetCDF_short(namefield,meteoname,var_global,GIMAX,ISMBEG,GJMAX, &
        JSMBEG,KMAX,nrec,nfetch,scalefactors,period_of_validity)
     else
        allocate(var_global(1,1,1)) !just to have the array defined
     endif

!note: var_global is defined only for me=0
     call global2local_short(var_global,var_local,MSG_READ4,GIMAX,GJMAX,&
                             KMAX,1,1)
     call gc_rbcast(200,2,0,NPROC,gc_info,scalefactors)
     call gc_bbcast(201,20,0,NPROC,gc_info,period_of_validity)
 
     deallocate(var_global)
 

     ijk=0
     do k=1,KMAX ! KMAX is =1 for 2D arrays
        do j=1,MAXLJMAX
           do i=1,MAXLIMAX
              ijk=ijk+1
              field(ijk)=var_local(i,j,k)*scalefactors(1)+scalefactors(2)
           enddo
        enddo
     enddo
     
     return
   end subroutine Getmeteofield

subroutine GetCDF_short(varname,fileName,var,GIMAX,ISMBEG,GJMAX,JSMBEG &
     ,KMAX,nstart,nfetch,scalefactors,period_of_validity)
  !
  ! open and reads CDF file
  !
  ! The nf90 are functions which return 0 if no error occur.
  ! check is only a subroutine which check wether the function returns zero
  !
  !
  use netcdf
  use Par_ml,                only : me,NPROC
  implicit none

  character (len=*),intent(in) :: fileName 

  character (len = *),intent(in) ::varname
  character (len = *),intent(out) ::period_of_validity
  real,intent(out) :: scalefactors(2)
  integer, intent(in) :: nstart,GIMAX,ISMBEG,GJMAX,JSMBEG,KMAX
  integer, intent(inout) ::  nfetch
  integer*2, dimension(GIMAX*GJMAX*KMAX*NFETCH),intent(out) :: var
  integer :: varID,ndims
  integer :: ncFileID,var_date,status
  real :: scale,offset
  character *100 :: period_read

  ndims=3
  if(KMAX==1)ndims=2
!  print *,'  reading ',trim(varname),'from',trim(fileName)
  !open an existing netcdf dataset
  call check(nf90_open(path=trim(fileName),mode=nf90_nowrite,ncid=ncFileID))

  !get varID:
  call check(nf90_inq_varid(ncid=ncFileID,name=trim(varname),varID=VarID))

  !get scale factors
  scalefactors(1) = 1.0 !default
  scalefactors(2) = 0.  !default

  status = nf90_get_att(ncFileID, VarID, "scale_factor", scale  )
  if(status == nf90_noerr) scalefactors(1) = scale
  status = nf90_get_att(ncFileID, VarID, "add_offset",  offset )
  if(status == nf90_noerr) scalefactors(2) = offset

  !find period_of_validity
  period_of_validity='                                     ' !initialisation
  period_read='                                     ' !initialisation
  status = nf90_get_att(ncFileID, VarID, "period_of_validity", period_read  )
  if(status == nf90_noerr)then
     period_of_validity  = trim(period_read)
  else
     period_of_validity='instantaneous' !default
  endif

 ! if(Nfetch<nrecords)then
 !    write(*,*)'Reading record',nstart,' to ',nstart+nfetch-1
 ! endif
 
  !get variable
  if(ndims==2)then
     call check(nf90_get_var(ncFileID, VarID, var,&
          start=(/ISMBEG,JSMBEG,nstart/),count=(/ GIMAX,GJMAX,nfetch /)))
  elseif(ndims==3)then
     call check(nf90_get_var(ncFileID, VarID, var,&
        start=(/ISMBEG,JSMBEG,1,nstart/),count=(/GIMAX,GJMAX,KMAX,nfetch /)))
  endif

 call check(nf90_close(ncFileID))

end subroutine GetCDF_short

   subroutine Getgridparams(meteoname,GRIDWIDTH_M,xp,yp,fi,xm,&
         ref_latitude,sigma_mid,Nhh,nyear,nmonth,nday,nhour,nhour_first)
!
! Get grid and time parameters as defined in the meteo file
! Do some checks on sizes and dates
!
! This routine is called only once (and is therefore not optimized for speed)
!

     use netcdf
     use Par_ml,                only : me,NPROC,ISMBEG,JSMBEG&
                                       ,GIMAX,GJMAX,MAXLIMAX,MAXLJMAX&
                                       ,gi0,gj0
     use ModelConstants_ml,     only : KMAX_MID

     implicit none

     character (len = *), intent(in) ::meteoname
     integer, intent(in):: nyear,nmonth,nday,nhour
     real, intent(out) :: GRIDWIDTH_M,xp,yp,fi, ref_latitude,&
                          xm(0:MAXLIMAX+1,0:MAXLJMAX+1),sigma_mid(KMAX_MID)
     integer, intent(out):: Nhh,nhour_first

     integer :: gc_info,nseconds(1),n1,i,j
     integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid
     integer :: GIMAX_file,GJMAX_file,KMAX_file,ihh,ndate(4)
     real ::xm_global(0:GIMAX+2,0:GJMAX+2)


  if(me==0)then
  print *,'  reading ',trim(meteoname)
  !open an existing netcdf dataset
  call check(nf90_open(path=trim(meteoname),mode=nf90_nowrite,ncid=ncFileID))


  !get dimensions id
  call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
  call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
  call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
  call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

  !get dimensions length
  call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_file))
  call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_file))
  call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_file))
  call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=Nhh))

  write(*,*)'dimensions meteo grid',GIMAX_file,GJMAX_file,KMAX_file,Nhh

  if(GIMAX+ISMBEG-1>GIMAX_file.or.GJMAX+JSMBEG-1>GJMAX_file)then
     write(*,*)'outside domain',GIMAX,GIMAX_file,GJMAX,GJMAX_file
     call gc_abort(me,NPROC,"error in NetCDF_ml")
  endif
  if(KMAX_MID>KMAX_file)then
     write(*,*)'wrong vertical dimension',KMAX_MID,KMAX_file
     call gc_abort(me,NPROC,"error in NetCDF_ml")
  endif
  if(24/Nhh/=METSTEP)then
     write(*,*)'ERROR: METSTEP and meteo step not equal',Nhh,METSTEP
     call gc_abort(me,NPROC,"error in NetCDF_ml")
  endif
  
  
  if(nhour/=0.and.nhour/=3)then
     write(*,*)'WARNING: must start at hour=0 or 3'
     write(*,*)nhour,' not tested'
     call gc_abort(me,NPROC,"error in NetCDF_ml")
  endif
  ihh=1
  n1=1
  call check(nf90_get_var(ncFileID,timeDimID,nseconds,&
             start=(/ihh/),count=(/n1 /)))
  call datefromsecondssince1970(ndate,nseconds(1),0)
  nhour_first=ndate(4)

  if(ndate(1)/=nyear.or.ndate(2)/=nmonth.or.ndate(3)/=nday)then
     write(*,*)'ERROR: wrong meteo date',ndate(1),nyear,ndate(2),&
          nmonth,ndate(3),nday
     call gc_abort(me,NPROC,"error in NetCDF_ml")
  endif

  do ihh=1,Nhh
     call check(nf90_get_var(ncFileID, timeDimID, nseconds,&
                start=(/ ihh /),count=(/ n1 /)))   
     call datefromsecondssince1970(ndate,nseconds(1),0)
     if(mod((ihh-1)*METSTEP+nhour_first,24)/=ndate(4))then
        write(*,*)'ERROR: wrong meteo hour',ihh,Nhh,ndate(4),nhour_first
        call gc_abort(me,NPROC,"error in NetCDF_ml")
     endif
  enddo

   
  !get global attributes
  call check(nf90_get_att(ncFileID,nf90_global,"Grid_resolution",GRIDWIDTH_M))
  call check(nf90_get_att(ncFileID,nf90_global,"ref_latitude",ref_latitude))
  call check(nf90_get_att(ncFileID, nf90_global, "xcoordinate_NorthPole",xp ))
  call check(nf90_get_att(ncFileID, nf90_global, "ycoordinate_NorthPole",yp ))
  call check(nf90_get_att(ncFileID, nf90_global, "fi",fi ))
  
  !get variables
  call check(nf90_inq_varid(ncid=ncFileID, name="map_factor", varID=varID))

  call check(nf90_get_var(ncFileID, varID, xm_global(1:GIMAX,1:GJMAX) &
          ,start=(/ ISMBEG,JSMBEG /),count=(/ GIMAX,GJMAX /)))

  call check(nf90_inq_varid(ncid = ncFileID, name = "k", varID = varID))
  call check(nf90_get_var(ncFileID, varID, sigma_mid ))


  endif !me=0



  call gc_ibcast(199,1,0,NPROC,gc_info,Nhh)
  call gc_rbcast(200,1,0,NPROC,gc_info,GRIDWIDTH_M)
  call gc_rbcast(198,1,0,NPROC,gc_info,ref_latitude)
  call gc_rbcast(201,1,0,NPROC,gc_info,xp)
  call gc_rbcast(202,1,0,NPROC,gc_info,yp)
  call gc_rbcast(203,1,0,NPROC,gc_info,fi)
  call gc_rbcast(204,KMAX_MID,0,NPROC,gc_info,sigma_mid)
  call gc_rbcast(205,GIMAX*GJMAX,0,NPROC,gc_info,xm_global(1:GIMAX,1:GJMAX))

!complete along the four lateral sides
  do i=1,GIMAX
     xm_global(i,0)=xm_global(i,1)
     xm_global(i,GJMAX+1)=xm_global(i,GJMAX)
     xm_global(i,GJMAX+2)=xm_global(i,GJMAX)
  enddo
  do j=0,GJMAX+2
     xm_global(0,j)=xm_global(1,j)
     xm_global(GIMAX+1,j)=xm_global(GIMAX,j)
     xm_global(GIMAX+2,j)=xm_global(GIMAX,j)
  enddo


!keep only part of xm relevant to the local domain
!note that xm has dimensions larger than local domain
  if(MAXLIMAX+1>limax+2.or.MAXLJMAX+1>ljmax+2)then
     call gc_abort(me,NPROC,"error in Met_ml sizes definitions")
  endif
  do j=0,MAXLJMAX+1
     do i=0,MAXLIMAX+1
        xm(i,j)=xm_global(gi0+i-1,gj0+j-1)
     enddo
  enddo


  end subroutine Getgridparams

  subroutine check(status)
    use Par_ml,                only : me,NPROC
    use netcdf
    implicit none
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      call gc_abort(me,NPROC,"error in NetCDF_ml")
    end if
  end subroutine check  
!_______________________________________________________________________

subroutine datefromsecondssince1970(ndate,nseconds,printdate)
  !calculate date from seconds that have passed since the start of the year 1970

!  use Dates_ml, only : nmdays
  implicit none

  integer, intent(out) :: ndate(4)
  integer, intent(in) :: nseconds
  integer,  intent(in) :: printdate

  integer :: n,nday,nmdays(12),nmdays2(13)
  nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/) 

  nmdays2(1:12)=nmdays
  nmdays2(13)=0
  ndate(1)=1969
  n=0
  do while(n<=nseconds)
     n=n+24*3600*365
     ndate(1)=ndate(1)+1
     if(mod(ndate(1),4)==0)n=n+24*3600
  enddo
  n=n-24*3600*365
  if(mod(ndate(1),4)==0)n=n-24*3600
  if(mod(ndate(1),4)==0)nmdays2(2)=29
  ndate(2)=0
  do while(n<=nseconds)
     ndate(2)=ndate(2)+1
     n=n+24*3600*nmdays2(ndate(2))
  enddo
  n=n-24*3600*nmdays2(ndate(2))
  ndate(3)=0
  do while(n<=nseconds)
     ndate(3)=ndate(3)+1
     n=n+24*3600
  enddo
  n=n-24*3600
  ndate(4)=-1
  do while(n<=nseconds)
     ndate(4)=ndate(4)+1
     n=n+3600
  enddo
  n=n-3600
  !    ndate(5)=nseconds-n
  if(printdate>0)then
  write(*,*)'year: ',ndate(1),', month: ',ndate(2),', day: ',&
       ndate(3),', hour: ',ndate(4),', seconds: ',nseconds-n
  endif
end subroutine datefromsecondssince1970

end module met_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
!  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



