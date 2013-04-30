implicit none
!  Reads from input file (e.g. model file sondes_2007.csv say), and provides 
!  various outputs for a chosen pollutant.
!
!   SONDES.vals - provides output for each hour (by row), with here NK=20
!                 columns representing 20 height levels.
!   SONDES.mmean - monthly means for each k level
!   SONDES.dmax  - mean of daily max for each k level
!


character(len=100) :: infile   ! for GetArg
integer, parameter :: IO_5=5, IO_IN = 10,  &
   IO_SNAM=12, IO_DMAX=14, IO_MMEAN=15, IO_VALS=16, IO_DATES= 17, &
   NK = 20, &   ! Number of vertical levels to be printed
   MAX_SITES= 100, MAX_SPEC = 200
character(len=35000) :: longline
character(len=20), dimension(MAX_SITES) :: sitename
character(len=20), dimension(MAX_SPEC) :: specname
real             , dimension(NK,MAX_SPEC) :: xin
real             , dimension(366,0:23,NK) :: data
integer          , dimension(366) :: cmonth, cday
real :: maxx = -9.99e19   ! Max xin value found in 1:NK range 
character(len=20) :: fmt
character(len=20) :: insite, wpoll, wsite, date, txt1, txt2
integer :: comma, comma2
integer :: yy, mm, dd, hh, freq, last_day, jd, old_dd
integer :: i, k, n, nspec, nsites, isite  !csv, ispec
integer :: wanted_site, wanted_spec
logical :: DEBUG = .false.
integer :: old_mm, nmm(12)
real :: mmean(12)

   if ( iargc() > 0 ) then
     call getarg(1,infile)
     call getarg(2,wsite)
     call getarg(3,wpoll)
     print *, "Input argument: ", infile
     if ( len_trim(wsite)>0 ) print *, "Site wanted ", wsite
     if ( len_trim(wpoll)>0 ) print *, "Poll wanted ", wpoll

   end if
   if ( infile == "-h" .or. iargc() == 0 ) then
     print *, "Usage: Rd_csvsondes sondes.mmyy  [site_name] [poll_name]"
     print *, " => "
     print *, "     SONDES.dmax      - daily-max values for NK (e.g. 15) levels"
     print *, "     SONDES.mmean     - monthly-mean values for NK (e.g. 15) levels"
     print *, "     SONDES.sname    - name of site chosen"
     print *, "     SONDES.vals     - all values, printed 24 hourly values per line"
     print *, "------------------------------"
     print *, "See code for details of formats"
     stop
   end if

  !In
   open(IO_IN,file=infile)

  !Out
   open(IO_VALS,file="SONDES.vals")
   open(IO_SNAM,file="SONDES.sname")
   open(IO_DMAX,file="SONDES.dmax")
   open(IO_MMEAN,file="SONDES.mmean")
   open(IO_DATES,file="SONDES.dates")

   read(IO_IN,*)  nsites
   read(IO_IN,*)  freq
!============================================================================
!IN: 42  sondes in domain   1 360   1 180
!IN:  1 Hours between outputs
!IN:Uccle                                             , 184, 141, 141
!IN:...
!IN:Land40N55E                                        , 235, 131, 131

   do isite = 1, nsites
     read(unit=io_in,fmt=*)   sitename(isite)
     if ( len_trim(wsite)>0 ) then
        if ( trim(wsite)==trim( sitename(isite))) wanted_site = isite
     else
        print *, isite, sitename(isite)
     end if
   end do
   if ( len_trim(wsite)==0 ) then
      print *, "Choose site"
      read(unit=5,fmt=*) wanted_site
   end if
   if ( wanted_site <1) then
     print *, "NO SITE OF THAT NAME, START AGAIN!"
     do isite = 1, nsites
        print *, isite, sitename(isite)
     end do
     stop
   end if
   print *, "Chose site", sitename(wanted_site)
!============================================================================
!IN: 16 Variables units: ppb
!IN:site,date,O3,NO2,NO,PAN,NO3_C,NO3_F,SO4,NH4_F,NH3,OH,OD,OP,NOy,z_mid,p_mid,th
  wanted_spec = -1
   read(unit=io_in,fmt=*) nspec
   read(unit=io_in,fmt=*) txt1, txt2, (specname(n), n = 1, nspec)
   do n = 1, nspec
     !csv     read(unit=io_in,fmt=*) n, specname(n)
     if ( len_trim(wpoll)>0 ) then
       if ( trim(wpoll)==trim( specname(n))) wanted_spec = n
     else
       print *,           n, specname(n)
     end if
   end do
   if ( len_trim(wpoll)==0 ) then
     print *, "Choose species"
     read(unit=5,fmt=*) wanted_spec
   end if
   if ( wanted_spec < 1 ) then
     print *, "NO POLL OF THAT NAME, START AGAIN!"
     print "(5(i4,a1,a))", (i, " ", trim(specname(i)), i=1,nspec)
     stop
   end if
   print *, "Site ", wanted_site, " Poll ", specname(wanted_spec)
!============================================================================
!IN:Uccle, 01/01/2006 01:00, 3.282E+01, 3.591E+01, 3.624E+01, 3.709E+01, 3.811E+01, 
!IN:..

   data(:,:,:) = 0.0
   jd = 0
   old_dd = 0


   do while (.true.)
     old_mm = -1
     do isite = 1, nsites
     ! ----- new csv -----------------------------
      read(io_in,"(a)",end=100) longline ! wsite, txt1, txt2 !, v1

       comma=index(longline, ",")
       wsite=longline(1:comma-1)      ! site name
       longline=longline(comma+1:)    !
       comma2=index(longline, ",")     ! end of date section
       date= longline(1:comma2-1)      ! date
       date=adjustl(date)
       read(unit=date(1:2),fmt="(i2)")  dd
       read(unit=date(4:5),fmt="(i2)")  mm
       read(unit=date(7:10),fmt="(i4)")  yy
       read(unit=date(12:13),fmt="(i2)")  hh
       if ( dd /= old_dd ) then
          jd = jd + 1
          cmonth(jd) = mm
          cday(jd) = dd
          old_dd = dd
       end if
       longline=longline(comma2+1:)
      !print *, "longline ", longline(1:40)
       read(longline,*) ( xin(1:NK,n), n=1, nspec )

      if( DEBUG ) print *, insite, yy, mm, dd, hh
!stop "HERE"

        if (isite == wanted_site ) then
           last_day = jd
           data(jd,hh,:) = xin(:,wanted_spec)
           maxx = max( maxval(xin(1:NK,wanted_spec)), maxx )
           if(DEBUG) then
             print "(a,a,6i5,3f12.3)", "SET:" , trim(date), &
              jd, last_day, dd, mm, yy, hh, &
               xin(1,wanted_spec), xin(20,wanted_spec), maxx
           end if
        end if
     end do
   end do
100 continue

    write(IO_SNAM,"(a20)") sitename(wanted_site)
    fmt = "(4i3,30f8.3)"    ! Allow up to 30 fields
    if ( maxx > 1.0e3 ) fmt = "(4i3,30f8.2)"
    if ( maxx > 1.0e4 ) fmt = "(4i3,30es10.2)"

   ! Usually the zeroth hour of the 1st month is not printed. We
   ! make an approximation here:
    data(1,0,:) =  data(1,1,:)

    do jd = 1, last_day

      do hh = 0, 23, freq
        write(IO_VALS,fmt=fmt) jd, cmonth(jd), cday(jd),hh, (data(jd,hh,k), k= 1, NK )
      end do

      write(IO_DMAX,"(3i4,30es10.3)") jd, cmonth(jd), cday(jd), (maxval(data(jd,:,k)), k=1,NK)


    end do ! dd

    do k = NK, 1, -1
      old_mm = -1
      mmean(:) = 0.0
      nmm(mm) = 0
      do jd = 1, last_day
        mm = cmonth(jd)
        if( mm /= old_mm ) then
           nmm(mm) = 0
           old_mm = mm
        end if
        nmm(mm) = nmm(mm) + 1
        mmean( mm ) =  mmean( mm )  + sum( data(jd,:, k)/24.0 )
       !print *, "Last ", jd, last_day
      end do
      write(IO_MMEAN,"(i4, 12f10.3)") k, (mmean(mm)/nmm(mm), mm= 1, 12) 
    end do

end
