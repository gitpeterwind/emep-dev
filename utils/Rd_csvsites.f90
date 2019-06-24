program Rdsites
implicit none
!  Reads from new-style sites_yyyy.csv nput file (e.g. sites_2009.csv), and provides various 
!  outputs for a chosen pollutant:
!

character(len=200) :: infile, prefix   ! for GetArg, output label
integer, parameter :: MAXDIM = 600   ! No. variables and /or sites 
character(len=25000) :: longline
character(len=20), dimension(MAXDIM) :: sitename
character(len=20), dimension(MAXDIM) :: specname
real             , dimension(MAXDIM) :: xin
real             , dimension(366,0:23) :: data
integer          , dimension(366) :: cmonth, cday, nhours = 0
character(len=20) :: wpoll, wsite, date, txt1, txt2, txthh, dbg_flag
integer :: i, yy, mm, dd, hh, isite, nsites, n, ispec, nspec, freq
integer :: wanted_site, wanted_spec, last_day, jd, old_dd, wanted_year
integer :: io_in, io_dmax, io_dmean, io_vals, io_hrly, io_dates
integer :: comma, comma2
logical :: first=.true., dbg = .false., force_output=.false.

   if ( iargc() > 0 ) then
     call getarg(1,infile)
     call getarg(2,wsite)
     call getarg(3,wpoll)
     call getarg(4,dbg_flag)  ! -d or -f
     print *, "Input argument: ", infile
     if ( len_trim(wsite)>0 ) print *, "Site wanted ", wsite
     if ( len_trim(wpoll)>0 ) print *, "Poll wanted ", wpoll
     if ( len_trim(dbg_flag)>0 ) then
       if ( dbg_flag == '-d' ) dbg = .true.
       if ( dbg_flag == '-f' ) then
          dbg = .true.
          force_output = .true.
       end if
       print *, "Debug wanted? ", trim(dbg_flag),  dbg
     end if
   end if
   if ( infile == "-h"  .or. iargc() == 0 ) then
     print *, "Usage: Rd_sites sites.mmyy  [site_name] [poll_name] {-d}"
     print *, " (-d triggers extra debug info)"
     print *, " => "
     print *, "     SITES.dmax     - dailymax values  "
     print *, "     SITES.dmean    - dailymean values  "
     print *, "     SITES.hrly     - all hourly values"
     print *, "     SITES.sname    - name of site chosen"
     print *, "     SITES.vals     - all values, printed 24 hourly values per line"
     stop
   end if

   io_in = 10
   open(io_in,file=infile)

   read(unit=io_in,fmt=*)  nsites
   read(unit=io_in,fmt=*)  freq
!============================================================================
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
print *, "NSITES ", nsites, "FREQ ", freq, "WANTED_SITEXX ", wanted_site, trim(wsite), trim(sitename(wanted_site))
   print *, "Chose site", sitename(wanted_site)
!stop
!============================================================================
   wanted_spec = -1
   read(unit=io_in,fmt=*) nspec
  !new read(unit=io_in,fmt=*) n, specname(n)
   read(unit=io_in,fmt=*) txt1, txt2, txthh, (specname(ispec), ispec=1, nspec)
   do n = 1, nspec
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
   prefix="SITES_" // trim(wsite) // "_" // trim(wpoll) 
   io_dmax = 11; open(io_dmax,file=trim(prefix) // ".dmax")
   io_vals = 12; open(io_vals,file=trim(prefix) // ".vals")
   io_hrly = 14; open(io_hrly,file=trim(prefix) // ".hrly")
   io_dmean = 15; open(io_dmean,file=trim(prefix) // ".dmean")
   io_dates = 17; open(io_dates,file=trim(prefix) // ".dates")
!============================================================================


!Data bits start with e.g.:
! AT02_Illmitz, 01/01/2009 01:00, 3.020E+01, 6.452E-03,....

   data(:,:) = 0.0
   jd = 0
   old_dd = 0
   do while (.true.)
     do isite = 1, nsites

       read(io_in,"(a)",err=100,end=100) longline 

       if (isite == wanted_site ) then !NEW HERE

        ! We process by splitting on commas
         comma=index(longline, ",")
         wsite=longline(1:comma-1)      ! site name
         longline=longline(comma+1:)    !
         comma2=index(longline, ",")     ! end of date section
         date= longline(1:comma2-1)      ! date  eg 20/03/16
         longline=longline(comma2+1:)    !
         comma2=index(longline, ",")     ! end of date section
         txthh = longline(1:comma2-1)      ! hour
         read(unit=date(2:3),fmt="(i2)")  dd
         read(unit=date(5:6),fmt="(i2)")  mm
         read(unit=date(8:11),fmt="(i4)")  yy
         read(unit=txthh,fmt="(i2)")  hh
!print *, "TMPA ", comma, wsite, date, txthh, yy, mm, dd, hh
         if ( dd /= old_dd ) then
!print *, "TMPB ", first, nspec, longline(1:20)
           if ( first ) then
             wanted_year = yy
             first = .false.
           end if
           if ( yy > wanted_year ) EXIT
           jd = jd + 1
           cmonth(jd) = mm
           cday(jd) = dd
           old_dd = dd
         end if
         nhours(jd) = nhours(jd)  + 1 
         longline=longline(comma2+1:)
         read(longline,*) xin(1:3)
!print *, "TMPC ", xin(1:3)
         read(longline,*) xin(1:nspec)
!print *, "TMP ", comma, wsite, date, txthh, yy, mm, dd, hh
!stop 'TMP'

     ! ----------------------------------
     !OLD  if (isite == wanted_site ) then
         last_day = jd
         if ( dd==15 .and. hh==12) write(*,*) "date ", date, "HH ", mm, dd, hh, " jd ", jd, last_day
         data(jd,hh) = xin(wanted_spec)
         if( jd == 1 .and. hh == freq) then
                 data(jd,0) = data(jd,hh)
                 nhours(jd) = 1  ! Initialise
         end if
      end if
     end do ! isite
   end do ! records
100 continue

    !write(unit=io_snam,fmt="(a20)") sitename(wanted_site)
    if ( nhours(last_day) < 24 ) then
       if ( dbg ) print *, "Incomplete last day ", last_day, nhours(last_day)
       if ( .not. force_output ) last_day = last_day  - 1
    end if
    do jd = 1, last_day

       ! The file SITES.vals has all values for one day in one record.
       ! We format with integers if possible, otherwise es

        if( maxval(data) < 10.0 ) then
          write(unit=io_vals,fmt="(24es9.2)") ( data(jd,hh), hh= 0,23)
        else if( maxval(data) < 100.0 ) then
          write(unit=io_vals,fmt="(24i3)") ( int(data(jd,hh)), hh= 0,23)
        else
          write(unit=io_vals,fmt="(24i5)") ( int(data(jd,hh)), hh= 0,23)

        end if

       ! Max, Mean
        print *, 'MaxDay', jd, nhours(jd), maxval(data(jd,:))
        write(unit=io_dmax,fmt="(g14.4)") maxval(data(jd,:))
        write(unit=io_dmean,fmt="(g14.4)") sum(data(jd,:))/24.0
        write(unit=io_dates,fmt="(4i4)") yy, cmonth(jd), cday(jd), jd
       ! SITES.hrly
        do hh = 0, 23, freq
           write(unit=io_hrly,fmt="(4i3,g14.5)") jd,cmonth(jd), cday(jd), hh,data(jd,hh)
        end do
    end do

end program Rdsites
