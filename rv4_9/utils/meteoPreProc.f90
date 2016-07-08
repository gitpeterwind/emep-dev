!=============================================================================
! Program which reads typical EMEP netcdf files, either output or
! inputs. Can cope with 1d or 2d lon/lat descriptions
!
! Needs  netcdf of course..
! Compile examples:
!ifort -o meteoproc -shared-intel -O3 -r8  -I/global/apps/netcdf/4.1.1/include  -L/global/apps/netcdf/4.1.1/lib meteoproc.f90 -lnetcdf -lnetcdff
!
!=============================================================================
! Last modified
! DS Nov 30th 2012 - modified from rd_longlat
!=============================================================================
module nxy              ! Keep the dimension variables here
  implicit none
  integer, public, save :: nx, ny, ntime
  real, public, allocatable, dimension(:,:,:) :: field ! =>field(nx,ny,ntime)
  ! For simple lat, lon based projections we have 1-d arrays
  real, public, allocatable, dimension(:) :: xlat , xlong
  ! For stereographic we will need 2-D, and need a copy of field
  real, public, allocatable, dimension(:,:) :: xlat2d , xlong2d
  real, public, allocatable, dimension(:,:,:) :: field2d
  logical, public, save :: longlat = .false.

  type :: gatt !global attributes from netcdf files
    character(len=50) :: in
    character(len=100) :: outtxt
    real              :: outnum
  end type gatt
  type(gatt), dimension(8) :: gatts
end module nxy
!=============================================================================
program meteoproc
use nxy  ! , only : nx, ny, ntime, field, xlat, xlong, longlat, xlong2d, xlat2d, field2d
use netcdf
implicit none
integer ::i,j,n, ncid
character(len=120) :: OutFileName = ''
character(len=120), dimension(12)  :: args
character(len=12)  :: proj  = 'PS' ! PS or lonlat
!HDD extra 
character(len=4) :: year, cTemp
character(len=120) :: metdir, fname, label
character(len=8):: txtdate
integer,parameter :: NREC=8
integer ::  Tbase   ! Std =18  ! 18.3Temperature threshold for degree-days, use WRI
integer :: jday,im,dd,tt, iyear, ndays, ioerror, status
integer :: idimID, jdimID, tdimID, lon_dimid, lat_dimid
integer :: ivarID, jvarID, tvarID, lon_varid, lat_varid, dgd_varid, temp_varid
real :: scalefac,addoffset

integer,dimension(12),save:: & ! do leaps later
     mdays=(/31,28,31,30,31,30,31,31,30,31,30,31/)

logical :: debug=.true., extra_debug=.false., firstcall=.true.
integer :: dbx = 60, dby = 59  ! debug coords
real, allocatable, dimension(:,:,:) :: DegreeDays, t24h
real, allocatable, dimension(:,:)   :: t2m
real :: x,y, maxdd, meandd
real :: dds = 0.1 ! INERIS shifting factor, B.Bessagnet email May 13th 2012

!---ARGS PROCESSING START --------------------------------------------------
do i = 1, iargc()
     call getarg(i,args(i))
end do

if ( args(1) == "-h" .or. iargc() == 0 ) then
 100 continue
     print *, "--------------------------------------------------------"
     print *, " "
     print *, "Usage : DegreeDayFactors  Tbase year  metdir label"
     print *, " e.g. : DegreeDayFactors  18    2008 /global/work/mifapw/.../metdata_EC/2008   EECCA-2008"
     stop
else
     cTemp = args(1)
     year = args(2)
     metdir = args(3)
     label  = args(4)
end if

if ( len_trim(label)  < 1 ) then
        print *, "Error - no label specified, or too few arguments?"
        goto 100 !oh no, the dreaded goto
end if
!---ARGS PROCESSING END   --------------------------------------------------

OutFileName= "HDD"//trim(cTemp) // "-" // trim(label) // ".nc"

jday = 0
read(unit=year,fmt="(i4)",iostat=ioerror) iyear
read(unit=cTemp,fmt="(i)",iostat=ioerror) Tbase
if( mod( iyear, 4 ) == 0 ) mdays(2) = 29
!/ Test for 29th Feb
if(mdays(2) == 29)then
   fname = trim(metdir)//"/meteo"//year//"0229.nc"
   status = nf90_open(fname,nf90_nowrite,ncid=ncid)
   if ( status /= 0 ) then
      write(*,*)'WARNING!'
      write(*,*)'Did not find 29th February. Will jump over that date'    
      mdays(2) = 28
   endif
endif
ndays = sum( mdays )
print *, "IYEAR year ", iyear, " No days ", ndays, "FILE OUT ", trim(OutFileName)

!/ Test for 1st and last files and 29th Feb

fname = trim(metdir)//"/meteo"//year//"0101.nc"
status = nf90_open(fname,nf90_nowrite,ncid=ncid)
call checkStop(status, "Missing 1st Jan? Or wrong directory? Or forgot to compile with netcdf4 :: "// fname )

fname = trim(metdir)//"/meteo"//year//"1231.nc"
status = nf90_open(fname,nf90_nowrite,ncid=ncid)
call checkStop(status,"Missing 31st Dec? :: "// fname )



do im=1, 12
  do dd=1, mdays(im)

      write(txtdate,'(a4,i2.2,i2.2)')year,im,dd
      fname = trim(metdir)//"/meteo"//txtdate//".nc"

      jday=jday+1

      call readnetcdfiles(fname, "temperature_2m",proj)

      if ( firstcall ) then
        allocate( t2m (nx,ny) )
        allocate( t24h(nx,ny,ndays) )
        allocate( DegreeDays(nx,ny,ndays) )
        if( proj == "PS" )  then
          x = xlong2d(dbx,dby); y = xlat2d (dbx,dby)
        else
          x = xlong(dbx);       y = xlat(dby)
        end if
      end if

      forall(i=1:nx,j=1:ny)

     ! Start collecting energy drivers:
     ! We enforce a min value of 1, and don't allow more than 40 (-22C), since
     ! there are limits as to how much wood can be burned. Also,
     ! might avoid numerical oddities, and since below a certain temperature
     ! emissions are likely supplemented by SNAP-1

        t2m(i,j) = sum( field(i,j,:))/8.0  ! Daily average
        t24h(i,j,jday) = t2m(i,j)

        DegreeDays(i,j,jday) = max( Tbase-t24h(i,j,jday), 1.0)
        DegreeDays(i,j,jday) = min( DegreeDays(i,j,jday) , 40.0 )

      end forall
      maxdd = max(  maxval( DegreeDays(:,:,jday)  ), maxdd )

      if ( debug ) then
        print "(a,a,i4,3f8.3,es12.2)", "DEBUG ", trim(proj), &
          jday, x, y, t2m(dbx,dby), DegreeDays(dbx,dby,jday)
      end if

      firstcall = .false.

end do !id
end do !im
status = nf90_close(ncid)  ! SAFETY?

!!!!!!!!!!!!!!   SCALE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j = 1, ny
do i = 1, nx
  meandd = sum(DegreeDays(i,j,:))/ real( jday )

  if ( meandd > 0.0 ) then
   ! Orig EMEP
   ! DegreeDays(i,j,:) = DegreeDays(i,j,:)/meandd

   ! INERIS/TFMM suggest a "shifting" factor (dds) to modulate the strength
   ! between cold and warm season
    DegreeDays(i,j,:) = ( DegreeDays(i,j,:) + dds*meandd)/&
                                ((dds+1.0)*meandd )
  else
    DegreeDays(i,j,:) = 1.0
  end if

end do 
end do ! 
do jday = 1, ndays
   print "(a,a,i4,3f8.3,es12.2)", "OUTDD ", trim(proj), &
          jday, x, y, t24h(dbx,dby,jday), DegreeDays(dbx,dby,jday)
end do
print "(a,2f12.4)", "OUTDD means ", sum(t24h(dbx,dby,:))/real(ndays), &
           sum(DegreeDays(dbx,dby,:))/real(ndays) 
!!!!!!!!!!!!!!!   WRITE NETCDF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

status = nf90_create(OutFileName,nf90_clobber,ncid=ncid)
!x,y would have been consistent with cdo, but EMEP code expects i,j for now
!status = nf90_def_dim(ncid,"x",nx,idimID)
!status = nf90_def_dim(ncid,"y",ny,jdimID)
status = nf90_def_dim(ncid,"i",nx,idimID)
status = nf90_def_dim(ncid,"j",ny,jdimID)
status = nf90_def_dim(ncid,"time",NF90_UNLIMITED,tdimID)
!WED status = nf90_def_dim(ncid,"daynumber",ndays,tdimID)
!        time = UNLIMITED ; // (8 currently)

!put the attribute
status = nf90_put_att(ncid,ivarID,"units","EMEP grid x coordinate")
status = nf90_put_att(ncid,jvarID,"units","EMEP grid y coordinate")
status = nf90_put_att(ncid,tvarID,"units","Day number")
  
!2d=(/ idimID, jdimID /)
if ( proj == "PS" ) then
 status = nf90_def_var(ncid,"lon",NF90_DOUBLE,(/ idimID,jdimID /),lon_varid)
 status = nf90_def_var(ncid,"lat",NF90_DOUBLE,(/ idimID,jdimID /),lat_varid)
else
 status = nf90_def_var(ncid,"lon",NF90_DOUBLE,(/ idimID /),lon_varid)
 status = nf90_def_var(ncid,"lat",NF90_DOUBLE,(/ jdimID /),lat_varid)
end if
status = nf90_put_att(ncid,lon_varid,"units","degrees_east")
status = nf90_put_att(ncid,lat_varid,"units","degrees_north")
  

!Out 3D
status = nf90_def_var(ncid,"HDD_Facs",NF90_FLOAT,(/ idimID,jdimID,tdimID /),dgd_varid)
status = nf90_put_att(ncid,dgd_varid,"units","Degree day Factors")

status = nf90_def_var(ncid,"temperature_24h",NF90_FLOAT,(/ idimID,jdimID,tdimID /),temp_varid)
status = nf90_put_att(ncid,temp_varid,"units","celsius")

status = nf90_put_att(ncid,NF90_GLOBAL,"title","Degree day factors, " // trim(label))
status = nf90_put_att(ncid,NF90_GLOBAL,"Tbase", trim(cTemp))
status = nf90_put_att(ncid,NF90_GLOBAL,"source","Met data from " // trim(metdir))

! re-use from met file:
do n = 1, size( gatts(:)%in ) 
   if ( gatts(n)%outtxt == "nd" ) then
     print *, "GLOBAL ATTS NUM ", n, trim( gatts(n)%in ),gatts(n)%outnum
     if  ( nint( gatts(n)%outnum ) /= -999 ) then
       status = nf90_put_att(ncid,NF90_GLOBAL,gatts(n)%in,gatts(n)%outnum)
     end if
   else
     print *, "GLOBAL ATTS TXT ", n, trim( gatts(n)%in ),trim( gatts(n)%outtxt )
     status = nf90_put_att(ncid,NF90_GLOBAL,gatts(n)%in,gatts(n)%outtxt)
   end if
end do

call date_and_time( date=txtdate )
status = nf90_put_att(ncid,NF90_GLOBAL,"created_date",txtdate)

!tell the program that the definiton of the variable finished
status = nf90_enddef(ncid)

status = nf90_put_var(ncid,ivarID,(/(i, i = 1, nx)/))
status = nf90_put_var(ncid,jvarID,(/(j,j = 1, ny)/))
! status = nf90_put_var(ncid,tvarID,time)
status = nf90_put_var(ncid,tvarID,(/ (n, n=1,ndays) /) )

if( proj == "PS" ) then
 status = nf90_put_var(ncid,lon_varid,xlong2d)
 status = nf90_put_var(ncid,lat_varid,xlat2d)
else
 status = nf90_put_var(ncid,lon_varid,xlong)
 status = nf90_put_var(ncid,lat_varid,xlat)
end if
status = nf90_put_var(ncid,temp_varid,t24h)
status = nf90_put_var(ncid,dgd_varid,DegreeDays)
status = nf90_close(ncid) 

end program meteoproc
!==============================================================================

subroutine readnetcdfiles(filename,cdfname,proj) !DS ,field)
use netcdf
use nxy !, only : nx, ny, ntime, field, xlong, xlat, xlat2d, xlong2d

!input
character(len=*),intent(in) :: filename, cdfname
character(len=*),intent(inout) :: proj

!local
integer            :: ncid, ivar, ivarLL, i, j
logical            :: fexist, myfirstcall=.true.
integer            :: status
integer :: idimID, timeID

!check whether file exists
inquire(file=filename,exist=fexist)
if (.not. fexist) then
   print *, "OPEN_FILE ::: missing file is :: ", filename 
endif

! apparently it exists
! open the file
status = nf90_open(path = trim(filename), mode = nf90_nowrite, ncid = ncid)
if (status /= nf90_noerr) call handle_cdferr(fname_scen,status)

! now read the desired parameter
! get the variable ID
status = nf90_inq_varid(ncid, trim(cdfname), ivar)
if (status /= nf90_noerr) call handle_cdferr(cdfname,status)

if ( myfirstcall ) then

 ! Try i, j coords
  status = nf90_inq_dimid(ncid, name="i",dimID=idimID)
  status = nf90_inquire_dimension(ncid, dimID=idimID, len=nx )
  status = nf90_inq_dimid(ncid, name="j",dimID=idimID)
  status = nf90_inquire_dimension(ncid, dimID=idimID, len=ny )
  print *, "DIMIDGET I,J ", idimID, nx, ny

  if ( idimID >0 ) then !  PS eg EECCA
   allocate(xlong2d(nx,ny),stat=status)
   allocate(xlat2d(nx,ny),stat=status)
   status = nf90_inq_varid(ncid, "lon",ivarLL)
   status = nf90_get_var(ncid, ivarLL, xlong2d(:,:) )
   status = nf90_inq_varid(ncid, "lat",ivarLL)
   status = nf90_get_var(ncid, ivarLL, xlat2d(:,:) ) 

 else ! Try lat long instead
   longlat = .true.
   print *, "I,J didn't work , try long, lat"
   ! get dimensions:
   status = nf90_inq_dimid(ncid, name="lon",dimID=idimID)
   status = nf90_inquire_dimension(ncid, dimID=idimID, len=nx )
   ! get coor variables :
   status = nf90_inq_varid(ncid, "lon",ivarLL)
   status = nf90_inquire_dimension(ncid, dimID=idimID, len=nx )
   allocate(xlong(nx),stat=status)
   status = nf90_get_var(ncid, ivarLL, xlong(:) )
   print *, "VAR GET lon ", nx, xlong(1), xlong(100), xlong(nx)

   status = nf90_inq_dimid(ncid, name="lat",dimID=idimID)
   status = nf90_inquire_dimension(ncid, dimID=idimID, len=ny )
   allocate(xlat(ny),stat=status)
   status = nf90_inq_varid(ncid, "lat",ivarLL)
   status = nf90_get_var(ncid, ivarLL, xlat(:) )
   print *, "VAR GET lat ", ny, xlat(1), xlat(100), xlat(ny), maxval(xlat)

   proj="lonlat"
  end if ! coord tests
  print *, "Deduced projection" // trim(proj)

  status = nf90_inq_dimid(ncid, name="time",dimID=idimID)
  status = nf90_inquire_dimension(ncid, dimID=idimID, len=ntime )
  print *, "DIMTIME", ntime

  ntime=max(1,ntime) ! can be zero
  allocate(field(nx,ny,ntime),stat=status)

 !/ Find global conventions 
  gatts(:)%in = (/ "Conventions","projection", "projection_params",&
               "xcoordinate_NorthPole", "ycoordinate_NorthPole", &
               "Grid_resolution", "fi", "ref_latitude" /)

  do n = 1, size( gatts(:)%in )
     status = nf90_get_att(ncid,NF90_GLOBAL,gatts(n)%in,gatts(n)%outtxt)
     print *, "READING GLOBAL", status, trim(gatts(n)%in),trim(gatts(n)%outtxt)
     if ( status /= 0 ) then
       gatts(n)%outtxt  = "nd" !needed to be reset
       status = nf90_get_att(ncid,NF90_GLOBAL,gatts(n)%in,gatts(n)%outnum)
       if ( status /= 0 ) gatts(n)%outnum = -999.
       print *, "RE-READING GLOBAL", status, trim(gatts(n)%in), gatts(n)%outnum
     end if
  end do

  myfirstcall = .false.
end if ! myfirstcall

! get the variable itself
status = nf90_get_var(ncid, ivar, field(:,:,:) )

!get attribute scale_factor and add_offset
status = nf90_get_att(ncid,ivar,"scale_factor",scalefac)
status = nf90_get_att(ncid,ivar,"add_offset",addoffset)

if( status == 0 ) then  
  !print "(a,a,3es12.3)", "SCALE ADD ", trim(cdfname),  scalefac, addoffset
  field(:,:,:) =  field(:,:,:)*scalefac + addoffset - 273.14
end if

!if (status /= nf90_noerr) call handle_cdferr(' ',status)

! close the file
status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_cdferr(' ',status)

end subroutine

subroutine handle_cdferr(name,status)

use netcdf
integer :: status
character(len=*) :: name

print *, trim(name),': ', trim(nf90_strerror(status))
! we stop if there is an error!!!!
stop

end subroutine

subroutine checkStop( status, txt )
  integer, intent(in) :: status
  character(len=*) :: txt
   if ( status /= 0 ) then
      print *, "Error, checkStop called ", status
      print *,  trim(txt)
      stop
   end if
end subroutine checkStop
