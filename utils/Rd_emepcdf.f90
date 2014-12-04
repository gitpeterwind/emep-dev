!=============================================================================
! Program which reads typical EMEP netcdf files, either output or
! inputs. Can cope with 1d or 2d lon/lat descriptions
!
! Needs  netcdf of course..
! Compile examples:
! Netcdf4 (stallo, Vilje):
! ifort -o Rd_emepcdf -shared-intel -O3 -r8  -I/global/apps/netcdf/4.1.1/include  -L/global/apps/netcdf/4.1.1/lib Rd_emepcdf.f90 -lnetcdf -lnetcdff

!=============================================================================
! Last modified
! DS Jan 24th 2013
!=============================================================================
module nxy
  implicit none
  ! Keep the dimension variables here
  integer, public, save :: nx, ny, ntime
  real, public, allocatable, dimension(:,:,:) :: field ! =>field(nx,ny,1)
  ! For simple lat, lon based projections we have 1-d arrays
  real, public, allocatable, dimension(:) :: xlat , xlong, var1d
  ! For stereographic we will need 2-D, and need a copy of field
  real, public, allocatable, dimension(:,:) :: xlat2d , xlong2d
  logical, public, save :: longlat = .false.
  character(len=20), public, save :: proj

 ! We may need to try several types of coordinate dimensions. These are tested
 ! so far:
  type, public :: coord
    character(len=10) :: x, y
  end type coord
  type(coord), public, dimension(3), parameter :: crds = (/ &
    coord( "i", "j" ),  &   ! 
    coord( "x", "y" ),  &   ! x,y often output after cdo processing
    coord( "lon", "lat" ) &  !
   /)
    
end module nxy

program readfield
use nxy, only : nx, ny, ntime, field, xlat, xlong, longlat, xlong2d, xlat2d
implicit none
integer ::i,j,n
character(len=120) :: innfil = '', utfil = ''
character(len=120), dimension(12)  :: args
character(len=80)  :: comp = '', fltfmt = "es15.7", ofmt  ! default outout format
character(len=20)  :: sitexy = '-'     ! default site txt (e.g. "52.0,11.2")
character(len=12)  :: proj  = 'NotSet' ! PS or lonlat
logical ::  print_indices   = .true.  ! Prints i, j with output
logical ::  print_lonlat    = .false. ! Prints lonlat with output
logical ::  skip_fillvalue  = .false. !Skip FillValues (assumed -1.0e10)
logical ::  first_record    = .true. 
integer :: status

do i = 1, iargc()
     call getarg(i,args(i))
end do

if ( args(1) == "-h" .or. iargc() == 0 ) then
 100 continue
     print *, "--------------------------------------------------------"
     print *, " "
     print *, "Usage : readfiles.exe -i infile -c comp [-o utfil] [-d]"
     print *, "e.g.  : readfiles.exe -i /global/../Base_fullrun.nc  -c D2_AOT40"
     print *, " "
     print *, "Required: "
     print *, "   -i  :  name of input file"
     print *, "   -c  :  name of variable, e.g. D2_SO4"
     print *, " "
     print *, "Optional: "
     print *, "   -d  :  only print  data in output, not i,j"
     print *, "   -f  :  format for output of numbers, e.g. es10.3, f8.3"
     print *, "   -o  :  name of output file, default is OUT_{comp}.txt"
     print *, "   -L  :  output lon, lat, data"
     print *, "   -v  :  print only valid values (assumed > 1-0e10 for now)"
     print *, "--------------------------------------------------------"
     stop
else
   do i = 1, iargc()
      if ( trim( args(i) ) == "-i" ) innfil = args(i+1)
      if ( args(i) == "-o" ) utfil  = args(i+1)
      if ( args(i) == "-c" ) comp   = args(i+1)
      if ( args(i) == "-d" ) print_indices = .false.
      if ( args(i) == "-f" ) fltfmt = trim( args(i+1) )
      if ( args(i) == "-L" ) print_lonlat  = .true.
      if ( args(i) == "-s" ) sitexy = trim( args(i+1) )
      if ( args(i) == "-v" ) skip_fillvalue  = .true.
      print *, "ARGS ", i, trim(args(i)), print_indices
   end do
end if

!Cope with missing flags or arguments
if( len_trim(comp) == 0 ) then
        print *, "Error - no component specified"
        goto 100 !oh no, the dreaded goto
end if
if( len_trim(utfil) == 0 ) utfil = "OUT_" // trim(comp) // ".txt"

write(*,*) innfil, utfil, print_indices

call readnetcdfiles(innfil,comp,sitexy)

print *, "PROJ PROJ " , proj

! Check for format errors, allowing for 1 space for +ve values, 2 for neg
ofmt = "(" // trim(fltfmt) // ")" 
write(ofmt,fmt=ofmt) max( 10*maxval(field), abs(100*minval(field)) )
if( index(ofmt,"*") > 0 ) then
  print *, "FMT ERROR: ", trim(fltfmt),  " too small for maxval+1: ", maxval(field)
  stop 'Use -f option to change fmt'
end if

ofmt = "(a," // trim(fltfmt) // ")" 
print ofmt, "MAX VALUE: ", maxval(field)
print ofmt, "MIN VALUE: ", minval(field)

!write(*,fmt=ofmt,iostat=status) "MAX VALUE: ", maxval(field)
write(fltfmt,"(i0.0,a)") ntime, trim( fltfmt )  ! Start to build outout format
ofmt = "(2f10.4," // trim(fltfmt) // ")"        ! default, reset if needed below

open (10,file=utfil)
do i=1,nx
do j=1,ny
    if ( skip_fillvalue .and. field(i,j,1)< -1.0e10 ) cycle ! crude FillValue
    !if( print_lonlat .and. proj =="NotSet"  ) then
    ! Use 2d lat/lon for all projs. Should be okay.
    if( print_lonlat ) then
        write (10,ofmt) xlong2d(i,j),xlat2d(i,j), (field(i,j,n), n=1, ntime)
    !else if( print_lonlat ) then
    !    !print *, i, j, xlong(i),xlat(j), (field(i,j,n), n=1, ntime)
    !    write (10,fmt=ofmt) xlong(i),xlat(j), (field(i,j,n), n=1, ntime)
    else if( print_indices ) then
        if(first_record) ofmt = "(2i4," // trim(fltfmt) // ")"
        write (10,fmt=ofmt) i,j, (field(i,j,n), n=1, ntime)
    else
        if(first_record) ofmt = "(" // trim(fltfmt) // ")"
        write (10,fmt=ofmt) (field(i,j,n), n=1, ntime)
    endif
    if( first_record) print *, "Output fmt:", trim(ofmt)
    first_record=.false.
enddo
enddo
close(10)

end

subroutine readnetcdfiles(filename,cdfname,sitexy) !DS ,field)
use netcdf
use nxy, only : nx, ny, ntime, &
                field, xlong, xlat, var1d, xlong2d, xlat2d, crds, proj

!input
character(len=*),intent(in) :: filename, cdfname, sitexy

!local
integer            :: ncid, ivar,crd
integer            :: xdimID, ydimID,xvarID, yvarID
real               :: xsite, ysite, rdist, rdist2
integer            :: xmin, ymin, io_s
logical            :: fexist
integer            :: status
real :: scalefac=1.0, addoffset=0.0
real :: tmp
!ROT TEST. Had 1-d i, j  float arrays with rot lat/lon. Not sure
!what that means though in real lat/long. Skip code for now.
integer :: nvar1d

!check whether file exists
inquire(file=filename,exist=fexist)
if (.not. fexist) then
   print *, "OPEN_FILE ::: missing file is :: ", filename 
endif

! apparently it exists
! open the file
status = nf90_open(path = trim(filename), mode = nf90_nowrite, ncid = ncid)
if (status /= nf90_noerr) call handle_cdferr(fname_scen,status)

print *, "OPENED - GOT TO HERE"
! now read the desired parameter
! get the variable ID
status = nf90_inq_varid(ncid, trim(cdfname), ivar)
if (status /= nf90_noerr) call handle_cdferr(cdfname,status)

do crd = 1, size(crds)
  status = nf90_inq_dimid(ncid, name=crds(crd)%x,dimID=xdimID)
  status = nf90_inquire_dimension(ncid, dimID=xdimID, len=nx )
  status = nf90_inq_dimid(ncid, name=crds(crd)%y,dimID=ydimID)
  status = nf90_inquire_dimension(ncid, dimID=ydimID, len=ny )
  print *, "DIMIDSEARCH  crd ",  crd, crds(crd), xdimID, ydimID, nx, ny
  if( ydimID >0  ) exit
end do

print *, "DIMIDFOUND  ", crds(crd), xdimID, nx, ny
 

! We assign a 2-D long/lat to all options for simplicity
allocate(xlong2d(nx,ny),stat=status)
allocate(xlat2d(nx,ny),stat=status)

status = nf90_inq_varid(ncid, "lon",xvarID)
status = nf90_inq_varid(ncid, "lat",yvarID)

proj="PS" ! default guess

if ( crds(crd)%y == "lat" ) then
   proj="lonlat"
   longlat = .true. 
   allocate(xlong(nx),stat=status)
   allocate(xlat(ny),stat=status)
   status = nf90_get_var(ncid, xvarID, xlong(:) )
   status = nf90_get_var(ncid, yvarID, xlat(:) )
   do i = 1, nx
     ! For plotting we usually want -180 to 180, whereas e.g. IPCC use 1-360
     !tmp = xlong(i)
     if( xlong(i) > 180.0 ) xlong(i) = xlong(i) - 360.0
     !print *, "LONG SWITCH", i, tmp, xlong(i)
     xlong2d(i,:) = xlong(i)
   end do
   do j = 1, ny
     xlat2d(:,j) = xlat(j)
   end do
   print *, "TEXT1D LL ", xlong2d(2,2), xlat2d(2,2)
else
   status = nf90_get_var(ncid, xvarID, xlong2d(:,:) )
   status = nf90_get_var(ncid, yvarID, xlat2d(:,:) )
   print *, "TEXT2D LL ", xlong2d(2,2), xlat2d(2,2)
end if


status = nf90_inq_dimid(ncid, name="time",dimID=idimID)
status = nf90_inquire_dimension(ncid, dimID=idimID, len=ntime )
print *, "DIMTIME", ntime
ntime=max(1,ntime) ! can be zero
allocate(field(nx,ny,ntime),stat=status)

! get the variable itself

!ROT LATLON TEST. Have to cope with i,j as coordinates (i=1..85 for RCA) and
! i, j as variables.
status = nf90_inquire_variable(ncid, ivar, ndims=nvar1d )
print *, "TEST VAR NDIMS = ", nvar1d
if( nvar1d == 1 ) then
   stop 'Outout of 1-D varoables not coded yet'
!  status = -1
!  status = nf90_inquire_dimension(ncid, dimID=ivar, len=nvar1d ) !Obs! reset nvar1d
!  allocate( var1d(nvar1d) )
!  status = nf90_get_var(ncid, ivar, var1d(:) )
!  if( status < 0 ) stop 'TOT ij failed'
!  print *, "TEXTVAR RCA?",  nvar1d, var1d(1), var1d(2), var1d(nvar1d)
!!ROT LATLON TEST
else
  status = nf90_get_var(ncid, ivar, field(:,:,:) )
end if
if (status /= nf90_noerr) call handle_cdferr(' ',status)
print *, "TESTVAR POINTS?",  nvar1d, field(1,1,1), field(2,2,1)

! get attribute scale_factor and add_offset
status = nf90_get_att(ncid,ivar,"scale_factor",scalefac)
if( status == 0 ) then
  tmp = field(2,2,1)
  field(:,:,:) =  field(:,:,:)*scalefac
  print "(a,3es12.3)", "SCALE FAC "//trim(cdfname),  tmp, scalefac, field(2,2,1)
end if

status = nf90_get_att(ncid,ivar,"add_offset",addoffset)
if( status == 0 ) then
   field(:,:,:) =  field(:,:,:) + addofset
   print "(a,3es12.3)", "SCALE ADD "//trim(cdfname),  tmp, addoffset, field(2,2,1)
end if


! close the file
status = nf90_close(ncid)
if (status /= nf90_noerr) call handle_cdferr(' ',status)

if( sitexy /= "-" )  then
   open(newunit=io_s,file="OUT_Site.txt")
   rdist = 1.0e10
   read( sitexy, * ) xsite, ysite
   print *, "SITEXY ", xsite, ysite , trim(cdfname)
   print *, "SITEXLON ", xlong2D(1,1), xlong2D(nx,ny)
   print *, "SITEXLAT ", xlat2D(1,1), xlat2D(nx,ny)
   do j = 1, ny
   do i = 1, nx
     if( abs(xlong2D(i,j) - xsite) > 5.0 ) cycle
     if( abs(xlat2D(i,j) - ysite) > 5.0 ) cycle
       rdist2 = (xlong2D(i,j) - xsite)**2 + (xlat2D(i,j) - ysite)**2
       if(rdist2 < rdist ) then
         rdist=rdist2 ! still squared
         xmin = i; ymin = j
         print "(a,2i4,5f8.3)", "Site searching ", i,j, xsite, ysite,&
             xlong2D(i,j), xlat2D(i,j), sqrt(rdist)
       end if
   end do
   end do
   i = xmin; j=ymin
   write (*,"(a,2f8.3,2i4,2f8.3)") "#Sites:", xsite, ysite, i,j, xlong2d(i,j),xlat2d(i,j)
   write (io_s,"(a,2f8.3,2i4,2f8.3)") "#Sites:", xsite, ysite, i,j, xlong2d(i,j),xlat2d(i,j)
   do n=1, ntime
     write (*,"(a,i5,f12.3)") "Sites  ", n,  field(i,j,n)
     write (io_s,"(a,i5,f12.3)") "Sites  ", n,  field(i,j,n)
   end do
end if

end subroutine

subroutine handle_cdferr(name,status)

use netcdf
integer :: status
character(len=*) :: name

print *, trim(name),': ', trim(nf90_strerror(status))
! we stop if there is an error!!!!
stop

end subroutine
