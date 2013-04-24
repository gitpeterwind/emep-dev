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
  real, public, allocatable, dimension(:) :: xlat , xlong
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
character(len=80)  :: comp = ''
character(len=12)  :: proj  = 'PS' ! PS or lonlat
logical ::  print_indices   = .true.  ! Prints i, j with output
logical ::  print_lonlat    = .false. ! Prints lonlat with output
logical ::  skip_fillvalue  = .false. !Skip FillValues (assumed -1.0e10)
integer            :: status

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
      if ( args(i) == "-L" ) print_lonlat  = .true.
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

call readnetcdfiles(innfil,comp)

print *, "PROJ PROJ " , proj
open (10,file=utfil)
do i=1,nx
do j=1,ny
    if ( skip_fillvalue .and. field(i,j,1)< -1.0e10 ) cycle ! crude FillValue
    if( print_lonlat .and. proj =="PS"  ) then ! Allow for 366 days of output
        write (10,"(2f10.4,400es15.7)") xlong2d(i,j),xlat2d(i,j), (field(i,j,n), n=1, ntime)
    else if( print_lonlat ) then ! Allow for 366 days of output
        !print *, i, j, xlong(i),xlat(j), (field(i,j,n), n=1, ntime)
        write (10,"(2f10.4,400es15.7)") xlong(i),xlat(j), (field(i,j,n), n=1, ntime)
    else if( print_indices ) then ! Allow for 366 days of output
        write (10,"(2i4,400es15.7)") i,j, (field(i,j,n), n=1, ntime)
    else
        write (10,"(400es15.7)") (field(i,j,n), n=1, ntime)
    endif
enddo
enddo
close(10)

end

subroutine readnetcdfiles(filename,cdfname) !DS ,field)
use netcdf
use nxy, only : nx, ny, ntime, &
                field, xlong, xlat, xlong2d, xlat2d, crds, proj

!input
character(len=*),intent(in) :: filename, cdfname

!local
integer            :: ncid, ivar,crd
integer            :: xdimID, ydimID,xvarID, yvarID
logical            :: fexist
integer            :: status
real :: tmp

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

status = nf90_get_var(ncid, ivar, field(:,:,:) )
if (status /= nf90_noerr) call handle_cdferr(' ',status)

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
