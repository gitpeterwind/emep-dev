!This program is a tool to make F90 routines for reading and processing of netCDF files
!
!for details see: 
!http://www.unidata.ucar.edu/packages/netcdf/f90/Documentation/f90-html-docs/
!
!
!compile with:
!f90 -o ReadCDF -L/home/u4/mifahik/netcdf/lib64 -I/home/u4/mifahik/netcdf/include -64 ReadCDF_ml.f90 -lnetcdf
!
!view results.nc with:
!xrdb -load /home/u4/mifahik/.app-defaults/Ncview  (once only)
!/home/u4/mifahik/bin/ncview results.nc
!or
!/home/u4/mifahik/bin/ncdump results.nc |less
!____________________________________________________________________________________________________________

   PROGRAM ReadCDF_ml
!
! Routines for reading and manipulating netCDF files
!
! Written by Peter january 2003
!

  character (len=*), parameter :: fileName1 = 'results.nc'
  character (len=*), parameter :: fileName2 = 'results00.nc'
  character*30  :: varname
  integer :: GIMAX,GJMAX,KMAX
  real, dimension(170,133,20) :: var1,var2
  real :: factor(170,133)

  real :: ddepsum1,wdepsum1,depsum1,ddepsum2,wdepsum2,depsum2

  varname='Out_2D0521'
  call GetCDF(varname,fileName1,var1,GIMAX,GJMAX,KMAX)

     call mgm2toton(factor,GIMAX,GJMAX)
     depsum1=0.
     do j=1,GJMAX
     do i=1,GIMAX
        ddepsum1=ddepsum1+var1(i,j,1)*factor(i,j)
     enddo
     enddo
     print *, 'sum dry dep in ',fileName1,ddepsum1
  varname='Out_2D0541'
  call GetCDF(varname,fileName1,var1,GIMAX,GJMAX,KMAX)
     depsum1=0.
     do j=1,GJMAX
     do i=1,GIMAX
        wdepsum1=wdepsum1+var1(i,j,1)*factor(i,j)
     enddo
     enddo
     print *, 'sum wet dep in ',fileName1,wdepsum1
     depsum1=ddepsum1+wdepsum1
     print *, 'sum wet+dry dep in ',fileName1,depsum1

  varname='Out_2D0521'
  call GetCDF(varname,fileName2,var2,GIMAX,GJMAX,KMAX)
     depsum2=0.
     do j=1,GJMAX
     do i=1,GIMAX
        ddepsum2=ddepsum2+var2(i,j,1)*factor(i,j)
     enddo
     enddo
     print *, 'sum dry dep in ',fileName2,ddepsum2
  varname='Out_2D0541'
  call GetCDF(varname,fileName2,var2,GIMAX,GJMAX,KMAX)
     depsum2=0.
     do j=1,GJMAX
     do i=1,GIMAX
        wdepsum2=wdepsum2+var2(i,j,1)*factor(i,j)
     enddo
     enddo
     print *, 'sum wet dep in ',fileName2,wdepsum2
     depsum2=ddepsum2+wdepsum2
     print *, 'sum wet+dry dep in ',fileName2,depsum2
     print *
     print *, 'diff of wet+dry dep (tons) in ',fileName1,' and ',fileName2,' is ',depsum2-depsum1


   end PROGRAM ReadCDF_ml



!_______________________________________________________________________

subroutine GetCDF(varname,fileName,var,GIMAX,GJMAX,kmax)
!
! open and reads CDF file
!
! The nf90 are functions which return 0 if no error occur.
! check is only a subroutine which check wether the function returns zero
!
!
  use typeSizes
  use netcdf
  implicit none

  character (len=*),intent(in) :: fileName 

  character (len = *),intent(in) ::varname
  integer, intent(out) :: GIMAX,GJMAX,kmax
!  real, dimension(:,:,:),intent(out) :: var
  real, dimension(170,133,20),intent(out) :: var

  integer :: varID,nrecords,status,ndims,alloc_err
  integer :: KMAX_MID  
  integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,VarID,iVarID,jVarID,kVarID,i,j,k

  real (kind = FourByteReal), allocatable,dimension(:,:,:,:)  :: values
  real ::depsum
  character*20::attribute,attribute2

!open an existing netcdf dataset
  call check(nf90_open(path = fileName, mode = nf90_nowrite, ncid = ncFileID))

!get global attributes

  !example:
  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_hour", attribute ))
  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_date", attribute2 ))
  print *,'file last modified (yyyymmdd hhmmss.sss) ',attribute2,' ',attribute

!test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)
  
  if(status == nf90_noerr) then     
     print *, 'variable exists: ',varname
  else
     print *, 'variable does not exist: ',varname,nf90_strerror(status)
     return
  endif

     !get dimensions id
     call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

     !get dimensions length
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=idimID,  len=GIMAX))
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=jdimID,  len=GJMAX))
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=kdimID,  len=KMAX_MID))
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=timedimID,  len=nrecords))

     print *, 'dimensions ',GIMAX,GJMAX,KMAX_MID,nrecords

     !get variable info
     call check(nf90_inquire_variable(ncFileID, varID, ndims=ndims))
!     print *, 'dimensions ',ndims
     if(ndims==3)then
        kmax=1
        !allocate a 2D array 
        allocate(values(GIMAX,GJMAX,1,nrecords), stat=alloc_err)
        if ( alloc_err /= 0 ) then
           print *, 'alloc failed in ReadCDF_ml: ',alloc_err,ndims
           stop
        endif
     elseif(ndims==4)then
        kmax=KMAX_MID
        !allocate a 3D array 
        allocate(values(GIMAX,GJMAX,KMAX_MID,nrecords), stat=alloc_err)
        if ( alloc_err /= 0 ) then
           print *, 'alloc failed in ReadCDF_ml: ',alloc_err,ndims
           stop
        endif
       
     else
        print *, 'unexpected number of dimensions: ',ndims
        stop
     endif

     !get variable attributes
     !example:
     attribute=''
     call check(nf90_get_att(ncFileID, VarID, "periodlength", attribute))
     print *,'periodlength ',attribute
     
     !get variable
     call check(nf90_get_var(ncFileID, VarID, values))
     depsum=0.
     do j=1,GJMAX
     do i=1,GIMAX
        var(i,j,1)=values(i,j,1,1)
        depsum=depsum+values(i,j,1,1)
     enddo
     enddo
!        print *, values(51,64,1,1), values(50,64,1,1), values(49,63,1,1), values(49,62,1,1)
!        print *, 'sum dep ',depsum
!close file
  call check(nf90_close(ncFileID))

   end subroutine GetCDF
!_______________________________________________________________________

  subroutine check(status)

  use typeSizes
  use netcdf
    implicit none
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      print *, "error in NetCDF_ml"
      stop
    end if
  end subroutine check  

  subroutine mgm2toton(factor,GIMAX,GJMAX)

! GRIDWIDTH_M = 50000.0 meter
! AN = 6.371e6*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M = 237.768957  
! xp,yp:   coord. of the polar point (read from fields)
!           xp = ident(15)/100. (= 43.00  ?)
!           yp = ident(16)/100. (= 121.00 ?)
! PI = 3.14159265358979323
! i,j map coordinates (1<= i <=170 , 1<= j <=133 ?)

implicit none

real :: factor(170,133)
real :: gridwidth,PI,AN,xp,yp,an2,y,x,rpol2
real :: xm(170,133)

integer ::i,j,GIMAX,GJMAX

if(GIMAX/=132 .or. GJMAX/=111)then
   print *,'WARNING: mapping factor probably wrong!'
endif

 !  map factor = xm
  gridwidth = 50000.
  PI = 3.14159265358979323
  AN = 237.768957
  xp = 43.00-35.  !assumes  a domain @smalldomain = (36, 167,  12, 122)
  yp = 121.00-11. !
  an2 = AN*AN
  do j=1,GJMAX        
     y = j - yp     
     do i=1,GIMAX
        x = i - xp   
        rpol2 = (x*x + y*y)/an2
        xm(i,j) = 0.5*(1.0+sin(PI/3.0))*(1.0 + rpol2)
!        write(*,*)i,j,xm(i,j)
     end do
  end do

!convert xm into a factor which converts mg/m^2 into tons
  do i=1,GIMAX
     do j=1,GJMAX
        factor(i,j)=gridwidth**2/xm(i,j)**2/1000000000
     end do
  end do

  end subroutine mgm2toton
