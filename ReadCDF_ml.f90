program ReadCDF
  !
  ! Routines for reading and netCDF files 
  ! and writing results in binary xfeltformat
  !
  ! Written by Peter August 2003
  !
  !compile with:
  !f90 -o ReadCDF -L/home/u4/mifahik/netcdf/lib64 -I/home/u4/mifahik/netcdf/include -64 ReadCDF.f90 -lnetcdf

  

  implicit none
!  character (len=*), parameter :: dirroot = '/home/u4/mifapw/nonlin/'
  character (len=*), parameter :: dirroot = '/work/mifapw/last/'
  character (len=*), parameter :: fileName0 =  dirroot//'out_year.nc'

  integer, parameter :: Numberofrecords=1 !increase for out_month or out_day
  character*30  :: varname
  character*30 ::  SIA,  SO2,  Dsox,  Dnox,  Drdn
  character*30 ::  WDsox,  WDnox,  WDrdn,  PM25,  PMco
  integer :: GIMAX,GJMAX,KMAX
  real, dimension(170,133,Numberofrecords) :: var1,var2,var3
  real, dimension(170,133,13) :: tofelt
  real :: factor(170,133)

  real :: ddepsum1,wdepsum1,depsum1,ddepsum2,wdepsum2,depsum2
  real :: x,xmax,xmin
  integer ::n,i,j,k,ij,ntofelt,nin,imax,jmax

  integer ::nrecords,ident(20),lpack,iscal,E,IO_num2=77
  integer*2::iutdpack(20),nyipack(170*133+3)
  real :: maxfield,locmaxfield,fac1

  write(*,*)fileName0

  Dnox='DDEP_OXN'
  WDnox='WDEP_OXN'

  varname=Dnox
  call GetCDF(varname,fileName0,var1,GIMAX,GJMAX,nrecords,ident)

  varname=WDnox
  call GetCDF(varname,fileName0,var2,GIMAX,GJMAX,nrecords,ident)

  var3=(var1+var2)

  nin=1
  tofelt(:,:,1)=var1(:,:,nin)
  tofelt(:,:,2)=var2(:,:,nin)
  tofelt(:,:,3)=var3(:,:,nin)

     n=3
     maxfield=0.
     do j=1,ident(11)
        do i = 1,ident(10)
           if(abs(tofelt(i,j,n))>maxfield)then
              imax=i
              jmax=j
           endif
           maxfield=max(maxfield,abs(tofelt(i,j,n)))
        enddo
     enddo
     write(*,*)'max ',imax,jmax,maxfield
stop
  ntofelt=3

  open(IO_num2,file='utfelt.flt',form='unformatted',status='unknown')

  !     iutdpack(1)=88
  !     iutdpack(2)=1841
  !     iutdpack(3)=2
  !     iutdpack(4)=6
  !     iutdpack(5)=2
  !
  !     iutdpack(6)=606  !parameter code
  !
  !     iutdpack(7)=1000
  !     iutdpack(8)=0
  !     iutdpack(9)=1
  !     iutdpack(10)=IMAX
  !     iutdpack(11)=JMAX
  !     iutdpack(12)=1999
  !     iutdpack(13)=101
  !     iutdpack(14)=600
  !     iutdpack(15)=100*xp1
  !     iutdpack(16)=100*yp1
  !     iutdpack(17)=gridwidth1/100.
  !     iutdpack(18)=fi2
  !     iutdpack(19)=200
  !     iutdpack(20)=0



  do n=1,ntofelt
     iutdpack=ident
     iutdpack(6)=601
     lpack=iutdpack(10)*iutdpack(11)

     if(n>31)stop
!     iutdpack(13)=100+n
     iutdpack(6)=iutdpack(6)+n-1
     maxfield=0.0
     do j=1,iutdpack(11)
        do i = 1,iutdpack(10)
           maxfield=max(maxfield,abs(tofelt(i,j,n)))
        enddo
     enddo
     E = 0
     if(maxfield.gt.(1.e-32))then
        maxfield = maxfield/3.2e4
        locmaxfield = log10(maxfield)
        if(locmaxfield.le.0.)E = -int(-locmaxfield)
        if(locmaxfield.gt.0.)E = int(locmaxfield)+1
     endif

     iutdpack(20)=E
     fac1=10.**(-E)
     write(*,*)fac1,E,maxfield*3.2e4,lpack,lpack/4.
!     lpack = (lpack+3)/4
!     do i=lpack+1,

     write(IO_num2,err=90) (iutdpack(i),i=1,20)
     write(*,*)'wrote ident'
     ij = 0
     do j=1,iutdpack(11)
        do i = 1,iutdpack(10)
           ij = ij+1
           nyipack(ij)=nint(tofelt(i,j,n)*fac1)
        enddo
     enddo

     write(IO_num2,err=90) (nyipack(i),i=1,lpack)
  enddo

  stop

90 continue
  write(*,*)'ERROR'


end program ReadCDF

!_______________________________________________________________________

subroutine GetCDF(varname,fileName,var,GIMAX,GJMAX,nrecords,xfelt_ident)
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
  integer, intent(out) :: GIMAX,GJMAX,nrecords,xfelt_ident(20)
!  real, dimension(:,:,:),intent(out) :: var
  real, dimension(170,133,12),intent(out) :: var

  integer :: varID,status,ndims,alloc_err
  integer :: KMAX_MID,n,KMAX
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
        allocate(values(GIMAX,GJMAX,nrecords,1), stat=alloc_err)
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
     call check(nf90_get_att(ncFileID, VarID, "long_name", attribute))
     print *,'long_name ',attribute
     
     call check(nf90_get_att(ncFileID, VarID, "xfelt_ident", xfelt_ident))
 

     !get variable
     call check(nf90_get_var(ncFileID, VarID, values))
     depsum=0.
     do n=1,nrecords
!        write(*,*)n,values(62,92,n,1)
     do j=1,GJMAX
     do i=1,GIMAX
        var(i,j,n)=values(i,j,n,1)
!        depsum=depsum+values(i,j,1,1)
     enddo
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
