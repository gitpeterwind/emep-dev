
       module NetCDF_ml
!
! Routines for netCDF output
!
! Written by Peter january 2003
!
!compile with options:
!f90 -L/home/u4/mifahik/netcdf/lib64 -I/home/u4/mifahik/netcdf/include -64 NetCDF_ml.f90 -lnetcdf
!
!view results.nc with:
!xrdb -load /home/u4/mifahik/.app-defaults/Ncview  (once only)
!/home/u4/mifahik/bin/ncview results.nc
!or
!
!/home/u4/mifahik/bin/ncdump results.nc |less
!
!for details see: 
!http://www.unidata.ucar.edu/packages/netcdf/f90/Documentation/f90-html-docs/
!
!
!

  use typeSizes
  use netcdf
  implicit none

  character (len=*), parameter :: fileName_day = 'out_day.nc'
  character (len=*), parameter :: fileName_month = 'out_month.nc'
  character (len=*), parameter :: fileName_year = 'out_year.nc'
  character (len=15) :: fileName ,period_type

  integer,parameter ::closedID=-999     !flag for showing that a file is closed
  integer,save :: ncFileID_day=closedID  
  integer,save :: ncFileID_month=closedID
  integer,save :: ncFileID_year=closedID

  public :: InitnetCDF
  public :: Out_netCDF
  public :: CloseNetCDF

  private :: CreatenetCDFfile
  private :: createnewvariable
  private :: check
  private :: secondssince1970

contains
!_______________________________________________________________________

subroutine InitnetCDF

fileName = fileName_year
period_type = 'yearly'
call CreatenetCDFfile
fileName = fileName_month
period_type = 'monthly'
call CreatenetCDFfile
fileName = fileName_day
period_type = 'daily'
call CreatenetCDFfile

end subroutine InitnetCDF

subroutine CreatenetCDFfile
  ! Create the netCDF file

use GridValues_ml,         only : GRIDWIDTH_M
use Par_ml,                only : GIMAX,GJMAX
use ModelConstants_ml,     only : KMAX_MID   
use GridValues_ml,         only : coordzero,sigma_mid
use My_Derived_ml,         only : model
  implicit none

character (len=*), parameter :: version='Unimod rv1.4'       
character (len=*), parameter :: author_of_run='Unimod group' 
character (len=*), parameter :: projection='Stereographic'
character (len=*), parameter :: projection_params='60.0 -32.0 1.0'
character (len=*), parameter :: vert_coord='vertical coordinates = (p-p(top))/(p(surf)-p(top))'

real (kind = FourByteReal) :: xcoord(GIMAX),ycoord(GJMAX),kcoord(KMAX_MID)

character*8 ::created_date,lastmodified_date
character*10 ::created_hour,lastmodified_hour
integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,VarID,iVarID,jVarID,kVarID,i,j,k
real :: izero,jzero

  ! fileName: Name of the new created file 
  ! nf90_clobber: protect existing datasets
  ! ncFileID: netcdf ID


  call check(nf90_create(path = fileName, cmode = nf90_clobber, ncid = ncFileID))

  ! Define the dimensions
  call check(nf90_def_dim(ncid = ncFileID, name = "i", len = GIMAX, dimid = iDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "j", len = GJMAX, dimid = jDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "k", len = KMAX_MID, dimid = kDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "time", len = nf90_unlimited, dimid = timeDimID))

  call Date_And_Time(date=created_date,time=created_hour)
     print *, 'created_date: ',created_date
     print *, 'created_hour: ',created_hour

  ! Write global attributes
  call check(nf90_put_att(ncFileID, nf90_global, "Conventions", "GDV" ))
  call check(nf90_put_att(ncFileID, nf90_global, "version", version ))
  call check(nf90_put_att(ncFileID, nf90_global, "model", model))
  call check(nf90_put_att(ncFileID, nf90_global, "author_of_run", author_of_run))
  call check(nf90_put_att(ncFileID, nf90_global, "created_date", created_date))
  call check(nf90_put_att(ncFileID, nf90_global, "created_hour", created_hour))
  lastmodified_date = created_date
  lastmodified_hour = created_hour
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_date", lastmodified_date))
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour))

  call check(nf90_put_att(ncFileID, nf90_global, "projection",projection))
  call check(nf90_put_att(ncFileID, nf90_global, "projection_params",projection_params))
  call check(nf90_put_att(ncFileID, nf90_global, "vert_coord", vert_coord))
  call check(nf90_put_att(ncFileID, nf90_global, "period_type", period_type))

! define coordinate variables
  call check(nf90_def_var(ncFileID, "i", nf90_float, dimids = iDimID, varID = iVarID) )
  call check(nf90_put_att(ncFileID, iVarID, "coord_axis", "x"))
  call check(nf90_put_att(ncFileID, iVarID, "long_name", "EMEP grid x coordinate"))
  call check(nf90_put_att(ncFileID, iVarID, "units", "km"))

  call check(nf90_def_var(ncFileID, "j", nf90_float, dimids = jDimID, varID = jVarID) )
  call check(nf90_put_att(ncFileID, jVarID, "coord_axis", "y"))
  call check(nf90_put_att(ncFileID, jVarID, "long_name", "EMEP grid y coordinate"))
  call check(nf90_put_att(ncFileID, jVarID, "units", "km"))

  call check(nf90_def_var(ncFileID, "k", nf90_float, dimids = kDimID, varID = kVarID) )
  call check(nf90_put_att(ncFileID, kVarID, "coord_alias", "level"))
  call check(nf90_put_att(ncFileID, kVarID, "long_name", "vertical eta coordinates"))
  call check(nf90_put_att(ncFileID, kVarID, "units", "eta_level"))
  call check(nf90_put_att(ncFileID, kVarID, "positive", "down"))

  call check(nf90_def_var(ncFileID, "time", nf90_int, dimids = timeDimID, varID = VarID) )
  call check(nf90_put_att(ncFileID, VarID, "long_name", "current time at end of period"))
  call check(nf90_put_att(ncFileID, VarID, "units", "seconds since 1970-1-1 00:00:00.0 +00:00"))



  ! Leave define mode
  call check(nf90_enddef(ncFileID))

  call check(nf90_open(path = fileName, mode = nf90_write, ncid = ncFileID))

  call coordzero(izero,jzero)
!  print *, 'izero,jzero ',izero,jzero
  xcoord(1)=-izero*GRIDWIDTH_M/1000.
  do i=2,GIMAX
     xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
!     print *, i,xcoord(i)
  enddo
  call check(nf90_put_var(ncFileID, iVarID, xcoord) )
  ycoord(1)=-jzero*GRIDWIDTH_M/1000.
  do j=2,GJMAX
     ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
  enddo
  call check(nf90_put_var(ncFileID, jVarID, ycoord) )
  do k=1,KMAX_MID
     kcoord(k)=sigma_mid(k)
  enddo
  call check(nf90_put_var(ncFileID, kVarID, kcoord) )

  call check(nf90_close(ncFileID))

end subroutine CreatenetCDFfile

!_______________________________________________________________________

subroutine Out_netCDF(iotyp,def1,ndim,ident,kmax,icmp,dat,dim,scale)

use Par_ml,                only : NPROC,me,GIMAX,GJMAX,tgi0,tgj0,tlimax,tljmax,MAXLIMAX, MAXLJMAX
use ModelConstants_ml,     only : KMAX_MID, current_date  
use My_Derived_ml, only : NDDEP, NWDEP, NDERIV_2D, NDERIV_3D ,Deriv &
                         ,IOU_INST, IOU_YEAR ,IOU_MON, IOU_DAY  


  implicit none

integer ,intent(in) :: icmp,ndim,kmax,dim
type(Deriv),     intent(in) :: def1 ! definition of fields
integer,                         intent(in) :: iotyp
integer, dimension(:) ,intent(in) ::  ident
real    ,intent(in) :: scale 
!real, dimension(:,:,:,:), intent(in) :: dat ! Data arrays
real, dimension(dim,MAXLIMAX,MAXLJMAX,KMAX), intent(in) :: dat ! Data arrays


character*18 :: varname
character*8 ::lastmodified_date
character*10 ::lastmodified_hour,lastmodified_hour0,created_hour
integer :: varID,new,nrecords,ncFileID
integer :: nyear,nmonth,nday,nhour,ndate(4)
integer :: info,d,alloc_err,ijk,itag,status,i,j,k,nseconds
!real*4 :: buff 
real :: buff(MAXLIMAX*MAXLJMAX*KMAX_MID) 
real (kind = FourByteReal), allocatable,dimension(:,:,:)  :: data3D

!rv1.4.15 changed:
!==================================================================
  return     !TEMP - to avoid slowdown - suggested by pw, 26 Feb 2003.
!==================================================================

if(iotyp==IOU_YEAR)then
   fileName = fileName_year
   ncFileID = ncFileID_year
elseif(iotyp==IOU_MON)then
   fileName = fileName_month
   ncFileID = ncFileID_month
elseif(iotyp==IOU_DAY)then
   fileName = fileName_day
   ncFileID = ncFileID_day
!   return
else
   return
endif 


if(ndim /= 2 .and. ndim /= 3 )then
   print *, 'error in NetCDF_ml ',ndim
   call gc_abort(me,NPROC,"error in NetCDF_ml")
endif


!buffer the wanted part of data
ijk=0
do k=1,kmax
   do j = 1,tljmax(me)
      do i = 1,tlimax(me)
         ijk=ijk+1
         buff(ijk)=dat(icmp,i,j,k)*scale
      enddo
   enddo
enddo

!send all data to me=0
itag=icmp
if(me.eq.0)then
   
   !allocate a large array (only on one processor)
   allocate(data3D(GIMAX,GJMAX,kmax), stat=alloc_err)
   if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "alloc failed in NetCDF_ml")
   
   !write own data in global array
   ijk=0
   do k=1,kmax
      do j = tgj0(me),tgj0(me)+tljmax(me)-1
         do i = tgi0(me),tgi0(me)+tlimax(me)-1
            ijk=ijk+1
            data3D(i,j,k)=buff(ijk)
         enddo
      enddo
   enddo

   do d = 1, NPROC-1
      call gc_rrecv(itag,tlimax(d)*tljmax(d)*kmax, d, info,buff, buff)

      !copy data to global buffer
      ijk=0
      do k=1,kmax
         do j = tgj0(d),tgj0(d)+tljmax(d)-1
            do i = tgi0(d),tgi0(d)+tlimax(d)-1
               ijk=ijk+1
               data3D(i,j,k)=buff(ijk)
            enddo
         enddo
      enddo
      
   enddo
else
   call gc_rsend(itag,tlimax(me)*tljmax(me)*kmax, 0, info, buff, buff)
endif


if(me==0)then

   ndate(1)  = current_date%year
   ndate(2)  = current_date%month
   ndate(3)  = current_date%day
   ndate(4)  = current_date%hour

!test if the file is already open
   if(ncFileID==closedID)then
!open an existing netcdf dataset
      call check(nf90_open(path = fileName, mode = nf90_write, ncid = ncFileID))
      if(iotyp==IOU_YEAR)then
         ncFileID_year = ncFileID
      elseif(iotyp==IOU_MON)then
         ncFileID_month = ncFileID
      elseif(iotyp==IOU_DAY)then
         ncFileID_day = ncFileID
      endif
   endif

!make variable name
  if(ndim==2)then
     write(varname,fmt='(A,''_2D'')')trim(def1%name)
  elseif(ndim==3)then
     write(varname,fmt='(A,''_3D'')')trim(def1%name)
  else
     write(varname,fmt='(A,''_XX'')')trim(def1%name)
  endif
  !test first if the variable is already defined:
  status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)
  
  if(status == nf90_noerr) then     
!     print *, 'variable exists: ',varname
  else
     print *, 'creating variable: ',varname!,nf90_strerror(status)
     call  createnewvariable(ncFileID,varname,ndim,ident,ndate,def1)
  endif


  !get variable id
  call check(nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID))


  !find the number of records already written
  call check(nf90_get_att(ncFileID, VarID, "numberofrecords",   nrecords))
  !increase the last coordinate by one, to define position of new data
  nrecords=nrecords+1
!  print *,'number of dataset saved: ',nrecords


 !append new values
  if(ndim==3)then
  do k=1,kmax
  call check(nf90_put_var(ncFileID, VarID, data3D(:, :, k), start = (/ 1, 1, k,nrecords /)) )
  enddo
  else 
  call check(nf90_put_var(ncFileID, VarID, data3D(:, :, 1), start = (/ 1, 1, nrecords /)) )
  endif

!  if(icmp == dim)then
     deallocate(data3D, stat=alloc_err)
     if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "dealloc failed in NetCDF_ml")
!  endif



  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour0  ))
  call check(nf90_get_att(ncFileID, nf90_global, "created_hour", created_hour  ))
  call Date_And_Time(date=lastmodified_date,time=lastmodified_hour)
!    print *, 'date now: ',lastmodified_hour,' date before ',lastmodified_hour0,' date start ', created_hour

!write or change attributes NB: strings must be of same length as originally

  call check(nf90_put_att(ncFileID, VarID, "numberofrecords",   nrecords))

!update dates
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_date", lastmodified_date))
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour))
  call check(nf90_put_att(ncFileID, VarID, "current_date_last",ndate ))

  !get variable id
  call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = VarID))
  call secondssince1970(ndate,nseconds)
  call check(nf90_put_var(ncFileID, VarID, nseconds, start = (/nrecords/) ) )

!!close file
!  call check(nf90_close(ncFileID))

endif !me=0

return
end subroutine Out_netCDF

!_______________________________________________________________________


subroutine  createnewvariable(ncFileID,varname,ndim,ident,ndate,def1)

  !create new netCDF variable

use My_Derived_ml, only : Deriv 

  implicit none

  type(Deriv),     intent(in) :: def1 ! definition of fields
  character (len = *),intent(in) ::varname
  integer ,intent(in) ::ndim,ncFileID
  integer, dimension(:) ,intent(in) ::  ident,ndate

  integer :: iDimID,jDimID,kDimID,timeDimID
  integer :: varID,nrecords
  real (kind = FourByteReal) :: FillValue,scale


     call check(nf90_redef(ncid = ncFileID))

     !get dimensions id
     call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

     !define new variable
     if(ndim==3)then
        call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = nf90_float,     &
             dimids = (/ iDimID, jDimID, kDimID , timeDimID/), varID=varID ) )
     elseif(ndim==2)then
        call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = nf90_float,     &
             dimids = (/ iDimID, jDimID , timeDimID/), varID=varID ) )
     else
         print *, 'createnewvariable: unexpected ndim ',ndim   
     endif
     FillValue=0.
     scale=1.
     !define attributes of new variable
     call check(nf90_put_att(ncFileID, varID, "long_name",  def1%name ))
     nrecords=0
     call check(nf90_put_att(ncFileID, varID, "numberofrecords", nrecords))

     call check(nf90_put_att(ncFileID, varID, "units",   def1%unit))
     call check(nf90_put_att(ncFileID, varID, "class",   def1%class))
     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))
     call check(nf90_put_att(ncFileID, varID, "_FillValue", FillValue  ))

!     call check(nf90_put_att(ncFileID, varID, "periodlength",   "yearly"))

     call check(nf90_put_att(ncFileID, varID, "xfelt_ident",ident ))
     call check(nf90_put_att(ncFileID, varID, "current_date_first",ndate ))
     call check(nf90_put_att(ncFileID, varID, "current_date_last",ndate ))
   
     call check(nf90_enddef(ncid = ncFileID))

end subroutine  createnewvariable
!_______________________________________________________________________

  subroutine check(status)
    use Par_ml,                only : me,NPROC
    implicit none
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      call gc_abort(me,NPROC,"error in NetCDF_ml")
    end if
  end subroutine check  

  subroutine CloseNetCDF
!close open files
!NB the data in a NetCDF file is not "safe" before the file
!is closed. The files are NOT automatically properly 
!closed after end of program, and data may be lost if the files are not
!closed explicitely.

integer :: ncFileID

    if(ncFileID_year/=closedID)then
       ncFileID = ncFileID_year
       call check(nf90_close(ncFileID))
       ncFileID_year=closedID
    endif
    if(ncFileID_month/=closedID)then
       ncFileID = ncFileID_month
       call check(nf90_close(ncFileID))
       ncFileID_month=closedID
    endif
    if(ncFileID_day/=closedID)then
       ncFileID = ncFileID_day
       call check(nf90_close(ncFileID))
       ncFileID_day=closedID
    endif

  end subroutine CloseNetCDF

  subroutine secondssince1970(ndate,nseconds)
    !calculate how many seconds have passed since the start of the year

    use Dates_ml, only: nmdays,is_leap 
    implicit none

    integer, intent(in) :: ndate(4)
    integer, intent(out) :: nseconds
    integer :: n,nday
    nday=0
    do n=1,ndate(2)-1
       nday=nday+nmdays(n)
    enddo
    nday=nday+ndate(3)

    nseconds=3600*(ndate(4)+24*(nday-1))

!add seconds from each year since 1970
    do n=1970,ndate(1)-1
       nseconds=nseconds+24*3600*365+24*3600*is_leap(n)
    enddo

  end subroutine secondssince1970

end module NetCDF_ml
