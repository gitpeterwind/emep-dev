
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
!To improve: When output is onto the same file, but with different positions for the 
!lower left corner, the coordinates i_EMEP j_EMEP and long lat will be wrong
!

  use typeSizes
  use netcdf
  implicit none

  character (len=25), save :: fileName_inst = 'out_inst.nc'
  character (len=25), save :: fileName_hour = 'out_hour.nc'
  character (len=25), save :: fileName_day = 'out_day.nc'
  character (len=25), save :: fileName_month = 'out_month.nc'
  character (len=25), save :: fileName_year = 'out_year.nc'
  character (len=25) :: fileName ,period_type

  integer,parameter ::closedID=-999     !flag for showing that a file is closed
  integer,save :: ncFileID_inst=closedID  
  integer,save :: ncFileID_hour=closedID  
  integer,save :: ncFileID_day=closedID  
  integer,save :: ncFileID_month=closedID
  integer,save :: ncFileID_year=closedID
  integer, public, parameter :: Int1=1,Int2=2,Int4=3,Real4=4,Real8=5 !CDF typr for output


  public :: InitnetCDF
  public :: Out_netCDF
  public :: CloseNetCDF
  public :: Init_new_netCDF

  private :: CreatenetCDFfile
  private :: createnewvariable
  private :: check
  private :: secondssince1970

contains
!_______________________________________________________________________

subroutine InitnetCDF

use Par_ml,           only : GIMAX,GJMAX,ISMBEG,JSMBEG
use ModelConstants_ml,only : KMAX_MID   
use My_Outputs_ml,    only : NHOURLY_OUT, &      ! No. outputs
                             Asc2D, hr_out      ! Required outputs

integer :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
integer :: ih

write(*,*)'initnetcdf'
fileName = fileName_year
period_type = 'yearly'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)
fileName = fileName_month
period_type = 'monthly'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)
fileName = fileName_day
period_type = 'daily'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)

fileName = fileName_hour
period_type = 'hourly'
ISMBEGcdf=GIMAX+ISMBEG-1; JSMBEGcdf=GJMAX+JSMBEG-1
GIMAXcdf=0; GJMAXcdf=0
KMAXcdf=1
do ih=1,NHOURLY_OUT
   ISMBEGcdf=min(ISMBEGcdf,hr_out(ih)%ix1)
   JSMBEGcdf=min(JSMBEGcdf,hr_out(ih)%iy1)
   GIMAXcdf=max(GIMAXcdf,hr_out(ih)%ix2-hr_out(ih)%ix1+1)
   GJMAXcdf=max(GJMAXcdf,hr_out(ih)%iy2-hr_out(ih)%iy1+1)
   KMAXcdf =max(KMAXcdf,hr_out(ih)%nk)
enddo
GIMAXcdf=min(GIMAXcdf,GIMAX)
GJMAXcdf=min(GJMAXcdf,GJMAX)

!write(*,*)GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
call CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf)

fileName = fileName_inst
period_type = 'instant'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)

end subroutine InitnetCDF

subroutine Init_new_netCDF(fileName,iotyp) 

use Par_ml,           only : GIMAX,GJMAX,ISMBEG,JSMBEG
use ModelConstants_ml,only : KMAX_MID   
use My_Outputs_ml,    only : NHOURLY_OUT, &      ! No. outputs
                             Asc2D, hr_out      ! Required outputs
use My_Derived_ml,    only :IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY  

integer,  intent(in) :: iotyp
  character(len=*),  intent(in)  :: fileName 

integer :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
integer :: ih

!write(*,*)'Init_new_netCDF ',fileName,iotyp
call CloseNetCDF

if(iotyp==IOU_YEAR)then

fileName_year = trim(fileName)
period_type = 'yearly'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)

elseif(iotyp==IOU_MON)then

fileName_month = trim(fileName)
period_type = 'monthly'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)

elseif(iotyp==IOU_DAY)then

fileName_day = trim(fileName)
period_type = 'daily'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)

elseif(iotyp==IOU_HOUR)then

fileName_hour = trim(fileName)
period_type = 'hourly'
ISMBEGcdf=GIMAX+ISMBEG-1; JSMBEGcdf=GJMAX+JSMBEG-1
GIMAXcdf=0; GJMAXcdf=0
KMAXcdf=1
do ih=1,NHOURLY_OUT
   ISMBEGcdf=min(ISMBEGcdf,hr_out(ih)%ix1)
   JSMBEGcdf=min(JSMBEGcdf,hr_out(ih)%iy1)
   GIMAXcdf=max(GIMAXcdf,hr_out(ih)%ix2-hr_out(ih)%ix1+1)
   GJMAXcdf=max(GJMAXcdf,hr_out(ih)%iy2-hr_out(ih)%iy1+1)
   KMAXcdf =max(KMAXcdf,hr_out(ih)%nk)
enddo
GIMAXcdf=min(GIMAXcdf,GIMAX)
GJMAXcdf=min(GJMAXcdf,GJMAX)
!write(*,*)'sizes CDF ',GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
call CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf)

elseif(iotyp==IOU_INST)then

fileName_inst = trim(fileName)
period_type = 'instant'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)

else
period_type = 'unknown'
call CreatenetCDFfile(fileName,GIMAX,GJMAX,ISMBEG,JSMBEG,KMAX_MID)
endif

end subroutine Init_new_netCDF

subroutine CreatenetCDFfile(fileName,GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf)
  ! Create the netCDF file

use GridValues_ml,         only : GRIDWIDTH_M,fi,xp,yp,xp_EMEP_official&
                                  ,yp_EMEP_official,fi_EMEP,GRIDWIDTH_M_EMEP&
                                  ,GlobalPosition,gb_glob,gl_glob
use Par_ml,                only : GIMAX,GJMAX,ISMBEG,JSMBEG,IILARDOM,JJLARDOM
use ModelConstants_ml,     only : KMAX_MID, runlabel1, runlabel2  
use GridValues_ml,         only : coordzero,sigma_mid
use My_Derived_ml,         only : model
use PhysicalConstants_ml,  only : PI       
   implicit none
integer, intent(in) :: GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
character(len=*),  intent(in)  :: fileName 

character (len=*), parameter :: version='Unimod rv1.9'       
character (len=*), parameter :: author_of_run='Unimod group' 
character (len=*), parameter :: projection='Stereographic'
character (len=*), parameter :: vert_coord='vertical coordinates = (p-p(top))/(p(surf)-p(top))'
character (len=19) :: projection_params='90.0 -32.0 0.933013' !set later on

real (kind = FourByteReal) :: xcoord(GIMAX),ycoord(GJMAX),kcoord(KMAX_MID)

character*8 ::created_date,lastmodified_date
character*10 ::created_hour,lastmodified_hour
integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,VarID,iVarID,jVarID,kVarID,i,j,k
integer :: iEMEPVarID,jEMEPVarID,latVarID,longVarID
real :: izero,jzero

  ! fileName: Name of the new created file 
  ! nf90_clobber: protect existing datasets
  ! ncFileID: netcdf ID

write(*,*)'create ',trim(fileName)
write(*,*)'with sizes (IMAX,JMAX,IBEG,JBEG,KMAX) ',GIMAXcdf,GJMAXcdf,ISMBEGcdf,JSMBEGcdf,KMAXcdf
  call check(nf90_create(path = trim(fileName), cmode = nf90_clobber, ncid = ncFileID))

  ! Define the dimensions
  call check(nf90_def_dim(ncid = ncFileID, name = "i", len = GIMAXcdf, dimid = iDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "j", len = GJMAXcdf, dimid = jDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "k", len = KMAXcdf, dimid = kDimID))
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
  write(projection_params,fmt='(''90.0 '',F5.1,F9.6)')fi,(1.+sin(PI/3.0))/2.
  call check(nf90_put_att(ncFileID, nf90_global, "projection_params",projection_params))
  call check(nf90_put_att(ncFileID, nf90_global, "vert_coord", vert_coord))
  call check(nf90_put_att(ncFileID, nf90_global, "period_type", trim(period_type)))
  call check(nf90_put_att(ncFileID, nf90_global, "run_label", trim(runlabel2)))

! define coordinate variables
  call check(nf90_def_var(ncFileID, "i", nf90_float, dimids = iDimID, varID = iVarID) )
  call check(nf90_put_att(ncFileID, iVarID, "coord_axis", "x"))
  call check(nf90_put_att(ncFileID, iVarID, "long_name", "EMEP grid x coordinate"))
  call check(nf90_put_att(ncFileID, iVarID, "units", "km"))

  call check(nf90_def_var(ncFileID, "i_EMEP", nf90_float, dimids = iDimID, varID = iEMEPVarID) )
  call check(nf90_put_att(ncFileID, iEMEPVarID, "long_name", "official EMEP grid coordinate i"))
  call check(nf90_put_att(ncFileID, iEMEPVarID, "units", "gridcells"))

  call check(nf90_def_var(ncFileID, "j", nf90_float, dimids = jDimID, varID = jVarID) )
  call check(nf90_put_att(ncFileID, jVarID, "coord_axis", "y"))
  call check(nf90_put_att(ncFileID, jVarID, "long_name", "EMEP grid y coordinate"))
  call check(nf90_put_att(ncFileID, jVarID, "units", "km"))

  call check(nf90_def_var(ncFileID, "j_EMEP", nf90_float, dimids = jDimID, varID = jEMEPVarID) )
  call check(nf90_put_att(ncFileID, jEMEPVarID, "long_name", "official EMEP grid coordinate j"))
  call check(nf90_put_att(ncFileID, jEMEPVarID, "units", "gridcells"))

  call check(nf90_def_var(ncFileID, "lat", nf90_float, dimids = (/ iDimID, jDimID/), varID = latVarID) )
  call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"))
  call check(nf90_put_att(ncFileID, latVarID, "units", "degrees"))

  call check(nf90_def_var(ncFileID, "long", nf90_float, dimids = (/ iDimID, jDimID/), varID = longVarID) )
  call check(nf90_put_att(ncFileID, longVarID, "long_name", "longitude"))
  call check(nf90_put_att(ncFileID, longVarID, "units", "degrees"))


  call check(nf90_def_var(ncFileID, "k", nf90_float, dimids = kDimID, varID = kVarID) )
  call check(nf90_put_att(ncFileID, kVarID, "coord_alias", "level"))
  call check(nf90_put_att(ncFileID, kVarID, "long_name", "vertical eta coordinates"))
  call check(nf90_put_att(ncFileID, kVarID, "units", "eta_level"))
  call check(nf90_put_att(ncFileID, kVarID, "positive", "down"))

  call check(nf90_def_var(ncFileID, "time", nf90_int, dimids = timeDimID, varID = VarID) )
  call check(nf90_put_att(ncFileID, VarID, "long_name", "time at middle of period"))
  call check(nf90_put_att(ncFileID, VarID, "units", "seconds since 1970-1-1 00:00:00.0 +00:00"))


  ! Leave define mode
  call check(nf90_enddef(ncFileID))

  call check(nf90_open(path = trim(fileName), mode = nf90_write, ncid = ncFileID))

! Define horizontal distances in GDV conventions

  xcoord(1)=(ISMBEGcdf-xp)*GRIDWIDTH_M/1000.
  do i=2,GIMAXcdf
     xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
!     print *, i,xcoord(i)
  enddo
  call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAXcdf)) )

  ycoord(1)=(JSMBEGcdf-yp)*GRIDWIDTH_M/1000.
  do j=2,GJMAXcdf
     ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
  enddo
  call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAXcdf)) )

! Define horizontal coordinates in the official EMEP grid
!  xp_EMEP_official=8.
!  yp_EMEP_official=110.
!  GRIDWIDTH_M_EMEP=50000.
!  fi_EMEP=-32.
  if(fi==fi_EMEP)then
! Implemented only if fi = fi_EMEP = -32 (Otherwise needs a 2-dimensional mapping)
! uses (i-xp)*GRIDWIDTH_M = (i_EMEP-xp_EMEP)*GRIDWIDTH_M_EMEP
  do i=1,GIMAXcdf
     xcoord(i)=(i+ISMBEGcdf-1-xp)*GRIDWIDTH_M/GRIDWIDTH_M_EMEP + xp_EMEP_official
!     print *, i,xcoord(i)
  enddo
  do j=1,GJMAXcdf
     ycoord(j)=(j+JSMBEGcdf-1-yp)*GRIDWIDTH_M/GRIDWIDTH_M_EMEP + yp_EMEP_official
!     print *, j,ycoord(j)
  enddo
  else
  do i=1,GIMAXcdf
     xcoord(i)=NF90_FILL_FLOAT
  enddo
  do j=1,GJMAXcdf
     ycoord(j)=NF90_FILL_FLOAT
  enddo
  endif
  call check(nf90_put_var(ncFileID, iEMEPVarID, xcoord(1:GIMAXcdf)) )
  call check(nf90_put_var(ncFileID, jEMEPVarID, ycoord(1:GJMAXcdf)) )

!Define longitude and latitude
  call GlobalPosition
  if(ISMBEGcdf+GIMAXcdf-1<=IILARDOM .and. JSMBEGcdf+GJMAXcdf-1<=JJLARDOM)then
  call check(nf90_put_var(ncFileID, latVarID, gb_glob(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)) )
  call check(nf90_put_var(ncFileID, longVarID, gl_glob(ISMBEGcdf:ISMBEGcdf+GIMAXcdf-1,JSMBEGcdf:JSMBEGcdf+GJMAXcdf-1)) )
  endif
!Define vertical levels
  if(KMAXcdf==KMAX_MID)then
     do k=1,KMAX_MID
        kcoord(k)=sigma_mid(k)
     enddo
  else
     do k=1,KMAXcdf
        kcoord(k)=sigma_mid(KMAX_MID-k+1) !REVERSE order of k !
     enddo
  endif
  call check(nf90_put_var(ncFileID, kVarID, kcoord(1:KMAXcdf)) )

  call check(nf90_close(ncFileID))

end subroutine CreatenetCDFfile

!_______________________________________________________________________

subroutine Out_netCDF(iotyp,def1,ndim,ident,kmax,icmp,dat,dim,scale,CDFtype,ist,jst,ien,jen,ik)

use Par_ml,                only : NPROC,me,GIMAX,GJMAX,tgi0,tgj0,tlimax,tljmax,MAXLIMAX, MAXLJMAX
use ModelConstants_ml,     only : KMAX_MID, current_date  
use My_Derived_ml, only : NDDEP, NWDEP, NDERIV_2D, NDERIV_3D ,Deriv &
                         ,IOU_INST,IOU_HOUR, IOU_YEAR ,IOU_MON, IOU_DAY  
use Dates_ml, only: nmdays,is_leap 


  implicit none

integer ,intent(in) :: icmp,ndim,kmax,dim
type(Deriv),     intent(in) :: def1 ! definition of fields
integer,                         intent(in) :: iotyp
integer, dimension(:) ,intent(in) ::  ident
real    ,intent(in) :: scale 
!real, dimension(:,:,:,:), intent(in) :: dat ! Data arrays
real, dimension(dim,MAXLIMAX,MAXLJMAX,KMAX), intent(in) :: dat ! Data arrays
integer, optional, intent(in) :: ist,jst,ien,jen,ik !start and end of saved area. Only level ik is written if defined
integer, optional, intent(in) :: CDFtype != OUTtype. output type (Integer*1, Integer*2,Integer*4, real*8 or real*4) 

character*18 :: varname
character*8 ::lastmodified_date
character*10 ::lastmodified_hour,lastmodified_hour0,created_hour
integer :: varID,new,nrecords,ncFileID
integer :: nyear,nmonth,nday,nhour,ndate(4)
integer :: info,d,alloc_err,ijk,itag,status,i,j,k,nseconds
integer :: i1,i2,j1,j2
!real*4 :: buff 
real :: buff(MAXLIMAX*MAXLJMAX*KMAX_MID) 
real*8 , allocatable,dimension(:,:,:)  :: Rdata3D
integer*4, allocatable,dimension(:,:,:)  :: Idata3D
integer, save ::SavedVar(KMAX_MID)=-999,SavedncFile(KMAX_MID)=-999
integer :: OUTtype !local version of CDFtype

!rv1.4.15 changed:
!==================================================================
!  return     !TEMP - to avoid slowdown - suggested by pw, 26 Feb 2003.
!==================================================================
!make variable name
  if(ndim==2)then
!     write(varname,fmt='(A,''_2D'')')trim(def1%name)
     write(varname,fmt='(A)')trim(def1%name)
  elseif(ndim==3)then
     write(varname,fmt='(A,''_3D'')')trim(def1%name)
  else
     write(varname,fmt='(A,''_XX'')')trim(def1%name)
  endif

!to shorten the output we can save only the components explicitely named here
!if(varname.ne.'D2_NO2'.and.varname.ne.'D2_O3' &
!                         .and.varname.ne.'D2_PM10')return

!do not write 3D fields (in order to shorten outputs)
!if(ndim==3)return

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
elseif(iotyp==IOU_HOUR)then
   fileName = fileName_hour
   ncFileID = ncFileID_hour
elseif(iotyp==IOU_INST)then
   fileName = fileName_inst
   ncFileID = ncFileID_inst
else
   return
endif 

if(ndim /= 2 .and. ndim /= 3 )then
   print *, 'error in NetCDF_ml ',ndim
   call gc_abort(me,NPROC,"error in NetCDF_ml")
endif

OUTtype=Real4  !default value
if(present(CDFtype))OUTtype=CDFtype

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
   if(OUTtype==Int1 .or. OUTtype==Int2 .or. OUTtype==Int4)then
      allocate(Idata3D(GIMAX,GJMAX,kmax), stat=alloc_err)
      if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "alloc failed in NetCDF_ml")
   elseif(OUTtype==Real4 .or. OUTtype==Real8)then
      allocate(Rdata3D(GIMAX,GJMAX,kmax), stat=alloc_err)
      if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "alloc failed in NetCDF_ml")
   else
      WRITE(*,*)'WARNING NetCDF:Data type not supported'
   endif

   !write own data in global array
   if(OUTtype==Int1 .or. OUTtype==Int2 .or. OUTtype==Int4)then
   ijk=0
   do k=1,kmax
      do j = tgj0(me),tgj0(me)+tljmax(me)-1
         do i = tgi0(me),tgi0(me)+tlimax(me)-1
            ijk=ijk+1
            Idata3D(i,j,k)=buff(ijk)
         enddo
      enddo
   enddo
   else
   ijk=0
   do k=1,kmax
      do j = tgj0(me),tgj0(me)+tljmax(me)-1
         do i = tgi0(me),tgi0(me)+tlimax(me)-1
            ijk=ijk+1
            Rdata3D(i,j,k)=buff(ijk)
         enddo
      enddo
   enddo
   endif

   do d = 1, NPROC-1
      call gc_rrecv(itag,tlimax(d)*tljmax(d)*kmax, d, info,buff, buff)

      !copy data to global buffer
      if(OUTtype==Int1 .or. OUTtype==Int2 .or. OUTtype==Int4)then
      ijk=0
      do k=1,kmax
         do j = tgj0(d),tgj0(d)+tljmax(d)-1
            do i = tgi0(d),tgi0(d)+tlimax(d)-1
               ijk=ijk+1
               Idata3D(i,j,k)=buff(ijk)
            enddo
         enddo
      enddo
      else
      ijk=0
      do k=1,kmax
         do j = tgj0(d),tgj0(d)+tljmax(d)-1
            do i = tgi0(d),tgi0(d)+tlimax(d)-1
               ijk=ijk+1
               Rdata3D(i,j,k)=buff(ijk)
            enddo
         enddo
      enddo
      endif
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
      call check(nf90_open(path = trim(fileName), mode = nf90_write, ncid = ncFileID))
      if(iotyp==IOU_YEAR)then
         ncFileID_year = ncFileID
      elseif(iotyp==IOU_MON)then
         ncFileID_month = ncFileID
      elseif(iotyp==IOU_DAY)then
         ncFileID_day = ncFileID
      elseif(iotyp==IOU_HOUR)then
         ncFileID_hour = ncFileID
      elseif(iotyp==IOU_INST)then
         ncFileID_inst = ncFileID
      endif
   endif

  !test first if the variable is already defined:
  status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)
  
  if(status == nf90_noerr) then     
!     print *, 'variable exists: ',varname
  else
     print *, 'creating variable: ',varname!,nf90_strerror(status)
     call  createnewvariable(ncFileID,varname,ndim,ident,ndate,def1,OUTtype)
  endif


  !get variable id
  call check(nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID))


  !find the number of records already written
  call check(nf90_get_att(ncFileID, VarID, "numberofrecords",   nrecords))
!  print *,'number of dataset saved: ',nrecords
  !increase the last coordinate by one, to define position of new data
  !test if new record is needed
  if(present(ik))then
     if((SavedVar(ik)==VarID.and.SavedncFile(ik)==ncFileID).or.nrecords==0)then
        nrecords=nrecords+1 !start a new record
        SavedVar(:)=-999; SavedncFile(:)=-999 !reset 
     endif
     !write which variable, file and level is written
     SavedVar(ik)=VarID;SavedncFile(ik)=ncFileID
  else
     nrecords=nrecords+1
  endif
!  print *,'writing on dataset: ',nrecords

  i1=1;i2=GIMAX;j1=1;j2=GJMAX  !start and end of saved area
  if(present(ist))i1=max(ist,i1)
  if(present(ien))i2=min(ien,i2)
  if(present(jst))j1=max(jst,j1)
  if(present(jen))j2=min(jen,j2)
 !append new values
  if(OUTtype==Int1 .or. OUTtype==Int2 .or. OUTtype==Int4)then

  if(ndim==3)then
     if(present(ik))then
!     print *, 'write: ',i1,i2, j1,j2,ik
        call check(nf90_put_var(ncFileID, VarID, &
        Idata3D(i1:i2, j1:j2, 1), start = (/ 1, 1, ik,nrecords /)) )
     else
        do k=1,kmax
           call check(nf90_put_var(ncFileID, VarID,&
           Idata3D(i1:i2, j1:j2, k), start = (/ 1, 1, k,nrecords /)) )
        enddo
     endif
  else 
     call check(nf90_put_var(ncFileID, VarID,&
     Idata3D(i1:i2, j1:j2, 1), start = (/ 1, 1, nrecords /)) )
  endif

!  if(icmp == dim)then
     deallocate(Idata3D, stat=alloc_err)
     if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "dealloc failed in NetCDF_ml")
!  endif

  else  
     !type Real
  if(ndim==3)then
     if(present(ik))then
!     print *, 'write: ',i1,i2, j1,j2,ik
        call check(nf90_put_var(ncFileID, VarID, &
        Rdata3D(i1:i2, j1:j2, 1), start = (/ 1, 1, ik,nrecords /)) )
     else
        do k=1,kmax
           call check(nf90_put_var(ncFileID, VarID,&
           Rdata3D(i1:i2, j1:j2, k), start = (/ 1, 1, k,nrecords /)) )
        enddo
     endif
  else 
     call check(nf90_put_var(ncFileID, VarID,&
     Rdata3D(i1:i2, j1:j2, 1), start = (/ 1, 1, nrecords /)) )
  endif

!  if(icmp == dim)then
     deallocate(Rdata3D, stat=alloc_err)
     if ( alloc_err /= 0 ) call gc_abort(me,NPROC, "dealloc failed in NetCDF_ml")
!  endif

  endif !type Real


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
  !middle of period: !NB WORKS ONLY FOR COMPLETE PERIODS
      if(iotyp==IOU_YEAR)then
         nseconds=nseconds-43200*365-43200*is_leap(ndate(1)-1)
      elseif(iotyp==IOU_MON)then
         nseconds=nseconds-43200*nmdays(max(ndate(2)-1,1))!nmdays(jan)=nmdays(dec)
      elseif(iotyp==IOU_DAY)then
         nseconds=nseconds-43200 !24*3600/2=43200
      elseif(iotyp==IOU_HOUR)then
         nseconds=nseconds-1800  !1800=half hour
      elseif(iotyp==IOU_INST)then
           nseconds=nseconds       
      endif

  call check(nf90_put_var(ncFileID, VarID, nseconds, start = (/nrecords/) ) )

!!close file
!  call check(nf90_close(ncFileID))

endif !me=0

return
end subroutine Out_netCDF

!_______________________________________________________________________


subroutine  createnewvariable(ncFileID,varname,ndim,ident,ndate,def1,OUTtype)

  !create new netCDF variable

use Par_ml,                only : NPROC,me
use My_Derived_ml, only : Deriv 

  implicit none

  type(Deriv),     intent(in) :: def1 ! definition of fields
  character (len = *),intent(in) ::varname
  integer ,intent(in) ::ndim,ncFileID,OUTtype
  integer, dimension(:) ,intent(in) ::  ident,ndate

  integer :: iDimID,jDimID,kDimID,timeDimID
  integer :: varID,nrecords
  real :: scale
  integer :: OUTtypeCDF !NetCDF code for type

  if(OUTtype==Int1)then
     OUTtypeCDF=nf90_byte
  elseif(OUTtype==Int2)then
     OUTtypeCDF=nf90_short
  elseif(OUTtype==Int4)then
     OUTtypeCDF=nf90_int
  elseif(OUTtype==Real4)then
     OUTtypeCDF=nf90_float
  elseif(OUTtype==Real8)then
     OUTtypeCDF=nf90_double
  else
     call gc_abort(me,NPROC,"NetCDF_ml: undefined datatype")
  endif

     call check(nf90_redef(ncid = ncFileID))

     !get dimensions id
     call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

     !define new variable
     if(ndim==3)then
        call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = OUTtypeCDF,     &
             dimids = (/ iDimID, jDimID, kDimID , timeDimID/), varID=varID ) )
     elseif(ndim==2)then
        call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = OUTtypeCDF,     &
             dimids = (/ iDimID, jDimID , timeDimID/), varID=varID ) )
     else
         print *, 'createnewvariable: unexpected ndim ',ndim   
     endif
!     FillValue=0.
     scale=1.
     !define attributes of new variable
     call check(nf90_put_att(ncFileID, varID, "long_name",  def1%name ))
     nrecords=0
     call check(nf90_put_att(ncFileID, varID, "numberofrecords", nrecords))

     call check(nf90_put_att(ncFileID, varID, "units",   def1%unit))
     call check(nf90_put_att(ncFileID, varID, "class",   def1%class))
     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))

  if(OUTtype==Int1)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_byte  ))
  elseif(OUTtype==Int2)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_short  ))
  elseif(OUTtype==Int4)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_int   ))
  elseif(OUTtype==Real4)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_float  ))
  elseif(OUTtype==Real8)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_double  ))
  endif
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
    if(ncFileID_hour/=closedID)then
       ncFileID = ncFileID_hour
       call check(nf90_close(ncFileID))
       ncFileID_hour=closedID
    endif
    if(ncFileID_inst/=closedID)then
       ncFileID = ncFileID_inst
       call check(nf90_close(ncFileID))
       ncFileID_inst=closedID
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
