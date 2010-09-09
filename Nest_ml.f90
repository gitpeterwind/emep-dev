! <Nest_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007 met.no
!*
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************!

module Nest_ml
  !
  !This module performs the reading or writing of data for nested runs
  !
  !To make a nested run:
  !1) run with MODE=1 to write out 3d BC
  !2) run (in a smaller domain) with MODE=2


  !
  !Set MODE and istart,jstart,iend,jend
  !Choose NHOURSAVE and NHOURREAD
  !


  !Grids may have any projection.
  !Horizontal interpolation uses a weighted average of the four closest points
  !This will work also if points in the present grid are not covered by the external grid.


  !To do:
  !At present the vertical coordinates cannot be interpolated and must be the same in both grid.
  !It should be possible to save only xn_adv_bnd if the inner grid is known for the outer grid.
  !The routines should be thought together with GlobalBC_ml (can it replace it?)


  !Peter May 2006

  use OwnDataTypes_ml, only :Deriv
  use TimeDate_ml,       only : date
  use GridValues_ml,  only : gl,gb
  use ChemChemicals_ml , only :species
  use ChemSpecs_shl_ml , only :NSPEC_SHL
  use ChemSpecs_adv_ml
  use ChemSpecs_tot_ml , only :NSPEC_TOT
  use netcdf
  use netcdf_ml,      only : GetCDF,Out_netCDF,Init_new_netCDF,&
       secondssince1970,dayssince1900,Int1,Int2,Int4,Real4,Real8
  use Functions_ml,    only : great_circle_distance
  use ModelConstants_ml,    only : KMAX_MID, NPROC &
          , FORECAST, DEBUG_ICBC=>DEBUG_NEST_ICBC &
          , IOU_INST,IOU_HOUR, IOU_YEAR,IOU_MON, IOU_DAY,RUNDOMAIN
  use Par_ml   ,      only : MAXLIMAX, MAXLJMAX, GIMAX,GJMAX,IRUNBEG,JRUNBEG &
       , me, li0,li1,lj0,lj1,limax,ljmax, tgi0, tgj0, tlimax, tljmax
  use Chemfields_ml,  only : xn_adv, xn_shl    ! emep model concs.

  use NetCDF_ml, only :WriteCDF

  implicit none

  INCLUDE 'mpif.h'
  INTEGER INFO

  integer,parameter ::MODE=0   !0=donothing , 1=write , 2=read , 3=read and write
  !10=write at end of run, 11=read at start , 12=read at start and write at end (BIC)


! Nested input/output on FORECAST mode
  integer, public, parameter :: FORECAST_NDUMP = 1  ! Number of nested output
  ! on FORECAST mode (1: starnt next forecast; 2: NMC statistics)
  type(date), public :: outdate(FORECAST_NDUMP)=date(-1,-1,-1,-1,-1)
  ! Nested output dates on FORECAST mode
! IFS-MOZART BC
  type, private :: icbc                 ! Inital (IC) & Boundary Conditions (BC)
    character(len=24) :: varname=""
    logical           :: wanted=.false.,found=.false.
  end type icbc
  type(icbc), dimension(NSPEC_ADV), private :: &
    adv_ic=icbc('',.false.,.false.), &  ! Initial 3D IC/CB
    adv_bc=icbc('',.false.,.false.)     ! Time dependent BC
  type, private :: adv_icbc             ! IC/BC Set, included intended ixadv
    integer           :: ixadv=-1
    type(icbc)        :: icbc
  end type adv_icbc
  type(adv_icbc), dimension(9), private, parameter :: &  ! BC from IFS-MOZART
    FORECAST_BC=(/adv_icbc(IXADV_O3    ,icbc('O3_VMR_inst'    ,.true.,.false.)), &
                  adv_icbc(IXADV_NO    ,icbc('NO_VMR_inst'    ,.true.,.false.)), &
                  adv_icbc(IXADV_NO2   ,icbc('NO2_VMR_inst'   ,.true.,.false.)), &
                  adv_icbc(IXADV_PAN   ,icbc('PAN_VMR_inst'   ,.true.,.false.)), &
                  adv_icbc(IXADV_HNO3  ,icbc('HNO3_VMR_inst'  ,.true.,.false.)), &
                  adv_icbc(IXADV_CO    ,icbc('CO_VMR_inst'    ,.true.,.false.)), &
                  adv_icbc(IXADV_C2H6  ,icbc('C2H6_VMR_inst'  ,.true.,.false.)), &
                  adv_icbc(IXADV_HCHO  ,icbc('CH2O_VMR_inst'  ,.true.,.false.)), &
                  adv_icbc(IXADV_CH3CHO,icbc('CH3CHO_VMR_inst',.true.,.false.))/)

  !coordinates of subdomain to write
  !coordinates relative to LARGE domain (only used in write mode)
  integer ::istart=60,jstart=11,iend=107,jend=58 !ENEA NB: version has changed, these numbers where for small domain!!!

  !/-- subroutines

  public  :: readxn
  public  :: wrtxn


  !  logical, save, public::Nest_BC,Nest_3D

  integer, public, parameter :: NHOURSAVE=3 !time between two saves. should be a fraction of 24
  integer, public, parameter :: NHOURREAD=1 !time between two reads. should be a fraction of 24
  !if(NHOURREAD<NHOURSAVE) the data is interpolated in time

  private

  !Use TOP BC on forecast mode
  logical, parameter :: TOP_BC=.false..or.FORECAST
  integer,save :: iw, ie, js, jn, kt ! i West/East bnd; j North/South bnd; k Top
  !BC values at boundaries in present grid
  real, save, allocatable, dimension(:,:,:,:) :: &
    xn_adv_bndw, xn_adv_bnde, & ! west and east
    xn_adv_bnds, xn_adv_bndn, & ! north and south
    xn_adv_bndt                 ! top

  !dimension of external grid for BC
  integer,save :: Next_BC,KMAX_ext_BC

  integer,save :: itime!itime_saved(2),
  real(kind=8),save :: rtime_saved(2)
  character(len=30),save  :: filename_read_BC='EMEP_IN.nc'
  character(len=30),save  :: filename_read_3D='EMEP_IN.nc'
  character(len=30),save  :: filename_write='EMEP_OUT.nc'
  real(kind=8), parameter :: halfsecond=1.0/(24.0*3600.0)!used to avoid rounding errors

contains

subroutine readxn(indate)
  implicit none
  type(date), intent(in) :: indate           ! Gives year..seconds
  integer,save  :: first_data=-1

  integer :: n,i,j,k !nseconds(1),n1,II,JJ
  integer :: ndate(4) !nstart,nfetch,nseconds_indate
  real(kind=8):: ndays_indate

  !    real , dimension(48,48,20) ::data
  real :: W1,W2
  logical, save :: first_call=.true.
  logical :: fexist=.false.

  if(MODE /= 2.and.MODE /= 3.and. MODE /= 11.and. MODE /= 12.and. .not.FORECAST)return

  ndate(1)  = indate%year
  ndate(2)  = indate%month
  ndate(3)  = indate%day
  ndate(4)  = indate%hour
 !call secondssince1970(ndate,nseconds_indate)
  call dayssince1900(ndate,ndays_indate)

  if(FORECAST)then ! FORECAST mode superseeds nest MODE
    filename_read_3D='EMEP_IN_IC.nc'          !IC file: dump/re-start
    filename_read_BC='EMEP_IN_BC_YYYYMMDD.nc' !BC file: 01,...,24 UTC rec for 1 day
    if(first_call)then
      first_call=.false.
      inquire(file=filename_read_3D,exist=fexist)
      if(.not.fexist)then
        if(me==0) print *,'No nest IC file found: ',trim(filename_read_3D)
      else
        if(me==0) print *,'RESET ALL XN 3D'
        call reset_3D(ndays_indate)
      endif
    endif
    if(mod(indate%hour,NHOURREAD)/=0.or.indate%seconds/=0) return
    n=index(filename_read_BC,'YYYYMMDD')
    if(n>0) write(filename_read_BC(n:n+7),"(I4.4,2I2.2)")indate%year,indate%month,indate%day
    inquire(file=filename_read_BC,exist=fexist)
    if(.not.fexist)then
      if(me==0) print *,'No nest BC file found: ',trim(filename_read_BC)
      return
    endif
  elseif(MODE == 11.or.MODE == 12)then
    if(.not. first_call)return
    first_call=.false.
    if(me==0)   print *,'RESET ALL XN 3D'
    call reset_3D(ndays_indate)
    return
  else
   !if(me==0)   print *,'call to READXN',indate%hour,indate%seconds
    if(mod(indate%hour,NHOURREAD)/=0.or.indate%seconds/=0)return
  endif

!never comes to this point if MODE=11 or 12

  if(me==0) print *,'NESTING'

  if(first_data==-1)then
    if(.not.FORECAST) call reset_3D(ndays_indate)
    call read_newdata_LATERAL(ndays_indate,1)
    call read_newdata_LATERAL(ndays_indate,2)
  endif

  if(ndays_indate-rtime_saved(2)>halfsecond)then
    !look for a new data set
    if(me==0)write(*,*)'NEST: READING NEW BC DATA'
    call read_newdata_LATERAL(ndays_indate,2)
  endif

!   make weights for time interpolation
  W1=1.0;  W2=0.0!default
  if(ndays_indate-rtime_saved(1)>halfsecond)then
    !interpolate
    W2=(ndays_indate-rtime_saved(1))/(rtime_saved(2)-rtime_saved(1))
    W1=1.0-W2
!   if(me==1)then
!   call datefromdayssince1900(ndate,ndays_indate,1)
!   write(*,*)'interpolating between'
!   call datefromdayssince1900(ndate,rtime_saved(1),1)
!   write(*,*)'and'
!   call datefromdayssince1900(ndate,rtime_saved(2),1)
!   write(*,*)'with weights : ',W1,W2
!   endif
  endif
! if(me==0)write(*,*)'weights : ',W1,W2,rtime_saved(1),rtime_saved(2)

  forall (n=1:NSPEC_ADV, adv_bc(n)%wanted.and.adv_bc(n)%found)
    forall (i=iw:iw, k=1:KMAX_ext_BC, j=1:ljmax, i>=1) &
      xn_adv(n,i,j,k)=W1*xn_adv_bndw(n,j,k,1)+W2*xn_adv_bndw(n,j,k,2)
    forall (i=ie:ie, k=1:KMAX_ext_BC, j=1:ljmax, i<=limax) &
      xn_adv(n,i,j,k)=W1*xn_adv_bnde(n,j,k,1)+W2*xn_adv_bnde(n,j,k,2)
    forall (j=js:js, k=1:KMAX_ext_BC, i=1:limax, j>=1) &
      xn_adv(n,i,j,k)=W1*xn_adv_bnds(n,i,k,1)+W2*xn_adv_bnds(n,i,k,2)
    forall (j=jn:jn, k=1:KMAX_ext_BC, i=1:limax, j<=ljmax) &
      xn_adv(n,i,j,k)=W1*xn_adv_bndn(n,i,k,1)+W2*xn_adv_bndn(n,i,k,2)
    forall (k=kt:kt, i=1:limax, j=1:ljmax, k>=1) &
      xn_adv(n,i,j,k)=W1*xn_adv_bndt(n,i,j,1)+W2*xn_adv_bndt(n,i,j,2)
  end forall

  first_data=0
  return
end subroutine readxn

subroutine wrtxn(indate,WriteNow)
  implicit none
  type(date), intent(in) :: indate
  logical, intent(in) :: WriteNow !Do not check indate value
  real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: dat ! Data arrays

  type(Deriv) :: def1 ! definition of fields
  integer :: n,iotyp,ndim,kmax
  real :: scale
  logical, save ::first_call=.true.

  if(MODE /= 1.and.MODE /= 3.and.MODE /= 10.and.MODE /= 12.and. .not.FORECAST)return

  if(FORECAST)then ! FORECAST mode superseeds nest MODE
    outdate(:)%seconds=0   ! output only at full hours
    if(.not.any(indate%year   ==outdate%year  .and.   &
                indate%month  ==outdate%month .and.   &
                indate%day    ==outdate%day   .and.   &
                indate%hour   ==outdate%hour  .and.   &
                indate%seconds==outdate%seconds))return
    if(me==0) print "(A,I5.4,2('-',I2.2),I3.2,2(':',I2.2))",&
      " Forecast nest/dump at",                             &
      indate%year,indate%month,indate%day,                  &
      indate%hour,indate%seconds/60,mod(indate%seconds,60)
    istart=RUNDOMAIN(1)
    jstart=RUNDOMAIN(3)
    iend=RUNDOMAIN(2)
    jend=RUNDOMAIN(4)
  elseif(MODE == 10.or.MODE == 12)then
    if(.not.WriteNow)return
    istart=RUNDOMAIN(1)
    jstart=RUNDOMAIN(3)
    iend=RUNDOMAIN(2)
    jend=RUNDOMAIN(4)
  else
    if(mod(indate%hour,NHOURSAVE)/=0.or.indate%seconds/=0)return
  endif


!222 FORMAT(A,I2.2,I4.4,A)
!    write(fileName_write,222)'EMEP_BC_',indate%month,indate%year,'.nc'!for different names each month
                                                                     !NB: readxn should have same name
  if(me==0)write(*,*)'write Nest data ',trim(fileName_write)

  iotyp=IOU_INST
  if(first_call)then
    if(me==0)then
      write(*,*)'Writing BC on ',trim(fileName_write)
     !write(command,*)'rm ',trim(fileName_write)
     !call system(command)
    endif
    first_call=.false.
  endif

  ndim=3 !3-dimensional
  kmax=KMAX_MID
  scale=1.0
  def1%class='Advected' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=scale      !not used
  def1%rho=.false.      !not used
  def1%inst=.true.      !not used
  def1%year=.false.     !not used
  def1%month=.false.    !not used
  def1%day=.false.      !not used
  def1%name=''        !written
  def1%unit='mix_ratio'       !written

  do n= 1, NSPEC_ADV
 !do n= 1, NSPEC_ADV-4  !ENEA
    def1%name= species(NSPEC_SHL+n)%name       !written
    dat=xn_adv(n,:,:,:)
    call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real4,&
        ist=istart,jst=jstart,ien=iend,jen=jend,fileName_given=fileName_write)
  enddo

  return
end subroutine wrtxn


  subroutine check(status)
    use netcdf
    implicit none
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
         WRITE(*,*) 'MPI_ABORT: ', "errorin NetCDF_ml"
         call  MPI_ABORT(MPI_COMM_WORLD,9,INFO)
    end if
  end subroutine check

subroutine datefromsecondssince1970(ndate,nseconds,printdate)
!calculate date from seconds that have passed since the start of the year 1970

  !use Dates_ml, only : nmdays
  implicit none
  integer, intent(out) :: ndate(4)
  integer, intent(in) :: nseconds
  integer,  intent(in) :: printdate

  integer :: n,nmdays(12),nmdays2(13) !,nday
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
 !ndate(5)=nseconds-n
  if(printdate>0) write(*,"(A,I5,3(A,I4),A,I10)")&
    'year: ',ndate(1),', month: ',ndate(2),', day: ',ndate(3),&
    ', hour: ',ndate(4),', seconds: ',nseconds-n
end subroutine datefromsecondssince1970

subroutine datefromdayssince1900(ndate,ndays,printdate)
!calculate date from seconds that have passed since the start of the year 1900

 !use Dates_ml, only : nmdays
  implicit none
  integer, intent(out) :: ndate(4)
  real(kind=8), intent(inout) :: ndays
  integer, intent(in) :: printdate
  real(kind=8) :: rn

  integer :: n,nmdays(12),nmdays2(13) !,nday
  nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)

!add 0.5 seconds to avoid numerical errors in (n<=ndays)
  ndays=ndays+halfsecond

  nmdays2(1:12)=nmdays
  nmdays2(13)=0
  ndate(1)=1899
  n=0
  do while(n<=ndays)
    n=n+365
    ndate(1)=ndate(1)+1
    if(mod(ndate(1),4)==0.and.ndate(1)/=1900)n=n+1
  enddo
  n=n-365
  if(mod(ndate(1),4)==0.and.ndate(1)/=1900)n=n-1
  if(mod(ndate(1),4)==0.and.ndate(1)/=1900)nmdays2(2)=29
  ndate(2)=0
  do while(n<=ndays)
    ndate(2)=ndate(2)+1
    n=n+nmdays2(ndate(2))
  enddo
  n=n-nmdays2(ndate(2))
  ndate(3)=0
  do while(n<=ndays)
    ndate(3)=ndate(3)+1
    n=n+1
  enddo
  rn=n-1
  ndate(4)=-1
  do while(rn<=ndays)
    ndate(4)=ndate(4)+1
    rn=rn+1/24.0
  enddo
  rn=rn-1/24.0

!correct for modification
  ndays=ndays-halfsecond
 !ndate(5)=(ndays-rn)*24*3600.0
  if(printdate>0) write(*,"(A,I5,3(A,I4),A,F10.2)")&
    'year: ',ndate(1),', month: ',ndate(2),', day: ',ndate(3),&
    ', hour: ',ndate(4),', seconds: ',(ndays-rn)*24*3600.0
end subroutine datefromdayssince1900

subroutine init_icbc()
  implicit none
  logical, save :: first_call=.true.
  integer :: n

  if(.not.first_call)return
  first_call=.false.

  if(all(adv_ic%varname==""))then
    adv_ic(:)%varname=species(NSPEC_SHL+1:NSPEC_SHL+NSPEC_ADV)%name
    adv_ic(:)%wanted=.true.
    adv_ic(:)%found=find_icbc(filename_read_3D,adv_ic%varname(:))
  endif
  if(all(adv_bc%varname==""))then
    if(FORECAST)then	! IFS-MOZART BC
      adv_bc(:)%wanted=.false.
      adv_bc(FORECAST_BC%ixadv)=FORECAST_BC%icbc
      adv_bc(:)%found=find_icbc(filename_read_bc,adv_bc%varname(:))
    else
      adv_bc(:)=adv_ic(:)
    endif
  endif

  if(DEBUG_ICBC.and.me==0)then
    print "(A)","DEBUG_ICBC Variebles:"
    print "(2(X,A,I3,'=',A24,2L2))",&
      ('ADV_IC',n,adv_ic(n),'ADV_BC',n,adv_bc(n),n=1,NSPEC_ADV)
  endif
contains
function find_icbc(filename_read,varname) result(found)
  implicit none
  character(len=*), intent(in)               :: filename_read
  character(len=*), dimension(:), intent(in) :: varname
  logical, dimension(size(varname))          :: found
  integer :: status,ncFileID,varID,n

  found(:)=.false.
  if(me==0)then
    status = nf90_open(path=trim(filename_read),mode=nf90_nowrite,ncid=ncFileID)
    if(status /= nf90_noerr) then
      print *,'not found ',trim(filename_read)
    else
      print *,'  reading ',trim(filename_read)
      do n=1,size(varname)
        if(varname(n)/="") &
          found(n)=(nf90_inq_varid(ncid=ncFileID,name=trim(varname(n)),varID=varID)==nf90_noerr)
      enddo
    endif
  endif
  CALL MPI_BCAST(found,size(found),MPI_LOGICAL,0,MPI_COMM_WORLD,INFO)
end function find_icbc
end subroutine init_icbc

subroutine init_nest(ndays_indate,filename_read,IIij,JJij,Weight,&
                      Next,KMAX_ext,GIMAX_ext,GJMAX_ext)

  implicit none
  character(len=*),intent(in) :: filename_read
  real ,intent(out):: Weight(MAXLIMAX,MAXLJMAX,4)
  integer ,intent(out)::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)
  integer ,intent(out)::Next,KMAX_ext,GIMAX_ext,GJMAX_ext
  real(kind=8) :: ndays_indate,ndays(1)
  integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,status !,timeVarID
  integer :: ndate(4) !nseconds_indate,
  real :: dist(0:4)
  integer :: i,j,II,JJ !nseconds(1),n,n1,k
  real, allocatable, dimension(:,:) ::lon_ext,lat_ext
  character(len=80) ::projection

  rtime_saved = -99999.9 !initialization

!Read dimensions (global)
  if(me==0)then
    status = nf90_open(path=trim(filename_read),mode=nf90_nowrite,ncid=ncFileID)
    if(status /= nf90_noerr) then
      print *,'not found',trim(filename_read)
      return
    else
      print *,'  reading ',trim(filename_read)
    endif

    projection=''
    call check(nf90_get_att(ncFileID,nf90_global,"projection",projection))
    write(*,*)'projection: ',trim(projection)
    !get dimensions id
    if(trim(projection)=='Stereographic') then
      call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
      call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
    elseif(trim(projection)==trim('lon lat').or. &
          trim(projection)==trim('lon_lat')) then
      projection='lon lat'
      call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
      call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
    else
     !write(*,*)'GENERAL PROJECTION ',trim(projection)
      call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
      call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
     !WRITE(*,*) 'MPI_ABORT: ', "PROJECTION NOT RECOGNIZED"
     !call  MPI_ABORT(MPI_COMM_WORLD,9,INFO)
    endif

    call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
    call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

    call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_ext))
    call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_ext))
    call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_ext))
    call check(nf90_inquire_dimension(ncid=ncFileID,dimID=timedimID,len=Next))

    write(*,*)'dimensions external grid',GIMAX_ext,GJMAX_ext,KMAX_ext,Next
  endif
  CALL MPI_BCAST(GIMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(GJMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(KMAX_ext,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(Next,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

  allocate(lon_ext(GIMAX_ext,GJMAX_ext))
  allocate(lat_ext(GIMAX_ext,GJMAX_ext))

  if(me==0)then
   !Read lon lat of the external grid (global)
    if(trim(projection)==trim('lon lat')) then
      call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
      call check(nf90_get_var(ncFileID, varID, lon_ext(:,1) ))
      do i=1,GJMAX_ext
        lon_ext(:,i)=lon_ext(:,1)
      enddo
      call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
      call check(nf90_get_var(ncFileID, varID, lat_ext(1,:) ))
      do i=1,GIMAX_ext
        lat_ext(i,:)=lat_ext(1,:)
      enddo
    else
      call check(nf90_inq_varid(ncid = ncFileID, name = "lon", varID = varID))
      call check(nf90_get_var(ncFileID, varID, lon_ext ))

      call check(nf90_inq_varid(ncid = ncFileID, name = "lat", varID = varID))
      call check(nf90_get_var(ncFileID, varID, lat_ext ))
    endif
    call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))

   !do n=1,Next
    call check(nf90_get_var(ncFileID, varID, ndays,start=(/ 1 /),count=(/ 1 /) ))

    if(ndays(1)-ndays_indate>halfsecond)then
      write(*,*)'WARNING: did not find BIC for date:'
      call datefromdayssince1900(ndate,ndays_indate,1)
      write(*,*)'first date found:'
      call datefromdayssince1900(ndate,ndays(1),1)
    endif
   !enddo

    call check(nf90_close(ncFileID))
  endif

  CALL MPI_BCAST(lon_ext,8*GIMAX_ext*GJMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)
  CALL MPI_BCAST(lat_ext,8*GIMAX_ext*GJMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

 !find interpolation constants
 !note that i,j are local
 !find the four closest points
  do j=1,ljmax
    do i=1,limax
      dist=1.0E40
      do JJ=1,GJMAX_ext
        do II=1,GIMAX_ext
         !distance between (i,j) and (II,JJ)
          dist(0)=great_circle_distance(lon_ext(II,JJ),lat_ext(II,JJ),gl(i,j),gb(i,j))
          if(dist(0)<dist(1))then
            dist(4)=dist(3)
            dist(3)=dist(2)
            dist(2)=dist(1)
            dist(1)=dist(0)
            IIij(i,j,4)=IIij(i,j,3)
            JJij(i,j,4)=JJij(i,j,3)
            IIij(i,j,3)=IIij(i,j,2)
            JJij(i,j,3)=JJij(i,j,2)
            IIij(i,j,2)=IIij(i,j,1)
            JJij(i,j,2)=JJij(i,j,1)
            IIij(i,j,1)=II
            JJij(i,j,1)=JJ
          elseif(dist(0)<dist(2))then
            dist(4)=dist(3)
            dist(3)=dist(2)
            dist(2)=dist(0)
            IIij(i,j,4)=IIij(i,j,3)
            JJij(i,j,4)=JJij(i,j,3)
            IIij(i,j,3)=IIij(i,j,2)
            JJij(i,j,3)=JJij(i,j,2)
            IIij(i,j,2)=II
            JJij(i,j,2)=JJ
          elseif(dist(0)<dist(3))then
            dist(4)=dist(3)
            dist(3)=dist(0)
            IIij(i,j,4)=IIij(i,j,3)
            JJij(i,j,4)=JJij(i,j,3)
            IIij(i,j,3)=II
            JJij(i,j,3)=JJ
          elseif(dist(0)<dist(4))then
            dist(4)=dist(0)
            IIij(i,j,4)=II
            JJij(i,j,4)=JJ
          endif
        enddo
      enddo

      dist(0)=(dist(1)+dist(2)+dist(3)+dist(4))
      Weight(i,j,1)=1.0-3.0*dist(1)/dist(0)
      dist(0)=(dist(2)+dist(3)+dist(4))
      Weight(i,j,2)=(1.0-Weight(i,j,1))*(1.0-2.0*dist(2)/dist(0))
      dist(0)=(dist(3)+dist(4))
      Weight(i,j,3)=(1.0-Weight(i,j,1)-Weight(i,j,2))*(1.0-dist(3)/dist(0))
      Weight(i,j,4)=1.0-Weight(i,j,1)-Weight(i,j,2)-Weight(i,j,3)
    enddo
  enddo

  deallocate(lon_ext,lat_ext)
end subroutine init_nest

subroutine read_newdata_LATERAL(ndays_indate,nr)
  implicit none
  real(kind=8), intent(in)::ndays_indate
  real, allocatable, dimension(:,:,:) ::data
  integer :: ncFileID,varid,status
  integer :: ndate(4),n,i,j,k,nr
  real(kind=8) :: ndays(1),ndays_old
  logical, save :: first_call=.true.

  !4 nearest points from external grid
  integer, save ::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)

  !weights of the 4 nearest points
  real, save :: Weight(MAXLIMAX,MAXLJMAX,4)

  !dimensions of external grid for BC
  integer, save ::GIMAX_ext,GJMAX_ext

  if(first_call)then
    first_call=.false.
    call init_icbc()
    call init_nest(ndays_indate,filename_read_BC,IIij,JJij,Weight,&
                  Next_BC,KMAX_ext_BC,GIMAX_ext,GJMAX_ext)

    !Define & allocate West/East/South/Nort Boundaries
    iw=li0-1;ie=li1+1   ! i West/East   boundaries
    js=lj0-1;jn=lj1+1   ! j South/North boundaries
    kt=0;if(TOP_BC)kt=1 ! k Top         boundary
    if(iw>=1    .and..not.allocated(xn_adv_bndw)) &
      allocate(xn_adv_bndw(NSPEC_ADV,MAXLJMAX,KMAX_MID,2)) ! West
    if(ie<=limax.and..not.allocated(xn_adv_bnde)) &
      allocate(xn_adv_bnde(NSPEC_ADV,MAXLJMAX,KMAX_MID,2)) ! East
    if(js>=1    .and..not.allocated(xn_adv_bnds)) &
      allocate(xn_adv_bnds(NSPEC_ADV,MAXLIMAX,KMAX_MID,2)) ! South
    if(jn<=ljmax.and..not.allocated(xn_adv_bndn)) &
      allocate(xn_adv_bndn(NSPEC_ADV,MAXLIMAX,KMAX_MID,2)) ! North
    if(kt>=1    .and..not.allocated(xn_adv_bndt)) &
      allocate(xn_adv_bndt(NSPEC_ADV,MAXLIMAX,MAXLJMAX,2)) ! Top
    if(DEBUG_ICBC)then
      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
      if(me==0) print "(A)","DEBUG_ICBC Boundaries:"
      print "(X,'me=',i3,5(X,A,I0,'=',L1))",&
        me,'W:i',iw,allocated(xn_adv_bndw),'E:i',ie,allocated(xn_adv_bnde),&
           'S:j',js,allocated(xn_adv_bnds),'N:j',jn,allocated(xn_adv_bndn),&
           'T:k',kt,allocated(xn_adv_bndt)
      call flush(6)
      CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
    endif
  endif

  ndays_old=rtime_saved(2)
  allocate(data(GIMAX_ext,GJMAX_ext,KMAX_ext_BC), stat=status)
  if(me==0)then
    call check(nf90_open(path = trim(fileName_read_BC), mode = nf90_nowrite, ncid = ncFileID))

    call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))
    do n=1,Next_BC
      call check(nf90_get_var(ncFileID, varID, ndays,start=(/ n /),count=(/ 1 /) ))
      if(ndays_indate-ndays(1)<halfsecond)then
        write(*,*)'Using date '
        call datefromdayssince1900(ndate,ndays(1),1)
        goto 876
      endif
    enddo
    write(*,*)'WARNING: did not find correct date'
    n=Next_BC
    write(*,*)'Using date '
    call datefromdayssince1900(ndate,ndays(1),1)
876 continue
    itime=n
    rtime_saved(2)=ndays(1)
  endif

  CALL MPI_BCAST(rtime_saved,8*2,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

  DO_SPEC: do n= 1, NSPEC_ADV
    if(.not.(adv_bc(n)%wanted.and.adv_bc(n)%found)) cycle DO_SPEC
    if(nr==2)then
      !store the old values in 1
      rtime_saved(1)=ndays_old
      if(iw>=1)     forall (k=1:KMAX_ext_BC, j=1:ljmax) &
        xn_adv_bndw(n,j,k,1)=xn_adv_bndw(n,j,k,2)
      if(ie<=limax) forall (k=1:KMAX_ext_BC, j=1:ljmax) &
        xn_adv_bnde(n,j,k,1)=xn_adv_bnde(n,j,k,2)
      if(js>=1)     forall (k=1:KMAX_ext_BC, i=1:limax) &
        xn_adv_bnds(n,i,k,1)=xn_adv_bnds(n,i,k,2)
      if(jn<=ljmax) forall (k=1:KMAX_ext_BC, i=1:limax) &
        xn_adv_bndn(n,i,k,1)=xn_adv_bndn(n,i,k,2)
      if(kt>=1)     forall (i=1:limax, j=1:ljmax) &
        xn_adv_bndt(n,i,j,1)=xn_adv_bndt(n,i,j,2)
    endif
    if(me==0)then
    !Could fetch one level at a time if sizes becomes too big
      call check(nf90_inq_varid(ncid=ncFileID, name=trim(adv_bc(n)%varname), varID=varID))

      call check(nf90_get_var(ncFileID, varID, data &
            ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext_BC,1 /) ))
    endif
    CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext_BC,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

   !overwrite Global Boundaries (lateral faces)
    if(iw>=1)     forall (k=1:KMAX_ext_BC, j=1:ljmax) &
      xn_adv_bndw(n,j,k,2)=Weight(iw,j,1)*data(IIij(iw,j,1),JJij(iw,j,1),k) &
                          +Weight(iw,j,2)*data(IIij(iw,j,2),JJij(iw,j,2),k) &
                          +Weight(iw,j,3)*data(IIij(iw,j,3),JJij(iw,j,3),k) &
                          +Weight(iw,j,4)*data(IIij(iw,j,4),JJij(iw,j,4),k)
    if(ie<=limax) forall (k=1:KMAX_ext_BC, j=1:ljmax) &
      xn_adv_bnde(n,j,k,2)=Weight(ie,j,1)*data(IIij(ie,j,1),JJij(ie,j,1),k) &
                          +Weight(ie,j,2)*data(IIij(ie,j,2),JJij(ie,j,2),k) &
                          +Weight(ie,j,3)*data(IIij(ie,j,3),JJij(ie,j,3),k) &
                          +Weight(ie,j,4)*data(IIij(ie,j,4),JJij(ie,j,4),k)
    if(js>=1)     forall (k=1:KMAX_ext_BC, i=1:limax) &
      xn_adv_bnds(n,i,k,2)=Weight(i,js,1)*data(IIij(i,js,1),JJij(i,js,1),k) &
                          +Weight(i,js,2)*data(IIij(i,js,2),JJij(i,js,2),k) &
                          +Weight(i,js,3)*data(IIij(i,js,3),JJij(i,js,3),k) &
                          +Weight(i,js,4)*data(IIij(i,js,4),JJij(i,js,4),k)
    if(jn<=ljmax) forall (k=1:KMAX_ext_BC, i=1:limax) &
      xn_adv_bndn(n,i,k,2)=Weight(i,jn,1)*data(IIij(i,jn,1),JJij(i,jn,1),k) &
                          +Weight(i,jn,2)*data(IIij(i,jn,2),JJij(i,jn,2),k) &
                          +Weight(i,jn,3)*data(IIij(i,jn,3),JJij(i,jn,3),k) &
                          +Weight(i,jn,4)*data(IIij(i,jn,4),JJij(i,jn,4),k)
    if(kt>=1)     forall (i=1:limax, j=1:ljmax) &
      xn_adv_bndt(n,i,j,2)=Weight(i,j,1)*data(IIij(i,j,1),JJij(i,j,1),kt) &
                          +Weight(i,j,2)*data(IIij(i,j,2),JJij(i,j,2),kt) &
                          +Weight(i,j,3)*data(IIij(i,j,3),JJij(i,j,3),kt) &
                          +Weight(i,j,4)*data(IIij(i,j,4),JJij(i,j,4),kt)
  enddo DO_SPEC

  deallocate(data)
  if(me==0) call check(nf90_close(ncFileID))
end subroutine read_newdata_LATERAL

subroutine reset_3D(ndays_indate)
  implicit none
  real(kind=8), intent(in)::ndays_indate
  real, allocatable, dimension(:,:,:) ::data
  integer :: ndate(4),n,i,j,k,itime=0,status
  integer :: ncFileID,varid
  real(kind=8) :: ndays(1)
  logical, save :: first_call=.true.

  !4 nearest points from external grid
  integer, save ::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)

  !weights of the 4 nearest points
  real, save :: Weight(MAXLIMAX,MAXLJMAX,4)

  !dimensions of external grid for 3D
  integer, save ::Next,KMAX_ext,GIMAX_ext,GJMAX_ext

  if(first_call)then
    first_call=.false.
    call init_icbc()
    call init_nest(ndays_indate,filename_read_3D,IIij,JJij,Weight,&
                  Next,KMAX_ext,GIMAX_ext,GJMAX_ext)
  endif
  allocate(data(GIMAX_ext,GJMAX_ext,KMAX_ext), stat=status)
  if(me==0)then
    call check(nf90_open(path = trim(fileName_read_3D), mode = nf90_nowrite, ncid = ncFileID))

    call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = varID))
    do n=1,Next
      call check(nf90_get_var(ncFileID, varID, ndays,start=(/ n /),count=(/ 1 /) ))
      if(ndays(1)>=ndays_indate)then
        write(*,*)'found date '
        call datefromdayssince1900(ndate,ndays(1),1)
        goto 876
      endif
    enddo
    write(*,*)'WARNING: did not find correct date'
    n=Next
876 continue
    itime=n
  endif

  DO_SPEC: do n= 1, NSPEC_ADV
    if(.not.(adv_ic(n)%wanted.and.adv_ic(n)%found)) cycle DO_SPEC
    if(me==0)then
    !Could fetch one level at a time if sizes becomes too big
      call check(nf90_inq_varid(ncid=ncFileID, name=trim(species(NSPEC_SHL+n)%name), varID=varID))

      call check(nf90_get_var(ncFileID, varID, data &
            ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext,1 /) ))
    endif
    CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

   !overwrite everything 3D (init)
    forall (k=1:KMAX_ext, j=1:ljmax, i=1:limax) &
      xn_adv(n,i,j,k)=Weight(i,j,1)*data(IIij(i,j,1),JJij(i,j,1),k) &
                      +Weight(i,j,2)*data(IIij(i,j,2),JJij(i,j,2),k) &
                      +Weight(i,j,3)*data(IIij(i,j,3),JJij(i,j,3),k) &
                      +Weight(i,j,4)*data(IIij(i,j,4),JJij(i,j,4),k)
  enddo DO_SPEC

  deallocate(data)
  if(me==0) call check(nf90_close(ncFileID))
end subroutine reset_3D

end module Nest_ml

