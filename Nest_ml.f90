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
  use ChemSpecs_adv_ml , only :NSPEC_ADV
  use ChemSpecs_tot_ml , only :NSPEC_TOT
  use netcdf
  use netcdf_ml,      only : GetCDF,Out_netCDF,Init_new_netCDF,&
       secondssince1970,dayssince1900,Int1,Int2,Int4,Real4,Real8
  use Functions_ml,    only : great_circle_distance
  use ModelConstants_ml,    only : KMAX_MID, NPROC &
          , FORECAST & ! AMVB 2009-11-06: nested input/output on FORECAST mode
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


! AMVB 2009-11-06: nested input/output on FORECAST mode
  integer, public, parameter :: FORECAST_NDUMP = 1  ! Number of nested output
  ! on FORECAST mode (1: starnt next forecast; 2: NMC statistics)
  type(date), public :: outdate(FORECAST_NDUMP)=date(-1,-1,-1,-1,-1)
  ! Nested output dates on FORECAST mode

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

  !BC values at boundaries in present grid
  real,save, dimension(NSPEC_ADV,MAXLIMAX,KMAX_MID,2) :: xn_adv_bnds,xn_adv_bndn !north and south
  real,save, dimension(NSPEC_ADV,MAXLJMAX,KMAX_MID,2) :: xn_adv_bndw,xn_adv_bnde !west and east

  !dimension of external grid for BC
  integer,save :: Next_BC,KMAX_ext_BC

  integer,save :: itime!itime_saved(2),
  real*8,save :: rtime_saved(2)
  character*30,save  :: filename_read_BC='EMEP_IN.nc'
  character*30,save  :: filename_read_3D='EMEP_IN.nc'
  character*30,save  :: filename_write='EMEP_OUT.nc'
  real*8, parameter :: halfsecond=1.0/(24.0*3600.0)!used to avoid rounding errors

contains

  subroutine readxn(indate)
    implicit none
    type(date), intent(in) :: indate           ! Gives year..seconds
    integer,save  :: first_data=-1


    integer :: nseconds(1),n1,n,i,j,k,II,JJ
    integer :: nstart,nfetch,ndate(4),nseconds_indate
    real*8:: ndays_indate

    !    real , dimension(48,48,20) ::data

    real :: W1,W2
    logical, save :: first_call=.true.
    logical :: fexist=.false. ! AMVB 2009-11-06: nested input/output on FORECAST mode

    if(MODE /= 2.and.MODE /= 3.and. MODE /= 11.and. MODE /= 12.and. .not.FORECAST)return ! AMVB 2009-11-06: nested input/output on FORECAST mode

    ndate(1)  = indate%year
    ndate(2)  = indate%month
    ndate(3)  = indate%day
    ndate(4)  = indate%hour
!    call secondssince1970(ndate,nseconds_indate)
    call dayssince1900(ndate,ndays_indate)

! AMVB 2009-11-06: nested input/output on FORECAST mode
    if(FORECAST)then ! FORECAST mode superseeds nest MODE
      if(.not. first_call)return
      first_call=.false.
      inquire(file=filename_read_3D,exist=fexist)
      if(.not. fexist)then
        if(me==0) print *,'No nest/dump file found: ',trim(filename_read_3D)
        return
      end if
      if(me==0)   print *,'RESET ALL XN 3D'
      call reset_3D(ndays_indate)
      return
    elseif(MODE == 11.or.MODE == 12)then
       if(.not. first_call)return
       first_call=.false.
       if(me==0)   print *,'RESET ALL XN 3D'
       call reset_3D(ndays_indate)
       return
    else
!    if(me==0)   print *,'call to READXN',indate%hour,indate%seconds
       if(mod(indate%hour,NHOURREAD)/=0.or.indate%seconds/=0)return
    endif

!never comes to this point if MODE=FORECAST, 11 or 12

    if(me==0)   print *,'NESTING'

    if(first_data==-1)then
       call reset_3D(ndays_indate)
       call read_newdata_LATERAL(ndays_indate,1)
       call read_newdata_LATERAL(ndays_indate,2)
    endif


    if(ndays_indate-rtime_saved(2)>halfsecond)then
       !look for a new data set
       if(me==0)write(*,*)'NEST: READING NEW BC DATA'
       call read_newdata_LATERAL(ndays_indate,2)
    endif


    !    make weights for time interpolation
    W1=1.0;  W2=0.0!default
    if(ndays_indate-rtime_saved(1)>halfsecond)then
       !interpolate
       W2=(ndays_indate-rtime_saved(1))/(rtime_saved(2)-rtime_saved(1))
       W1=1.0-W2
!       if(me==1)then
!       call datefromdayssince1900(ndate,ndays_indate,1)
!       write(*,*)'interpolating between'
!       call datefromdayssince1900(ndate,rtime_saved(1),1)
!       write(*,*)'and'
!       call datefromdayssince1900(ndate,rtime_saved(2),1)
!       write(*,*)'with weights : ',W1,W2
!       endif
    endif
!    if(me==0)write(*,*)'weights : ',W1,W2,rtime_saved(1),rtime_saved(2)

    do n=1,NSPEC_ADV
       do k=1,KMAX_ext_BC
          do i=1,li0-1
             do j=1,ljmax
                xn_adv(n,i,j,k)=W1*xn_adv_bndw(n,j,k,1)+W2*xn_adv_bndw(n,j,k,2)
             enddo
          enddo
          do j=1,lj0-1
             do i=1,limax
                xn_adv(n,i,j,k)=W1*xn_adv_bnds(n,i,k,1)+W2*xn_adv_bnds(n,i,k,2)
             enddo
          enddo
          do i=li1+1,limax
             do j=1,ljmax
                xn_adv(n,i,j,k)=W1*xn_adv_bnde(n,j,k,1)+W2*xn_adv_bnde(n,j,k,2)
             enddo
          enddo
          do j=lj1+1,ljmax
             do i=1,limax
                xn_adv(n,i,j,k)=W1*xn_adv_bndn(n,i,k,1)+W2*xn_adv_bndn(n,i,k,2)
             enddo
          enddo
       enddo
    enddo


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
    character *40:: command

    if(MODE /= 1.and.MODE /= 3.and.MODE /= 10.and.MODE /= 12.and. .not.FORECAST)return ! AMVB 2009-11-06: nested input/output on FORECAST mode

 ! AMVB 2009-11-06: nested input/output on FORECAST mode
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


    222 FORMAT(A,I2.2,I4.4,A)
!    write(fileName_write,222)'EMEP_BC_',indate%month,indate%year,'.nc'!for different names each month
                                                                      !NB: readxn should have same name
    if(me==0)write(*,*)'write Nest data ',trim(fileName_write)



    iotyp=IOU_INST
    if(first_call)then
       if(me==0)then
          write(*,*)'Writing BC on ',trim(fileName_write)
!          write(command,*)'rm ',trim(fileName_write)
!          call system(command)
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
!    do n= 1, NSPEC_ADV-4  !ENEA
       def1%name= species(NSPEC_SHL+n)%name       !written
       dat=xn_adv(n,:,:,:)
       call Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype=Real4,&
            ist=istart,jst=jstart,ien=iend,jen=jend&
            ,fileName_given=fileName_write)
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

    !  use Dates_ml, only : nmdays
    implicit none

    integer, intent(out) :: ndate(4)
    integer, intent(in) :: nseconds
    integer,  intent(in) :: printdate

    integer :: n,nday,nmdays(12),nmdays2(13)
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
    !    ndate(5)=nseconds-n
    if(printdate>0)then
       write(*,55)'year: ',ndate(1),', month: ',ndate(2),', day: ',&
            ndate(3),', hour: ',ndate(4),', seconds: ',nseconds-n
    endif
55  format(A,I5,A,I4,A,I4,A,I4,A,I10)
  end subroutine datefromsecondssince1970

  subroutine datefromdayssince1900(ndate,ndays,printdate)
    !calculate date from seconds that have passed since the start of the year 1900

    !  use Dates_ml, only : nmdays
    implicit none

    integer, intent(out) :: ndate(4)
    real*8, intent(inout) :: ndays
    integer,  intent(in) :: printdate
    real*8 :: rn

    integer :: n,nday,nmdays(12),nmdays2(13)
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
    !    ndate(5)=(ndays-rn)*24*3600.0
    if(printdate>0)then
       write(*,55)'year: ',ndate(1),', month: ',ndate(2),', day: ',&
            ndate(3),', hour: ',ndate(4),', seconds: ',(ndays-rn)*24*3600.0
    endif
55  format(A,I5,A,I4,A,I4,A,I4,A,F10.2)
  end subroutine datefromdayssince1900

  subroutine init_nest(ndays_indate,filename_read,IIij,JJij,&
       Weight1,Weight2,Weight3,Weight4,Next,KMAX_ext,GIMAX_ext,GJMAX_ext)

    implicit none
    character(len=*),intent(in) :: filename_read
    real ,intent(out):: Weight1(MAXLIMAX,MAXLJMAX),Weight2(MAXLIMAX,MAXLJMAX)
    real ,intent(out):: Weight3(MAXLIMAX,MAXLJMAX),Weight4(MAXLIMAX,MAXLJMAX)
    integer ,intent(out)::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)
    integer ,intent(out)::Next,KMAX_ext,GIMAX_ext,GJMAX_ext
    real*8 :: ndays_indate,ndays(1)
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID,status
    integer :: nseconds_indate,ndate(4)
    real :: dist(0:4)
    integer :: nseconds(1),n1,n,i,j,k,II,JJ
    real, allocatable, dimension(:,:) ::lon_ext,lat_ext
    character*80 ::projection

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
       elseif(trim(projection)==trim('lon lat')) then
          call check(nf90_inq_dimid(ncid = ncFileID, name = "lon", dimID = idimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "lat", dimID = jdimID))
       else
          !     write(*,*)'GENERAL PROJECTION ',trim(projection)
          call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
          call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
          !       WRITE(*,*) 'MPI_ABORT: ', "PROJECTION NOT RECOGNIZED"
            !     call  MPI_ABORT(MPI_COMM_WORLD,9,INFO)
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

       !          do n=1,Next
       call check(nf90_get_var(ncFileID, varID, ndays,start=(/ 1 /),count=(/ 1 /) ))

       if(ndays(1)-ndays_indate>halfsecond)then
          write(*,*)'WARNING: did not find BIC for date:'
          call datefromdayssince1900(ndate,ndays_indate,1)
          write(*,*)'first date found:'
          call datefromdayssince1900(ndate,ndays(1),1)
       endif
       !          enddo

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
          Weight1(i,j)=1.0-3.0*dist(1)/dist(0)
          dist(0)=(dist(2)+dist(3)+dist(4))
          Weight2(i,j)=(1.0-Weight1(i,j))*(1.0-2.0*dist(2)/dist(0))
          dist(0)=(dist(3)+dist(4))
          Weight3(i,j)=(1.0-Weight1(i,j)-Weight2(i,j))*(1.0-dist(3)/dist(0))
          Weight4(i,j)=1.0-Weight1(i,j)-Weight2(i,j)-Weight3(i,j)

       enddo
    enddo

    deallocate(lon_ext)
    deallocate(lat_ext)

  end subroutine init_nest

  subroutine read_newdata_LATERAL(ndays_indate,nr)

    implicit none
    real*8, intent(in)::ndays_indate
    real, allocatable, dimension(:,:,:) ::data
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID,status
    integer :: nseconds(1),ndate(4),n1,n,i,j,k,II,JJ,nseconds_indate,nr
    integer :: nseconds_old
    real*8:: ndays(1),ndays_old
    logical, save :: first_call=.true.

    !BC values at boundaries in present grid
    real,save, dimension(NSPEC_ADV,MAXLIMAX,KMAX_MID,2) :: xn_adv_bnds,xn_adv_bndn !north and south
    real,save, dimension(NSPEC_ADV,MAXLJMAX,KMAX_MID,2) :: xn_adv_bndw,xn_adv_bnde !west and east

    !4 nearest points from external grid
    integer, save ::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)

    !weights of the 4 nearest points
    real, save :: Weight1(MAXLIMAX,MAXLJMAX),Weight2(MAXLIMAX,MAXLJMAX)
    real, save :: Weight3(MAXLIMAX,MAXLJMAX),Weight4(MAXLIMAX,MAXLJMAX)

    !dimensions of external grid for BC
    integer, save ::GIMAX_ext,GJMAX_ext

    if(first_call)then
       first_call=.false.
       call init_nest(ndays_indate,filename_read_BC,IIij,JJij,&
            Weight1,Weight2,Weight3,Weight4,Next_BC,KMAX_ext_BC,GIMAX_ext,GJMAX_ext)
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
       itime=Next_BC
       write(*,*)'Using date '
       call datefromdayssince1900(ndate,ndays(1),1)
876    continue
       itime=n
       rtime_saved(2)=ndays(1)
    endif

      CALL MPI_BCAST(rtime_saved,8*2,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

    do n= 1, NSPEC_ADV
       if(nr==2)then
          !store the old values in 1
          rtime_saved(1)=ndays_old
          do k=1,KMAX_ext_BC
             do i=1,li0-1
                do j=1,ljmax
                   xn_adv_bndw(n,j,k,1)=xn_adv_bndw(n,j,k,2)
                enddo
             enddo
             do j=1,lj0-1
                do i=1,limax
                   xn_adv_bnds(n,i,k,1)=xn_adv_bnds(n,i,k,2)
                enddo
             enddo
             do i=li1+1,limax
                do j=1,ljmax
                   xn_adv_bnde(n,j,k,1)=xn_adv_bnde(n,j,k,2)
                enddo
             enddo
             do j=lj1+1,ljmax
                do i=1,limax
                   xn_adv_bndn(n,i,k,1)=xn_adv_bndn(n,i,k,2)
                enddo
             enddo
          enddo
       endif
       if(me==0)then
          !Could fetch one level at a time if sizes becomes too big

          call check(nf90_inq_varid(ncid=ncFileID, name=trim(species(NSPEC_SHL+n)%name), varID=varID))

          call check(nf90_get_var(ncFileID, varID, data &
               ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext_BC,1 /) ))

       endif
         CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext_BC,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

       !overwrite Global Boundaries (lateral faces)

       do k=1,KMAX_ext_BC
          do i=1,li0-1
             do j=1,ljmax
                xn_adv_bndw(n,j,k,2)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)

             enddo
          enddo
          do j=1,lj0-1
             do i=1,limax
                xn_adv_bnds(n,i,k,2)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)
             enddo
          enddo
          do i=li1+1,limax
             do j=1,ljmax
                xn_adv_bnde(n,j,k,2)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)
             enddo
          enddo
          do j=lj1+1,ljmax
             do i=1,limax
                xn_adv_bndn(n,i,k,2)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)
             enddo
          enddo
       enddo
    enddo

    deallocate(data)
    if(me==0)then
       call check(nf90_close(ncFileID))
    endif

  end subroutine read_newdata_LATERAL

  subroutine reset_3D(ndays_indate)
    implicit none
    real*8, intent(in)::ndays_indate
    real, allocatable, dimension(:,:,:) ::data
    integer :: nseconds(1),ndate(4),n1,n,i,j,k,II,JJ,itime,status
    integer :: nseconds_indate
    integer :: ncFileID,idimID,jdimID, kdimID,timeDimID,varid,timeVarID
    real*8 :: ndays(1)
    logical, save :: first_call=.true.

    !BC values at boundaries in present grid
    real,save, dimension(NSPEC_ADV,MAXLIMAX,KMAX_MID,2) :: xn_adv_bnds,xn_adv_bndn !north and south
    real,save, dimension(NSPEC_ADV,MAXLJMAX,KMAX_MID,2) :: xn_adv_bndw,xn_adv_bnde !west and east

    !4 nearest points from external grid
    integer, save ::IIij(MAXLIMAX,MAXLJMAX,4),JJij(MAXLIMAX,MAXLJMAX,4)

    !weights of the 4 nearest points
    real, save :: Weight1(MAXLIMAX,MAXLJMAX),Weight2(MAXLIMAX,MAXLJMAX)
    real, save :: Weight3(MAXLIMAX,MAXLJMAX),Weight4(MAXLIMAX,MAXLJMAX)

    !dimensions of external grid for 3D
    integer, save ::Next,KMAX_ext,GIMAX_ext,GJMAX_ext

    if(first_call)then
       first_call=.false.
       call init_nest(ndays_indate,filename_read_3D,IIij,JJij,&
            Weight1,Weight2,Weight3,Weight4,Next,KMAX_ext,GIMAX_ext,GJMAX_ext)
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
       itime=Next
876    continue
       itime=n
    endif

    do n= 1, NSPEC_ADV
       if(me==0)then
          !Could fetch one level at a time if sizes becomes too big

          call check(nf90_inq_varid(ncid=ncFileID, name=trim(species(NSPEC_SHL+n)%name), varID=varID))

          call check(nf90_get_var(ncFileID, varID, data &
               ,start=(/ 1,1,1,itime /),count=(/ GIMAX_ext,GJMAX_ext,KMAX_ext,1 /) ))

       endif
         CALL MPI_BCAST(data,8*GIMAX_ext*GJMAX_ext*KMAX_ext,MPI_BYTE,0,MPI_COMM_WORLD,INFO)

       ! overwrite everything 3D (init)
       do k=1,KMAX_ext
          do j=1,ljmax
             do i=1,limax
                xn_adv(n,i,j,k)=Weight1(i,j)*data(IIij(i,j,1),JJij(i,j,1),k)+&
                     Weight2(i,j)*data(IIij(i,j,2),JJij(i,j,2),k)+&
                     Weight3(i,j)*data(IIij(i,j,3),JJij(i,j,3),k)+&
                     Weight4(i,j)*data(IIij(i,j,4),JJij(i,j,4),k)
             enddo
          enddo
       enddo

    enddo

    deallocate(data)
    if(me==0)then
       call check(nf90_close(ncFileID))
    endif

  end subroutine reset_3D

end module Nest_ml

