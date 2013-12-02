module Util_ml
use netcdf
use ModelConstants_ml,      only: MasterProc
use TimeDate_ml,            only: date
use TimeDate_ExtraUtil_ml,  only: date2nctime,date2string,nctime2string
use CheckStop_ml,           only: CheckStop
implicit none
integer(4) :: numLev,numLon,numLat,numRec
public  :: norm
private :: norm_r8d1,norm_r8d2,norm_r8d3,norm_r8d4,&
           norm_c8d1,norm_c8d2,norm_c8d3,norm_c8d4
interface norm
  module procedure norm_r8d1,norm_r8d2,norm_r8d3,norm_r8d4,&
                   norm_c8d1,norm_c8d2,norm_c8d3,norm_c8d4
end interface norm
public  :: infovar
private :: info_1var2d,info_2var2d,info_3var2d
interface infovar
  module procedure info_1var2d,info_2var2d,info_3var2d
end interface infovar
public  :: allocvar
private :: allocvar_int1d,allocvar_real1d,allocvar_check
interface allocvar
  module procedure allocvar_int1d,allocvar_real1d
end interface allocvar
contains
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
subroutine allocvar_check(assoc,iaux,vsize,vrange,vname,mfmt)
  logical, intent(in) :: assoc
  integer, intent(in) :: iaux,vsize,vrange(0:1)
  character(len=*), intent(in) :: vname,mfmt
  call CheckStop(vsize<vrange(0).or.vsize>vrange(1),&
                 'Out or range size '//trim(vname))
  if(assoc)then
    call CheckStop(iaux/=vsize,'Wrong size '//trim(vname))
  else
    if(vsize<=0) write(*,mfmt)'WARNING size<=0 '//trim(vname)
  endif
endsubroutine allocvar_check
subroutine allocvar_int1d(var,vsize,vrange,vname,mfmt)
  integer, intent(out), dimension(:), pointer :: var
  integer, intent(in) :: vsize,vrange(0:1)
  character(len=*), intent(in) :: vname,mfmt
  integer :: ierr
  if(associated(var))then
    call allocvar_check(.true.,size(var),vsize,vrange,vname,mfmt)
  else
    allocate(var(max(vsize,0)),stat=ierr)
    call allocvar_check(.false.,ierr,vsize,vrange,vname,mfmt)
    var=0
  endif
endsubroutine allocvar_int1d
subroutine allocvar_real1d(var,vsize,vrange,vname,mfmt)
  real, intent(out), dimension(:), pointer :: var
  integer, intent(in) :: vsize,vrange(0:1)
  character(len=*), intent(in) :: vname,mfmt
  integer :: ierr
  if(associated(var))then
    call allocvar_check(.true.,size(var),vsize,vrange,vname,mfmt)
  else
    allocate(var(max(vsize,0)),stat=ierr)
    call allocvar_check(.false.,ierr,vsize,vrange,vname,mfmt)
    var=0
  endif
endsubroutine allocvar_real1d
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
function info_1var2d(var,trimdomain) result(msg)
  implicit none
  real, intent(in), dimension(:,:) :: var
  logical, intent(in), optional    :: trimdomain
  character(len=35) :: msg
  character(len=*), parameter :: msgfmt="(2(I3,'x',I3,':',G9.2,:,1X))"
  write(msg,msgfmt) minloc(var),minval(var),maxloc(var),maxval(var)
end function info_1var2d
function info_2var2d(var1,var2,trimdomain) result(msg)
  implicit none
  real, intent(in), dimension(:,:) :: var1,var2
  logical, intent(in), optional    :: trimdomain
  character(len=(35+1)*2-1) :: msg
  character(len=*), parameter :: msgfmt="(2(A35,:,1X))"
  logical :: td
  integer, dimension(2) :: lb,ub
  td=.false.;if(present(trimdomain))td=trimdomain
  if(td)then
    lb=max(lbound(var1),lbound(var2))
    ub=min(ubound(var1),ubound(var2))
    write(msg,msgfmt) info_1var2d(var1(lb(1):ub(1),lb(2):ub(2))),&
                      info_1var2d(var2(lb(1):ub(1),lb(2):ub(2)))
  else
    write(msg,msgfmt) info_1var2d(var1),info_1var2d(var2)
  endif
end function info_2var2d
function info_3var2d(var1,var2,var3,trimdomain) result(msg)
  implicit none
  real, intent(in), dimension(:,:) :: var1,var2,var3
  logical, intent(in), optional    :: trimdomain
  character(len=(35+1)*3-1) :: msg
  character(len=*), parameter :: msgfmt="(3(A35,:,1X))"
  logical :: td
  integer, dimension(2) :: lb,ub
  td=.false.;if(present(trimdomain))td=trimdomain
  if(td)then
    lb=max(lbound(var1),lbound(var2),lbound(var3))
    ub=min(ubound(var1),ubound(var2),ubound(var3))
    write(msg,msgfmt) info_1var2d(var1(lb(1):ub(1),lb(2):ub(2))),&
                      info_1var2d(var2(lb(1):ub(1),lb(2):ub(2))),&
                      info_1var2d(var3(lb(1):ub(1),lb(2):ub(2)))
  else
    write(msg,msgfmt) info_1var2d(var1),info_1var2d(var2),info_1var2d(var3)
  endif
end function info_3var2d
!+------------------------------------------------------------------
! nc_check (CDF_Utils.check)
!   check wether the function returns nf90_noerr
! io_check
!   check wether the function returns zero
!+------------------------------------------------------------------
subroutine io_check(status,msg)
  implicit none
  integer, intent (in) :: status
  character(len=*), intent (in) :: msg
  select case (status)
    case(0) ! correct IO
   !case()  ! recoverable errors
   !  if(MasterProc) print *,"RecoverableError in I/O call: "//trim(msg))
    case default
      call CheckStop("Error in I/O call: "//trim(msg))
   endselect
end subroutine io_check
subroutine nc_check(status)
  use netcdf
  implicit none
  integer, intent (in) :: status
  call CheckStop(status,nf90_noerr,"Error in NetCDF call: "//trim(nf90_strerror(status)))
end subroutine nc_check
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
function norm_r8d1(a,squared) result(norm)
  implicit none
  intent(in)      :: a,squared
  optional        :: squared
  logical         :: squared,sq
  real(kind=8)    :: a(:),norm
  sq=.false.;if(present(squared))sq=squared
  norm=sum(a**2)
  if(.not.sq)norm=sqrt(norm)
end function norm_r8d1
function norm_c8d1(a,squared) result(norm)
  implicit none
  intent(in)      :: a,squared
  optional        :: squared
  logical         :: squared,sq
  complex(kind=8) :: a(:)
  real(kind=8)    :: norm
  sq=.false.;if(present(squared))sq=squared
  norm=sum(conjg(a)*a)
  if(.not.sq)norm=sqrt(norm)
end function norm_c8d1
function norm_r8d2(a,squared) result(norm)
  implicit none
  intent(in)      :: a,squared
  optional        :: squared
  logical         :: squared,sq
  real(kind=8)    :: a(:,:),norm
  sq=.false.;if(present(squared))sq=squared
  norm=sum(a**2)
  if(.not.sq)norm=sqrt(norm)
end function norm_r8d2
function norm_c8d2(a,squared) result(norm)
  implicit none
  intent(in)      :: a,squared
  optional        :: squared
  logical         :: squared,sq
  complex(kind=8) :: a(:,:)
  real(kind=8)    :: norm
  sq=.false.;if(present(squared))sq=squared
  norm=sum(conjg(a)*a)
  if(.not.sq)norm=sqrt(norm)
end function norm_c8d2
function norm_r8d3(a,squared) result(norm)
  implicit none
  intent(in)      :: a,squared
  optional        :: squared
  logical         :: squared,sq
  real(kind=8)    :: a(:,:,:),norm
  sq=.false.;if(present(squared))sq=squared
  norm=sum(a**2)
  if(.not.sq)norm=sqrt(norm)
end function norm_r8d3
function norm_c8d3(a,squared) result(norm)
  implicit none
  intent(in)      :: a,squared
  optional        :: squared
  logical         :: squared,sq
  complex(kind=8) :: a(:,:,:)
  real(kind=8)    :: norm
  sq=.false.;if(present(squared))sq=squared
  norm=sum(conjg(a)*a)
  if(.not.sq)norm=sqrt(norm)
end function norm_c8d3
function norm_r8d4(a,squared) result(norm)
  implicit none
  intent(in)      :: a,squared
  optional        :: squared
  logical         :: squared,sq
  real(kind=8)    :: a(:,:,:,:),norm
  sq=.false.;if(present(squared))sq=squared
  norm=sum(a**2)
  if(.not.sq)norm=sqrt(norm)
end function norm_r8d4
function norm_c8d4(a,squared) result(norm)
  implicit none
  intent(in)      :: a,squared
  optional        :: squared
  logical         :: squared,sq
  complex(kind=8) :: a(:,:,:,:)
  real(kind=8)    :: norm
  sq=.false.;if(present(squared))sq=squared
  norm=sum(conjg(a)*a)
  if(.not.sq)norm=sqrt(norm)
end function norm_c8d4
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
function compare_date(n,dateA,dateB,wildcard) result(equal)
  implicit none
  integer,   intent(in)           :: n
  type(date),intent(in)           :: dateA,dateB(n)
  integer,   intent(in), optional :: wildcard
  logical :: equal
  integer :: dA(5),dB(5),i
  equal=.false.
  dA=(/dateA%year,dateA%month,dateA%day,dateA%hour,dateA%seconds/)
  do i=1,n
    dB=(/dateB(i)%year,dateB(i)%month,dateB(i)%day,&
         dateB(i)%hour,dateB(i)%seconds/)
    if(present(wildcard))then
      equal=equal.or.all((dA==dB).or.(dA==wildcard).or.(dB==wildcard))
    else
      equal=equal.or.all(dA==dB)
    endif
  enddo
end function compare_date
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
subroutine GetNCDim(ncFileID,name,num,varInt,varReal)
  implicit none
  integer(4), intent(in) :: ncFileID
  character(len=*), intent(in):: name
  integer, intent(inout) :: num
!   integer(4), dimension(:), pointer, intent(out), optional :: varInt
!   real(8),    dimension(:), pointer, intent(out), optional :: varReal
  integer(4), dimension(:), pointer, optional :: varInt
  real(8),    dimension(:), pointer, optional :: varReal
  integer(4) :: varID, dimID, ierr
  call nc_check(nf90_inq_dimid(ncFileID,trim(name),dimID))
  call nc_check(nf90_inquire_dimension(ncFileID,dimID,len=num))
  call nc_check(nf90_inq_varid(ncFileID,trim(name),varID))
  if(present(varInt))then
    if(associated(varInt))deallocate(varInt)
    allocate(varInt(num),stat=ierr)
    call nc_check(nf90_get_var(ncFileID,varID,varInt))
  endif
  if(present(varReal))then
    if(associated(varReal))deallocate(varReal)
    allocate(varReal(num),stat=ierr)
    call nc_check(nf90_get_var(ncFileID,varID,varReal))
  endif
end subroutine GetNCDim
subroutine PrintNCDim(ncFileID,dOut)
  implicit none
  integer(4), intent(in) :: ncFileID, dOut
  integer(4) :: i, j, ind
  real(8), dimension(:), pointer :: lon=>null(), lat=>null()
  character(len=64) :: fmtll
  call GetNCDim(ncFileID,'lon' ,numLon,varReal=lon)
  call GetNCDim(ncFileID,'lat' ,numLat,varReal=lat)
  call GetNCDim(ncFileID,'k'   ,numLev)
  call GetNCDim(ncFileID,'time',numRec)
  print "(/A,':',4(1X,A,':',I0))",'NC File: Dimensions',&
      'numLon',numLon,'numLat',numLat,'numLev',numLev,'numRec',numRec
  print "(/A,':')",'NC File: Lon/Lat Coordinates'
  fmtll="(LON(1X,I3,':',F6.2,LAT(1X,F6.2,:),/))"
  ind=index(fmtll,'LON' );write(fmtll(ind:ind+2),"(I3.3)")numLon/dOut+1
  ind=index(fmtll,'LAT' );write(fmtll(ind:ind+2),"(I3.3)")numLat/dOut+1
  print "((1X,A3,':',A6,200(1X,A3,I3.3,:)))",'N',"Lon",("Lat",j,j=1,numLat,dOut)
  print fmtll,(i,lon(i),(lat(j),j=1,numLat,dOut),i=1,numLon,dOut)
end subroutine PrintNCDim
subroutine AllocateNCVar(var2D,var3D,rec)
  implicit none
!   real(8), dimension(:,:,:),   pointer, intent(out), optional :: var2D
!   real(8), dimension(:,:,:,:), pointer, intent(out), optional :: var3D
  real(8), dimension(:,:,:),   pointer, optional :: var2D
  real(8), dimension(:,:,:,:), pointer, optional :: var3D
  integer(4), optional :: rec
  integer(4) :: rec0,rec1,ierr
  rec0=1;rec1=numRec
  if(present(rec))then
   rec0=rec;rec1=rec
  endif
!e.g. ps(lon,lat,time) ;
  if(present(var2D))then
    if(associated(var2D))deallocate(var2D)
    allocate(var2D(numLon,numLat,rec0:rec1),stat=ierr)
    var2D=-999.0
  endif
!e.g. no2(lon,lat,k,time) ;
  if(present(var3D))then
    if(associated(var3D))deallocate(var3D)
    allocate(var3D(numLon,numLat,numLev,rec0:rec1),stat=ierr)
    var3D=-999.0
  endif
end subroutine AllocateNCVar
subroutine GetNCVar(ncFileID,name,var2D,var3D,rec,unitconv)
  implicit none
  integer(4), intent(in) :: ncFileID
  character(len=*), intent(in):: name
!   real(8), dimension(:,:,:),   pointer, intent(out), optional :: var2D
!   real(8), dimension(:,:,:,:), pointer, intent(out), optional :: var3D
  real(8), dimension(:,:,:),   pointer, optional :: var2D
  real(8), dimension(:,:,:,:), pointer, optional :: var3D
  integer(4), intent(in), optional :: rec
  real(8),    intent(in), optional :: unitconv
  integer(4) :: varID, ierr, beg(4), cnt(4)
  real(8)    :: add_offset, scale_factor
  call AllocateNCVar(var2D=var2D,var3D=var3D,rec=rec)
  ierr=nf90_inq_varid(ncFileID,trim(name),varID)
  if(ierr/=nf90_noerr)then
    print *, "variable '",trim(name),"' does not exist: ",nf90_strerror(ierr)
    return
  endif
  ierr=nf90_get_att(ncFileID,VarID,'scale_factor',scale_factor)
  if(ierr/=nf90_noerr) scale_factor=1.0
  ierr=nf90_get_att(ncFileID,VarID,'add_offset',add_offset)
  if(ierr/=nf90_noerr) add_offset=0.0
! number of records to fetch
  beg=(/1,1,1,1/);cnt=(/numLon,numLat,numLev,numRec/)
  if(present(rec))then
    beg(4)=rec;cnt(4)=1
  endif
!e.g. ps(lon,lat,time) ;
  if(present(var2D))then
    call nc_check(nf90_get_var(ncFileID,varID,var2D,start=beg,count=cnt))
    var2D=var2D*scale_factor+add_offset
    if(present(unitconv))var2D=var2D*unitconv
  endif
!e.g. no2(lon,lat,k,time) ;
  if(present(var3D))then
    call nc_check(nf90_get_var(ncFileID,varID,var3D,start=beg,count=cnt))
    var3D=var3D*scale_factor+add_offset
    if(present(unitconv))var3D=var3D*unitconv
  endif
end subroutine GetNCVar
subroutine GetNCRec(ncFileID,ncDate,rec,exact)
  implicit none
  integer(4), intent(in) :: ncFileID
  type(date), intent(in) :: ncDate
  integer(4), intent(out):: rec
  logical,    intent(in),optional :: exact
  logical, parameter :: debug=.false.
! integer(4) :: ncTime
! integer(4), dimension(:), pointer :: time
  real(8) :: ncTime
  real(8), dimension(:), pointer :: time
  real,  parameter :: halfsec=1e0/(24*3600*2e0)
! call GetNCDim(ncFileID,'time',numRec,varInt=time)
  call GetNCDim(ncFileID,'time',numRec,varReal=time)
  call date2nctime(ncDate,ncTime)
  if(debug)print "(/A)","Looking for "//date2string("YYYY-MM-DD hh:mm",ncDate)
  recloop: do rec=1,numRec
    if(debug)print "(I4.4,': ',A)",rec,nctime2string("YYYY-MM-DD hh:mm",time(rec))
    if(abs(ncTime-time(rec))<halfsec) exit recloop
  enddo recloop
  if(rec>numRec.or.abs(ncTime-time(rec))>=halfsec)then
    if(present(exact)) &
      call CheckStop(exact,date2string('Date '//"YYYY-MM-DD hh:mm",ncDate)//&
                           ' not found on File!')
    rec=numRec
    print*,'WARNING: date '//date2string("YYYY-MM-DD hh:mm",ncDate)//&
           ' not found on File! Using lastrecord:',rec
  endif
endsubroutine GetNCRec
endmodule Util_ml

