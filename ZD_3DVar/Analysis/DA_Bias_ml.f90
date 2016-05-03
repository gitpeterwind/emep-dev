#define STRING2(x) #x
#define STRING(x) STRING2(x)
#define HERE(MSG) MSG//" ("//__FILE__//":"//STRING(__LINE__)//")."
module DA_Bias_ml
use CheckStop_ml,     only: CheckStop
use Io_ml,            only: IO_TMP,PrintLog
use SmallUtils_ml,    only: find_index
use TimeDate_ml,      only: date,current_date
use TimeDate_ExtraUtil_ml, only: date2string
use Units_ml,         only: Units_Scale,Group_Units
use DA_ml,            only: debug=>DEBUG_DA,dafmt=>da_fmt_msg,damsg=>da_msg
use DA_Obs_ml,        only: nobsData,nobsDataMax,obsData
!-----------------------------------------------------------------------
implicit none
integer, parameter :: NBIAS_PREDICTORS=0
!-----------------------------------------------------------------------
type, private :: bias_data
  real :: beta(0:NBIAS_PREDICTORS)=0e0,stddev=1e3,obsFrac=1e2,weight=1e6
endtype
integer :: nbiasData
type(bias_data), target, save :: biasData(nobsDataMax)
namelist /BIAS_PREDICTOR/nbiasData,biasData
!-----------------------------------------------------------------------
real, dimension(:), pointer :: bias=>null(),biasStdDev=>null(),&
                               biasOFrac=>null(),biasWeight=>null()
real, dimension(:), allocatable, save :: dbias
!-----------------------------------------------------------------------
contains
subroutine allocate_bias()
  integer :: ierr,ipar
  call CheckStop(nbiasData,nobsData,HERE('Inconsistent nbiasData'))
  call CheckStop(any(.not.obsData(:nobsData)%set),HERE('.not.obsData%set'))
  if(.not.associated(bias))       &
    bias=>biasData(:nobsData)%beta(0)
  if(.not.associated(biasStdDev)) &
    biasStdDev=>biasData(:nobsData)%stddev
  if(.not.associated(biasOFrac))      &
    biasOFrac=>biasData(:nobsData)%obsFrac
  if(.not.associated(biasWeight))then  
    biasWeight=>biasData(:nobsData)%weight
    biasWeight=obsData(:nobsData)%nobs*biasOFrac/100
  endif
  if(.not.allocated(dbias))then
    allocate(dbias(nobsData),stat=ierr)
    call CheckStop(ierr,HERE('Allocate DBIAS'))
    dbias=0.0
  endif
  do ipar=1,nobsData
    if(obsData(ipar)%set)then
      print dafmt,"Bias-info for "//trim(obsData(ipar)%tag)
    endif
  enddo
endsubroutine
subroutine deallocate_bias()
  if(associated(bias))      nullify(bias)
  if(associated(biasStdDev))nullify(biasStdDev)
  if(associated(biasOFrac)) nullify(biasOFrac)
  if(associated(biasWeight)) nullify(biasWeight)
  if(allocated(dbias))      deallocate(dbias)
endsubroutine
subroutine write_bias(fname)
  character(len=*), intent(in), optional :: fname
  integer :: ierr
  if(present(fname))then
    open(unit=IO_TMP,file=date2string(fname,current_date),&
      status='REPLACE',action='WRITE',form='FORMATTED',iostat=ierr)
  else
    open(unit=IO_TMP,file='namelist.nml',&
      status='OLD',action='WRITE',position='APPEND',form='FORMATTED',iostat=ierr)
  endif
  call CheckStop(ierr,HERE('open namelist: BIAS_PREDICTOR'))
  write(unit=IO_TMP,nml=BIAS_PREDICTOR,iostat=ierr)
  call CheckStop(ierr,HERE('write namelist: BIAS_PREDICTOR'))
  close(IO_TMP)
endsubroutine
endmodule DA_Bias_ml
