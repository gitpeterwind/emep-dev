module CheckStop_ml
! Provides routines to check for errors and if necessary close
! down the code neatly (** all processors **). 

! The generic routine CheckStopAll is defined, so that the code may be 
! stopped if:
!   (a)  errmsg   /= ok
!   (b)  int      /= 0               (e.g. iostat index after read)
!   (c)  int1     /= int2
!   (d)  string1  /= string2
!   (e)  logical  expression = true  (e.g. lu < 0 for landuse index)
!   (f)  range    outside [range(0)..range(1)]

implicit none
INCLUDE 'mpif.h'
INTEGER, private :: STATUS(MPI_STATUS_SIZE),INFO

public  :: StopAll, CheckStop
private :: CheckStop_ok, CheckStop_okinfo, CheckStop_int1, CheckStop_int2, &
           CheckStop_str2, CheckStop_TF, CheckStop_range

interface CheckStop
   module procedure CheckStop_ok
   module procedure CheckStop_okinfo
   module procedure CheckStop_int1
   module procedure CheckStop_int2
   module procedure CheckStop_str2
   module procedure CheckStop_TF
   module procedure CheckStop_range
end interface CheckStop

contains

subroutine StopAll(errmsg)
  character(len=*), intent(in) :: errmsg
  ! Stops all processors.
  ! MPI_COMM_WORLD indicates all processors, in other programs you could have
  ! different groups of processes.
  ! INFO is the error message from MPI

  if(errmsg/="ok") then
    write(*,*) "STOP-ALL ERROR: ", trim(errmsg)
    call MPI_ABORT(MPI_COMM_WORLD,9,INFO)
  endif
endsubroutine StopAll

!---- Variations on CheckStop:
subroutine CheckStop_ok(errmsg)                 ! Test if errmsg /= "ok"
  character(len=*), intent(in) :: errmsg

  if(errmsg/="ok") then
   !write(*,*) "CheckStop_ok Called with:  errmsg ", errmsg
    call StopAll(errmsg)
  endif
endsubroutine CheckStop_ok

subroutine CheckStop_okinfo(errmsg,infomsg)     ! Test if errmsg /= "ok"
  character(len=*), intent(in) :: errmsg
  character(len=*), intent(in) :: infomsg

  if(errmsg/="ok") then
   !write(*,*) "CheckStop_ok Called with:  errmsg ", errmsg
    write(*,*) "                          infomsg ", infomsg
    call StopAll(errmsg)
  endif
endsubroutine CheckStop_okinfo

subroutine CheckStop_int1(int1,infomsg)         ! Test if int1 /= 0
  integer, intent(in)          :: int1
  character(len=*), intent(in) :: infomsg

  if(int1/=0) then
    write(*,*) "CheckStopl_int1 Called with:    int1 ", int1
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  endif
endsubroutine CheckStop_int1

subroutine CheckStop_int2(int1,int2, infomsg)   ! Test if int1 /= int2
  integer, intent(in)          :: int1, int2
  character(len=*), intent(in) :: infomsg

  if(int1/=int2) then
    write(*,*) "CheckStopl_int2 Called with: int1 ", int1, " int2 ", int2
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  endif
endsubroutine CheckStop_int2

subroutine CheckStop_str2(str1,str2, infomsg)   ! Test if str1 /= str2
  character(len=*), intent(in) :: str1, str2, infomsg

  if(trim(str1)/=trim(str2)) then
    write(*,*) "CheckStopl_str2 Called with: str1 ", str1, " str2 ", str2
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  endif
endsubroutine CheckStop_str2

subroutine CheckStop_TF(is_error, infomsg)   ! Test expression, e.g. lu<0
  logical, intent(in)          :: is_error  
  character(len=*), intent(in) :: infomsg

  if(is_error) then
   !write(*,*) "CheckStopl_TF   Called with: logical ", is_error
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  endif
endsubroutine CheckStop_TF

subroutine CheckStop_range(var,vrange,infomsg)  ! test .not.(vrange(0)<=var<=vrange(1))
  real, intent(in) :: var,vrange(0:1)
  character(len=*), intent(in) :: infomsg
  character(len=*), parameter :: &
    errfmt="(A,'=',ES12.3,' is out of range ',ES12.3,'..',F6.2)"

  if(var<vrange(0).or.var>vrange(1))then
    write(*,errfmt) "CheckStopl_range: variable",var,vrange
   !write(*,*) "                             infomsg ", infomsg
    call StopAll(infomsg)
  endif
endsubroutine CheckStop_range

endmodule CheckStop_ml

