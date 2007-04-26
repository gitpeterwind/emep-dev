
module CheckStop_ml

! Provides routines to check for errors and if necessary close
! down the code neatly (** all processors **). 

! The generic routine CheckStopAll is defined, so that the code may be 
! stopped if:
!   (a)  errmsg  /= ok
!   (b)  int     /= 0               (e.g. iostat index after read)
!   (c)  int1    /= int2
!   (d)  logical expression = true  (e.g. lu < 0 for landuse index)

! Dave Simpson, 25 April 2007
!

 implicit none
 INCLUDE 'mpif.h'
 INTEGER, private :: STATUS(MPI_STATUS_SIZE),INFO

  public :: StopAll
  public :: CheckStop
  private :: CheckStop_ok, CheckStop_int1, CheckStop_int2, CheckStop_TF

  interface CheckStop
     module procedure CheckStop_ok
     module procedure CheckStop_okinfo
     module procedure CheckStop_int1
     module procedure CheckStop_int2
     module procedure CheckStop_TF   
  end interface CheckStop

 contains

   subroutine StopAll(errmsg)
      character(len=*), intent(in) :: errmsg

    ! Stops all processors.
    ! MPI_COMM_WORLD indicates all processors, in other programs you could have
    ! different groups of processes.
    ! INFO is the error message from MPI

      if ( errmsg /= "ok" ) then
        write(*,*) "StopAll Called with", errmsg
        write(*,*) "MPI_ABORT!!"
        call MPI_ABORT(MPI_COMM_WORLD,9,INFO)
      end if
  end subroutine StopAll


  !---- Four variations on CheckStop:

  subroutine CheckStop_ok(errmsg)    ! Test if errmsg /= "ok"
      character(len=*), intent(in) :: errmsg

      if ( errmsg /= "ok" ) then
        write(*,*) "CheckStop_ok Called with:  errmsg ", errmsg
        call StopAll(errmsg)
      end if
  end subroutine CheckStop_ok

  subroutine CheckStop_okinfo(errmsg,infomsg)    ! Test if errmsg /= "ok"
      character(len=*), intent(in) :: errmsg
      character(len=*), intent(in) :: infomsg

      if ( errmsg /= "ok" ) then
        write(*,*) "CheckStop_ok Called with:  errmsg ", errmsg
        write(*,*) "                          infomsg ", infomsg
        call StopAll(errmsg)
      end if
  end subroutine CheckStop_okinfo

  subroutine CheckStop_int1(int1,infomsg)    ! Test if int1 /= 0
      integer, intent(in)          :: int1
      character(len=*), intent(in) :: infomsg

      if ( int1 /= 0 ) then
        write(*,*) "CheckStopl_int1 Called with:    int1 ", int1
        write(*,*) "                             infomsg ", infomsg
        call StopAll(infomsg)
      end if
  end subroutine CheckStop_int1

  subroutine CheckStop_int2(int1,int2, infomsg)   ! Test if int1 /= int2
      integer, intent(in)          :: int1, int2
      character(len=*), intent(in) :: infomsg

      if ( int1 /= int2 ) then
        write(*,*) "CheckStopl_int2 Called with: int1 ", int1, " int2 ", int2
        write(*,*) "                             infomsg ", infomsg
        call StopAll(infomsg)
      end if
  end subroutine CheckStop_int2

  subroutine CheckStop_TF(is_true, infomsg)   ! Test expression, e.g. lu<0
      logical, intent(in)          :: is_true  
      character(len=*), intent(in) :: infomsg

      if ( is_true ) then
        write(*,*) "CheckStopl_TF   Called with: logical ", is_true
        write(*,*) "                             infomsg ", infomsg
        call StopAll(infomsg)
      end if
  end subroutine CheckStop_TF

end module CheckStop_ml

