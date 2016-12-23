!###############################################################################
!
! Fortran file units.
!
!###############################################################################
!
#include "go.inc"
!
!###############################################################################

module GO_FU

  implicit none
  
  
  ! --- in/out ------------------------------
  
  private
  
  public  ::  goStdIn, goStdOut, goStdErr
  public  ::  goFuRange
  
  public  ::  goGetFU

  
  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'GO_FU'

  ! standard file units
  integer, parameter    ::  goStdErr = 0
  integer, parameter    ::  goStdIn  = 5
  integer, parameter    ::  goStdOut = 6
  
  ! range of file units that might be used by this program:
  integer, parameter    ::  goFuRange(2) = (/200,999/)


contains

 
  ! ====================================================================
  
  
  ! Return the first free available file unit number.

  subroutine goGetFU( fu, status )
  
    ! --- in/out --------------------------

    integer, intent(out)      ::  fu
    integer, intent(out)      ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/goGetFU'
    
    ! --- local --------------------------

    integer               ::  i
    character(len=256)    ::  fname
    logical               ::  opened

    ! --- local ---------------------------
    
    ! start with lowest possible unit:
    fu = goFuRange(1) - 1
    
    ! loop until unopned unit is found:
    do

      ! try next file unit:
      fu = fu + 1

      ! too large ?
      if ( fu > goFuRange(2) ) then
        write (*,'("unable to select free file unit within allowed range ...")')
        write (*,'("close some files or increase goFuRange in module GO_FU")')
        write (*,'("current goFuRange : ",i6," .. ",i6)') goFuRange
        write (*,'("open files:")')
        do i = goFuRange(1), goFuRange(2)
          inquire( unit=i, name=fname )
          write (*,'(i6," : ",a)') i, trim(fname)
        end do
        write (*,'("in ",a)') rname; status=1; return
      end if
      
      ! skip ?
      if ( fu==goStdIn  ) cycle
      if ( fu==goStdOut ) cycle
      if ( fu==goStdErr ) cycle

      ! free available unit ? then ok
      inquire( unit=fu, opened=opened )
      if ( .not. opened ) exit
      
    end do

    ! ok
    status = 0

  end subroutine goGetFU


end module GO_FU


!*********************************************
! Test program to identify unit numbers for standard input etc.
! Compile:
!   f90 -o test.exe test.f90
! Run and see what comes to the terminal:
!   test.exe
! Files named 'fort.<fu>' will be created for not-special units.
! Error message when writing to standard input;
! change the range of file units, uncomment the test lines,
! and execute:
!   echo 'hello' | test.exe
! Change the range of file units if nothing happends.
!
!program test_fu
!
!  integer            ::  fu
!  character(len=10)  ::  s
!  
!  do fu = 0, 10
!    write (*,*) 'try to write to file unit ', fu
!    write (fu,*) 'THIS IS FILE UNIT ', fu
!  end do
!
!  ! uncoment following line to check standard input:
!  !read (5,*) s
!  !print *, 'read from standard input: ', s
!  
!end program test_fu
!
!*********************************************
