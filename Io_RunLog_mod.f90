module Io_RunLog_mod
!_____________________________________________________________________________
! -- routines to write out to both screen and RunLog file.
! -- keep simple to avoid circularity
!_____________________________________________________________________________
use Io_Nums_mod,             only: IO_LOG
implicit none


 public :: PrintLog                   !writes message to both RunLog and unit 6

 character(len=200), public :: logtxt  !long text string used to convert mixes
                                      ! of str, float to str

contains
!-------------------------------------------------------------------------
subroutine PrintLog(txt,OutputProc,ioOption)
  character(len=*), intent(in) :: txt
  logical, intent(in), optional :: OutputProc  !typically MasterProc, me==0
  integer, intent(in), optional :: ioOption    !use for other files
  logical :: ok2print
  integer :: io
  ok2print = .true.
  if ( present(OutputProc) ) ok2print = OutputProc
  if ( ok2print) then
    io = IO_LOG
    if ( present(ioOption) ) io = ioOption
    write(*,fmt='(A)')   trim(txt)
    write(io,fmt='(A)')  trim(txt)
  end if
end subroutine PrintLog
!-------------------------------------------------------------------------

end module Io_RunLog_mod
