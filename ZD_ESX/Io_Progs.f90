
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!>MODULE  Io_Progs   - small help routines
!! Mini-veresion of code from EMEP model.

   module Io_Progs
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

  use Io_ml, only : IO_LOG 
  implicit none
  private

  public :: printlog     !> writes to both io number and screen 

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    write(*,*)  trim(txt)
    write(io,*)  trim(txt)
  end if
end subroutine PrintLog
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
                        end module Io_Progs
! __________________________________________________________________________ !
