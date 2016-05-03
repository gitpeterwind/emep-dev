module DA_ml
implicit none
logical, parameter ::     &
  DEBUG_DA_1STEP=.false., &   ! run only 1 DA step (no adv/chem)
  DEBUG_DA_OUTPUT=.false.     ! hourly output before/after DA step
endmodule DA_ml
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module DA_3DVar_ml
use CheckStop_ml,     only: CheckStop
use ModelConstants_ml,only: ANALYSIS
implicit none
integer, parameter :: NTIMING_3DVAR=0, T_3DVAR=0
contains
subroutine main_3dvar()
!-----------------------------------------------------------------------
! Empty call to 3dvar, for "standrd" model compilation
!-----------------------------------------------------------------------
implicit none
logical, save :: first_call=.true.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  if(.not.first_call)return
  call CheckStop(ANALYSIS,&
    "No 3DVar available. Need to recompile, e.g. make MACC-3DVar")
  first_call=.false.
endsubroutine main_3dvar
endmodule DA_3DVar_ml
