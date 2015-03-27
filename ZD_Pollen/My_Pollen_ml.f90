!-----------------------------------------------------------------------!
! Empty Pollen rountines for "standrd" model compilation
!-----------------------------------------------------------------------!
! Birch pollen emission calculation based on
! M. Sofiev et al. 2006, doi:10.1007/s00484-006-0027-x
!
! Pollen emission based upon meteorology paparameters, and heatsum.
! Pollen particles are assumed of 22 um diameter and 800 kg/m3 density. 
!-----------------------------------------------------------------------!
module Pollen_const_ml
use PhysicalConstants_ml, only: PI
use ModelConstants_ml,    only: USE_POLLEN,DEBUG=>DEBUG_POLLEN
use ChemSpecs,            only: species_adv
use CheckStop_ml,         only: CheckStop
implicit none
public

real, parameter :: &
  D_POLL   = 22e-6,     & ! Pollen grain diameter [m]
  POLL_DENS= 800e3        ! Pollen density [g/m3]

real, parameter :: &
  grain_wt = POLL_DENS*PI*D_POLL**3/6.0, &  ! 1 grain weight [g]
  ug2grains= 1e-6/grain_wt                  ! # grains in 1 ug

contains

subroutine pollen_check(igrp)
  integer, intent(out), optional :: igrp
  logical,save :: first_call=.true.
  if(present(igrp))igrp=-1
  if(.not.first_call)return
  first_call=.false.
  call CheckStop(USE_POLLEN.or.DEBUG,&
    "USE_POLLEN/DEBUG_POLLEN on model compiled without pollen modules")
endsubroutine pollen_check
endmodule Pollen_const_ml
!-----------------------------------------------------------------------!
! Empty Pollen rountines for "standrd" model compilation
!-----------------------------------------------------------------------!
! Birch pollen emission calculation based on
! M. Sofiev et al. 2006, doi:10.1007/s00484-006-0027-x
!
! Pollen emission based upon meteorology paparameters, and heatsum.
! Pollen particles are assumed of 22 um diameter and 800 kg/m3 density. 
!-----------------------------------------------------------------------!
module Pollen_ml
  use Pollen_const_ml
  implicit none
  public:: pollen_flux,pollen_dump,pollen_read,pollen_check

  real,public,save, allocatable,dimension(:,:,:) :: &
    AreaPOLL,     & ! emission of pollen 
    heatsum,      & ! heatsum, needs to be remembered for forecast
    pollen_left     ! amount of pollen left in catkins, relative amount... 0:sr 1:end 

contains

subroutine pollen_flux(i,j,debug_flag) 
  implicit none
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug_flag
  call pollen_check()
endsubroutine pollen_flux

subroutine pollen_read()
  call pollen_check()
endsubroutine pollen_read

subroutine pollen_dump()
  call pollen_check()
endsubroutine pollen_dump
endmodule Pollen_ml
