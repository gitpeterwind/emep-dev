module AeroConstants_mod

  ! BoxChem/ESX/EMEP need some specific calculations of aerosol
  ! surface area, and we define 7 aerosol types
  ! We end up with variables here to avoid circularity in
  ! Makefile dependencies
   
  implicit none

  private 

  !A2018 - added to allow emepctm-like aerosol reactions
  !        so we can refer to eg AERO%PM_F

   integer, parameter, public :: NSAREA_DEF = 6 ! skip SIA_F - not needed!

   type, public :: aero_t
     ! EMEP only
     character(len=15) :: EQUILIB  ='EQSAM'  !'MARS ' !aerosol themodynamics 
     logical          :: DYNAMICS = .false.
     integer          :: NSIZE    = 7
     integer :: PM_F=1,SS_F=2,DU_F=3,SS_C=4,DU_C=5,PM=6  ! Will be set in GasParticleCoeffs_mod
   end type aero_t
   type(aero_t), public, save :: AERO = aero_t()

!M24 integer, parameter, public :: NSAREA_DEF = 7 ! skip ORIG=Riemer
!M24       SIA_F=1,PM_F=2,SS_F=3,DU_F=4,SS_C=5,DU_C=6,PM=7,NSAREA=NSAREA_DEF

end module AeroConstants_mod
