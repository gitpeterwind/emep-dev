module AeroConstants_mod

  ! BoxChem/ESX/EMEP need some specific calculations of aerosol
  ! surface area, and we define 7 aerosol types
  ! We end up with variables here to avoid circularity in
  ! Makefile dependencies

  ! To allow emepctm-like aerosol reactions so we can refer to eg AERO%PM_F:
   
  implicit none

  private 

   integer, parameter, public :: NSAREA_DEF = 9 ! skip SIA_F - not needed!

   type, public :: aero_t
     ! EMEP only
     character(len=15) :: EQUILIB  ='MARS'        ! 'ISORROPIA', 'EQSAM' or 'MARS' !aerosol thermodynamics 
     character(len=15) :: EQUILIB_WATER  = 'MARS' ! 'ISORROPIA', 'MARS' or 'EQSAM' !aerosol thermodynamics for PM water
     logical           :: DYNAMICS = .false.
     logical           :: INTERNALMIXED = .true.  ! sea salt assumption, only used by ISORROPIA and EQSAM
     logical           :: CATIONS = .true.        ! dust cat assumption, now only used by ISORROPIA
     real              :: RH_UPLIM_AERO = 0.98    ! RH upper limit used in thermodynamic equilibrium calls
     real              :: RH_LOLIM_AERO = 0.15    ! RH lower limit used in thermodynamic equilibrium calls to avoid div. by zero
     logical           :: ORGANIC_WATER = .false.  ! add organic matter water uptake to PM25 aerosol water
     logical           :: ThermoH2OSurfArea = .false.  ! calculate aerosol surf. area based on thermodynamics water uptake
     real              :: OM_KAPPA = 0.15         ! OM kappa hygroscopicity factor, 0.15 default from ISORROPIA II
     real              :: OM_RHO = 1400           ! aerosol density kg/m3; based on observations (Kakavas, 2023) & florou et al., 2014  
     integer           :: NSIZE = 7               ! can be removed?
     integer :: PM_F=1,SS_F=2,DU_F=3,SS_C=4,DU_C=5,SS_F_LS=6,SS_C_LS=7,PM_F_EQUI=8,PM=9  ! Will be set in GasParticleCoeffs_mod
     logical :: JUN21AERO = .false.   ! Flag to trigger ST's 2021 EQSAM and Aero tests
   end type aero_t
   type(aero_t), public, save :: AERO = aero_t()

end module AeroConstants_mod
