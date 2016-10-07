module Pollen_const_ml
!-----------------------------------------------------------------------!
! Birch pollen emission calculation based on
! M. Sofiev et al. 2006, doi:10.1007/s00484-006-0027-x
!
! Pollen emission based upon meteorology paparameters, and heatsum.
! Pollen particles are assumed of 22 um diameter and 800 kg/m3 density. 
!-----------------------------------------------------------------------!
use PhysicalConstants_ml, only: PI,ATWAIR,AVOG
use ModelConstants_ml,    only: USE_POLLEN,DEBUG=>DEBUG_POLLEN
use CheckStop_ml,         only: CheckStop
use ChemSpecs,            only: NSPEC_ADV,NSPEC_SHL,species
use ChemGroups_ml,        only: chemgroups
use SmallUtils_ml,        only: find_index
implicit none
public

real, parameter  :: &
  T_cutoff_birch = 273.2+3.5, & ! Cut-off temperature [K] birch
  T_cutoff_olive = 273.2,     & ! Cut-off temperature [K] olive
  dH_d_birch    =  50,        & ! Flowering period [degree days] birch
  dH_d_olive    = 275,        & ! Flowering period [degree days] olive
  dH_birch=dH_d_birch*24*3600,& ! Flowering period [degree seconds] birch
  dH_olive=dH_d_olive*24*3600,& ! Flowering period [degree seconds] olive
  PREC_MIN = 0.0,             & ! Min cut-off precipitation [mm/h]
  PREC_MAX = 0.5,             & ! Max cut-off precipitation [mm/h]
  N_TOT_birch = 1.0e9,        & ! Available pollen [grains/m2] birch
  N_TOT_olive = 3.9e8,        & ! Available pollen [grains/m2] olive
  N_TOT_grass = 2.0e8,        & ! Available pollen [grains/m2] grass
  RH_LOW   = 0.50,            & ! Min cut-off relative humidity [1]
  RH_HIGH  = 0.80,            & ! Max cut-off relative humidity [1]
  PROB_IN_birch  = 0.2,       & ! Probability for flowering to start
  PROB_OUT_birch = 0.2,       & ! Probability for flowering to end 
  PROB_IN_olive  = 0.1,       & ! Probability for flowering to start
  PROB_OUT_olive = 0.1,       & ! Probability for flowering to end 
                                ! (could be assumed to be larger than PROB_IN)
  uncert_grass_day = 7,       &
  uncert_tot_grass_poll = 0.2,& ! end uncertainty for linear releases
! uncert_tot_grass_poll = 0.6,& ! end uncertainty for linear releases (new value)
  D_POLL_birch = 22.0,        & ! Pollen grain diameter [um] birch
  D_POLL_olive = 28.0,        & ! Pollen grain diameter [um] olive
  D_POLL_grass = 32.0,        & ! Pollen grain diameter [um] grass
  POLL_DENS= 800e3              ! Pollen density [g/m3]

! pollen arrays indexing, order must match with POLLEN_GROUP: birch,olive,grass
character(len=*), parameter :: &
  BIRCH = "POLLEN_BIRCH",&
  OLIVE = "POLLEN_OLIVE",&
  GRASS = "POLLEN_GRASS",&
  POLLEN_GROUP(3)=[BIRCH,OLIVE,GRASS]
integer, parameter :: &
  POLLEN_NUM=size(POLLEN_GROUP)
real, parameter  :: &
  N_TOT(POLLEN_NUM)=[N_TOT_birch,N_TOT_olive,N_TOT_grass]

real, parameter  :: &
  D_POLL(POLLEN_NUM)=[D_POLL_birch,D_POLL_olive,D_POLL_grass], & ! pollen diameter
  grain_wt(POLLEN_NUM) = POLL_DENS*PI*(D_POLL*1e-6)**3/6.0       ! 1 grain weight [g]
! weight 1 grain [ug], 1 mol of grains (AVOG*grain_wt) [Tonne=1e3 kg]
! BIRCH: 4.460e-3, 2686e6
! OLIVE: 9.195e-3, 5538e6
! GRASS: 13.73e-3, 8267e6

private :: N_TOT_birch,N_TOT_olive,N_TOT_grass,&
           D_POLL_birch,D_POLL_olive,D_POLL_grass

contains
subroutine pollen_check(igrp,uconv_adv)
  integer, intent(inout), optional :: igrp
  real, dimension(NSPEC_ADV), intent(inout), optional :: uconv_adv
  integer :: poll
  logical,save :: first_call=.true.
  poll=find_index("POLLEN",chemgroups(:)%name)
  if(present(igrp))igrp=poll
  if(.not.first_call)return
  first_call=.false.
  call CheckStop(USE_POLLEN.and.(poll<1),&
    "USE_POLLEN on model compiled without pollen")
  call CheckStop(DEBUG.and..not.USE_POLLEN,&
    "DEBUG_POLLEN on run without USE_POLLEN")
  call CheckStop(size(chemgroups(poll)%specs),size(POLLEN_GROUP),&
    "pollen_check: Inconsistent POLLEN group size")
  call CheckStop(any(species(chemgroups(poll)%specs)%name/=POLLEN_GROUP),&
    "pollen_check: Inconsistent POLLEN group species")
  if(present(uconv_adv))&
    uconv_adv(chemgroups(poll)%specs-NSPEC_SHL)=&
      uconv_adv(chemgroups(poll)%specs-NSPEC_SHL)/grain_wt
end subroutine pollen_check
end module Pollen_const_ml
