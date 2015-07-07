module Chemfields_ml
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
use AllocInits,           only: AllocInit
use ChemSpecs,            only: NSPEC_ADV, NSPEC_SHL, NSPEC_TOT ! => No. species 
use ModelConstants_ml,    only: KMAX_MID, KCHEMTOP, AERO        ! =>  z dimension
use NumberConstants,      only: UNDEF_R
use Par_ml,               only: MAXLIMAX,MAXLJMAX   ! => x, y dimensions
use Setup_1dfields_ml
implicit none
private

!-------- this snipppet was from older GenSpec_bgn_ml. ------
  ! PRETTY MUCH FAKED FOR NOW. CAN BE DELETED SOON IN HOPE!
  !+ Defines indices and NSPEC for bgn : Background species

   ! Species which can be specified simply for each column, e.g.
   ! as function of local meteorology or zenith angle
   !   o2, m,  and for MADE-like, oh, ch3coo2

   integer, public, parameter ::  NSPEC_BGN = 0 ! No. 3D bgn species
   integer, public, parameter ::  NSPEC_COL = 0 ! total no. prescribed specs

    !/ define xn_2d_bgn here.
     real, public, save, allocatable, dimension(:,:) :: xn_2d_bgn

!-------- end of this snipppet from older GenSpec_bgn_ml. ------

    !----------------- basic chemical fields ----------------------------------!
    !  Here we declare and initialise to zero the chemical fields used in the  !
    !  model, as well as cfac (converts from 50m to 1m/3m output)         ! 
    !---------------------------------------------------------------------!

  real, save, allocatable, public :: &
     xn_adv(:,:,:,:)  &
    ,xn_shl(:,:,:,:)  &
    ,xn_bgn(:,:,:,:) &
    ,PM25_water(:,:,:) &  !3D PM water
    ,PM25_water_rh50(:,:) &   !gravimetric PM water
    ,Gerber_water(:,:,:)   !3D PM water from GERBER

  real, save, allocatable, public :: &
    SurfArea_um2cm3(:,:,:)  !2D  aerosol surface area, um2/cm3  (n,i,j)

  real, public, save, allocatable:: Fgas3d (:,:,:,:)  ! for SOA

  real, public, save, allocatable :: AOD(:,:,:,:),Extin_coeff(:,:,:,:,:)

  real, save, allocatable, public :: &
     cfac   (:,:,:)   

  real, save, allocatable, public :: &
     so2nh3_24hr(:,:)!hf CoDep

  real, save, allocatable, public :: &
     Grid_snow(:,:) !snow_flag fraction in grid

  public ::alloc_ChemFields

contains

  subroutine alloc_ChemFields

    implicit none
    integer :: nk

    allocate(xn_adv(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID))
    xn_adv=0.0
    allocate(xn_shl(NSPEC_SHL,MAXLIMAX,MAXLJMAX,KMAX_MID))
    xn_shl=0.0
    allocate(xn_bgn(NSPEC_BGN,MAXLIMAX,MAXLJMAX,KMAX_MID))
    xn_bgn=0.0
    allocate(PM25_water(MAXLIMAX,MAXLJMAX,KMAX_MID))
    PM25_water=0.0
    allocate(PM25_water_rh50(MAXLIMAX,MAXLJMAX))
    PM25_water_rh50=0.0
!   allocate(AOD(MAXLIMAX,MAXLJMAX))
!   AOD=0.0
!   allocate(Extin_coeff(MAXLIMAX,MAXLJMAX,KMAX_MID))
!   Extin_coeff=0.0
    allocate(cfac(NSPEC_ADV,MAXLIMAX,MAXLJMAX))
    cfac=1.0
    allocate(so2nh3_24hr(MAXLIMAX,MAXLJMAX))
    so2nh3_24hr=0.0
    allocate(Grid_snow(MAXLIMAX,MAXLJMAX))
    Grid_snow=0.0
    allocate(xn_2d_bgn(1,KCHEMTOP:KMAX_MID))

    allocate(xn_2d(NSPEC_TOT,KCHEMTOP:KMAX_MID))
  xn_2d = 0.0
 !   nk = KMAX_MID-KCHEMTOP+1   ! number of levels used in column chemistry
 !   call AllocInit(xn_2d,0.0, NSPEC_TOT, nk)
    allocate(Fgas(NSPEC_TOT,KCHEMTOP:KMAX_MID),Fpart(NSPEC_TOT,KCHEMTOP:KMAX_MID))
    Fgas  = 1.0! Fraction as gas-phase
    Fpart = 0.0
    allocate(rcemis(NSPEC_SHL+1:NSPEC_TOT,KCHEMTOP:KMAX_MID))
    rcemis = 0.0
    allocate(rh(KCHEMTOP:KMAX_MID),amk(KCHEMTOP:KMAX_MID),o2(KCHEMTOP:KMAX_MID))
    allocate(n2(KCHEMTOP:KMAX_MID),h2o(KCHEMTOP:KMAX_MID),temp(KCHEMTOP:KMAX_MID))
    allocate(tinv(KCHEMTOP:KMAX_MID),pp(KCHEMTOP:KMAX_MID))
    allocate(itemp(KCHEMTOP:KMAX_MID))
    CHEMSIZE = KMAX_MID-KCHEMTOP+1

   ! Surface area and water
!if( USES%SURF_AREA) then
    allocate(surfarea_um2cm3(AERO%NSAREA,MAXLIMAX,MAXLJMAX))
    SurfArea_um2cm3=0.0
    allocate(Gerber_water(MAXLIMAX,MAXLJMAX,KMAX_MID))
    Gerber_water=0.0

   ! wet DpgN and defaults from dry values
    allocate(DpgNw(AERO%NSAREA, KCHEMTOP:KMAX_MID))

    allocate(S_m2m3(AERO%NSAREA, KCHEMTOP:KMAX_MID)) ! GERBER
    S_m2m3=0.0

    ! Mol speeds
    allocate(cn2o5(KCHEMTOP:KMAX_MID),chno3(KCHEMTOP:KMAX_MID),&
              cho2(KCHEMTOP:KMAX_MID),co3(KCHEMTOP:KMAX_MID))
    cn2o5=UNDEF_R
    chno3=UNDEF_R
    cho2= UNDEF_R
    co3=  UNDEF_R
    allocate(aero_fom(KCHEMTOP:KMAX_MID),aero_fdust(KCHEMTOP:KMAX_MID),&
              aero_fss(KCHEMTOP:KMAX_MID))
    aero_fom    = UNDEF_R
    aero_fss    = UNDEF_R
    aero_fdust  = UNDEF_R
  
!end if


  end subroutine alloc_ChemFields


!_____________________________________________________________________________
endmodule Chemfields_ml
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
