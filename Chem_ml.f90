!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module Chemfields_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
use Par_ml               , only: MAXLIMAX,MAXLJMAX   ! => x, y dimensions
use ModelConstants_ml    , only: KMAX_MID     ! =>  z dimension
use GenSpec_adv_ml,  only: NSPEC_ADV         ! => No. species 
use GenSpec_shl_ml,  only: NSPEC_SHL         ! => No. species 
!u2 !hf
!u2 use GenSpec_bgn_ml,  only: NSPEC_BGN,NSPEC_UIO         ! => No. species 
use GenSpec_bgn_ml,  only: NSPEC_BGN         ! => No. species 
implicit none
private

    !----------------- basic chemical fields ----------------------------------!
    !  Here we declare and initialise to zero the chemical fields used in the  !
    !  model                                                                   !
    !  and cfac which is ....
    !  From eulmc.inc                     
    !---------------------------------------------------------------------!

  real, save, public :: &
     xn_adv (NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0     &
    ,xn_shl (NSPEC_SHL,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0     &
    ,xn_bgn (NSPEC_BGN,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0   
!u2    ,xn_bgn (NSPEC_UIO,MAXLIMAX,MAXLJMAX,KMAX_MID)   = 0.0   

  real, save, public :: &
     cfac   (NSPEC_ADV,MAXLIMAX,MAXLJMAX) = 1.0    

!_____________________________________________________________________________
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
                     end module Chemfields_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!_____________________________________________________________________________
