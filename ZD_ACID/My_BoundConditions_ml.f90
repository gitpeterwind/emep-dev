!############ ACID model ###################################################
module My_BoundConditions_ml
!____________________________________________________________________________
! This module specifies model-dependant:
!  (A) boundary-condition (bc) indices and 
!  (B) mapping arrays "bc2xn_adv" and "bc2xn_bgn" 
!       - which tell how  boundary conditions (bcs) are to be assigned to emep 
!         concentrations (xn_adv, xn_bgn). 
!
! THIS FILE WILL CHANGE FOR DIFFERENT CHEMISTRIES - MUST BE SUPPLIED BY USER. 
! AS A FIRST INDICATION OF THIS I HAVE SUPPLIED A "MY_MODEL" LABEL BELOW.
!
! The module BoundaryConditions_ml calls up this module with just:
!
!  call My_bcmap() 
!
! So far this module copes only with boundary condistion supplied either
! by the UiO global model (through the UiO_ml), or defined here as
! constant mixing ratios.
!
! Language: F
! History :
! ds - December 2000-January 2001
! hf - september 2001 Misc BC's added as a function of sigma
!____________________________________________________________________________
! IMPORTANT NOTE:
! The routines given here are constructed around the global model fields from 
! the University of Oslo (T21) global model. In order to use other models as 
! bcs  then usually these routines will have to be replaced by model-specific 
! routines. The important thing is that the inputs and outputs from the routine
! are independant of the global module used.
!u2 - "use" statements moved to top, stop_test replaced by gc_abort
!hf MADE comment:
! **/ Background species prescribed from solar zenith angle
! are set in Setup_1d_ml since it has to be reset each ? timestep anyway
! and is i,j dependendant
!_____________________________________________________________________________
!hf
  use GenSpec_bgn_ml, only: NSPEC_BGN,IXBGN_O3  !u3 ,IXBGN_H2O2
  use GenSpec_adv_ml, only: NSPEC_ADV   & ! No. advected species
                           ,IXADV_HNO3,IXADV_SO4,IXADV_PAN  &
                           ,IXADV_NO,IXADV_NO2, IXADV_SO2
  use GridValues_ml,   only : sigma_mid    !sigma layer midpoint
  use Met_ml         , only : z_mid        ! height of half layers
  use ModelConstants_ml, only: KMAX_MID    !hf Number of levels in vertical
  use Par_ml, only: NPROC,me
  !u3 use UiO_ml , only:   NGLOB_BC  & !  - indices from UiO model
  use GlobalBCs_ml , only:   NGLOB_BC  & !  - indices from global model
                      ,IBC_O3          & ! u3 - tmp
                      ,IBC_HNO3,IBC_PAN   &
                      ,IBC_SO2, IBC_SO4   &
                      ,IBC_NO,IBC_NO2    !u3   ,IBC_H2O2
 implicit none
 private

 !/-- subroutines
!hf h2o2
 public :: My_bcmap            ! sets bc2xn_adv, bc2xn_bc, and  misc_bc
           !u3 set_daily           !set daily bgn conc of h2o2
 
 !/-- model-type
 !    for consistency checks, possibly to match label in My_model_ml??
   character(len=12), public :: MY_MODEL = "uni-made"

   logical, public, parameter  :: BGN_2D = .true. !2d bgn species
   logical, private, parameter :: DEBUG_MYBC = .false.


 ! A. Set indices
 ! ===========================================================================
 ! For species which have constant mixing ratios:

!hf integer, public, parameter ::  NMISC_BC  =  1                  ! H2, OC
 integer, public, parameter ::  NMISC_BC  =  0
!hf integer, public, parameter :: IBC_H2  = NGLOB_BC + 1   !6s &
                              !6s ,IBC_OC  = NGLOB_BC + 2

 integer, public, parameter ::  NTOT_BC  = NGLOB_BC + NMISC_BC

 ! We also need the array misc_bc to specify concentrations of these species:

 real, public, save, dimension(NGLOB_BC+1:NTOT_BC,KMAX_MID) :: misc_bc 
!real, public, save, dimension(NGLOB_BC+1:NTOT_BC) :: misc_bc 

 ! B. Define mapping arrays
 ! ===========================================================================
 ! The mapping is done through the arrays bc2xn_adv and bc2xn_bgn, such that 
 ! the emep species are given along the x-dimension and the bc species along 
 ! the y.  e.g., the statement
 !
 !     bc2xn_adv(IBC_NOX,IXADV_NO2) = 0.55 
 !
 ! would assign the BC concentration of NOX  to the EMEP model concentration
 ! of NO2 after multiplication with a factor 0.55.
 !(The CTM2 concentration of NOx used as BC after a multilication 
 ! with a factor 0.55.)
 !_______________________________________________________________________

    real, public, save, dimension(NTOT_BC,NSPEC_ADV) :: bc2xn_adv  ! see above
    real, public, save, dimension(NTOT_BC,NSPEC_BGN) :: bc2xn_bgn ! see above
!u3 real, public, save ::  h2o2conc
 !-------
 contains
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 subroutine My_bcmap()    ! sets bc2xn_adv, bc2xn_bc, and  misc_bc
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    real    :: ppb = 1.0e-9
    integer :: ii,i,j,k 
    real :: decrease_factor(NGLOB_BC+1:NTOT_BC) ! Decrease factor for misc bc's
                                      ! Gives the factor for how much of 
                                      ! the top-layer conc. that is left 
                                      ! at bottom layer

    real :: top_misc_bc(NGLOB_BC+1:NTOT_BC) ! Conc. at top of misc bc
!    real :: ratio_length(KMAX_MID)    ! Vertical length of the actual layer
                                      ! divided by length from midpoint of 
                                      ! layer 1 to layer KMAX_MID

!hf MADE    misc_bc            = 0.0  ! Initialise
    bc2xn_adv          = 0.0  ! Initialise
    bc2xn_bgn          = 0.0  ! Initialise

    ! Own (constant mixing ratio) boundary conditions **********

    ! NOTE - these species have to have the bc2xn_ indices set to 1.0 for either
    ! the advected or the background concentrations, in order that the 
    ! concentrations specified in misc_bc are transferred correctly into the 
    ! boundary conditions.
    !
    ! 18.09.01 -hf- misc bc'c as function of sigma possible 

        !! misc_bc(IBC_CH4) = 1760.0 * ppb   !! used for MACHO-jej chem
!        top_misc_bc(IBC_H2) = 600.0 * ppb
        !! top_misc_bc(IBC_OC) =   0.0 !!! 1.0 * ppb

!        decrease_factor(IBC_H2)=0.0 !TESTING
        !! decrease_factor(IBC_OC)=0.0 !No increase/decrease with height


!a) Function of height(not included yet!):
!ratio_length(i,j,k)=(z_mid(i,j,1)-z_mid(i,j,k))/(z_mid(i,j,1)-z_mid(i,j,KMAX_MID))
!Replace sigma_mid with ratio_length and make misc_bc 4 dimentional
!
!b)Function of sigma_mid (I assume that top_misc_bc is given for
!top(sigma_bnd(1)=0.), and that at ground (sigma_bnd(KMAX_BND)=1.) the conc.
!is top_misc_bc -decrease_factor*top_misc_bc. Since I choose to set the concentration
! as a factor of sigma_mid, the concentration in the lowest grid cell will not be 
! excactly  top_misc_bc -decrease_factor*top_misc_bc, but close.

        do ii=NGLOB_BC+1,NTOT_BC
           do k=1,KMAX_MID
!              misc_bc(ii,k)=top_misc_bc(ii) - &
!                  top_misc_bc(ii)*decrease_factor(ii)*sigma_mid(k) 
!              if (me.eq.0) then
!              write(*,*)'height,misc_vert,k', sigma_mid(k),misc_bc(ii,k),k
!              endif
           enddo
        enddo

!        misc_bc(IBC_H2)     = 600.0 * ppb
!        misc_bc(IBC_OC)     =   1.0*ppb !!! 0.0

!        bc2xn_adv(IBC_H2,  IXADV_H2)    = 1.0
        !! bc2xn_adv(IBC_OC,  IXADV_OC)    = 1.0

        !/-- check, just in case we forgot something...!

        if ( DEBUG_MYBC  ) then
             print *, "In My_bcmap, NGLOB_BC is", NGLOB_BC
             print *, "In My_bcmap, NTOT_BC  is", NTOT_BC
             do i = NGLOB_BC+1 , NTOT_BC
                print *, "In My_bcmap, sum-adv", i, " is", sum(bc2xn_adv(i,:))
                print *, "In My_bcmap, sum-bgn", i, " is", sum(bc2xn_bgn(i,:))
             end do
        end if ! DEBUG

        do i = NGLOB_BC+1 , NTOT_BC
           if ( sum(bc2xn_adv(i,:))  + &
                sum(bc2xn_bgn(i,:))      /= 1.0 )          &
                !u2 call stop_test(.true.,me,NPROC,99,"BC problem - my")
                call gc_abort(me,NPROC,"BC problem - my")
        end do


 ! mappings for species from global model ***********

  bc2xn_adv(IBC_HNO3    ,IXADV_HNO3    )   =   1.0 
!u3 bc2xn_adv(IBC_PANX    ,IXADV_PAN     )   =   1.0   ! ds
  bc2xn_adv(IBC_PAN     ,IXADV_PAN     )   =   1.0 
  bc2xn_adv(IBC_NO      ,IXADV_NO      )   =   1.0
  bc2xn_adv(IBC_NO2     ,IXADV_NO2     )   =   1.0
  bc2xn_adv(IBC_SO2     ,IXADV_SO2     )   =   1.0
  bc2xn_adv(IBC_SO4     ,IXADV_SO4     )   =   0.0

  
! The following species are excluded either because they have no corresponding
! species in the emep model, or because they have lifetimes which are so
! short that initialisation is uncessary.
!-----------------------------------------------------------------------------
!!bc2xn_adv(IBC_NOX     ,IXADV_NOX     )   =  -1.0   ! Excluded, we have NO and NO2


!/** mappings for bgn species from global model ***********
  bc2xn_bgn(IBC_O3      ,IXBGN_O3      )   =   1.0
!hf h2o2  bc2xn_bgn(IBC_H2O2    ,IXBGN_H2O2    )   =   1.0

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 end subroutine My_bcmap
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !u3 set_daily replaced by Daily_halfsine function in Tabulations and
 !u3 Init_mychem in My_Chem_ml.
 !u3 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !u3 subroutine set_daily
 !u3 !hf Set daily bgn conc. of h2o2
 !u3
 !u3 use ModelConstants_ml,    only : current_date
 !u3 use PhysicalConstants_ml, only : PI
 !u3
 !u3
 !u3 real dd 
 !u3
 !u3 dd = current_date%day
 !u3 h2o2conc = 0.35e-9 + 0.3e-9 * sin (PI * ((365. - dd)/365.))
 !u3 end subroutine set_daily
 !u3 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module My_BoundConditions_ml
!_____________________________________________________________________________












