module My_UKDep_ml    ! DryDep_ml
!+
! Module to define the dry deposition components and rates. We
! define the min (vd_min) and max dep. velocity (Vg) first and then derive the
! daytime addition (vd_day).
!/**************************************************************************
!  Specifies which of the possible species (from Wesely's list)
!  are required in the current air pollution model   
!/**************************************************************************

 use My_Derived_ml , only : DDEP_SOX,DDEP_OXN,DDEP_RDN, &
                             IOU_INST    &!updates inst. dep. fields
                           , ddep         ! 2d fields
 use GenSpec_adv_ml , only: NSPEC_ADV &
                   ,IXADV_SO2,IXADV_NO2  &
                   ,IXADV_HNO3 &  !! ,IXADV_EC  ,IXADV_OC & 
    ,  IXADV_PAN ,IXADV_SO4,IXADV_NH3,IXADV_AMNI, IXADV_AMSU
 use ModelConstants_ml , only : atwS, atwN
 use Wesely_ml
 implicit none
 private

  public :: Init_vd
  public :: Add_ddep


  !/** Variables used in deposition calculations
 
  ! DDEP_xx gives the index that will be used in the EMEP model
  ! WES_xx gives the index of the Wesely gas to which this corresponds


  ! Here we define the minimum set of species which has different
  ! deposition velocities. We calculate Vg for these, and then
  ! can use the rates for other similar species. (e.g. AMSU can use
  ! the Vg for SO4.  Must set NDRYDEP_CALC species

  !/** IMPORTANT: the variables below must match up in the sense that, for 
  ! example, if DDEP_NH3=4 then the 4th element of DRYDEP must be WES_NH3.

  integer, public, parameter :: NDRYDEP_CALC = 6

  integer, public, parameter :: &
       CDEP_HNO3 = 1, CDEP_O3  = 2, CDEP_SO2 = 3  &
      ,CDEP_NH3  = 4, CDEP_NO2 = 5, CDEP_PAN = 6

  integer, public, parameter :: CDEP_SET = -99    

!also have CDEP_H2O2=, CDEP_ALD, CDEP_HCHO, CDEP
 
  integer, public, parameter, dimension(NDRYDEP_CALC) :: &
    DRYDEP_CALC = (/ WES_HNO3, WES_O3,   WES_SO2, &
                     WES_NH3,  WES_NO2 , WES_PAN /)

  !/** Compensation pount approach from CEH used?:

  logical, public, parameter :: COMPENSATION_PT = .false. 


! **sc: characters for printout which is generated directly from Rsurface 
!        module or the program test_dep

!  character(len=6), public, parameter, dimension(NDRYDEP_CALC) :: & 
!      GASNAME = (/ "  HNO3", "    O3", "   SO2", &
!                   "   NH3", "   NO2", "PAN"  /)
  
      

  ! We define also the number of species which will be deposited in
  ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_CALC
  ! The actual species used and their relation to the CDEP_ indices
  ! above will be defined in Init_Vg

  integer, public, parameter ::  NDRYDEP_ADV  = 8

  !/-- we define a type to map indices of species to be deposited
  !   to the lesser number of species where Vg is calculated

   type, public :: depmap
      integer :: adv   ! Index of species in IXADV_ arrays
      integer :: calc ! Index of species in  calculated dep arrays
      real    :: vg   ! if CDEP_SET, give vg in m/s
   end type depmap

   type(depmap), public, dimension(NDRYDEP_ADV):: Dep

   real, public, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost


   logical, private, parameter :: MY_DEBUG = .false.

contains
  subroutine Init_vd
   real :: cms = 0.01     ! Convert to m/s

 ! .... Define the mapping between the advected species and
 !      the specied for which the calculation needs to be done.

   Dep(1) =  depmap( IXADV_HNO3 , CDEP_HNO3, -1.)
   Dep(2) =  depmap( IXADV_PAN,   CDEP_PAN, -1. ) 
   Dep(3) =  depmap( IXADV_NO2,   CDEP_NO2, -1. )
   Dep(4) =  depmap( IXADV_SO2,   CDEP_SO2, -1. )
   Dep(5) =  depmap( IXADV_SO4,   CDEP_SET,  0.1 * cms )
   Dep(6) =  depmap( IXADV_NH3,   CDEP_NH3, -1. )
   Dep(7) =  depmap( IXADV_AMSU,  CDEP_SET,  0.1 * cms  )
   Dep(8) =  depmap( IXADV_AMNI,  CDEP_SET,  0.1 * cms  )

  end subroutine Init_vd

  subroutine Add_ddep(i,j,convfac)
     ! Adds deposition losses to ddep arrays
     integer, intent(in) :: i,j             ! coordinates
     real,    intent(in) ::  convfac   !

     ddep(DDEP_SOX,i,j,IOU_INST) = (  &
          DepLoss(IXADV_SO2) + &
          DepLoss(IXADV_SO4) + &
          DepLoss(IXADV_AMSU)  &
                                    ) * convfac * atwS

     ddep(DDEP_OXN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_HNO3) + &
          DepLoss(IXADV_PAN) + &
          DepLoss(IXADV_NO2) + &
          DepLoss(IXADV_AMNI)  &
                                    ) * convfac * atwN

     ddep(DDEP_RDN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_NH3) + &
    1.5 * DepLoss(IXADV_AMSU) + &
          DepLoss(IXADV_AMNI)  &
                                    ) * convfac * atwN
  end subroutine  Add_ddep

  end module My_UKDep_ml

