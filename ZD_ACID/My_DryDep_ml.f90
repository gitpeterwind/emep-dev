module My_UKDep_ml    ! DryDep_ml
!+
! Module to define the dry deposition components and rates. We
! define the min (vd_min) and max dep. velocity (Vg) first and then derive the
! daytime addition (vd_day).
!/**************************************************************************
!  Specifies which of the possible species (from Wesely's list)
!  are required in the current air pollution model   
!/**************************************************************************

 use DepVariables_ml, only : &
            ECO_WATER,ECO_CONIF_FOREST,ECO_DECID_FOREST, & !ds rv1.6.12
            ECO_CROP,ECO_SEMINAT,ECO_WETLAND               !ds rv1.6.12

 use Derived_ml,    only : f_ddep, ddep,  &   !ds NEW system 16/12/2003
                           f_2d,   d_2d,  &
                           find_one_index, IOU_INST

 use GenSpec_adv_ml               !   e.g. NSPEC_ADV,IXADV_O3,IXADV_H2O2,
 use ModelConstants_ml , only : atwS, atwN &
                              , current_date  !ds rv1_9_17
 use PhysicalConstants_ml, only : AVOG
 !ds mar2005 use Radiation_ml,  only :  zen
 use Wesely_ml
 implicit none
 private

  public :: Init_DepMap
  public :: Add_ddep


  !/** Variables used in deposition calculations
 
  ! DDEP_xx gives the index that will be used in the EMEP model
  ! WES_xx gives the index of the Wesely gas to which this corresponds

  ! ds - new system: 16/12/2003
  ! copied to ACID 16/1/2004

  integer, private, save :: &
    DDEP_SOX,   DDEP_OXN,   DDEP_RDN   &
   ,DDEP_OXNSW, DDEP_OXNCF, DDEP_OXNDF &
   ,DDEP_RDNSW, DDEP_RDNCF, DDEP_RDNDF 


  ! Here we define the minimum set of species which has different
  ! deposition velocities. We calculate Vg for these, and then
  ! can use the rates for other similar species. (e.g. AMSU can use
  ! the Vg for SO4.  Must set NDRYDEP_CALC species

  !/** IMPORTANT: the variables below must match up in the sense that, for 
  ! example, if DDEP_NH3=4 then the 4th element of DRYDEP must be WES_NH3.

  integer, public, parameter :: NDRYDEP_CALC = 5
  integer, public, parameter :: NDRYDEP_AER = 2             !stDep
  integer, public, parameter :: NDRYDEP_TOT = NDRYDEP_CALC + NDRYDEP_AER


  integer, public, parameter :: &
       CDEP_HNO3 = 1, CDEP_NO2 = 2, CDEP_SO2 = 3  &
      ,CDEP_NH3  = 4, CDEP_PAN = 5 &
      ,CDEP_FIN  = 6, CDEP_COA = 7  ! stDep

  integer, public, parameter :: CDEP_SET = -99    



 ! WE NEED A FLUX_CDEP, FLUX_ADV FOR OZONE;
 ! (set to one for non-ozone models)

  logical, public, parameter :: STO_FLUXES = .false.
  integer, public, parameter :: FLUX_CDEP  = 1
  integer, public, parameter :: FLUX_ADV   = 1

 
  integer, public, parameter, dimension(NDRYDEP_CALC) :: &
    DRYDEP_CALC = (/ WES_HNO3, WES_NO2,  WES_SO2, &
                     WES_NH3,  WES_PAN  /)

  !/** Compensation pount approach from CEH used?:

  logical, public, parameter :: COMPENSATION_PT = .false. 



  ! We define also the number of species which will be deposited in
  ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_CALC
  ! The actual species used and their relation to the CDEP_ indices
  ! above will be defined in Init_DepMap

  integer, public, parameter ::  NDRYDEP_ADV  = 13  !SeaS

  !/-- we define a type to map indices of species to be deposited
  !   to the lesser number of species where Vg is calculated

   type, public :: depmap
      integer :: adv   ! Index of species in IXADV_ arrays
      integer :: calc  ! Index of species in  calculated dep arrays
      real    :: vg    ! if CDEP_SET, give vg in m/s
   end type depmap

   type(depmap), public, dimension(NDRYDEP_ADV):: Dep

   real, public, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost


   logical, private, parameter :: MY_DEBUG = .false.

contains
  subroutine Init_DepMap
   real :: cms = 0.01     ! Convert to m/s

 ! .... Define the mapping between the advected species and
 !      the specied for which the calculation needs to be done.

  Dep(1) =  depmap( IXADV_HNO3 , CDEP_HNO3, -1.)
   Dep(2) =  depmap( IXADV_PAN,   CDEP_PAN, -1. )
   Dep(3) =  depmap( IXADV_NO2,   CDEP_NO2, -1. )
   Dep(4) =  depmap( IXADV_SO2,   CDEP_SO2, -1. )
   Dep(5) =  depmap( IXADV_SO4,   CDEP_FIN, -1. )
   Dep(6) =  depmap( IXADV_NH3,   CDEP_NH3, -1. )
   Dep(7) =  depmap( IXADV_aNH4,  CDEP_FIN, -1. )
   Dep(8) =  depmap( IXADV_aNO3,  CDEP_FIN, -1. )
   Dep(9) =  depmap( IXADV_PM25,  CDEP_FIN, -1. )
   Dep(10)=  depmap( IXADV_PMco,  CDEP_COA, -1. )
   Dep(11)=  depmap( IXADV_pNO3,  CDEP_COA, -1.)
   Dep(12)=  depmap( IXADV_SSfi,  CDEP_FIN, -1. )   !SeaS
   Dep(13)=  depmap( IXADV_SSco,  CDEP_COA, -1. )   !SeaS


!####################### ds NEW define indices here #######################

DDEP_SOX   = find_one_index("DDEP_SOX",f_ddep(:)%name)
DDEP_OXN   = find_one_index("DDEP_OXN",f_ddep(:)%name)
DDEP_RDN   = find_one_index("DDEP_RDN",f_ddep(:)%name)

!ds 25/3/2004: General change: waters and wetland categories removed after
!    discussions with CCE/IIASA
!
!ds DDEP_OXNSW = find_one_index("DDEP_OXNSW",f_ddep(:)%name)
DDEP_OXNCF = find_one_index("DDEP_OXNCF",f_ddep(:)%name)
DDEP_OXNDF = find_one_index("DDEP_OXNDF",f_ddep(:)%name)

!ds DDEP_RDNSW = find_one_index("DDEP_RDNSW",f_ddep(:)%name)
DDEP_RDNCF = find_one_index("DDEP_RDNCF",f_ddep(:)%name)
DDEP_RDNDF = find_one_index("DDEP_RDNDF",f_ddep(:)%name)

!####################### ds END of define indices #######################

  end subroutine Init_DepMap

  !<==========================================================================
  subroutine Add_ddep(debug_flag,dt,i,j,convfac,lossfrac,fluxfrac,c_hvegppb)

  !<==========================================================================
     ! Adds deposition losses to ddep arrays
     logical, intent(in) :: debug_flag
     real,    intent(in) :: dt              ! time-step
     integer, intent(in) :: i,j             ! coordinates
     real,    intent(in) ::  convfac, lossfrac
     !ds real,    intent(in) ::  convfac, convfaco3   !
     real, dimension(:,:), intent(in) ::  fluxfrac   ! dim (NADV, NLANDUSE)
     real, dimension(:), intent(in) ::  c_hvegppb   ! dim (NLANDUSE)
     integer :: n, nadv, ihh, idd
     !ds mar2005 integer :: izen                    ! integer of zenith angle
     logical, parameter :: DEBUG_ECO = .false.

     integer, parameter :: N_OXS = 2        ! Number in ox. sulphur family
     real, parameter, dimension(N_OXS) :: OXS = &
             (/ IXADV_SO2, IXADV_SO4 /)
     integer, parameter :: N_OXN = 5        ! Number in ox. nitrogen family
     real, parameter, dimension(N_OXN) :: OXN = &
             (/ IXADV_HNO3, IXADV_PAN, IXADV_NO2, IXADV_aNO3, IXADV_pNO3 /)
     integer, parameter :: N_RDN = 2        ! Number in red. nitrogen family
     real, parameter, dimension(N_RDN) :: RDN = &
             (/ IXADV_NH3, IXADV_aNH4 /)

     real, parameter  :: NMOLE_M3 = 1.0e6*1.0e9/AVOG  ! Converts from 
                                                      ! mol/cm3 to nmole/m3

     real ::  to_nmole, timefrac 
     to_nmole =  NMOLE_M3
     timefrac = dt/3600.0

!! OXIDIZED SULPHUR
!!-----------------------

     ddep(DDEP_SOX,i,j,IOU_INST) = (  &
          DepLoss(IXADV_SO2) + DepLoss(IXADV_SO4) ) * convfac * atwS

!BUG FIX!!!!??? 17/1/2004!!!!

    ! Not used here anyway:;
    ! ddep(DDEP_OXSSW,i,j,IOU_INST) = 0.0
    ! ddep(DDEP_OXSCF,i,j,IOU_INST) = 0.0
    ! ddep(DDEP_OXSDF,i,j,IOU_INST) = 0.0
    ! ddep(DDEP_OXSCR,i,j,IOU_INST) = 0.0
    ! ddep(DDEP_OXSSN,i,j,IOU_INST) = 0.0
    ! ddep(DDEP_OXSWE,i,j,IOU_INST) = 0.0

     do n = 1,N_OXS
       nadv = OXS(n)

      ! == ds make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !

      !ddep(DDEP_OXSSW,i,j,IOU_INST) = ddep(DDEP_OXSSW,i,j,IOU_INST) +  &
      !     sum( fluxfrac(nadv,ECO_WATER) ) * DepLoss(nadv)

      !ddep(DDEP_OXSCF,i,j,IOU_INST) = ddep(DDEP_OXSCF,i,j,IOU_INST) +  &
      !     sum( fluxfrac(nadv,ECO_CONIF_FOREST) ) * DepLoss(nadv)


      !ddep(DDEP_OXSDF,i,j,IOU_INST) = ddep(DDEP_OXSDF,i,j,IOU_INST) +  &
      !     sum( fluxfrac(nadv,ECO_DECID_FOREST) ) * DepLoss(nadv)

      !ddep(DDEP_OXSCR,i,j,IOU_INST) = ddep(DDEP_OXSCR,i,j,IOU_INST) +  &
      !     sum( fluxfrac(nadv,ECO_CROP) ) * DepLoss(nadv)


      !ddep(DDEP_OXSSN,i,j,IOU_INST) = ddep(DDEP_OXSSN,i,j,IOU_INST) +  &
      !     sum( fluxfrac(nadv,ECO_SEMINAT) ) * DepLoss(nadv)

      !ddep(DDEP_OXSWE,i,j,IOU_INST) = ddep(DDEP_OXSWE,i,j,IOU_INST) +  &
      !     sum( fluxfrac(nadv,ECO_WETLAND) ) * DepLoss(nadv)
      ! ==                                                            ==== !

     end do

    !ddep(DDEP_OXSSW,i,j,IOU_INST) = ddep(DDEP_OXSSW,i,j,IOU_INST)*convfac*atwS
    !ddep(DDEP_OXSCF,i,j,IOU_INST) = ddep(DDEP_OXSCF,i,j,IOU_INST)*convfac*atwS
    !ddep(DDEP_OXSDF,i,j,IOU_INST) = ddep(DDEP_OXSDF,i,j,IOU_INST)*convfac*atwS
    !ddep(DDEP_OXSCR,i,j,IOU_INST) = ddep(DDEP_OXSCR,i,j,IOU_INST)*convfac*atwS
    !ddep(DDEP_OXSSN,i,j,IOU_INST) = ddep(DDEP_OXSSN,i,j,IOU_INST)*convfac*atwS
    !ddep(DDEP_OXSWE,i,j,IOU_INST) = ddep(DDEP_OXSWE,i,j,IOU_INST)*convfac*atwS



!! OXIDIZED NITROGEN
!!-----------------------

     ddep(DDEP_OXN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_HNO3) +  DepLoss(IXADV_PAN) +  DepLoss(IXADV_NO2) + &
          DepLoss(IXADV_aNO3)+  DepLoss(IXADV_pNO3)  ) * convfac * atwN

!BUG FIX!!!!??? 16/1/2004!!!!
!ds 25/3/2004: General change: waters and wetland categories removed after
!    discussions with CCE/IIASA

     !ds ddep(DDEP_OXNSW,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXNCF,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXNDF,i,j,IOU_INST) = 0.0

     do n = 1, N_OXN
       nadv = OXN(n)
         
      ! == ds make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !

       !ds ddep(DDEP_OXNSW,i,j,IOU_INST) = ddep(DDEP_OXNSW,i,j,IOU_INST) +  &
       !ds      sum( fluxfrac(nadv,ECO_WATER) ) * DepLoss(nadv)

       ddep(DDEP_OXNCF,i,j,IOU_INST) = ddep(DDEP_OXNCF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_CONIF_FOREST) ) * DepLoss(nadv)


       ddep(DDEP_OXNDF,i,j,IOU_INST) = ddep(DDEP_OXNDF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_DECID_FOREST) ) * DepLoss(nadv)

      ! ==                                                            ==== !


     end do

     !ds ddep(DDEP_OXNSW,i,j,IOU_INST) = ddep(DDEP_OXNSW,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNCF,i,j,IOU_INST) = ddep(DDEP_OXNCF,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNDF,i,j,IOU_INST) = ddep(DDEP_OXNDF,i,j,IOU_INST)*convfac*atwN



!! REDUCED NITROGEN
!!-----------------------

     ddep(DDEP_RDN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_NH3) +  DepLoss(IXADV_aNH4)  ) * convfac * atwN

!BUG FIX!!!!??? 16/1/2004!!!!

     !ds ddep(DDEP_RDNSW,i,j,IOU_INST) = 0.0
     ddep(DDEP_RDNCF,i,j,IOU_INST) = 0.0
     ddep(DDEP_RDNDF,i,j,IOU_INST) = 0.0

     do n = 1, N_RDN
       nadv = RDN(n)
         
      ! == ds make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !

      !ds  ddep(DDEP_RDNSW,i,j,IOU_INST) = ddep(DDEP_RDNSW,i,j,IOU_INST) +  &
      !ds       sum( fluxfrac(nadv,ECO_WATER) ) * DepLoss(nadv)

       ddep(DDEP_RDNCF,i,j,IOU_INST) = ddep(DDEP_RDNCF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_CONIF_FOREST) ) * DepLoss(nadv)


       ddep(DDEP_RDNDF,i,j,IOU_INST) = ddep(DDEP_RDNDF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_DECID_FOREST) ) * DepLoss(nadv)

      ! ==                                                            ==== !

     end do

    !ds  ddep(DDEP_RDNSW,i,j,IOU_INST) = ddep(DDEP_RDNSW,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNCF,i,j,IOU_INST) = ddep(DDEP_RDNCF,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNDF,i,j,IOU_INST) = ddep(DDEP_RDNDF,i,j,IOU_INST)*convfac*atwN
    !ddep(DDEP_RDNCR,i,j,IOU_INST) = ddep(DDEP_RDNCR,i,j,IOU_INST)*convfac*atwN
    !ddep(DDEP_RDNSN,i,j,IOU_INST) = ddep(DDEP_RDNSN,i,j,IOU_INST)*convfac*atwN
    !ddep(DDEP_RDNWE,i,j,IOU_INST) = ddep(DDEP_RDNWE,i,j,IOU_INST)*convfac*atwN

   !---- end ecosystem specific ----------------------------------------------

  end subroutine  Add_ddep

  end module My_UKDep_ml

