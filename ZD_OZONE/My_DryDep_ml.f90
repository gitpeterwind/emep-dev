module My_UKDep_ml    ! DryDep_ml
!+
! Module to define the dry deposition components and rates. We
! define the min (vd_min) and max dep. velocity (Vg) first and then derive the
! daytime addition (vd_day).
!/**************************************************************************
!  Specifies which of the possible species (from Wesely's list)
!  are required in the current air pollution model   
!/**************************************************************************

 use DepVariables_ml, only : unit_flux, lai_flux, leaf_flux, & 
            WHEAT, BEECH, & !ds rv1_9_15
            ECO_WATER,ECO_CONIF_FOREST,ECO_DECID_FOREST, & !ds rv1.6.12
            ECO_CROP,ECO_SEMINAT,ECO_WETLAND               !ds rv1.6.12

 use Derived_ml,    only : f_ddep, ddep,  &   !ds NEW system 16/12/2003
                           f_2d,   d_2d,  &
                           find_one_index, IOU_INST

 use GenSpec_adv_ml               !   e.g. NSPEC_ADV,IXADV_O3,IXADV_H2O2,
 use ModelConstants_ml , only : atwS, atwN &
                              , current_date, AOT_HORIZON  !ds rv1_9_17
 use PhysicalConstants_ml, only : AVOG
 use Radiation_ml,  only :  zen               !ds rv1_9_17
 use Wesely_ml
 implicit none
 private

  public :: Init_DepMap
  public :: Add_ddep


  !/** Variables used in deposition calculations
 
  ! DDEP_xx gives the index that will be used in the EMEP model
  ! WES_xx gives the index of the Wesely gas to which this corresponds

  ! ds - new system: 16/12/2003

  integer, private, save :: &
    DDEP_SOX,   DDEP_OXN,   DDEP_RDN,  &
    DDEP_OXSSW, DDEP_OXSCF, DDEP_OXSDF, DDEP_OXSCR, DDEP_OXSSN, DDEP_OXSWE, &
    DDEP_OXNSW, DDEP_OXNCF, DDEP_OXNDF, DDEP_OXNCR, DDEP_OXNSN, DDEP_OXNWE, &
    DDEP_RDNSW, DDEP_RDNCF, DDEP_RDNDF, DDEP_RDNCR, DDEP_RDNSN, DDEP_RDNWE, &
    D2_FSTDF00, D2_FSTDF08, D2_FSTDF16, D2_FSTWH00, D2_FSTWH20, D2_FSTWH40, &
    D2_FSTWH60,&
    D2_O3DF,    D2_O3WH, &
    D2_EUAOT30WH, D2_EUAOT40WH, & !ds rv1_9_17
    D2_EUAOT30DF, D2_EUAOT40DF, &
    D2_UNAOT30WH, D2_UNAOT40WH, & !ds rv1_9_17
    D2_UNAOT30DF, D2_UNAOT40DF


  ! Here we define the minimum set of species which has different
  ! deposition velocities. We calculate Vg for these, and then
  ! can use the rates for other similar species. (e.g. AMSU can use
  ! the Vg for SO4.  Must set NDRYDEP_CALC species

  !/** IMPORTANT: the variables below must match up in the sense that, for 
  ! example, if DDEP_NH3=4 then the 4th element of DRYDEP must be WES_NH3.

  integer, public, parameter :: NDRYDEP_CALC = 10
  integer, public, parameter :: NDRYDEP_AER = 2             !stDep
  integer, public, parameter :: NDRYDEP_TOT = NDRYDEP_CALC + NDRYDEP_AER


  integer, public, parameter :: &
       CDEP_HNO3 = 1, CDEP_O3  = 2, CDEP_SO2 = 3  &
      ,CDEP_NH3  = 4, CDEP_NO2 = 5, CDEP_PAN  = 6 &
      ,CDEP_H2O2 = 7, CDEP_ALD = 8, CDEP_HCHO = 9, &
       CDEP_OP = 10, &
       CDEP_FIN = 11, CDEP_COA = 12 !stDep

  integer, public, parameter :: CDEP_SET = -99    



 ! WE NEED A FLUX_CDEP, FLUX_ADV FOR OZONE;
 ! (set to one for non-ozone models)

  logical, public, parameter :: STO_FLUXES = .true.
  integer, public, parameter :: FLUX_CDEP  = CDEP_O3
  integer, public, parameter :: FLUX_ADV   = IXADV_O3

 
  integer, public, parameter, dimension(NDRYDEP_CALC) :: &
    DRYDEP_CALC = (/ WES_HNO3, WES_O3,   WES_SO2, &
                     WES_NH3,  WES_NO2 , WES_PAN, &
                     WES_H2O2, WES_ALD, WES_HCHO, WES_OP    /)

  !/** Compensation pount approach from CEH used?:

  logical, public, parameter :: COMPENSATION_PT = .false. 



  ! We define also the number of species which will be deposited in
  ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_CALC
  ! The actual species used and their relation to the CDEP_ indices
  ! above will be defined in Init_DepMap

  integer, public, parameter ::  NDRYDEP_ADV  = 19

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
   Dep(5) =  depmap( IXADV_SO4,   CDEP_FIN,  -1) !ds 0.1 * cms )
   Dep(6) =  depmap( IXADV_NH3,   CDEP_NH3, -1. )
   Dep(7) =  depmap( IXADV_aNH4,  CDEP_FIN,  -1) !ds 0.1 * cms  )
   Dep(8) =  depmap( IXADV_aNO3,  CDEP_FIN,  -1) !ds 0.1 * cms  )
   Dep(9) =  depmap( IXADV_O3   , CDEP_O3  , -1.)
   Dep(10) =  depmap( IXADV_H2O2 , CDEP_H2O2, -1.)
   Dep(11) =  depmap( IXADV_MPAN , CDEP_PAN , -1.)
   Dep(12) =  depmap( IXADV_HCHO , CDEP_HCHO, -1.)
   Dep(13) =  depmap( IXADV_CH3CHO,CDEP_ALD , -1.)
   Dep(14) =  depmap( IXADV_MAL   ,CDEP_ALD , -1.)
   Dep(15) =  depmap( IXADV_CH3O2H,CDEP_OP  , -1.)
   Dep(16) =  depmap( IXADV_C2H5OOH,CDEP_OP  , -1.)
   Dep(17) =  depmap( IXADV_pNO3,  CDEP_COA, -1.)
   Dep(18) =  depmap( IXADV_PM25,  CDEP_FIN, -1. )
   Dep(19) =  depmap( IXADV_PMco,  CDEP_COA, -1. )


!####################### ds NEW define indices here #######################

DDEP_SOX   = find_one_index("DDEP_SOX",f_ddep(:)%name)
DDEP_OXN   = find_one_index("DDEP_OXN",f_ddep(:)%name)
DDEP_RDN   = find_one_index("DDEP_RDN",f_ddep(:)%name)

DDEP_OXSSW = find_one_index("DDEP_OXSSW",f_ddep(:)%name)
DDEP_OXSCF = find_one_index("DDEP_OXSCF",f_ddep(:)%name)
DDEP_OXSDF = find_one_index("DDEP_OXSDF",f_ddep(:)%name)
DDEP_OXSCR = find_one_index("DDEP_OXSCR",f_ddep(:)%name)
DDEP_OXSSN = find_one_index("DDEP_OXSSN",f_ddep(:)%name)
DDEP_OXSWE = find_one_index("DDEP_OXSWE",f_ddep(:)%name)

DDEP_OXNSW = find_one_index("DDEP_OXNSW",f_ddep(:)%name)
DDEP_OXNCF = find_one_index("DDEP_OXNCF",f_ddep(:)%name)
DDEP_OXNDF = find_one_index("DDEP_OXNDF",f_ddep(:)%name)
DDEP_OXNCR = find_one_index("DDEP_OXNCR",f_ddep(:)%name)
DDEP_OXNSN = find_one_index("DDEP_OXNSN",f_ddep(:)%name)
DDEP_OXNWE = find_one_index("DDEP_OXNWE",f_ddep(:)%name)

DDEP_RDNSW = find_one_index("DDEP_RDNSW",f_ddep(:)%name)
DDEP_RDNCF = find_one_index("DDEP_RDNCF",f_ddep(:)%name)
DDEP_RDNDF = find_one_index("DDEP_RDNDF",f_ddep(:)%name)
DDEP_RDNCR = find_one_index("DDEP_RDNCR",f_ddep(:)%name)
DDEP_RDNSN = find_one_index("DDEP_RDNSN",f_ddep(:)%name)
DDEP_RDNWE = find_one_index("DDEP_RDNWE",f_ddep(:)%name)

!ds rv1_9_15 changes:
D2_FSTDF00 = find_one_index("D2_FSTDF00",f_2d(:)%name)
D2_FSTDF08 = find_one_index("D2_FSTDF08",f_2d(:)%name)
D2_FSTDF16 = find_one_index("D2_FSTDF16",f_2d(:)%name)

D2_FSTWH00 = find_one_index("D2_FSTWH00",f_2d(:)%name)
D2_FSTWH20 = find_one_index("D2_FSTWH20",f_2d(:)%name)
D2_FSTWH40 = find_one_index("D2_FSTWH40",f_2d(:)%name)
D2_FSTWH60 = find_one_index("D2_FSTWH60",f_2d(:)%name)

D2_O3DF    = find_one_index("D2_O3DF   ",f_2d(:)%name)
D2_O3WH    = find_one_index("D2_O3WH   ",f_2d(:)%name)

D2_EUAOT30WH    = find_one_index("D2_EUAOT30WH",f_2d(:)%name)
D2_EUAOT40WH    = find_one_index("D2_EUAOT40WH",f_2d(:)%name)
D2_EUAOT30DF    = find_one_index("D2_EUAOT30DF",f_2d(:)%name)
D2_EUAOT40DF    = find_one_index("D2_EUAOT40DF",f_2d(:)%name)

D2_UNAOT30WH    = find_one_index("D2_UNAOT30WH",f_2d(:)%name)
D2_UNAOT40WH    = find_one_index("D2_UNAOT40WH",f_2d(:)%name)
D2_UNAOT30DF    = find_one_index("D2_UNAOT30DF",f_2d(:)%name)
D2_UNAOT40DF    = find_one_index("D2_UNAOT40DF",f_2d(:)%name)

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
     real :: o3WH, o3DF   ! O3 over wheat, decid forest
     integer :: izen                    ! integer of zenith angle
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

   !BUG FIX!!!!??? 16/1/2004!!!!

     ddep(DDEP_OXSSW,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXSCF,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXSDF,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXSCR,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXSSN,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXSWE,i,j,IOU_INST) = 0.0

     do n = 1,N_OXS
       nadv = OXS(n)

      ! == ds make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !

       ddep(DDEP_OXSSW,i,j,IOU_INST) = ddep(DDEP_OXSSW,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_WATER) ) * DepLoss(nadv)

       ddep(DDEP_OXSCF,i,j,IOU_INST) = ddep(DDEP_OXSCF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_CONIF_FOREST) ) * DepLoss(nadv)


       ddep(DDEP_OXSDF,i,j,IOU_INST) = ddep(DDEP_OXSDF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_DECID_FOREST) ) * DepLoss(nadv)

       ddep(DDEP_OXSCR,i,j,IOU_INST) = ddep(DDEP_OXSCR,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_CROP) ) * DepLoss(nadv)


       ddep(DDEP_OXSSN,i,j,IOU_INST) = ddep(DDEP_OXSSN,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_SEMINAT) ) * DepLoss(nadv)

       ddep(DDEP_OXSWE,i,j,IOU_INST) = ddep(DDEP_OXSWE,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_WETLAND) ) * DepLoss(nadv)
      ! ==                                                            ==== !

     end do

     ddep(DDEP_OXSSW,i,j,IOU_INST) = ddep(DDEP_OXSSW,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSCF,i,j,IOU_INST) = ddep(DDEP_OXSCF,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSDF,i,j,IOU_INST) = ddep(DDEP_OXSDF,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSCR,i,j,IOU_INST) = ddep(DDEP_OXSCR,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSSN,i,j,IOU_INST) = ddep(DDEP_OXSSN,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSWE,i,j,IOU_INST) = ddep(DDEP_OXSWE,i,j,IOU_INST)*convfac*atwS



!! OXIDIZED NITROGEN
!!-----------------------


     ddep(DDEP_OXN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_HNO3) +  DepLoss(IXADV_PAN) +  DepLoss(IXADV_NO2) + &
          DepLoss(IXADV_aNO3)+  DepLoss(IXADV_pNO3)  ) * convfac * atwN

   !BUG FIX!!!!??? 17/1/2004!!!!

     ddep(DDEP_OXNSW,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXNCF,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXNDF,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXNCR,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXNSN,i,j,IOU_INST) = 0.0
     ddep(DDEP_OXNWE,i,j,IOU_INST) = 0.0

     do n = 1, N_OXN
       nadv = OXN(n)
         
      ! == ds make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !

       ddep(DDEP_OXNSW,i,j,IOU_INST) = ddep(DDEP_OXNSW,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_WATER) ) * DepLoss(nadv)

       ddep(DDEP_OXNCF,i,j,IOU_INST) = ddep(DDEP_OXNCF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_CONIF_FOREST) ) * DepLoss(nadv)


       ddep(DDEP_OXNDF,i,j,IOU_INST) = ddep(DDEP_OXNDF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_DECID_FOREST) ) * DepLoss(nadv)

       ddep(DDEP_OXNCR,i,j,IOU_INST) = ddep(DDEP_OXNCR,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_CROP) ) * DepLoss(nadv)


       ddep(DDEP_OXNSN,i,j,IOU_INST) = ddep(DDEP_OXNSN,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_SEMINAT) ) * DepLoss(nadv)

       ddep(DDEP_OXNWE,i,j,IOU_INST) = ddep(DDEP_OXNWE,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_WETLAND) ) * DepLoss(nadv)
      ! ==                                                            ==== !


     end do

     ddep(DDEP_OXNSW,i,j,IOU_INST) = ddep(DDEP_OXNSW,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNCF,i,j,IOU_INST) = ddep(DDEP_OXNCF,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNDF,i,j,IOU_INST) = ddep(DDEP_OXNDF,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNCR,i,j,IOU_INST) = ddep(DDEP_OXNCR,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNSN,i,j,IOU_INST) = ddep(DDEP_OXNSN,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNWE,i,j,IOU_INST) = ddep(DDEP_OXNWE,i,j,IOU_INST)*convfac*atwN



!! REDUCED NITROGEN
!!-----------------------

     ddep(DDEP_RDN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_NH3) +  DepLoss(IXADV_aNH4)  ) * convfac * atwN


   !BUG FIX!!!!??? 17/1/2004!!!!

     ddep(DDEP_RDNSW,i,j,IOU_INST) = 0.0
     ddep(DDEP_RDNCF,i,j,IOU_INST) = 0.0
     ddep(DDEP_RDNDF,i,j,IOU_INST) = 0.0
     ddep(DDEP_RDNCR,i,j,IOU_INST) = 0.0
     ddep(DDEP_RDNSN,i,j,IOU_INST) = 0.0
     ddep(DDEP_RDNWE,i,j,IOU_INST) = 0.0

     do n = 1, N_RDN
       nadv = RDN(n)
         
      ! == ds make use of ECO_ arrays from DepVariables - specifies   ==== !
      !    which landuse is in which category                         ==== !

       ddep(DDEP_RDNSW,i,j,IOU_INST) = ddep(DDEP_RDNSW,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_WATER) ) * DepLoss(nadv)

       ddep(DDEP_RDNCF,i,j,IOU_INST) = ddep(DDEP_RDNCF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_CONIF_FOREST) ) * DepLoss(nadv)


       ddep(DDEP_RDNDF,i,j,IOU_INST) = ddep(DDEP_RDNDF,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_DECID_FOREST) ) * DepLoss(nadv)

       ddep(DDEP_RDNCR,i,j,IOU_INST) = ddep(DDEP_RDNCR,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_CROP) ) * DepLoss(nadv)


       ddep(DDEP_RDNSN,i,j,IOU_INST) = ddep(DDEP_RDNSN,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_SEMINAT) ) * DepLoss(nadv)

       ddep(DDEP_RDNWE,i,j,IOU_INST) = ddep(DDEP_RDNWE,i,j,IOU_INST) +  &
            sum( fluxfrac(nadv,ECO_WETLAND) ) * DepLoss(nadv)
      ! ==                                                            ==== !

     end do

     ddep(DDEP_RDNSW,i,j,IOU_INST) = ddep(DDEP_RDNSW,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNCF,i,j,IOU_INST) = ddep(DDEP_RDNCF,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNDF,i,j,IOU_INST) = ddep(DDEP_RDNDF,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNCR,i,j,IOU_INST) = ddep(DDEP_RDNCR,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNSN,i,j,IOU_INST) = ddep(DDEP_RDNSN,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNWE,i,j,IOU_INST) = ddep(DDEP_RDNWE,i,j,IOU_INST)*convfac*atwN

   !ds  d_2d(D2_FSTCF0,i,j,IOU_INST) =  leaf_flux(1) !* to_nmole
   !ds  d_2d(D2_FSTDF0,i,j,IOU_INST) =  leaf_flux(2) !* to_nmole
   !ds  d_2d(D2_FSTTC0,i,j,IOU_INST) =  leaf_flux(5) !* to_nmole
   !ds  d_2d(D2_FSTMC0,i,j,IOU_INST) =  leaf_flux(6) !* to_nmole
   !ds  d_2d(D2_FSTGR0,i,j,IOU_INST) =  leaf_flux(10)! * to_nmole
   !ds  d_2d(D2_FSTWH0,i,j,IOU_INST) =  leaf_flux(9)! * to_nmole


!MAPPING_MANUAL CHANGES: !ds rv1_9_15

!Beech:
     d_2d(D2_FSTDF00,i,j,IOU_INST) =  leaf_flux(BEECH)
     d_2d(D2_FSTDF08,i,j,IOU_INST) =  max(leaf_flux(BEECH)-0.8,0.0)
     d_2d(D2_FSTDF16,i,j,IOU_INST) =  max(leaf_flux(BEECH)-1.6,0.0)
!Wheat
     d_2d(D2_FSTWH00,i,j,IOU_INST) =  leaf_flux(WHEAT)
     d_2d(D2_FSTWH20,i,j,IOU_INST) =  max(leaf_flux(WHEAT)-2.0,0.0)
     d_2d(D2_FSTWH40,i,j,IOU_INST) =  max(leaf_flux(WHEAT)-4.0,0.0)
     d_2d(D2_FSTWH60,i,j,IOU_INST) =  max(leaf_flux(WHEAT)-6.0,0.0)

   !--- ecosystem specific concentrations..
   !ds - use Conif forest for forests - safer for growing seasons

     idd      =    current_date%day              ! for debugging
     ihh      =    current_date%hour             ! for debugging


     o3WH = c_hvegppb(WHEAT)* lossfrac
     o3DF = c_hvegppb(BEECH)* lossfrac

     d_2d(D2_O3DF,i,j,IOU_INST) =   o3DF
     d_2d(D2_O3WH,i,j,IOU_INST) =   o3WH

     if ( ihh >= 9 .and. ihh <= 21 ) then ! 8-20 CET, assuming summertime

        d_2d(D2_EUAOT30WH,i,j,IOU_INST) =  max(o3WH-30.0,0.0) * timefrac
        d_2d(D2_EUAOT40WH,i,j,IOU_INST) =  max(o3WH-40.0,0.0) * timefrac
        d_2d(D2_EUAOT30DF,i,j,IOU_INST) =  max(o3DF-30.0,0.0) * timefrac
        d_2d(D2_EUAOT40DF,i,j,IOU_INST) =  max(o3DF-40.0,0.0) * timefrac
     else ! rv1_9_19 bug fix!!!!
        d_2d(D2_EUAOT30WH,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT40WH,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT30DF,i,j,IOU_INST) =  0.0
        d_2d(D2_EUAOT40DF,i,j,IOU_INST) =  0.0
     end if


    !/-- Calcuates AOT values for specific veg. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )

           izen = max(1,int( zen(i,j) + 0.5))

           if ( izen < AOT_HORIZON ) then

             d_2d(D2_UNAOT30WH,i,j,IOU_INST) =  max(o3WH-30.0,0.0) * timefrac
             d_2d(D2_UNAOT40WH,i,j,IOU_INST) =  max(o3WH-40.0,0.0) * timefrac
             d_2d(D2_UNAOT30DF,i,j,IOU_INST) =  max(o3DF-30.0,0.0) * timefrac
             d_2d(D2_UNAOT40DF,i,j,IOU_INST) =  max(o3DF-40.0,0.0) * timefrac

             !Derivd had: d_2d(n, i,j,IOU_INST ) = o3 * timefrac  
           else ! rv1_9_19 bug fix!!!!
             d_2d(D2_UNAOT30WH,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT40WH,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT30DF,i,j,IOU_INST) =  0.0
             d_2d(D2_UNAOT40DF,i,j,IOU_INST) =  0.0
           end if

    if ( DEBUG_ECO .and. debug_flag ) then
          write(6,"(a12,i5,f7.2,2es12.3)") "DEBUG_ECO ", ihh, o3DF, &
             leaf_flux(BEECH), d_2d(D2_FSTDF00,i,j,IOU_INST)
    end if ! DEBUG

   !---- end ecosystem specific ----------------------------------------------

  end subroutine  Add_ddep

  end module My_UKDep_ml

