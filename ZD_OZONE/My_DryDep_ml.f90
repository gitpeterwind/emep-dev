module My_UKDep_ml    ! DryDep_ml
!+
! Module to define the dry deposition components and rates. We
! define the min (vd_min) and max dep. velocity (Vg) first and then derive the
! daytime addition (vd_day).
!/**************************************************************************
!  Specifies which of the possible species (from Wesely's list)
!  are required in the current air pollution model   
!/**************************************************************************

 use DepVariables_ml, only : unit_flux, lai_flux
 use My_Derived_ml , only : DDEP_SOX,DDEP_OXN,DDEP_RDN, &
                          !!  ECOSYSTEM Specific
                            DDEP_OXSW, DDEP_OXSF, DDEP_OXSC, DDEP_OXSM, & 
                            DDEP_OXNW, DDEP_OXNF, DDEP_OXNC, DDEP_OXNM, & 
                            DDEP_RDNW, DDEP_RDNF, DDEP_RDNC, DDEP_RDNM, & 
                            DDEP_O3W,  DDEP_O3F,  DDEP_O3C,  DDEP_O3M,  & 
                            DDEP_OXSR, DDEP_OXNR, DDEP_RDNR,            & 
                            DDEP_STOCFU, DDEP_STODFU, DDEP_STOTCU,      &
                            DDEP_STOMCU, DDEP_STOGRU,                   &
                            DDEP_STOCFL, DDEP_STODFL, DDEP_STOTCL,      &
                            DDEP_STOMCL, DDEP_STOGRL, &
                             IOU_INST    &!updates inst. dep. fields
                           , ddep         ! 2d fields
 use GenSpec_adv_ml !, only: NSPEC_ADV &
                   !,IXADV_O3,IXADV_H2O2,
                  ! ,IXADV_SO2,IXADV_NO2  &
                  ! ,IXADV_HNO3 &  !! ,IXADV_EC  ,IXADV_OC & 
   ! ,  IXADV_PAN ,IXADV_SO4,IXADV_NH3,IXADV_aNO3, IXADV_aNH4
 use ModelConstants_ml , only : atwS, atwN
 use PhysicalConstants_ml, only : AVOG
 use Wesely_ml
 implicit none
 private

  public :: Init_DepMap
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

  integer, public, parameter :: NDRYDEP_CALC = 10
!stDep
  integer, public, parameter :: NDRYDEP_AER = 2
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


! **sc: characters for printout which is generated directly from Rsurface 
!        module or the program test_dep

!  character(len=6), public, parameter, dimension(NDRYDEP_CALC) :: & 
!      GASNAME = (/ "  HNO3", "    O3", "   SO2", &
!                   "   NH3", "   NO2", "PAN"  /)
  
      

  ! We define also the number of species which will be deposited in
  ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_CALC
  ! The actual species used and their relation to the CDEP_ indices
  ! above will be defined in Init_DepMap

  integer, public, parameter ::  NDRYDEP_ADV  = 17

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
   Dep(17) =  depmap( IXADV_pNO3,CDEP_COA, -1.)

  end subroutine Init_DepMap

  subroutine Add_ddep(i,j,convfac,convfaco3,fluxfrac)
     ! Adds deposition losses to ddep arrays
     integer, intent(in) :: i,j             ! coordinates
     real,    intent(in) ::  convfac, convfaco3   !
     real, dimension(:,:), intent(in) ::  fluxfrac   ! dim (NADV, NLANDUSE)
     integer :: n, nadv

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

     real to_nmole 
     to_nmole =  NMOLE_M3

!! OXIDIZED SULPHUR
!!-----------------------

     ddep(DDEP_SOX,i,j,IOU_INST) = (  &
          DepLoss(IXADV_SO2) + &
          DepLoss(IXADV_SO4)&!hf + &
!hf          DepLoss(IXADV_AMSU)  &
                                    ) * convfac * atwS

     do n = 1,N_OXS
       nadv = OXS(n)
         
       ddep(DDEP_OXSW,i,j,IOU_INST) = ddep(DDEP_OXSW,i,j,IOU_INST) +  &
            fluxfrac(nadv,15) * DepLoss(nadv)  !CRUDE, 15=water for now

       ddep(DDEP_OXSF,i,j,IOU_INST) = ddep(DDEP_OXSF,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,1) + fluxfrac(nadv,2) + &
              fluxfrac(nadv,3) + fluxfrac(nadv,4)   ) * DepLoss(nadv) 

       ddep(DDEP_OXSC,i,j,IOU_INST) = ddep(DDEP_OXSC,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,5) + fluxfrac(nadv,6) ) * DepLoss(nadv)

       ddep(DDEP_OXSM,i,j,IOU_INST) = ddep(DDEP_OXSM,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,8) + fluxfrac(nadv,11) ) * DepLoss(nadv)

       ddep(DDEP_OXSR,i,j,IOU_INST) = ddep(DDEP_OXSR,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,7) + fluxfrac(nadv,9) + &
              fluxfrac(nadv,10) + fluxfrac(nadv,12) + &
              fluxfrac(nadv,13) + fluxfrac(nadv,14) + &
              fluxfrac(nadv,16) + fluxfrac(nadv,17)   ) * DepLoss(nadv) 

     end do
     ddep(DDEP_OXSW,i,j,IOU_INST) = ddep(DDEP_OXSW,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSF,i,j,IOU_INST) = ddep(DDEP_OXSF,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSC,i,j,IOU_INST) = ddep(DDEP_OXSC,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSM,i,j,IOU_INST) = ddep(DDEP_OXSM,i,j,IOU_INST)*convfac*atwS
     ddep(DDEP_OXSR,i,j,IOU_INST) = ddep(DDEP_OXSR,i,j,IOU_INST)*convfac*atwS



!! OXIDIZED NITROGEN
!!-----------------------

     ddep(DDEP_OXN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_HNO3) + &
          DepLoss(IXADV_PAN) + &
          DepLoss(IXADV_NO2) + &
          DepLoss(IXADV_aNO3)+ &
          DepLoss(IXADV_pNO3)  &
                                    ) * convfac * atwN

     do n = 1, N_OXN
       nadv = OXN(n)
         
       ddep(DDEP_OXNW,i,j,IOU_INST) = ddep(DDEP_OXNW,i,j,IOU_INST) +  &
            fluxfrac(nadv,15) * DepLoss(nadv)  !CRUDE, 15=water for now

       ddep(DDEP_OXNF,i,j,IOU_INST) = ddep(DDEP_OXNF,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,1) + fluxfrac(nadv,2) + &
              fluxfrac(nadv,3) + fluxfrac(nadv,4)   ) * DepLoss(nadv) 

       ddep(DDEP_OXNC,i,j,IOU_INST) = ddep(DDEP_OXNC,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,5) + fluxfrac(nadv,6) ) * DepLoss(nadv)

       ddep(DDEP_OXNM,i,j,IOU_INST) = ddep(DDEP_OXNM,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,8) + fluxfrac(nadv,11) ) * DepLoss(nadv)

       ddep(DDEP_OXNR,i,j,IOU_INST) = ddep(DDEP_OXNR,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,7) + fluxfrac(nadv,9) + &
              fluxfrac(nadv,10) + fluxfrac(nadv,12) + &
              fluxfrac(nadv,13) + fluxfrac(nadv,14) + &
              fluxfrac(nadv,16) + fluxfrac(nadv,17)   ) * DepLoss(nadv) 


     end do
     ddep(DDEP_OXNW,i,j,IOU_INST) = ddep(DDEP_OXNW,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNF,i,j,IOU_INST) = ddep(DDEP_OXNF,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNC,i,j,IOU_INST) = ddep(DDEP_OXNC,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNM,i,j,IOU_INST) = ddep(DDEP_OXNM,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_OXNR,i,j,IOU_INST) = ddep(DDEP_OXNR,i,j,IOU_INST)*convfac*atwN



!! REDUCED NITROGEN
!!-----------------------

     ddep(DDEP_RDN,i,j,IOU_INST) = ( &
          DepLoss(IXADV_NH3) + &
!hf    1.5 * DepLoss(IXADV_AMSU) + &
          DepLoss(IXADV_aNH4)  &
                                    ) * convfac * atwN

     do n = 1, N_RDN
       nadv = RDN(n)
         
       ddep(DDEP_RDNW,i,j,IOU_INST) = ddep(DDEP_RDNW,i,j,IOU_INST) +  &
            fluxfrac(nadv,15) * DepLoss(nadv)  !CRUDE, 15=water for now

       ddep(DDEP_RDNF,i,j,IOU_INST) = ddep(DDEP_RDNF,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,1) + fluxfrac(nadv,2) + &
              fluxfrac(nadv,3) + fluxfrac(nadv,4)   ) * DepLoss(nadv) 

       ddep(DDEP_RDNC,i,j,IOU_INST) = ddep(DDEP_RDNC,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,5) + fluxfrac(nadv,6) ) * DepLoss(nadv)

       ddep(DDEP_RDNM,i,j,IOU_INST) = ddep(DDEP_RDNM,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,8) + fluxfrac(nadv,11) ) * DepLoss(nadv)

       ddep(DDEP_RDNR,i,j,IOU_INST) = ddep(DDEP_RDNR,i,j,IOU_INST) +  &
            ( fluxfrac(nadv,7) + fluxfrac(nadv,9) + &
              fluxfrac(nadv,10) + fluxfrac(nadv,12) + &
              fluxfrac(nadv,13) + fluxfrac(nadv,14) + &
              fluxfrac(nadv,16) + fluxfrac(nadv,17)   ) * DepLoss(nadv) 

     end do
     ddep(DDEP_RDNW,i,j,IOU_INST) = ddep(DDEP_RDNW,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNF,i,j,IOU_INST) = ddep(DDEP_RDNF,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNC,i,j,IOU_INST) = ddep(DDEP_RDNC,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNM,i,j,IOU_INST) = ddep(DDEP_RDNM,i,j,IOU_INST)*convfac*atwN
     ddep(DDEP_RDNR,i,j,IOU_INST) = ddep(DDEP_RDNR,i,j,IOU_INST)*convfac*atwN

     ddep(DDEP_STOCFU,i,j,IOU_INST) = unit_flux(1) ! * to_nmole
     ddep(DDEP_STODFU,i,j,IOU_INST) = unit_flux(2) !* to_nmole
     ddep(DDEP_STOTCU,i,j,IOU_INST) = unit_flux(5) !* to_nmole
     ddep(DDEP_STOMCU,i,j,IOU_INST) = unit_flux(6) !* to_nmole
     ddep(DDEP_STOGRU,i,j,IOU_INST) = unit_flux(10) !* to_nmole


     ddep(DDEP_STOCFL,i,j,IOU_INST) = lai_flux(1) !* to_nmole
     ddep(DDEP_STODFL,i,j,IOU_INST) =  lai_flux(2) !* to_nmole
     ddep(DDEP_STOTCL,i,j,IOU_INST) =  lai_flux(5) !* to_nmole
     ddep(DDEP_STOMCL,i,j,IOU_INST) =  lai_flux(6) !* to_nmole
     ddep(DDEP_STOGRL,i,j,IOU_INST) =  lai_flux(10)! * to_nmole

   !---- end ecosystem specific ----------------------------------------------

  end subroutine  Add_ddep

  end module My_UKDep_ml

