module My_WetDep_ml
  use MassBudget_ml,     only : totwdep
  use ModelConstants_ml, only : atwS, atwN, atwPM
  use My_Derived_ml    , only : NWDEP, WDEP_SOX, WDEP_OXN, WDEP_RDN &
                                ,IOU_INST &   ! Index: instantaneous values
                                ,wdep         ! Wet deposition fields
  use GenSpec_tot_ml          ! SO2, SO4, etc.
  use GenSpec_adv_ml          ! IXADV_SO2, IXADV_SO4, etc.
  implicit none
  private

  public :: Init_WetDep       ! Call from Unimod
  public :: WetDep_Budget     ! Call from Aqueous_ml

  type, public :: WScav
     integer  :: adv
     real  :: W_sca       ! Scavenging ratio/z_Sca/rho = W_sca/1.0e6
     real  :: W_sub       ! same for subcloud
  end type WScav
  

  integer, public, parameter :: NWETDEP = 9  ! Number of solublity classes
  type(WScav), public, dimension(NWETDEP), save  :: WetDep
  
contains

  subroutine Init_WetDep()

  !/ INCLOUDFAC is A/v where A is 5.2 m3 kg-1 s-1,   and v is the 
  !                                            fallspeed (5 m/s). 
    real, parameter ::  FALLSPEED = 5.0                   ! m/s 
    real, parameter ::  INCLOUDFAC = 5.2 / FALLSPEED

  !/ e is the scavenging efficiency (0.1 for fine particles, 0.4 for course)

    real, parameter ::  EFF25 = 0.1*INCLOUDFAC  & 
                      , EFFCO = 0.4*INCLOUDFAC  ! collection efficiency
                                                ! below clouds - coarse

   !/.. setup the scavenging ratios for in-cloud and sub-cloud. For
   !    gases, sub-cloud = 0.5 * incloud. For particles, sub-cloud=
   !    efficiency * INCLOUDFAC. See also notes in Aqueous_ml.

   !/..                        W_Sca  W_sub
  
    WetDep(1)   = WScav(SO2,    0.3,  0.15)   ! Berge+Jakobsen, issh
    WetDep(2)   = WScav(SO4,    0.7,  EFF25)  ! Berge+Jokobsen, jej
!hf amsu    WetDep(3)   = WScav(AMSU,   0.7,  EFF25)
    WetDep(3)   = WScav(aNH4,   0.7,  EFF25)
    WetDep(4)   = WScav(NH3,    0.7,  0.35)
!hf amsu    WetDep(5)   = WScav(AMNI,   0.7,  EFF25)  
    WetDep(5)   = WScav(aNO3,   0.7,  EFF25)  
    WetDep(6)   = WScav(HNO3,   0.7,  0.35)
!OZ WetDep(7)   = WScav(H2O2,   0.6,  0.3)   ! jej, maybe should be 0.6*0.7??
!OZ WetDep(8)   = WScav(HCHO,   0.1,  0.05)  ! jej
    WetDep(7)   = WScav(PM25,   0.7,  EFF25)
    WetDep(8)   = WScav(PMCO,   0.7,  EFFCO)
    WetDep(9)   = WScav(pNO3,   0.7,  EFFCO)

  end subroutine Init_WetDep

  subroutine WetDep_Budget(i,j,sumloss,invgridarea)
     integer,            intent(in) ::  i,j
     real, dimension(:), intent(in) :: sumloss
     real, intent(in)  :: invgridarea
     real :: wdeps, wdepox, wdepred


      !wdeps = sumloss(SO2) + sumloss(SO4) + sumloss(AMSU)
       wdeps = sumloss(1) + sumloss(2)  !hf amsu + sumloss(3)

      !wdepred = sumloss(NH3)  + sumloss(AMNI) & 
      !                    +1.5 * sumloss(AMSU) !!! (NH4)1.5 SO4
!hf amsu       wdepred = sumloss(4)  + sumloss(5) & 
!hf amsu                           +1.5 * sumloss(3) !!! (NH4)1.5 SO4
         wdepred = sumloss(3)  + sumloss(4) 

      !wdepox  = sumloss(HNO3) + sumloss(AMNI) + pNO3
       wdepox  = sumloss(6) + sumloss(5) + sumloss(9)


       totwdep(IXADV_SO4)  = totwdep(IXADV_SO4) + wdeps
       totwdep(IXADV_HNO3) = totwdep(IXADV_HNO3) + wdepox
       totwdep(IXADV_NH3)  = totwdep(IXADV_NH3)  + wdepred


       wdep(WDEP_SOX,i,j,IOU_INST) = wdeps * atwS * invgridarea 
       wdep(WDEP_OXN,i,j,IOU_INST) = wdepox * atwN * invgridarea 
       wdep(WDEP_RDN,i,j,IOU_INST) = wdepred * atwN * invgridarea 

  end subroutine WetDep_Budget
end module My_WetDep_ml
