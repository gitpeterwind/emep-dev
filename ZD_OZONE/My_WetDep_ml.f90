module My_WetDep_ml
  use MassBudget_ml,     only : totwdep
  use ModelConstants_ml, only : atwS, atwN, atwPM
!ds  use My_Derived_ml    , only : NWDEP, WDEP_SOX, WDEP_OXN, WDEP_RDN &
!ds                                ,IOU_INST &   ! Index: instantaneous values
!ds                                ,wdep         ! Wet deposition fields

 !ds use My_Derived_ml, only : WDEP_USED,D2_USED   !ds NEW system 16/12/2003
 use Derived_ml,    only : f_wdep, wdep,  &   !ds NEW system 16/12/2003
                           f_2d,   d_2d,  &
                           find_one_index, IOU_INST


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
  

  integer, public, parameter :: NWETDEP =  13  !SeaS 11  ! Number of solublity classes
  type(WScav), public, dimension(NWETDEP), save  :: WetDep
  
 !ds NEW 16/12/2003:

  integer, public, save  :: WDEP_PREC   ! Used in Aqueous_ml
  integer, private, save :: WDEP_SOX, WDEP_OXN, WDEP_RDN

contains

  subroutine Init_WetDep()

  !/ INCLOUDFAC is A/v where A is 5.2 m3 kg-1 s-1, !  and v is the fallspeed (5 m/s). 
    real, parameter ::  FALLSPEED = 5.0                            ! m/s 
    real, parameter ::  SUBCLFAC = 5.2 / FALLSPEED

  !/ e is the scavenging efficiency (0.1 for fine particles, 0.4 for course)

    real, parameter ::  EFF25 = 0.1*SUBCLFAC  & 
                      , EFFCO = 0.4*SUBCLFAC  ! collection efficiency b/clouds - coarse

   !/.. setup the scavenging ratios for in-cloud and sub-cloud. For
   !    gases, sub-cloud = 0.5 * incloud. For particles, sub-cloud=
   !    efficiency * INCLOUDFAC. See also notes in Aqueous_ml.

   !/..                        W_Sca  W_sub
  
    WetDep(1)   = WScav(SO2,    0.3,  0.15)   ! Berge+Jakobsen, issh
    WetDep(2)   = WScav(SO4,    1.0,  EFF25)  ! Berge+Jokobsen, jej
    WetDep(3)   = WScav(aNH4,   1.0,  EFF25)
    WetDep(4)   = WScav(NH3,    1.4,  0.5 )  ! ds, subcloud = 1/3 of cloud for gases
    WetDep(5)   = WScav(aNO3,   1.0,  EFF25)  
    WetDep(6)   = WScav(HNO3,   1.4,  0.5)   ! ds
    WetDep(7)   = WScav(H2O2,   1.4,  0.5)   ! jej, maybe should be 0.6*0.7??
    WetDep(8)   = WScav(HCHO,   0.1,  0.03)  ! jej, plus ds 1/3 rule
    WetDep(9)   = WScav(pNO3,   1.0,  EFFCO) !!dstest rv1_8b  ds, from Svetlana's PMco stuff
    WetDep(10)  = WScav(PM25,   1.0,  EFF25)
    WetDep(11)  = WScav(PMCO,   1.0,  EFFCO)
    WetDep(12)  = WScav(SSFI,   1.0,  EFF25)   !SeaS
    WetDep(13)  = WScav(SSCO,   1.0,  EFFCO)   !SeaS

   !####################### ds NEW define indices here ##########

     WDEP_PREC= find_one_index("WDEP_PREC",f_wdep(:)%name)
     WDEP_SOX = find_one_index("WDEP_SOX",f_wdep(:)%name)
     WDEP_OXN = find_one_index("WDEP_OXN",f_wdep(:)%name)
     WDEP_RDN = find_one_index("WDEP_RDN",f_wdep(:)%name)
   !####################### ds END define indices here ##########

  end subroutine Init_WetDep

  subroutine WetDep_Budget(i,j,sumloss,invgridarea)
     integer,            intent(in) ::  i,j
     real, dimension(:), intent(in) :: sumloss
     real, intent(in)  :: invgridarea
     real :: wdeps, wdepox, wdepred, wdeppm25, wdeppmco


      !wdeps = sumloss(SO2) + sumloss(SO4) + sumloss(AMSU)
       wdeps = sumloss(1) + sumloss(2) !hf+ sumloss(3)

      !wdepred = sumloss(NH3)  + sumloss(AMNI) & 
      !                    +1.5 * sumloss(AMSU) !!! (NH4)1.5 SO4
       wdepred = sumloss(4)  + sumloss(3) !hf& 
                          !hf +1.5 * sumloss(3) !!! (NH4)1.5 SO4
  
      !wdepox  = sumloss(HNO3) + sumloss(AMNI) + pNO3
       wdepox  = sumloss(6) + sumloss(5) + sumloss(9)
      wdeppm25= sumloss(7) 
      wdeppmco= sumloss(8) 

       totwdep(IXADV_SO4)  = totwdep(IXADV_SO4) + wdeps
       totwdep(IXADV_HNO3) = totwdep(IXADV_HNO3) + wdepox
       totwdep(IXADV_NH3)  = totwdep(IXADV_NH3)  + wdepred
       totwdep(IXADV_PM25)  = totwdep(IXADV_PM25)  + wdeppm25
       totwdep(IXADV_PMco)  = totwdep(IXADV_PMco)  + wdeppmco


       wdep(WDEP_SOX,i,j,IOU_INST) = wdeps * atwS * invgridarea 
       wdep(WDEP_OXN,i,j,IOU_INST) = wdepox * atwN * invgridarea 
       wdep(WDEP_RDN,i,j,IOU_INST) = wdepred * atwN * invgridarea 


  end subroutine WetDep_Budget
end module My_WetDep_ml
