module My_DryDep_ml
!+
! Module to define the dry deposition components and rates. We
! define the min (vd_min) and max dep. velocity (Vg) first and then derive the
! daytime addition (vd_day).

!Bugs found:
! PAN Vg not defined
! coszen used for day/night, but this gives v. low values for high latitudes
! z_mid used instead of z_bnd for height of lowest layer
! 


  use My_Derived_ml , only : DDEP_SOX,DDEP_OXN,DDEP_RDN, &
                             IOU_INST    &!updates inst. dep. fields
                           , ddep         ! 2d fields
  use GenSpec_adv_ml , only: NSPEC_ADV &
                   ,IXADV_SO2,IXADV_NO2  &
                   ,IXADV_HNO3 &  !! ,IXADV_EC  ,IXADV_OC & 
    ,  IXADV_PAN ,IXADV_SO4,IXADV_NH3,IXADV_AMNI, IXADV_AMSU
 use ModelConstants_ml , only : atwS, atwN
 implicit none
 private

  public :: Init_vd
  public :: Add_ddep

  ! Indices from land-use database:

   integer, public, parameter :: SEA=0, ICE=1, TUNDRA=2, AGRIC=3, &
                         UNDEF=4, FOREST=5, DESERT=6

  ! here we define the minimum set of species which has different
  ! deposition velocities. We calculate Vg for these, and then
  ! can use the rates for other similar species. (e.g. AMSU can use
  ! the Vg for SO4.  Must set NDRYDEP species

  integer, public, parameter ::  NDRYDEP_CALC = 5

  ! .. CDEP = calculated dep:

  integer, public, parameter ::  &
    CDEP_HNO3 = 1, CDEP_NO2 = 2, &
    CDEP_SO2= 3, CDEP_SO4  = 4, CDEP_PAN = 5

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
   end type depmap

   type(depmap), public, dimension(NDRYDEP_ADV):: Dep

   real, public, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost


  !/-- we define two velocities, vd_day and vd_min, in order to allow
  !    diurnal variations (Vd). In subroutine Init_vd an
  !    additional set of values, vd_max is used. vd_day is the
  !    difference between vd_max and vd_min. For some components the max 
  !    and min will be identical, whereas for others (e.g. ozone) the
  !    max will be significantly higher than the min.

   real, public, save, dimension(NDRYDEP_CALC,0:6) :: &
       vd_day,      & ! Extra (noontime-ish)  Vd above the min Vd
       vd_min         ! min. (nightime-ish)  Vd.

   logical, private, parameter :: MY_DEBUG = .false.

contains
  subroutine Init_vd

!------------------------------------------------------------------------------
!------- DEPOSITION VELOCITIES, FROM HOUGH JGR,96, D4, pp 7325-7362, 1991 -----
!------- OR SEINFELD+PANDIS OR EMEP REPORTS 2/93 OR 2/2001 
!------- 6 landclasses in the model: 0 = ocean, 1 = ice, 2 = tundra, 
!        3 = agriculture, 4 = undef, 5 = forest, 6 = desert 
!------------------------------------------------------------------------------
   real    :: cm_s = 0.01   !  Converts cm/s to m/s
   real, dimension(NDRYDEP_CALC,0:6) :: vd_max

   vd_max(:,:) = 0.0
   vd_min(:,:) = 0.0
!
!  The deposition array is declared as vd(NSPEC_ADV,0:6)


!..2) HNO3
        vd_max(CDEP_HNO3,SEA)          =  1.0 * cm_s    ! Seinfeld, p969
        vd_max(CDEP_HNO3,ICE:TUNDRA)   =  0.5 * cm_s    !  " ", guessed for tundra
        vd_max(CDEP_HNO3,AGRIC:FOREST) =  4.0 * cm_s    ! Seinfeld, p969
        vd_max(CDEP_HNO3,DESERT)       =  2.0 * cm_s    ! ?

        vd_min (CDEP_HNO3,:)           = vd_max(CDEP_HNO3,:)  ! default

!..5) NO2

        vd_max(CDEP_NO2,SEA)          =  0.02 * cm_s  ! Seinfeld, p969
        vd_max(CDEP_NO2,ICE)          = 0.002 * cm_s  ! Seinfeld, p969
        vd_max(CDEP_NO2,TUNDRA)       = 0.2  * cm_s   ! TROTREP1c
        vd_max(CDEP_NO2,AGRIC:FOREST) = 0.5  * cm_s   ! TROTREP1a EMEP 2/93
        vd_max(CDEP_NO2,DESERT)       = 0.02 * cm_s      

        vd_min(CDEP_NO2,:)            = vd_max(CDEP_NO2,:)  ! default
        !u6 ERROR vd_min(CDEP_NO2,AGRIC:FOREST) = 0.5 * cm_s
        vd_min(CDEP_NO2,TUNDRA:FOREST) = 0.06 * cm_s

!..6) SO2
        vd_max(CDEP_SO2,:)      = 1.0 * cm_s       
        vd_max(CDEP_SO2,ICE)    = 0.1 * cm_s       ! Guess....
        vd_max(CDEP_SO2,DESERT) = 0.1 * cm_s       ! Guess....

        vd_min(CDEP_SO2,:)            = vd_max(CDEP_SO2,:)  ! default
        vd_min(CDEP_SO2,TUNDRA:FOREST) = 0.3 *  cm_s

!..7) SO4
	vd_max(CDEP_SO4,:)            =  0.1 * cm_s   ! Small for particles
        vd_min(CDEP_SO4,:)            =  vd_max(CDEP_SO4,:)  ! default

!..8) PAN
        vd_max(CDEP_PAN,:) = 0.5 * vd_max(CDEP_NO2,:)
        vd_min(CDEP_PAN,:) = 0.5 * vd_min(CDEP_NO2,:)

! ... Now find daytime extra by subtraction:

       vd_day(:,:) = vd_max(:,:) - vd_min(:,:)

       if ( any (vd_min<0.0) .or. any(vd_day<0.0)  ) then
          print *, "ERROR IN VDINIT!!!! VD_MAX: ", vd_max
          print *, "ERROR IN VDINIT!!!! VD_MIN: ", vd_min
          print *, "ERROR IN VDINIT!!!! VD_DAY: ", vd_day
       end if

 ! .... Now, we define the mapping between the advected species and
 !      the specied for which the calculation needs to be done.

   Dep(1) =  depmap( IXADV_HNO3 , CDEP_HNO3)
   Dep(2) =  depmap( IXADV_PAN,   CDEP_PAN)   !u7
   Dep(3) =  depmap( IXADV_NO2,   CDEP_NO2)
   Dep(4) =  depmap( IXADV_SO2,   CDEP_SO2)
   Dep(5) =  depmap( IXADV_SO4,   CDEP_SO4)
   Dep(6) =  depmap( IXADV_NH3,   CDEP_SO2)
   Dep(7)=  depmap( IXADV_AMSU,  CDEP_SO4)
   Dep(8)=  depmap( IXADV_AMNI,  CDEP_SO4)

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

end module My_DryDep_ml

