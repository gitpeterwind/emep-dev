!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module My_Emis_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
implicit none

   !-----------------  "my" emissions choices    ---------------------------!
   !  Here we define the emissions stuff likely to change between different !
   !  model versions - e.g. the number of emission files. The consistency   !
   !  of some of these choices will be checked in the Emissions_ml module   !
   !  Includes setting of biogenic emissions and AIRNOX (lightning+aircraft)!
   !------------------------------------------------------------------------!

   !-- contains subroutines:

      public :: set_molwts     ! Called from Emissions_ml

   !/-- model-type
   !    for consistency checks, possibly to match label in My_model_ml??

       character(len=12), public :: MY_MODEL = "emepds"

   !/** emissions. 
   !
   !   NEMIS_PLAIN  gives the number of emission file used for compounds which
   !   do not have speciation (typically, so2, nox).
   !   NEMIS_SPLIT gives the number of emission file used for compounds which
   !   do have speciation (typically, voc, pm).
   !
   !   NRCEMIS: the total number of species where emissions enter into the 
   !   rate-coefficients. This should be equal to the sum of
   !   NEMIS_PLAIN + sum(EMIS_SPLIT). ( Checked in consistency check )
   !   ------------------------------------------------------------------------

   integer, public, parameter :: &
         NEMIS_PLAIN =  6   & ! No. emission files to be read for non-speciated
       , NEMIS_SPLIT =  1   & ! No. emission files to be read for speciated
       , NRCEMIS     = 16     ! No. chemical species with emissions 

   integer, public, parameter :: &   ! ** derived ** shouldn't need to change:
         NEMIS       =      & ! Sum of the above - all emissions
           NEMIS_PLAIN + NEMIS_SPLIT   &
       , NRCSPLIT    =      & ! No. species from speciated (split) compounds
           NRCEMIS - NEMIS_PLAIN  

   !/** The names used below must have length of 6 characters and must
   !    belong to the full list given in Emissions_ml, as they will
   !    be used together with information in that module. An error
   !    message (fatal)  will be produced if this is not the case.
   !-----------------------------------------------------------------

    character(len=6), public, save, dimension(NEMIS) :: &
      EMIS_NAME  = &
      (/ "sox   ", "nox   ", "co    "   &   ! =non-split first
       , "nh3   ", "pm25  ", "pmco  "   &
       , "voc   "  /)                       ! =to be split

    character(len=6), public, save, dimension(NEMIS_SPLIT) :: &
      SPLIT_NAME = &
       (/ "voc   "  /)
 !! for SOA      (/ "voc   ", "pm25  "  /)

    integer, public, save, dimension(NEMIS_SPLIT) :: &
      EMIS_NSPLIT  = &
       (/  10    /)    !!!! (check - excluding bio?)
 !! for SOA       (/  10   ,      3  /)    !!!! (check - excluding bio?)

    !/-- and now  join the above name arrays  to make the complete list:

!    character(len=6), public, save, dimension(NEMIS) :: &
!      EMIS_NAME = &
!      (/ (EMIS_PLAIN(1:NEMIS_PLAIN)), &
!         (EMIS_SPLIT(1:NEMIS_SPLIT))     /)


!... define integers for emitted species which are entered into the chemistry 
!    calculation. The index Qxxx is used for the species which have gridded 
!    emissions, e.g. so2 or total VOC. The index QRCxxx is then used to allow 
!    for emissions of VOC-split species,e.g.  C2H6.
!
!    nb species such as so4 and no2 are not included here as they are
!       derived in setup as simple fractions of the so2 and no emissions
!       (Of course, they could be defined as SPLIT above and then they
!       should be included).

   integer, public, parameter ::   &
           QRCSO2 =   1      & ! IQSO2   &   ! 1
         , QRCNO  =   2      & ! IQNOX   &   ! 2
         , QRCCO  =   3      & ! IQCO        ! 4
         , QRCNH3 =   4      & ! IQCO        ! 4
         , QRCPM25=   5      & ! IQSO2   &   ! 1
         , QRCPMCO=   6      &
      !/**now we deal with the emissions which are split,e.g.VOC
      !  ******************************************************
      !  **** must be in same order as EMIS_SPLIT array **** **
      !  **** ds rv1.8.4 bug-fix:
      !  **** AND vocsplit.defaults file !!!!!!!  ***** **** **
      !  ******************************************************
         , QRCC2H6    = 7      & 
         , QRCNC4H10 =  8    & 
         , QRCC2H4    = 9      & 
         , QRCC3H6    =10      & 
         , QRCOXYL   = 11    & 
         , QRCHCHO   = 12    & 
         , QRCCH3CHO = 13   & 
         , QRCMEK    = 14   & 
         , QRCC2H5OH = 15    & 
         , QRCCH3OH  = 16

      !ds CHanged from:
      !  , QRCC2H4  = 7      & ! MACHDS 
      !  , QRCC2H6  = 8      & ! MACHDS
      !  , QRCC3H6  = 9      & ! MACHDS
      !  , QRCNC4H10 = 10     &  ! MACHDS
      !  , QRCOXYL   = 11    &  ! MACHDS
      !  , QRCC2H5OH = 12    &  ! MACHDS
      !  , QRCHCHO   = 13    &  ! MACHDS
      !  , QRCCH3CHO  = 14    &  ! MACHDS
      !  , QRCCH3OH   = 15    &  ! MACHDS
      !  , QRCMEK     = 16      ! MACHDS

 !! for SOA          , QRCOC      = 15    &  ! MACHDS
 !! for SOA          , QRCEC      = 16    &  ! MACHDS
 !! for SOA          , QRCINORG   = 17       ! MACHDS

  ! Biogenics

   integer, public, parameter ::   NFORESTVOC = 2   
   character(len=8),public, save, dimension(NFORESTVOC) :: &
                                   FORESTVOC = (/ "isoprene","terpene "/)   
   integer, public, parameter ::   & 
               QRCISOP    = 18     &
              ,QRCTERP    = 19
!SeaS
   integer, public, parameter ::  NSS   = 2 &   ! number of sea salt size modes
                                 ,QSSFI = 1 &   ! production of fine SS
                                 ,QSSCO = 2     ! production of coarse SS  

    real, public, dimension(NRCEMIS), save  :: molwt ! Molecular weights                                                             

   !/** Lightning and aircraft NOx. QRCAIR is set equal to QRCNO
   ! if AIRNOX is true, otherwise to one. Avoids problems with
   ! dimensions.

    logical, public, parameter :: AIRNOX = .true.   ! Gives NOx emission
    integer, public, parameter :: QRCAIR = QRCNO    ! 
 
   !hf u2
   !/** Volcanos. QRCVOL is set equal to QRCSO2
   ! if VOLCANOES is true, otherwise to one. Avoids problems with 
   ! dimensions

    logical, public, parameter :: VOLCANOES  = .true.  ! Gives Volcanos
    integer, public, parameter :: QRCVOL     = QRCSO2   

  contains
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  !+ set molecular weights for the QRCxx species . 
  subroutine set_molwts()
  !---------------------------------------------------------------------  

  !  
  ! ACID and OZONE ..
        molwt(QRCSO2)   = 32.0  ! Emissions as S
        molwt(QRCNO)    = 14.0  ! Emissions as N
        molwt(QRCNH3)   = 14.0  ! Emissions as N
        molwt(QRCPM25)  = 100.0  !  Fake for PM2.5
        molwt(QRCPMCO)  = 100.0  !  Fake for PM2.5

        molwt(QRCCO )    = 28.0  ! Emissions as N      
        molwt(QRCC2H4)   = 24.0 ! Emissions as C ???   
        molwt(QRCC2H6)   = 24.0 ! Emissions as C ???   
        molwt(QRCC3H6)   = 36.0 ! Emissions as C ???   
        molwt(QRCNC4H10) = 48.0 ! Emissions as C ???   
        molwt(QRCOXYL)   = 106.0 ! Emissions as C ???   
        molwt(QRCC2H5OH) = 46.0 
        molwt(QRCHCHO)   = 30.0 
        molwt(QRCCH3CHO) = 44.0 
        molwt(QRCCH3OH)  = 32.0 
        molwt(QRCMEK)    = 72.0 
!        molwt(QRCOC)    = 150.0 
!        molwt(QRCEC)    = 150.0 
!        molwt(QRCINORG)  = 150.0 
  end subroutine set_molwts
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end module My_Emis_ml
