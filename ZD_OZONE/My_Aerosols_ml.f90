 module My_Aerosols_ml

   !----------------------------------------------------------------------
   ! Allows to select aerosol types for the model run:
   ! 1. AERO_DYNAMICS - for running UNI-AERO 
   ! 2. INORGANIC_AEROSOLS - for run old Ammonium routine
   ! 3. RUN_MARS - run MARS eq model
   ! 4. RUN_EQSAM - run EQSAM eq model ! One of 2.,3. or 4. must be true
   ! 5. ORGANIC_AEROSOLS - for including Secondary Organic Aerosol
   !----------------------------------------------------------------------

!st.. Made from the previous My_Aerosols, where Ammonium module is moved to
!st.. Ammonium_ml.f90, and OrganicAerosol module is moved to SOA_ml.f90

   implicit none

   !/-- public           !!  true if wanted
                    
    logical, public, parameter :: AERO_DYNAMICS      = .false.   &  
                                , INORGANIC_AEROSOLS = .false.  & !old Ammonium stuff
                                , RUN_MARS           = .false. & !MARS
                                , RUN_EQSAM          = .true. & !EQSAM
                                , ORGANIC_AEROSOLS   = .false.   
!stDep
    integer, public, parameter :: NSIZE = 2    ! number of aerosol sizes (1-fine, 2-coarse)

contains

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine My_MARS(deb)

    !..................................................................
    ! Different pre-processing (and for aerosol dynamisc: fractioning) 
    ! is needed for UNI-ACID,UNI-OZONE and UNI-AERO. 
    !
    ! hf dec-2002
    ! changed eg aNH4->aNH4out in call to rpmares
    ! Should put k-loop into rpmares
    ! Should change name rpmares to MARS
    !..................................................................

 use Setup_1dfields_ml,  only :  xn_2d     ! SIA concentration 
 use GenSpec_tot_ml,     only :  NH3, HNO3, SO4, aNO3, aNH4, NO3
 use Setup_1dfields_ml,  only :  temp, rh
 use ModelConstants_ml,  only :  KMAX_MID, KCHEMTOP   
 use GenChemicals_ml,    only :  species
 use PhysicalConstants_ml, only : AVOG
 use MARS_ml, only: rpmares

 implicit none
 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  

 logical, intent(in)  :: deb


 !.. local
  real    :: so4in, no3in, nh4in, hno3in, nh3in,   &
             aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out,   &
             coef
  integer :: k, ic, bin, spec, errmark
  logical :: debsub
 !-----------------------------------

   coef = 1.e12 / AVOG

   do k = KCHEMTOP, KMAX_MID
  
!//.... molec/cm3 -> ug/m3
      so4in  = xn_2d(SO4,k) * species(SO4)%molwt  *coef
      hno3in = xn_2d(HNO3,k)* species(HNO3)%molwt *coef 
      nh3in  = xn_2d(NH3,k) * species(NH3)%molwt  *coef
      no3in  = xn_2d(aNO3,k) * species(aNO3)%molwt  *coef ! + &
!               xn_2d(NO3,k) * species(NO3)%molwt  *coef   !different for UNI-OZONE/ACID
      nh4in  = xn_2d(aNH4,k) * species(aNH4)%molwt  *coef

 !--------------------------------------------------------------------------                
    call rpmares (so4in, hno3in,no3in ,nh3in, nh4in , rh(k), temp(k),   &
                  aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, ERRMARK,debsub) 
 !--------------------------------------------------------------------------
 ! SO4 is not changed so do not need to be reset
     
!      xn_2d(NO3,k)  = FLOOR !different for UNI-OZONE/ACID

      xn_2d(HNO3,k)  = max (FLOOR, gNO3out / (species(HNO3)%molwt *coef) )
      xn_2d(NH3,k)   = max (FLOOR, gNH3out / (species(NH3)%molwt  *coef) )
      xn_2d(aNO3,k)  = max (FLOOR, aNO3out / (species(aNO3)%molwt  *coef) )
      xn_2d(aNH4,k)  = max (FLOOR, aNH4out / (species(aNH4)%molwt  *coef) )

   enddo  ! K-levels

 end subroutine My_MARS

!.............................

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine My_EQSAM(debug_cell)

    !..................................................................
    ! Different pre-processing (and for aerosol dynamisc: fractioning) 
    ! is needed for UNI-ACID and UNI-AERO. UNI-OZONE and UNI-ACID can 
    ! use the same.
    !
    ! hf dec-2002
    ! Version aero_martijn included (same as in LOTUS)
    ! 
    ! 
    !..................................................................

 use EQSAM_v03d_ml,      only :  eqsam_v03d
 use Setup_1dfields_ml,  only :  xn_2d     ! SIA concentration 
 use GenSpec_tot_ml,     only :  NH3, HNO3, SO4, aNO3, aNH4,NO3
 use Setup_1dfields_ml,  only :  temp, rh,pp
 use ModelConstants_ml,  only :  KMAX_MID, KCHEMTOP   
 use PhysicalConstants_ml, only : AVOG
! use EQSAM_ml, only: aero

 implicit none
 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  

 logical, intent(in)  :: debug_cell


 !.. local
  real    :: so4in(KCHEMTOP:KMAX_MID), &
             no3in(KCHEMTOP:KMAX_MID), &
             nh4in(KCHEMTOP:KMAX_MID), &
             hno3in(KCHEMTOP:KMAX_MID), &
             nh3in(KCHEMTOP:KMAX_MID),   &
             aSO4out(KCHEMTOP:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), & 
             gNH3out(KCHEMTOP:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID), &
!HF tmp fix
             NAin(KCHEMTOP:KMAX_MID)  , &
             CLin(KCHEMTOP:KMAX_MID) , &
             aNAout(KCHEMTOP:KMAX_MID),&
             aCLout(KCHEMTOP:KMAX_MID),&
             gCLout(KCHEMTOP:KMAX_MID),&
             gSO4out(KCHEMTOP:KMAX_MID)

  integer :: i,j,k, errmark
  logical :: debsub = .false.
 !-----------------------------------


  if ( debsub .and. debug_cell ) then ! Selected debug cell
    write(*,*)'Before EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(aNO3,20),xn_2d(aNH4,20)
  endif

!//.... molec/cm3 -> micromoles/m**3
      so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
      no3in(KCHEMTOP:KMAX_MID)  = xn_2d(aNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG !+&
!                                  xn_2d(NO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG !different for ACID/OZONE
      nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(aNH4,KCHEMTOP:KMAX_MID)*1.e12/AVOG

      NAin(KCHEMTOP:KMAX_MID)  = 0.
      CLin(KCHEMTOP:KMAX_MID)  = 0.

 !--------------------------------------------------------------------------                
 !   call aero (so4in, hno3in,no3in ,nh3in, nh4in , rh, temp,   &
 !                 aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, ERRMARK,debsub) 
 !--------------------------------------------------------------------------
 !--------------------------------------------------------------------------                
  
   call eqsam_v03d (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh,temp,pp,   &
                    aSO4out, aNO3out, aNH4out, aNAout, aCLout,               &
                    gSO4out, gNH3out, gNO3out, gClout, aH2Oout)
 
 !--------------------------------------------------------------------------
 ! SO4 is not changed so do not need to be reset 

!      if (so4in .ne. aSO4out) then
!      write(*,*) 'STARNGE EQSAM RESULTS'
!      endif

!//.... micromoles/m**3  -> molec/cm3 
!      xn_2d(NO3,KCHEMTOP:KMAX_MID)  = FLOOR !different for ACID/OZONE

      xn_2d(HNO3,KCHEMTOP:KMAX_MID)  = max(FLOOR,gNO3out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )
      xn_2d(NH3,KCHEMTOP:KMAX_MID)   = max(FLOOR,gNH3out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )
      xn_2d(aNO3,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNO3out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 ) 
      xn_2d(aNH4,KCHEMTOP:KMAX_MID)  = max(FLOOR,aNH4out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )
      xn_2d(SO4,KCHEMTOP:KMAX_MID)  = max(FLOOR,aSO4out(KCHEMTOP:KMAX_MID)*AVOG/1.e12 )

 if ( debsub .and. debug_cell ) then ! Selected debug cell
    write(*,*)'After EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
               xn_2d(NH3,20),xn_2d(aNO3,20),xn_2d(aNH4,20)
  endif

 end subroutine My_EQSAM

!.............................

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


 !water 
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      subroutine Aero_water(i,j)

    !..................................................................
    !EQSAM is called just before writing out to calculate aerosol water
    !simulating "PM samples equilibration" at T=20C and Rh = 50%
    !..................................................................
 use EQSAM_v03d_ml,      only :  eqsam_v03d
 use Setup_1dfields_ml,  only :  xn_2d      ! SIA concentration 
 use Chemfields_ml,      only :  PM_water               !PMwater  
 use GenSpec_tot_ml,     only :  NH3, HNO3, SO4, aNO3, aNH4
 use Setup_1dfields_ml,  only :  temp, rh,pp
 use ModelConstants_ml,  only :  KMAX_MID, KCHEMTOP   
 use PhysicalConstants_ml, only : AVOG
! use EQSAM_ml, only: aero

 implicit none
 real, parameter ::    FLOOR = 1.0E-30         ! minimum concentration  

 integer, intent(in)  :: i, j

 !.. local
  real    :: so4in(KCHEMTOP:KMAX_MID), &
             no3in(KCHEMTOP:KMAX_MID), &
             nh4in(KCHEMTOP:KMAX_MID), &
             hno3in(KCHEMTOP:KMAX_MID), &
             nh3in(KCHEMTOP:KMAX_MID),   &
             aSO4out(KCHEMTOP:KMAX_MID), &
             aNO3out(KCHEMTOP:KMAX_MID), &
             aH2Oout(KCHEMTOP:KMAX_MID), &
             aNH4out(KCHEMTOP:KMAX_MID), & 
             gNH3out(KCHEMTOP:KMAX_MID), &
             gNO3out(KCHEMTOP:KMAX_MID), &
!HF tmp fix
             NAin(KCHEMTOP:KMAX_MID)  , &
             CLin(KCHEMTOP:KMAX_MID) , &
             aNAout(KCHEMTOP:KMAX_MID),&
             aCLout(KCHEMTOP:KMAX_MID),&
             gCLout(KCHEMTOP:KMAX_MID),&
             gSO4out(KCHEMTOP:KMAX_MID), &

             rh50(KCHEMTOP:KMAX_MID),t20(KCHEMTOP:KMAX_MID)

  integer :: i,j,k, errmark
  logical :: debsub = .false.
 !-----------------------------------


!  if ( debsub .and. debug_cell ) then ! Selected debug cell
!    write(*,*)'Before EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
!               xn_2d(NH3,20),xn_2d(aNO3,20),xn_2d(aNH4,20)
!  endif

!//.... molec/cm3 -> micromoles/m**3
      so4in(KCHEMTOP:KMAX_MID)  = xn_2d(SO4,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      hno3in(KCHEMTOP:KMAX_MID) = xn_2d(HNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      nh3in(KCHEMTOP:KMAX_MID)  = xn_2d(NH3,KCHEMTOP:KMAX_MID)*1.e12/AVOG 
      no3in(KCHEMTOP:KMAX_MID)  = xn_2d(aNO3,KCHEMTOP:KMAX_MID)*1.e12/AVOG
      nh4in(KCHEMTOP:KMAX_MID)  = xn_2d(aNH4,KCHEMTOP:KMAX_MID)*1.e12/AVOG

      NAin(KCHEMTOP:KMAX_MID)  = 0.
      CLin(KCHEMTOP:KMAX_MID)  = 0.

      rh50(:) = 0.5
      t20(:)  = 283.

 !--------------------------------------------------------------------------                
 !   call aero (so4in, hno3in,no3in ,nh3in, nh4in , rh50, temp20,   &
 !                 aSO4out, aNO3out, aH2Oout, aNH4out, gNH3out, gNO3out, ERRMARK,debsub) 
 !--------------------------------------------------------------------------
 !--------------------------------------------------------------------------                
  
   call eqsam_v03d (so4in, hno3in,no3in,nh3in,nh4in,NAin,CLin, rh50,t20,pp,   &
                    aSO4out, aNO3out, aNH4out, aNAout, aCLout,               &
                    gSO4out, gNH3out, gNO3out, gClout, aH2Oout)
 
 !--------------------------------------------------------------------------

!//....aerosol water (ug/m**3 

      PM_water(i,j,KCHEMTOP:KMAX_MID) = max(0., aH2Oout(KCHEMTOP:KMAX_MID) )

! if ( debsub .and. debug_cell ) then ! Selected debug cell
!    write(*,*)'After EQSAM',xn_2d(SO4,20),xn_2d(HNO3,20),&
!               xn_2d(NH3,20),xn_2d(aNO3,20),xn_2d(aNH4,20)
!  endif

 end subroutine  Aero_water
!.............................

 end module My_Aerosols_ml






