!>____________________________________________________________________________<
module Aqueous_ml
!>____________________________________________________________________________<
  !+ aqueous scavenging routines.
  !-------------------------------------------------------------------------
  !
  ! ds May 2002
  !   - a "minimal" version of the aqueous reactions.
  !  The scavenging of soluble compounds is based upon the work of Berge 
  !   and Jakobsen, Eliassen+Saltbones.  Simple scavenging coefficients
  !   are used, in order to avoid the supsect HIRLAM cloud-water values.
  !   A distinction is made between in-cloud and sub-cloud scavenging.
  !
  !
  !**USAGE
  !   Uses My_WetDep_ml, where the scavenging coefficients, divided by
  !   relevant factors, are specified.
  !
  !   call Setup_Clouds(i,j)     from Runchem_ml
  !   call WetDeposition(i,j)    from Runchem_ml if prec. clouds_present
  !
  use My_WetDep_ml, only : WetDep, NWETDEP, WetDep_Budget
  use My_Derived_ml       , only : NWDEP, WDEP_PREC &
                                    ,IOU_INST &   ! Index: instantaneous values
                                    ,wdep         ! Wet deposition fields
  use GridValues_ml       , only : gridwidth_m,xm2,xmd,carea
  use ModelConstants_ml, only: &
      CHEMTMIN, CHEMTMAX       & ! Range of temperature 
     ,KMAX_MID                 & ! ground, k=20
     ,KUPPER                   & ! top of cloud-chemistry, now k=6
     ,KCHEMTOP                 & ! top of chemistry, now k=2
     , dt => dt_advec          & ! model timestep
     ,PT, ATWAIR                 ! Pressure at top, atw. air
  use Met_ml,               only :  pr, roa, z_bnd, cc3d, ps&
!hf cum
                                    ,lwc
  use Par_ml              , only : me   ! for DEBUG
  use Setup_1dfields_ml,    only : xn_2d, amk
  implicit none
  private

  !/-- subroutines
  public ::  init_aqueous
  public :: Setup_Clouds   ! characterises clouds and calls WetDeposition if rain
  public :: WetDeposition   !u7.2, ds -  simplified setup_wetdep
  private ::  tabulate_aqueous
  private ::  get_frac
  private ::  setup_aqurates


 !/ Outputs:
 !-------------------------------------------------------------------
   logical, public,save :: prclouds_present ! true if precipitating clouds

   logical, public, save, dimension(KUPPER:KMAX_MID) :: &
           incloud               ! True for in-cloud k values 


!/-- variables used in module:
   logical, private, save :: DEBUG_AQ = .false.

   real, private, save, dimension(KUPPER:KMAX_MID) :: &
         pr_acc                 ! Accumulated precipitation

   integer, private, save  :: kcloudtop   ! k-level of highest-cloud
   integer, private, save  :: ksubcloud   ! k-level just below cloud

   real, private, parameter :: &      ! Define limits for "cloud"
       PR_LIMIT = 1.0e-7  &   ! for accumulated precipitation
       !hf ,CW_LIMIT = 1.0e-7      ! for cloud water
       !hf lwc is another unit than cw(LWC g/m3, cw kg(H2o)/kg(air))
       ,CW_LIMIT = 1.0e-10   &   ! for cloud water
       ,B_LIMIT = 1.0e-3        ! for cloud cover (fraction)
   real, private, save  :: & !  Set in init below (F wouldn't accept **0.4 in param
        INV_Hplus             & ! = 1.0/Hplus      1/H+,         was:HINRAT
       ,INV_Hplus0p4            ! =INV_Hplus**0.4  (1/H+)**0.4

 ! The Henry's law coefficients, K, given in units of M or M atm-1,
 ! are calculated as effective. A factor K1fac = 1+K1/H+ is defined
 ! here also. 
 !-------------------------------------------------------------------
  integer, public, parameter :: &
  NHENRY  = 3, &       ! No.species with Henry's law applied
  NK1     = 1, &       ! No.species needing effective Henry's calc.
  IH_SO2  = 1, &
  IH_H2O2 = 2, &
  IH_O3   = 3

   ! Aqueous fractions:
    real, public, dimension ( NHENRY, KUPPER:KMAX_MID ), save :: frac_aq

   real, private, dimension (NHENRY, CHEMTMIN:CHEMTMAX), save  :: H
   real, private, dimension (NK1,CHEMTMIN:CHEMTMAX), save  :: K1fac  ! old Hso2eff

 !/ Outputs:
 ! Aqueous reaction rates for usage in gas-phase chemistry
 !-------------------------------------------------------------------
   integer, private, parameter :: NAQUEOUS = 4  ! No. aqueous rates
   integer, private, parameter :: NAQRC    = 3   ! No. constant rates

   real, public, dimension(NAQUEOUS,KCHEMTOP:KMAX_MID), save :: aqrck
   real, private, dimension(NAQRC),        save :: aqrc ! constant rates for
                                                       !  so2 oxidn.
   real, private, dimension(2),        save :: vw ! constant rates for
   logical, public,save :: prclouds_present ! true if precipitating clouds

   integer, public, parameter :: &
      ICLOHSO2  = 1    &   ! for oh + so2 ->                (was IC4055)
     ,ICLRC1    = 2    &   ! for  [h2o2] + [so2]            (was IAQRC1)
     ,ICLRC2    = 3    &   ! for  [o3] + [so2]              (was IAQRC2)
     ,ICLRC3    = 4        ! for [o3] + [o2] (Fe catalytic) (was IAQRC3)

                                                       !  so2 oxidn.


!==============================================================================!

!==============================================================================!
! Notes - Scavenging for Dummies (me - ds!)
!==============================================================================!
! partly re-done:

!  Incloud scavenging: (only dependant on precipitation (not cloud water)
!----------------------

        ! The parameterization of the scavenging of soluble chemical 
        ! components is scaled to the precipitation in each layer. For 
        ! incloud scavenging it is based on the parameterization described 
        ! in \cite{Berge1998}. The incloud scavenging of a soluble component 
        ! X is given by the expression:
        !
        ! Q = X x f x pr_acc x W_sca
        !     ----------------------
        !          Z_sca x rho_water
        !
        ! where pr_acc is the accumulated precipitation in the layer, 
        ! Z_sca is the scavenging depth (scale=1000m) and rho_water is the density of water. 
        ! f is an efficiency  parameter. In the preiovus Aqueous_ml then f was given
        ! by the fractional solubility factor frac_aq. In this version
        ! the module My_WetDep_ml should return the user-defined array 
        ! WetDep%W_sca, with  of W_Sca representing the value 
        ! of f.W_Sca/Z_sca/rho_water


!  Sub cloud scavenging:
!--------------------------

        ! the sub-cloud scavenging distinguihses between particulate and
        ! gas-phase components. The scavenging of gases is calculated
        ! as Q = vw x P, where P is the accumulated precipitation (pr_acc, m)
        ! from all the layers above. (From setup_1d)
  
        ! For particles the scavenging is believed to much less effective,
        ! as they follow the air-current around the droplets (Berge, 1993).
        ! Scavenging for particles is calculated as 
        !    Q = A.e.P/v   
        ! where A is 5.2 m3 kg-1 s-1, v is the fall-speed of the droplets,
        ! 5 m s-1, and eff is the scavenging efficiency, 0.1.

contains
!------------------------------------------------------------------------------
   ! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    subroutine Setup_Clouds(i,j)
   !
   !   Define incloud and precipitating clouds.
   !   The layer must contain at least 1.e-7 kgwater/kg air to be considered 
   !   a cloud. 
   !
   !   Also calculates 
   !    pr_acc     - the accumulated precipitation for each k
   !    b          - fractional cloud cover for each k
   !--------------------------------------------------------------------------

    integer, intent(in) ::  i,j

    !/-local
!u7.2   real, dimension(KUPPER:KMAX_MID) :: &
!u7.2          b                     !  Cloud-area (fraction)
!u7.2    ,cloudwater            !  Cloud-water (volume mixing ratio) 
!u7.2                                !  cloudwater = 1.e-6 same as 1.g m^-3
!hf cum
   real, dimension(KUPPER:KMAX_MID) :: &
          b           &          !  Cloud-area (fraction)
         ,cloudwater             !  Cloud-water (volume mixing ratio) 

    integer         ::  k


! Add up the precipitation in the column
    pr_acc(KUPPER) = sum ( pr(i,j,1:KUPPER) )   ! prec inc. from above 
    do k= KUPPER+1, KMAX_MID
      pr_acc(k) = pr_acc(k-1) + pr(i,j,k)
      pr_acc(k) = max( pr_acc(k), 0.0 ) !u7.2 ds - FIX
    end do

   prclouds_present = .false.  
   if ( pr_acc(KMAX_MID) > PR_LIMIT )  prclouds_present = .true. ! Precipitation
                                                            ! at the surface


  !su  initialise with .false. and 0, 

    incloud(:)  = .false.
!hf SO2aq
    cloudwater(:) = 0.



  !  Loop starting at surface finding the cloud base

    ksubcloud = KMAX_MID+1    ! k-coordinate of sub-cloud limit

    do  k = KMAX_MID, KUPPER, -1  
        !hf if ( cw(i,j,k,1) >  CW_LIMIT ) exit  !  out of loop
        if ( lwc(i,j,k) >  CW_LIMIT ) exit  !  out of loop
        ksubcloud = k
    end do
    if ( ksubcloud == 0 ) return  ! No cloud water found below level 6
                             ! Cloud above level 6 are likely thin cirrus 
                             ! clouds, and if included may need special 
                             ! treatment.     
                             ! ==>  assume no cloud 

     !ds -query - this excludes some pr values from the summation into
     !       wdep(WDEP_PREC). Is this intended??


 !  Define incloud part of the column requiring that both cloud water
 !  and cloud fractions are above limit values 
 !u7.2.  Skip test for cloud_fraction

    kcloudtop = -1              ! k-level of cloud top
    do k = KUPPER, ksubcloud-1
!hf cum
        b(k) = cc3d(i,j,k)
        !u7.2 b(k) = cc3d(i,j,k)  ! ds simplify tests
        !u7.2 if ( cw(i,j,k,1) > CW_LIMIT .and. b(k) > B_LIMIT ) then
        !  Units: kg(w)/kg(air) * kg(air(m^3) / density of water 10^3 kg/m^3
        !  ==>  cloudwater (volume mixing ratio of water to air in cloud 
        !  (when devided bu cloud fraction b )
        !u7.2   cloudwater(k) = 1.0e-3 * cw(i,j,k,1) * roa(i,j,k,1) / b(k)   

 
!hf cum        if ( cw(i,j,k,1) > CW_LIMIT ) then
!hf old HIRLAM cw used        if ( cw(i,j,k,1) > 1.0e-7   .and. b(k) > B_LIMIT ) then
!hf                           cloudwater(k) = 1.0e-3 * cw(i,j,k,1) * roa(i,j,k,1) / b(k)  
        if ( lwc(i,j,k) >  CW_LIMIT ) then
                cloudwater(k) = lwc(i,j,k) / b(k) !value of cloudwater in the cloud fraction of the grid
                incloud(k) = .true.
                if ( kcloudtop < 0 ) kcloudtop   = k
        end if

    end do

    
    if ( prclouds_present .and. kcloudtop == -1 ) then
           if ( DEBUG_AQ ) write(6,"(a20,i3,2i4,3es12.4)") &
               "ERROR prclouds sum_cw", &
                me, i,j, maxval(lwc(i,j,KUPPER:KMAX_MID),1) , &
                 maxval(pr(i,j,:)), pr_acc(KMAX_MID)
           kcloudtop = KUPPER ! for safety
    end if

!  sets up the aqueous phase reaction rates (SO2 oxidation) and the 
!  fractional solubility 

    call setup_aqurates(b ,cloudwater,incloud)
   !u7.2 if ( prclouds_present )  then
   !u7.2        call WetDeposition(i,j)
   !u7.2 end if

      if ( DEBUG_AQ .and. i == 3 .and. j == 3 ) then
        print *, "DEBUG_AQ me presetn", me, prclouds_present
        print "(a15,i3,2i4,es14.4)", "DEBUG_AQ me ", me, &
                      kcloudtop, ksubcloud, pr_acc(KMAX_MID)
      end if
           
          


   end subroutine Setup_Clouds
   !--------------------------------------------------------------------------
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine init_aqueous() 
  use PhysicalConstants_ml, only: AVOG    ! Gas-constant, Avo's No.
   
   !** DESCRIPTION
   ! Calls initial tabulations, sets frac_aq to zero above cloud level, and
   ! sets constant rates.

   !   MTRLIM represents mass transport limitations between the clouds
   !   and the reminder of the grid-box volume. (so2 will be rapidly depleted 
   !   within the clouds, and must be replenished from the surounding cloud 
   !   free volume.
   !     old tests:    MTRLIM = clfr_loc, 0.8, 1.0

   real, parameter ::           &! H+ stuff. Assume ph approx 4.5, then:
           Hplus     = 5.0e-5          ! Hydrogen ion concentration
   real, parameter :: MASSTRLIM = 1.0  ! Mass transport limitation

   ! H+ stuff:

      INV_Hplus     = 1.0/Hplus      ! 1/H+,         was:HINRAT
      INV_Hplus0p4  =INV_Hplus**0.4  ! (1/H+)**0.4


   ! tabulations

     !======================
     call tabulate_aqueous()
     !======================


   ! Constant rates: The rates given in Berge (1993) are in mol-1 l.
   ! These need to be multiplied by 1.0e3/AVOG/Vf,so we perform the
   ! 1.0e3/AVOG scaling here.

   !so2aq + h2o2   ---> so4,             ref: B93, M\"oller 1980

      aqrc(1) = 8.3e5 * 1.0e3/AVOG * MASSTRLIM


   ! (so2aq + hso3-) + H+ + o3 ---> so4, ref: B93, from Martin&Damschen 1981

      aqrc(2) = 1.8e4 * 1.0e3/AVOG * MASSTRLIM


   ! (so2aq + hso3-) + o2 ( + Fe ) --> so4,  ! See documentation below

      aqrc(3) = 3.3e-10  * MASSTRLIM  

   ! Regarding aqrc(3):
   ! catalytic oxidation with Fe. The assumption is that 2% of SIV
   ! is oxidised per hour inside the droplets, corresponding to a conversion
   ! rate of 5.6^-6 (units s^-1  -- Therfore no conversion from mol l^-1)

   ! Ref: Seland, \O. and T. Iversen (1999) A scheme for black carbon and
   ! sulphate aerosols tested in a hemispheric scale, Eulerian dispersion
   ! model. Atm.  Env. Vol. 33, pp 2853 -- 2879.

   ! !   5.6e-6 * 0.5e-6 (liquid water fraction) /8.5e-3 (fso2 at 10 deg C)

   ! I multiply with the assumed liquid water fraction from Seland and Iversen
   ! (0.5e-6) and with an assumed fso2 since the reaction is scaled by the
   ! calculated value for these parameters later.



end subroutine init_aqueous

!------------------------------------------------------------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine tabulate_aqueous()

  !** DESCRIPTION
  ! Tabulates Henry's law coefficients over the temperature range
  ! defined in Tabulations_ml.
  !   For SO2, the effective Henry's law is given by
  ! Heff = H * ( 1 + K1/H+ )  
  ! where k2 is omitted as it is significant only at high pH.
  ! We tabulate also the factor 1+K1/H+ as K1fac.

  real, dimension(CHEMTMIN:CHEMTMAX)  :: t, tfac   ! Temperature, K, & factor
  integer :: i

     t(:)           = (/ ( real(i), i=CHEMTMIN, CHEMTMAX ) /)
     tfac(:)        = 1.0/t(:) -  1.0/298.0

     H (IH_SO2 ,:)  = 1.23    * exp(3020.0*tfac(:) ) 
     H (IH_H2O2,:)  = 7.1e4   * exp(6800.0*tfac(:) )
     H (IH_O3  ,:)  = 1.13e-2 * exp(2300.0*tfac(:) )
!hf     H (IH_HCHO,:)  = 2.97e3  * exp(7194.0*tfac(:) )


     ! Need  effective Henry's coefficient for SO2:
     K1fac(IH_SO2  ,:)  =  &
            ( 1.0 + 1.23e-2 * exp(2010.0*tfac(:) ) * INV_Hplus)


     H (IH_SO2 ,:)  = H(IH_SO2,:) * K1fac(IH_SO2,:)


  end subroutine tabulate_aqueous
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine setup_aqurates(b ,cloudwater,incloud)
  use Setup_1dfields_ml, only : &
      itemp          ! temperature (K)
   
   !**DESCRIPTION
   !  sets the rate-coefficients for thr aqueous-phase reactions

   real, dimension(KUPPER:KMAX_MID) :: &
          b                   & !  Cloud-aread (fraction)
         ,cloudwater            !  Cloud-water
   logical, dimension(KUPPER:KMAX_MID) :: &
           incloud               ! True for in-cloud k values 
   !/Outputs  -> aqurates

   !/-- local
   real, dimension(KUPPER:KMAX_MID) :: &
                      fso2grid  & ! f_aq * b = {f_aq}      (was: fso2loc)
                     ,fso2aq    & ! only so2.h2o part (not hso4-)
                     ,caqh2o2   & ! rate of oxidation of so2 with H2O2
                     ,caqo3     & ! rate of oxidation of so2 with H2O2
                     ,caqsx       ! rate of oxidation of so2 with o2 ( Fe )
     
	integer k


   call get_frac(cloudwater,incloud)     ! => frac_aq
    
!hf initialize
       aqrck(:,:)=0.


!hf Gas phase ox. of SO2 is "default"
!in cloudy air, only the part remaining in gas phase (not dissolved)
!is oxidized 

    aqrck(ICLOHSO2,:) = 1.0 
!    prclouds_present = .false.  

   do k = KUPPER,KMAX_MID
      if ( incloud(k) ) then       !!! Vf > 1.0e-10) !lwc>CW_limit
       fso2grid(k) = b(k) * frac_aq(IH_SO2,k) 
       fso2aq  (k) = fso2grid(k) / K1fac(IH_SO2,itemp(k)) 
       caqh2o2 (k) = aqrc(1) * frac_aq(IH_H2O2,k)/cloudwater(k)
       caqo3   (k) = aqrc(2) * frac_aq(IH_O3,k)  /cloudwater(k)
   !6z - error spotted, jej, caqsx(k) = aqrc(3) * fso2grid(k)/cloudwater(k)  
       caqsx   (k) = aqrc(3) /cloudwater(k)  


    
      ! oh + so2 gas-phase
       aqrck(ICLOHSO2,k) = ( 1.0-fso2grid(k) )   ! Now correction factor!
    
       aqrck(ICLRC1,k)   = caqh2o2(k) * fso2aq(k)
    
       aqrck(ICLRC2,k)   = caqo3(k) * INV_Hplus0p4 * fso2grid(k)
    
       aqrck(ICLRC3,k)   = caqsx(k) *  fso2grid(k) 

     end if

   enddo

 end subroutine setup_aqurates
!------------------------------------------------------------------------------
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
subroutine get_frac(cloudwater,incloud)
  use Setup_1dfields_ml, only : &
      temp           &! temperature (K)
     ,itemp          ! temperature (K)
  use PhysicalConstants_ml, only: RGAS_ATML    ! Gas-constant, Avo's No.

  !** DESCRIPTION
  ! Calculating pH dependant solubility fractions:
  ! Calculates the fraction of each soluble gas in the aqueous phase, frac_aq
  !

  !/-- in from used modules :  cloudwater and logical incloud
  !/-- out to rest of module: frac_aq
  !/-- local
   real, dimension (KUPPER:KMAX_MID) :: &
                      cloudwater         ! Volume fraction  - see notes above. Average for grid or value for cloud part?
   logical, dimension(KUPPER:KMAX_MID) :: &
          incloud               ! True for in-cloud k values 
  real, dimension (KUPPER:KMAX_MID) :: VfRT       ! Vf * Rgas * Temp

  integer :: ih, k   ! index over species with Henry's law, vertical level k


  ! Make sure frac_aq is zero outside clouds

  frac_aq(:,:) = 0.

  do k = KUPPER, KMAX_MID

     ! the old test was clw=1.0e-3cw/clfr_loc > 0.05e-7
     ! ie  cw > 0.05e-4.clfr_loc, dvs cw > 0.05e-7 for b=0.001

     if  ( incloud(k) ) then  


         VfRT(k) = cloudwater(k) * RGAS_ATML * temp(k)

         ! Get aqueous fractions:
         do ih = 1, NHENRY

            frac_aq(ih,k) = 1.0/ ( 1.0+1.0/( H(ih,itemp(k)) * VfRT(k) ) )

         end do

     end if
  end do

end subroutine get_frac
   
!------------------------------------------------------------------------------

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine WetDeposition(i,j)

  !+ Calculates wet deposition and changes in xn concentrations
  !  WetDeposition called from RunChem if precipitation reach the surface

  !  input
    integer, intent(in) ::  i,j

  !  local


    integer :: adv     ! index in IXADV_  arrays
    integer :: spec    !  species index from WetDep array
    integer :: n,k

    real    :: invgridarea     !   xm2/(h*h)
    real    :: f_rho           !   Factors in rho calculation
    real    :: rho(KUPPER:KMAX_MID)
    real, dimension(KUPPER:KMAX_MID) :: vw  ! Svavenging rates (tmp. array)

    real                    :: loss    ! conc. loss due to scavenging
    real,dimension(NWETDEP) :: sumloss ! sum conc. loss due to scavenging




!  Loop starting from above

    f_rho  = xmd(i,j)*(ps(i,j,1) - PT)/ATWAIR

    do k=kcloudtop, KMAX_MID    ! No need to go above cloudtop
        rho(k)  = f_rho * carea(k)/amk(k)
    end do


    sumloss(:) = 0.0

    invgridarea = xm2(i,j)/( gridwidth_m*gridwidth_m )


  ! calculate concentration after wet deposition and sum up the vertical
  ! column of the depositions for the fully soluble species.

    if ( DEBUG_AQ .and. i == 3 .and. j == 3 ) then
        print "(a15,i3,2i4,es14.4)", "DEBUG_WDEP2", me, &
                      kcloudtop, ksubcloud, pr_acc(KMAX_MID)
    end if ! DEBUG

    do spec = 1, NWETDEP

        ! Put both in- and sub-cloud scavenging ratios in the array vw:

         vw(kcloudtop:ksubcloud-1) = WetDep(spec)%W_sca ! Scav. for incloud
         vw(ksubcloud:KMAX_MID  )  = WetDep(spec)%W_sub ! Scav. for sub-cloud

         if ( DEBUG_AQ .and. i == 3 .and. j == 3 ) then
             print "(a15,i3,2es14.4)", "DEBUG_WDEP Sca", &
                   spec, WetDep(spec)%W_sca, WetDep(spec)%W_sub
         end if ! DEBUG

         adv = WetDep(spec)%adv

         do k = kcloudtop, KMAX_MID
              loss =  xn_2d(adv,k) * ( 1.0 - exp( -vw(k)*pr_acc(k)*dt ) )
              xn_2d(adv,k) = xn_2d(adv,k) - loss
              sumloss(spec) = sumloss(spec) + loss * rho(k)

             if ( DEBUG_AQ .and. i == 3 .and. j == 3 ) then
                print "(a30,4i4)", "DEBUG_WDEP me, k, adv, spec", me, k, adv, spec
                print "(a30,i4, 2es12.4, f8.2)", "DEBUG_WDEP me, vw, pr, dt ", me,  &
                       vw(k), pr_acc(k), dt
                print "(a30,3es12.4)", "DEBUG_WDEP loss, xn, sumloss", &
                         loss, xn_2d(adv,k), sumloss(spec)
             end if ! DEBUG
         end do   ! k loop
     end do  ! spec loop



    !hf add up precipitation release, omit divt, ds- sum 1:KMAX_MID
    !ds-query .... isn't this the same as pr_acc(KMAX_MID) ???
    !ds-answer - no, since we had negative pr_acc values !!! 

      wdep(WDEP_PREC,i,j,IOU_INST) = sum ( pr(i,j,:) ) * dt   ! Same for all models


    !/.. add other losses into twetdep and wdep arrays:

      call WetDep_Budget(i,j,sumloss,invgridarea) ! Model-specific


  end subroutine WetDeposition

!------------------------------------------------------------------------------
end module Aqueous_ml





