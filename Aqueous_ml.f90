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
  use Met_ml,               only :  pr, cw, roa, z_bnd, cc3d, ps
  use Par_ml              , only : me   ! for DEBUG
  use Setup_1dfields_ml,    only : xn_2d, amk
  implicit none
  private

  !/-- subroutines

  public :: Setup_Clouds   ! characterises clouds and calls WetDeposition if rain

  public :: WetDeposition   !u7.2, ds -  simplified setup_wetdep


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
      ,CW_LIMIT = 1.0e-7      ! for cloud water
!u7.2 , B_LIMIT = 1.0e-3      ! for cloud cover (fraction)

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
      !u7.2 cloudwater(:) = 0.



  !  Loop starting at surface finding the cloud base

    ksubcloud = KMAX_MID+1    ! k-coordinate of sub-cloud limit

    do  k = KMAX_MID, KUPPER, -1  
        if ( cw(i,j,k,1) >  CW_LIMIT ) exit  !  out of loop
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

        !u7.2 b(k) = cc3d(i,j,k)  ! ds simplify tests
        !u7.2 if ( cw(i,j,k,1) > CW_LIMIT .and. b(k) > B_LIMIT ) then
        !  Units: kg(w)/kg(air) * kg(air(m^3) / density of water 10^3 kg/m^3
        !  ==>  cloudwater (volume mixing ratio of water to air in cloud 
        !  (when devided bu cloud fraction b )
        !u7.2   cloudwater(k) = 1.0e-3 * cw(i,j,k,1) * roa(i,j,k,1) / b(k)   

        if ( cw(i,j,k,1) > CW_LIMIT ) then

                incloud(k) = .true.
                if ( kcloudtop < 0 ) kcloudtop   = k
        end if

    end do
    
    if ( prclouds_present .and. kcloudtop == -1 ) then
           if ( DEBUG_AQ ) write(6,"(a20,i3,2i4,3es12.4)") &
               "ERROR prclouds sum_cw", &
                me, i,j, maxval(cw(i,j,KUPPER:KMAX_MID,1) ) , &
                 maxval(pr(i,j,:)), pr_acc(KMAX_MID)
           kcloudtop = KUPPER ! for safety
    end if

!  sets up the aqueous phase reaction rates (SO2 oxidation) and the 
!  fractional solubility 
!u7.2    call setup_aqurates(b ,cloudwater,incloud)



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
!u7.2    do k = KUPPER, KMAX_MID
!u7.2       if ( cloudwater(k) .gt.1.e-8 ) then ! in cloud scavenging only for 
                                           ! cloud water above 10^-2 g m^-3 
!u7.2       end if
!u7.2    end do
!u7.2 roak1 = xmd(i,j)*(ps(i,j,1) - PT)*carea(k)/amk(k)/ATWAIR


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
              loss =   xn_2d(adv,k) * ( 1.0 - exp( -vw(k)*pr_acc(k)*dt ) )
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


