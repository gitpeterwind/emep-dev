
module DryDep_ml

  ! Module contains relevant parts of old readpar_mach and drydep.f
  ! MACHO method. This is the drag-coefficient based approach of
  ! BJ98:
  ! Berge, E. and  Jakobsen, H.A., A regional scale multi-layer model
  ! for the calculation of long-term transport and deposition of air
  ! pollution in Europe, Tellus B (1998), 50, 105-223.

 !-- model specific dry dep values are set in My_UKDep_ml:
 !   (in file My_DryDep_ml)
 !st mar-03 Particle stuff introduced
 !ds jan-03 STO_FLUX stuff re-introduced
 !hf dec-02 Structure changed so that DryDep can go inside loop in Runchem

 use My_UKDep_ml, only : Init_DepMap,  &   ! Maps indices between Vg-calculated (CDEP..
                                       !&   !and advected  (IXADV_..)
                          NDRYDEP_CALC, &  ! No. Vd values calculated 
                          NDRYDEP_ADV, &   ! No. advected species affected
                          NDRYDEP_AER, &   !st No. Vd values for aerosol size modes
                          NDRYDEP_TOT, &   !st Total No. of  Vd values
                          STO_FLUXES,  &   ! true if fluxes wanted.
                          CDEP_SET,    &   ! for so4
                          CDEP_NO2,    &   ! for NO2 comp pt. approach
                          FLUX_CDEP,   &   ! index O3 in CALC array, for STO_FLUXES
                          FLUX_ADV ,   &   ! index O3 in ADV  array, for STO_FLUXES
                          DepLoss, Add_ddep, &
                          Dep        ! Mapping (type = depmap)


 use My_Derived_ml      ! -> OUTPUT_ABS_HEIGHTS, d_2d, IOU_INST, D2_VG etc...

 use Dates_ml,       only : daynumber
 use DepVariables_ml,only : NLANDUSE,  & !ds jan2003 LU_WATER, &
                            forest, water, g_pot, g_temp,g_vpd, &
                            crops,         & !rv1_7_4 for SAIadd
                            luflux_wanted, & ! true if fluxes wanted for landuse lu
                            g_light,g_swp, &
                            STUBBLE,       & ! Small ht., 1 cm
                            WHEAT,         & ! Used for Wheat Fst calculations, rv1_9_7
                            SAIadd,        &
                            SLAIlen,       & !rv1_7_4 for SAIadd
                            leaf_flux,   &! = flag-leaf sto. flux per m2
                            unit_flux,   &! = sto. flux per m2
                            lai_flux      ! = lai * unit_flux
 use Chemfields_ml , only : cfac!,xn_adv
 use GenSpec_adv_ml, only : NSPEC_ADV, IXADV_NO2, IXADV_SO2, IXADV_NH3

 use Functions_ml,   only : AerRes,PsiM
 use GridValues_ml , only : GRIDWIDTH_M,xmd,xm2,carea, gb, &
                            i_glob, j_glob   ! for testing
 use MassBudget_ml,  only : totddep,DryDep_Budget !hf
 use Met_ml,         only : roa,fm,fh,z_bnd,th2m,th,ps,u,v,z_mid&
                            ,snow, pr, psurf, cc3dmax, t2, q &
                            ,surface_precip & ! ds rv1.6.2 
                            ,iclass   &
                            ,ustar,foundustar, fl &
                            ,pzpbl  !stDep
 use ModelConstants_ml,    only : dt_advec,PT,KMAX_MID, KMAX_BND ,&
                                  current_date,     &  ! u7.lu
                                  DEBUG_i, DEBUG_j, &
                                  ATWAIR, atwS, atwN, PPBINV,&
                                  KUPPER     !hf ddep
 use Par_ml,               only : me,NPROC,li0,li1,lj0,lj1
 use PhysicalConstants_ml, only : XKAP, PI, KARMAN, GRAV, RGAS_KG, CP, AVOG
 use Radiation_ml,         only : zen         &! zenith angle (degrees)
                                 ,SolBio      &!u7.lu extras
                                 ,coszen
 
 use SubMet_ml,        only: Get_Submet
 use UKdep_ml,         only : Init_ukdep, ReadLanduse, SetLandUse  & 
                              ,NLUMAX &  ! Max. no countries per grid
                              ,landuse_ncodes, landuse_codes, landuse_data  &
                              ,landuse_SGS, landuse_EGS &
                              ,landuse_LAI,    landuse_hveg , landuse_gpot

 use Rsurface_ml
 use SoilWater_ml, only : SWP ! = 0.0 always for now!
 use Wesely_ml,    only : Init_GasCoeff !  Wesely stuff, DRx, Rb_Cor, ...
 use Setup_1dfields_ml,    only : xn_2d,amk,Idrctt,Idfuse
 use GenSpec_shl_ml,        only :  NSPEC_SHL
!stDep
 use Aero_DryDep_ml,        only : Aero_Rb
 use My_Aerosols_ml,        only : NSIZE

 implicit none
 private

 public :: drydep, init_drydep

  logical, public, dimension(NDRYDEP_ADV), save :: vg_set 

  logical, private, save :: my_first_call = .true.
  logical, private, parameter :: DEBUG_VG = .false.
  logical, private, parameter :: DEBUG_UK = .false.
  logical, private, parameter :: DEBUG_WET = .false.
  logical, private, parameter :: DEBUG_FLUX = .false.
  logical, private, parameter :: DEBUG_NO2 = .false.
  logical, private, parameter :: DEBUG_AERO = .false.


 contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     subroutine init_drydep


     integer, save ::  old_daynumber = -99
     integer ::nadv,n

  if ( my_first_call ) then 

     call Init_DepMap()               ! Maps CDEP to IXADV


     call Init_ukdep()                ! reads ukdep_biomass, etc.
     call Init_GasCoeff()             ! Sets Wesely coeffs.
     call ReadLanduse()

     nadv = 0
     do n = 1, NDRYDEP_ADV  
         nadv       = max( Dep(n)%adv, nadv )  ! Looking for highest IXADV
         vg_set(n)  = ( Dep(n)%calc == CDEP_SET ) ! for set vg
         if ( DEBUG_UK .and. me == 0 ) print *, "VGSET ", n, nadv, vg_set(n)
     end do

     my_first_call = .false.

  end if !  my_first_call

  if ( old_daynumber /= daynumber ) then

      call SetLandUse()         ! Sets landuse_LAI, landuse_hveg 
      old_daynumber = daynumber

  end if

  end subroutine init_drydep

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine drydep(i,j)

 real, dimension(NDRYDEP_CALC) :: &
          Rb           & ! Quasi-boundary layer rsis.
         ,Rsur         & ! 
         ,Rsur_wet       ! rv1.4.8 - Rsur over wet surface
!stDry
 real, dimension(NDRYDEP_TOT) :: &
          gradient_fac & ! Ratio of conc. at zref (ca. 50m) and 3m
         ,vg_fac       & ! Loss factor due to dry dep.
         ,Vg_ref       & ! Vg at ref ht.
         ,Vg_3m        & ! Vg at  3m
         ,Grid_Vg_ref  & ! Grid average of Vg at ref ht. (effective Vg for cell)
         ,Grid_Vg_3m   & ! Grid average Vg at  3m (or tree height)
         ,Vg_ratio     & ! Ratio Vg_ref/Vg_3m = ratio C(3m)/C(ref), over land
         ,sea_ratio      ! Ratio Vg_ref/Vg_3m = ratio C(3m)/C(ref), over sea

 integer, intent(in):: i,j
 integer n, ilu, lu, nlu, ncalc, nadv, ind, err,k   ! help indexes
 integer :: imm, idd, ihh     ! date
 integer :: nadv2d !index of adv species in xn_2d array
 real ustar_sq, & ! u star squared
      abshd,    & ! abs value of surfae flux of sens. heat W/m^2
      u_ref,    & ! horizontal vind   !ds - was uvhs
      z_ref       !  reference level - middle of cell, ca. 45m

 real :: no2fac  ! Reduces Vg for NO2 in ration (NO2-4ppb)/NO2
 real :: RaVs    ! Ra_ref *Vs for particles
 real, dimension(NOUTPUT_ABS_HEIGHTS) :: htz, Raz  ! Heights and Ra(z)

 real convfac,  & ! rescaling to different units
      convfac2, & ! rescaling to different units
      lossfrac,  & !  If needed in My_DryDep - not used now.
      !ds convfaco3,& ! rescaling to different units
      dtz         ! scaling factor for veff ( = dt/z, where dt=timestep and 
                  ! z = height of layer)
 logical :: is_wet  ! set true if surface_precip>0
 logical :: hirlam_sea   !ds .. set true if iclass == 0

!stDep
 integer :: nae
 real, dimension(NSIZE):: aeRb, aeRbw , Vs
 real :: wstar, convec   


  real, save :: inv_gridarea  ! inverse of grid area, m2

! Landuse, UKDEP
   real :: cover , Sumcover, Sumland   ! Land-coverage
   real :: z0, g_sto   &
         ,Ra_ref       & ! Ra from ref ht.  to z0
         ,Ra_3m          ! Ra from 3m over veg. to z0
         !ds ,Idrctt, Idfuse   ! Direct-total and diffuse radiation
   real :: wetarea         ! Fraction of grid square assumed wet 
   real :: ustar_nwp, ustar_loc, vpd, invL, rho_surf, invL_nwp, d, rh, Ts_C
   real :: lai, hveg        ! For convenience  
   real ::  Hd   ! Heat flux, into surface (opp. sign to fh!)
   real ::  LE   ! Latent Heat flux, into surface (opp. sign to fh!)
   logical :: debug_flag   ! set true when i,j match DEBUG_i, DEBUG_j
!ICP: re-instated:
   real, parameter  :: NMOLE_M3 = 1.0e6*1.0e9/AVOG   ! Converts from mol/cm3 to nmole/m3
   real :: nmole_o3, ppb_o3    ! O3 in nmole/m3, ppb
   real :: loss,sto_frac,Vg_scale
   real :: gsun        ! g_sto for sunlit upper-canopy (flag) leaves
   real :: lat_factor   ! latitide-based correction for lai, hveg
   real :: so2nh3ratio  ! So2/NH3 ratio, for Rsur calc

 ! Ecosystem specific deposition requires the fraction of dep in each landuse, lu:

!stDep   real, dimension(NDRYDEP_CALC,NLUMAX):: Vg_ref_lu
   real, dimension(NDRYDEP_TOT,NLUMAX):: Vg_ref_lu
   real, dimension(NSPEC_ADV ,NLANDUSE):: fluxfrac_adv
   integer :: lu_used(NLUMAX), nlu_used
   real    :: lu_cover(NLUMAX)
   real    :: lu_lai(NLUMAX)  ! TMP for debug
!ICP:
   real :: c_hveg, u_hveg ! Values at canopy top, for fluxes
   real    :: c_hvegppb(NLANDUSE)
   real ::  Ra_diff       ! Ra from z_ref to hveg
   real :: gext_leaf, rc_leaf, rb_leaf, Fst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!   real, save, dimension ::: leafw
!if( lu ==  9 ) g_pot = 0.8  !!! TFMM FOR CLe wheat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     first calculate the 3m deposition velocity following the same
!     procedure as in the lagmod. second the flux and the accumulated 
!     deposition is calculated.
!
!     effective dry deposition velocity applied to the model concentration
!     at the top of the constant flux layer, zdep 
!     Dry deposion rates are specified in subroutine readpar
!

!hf  do j = lj0,lj1
!hf    do i = li0,li1

!hf FOR DEBUGGING s
  ustar_sq = -99.99                           ! for debugging
  imm      =    current_date%month            ! for debugging
  idd      =    current_date%day              ! for debugging
  ihh      =    current_date%hour             ! for debugging

     inv_gridarea = 1.0/(GRIDWIDTH_M*GRIDWIDTH_M) 

  !ds   We assume that the area of grid which is wet is
  !     proportional to cloud-cover:

     is_wet = ( surface_precip(i,j) > 0.0 )

     if ( is_wet ) then
         wetarea = cc3dmax(i,j,KMAX_MID) ! cloud cover above surface
     else
         wetarea = 0.0
     end if


     ! - Set up debugging coordinates first. ---------------------------!
     ! If location matches debug i,j value, set debug_flag. Also passed
     ! to Rsurface_ml

      debug_flag = .false. 
      if ( i_glob(i)==DEBUG_i .and. j_glob(j)==DEBUG_j) debug_flag = .true.


      ind = iclass(i,j)   ! nb 0 = sea, 1= ice, 2=tundra

      hirlam_sea  = ( iclass(i,j) == 0 )   ! nb 0 = sea, 1= ice, 2=tundra

     ! -----------------------------------------------------------------!
     !.and conversion facrtor,  convfac (( ps-pt)/grav... )  ===> 
     !      pressure in kg m-1 s-2

      convfac = (ps(i,j,1) - PT)*carea(KMAX_MID)*xmd(i,j)/ATWAIR
     ! -----------------------------------------------------------------!

!JEJ  to be implemented
!
!     !.and factor,  kg_air_ij (( ps-pt)/grav... )  ===> 
!     !      pressure in kg m-1 s
!     !      used for converting from mixing ratio to kg
!
!      kg_air_ij = (ps(i,j,1) - PT)*carea(KMAX_MID)
!     ! -----------------------------------------------------------------!
!
      lossfrac = 1.0 !  Ratio of xn before and after deposition

     !ds - I prefer micromet terminology here:

      Hd = -fh(i,j,1)       ! Heat flux, *away from* surface
      LE = -fl(i,j,1)       ! Heat flux, *away from* surface


      if(foundustar)then                      !pw u3
         ustar_sq = ustar(i,j,1)*ustar(i,j,1)
      else
         ustar_sq = fm(i,j,1)/roa(i,j,KMAX_MID,1)
         ustar(i,j,1) = sqrt(ustar_sq)
         ustar(i,j,1) = max( ustar(i,j,1), 1.0e-5 )
      endif
      ustar_sq = max(ustar_sq, 1.e-10)

     ! wind-speed at reference height, which we take as centre of
     ! lowest layer:

      u_ref = 0.5*sqrt( (u(i,j,KMAX_MID,1)+u(i-1,j,KMAX_MID,1))**2  &
                     + (v(i,j,KMAX_MID,1)+v(i,j-1,KMAX_MID,1))**2 )

      z_ref = z_mid(i,j,KMAX_MID)

      dtz      = dt_advec/z_bnd(i,j,KMAX_BND-1)

     ! wstar for particle deposition
        wstar = 0.                   ! w* ......Based on Wesely
        if(Hd <  0.0 )   &          ! for unstable stratification
        wstar = (-GRAV * pzpbl(i,j) * Hd /      &
        (roa(i,j,KMAX_MID,1)*CP*th(i,j,KMAX_MID,1))) ** (1./3.)


    call SolBio(daynumber, coszen(i,j), &
               cc3dmax(i,j,KMAX_MID)  & ! cloud cover above surface
               ,psurf(i,j)            &
               ,Idfuse, Idrctt)   ! output radiation
      

    !   we must use L (the Monin-Obukhov length) to calculate deposition,
    !   therefore we calculate u*, t* from NWP-model data. 

      rho_surf       = psurf(i,j)/(RGAS_KG * t2(i,j) )   ! at surface

      ustar_nwp = ustar(i,j,1)    ! First guess = NWP value
      !tstar_nwp = -fh(i,j,1)/ ( CP*rho_surf*ustar_nwp )
      invL_nwp  = -KARMAN * GRAV * Hd/ &
         (CP*rho_surf* ustar_nwp*ustar_nwp*ustar_nwp * th2m(i,j,1))

      !.. we limit the range of 1/L to prevent numerical and printout problems
      !.. and because I don't trust HIRLAM enough.
      !   This range is very wide anyway.

      invL_nwp  = max( -1.0, invL_nwp ) !! limit very unstable
      invL_nwp  = min(  1.0, invL_nwp ) !! limit very stable

    if ( DEBUG_UK .and. debug_flag ) then
          print "(a26,4i4)", "UKDEP DryDep me, i,j, ind ", me, i,j, ind
          print "(a10,i4,3i3,i6,f6.1,f8.3,4f7.2)", "UKDEP SOL", &
                  daynumber, imm, idd, ihh, current_date%seconds, &
                   zen(i,j), coszen(i,j), cc3dmax(i,j,KMAX_MID), &
                     1.0e-5*psurf(i,j), Idfuse, Idrctt
          print "(a10,i4,3i3,2f8.3,es12.4,f8.4)", "UKDEP NWP", &
                  daynumber, imm, idd, ihh,  &
                  Hd, LE, invL_nwp, ustar_nwp
    end if
                 

 
    !/ Initialise Grid-avg Vg for this grid square:

    Grid_Vg_ref(:) = 0.0
    Grid_Vg_3m(:) = 0.0
    Vg_ref_lu(:,:) = 0.0
    Vg_ratio(:) = 0.0
    Sumcover = 0.0
    Sumland  = 0.0
    fluxfrac_adv (:,:) = 0.0
    c_hvegppb(:) = 0.0   !ds rv1_9_16 needed for safety
    leaf_flux(:) = 0.0   !ds rv1_9_17 bug-fix - 26/12/2003

    !/ SO2/NH3 for Rsur calc
    so2nh3ratio = &
               xn_2d(NSPEC_SHL+IXADV_SO2,KMAX_MID) / & !SO2/NH3
               max(1.0,xn_2d(NSPEC_SHL+IXADV_NH3,KMAX_MID))



    if ( STO_FLUXES ) then
       !dsDEC unit_flux(:) = 0.0
       !dsDEC leaf_flux(:) = 0.0

    ! --- ICP -----
       !ppb_o3   =  xn_2d(NSPEC_SHL+IXADV_O3,KMAX_MID) * 4.0e-11  !Crude for now
       !nmole_o3 =  xn_2d(NSPEC_SHL+IXADV_O3,KMAX_MID) * NMOLE_M3
       !ds ppb_o3   =  xn_2d(NSPEC_SHL+FLUX_ADV,KMAX_MID) * 4.0e-11  !Crude for now

    ! xn_2d is in #/cm3. x2_2d/amk gives mixing ratio, and 1.0e9*xn_2d/amk = ppb

       ppb_o3   =  1.0e9* xn_2d(NSPEC_SHL+FLUX_ADV,KMAX_MID)/amk(KMAX_MID)
       nmole_o3 =  xn_2d(NSPEC_SHL+FLUX_ADV,KMAX_MID) * NMOLE_M3

    end if ! STO_FLUXES


    !/ And start the sub-grid stuff over different landuse (lu)

    nlu = landuse_ncodes(i,j)
    LULOOP: do ilu= 1, nlu
        lu      = landuse_codes(i,j,ilu)
        cover   = landuse_data (i,j,ilu)

        lu_used (ilu) = lu    ! for eco dep
        lu_cover(ilu) = cover

        Ts_C    = t2(i,j)-273.15
        lai     = landuse_LAI(i,j,ilu)
        hveg    = landuse_hveg(i,j,ilu)
        g_pot   = landuse_gpot(i,j,ilu)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if( lu ==  WHEAT ) g_pot = 1.0  !!!ds New Mapping Manual IAM suggestion, March 2004
        !ds rv1_9_28 if( lu ==  WHEAT ) g_pot = 0.8  !!! TFMM FOR CLe wheat
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if ( DEBUG_UK .and. water(lu) .and. hveg > 0.0 ) then
            print *, "HIGHW!!! h,lai,cov", i,j,ilu, hveg, lai,cover
            call gc_abort(me,NPROC,"WATER!")
        end if


        if ( forest(lu) ) then  !! More realistic forests, test:



           ! Just used reduced LAI for high latitudes for now, because of tests
           ! which suggest that the big-leaf model as coded will overestimate
           ! Gsto if we allow higher LAI in central Europe.

            if ( gb(i,j) >= 60.0 ) then
               lat_factor  = max(0.3, ( 1.0 - 0.05* (gb(i,j)-60.0)) )
               lai   = lai  *  lat_factor
               hveg  = hveg *  lat_factor
            end if
        end if ! forest

      !ds code moved here from UK_ml to fix bug spotted by Peter, 19/8/2003
      !SAIadd for other vegetation still defined in UK_ml
       if (  crops(lu) ) then

            if ( daynumber < landuse_SGS(i,j,ilu) .or. &
                   daynumber > landuse_EGS(i,j,ilu)  ) then
                     SAIadd(lu) = 0.0
                else if ( daynumber < &
                          (landuse_SGS(i,j,ilu) + SLAIlen(lu)) ) then
                     SAIadd(lu) = ( 5.0/3.5 - 1.0) * lai
            !pw else if ( daynumber < landuse_EGS(i,j,lu) ) then
!rv2_0_4pw      else if ( daynumber <= landuse_EGS(i,j,lu) ) then
                else if ( daynumber <= landuse_EGS(i,j,ilu) ) then
                     SAIadd(lu) = 1.5   ! Sensescent
            end if
        end if ! crops




        lu_lai(ilu) = lai  ! tmp for debug
        ustar_loc = ustar_nwp       ! First guess = NWP value
        invL      = invL_nwp        ! First guess = NWP value

        !ds - rv1.6.xx coast
        ! If HIRLAM thinks this is a sea-square, but we anyway have land,
        ! the surface temps will be wrong and so will stability gradients.
        ! As a simple subsistute, we assume neutral condiutions for these
        ! situations.

        if( .not. water(lu) .and. hirlam_sea  ) then
             invL = 0.0
             Hd = 0.0
        end if

        call Get_Submet(hveg,t2(i,j),Hd,LE, psurf(i,j), &
               z_ref, u_ref, q(i,j,KMAX_MID,1), & ! qw_ref
                       debug_flag,                        &    ! in
                       ustar_loc, invL,                   &    ! in-out
                       z0,d, Ra_ref,Ra_3m,rh,vpd)                    ! out

             if ( DEBUG_UK .and. debug_flag ) then
                write(6,"(a40,4i3,f6.1,2i4,3f7.3,2i4,2f6.2)") &
                    "DEBUG_veg: me,nlu,ilu,lu, lat, SGS, EGS ", &
                    me,nlu,ilu, lu, gb(i,j), &
                   landuse_SGS(i,j,ilu), landuse_EGS(i,j,ilu), &
                   cover, lai, hveg,daynumber, snow(i,j), SWP(daynumber),Ts_C

                write(6,"(a10,2i4,2f7.2,2es12.3,3f8.3)") "UKDEP SUB", me, &
                  lu, ustar_nwp, ustar_loc, invL_nwp, invL, Ra_ref, Ra_3m,rh

             end if


         call Rsurface(rh,lu,debug_flag, &
                   lai, &   ! tmp
                   hveg, &  ! tmp
                   z0,ustar_loc,       &
                   Ts_C,               & ! Ts_C
                   vpd,                & ! CHECK 
                   SWP(daynumber),     & ! NEEDS i,j later
                   psurf(i,j),         &
                   is_wet,             &  !ds rv1.6.2 true if wet
                   !ds pr_acc(KMAX_MID),   &  ! Bug fixed by hf - was pr
                   coszen(i,j),        &  ! CHECK   
                   Idfuse,             &  ! CHECK   
                   Idrctt,             &  ! CHECK   
                   snow(i,j),          &
                   so2nh3ratio,        & !SO2/NH3
                   g_sto,  &
                   gsun,   &        ! sun-leaf g_sto, TFMM
                   Rsur, &
                   Rsur_wet, &
                   Rb)


       !===================
       !stDep
       !// calculate dry deposition velocities for fine/coarse particles

         convec = wstar*wstar/(ustar_loc*ustar_loc)     ! Convection velocity scale  

        call Aero_Rb ( ustar_loc, convec, roa(i,j,KMAX_MID,1)     &
                     , u_ref, lu, snow(i,j), wetarea, t2(i,j)      &   
                     , Vs, aeRb, aeRbw )
       !===================

         if( DEBUG_AERO .and. debug_flag) then
             write(6,*)  '  -- Gravitational settling --'
             write(6,*) (Vs(n), n=1,NDRYDEP_AER)
             write(6,*)  '  --Rb  --'
             write(6,*) (aeRb(n), n=1,NDRYDEP_AER )
         endif



       !/... add to grid-average Vg:


         do n = 1, NDRYDEP_TOT  !stDep  NDRYDEP_CALC

            !stDep
            if ( n > NDRYDEP_CALC)  then    ! particles

                nae = n - NDRYDEP_CALC
                RaVs = Ra_ref * Vs(nae)  ! Mainly to keep lines <78 chars, ds!

                Vg_ref(n) = Vs(nae) +      &
                  (1.0-wetarea) / (Ra_ref + aeRb(nae)  + RaVs  *aeRb(nae)  ) &
                +      wetarea  / (Ra_ref + aeRbw(nae) + RaVs  *aeRbw(nae) )
                    
                RaVs = Ra_3m  * Vs(nae)  ! Mainly to keep the equations short, ds!

                Vg_3m(n)  = Vs(nae) +       &
                  (1.0-wetarea) / (Ra_3m + aeRb(nae)  + RaVs  *aeRb(nae)  ) &
                +      wetarea  / (Ra_3m + aeRbw(nae) + RaVs  *aeRbw(nae) )

            else                           ! gases
               Vg_ref(n) = (1.0-wetarea) / ( Ra_ref + Rb(n) + Rsur(n) ) &
                     +      wetarea  / ( Ra_ref + Rb(n) + Rsur_wet(n) )

               Vg_3m (n) = (1.0-wetarea) / ( Ra_3m + Rb(n) + Rsur(n) ) &
                     +      wetarea  / ( Ra_3m + Rb(n) + Rsur_wet(n) )

            endif

         ! NO2 compensation point approach:        

            if ( n == CDEP_NO2 ) then

            no2fac = xn_2d(NSPEC_SHL+IXADV_NO2,KMAX_MID)   ! Here we have no2 in cm-3
            no2fac = max(1.0, no2fac)
            no2fac = max(0.00001,  (no2fac-1.0e11)/no2fac) ! Comp. point of 4 ppb

            Vg_ref(CDEP_NO2) = Vg_ref(CDEP_NO2) * no2fac
            Vg_3m (CDEP_NO2) = Vg_3m (CDEP_NO2) * no2fac
          end if ! CDEP_NO2

           Vg_ref_lu(n,ilu) = Vg_ref(n)
           Grid_Vg_ref(n) = Grid_Vg_ref(n) + cover * Vg_ref(n)
           Grid_Vg_3m(n)  = Grid_Vg_3m(n)  + cover * Vg_3m(n)

         end do


         Sumcover = Sumcover + cover


        !/-- only grab gradients over land-areas
        !ds - Using water() is better than LU_WATER

         if ( water(lu) ) then
            do n = 1, NDRYDEP_TOT  !stDep NDRYDEP_CALC
               sea_ratio(n) =  Vg_ref(n)/Vg_3m(n)
            end do
         else
            Sumland = Sumland + cover
            do n = 1, NDRYDEP_TOT  !stDep  NDRYDEP_CALC
                Vg_ratio(n) =  Vg_ratio(n) + cover * Vg_ref(n)/Vg_3m(n)
            end do
         end if

        if( DEBUG_AERO .and. debug_flag ) then
           write(6,*) 
           write(6,*) ' >=>=>=>=>=>=>  Dry deposition velocity at :', &
              i_glob(i), j_glob(j)
           write(6,'(7(i3,f8.3))') &
              (n, 100.*Grid_Vg_3m(n), n = 1,NDRYDEP_TOT)
           write(6,'(7(i3,f8.3))') &
              (n, 100.*Grid_Vg_ref(n), n = 1,NDRYDEP_TOT)
           write(6,*)' >=>=>=>=>=>=>>=>=>=>=>=>=>>=>=>=>=>=>=>=>=>=>=>=>=>'
        endif

        if ( DEBUG_UK .and. debug_flag ) then
            print "(a14,2i4,f8.3,5f7.2,f12.3)", "UKDEP RATVG",ilu,lu,cover,&
                100.0*Vg_3m(4), 100.0*Vg_ref(4), Sumland, Sumcover, &
                   Vg_ratio(4)   !helcom, 4=NH3
            print "(a14,2i4,f8.3,5f7.2,f12.3)", "UKDEP FINVG", ilu, lu, cover, &
                   lai,100.0*g_sto, &  ! tmp, in cm/s 
                   Ra_ref, Rb(4), Rsur(4), 100.0/(Ra_ref + Rb(4) + Rsur(4) )
        end if
         
       !=======================
       ! The fraction going to the stomata = g_sto/g_sur = g_sto * R_sur.
       ! Vg*nmole_o3 is the total deposition flux of ozone, but
       ! we calculate the actual flux later (once we know DepLoss(O3)).
       ! For now we just calculate the g_sto*R_sur bit:
       ! (Caution - g_sto is for O3 only)

        if ( STO_FLUXES .and. luflux_wanted(lu) ) then
          !dsDEC2003 - following stuff not needed any more. leaf_flux
          !dsDEC2003   used instead
          !dsDEC unit_flux(lu) = g_sto * Rsur(FLUX_CDEP)
          !dsDEC lai_flux(lu)  = lai * unit_flux(lu)

          !ds_sep27
          if ( hveg < 1.1 * z0 ) then   ! needed for temp crops, outside growing season
                leaf_flux(lu) = 0.0
          else 

          ! ICP method for flag-leaf

          Ra_diff = AerRes(max(hveg-d, STUBBLE) ,z_ref-d,ustar_loc,invL,KARMAN)
          c_hveg         = nmole_o3 * ( 1.0 - Ra_diff * Vg_ref(FLUX_CDEP) )
          c_hvegppb(lu)  = ppb_o3   * ( 1.0 - Ra_diff * Vg_ref(FLUX_CDEP) )


          if ( DEBUG_FLUX .and. c_hveg <= 1.0e-19 ) then
            write(6,*) "FLUX ZERO ", lu, nmole_o3, FLUX_CDEP, Vg_ref(FLUX_CDEP)
          end if

          u_hveg  = u_ref *  &
             ( log((hveg-d)/z0)  -PsiM((hveg-d)*invL) + PsiM(z0*invL)  )/ &
             ( log((z_ref-d)/z0) -PsiM((z_ref-d)*invL) + PsiM(z0*invL))

          if ( DEBUG_FLUX .and. u_hveg <= 1.0e-19 ) then
            print *, "ERRROR UHVEG", u_ref, u_hveg, z_ref,d,z0,hveg
            print *, "ERRROR UHVEG DATE ",  imm, idd, ihh
            print *, "ERRROR UHVEG HTS ", lu, hveg, d, z0
            print *, "ERRROR UHVEG TERMS", &
                             invL, &
                             log((hveg-d)/z0), &  
                             PsiM((hveg-d)*invL), &
                             log((z_ref-d)/z0), &
                             PsiM((z_ref-d)*invL), &
                             PsiM(z0*invL) 
            call gc_abort(me,NPROC,"UHVEG!")
          end if

          gext_leaf = 1.0/2500.0
          rc_leaf = 1.0/(g_sto+ gext_leaf)

          !McNaughton + van den Hurk:
          ! 2cm leaf for wheat - too large for trees?

	if ( forest(lu)) then  !BEECH
          rb_leaf = 1.3 * 150.0 * sqrt(0.04/u_hveg)  ! 2cm leaves? 
	else ! WHEAT
          rb_leaf = 1.3 * 150.0 * sqrt(0.02/u_hveg)  ! 2cm leaves? 
        endif

         ! Flux in nmole/m2/s:
          leaf_flux(lu) = c_hveg * rc_leaf/(rb_leaf+rc_leaf) * gsun 

          if ( DEBUG_FLUX .and. STO_FLUXES .and. debug_flag ) then 
            write(6,"(a8,3i3,i4,2f6.2,2f8.1,2es10.2,f6.2,es12.3)") &
                "FST ", lu, imm, idd, ihh, lai, SAIadd(lu), &
                 nmole_o3, c_hveg, g_sto, gsun, u_hveg,leaf_flux(lu)
          end if

          end if !ds_sep27-----------
        end if ! STO_FLUXES

       !=======================
        end do LULOOP
       !=======================
       !=======================

        if ( DEBUG_UK .and. Sumland > 1.011  ) then
            print *, "SUMLAND ", me, nlu, i,j,i_glob(i), j_glob(j), Sumland
            call gc_abort(me,NPROC,"SUMLAND!")
        end if


        if ( Sumland > 0.01 ) then
            gradient_fac(:) = Vg_ratio(:) / Sumland
        else
            gradient_fac(:) = sea_ratio(:)
        end if

        if ( DEBUG_UK .and. debug_flag ) then
            print "(a14,i2,3i3,10f6.2)", "UKDEP VG_UKR", &
                   snow(i,j), imm, idd, ihh,  &
                 (100.0*Grid_Vg_ref(n), n = 1, min(5,NDRYDEP_CALC)), &
                 (100.0*Grid_Vg_3m(n), n = 1, min(5,NDRYDEP_CALC))
        end if


!-- loop through all affected advected species to calculate changes in
!   concentration (xn_adv), the conc. ratios (cfac), and deposition 

    do ncalc = 1, NDRYDEP_TOT  !stDep NDRYDEP_CALC

        vg_fac (ncalc) = 1.0 - exp ( -Grid_Vg_ref(ncalc) * dtz ) 

    end do ! n

      !hf - replaced xn_adv by xn_2d below, since we are now working within
      !     an i,j loop

      do n = 1, NDRYDEP_ADV 
         nadv    = Dep(n)%adv
         nadv2d  = NSPEC_SHL + Dep(n)%adv

         ncalc   = Dep(n)%calc

         if ( vg_set(n) ) then

             DepLoss(nadv) =   & ! Use directly set Vg
                 ( 1.0 - exp ( -Dep(n)%vg * dtz ) ) * xn_2d( nadv2d,KMAX_MID)
             cfac(nadv, i,j) = 1.0   ! Crude, for now.
  
         else
             DepLoss(nadv) =   vg_fac( ncalc )  * xn_2d( nadv2d,KMAX_MID)
             cfac(nadv, i,j) = gradient_fac( ncalc )
         end if

         if ( DepLoss(nadv) < 0.0 .or. &
              DepLoss(nadv)>xn_2d(nadv2d,KMAX_MID) ) then
           call gc_abort(me,NPROC,"NEG XN DEPLOSS")
         end if


        xn_2d( nadv2d,KMAX_MID) = &
             xn_2d( nadv2d,KMAX_MID) - DepLoss(nadv)



        if ( STO_FLUXES .and. nadv == FLUX_ADV ) then
           ! fraction by which xn is reduced - used in
           !ds rv1_9_17: safety measure:
             
              if( xn_2d( nadv2d,KMAX_MID)  > 1.0e-30 ) then
                  lossfrac = ( 1.0 - DepLoss(nadv)/ &
                                (DepLoss(nadv)+xn_2d( nadv2d,KMAX_MID)))
              end if
              if ( DEBUG_FLUX .and. lossfrac < 0.1 ) then
                  print *, "ERROR: LOSSFRAC ", lossfrac, nadv, nadv2d
                  call gc_abort(me,NPROC,"LOSSFRAC!")
              end if
        end if


      !.. ecosystem specific deposition - translate from calc to adv and normalise

         do ilu = 1, nlu
            lu      = lu_used(ilu)    ! faster than using landuse_codes(i,j,ilu)????

            if ( vg_set(n) )  then
               fluxfrac_adv(nadv,lu) = lu_cover(ilu)  ! Since all vg_set equal
            else
               Vg_scale = Vg_ref_lu(ncalc,ilu)/ Grid_Vg_ref(ncalc)
               fluxfrac_adv(nadv,lu) = lu_cover(ilu)*Vg_scale
            end if


          !=======================
          ! The fraction going to the stomata = g_sto/g_sur = g_sto * R_sur.
          ! Vg*nmole_o3 is the instantaneous deposition flux of ozone, and
          ! the actual amount "deposited" is obtained from DepLoss(O3) using
          ! fluxfrac as obtained above.

          ! Now, DepLoss is loss of molecules/cm3 over time-period dt_advec
          ! and depth z_bnd over 1m2 of grid. For sto fluxes we need to
          ! find values over 1m2 of vegeation (regardless of how much veg
          ! is in grid, so we don't need lu_cover. Instead:

          !dsDEC2003 - Most of the following stuff not needed any more. leaf_flux
          !            used instead

          !dsDEC2003 if ( STO_FLUXES .and. nadv == FLUX_ADV ) then
!if(debug_flag) then
!  print "(a15, 3i3,2es12.3,f8.4)","CHVEG ", i,j, lu, c_hveg, c_hvegppb(lu), lossfrac
!end if

             !dsDEC loss  = Vg_scale &
             !dsDEC       * DepLoss(FLUX_ADV)/AVOG*1.e6  &  ! From mol/cm3 to mole/m3
             !dsDEC       * z_bnd(i,j,KMAX_BND-1) * 1.e9     ! mole/m3 -> nmole/m2

             ! dt_advec or dt_drydep ????? I guess ddep takes care of this
             ! for accumulated?

             !dsDEC unit_flux(lu) = unit_flux(lu) * loss

             !dsDEC sto_frac      = lai_flux(lu)       ! eqv. to  LAI*gsto/Gsur
             !dsDEC lai_flux(lu)  = lai_flux(lu)  * loss


             !dsDEC if ( DEBUG_UK   .and. debug_flag) then
             !dsDEC    write (6,"(a5,i4,3i4,f6.1,f8.4,es10.2,3f8.4,3es10.2)") &
             !dsDEC      "OSTAD", lu, imm, idd, ihh, lu_lai(ilu), Vg_scale, &
             !dsDEC         DepLoss(FLUX_ADV), sto_frac, &
             !dsDEC            4.0e-11*xn_2d(nadv2d,KMAX_MID),& ! O3 in ppb
             !dsDEC           cfac(nadv,i,j), &
             !dsDEC           unit_flux(lu)/dt_advec, lai_flux(lu)/dt_advec,&
             !dsDEC           leaf_flux(lu)
             !dsDEC end if
           !dsDEC2003 end if

            if ( DEBUG_UK .and. debug_flag ) then

               if ( vg_set(n) )  then
                 write(6,"(a12,3i3,3f12.3)") "FLUXSET  ", ilu, lu, nadv, &
                     100.0*Dep(n)%vg, lu_cover(ilu), fluxfrac_adv(nadv,lu)
               else
                 write(6,"(a12,3i3,f6.3,4f8.3)") "FLUXFRAC ", ilu, lu, nadv, &
            lu_cover(ilu), &
            100.0*Grid_Vg_ref(ncalc), 100.0*Vg_ref_lu(ncalc,ilu), &
            100.0*lu_cover(ilu)*Vg_ref_lu(ncalc,ilu), fluxfrac_adv(nadv,lu)
               end if
            end if
         end do
             

          

      !..accumulated dry deposition per grid square and summed over the whole
      !  domain

!         totddep( nadv ) = totddep (nadv) + Deploss(nadv) * convfac


        if ( DEBUG_UK .and. debug_flag ) then
          if ( vg_set(n) ) then
              print "(a30,2i4,f8.3)", "DEBUG DryDep SET ",  n,nadv, Dep(n)%vg
          else
              print "(a30,3i4,f12.5)", &
                  "DEBUG DryDep n, adv, calc, fac ",  n,nadv, ncalc, gradient_fac( ncalc)
              print "(a20,2e12.4)", &
                "DEBUG xn, DepLoss ", xn_2d(nadv2d,KMAX_MID), DepLoss(nadv)
              print "(a20,2f8.4)", "DEBUG gv_fac( ncalc)", &
                 vg_fac(ncalc), 1.0-vg_fac(ncalc)
          end if
        end if

       end do ! n

!hf Not needed inside IXADV_ loop
!ds - Could use DEBUG_i, DEBUG_j here also. I added DEBUG_VG
!     since the compiler will ignor tis if-test and hence be faster
!     unless DEBUG_VG is set.

        if (DEBUG_VG .and. i==2 .and. j==2 .and. nadv==49)then
            write(*,*)'nadv, Deploss',49,DepLoss(nadv)
        endif
        call DryDep_Budget(i,j,Deploss,convfac)

       ! inv_gridarea = xm2(i,j)/(GRIDWIDTH_M*GRIDWIDTH_M)
       convfac2 = convfac * xm2(i,j) * inv_gridarea/amk(KMAX_MID)
       !ds convfaco3 = convfac / (amk(KMAX_MID)*dtz)


      !.. Add DepLoss to budgets if needed:

       !ds call Add_ddep(i,j,convfac2,convfaco3,fluxfrac_adv,c_hvegppb)
       call Add_ddep(debug_flag,dt_advec,i,j,convfac2,lossfrac,fluxfrac_adv,c_hvegppb)

!hf     enddo   !j
!hf   enddo    !i
 end subroutine drydep

end module DryDep_ml
