
module DryDep_ml

  ! Module contains relevant parts of old readpar_mach and drydep.f
  ! MACHO method. This is the drag-coefficient based approach of
  ! BJ98:
  ! Berge, E. and  Jakobsen, H.A., A regional scale multi-layer model
  ! for the calculation of long-term transport and deposition of air
  ! pollution in Europe, Tellus B (1998), 50, 105-223.

 !-- model specific dry dep values are set in My_UKDep_ml:
 !   (in file My_DryDep_ml)
 !ds jan-03 STO_FLUX stuff re-introduced
 !hf dec-02 Structure changed so that DryDep can go inside loop in Runchem

 use My_UKDep_ml, only : Init_DepMap,  &   ! Maps indices between Vg-calculated (CDEP..
                                       &   !and advected  (IXADV_..)
                          NDRYDEP_CALC, &  ! No. Vd values calculated 
                          NDRYDEP_ADV, &   ! No. advected species affected
                          STO_FLUXES,  &   ! true if fluxes wanted.  !got
                          CDEP_SET,    &   ! for so4
                          CDEP_NO2,    &   ! for NO2 comp pt. approach
                          FLUX_CDEP,   &   ! index O3 in CALC array, for STO_FLUXES
                          FLUX_ADV ,   &   ! index O3 in ADV  array, for STO_FLUXES
                          DepLoss, Add_ddep, &
                          Dep        ! Mapping (type = depmap)


 use My_Derived_ml                ! -> d_2d, IOU_INST, D2_VG etc...

 use Dates_ml,       only : daynumber
 use DepVariables_ml,only : NLANDUSE,  & !ds jan2003 LU_WATER, &
                            forest, water, g_pot, g_temp,g_vpd, &
                            g_light,g_swp, &
                            unit_flux,   &! = sto. flux per m2
                            lai_flux      ! = lai * unit_flux
 use Chemfields_ml , only : cfac!,xn_adv
 use GenSpec_adv_ml, only : NSPEC_ADV, IXADV_NO2, IXADV_SO2, IXADV_NH3

 use GridValues_ml , only : GRIDWIDTH_M,xmd,xm2,carea, gb, &
                            i_glob, j_glob   !u1 for testing
 use MassBudget_ml,  only : totddep,DryDep_Budget !hf
 use Met_ml,         only : roa,fm,fh,z_bnd,th2m,th,ps,u,v,z_mid&
                            ,snow, pr, psurf, cc3dmax, t2, q    & ! u7.lu
                            ,iclass   & ! ds rv1.2
                            ,ustar,foundustar, fl
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
 use Setup_1dfields_ml,    only : xn_2d,amk
 use GenSpec_shl_ml,        only :  NSPEC_SHL

 implicit none
 private

 public :: drydep, init_drydep

!hf ddep
   real, private, save, dimension(KUPPER:KMAX_MID) :: &
         pr_acc                 ! Accumulated precipitation
     logical, public, dimension(NDRYDEP_ADV), save :: vg_set 

  logical, private, save :: my_first_call = .true.
  logical, private, parameter :: DEBUG_VG = .false.
  logical, private, parameter :: DEBUG_UK = .false.
  logical, private, parameter :: DEBUG_WET = .false.
  logical, private, parameter :: DEBUG_FLUX = .false.
  logical, private, parameter :: DEBUG_NO2 = .false.


 contains

!hf  subroutine drydep

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

!hf  ustar_sq = -99.99                           ! for debugging
!hf  imm      =    current_date%month            ! for debugging
!hf  idd      =    current_date%day              ! for debugging
!hf  ihh      =    current_date%hour             ! for debugging

  end subroutine init_drydep

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine drydep(i,j)

 real, dimension(NDRYDEP_CALC) :: &
          gradient_fac & ! Ratio of conc. at zref (ca. 50m) and 1m
         ,vg_fac       & ! Loss factor due to dry dep.
         ,Rb           & ! Quasi-boundary layer rsis.
         ,Rsur         & ! 
         ,Rsur_wet     & ! rv1.4.8 - Rsur over wet surface
         ,Vg_ref       & ! Vg at ref ht.
         ,Vg_1m        & ! Vg at  1m over veg.
         ,Grid_Vg_ref  & ! Grid average of Vg at ref ht. (effective Vg for cell)
         ,Grid_Vg_1m   & ! Grid average Vg at  1m over veg.
         ,Vg_ratio     & ! Ratio Vg_ref/Vg_1m = ratio C(1m)/C(ref), over land
         ,sea_ratio      ! Ratio Vg_ref/Vg_1m = ratio C(1m)/C(ref), over sea

!hf logical, dimension(NDRYDEP_ADV), save :: vg_set 
 integer, intent(in):: i,j
 integer n, ilu, lu, nlu, ncalc, nadv, ind, err,k   ! help indexes
 integer :: imm, idd, ihh     ! date
 integer :: nadv2d !index of adv species in xn_2d array
 real ustar_sq, & ! u star squared
      abshd,    & ! abs value of surfae flux of sens. heat W/m^2
      u_ref       ! horizontal vind   !ds - was uvhs

 real :: no2fac  ! Reduces Vg for NO2 in ration (NO2-4ppb)/NO2

 real convfac,  & ! rescaling to different units
      convfac2, & ! rescaling to different units
      convfaco3,& ! rescaling to different units
      dtz         ! scaling factor for veff ( = dt/z, where dt=timestep and 
                  ! z = height of layer)

  real, save :: inv_gridarea  ! inverse of grid area, m2
!hf  integer, save ::  old_daynumber = -99

! Landuse, UKDEP
   real :: cover , Sumcover, Sumland   ! Land-coverage
   real :: z0, g_sto   &
         ,Ra_ref       & ! Ra from ref ht.  to z0
         ,Ra_1m        & ! Ra from 1m over veg. to z0
         ,Idrctt, Idfuse   ! Direct-total and diffuse radiation
   real :: wetarea         ! Fraction of grid square assumed wet 
   real :: ustar_nwp, ustar_loc, vpd, invL, rho_surf, invL_nwp, d, rh, Ts_C
   real :: lai, hveg        ! For convenience  
   real ::  Hd   ! Heat flux, into surface (opp. sign to fh!)
   real ::  LE   ! Latent Heat flux, into surface (opp. sign to fh!)
   logical :: debug_flag   ! set true when i,j match DEBUG_i, DEBUG_j
!   real, parameter  :: NMOLE_M3 = 1.0e6*1.0e9/AVOG   ! Converts from mol/cm3 to nmole/m3
   real :: nmole_o3    ! O3 in nmole/m3
   real :: loss,sto_frac,Vg_scale
   real :: lat_factor   ! latitide-based correction for lai, hveg

 ! Ecosystem specific deposition requires the fraction of dep in each landuse, lu:

   real, dimension(NDRYDEP_CALC,NLUMAX):: Vg_ref_lu
   real, dimension(NSPEC_ADV ,NLANDUSE):: fluxfrac_adv
   integer :: lu_used(NLUMAX), nlu_used
   real    :: lu_cover(NLUMAX)
   real    :: lu_lai(NLUMAX)  ! TMP for debug

!     first calculate the 1m deposition velocity following the same
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
!hf
! Need pr_acc for wet surfaces
! Add up the precipitation in the column
! Note that this is accumulated prec per second
! but this is the same as used in Aqueues 
! as threshold for precipitating clouds
! ds -  in the longer term we might save pr_acc in Met_ml, or even the
!  original "pr" as read in from HIRLAM, but this can wait until we
!  settle on a final scheme for RextS in Rsurface_ml.

    pr_acc(KUPPER) = sum ( pr(i,j,1:KUPPER) )   ! prec inc. from above 
    do k= KUPPER+1, KMAX_MID
      pr_acc(k) = pr_acc(k-1) + pr(i,j,k)
      pr_acc(k) = max( pr_acc(k), 0.0 ) !u7.2 ds - FIX
    end do


     ! - Set up debugging coordinates first. ---------------------------!
     ! If location matches debug i,j value, set debug_flag. Also passed
     ! to Rsurface_ml

      debug_flag = .false. 
      if ( i_glob(i)==DEBUG_i .and. j_glob(j)==DEBUG_j) debug_flag = .true.

      ind = iclass(i,j)   ! nb 0 = sea, 1= ice, 2=tundra
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


      dtz      = dt_advec/z_bnd(i,j,KMAX_BND-1)

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
    Grid_Vg_1m(:) = 0.0
    Vg_ref_lu(:,:) = 0.0
    Vg_ratio(:) = 0.0
    Sumcover = 0.0
    Sumland  = 0.0
    fluxfrac_adv (:,:) = 0.0


    if ( STO_FLUXES ) then
       unit_flux(:) = 0.0
       !ds lai_flux(:) = 0.0
    end if


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

        !ds if ( DEBUG_UK .and. lu == LU_WATER .and. hveg > 0.0 ) then
        if ( DEBUG_UK .and. water(lu) .and. hveg > 0.0 ) then
            print *, "HIGHW!!! h,lai,cov", i,j,ilu, hveg, lai,cover
            call gc_abort(me,NPROC,"WATER!")
        end if


        if ( forest(lu) ) then  !! More realistic forests, test:


           ! if ( gb(i,j) < 50.0 .and. gb(i,j) > 42.5  ) then
           !    lai = 1.5 * lai
           ! else if ( gb(i,j) >= 62.0 ) then

           ! Just used reduced LAI for high latitudes for now, because of tests
           ! which suggest that the big-leaf model as coded will overestimate
           ! Gsto if we allow higher LAI in central Europe.

            if ( gb(i,j) >= 60.0 ) then
               lat_factor  = max(0.3, ( 1.0 - 0.05* (gb(i,j)-60.0)) )
               lai   = lai  *  lat_factor
               hveg  = hveg *  lat_factor
            end if
        end if ! forest



        lu_lai(ilu) = lai  ! tmp for debug
        ustar_loc = ustar_nwp       ! First guess = NWP value
        invL      = invL_nwp        ! First guess = NWP value

        call Get_Submet(hveg,t2(i,j),Hd,LE, psurf(i,j), &
               z_mid(i,j,KMAX_MID), u_ref, q(i,j,KMAX_MID,1), & ! qw_ref
                       debug_flag,                        &    ! in
                       ustar_loc, invL,                   &    ! in-out
                       z0,d, Ra_ref,Ra_1m,rh,vpd)                    ! out

             if ( DEBUG_UK .and. debug_flag ) then
                write(6,"(a40,4i3,f6.1,2i4,3f7.3,2i4,2f6.2)") &
                    "DEBUG_veg: me,nlu,ilu,lu, lat, SGS, EGS ", &
                    me,nlu,ilu, lu, gb(i,j), &
                   landuse_SGS(i,j,ilu), landuse_EGS(i,j,ilu), &
                   cover, lai, hveg,daynumber, snow(i,j), SWP(daynumber),Ts_C

                write(6,"(a10,2i4,2f7.2,2es12.3,3f8.3)") "UKDEP SUB", me, &
                  lu, ustar_nwp, ustar_loc, invL_nwp, invL, Ra_ref, Ra_1m,rh

             end if


         call Rsurface(rh,lu,debug_flag, &
                   lai, &   ! tmp
                   hveg, &  ! tmp
                   z0,ustar_loc,       &
                   Ts_C,               & ! Ts_C
                   vpd,                & ! CHECK 
                   SWP(daynumber),     & ! NEEDS i,j later
                   psurf(i,j),         &
                   pr_acc(KMAX_MID),   &  ! Bug fixed by hf - was pr
                   coszen(i,j),        &  ! CHECK   
                   Idfuse,             &  ! CHECK   
                   Idrctt,             &  ! CHECK   
                   snow(i,j),          &
                   xn_2d(NSPEC_SHL+IXADV_SO2,KMAX_MID) / & !SO2/NH3
                   xn_2d(NSPEC_SHL+IXADV_NH3,KMAX_MID), &
                   g_sto,  &
                   Rsur, &
                   Rsur_wet, &
                   Rb)

       !ds rv1.4.8 - we now allow for areas which are treated
       !   as wet and use Rsur_wet.
       !
       !   We assume that the area of grid which is wet is
       !   proportional to cloud-cover:

         if ( pr_acc(KMAX_MID) > 0.0 ) then
             wetarea = cc3dmax(i,j,KMAX_MID) ! cloud cover above surface
         else
             wetarea = 0.0
         end if


       !/... add to grid-average Vg:


         do n = 1, NDRYDEP_CALC

           Vg_ref(n) = (1.0-wetarea) / ( Ra_ref + Rb(n) + Rsur(n) ) &
                     +      wetarea  / ( Ra_ref + Rb(n) + Rsur_wet(n) )

           Vg_1m (n) = (1.0-wetarea) / ( Ra_1m + Rb(n) + Rsur(n) ) &
                     +      wetarea  / ( Ra_1m + Rb(n) + Rsur_wet(n) )

         ! NO2 compensation point approach:        

         no2fac = xn_2d(NSPEC_SHL+IXADV_NO2,KMAX_MID)   ! Here we have no2 in cm-3
         no2fac = max(1.0, no2fac)
         !   print *, "NO2FAC pre", xn_2d(NSPEC_SHL+IXADV_NO2 ,KMAX_MID)
         no2fac = max(0.00001,  (no2fac-1.0e11)/no2fac)      ! Comp. point of 4 ppb

         if (DEBUG_NO2 .and. debug_flag .and.  &
                ( lu == 1 .or. lu == 10 ) .and. & 
                   (current_date%seconds == 0)  ) then
            if (lu==1) then 
            write(6,"(a10,3i3,i5,2f12.5,2f8.3)") "CONIF-NO2", &
                  imm, idd, ihh, lu, &
                     xn_2d(NSPEC_SHL+IXADV_NO2 ,KMAX_MID)*4.0e-11, no2fac,&
                      100.0*Vg_ref(CDEP_NO2), 100.0*Vg_ref(CDEP_NO2)*no2fac
            else
            write(6,"(a10,3i3,i5,2f12.5,2f8.3)") "GRASS-NO2", &
                  imm, idd, ihh, lu, &
                     xn_2d(NSPEC_SHL+IXADV_NO2 ,KMAX_MID)*4.0e-11, no2fac,&
                      100.0*Vg_ref(CDEP_NO2), 100.0*Vg_ref(CDEP_NO2)*no2fac
            end if
         end if

         Vg_ref(CDEP_NO2) = Vg_ref(CDEP_NO2) * no2fac
         Vg_1m (CDEP_NO2) = Vg_1m (CDEP_NO2) * no2fac

           Vg_ref_lu(n,ilu) = Vg_ref(n)
           Grid_Vg_ref(n) = Grid_Vg_ref(n) + cover * Vg_ref(n)
           Grid_Vg_1m(n)  = Grid_Vg_1m(n)  + cover * Vg_1m(n)

           if ( DEBUG_VG .and.  Vg_ref(n) < 0.0 .or. Vg_1m(n)<Vg_ref(n) ) then
               print *, "VGREF ERROR", me, n, Vg_ref(n), " Ras ", &
                  Ra_ref, Rb(n), Rsur(n)
               print *, "VGREF ERROR stab", z0, d, ustar_loc, invL, g_sto
               call gc_abort(me, NPROC, "VGREF ERROR")
           end if

         end do


         Sumcover = Sumcover + cover


        !/-- only grab gradients over land-areas
        !ds - Using water() is better than LU_WATER

         if ( water(lu) ) then
            do n = 1, NDRYDEP_CALC
               sea_ratio(n) =  Vg_ref(n)/Vg_1m(n)
            end do
         else
            Sumland = Sumland + cover
            do n = 1, NDRYDEP_CALC
                Vg_ratio(n) =  Vg_ratio(n) + cover * Vg_ref(n)/Vg_1m(n)
            end do
         end if


        if ( DEBUG_UK .and. debug_flag ) then
            print "(a14,2i4,f8.3,5f7.2,f12.3)", "UKDEP RATVG",ilu,lu,cover,&
                100.0*Vg_1m(4), 100.0*Vg_ref(4), Sumland, Sumcover, &
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

        if ( STO_FLUXES ) then
          unit_flux(lu) = g_sto * Rsur(FLUX_CDEP)
          lai_flux(lu)  = lai * unit_flux(lu)
        end if        

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
                 (100.0*Grid_Vg_1m(n), n = 1, min(5,NDRYDEP_CALC))
        end if


!-- loop through all affected advected species to calculate changes in
!   concentration (xn_adv), the conc. ratios (cfac), and deposition 

    do ncalc = 1, NDRYDEP_CALC

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

          if ( STO_FLUXES .and. nadv == FLUX_ADV ) then

             loss  = Vg_scale &
                    * DepLoss(FLUX_ADV)/AVOG*1.e6  &  ! From mol/cm3 to mole/m3
                    * z_bnd(i,j,KMAX_BND-1) * 1.e9     ! mole/m3 -> nmole/m2

             ! dt_advec or dt_drydep ????? I guess ddep takes care of this
             ! for accumulated?

             unit_flux(lu) = unit_flux(lu) * loss

             sto_frac      = lai_flux(lu)       ! eqv. to  LAI*gsto/Gsur
             lai_flux(lu)  = lai_flux(lu)  * loss


             if ( DEBUG_FLUX .and. debug_flag) then
                 write (6,"(a5,i4,3i4,f6.1,f8.4,es10.2,f7.2,3f8.4,2es10.2)") &
                   "OSTAD", lu, imm, idd, ihh, &
                        lu_lai(ilu), &
                        Vg_scale, &
                      DepLoss(FLUX_ADV), & 
                      z_bnd(i,j,KMAX_BND-1), &
                       sto_frac, &
                         4.0e-11*xn_2d(nadv2d,KMAX_MID),& ! O3 in ppb
                        cfac(nadv,i,j), &
                        unit_flux(lu)/dt_advec, lai_flux(lu)/dt_advec
             end if
           end if

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
       convfaco3 = convfac / (amk(KMAX_MID)*dtz)


      !.. Add DepLoss to budgets if needed:

       !if ( DEBUG_VG .and. debug_flag ) then
       !  write(6,"(a7,5i5,3i3,i5,2f8.4,3f10.4)") "FLUX ",  me, i, j, nlu, &
       !      lu_used(1), imm, idd, ihh, current_date%seconds, &
       !      Sumcover, Sumland, fluxfrac_adv(7,15), fluxfrac_adv(8,15), &
       !      fluxfrac_adv(9,15)
       !end if 

       call Add_ddep(i,j,convfac2,convfaco3,fluxfrac_adv)

!hf     enddo   !j
!hf   enddo    !i
 end subroutine drydep

end module DryDep_ml
