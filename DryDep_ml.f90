
module DryDep_ml

  ! Module contains relevant parts of old readpar_mach and drydep.f
  ! MACHO method. This is the drag-coefficient based approach of
  ! BJ98:
  ! Berge, E. and  Jakobsen, H.A., A regional scale multi-layer model
  ! for the calculation of long-term transport and deposition of air
  ! pollution in Europe, Tellus B (1998), 50, 105-223.

 !-- model specific dry dep values are set in My_UKDep_ml:
 !   (in file My_DryDep_ml)

 use My_UKDep_ml, only : Init_DepMap,  &   ! Maps indices between Vg-calculated (CDEP..
                                       &   !and advected  (IXADV_..)
                          NDRYDEP_CALC, &  ! No. Vd values calculated 
                          NDRYDEP_ADV, &   ! No. advected species affected
                          CDEP_SET,    &   ! for so4
                          FLUX_CDEP,   &   ! for stom. fluxes if wanted
                          DepLoss, Add_ddep, &
                          Dep        ! Mapping (type = depmap)

 use My_Derived_ml                ! -> d_2d, IOU_INST, D2_VG etc...

 use Dates_ml,       only : daynumber !u7.lu
 use DepVariables_ml,only : NLANDUSE,  LU_WATER, &
                            forest, g_pot, g_temp,g_vpd, &
                            g_light,g_swp
 use Chemfields_ml , only : xn_adv,cfac
 use GenSpec_adv_ml, only : NSPEC_ADV
 use GridValues_ml , only : GRIDWIDTH_M,xmd,xm2,carea, gb, &
                            i_glob, j_glob   !u1 for testing
 use MassBudget_ml,  only : totddep
 use Met_ml,         only : roa,fm,fh,z_bnd,th2m,th,ps,u,v,z_mid&
                            ,snow, pr, psurf, cc3dmax, t2, q    & ! u7.lu
                            ,iclass   & ! ds rv1.2
                            ,ustar,foundustar, fl
 use ModelConstants_ml,    only : dt_advec,PT,KMAX_MID, KMAX_BND ,&
                                  current_date,     &  ! u7.lu
                                  DEBUG_i, DEBUG_j, &
                                  ATWAIR, atwS, atwN, PPBINV
 use Par_ml,               only : me,NPROC,li0,li1,lj0,lj1
 use PhysicalConstants_ml, only : XKAP, PI, KARMAN, GRAV, RGAS_KG, CP
 use Radiation_ml,         only : zen         &! zenith angle (degrees)
                                 ,SolBio      &!u7.lu extras
                                 ,coszen
!u7.4vg -
!rv1.2 use My_Derived_ml,    only : d_2d, IOU_INST , &  ! Store results here
!rv1.2                             D2_VG_REF, D2_VG_1M, D2_VG_STO, &
!rv1.2                             D2_FX_REF, D2_FX_STO
 
 use SubMet_ml,        only: Get_Submet
 use UKdep_ml,         only : Init_ukdep, ReadLanduse, SetLandUse  & 
                              ,NLUMAX &  ! Max. no countries per grid
                              ,landuse_ncodes, landuse_codes, landuse_data  &
                              ,landuse_SGS, landuse_EGS &
                              ,landuse_LAI,    landuse_hveg , landuse_gpot
!rv1.2                        ,dep_flux      &
!rv1.2                       ,DEP_VG_REF, DEP_VG_1M, DEP_VG_STO &
!rv1.2                       ,DEP_FL_REF, DEP_FL_STO

 use Rsurface_ml
 use SoilWater_ml, only : SWP ! = 0.0 always for now!
 use Wesely_ml,    only : Init_GasCoeff !  Wesely stuff, DRx, Rb_Cor, ...
 implicit none
 private

 public :: drydep



  logical, private, save :: my_first_call = .true.
  logical, private, parameter :: DEBUG_VG = .false.
  logical, private, parameter :: DEBUG_UK = .false.


 contains
  subroutine drydep

 real, dimension(NDRYDEP_CALC) :: &
          gradient_fac & ! Ratio of conc. at zref (ca. 50m) and 1m
         ,vg_fac       & ! Loss factor due to dry dep.
         ,Rb           & ! Quasi-boundary layer rsis.
         ,Rsur         & ! 
         ,Vg_ref       & ! Vg at ref ht.
         ,Vg_1m        & ! Vg at  1m over veg.
         ,Grid_Vg_ref  & ! Grid average of Vg at ref ht. (effective Vg for cell)
         ,Grid_Vg_1m   & ! Grid average Vg at  1m over veg.
         ,Vg_ratio     & ! Ratio Vg_ref/Vg_1m = ratio C(1m)/C(ref), over land
         ,sea_ratio      ! Ratio Vg_ref/Vg_1m = ratio C(1m)/C(ref), over sea

 logical, dimension(NDRYDEP_ADV), save :: vg_set 

 integer i, j, n, ilu, lu, nlu, ncalc, nadv, ind, err   ! help indexes
 integer :: imm, idd, ihh     ! date

 real ustar_sq, & ! u star squared
      abshd,    & ! abs value of surfae flux of sens. heat W/m^2
      u_ref       ! horizontal vind   !ds - was uvhs

 real convfac,  & ! rescaling to different units
      convfac2, & ! rescaling to different units
      dtz         ! scaling factor for veff ( = dt/z, where dt=timestep and 
                  ! z = height of layer)

  real, save :: inv_gridarea  ! inverse of grid area, m2
  integer, save ::  old_daynumber = -99

! Landuse, UKDEP
   real :: cover , Sumcover, Sumland   ! Land-coverage
   real :: z0, g_sto   &
         ,Ra_ref       & ! Ra from ref ht.  to z0
         ,Ra_1m        & ! Ra from 1m over veg. to z0
         ,Idrctt, Idfuse   ! Direct-total and diffuse radiation
   real :: ustar_nwp, ustar_loc, vpd, invL, rho_surf, invL_nwp, d, rh, Ts_C
   real :: lai, hveg        ! For convenience  
   real ::  Hd   ! Heat flux, into surface (opp. sign to fh!)
   real ::  LE   ! Latent Heat flux, into surface (opp. sign to fh!)
   logical :: debug_flag   ! set true when i,j match DEBUG_i, DEBUG_j
   real :: nmole_o3    ! O3 in nmole/m3
   real :: lat_factor   ! latitide-based correction for lai, hveg

 ! Ecosystem specific deposition requires the fraction of dep in each landuse, lu:

   real, dimension(NDRYDEP_CALC,NLANDUSE):: fluxfrac_calc
   real, dimension(NSPEC_ADV ,NLANDUSE):: fluxfrac_adv
   integer :: lu_used(NLUMAX), nlu_used
   real    :: lu_cover(NLUMAX)



!     first calculate the 1m deposition velocity following the same
!     procedure as in the lagmod. second the flux and the accumulated 
!     deposition is calculated.
!
!     effective dry deposition velocity applied to the model concentration
!     at the top of the constant flux layer, zdep 
!     Dry deposion rates are specified in subroutine readpar
!

  if ( my_first_call ) then 

     call Init_DepMap()                          ! Maps CDEP to IXADV

     inv_gridarea = 1.0/(GRIDWIDTH_M*GRIDWIDTH_M) 

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

  ustar_sq = -99.99                           ! for debugging
  imm      =    current_date%month            ! for debugging
  idd      =    current_date%day              ! for debugging
  ihh      =    current_date%hour             ! for debugging

  do j = lj0,lj1
    do i = li0,li1

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


      !BUG??  dtz      = dt_advec/z_mid(i,j,KMAX_MID)

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
    Vg_ratio(:) = 0.0
    Sumcover = 0.0
    Sumland  = 0.0
    !! fluxfrac_calc(:,:) = 0.0
    fluxfrac_adv (:,:) = 0.0

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

        if ( DEBUG_UK .and. lu == LU_WATER .and. hveg > 0.0 ) then
            print *, "HIGHW!!! h,lai,cov", i,j,ilu, hveg, lai,cover
            call gc_abort(me,NPROC,"WATER!")
        end if

        !rv1.2 if ( lu  <= 4 ) then  !! More realistic forests, test:
        !rv1.2                       !! Hard-coded 4  remove later!!

        if ( forest(lu) ) then  !! More realistic forests, test:
                              !! Hard-coded 4  remove later!!

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


         call Rsurface(lu,debug_flag, &
                   lai, &   ! tmp
                   hveg, &  ! tmp
                   z0,ustar_loc,       &
                   Ts_C,               & ! Ts_C
                   vpd,                & ! CHECK 
                   SWP(daynumber),     & ! NEEDS i,j later
                   psurf(i,j),         &
                   pr(i,j,KMAX_MID),   &  ! CHECK   
                   coszen(i,j),        &  ! CHECK   
                   Idfuse,             &  ! CHECK   
                   Idrctt,             &  ! CHECK   
                   snow(i,j),          &
                   g_sto,  &
                   Rsur, &
                   Rb)

         !/... add to grid-average Vg:


         do n = 1, NDRYDEP_CALC
           Vg_ref(n) = 1.0 / ( Ra_ref + Rb(n) + Rsur(n) )
           Vg_1m(n)  = 1.0 / ( Ra_1m  + Rb(n) + Rsur(n) )
           ! fluxfrac_calc(n,lu) = cover * Vg_ref(n)
           Grid_Vg_ref(n) = Grid_Vg_ref(n) + cover * Vg_ref(n)
           Grid_Vg_1m(n)  = Grid_Vg_1m(n)  + cover * Vg_1m(n)

           if ( Vg_ref(n) < 0.0 .or. Vg_ref(n) > 0.1 .or. Vg_1m(n)<Vg_ref(n) ) then
               print *, "VGREF ERROR", me, n, Vg_ref(n), " Ras ", Ra_ref, Rb(n), Rsur(n)
               print *, "VGREF ERROR stab", z0, d, ustar_loc, invL, g_sto
               call gc_abort(me, NPROC, "VGREF ERROR")
           end if

         end do

         if ( FLUX_CDEP > 0 .and. Vg_ref(1) > 0.04 ) then !CRUDE FIX WITH (1)
               print *, "VGREFO3 ERROR", me, n, Vg_ref(n)
               print *, "VGREFO3 ERROR Ras", Ra_ref, Rb(n), Rsur(n)
               print *, "VGREFO3 ERROR stab", z0, d, ustar_loc, invL, g_sto
               call gc_abort(me, NPROC, "VGREFO3 ERROR")
         end if

         Sumcover = Sumcover + cover


        !/-- only grab gradients over land-areas

         if ( lu /= LU_WATER ) then
            Sumland = Sumland + cover
            do n = 1, NDRYDEP_CALC 
                Vg_ratio(n) =  Vg_ratio(n) + cover * Vg_ref(n)/Vg_1m(n)
            end do
         else
            do n = 1, NDRYDEP_CALC 
               sea_ratio(n) =  Vg_ref(n)/Vg_1m(n)
            end do
         end if

        if ( DEBUG_UK .and. debug_flag ) then
            print "(a14,2i4,f8.3,5f7.2,f12.3)", "UKDEP RATVG", ilu, lu, cover, &
                Vg_1m(2), Vg_ref(2), Sumland, Sumcover, Vg_ratio(2)
            print "(a14,2i4,f8.3,5f7.2,f12.3)", "UKDEP FINVG", ilu, lu, cover, &
                   lai,100.0*g_sto, &  ! tmp, in cm/s 
                   Ra_ref, Rb(2), Rsur(2), 100.0/(Ra_ref + Rb(2) + Rsur(2) )
        end if
         
       !=======================
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

      do n = 1, NDRYDEP_ADV 
         nadv    = Dep(n)%adv
         ncalc   = Dep(n)%calc

         if ( vg_set(n) ) then

             DepLoss(nadv) =   & ! Use directly set Vg
                 ( 1.0 - exp ( -Dep(n)%vg * dtz ) ) * xn_adv( nadv,i,j,KMAX_MID)
             cfac(nadv, i,j) = 1.0   ! Crude, for now.
  
         else
             DepLoss(nadv) =   vg_fac( ncalc )  * xn_adv( nadv,i,j,KMAX_MID)
             cfac(nadv, i,j) = gradient_fac( ncalc )
         end if

         if ( DepLoss(nadv) < 0.0 .or. &
              DepLoss(nadv)>xn_adv(nadv,i,j,KMAX_MID) ) then
           print *,"NEG XN DEPLOSS!!! ", DepLoss(nadv), xn_adv(nadv,i,j,KMAX_MID)
           call gc_abort(me,NPROC,"NEG XN DEPLOSS")
         end if

         xn_adv( nadv,i,j,KMAX_MID) = &
             xn_adv( nadv,i,j,KMAX_MID) - DepLoss(nadv)

      !.. ecosystem specific deposition - translate from calc to adv and normalise

         do ilu = 1, nlu
            lu      = lu_used(ilu)    ! faster than using landuse_codes(i,j,ilu)????

            if ( vg_set(n) )  then
               fluxfrac_adv(nadv,lu) = lu_cover(ilu)  ! Since all vg_set equal
            else
               fluxfrac_adv(nadv,lu) = lu_cover(ilu)*Vg_ref(ncalc)/ Grid_Vg_ref(ncalc)
            end if

            if ( DEBUG_UK .and. debug_flag ) then

               if ( vg_set(n) )  then
                 write(6,"(a12,3i3,3f12.3)") "FLUXSET  ", ilu, lu, nadv, &
                           Dep(n)%vg, lu_cover(ilu), fluxfrac_adv(nadv,lu)
               else
                 write(6,"(a12,3i3,3f12.3)") "FLUXFRAC ", ilu, lu, nadv, &
                  Grid_Vg_ref(ncalc), lu_cover(ilu)*Vg_ref(ncalc), fluxfrac_adv(nadv,lu)
               end if
            end if
         end do
             

          

      !..accumulated dry deposition per grid square and summed over the whole
      !  domain

         totddep( nadv ) = totddep (nadv) + Deploss(nadv) * convfac

        if ( DEBUG_UK .and. debug_flag ) then
          if ( vg_set(n) ) then
              print "(a30,2i4,f8.3)", "DEBUG DryDep SET ",  n,nadv, Dep(n)%vg
          else
              print "(a30,3i4,f12.5)", &
                  "DEBUG DryDep n, adv, calc, fac ",  n,nadv, ncalc, gradient_fac( ncalc)
              print "(a20,2e12.4)", &
                "DEBUG xn, DepLoss ", xn_adv(nadv, i,j,KMAX_MID), DepLoss(nadv)
              print "(a20,2f8.4)", "DEBUG gv_fac( ncalc)", &
                 vg_fac(ncalc), 1.0-vg_fac(ncalc)
          end if
        end if

       end do ! n

       ! inv_gridarea = xm2(i,j)/(GRIDWIDTH_M*GRIDWIDTH_M)
       convfac2 = convfac * xm2(i,j) * inv_gridarea


      !.. Add DepLoss to budgets if needed:

       call Add_ddep(i,j,convfac2,fluxfrac_adv)

     enddo   !j
   enddo    !i
 end subroutine drydep

end module DryDep_ml
