
module DryDep_ml

  ! Module contains relevant parts of old readpar_mach and drydep.f
  ! MACHO method. This is the drag-coefficient based approach of
  ! BJ98:
  ! Berge, E. and  Jakobsen, H.A., A regional scale multi-layer model
  ! for the calculation of long-term transport and deposition of air
  ! pollution in Europe, Tellus B (1998), 50, 105-223.

 !-- model specific dry dep values are set in My_UKDep_ml:
 !   (in file My_DryDep_ml)

 use My_UKDep_ml, only : Init_vd, & ! Initialisation
                          NDRYDEP_CALC, &  ! No. Vd values calculated 
                          NDRYDEP_ADV, &   ! No. advected species affected
                          CDEP_SET,    &   ! for so4
                          CDEP_O3,     &   ! for fluxes - hardcoded :-(
                          DepLoss, Add_ddep, &
                          Dep        ! Mapping (type = depmap)

 use Dates_ml,       only : daynumber !u7.lu
 use Chemfields_ml , only : xn_adv,cfac
 use GenSpec_adv_ml, only : NSPEC_ADV  &
  ,IXADV_O3    ! TMP FIX UKDEP ONLY
 use GridValues_ml , only : GRIDWIDTH_M,xmd,xm2,carea,iclass, gb, &
                            i_glob, j_glob   !u1 for testing
 use MassBudget_ml,  only : totddep
 use Met_ml,         only : roa,fm,fh,z_bnd,th2m,th,ps,u,v,z_mid&
                            ,snow, pr, psurf, cc3dmax, t2, q    & ! u7.lu
                            ,ustar,foundustar, fl
 use ModelConstants_ml,    only : dt_advec,PT,KMAX_MID, KMAX_BND ,&
                                  !u7.lu  daynumber,  &  !u3 for seasonal  Vg
                                  current_date, &  ! u7.lu
                                  ATWAIR, atwS, atwN, PPBINV
 use Par_ml,               only : me,NPROC,li0,li1,lj0,lj1
 use PhysicalConstants_ml, only : XKAP, PI, KARMAN, GRAV, RGAS_KG, CP
 use Radiation_ml,         only : zen         &! zenith angle (degrees)
                                 ,SolBio      &!u7.lu extras
                                 ,coszen
!u7.4vg -
 use My_Derived_ml,    only : d_2d, IOU_INST , &  ! Store results here
                             D2_VG_REF, D2_VG_1M, D2_VG_STO, &
                             D2_FX_REF, D2_FX_STO
 
 use DepVariables_ml,  only: LU_WATER, &
             g_pot, g_temp,g_vpd,g_light,g_swp  !u7.4 for possible outputs
 use SubMet_ml,        only: Get_Submet
 use UKdep_ml,         only : Init_ukdep, ReadLanduse, SetLandUse  & 
                              ,landuse_ncodes, landuse_codes, landuse_data  &
                              ,landuse_SGS, landuse_EGS &
                              ,landuse_LAI,    landuse_hveg &
                              ,dep_flux      &
                             ,DEP_VG_REF, DEP_VG_1M, DEP_VG_STO &
                             ,DEP_FL_REF, DEP_FL_STO

 use Rsurface_ml
 use SoilWater_ml, only : SWP ! = 0.0 always for now!
 use Wesely_ml,    only : Init_GasCoeff !  Wesely stuff, DRx, Rb_Cor, ...
 implicit none
 private

 public :: drydep



  logical, private, save :: my_first_call = .true.
  logical, private, parameter :: DEBUG_VG = .false.
  logical, private, parameter :: DEBUG_UK = .false.
 !integer, private, parameter :: debug_i=79, debug_j=56 ! Eskdalemuir
 !integer, private, parameter :: debug_i=97, debug_j=62 !  Waldhof
 !integer, private, parameter :: debug_i=73, debug_j=48 ! Mace Head
 !integer, private, parameter :: debug_i=91, debug_j=71 ! Rorvik
 !integer, private, parameter :: debug_i=82, debug_j=72 !  Voss, has some snow
 integer, private, parameter :: debug_i=101, debug_j=51 !  Schauinsland


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

 integer i, j, n, ilu, lu, nlu, ncalc, nadv, ind   ! help indexes

 real ustar_sq, & ! u star squared
      abshd,    & ! abs value of surfae flux of sens. heat W/m^2
      u_ref       ! horizontal vind   !ds - was uvhs

 real convfac,  & ! rescaling to different units
      convfac2, & ! rescaling to different units
      dtz         ! scaling factor for veff ( = dt/z, where dt=timestep and 
                  ! z = height of layer)

  real :: sin2, cos2, seasonfac   !u3 For seasonal variation

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
   logical :: debug_flag   ! set true when i,j match debug_i, debug_j
   real :: nmole_o3    ! O3 in nmole/m3
   real :: lat_factor   ! latitide-based correction for lai, hveg



!     first calculate the 1m deposition velocity following the same
!     procedure as in the lagmod. second the flux and the accumulated 
!     deposition is calculated.
!
!     effective dry deposition velocity applied to the model concentration
!     at the top of the constant flux layer, zdep 
!     Dry deposion rates are specified in subroutine readpar
!

  if ( my_first_call ) then 

     call Init_vd()
     inv_gridarea = 1.0/(GRIDWIDTH_M*GRIDWIDTH_M) 

     call Init_ukdep()
     call Init_GasCoeff()
     call ReadLanduse(debug_i,debug_j)

     do n = 1, NDRYDEP_ADV  
         nadv    = Dep(n)%adv
         vg_set(n)  = ( Dep(n)%calc == CDEP_SET ) ! for set vg
         if ( DEBUG_UK .and. me == 0 ) print *, "VGSET ", n, nadv, vg_set(n)
     end do

     my_first_call = .false.

  end if !  my_first_call

  if ( old_daynumber /= daynumber ) then

      call SetLandUse()         ! Sets landuse_LAI, landuse_hveg 
      old_daynumber = daynumber

  end if

!1) return
  ustar_sq = -99.99   ! for debugging

  do j = lj0,lj1
    do i = li0,li1

     ! - Set up debugging coordinates first. ---------------------------!
     ! If location matches debug i,j value, set debug_flag. Also passed
     ! to Rsurface_ml

      debug_flag = .false. 
      if ( i_glob(i)==debug_i .and. j_glob(j)==debug_j) debug_flag = .true.

      ind = iclass(i,j)   ! nb 0 = sea, 1= ice, 2=tundra
      if ( DEBUG_UK .and. debug_flag ) then
          print "(a26,4i4)", "DEBUG DryDep me, i,j, ind ", me, i,j, ind
      end if
     ! -----------------------------------------------------------------!



     !ds - I prefer micromet terminology here:

      Hd = -fh(i,j,1)       ! Heat flux, *into* surface
      LE = -fl(i,j,1)       ! Heat flux, *into* surface


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


      convfac = (ps(i,j,1) - PT)*carea(KMAX_MID)*xmd(i,j)/ATWAIR

      !BUG??  dtz      = dt_advec/z_mid(i,j,KMAX_MID)

      dtz      = dt_advec/z_bnd(i,j,KMAX_BND-1)

    call SolBio(daynumber, coszen(i,j), &
               cc3dmax(i,j,KMAX_MID)  & ! cloud cover above surface
               ,psurf(i,j)            &
               ,Idfuse, Idrctt)   ! output radiation
      
    if ( DEBUG_UK .and. debug_flag ) then
             print "(a10,i4,3i3,f6.1,f8.3,4f7.2)", "UKDEP SOL", &
                  daynumber,&
                  current_date%month, &
                  current_date%day, &
                  current_date%hour, &
                 zen(i,j), coszen(i,j), &
                 cc3dmax(i,j,KMAX_MID), 1.0e-5*psurf(i,j), Idfuse, Idrctt
             print "(a10,i4,3i3,i6,2f8.3)", "UKDEP Hs ",  &
                  daynumber,&
                  current_date%month, &
                  current_date%day, &
                  current_date%hour,  current_date%seconds, Hd, LE
    end if

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
             print "(a10,i4,3i3,2f8.3,es12.4,f8.4)", "UKDEP NWP", &
                  daynumber,&
                  current_date%month, &
                  current_date%day, &
                  current_date%hour, &
                  Hd, LE, invL_nwp, ustar_nwp
    end if
                 

 
    !/ Initialise Grid-avg Vg for this grid square:

    Grid_Vg_ref(:) = 0.0
    Grid_Vg_1m(:) = 0.0
    Vg_ratio(:) = 0.0
    Sumcover = 0.0
    Sumland  = 0.0

    !/ And start the sub-grid stuff over different landuse (lu)

    nlu = landuse_ncodes(i,j)
    LULOOP: do ilu= 1, nlu
        lu      = landuse_codes(i,j,ilu)

        cover   = landuse_data (i,j,ilu)
        Ts_C    = t2(i,j)-273.15
        lai     = landuse_LAI(i,j,ilu)
        hveg    = landuse_hveg(i,j,ilu)
        if ( lu == 15 .and. hveg > 0.0 ) then
            print *, "HIGHW!!! is", i,j,ilu
            print *, "HIGHW!!! h,lai,cov",  hveg, lai,cover
            call gc_abort(me,NPROC,"WATER!")
        end if

        if ( lu  <= 4 ) then  !! More realistic forests, test:
                              !! Hard-coded 4  remove later!!

           ! if ( gb(i,j) < 50.0 .and. gb(i,j) > 42.5  ) then
           !    lai = 1.5 * lai
           ! else if ( gb(i,j) >= 62.0 ) then

           ! Just used reduced LAI for high latitudes for now, because of tests
           ! which suggest that the big-leaf model as coded will overestimate
           ! Gsto if we allow higher LAI in central Europe.

           !report if ( gb(i,j) >= 62.0 ) then
           !report    lai  = 0.5 * lai
           !report    hveg = 0.5 * hveg
           !report end if
            if ( gb(i,j) >= 60.0 ) then
               lat_factor  = max(0.3, ( 1.0 - 0.05* (gb(i,j)-60.0)) )
               lai   = lai  *  lat_factor
               hveg  = hveg *  lat_factor
            end if
        end if


        if ( DEBUG_VG )then
           if ( lu <= 0 ) call gc_abort(me,NPROC,"LU ERORR")
           if (  DEBUG_UK .and. debug_flag ) then
                write(6,"(a40,4i3,f6.1,2i4)") &
                    "DEBUG_veg: me,nlu,ilu,lu, lat, SGS, EGS ", &
                          me,nlu,ilu, lu, gb(i,j), &
                         landuse_SGS(i,j,ilu), landuse_EGS(i,j,ilu)
                write(6,"(a44,3f7.3,2i4,2f6.2)") &
                  "DEBUG_veg: cover, lai, hveg,jd,snow,SWP,Ts", &
                  cover, lai, hveg, &  ! tmp
                  daynumber, snow(i,j), SWP(daynumber), &
                  Ts_C
           end if ! debug_flag
        end if ! DEBUG

        ustar_loc = ustar_nwp       ! First guess = NWP value
        invL      = invL_nwp        ! First guess = NWP value

       call Get_Submet(hveg,t2(i,j),Hd,LE, psurf(i,j), &
               z_mid(i,j,KMAX_MID), u_ref, q(i,j,KMAX_MID,1), & ! qw_ref
                       debug_flag,                        &    ! in
                       ustar_loc, invL,                   &    ! in-out
                       z0,d, Ra_ref,Ra_1m,rh,vpd)                    ! out

    if ( DEBUG_UK .and. debug_flag ) then
             print "(a10,2i4,2f7.2,2es12.3,3f8.3)", &
             "UKDEP SUB", me, &
             lu, ustar_nwp, ustar_loc, invL_nwp, invL, Ra_ref, Ra_1m,rh
    end if


         call Rsurface(lu,debug_flag, &
                   landuse_SGS(i,j,ilu),&
                   landuse_EGS(i,j,ilu), &
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
                   current_date%month, &
                   current_date%day, &
                   current_date%hour, &
                   Ra_ref,   &
                   g_sto,  &
                   Rsur, &
                   Rb)

         !/... add to grid-average Vg:


         do n = 1, NDRYDEP_CALC
           Vg_ref(n) = 1.0 / ( Ra_ref + Rb(n) + Rsur(n) )
           Vg_1m(n)  = 1.0 / ( Ra_1m  + Rb(n) + Rsur(n) )
           Grid_Vg_ref(n) = Grid_Vg_ref(n) + cover * Vg_ref(n)
           Grid_Vg_1m(n)  = Grid_Vg_1m(n)  + cover * Vg_1m(n)
           if ( Vg_ref(n) < 0.0 .or. Vg_ref(n) > 0.1 ) then
               print *, "VGREF ERROR", me, n, Vg_ref(n)
               print *, "VGREF ERROR Ras", Ra_ref, Rb(n), Rsur(n)
               print *, "VGREF ERROR stab", z0, d, ustar_loc, invL, g_sto
               call gc_abort(me, NPROC, "VGREF ERROR")
           end if
         end do
           if ( Vg_ref(CDEP_O3) > 0.04 ) then
               print *, "VGREFO3 ERROR", me, n, Vg_ref(n)
               print *, "VGREFO3 ERROR Ras", Ra_ref, Rb(n), Rsur(n)
               print *, "VGREFO3 ERROR stab", z0, d, ustar_loc, invL, g_sto
               call gc_abort(me, NPROC, "VGREFO3 ERROR")
           end if

         Sumcover = Sumcover + cover


       ! Stomatal flux = G_sto/G_Sur * Vg * O3
       !               = ( LAI * g_Sto) / Gsur * Vg * O3
       !               = LAI*g_Sto*Rsur * Vg * O3

    !!real, parameter ::  MOLARVOL = 1000.0/22.41  = 44.6
          !Vg_sto_ref(lu) =  lai * g_sto*Rsur(CDEP_O3) * Vg_ref(CDEP_O3)

          if ( lu == 1 ) then   ! Conif. forest

             nmole_o3 = xn_adv(IXADV_O3,i,j,KMAX_MID) * 41.0 * PPBINV


             dep_flux(DEP_VG_REF,i,j,1) = Vg_ref(CDEP_O3)
             dep_flux(DEP_VG_1M ,i,j,1) = Vg_1m (CDEP_O3)
             dep_flux(DEP_VG_STO,i,j,1) = lai * g_sto*Rsur(CDEP_O3) * &
                                             Vg_ref(CDEP_O3)
             dep_flux(DEP_FL_REF,i,j,1) = Vg_ref(CDEP_O3) * nmole_o3
             dep_flux(DEP_FL_STO,i,j,1) = lai * g_sto*Rsur(CDEP_O3) * &
                                             Vg_ref(CDEP_O3) * nmole_o3

             d_2d(D2_VG_REF,i,j,IOU_INST) = Vg_ref(CDEP_O3)
             d_2d(D2_VG_1M ,i,j,IOU_INST) = Vg_1m (CDEP_O3)
             d_2d(D2_VG_STO,i,j,IOU_INST) = lai * g_sto*Rsur(CDEP_O3) * &
                                             Vg_ref(CDEP_O3)
             d_2d(D2_FX_REF,i,j,IOU_INST) = Vg_ref(CDEP_O3) * nmole_o3
             d_2d(D2_FX_STO,i,j,IOU_INST) = lai * g_sto*Rsur(CDEP_O3) * &
                                             Vg_ref(CDEP_O3) * nmole_o3

             if ( DEBUG_UK .and. debug_flag ) then
                print "(a10,4es12.3)","OZONE", &
                    nmole_o3, &
                    PPBINV,   &
                   xn_adv(IXADV_O3,i,j,KMAX_MID), &
                   xn_adv(IXADV_O3,i,j,KMAX_MID) * PPBINV

                print "(a10,3i3,3f8.5, f9.2, 2f10.3,)", "FLUXES: ", &
                  current_date%month, &
                  current_date%day, &
                  current_date%hour, &
                  Vg_ref(CDEP_O3), Vg_1m(CDEP_O3), &
                   lai * g_sto*Rsur(CDEP_O3) * Vg_ref(CDEP_O3), &
                  nmole_o3, &
                  dep_flux(DEP_FL_REF,i,j,1), dep_flux(DEP_FL_STO,i,j,1)
             end if
          end if

          !flux(lu) =  lai * g_sto*Rsur(CDEP_O3) * Vg_ref(CDEP_O3)

!IPM HERE
        !/-- only grab gradients over land-areas

         if ( lu /= LU_WATER ) then
            Sumland = Sumland + cover
            do n = 1, NDRYDEP_CALC 
                if( Vg_1m(n) == 0.0 .or. Vg_1m(n)<Vg_ref(n) ) then
                     print  "(a15,i3,2i4, i5,i4)", "Vg LANDcoorsd ", &
                         me, i, j, i_glob(i), j_glob(j)
                     print  "(a15,i3,2f12.4)", "Vg LANDVgs ", &
                              n, Vg_1m(n), Vg_ref(n)
                     print  "(a15,3f12.4)", "Vg LANDcovers ",  &
                               cover, Sumland, Sumcover
                     print  "(a15,4f12.4)", "Vg LANDRs ",  &
                               Ra_ref, Ra_1m, Rb(n), Rsur(n)
                     call gc_abort(me,NPROC,"Vg PROBLEMS!")
                end if
                Vg_ratio(n) =  Vg_ratio(n) + cover * Vg_ref(n)/Vg_1m(n)
            end do
         else
            do n = 1, NDRYDEP_CALC 
               sea_ratio(n) =  Vg_ref(n)/Vg_1m(n)
               if( Vg_1m(n) == 0.0 .or. Vg_1m(n)<Vg_ref(n) ) then
                    call gc_abort(me,NPROC,"Vg Sea PROBLEMS!")
               end if
            end do
         end if

        if ( DEBUG_UK .and. debug_flag ) then
            print "(a14,2i4,f8.3,5f7.2,f12.3)", "UKDEP RATVG", ilu, lu, cover, &
                Vg_1m(2), Vg_ref(2), Sumland, Sumcover, Vg_ratio(2)
            print "(a14,2i4,f8.3,5f7.2,f12.3)", "UKDEP FINVG", ilu, lu, cover, &
                   lai, &   ! tmp
                   100.0*g_sto, &  ! tmp, in cm/s 
                Ra_ref, Rb(2), Rsur(2), 100.0/(Ra_ref + Rb(2) + Rsur(2) )
        end if
         
       !=======================
       !=======================
        end do LULOOP
       !=======================
       !=======================

        if ( DEBUG_UK .and. Sumland > 1.011  ) then
            print *, "SUMLAND ", me, nlu, i,j,i_glob(i), j_glob(j)
            print "(a10,f12.6)", "SUMLAND ", Sumland
            cover = 0.0
            do ilu= 1, nlu
              lu  = landuse_codes(i,j,ilu)
              cover = cover + landuse_data (i,j,ilu)
              print "(a10,i4,f8.3,f12.4)", "SUMLAND",lu,&
                  landuse_data (i,j,ilu), cover
             end do
             call gc_abort(me,NPROC,"SUMLAND!")
        end if


        if ( Sumland > 0.01 ) then
            gradient_fac(:) = Vg_ratio(:) / Sumland
        else
            gradient_fac(:) = sea_ratio(:)
        end if

        if ( DEBUG_UK .and. debug_flag ) then
            print "(a14,i2,3i3,9f8.3)", "UKDEP VG_UKR", &
                   snow(i,j),          &
                   current_date%month, &
                   current_date%day, &
                   current_date%hour, &
                 (100.0*Grid_Vg_ref(n), n = 1, min(5,NDRYDEP_CALC))
            print "(a14,i3,3i3,9f8.3)", "UKDEP VG_UK1", &
                   snow(i,j),          &
                   current_date%month, &
                   current_date%day, &
                   current_date%hour, &
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

   !WHY??? DepLoss(:) = 0.0 !FIXXX
             DepLoss(nadv) =   & ! Use directly set Vg
                 ( 1.0 - exp ( -Dep(n)%vg * dtz ) ) * xn_adv( nadv,i,j,KMAX_MID)
             cfac(nadv, i,j) = 1.0   ! Crude, for now.
  
         else
             DepLoss(nadv) =   vg_fac( ncalc )  * xn_adv( nadv,i,j,KMAX_MID)
             cfac(nadv, i,j) = gradient_fac( ncalc )
         end if
         if ( DepLoss(nadv) < 0.0 ) then
           call gc_abort(me,NPROC,"NEG DEPLOSS")
         end if
         if ( DepLoss(nadv) > xn_adv( nadv,i,j,KMAX_MID) )  then
           call gc_abort(me,NPROC,"NEG XN DEPLOSS")
         end if

         xn_adv( nadv,i,j,KMAX_MID) = &
             xn_adv( nadv,i,j,KMAX_MID) - DepLoss(nadv)


      !..accumulated dry deposition per grid square and summed over the whole
      !  domain

         totddep( nadv ) = totddep (nadv) + Deploss(nadv) * convfac

        if ( DEBUG_UK .and. debug_flag ) then
          if ( vg_set(n) ) then
              print "(a30,2i4,f8.3)", "DEBUG DryDep SET ",  n,nadv, Dep(n)%vg
          else
              print "(a30,3i4)", "DEBUG DryDep n, adv, calc ",  n,nadv, ncalc
              print "(a20,f12.5)", "DEBUG grad fac", gradient_fac( ncalc)
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

       call Add_ddep(i,j,convfac2)

     enddo   !j
   enddo    !i
 end subroutine drydep

end module DryDep_ml
