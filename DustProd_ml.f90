! <DustProd_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
!* 
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*  
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!* 
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!* 
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************! 


 module DustProd_ml

  !-- Calculates wind blown dust production based on works by
  !   Zender et al. (2003). JGR, 108, 4416-4437
  !   Alfaro & Gomes (2001). JGR, 106, 18075-18084 
  !   Gomes et al. (2003). Catena, 52, 257-271.

 use CheckStop_ml,         only : CheckStop
 use EmisDef_ml,           only : NDU, QDUFI, QDUCO
 use Functions_ml,         only : ERF, ERFfunc
 use ChemChemicals_ml,     only : species
 use ChemSpecs_tot_ml,     only : DUST_NAT_F
 use GridValues_ml,        only : glat, glon, glat_fdom, glon_fdom, i_fdom, j_fdom 
 use Io_ml,                only : PrintLog
 use Landuse_ml,           only : LandCover, NLUMAX 
 use LandDefs_ml,          only:  LandType
 use LocalVariables_ml,    only : Sub, Grid
 use MetFields_ml,         only : z_bnd, z_mid, u_ref, ustar_nwp, roa,  &
                                  t2_nwp, sdepth, fh, ps, surface_precip,  &
                                  rho_surf, SoilWater, foundSoilWater,  &
                                  clay_frac, sand_frac !ACB snow_flag
 use ModelConstants_ml,    only : KMAX_MID, KMAX_BND, dt_advec, METSTEP, &
    NPROC, MasterProc, USE_DUST
 use MicroMet_ml,          only : Wind_at_h
 use Par_ml,               only : me,MAXLIMAX,MAXLJMAX
 use PhysicalConstants_ml, only : GRAV,  AVOG, PI, KARMAN, RGAS_KG, CP
                                  !! ECO_CROP, ECO_SEMINAT, Desert=13, Urban=17
 use TimeDate_ml,          only : daynumber
 use Setup_1dfields_ml,    only : rh

!-----------------------------------------------
  implicit none
  private

  public ::  WindDust       

  real, public           :: DU_prod (NDU,MAXLIMAX,MAXLJMAX)  
  real, public           :: DUST_flux (MAXLIMAX,MAXLJMAX)
  real, public, save     :: kg2molecDU, m_to_nDU, frac_fine, frac_coar,  &
                            soil_dns_dry, help_ustar_th
  real, parameter        :: soil_dens = 2650.0  ! [kg/m3]
  logical, private, save :: my_first_call = .true.
  integer, save          :: dry_period(MAXLIMAX, MAXLJMAX) = 72
  logical, parameter, private :: DEBUG_DUST = .false.
  character(len=20)           :: soil_type
  contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   subroutine WindDust (i,j,debug_flag)

   implicit none

   integer, intent(in) :: i,j            ! coordinates of column
   logical, intent(in) :: debug_flag

   integer, parameter  :: Ndust = 2      ! number of size classes 
   integer, parameter  :: LU_DESERT = 13
   real   , parameter  :: Ro_water = 1000.0 
   real, parameter, dimension(Ndust) ::                    &
                    dsoil = (/ 1.5, 6.7/)   & ! diameter of dust particles [mkm]
                   ,mfrac = (/0.05, 0.45/)    ! mass fraction of the total mass

   real, parameter:: D_opt = 75.e-6,     &  ! [m]
                     Dm_soil = 210.0e-6,  & ! [m] MMD of the coarsest soil  (100)
                     Z0s = Dm_soil/30.0 , & ! [m] Smooth roughness length MaB95 p.16426, 
                                            !     MaB97 p.4392, GMB98 p.6207
                     !z0 = 0.5e-3,  & !(for desert..) 1.e-4 saltation roughness length 
                     z10 = 10.0             ! Z=10m
   real ::  Mflux = 0.0 , Nflux = 0.0
   real ::  delz, cover, z0, vh2o_sat, gr_h2o, v_h2o, ustar_moist_cor   &
          , gwc_thr, dust_lim, soil_dns_dry, ustar_z0_cor, u_th10, u10  &
          , ustar_salt, alfa, ustar_th, uratio, ustar, invL, ln10       &
          , owens, clay, frac_fin, frac_coa, flx_hrz_slt,  flx_vrt_dst

   logical :: arable, dust_prod = .false.
   integer :: n, ii, jj,nlu, ilu, lu

!_______________________________________________________
    if ( USE_DUST .eqv. .false. ) then
        call PrintLog("Skipping soil dust")
        return
    end if

    if ( my_first_call.and.MasterProc ) then 
     write(6,*)'***    Call for init_dust     ***   '

        call init_dust

          kg2molecDU = 1.0e-3 * AVOG / species(DUST_NAT_F)%molwt 
          my_first_call = .false.

    end if !  my_first_call
!_______________________________________________________


 if(DEBUG_DUST .and. debug_flag) write(6,*)'***   WINDBLOWN DUST',  &
                                 i_fdom(i), j_fdom(j)

 if ((glat(i,j)>65.0 .and. glon(i,j)>50.0)) return ! Avoid dust production in N. Siberia

!/..  landuse types ..........

  DU_prod (:,i,j) = 0.0  
  DUST_flux (i,j) = 0.0
  flx_vrt_dst = 0.0   ! vertical dust flux 
  Mflux = 0.0  


!/.. Crude assumption: dust generation if Pr < 2.4 mm per day (surpace_prec[mm/h])
 NO_PRECIP:  if (surface_precip(i,j) < 0.1) then
                  dry_period(i,j) = dry_period(i,j) + 1

   if(DEBUG_DUST .and. debug_flag) write(6,'(a30,i5,es12.3)')'>> NO RAIN >>', &
                                   dry_period(i,j),surface_precip(i,j)

!/.. Crude assumption: Dry soil after 48 hours since precipitation
  DRY:  if (dry_period(i,j)*dt_advec/3600.0 >= 48.0) then

     if(DEBUG_DUST .and. debug_flag) &
     write(6,'(a30,f10.4)')'>> DRY >>', dry_period(i,j)*dt_advec/3600.0

!/.. No dust production when soil is frozen, or covered by snow,
!/.. or wet (crude approximation by surface Rh)
!ACB Grid%snow == 1
  FROST: if ( Grid%t2C > 0.0 .and.  Grid%sdepth > 0 .and.  & 
                                    rh(KMAX_MID) < 0.85)  then

    if(DEBUG_DUST .and. debug_flag) write(6,'(a25,2f10.2,i4)')   &
          '>>      FAVOURABLE >>', Grid%t2C, rh(KMAX_MID), Grid%sdepth 
!ACB Grid%snow


!==  Land-use loop ===============================================

 nlu = LandCover(i,j)%ncodes
 LUC: do ilu= 1, nlu
        lu      = LandCover(i,j)%codes(ilu)
        cover   = LandCover(i,j)%fraction(ilu)
        arable = .false.

!/.. Dust production from arables (crops outside the growing period)
   if ( LandType(lu)%is_crop) then  
     if (daynumber <  LandCover(i,j)%SGS(ilu) .or.   &
         daynumber >  LandCover(i,j)%EGS(ilu) )      &
             arable = .true.
   endif

!/.. Dust erosion from Crops (Arable) and Desert (Mediterr.Scrubs lu==11 ???)
!/.. Creates problems on e.g. Greenland, as Bare land is included in Desert

   DUST: if (arable .or. lu == LU_DESERT) then

     if(DEBUG_DUST .and. debug_flag) write(6,'(a15,i5,f10.3)')   &
                                     '-----> Landuse', lu, cover

!//  Sandblasting coefficient  Alfa [1/m] = Flux_vert/Flux_horiz
!--------------------------------------------------------------
!   Dust production is EXTREMELY sensitive to this parameter, which changes 
!   flux by 3 orders of magnitude in 0.0 < fr_clay < 0.20
!   However, larger clay content decreases sandblasting effect due to coalescing
!   crust formation, which traps loose soil aggregates (Gomes etal. Catena, 2003)  
!   (CHIMERE uses C = alfa*crust_fac = 5e-5 * 4e-3 = 2e-7)
!--------------------------------------------------------------

!//__ Clay content in soil (limited to 0.2 ___
     clay = min (clay_frac(i,j), 0.20) ! [frc]

!== Testing ====
!    ln10 = log(10.0) 
!    alfa = 100.0 * exp(ln10 * (13.4 * clay - 6.0))  ! [m-1]     MaB95 p. 16423 (47)
!    if (DEBUG_DUST .and. debug_flag)  &
!    write(6,'(a15,2f10.5,a10,e10.3)') '>>  CLAY  >>',clay_frac(i,j),  &
!                                       clay, '  Alfa: ', alfa

!! >>>>>  For now values to ALFA are assigned below 


!//__ U* from NWP model ____
     ustar = Grid%ustar 

!//__ Threshold U* for saltation for particle (D_opt=70um) - dry ground

     ustar_th = help_ustar_th / sqrt(roa(i,j,KMAX_MID,1)) ! [m s-1] 

     if(DEBUG_DUST .and. debug_flag) write(6,'(a20,3f15.5)')  &
                '>> U*/U*th/ro >>', ustar,ustar_th,roa(i,j,KMAX_MID,1)

!//___  Inhibition of saltation by soil moisture

   ustar_moist_cor = 1.0

   if (foundSoilWater) then

!__ Saturated volumetric water content (sand-dependent)
     vh2o_sat = 0.489-0.126 * sand_frac(i,j)       ! [m3 m-3] 
!__ Bulk density of dry surface soil [Zender, (8)]
     soil_dns_dry = soil_dens * (1.0 - vh2o_sat)    ! [kg m-3]
!__ Volumetric soil water content [m3/m3] SoilWater[m/0.072m]
     v_h2o = SoilWater(i,j,1)/0.072
!__ Gravimetric soil water content [kg kg-1]  
     gr_h2o = v_h2o * Ro_water/soil_dns_dry       

     if(DEBUG_DUST .and. debug_flag) then
     write(6,'(a30,f8.2,2f12.4)') 'Sand/Water_sat/soil_dens',  &
                                   sand_frac(i,j),vh2o_sat,soil_dns_dry
     write(6,'(a30,3f15.5)') ' SW/VolW/GrW/ ',SoilWater(i,j,1),v_h2o,gr_h2o
     endif

!//__ Minimal soil moisture at which U*_thresh icreases
     gwc_thr = (0.17 + 0.14 * clay)* clay  *5.   ! [kg/kg] [m3 m-3]  
 
   else  !.. No SoilWater in met.input
     gwc_thr = 0.1 
  endif

     if (gr_h2o > gwc_thr) &
       ustar_moist_cor = sqrt(1.0 + 1.21 *  &
                 (100.*(gr_h2o - gwc_thr))**0.68) ! [frc] FMB99 p.155(15) 

     if(DEBUG_DUST .and. debug_flag) write(6,'(a25,e12.4,f10.4)')   &
                     '>> U*_moist_corr  >>',gwc_thr, ustar_moist_cor

! ===================================

  
!// Limitation of available erodible elements and aerodynamic roughness length
  if (lu == LU_DESERT)       then      ! ---------  desert -----
        soil_type = 'Saharan desert'
        z0 = 0.5e-4        !TEST
        dust_lim = 0.5     ! 1.0
        alfa = 2.0e-5 
 !//  European deserts: lat(gb), long (gl)
     if ( (glat(i,j) > 36.0 .and. glon(i,j) < 0.0) .or.   &  
          (glat(i,j) > 37.0 .and. glon(i,j) < 45.0) )    then 
        soil_type = 'European Arid'
        dust_lim = 0.05
        alfa = 1.3e-5 
!        alfa = 1.5e-5    ! As for TFMM spring 2005
     endif

!// limit emissions in the Spanish desert grid (covered with greenhouses??)
!     if ( (i_glob(i) == 102.0 .and. j_glob(j) == 18.0) )  & 
!         dust_lim = 0.02
!         alfa = 1.3e-5 
!     elseif (lu == 6)  then           ! ------  Mediter. crops --

 !// Crops
  else                                 ! ------  Crops -------  
        soil_type = 'Crops'
        z0  = max (0.1 * LandCover(i,j)%hveg(ilu), 0.001)
        dust_lim = 0.02
        alfa = 1.0e-5  !1.e-5 
  endif

!*** TEST as in CHIMERE
!        alfa = 5.0e-5  !1.e-5 
!*** TEST as in CHIMERE

!  else                                 ! ---------  temp/root crops ------
!        z0  = max (0.1 * landuse_hveg(i,j,ilu), 0.001)
!        dust_lim = 0.01
!        alfa = 0.8e-5  !1.e-5 


!//__ Inhibition of threshold velocity by roughness elements
!//   Roughness length momentum for erodible surfaces MaB95 p.16420, GMB98 p.6205

   ustar_z0_cor = 1.0 -          & ! [frc] MaB95 p. 16420, GMB98 p. 6207
                  log(z0 / Z0s) / log( 0.35 * (0.1/Z0s)**0.8 )

   if(DEBUG_DUST .and. debug_flag) write(6,'(a20,es12.3)') &
                            '>> ustar_zo_corr  >>', ustar_z0_cor

   ustar_z0_cor = min ( 1.0, max(0.0001,ustar_z0_cor) )

   call CheckStop(ustar_z0_cor <= 0.0 .or. ustar_z0_cor > 1.0,  &
                        "DUST ERROR: ustar_z0_cor out of range")

   ustar_z0_cor = 1.0 / ustar_z0_cor   ! [frc] 

   if(DEBUG_DUST .and. debug_flag) write(6,'(a25,3es12.3)')   &
                         '>> U*_zo_corr  >>',z0, Z0s, ustar_z0_cor


!//___   Final threshold friction velocity


!TEST       ustar_th =0.25  ! TEST

   ustar_th =                 & ! [m s-1] Threshold friction velocity for saltation
            ustar_th        * & ! [m s-1] Threshold for dry, flat ground
            ustar_moist_cor * & ! [frc] Adjustment for moisture
            ustar_z0_cor        ! [frc] Adjustment for roughness

   if (DEBUG_DUST .and. debug_flag)  then
      write(6,'(3es12.3)') ustar_th, ustar_moist_cor, ustar_z0_cor
      write(6,'(a35,f8.3,2(a8,f6.3))')  &
        'FINAL U*_th= ',ustar_th,'U*=',ustar,'Ratio=',ustar_th/ustar
   endif


    flx_hrz_slt = 0.0
!   flx_vrt_dst = 0.0


! >>>>>  Check for saltation to occur [Whi79 p.4648(19), MaB97 p.16422(28)]        

    if (ustar > ustar_th .and. ustar > 0.35) then   
         dust_prod = .true.


    if(DEBUG_DUST .and. debug_flag) write(6,'(a30,f8.2)')  &
         ' Saltation occur  => ',   ustar_th/ustar

!//__ Increase of friction velocity by surface roughness length
!//__                                 and friction speeds (Owens effect)

!//__ wind speed at 10m
    invL = -KARMAN * GRAV * Grid%Hd / ( CP*Grid%rho_s   &
                                       *ustar*ustar*ustar *Grid%t2)
!  u10 = u_ref(i,j) + ustar/KARMAN * (log(z10/z_ref) - PsiM(z10*invL) + PsiM(z_ref*invL))

    u10 = Wind_at_h (Grid%u_ref, Grid%z_ref, Z10, Sub(lu)%d,   &
                     Sub(lu)%z0,  Sub(lu)%invL )

!//__ Threshold 10 m wind speed [m s-1] for saltation
    u_th10 = ustar_th/KARMAN * log(z10/z0) 
!..or    u_th10 = u10 * ustar_th / ustar ! [m s-1]

    if(DEBUG_DUST .and. debug_flag)  write(6,'(a20,2f12.4,2e12.4)')   &
            'Ro/T2/Hd/invL', Grid%rho_s, Grid%t2, Grid%Hd, Grid%invL
    if(DEBUG_DUST .and. debug_flag)  write(6,'(a20,4f12.4)')          &
            'Uref/U10/ratio/Uth10', Grid%u_ref, u10, ustar_th/ustar, u_th10


  !== Saltation roughens the boundary layer, AKA "Owen's effect"
  !  GMB98 p. 6206 Fig. 1 shows observed/computed u* dependence on observed U(1 m)
  !  GMB98 p. 6209 (12) has u* in cm s-1 and U, Ut in m s-1, personal communication, 
  !  D. Gillette, 19990529
  !  With everything in MKS, the 0.3 coefficient in GMB98 (12) becomes 0.003 
  !  Increase in friction velocity due to saltation varies as square of 
  !  difference between reference wind speed and reference threshold speed
 
!//__  Friction velocity increase from saltation GMB98 p. 6209

    owens = 0.003* (u10 - u_th10)*(u10 - u_th10)   !  [m s-1]
!      ustar = ustar + owens ! [m s-1] Saltating friction velocity

    if(DEBUG_DUST .and. debug_flag)  write(6,'(a20,3f10.3)')   &
             'Owens effect ',ustar-owens, owens , ustar

!//__ Calculate U*th / U* ratio
    uratio = ustar_th / ustar 

!//___ Horizontal saltation flux [kg m-1 s-1] 

    flx_hrz_slt = ustar**3 * (1.0-uratio)*(1.0+uratio)*(1.0+uratio) &   
                    * Grid%rho_ref / GRAV 

    if (DEBUG_DUST .and. debug_flag) then
      write(6,*)' '
      write(6,'(a35,es12.3)')  ' Horizontal Flux => ',   flx_hrz_slt
    endif

    if(DEBUG_DUST .and. debug_flag) then
      write(6,'(/a25,4f8.2)') soil_type,glat(i,j),glon(i,j),glat_fdom(i,j),glon_fdom(i,j)
      write(6,'(a25,3f10.3,es10.2)') '>>  U*/U*t/Klim,alfa  >>',  &
            ustar,ustar_th, dust_lim, alfa
      write(6,'(a15,f10.3,2es12.3)') 'FLUXES:',uratio, flx_hrz_slt*1000.0,  &
            flx_hrz_slt*dust_lim*alfa *1000.0
    endif

!TEST  to limit the dust production
!     if (lu == LU_DESERT)  then 
!        flx_hrz_slt  = min(10.e-3, flx_hrz_slt* dust_lim )
!     else
!        flx_hrz_slt  = min(2.e-3, flx_hrz_slt* dust_lim  )
!     endif
!        flx_vrt_dst =  flx_vrt_dst + flx_hrz_slt * alfa * cover

!// ___  Vertical dust flux [kg/m2/s], scaled with area fraction and
!//       added for all landuses
!TEST         flx_vrt_dst =  flx_vrt_dst + flx_hrz_slt * alfa * dust_lim
  

!TEST..  to limit the dust production to reasonable (estimated) values
   if (lu == LU_DESERT)  then 
      flx_vrt_dst  =  flx_vrt_dst + min(1.e-7, flx_hrz_slt * alfa * dust_lim)
   else
      flx_vrt_dst  =  flx_vrt_dst + min(1.e-8, flx_hrz_slt * alfa * dust_lim)
   endif

   flx_vrt_dst =  flx_vrt_dst * cover

!.. Mass flus for test output
   Mflux = Mflux + flx_hrz_slt * alfa * dust_lim

   if(DEBUG_DUST .and. debug_flag)  then
    write(6,'(a35,es12.3/)')  ' Vertical Flux   => ',  Mflux
    write(6,'(a35,es12.3,i4,f8.3)')  'DUST Flux   => ',   flx_vrt_dst, lu, cover
    write(6,'(a15,f10.3,2es12.3)') 'FLUXES:',uratio, flx_hrz_slt*1000., flx_vrt_dst*1000.
   endif

   endif  ! U* > U*_threshold

   endif DUST

   enddo LUC ! landuse

  endif FROST
  endif DRY 

 else  ! PREC
    dry_period(i,j) = 0

    if(DEBUG_DUST .and. debug_flag) write(6,'(a30,i5,es12.3)')   &
        '>> RAIN-RAIN >>', dry_period(i,j),surface_precip(i,j)
 endif NO_PRECIP

!//__  N production [ 1/m2/s]:   d3(mkm->m) * 1e-18
!TEST       Nflux(n) = Mflux(n) *m_to_nDU / dsoil(n)**3 *1.e18

  if (dust_prod) then

!!.. vertical dust mass flux: fine and coarse [kg/m2/s] -> [ng/m3/s] 
! fine  (0.10-2006; 0.05-2007)
     DU_prod(1,i,j) = 0.05 * flx_vrt_dst * kg2molecDU /Grid%DeltaZ 
! coarse (0.23-2006; 0.15-2007)
     DU_prod(2,i,j) = 0.17 * flx_vrt_dst * kg2molecDU /Grid%DeltaZ   

! TEST: stuff below needs to be tested
!   call get_dustfrac(frac_fin, frac_coa, ustar) 
!!.. vertical dust flux [kg/m2/s] -> [kg/m3/s]*AVOG/M e-3 -> [molec/cm3/s] 
!   DU_prod(1,i,j) = frac_fin * flx_vrt_dst * kg2molecDU / delz   
!   DU_prod(2,i,j) = frac_coa * flx_vrt_dst * kg2molecDU / delz  

!//__Dust flux [kg/m2/s] for check
    DUST_flux(i,j) =   flx_vrt_dst ! * (frac_fin + frac_coar) 
 
    dust_prod = .false.   ! Zero-setting

    if(DEBUG_DUST .and. debug_flag) &
      write(6,'(//a15,2es12.4,a15,e12.4,2f6.3)') '>>M-DU prod >>',   &
      DU_prod(1,i,j),  DU_prod(2,i,j), ' > TOTAL >', DUST_flux(i,j), &
      frac_fin, frac_coa
  endif  ! dust_prod

  if(DEBUG_DUST .and. debug_flag) write(6,'(//a15,2es12.4,a15,e12.4,2f6.3)') &
     '<< DUST OUT>>', DU_prod(1,i,j),  DU_prod(2,i,j), ' > TOTAL >',         &
     DUST_flux(i,j), frac_fin, frac_coa

   end subroutine WindDust
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine init_dust

!_____ Initialization, assignments etc.
!  Calculation of fine/coarse dust fractions - not used yet
 
  implicit none

  integer :: i, k
  real    :: Re_opt, k_help1, k_help2, help_ust
  real, parameter   :: D_opt = 75.0e-6  ! [um] Optimal particle size for uplift
!......................................
  integer, parameter  :: Nsoil = 3, Ndust = 4
  real    :: lnsig, sqrt2, dead_erf, ERFfun, x1,x2, y, tot_soil
  integer :: isoil, idu
  real, dimension(Nsoil,Ndust)  :: dif1, dif2
  real, dimension(Ndust)    :: sum_soil
  real, dimension (Ndust) :: d1 = (/0.1, 1., 2.5, 5. /)
  real, dimension (Ndust) :: d2 = (/1., 2.5, 5., 10. /)
!//__ different soil size distributions
!  real, dimension (3) :: mass_fr = (/ 2.6e-6,    0.781,    0.219/) 
!  real, dimension (3) :: dsoil   = (/ 0.0111, 2.524, 42.10 /)  ! [um] MMD 
!  real, dimension (3) :: sig_soil= (/ 1.89 ,     2.0 ,      2.13    /)  ! [frc] GSD
!  real, dimension (3) :: dsoil   = (/  1.5 ,   6.7,   14.2 /)  ! [um] MMD 
!  real, dimension (3) :: sig_soil= (/  1.7 ,  1.6 ,    1.5 /)  ! [frc] GSD
  real, dimension (Nsoil) :: dsoil   = (/ 0.832 ,  4.82 , 19.38 /)
  real, dimension (Nsoil) :: mass_fr = (/ 0.036,  0.957,  0.007 /) 
  real, dimension (Nsoil) :: sig_soil= (/ 2.10 ,  1.9 ,   1.60  /)
!---------------------------------------------

  if (me == 0) then
   write(6,*)
   write(6,*) ' >> DUST init <<',soil_dens
  endif

!//__ Reynold's number ( uses D_opt[cm] ) 
    Re_opt=0.38 + 1331. *(100.*D_opt)**1.56  ! [frc] "B" MaB95 p. 16417 (5)

 !.. Given Re*t(D_opt), compute time independent factors contributing to u*t

!//__ [frc] IvW82 p. 115 (6) MaB95 p. 16417 (4) Inter-particle cohesive forces
    k_help1 = 1.0 + 6.e-7 / ( soil_dens *GRAV *(D_opt**2.5) )      !   SQUARED      
    k_help2 = soil_dens * GRAV * D_opt
                             !   SQUARED
!   if (me == 0) write(6,'(a25,f5.1,3e12.4)')  &
!               'ROsoil/Re_opt/K1/K2 ',soil_dens,Re_opt, k_help1, k_help2

!//__ U_star_threshold
    if ( Re_opt < 0.03 ) then
      stop 'Dust: Reynolds < 0.03'

    else if ( Re_opt < 10.0 ) then
      help_ust = 1.928 * Re_opt**0.0922 - 1.0 ! [frc] IvW82 p. 114 (3), MaB95 p. 16417 (6)
      help_ust = 0.129 * 0.129 /  help_ust ! [frc]                       SQUARED

     if (me == 0)  write(6,'(a20,e12.4)') 'U* =', help_ust

    else 
      help_ust = 1.0- 0.0858 * exp(-0.0617 *(Re_opt-10.0)) ! [frc] IvW82 p. 114(3), 
                                                           ! MaB95 p. 16417 (7)
      help_ust = 0.12*0.12 * help_ust*help_ust  ! [frc]                  SQUARED

     if (me == 0) write(6,'(a20,e12.4)') 'U* =', help_ust

    endif     ! Re_opt < 0.03
         
!//__ This method minimizes the number of square root computations performed

    help_ustar_th = sqrt (k_help1 * k_help2 * help_ust)

    if (me == 0) write(6,'(a20,e12.4)') 'U*t help =', help_ustar_th

!// =======================================
!TEST       sand_frac = 0.6   ! sand fraction in soil
!// Saturated volumetric water content (sand-dependent)
!     vh2o_sat = 0.489-0.126 * sand_frac(i,j)       ! [m3 m-3] 
!// Bulk density of dry surface soil [Zender, (8)]
!     soil_dns_dry = soil_dens * (1.0 - vh2o_sat)    ! [kg m-3]

!//__ Calculate mass dust fractions : fine and coarse
    if (me == 0) then
      write(6,*)
      write(6,*) ' >> DUST fractions <<', Nsoil, Ndust
      write(6,'(a15,3e12.4)') 'Sigma =', (sig_soil(i),i=1,Nsoil)
    endif

    sum_soil(:) = 0.0
    tot_soil    = 0.0

    sqrt2 = sqrt (2.0)

    do isoil = 1, Nsoil
       y = log (sig_soil(isoil) ) * sqrt2

!      if (me == 0) write (6,'(a15,i4,2e12.4)') &
!                   'TEST 1', isoil, log(sig_soil(isoil)),y

    do idu = 1, Ndust 

!      if (me == 0) write (6,'(a15,2i3,3f8.3,e12.4)')  &
!                  'TEST 2', isoil,idu,d1(idu),d2(idu),dsoil(isoil),y

       x1 = log ( d1(idu) / dsoil(isoil) ) / y
       x2 = log ( d2(idu) / dsoil(isoil) ) / y

       if (MasterProc) write (6,'(a,4e12.4)') 'DUST TEST 3', &
                               x1,x2,ERFfunc(x1),ERFfunc(x2)

       dif2(isoil,idu) = 0.5 * ( ERFfunc(x2) - ERFfunc(x1) )
       dif2(isoil,idu) =  dif2(isoil,idu) * mass_fr(isoil)

       sum_soil(idu) =  sum_soil(idu) + dif2(isoil,idu)
       tot_soil = tot_soil + dif2(isoil,idu)

     enddo
   enddo

   frac_fine =  sum_soil(1) + sum_soil(2) +  sum_soil(3) 
   frac_coar =  sum_soil(4)

!   if (me == 0)  then
!     do idu = 1,Ndust 
!      write (6,'(a25,2f8.4,3(f8.3),2f12.3)') ' Dust frac in bins:',  &
!             d1(idu), d2(idu), (dif2(isoil,idu), isoil=1,3),         &
!             sum_soil(idu),sum_soil(idu)/tot_soil
!     enddo
!     write (6,'(a30,2f8.4)') ' ** FINE / COARSE **',frac_fine, frac_coar
!   endif

  end subroutine init_dust
!>------------------------------------------------------<

!<------------------------------------------------------<
  subroutine get_dustfrac(frac_fine, frac_coarse, ustar)
!== Calculates fine/coarse dust fractions dependenden on U*,
!   based Alfaro & Gomes (2001) - under testing

   real :: frac_fine, frac_coarse, ustar, a
   integer :: iclass

   if (ustar < 0.35) then
      frac_fine   = 0.02
      frac_coarse = 0.09
   else if (ustar < 0.40) then
      a           = (ustar-0.35)/0.05
      frac_fine   = (1.0-a)*0.02 + a*0.04
      frac_coarse = (1.0-a)*0.09 + a*0.11
   else if (ustar < 0.55) then
      a           = (ustar-0.40)/0.15
      frac_fine   = (1.0-a)*0.04 + a*0.26
      frac_coarse = (1.0-a)*0.11 + a*0.30
   else if (ustar < 0.80) then
      a           = (ustar-0.55)/0.25
      frac_fine   = (1.0-a)*0.26 + a*0.35
      frac_coarse = (1.0-a)*0.30 + a*0.11
   else
      frac_fine   = 0.35
      frac_coarse = 0.11
   endif

 end subroutine get_dustfrac

! >----------------------------------------------------------<

 end module DustProd_ml

