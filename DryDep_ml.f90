! <DryDep_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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

module DryDep_ml

  ! Module started from the drag-coefficient based approach of BJ98:
  ! Berge, E. and  Jakobsen, H.A., A regional scale multi-layer model
  ! for the calculation of long-term transport and deposition of air
  ! pollution in Europe, Tellus B (1998), 50, 105-223.

  ! but has been extensively re-written since. See ....
  ! Emberson, L.,Simpson, D.,Tuovinen, J.-P.,Ashmore, M.R., Cambridge, H.M.",
  !  2000, Towards a model of ozone deposition and stomatal uptake over 
  !  Europe, EMEP MSC-W Note 6/2000,
  ! Simpson, D.,Tuovinen, J.-P.,Emberson, L.D.,Ashmore, M.R.,2001,
  !  "Characteristics of an ozone deposition module",WASP:Focus,1,253-262
  ! Simpson, D.,Tuovinen, J.-P.,Emberson, L.D.,Ashmore, M.R.,2003,
  !  "Characteristics of an ozone deposition module II: sensitivity analysis",
  !   WASP, 143, 123-137
  ! Tuovinen, J.-P.,Ashmore, M.R.,Emberson, L.D.,Simpson, D., 2004, "Testing
  !   and improving the EMEP ozone deposition module", Atmos.Env.,38,2373-2385

 !-- model specific dry dep values are set in My_DryDep_ml


 use My_DryDep_ml, only : Init_DepMap, &  ! Maps indices between 
                          ! Vg-calculated (CDDEP..) and advected  (IXADV_..)
                          NDRYDEP_GASES , &  ! No. Vd values calculated 
                          NDRYDEP_ADV, &   ! No. advected species affected
                          NDRYDEP_AER, &   ! No. aerosol size modes for Vd
                          NDRYDEP_CALC, &  ! Total No. of  Vd values
                          DRYDEP_GASES, &  ! Wesely Index Vd values calculated 
                          CDDEP_SET,    &   ! for so4
                          CDDEP_NO2,CDDEP_O3,    &   ! for NO2 comp pt. approach
                          CDDEP_SO2,CDDEP_NH3, & ! CoDep extra
                          CDDEP_PMf,   &   !
                          FLUX_CDDEP,   &   ! index O3 in CALC array, for STO_FLUXES
                          FLUX_ADV ,   &   ! index O3 in ADV  array, for STO_FLUXES
                          DepLoss, Add_MosaicOutput, &
                          DDepMap        ! Mapping (type = depmap)


use LandDefs_ml, only : LandDefs
use My_Derived_ml, only : METCONC_PARAMS &    ! ->  d_2d, IOU_INST, D2_VG etc...
   ,MMC_USTAR, MMC_RH, MMC_INVL  !

 use Aero_Vds_ml,  only : SettlingVelocity, GPF_Vds300, Wesely300
 use CheckStop_ml, only: CheckStop
 use Chemfields_ml , only : cfac, so2nh3_24hr,Grid_snow !hf CoDep!,xn_adv


 use ChemSpecs_adv_ml, only : NSPEC_ADV, IXADV_NO2, IXADV_SO2, IXADV_NH3
 use ChemSpecs_tot_ml, only : NSPEC_TOT
!LPJ  use DO3SE_ml,       only : Init_DO3SE, do3se, f_phen
 use DO3SE_ml,       only : do3se, f_phen
 use EcoSystem_ml,   only : EcoSystemFrac, Is_EcoSystem,  &
                             NDEF_ECOSYSTEMS, DEF_ECOSYSTEMS
 use GridValues_ml , only : GRIDWIDTH_M,xmd,xm2, gb,dA,dB, &
          debug_proc, debug_li, debug_lj, i_fdom, j_fdom   ! for testing
!LPJ use Io_Nums_ml,     only : IO_DO3SE
 use Landuse_ml,     only : Land_codes,LU_cdf, LandCover
 use LandDefs_ml,    only : LandType
 use LocalVariables_ml, only : Grid, Sub, L, iL ! Grid and sub-scale Met/Veg data
 use MassBudget_ml,  only : totddep,DryDep_Budget
 use MicroMet_ml,   only : AerRes, Wind_at_h
 use ModelConstants_ml,    only : dt_advec,PT,KMAX_MID, KMAX_BND ,&
                                  DEBUG_i, DEBUG_j, NPROC, &
                                  DEBUG_DRYDEP, DEBUG_VDS, MasterProc, &
                                  ATWAIR, atwS, atwN, PPBINV,&
                                  KUPPER, NLANDUSEMAX
 use Par_ml,               only : li0,li1,lj0,lj1, me
 use PhysicalConstants_ml, only : PI, KARMAN, GRAV, RGAS_KG, CP, AVOG
 
 use Landuse_ml,       only : SetLandUse  & 
                              ,NLUMAX &  ! Max. no countries per grid
                              ,LandCover   ! Provides codes, SGS, LAI, etc,

 use Rb_ml,        only: Rb_gas
 use Rsurface_ml
 use SoilWater_ml, only : SWP ! = 0.0 always for now!
 use Wesely_ml,    only : Init_GasCoeff !  Wesely stuff, DRx, Rb_Cor, ...
 use Setup_1dfields_ml, only : xn_2d,amk
 use StoFlux_ml,  only:   STO_FLUXES,  &   ! true if fluxes wanted.
                            unit_flux, &! = sto. flux per m2
                            lai_flux,  &! = lai * unit_flux
                            Setup_StoFlux, Calc_StoFlux  ! subs
                            !Init_StoFlux, Setup_StoFlux, Calc_StoFlux  ! subs
 use ChemSpecs_shl_ml,    only :  NSPEC_SHL
 use My_Aerosols_ml,    only : NSIZE
 use TimeDate_ml,       only : daynumber, current_date

!FEB2009 - Met_ml needed for temp 
 use MetFields_ml, only: u_ref
 use MetFields_ml, only: tau, sdepth, SoilWater, SoilWater_deep, rh2m, th,pzpbl
! use Sites_ml,         only : nlocal_sites,site_x,site_y,site_n,site_name  ! JP-DEPSITES

 implicit none
 private

 public :: drydep, init_drydep
 
 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO

  logical, public, dimension(NDRYDEP_ADV), save :: vg_set 

  logical, private, save :: my_first_call = .true.
  character(len=30),private, save :: errmsg = "ok"


 contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     subroutine init_drydep


     integer, save ::  old_daynumber = -99
     integer ::nadv,n
     integer :: i, j, lc, ilc, nlc, iEco   ! for EcoSystem stuff
     logical :: debug_flag  ! for EcoSystem stuff
     real    :: coverage    ! for EcoSystem stuff

  if ( my_first_call ) then 

     call Init_DepMap()               ! Maps CDDEP to IXADV
     call Init_GasCoeff()             ! Sets Wesely coeffs.

! Read data for DO3SE (deposition O3 and  stomatal exchange) module
! (also used for other gases!)
!LPJ     if(.not.LU_cdf)then
!LPJ        call Init_DO3SE(IO_DO3SE,"Inputs_DO3SE.csv",Land_codes, errmsg)
!LPJ        call CheckStop(errmsg, "Reading DO3SE ")
!LPJ     endif
!     call Init_StoFlux()              

     nadv = 0
     do n = 1, NDRYDEP_ADV  
         nadv       = max( DDepMap(n)%ind, nadv )  ! Looking for highest IXADV
         vg_set(n)  = ( DDepMap(n)%calc == CDDEP_SET ) ! for set vg
         !if ( DEBUG_DRYDEP .and. MasterProc ) write(*,*) "VGSET ", n, nadv, vg_set(n)
     end do

     my_first_call = .false.
     if( MasterProc  .and. DEBUG_DRYDEP) write(*,*) "INIT_DRYDEP day ", &
           daynumber, old_daynumber

!=============================================================================
! From EcoSystems, but caused some circularity problem
! use EcoSystem_ml, only :: EcoSystemFrac, Is_EcoSystem

   EcoSystemFrac(:,:,:) = 0.0
   do j = lj0, lj1
     do i = li0, li1
       debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj )

       nlc = LandCover(i,j)%ncodes
       LCLOOP: do ilc= 1, nlc
           lc       = LandCover(i,j)%codes(ilc)
           coverage = LandCover(i,j)%fraction(ilc)
           ECOLOOP: do iEco= 1, NDEF_ECOSYSTEMS
              if( Is_EcoSystem(iEco,lc) ) then
                 EcoSystemFrac(iEco,i,j) = EcoSystemFrac(iEco,i,j) + coverage
                if( debug_flag ) then
                     write(6,"(a,2i4,a12,3f10.4)") "ECOSYS AREA ",&
                     ilc, lc, "=> "//trim(DEF_ECOSYSTEMS(iEco)), &
                         coverage, EcoSystemFrac(iEco,i,j)
                end if
             end if
          end do ECOLOOP
        end do LCLOOP
      end do ! i
    end do ! j

!     invEcoFrac(:) = 0.0
!
!     do n = 0, size(DEF_ECOSYSTEMS)-1
!        if ( EcoFrac(n) > 1.0e-39 ) invEcoFrac(n) = 1.0/EcoFrac(n)
!     end do
!=============================================================================
  end if !  my_first_call

  end subroutine init_drydep

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine DryDep(i,j)
    integer, intent(in):: i,j

    real, dimension(NDRYDEP_GASES ) :: &
          Rb           & ! Quasi-boundary layer rsis.
         ,Rsur    &   ! Surface Resistance (s/m) 
         ,Gns !Surface conductance
         
    real, dimension(NDRYDEP_CALC) :: &
          gradient_fac & ! Ratio of conc. at zref (ca. 50m) and 3m
         ,vg_fac       & ! Loss factor due to dry dep.
         ,Vg_ref       & ! Vg at ref ht.
         ,Vg_3m        & ! Vg at  3m
         ,Vg_ratio     & ! Ratio Vg_ref/Vg_3m = ratio C(3m)/C(ref), over land
         ,sea_ratio     ! Ratio Vg_ref/Vg_3m = ratio C(3m)/C(ref), over sea

    real, dimension(NDRYDEP_GASES ) :: &
          Grid_Rs      & ! Grid average of Rsurface
         ,Grid_Gns       ! Grid average of Gns

    integer n, iiL, nlu, ncalc, nadv  ! help indexes
    integer :: imm, idd, ihh, iss     ! date
    integer :: nadv2d !index of adv species in xn_2d array

    real :: no2fac  ! Reduces Vg for NO2 in ration (NO2-4ppb)/NO2
!dsVDS    real :: RaVs    ! Ra_ref *Vs for particles

    real convfac,  & ! rescaling to different units
         convfac2, & ! rescaling to different units
         lossfrac,  & !  If needed in My_DryDep - not used now.
         dtz         ! scaling factor for veff ( = dt/z, where dt=timestep and 
                     ! z = height of layer)

    integer :: nae
    !dsVDS real, dimension(NSIZE):: aeRb, aeRbw , Vs
    !dsVDS real :: convec
    real, dimension(NSIZE)::  Vs


     real, save :: inv_gridarea  ! inverse of grid area, m2

      real ::  Sumcover, Sumland   ! Land-coverage
      logical :: debug_flag        ! set true when i,j match DEBUG_i, DEBUG_j
      real :: Vg_scale

 ! Ecosystem specific deposition requires the fraction of dep in each 
 !  landuse mosaic, iL:

      real, dimension(NDRYDEP_CALC,0:NLANDUSEMAX):: &
              Mosaic_Gsur  & ! Gsur for output, LC specific, 0=GRID
            , Mosaic_Gns & ! Gns for output, LC specific
            , Mosaic_VgRef & ! Vg at ref height, e.g. 45m
            , Mosaic_Vg3m  ! Vg at 3m

      real, dimension(size(METCONC_PARAMS),0:NLANDUSEMAX):: &
            Mosaic_Met  ! met, just 

      real, dimension(NSPEC_ADV ,NLANDUSEMAX):: fluxfrac_adv
      integer :: iL_used(NLUMAX)
      real :: wet, dry    ! Fractions
      real :: snow_iL !snow fraction for one landuse
      real :: Vds     ! FEB2009 Vds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Extra outputs sometime used. Important that this 
!! line is kept at the end of the variable definitions and the start of real
!!  code - allows both in .inc file
!! Uncomment and make .inc file as required
!   include 'EXTRA_LU_Setup.inc'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     first calculate the 3m deposition velocity following the same
!     procedure as in the lagmod. second the flux and the accumulated 
!     deposition is calculated.
!
!     effective dry deposition velocity applied to the model concentration
!     at the top of the constant flux layer, zdep 
!     Dry deposion rates are specified in subroutine readpar
!

! FOR DEBUGGING
  imm      =    current_date%month            ! for debugging
  idd      =    current_date%day              ! for debugging
  ihh      =    current_date%hour             ! for debugging
  iss      =    current_date%seconds          ! for debugging

     inv_gridarea = 1.0/(GRIDWIDTH_M*GRIDWIDTH_M) 


     ! - Set up debugging coordinates first. ---------------------------!
     ! If location matches debug i,j value, set debug_flag. Also passed
     ! to Rsurface_ml

      debug_flag= ( debug_proc .and. i == debug_li .and. j == debug_lj) 


     ! -----------------------------------------------------------------!
     !.and conversion factor,  convfac (( ps-pt)/grav... )  ===> 
     !      pressure in kg m-1 s-2
     ! 
      convfac = (dA(KMAX_MID) + dB(KMAX_MID)*Grid%psurf)&!dP
                 *xmd(i,j)/(ATWAIR*GRAV*inv_gridarea)

     ! -----------------------------------------------------------------!

!     !.and factor,  kg_air_ij (( ps-pt)/grav... )  ===> 
!     !      pressure in kg m-1 s
!     !      used for converting from mixing ratio to kg
!
!      kg_air_ij = (ps(i,j,1) - PT)*carea(KMAX_MID) = dP*dx**2/g
!     ! -----------------------------------------------------------------!
!
      lossfrac = 1.0 !  Ratio of xn before and after deposition


      dtz      = dt_advec/Grid%DeltaZ

      if ( DEBUG_DRYDEP .and. debug_flag ) then
         write(*,"(a26,4i4)") "UKDEP DryDep me, i,j ", me, i,j
         write(*,"(a10,i4,3i3,i6,10f10.3)") "UKDEP SOL", &
              daynumber, imm, idd, ihh, current_date%seconds, &
              Grid%zen, Grid%coszen, Grid%wetarea, &
              1.0e-5*Grid%psurf, Grid%Idiffuse, Grid%Idirect
         write(*,"(a10,i4,3i3,2f8.3,es12.4,f8.4)") "UKDEP NWP", &
              daynumber, imm, idd, ihh,  &
              Grid%Hd, Grid%LE, Grid%invL, Grid%ustar
      end if
      
      
    !/ Initialise Grid-avg Vg for this grid square:


    Mosaic_VgRef(:,:) = 0.0
    Mosaic_Vg3m(:,:)  = 0.0
    Mosaic_Gsur(:,:)  = 0.0
    Mosaic_Gns(:,:)   = 0.0
    Mosaic_Met(:,:)   = 0.0
    if(MMC_RH >0) Mosaic_Met(MMC_RH,0) = rh2m(i,j,1)  ! NWs 3 P output

    Vg_ratio(:) = 0.0
    Sumcover = 0.0
    Sumland  = 0.0
    fluxfrac_adv (:,:) = 0.0
    Grid_snow(i,j)=0.0 
    Grid_Rs(:) = 0.0
    Grid_Gns(:)  = 0.0
 
    !/ SO2/NH3 for Rsur calc
    Grid%so2nh3ratio = &
               xn_2d(NSPEC_SHL+IXADV_SO2,KMAX_MID) / & 
               max(1.0,xn_2d(NSPEC_SHL+IXADV_NH3,KMAX_MID))

    Grid%so2nh3ratio24hr = so2nh3_24hr(i,j)

    if ( DEBUG_DRYDEP .and. debug_flag ) then
         write(*,"(a,2i4,2es15.4)") "DRYDEP CONCS ", i,j, &
          xn_2d(NSPEC_SHL+IXADV_SO2,KMAX_MID), xn_2d(NSPEC_SHL+IXADV_NH3,KMAX_MID)
    end if

    if ( STO_FLUXES ) call Setup_StoFlux(daynumber, &
         xn_2d(NSPEC_SHL+FLUX_ADV,KMAX_MID),amk(KMAX_MID))

!dsVDS - can set settling velcoty here since not landuse dependent

        Vs = SettlingVelocity( Grid%t2, Grid%rho_ref )
        ! Restrict settling velocity to 2cm/s. Seems
        ! very high otherwise,  e.g. see Fig. 4, Petroff..

         Vs(:) = min( Vs(:), 0.02) 

    !/ And start the sub-grid stuff over different landuse (iL)

    nlu = LandCover(i,j)%ncodes
    LULOOP: do iiL= 1, nlu
        iL      = LandCover(i,j)%codes(iiL)

        iL_used (iiL) = iL    ! for eco dep

        f_phen   = LandCover(i,j)%fphen(iiL)


        L = Sub(iL)    ! ! Assign e.g. Sub(iL)ustar to ustar
        L%SGS = LandCover(i,j)%SGS(iiL)   !NOT NEEDED???
        L%EGS = LandCover(i,j)%EGS(iiL)


             if ( DEBUG_DRYDEP .and. debug_flag ) then
                write(6,"(a40,4i3,f6.1,2i4,3f7.3,2i4,2f6.2)") &
                    "DEBUG_veg: me,nlu,iiL,iL, lat, SGS, EGS ", &
                    me,nlu,iiL, iL, gb(i,j), L%SGS, L%EGS, &
                   L%coverage, L%LAI, L%hveg,daynumber, &
                   Grid%snow, SWP(daynumber),L%t2C

                write(6,"(a10,2i4,2f7.2,2es12.3,3f8.3)") "UKDEP SUB", me, &
                  iL, Grid%ustar, L%ustar, Grid%invL, &
                  L%invL, L%Ra_ref, L%Ra_3m,L%rh

             end if


         call Rb_gas(L%is_water, L%ustar, L%z0, DRYDEP_GASES ,Rb)

         call Rsurface(i,j,DRYDEP_GASES ,Gns,Rsur,errmsg,debug_flag,snow_iL)

         Grid_snow(i,j) = Grid_snow(i,j) +  L%coverage * snow_iL 

         if( MMC_USTAR >0) Mosaic_Met(MMC_USTAR,iL) = L%ustar !ds was 1
         if( MMC_INVL  >0) Mosaic_Met(MMC_INVL ,iL) = L%invL  !ds was 2
         if( MMC_RH    >0) Mosaic_Met(MMC_RH   ,iL) = L%rh    !ds was 3
!if( debug_flag ) then
!   write(6,"(a,i4,a,2f12.3)") "MOSAIC ", iL, &
!     trim(METCONC_PARAMS(3)), Mosaic_Met(3,iL), Mosaic_Met(3,0)
!end if


       !===================
       !// calculate dry deposition velocities for fine/coarse particles

        !dsVDS - set Vs earlier
        !dsVDS convec = Grid%wstar/L%ustar     ! Convection velocity scale  
        !dsVDS convec = convec * convec
        !dsVDS
        !dsVDS call Aero_Rb ( L%ustar, convec, Grid%rho_ref &
        !dsVDS              , Grid%u_ref, iL, Grid%snow, Grid%wetarea, L%t2   &   
        !dsVDS              , Vs, aeRb, aeRbw )
       !===================
       !Vds changes
       !Could be coded faster with Ra....
       ! u_hveg  = Wind_at_h( Grid%u_ref, Grid%z_ref, L%hveg, &
       !              L%d, L%z0, L%invL )
       !dsVDS   u_hveg = 1.0 ! FAKE for Wesely
       !dsVds call Aero_Vds( L%ustar,  L%invL, u_hveg, pzpbl(i,j), Vds )

       !===================


       !/... add to grid-average Vg:


         wet =   Grid%wetarea
         dry =   1.0 - wet


         do n = 1, NDRYDEP_CALC

            if ( n > NDRYDEP_GASES )  then    ! particles

                nae = n - NDRYDEP_GASES 

!dsVDS                RaVs = L%Ra_ref * Vs(nae)
!dsVDS
!dsVDS                Vg_ref(n) = Vs(nae) +      &
!dsVDS                  dry / (L%Ra_ref + aeRb(nae)  + RaVs  *aeRb(nae)  ) &
!dsVDS                + wet / (L%Ra_ref + aeRbw(nae) + RaVs  *aeRbw(nae) )
!dsVDS                    
!dsVDS                RaVs = L%Ra_3m  * Vs(nae)
!dsVDS
!dsVDS                Vg_3m(n)  = Vs(nae) +       &
!dsVDS                  dry / (L%Ra_3m + aeRb(nae)  + RaVs  *aeRb(nae)  ) &
!dsVDS                + wet / (L%Ra_3m + aeRbw(nae) + RaVs  *aeRbw(nae) )
!dsVDS
!dsVDS
!dsVDS                tmpref = Vg_ref(n) ! Store old values for comparison
!dsVDS                tmp3m  = Vg_3m (n)

!dsVDS        !FEB2009 Vds ================================================

              if ( LandType(iL)%is_forest  ) then ! Vds NOV08

                 !/ Use eqn *loosely* derived from Petroff results

                  Vds = GPF_Vds300(L%ustar,L%invL, L%SAI )

              else !!!  Vds NOV08

                !/  Use Wesely et al  for other veg & sea

                 ! Vds = Nemitz2004( 0.4, L%ustar, L%invL )
                 Vds = Wesely300( L%ustar, L%invL )

              end if

                !Vg_ref(n) = Vs(nae) + 1.0/ ( L%Ra_ref + 1.0/Vds  )
                !Vg_3m (n) = Vs(nae) + 1.0/ ( L%Ra_3m  + 1.0/Vds  )

                ! Use non-electrical-analogy version of Venkatram+Pleim (AE,1999)
                 Vg_ref(n) =  Vs(nae)/ ( 1.0 - exp( -( L%Ra_ref + 1.0/Vds)* Vs(nae)))
                 Vg_3m (n) =  Vs(nae)/ ( 1.0 - exp( -( L%Ra_3m  + 1.0/Vds)* Vs(nae)))


        if ( DEBUG_VDS ) then
            if ( debug_flag .or. (Vg_3m(n)>0.50 .or. Vg_ref(n)>0.50 )) then
              print "(a,5i3,2i4,2f7.3,f8.2,20f7.2)", "AEROCHECK", imm, idd, ihh, &
              n, iL, i_fdom(i), j_fdom(j), L%ustar,  L%invL, pzpbl(i,j), &
                Grid%ustar, Grid%Hd,  100.0*Vds, &
              100.0*Vs(nae), 100.0*Vg_ref(n),  100.0*Vg_3m (n) &
             , Grid%t2, Grid%rho_ref 
            
            call CheckStop((Vg_3m(n)>0.50 .or. Vg_ref(n)>0.50 ), "AEROSTOP")
          end if
        end if

        !FEB2009 Vds ================================================


            else                           ! gases

           ! After's Hilde's changes, no wet-dry needed

              Vg_ref(n) = 1. / ( L%Ra_ref + Rb(n) + Rsur(n) ) 

              Vg_3m (n) = 1. / ( L%Ra_3m + Rb(n) + Rsur(n) ) 


            endif

         ! Surrogate for NO2 compensation point approach, 
         ! assuming c.p.=4 ppb (ca. 1.0e11 #/cm3):        
         ! Note, xn_2d has no2 in #/cm-3

           if ( n == CDDEP_NO2 ) then

            no2fac = xn_2d(NSPEC_SHL+IXADV_NO2,KMAX_MID)   
            no2fac = max(1.0, no2fac)
            no2fac = max(0.00001,  (no2fac-1.0e11)/no2fac)

            Vg_ref(CDDEP_NO2) = Vg_ref(CDDEP_NO2) * no2fac
            Vg_3m (CDDEP_NO2) = Vg_3m (CDDEP_NO2) * no2fac
          end if ! CDDEP_NO2

           Mosaic_VgRef(n,iL) = Vg_ref(n)  ! Note iL, not iiL 
           Mosaic_Vg3m(n,iL) = Vg_3m(n)  ! Note iL, not iiL 
           Mosaic_Vg3m(n,0)  = Mosaic_Vg3m(n,0)  + L%coverage * Vg_3m(n)
           Mosaic_VgRef(n,0) = Mosaic_VgRef(n,0) + L%coverage * Vg_ref(n)

!hf rs
         if ( n <= NDRYDEP_GASES )  then    ! gases
             Mosaic_Gsur(n,iL) = 1.0/Rsur(n) ! Note iL, not iiL 
             Mosaic_Gns(n,iL) = Gns(n)  ! Note iL, not iiL 
          ! Grid averages. Wait with Rs
             Mosaic_Gsur (n,0) =  Mosaic_Gsur(n,0) + L%coverage / Rsur(n)
             Mosaic_Gns(n,0) =  Mosaic_Gns(n,0)+ L%coverage * Gns(n)
         endif
         end do !species loop

         Sumcover = Sumcover + L%coverage


        !/-- only grab gradients over land-areas

         if ( L%is_water ) then
            do n = 1, NDRYDEP_CALC
               sea_ratio(n) =  Vg_ref(n)/Vg_3m(n)
            end do
         else
            Sumland = Sumland + L%coverage
            do n = 1, NDRYDEP_CALC
                Vg_ratio(n) =  Vg_ratio(n) + L%coverage * Vg_ref(n)/Vg_3m(n)
            end do
         end if


        if ( DEBUG_DRYDEP .and. debug_flag ) then
            do n = 1 , NDRYDEP_GASES 
               write(*,"(a14,2i4,f7.3,i3,2f10.2,es12.2,2f8.2,a5,f8.3,2es18.6)") &
                  "UKDEP EXT: ", iiL, iL, L%coverage, n,&
                   L%LAI,100.0*L%g_sto, &  ! tmp, in cm/s 
                   L%Ra_ref, Rb(n), min( 999.0,Rsur(n) ),  &
                  " Vg: ", 100.0*Vg_3m(n), 100.0*Vg_ref(n), Vg_ratio(n)
            end do

        end if

       !
         
       !=======================

          if ( DEBUG_DRYDEP .and. debug_flag ) then
               write(*,*) "DEPSTO IN?: ", iL, LandType(iL)%flux_wanted , Sub(iL)%cano3_ppb
          end if
        if ( STO_FLUXES .and. LandType(iL)%flux_wanted ) then
             call Calc_StoFlux(iL,  Vg_ref(FLUX_CDDEP), debug_flag )
        end if ! STO_FLUXES
          if ( DEBUG_DRYDEP .and. debug_flag ) then
               write(*,*) "DEPSTO OUT?: ", iL, LandType(iL)%flux_wanted , Sub(iL)%cano3_ppb
          end if

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! Extra outputs sometime used for Sweden/IVL/SEI/CEH
          !! Uncomment and make .inc file as required

!            include 'EXTRA_LU_Outputs.inc'
!hf CoDep extra

!if  (debug_flag )write(6,*) "LANDCOVER ", iL, iiL,"!",LandDefs(iL)%code,"!"
       !=======================
        end do LULOOP
       !=======================
       !=======================



        if ( DEBUG_DRYDEP .and. Sumland > 1.011  ) then
            print *, "SUMLAND ", me, nlu, i,j,i_fdom(i), j_fdom(j), Sumland
            call CheckStop( "SUMLAND TOO BUG")
        end if


        if ( Sumland > 0.01 ) then
            gradient_fac(:) = Vg_ratio(:) / Sumland
        else
            gradient_fac(:) = sea_ratio(:)
        end if

        if ( DEBUG_DRYDEP .and. debug_flag ) then
            write(*, "(a14,i2,3i3,10f6.2)") "UKDEP VG_UKR", &
                   Grid%snow, imm, idd, ihh,  &
                 (100.0*Mosaic_VgRef(n,0), n = 1, min(5,NDRYDEP_GASES )), &
                 (100.0*Mosaic_Vg3m(n,0), n = 1, min(5,NDRYDEP_GASES ))
        end if


!-- loop through all affected advected species to calculate changes in
!   concentration (xn_adv), the conc. ratios (cfac), and deposition 

    do ncalc = 1, NDRYDEP_CALC

        vg_fac (ncalc) = 1.0 - exp ( -Mosaic_VgRef(ncalc,0) * dtz ) 

    end do ! n

      do n = 1, NDRYDEP_ADV 
         nadv    = DDepMap(n)%ind
         nadv2d  = NSPEC_SHL + DDepMap(n)%ind

         ncalc   = DDepMap(n)%calc

         if ( vg_set(n) ) then

             DepLoss(nadv) =   & ! Use directly set Vg
                 ( 1.0 - exp ( -DDepMap(n)%vg * dtz ) ) * xn_2d( nadv2d,KMAX_MID)
             cfac(nadv, i,j) = 1.0   ! Crude, for now.
  
         else
             DepLoss(nadv) =   vg_fac( ncalc )  * xn_2d( nadv2d,KMAX_MID)
             cfac(nadv, i,j) = gradient_fac( ncalc )
         end if

         if ( DepLoss(nadv) < 0.0 .or. &
              DepLoss(nadv)>xn_2d(nadv2d,KMAX_MID) ) then
             call CheckStop("NEGXN DEPLOSS" )
         end if


        xn_2d( nadv2d,KMAX_MID) = &
             xn_2d( nadv2d,KMAX_MID) - DepLoss(nadv)



        if ( STO_FLUXES .and. nadv == FLUX_ADV ) then
           ! fraction by which xn is reduced - used in
           ! safety measure:
             
              if( xn_2d( nadv2d,KMAX_MID)  > 1.0e-30 ) then
                  lossfrac = ( 1.0 - DepLoss(nadv)/ &
                                (DepLoss(nadv)+xn_2d( nadv2d,KMAX_MID)))
              end if
              if ( DEBUG_DRYDEP .and. lossfrac < 0.1 ) then
                  call CheckStop( lossfrac < 0.1, "ERROR: LOSSFRAC " )
                  !print *, "ERROR: LOSSFRAC ", lossfrac, nadv, nadv2d
              end if
        end if


      !.. ecosystem specific deposition - translate from calc to adv 
      !  and normalise

         do iiL = 1, nlu
            iL      = iL_used(iiL)

            if ( vg_set(n) )  then
               fluxfrac_adv(nadv,iL) = Sub(iL)%coverage  ! Since all vg_set equal
            else
               Vg_scale = Mosaic_VgRef(ncalc,iL)/ Mosaic_VgRef(ncalc,0)
               fluxfrac_adv(nadv,iL) = Sub(iL)%coverage*Vg_scale
            end if


          !=======================
          ! The fraction going to the stomata = g_sto/g_sur = g_sto * R_sur.
          ! Vg*nmole_o3 is the instantaneous deposition flux of ozone, and
          ! the actual amount "deposited" is obtained from DepLoss(O3) using
          ! fluxfrac as obtained above.

          ! Now, DepLoss is loss of molecules/cm3 over time-period dt_advec
          ! and depth z_bnd over 1m2 of grid. For sto fluxes we need to
          ! find values over 1m2 of vegeation (regardless of how much veg
          ! is in grid, so we don't need cover. Instead:


            if ( DEBUG_DRYDEP .and. debug_flag ) then
            if ( n == CDDEP_SO2 .and. iL == 1 ) then ! SO2, CF

               if ( vg_set(n) )  then
                 write(6,"(a12,3i3,3f12.3)") "FLUXSET  ", iiL, iL, nadv, &
                     100.0*DDepMap(n)%vg, Sub(iL)%coverage, fluxfrac_adv(nadv,iL)
               else
                 write(6,"(a12,3i3,f8.5,5f8.3)") "FLUXFRAC ", iiL, iL, nadv, &
                  Sub(iL)%coverage, &
                  100.0*Mosaic_VgRef(ncalc,0), & ! GRID
                  100.0*Mosaic_VgRef(ncalc,iL), &
                  100.0*Sub(iL)%coverage*Mosaic_VgRef(ncalc,iL), &
                   fluxfrac_adv(nadv,iL)
               end if
            end if !SO2 CF
            end if
         end do
          

      !..accumulated dry deposition per grid square and summed over the whole
      !  domain

!         totddep( nadv ) = totddep (nadv) + Deploss(nadv) * convfac


        if ( DEBUG_DRYDEP .and. debug_flag ) then
          if ( vg_set(n) ) then
              write(*, "(a30,2i4,f8.3)") "DEBUG DryDep SET ",  n,nadv, DDepMap(n)%vg
          else
              write(*, "(a30,3i4,f12.5)") &
                  "DEBUG DryDep n, adv, calc, fac ",  n,nadv, ncalc, gradient_fac( ncalc)
              write(*, "(a20,2e12.4)") &
                "DEBUG xn, DepLoss ", xn_2d(nadv2d,KMAX_MID), DepLoss(nadv)
              write(*, "(a20,2f8.4)") "DEBUG gv_fac( ncalc)", &
                 vg_fac(ncalc), 1.0-vg_fac(ncalc)
          end if
          !write(*,*) "XNSPEC DATES ", current_date
          !do ispec = 1, NSPEC_TOT
          !   write(*,"(a7,i3,es15.8)")  "XNSPEC ", ispec, xn_2d(ispec,20)
          !end do
        end if

       end do ! n

        call DryDep_Budget(i,j,Deploss,convfac)

       ! inv_gridarea = xm2(i,j)/(GRIDWIDTH_M*GRIDWIDTH_M)
       convfac2 = convfac * xm2(i,j) * inv_gridarea/amk(KMAX_MID)


      !.. Add DepLoss to budgets if needed:

       call Add_MosaicOutput(debug_flag,dt_advec,i,j,convfac2,&
           fluxfrac_adv, Mosaic_Vg3m,Mosaic_Gsur, Mosaic_Gns, Mosaic_Met)

       !DSX call Add_Vg(debug_flag,i,j, Mosaic_Vg3m)
       !DSX call Add_RG(debug_flag,i,j, Mosaic_Gsur, Mosaic_Gns)
       !DSX call Add_LCC_Met(debug_flag,i,j, Mosaic_Met)

 end subroutine drydep

end module DryDep_ml
