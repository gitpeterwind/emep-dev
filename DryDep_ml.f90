!*****************************************************************************! 
! <DryDep_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007-2011 met.no
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

  ! Dry deposition scheme uses a mosaic approach. Updated extensively
  ! in 2010-2011.
  ! History
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
  !
  ! Also, handling of dry/wet and co-dep procedure changed following discussions
  ! with CEH: Fagerli et al., in preperation....

 use My_Aerosols_ml,    only : NSIZE  

 use Aero_Vds_ml,  only : SettlingVelocity, GPF_Vds300, Wesely300
 use CheckStop_ml, only: CheckStop
 use Chemfields_ml , only : cfac, so2nh3_24hr,Grid_snow 


 use ChemChemicals_ml, only : species
 use ChemSpecs_adv_ml         ! several species needed
 use ChemSpecs_tot_ml, only : NSPEC_TOT, FIRST_SEMIVOL, LAST_SEMIVOL, &
                               NO2, SO2, NH3, O3
 use DO3SE_ml,         only : do3se, f_phen
 use EcoSystem_ml,     only : EcoSystemFrac, Is_EcoSystem,  &
                             NDEF_ECOSYSTEMS, DEF_ECOSYSTEMS
 use GridValues_ml ,   only : GRIDWIDTH_M,xmd,xm2, glat,dA,dB, &
          debug_proc, debug_li, debug_lj, i_fdom, j_fdom   ! for testing
 use Io_Progs_ml,      only : datewrite
 use Landuse_ml,       only : SetLandUse, Land_codes  & 
                              ,NLUMAX &  ! Max. no countries per grid
                              ,LandCover   ! Provides codes, SGS, LAI, etc,
 use LandDefs_ml,      only : LandType, LandDefs, STUBBLE
 use LocalVariables_ml,only : Grid, Sub, L, iL ! Grid and sub-scale Met/Veg data
 use MassBudget_ml,    only : totddep
 use MetFields_ml,     only : u_ref, rh2m
 use MetFields_ml,     only : tau, sdepth, SoilWater, SoilWater_deep, th,pzpbl
 use MicroMet_ml,      only : AerRes, Wind_at_h
 use ModelConstants_ml,only : dt_advec,PT,KMAX_MID, KMAX_BND ,&
                                  DEBUG_i, DEBUG_j, NPROC, &
                                  DEBUG_DRYDEP, DEBUG_ECOSYSTEMS, DEBUG_VDS,&
                                  MasterProc, &
                                  DEBUG_AOT, & !JUST TESTING
                                  ATWAIR, atwS, atwN, PPBINV,&
                                  KUPPER, NLANDUSEMAX

 use MosaicOutputs_ml,     only : Add_MosaicOutput, MMC_RH
 use OwnDataTypes_ml,      only : depmap
 use Par_ml,               only : li0,li1,lj0,lj1, me
 use PhysicalConstants_ml, only : PI, KARMAN, GRAV, RGAS_KG, CP, AVOG, NMOLE_M3
 use Rb_ml,                only : Rb_gas
 use Rsurface_ml
 use Setup_1dfields_ml,    only : xn_2d,amk, Fpart, Fgas
 use SoilWater_ml,         only : fSW !  =1.0 always for now with EC metdata. Can change with PARLAM
 use StoFlux_ml,  only:   unit_flux, &! = sto. flux per m2
                          lai_flux,  &! = lai * unit_flux
                          Setup_StoFlux, Calc_StoFlux  ! subs
 use ChemSpecs_shl_ml,  only :  NSPEC_SHL
 use TimeDate_ml,       only : daynumber, current_date
 use Wesely_ml         ! ... Init_GasCoeff, DRx, Rb_Cor, ...


 implicit none
 private

 public  :: DryDep, init_drydep
 private :: Init_DepMap
 
 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO

  ! Maps from adv index to one of calc indices
  integer, public, save, dimension(NSPEC_ADV) :: DepAdv2Calc 

  logical, private, save :: my_first_call = .true.
  character(len=30),private, save :: errmsg = "ok"

 ! WE NEED A FLUX_CDDEP, FLUX_ADV FOR OZONE;
 ! (set to one for non-ozone models)

  integer, public, parameter :: FLUX_CDDEP  = CDDEP_O3
  integer, public, parameter :: FLUX_ADV   = IXADV_O3
  integer, public, parameter :: FLUX_TOT   = O3


  logical, public, parameter :: COMPENSATION_PT = .false. 



!/**************************************************************************
!  Specifies which of the possible species (from Wesely's list)
!  are required in the current air pollution model   
!/**************************************************************************
   ! .... Define the mapping between the advected species and
   !      the specied for which the calculation needs to be done.
   !  We also define the number of species which will be deposited in
   ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_GASES
   ! The actual species used and their relation to the CDDEP_ indices
   ! above will be defined in Init_DepMap

       include 'CM_DryDep.inc'

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   logical, public, dimension(NDRYDEP_ADV), save :: vg_set 


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

     nadv = 0
     do n = 1, NDRYDEP_ADV  
         nadv       = max( DDepMap(n)%ind, nadv )  ! Looking for highest IXADV
         vg_set(n)  = ( DDepMap(n)%calc == CDDEP_SET ) ! for set vg
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
       debug_flag = ( DEBUG_ECOSYSTEMS .and. debug_proc .and. &
          i == debug_li .and. j == debug_lj )

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
  subroutine Init_DepMap
   integer :: iadv, i

     do i = 1, NDRYDEP_ADV  ! 22
      iadv = DDepMap(i)%ind
      if(DEBUG_DRYDEP .and. MasterProc) &
         write(6,*) "DEPMAP   ", DDepMap(i)%ind, DDepMap(i)%calc
      call CheckStop( iadv < 1, "ERROR: Negative iadv" )
      DepAdv2Calc(iadv) = DDepMap(i)%calc
    end do
  
   ! We process the various combinations of gas-species and ecosystem:
   ! starting with DryDep, e.g. DDEP_SO2_m2CF
 
     if(MasterProc.and.DEBUG_DRYDEP) write(6,*) "Init_DepMap D2D FINISHED"

  end subroutine Init_DepMap


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine DryDep(i,j)
    integer, intent(in):: i,j
    real, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost

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

    integer n, iiL, nlu, ncalc, nadv, nFlux  ! help indexes
    integer :: imm, idd, ihh, iss     ! date
    integer :: ntot !index of adv species in xn_2d array

    real :: no2fac  ! Reduces Vg for NO2 in ration (NO2-4ppb)/NO2

    real convfac,  & ! rescaling to different units
         convfac2, & ! rescaling to different units
         lossfrac,  & !  If needed in My_DryDep - not used now.
         dtz         ! scaling factor for veff ( = dt/z, where dt=timestep and 
                     ! z = height of layer)

    integer :: nae
    real, dimension(NSIZE)::  Vs


     real, save :: inv_gridarea  ! inverse of grid area, m2

      real ::  Sumcover, Sumland   ! Land-coverage
      logical :: debug_flag        ! set true when i,j match DEBUG_i, DEBUG_j
      real :: Vg_scale

      real, dimension(NSPEC_ADV ,NLANDUSEMAX):: fluxfrac_adv
      integer, dimension(NLUMAX)  :: iL_used, iL_fluxes
      real :: wet, dry         ! Fractions
      real :: snow_iL          !snow_flag fraction for one landuse
      real :: Vds              ! Aerosol
      real :: Vg_3mN           ! Crude nitrate correction

      real :: c_hveg, Ra_diff, surf_ppb  ! for O3 fluxes and Fst where needed
      real :: c_hveg3m  !TESTS ONLY
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
     ! conver molecules/cm3 to ppb for surface:
      surf_ppb   = PPBINV /amk(KMAX_MID)
      if ( DEBUG_AOT .and. debug_flag ) write(*,"(a,es12.4)") "CHAMK", surf_ppb

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
         write(*,"(a,4i4)") "DMET me, i,j ", me, i,j
         write(*,"(a,i4,3i3,i6,4f8.3,10f8.2)") "DMET SOL", &
              daynumber, imm, idd, ihh, current_date%seconds, &
              Grid%zen, Grid%coszen, Grid%wetarea, &
              1.0e-5*Grid%psurf, Grid%Idiffuse, Grid%Idirect
         write(*,"(a,i4,3i3,i6,2f8.2,es10.2,f8.3)") "DMET NWP", &
              daynumber, imm, idd, ihh, current_date%seconds, &
              Grid%Hd, Grid%LE, Grid%invL, Grid%ustar
      end if
      
      
    !/ Initialise Grid-avg Vg for this grid square:
      
      !call Init_GridMosaic(i,j)
      if( MMC_RH > 0 ) Grid%rh2m = rh2m(i,j,1) !NWP value


    Vg_ratio(:) = 0.0
    Sumcover = 0.0
    Sumland  = 0.0
    fluxfrac_adv (:,:) = 0.0
    Grid_snow(i,j)=0.0 
    !SUB0 Grid%Gsur = 0.0
    !SUB0 Grid%Gns  = 0.0
    !SUB0 Grid%Vg_Ref  = 0.0
    !SUB0 Grid%Vg_3m   = 0.0
    Sub(0)%Gsur = 0.0
    Sub(0)%Gns  = 0.0
    Sub(0)%Vg_Ref  = 0.0
    Sub(0)%Vg_3m   = 0.0
 
    !/ SO2/NH3 for Rsur calc
    Grid%so2nh3ratio = &
               xn_2d(SO2,KMAX_MID) / max(1.0,xn_2d(NH3,KMAX_MID))

    Grid%so2nh3ratio24hr = so2nh3_24hr(i,j)

  ! Surrogate for NO2 compensation point approach, assuming 
  ! c.p.=4 ppb (actually use 1.0e11 #/cm3):        
  ! Aiming at dC/dt = -Vg.(C-4)/dz instead of -Vg.C/dz
  ! factor difference is then:   C-4/C in ppb units
  ! Note, xn_2d has no2 in #/cm-3

    no2fac = 1.0 ! QUERY!!!!
    no2fac = max( 1.0, xn_2d(NO2,KMAX_MID) )
    no2fac = max(0.00001,  (no2fac-1.0e11)/no2fac)


    if ( DEBUG_DRYDEP .and. debug_flag ) then
         write(*,"(a,2i4,4es12.4)") "DRYDEP CONCS SO2,NH3,O3 (ppb) ", i,j, &
          xn_2d(SO2,KMAX_MID)*surf_ppb, xn_2d(NH3,KMAX_MID)*surf_ppb, &
            xn_2d(O3,KMAX_MID)*surf_ppb, no2fac
    end if

    call Setup_StoFlux( daynumber )

! - can set settling velcoty here since not landuse dependent

    Vs = SettlingVelocity( Grid%t2, Grid%rho_ref )

  ! Restrict settling velocity to 2cm/s. Seems
  ! very high otherwise,  e.g. see Fig. 4, Petroff et al., 2008 (Part I), where
  ! observed Vg for forests is usually < 2cm/s.

    if ( DEBUG_DRYDEP .and. debug_flag ) then
         call datewrite("DRYDEP VS",NSIZE,(/ Grid%t2, Grid%rho_ref, Vs /) )
   !      if( maxval(Vs) > 0.02 ) write(*,*) "DRYDEP LIM!"
    end if
   !  Vs(:) = min( Vs(:), 0.02) 

    !/ And start the sub-grid stuff over different landuse (iL)

    nlu = LandCover(i,j)%ncodes
    nFlux = 0           ! Numberof LC needing flux calcs
    LULOOP: do iiL= 1, nlu
        iL      = LandCover(i,j)%codes(iiL)

        iL_used (iiL) = iL    ! for eco dep

        if ( LandType(iL)%flux_wanted ) then
          nFlux = nFlux + 1
          iL_fluxes (nFlux) = iL    ! for eco dep
        end if

        f_phen   = LandCover(i,j)%fphen(iiL)

        Sub(iL)%SGS = LandCover(i,j)%SGS(iiL)   !used for AOT CHECK?
        Sub(iL)%EGS = LandCover(i,j)%EGS(iiL)

        L = Sub(iL)    ! ! Assign e.g. Sub(iL)ustar to ustar


             if ( DEBUG_DRYDEP .and. debug_flag ) then
                write(6,"(a,3i3,f6.1,2i4,3f7.3,i4,i2,2f6.2)") "DVEG: ", &
                    nlu,iiL, iL, glat(i,j), L%SGS, L%EGS, &
                   L%coverage, L%LAI, L%hveg,daynumber, &
                   Grid%sdepth, fSW(i,j),L%t2C !ACB Grid%snow_flag

                write(6,"(a,i4,2f7.2,2es10.2,3f8.3)") "DMET SUB", &
                  iL, Grid%ustar, L%ustar, Grid%invL, &
                  L%invL, L%Ra_ref, L%Ra_3m,L%rh

             end if


         call Rb_gas(L%is_water, L%ustar, L%z0, DRYDEP_GASES ,Rb)

         call Rsurface(i,j,DRYDEP_GASES ,Gns,Rsur,errmsg,debug_flag,snow_iL)

           Sub(iL)%g_sto = L%g_sto   ! needed elsewhere
           Sub(iL)%g_sun = L%g_sun

         Grid_snow(i,j) = Grid_snow(i,j) +  L%coverage * snow_iL 


       !/... add to grid-average Vg:


         wet =   Grid%wetarea
         dry =   1.0 - wet


         do n = 1, NDRYDEP_CALC

            if ( n > NDRYDEP_GASES )  then    ! particles

                !nae = n - NDRYDEP_GASES 
                nae = AERO_SIZE(n)


              if ( LandType(iL)%is_forest  ) then ! Vds NOV08

                 !/ Use eqn *loosely* derived from Petroff results

                  Vds = GPF_Vds300(L%ustar,L%invL, L%SAI )
                  if (n==CDDEP_PMfN .and. L%invL<0.0 ) then ! We allow nitrate to deposit x 2
                       Vds = Vds * 2.0 ! for nitrate-like
                  end if

              else !!!

                !/  Use Wesely et al  for other veg & sea

                 ! Vds = Nemitz2004( 0.4, L%ustar, L%invL )
                 Vds = Wesely300( L%ustar, L%invL )
                  ! We allow nitrate to deposit x 2
                  if (n==CDDEP_PMfN .and. L%invL<0.0 ) then 
                       Vds = Vds * 2.0 ! for nitrate-like
                  end if

              end if

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
               print "(a,2i4,3es10.3)", "VDS CHECK ",n, nae, Vds, Vg_ref(n)
            
            call CheckStop((Vg_3m(n)>0.50 .or. Vg_ref(n)>0.50 ), "AEROSTOP")
          end if
        end if

         ! ================================================


            else                           ! gases

           ! NB no wet-dry difference needed here
 
              Vg_ref(n) = 1. / ( L%Ra_ref + Rb(n) + Rsur(n) ) 

              Vg_3m (n) = 1. / ( L%Ra_3m + Rb(n) + Rsur(n) ) 


            endif

         ! Surrogate for NO2 compensation point approach, 
         ! assuming c.p.=4 ppb (ca. 1.0e11 #/cm3):        
         ! Note, xn_2d has no2 in #/cm-3

          if ( n == CDDEP_NO2 ) then

            Vg_ref(CDDEP_NO2) = Vg_ref(CDDEP_NO2) * no2fac
            Vg_3m (CDDEP_NO2) = Vg_3m (CDDEP_NO2) * no2fac
          end if ! CDDEP_NO2

          !SUB0   Grid%Vg_ref(n) = Grid%Vg_ref(n) + L%coverage * Vg_ref(n)
          !SUB0   Grid%Vg_3m(n)  = Grid%Vg_3m(n)  + L%coverage * Vg_3m(n)
            Sub(0)%Vg_ref(n) = Sub(0)%Vg_ref(n) + L%coverage * Vg_ref(n)
            Sub(0)%Vg_3m(n)  = Sub(0)%Vg_3m(n)  + L%coverage * Vg_3m(n)
            Sub(iL)%Vg_ref(n) = Vg_ref(n)
            Sub(iL)%Vg_3m(n) = Vg_3m(n)


         if ( n <= NDRYDEP_GASES )  then    ! gases
             !QUERY - do we need Gsur for anything now?!

             Sub(iL)%Gsur(n) = 1.0/Rsur(n) ! Note iL, not iiL 
             Sub(iL)%Gns(n)  = Gns(n)  ! Note iL, not iiL 

             !SUB0 Grid%Gsur(n)  =  Grid%Gsur(n) + L%coverage / Rsur(n)
             !SUB0 Grid%Gns(n)  =  Grid%Gns(n)+ L%coverage * Gns(n)
             Sub(0)%Gsur(n)  =  Sub(0)%Gsur(n) + L%coverage / Rsur(n)
             Sub(0)%Gns(n)   =  Sub(0)%Gns(n)  + L%coverage * Gns(n)
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
            do n = CDDEP_O3 , CDDEP_O3 !!! 1,NDRYDEP_GASES 
               call datewrite("DEPO3 ", iL, &
                   (/ Vg_ref(n), Sub(iL)%Vg_ref(n) /) )
                   !(/ Mosaic_VgRef(n,iL) , Vg_ref(n), Sub(iL)%Vg_ref(n) /) )
               call datewrite("DEPDVG", iL, (/ L%coverage, 1.0*n,& ! gs in cm/s :
                 L%LAI,100.0*L%g_sto, L%Ra_ref, Rb(n), min( 999.0,Rsur(n) ),  &
                100.0*Vg_3m(n), 100.0*Vg_ref(n), Vg_ratio(n) /) )
            end do

        end if

       !
         
       !=======================

        if (  LandType(iL)%flux_wanted ) then

            n = CDDEP_O3

           Ra_diff = L%Ra_ref - L%Ra_3m   
           c_hveg3m = xn_2d(FLUX_TOT,KMAX_MID)  &     ! #/cm3 units
                        * ( 1.0-Ra_diff*Vg_ref(n) )

         ! Flux = Vg_ref*c_ref = Vg_h * c_h = (c_ref-c_h)/Ra(z_ref,z_h)
         ! which gives:  c_h = c_ref * [ 1-Ra(z_ref,z_h)*Vg_ref ]
         ! Resistance Ra from z_ref to top of canopy:

           Ra_diff  = AerRes(max( L%hveg-L%d, STUBBLE) , Grid%z_ref-L%d,&
                       L%ustar,L%invL,KARMAN)

           c_hveg = xn_2d(FLUX_TOT,KMAX_MID)  &     ! #/cm3 units
                        * ( 1.0-Ra_diff*Vg_ref(n) )

          if ( DEBUG_AOT .and. debug_flag ) then
              write(*, "(a,3i3,i5,i3, 3f9.3,f5.2,9es10.3)") &
               "CHVEG ", imm, idd, ihh, current_date%seconds,  iL, &
                  xn_2d(FLUX_TOT,KMAX_MID)*surf_ppb, c_hveg*surf_ppb,&
                   c_hveg3m * surf_ppb, &
                   100.0*Vg_ref(n), L%Ra_ref, (L%Ra_ref-L%Ra_3m), Ra_diff, Rb(n),Rsur(n)
          end if
           

         ! Need to be careful with scope. L is within iL loop, whereas Sub
         ! will be kept throughout i,j calculations:

          Sub(iL)%cano3_ppb   = c_hveg * surf_ppb  ! change units
          Sub(iL)%cano3_nmole = c_hveg * NMOLE_M3  ! units of nmole/m3

        end if !

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! Extra outputs sometime used for Sweden/IVL/SEI/CEH
          !! Uncomment and make .inc file as required

!            include 'EXTRA_LU_Outputs.inc'

       !=======================
        end do LULOOP
       !=======================
       !=======================

         call Calc_StoFlux(nFlux, iL_fluxes(1:nFlux), debug_flag )


        if ( DEBUG_DRYDEP .and. Sumland > 1.011  ) then
            call CheckStop( "DryDep:SUMLAND TOO BIG")
        end if


        if ( Sumland > 0.01 ) then
            gradient_fac(:) = Vg_ratio(:) / Sumland
        else
            gradient_fac(:) = sea_ratio(:)
        end if

        if ( DEBUG_DRYDEP .and. debug_flag ) then
            call datewrite("DEP VGR snow_flag Vg", (/ Grid%sdepth, & !ACB Grid%snow_flag
             !SUB0 (/ (100.0*Grid%Vg_Ref(n), n = 1, min(4,NDRYDEP_GASES )) , &
                (100.0*Sub(0)%Vg_Ref(n), n = 1, min(4,NDRYDEP_GASES )) , &
                (100.0*Sub(iL)%Vg_Ref(n), n = 1, min(4,NDRYDEP_GASES )) /) )
        end if


!-- loop through all affected advected species to calculate changes in
!   concentration (xn_adv), the conc. ratios (cfac), and deposition 

    do ncalc = 1, NDRYDEP_CALC

        !SUB0 vg_fac (ncalc) = 1.0 - exp ( -Grid%Vg_Ref(ncalc) * dtz ) 
        vg_fac (ncalc) = 1.0 - exp ( -Sub(0)%Vg_Ref(ncalc) * dtz ) 

    end do ! n


    GASLOOP2 :  do n = 1, NDRYDEP_ADV 
         nadv    = DDepMap(n)%ind
         ntot  = NSPEC_SHL + DDepMap(n)%ind
         ncalc   = DDepMap(n)%calc


         if ( vg_set(n) ) then

             DepLoss(nadv) =   & ! Use directly set Vg
                 ( 1.0 - exp ( -DDepMap(n)%vg * dtz ) ) * xn_2d( ntot,KMAX_MID)
             cfac(nadv, i,j) = 1.0   ! Crude, for now.
  
         else
           if ( ntot >= FIRST_SEMIVOL .and. ntot <= LAST_SEMIVOL ) THEN
              ! Assuming dry deposition of particulate part of
              ! semi-volatile components as PMfS and the gaseous part as
              ! specified in GenIn.species.

             DepLoss(nadv) =  &
              Fgas(ntot,KMAX_MID)*vg_fac( ncalc ) * xn_2d(ntot,KMAX_MID) + &
              Fpart(ntot,KMAX_MID)*vg_fac( CDDEP_PMfS ) * xn_2d(ntot,KMAX_MID)

              cfac(nadv, i,j) = Fgas(ntot,KMAX_MID)*gradient_fac(ncalc) + &
                   Fpart(ntot,KMAX_MID)*gradient_fac( CDDEP_PMfS )
            else
               DepLoss(nadv) =   vg_fac( ncalc )  * xn_2d( ntot,KMAX_MID)
               cfac(nadv, i,j) = gradient_fac( ncalc )
            endif
         end if

         if ( DepLoss(nadv) < 0.0 .or. &
              DepLoss(nadv)>xn_2d(ntot,KMAX_MID) ) then
             call CheckStop("NEGXN DEPLOSS" )
         end if


         if ( ntot == O3 ) then 

           Grid%surf_o3_ppb = xn_2d(O3,KMAX_MID)*gradient_fac( ncalc )*surf_ppb

           Grid%O3factor = vg_fac(ncalc)

           if ( DEBUG_DRYDEP .and. debug_flag ) then
              call datewrite("O3_ppb_ratios ", n, (/ &
                 Grid%surf_o3_ppb, Grid%O3factor /) )
           end if ! DEBUG
         end if

        xn_2d( ntot,KMAX_MID) = &
             xn_2d( ntot,KMAX_MID) - DepLoss(nadv)



        if ( ntot == FLUX_TOT ) then
           ! fraction by which xn is reduced - safety measure:
             
              if( xn_2d( ntot,KMAX_MID)  > 1.0e-30 ) then
                  lossfrac = ( 1.0 - DepLoss(nadv)/ &
                                (DepLoss(nadv)+xn_2d( ntot,KMAX_MID)))
              end if
              if ( DEBUG_DRYDEP .and. lossfrac < 0.1 ) then
                  call datewrite( "LOSSFRACING ", nadv, (/ 1.0*iL, &
                  !SUB0 Grid%Vg_Ref(n), DepLoss(nadv), vg_fac(ncalc), lossfrac /) )
                    Sub(0)%Vg_Ref(n), DepLoss(nadv), vg_fac(ncalc), lossfrac /) )
                  call CheckStop( lossfrac < 0.1, "ERROR: LOSSFRAC " )
              end if
        end if


      !.. ecosystem specific deposition - translate from calc to adv 
      !  and normalise

         IIL_LOOP : do iiL = 1, nlu
            iL      = iL_used(iiL)

            if ( vg_set(n) )  then
               fluxfrac_adv(nadv,iL) = Sub(iL)%coverage  ! Since all vg_set equal
            else
               !Vg_scale = Mosaic_VgRef(ncalc,iL)/ Mosaic_VgRef(ncalc,0)
               Vg_scale = Sub(iL)%Vg_Ref(ncalc)/ Sub(0)%Vg_Ref(ncalc)
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
            if ( iL == 1 ) then ! SO2, CF

               if ( vg_set(n) )  then
                 write(6,"(a,3i3,3f12.3)") "FLUXSET  ", iiL, iL, nadv, &
                     100.0*DDepMap(n)%vg, Sub(iL)%coverage, fluxfrac_adv(nadv,iL)
               else
                 write(6,"(a,3i3,f8.5,5f8.3)") "FLUXFRAC ", iiL, iL, nadv, &
                  Sub(iL)%coverage, &
                  100.0*Sub(0)%Vg_Ref(ncalc), &  ! GRID
                  100.0*Sub(iL)%Vg_Ref(ncalc), & ! Mosaic VgRef &
                  100.0*Sub(iL)%coverage*Sub(iL)%Vg_Ref(ncalc), & 
                   fluxfrac_adv(nadv,iL)
               end if
            end if !SO2 CF
            end if
         end do IIL_LOOP
          

      !..accumulated dry deposition per grid square and summed over the whole
      !  domain

!         totddep( nadv ) = totddep (nadv) + Deploss(nadv) * convfac


        if ( DEBUG_DRYDEP .and. debug_flag ) then
          if ( vg_set(n) ) then
              write(*, "(a,2i4,f8.3)") "DEBUG DryDep SET ", &
                   n,nadv, DDepMap(n)%vg
          else
              if( n == 1 ) & ! O3
              call datewrite( "DEBUG DDEPxnd: "// trim(species(ntot)%name), &
                n, (/ real(nadv),real(ncalc), gradient_fac( ncalc),&
                   xn_2d(ntot,KMAX_MID), vg_fac(ncalc) /) )
          end if
        end if

          if ( DEBUG_AOT .and. debug_flag .and. ntot == FLUX_TOT  ) then
              write(*, "(a,3i3,i5,i3,2f9.4,f7.3)") &
               "AOTCHXN ", imm, idd, ihh, current_date%seconds, &
                   iL, xn_2d(FLUX_TOT,KMAX_MID)*surf_ppb, &
                    (xn_2d( FLUX_TOT,KMAX_MID) + DepLoss(nadv) )*surf_ppb, &
                     gradient_fac( ncalc)
          end if
       end do GASLOOP2 ! n

    !  DryDep Budget terms

      convfac =  convfac/amk(KMAX_MID)
      do n = 1, NDRYDEP_ADV
         nadv    = DDepMap(n)%ind
         totddep( nadv ) = totddep (nadv) + DepLoss(nadv)*convfac
      enddo

       convfac2 = convfac * xm2(i,j) * inv_gridarea


      !.. Add DepLoss to budgets if needed:

       call Add_MosaicOutput(debug_flag,i,j,convfac2,&
           DepAdv2Calc, fluxfrac_adv, Deploss ) 

 end subroutine drydep

end module DryDep_ml
