!> <PBAP_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!  **************************************************************************! 

module PBAP_mod

  !/-- Module to deal with primary biological aerosol pollutants (PBAPs).
  !
  !    Gunnar Felix Lange 2023
  !---------------------------------------------------------------------------

  use CheckStop_mod,      only: CheckStop, StopAll
  use ChemSpecs_mod,         only : species
  use Config_module, only : NPROC, MasterProc, TINY, &
                           NLANDUSEMAX, IOU_INST, & 
                           KT => KCHEMTOP, KG => KMAX_MID, & 
                           EURO_SOILNOX_DEPSCALE, & 
                           MasterProc, &
                           USES, &
                           NATBIO, EmBio, EMEP_EuroBVOCFile
  use Debug_module,       only: DebugCell, DEBUG
  use GridValues_mod,     only: i_fdom,j_fdom, debug_proc,debug_li,debug_lj
  use Io_mod,             only: IO_FORES, open_file, ios, datewrite
  use Io_RunLog_mod,      only: PrintLog
  use KeyValueTypes,      only: KeyVal,KeyValue
  use LandDefs_mod,       only: LandType, LandDefs
  use LandPFT_mod,        only: MapPFT_LAI, pft_lai
  use Landuse_mod,        only: LandCover
  use LocalVariables_mod, only: Grid  ! -> izen, DeltaZ
  use MetFields_mod,      only: t2_nwp, q
  use MetFields_mod,      only: PARdbh, PARdif !WN17, in W/m2
  use NetCDF_mod,         only: ReadField_CDF, printCDF
  use OwnDataTypes_mod,   only: Deriv, TXTLEN_SHORT
!  use Paleo_mod, only : PALEO_modai, PALEO_miso, PALEO_mmon
  use Par_mod,            only: MSG_READ1,me, limax, ljmax
  use PhysicalConstants_mod,  only:  AVOG, GRAV, PI
  use Radiation_mod,      only: PARfrac, Wm2_uE
  use SmallUtils_mod,     only: find_index
  use TimeDate_mod,       only: current_date, daynumber
  use ZchemData_mod,      only: rcemis, rcbio
  implicit none
  private

  !/-- subroutines for PBAPs
  public ::  Init_PBAPs
  !/-- subroutines for Fungal Spores
  private :: Set_FungalSpores


  integer, public, parameter ::   NPBAP = 1 !Number of PBAPs currently implemented
                                            !(only fungal spores at the moment)

  real*8, DIMENSION(3), parameter  ::  &
  FUNG_PARAMS = [20.426d0, 275.82d0, &
         39300.0d0] !From Fungal paramterization, Eq. (2) of
                    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                    !DOI 10.1007/978-3-319-35095-0_121

  ! - allows rcbio in CM_Reactions, but we access elements with
  ! the natbio indices here. These much match the indices used in rcbio
  ! We only use rcbio for isoprene and terpenes so far,  since
  ! soil NO, NH3 emissions etc are dealt with through rcemis.

  ! We hard-code these indices, but only calculate emissions if needed
  ! Must match order of NATBIO to start with 
  character(len=16), save, dimension(NPBAP), public:: &
      EMIS_PBAP = [character(len=16):: &
             "FUNGAL_SPORES"]

  integer, public, parameter :: &
      N_ECF=2, ECF_ISOP=1, ECF_TERP=2   &! canopy factors, BVOC
     ,BIO_ISOP=1, BIO_MTP=2, BIO_MTL=3  &! BIO_SOILNO=4, BIO_SOILNH3=5
     ,BIO_TERP=2 ! Used for final emis, sum of MTP+MTL
  integer, public, save ::  last_bvoc_LC   !max index land-cover with BVOC (min 4)
                                                        
  ! Soil NOx
   real,public, save, allocatable, dimension(:,:) :: &
      AnnualNdep, &  ! N-dep in mgN/m2/
      SoilNOx, SoilNH3
   real,public, save, allocatable, dimension(:,:,:) :: SoilNOx3d 

 ! Set true if LCC read from e.g. EMEP_EuroBVOC.nc:
 ! (Currently for 1st four LCC, CF, DF, BF, NF)
  logical, private, dimension(NLANDUSEMAX), save :: HaveLocalEF 

  ! EmisNat is used for BVOC; soil-NO, also in futur for sea-salt etc.
  ! Main criteria is not provided in gridded data-bases, often land-use
  ! dependent.

  real, public, save, allocatable, dimension(:,:,:) :: &
     EmisNat       !  will be transferred to d_2d emis sums


  !standard emission factors (EFs) per LC
  !Need to dimension later for Emtp, Emtl, last_bvoc_LC
  real, private, save, allocatable, dimension(:,:,:,:) :: &
     bvocEF       !  Gridded std. emissions per PFT

  !standard emission factors per LC for daily LAI
  real, private, save, allocatable, dimension(:,:,:) :: &
     day_embvoc    !  emissions scaled by daily LAI

  logical, private, save, allocatable, dimension(:,:) :: EuroMask

  !/-- Canopy environmental correction factors-----------------------------
  !
  !    - to correct for temperature and light of the canopy 
  !    - from Guenther's papers. (Limit to 0<T<40 deg C.)

  real, public, save, dimension(N_ECF,40) :: canopy_ecf  ! Canopy env. factors

 ! Indices for the species defined in this routine. Only set if found
  integer, private, save :: itot_C5H8,  itot_TERP,  itot_NO , itot_NH3

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine Init_PBAPs()
    allocate(FSpores(LIMAX,LJMAX)) !Spatial distribution of fungal spores
  
  end subroutine Init_PBAPs

  

  subroutine Set_FungalSpores()
    !!!!!!!
    !Fills FungalSpores(i,j)
    !Based on parameterization from
    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
    !DOI 10.1007/978-3-319-35095-0_121
    integer :: i, j, nLC, iLC, LC
    logical :: my_first_call = .true.
    real    :: f, ft, fn, ftn
    real    :: eFGS !, emissions Fungal spores TODO: Fix units
    real :: beta, bmin, bmax, bx, by ! for beta function
    real :: hfac


    if ( .not. USES%FUNGAL_SPORES) return ! TODO: Add the spores to this module

    if( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
       write(*,*)"PBAP_mod DEBUG FUNGAL_SPORES: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%FUNGAL_SPORES
    end if


      do j = 1, ljmax
         do i =  1, limax

           nlc = LandCover(i,j)%ncodes

           F_FNG =  FUNG_PARAMS(1)*(t2_nwp(i,j,1)-FUNG_PARAMS(2))+FUNG_PARAMS(3)*q(i,j,1)*LAIfac
           !Fungal spores flux, Eq.(2) of S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)


   LCLOOP: do ilc= 1, nLC

               LC = LandCover(i,j)%codes(ilc)
               if ( LandType(LC)%is_water ) cycle
               if ( LandType(LC)%is_ice   ) cycle
               if ( LandType(LC)%is_iam   ) cycle

             ! Soil NO
             ! for 1 ugN/m2/hr, the temp funct alone would give
             ! ca. 6 mgN/m2/a in Germany, where dep is about 5000 mgN/m2 max ca. 9
             ! Conif Forests in Germany  

               f  = LandCover(i,j)%fraction(ilc) 
               beta = 0.0

               if ( LandType(LC)%is_conif ) then
                  enox = enox + f*ftn*150.0
                  !enh3 = enh3 + f*ftn*1500.0 ! Huge?! W+E ca. 600 ngNH3/m2/s -> 1800 ugN/m2/h
               else if ( LandType(LC)%is_decid ) then
                  enox = enox + f*ftn* 50.0
                  !enh3 = enh3 + f*ftn*500.0 !  Just guessing
               else if ( LandType(LC)%is_seminat ) then
                  enox = enox + f*ftn* 50.0
                  !enh3 = enh3 + f * ftn  *20.0 !mg/m2/h approx from US report 1 ng/m2/s

               else if ( LandType(LC)%is_crop    ) then ! emissions in 1st 70 days

                  bmin = Landcover(i,j)%SGS(iLC) -30 ! !st March LandCover(i,j)%SGS(iLC) - 30 
                  bmax = Landcover(i,j)%SGS(iLC) +30 ! End April  LandCover(i,j)%SGS(iLC) + 40 

                  ! US p.29, Suttn had ca. 20 ng/m2/s  = 60ugN/m2/hfor crops
                  ! throughout growing season
                  if ( daynumber >= Landcover(i,j)%SGS(iLC) .and. &
                       daynumber <= Landcover(i,j)%EGS(iLC) ) then
                       enox = enox + f* 1.0
                       !enh3 = enh3 + f * 60.0
                  end if

                  ! CRUDE - just playing for NH3.
                  ! NH3 from fertilizer? Assume e.g. 120 kg/ha over 1 month
                  ! with 10% giving emission, i.e. 10 kg/ha
                  ! 10 kg/ha/month =  ca. 1000 ugN/m2/h

                  ! For NO, numbers based upon papers by e.g. Rolland,
                  ! Butterbach, etc.
                  if ( daynumber >= bmin .and. daynumber <= bmax ) then

                       bx = (daynumber-bmin)/( bmax-bmin)
                       bx = max(bx,0.0)
                       by = 1.0 - bx
                       beta =  ( bx*by *4.0) 
                       enox = enox + f*80.0*ft* beta 
                       !enh3 = enh3 + f * 1000.0*ft * beta
                  end if

                  
               end if
               if (  DEBUG%SOILNOX .and. debug_proc .and. &
                   i == debug_li .and. j == debug_lj ) then
                 write(*, "(a,4i4,f7.2,9g12.3)") "LOOPING SOIL", daynumber, &
                 iLC, LC, LandCover(i,j)%SGS(iLC), t2_nwp(i,j,1)-273.15, &
                    f, ft, fn, ftn,  beta, enox!, enh3
                 if(iLC==1) &
                   call datewrite("HFAC SOIL", (/ 1.0*daynumber,hfac /) )
               end if
               enox = max( 0.001, enox ) ! Just to stop negatives while testing
  
             ! Soil NH3
         end do LCLOOP


   ! And we scale EmisNat to get units kg/m2 consistent with
   ! Emissions_mod (snapemis).  ug/m2/h -> kg/m2/s needs 1.0-9/3600.0. 

         SoilNOx(i,j) = enox

           !enh3 = 0.0 ! BIDIR SOON .... we don't want enh3
           !SoilNH3(i,j) = enh3

       end do
    end do

    if ( DEBUG%SOILNOX .and. debug_proc ) then
       i = debug_li
       j = debug_lj
       write(*,"(a,4i4)") "RESET_SOILNOX: ",  1, limax, 1, ljmax
       write(*,"(a,2i4,2f12.4,es12.4)") "RESET_SOILNOX: ", &
               daynumber, current_date%hour, t2_nwp(i,j,1), SoilNOx(i,j), AnnualNdep(i,j)
    end if

    my_first_call = .false.

  end subroutine Set_FungalSpores

  subroutine setup_PBAPs(i,j)
  !
  !---- assign fungal rates  ------------------------------------------------
  !
  !  So far, assigns isoprene using surface (2m) temperature, and for all
  !  zenith angles <90. Should include light dependance at some stage
  !
  !  Output : rcbio added to rcemis - isoprene emissions for 1d column
  !
  !  Called from setup_mod, every  advection step.
  !----------------------------------------------------------------------------

  integer, intent(in) ::  i,j

  character(len=*), parameter :: dtxt='BioModSetup:' 
  integer :: it2m, gmt_3hour
  real    :: E_ISOP, E_MTP, E_MTL, sFac

! To get from ug/m2/h to molec/cm3/s
! ug -> g  1.0e-6; m2-> cm2 1e-4, g -> mole / MW; x AVOG
! will use /Grid%DeltaZ, which is in m, so anoter 1e-2  tp et cm-3
  real, parameter :: & 
        biofac_ISOP   = 1.0e-12*AVOG/68.0 /3600.0  &
       ,biofac_TERP   = 1.0e-12*AVOG/136.0/3600.0  &
       ,biofac_SOILNO = 1.0e-12*AVOG/14.0 /3600.0  &
       ,biofac_SOILNH3= 1.0e-12*AVOG/14.0 /3600.0  
  logical :: dbg

 ! Light effects added for isoprene emissions

  real            :: par   ! Photosynthetically active radiation
  real            :: cL    ! Factor for light effects
  real, parameter :: &
      CL1 = 1.066  , &    ! Guenther et al's G93/G95 params
      ALPHA = 0.0027!,&   ! Guenther et al's G93/G95 params
!     AG99  = 0.001 * 1.42  ! " Warneke update, but not same as G99?

  if ( size(BVOC_USED) == 0  ) return   ! e.g. for ACID only

  dbg = ( DEBUG%BIO .and. debug_proc .and. &
          i==debug_li .and. j==debug_lj .and. current_date%seconds == 0 )

  !TLEAF: prepping for future Tleaf calculation.
  ! We have a mix of forest species for any emitting PFT, so we just
  ! use the one Grid based value, which is appropriate for trees.
  ! Oct 2021: Grid%dTleaf is just zero
  it2m = nint( Grid%t2C - TINY )
  if ( dbg ) write(*,*)'DBGITA',  Grid%t2C, it2m, canopy_ecf(BIO_ISOP,it2m)
  it2m = nint( Grid%t2C + Grid%dTleaf  - TINY )
  it2m = max(it2m,1)
  it2m = min(it2m,40)

  !ASSUME C5H8 FOR NOW if ( ibn_C5H8 > 0 ) then
    if ( Grid%izen <= 90) then ! Isoprene in daytime only:

     ! Light effects from Guenther G93. Need uE:

      par = ( PARdbh(i,j) + PARdif(i,j)  ) * Wm2_uE

      cL = ALPHA * CL1 * par/ sqrt( 1 + ALPHA*ALPHA * par*par)

     ! E in ug/m2/h

       E_ISOP = day_embvoc(i,j,BIO_ISOP)*canopy_ecf(BIO_ISOP,it2m) * cL &
                  * EmBio%IsopFac
       if ( dbg ) write(*,*)'DBGITB',  Grid%dTleaf, it2m, canopy_ecf(BIO_ISOP,it2m), cL, E_ISOP

      ! Add light-dependent terpenes to pool-only
      if(BIO_TERP > 0) E_MTL = &
             day_embvoc(i,j,BIO_MTL)*canopy_ecf(ECF_TERP,it2m)*cL * EmBio%TerpFac

     !  molecules/cm3/s
     ! And we scale EmisNat to get units kg/m2 consistent with
     ! Emissions_mod (snapemis).  ug/m2/h -> kg/m2/s needs 1.0-9/3600.0. 


      rcbio(NATBIO%C5H8,KG)   = E_ISOP * biofac_ISOP/Grid%DeltaZ
      EmisNat(NATBIO%C5H8,i,j)= E_ISOP * 1.0e-9/3600.0

  else ! night
     rcbio(NATBIO%C5H8,KG)    = 0.0
     EmisNat(NATBIO%C5H8,i,j) = 0.0
     E_MTL = 0.0
     E_ISOP = 0.0
     par = 0.0   ! just for printout
     cL  = 0.0   ! just for printout
  end if ! daytime

 ! add pool-only terpenes rate;
  E_MTP = day_embvoc(i,j,BIO_MTP)*canopy_ecf(ECF_TERP,it2m) * EmBio%TerpFac
  rcbio(NATBIO%TERP,KG)    = (E_MTL+E_MTP) * biofac_TERP/Grid%DeltaZ
  EmisNat(NATBIO%TERP,i,j) = (E_MTL+E_MTP) * 1.0e-9/3600.0

  if ( USES%SOILNOX ) then
    if ( USES%SOILNOX_METHOD == 'OLD_EURO' ) then
      rcemis(itot_NO,KG)    = rcemis(itot_NO,KG) + &
           SoilNOx(i,j) * biofac_SOILNO/Grid%DeltaZ
      EmisNat(NATBIO%NO,i,j) =  SoilNOx(i,j) * 1.0e-9/3600.0

    else !GLOBAL emissions should be in molecules/m2/s (NB: not molecules/cm3/s!)
      EmisNat(NATBIO%NO,i,j) =  SoilNOx(i,j)/biofac_SOILNO * 1.0e-15/3600.0 !molecules/m2/s -> kg/m2/h ?
      gmt_3hour = 1 + int(current_date%hour/3)
      !molecules/m2/s -> molecules/cm3/s:
      rcemis(itot_NO,KG)    = rcemis(itot_NO,KG) + &
         SoilNOx3D(i,j,gmt_3hour)/Grid%DeltaZ * 1.0e-6

    end if ! USES%SOILNOX_METHOD
  end if ! USES%SOILNOX

    !EXPERIMENTAL
    !if ( USES%SOILNH3 ) then
    if ( USES%BIDIR ) then
       !FUTURE? rcbio(NATBIO%NH3,KG)    = &
       rcemis(itot_NH3,KG)    = rcemis(itot_NH3,KG) + &
           SoilNH3(i,j) * biofac_SOILNH3/Grid%DeltaZ
        if(NATBIO%NH3>0)EmisNat(NATBIO%NH3,i,j) =  SoilNH3(i,j) * 1.0e-9/3600.0
    else
        if(NATBIO%NH3>0)EmisNat(NATBIO%NH3,i,j) = 0.0
    end if
     
 
    if ( dbg .and. current_date%seconds==0 ) then 

      !molecules/m2/s -> molecules/cm3/s
      sFac = 1/Grid%DeltaZ * 1.0e-6

      call datewrite(dtxt//" env ", it2m, (/ max(par,0.0), max(cL,0.0), &
            canopy_ecf(BIO_ISOP,it2m),canopy_ecf(BIO_TERP,it2m) /) )
      write(*,*) dtxt//" EISOP RAW ",  gmt_3hour, E_ISOP
      call datewrite(dtxt//" E_SOI ", [  gmt_3hour ], [ SoilNOx(i,j) ] )
      call datewrite(dtxt//" EISOP EMTP EMTL ESOIL-N ", [  gmt_3hour ], &
       [ E_ISOP, E_MTP, E_MTL, SoilNOx(i,j) * sFac, &
        SoilNOx3D(i,j,gmt_3hour) * sFac ] ) 

      if (USES%BIDIR) call datewrite(dtxt//" BIDIR ", (/  SoilNOx(i,j) * sFac,&
                                      SoilNH3(i,j), rcbio(NATBIO%NH3,KG) /) ) 
      call datewrite(dtxt//" rcemisL ", (/ Grid%t2C , Grid%dTleaf, &
            rcbio(NATBIO%C5H8,KG), rcbio(NATBIO%TERP,KG) /))
      call datewrite(dtxt//" EmisNat1:8 ", EmisNat(1:8,i,j) )

     end if


  end subroutine setup_PBAPs

  !----------------------------------------------------------------------------
end module PBAP_mod
