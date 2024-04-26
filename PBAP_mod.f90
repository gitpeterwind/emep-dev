!> <PBAP_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!  **************************************************************************!

module PBAP_mod

  !/-- Module to deal with primary biological aerosol pollutants (PBAPs).
  !    Currently implemented: Fungal spores and bacteria
  !    TODO: Solubility of Fungal Spores
  !
  !    Gunnar Felix Lange 2024
  !---------------------------------------------------------------------------

  use CheckStop_mod,      only: CheckStop, StopAll
  use ChemSpecs_mod,         only : species
  use Config_module, only : NPROC, MasterProc,&
                           KT => KCHEMTOP, KG => KMAX_MID,&
                           USES, &
                           NATBIO, EmBio,OceanChlorophyll_File
  use Debug_module,       only: DebugCell, DEBUG
  use GridValues_mod,     only: debug_proc,debug_li,debug_lj
  use LandDefs_mod,       only: LandType, LandDefs
  use Landuse_mod,        only: LandCover, likely_coastal
  use LocalVariables_mod, only: Grid  ! -> izen, DeltaZ
  use MetFields_mod,      only: t2_nwp, q,ustar_nwp
  use SeaSalt_mod,        only: SeaSalt_flux
  use Par_mod,            only: MSG_READ1,me, limax, ljmax
  use PhysicalConstants_mod,  only:  AVOG, GRAV, PI
  use SmallUtils_mod,     only: find_index
  use TimeDate_mod,       only: current_date, daynumber
  use ZchemData_mod,      only: rcemis, rcbio
  use Biogenics_mod,      only: NEMIS_BioNat,EMIS_BioNat,EmisNat

  implicit none
  private

  !/-- subroutines for PBAPs
  public ::  init_PBAPs,set_PBAPs
  !/-- subroutines for Fungal Spores
  private :: Set_FungalSpores, Set_Bacteria, Set_MarineOA

  integer, public, save ::   NPBAP !Number of PBAPs!(only fungal spores bacteria and Marine OA implemented at the moment)

  real,public, save, allocatable, dimension(:,:,:) :: PBAP_flux !Dim: i,j,NPBAP

  real,public, save, allocatable, dimension(:,:) :: O_Chlorophyll !Amount of Chlorophyll in the ocean. Dim: i,j
                                                                  !Used for MarineOA


  real,public, save, allocatable, dimension(:) :: WEIGHTS ,DIAMETERS, DENSITIES !Physical parameters. Dim: NPBAP
  real,private, save, allocatable, dimension(:) :: n2m, kgm2h !Conversion factors. Dim: NPBAP
  real,private, save, allocatable, dimension(:) :: itot, inat !indices in species and EMIS_Bionat. Dim: NPBAP
  real, public :: FUNGAL_DIAMETER, FUNGAL_WEIGHT !Diameter depends on chosen fungal parameterization scheme

  character(len=20),private,save,allocatable,dimension(:) ::PBAP_names !For debugging only


  integer, private, save :: iint_FungalSpores, itot_FungalSpores, inat_FungalSpores !Index of fungal spores internall, in Species and in EMIS_BioNat
  integer, private, save :: iint_Bacteria,itot_Bacteria, inat_Bacteria !Index of bacteria spores internally, in Species and in EMIS_BioNat
  integer, private, save :: iint_MarineOA,itot_MarineOA, inat_MarineOA !Index of Marine OA internally, in Species and in EMIS_BioNat



  !Choice of fungal flux parameterization:
  !HM: Hummel(2015) DOI: 10.5194/acp-15-6127-2015
  !SD: Sesartic and Dallafior (2011) DOI: 10.5194/bg-8-1181-2011
  !HS: Heald and Spracklen (2009) DOI:10.1029/2009GL037493
  !HO: Hoose et al (2010) DOI: 10.1088/1748-9326/5/2/024009, based on HS
  !JS: Janssen et al. (2021) DOI: 10.5194/acp-21-4381-2021

  !(GFL Feb 2024): Hummel parameterization uses different settling scheme than EMEP, Sesartic and Dallfior uses
  !no settling scheme at all, so currently seems that Heald and Spracken parameterization modified by Hoose
  !(parameterization choice = "HO") gives best results for fungal spores

  real*8, DIMENSION(6), parameter  ::  &
  BACTERIA_PARAMS = [900.0,704.0,648.0,7.7,502.0,196.0] !From bacteria paramterization, Eq. (1) of
                    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                    !DOI 10.1007/978-3-319-35095-0_121
                    !Note that most of these are set in the LandInput file, except for the coastal parameter (as coastal is not a LandType)
  real*8, DIMENSION(3), parameter  ::  &
  FUNG_PARAMS_HM = [20.426, 275.82, 39300.0] !From Fungal paramterization, Eq. (2) of
  !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
  !DOI 10.1007/978-3-319-35095-0_121, based on Hummel et al. Atmos. Chem. Phys., 15, 6127–6146,
                                        !       https://doi.org/10.5194/acp-15-6127-2015, 2015

  !real*8, DIMENSION(3), parameter  ::  &
  !FUNG_PARAMS_HS_LARGE = [500.0, 5.0,0.015] !From Fungal parameterization, Heald and Spracken for 5um
                                       !fungal spores. DOI:10.1029/2009GL037493

  real*8, DIMENSION(3)  ::  &
  FUNG_PARAMS_HS = [2315.0, 5.0,0.015] !From Fungal parameterization, Heald and Spracken for 3um
                                       !modifiend in Hoose et al. (2010),DOI:


  real*8, DIMENSION(4)  ::  &
  FUNG_PARAMS_JS = [2.63*1.0e-5, 6.10*1.0e3,46.7,59.0] !From Fungal parameterization, Janssen for North America.
                                       !https://doi.org/10.5194/acp-21-4381-2021 (2021) [Has not been tested!]



  real, parameter :: FUNGAL_DENS = 1.0e6 !Fungal density [g/m3]
                                         !From Hummel et al. Atmos. Chem. Phys., 15, 6127–6146,
                                         !https://doi.org/10.5194/acp-15-6127-2015, 2015


  real, parameter :: BACTERIA_DIAMETER = 1.0 !Bacteria diameter [um]
                                             !From S. M. Burrows et al. Atmos. Chem. Phys., 9, 9281–9297
                                             !https://doi.org/10.5194/acp-9-9281-2009
  real, parameter :: BACTERIA_WEIGHT = 0.52*1.0e-12 !Bacterial weight [g]  (ibid)
  real, parameter :: BACTERIA_DENS   = BACTERIA_WEIGHT/((4/3.0)*PI*(0.5*BACTERIA_DIAMETER*1e-6)**3) !Bacteria density [g/m3]


  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine init_PBAPs()
    logical, save ::  my_first_call = .true.

    !Initializes the PBAP calculation
    !sets up number of PBAPs (NPBAP)
    !And fills array accordingly.
    !GFL Jan 2024: Is this parallelization safe?

    if (my_first_call) then
      NPBAP = 0
      iint_FungalSpores = -1
      itot_FungalSpores = -1
      iint_Bacteria = -1
      itot_Bacteria = -1
      iint_MarineOA = -1
      itot_MarineOA = -1


      if (USES%FUNGAL_SPORES) then
          if (USES%FUNGAL_METHOD == "HS") then
            itot_FungalSpores = find_index( "FUNGAL_SPORES_5", species(:)%name)
            inat_FungalSpores = find_index( "FUNGAL_SPORES_5", EMIS_BioNat(:))
            FUNGAL_DIAMETER = 5.0 !um

          else
            itot_FungalSpores = find_index( "FUNGAL_SPORES_3", species(:)%name)
            inat_FungalSpores = find_index( "FUNGAL_SPORES_3", EMIS_BioNat(:))
            FUNGAL_DIAMETER = 3.0 !um
          end if 

          if (itot_FungalSpores < 0 ) then
            if(MasterProc)  write(*,*) "WARNING: No fungal spores found in species, not including fungal spores!"
          else
            FUNGAL_WEIGHT = (4/3.0)*FUNGAL_DENS*PI*(0.5*FUNGAL_DIAMETER*1e-6)**3!Spore weight [g]
            NPBAP = NPBAP + 1
            iint_FungalSpores = NPBAP
            if(MasterProc) write(*,*) "USING FUNGAL_METHOD:",USES%FUNGAL_METHOD
          end if
       end if


       if (USES%BACTERIA) then
          itot_Bacteria = find_index( "BACTERIA", species(:)%name)
          if (itot_Bacteria< 0) then
            if(MasterProc)  write(*,*) "WARNING: No bacteria found in species, not including bacteria!"
          else
            NPBAP = NPBAP + 1
            iint_Bacteria = NPBAP
            inat_Bacteria = find_index( "BACTERIA", EMIS_BioNat(:))
          end if
       end if

       if (USES%MARINE_OA) then
        itot_MarineOA = find_index( "MARINE_OA_NEW", species(:)%name)
        if (itot_MarineOA< 0) then
          if(MasterProc)  write(*,*) "WARNING: No Marine OA found in species, not including Marine OA!"
        else
          allocate(O_Chlorophyll(LIMAX,LJMAX))
          if(MasterProc) write(*,*)'Reading Ocean Chlorophyll'
          !call ReadField_CDF(trim(OceanChlorophyll_File),'chlor_a',O_Chlorophyll,&
          !      nstart=current_date%month+12*3,interpol='conservative',known_projection="lon lat",&
          !      needed=.true.,debug_flag=.false.,UnDef=0.0) !In mg/m3
            NPBAP = NPBAP + 1
            iint_MarineOA = NPBAP
            inat_MarineOA = find_index( "MarineOA_NEW", EMIS_BioNat(:))
        end if
     end if


       if (NPBAP > 0) then
          allocate(PBAP_Flux(LIMAX,LJMAX,NPBAP))
          PBAP_Flux = 0.0

          allocate(WEIGHTS(NPBAP))
          allocate(DIAMETERS(NPBAP))
          allocate(DENSITIES(NPBAP))

          allocate(n2m(NPBAP))
          n2m = 0.0

          allocate(kgm2h(NPBAP))
          kgm2h = 0.0

          allocate(inat(NPBAP))
          allocate(itot(NPBAP))
          allocate(PBAP_names(NPBAP))

          if (iint_FungalSpores > 0) then
              WEIGHTS(iint_FungalSpores) = FUNGAL_WEIGHT
              DIAMETERS(iint_FungalSpores) = FUNGAL_DIAMETER
              DENSITIES(iint_FungalSpores) = FUNGAL_DENS
              inat(iint_FungalSpores) = inat_FungalSpores
              itot(iint_FungalSpores) = itot_FungalSpores
              PBAP_names(iint_FungalSpores) = "FUNGAL_SPORES"
              if (FUNGAL_DIAMETER > 3.5) then
                FUNG_PARAMS_HS(1) = 500 !New parameterization (see Ref. above)
                                        !for larger particles according to
                                        !Heald and Spracken
              end if
          end if

          if (iint_Bacteria > 0) then
              WEIGHTS(iint_Bacteria) = BACTERIA_WEIGHT
              DIAMETERS(iint_Bacteria) = BACTERIA_DIAMETER
              DENSITIES(iint_Bacteria) = BACTERIA_DENS
              inat(iint_Bacteria) = inat_Bacteria
              itot(iint_Bacteria) = itot_Bacteria
              PBAP_names(iint_Bacteria) = "BACTERIA"
          end if

          if (iint_MarineOA > 0) then
              inat(iint_MarineOA) = inat_MarineOA
              itot(iint_MarineOA) = itot_MarineOA
              PBAP_names(iint_MarineOA) = "MARINE_OA"
              n2m(iint_MarineOA) = 1
              kgm2h(iint_MarineOA) = 1
          end if
      end if !NPBAP > 0

     my_first_call = .false.

     if( DEBUG%PBAP .and. debug_proc ) then
         write(*,*)"INIT PBAPs (should only happen once!). NPBAP = ",NPBAP
     end if

    end if

  end subroutine init_PBAPs



  subroutine Set_FungalSpores(i,j)
    !!!!!!!
    !Fills PBAP_Flux(i,j,iint_FugalSpores)
    !Currently three possible choices of parameterization

    integer, intent(in) ::  i,j
    integer :: nlu,iiL,LC,i_d,j_d
    real    :: F_FNG, temp_val,sum_LC
    logical, save ::  my_first_call = .true.

    if( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
      if (i .eq. debug_li .and. j .eq. debug_lj) then
       write(*,*)"PBAP_mod DEBUG FUNGAL_SPORES: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%FUNGAL_SPORES, itot(iint_FungalSpores),inat(iint_FungalSpores)
      end if
    end if

    if (my_first_call) then
      n2m(iint_FungalSpores) = (1e-6/Grid%DeltaZ)*WEIGHTS(iint_FungalSpores)*AVOG/species(itot(iint_FungalSpores))%molwt
      !Converts from number of spores/m^2/s -> mol/cm^3/s
      kgm2h(iint_FungalSpores) = WEIGHTS(iint_FungalSpores)*1e-6*3600 !num/m2/s -> kg/m2/h
      my_first_call = .false.
    end if

    F_FNG = 0.0 !Fungal spores flux

    nlu = LandCover(i,j)%ncodes
    sum_LC = 0.0

    if (USES%FUNGAL_METHOD=="HM") then
      do iiL = 1,nlu
        LC = LandCover(i,j)%codes(iiL)
        sum_LC = sum_LC + LandCover(i,j)%fraction(iiL)
        if (LandCover(i,j)%fraction(iiL)<3e-4) then !Avoid integrated assesment (IAM) landtypes as no LAI
          cycle
        else if ( LandType(LC)%is_water) then
            cycle
        else if ( LandType(LC)%is_ice) then
            cycle
        else
          temp_val =  FUNG_PARAMS_HM(1)*(t2_nwp(i,j,1)-FUNG_PARAMS_HM(2))+FUNG_PARAMS_HM(3)*q(i,j,KG,1)*LandCover(i,j)%LAI(iiL)
          temp_val = LandCover(i,j)%fraction(iiL)*temp_val
          F_FNG = F_FNG + max(0.0, temp_val)
          !Eq.(2) of S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
          !DOI 10.1007/978-3-319-35095-0_121, scaled by fraction
          !Leaf-area index (LAI) should be in m2/m2
          !Specific humidity q should be in kg/kg
          !Temperature at 2m (t2_nwp) should be in K
          !Based on Hummel et al. Atmos. Chem. Phys., 15, 6127–6146,
          !https://doi.org/10.5194/acp-15-6127-2015, 2015

        end if
      end do !iiL

    else if (USES%FUNGAL_METHOD=="SD") then
      do iiL = 1,nlu
          sum_LC = sum_LC + LandCover(i,j)%fraction(iiL)
          if (LandDefs(iil)%FungalFlux > 0.0) then
            temp_val = LandDefs(iiL)%FungalFlux
          end if

          F_FNG = F_FNG + LandCover(i,j)%fraction(iiL)*temp_val

          !Flux from S. Sesartic and T.N Dallafiro (2011)
          !DOI 10.5194/bg-8-1181-2011
      end do !iiL

  else if (USES%FUNGAL_METHOD == "HS" .or. USES%FUNGAL_METHOD == "HO") then
    do iiL = 1,nlu
      LC = LandCover(i,j)%codes(iiL)
      sum_LC = sum_LC + LandCover(i,j)%fraction(iiL)
      if (LandCover(i,j)%fraction(iiL)<3e-4) then !Avoid integrated assesment (IAM) landtypes as no LAI
        cycle
      else if ( LandType(LC)%is_water) then
          cycle
      else if ( LandType(LC)%is_ice) then
          cycle
      else
        temp_val =  LandCover(i,j)%fraction(iiL)*(FUNG_PARAMS_HS(1)*(q(i,j,KG,1)*LandCover(i,j)%LAI(iiL))/(FUNG_PARAMS_HS(2)*FUNG_PARAMS_HS(3)))
        F_FNG = F_FNG + max(0.0, temp_val)
        !Flux from  Heald and Spracklen (2009) DOI:10.1029/2009GL037493,
        !as specified in Hoose et al (2010) DOI: 10.1088/1748-9326/5/2/024009
        !for spores of size 5um and in Hummel (2015) DOI: 10.5194/acp-15-6127-2015
        !for spores of size 3um.
      end if
    end do !iiL
  
  else if (USES%FUNGAL_METHOD == "JS") then !N.B: Has not been tested (April 2024)
    do iiL = 1,nlu
      LC = LandCover(i,j)%codes(iiL)
      sum_LC = sum_LC + LandCover(i,j)%fraction(iiL)
      if (LandCover(i,j)%fraction(iiL)<3e-4) then !Avoid integrated assesment (IAM) landtypes as no LAI
        cycle
      else if ( LandType(LC)%is_water) then
          cycle
      else if ( LandType(LC)%is_ice) then
          cycle
      else
        temp_val = FUNG_PARAMS_JS(1)+FUNG_PARAMS_JS(2)*q(i,j,KG,1)
        temp_val = temp_val+FUNG_PARAMS_JS(3)*LandCover(i,j)%LAI(iiL)
        temp_val = temp_val+FUNG_PARAMS_JS(4)*ustar_nwp(i,j)
        temp_val =  LandCover(i,j)%fraction(iiL)*(temp_val)
        F_FNG = F_FNG + max(0.0, temp_val)
        !Eq.(2) of Janssen et. al. (2021) 
        !Atmos. Chem. Phys., 21, 4381–4401,
        !https://doi.org/10.5194/acp-21-4381-2021


      end if
    end do !iiL  
  else
      call StopAll('Unknown FUNGAL_METHOD chosen! Valid options are HM, HS, SD and JS. ')
  end if

    PBAP_flux(i,j,iint_FungalSpores) = F_FNG

    if ( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
       if (i .eq. debug_li .and. j .eq. debug_lj) then
        write(*,"(a,4i4)") "FUNGAL_SPORES i,j: ",  1, limax, 1, ljmax
        write(*,"(a,2f12.4)") "Unit conversions: ", n2m(iint_FungalSpores), kgm2h(iint_FungalSpores)
        write(*,*) "Fungal flux:", F_FNG
        if (abs(sum_LC-1)>1e-4) then
          write(*,*) "WARNING: Land Cover classes fraction do not sum to 1", sum_LC
        end if
      end if
    end if

  end subroutine Set_FungalSpores

  subroutine Set_Bacteria(i,j)
    !NOT YET TESTED
    !!!!!!!
    !Fills PBAP_Flux(i,j,iint_Bacteria)
    !Based on parameterization from
    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
    !DOI 10.1007/978-3-319-35095-0_121
    integer, intent(in) ::  i,j
    integer :: nlu,iiL,LC,i_d,j_d
    real    :: F_Bacteria, temp_val
    logical, save ::  my_first_call = .true.


    if( DEBUG%BACTERIA .and. debug_proc ) then
      if (i .eq. debug_li .and. j .eq. debug_lj) then
       write(*,*)"PBAP_mod DEBUG BACTERIA: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%BACTERIA, itot(iint_Bacteria),inat(iint_Bacteria)
      end if
    end if

    if (my_first_call) then
      n2m(iint_Bacteria) = (1e-6/Grid%DeltaZ)*WEIGHTS(iint_Bacteria)*AVOG/species(itot(iint_Bacteria))%molwt
      !Converts from number of spores/m^2/s -> mol/cm^3/s
      kgm2h(iint_Bacteria) = WEIGHTS(iint_Bacteria)*1e-6*3600 !num/m2/s -> kg/m2/h
      my_first_call = .false.
    end if

    F_Bacteria = 0.0 !Bacteria flux
    temp_val = 0.0

    if (likely_coastal(i,j)) then
      F_Bacteria = BACTERIA_PARAMS(1) !Coastal not LandType, therefore treated
                                      !differentely, following parameterization from
                                      !Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                                      !DOI 10.1007/978-3-319-35095-0_121

    else
      nlu = LandCover(i,j)%ncodes

      do iiL = 1,nlu
          if (LandDefs(iil)%BacteriaFlux > 0.0) then
            temp_val = LandDefs(iiL)%BacteriaFlux
          end if

          F_Bacteria = F_Bacteria + LandCover(i,j)%fraction(iiL)*temp_val

          !Eq.(1) of S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
          !DOI 10.1007/978-3-319-35095-0_121, scaled by fraction
          !Note that most of these valuse are set in the LandUse input file
      end do !iiL
    end if


    PBAP_flux(i,j,iint_Bacteria) = F_Bacteria

    if ( DEBUG%BACTERIA .and. debug_proc ) then
       if (i .eq. debug_li .and. j .eq. debug_lj) then
        write(*,"(a,4i4)") "Bacteria i,j: ",  1, limax, 1, ljmax
        write(*,"(a,2f12.4)") "Unit conversions: ", n2m(iint_Bacteria), kgm2h(iint_Bacteria)
       end if
    end if

  end subroutine Set_Bacteria

  subroutine Set_MarineOA(i,j)
    !NOT YET TESTED
    !!!!!!!
    !Fills PBAP_Flux(i,j,iint_MarineOA)
    !Based on parameterization from

    integer, intent(in) ::  i,j
    integer :: nlu, iiL,lu
    real    :: F_MarineOA, temp_val


    if( DEBUG%MARINE_OA .and. debug_proc ) then
      if (i .eq. debug_li .and. j .eq. debug_lj) then
       write(*,*)"PBAP_mod DEBUG Marine OA: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%MARINE_OA, itot(iint_MarineOA),inat(iint_MarineOA)
      end if
    end if


    F_MarineOA = 0.0 !MarineOA flux
    temp_val = 0.0

     !if (.not. USES%SEASALT) then
     !  call SeaSalt_flux(i,j,DEBUG%MARINE_OA) ! sets rcemis(SEASALT_...)
     !end if


      !nlu = LandCover(i,j)%ncodes

      !do iiL = 1,nlu
      !  lu =  LandCover(i,j)%codes(iiL)
      !  if ( Sub(lu)%is_water ) then
      !    F_MarineOA = F_MarineOA + LandCover(i,j)%fraction(iiL)*0.4*O_Chlorophyll(i,j)
      !  end if

      !end do !iiL


    PBAP_flux(i,j,iint_MarineOA) = F_MarineOA

    if ( DEBUG%MARINE_OA .and. debug_proc ) then
       if (i .eq. debug_li .and. j .eq. debug_lj) then
        write(*,"(a,4i4)") "MarineOA i,j: ",  1, limax, 1, ljmax
        write(*,"(a,2f12.4)") "Flux, Cholorphyll: ", F_MarineOA, O_Chlorophyll(i,j)
      end if
    end if

  end subroutine Set_MarineOA


  subroutine set_PBAPs(i,j)
  !
  !---- Adds PBAPs to rcemis  and EmisNat--------------------------------------
  !
  !  So far, only adds Fungal Spores and Bacteria to rcemis and EmisNat
  !
  !  Called from setup_1d_mod, every  advection step.
  !----------------------------------------------------------------------------

  integer, intent(in) ::  i,j

  character(len=*), parameter :: dtxt='PBAPModSetup:'

  logical :: dbg

  integer :: i_PBAP
  if ( NPBAP == 0  ) return   ! Number of PBAPs

  dbg = ( DEBUG%PBAP .and. debug_proc .and. &
          i==debug_li .and. j==debug_lj .and. current_date%seconds == 0 )


  if ( iint_FungalSpores > 0 ) then
      call Set_FungalSpores(i,j)
  end if

  if ( iint_Bacteria > 0 ) then
    call Set_Bacteria(i,j)
  end if

  if (iint_MarineOA > 0) then
    call Set_MarineOA(i,j)
  end if

  do i_PBAP = 1,NPBAP
    rcemis(itot(i_PBAP),KG) = rcemis(itot(i_PBAP),KG)+n2m(i_PBAP)*PBAP_flux(i,j,i_PBAP)![mol/cm3/s]
    if (dbg) then
          write(*,"(2a,f12.5)") PBAP_names(i_PBAP),": rcemis ",rcemis(itot(i_PBAP),KG)
    end if

    if (inat(i_PBAP) > 0) then
        EmisNat(inat(i_PBAP),i,j) = kgm2h(i_PBAP)*PBAP_flux(i,j,i_PBAP) !Emissions in molec/m2/s
      if (dbg) then
        write(*,"(2a,f12.5)") PBAP_names(i_PBAP),": EmisNat ",EmisNat(inat(i_PBAP),i,j)
      end if
    end if
  end do

  end subroutine set_PBAPs

end module PBAP_mod
