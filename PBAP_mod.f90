!> <PBAP_mod.f90 - A component of the EMEP MSC-W Chemical transport Model>
!  **************************************************************************!

module PBAP_mod

  !/-- Module to deal with primary biological aerosol particles (PBAPs).
  !    IMPLEMENTED (Aug 2024):
  !    - Fungal parameterization schemes (JS parameterization has not been tested)
  !    - Bacteria based on LandUse classes (NOT TESTED)
  !
  !    TODO:
  !    - Test bacteria
  !    - Marine OA (most of the code is there, but not yet tested)
  !    - Aging of fungal spores
  !    - Water uptake of fungal Spores
  !
  !    Gunnar Felix Lange 2024
  !---------------------------------------------------------------------------

  use CheckStop_mod,      only: CheckStop, StopAll
  use ChemSpecs_mod,         only : species
  use Config_module, only : NPROC, MasterProc,&
                           KT => KCHEMTOP, KG => KMAX_MID,&
                           USES, OceanChlorophyll_File
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
  use ZchemData_mod,      only: rcemis

  implicit none
  private

  !/-- subroutines for all PBAPs
  public ::  init_PBAPs,set_PBAPs
  !/-- subroutines for individual PBAPs
  private :: set_fungal_spores, set_bacteria, set_marineOA

  integer, public, save ::   NPBAP !Number of PBAPs!(only fungal spores, bacteria and marine organic aerosol (OA) implemented at the moment)
                                   !Only fungal spores have been properly tested!

  real,public, save, allocatable, dimension(:,:,:) :: PBAP_flux !Dim: i,j,NPBAP

  real,public, save, allocatable, dimension(:,:) :: O_Chlorophyll !Amount of Chlorophyll in the ocean. Dim: i,j
                                                                  !Used for marineOA (not yet implemented)


  real,public, save, allocatable, dimension(:) :: WEIGHTS ,DIAMETERS, DENSITIES !Physical parameters. Dim: NPBAP
  real,public, save, allocatable, dimension(:) :: num2emis  !Conversion factor. Dim: NPBAP.
                                                            !num2emis: number emission flux [num/m2/s] to units of rcemis [molec/cm3/s]
                                                            !The area -> volume conversion is done by multiplying by height of lowest gridbox
                                                            !The molar weight of the PBAPs is just a dummy value (set to 1), as it gets divided
                                                            !out later.

  real,private, save, allocatable, dimension(:) :: itot !indices in species. Dim: NPBAP

  character(len=20),private,save,allocatable,dimension(:) ::PBAP_names !For debugging only


  integer, private, save :: iint_fungal_spores, itot_fungal_spores !Index of fungal spores internally/in species
  integer, private, save :: iint_bacteria,itot_bacteria !Index of bacteria internally/in species
  integer, private, save :: iint_marineOA,itot_marineOA !Index of marine OA internally/in species

  !!!!!!!!!!!!!!!!!FUNGAL PARAMETERIZATION CHOICES!!!!!!!!!!!!!!!!
  !Choice of fungal flux parameterization:
  !HM: Hummel(2015) DOI: 10.5194/acp-15-6127-2015
  !SD: Sesartic and Dallafior (2011) DOI: 10.5194/bg-8-1181-2011
  !HS_3: Heald and Spracklen (2009) DOI:10.1029/2009GL037493 for 3um spores [see Hummel(2015) DOI: 10.5194/acp-15-6127-2015]
  !HS_5: Heald and Spracklen (2009) DOI:10.1029/2009GL037493 for 5um spores [see Hoose et al (2010) DOI: 10.1088/1748-9326/5/2/024009]
  !JS: Janssen et al. (2021) DOI: 10.5194/acp-21-4381-2021 [NOT TESTED!]
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !For similar implementation in CHIMERE,see Vida et al. DOI: 10.5194/egusphere-2024-698

  !(GFL Feb 2024): Hummel parameterization uses different settling scheme than EMEP, Sesartic and Dallfior uses
  !no settling scheme at all, so currently seems that Heald and Spracken parameterization modified by Hoose
  !(parameterization choice = "HS_5") gives best results for fungal spores


  real, public :: FUNGAL_DIAMETER, FUNGAL_WEIGHT !Diameter depends on chosen fungal parameterization scheme so allocated dynamically
  real*8, DIMENSION(3), parameter  ::  &
  FUNG_PARAMS_HM = [20.426, 275.82, 39300.0] !From Fungal paramterization, Eq. (2) of
  !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
  !DOI 10.1007/978-3-319-35095-0_121, based on Hummel et al. Atmos. Chem. Phys., 15, 6127–6146,
                                        !       https://doi.org/10.5194/acp-15-6127-2015, 2015

  !real*8, DIMENSION(3), parameter  ::  &
  !FUNG_PARAMS_HS_LARGE = [500.0, 5.0,0.015] !From Fungal parameterization, Heald and Spracken
                                       !fungal spores. DOI:10.1029/2009GL037493, for 5um spores
                                       !described by Hoose et al. (2010), DOI: 10.1088/1748-9326/5/2/024009.
                                       !Gets updated automatically dependent on parameterization choice

  real*8, DIMENSION(3)  ::  &
  FUNG_PARAMS_HS = [2315.0, 5.0,0.015] !From Fungal parameterization, Heald and Spracken for 3um
                                       !as described in Hummel et al. (2015) DOI: https://doi.org/10.5194/acp-15-6127-2015


  real*8, DIMENSION(4)  ::  &
  FUNG_PARAMS_JS = [2.63*1.0e-5, 6.10*1.0e3,46.7,59.0] !From Fungal parameterization, Janssen (2021) for North America.
                                       !DOI: 10.5194/acp-21-4381-2021 [N.B: Has not been tested!]



  real, parameter :: FUNGAL_DENS = 1.0e6 !Fungal density [g/m3]
                                         !From Hummel et al. Atmos. Chem. Phys., 15, 6127–6146,
                                         !https://doi.org/10.5194/acp-15-6127-2015, 2015

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Bacteria PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8, DIMENSION(6), parameter  ::  &
  Bacteria_PARAMS = [900.0,704.0,648.0,7.7,502.0,196.0] !From bacteria paramterization, Eq. (1) of
                    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                    !DOI 10.1007/978-3-319-35095-0_121
                    !Note that most of these are set in the LandInput file, except for the coastal parameter (as coastal is not a LandType)
                    !Has not been tested extensively.


  real, parameter :: Bacteria_DIAMETER = 1.0 !bacteria diameter [um]
                                             !From S. M. Burrows et al. Atmos. Chem. Phys., 9, 9281–9297
                                             !https://doi.org/10.5194/acp-9-9281-2009

  real, parameter :: Bacteria_WEIGHT = 0.52*1.0e-12 !bacterial weight [g]  (ibid)
  real, parameter :: Bacteria_DENS   = bacteria_WEIGHT/((4/3.0)*PI*(0.5*bacteria_DIAMETER*1e-6)**3) !bacteria density [g/m3]


  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine init_PBAPs()
    logical, save ::  my_first_call = .true.

    !Initializes the PBAP calculation
    !sets up number of PBAPs (NPBAP)
    !And fills array accordingly.
    !GFL Jul 2024: Is this parallelization safe?

    if (my_first_call) then
      NPBAP = 0
      iint_fungal_spores = -1
      itot_fungal_spores = -1
      iint_bacteria = -1
      itot_bacteria = -1
      iint_marineOA = -1
      itot_marineOA = -1


      if (USES%FUNGAL_SPORES) then
          if (USES%FUNGAL_METHOD == "HS_5") then !Only this method uses 5um spores
            itot_fungal_spores = find_index( "FUNGAL_SPORES_5", species(:)%name)
            FUNGAL_DIAMETER = 5.0 !um

          else
            itot_fungal_spores = find_index( "FUNGAL_SPORES_3", species(:)%name)
            FUNGAL_DIAMETER = 3.0 !um
          end if 

          if (itot_fungal_spores < 0 ) then
            if(MasterProc)  write(*,*) "WARNING: No fungal spores found in species, not including fungal spores!"
          else
            FUNGAL_WEIGHT = (4/3.0)*FUNGAL_DENS*PI*(0.5*FUNGAL_DIAMETER*1e-6)**3 !Spore weight [g]
            NPBAP = NPBAP + 1
            iint_fungal_spores = NPBAP
            if(MasterProc) write(*,*) "USING FUNGAL_METHOD:",USES%FUNGAL_METHOD
          end if
       end if


       if (USES%bacteria) then
          itot_bacteria = find_index( "Bacteria", species(:)%name)
          if (itot_bacteria< 0) then
            if(MasterProc)  write(*,*) "WARNING: No Bacteria found in species, not including bacteria!"
          else
            NPBAP = NPBAP + 1
            iint_bacteria = NPBAP
          end if
       end if

       if (USES%MARINE_OA) then
        itot_marineOA = find_index( "MARINE_OA_NEW", species(:)%name)
        if (itot_marineOA< 0) then
          if(MasterProc)  write(*,*) "WARNING: No Marine OA found in species, not including Marine OA!"
        else
          allocate(O_Chlorophyll(LIMAX,LJMAX))
          if(MasterProc) write(*,*)'Reading Ocean Chlorophyll'
            !call ReadField_CDF(trim(OceanChlorophyll_File),'chlor_a',O_Chlorophyll,&
            !    nstart=current_date%month+12*3,interpol='conservative',known_projection="lon lat",&
            !    needed=.true.,debug_flag=.false.,UnDef=0.0) !In mg/m3. Should call this and add to fraction of seasalt.
            NPBAP = NPBAP + 1
            iint_marineOA = NPBAP
        end if
     end if


       if (NPBAP > 0) then
          allocate(PBAP_Flux(LIMAX,LJMAX,NPBAP))
          PBAP_Flux = 0.0

          allocate(WEIGHTS(NPBAP))
          allocate(DIAMETERS(NPBAP))
          allocate(DENSITIES(NPBAP))

          allocate(num2emis(NPBAP))
          num2emis = 0.0

          allocate(itot(NPBAP)) !Indices in species
          allocate(PBAP_names(NPBAP)) !For debugging

          if (iint_fungal_spores > 0) then
              WEIGHTS(iint_fungal_spores) = FUNGAL_WEIGHT
              DIAMETERS(iint_fungal_spores) = FUNGAL_DIAMETER
              DENSITIES(iint_fungal_spores) = FUNGAL_DENS
              itot(iint_fungal_spores) = itot_fungal_spores
              PBAP_names(iint_fungal_spores) = "FUNGAL_SPORES"
              if (FUNGAL_DIAMETER > 3.5) then
                FUNG_PARAMS_HS(1) = 500 !Parameterization (see Ref. above)
                                        !for larger spores according to
                                        !Heald and Spracken
              end if
              num2emis(iint_fungal_spores) = WEIGHTS(iint_fungal_spores)*AVOG/species(itot(iint_fungal_spores))%molwt
              !Converts from number of spores/m2/s -> mol/m2/s (we eventually output cm3 by dividing by grid box further down)
          end if

          if (iint_bacteria > 0) then
              WEIGHTS(iint_bacteria) = bacteria_WEIGHT
              DIAMETERS(iint_bacteria) = bacteria_DIAMETER
              DENSITIES(iint_bacteria) = bacteria_DENS
              itot(iint_bacteria) = itot_bacteria
              PBAP_names(iint_bacteria) = "bacteria"
              num2emis(iint_bacteria) = (1e-6)*WEIGHTS(iint_bacteria)*AVOG/species(itot(iint_bacteria))%molwt
              !Converts from number of spores/m2/s -> mol/cm2/s (we eventually output cm3 by dividing by grid box further down)
            end if

          if (iint_marineOA > 0) then
              itot(iint_marineOA) = itot_marineOA
              PBAP_names(iint_marineOA) = "MARINE_OA"
              num2emis(iint_marineOA) = 1 !Not yet implemented!
          end if
      end if !NPBAP > 0

     my_first_call = .false.

     if( (DEBUG%FUNGAL_SPORES .or. DEBUG%BACTERIA .or. DEBUG%MARINE_OA) .and. debug_proc ) then
         write(*,*)"INIT PBAPs (should only happen once!). NPBAP = ",NPBAP
     end if

    end if

  end subroutine init_PBAPs



  subroutine set_fungal_spores(i,j)
    !!!!!!!
    !Fills PBAP_Flux(i,j,iint_FugalSpores)
    !!!!!!!

    integer, intent(in) ::  i,j
    integer :: nlu,iiL,LC,i_d,j_d
    real    :: F_FNG, temp_val,sum_LC

    if( DEBUG%FUNGAL_SPORES .and. debug_proc ) then
      if (i .eq. debug_li .and. j .eq. debug_lj) then
       write(*,*)"PBAP_mod DEBUG FUNGAL_SPORES: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%FUNGAL_SPORES, itot(iint_fungal_spores)
      end if
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

  else if (USES%FUNGAL_METHOD == "HS_3" .or. USES%FUNGAL_METHOD == "HS_5") then
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
  
  else if (USES%FUNGAL_METHOD == "JS") then !N.B: Has not been tested (July 2024)
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
      call StopAll('Unknown FUNGAL_METHOD chosen! Valid options are HM, SD, HS_3, HS_5 and JS. ')
  end if

    PBAP_flux(i,j,iint_fungal_spores) = F_FNG

    if ( DEBUG%FUNGAL_SPORES .and. debug_proc) then
       if (i .eq. debug_li .and. j .eq. debug_lj) then
        write(*,"(a,4i4)") "FUNGAL_SPORES i,j: ",  1, limax, 1, ljmax
        write(*,"(a,1f12.4)") "Unit conversion: ", num2emis(iint_fungal_spores)
        write(*,*) "Fungal flux:", F_FNG
        if (abs(sum_LC-1)>1e-4) then
          write(*,*) "WARNING: Land Cover classes fraction do not sum to 1", sum_LC
        end if
      end if
    end if

  end subroutine set_fungal_spores

  subroutine set_bacteria(i,j)
    !NOT YET TESTED
    !!!!!!!
    !Fills PBAP_Flux(i,j,iint_bacteria)
    !Based on parameterization from
    !S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
    !DOI 10.1007/978-3-319-35095-0_121
    integer, intent(in) ::  i,j
    integer :: nlu,iiL,LC,i_d,j_d
    real    :: F_bacteria, temp_val

    if( DEBUG%bacteria .and. debug_proc ) then
      if (i .eq. debug_li .and. j .eq. debug_lj) then
       write(*,*)"PBAP_mod DEBUG bacteria: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%bacteria, itot(iint_bacteria)
      end if
    end if

    F_bacteria = 0.0 !bacteria flux
    temp_val = 0.0

    if (likely_coastal(i,j)) then
      F_bacteria = bacteria_PARAMS(1) !Coastal not LandType, therefore treated
                                      !differentely, following parameterization from
                                      !Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
                                      !DOI 10.1007/978-3-319-35095-0_121

    else
      nlu = LandCover(i,j)%ncodes

      do iiL = 1,nlu
          if (LandDefs(iil)%bacteriaFlux > 0.0) then
            temp_val = LandDefs(iiL)%bacteriaFlux
          end if

          F_bacteria = F_bacteria + LandCover(i,j)%fraction(iiL)*temp_val

          !Eq.(1) of S. Myriokefalitakis, G. Fanourgakis and M. Kanakidou (2017)
          !DOI 10.1007/978-3-319-35095-0_121, scaled by fraction
          !Note that most of these valuse are set in the LandUse input file
      end do !iiL
    end if

    PBAP_flux(i,j,iint_bacteria) = F_bacteria

    if ( DEBUG%bacteria .and. debug_proc ) then
       if (i .eq. debug_li .and. j .eq. debug_lj) then
        write(*,"(a,4i4)") "bacteria i,j: ",  1, limax, 1, ljmax
        write(*,"(a,1f12.4)") "Unit conversion: ", num2emis(iint_bacteria)
       end if
    end if

  end subroutine set_bacteria

  subroutine set_marineOA(i,j)
    !NOT YET IMPLEMENTED!!
    !!!!!!!
    !Fills PBAP_Flux(i,j,iint_marineOA)


    integer, intent(in) ::  i,j
    integer :: nlu, iiL,lu
    real    :: F_marineOA, temp_val


    if( DEBUG%MARINE_OA .and. debug_proc ) then
      if (i .eq. debug_li .and. j .eq. debug_lj) then
       write(*,*)"PBAP_mod DEBUG Marine OA: ",&
        current_date%day, current_date%hour, current_date%seconds,&
        USES%MARINE_OA, itot(iint_marineOA)
      end if
    end if


    F_marineOA = 0.0 !marineOA flux
    temp_val = 0.0

     !if (.not. USES%SEASALT) then
     !  call SeaSalt_flux(i,j,DEBUG%MARINE_OA) ! sets rcemis(SEASALT_...)
     !end if


      !nlu = LandCover(i,j)%ncodes

      !do iiL = 1,nlu
      !  lu =  LandCover(i,j)%codes(iiL)
      !  if ( Sub(lu)%is_water ) then
      !    F_marineOA = F_marineOA + LandCover(i,j)%fraction(iiL)*0.4*O_Chlorophyll(i,j)
      !  end if

      !end do !iiL


    PBAP_flux(i,j,iint_marineOA) = F_marineOA

    if ( DEBUG%MARINE_OA .and. debug_proc ) then
       if (i .eq. debug_li .and. j .eq. debug_lj) then
        write(*,"(a,4i4)") "marineOA i,j: ",  1, limax, 1, ljmax
        write(*,"(a,2f12.4)") "Flux, Cholorphyll: ", F_marineOA, O_Chlorophyll(i,j)
      end if
    end if

  end subroutine set_marineOA


  subroutine set_PBAPs(i,j)
  !
  !---- Adds PBAPs to rgcemis--------------------------------------
  !
  !
  !  Called from setup_1d_mod, every  advection step.
  !----------------------------------------------------------------------------

  integer, intent(in) ::  i,j

  character(len=*), parameter :: dtxt='PBAPModSetup:'

  logical :: dbg

  integer :: i_PBAP
  if ( NPBAP == 0  ) return   ! Number of PBAPs

  dbg = ((DEBUG%FUNGAL_SPORES .or. DEBUG%BACTERIA .or. DEBUG%MARINE_OA) .and. debug_proc .and. &
        i == debug_li .and. j == debug_lj)

  if ( iint_fungal_spores > 0 ) then
      call set_fungal_spores(i,j)
  end if

  if ( iint_bacteria > 0 ) then
    call set_bacteria(i,j)
  end if

  if (iint_marineOA > 0) then
    call set_marineOA(i,j)
  end if

  do i_PBAP = 1,NPBAP
    rcemis(itot(i_PBAP),KG) = rcemis(itot(i_PBAP),KG)+(1e-6/Grid%DeltaZ)*num2emis(i_PBAP)*PBAP_flux(i,j,i_PBAP)!Emission in [molec/cm3/s]
    if (dbg) then
          write(*,*) PBAP_names(i_PBAP),": rcemis, deltaZ: ",rcemis(itot(i_PBAP),KG), Grid%DeltaZ
          write(*,*) PBAP_names(i_PBAP),": num2emis, flux: ",num2emis(i_PBAP),PBAP_flux(i,j,i_PBAP)
    end if

  end do

  end subroutine set_PBAPs

end module PBAP_mod
