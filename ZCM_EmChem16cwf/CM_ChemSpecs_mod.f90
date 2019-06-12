! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem16x BVOC_EmChem16x Aqueous_EmChem16x Aero2017nx VBS_EmChem16x Ash_EmChem16x ShipNOx FFireInert SeaSalt DustExtended Pollen
module ChemSpecs_mod

  use ChemDims_mod      ! => NSPEC_TOT, NCHEMRATES, ....
  
  implicit none
  private
  character(len=*),parameter, public :: CM_schemes_ChemSpecs = " EmChem16x BVOC_EmChem16x Aqueous_EmChem16x Aero2017nx VBS_EmChem16x Ash_EmChem16x ShipNOx FFireInert SeaSalt DustExtended Pollen"
  
  
  integer, public, parameter :: &
      RO2POOL     =   1  &
    , OD          =   2  &
    , OP          =   3  &
    , OH          =   4  &
    , HO2         =   5  &
    , CH3O2       =   6  &
    , C2H5O2      =   7  &
    , SECC4H9O2   =   8  &
    , ISRO2       =   9  &
    , ETRO2       =  10
  
  integer, public, parameter :: &
      PRRO2       =  11  &
    , OXYO2       =  12  &
    , MEKO2       =  13  &
    , MALO2       =  14  &
    , MACRO2      =  15  &
    , CH3CO3      =  16  &
    , TERPO2      =  17  &
    , XMTO3_RO2   =  18  &
    , NMVOC       =  19  &
    , O3          =  20
  
  integer, public, parameter :: &
      NO          =  21  &
    , NO2         =  22  &
    , PAN         =  23  &
    , MPAN        =  24  &
    , NO3         =  25  &
    , N2O5        =  26  &
    , HNO3        =  27  &
    , HONO        =  28  &
    , HO2NO2      =  29  &
    , MACR        =  30
  
  integer, public, parameter :: &
      ISON        =  31  &
    , GLYOX       =  32  &
    , MGLYOX      =  33  &
    , MAL         =  34  &
    , MEK         =  35  &
    , HCHO        =  36  &
    , CH3CHO      =  37  &
    , C2H6        =  38  &
    , NC4H10      =  39  &
    , C2H4        =  40
  
  integer, public, parameter :: &
      C3H6        =  41  &
    , OXYL        =  42  &
    , C5H8        =  43  &
    , APINENE     =  44  &
    , CH3OOH      =  45  &
    , C2H5OOH     =  46  &
    , BURO2H      =  47  &
    , ETRO2H      =  48  &
    , PRRO2H      =  49  &
    , OXYO2H      =  50
  
  integer, public, parameter :: &
      MEKO2H      =  51  &
    , MALO2H      =  52  &
    , MACROOH     =  53  &
    , ISRO2H      =  54  &
    , H2O2        =  55  &
    , CH3CO3H     =  56  &
    , CH3OH       =  57  &
    , C2H5OH      =  58  &
    , ACETOL      =  59  &
    , NALD        =  60
  
  integer, public, parameter :: &
      HPALD       =  61  &
    , PACALD      =  62  &
    , IEPOX       =  63  &
    , H2          =  64  &
    , CO          =  65  &
    , CH4         =  66  &
    , SO2         =  67  &
    , MVK         =  68  &
    , BPINENE     =  69  &
    , XTERP       =  70
  
  integer, public, parameter :: &
      SQT_SOA_NV  =  71  &
    , TERPOOH     =  72  &
    , TERPPeroxy  =  73  &
    , VBS_TEST    =  74  &
    , SHIPNOX     =  75  &
    , Dust_ROAD_f =  76  &
    , Dust_ROAD_c =  77  &
    , Dust_WB_f   =  78  &
    , Dust_WB_c   =  79  &
    , Dust_SAH_f  =  80
  
  integer, public, parameter :: &
      Dust_SAH_c  =  81  &
    , POLLEN_BIRCH=  82  &
    , POLLEN_OLIVE=  83  &
    , POLLEN_RWEED=  84  &
    , POLLEN_GRASS=  85  &
    , ASOC_ng100  =  86  &
    , ASOC_ug1    =  87  &
    , ASOC_ug10   =  88  &
    , ASOC_ug1e2  =  89  &
    , ASOC_ug1e3  =  90
  
  integer, public, parameter :: &
      NON_C_ASOA_ng100=  91  &
    , NON_C_ASOA_ug1=  92  &
    , NON_C_ASOA_ug10=  93  &
    , NON_C_ASOA_ug1e2=  94  &
    , NON_C_ASOA_ug1e3=  95  &
    , BSOC_ng100  =  96  &
    , BSOC_ug1    =  97  &
    , BSOC_ug10   =  98  &
    , BSOC_ug1e2  =  99  &
    , BSOC_ug1e3  = 100
  
  integer, public, parameter :: &
      NON_C_BSOA_ng100= 101  &
    , NON_C_BSOA_ug1= 102  &
    , NON_C_BSOA_ug10= 103  &
    , NON_C_BSOA_ug1e2= 104  &
    , NON_C_BSOA_ug1e3= 105  &
    , FFFUEL_ng10 = 106  &
    , WOODOA_ng10 = 107  &
    , SO4         = 108  &
    , NH3         = 109  &
    , NO3_f       = 110
  
  integer, public, parameter :: &
      NO3_c       = 111  &
    , NH4_f       = 112  &
    , POM_f_WOOD  = 113  &
    , POM_f_FFUEL = 114  &
    , POM_c_FFUEL = 115  &
    , EC_f_WOOD_new= 116  &
    , EC_f_WOOD_age= 117  &
    , EC_c_WOOD   = 118  &
    , EC_f_FFUEL_new= 119  &
    , EC_f_FFUEL_age= 120
  
  integer, public, parameter :: &
      EC_c_FFUEL  = 121  &
    , REMPPM25    = 122  &
    , REMPPM_c    = 123  &
    , OM25_BGND   = 124  &
    , OM25_p      = 125  &
    , Ash_f       = 126  &
    , Ash_c       = 127  &
    , FFIRE_CO    = 128  &
    , FFIRE_OM    = 129  &
    , FFIRE_BC    = 130
  
  integer, public, parameter :: &
      FFIRE_REMPPM25= 131  &
    , SeaSalt_f   = 132  &
    , SeaSalt_c   = 133
  
  !+ Defines indices for ADV : Advected species
  integer, public, parameter :: FIRST_ADV=19, &
                                 LAST_ADV=133
  
  integer, public, parameter :: &
      IXADV_NMVOC       =   1  &
    , IXADV_O3          =   2  &
    , IXADV_NO          =   3  &
    , IXADV_NO2         =   4  &
    , IXADV_PAN         =   5  &
    , IXADV_MPAN        =   6  &
    , IXADV_NO3         =   7  &
    , IXADV_N2O5        =   8  &
    , IXADV_HNO3        =   9  &
    , IXADV_HONO        =  10
  
  integer, public, parameter :: &
      IXADV_HO2NO2      =  11  &
    , IXADV_MACR        =  12  &
    , IXADV_ISON        =  13  &
    , IXADV_GLYOX       =  14  &
    , IXADV_MGLYOX      =  15  &
    , IXADV_MAL         =  16  &
    , IXADV_MEK         =  17  &
    , IXADV_HCHO        =  18  &
    , IXADV_CH3CHO      =  19  &
    , IXADV_C2H6        =  20
  
  integer, public, parameter :: &
      IXADV_NC4H10      =  21  &
    , IXADV_C2H4        =  22  &
    , IXADV_C3H6        =  23  &
    , IXADV_OXYL        =  24  &
    , IXADV_C5H8        =  25  &
    , IXADV_APINENE     =  26  &
    , IXADV_CH3OOH      =  27  &
    , IXADV_C2H5OOH     =  28  &
    , IXADV_BURO2H      =  29  &
    , IXADV_ETRO2H      =  30
  
  integer, public, parameter :: &
      IXADV_PRRO2H      =  31  &
    , IXADV_OXYO2H      =  32  &
    , IXADV_MEKO2H      =  33  &
    , IXADV_MALO2H      =  34  &
    , IXADV_MACROOH     =  35  &
    , IXADV_ISRO2H      =  36  &
    , IXADV_H2O2        =  37  &
    , IXADV_CH3CO3H     =  38  &
    , IXADV_CH3OH       =  39  &
    , IXADV_C2H5OH      =  40
  
  integer, public, parameter :: &
      IXADV_ACETOL      =  41  &
    , IXADV_NALD        =  42  &
    , IXADV_HPALD       =  43  &
    , IXADV_PACALD      =  44  &
    , IXADV_IEPOX       =  45  &
    , IXADV_H2          =  46  &
    , IXADV_CO          =  47  &
    , IXADV_CH4         =  48  &
    , IXADV_SO2         =  49  &
    , IXADV_MVK         =  50
  
  integer, public, parameter :: &
      IXADV_BPINENE     =  51  &
    , IXADV_XTERP       =  52  &
    , IXADV_SQT_SOA_NV  =  53  &
    , IXADV_TERPOOH     =  54  &
    , IXADV_TERPPeroxy  =  55  &
    , IXADV_VBS_TEST    =  56  &
    , IXADV_SHIPNOX     =  57  &
    , IXADV_Dust_ROAD_f =  58  &
    , IXADV_Dust_ROAD_c =  59  &
    , IXADV_Dust_WB_f   =  60
  
  integer, public, parameter :: &
      IXADV_Dust_WB_c   =  61  &
    , IXADV_Dust_SAH_f  =  62  &
    , IXADV_Dust_SAH_c  =  63  &
    , IXADV_POLLEN_BIRCH=  64  &
    , IXADV_POLLEN_OLIVE=  65  &
    , IXADV_POLLEN_RWEED=  66  &
    , IXADV_POLLEN_GRASS=  67  &
    , IXADV_ASOC_ng100  =  68  &
    , IXADV_ASOC_ug1    =  69  &
    , IXADV_ASOC_ug10   =  70
  
  integer, public, parameter :: &
      IXADV_ASOC_ug1e2  =  71  &
    , IXADV_ASOC_ug1e3  =  72  &
    , IXADV_NON_C_ASOA_ng100=  73  &
    , IXADV_NON_C_ASOA_ug1=  74  &
    , IXADV_NON_C_ASOA_ug10=  75  &
    , IXADV_NON_C_ASOA_ug1e2=  76  &
    , IXADV_NON_C_ASOA_ug1e3=  77  &
    , IXADV_BSOC_ng100  =  78  &
    , IXADV_BSOC_ug1    =  79  &
    , IXADV_BSOC_ug10   =  80
  
  integer, public, parameter :: &
      IXADV_BSOC_ug1e2  =  81  &
    , IXADV_BSOC_ug1e3  =  82  &
    , IXADV_NON_C_BSOA_ng100=  83  &
    , IXADV_NON_C_BSOA_ug1=  84  &
    , IXADV_NON_C_BSOA_ug10=  85  &
    , IXADV_NON_C_BSOA_ug1e2=  86  &
    , IXADV_NON_C_BSOA_ug1e3=  87  &
    , IXADV_FFFUEL_ng10 =  88  &
    , IXADV_WOODOA_ng10 =  89  &
    , IXADV_SO4         =  90
  
  integer, public, parameter :: &
      IXADV_NH3         =  91  &
    , IXADV_NO3_f       =  92  &
    , IXADV_NO3_c       =  93  &
    , IXADV_NH4_f       =  94  &
    , IXADV_POM_f_WOOD  =  95  &
    , IXADV_POM_f_FFUEL =  96  &
    , IXADV_POM_c_FFUEL =  97  &
    , IXADV_EC_f_WOOD_new=  98  &
    , IXADV_EC_f_WOOD_age=  99  &
    , IXADV_EC_c_WOOD   = 100
  
  integer, public, parameter :: &
      IXADV_EC_f_FFUEL_new= 101  &
    , IXADV_EC_f_FFUEL_age= 102  &
    , IXADV_EC_c_FFUEL  = 103  &
    , IXADV_REMPPM25    = 104  &
    , IXADV_REMPPM_c    = 105  &
    , IXADV_OM25_BGND   = 106  &
    , IXADV_OM25_p      = 107  &
    , IXADV_Ash_f       = 108  &
    , IXADV_Ash_c       = 109  &
    , IXADV_FFIRE_CO    = 110
  
  integer, public, parameter :: &
      IXADV_FFIRE_OM    = 111  &
    , IXADV_FFIRE_BC    = 112  &
    , IXADV_FFIRE_REMPPM25= 113  &
    , IXADV_SeaSalt_f   = 114  &
    , IXADV_SeaSalt_c   = 115
  
  !+ Defines indices for SHL : Short-lived (non-advected) species
  integer, public, parameter :: FIRST_SHL=1, &
                                 LAST_SHL=18
  
  integer, public, parameter :: &
      IXSHL_RO2POOL     =   1  &
    , IXSHL_OD          =   2  &
    , IXSHL_OP          =   3  &
    , IXSHL_OH          =   4  &
    , IXSHL_HO2         =   5  &
    , IXSHL_CH3O2       =   6  &
    , IXSHL_C2H5O2      =   7  &
    , IXSHL_SECC4H9O2   =   8  &
    , IXSHL_ISRO2       =   9  &
    , IXSHL_ETRO2       =  10
  
  integer, public, parameter :: &
      IXSHL_PRRO2       =  11  &
    , IXSHL_OXYO2       =  12  &
    , IXSHL_MEKO2       =  13  &
    , IXSHL_MALO2       =  14  &
    , IXSHL_MACRO2      =  15  &
    , IXSHL_CH3CO3      =  16  &
    , IXSHL_TERPO2      =  17  &
    , IXSHL_XMTO3_RO2   =  18
  
  !+ Defines indices for SEMIVOL : Semi-volatile organic aerosols
  integer, public, parameter :: FIRST_SEMIVOL=86, &
                                 LAST_SEMIVOL=107
  
  integer, public, parameter :: &
      IXSOA_ASOC_ng100  =   1  &
    , IXSOA_ASOC_ug1    =   2  &
    , IXSOA_ASOC_ug10   =   3  &
    , IXSOA_ASOC_ug1e2  =   4  &
    , IXSOA_ASOC_ug1e3  =   5  &
    , IXSOA_NON_C_ASOA_ng100=   6  &
    , IXSOA_NON_C_ASOA_ug1=   7  &
    , IXSOA_NON_C_ASOA_ug10=   8  &
    , IXSOA_NON_C_ASOA_ug1e2=   9  &
    , IXSOA_NON_C_ASOA_ug1e3=  10
  
  integer, public, parameter :: &
      IXSOA_BSOC_ng100  =  11  &
    , IXSOA_BSOC_ug1    =  12  &
    , IXSOA_BSOC_ug10   =  13  &
    , IXSOA_BSOC_ug1e2  =  14  &
    , IXSOA_BSOC_ug1e3  =  15  &
    , IXSOA_NON_C_BSOA_ng100=  16  &
    , IXSOA_NON_C_BSOA_ug1=  17  &
    , IXSOA_NON_C_BSOA_ug10=  18  &
    , IXSOA_NON_C_BSOA_ug1e2=  19  &
    , IXSOA_NON_C_BSOA_ug1e3=  20
  
  integer, public, parameter :: &
      IXSOA_FFFUEL_ng10 =  21  &
    , IXSOA_WOODOA_ng10 =  22
  
  !/--   Characteristics of species:
  !/--   Number, name, molwt, carbon num, nmhc (1) or not(0)
  
  public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc
  
  type, public :: Chemical
       character(len=20) :: name
       real              :: molwt
       integer           :: nmhc      ! nmhc (1) or not(0)
       real              :: carbons   ! Carbon-number
       real              :: nitrogens ! Nitrogen-number
       real              :: sulphurs  ! Sulphur-number
  endtype Chemical
  type(Chemical), public, dimension(NSPEC_TOT), target :: species
  
  ! Pointers to parts of species (e.g. short-lived, advected)
  type(Chemical), public, dimension(:), pointer :: species_adv=>null()
  type(Chemical), public, dimension(:), pointer :: species_shl=>null()
  type(Chemical), public, dimension(:), pointer :: species_semivol=>null()

contains
  subroutine define_chemicals()
    integer :: istart ! For NAG compliance
    !+
    ! Pointers to parts of species (e.g. short-lived, advected), only assigned if
    ! non-empty.
    !
    istart = 1
    if (NSPEC_ADV > 0) then
      if( FIRST_ADV > 0 ) istart = FIRST_ADV
      species_adv => species(istart:LAST_ADV)
    end if
    istart = 1
    if (NSPEC_SHL > 0) then
      if( FIRST_SHL > 0 ) istart = FIRST_SHL
      species_shl => species(istart:LAST_SHL)
    end if
    istart = 1
    if (NSPEC_SEMIVOL > 0) then
      if( FIRST_SEMIVOL > 0 ) istart = FIRST_SEMIVOL
      species_semivol => species(istart:LAST_SEMIVOL)
    end if
    
    !+
    ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
    ! array, using indices from total list of species (advected + short-lived).
    !                                                      MW  NM   C    N   S
    species(RO2POOL     ) = Chemical("RO2POOL     ",  64.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(OD          ) = Chemical("OD          ",  16.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(OP          ) = Chemical("OP          ",  16.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(OH          ) = Chemical("OH          ",  17.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(HO2         ) = Chemical("HO2         ",  33.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(CH3O2       ) = Chemical("CH3O2       ",  47.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(C2H5O2      ) = Chemical("C2H5O2      ",  61.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(SECC4H9O2   ) = Chemical("SECC4H9O2   ",  89.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(ISRO2       ) = Chemical("ISRO2       ", 101.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(ETRO2       ) = Chemical("ETRO2       ",  77.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(PRRO2       ) = Chemical("PRRO2       ",  91.0000,  0,   3.0000,   0.0000,   0.0000 )
    species(OXYO2       ) = Chemical("OXYO2       ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(MEKO2       ) = Chemical("MEKO2       ", 103.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(MALO2       ) = Chemical("MALO2       ", 147.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(MACRO2      ) = Chemical("MACRO2      ", 119.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(CH3CO3      ) = Chemical("CH3CO3      ",  75.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(TERPO2      ) = Chemical("TERPO2      ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(XMTO3_RO2   ) = Chemical("XMTO3_RO2   ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NMVOC       ) = Chemical("NMVOC       ",  36.7670,  0,   0.0000,   0.0000,   0.0000 )
    species(O3          ) = Chemical("O3          ",  48.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NO          ) = Chemical("NO          ",  30.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NO2         ) = Chemical("NO2         ",  46.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(PAN         ) = Chemical("PAN         ", 121.0000,  0,   2.0000,   1.0000,   0.0000 )
    species(MPAN        ) = Chemical("MPAN        ", 132.0000,  0,   4.0000,   1.0000,   0.0000 )
    species(NO3         ) = Chemical("NO3         ",  62.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(N2O5        ) = Chemical("N2O5        ", 108.0000,  0,   0.0000,   2.0000,   0.0000 )
    species(HNO3        ) = Chemical("HNO3        ",  63.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(HONO        ) = Chemical("HONO        ",  47.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(HO2NO2      ) = Chemical("HO2NO2      ",  79.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(MACR        ) = Chemical("MACR        ",  70.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(ISON        ) = Chemical("ISON        ", 147.0000,  0,   5.0000,   1.0000,   0.0000 )
    species(GLYOX       ) = Chemical("GLYOX       ",  58.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(MGLYOX      ) = Chemical("MGLYOX      ",  72.0000,  0,   3.0000,   0.0000,   0.0000 )
    species(MAL         ) = Chemical("MAL         ",  98.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(MEK         ) = Chemical("MEK         ",  72.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(HCHO        ) = Chemical("HCHO        ",  30.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(CH3CHO      ) = Chemical("CH3CHO      ",  44.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(C2H6        ) = Chemical("C2H6        ",  30.0000,  1,   2.0000,   0.0000,   0.0000 )
    species(NC4H10      ) = Chemical("NC4H10      ",  58.0000,  1,   4.0000,   0.0000,   0.0000 )
    species(C2H4        ) = Chemical("C2H4        ",  28.0000,  1,   2.0000,   0.0000,   0.0000 )
    species(C3H6        ) = Chemical("C3H6        ",  42.0000,  1,   3.0000,   0.0000,   0.0000 )
    species(OXYL        ) = Chemical("OXYL        ", 106.0000,  1,   8.0000,   0.0000,   0.0000 )
    species(C5H8        ) = Chemical("C5H8        ",  68.0000,  1,   5.0000,   0.0000,   0.0000 )
    species(APINENE     ) = Chemical("APINENE     ", 136.0000,  1,  10.0000,   0.0000,   0.0000 )
    species(CH3OOH      ) = Chemical("CH3OOH      ",  48.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(C2H5OOH     ) = Chemical("C2H5OOH     ",  62.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(BURO2H      ) = Chemical("BURO2H      ",  90.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(ETRO2H      ) = Chemical("ETRO2H      ",  78.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(PRRO2H      ) = Chemical("PRRO2H      ",  92.0000,  0,   3.0000,   0.0000,   0.0000 )
    species(OXYO2H      ) = Chemical("OXYO2H      ",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(MEKO2H      ) = Chemical("MEKO2H      ", 104.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(MALO2H      ) = Chemical("MALO2H      ", 147.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(MACROOH     ) = Chemical("MACROOH     ", 120.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(ISRO2H      ) = Chemical("ISRO2H      ",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(H2O2        ) = Chemical("H2O2        ",  34.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(CH3CO3H     ) = Chemical("CH3CO3H     ",  76.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(CH3OH       ) = Chemical("CH3OH       ",  32.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(C2H5OH      ) = Chemical("C2H5OH      ",  46.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(ACETOL      ) = Chemical("ACETOL      ",  74.0000,  0,   3.0000,   0.0000,   0.0000 )
    species(NALD        ) = Chemical("NALD        ", 105.0000,  0,   2.0000,   1.0000,   0.0000 )
    species(HPALD       ) = Chemical("HPALD       ", 116.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(PACALD      ) = Chemical("PACALD      ", 130.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(IEPOX       ) = Chemical("IEPOX       ", 118.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(H2          ) = Chemical("H2          ",   2.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(CO          ) = Chemical("CO          ",  28.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(CH4         ) = Chemical("CH4         ",  16.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(SO2         ) = Chemical("SO2         ",  64.0000,  0,   0.0000,   0.0000,   1.0000 )
    species(MVK         ) = Chemical("MVK         ",  70.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(BPINENE     ) = Chemical("BPINENE     ", 136.0000,  1,  10.0000,   0.0000,   0.0000 )
    species(XTERP       ) = Chemical("XTERP       ", 136.0000,  1,  10.0000,   0.0000,   0.0000 )
    species(SQT_SOA_NV  ) = Chemical("SQT_SOA_NV  ", 302.0000,  0,  14.0000,   0.0000,   0.0000 )
    species(TERPOOH     ) = Chemical("TERPOOH     ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(TERPPeroxy  ) = Chemical("TERPPeroxy  ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(VBS_TEST    ) = Chemical("VBS_TEST    ",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SHIPNOX     ) = Chemical("SHIPNOX     ",  46.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(Dust_ROAD_f ) = Chemical("Dust_ROAD_f ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_ROAD_c ) = Chemical("Dust_ROAD_c ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_WB_f   ) = Chemical("Dust_WB_f   ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_WB_c   ) = Chemical("Dust_WB_c   ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_SAH_f  ) = Chemical("Dust_SAH_f  ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_SAH_c  ) = Chemical("Dust_SAH_c  ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_BIRCH) = Chemical("POLLEN_BIRCH",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_OLIVE) = Chemical("POLLEN_OLIVE",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_RWEED) = Chemical("POLLEN_RWEED",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_GRASS) = Chemical("POLLEN_GRASS",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(ASOC_ng100  ) = Chemical("ASOC_ng100  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ASOC_ug1    ) = Chemical("ASOC_ug1    ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ASOC_ug10   ) = Chemical("ASOC_ug10   ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ASOC_ug1e2  ) = Chemical("ASOC_ug1e2  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ASOC_ug1e3  ) = Chemical("ASOC_ug1e3  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(NON_C_ASOA_ng100) = Chemical("NON_C_ASOA_ng100",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NON_C_ASOA_ug1) = Chemical("NON_C_ASOA_ug1",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NON_C_ASOA_ug10) = Chemical("NON_C_ASOA_ug10",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NON_C_ASOA_ug1e2) = Chemical("NON_C_ASOA_ug1e2",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NON_C_ASOA_ug1e3) = Chemical("NON_C_ASOA_ug1e3",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(BSOC_ng100  ) = Chemical("BSOC_ng100  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(BSOC_ug1    ) = Chemical("BSOC_ug1    ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(BSOC_ug10   ) = Chemical("BSOC_ug10   ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(BSOC_ug1e2  ) = Chemical("BSOC_ug1e2  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(BSOC_ug1e3  ) = Chemical("BSOC_ug1e3  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(NON_C_BSOA_ng100) = Chemical("NON_C_BSOA_ng100",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NON_C_BSOA_ug1) = Chemical("NON_C_BSOA_ug1",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NON_C_BSOA_ug10) = Chemical("NON_C_BSOA_ug10",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NON_C_BSOA_ug1e2) = Chemical("NON_C_BSOA_ug1e2",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NON_C_BSOA_ug1e3) = Chemical("NON_C_BSOA_ug1e3",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(FFFUEL_ng10 ) = Chemical("FFFUEL_ng10 ",  15.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(WOODOA_ng10 ) = Chemical("WOODOA_ng10 ",  20.4000,  0,   1.0000,   0.0000,   0.0000 )
    species(SO4         ) = Chemical("SO4         ",  96.0000,  0,   0.0000,   0.0000,   1.0000 )
    species(NH3         ) = Chemical("NH3         ",  17.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NO3_f       ) = Chemical("NO3_f       ",  62.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NO3_c       ) = Chemical("NO3_c       ",  62.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NH4_f       ) = Chemical("NH4_f       ",  18.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(POM_f_WOOD  ) = Chemical("POM_f_WOOD  ",  20.4000,  0,   1.0000,   0.0000,   0.0000 )
    species(POM_f_FFUEL ) = Chemical("POM_f_FFUEL ",  15.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(POM_c_FFUEL ) = Chemical("POM_c_FFUEL ",  15.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_f_WOOD_new) = Chemical("EC_f_WOOD_new",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_f_WOOD_age) = Chemical("EC_f_WOOD_age",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_c_WOOD   ) = Chemical("EC_c_WOOD   ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_f_FFUEL_new) = Chemical("EC_f_FFUEL_new",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_f_FFUEL_age) = Chemical("EC_f_FFUEL_age",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_c_FFUEL  ) = Chemical("EC_c_FFUEL  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(REMPPM25    ) = Chemical("REMPPM25    ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(REMPPM_c    ) = Chemical("REMPPM_c    ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(OM25_BGND   ) = Chemical("OM25_BGND   ",  24.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(OM25_p      ) = Chemical("OM25_p      ",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Ash_f       ) = Chemical("Ash_f       ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Ash_c       ) = Chemical("Ash_c       ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(FFIRE_CO    ) = Chemical("FFIRE_CO    ",  28.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(FFIRE_OM    ) = Chemical("FFIRE_OM    ",  20.4000,  0,   1.0000,   0.0000,   0.0000 )
    species(FFIRE_BC    ) = Chemical("FFIRE_BC    ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(FFIRE_REMPPM25) = Chemical("FFIRE_REMPPM25",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SeaSalt_f   ) = Chemical("SeaSalt_f   ",  58.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SeaSalt_c   ) = Chemical("SeaSalt_c   ",  58.0000,  0,   0.0000,   0.0000,   0.0000 )
  end subroutine define_chemicals

end module ChemSpecs_mod
