! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem16x ShipNOx BVOC_EmChem16x Aqueous_EmChem16x Aero2017nx VBS_acp2012 FFireInert SeaSalt DustExtended Ash_EmChem16x
module ChemSpecs_mod

  use ChemDims_mod      ! => NSPEC_TOT, NCHEMRATES, ....
  
  implicit none
  private
  character(len=*),parameter, public :: CM_schemes_ChemSpecs = " EmChem16x ShipNOx BVOC_EmChem16x Aqueous_EmChem16x Aero2017nx VBS_acp2012 FFireInert SeaSalt DustExtended Ash_EmChem16x"
  
  
  integer, public, parameter :: &
      OD          =   1  &
    , OP          =   2  &
    , OH          =   3  &
    , HO2         =   4  &
    , CH3O2       =   5  &
    , C2H5O2      =   6  &
    , SECC4H9O2   =   7  &
    , ISRO2       =   8  &
    , ETRO2       =   9  &
    , PRRO2       =  10
  
  integer, public, parameter :: &
      OXYO2       =  11  &
    , MEKO2       =  12  &
    , MALO2       =  13  &
    , MACRO2      =  14  &
    , CH3CO3      =  15  &
    , TERPO2      =  16  &
    , XMTO3_RO2   =  17  &
    , NMVOC       =  18  &
    , O3          =  19  &
    , NO          =  20
  
  integer, public, parameter :: &
      NO2         =  21  &
    , PAN         =  22  &
    , MPAN        =  23  &
    , NO3         =  24  &
    , N2O5        =  25  &
    , HNO3        =  26  &
    , HONO        =  27  &
    , HO2NO2      =  28  &
    , MACR        =  29  &
    , ISON        =  30
  
  integer, public, parameter :: &
      GLYOX       =  31  &
    , MGLYOX      =  32  &
    , MAL         =  33  &
    , MEK         =  34  &
    , HCHO        =  35  &
    , CH3CHO      =  36  &
    , C2H6        =  37  &
    , NC4H10      =  38  &
    , C2H4        =  39  &
    , C3H6        =  40
  
  integer, public, parameter :: &
      OXYL        =  41  &
    , C5H8        =  42  &
    , APINENE     =  43  &
    , CH3OOH      =  44  &
    , C2H5OOH     =  45  &
    , BURO2H      =  46  &
    , ETRO2H      =  47  &
    , PRRO2H      =  48  &
    , OXYO2H      =  49  &
    , MEKO2H      =  50
  
  integer, public, parameter :: &
      MALO2H      =  51  &
    , MACROOH     =  52  &
    , ISRO2H      =  53  &
    , H2O2        =  54  &
    , CH3CO3H     =  55  &
    , CH3OH       =  56  &
    , C2H5OH      =  57  &
    , ACETOL      =  58  &
    , NALD        =  59  &
    , HPALD       =  60
  
  integer, public, parameter :: &
      PACALD      =  61  &
    , IEPOX       =  62  &
    , H2          =  63  &
    , CO          =  64  &
    , CH4         =  65  &
    , SO2         =  66  &
    , SHIPNOX     =  67  &
    , MVK         =  68  &
    , BPINENE     =  69  &
    , XTERP       =  70
  
  integer, public, parameter :: &
      SQT_SOA_NV  =  71  &
    , TERPOOH     =  72  &
    , TERPPeroxy  =  73  &
    , VBS_TEST    =  74  &
    , Dust_ROAD_f =  75  &
    , Dust_ROAD_c =  76  &
    , Dust_WB_f   =  77  &
    , Dust_WB_c   =  78  &
    , Dust_SAH_f  =  79  &
    , Dust_SAH_c  =  80
  
  integer, public, parameter :: &
      ASOC_ng100  =  81  &
    , ASOC_ug1    =  82  &
    , ASOC_ug10   =  83  &
    , ASOC_ug1e2  =  84  &
    , ASOC_ug1e3  =  85  &
    , NON_C_ASOA_ng100=  86  &
    , NON_C_ASOA_ug1=  87  &
    , NON_C_ASOA_ug10=  88  &
    , NON_C_ASOA_ug1e2=  89  &
    , NON_C_ASOA_ug1e3=  90
  
  integer, public, parameter :: &
      BSOC_ng100  =  91  &
    , BSOC_ug1    =  92  &
    , BSOC_ug10   =  93  &
    , BSOC_ug1e2  =  94  &
    , BSOC_ug1e3  =  95  &
    , NON_C_BSOA_ng100=  96  &
    , NON_C_BSOA_ug1=  97  &
    , NON_C_BSOA_ug10=  98  &
    , NON_C_BSOA_ug1e2=  99  &
    , NON_C_BSOA_ug1e3= 100
  
  integer, public, parameter :: &
      FFFUEL_ng10 = 101  &
    , WOODOA_ng10 = 102  &
    , SO4         = 103  &
    , NH3         = 104  &
    , NO3_f       = 105  &
    , NO3_c       = 106  &
    , NH4_f       = 107  &
    , POM_f_WOOD  = 108  &
    , POM_f_FFUEL = 109  &
    , POM_c_FFUEL = 110
  
  integer, public, parameter :: &
      EC_f_WOOD_new= 111  &
    , EC_f_WOOD_age= 112  &
    , EC_c_WOOD   = 113  &
    , EC_f_FFUEL_new= 114  &
    , EC_f_FFUEL_age= 115  &
    , EC_c_FFUEL  = 116  &
    , REMPPM25    = 117  &
    , REMPPM_c    = 118  &
    , OM25_BGND   = 119  &
    , OM25_p      = 120
  
  integer, public, parameter :: &
      FFIRE_CO    = 121  &
    , FFIRE_OM    = 122  &
    , FFIRE_BC    = 123  &
    , FFIRE_REMPPM25= 124  &
    , SeaSalt_f   = 125  &
    , SeaSalt_c   = 126  &
    , Ash_f       = 127  &
    , Ash_c       = 128
  
  !+ Defines indices for ADV : Advected species
  integer, public, parameter :: FIRST_ADV=18, &
                                 LAST_ADV=128
  
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
    , IXADV_SHIPNOX     =  50
  
  integer, public, parameter :: &
      IXADV_MVK         =  51  &
    , IXADV_BPINENE     =  52  &
    , IXADV_XTERP       =  53  &
    , IXADV_SQT_SOA_NV  =  54  &
    , IXADV_TERPOOH     =  55  &
    , IXADV_TERPPeroxy  =  56  &
    , IXADV_VBS_TEST    =  57  &
    , IXADV_Dust_ROAD_f =  58  &
    , IXADV_Dust_ROAD_c =  59  &
    , IXADV_Dust_WB_f   =  60
  
  integer, public, parameter :: &
      IXADV_Dust_WB_c   =  61  &
    , IXADV_Dust_SAH_f  =  62  &
    , IXADV_Dust_SAH_c  =  63  &
    , IXADV_ASOC_ng100  =  64  &
    , IXADV_ASOC_ug1    =  65  &
    , IXADV_ASOC_ug10   =  66  &
    , IXADV_ASOC_ug1e2  =  67  &
    , IXADV_ASOC_ug1e3  =  68  &
    , IXADV_NON_C_ASOA_ng100=  69  &
    , IXADV_NON_C_ASOA_ug1=  70
  
  integer, public, parameter :: &
      IXADV_NON_C_ASOA_ug10=  71  &
    , IXADV_NON_C_ASOA_ug1e2=  72  &
    , IXADV_NON_C_ASOA_ug1e3=  73  &
    , IXADV_BSOC_ng100  =  74  &
    , IXADV_BSOC_ug1    =  75  &
    , IXADV_BSOC_ug10   =  76  &
    , IXADV_BSOC_ug1e2  =  77  &
    , IXADV_BSOC_ug1e3  =  78  &
    , IXADV_NON_C_BSOA_ng100=  79  &
    , IXADV_NON_C_BSOA_ug1=  80
  
  integer, public, parameter :: &
      IXADV_NON_C_BSOA_ug10=  81  &
    , IXADV_NON_C_BSOA_ug1e2=  82  &
    , IXADV_NON_C_BSOA_ug1e3=  83  &
    , IXADV_FFFUEL_ng10 =  84  &
    , IXADV_WOODOA_ng10 =  85  &
    , IXADV_SO4         =  86  &
    , IXADV_NH3         =  87  &
    , IXADV_NO3_f       =  88  &
    , IXADV_NO3_c       =  89  &
    , IXADV_NH4_f       =  90
  
  integer, public, parameter :: &
      IXADV_POM_f_WOOD  =  91  &
    , IXADV_POM_f_FFUEL =  92  &
    , IXADV_POM_c_FFUEL =  93  &
    , IXADV_EC_f_WOOD_new=  94  &
    , IXADV_EC_f_WOOD_age=  95  &
    , IXADV_EC_c_WOOD   =  96  &
    , IXADV_EC_f_FFUEL_new=  97  &
    , IXADV_EC_f_FFUEL_age=  98  &
    , IXADV_EC_c_FFUEL  =  99  &
    , IXADV_REMPPM25    = 100
  
  integer, public, parameter :: &
      IXADV_REMPPM_c    = 101  &
    , IXADV_OM25_BGND   = 102  &
    , IXADV_OM25_p      = 103  &
    , IXADV_FFIRE_CO    = 104  &
    , IXADV_FFIRE_OM    = 105  &
    , IXADV_FFIRE_BC    = 106  &
    , IXADV_FFIRE_REMPPM25= 107  &
    , IXADV_SeaSalt_f   = 108  &
    , IXADV_SeaSalt_c   = 109  &
    , IXADV_Ash_f       = 110
  
  integer, public, parameter :: &
      IXADV_Ash_c       = 111
  
  !+ Defines indices for SHL : Short-lived (non-advected) species
  integer, public, parameter :: FIRST_SHL=1, &
                                 LAST_SHL=17
  
  integer, public, parameter :: &
      IXSHL_OD          =   1  &
    , IXSHL_OP          =   2  &
    , IXSHL_OH          =   3  &
    , IXSHL_HO2         =   4  &
    , IXSHL_CH3O2       =   5  &
    , IXSHL_C2H5O2      =   6  &
    , IXSHL_SECC4H9O2   =   7  &
    , IXSHL_ISRO2       =   8  &
    , IXSHL_ETRO2       =   9  &
    , IXSHL_PRRO2       =  10
  
  integer, public, parameter :: &
      IXSHL_OXYO2       =  11  &
    , IXSHL_MEKO2       =  12  &
    , IXSHL_MALO2       =  13  &
    , IXSHL_MACRO2      =  14  &
    , IXSHL_CH3CO3      =  15  &
    , IXSHL_TERPO2      =  16  &
    , IXSHL_XMTO3_RO2   =  17
  
  !+ Defines indices for SEMIVOL : Semi-volatile organic aerosols
  integer, public, parameter :: FIRST_SEMIVOL=81, &
                                 LAST_SEMIVOL=102
  
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
    species(SHIPNOX     ) = Chemical("SHIPNOX     ",  46.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(MVK         ) = Chemical("MVK         ",  70.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(BPINENE     ) = Chemical("BPINENE     ", 136.0000,  1,  10.0000,   0.0000,   0.0000 )
    species(XTERP       ) = Chemical("XTERP       ", 136.0000,  1,  10.0000,   0.0000,   0.0000 )
    species(SQT_SOA_NV  ) = Chemical("SQT_SOA_NV  ", 202.0000,  1,  15.0000,   0.0000,   0.0000 )
    species(TERPOOH     ) = Chemical("TERPOOH     ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(TERPPeroxy  ) = Chemical("TERPPeroxy  ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(VBS_TEST    ) = Chemical("VBS_TEST    ",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_ROAD_f ) = Chemical("Dust_ROAD_f ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_ROAD_c ) = Chemical("Dust_ROAD_c ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_WB_f   ) = Chemical("Dust_WB_f   ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_WB_c   ) = Chemical("Dust_WB_c   ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_SAH_f  ) = Chemical("Dust_SAH_f  ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_SAH_c  ) = Chemical("Dust_SAH_c  ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
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
    species(FFIRE_CO    ) = Chemical("FFIRE_CO    ",  28.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(FFIRE_OM    ) = Chemical("FFIRE_OM    ",  20.4000,  0,   1.0000,   0.0000,   0.0000 )
    species(FFIRE_BC    ) = Chemical("FFIRE_BC    ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(FFIRE_REMPPM25) = Chemical("FFIRE_REMPPM25",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SeaSalt_f   ) = Chemical("SeaSalt_f   ",  58.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SeaSalt_c   ) = Chemical("SeaSalt_c   ",  58.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Ash_f       ) = Chemical("Ash_f       ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Ash_c       ) = Chemical("Ash_c       ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
  end subroutine define_chemicals

end module ChemSpecs_mod
