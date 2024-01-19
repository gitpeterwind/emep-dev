! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem19cj PM_VBS_EmChem19 Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt DustExtended Ash PM_ResNonResInert BVOC_SQT_NV BVOC_IsoMT1_emis Pollen FungalSpores
module ChemSpecs_mod

  use ChemDims_mod      ! => NSPEC_TOT, NCHEMRATES, ....
  
  implicit none
  private
  character(len=*),parameter, public :: CM_schemes_ChemSpecs = " EmChem19cj PM_VBS_EmChem19 Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt DustExtended Ash PM_ResNonResInert BVOC_SQT_NV BVOC_IsoMT1_emis Pollen FungalSpores"
  
  
  integer, public, parameter :: &
      OD          =   1  &
    , OP          =   2  &
    , OH          =   3  &
    , HO2         =   4  &
    , CH3O2       =   5  &
    , C2H5O2      =   6  &
    , C4H9O2      =   7  &
    , ISRO2       =   8  &
    , ETRO2       =   9  &
    , PRRO2       =  10
  
  integer, public, parameter :: &
      OXYO2       =  11  &
    , MEKO2       =  12  &
    , C5DICARBO2  =  13  &
    , MACRO2      =  14  &
    , TERPO2      =  15  &
    , RO2POOL     =  16  &
    , CH3CO3      =  17  &
    , O3          =  18  &
    , NO          =  19  &
    , NO2         =  20
  
  integer, public, parameter :: &
      NO3         =  21  &
    , N2O5        =  22  &
    , H2          =  23  &
    , H2O2        =  24  &
    , HONO        =  25  &
    , HNO3        =  26  &
    , HO2NO2      =  27  &
    , CO          =  28  &
    , CH4         =  29  &
    , C2H6        =  30
  
  integer, public, parameter :: &
      NC4H10      =  31  &
    , C2H4        =  32  &
    , C3H6        =  33  &
    , BENZENE     =  34  &
    , TOLUENE     =  35  &
    , OXYL        =  36  &
    , C5H8        =  37  &
    , CH3OH       =  38  &
    , C2H5OH      =  39  &
    , HCHO        =  40
  
  integer, public, parameter :: &
      CH3CHO      =  41  &
    , MACR        =  42  &
    , MEK         =  43  &
    , ACETOL      =  44  &
    , GLYOX       =  45  &
    , MGLYOX      =  46  &
    , BIACET      =  47  &
    , C5DICARB    =  48  &
    , CH3OOH      =  49  &
    , C2H5OOH     =  50
  
  integer, public, parameter :: &
      BURO2H      =  51  &
    , ETRO2H      =  52  &
    , MEKO2H      =  53  &
    , ISRO2H      =  54  &
    , C5DICAROOH  =  55  &
    , HPALD       =  56  &
    , MACROOH     =  57  &
    , OXYO2H      =  58  &
    , CH3CO3H     =  59  &
    , PACALD      =  60
  
  integer, public, parameter :: &
      IEPOX       =  61  &
    , SC4H9NO3    =  62  &
    , NALD        =  63  &
    , ISON        =  64  &
    , PAN         =  65  &
    , MPAN        =  66  &
    , APINENE     =  67  &
    , TERPOOH     =  68  &
    , SO2         =  69  &
    , shipNOx     =  70
  
  integer, public, parameter :: &
      Dust_road_f =  71  &
    , Dust_road_c =  72  &
    , Dust_wb_f   =  73  &
    , Dust_wb_c   =  74  &
    , Dust_sah_f  =  75  &
    , Dust_sah_c  =  76  &
    , SQT_SOA_NV  =  77  &
    , POLLEN_BIRCH=  78  &
    , POLLEN_OLIVE=  79  &
    , POLLEN_ALDER=  80
  
  integer, public, parameter :: &
      POLLEN_RWEED=  81  &
    , POLLEN_GRASS=  82  &
    , POLLEN_MUGWORT1=  83  &
    , POLLEN_MUGWORT2=  84  &
    , POLLEN_MUGWORT3=  85  &
    , POLLEN_MUGWORT4=  86  &
    , POLLEN_MUGWORT5=  87  &
    , ASOC_ng1e2  =  88  &
    , ASOC_ug1    =  89  &
    , ASOC_ug10   =  90
  
  integer, public, parameter :: &
      ASOC_ug1e2  =  91  &
    , ASOC_ug1e3  =  92  &
    , non_C_ASOA_ng1e2=  93  &
    , non_C_ASOA_ug1=  94  &
    , non_C_ASOA_ug10=  95  &
    , non_C_ASOA_ug1e2=  96  &
    , non_C_ASOA_ug1e3=  97  &
    , BSOC_ng1e2  =  98  &
    , BSOC_ug1    =  99  &
    , BSOC_ug10   = 100
  
  integer, public, parameter :: &
      BSOC_ug1e2  = 101  &
    , BSOC_ug1e3  = 102  &
    , non_C_BSOA_ng1e2= 103  &
    , non_C_BSOA_ug1= 104  &
    , non_C_BSOA_ug10= 105  &
    , non_C_BSOA_ug1e2= 106  &
    , non_C_BSOA_ug1e3= 107  &
    , SO4         = 108  &
    , NH3         = 109  &
    , NO3_f       = 110
  
  integer, public, parameter :: &
      NO3_c       = 111  &
    , NH4_f       = 112  &
    , OM25_bgnd   = 113  &
    , OM25_p      = 114  &
    , ffire_OM    = 115  &
    , ffire_BC    = 116  &
    , ffire_remPPM25= 117  &
    , ffire_c     = 118  &
    , SeaSalt_f   = 119  &
    , SeaSalt_c   = 120
  
  integer, public, parameter :: &
      Ash_f       = 121  &
    , Ash_c       = 122  &
    , POM_f_Res   = 123  &
    , POM_c_Res   = 124  &
    , POM_f_nonRes= 125  &
    , POM_c_nonRes= 126  &
    , EC_f_Res_new= 127  &
    , EC_f_Res_age= 128  &
    , EC_c_Res    = 129  &
    , EC_f_nonRes_new= 130
  
  integer, public, parameter :: &
      EC_f_nonRes_age= 131  &
    , EC_c_nonRes = 132  &
    , remPPM25_nonRes= 133  &
    , remPPM25_Res= 134  &
    , remPPM_c_nonRes= 135  &
    , remPPM_c_Res= 136  &
    , FUNGAL_SPORES= 137
  
  !+ Defines indices for ADV : Advected species
  integer, public, parameter :: FIRST_ADV=16, &
                                 LAST_ADV=137
  
  integer, public, parameter :: &
      IXADV_RO2POOL     =   1  &
    , IXADV_CH3CO3      =   2  &
    , IXADV_O3          =   3  &
    , IXADV_NO          =   4  &
    , IXADV_NO2         =   5  &
    , IXADV_NO3         =   6  &
    , IXADV_N2O5        =   7  &
    , IXADV_H2          =   8  &
    , IXADV_H2O2        =   9  &
    , IXADV_HONO        =  10
  
  integer, public, parameter :: &
      IXADV_HNO3        =  11  &
    , IXADV_HO2NO2      =  12  &
    , IXADV_CO          =  13  &
    , IXADV_CH4         =  14  &
    , IXADV_C2H6        =  15  &
    , IXADV_NC4H10      =  16  &
    , IXADV_C2H4        =  17  &
    , IXADV_C3H6        =  18  &
    , IXADV_BENZENE     =  19  &
    , IXADV_TOLUENE     =  20
  
  integer, public, parameter :: &
      IXADV_OXYL        =  21  &
    , IXADV_C5H8        =  22  &
    , IXADV_CH3OH       =  23  &
    , IXADV_C2H5OH      =  24  &
    , IXADV_HCHO        =  25  &
    , IXADV_CH3CHO      =  26  &
    , IXADV_MACR        =  27  &
    , IXADV_MEK         =  28  &
    , IXADV_ACETOL      =  29  &
    , IXADV_GLYOX       =  30
  
  integer, public, parameter :: &
      IXADV_MGLYOX      =  31  &
    , IXADV_BIACET      =  32  &
    , IXADV_C5DICARB    =  33  &
    , IXADV_CH3OOH      =  34  &
    , IXADV_C2H5OOH     =  35  &
    , IXADV_BURO2H      =  36  &
    , IXADV_ETRO2H      =  37  &
    , IXADV_MEKO2H      =  38  &
    , IXADV_ISRO2H      =  39  &
    , IXADV_C5DICAROOH  =  40
  
  integer, public, parameter :: &
      IXADV_HPALD       =  41  &
    , IXADV_MACROOH     =  42  &
    , IXADV_OXYO2H      =  43  &
    , IXADV_CH3CO3H     =  44  &
    , IXADV_PACALD      =  45  &
    , IXADV_IEPOX       =  46  &
    , IXADV_SC4H9NO3    =  47  &
    , IXADV_NALD        =  48  &
    , IXADV_ISON        =  49  &
    , IXADV_PAN         =  50
  
  integer, public, parameter :: &
      IXADV_MPAN        =  51  &
    , IXADV_APINENE     =  52  &
    , IXADV_TERPOOH     =  53  &
    , IXADV_SO2         =  54  &
    , IXADV_shipNOx     =  55  &
    , IXADV_Dust_road_f =  56  &
    , IXADV_Dust_road_c =  57  &
    , IXADV_Dust_wb_f   =  58  &
    , IXADV_Dust_wb_c   =  59  &
    , IXADV_Dust_sah_f  =  60
  
  integer, public, parameter :: &
      IXADV_Dust_sah_c  =  61  &
    , IXADV_SQT_SOA_NV  =  62  &
    , IXADV_POLLEN_BIRCH=  63  &
    , IXADV_POLLEN_OLIVE=  64  &
    , IXADV_POLLEN_ALDER=  65  &
    , IXADV_POLLEN_RWEED=  66  &
    , IXADV_POLLEN_GRASS=  67  &
    , IXADV_POLLEN_MUGWORT1=  68  &
    , IXADV_POLLEN_MUGWORT2=  69  &
    , IXADV_POLLEN_MUGWORT3=  70
  
  integer, public, parameter :: &
      IXADV_POLLEN_MUGWORT4=  71  &
    , IXADV_POLLEN_MUGWORT5=  72  &
    , IXADV_ASOC_ng1e2  =  73  &
    , IXADV_ASOC_ug1    =  74  &
    , IXADV_ASOC_ug10   =  75  &
    , IXADV_ASOC_ug1e2  =  76  &
    , IXADV_ASOC_ug1e3  =  77  &
    , IXADV_non_C_ASOA_ng1e2=  78  &
    , IXADV_non_C_ASOA_ug1=  79  &
    , IXADV_non_C_ASOA_ug10=  80
  
  integer, public, parameter :: &
      IXADV_non_C_ASOA_ug1e2=  81  &
    , IXADV_non_C_ASOA_ug1e3=  82  &
    , IXADV_BSOC_ng1e2  =  83  &
    , IXADV_BSOC_ug1    =  84  &
    , IXADV_BSOC_ug10   =  85  &
    , IXADV_BSOC_ug1e2  =  86  &
    , IXADV_BSOC_ug1e3  =  87  &
    , IXADV_non_C_BSOA_ng1e2=  88  &
    , IXADV_non_C_BSOA_ug1=  89  &
    , IXADV_non_C_BSOA_ug10=  90
  
  integer, public, parameter :: &
      IXADV_non_C_BSOA_ug1e2=  91  &
    , IXADV_non_C_BSOA_ug1e3=  92  &
    , IXADV_SO4         =  93  &
    , IXADV_NH3         =  94  &
    , IXADV_NO3_f       =  95  &
    , IXADV_NO3_c       =  96  &
    , IXADV_NH4_f       =  97  &
    , IXADV_OM25_bgnd   =  98  &
    , IXADV_OM25_p      =  99  &
    , IXADV_ffire_OM    = 100
  
  integer, public, parameter :: &
      IXADV_ffire_BC    = 101  &
    , IXADV_ffire_remPPM25= 102  &
    , IXADV_ffire_c     = 103  &
    , IXADV_SeaSalt_f   = 104  &
    , IXADV_SeaSalt_c   = 105  &
    , IXADV_Ash_f       = 106  &
    , IXADV_Ash_c       = 107  &
    , IXADV_POM_f_Res   = 108  &
    , IXADV_POM_c_Res   = 109  &
    , IXADV_POM_f_nonRes= 110
  
  integer, public, parameter :: &
      IXADV_POM_c_nonRes= 111  &
    , IXADV_EC_f_Res_new= 112  &
    , IXADV_EC_f_Res_age= 113  &
    , IXADV_EC_c_Res    = 114  &
    , IXADV_EC_f_nonRes_new= 115  &
    , IXADV_EC_f_nonRes_age= 116  &
    , IXADV_EC_c_nonRes = 117  &
    , IXADV_remPPM25_nonRes= 118  &
    , IXADV_remPPM25_Res= 119  &
    , IXADV_remPPM_c_nonRes= 120
  
  integer, public, parameter :: &
      IXADV_remPPM_c_Res= 121  &
    , IXADV_FUNGAL_SPORES= 122
  
  !+ Defines indices for SHL : Short-lived (non-advected) species
  integer, public, parameter :: FIRST_SHL=1, &
                                 LAST_SHL=15
  
  integer, public, parameter :: &
      IXSHL_OD          =   1  &
    , IXSHL_OP          =   2  &
    , IXSHL_OH          =   3  &
    , IXSHL_HO2         =   4  &
    , IXSHL_CH3O2       =   5  &
    , IXSHL_C2H5O2      =   6  &
    , IXSHL_C4H9O2      =   7  &
    , IXSHL_ISRO2       =   8  &
    , IXSHL_ETRO2       =   9  &
    , IXSHL_PRRO2       =  10
  
  integer, public, parameter :: &
      IXSHL_OXYO2       =  11  &
    , IXSHL_MEKO2       =  12  &
    , IXSHL_C5DICARBO2  =  13  &
    , IXSHL_MACRO2      =  14  &
    , IXSHL_TERPO2      =  15
  
  !+ Defines indices for SEMIVOL : Semi-volatile organic aerosols
  integer, public, parameter :: FIRST_SEMIVOL=88, &
                                 LAST_SEMIVOL=107
  
  integer, public, parameter :: &
      IXSOA_ASOC_ng1e2  =   1  &
    , IXSOA_ASOC_ug1    =   2  &
    , IXSOA_ASOC_ug10   =   3  &
    , IXSOA_ASOC_ug1e2  =   4  &
    , IXSOA_ASOC_ug1e3  =   5  &
    , IXSOA_non_C_ASOA_ng1e2=   6  &
    , IXSOA_non_C_ASOA_ug1=   7  &
    , IXSOA_non_C_ASOA_ug10=   8  &
    , IXSOA_non_C_ASOA_ug1e2=   9  &
    , IXSOA_non_C_ASOA_ug1e3=  10
  
  integer, public, parameter :: &
      IXSOA_BSOC_ng1e2  =  11  &
    , IXSOA_BSOC_ug1    =  12  &
    , IXSOA_BSOC_ug10   =  13  &
    , IXSOA_BSOC_ug1e2  =  14  &
    , IXSOA_BSOC_ug1e3  =  15  &
    , IXSOA_non_C_BSOA_ng1e2=  16  &
    , IXSOA_non_C_BSOA_ug1=  17  &
    , IXSOA_non_C_BSOA_ug10=  18  &
    , IXSOA_non_C_BSOA_ug1e2=  19  &
    , IXSOA_non_C_BSOA_ug1e3=  20
  
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
    species(C4H9O2      ) = Chemical("C4H9O2      ",  89.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(ISRO2       ) = Chemical("ISRO2       ", 101.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(ETRO2       ) = Chemical("ETRO2       ",  77.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(PRRO2       ) = Chemical("PRRO2       ",  91.0000,  0,   3.0000,   0.0000,   0.0000 )
    species(OXYO2       ) = Chemical("OXYO2       ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(MEKO2       ) = Chemical("MEKO2       ", 103.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(C5DICARBO2  ) = Chemical("C5DICARBO2  ", 147.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(MACRO2      ) = Chemical("MACRO2      ", 119.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(TERPO2      ) = Chemical("TERPO2      ",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(RO2POOL     ) = Chemical("RO2POOL     ",  64.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(CH3CO3      ) = Chemical("CH3CO3      ",  75.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(O3          ) = Chemical("O3          ",  48.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(NO          ) = Chemical("NO          ",  30.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NO2         ) = Chemical("NO2         ",  46.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NO3         ) = Chemical("NO3         ",  62.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(N2O5        ) = Chemical("N2O5        ", 108.0000,  0,   0.0000,   2.0000,   0.0000 )
    species(H2          ) = Chemical("H2          ",   2.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(H2O2        ) = Chemical("H2O2        ",  34.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(HONO        ) = Chemical("HONO        ",  47.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(HNO3        ) = Chemical("HNO3        ",  63.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(HO2NO2      ) = Chemical("HO2NO2      ",  79.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(CO          ) = Chemical("CO          ",  28.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(CH4         ) = Chemical("CH4         ",  16.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(C2H6        ) = Chemical("C2H6        ",  30.0000,  1,   2.0000,   0.0000,   0.0000 )
    species(NC4H10      ) = Chemical("NC4H10      ",  58.0000,  1,   4.0000,   0.0000,   0.0000 )
    species(C2H4        ) = Chemical("C2H4        ",  28.0000,  1,   2.0000,   0.0000,   0.0000 )
    species(C3H6        ) = Chemical("C3H6        ",  42.0000,  1,   3.0000,   0.0000,   0.0000 )
    species(BENZENE     ) = Chemical("BENZENE     ",  78.0000,  1,   6.0000,   0.0000,   0.0000 )
    species(TOLUENE     ) = Chemical("TOLUENE     ",  92.0000,  1,   7.0000,   0.0000,   0.0000 )
    species(OXYL        ) = Chemical("OXYL        ", 106.0000,  1,   8.0000,   0.0000,   0.0000 )
    species(C5H8        ) = Chemical("C5H8        ",  68.0000,  1,   5.0000,   0.0000,   0.0000 )
    species(CH3OH       ) = Chemical("CH3OH       ",  32.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(C2H5OH      ) = Chemical("C2H5OH      ",  46.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(HCHO        ) = Chemical("HCHO        ",  30.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(CH3CHO      ) = Chemical("CH3CHO      ",  44.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(MACR        ) = Chemical("MACR        ",  70.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(MEK         ) = Chemical("MEK         ",  72.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(ACETOL      ) = Chemical("ACETOL      ",  74.0000,  0,   3.0000,   0.0000,   0.0000 )
    species(GLYOX       ) = Chemical("GLYOX       ",  58.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(MGLYOX      ) = Chemical("MGLYOX      ",  72.0000,  0,   3.0000,   0.0000,   0.0000 )
    species(BIACET      ) = Chemical("BIACET      ",  86.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(C5DICARB    ) = Chemical("C5DICARB    ",  98.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(CH3OOH      ) = Chemical("CH3OOH      ",  48.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(C2H5OOH     ) = Chemical("C2H5OOH     ",  62.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(BURO2H      ) = Chemical("BURO2H      ",  90.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(ETRO2H      ) = Chemical("ETRO2H      ",  78.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(MEKO2H      ) = Chemical("MEKO2H      ", 104.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(ISRO2H      ) = Chemical("ISRO2H      ", 118.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(C5DICAROOH  ) = Chemical("C5DICAROOH  ", 147.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(HPALD       ) = Chemical("HPALD       ", 116.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(MACROOH     ) = Chemical("MACROOH     ", 120.0000,  0,   4.0000,   0.0000,   0.0000 )
    species(OXYO2H      ) = Chemical("OXYO2H      ", 188.0000,  0,   8.0000,   0.0000,   0.0000 )
    species(CH3CO3H     ) = Chemical("CH3CO3H     ",  76.0000,  0,   2.0000,   0.0000,   0.0000 )
    species(PACALD      ) = Chemical("PACALD      ", 130.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(IEPOX       ) = Chemical("IEPOX       ", 118.0000,  0,   5.0000,   0.0000,   0.0000 )
    species(SC4H9NO3    ) = Chemical("SC4H9NO3    ", 119.0000,  0,   4.0000,   1.0000,   0.0000 )
    species(NALD        ) = Chemical("NALD        ", 105.0000,  0,   2.0000,   1.0000,   0.0000 )
    species(ISON        ) = Chemical("ISON        ", 147.0000,  0,   5.0000,   1.0000,   0.0000 )
    species(PAN         ) = Chemical("PAN         ", 121.0000,  0,   2.0000,   1.0000,   0.0000 )
    species(MPAN        ) = Chemical("MPAN        ", 132.0000,  0,   4.0000,   1.0000,   0.0000 )
    species(APINENE     ) = Chemical("APINENE     ", 136.0000,  1,  10.0000,   0.0000,   0.0000 )
    species(TERPOOH     ) = Chemical("TERPOOH     ", 186.0000,  0,  10.0000,   0.0000,   0.0000 )
    species(SO2         ) = Chemical("SO2         ",  64.0000,  0,   0.0000,   0.0000,   1.0000 )
    species(shipNOx     ) = Chemical("shipNOx     ",  46.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(Dust_road_f ) = Chemical("Dust_road_f ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_road_c ) = Chemical("Dust_road_c ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_wb_f   ) = Chemical("Dust_wb_f   ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_wb_c   ) = Chemical("Dust_wb_c   ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_sah_f  ) = Chemical("Dust_sah_f  ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Dust_sah_c  ) = Chemical("Dust_sah_c  ", 200.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SQT_SOA_NV  ) = Chemical("SQT_SOA_NV  ", 302.0000,  0,  14.0000,   0.0000,   0.0000 )
    species(POLLEN_BIRCH) = Chemical("POLLEN_BIRCH",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_OLIVE) = Chemical("POLLEN_OLIVE",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_ALDER) = Chemical("POLLEN_ALDER",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_RWEED) = Chemical("POLLEN_RWEED",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_GRASS) = Chemical("POLLEN_GRASS",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_MUGWORT1) = Chemical("POLLEN_MUGWORT1",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_MUGWORT2) = Chemical("POLLEN_MUGWORT2",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_MUGWORT3) = Chemical("POLLEN_MUGWORT3",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_MUGWORT4) = Chemical("POLLEN_MUGWORT4",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POLLEN_MUGWORT5) = Chemical("POLLEN_MUGWORT5",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(ASOC_ng1e2  ) = Chemical("ASOC_ng1e2  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ASOC_ug1    ) = Chemical("ASOC_ug1    ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ASOC_ug10   ) = Chemical("ASOC_ug10   ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ASOC_ug1e2  ) = Chemical("ASOC_ug1e2  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ASOC_ug1e3  ) = Chemical("ASOC_ug1e3  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(non_C_ASOA_ng1e2) = Chemical("non_C_ASOA_ng1e2",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(non_C_ASOA_ug1) = Chemical("non_C_ASOA_ug1",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(non_C_ASOA_ug10) = Chemical("non_C_ASOA_ug10",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(non_C_ASOA_ug1e2) = Chemical("non_C_ASOA_ug1e2",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(non_C_ASOA_ug1e3) = Chemical("non_C_ASOA_ug1e3",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(BSOC_ng1e2  ) = Chemical("BSOC_ng1e2  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(BSOC_ug1    ) = Chemical("BSOC_ug1    ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(BSOC_ug10   ) = Chemical("BSOC_ug10   ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(BSOC_ug1e2  ) = Chemical("BSOC_ug1e2  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(BSOC_ug1e3  ) = Chemical("BSOC_ug1e3  ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(non_C_BSOA_ng1e2) = Chemical("non_C_BSOA_ng1e2",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(non_C_BSOA_ug1) = Chemical("non_C_BSOA_ug1",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(non_C_BSOA_ug10) = Chemical("non_C_BSOA_ug10",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(non_C_BSOA_ug1e2) = Chemical("non_C_BSOA_ug1e2",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(non_C_BSOA_ug1e3) = Chemical("non_C_BSOA_ug1e3",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SO4         ) = Chemical("SO4         ",  96.0000,  0,   0.0000,   0.0000,   1.0000 )
    species(NH3         ) = Chemical("NH3         ",  17.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NO3_f       ) = Chemical("NO3_f       ",  62.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NO3_c       ) = Chemical("NO3_c       ",  62.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(NH4_f       ) = Chemical("NH4_f       ",  18.0000,  0,   0.0000,   1.0000,   0.0000 )
    species(OM25_bgnd   ) = Chemical("OM25_bgnd   ",  24.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(OM25_p      ) = Chemical("OM25_p      ",   1.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(ffire_OM    ) = Chemical("ffire_OM    ",  20.4000,  0,   1.0000,   0.0000,   0.0000 )
    species(ffire_BC    ) = Chemical("ffire_BC    ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(ffire_remPPM25) = Chemical("ffire_remPPM25",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(ffire_c     ) = Chemical("ffire_c     ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SeaSalt_f   ) = Chemical("SeaSalt_f   ",  58.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(SeaSalt_c   ) = Chemical("SeaSalt_c   ",  58.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Ash_f       ) = Chemical("Ash_f       ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(Ash_c       ) = Chemical("Ash_c       ",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(POM_f_Res   ) = Chemical("POM_f_Res   ",  20.4000,  0,   1.0000,   0.0000,   0.0000 )
    species(POM_c_Res   ) = Chemical("POM_c_Res   ",  20.4000,  0,   1.0000,   0.0000,   0.0000 )
    species(POM_f_nonRes) = Chemical("POM_f_nonRes",  15.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(POM_c_nonRes) = Chemical("POM_c_nonRes",  15.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_f_Res_new) = Chemical("EC_f_Res_new",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_f_Res_age) = Chemical("EC_f_Res_age",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_c_Res    ) = Chemical("EC_c_Res    ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_f_nonRes_new) = Chemical("EC_f_nonRes_new",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_f_nonRes_age) = Chemical("EC_f_nonRes_age",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(EC_c_nonRes ) = Chemical("EC_c_nonRes ",  12.0000,  0,   1.0000,   0.0000,   0.0000 )
    species(remPPM25_nonRes) = Chemical("remPPM25_nonRes",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(remPPM25_Res) = Chemical("remPPM25_Res",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(remPPM_c_nonRes) = Chemical("remPPM_c_nonRes",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(remPPM_c_Res) = Chemical("remPPM_c_Res",  12.0000,  0,   0.0000,   0.0000,   0.0000 )
    species(FUNGAL_SPORES) = Chemical("FUNGAL_SPORES",   0.0000,  0,   0.0000,   0.0000,   0.0000 )
  end subroutine define_chemicals

end module ChemSpecs_mod
