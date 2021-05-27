! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem19a PM_VBS_EmChem19 Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt DustExtended Ash PM_ResNonResInert BVOC_SQT_NV BVOC_IsoMT1_emis
module ChemSpecs_mod

  use ChemDims_mod      ! => NSPEC_TOT, NCHEMRATES, ....
  
  implicit none
  private
  character(len=*),parameter, public :: CM_schemes_ChemSpecs = " EmChem19a PM_VBS_EmChem19 Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt DustExtended Ash PM_ResNonResInert BVOC_SQT_NV BVOC_IsoMT1_emis"
  
  
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
    , NMVOC       =  18  &
    , O3          =  19  &
    , NO          =  20
  
  integer, public, parameter :: &
      NO2         =  21  &
    , NO3         =  22  &
    , N2O5        =  23  &
    , H2          =  24  &
    , H2O2        =  25  &
    , HONO        =  26  &
    , HNO3        =  27  &
    , HO2NO2      =  28  &
    , CO          =  29  &
    , CH4         =  30
  
  integer, public, parameter :: &
      C2H6        =  31  &
    , NC4H10      =  32  &
    , C2H4        =  33  &
    , C3H6        =  34  &
    , BENZENE     =  35  &
    , TOLUENE     =  36  &
    , OXYL        =  37  &
    , C5H8        =  38  &
    , CH3OH       =  39  &
    , C2H5OH      =  40
  
  integer, public, parameter :: &
      HCHO        =  41  &
    , CH3CHO      =  42  &
    , MACR        =  43  &
    , MEK         =  44  &
    , ACETOL      =  45  &
    , GLYOX       =  46  &
    , MGLYOX      =  47  &
    , BIACET      =  48  &
    , C5DICARB    =  49  &
    , CH3OOH      =  50
  
  integer, public, parameter :: &
      C2H5OOH     =  51  &
    , BURO2H      =  52  &
    , ETRO2H      =  53  &
    , MEKO2H      =  54  &
    , ISRO2H      =  55  &
    , C5DICAROOH  =  56  &
    , HPALD       =  57  &
    , MACROOH     =  58  &
    , OXYO2H      =  59  &
    , CH3CO3H     =  60
  
  integer, public, parameter :: &
      PACALD      =  61  &
    , IEPOX       =  62  &
    , SC4H9NO3    =  63  &
    , NALD        =  64  &
    , ISON        =  65  &
    , PAN         =  66  &
    , MPAN        =  67  &
    , APINENE     =  68  &
    , TERPOOH     =  69  &
    , SO2         =  70
  
  integer, public, parameter :: &
      shipNOx     =  71  &
    , Dust_road_f =  72  &
    , Dust_road_c =  73  &
    , Dust_wb_f   =  74  &
    , Dust_wb_c   =  75  &
    , Dust_sah_f  =  76  &
    , Dust_sah_c  =  77  &
    , SQT_SOA_NV  =  78  &
    , ASOC_ng1e2  =  79  &
    , ASOC_ug1    =  80
  
  integer, public, parameter :: &
      ASOC_ug10   =  81  &
    , ASOC_ug1e2  =  82  &
    , ASOC_ug1e3  =  83  &
    , non_C_ASOA_ng1e2=  84  &
    , non_C_ASOA_ug1=  85  &
    , non_C_ASOA_ug10=  86  &
    , non_C_ASOA_ug1e2=  87  &
    , non_C_ASOA_ug1e3=  88  &
    , BSOC_ng1e2  =  89  &
    , BSOC_ug1    =  90
  
  integer, public, parameter :: &
      BSOC_ug10   =  91  &
    , BSOC_ug1e2  =  92  &
    , BSOC_ug1e3  =  93  &
    , non_C_BSOA_ng1e2=  94  &
    , non_C_BSOA_ug1=  95  &
    , non_C_BSOA_ug10=  96  &
    , non_C_BSOA_ug1e2=  97  &
    , non_C_BSOA_ug1e3=  98  &
    , SO4         =  99  &
    , NH3         = 100
  
  integer, public, parameter :: &
      NO3_f       = 101  &
    , NO3_c       = 102  &
    , NH4_f       = 103  &
    , OM25_bgnd   = 104  &
    , OM25_p      = 105  &
    , ffire_OM    = 106  &
    , ffire_BC    = 107  &
    , ffire_remPPM25= 108  &
    , ffire_c     = 109  &
    , SeaSalt_f   = 110
  
  integer, public, parameter :: &
      SeaSalt_c   = 111  &
    , Ash_f       = 112  &
    , Ash_c       = 113
  
  !+ Defines indices for ADV : Advected species
  integer, public, parameter :: FIRST_ADV=16, &
                                 LAST_ADV=113
  
  integer, public, parameter :: &
      IXADV_RO2POOL     =   1  &
    , IXADV_CH3CO3      =   2  &
    , IXADV_NMVOC       =   3  &
    , IXADV_O3          =   4  &
    , IXADV_NO          =   5  &
    , IXADV_NO2         =   6  &
    , IXADV_NO3         =   7  &
    , IXADV_N2O5        =   8  &
    , IXADV_H2          =   9  &
    , IXADV_H2O2        =  10
  
  integer, public, parameter :: &
      IXADV_HONO        =  11  &
    , IXADV_HNO3        =  12  &
    , IXADV_HO2NO2      =  13  &
    , IXADV_CO          =  14  &
    , IXADV_CH4         =  15  &
    , IXADV_C2H6        =  16  &
    , IXADV_NC4H10      =  17  &
    , IXADV_C2H4        =  18  &
    , IXADV_C3H6        =  19  &
    , IXADV_BENZENE     =  20
  
  integer, public, parameter :: &
      IXADV_TOLUENE     =  21  &
    , IXADV_OXYL        =  22  &
    , IXADV_C5H8        =  23  &
    , IXADV_CH3OH       =  24  &
    , IXADV_C2H5OH      =  25  &
    , IXADV_HCHO        =  26  &
    , IXADV_CH3CHO      =  27  &
    , IXADV_MACR        =  28  &
    , IXADV_MEK         =  29  &
    , IXADV_ACETOL      =  30
  
  integer, public, parameter :: &
      IXADV_GLYOX       =  31  &
    , IXADV_MGLYOX      =  32  &
    , IXADV_BIACET      =  33  &
    , IXADV_C5DICARB    =  34  &
    , IXADV_CH3OOH      =  35  &
    , IXADV_C2H5OOH     =  36  &
    , IXADV_BURO2H      =  37  &
    , IXADV_ETRO2H      =  38  &
    , IXADV_MEKO2H      =  39  &
    , IXADV_ISRO2H      =  40
  
  integer, public, parameter :: &
      IXADV_C5DICAROOH  =  41  &
    , IXADV_HPALD       =  42  &
    , IXADV_MACROOH     =  43  &
    , IXADV_OXYO2H      =  44  &
    , IXADV_CH3CO3H     =  45  &
    , IXADV_PACALD      =  46  &
    , IXADV_IEPOX       =  47  &
    , IXADV_SC4H9NO3    =  48  &
    , IXADV_NALD        =  49  &
    , IXADV_ISON        =  50
  
  integer, public, parameter :: &
      IXADV_PAN         =  51  &
    , IXADV_MPAN        =  52  &
    , IXADV_APINENE     =  53  &
    , IXADV_TERPOOH     =  54  &
    , IXADV_SO2         =  55  &
    , IXADV_shipNOx     =  56  &
    , IXADV_Dust_road_f =  57  &
    , IXADV_Dust_road_c =  58  &
    , IXADV_Dust_wb_f   =  59  &
    , IXADV_Dust_wb_c   =  60
  
  integer, public, parameter :: &
      IXADV_Dust_sah_f  =  61  &
    , IXADV_Dust_sah_c  =  62  &
    , IXADV_SQT_SOA_NV  =  63  &
    , IXADV_ASOC_ng1e2  =  64  &
    , IXADV_ASOC_ug1    =  65  &
    , IXADV_ASOC_ug10   =  66  &
    , IXADV_ASOC_ug1e2  =  67  &
    , IXADV_ASOC_ug1e3  =  68  &
    , IXADV_non_C_ASOA_ng1e2=  69  &
    , IXADV_non_C_ASOA_ug1=  70
  
  integer, public, parameter :: &
      IXADV_non_C_ASOA_ug10=  71  &
    , IXADV_non_C_ASOA_ug1e2=  72  &
    , IXADV_non_C_ASOA_ug1e3=  73  &
    , IXADV_BSOC_ng1e2  =  74  &
    , IXADV_BSOC_ug1    =  75  &
    , IXADV_BSOC_ug10   =  76  &
    , IXADV_BSOC_ug1e2  =  77  &
    , IXADV_BSOC_ug1e3  =  78  &
    , IXADV_non_C_BSOA_ng1e2=  79  &
    , IXADV_non_C_BSOA_ug1=  80
  
  integer, public, parameter :: &
      IXADV_non_C_BSOA_ug10=  81  &
    , IXADV_non_C_BSOA_ug1e2=  82  &
    , IXADV_non_C_BSOA_ug1e3=  83  &
    , IXADV_SO4         =  84  &
    , IXADV_NH3         =  85  &
    , IXADV_NO3_f       =  86  &
    , IXADV_NO3_c       =  87  &
    , IXADV_NH4_f       =  88  &
    , IXADV_OM25_bgnd   =  89  &
    , IXADV_OM25_p      =  90
  
  integer, public, parameter :: &
      IXADV_ffire_OM    =  91  &
    , IXADV_ffire_BC    =  92  &
    , IXADV_ffire_remPPM25=  93  &
    , IXADV_ffire_c     =  94  &
    , IXADV_SeaSalt_f   =  95  &
    , IXADV_SeaSalt_c   =  96  &
    , IXADV_Ash_f       =  97  &
    , IXADV_Ash_c       =  98
  
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
  integer, public, parameter :: FIRST_SEMIVOL=79, &
                                 LAST_SEMIVOL=98
  
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
    species(NMVOC       ) = Chemical("NMVOC       ",  36.7670,  0,   0.0000,   0.0000,   0.0000 )
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
  end subroutine define_chemicals

end module ChemSpecs_mod
