! Generated by GenChem.py - DO NOT EDIT
module ChemSpecs

  implicit none
  
  !+ Defines indices and NSPEC for TOT : All reacting species
  
  integer, public, parameter :: NSPEC_TOT=75
  
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
    , MVKO2       =  14  &
    , MACRO2      =  15  &
    , MACO3       =  16  &
    , TRACER1     =  17  &
    , TRACER2     =  18  &
    , O3          =  19  &
    , NO          =  20
  
  integer, public, parameter :: &
      NO2         =  21  &
    , PAN         =  22  &
    , MPAN        =  23  &
    , NO3         =  24  &
    , N2O5        =  25  &
    , ISONO3      =  26  &
    , HNO3        =  27  &
    , HONO        =  28  &
    , CH3COO2     =  29  &
    , MACR        =  30
  
  integer, public, parameter :: &
      ISNI        =  31  &
    , ISNIR       =  32  &
    , GLYOX       =  33  &
    , MGLYOX      =  34  &
    , MAL         =  35  &
    , MEK         =  36  &
    , MVK         =  37  &
    , HCHO        =  38  &
    , CH3CHO      =  39  &
    , C2H6        =  40
  
  integer, public, parameter :: &
      NC4H10      =  41  &
    , C2H4        =  42  &
    , C3H6        =  43  &
    , OXYL        =  44  &
    , C5H8        =  45  &
    , APINENE     =  46  &
    , CH3O2H      =  47  &
    , C2H5OOH     =  48  &
    , BURO2H      =  49  &
    , ETRO2H      =  50
  
  integer, public, parameter :: &
      PRRO2H      =  51  &
    , OXYO2H      =  52  &
    , MEKO2H      =  53  &
    , MALO2H      =  54  &
    , MVKO2H      =  55  &
    , MACROOH     =  56  &
    , MACO3H      =  57  &
    , MACO2H      =  58  &
    , ISRO2H      =  59  &
    , H2O2        =  60
  
  integer, public, parameter :: &
      CH3COO2H    =  61  &
    , ISONO3H     =  62  &
    , ISNIRH      =  63  &
    , CH3OH       =  64  &
    , C2H5OH      =  65  &
    , ACETOL      =  66  &
    , H2          =  67  &
    , CO          =  68  &
    , CH4         =  69  &
    , SO2         =  70
  
  integer, public, parameter :: &
      SO4         =  71  &
    , NH3         =  72  &
    , NO3_f       =  73  &
    , NO3_c       =  74  &
    , NH4_f       =  75
  
  !+ Defines indices and NSPEC for ADV : Advected species
  
  integer, public, parameter :: NSPEC_ADV=59
  integer, public, parameter :: FIRST_ADV=17, &
                                 LAST_ADV=75
  
  integer, public, parameter :: &
      IXADV_TRACER1     =   1  &
    , IXADV_TRACER2     =   2  &
    , IXADV_O3          =   3  &
    , IXADV_NO          =   4  &
    , IXADV_NO2         =   5  &
    , IXADV_PAN         =   6  &
    , IXADV_MPAN        =   7  &
    , IXADV_NO3         =   8  &
    , IXADV_N2O5        =   9  &
    , IXADV_ISONO3      =  10
  
  integer, public, parameter :: &
      IXADV_HNO3        =  11  &
    , IXADV_HONO        =  12  &
    , IXADV_CH3COO2     =  13  &
    , IXADV_MACR        =  14  &
    , IXADV_ISNI        =  15  &
    , IXADV_ISNIR       =  16  &
    , IXADV_GLYOX       =  17  &
    , IXADV_MGLYOX      =  18  &
    , IXADV_MAL         =  19  &
    , IXADV_MEK         =  20
  
  integer, public, parameter :: &
      IXADV_MVK         =  21  &
    , IXADV_HCHO        =  22  &
    , IXADV_CH3CHO      =  23  &
    , IXADV_C2H6        =  24  &
    , IXADV_NC4H10      =  25  &
    , IXADV_C2H4        =  26  &
    , IXADV_C3H6        =  27  &
    , IXADV_OXYL        =  28  &
    , IXADV_C5H8        =  29  &
    , IXADV_APINENE     =  30
  
  integer, public, parameter :: &
      IXADV_CH3O2H      =  31  &
    , IXADV_C2H5OOH     =  32  &
    , IXADV_BURO2H      =  33  &
    , IXADV_ETRO2H      =  34  &
    , IXADV_PRRO2H      =  35  &
    , IXADV_OXYO2H      =  36  &
    , IXADV_MEKO2H      =  37  &
    , IXADV_MALO2H      =  38  &
    , IXADV_MVKO2H      =  39  &
    , IXADV_MACROOH     =  40
  
  integer, public, parameter :: &
      IXADV_MACO3H      =  41  &
    , IXADV_MACO2H      =  42  &
    , IXADV_ISRO2H      =  43  &
    , IXADV_H2O2        =  44  &
    , IXADV_CH3COO2H    =  45  &
    , IXADV_ISONO3H     =  46  &
    , IXADV_ISNIRH      =  47  &
    , IXADV_CH3OH       =  48  &
    , IXADV_C2H5OH      =  49  &
    , IXADV_ACETOL      =  50
  
  integer, public, parameter :: &
      IXADV_H2          =  51  &
    , IXADV_CO          =  52  &
    , IXADV_CH4         =  53  &
    , IXADV_SO2         =  54  &
    , IXADV_SO4         =  55  &
    , IXADV_NH3         =  56  &
    , IXADV_NO3_f       =  57  &
    , IXADV_NO3_c       =  58  &
    , IXADV_NH4_f       =  59
  
  !+ Defines indices and NSPEC for SHL : Short-lived (non-advected) species
  
  integer, public, parameter :: NSPEC_SHL=16
  integer, public, parameter :: FIRST_SHL=1, &
                                 LAST_SHL=16
  
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
    , IXSHL_MVKO2       =  14  &
    , IXSHL_MACRO2      =  15  &
    , IXSHL_MACO3       =  16
  
  !+ Defines indices and NSPEC for SEMIVOL : Semi-volatile organic aerosols
  
  integer, public, parameter :: NSPEC_SEMIVOL=0
  integer, public, parameter :: FIRST_SEMIVOL=-999, &
                                 LAST_SEMIVOL=-999
  ! Compatibility with existing code
  integer, public, parameter :: NAEROSOL=NSPEC_SEMIVOL
  
  !/--   Characteristics of species:
  !/--   Number, name, molwt, carbon num, nmhc (1) or not(0)
  
  public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc
  
  type, public :: Chemical
       character(len=20) :: name
       real              :: molwt
       integer           :: nmhc      ! nmhc (1) or not(0)
       integer           :: carbons   ! Carbon-number
       real              :: nitrogens ! Nitrogen-number
       integer           :: sulphurs  ! Sulphur-number
       real              :: CiStar     ! VBS param
       real              :: DeltaH    ! VBS param
  endtype Chemical
  type(Chemical), public, dimension(NSPEC_TOT), target :: species
  
  ! Pointers to parts of species (e.g. short-lived, advected)
  type(Chemical), public, dimension(:), pointer :: species_adv=>null()
  type(Chemical), public, dimension(:), pointer :: species_shl=>null()
  type(Chemical), public, dimension(:), pointer :: species_semivol=>null()

contains
  subroutine define_chemicals()
    !+
    ! Pointers to parts of species (e.g. short-lived, advected), only assigned if
    ! non-empty.
    !
    if (NSPEC_ADV > 0) then
      species_adv => species(FIRST_ADV:LAST_ADV)
    end if
    if (NSPEC_SHL > 0) then
      species_shl => species(FIRST_SHL:LAST_SHL)
    end if
    if (NSPEC_SEMIVOL > 0) then
      species_semivol => species(FIRST_SEMIVOL:LAST_SEMIVOL)
    end if
    
    !+
    ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
    ! array, using indices from total list of species (advected + short-lived).
    !                                                        MW  NM   C    N   S       C*      dH
    species(OD          ) = Chemical("OD          ",  16.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(OP          ) = Chemical("OP          ",  16.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(OH          ) = Chemical("OH          ",  17.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(HO2         ) = Chemical("HO2         ",  33.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(CH3O2       ) = Chemical("CH3O2       ",  47.0000,  0,  1,   0,  0,  0.0000,    0.0 )
    species(C2H5O2      ) = Chemical("C2H5O2      ",  61.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(SECC4H9O2   ) = Chemical("SECC4H9O2   ",  89.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(ISRO2       ) = Chemical("ISRO2       ", 101.0000,  0,  5,   0,  0,  0.0000,    0.0 )
    species(ETRO2       ) = Chemical("ETRO2       ",  77.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(PRRO2       ) = Chemical("PRRO2       ",  91.0000,  0,  3,   0,  0,  0.0000,    0.0 )
    species(OXYO2       ) = Chemical("OXYO2       ",   0.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(MEKO2       ) = Chemical("MEKO2       ", 103.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(MALO2       ) = Chemical("MALO2       ", 147.0000,  0,  5,   0,  0,  0.0000,    0.0 )
    species(MVKO2       ) = Chemical("MVKO2       ", 119.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(MACRO2      ) = Chemical("MACRO2      ", 119.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(MACO3       ) = Chemical("MACO3       ", 101.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(TRACER1     ) = Chemical("TRACER1     ",  14.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(TRACER2     ) = Chemical("TRACER2     ",  14.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(O3          ) = Chemical("O3          ",  48.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(NO          ) = Chemical("NO          ",  30.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(NO2         ) = Chemical("NO2         ",  46.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(PAN         ) = Chemical("PAN         ", 121.0000,  0,  2,   1,  0,  0.0000,    0.0 )
    species(MPAN        ) = Chemical("MPAN        ", 132.0000,  0,  4,   1,  0,  0.0000,    0.0 )
    species(NO3         ) = Chemical("NO3         ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(N2O5        ) = Chemical("N2O5        ", 108.0000,  0,  0,   2,  0,  0.0000,    0.0 )
    species(ISONO3      ) = Chemical("ISONO3      ",   0.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(HNO3        ) = Chemical("HNO3        ",  63.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(HONO        ) = Chemical("HONO        ",  47.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(CH3COO2     ) = Chemical("CH3COO2     ",  75.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(MACR        ) = Chemical("MACR        ",  70.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(ISNI        ) = Chemical("ISNI        ",   0.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(ISNIR       ) = Chemical("ISNIR       ",   0.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(GLYOX       ) = Chemical("GLYOX       ",  58.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(MGLYOX      ) = Chemical("MGLYOX      ",  72.0000,  0,  3,   0,  0,  0.0000,    0.0 )
    species(MAL         ) = Chemical("MAL         ",  98.0000,  0,  5,   0,  0,  0.0000,    0.0 )
    species(MEK         ) = Chemical("MEK         ",  72.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(MVK         ) = Chemical("MVK         ",  70.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(HCHO        ) = Chemical("HCHO        ",  30.0000,  0,  1,   0,  0,  0.0000,    0.0 )
    species(CH3CHO      ) = Chemical("CH3CHO      ",  44.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(C2H6        ) = Chemical("C2H6        ",  30.0000,  1,  2,   0,  0,  0.0000,    0.0 )
    species(NC4H10      ) = Chemical("NC4H10      ",  58.0000,  1,  4,   0,  0,  0.0000,    0.0 )
    species(C2H4        ) = Chemical("C2H4        ",  28.0000,  1,  2,   0,  0,  0.0000,    0.0 )
    species(C3H6        ) = Chemical("C3H6        ",  42.0000,  1,  3,   0,  0,  0.0000,    0.0 )
    species(OXYL        ) = Chemical("OXYL        ", 106.0000,  1,  8,   0,  0,  0.0000,    0.0 )
    species(C5H8        ) = Chemical("C5H8        ",  68.0000,  1,  5,   0,  0,  0.0000,    0.0 )
    species(APINENE     ) = Chemical("APINENE     ", 136.0000,  1, 10,   0,  0,  0.0000,    0.0 )
    species(CH3O2H      ) = Chemical("CH3O2H      ",  48.0000,  0,  1,   0,  0,  0.0000,    0.0 )
    species(C2H5OOH     ) = Chemical("C2H5OOH     ",  62.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(BURO2H      ) = Chemical("BURO2H      ",  90.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(ETRO2H      ) = Chemical("ETRO2H      ",  78.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(PRRO2H      ) = Chemical("PRRO2H      ",  92.0000,  0,  3,   0,  0,  0.0000,    0.0 )
    species(OXYO2H      ) = Chemical("OXYO2H      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(MEKO2H      ) = Chemical("MEKO2H      ", 104.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(MALO2H      ) = Chemical("MALO2H      ", 147.0000,  0,  5,   0,  0,  0.0000,    0.0 )
    species(MVKO2H      ) = Chemical("MVKO2H      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(MACROOH     ) = Chemical("MACROOH     ", 120.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(MACO3H      ) = Chemical("MACO3H      ", 102.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(MACO2H      ) = Chemical("MACO2H      ",  86.0000,  0,  4,   0,  0,  0.0000,    0.0 )
    species(ISRO2H      ) = Chemical("ISRO2H      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(H2O2        ) = Chemical("H2O2        ",  34.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(CH3COO2H    ) = Chemical("CH3COO2H    ",  76.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(ISONO3H     ) = Chemical("ISONO3H     ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(ISNIRH      ) = Chemical("ISNIRH      ",   1.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(CH3OH       ) = Chemical("CH3OH       ",  32.0000,  0,  1,   0,  0,  0.0000,    0.0 )
    species(C2H5OH      ) = Chemical("C2H5OH      ",  46.0000,  0,  2,   0,  0,  0.0000,    0.0 )
    species(ACETOL      ) = Chemical("ACETOL      ",  74.0000,  0,  3,   0,  0,  0.0000,    0.0 )
    species(H2          ) = Chemical("H2          ",   2.0000,  0,  0,   0,  0,  0.0000,    0.0 )
    species(CO          ) = Chemical("CO          ",  28.0000,  0,  1,   0,  0,  0.0000,    0.0 )
    species(CH4         ) = Chemical("CH4         ",  16.0000,  0,  1,   0,  0,  0.0000,    0.0 )
    species(SO2         ) = Chemical("SO2         ",  64.0000,  0,  0,   0,  1,  0.0000,    0.0 )
    species(SO4         ) = Chemical("SO4         ",  96.0000,  0,  0,   0,  1,  0.0000,    0.0 )
    species(NH3         ) = Chemical("NH3         ",  17.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(NO3_f       ) = Chemical("NO3_f       ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(NO3_c       ) = Chemical("NO3_c       ",  62.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(NH4_f       ) = Chemical("NH4_f       ",  18.0000,  0,  0,   1,  0,  0.0000,    0.0 )
    species(CO          ) = Chemical("CO          ",  28.0000,  0,  1,   0,  0,  0.0000,    0.0 ) 
  end subroutine define_chemicals

end module ChemSpecs
