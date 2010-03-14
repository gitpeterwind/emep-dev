!>_________________________________________________________<

  module  ChemSpecs_adv_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 63 
 


   integer, public, parameter ::   & 
     IXADV_O3          =   1   &
  ,  IXADV_NO          =   2   &
  ,  IXADV_NO2         =   3   &
  ,  IXADV_PAN         =   4   &
  ,  IXADV_MPAN        =   5   &
  ,  IXADV_NO3         =   6   &
  ,  IXADV_N2O5        =   7   &
  ,  IXADV_ISONO3      =   8   &
  ,  IXADV_HNO3        =   9

   integer, public, parameter ::   & 
     IXADV_HONO        =  10   &
  ,  IXADV_CH3COO2     =  11   &
  ,  IXADV_MACR        =  12   &
  ,  IXADV_ISNI        =  13   &
  ,  IXADV_ISNIR       =  14   &
  ,  IXADV_GLYOX       =  15   &
  ,  IXADV_MGLYOX      =  16   &
  ,  IXADV_MAL         =  17   &
  ,  IXADV_MEK         =  18   &
  ,  IXADV_MVK         =  19

   integer, public, parameter ::   & 
     IXADV_HCHO        =  20   &
  ,  IXADV_CH3CHO      =  21   &
  ,  IXADV_C2H6        =  22   &
  ,  IXADV_NC4H10      =  23   &
  ,  IXADV_C2H4        =  24   &
  ,  IXADV_C3H6        =  25   &
  ,  IXADV_OXYL        =  26   &
  ,  IXADV_C5H8        =  27   &
  ,  IXADV_CH3O2H      =  28   &
  ,  IXADV_C2H5OOH     =  29

   integer, public, parameter ::   & 
     IXADV_BURO2H      =  30   &
  ,  IXADV_ETRO2H      =  31   &
  ,  IXADV_PRRO2H      =  32   &
  ,  IXADV_OXYO2H      =  33   &
  ,  IXADV_MEKO2H      =  34   &
  ,  IXADV_MALO2H      =  35   &
  ,  IXADV_MVKO2H      =  36   &
  ,  IXADV_MACROOH     =  37   &
  ,  IXADV_MACO3H      =  38   &
  ,  IXADV_MACO2H      =  39

   integer, public, parameter ::   & 
     IXADV_ISRO2H      =  40   &
  ,  IXADV_H2O2        =  41   &
  ,  IXADV_CH3COO2H    =  42   &
  ,  IXADV_ISONO3H     =  43   &
  ,  IXADV_ISNIRH      =  44   &
  ,  IXADV_CH3OH       =  45   &
  ,  IXADV_C2H5OH      =  46   &
  ,  IXADV_ACETOL      =  47   &
  ,  IXADV_H2          =  48   &
  ,  IXADV_CO          =  49

   integer, public, parameter ::   & 
     IXADV_CH4         =  50   &
  ,  IXADV_SO2         =  51   &
  ,  IXADV_SO4         =  52   &
  ,  IXADV_PNO3        =  53   &
  ,  IXADV_NH3         =  54   &
  ,  IXADV_ANH4        =  55   &
  ,  IXADV_ANO3        =  56   &
  ,  IXADV_PPM25       =  57   &
  ,  IXADV_PPM25_FIRE  =  58   &
  ,  IXADV_PPMCO       =  59

   integer, public, parameter ::   & 
     IXADV_SSFI        =  60   &
  ,  IXADV_SSCO        =  61   &
  ,  IXADV_RN222       =  62   &
  ,  IXADV_PB210       =  63

 !-----------------------------------------------------------
  end module ChemSpecs_adv_ml
!>_________________________________________________________<

  module  ChemSpecs_shl_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 16 
 


   integer, public, parameter ::   & 
     IXSHL_OD          =   1   &
  ,  IXSHL_OP          =   2   &
  ,  IXSHL_OH          =   3   &
  ,  IXSHL_HO2         =   4   &
  ,  IXSHL_CH3O2       =   5   &
  ,  IXSHL_C2H5O2      =   6   &
  ,  IXSHL_SECC4H9O2   =   7   &
  ,  IXSHL_ISRO2       =   8   &
  ,  IXSHL_ETRO2       =   9

   integer, public, parameter ::   & 
     IXSHL_PRRO2       =  10   &
  ,  IXSHL_OXYO2       =  11   &
  ,  IXSHL_MEKO2       =  12   &
  ,  IXSHL_MALO2       =  13   &
  ,  IXSHL_MVKO2       =  14   &
  ,  IXSHL_MACRO2      =  15   &
  ,  IXSHL_MACO3       =  16

 !-----------------------------------------------------------
  end module ChemSpecs_shl_ml
!>_________________________________________________________<

  module  ChemSpecs_tot_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_TOT = 79 
 
  ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL=0,   &!   Number of aerosol species
                FIRST_SOA=-999, &!   First aerosol species
                LAST_SOA=-999     !   Last  aerosol species  



   integer, public, parameter ::   & 
     OD          =   1   &
  ,  OP          =   2   &
  ,  OH          =   3   &
  ,  HO2         =   4   &
  ,  CH3O2       =   5   &
  ,  C2H5O2      =   6   &
  ,  SECC4H9O2   =   7   &
  ,  ISRO2       =   8   &
  ,  ETRO2       =   9

   integer, public, parameter ::   & 
     PRRO2       =  10   &
  ,  OXYO2       =  11   &
  ,  MEKO2       =  12   &
  ,  MALO2       =  13   &
  ,  MVKO2       =  14   &
  ,  MACRO2      =  15   &
  ,  MACO3       =  16   &
  ,  O3          =  17   &
  ,  NO          =  18   &
  ,  NO2         =  19

   integer, public, parameter ::   & 
     PAN         =  20   &
  ,  MPAN        =  21   &
  ,  NO3         =  22   &
  ,  N2O5        =  23   &
  ,  ISONO3      =  24   &
  ,  HNO3        =  25   &
  ,  HONO        =  26   &
  ,  CH3COO2     =  27   &
  ,  MACR        =  28   &
  ,  ISNI        =  29

   integer, public, parameter ::   & 
     ISNIR       =  30   &
  ,  GLYOX       =  31   &
  ,  MGLYOX      =  32   &
  ,  MAL         =  33   &
  ,  MEK         =  34   &
  ,  MVK         =  35   &
  ,  HCHO        =  36   &
  ,  CH3CHO      =  37   &
  ,  C2H6        =  38   &
  ,  NC4H10      =  39

   integer, public, parameter ::   & 
     C2H4        =  40   &
  ,  C3H6        =  41   &
  ,  OXYL        =  42   &
  ,  C5H8        =  43   &
  ,  CH3O2H      =  44   &
  ,  C2H5OOH     =  45   &
  ,  BURO2H      =  46   &
  ,  ETRO2H      =  47   &
  ,  PRRO2H      =  48   &
  ,  OXYO2H      =  49

   integer, public, parameter ::   & 
     MEKO2H      =  50   &
  ,  MALO2H      =  51   &
  ,  MVKO2H      =  52   &
  ,  MACROOH     =  53   &
  ,  MACO3H      =  54   &
  ,  MACO2H      =  55   &
  ,  ISRO2H      =  56   &
  ,  H2O2        =  57   &
  ,  CH3COO2H    =  58   &
  ,  ISONO3H     =  59

   integer, public, parameter ::   & 
     ISNIRH      =  60   &
  ,  CH3OH       =  61   &
  ,  C2H5OH      =  62   &
  ,  ACETOL      =  63   &
  ,  H2          =  64   &
  ,  CO          =  65   &
  ,  CH4         =  66   &
  ,  SO2         =  67   &
  ,  SO4         =  68   &
  ,  PNO3        =  69

   integer, public, parameter ::   & 
     NH3         =  70   &
  ,  ANH4        =  71   &
  ,  ANO3        =  72   &
  ,  PPM25       =  73   &
  ,  PPM25_FIRE  =  74   &
  ,  PPMCO       =  75   &
  ,  SSFI        =  76   &
  ,  SSCO        =  77   &
  ,  RN222       =  78   &
  ,  PB210       =  79

 !-----------------------------------------------------------
  end module ChemSpecs_tot_ml
!>_________________________________________________________<

  module  ChemChemicals_ml
!-----------------------------------------------------------

   use ChemSpecs_tot_ml  ! => NSPEC_TOT, species indices
  implicit none
  private

  !/--   Characteristics of species: 
  !/--   Number, name, molwt, carbon num, nmhc (1) or not(0)
 
  public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc

  type, public :: Chemical 
       character(len=12) :: name
       integer           :: molwt
       integer           :: nmhc      ! nmhc (1) or not(0)
       integer           :: carbons   ! Carbon-number
       real              :: nitrogens ! Nitrogen-number
       integer           :: sulphurs  ! Sulphur-number
  end type Chemical
  type(Chemical), public, dimension(NSPEC_TOT) :: species

  contains
    subroutine define_chemicals()
    !+
    ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
    ! array, using indices from total list of species (advected + short-lived).
    !                                           MW  NM   C    N   S
       species(OD) = Chemical("OD          ",  16,  0,  0,   0,  0 ) 
       species(OP) = Chemical("OP          ",  16,  0,  0,   0,  0 ) 
       species(OH) = Chemical("OH          ",  17,  0,  0,   0,  0 ) 
       species(HO2) = Chemical("HO2         ",  33,  0,  0,   0,  0 ) 
       species(CH3O2) = Chemical("CH3O2       ",  47,  0,  1,   0,  0 ) 
       species(C2H5O2) = Chemical("C2H5O2      ",  61,  0,  2,   0,  0 ) 
       species(SECC4H9O2) = Chemical("SECC4H9O2   ",  89,  0,  4,   0,  0 ) 
       species(ISRO2) = Chemical("ISRO2       ", 101,  0,  5,   0,  0 ) 
       species(ETRO2) = Chemical("ETRO2       ",  77,  0,  2,   0,  0 ) 
       species(PRRO2) = Chemical("PRRO2       ",  91,  0,  3,   0,  0 ) 
       species(OXYO2) = Chemical("OXYO2       ",   1,  0,  0,   0,  0 ) 
       species(MEKO2) = Chemical("MEKO2       ", 103,  0,  4,   0,  0 ) 
       species(MALO2) = Chemical("MALO2       ", 147,  0,  5,   0,  0 ) 
       species(MVKO2) = Chemical("MVKO2       ", 119,  0,  4,   0,  0 ) 
       species(MACRO2) = Chemical("MACRO2      ", 119,  0,  4,   0,  0 ) 
       species(MACO3) = Chemical("MACO3       ", 101,  0,  4,   0,  0 ) 
       species(O3) = Chemical("O3          ",  48,  0,  0,   0,  0 ) 
       species(NO) = Chemical("NO          ",  30,  0,  0,   1,  0 ) 
       species(NO2) = Chemical("NO2         ",  46,  0,  0,   1,  0 ) 
       species(PAN) = Chemical("PAN         ", 121,  0,  2,   1,  0 ) 
       species(MPAN) = Chemical("MPAN        ", 132,  0,  4,   1,  0 ) 
       species(NO3) = Chemical("NO3         ",  62,  0,  0,   1,  0 ) 
       species(N2O5) = Chemical("N2O5        ", 108,  0,  0,   2,  0 ) 
       species(ISONO3) = Chemical("ISONO3      ",   1,  0,  0,   0,  0 ) 
       species(HNO3) = Chemical("HNO3        ",  63,  0,  0,   1,  0 ) 
       species(HONO) = Chemical("HONO        ",  47,  0,  0,   1,  0 ) 
       species(CH3COO2) = Chemical("CH3COO2     ",  75,  0,  2,   0,  0 ) 
       species(MACR) = Chemical("MACR        ",  70,  0,  4,   0,  0 ) 
       species(ISNI) = Chemical("ISNI        ",   1,  0,  0,   0,  0 ) 
       species(ISNIR) = Chemical("ISNIR       ",   1,  0,  0,   0,  0 ) 
       species(GLYOX) = Chemical("GLYOX       ",  58,  0,  2,   0,  0 ) 
       species(MGLYOX) = Chemical("MGLYOX      ",  72,  0,  3,   0,  0 ) 
       species(MAL) = Chemical("MAL         ",  98,  0,  5,   0,  0 ) 
       species(MEK) = Chemical("MEK         ",  72,  0,  4,   0,  0 ) 
       species(MVK) = Chemical("MVK         ",  70,  0,  4,   0,  0 ) 
       species(HCHO) = Chemical("HCHO        ",  30,  0,  1,   0,  0 ) 
       species(CH3CHO) = Chemical("CH3CHO      ",  44,  0,  2,   0,  0 ) 
       species(C2H6) = Chemical("C2H6        ",  30,  1,  2,   0,  0 ) 
       species(NC4H10) = Chemical("NC4H10      ",  58,  1,  4,   0,  0 ) 
       species(C2H4) = Chemical("C2H4        ",  28,  1,  2,   0,  0 ) 
       species(C3H6) = Chemical("C3H6        ",  42,  1,  3,   0,  0 ) 
       species(OXYL) = Chemical("OXYL        ", 106,  1,  8,   0,  0 ) 
       species(C5H8) = Chemical("C5H8        ",  68,  1,  5,   0,  0 ) 
       species(CH3O2H) = Chemical("CH3O2H      ",  48,  0,  1,   0,  0 ) 
       species(C2H5OOH) = Chemical("C2H5OOH     ",  62,  0,  2,   0,  0 ) 
       species(BURO2H) = Chemical("BURO2H      ",  90,  0,  4,   0,  0 ) 
       species(ETRO2H) = Chemical("ETRO2H      ",  78,  0,  2,   0,  0 ) 
       species(PRRO2H) = Chemical("PRRO2H      ",  92,  0,  3,   0,  0 ) 
       species(OXYO2H) = Chemical("OXYO2H      ",   1,  0,  0,   0,  0 ) 
       species(MEKO2H) = Chemical("MEKO2H      ", 104,  0,  4,   0,  0 ) 
       species(MALO2H) = Chemical("MALO2H      ", 147,  0,  5,   0,  0 ) 
       species(MVKO2H) = Chemical("MVKO2H      ",   1,  0,  0,   0,  0 ) 
       species(MACROOH) = Chemical("MACROOH     ", 120,  0,  4,   0,  0 ) 
       species(MACO3H) = Chemical("MACO3H      ", 102,  0,  4,   0,  0 ) 
       species(MACO2H) = Chemical("MACO2H      ",  86,  0,  4,   0,  0 ) 
       species(ISRO2H) = Chemical("ISRO2H      ",   1,  0,  0,   0,  0 ) 
       species(H2O2) = Chemical("H2O2        ",  34,  0,  0,   0,  0 ) 
       species(CH3COO2H) = Chemical("CH3COO2H    ",  76,  0,  2,   0,  0 ) 
       species(ISONO3H) = Chemical("ISONO3H     ",   1,  0,  0,   0,  0 ) 
       species(ISNIRH) = Chemical("ISNIRH      ",   1,  0,  0,   0,  0 ) 
       species(CH3OH) = Chemical("CH3OH       ",  32,  0,  1,   0,  0 ) 
       species(C2H5OH) = Chemical("C2H5OH      ",  46,  0,  2,   0,  0 ) 
       species(ACETOL) = Chemical("ACETOL      ",  74,  0,  3,   0,  0 ) 
       species(H2) = Chemical("H2          ",   2,  0,  0,   0,  0 ) 
       species(CO) = Chemical("CO          ",  28,  0,  1,   0,  0 ) 
       species(CH4) = Chemical("CH4         ",  16,  0,  1,   0,  0 ) 
       species(SO2) = Chemical("SO2         ",  64,  0,  0,   0,  1 ) 
       species(SO4) = Chemical("SO4         ",  96,  0,  0,   0,  1 ) 
       species(PNO3) = Chemical("PNO3        ",  62,  0,  0,   1,  0 ) 
       species(NH3) = Chemical("NH3         ",  17,  0,  0,   1,  0 ) 
       species(ANH4) = Chemical("ANH4        ",  18,  0,  0,   1,  0 ) 
       species(ANO3) = Chemical("ANO3        ",  62,  0,  0,   1,  0 ) 
       species(PPM25) = Chemical("PPM25       ",  12,  0,  0,   0,  0 ) 
       species(PPM25_FIRE) = Chemical("PPM25_FIRE  ",  12,  0,  0,   0,  0 ) 
       species(PPMCO) = Chemical("PPMCO       ",   1,  0,  0,   0,  0 ) 
       species(SSFI) = Chemical("SSFI        ",  58,  0,  0,   0,  0 ) 
       species(SSCO) = Chemical("SSCO        ",  58,  0,  0,   0,  0 ) 
       species(RN222) = Chemical("RN222       ", 222,  0,  0,   0,  0 ) 
       species(PB210) = Chemical("PB210       ", 210,  0,  0,   0,  0 ) 
   end subroutine define_chemicals
 end module ChemChemicals_ml
 !-----------------------------------------------------------
!>_________________________________________________________<

  module  ChemGroups_ml
!-----------------------------------------------------------

  use ChemSpecs_tot_ml  ! => species indices
  implicit none
  private
! Assignment of groups from GenIn.species:

! ------- Gas/particle species ------------------

  integer, public, parameter, dimension(2) :: &
             PMCO_GROUP     = (/ PNO3,PPMCO /)

  integer, public, parameter, dimension(2) :: &
             SOX_GROUP     = (/ SO2,SO4 /)

  integer, public, parameter, dimension(2) :: &
             RDN_GROUP     = (/ NH3,ANH4 /)

  integer, public, parameter, dimension(1) :: &
             BVOC_GROUP     = (/ C5H8 /)

  integer, public, parameter, dimension(2) :: &
             OX_GROUP     = (/ O3,NO2 /)

  integer, public, parameter, dimension(4) :: &
             SIA_GROUP     = (/ SO4,PNO3,ANH4,ANO3 /)

  integer, public, parameter, dimension(6) :: &
             PM25_GROUP     = (/ SO4,ANH4,ANO3,PPM25,PPM25_FIRE,SSFI /)

  integer, public, parameter, dimension(2) :: &
             NOX_GROUP     = (/ NO,NO2 /)

  integer, public, parameter, dimension(2) :: &
             SS_GROUP     = (/ PPMCO,SSFI /)

  integer, public, parameter, dimension(13) :: &
             OXN_GROUP     = (/ NO,NO2,PAN,MPAN,NO3,N2O5,ISONO3,HNO3,HONO,ISNI,ISNIR,PNO3,ANO3 /)

  integer, public, parameter, dimension(2) :: &
             TNO3_GROUP     = (/ PNO3,ANO3 /)

! ------- Dry dep      species ------------------
  integer, public, parameter, dimension(7) :: &
               DDEP_OXNGROUP = (/ HNO3,HONO,PAN,NO2,ANO3,MPAN,PNO3 /)
  integer, public, parameter, dimension(2) :: &
               DDEP_SOXGROUP = (/ SO2,SO4 /)
  integer, public, parameter, dimension(2) :: &
               DDEP_RDNGROUP = (/ NH3,ANH4 /)
  integer, public, dimension(7) :: DDEP_GROUP

! ------- RO2 Pool     species ------------------
  integer, public, parameter :: SIZE_RO2_POOL      = 1
  integer, public, parameter, dimension(1) :: &
     RO2_POOL      = (/ -99 /)


 end module ChemGroups_ml
 !-----------------------------------------------------------
 module ChemSpecs_bgn_ml
!-----------------------------------------------------------
! PRETTY MUCH FAKED FOR NOW. CAN BE DELETED SOON IN HOPE!
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP
  implicit none
  private

!+ Defines indices and NSPEC for bgn : Background species

 ! Species which can be specified simply for each column, e.g.
 ! as function of local meteorology or zenith angle
 !   o2, m,  and for MADE-like, oh, ch3coo2

   integer, public, parameter ::  NSPEC_BGN = 0 ! No. 3D bgn species
   integer, public, parameter ::  NSPEC_COL = 0 ! total no. prescribed specs

  !/ define xn_2d_bgn here.
   real, public, save, dimension(1,KCHEMTOP:KMAX_MID) :: xn_2d_bgn

!-----------------------------------------------------------
  end module ChemSpecs_bgn_ml
