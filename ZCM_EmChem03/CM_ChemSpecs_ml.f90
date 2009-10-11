!>_________________________________________________________<

  module  ChemSpecs_adv_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 60 
 


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
     IXADV_CH2CCH3     =  10   &
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
  ,  IXADV_ISOP        =  27   &
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
  ,  IXADV_MARO2H      =  37   &
  ,  IXADV_ISRO2H      =  38   &
  ,  IXADV_H2O2        =  39

   integer, public, parameter ::   & 
     IXADV_CH3COO2H    =  40   &
  ,  IXADV_CH2CO2HCH3  =  41   &
  ,  IXADV_ISONO3H     =  42   &
  ,  IXADV_ISNIRH      =  43   &
  ,  IXADV_CH3OH       =  44   &
  ,  IXADV_C2H5OH      =  45   &
  ,  IXADV_H2          =  46   &
  ,  IXADV_CO          =  47   &
  ,  IXADV_CH4         =  48   &
  ,  IXADV_SO2         =  49

   integer, public, parameter ::   & 
     IXADV_SO4         =  50   &
  ,  IXADV_PNO3        =  51   &
  ,  IXADV_NH3         =  52   &
  ,  IXADV_ANH4        =  53   &
  ,  IXADV_ANO3        =  54   &
  ,  IXADV_PPM25       =  55   &
  ,  IXADV_PPMCO       =  56   &
  ,  IXADV_SSFI        =  57   &
  ,  IXADV_SSCO        =  58   &
  ,  IXADV_RN222       =  59

   integer, public, parameter ::   & 
     IXADV_PB210       =  60

 !-----------------------------------------------------------
  end module ChemSpecs_adv_ml
!>_________________________________________________________<

  module  ChemSpecs_shl_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 15 
 


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
  ,  IXSHL_MACRO2      =  15

 !-----------------------------------------------------------
  end module ChemSpecs_shl_ml
!>_________________________________________________________<

  module  ChemSpecs_tot_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_TOT = 75 
 
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
  ,  O3          =  16   &
  ,  NO          =  17   &
  ,  NO2         =  18   &
  ,  PAN         =  19

   integer, public, parameter ::   & 
     MPAN        =  20   &
  ,  NO3         =  21   &
  ,  N2O5        =  22   &
  ,  ISONO3      =  23   &
  ,  HNO3        =  24   &
  ,  CH2CCH3     =  25   &
  ,  CH3COO2     =  26   &
  ,  MACR        =  27   &
  ,  ISNI        =  28   &
  ,  ISNIR       =  29

   integer, public, parameter ::   & 
     GLYOX       =  30   &
  ,  MGLYOX      =  31   &
  ,  MAL         =  32   &
  ,  MEK         =  33   &
  ,  MVK         =  34   &
  ,  HCHO        =  35   &
  ,  CH3CHO      =  36   &
  ,  C2H6        =  37   &
  ,  NC4H10      =  38   &
  ,  C2H4        =  39

   integer, public, parameter ::   & 
     C3H6        =  40   &
  ,  OXYL        =  41   &
  ,  ISOP        =  42   &
  ,  CH3O2H      =  43   &
  ,  C2H5OOH     =  44   &
  ,  BURO2H      =  45   &
  ,  ETRO2H      =  46   &
  ,  PRRO2H      =  47   &
  ,  OXYO2H      =  48   &
  ,  MEKO2H      =  49

   integer, public, parameter ::   & 
     MALO2H      =  50   &
  ,  MVKO2H      =  51   &
  ,  MARO2H      =  52   &
  ,  ISRO2H      =  53   &
  ,  H2O2        =  54   &
  ,  CH3COO2H    =  55   &
  ,  CH2CO2HCH3  =  56   &
  ,  ISONO3H     =  57   &
  ,  ISNIRH      =  58   &
  ,  CH3OH       =  59

   integer, public, parameter ::   & 
     C2H5OH      =  60   &
  ,  H2          =  61   &
  ,  CO          =  62   &
  ,  CH4         =  63   &
  ,  SO2         =  64   &
  ,  SO4         =  65   &
  ,  PNO3        =  66   &
  ,  NH3         =  67   &
  ,  ANH4        =  68   &
  ,  ANO3        =  69

   integer, public, parameter ::   & 
     PPM25       =  70   &
  ,  PPMCO       =  71   &
  ,  SSFI        =  72   &
  ,  SSCO        =  73   &
  ,  RN222       =  74   &
  ,  PB210       =  75

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
       species(MACRO2) = Chemical("MACRO2      ", 102,  0,  4,   0,  0 ) 
       species(O3) = Chemical("O3          ",  48,  0,  0,   0,  0 ) 
       species(NO) = Chemical("NO          ",  30,  0,  0,   1,  0 ) 
       species(NO2) = Chemical("NO2         ",  46,  0,  0,   1,  0 ) 
       species(PAN) = Chemical("PAN         ", 121,  0,  2,   1,  0 ) 
       species(MPAN) = Chemical("MPAN        ", 132,  0,  4,   1,  0 ) 
       species(NO3) = Chemical("NO3         ",  62,  0,  0,   1,  0 ) 
       species(N2O5) = Chemical("N2O5        ", 108,  0,  0,   2,  0 ) 
       species(ISONO3) = Chemical("ISONO3      ",   1,  0,  0,   0,  0 ) 
       species(HNO3) = Chemical("HNO3        ",  63,  0,  0,   1,  0 ) 
       species(CH2CCH3) = Chemical("CH2CCH3     ",  73,  0,  3,   0,  0 ) 
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
       species(ISOP) = Chemical("ISOP        ",  68,  1,  5,   0,  0 ) 
       species(CH3O2H) = Chemical("CH3O2H      ",  48,  0,  1,   0,  0 ) 
       species(C2H5OOH) = Chemical("C2H5OOH     ",  62,  0,  2,   0,  0 ) 
       species(BURO2H) = Chemical("BURO2H      ",  90,  0,  4,   0,  0 ) 
       species(ETRO2H) = Chemical("ETRO2H      ",  78,  0,  2,   0,  0 ) 
       species(PRRO2H) = Chemical("PRRO2H      ",  92,  0,  3,   0,  0 ) 
       species(OXYO2H) = Chemical("OXYO2H      ",   1,  0,  0,   0,  0 ) 
       species(MEKO2H) = Chemical("MEKO2H      ", 104,  0,  4,   0,  0 ) 
       species(MALO2H) = Chemical("MALO2H      ", 147,  0,  5,   0,  0 ) 
       species(MVKO2H) = Chemical("MVKO2H      ",   1,  0,  0,   0,  0 ) 
       species(MARO2H) = Chemical("MARO2H      ",   1,  0,  0,   0,  0 ) 
       species(ISRO2H) = Chemical("ISRO2H      ",   1,  0,  0,   0,  0 ) 
       species(H2O2) = Chemical("H2O2        ",  34,  0,  0,   0,  0 ) 
       species(CH3COO2H) = Chemical("CH3COO2H    ",  76,  0,  2,   0,  0 ) 
       species(CH2CO2HCH3) = Chemical("CH2CO2HCH3  ",  74,  0,  3,   0,  0 ) 
       species(ISONO3H) = Chemical("ISONO3H     ",   1,  0,  0,   0,  0 ) 
       species(ISNIRH) = Chemical("ISNIRH      ",   1,  0,  0,   0,  0 ) 
       species(CH3OH) = Chemical("CH3OH       ",  32,  0,  1,   0,  0 ) 
       species(C2H5OH) = Chemical("C2H5OH      ",  46,  0,  2,   0,  0 ) 
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
! Automate assignment of OXN species:

! ------- Gas/particle species ------------------
  integer, public, parameter :: SIZE_OXNGROUP = 9
  integer, public, parameter, dimension(9) :: &
     OXNGROUP = (/ NO,NO2,PAN,MPAN,NO3,N2O5,HNO3,PNO3,ANO3 /)

! ------- Dry dep      species ------------------
  integer, public, parameter :: SIZE_DDEP_OXNGROUP = 6
  integer, public, parameter, dimension(6) :: &
     DDEP_OXNGROUP = (/ HNO3,PAN,NO2,ANO3,MPAN,PNO3 /)

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
