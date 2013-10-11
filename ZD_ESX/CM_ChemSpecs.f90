!>_________________________________________________________<

  module  ChemSpecs
!-----------------------------------------------------------

  
  implicit none
!ESXOUT TESTING 
!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_TOT = 20 
 
  ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL=0,   &!   Number of aerosol species
                FIRST_SEMIVOL=-999, &!   First aerosol species
                LAST_SEMIVOL=-999     !   Last  aerosol species  



   integer, public, parameter ::   & 
     OD          =   1   &
  ,  OP          =   2   &
  ,  OH          =   3   &
  ,  HO2         =   4   &
  ,  RO2         =   5   &
  ,  TRACER      =   6   &
  ,  O3          =   7   &
  ,  NO          =   8   &
  ,  NO2         =   9

   integer, public, parameter ::   & 
     HNO3        =  10   &
  ,  VOC         =  11   &
  ,  CH3CHO      =  12   &
  ,  PAN         =  13   &
  ,  CH3COO2     =  14   &
  ,  SO2         =  15   &
  ,  SO4         =  16   &
  ,  NH3         =  17   &
  ,  NO3_F       =  18   &
  ,  NH4_F       =  19

   integer, public, parameter ::   & 
     CO          =  20

 !-----------------------------------------------------------
 
!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 15 
 


   integer, public, parameter ::   & 
     IXADV_TRACER      =   1   &
  ,  IXADV_O3          =   2   &
  ,  IXADV_NO          =   3   &
  ,  IXADV_NO2         =   4   &
  ,  IXADV_HNO3        =   5   &
  ,  IXADV_VOC         =   6   &
  ,  IXADV_CH3CHO      =   7   &
  ,  IXADV_PAN         =   8   &
  ,  IXADV_CH3COO2     =   9

   integer, public, parameter ::   & 
     IXADV_SO2         =  10   &
  ,  IXADV_SO4         =  11   &
  ,  IXADV_NH3         =  12   &
  ,  IXADV_NO3_F       =  13   &
  ,  IXADV_NH4_F       =  14   &
  ,  IXADV_CO          =  15

 !-----------------------------------------------------------
 
!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 5 
 


   integer, public, parameter ::   & 
     IXSHL_OD          =   1   &
  ,  IXSHL_OP          =   2   &
  ,  IXSHL_OH          =   3   &
  ,  IXSHL_HO2         =   4   &
  ,  IXSHL_RO2         =   5

 !-----------------------------------------------------------
 
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
       real              :: ExtC      ! Extinction coef (aerosols)
       real              :: CiStar     ! VBS param
       real              :: DeltaH    ! VBS param
  end type Chemical
  type(Chemical), public, dimension(NSPEC_TOT), target :: species
  type(Chemical), public, dimension(:), pointer :: &
    species_shl=>null(),&             ! => species(..short lived..)
    species_adv=>null()               ! => species(..advected..)

  contains
  subroutine define_chemicals()
  !+
  ! Pointers to short lived and advected portions of species
  !
    species_shl=>species(1:NSPEC_SHL)
    species_adv=>species(NSPEC_SHL+1:NSPEC_SHL+NSPEC_ADV)
  !+
  ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
  ! array, using indices from total list of species (advected + short-lived).
  !                                           MW  NM   C    N   S  ExtC C*  dH
    species(OD          ) = Chemical("OD          ",  16.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(OP          ) = Chemical("OP          ",  16.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(OH          ) = Chemical("OH          ",  17.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(HO2         ) = Chemical("HO2         ",  33.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(RO2         ) = Chemical("RO2         ",  47.0000,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(TRACER      ) = Chemical("TRACER      ",  14.0000,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(O3          ) = Chemical("O3          ",  48.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(NO          ) = Chemical("NO          ",  30.0000,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(NO2         ) = Chemical("NO2         ",  46.0000,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(HNO3        ) = Chemical("HNO3        ",  63.0000,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(VOC         ) = Chemical("VOC         ",  58.0000,  1,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(CH3CHO      ) = Chemical("CH3CHO      ",  44.0000,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(PAN         ) = Chemical("PAN         ", 121.0000,  0,  2,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(CH3COO2     ) = Chemical("CH3COO2     ",  75.0000,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(SO2         ) = Chemical("SO2         ",  64.0000,  0,  0,   0,  1,  0.0,  0.0000,    0.0 ) 
    species(SO4         ) = Chemical("SO4         ",  96.0000,  0,  0,   0,  1,  8.5,  0.0000,    0.0 ) 
    species(NH3         ) = Chemical("NH3         ",  17.0000,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(NO3_F       ) = Chemical("NO3_F       ",  62.0000,  0,  0,   1,  0,  8.5,  0.0000,    0.0 ) 
    species(NH4_F       ) = Chemical("NH4_F       ",  18.0000,  0,  0,   1,  0,  8.5,  0.0000,    0.0 ) 
    species(CO          ) = Chemical("CO          ",  28.0000,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
  end subroutine define_chemicals
end module ChemSpecs
 !-----------------------------------------------------------
