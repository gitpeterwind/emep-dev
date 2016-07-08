ESXOUT TESTING 
!>_________________________________________________________<

  module  ChemSpecs_tot_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_TOT = 13 
 
  ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL=0,   &!   Number of aerosol species
                FIRST_SEMIVOL=-999, &!   First aerosol species
                LAST_SEMIVOL=-999     !   Last  aerosol species  



   integer, public, parameter ::   & 
     OP          =   1   &
  ,  OH          =   2   &
  ,  HO2         =   3   &
  ,  RO2         =   4   &
  ,  O3          =   5   &
  ,  NO          =   6   &
  ,  NO2         =   7   &
  ,  HNO3        =   8   &
  ,  VOC         =   9

   integer, public, parameter ::   & 
     CH3CHO      =  10   &
  ,  PAN         =  11   &
  ,  CH3COO2     =  12   &
  ,  CO          =  13

 !-----------------------------------------------------------
  end module ChemSpecs_tot_ml
!>_________________________________________________________<

  module  ChemSpecs_adv_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 9 
 


   integer, public, parameter ::   & 
     IXADV_O3          =   1   &
  ,  IXADV_NO          =   2   &
  ,  IXADV_NO2         =   3   &
  ,  IXADV_HNO3        =   4   &
  ,  IXADV_VOC         =   5   &
  ,  IXADV_CH3CHO      =   6   &
  ,  IXADV_PAN         =   7   &
  ,  IXADV_CH3COO2     =   8   &
  ,  IXADV_CO          =   9

 !-----------------------------------------------------------
  end module ChemSpecs_adv_ml
!>_________________________________________________________<

  module  ChemSpecs_shl_ml
!-----------------------------------------------------------

  
  implicit none
!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 4 
 


   integer, public, parameter ::   & 
     IXSHL_OP          =   1   &
  ,  IXSHL_OH          =   2   &
  ,  IXSHL_HO2         =   3   &
  ,  IXSHL_RO2         =   4

 !-----------------------------------------------------------
  end module ChemSpecs_shl_ml

  module ChemChemicals_ml 
  
  use ChemSpecs_tot_ml  ! => NSPEC_TOT, species indices
  use ChemSpecs_shl_ml, only: NSPEC_SHL
  use ChemSpecs_adv_ml, only: NSPEC_ADV
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
    species(OP          ) = Chemical("OP          ",  16.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(OH          ) = Chemical("OH          ",  17.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(HO2         ) = Chemical("HO2         ",  33.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(RO2         ) = Chemical("RO2         ",  47.0000,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(O3          ) = Chemical("O3          ",  48.0000,  0,  0,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(NO          ) = Chemical("NO          ",  30.0000,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(NO2         ) = Chemical("NO2         ",  46.0000,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(HNO3        ) = Chemical("HNO3        ",  63.0000,  0,  0,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(VOC         ) = Chemical("VOC         ",  58.0000,  1,  4,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(CH3CHO      ) = Chemical("CH3CHO      ",  44.0000,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(PAN         ) = Chemical("PAN         ", 121.0000,  0,  2,   1,  0,  0.0,  0.0000,    0.0 ) 
    species(CH3COO2     ) = Chemical("CH3COO2     ",  75.0000,  0,  2,   0,  0,  0.0,  0.0000,    0.0 ) 
    species(CO          ) = Chemical("CO          ",  28.0000,  0,  1,   0,  0,  0.0,  0.0000,    0.0 ) 
  end subroutine define_chemicals
end module ChemChemicals_ml
 !-----------------------------------------------------------
