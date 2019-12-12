! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem19a PM_VBS_EmChem19 BVOC_MTERP1_EmChem19 Aqueous_EmChem16x Aero2017nx ShipNOx FFireInert SeaSalt DustExtended Ash BVOC_SQT_NV
module ChemDims_mod

  implicit none
  character(len=*),parameter, public :: CM_schemes_ChemDims = " EmChem19a PM_VBS_EmChem19 BVOC_MTERP1_EmChem19 Aqueous_EmChem16x Aero2017nx ShipNOx FFireInert SeaSalt DustExtended Ash BVOC_SQT_NV"
  

    ! NSPEC for TOT : All reacting species
    integer, public, parameter :: NSPEC_TOT=131
    
    ! NSPEC for ADV : Advected species
    integer, public, parameter :: NSPEC_ADV=116
    
    ! NSPEC for SHL : Short-lived (non-advected) species
    integer, public, parameter :: NSPEC_SHL=15
    
    ! NSPEC for SEMIVOL : Semi-volatile organic aerosols
    integer, public, parameter :: NSPEC_SEMIVOL=20
        
    ! No. DRY deposition species
    integer, public, parameter :: NDRYDEP_ADV = 93
    
    ! No. WET deposition species
    integer, public, parameter :: NWETDEP_ADV = 80
    
    ! No. rate coefficients
    integer, parameter, public :: NCHEMRATES = 101
    
    ! No. photolysis rates used
    integer, parameter, public :: NPHOTOLRATES = 16
    
    ! No. emission Files
    integer, parameter, public :: NEMIS_File = 7
    
    ! No. emission Specss
    integer, parameter, public :: NEMIS_Specs = 47

end module ChemDims_mod
