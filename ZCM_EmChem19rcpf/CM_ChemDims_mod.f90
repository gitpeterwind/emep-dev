! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem19c PM_VBS_EmChem19 Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt DustExtended Ash PM_ResNonResInert BVOC_SQT_NV BVOC_IsoMT1_emis Pollen PBAP
module ChemDims_mod

  implicit none
  character(len=*),parameter, public :: CM_schemes_ChemDims = " EmChem19c PM_VBS_EmChem19 Aqueous_EmChem16x Aero2017nx ShipNOx PM_FFireInert SeaSalt DustExtended Ash PM_ResNonResInert BVOC_SQT_NV BVOC_IsoMT1_emis Pollen PBAP"
  

    ! NSPEC for TOT : All reacting species
    integer, public, parameter :: NSPEC_TOT=141
    
    ! NSPEC for ADV : Advected species
    integer, public, parameter :: NSPEC_ADV=126
    
    ! NSPEC for SHL : Short-lived (non-advected) species
    integer, public, parameter :: NSPEC_SHL=15
    
    ! NSPEC for SEMIVOL : Semi-volatile organic aerosols
    integer, public, parameter :: NSPEC_SEMIVOL=20
        
    ! No. DRY deposition species
    integer, public, parameter :: NDRYDEP_ADV = 106
    
    ! No. WET deposition species
    integer, public, parameter :: NWETDEP_ADV = 93
    
    ! No. rate coefficients
    integer, parameter, public :: NCHEMRATES = 101
    
    ! No. photolysis rates used
    integer, parameter, public :: NPHOTOLRATES = 18
    
    ! No. emission Files
    integer, parameter, public :: NEMIS_File = 7
    
    ! No. emission Specss
    integer, parameter, public :: NEMIS_Specs = 64

end module ChemDims_mod
