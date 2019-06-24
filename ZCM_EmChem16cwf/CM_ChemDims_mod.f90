! Generated by GenChem.py - DO NOT EDIT
! scheme(s)  EmChem16x ShipNOx BVOC_EmChem16x Aqueous_EmChem16x Aero2017nx VBS_EmChem16x FFireInert SeaSalt DustExtended Ash_EmChem16x Pollen
module ChemDims_mod

  implicit none
  character(len=*),parameter, public :: CM_schemes_ChemDims = " EmChem16x ShipNOx BVOC_EmChem16x Aqueous_EmChem16x Aero2017nx VBS_EmChem16x FFireInert SeaSalt DustExtended Ash_EmChem16x Pollen"
  

    ! NSPEC for TOT : All reacting species
    integer, public, parameter :: NSPEC_TOT=133
    
    ! NSPEC for ADV : Advected species
    integer, public, parameter :: NSPEC_ADV=115
    
    ! NSPEC for SHL : Short-lived (non-advected) species
    integer, public, parameter :: NSPEC_SHL=18
    
    ! NSPEC for SEMIVOL : Semi-volatile organic aerosols
    integer, public, parameter :: NSPEC_SEMIVOL=22
        
    ! No. DRY deposition species
    integer, public, parameter :: NDRYDEP_ADV = 90
    
    ! No. WET deposition species
    integer, public, parameter :: NWETDEP_ADV = 79
    
    ! No. rate coefficients
    integer, parameter, public :: NCHEMRATES = 117
    
    ! No. photolysis rates used
    integer, parameter, public :: NPHOTOLRATES = 14
    
    ! No. emission Files
    integer, parameter, public :: NEMIS_File = 7
    
    ! No. emission Specss
    integer, parameter, public :: NEMIS_Specs = 48

end module ChemDims_mod
