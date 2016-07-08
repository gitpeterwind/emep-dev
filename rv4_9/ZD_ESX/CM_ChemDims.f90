module ChemDims
!! Basic dimensions of the chemical scheme
  integer, parameter, public  :: &
      NRCT = 9             &!! No. rate coefficients
     ,NSPEC_TOT = 14       &!! No. species (total)
     ,NSPEC_ADV = 9        &!! No. advected species
     ,NSPEC_SHL = 5         !! No. short-lived species
end module ChemDims

