# EMEP MSC-W open-source release rv4.17


The EMEP MSC-W code released as rv4.15 (and with an emep-internal update to rv4.16)
was the result of many changes compared to
previous versions, and partly as a result of this we have found both bugs and
areas of improvement in the model. We summarise below the main changes


# Chemistry - EmChem16a

EmChem16a is an update of the EmChem16 scheme presented in EMEP report 1/2017.

## Modifications

  Some of the gas-aerosol reactions introduced in rv4.15 have now been removed, either
  as they were negligable  (see Stadtler et al., 2017) or having a significant
  influence but causing model-performance degredation and being too uncertain to
  be reliable  (NO2+ .., ibid)

## bug fixes

  OH + HONO - rate coefficient had incorrect negative temperature dependence (small change)
  OD + H2O  - updated rate coefficient from IUPAC 2007

# Timefactors

  For EECCA runs, the monthly time-factors from the LOTOS-EUROS model can now be
  activated, and seem to generate betetr results. 

  To use, set MonthlyNH3 = 'LOTOS' in config_emep.nml


# Landuse_ml

 Fixed varius bugs  whereby ...

# Global landcover

  NEW landcover file: glc2000xCLMf18.nc - set this in config_emep.nml.

  In the rv4.15 glc2000mCLM.nc file, deserts were included under the
  category 'BARE' instead of 'DE'. In principal this could have been
  fine, except that some 'BARE' land was also assigned where the GLC2000
  database had sparse vegetation. The EMEP model´s calculation of dust
  production was taking place only over the 'DE' regions, which were
  only found in the inner EECCA domain.


# Tidy-ups

  A large number of tidy-ups have been done....

# Preliminary

  BiDir ...

