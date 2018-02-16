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

# config_emep.nml
  
    Pathes to all input files can now be steered by configuration settings.
    A short explanation for making own output fileds using 'USET' is given in the user guide.

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

# snow depth
  
  Modification of snow depth defintions: By default snow depth in the meteo fields is assumed to be in water equivalent, and multiplied by 5 to give snow depth. For WRF metdata, the field named SNOWH is used (was SNOWNC in earlier versions); the WRF snowdepth is assumed to be directly in snow meters (i.e. not multiplied by 5).

# femis.dat

  "lonlat" type reductions has an additional country flag.
  The "lonlat" reductions can also be dectivated for individual emission files, using the configuration flag emis_inputlist(1)%use_lonlat_femis = F
  
# Photolysis rates

  A height instead of level defintion is used. Should be more correct over mountainous terrain. Also in Southern hemisphere the wrong latitude was taken, now corrected.

# Preliminary

  BiDir ...

