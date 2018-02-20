# EMEP MSC-W open-source release rv4.17


The EMEP MSC-W code released as rv4.15 (and with an emep-internal update
to rv4.16) was the result of many changes compared to previous versions,
and partly as a result of this we have found both bugs and areas of
improvement in the model. We summarise below the main changes made in
the rv4.17 (and previous interim rv4.16) versions.


# Land-cover issues

## Landuse_ml

 The rv4.15 code had several issues, especially when trying to run
 the  EMEP model in domains away from Europe (e.g. Asia). The new code
 improves initialisations, and can cope with the omission of the Euriopean
 landcover data.

 One issue remains: the code still expects two landcover files to be
 named in the config settings. If not using the default:

```
 LandCoverInputs%MapFile =
    'DataDir/Landuse/Landuse_PS_5km_LC.nc', 
    'DataDir/LandInputs_Feb2018/glc2000xCLMf18.nc',
```

One can use a dummy 2nd file:

```
 LandCoverInputs%MapFile =
    'DataDir/LandInputs_Feb2018/glc2000xCLMf18.nc',
    'DataDir/LandInputs_Feb2018/glc2000xCLMf18.nc',
```




## Global landcover

  NEW landcover file: glc2000xCLMf18.nc - set this in config_emep.nml.

  In the rv4.15 glc2000mCLM.nc file, deserts were included under the
  category 'BARE' instead of 'DE'. In principal this could have been
  fine, except that some 'BARE' land was also assigned where the GLC2000
  database had sparse vegetation. The EMEP model´s calculation of dust
  production was taking place only over the 'DE' regions, which were
  only found in the inner EECCA domain.

## Biogenics_ml

  The GetEuroBVOC routine in Biogenics_ml.f90 needed the landcover codes
   CF,DF,NF and BF to be defined. Those are defined by the European land-cover
   file Landuse_PS_5km_LC.nc. Users running for e.g. Asian domains use only
   the global land-cover file, which does not have these codes. The rv4.17 
   update fixes this issue.


# Radiation_ml

  As of version rv4.16 of the model, the radiation scheme of Weiss & Norman (1985)
  is used to calculate PAR levels. This gives a better treatment of diffuse and
  direct radiation compared to previous versions. 
  

# Chemistry - EmChem16a

EmChem16a is an update of the EmChem16 scheme presented in EMEP report 1/2017.

## Modifications

  Version rv4.15 of the EMEP model introduced a number of gas-aerosol
  reactions as part of the global study of Stadtler et al., 2018.
  Some of these gas-aerosol reactions have now been removed, either
  as they were negligable  (see Stadtler et al., 2018) or having a
  significant influence but causing model-performance degredation and
  being too uncertain to be reliable (NO2+aerosol). The main impact
  compared to rv4.16 is to raise NO2 concentrations to some extent.

  As part of the rv4.16 update, dry and wet deposition of N2O5, using
  the same rate as HNO3, was added to the chemical mechanism.

## bug fixes

  OH + HONO - rate coefficient had incorrect negative temperature dependence (small change)
  OD + H2O  - updated rate coefficient from IUPAC 2007

# config_emep.nml
  
    Pathes to all input files can now be steered by configuration settings.
    A short explanation for making own output fileds using 'USET' is given in the user guide.

    The n2o5Hydrolysi  method was incorrectly set in the 2017 open-source. It should be 'Smix',
    which is now the default in Config_module.f90.

# Timefactors

  For EECCA runs, the monthly time-factors from the LOTOS-EUROS model can now be
  activated, and seem to generate betetr results. 

  To use, set MonthlyNH3 = 'LOTOS' in config_emep.nml

# Tidy-ups

  A large number of tidy-ups have been done, and are underway. For example, the code
  used a lot of configutation variables (set via config_emep.nml) with the prefix
  USES\_. Each of these had to be named in the namelist settings. Now, the more 
  general user-defined variable USES covers many of these cases, with e.g.
  USES%FOREST_FIRES, USES%DEGREEDAY_FACTORS

# snow depth
  
  Modification of snow depth defintions: By default snow depth in the
  meteo fields is assumed to be in water equivalent, and multiplied by
  5 to give snow depth. For WRF metdata, the field named SNOWH is used
  (was SNOWNC in earlier versions); the WRF snowdepth is assumed to be
  directly in snow meters (i.e. not multiplied by 5).

# femis.dat

  "lonlat" type reductions has an additional country flag.

  The "lonlat" reductions can also be dectivated for individual emission
  files, using the configuration flag emis_inputlist(1)%use_lonlat_femis
  = F

# Radiation

  A height instead of level defintion is used. Should be more correct
  over mountainous terrain. Also in Southern hemisphere the wrong latitude
  was taken, now corrected.

# Acknowledgements

Some of the improvements made in the code stem from issues raised by EMEP users, either
via the github issues sides or personal communication. Thanks are due to 
John Johansson and Robert Bergström (Chalmers, Sweden), Roy Wichink Kruit (RIVM),
Paul Hamer, and Massimo Vieno (CEH, Scotland) for such feedback. 


# References

Stadtler, S., Simpson, D., Schröder, S., Taraborrelli, D., Bott, A., and Schultz, M.: Ozone Impacts of Gas-Aerosol Uptake in Global Chemistry Transport Models, Atmos. Chem. Phys., in press, 2018

Weiss, A. and Norman, J. M.: Partitioning Solar-radiation into Direct and Diffuse, Visible and Near-infrared Components, Agricultural and Forest Meteorology, 34, 205–213, doi:10.1016/0168-1923(85)90020-6, 1985.


