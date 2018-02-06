# EMEP MSC-W open-source release rv4.17


The EMEP MSC-W code released as rv4.15 (and with an emep-internal update to rv4.16)
was the result of many changes compared to
previous versions, and partly as a result of this we have found both bugs and
areas of improvement in the model. We summarise below the main changes


# Chemistry

## Modifications

  Some of the gas-aerosol reactions introduced in rv4.15 have now been removed, either
  as they were negligable  (see Stadtler et al., 2017) or having a significant
  influence but causing model-performance degredation and being too uncertain to
  be reliable  (NO2+ .., ibid)

## bug fix

  OH + HONO  - rate coefficient had incorrect negative temperature dependence (small change)

# Timefactors

  Dave will likely use LOTOS for NH3 ... gives better results ;-)

# Landuse_ml

 Fixed a bug whereby ...

# Global landcover

  Deserts were included under the category 'BARE' instead of 'DE'. In principal this could
  have been fine, except that some 'BARE' land was also assigned where tew GLC2000 database
  had sparse vegetation. 

  To correct this, we have temporalily renamed 'BARE' from CLM to be 'DE', and reserved
  BARE for bare-soil assumed to be associated with vegetation 

# Tidy-ups

  A large number of tidy-ups have been done....

# Preliminary

  BiDir ...

