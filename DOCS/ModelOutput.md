# Output files

Most of the model output files are writte in netCDF format. 
The sellection of output variables and averaging/accumulation intervals
are set on different namelists in the configuration file: `config_emep.nml`.

## Output data files
- `*_fullrun.nc`: Gridded yearly  mean/total values.
- `*_month.nc`:   Gridded monthly mean/total values.
- `*_day.nc`:     Gridded daily   mean/total values.
- `*_hour.nc`:    Gridded hourly  mean/total values.
- `*_hourInst.nc`:  Gridded hourly instantaneous values.
- `*_hourExtra.nc`: Gridded legacy hourly outputs.
- `sites_*.nc`:  Hourly instantaneous surface values    on sellected locations.
- `sondes_*.nc`: Hourly instantaneous vertical profiles on sellected locations.

Gridded  files are named after the `runlabel1` parameter in `INPUT_PARA` namelist.
Location files are named after the year in the `startdate` parameter:
```Fortran
&INPUT_PARA
  GRID      = 'MACC14',
  iyr_trend = 2012,
  runlabel1 = 'EVA00AN-2012',
  runlabel2 = 'BM-EmChem09soa_svn3164_EVA00AN-201207_Trend2012',
  startdate = 2012,07,01,000000,
  enddate   = 2012,07,31,000000,
&end
```

## Output parameters

### Gridded output
```Fortran
&OutputConcs_config
OutputConcs=
  'O3'            ,'ug'  ,'2d' ,'AIR_CONCS','SPEC' ,'YMDI',
  'NO2'           ,'ug'  ,'2d' ,'AIR_CONCS','SPEC' ,'YMDI',
  'PMFINE'        ,'ug'  ,'2d' ,'AIR_CONCS','GROUP','',
  'PM10'          ,'ug'  ,'2d' ,'AIR_CONCS','GROUP','',
  'SURF_PM25water','ug'  ,'2d' ,'PM25water','MISC' ,'',
  'SURF_ug_PM25_rh50','ug','2d','PM25_rh50','MISC' ,'YMDI',
  'SURF_ug_PM10_rh50','ug','2d','PM10_rh50','MISC' ,'YMDI',
  'NO'            ,'ug'  ,'2d' ,'AIR_CONCS','SPEC' ,'YMDI',
  'SO2'           ,'ug'  ,'2d' ,'AIR_CONCS','SPEC' ,'YMDI',
  'CO'            ,'ug'  ,'2d' ,'AIR_CONCS','SPEC' ,'YMDI',
  'NH3'           ,'ug'  ,'2d' ,'AIR_CONCS','SPEC' ,'YMDI',
  'PAN'           ,'ug'  ,'2d' ,'AIR_CONCS','SPEC' ,'I',
  'NMVOC'         ,'ug'  ,'2d' ,'AIR_CONCS','GROUP','I',
  'POLLEN_BIRCH','grains','2d' ,'AIR_CONCS','SPEC' ,'I',
  'POLLEN_OLIVE','grains','2d' ,'AIR_CONCS','SPEC' ,'I',
  'POLLEN_GRASS','grains','2d' ,'AIR_CONCS','SPEC' ,'I',
  'POLLEN_RWEED','grains','2d' ,'AIR_CONCS','SPEC' ,'I',
  'O3'            ,'ug'  ,'3d' ,'AIR_CONCS','SPEC' ,'I',
  'NO2'           ,'ug'  ,'3d' ,'AIR_CONCS','SPEC' ,'I',
  'PMFINE'        ,'ug'  ,'3d' ,'AIR_CONCS','GROUP','',
  'PM10'          ,'ug'  ,'3d' ,'AIR_CONCS','GROUP','',
  'D3_PM25water'  ,'ug'  ,'3d' ,'PM25water','MISC' ,'',
  'D3_ug_PM25_wet','ug'  ,'3d' ,'PM25_wet' ,'MISC' ,'I',
  'D3_ug_PM10_wet','ug'  ,'3d' ,'PM10_wet' ,'MISC' ,'I',
  'NO'            ,'ug'  ,'3d' ,'AIR_CONCS','SPEC' ,'I',
  'SO2'           ,'ug'  ,'3d' ,'AIR_CONCS','SPEC' ,'I',
  'CO'            ,'ug'  ,'3d' ,'AIR_CONCS','SPEC' ,'I',
  'NH3'           ,'ug'  ,'3d' ,'AIR_CONCS','SPEC' ,'I',
  'PAN'           ,'ug'  ,'3d' ,'AIR_CONCS','SPEC' ,'I',
  'NMVOC'         ,'ug'  ,'3d' ,'AIR_CONCS','GROUP','I',
  'RN222'         ,'ug'  ,'2d' ,'AIR_CONCS','SPEC' ,'YMD',
  'NO2'           ,'ugm2','k20','COLUMN'   ,'MISC' ,'H',
  'O3'            ,'ugm2','k20','COLUMN'   ,'MISC' ,'H',
  'CO'            ,'ugm2','k20','COLUMN'   ,'MISC' ,'H',
  'HCHO'          ,'ugm2','k20','COLUMN'   ,'MISC' ,'H',
  'AOD'           ,' ' ,'550nm','AOD:GROUP','MISC' ,'H',
&end
&OutputDep_config
  DDEP_ECOS =
   'Grid'     , 'YMD',
  DDEP_WANTED =
   'SOX'      ,'GROUP','mgS',
   'OXN'      ,'GROUP','mgN',
   'RDN'      ,'GROUP','mgN',
   'O3'       ,'SPEC' , 'mg',
   'STO_O3'   ,'SPEC' , 'mg',
  WDEP_WANTED =
   'PREC'     ,'PREC' , 'mm',
   'SOX'      ,'GROUP','mgS',
   'OXN'      ,'GROUP','mgN',
   'RDN'      ,'GROUP','mgN',
   'SO2'      ,'SPEC' ,'mgS',
   'HNO3'     ,'SPEC' ,'mgN',
&end
&OutputSize_config
  num_lev3d =8,
  lev3d=1,2,3,4,6,9,10,12,
  lev3d_from_surface=T,
&end
```

### Legacy hourly output
Different sellection of legacy/deprecated hourly outputs (`*_hourExtra.nc`)
are preset `My_Derived_ml.f90`. The use of this outputs is disencuraged,
as most of the functionality can be found on the gridded hourly output.
If requred, one of these preset output sets can sellected
with the `MY_OUTPUTS` parameter in the `ModelConstants_config` namelist:
```Fortran
&ModelConstants_config
  EXP_NAME   = 'FORECAST',
  MY_OUTPUTS = 'MACC_ENS',
&end
```

### Sites and Sondes

Site and sonde outputs are also available on ASCII format with `*.csv` extension.
This format is easier to handle than netCDF, but it is unsuitable for large number of sites.
