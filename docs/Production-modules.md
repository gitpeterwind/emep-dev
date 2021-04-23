Projects like *CAMS50* and *CAMS71* require (daily/automated) operational forecast, analysis or source-receptor ensemble runs. For this the model binaries and configuration scripts are available on the target HPC as environment modules.

## *stratus.nsc.liu.se* suites
Setup files are tracked on the [metno_op@stratus_emep][] branch of the [modulefiles][] repository.

[modulefiles]: https://github.com/metno/modulefiles
[metno_op@stratus_emep]: https://github.com/metno/modulefiles/tree/metno_op%40stratus_emep

### *CAMS50* suites

Usage example for operational user (*metno_op*):

```bash
# load latest experimental suite
module purge
module use ~metno_op/emep/modulefiles
module load cams50/e-suite

# set up today's forecast, analysis and reanalysis runs
emep-config.py forecast analysis reanalysis

# submit the forecast
RUNDIR=`date +~/work/emep/cwf/CWF_12FC-%Y%m%d -d "today"`
$RUNDIR/run.sh

# submit the analysis
RUNDIR=`date +~/work/emep/cwf/CWF_00AN-%Y%m%d -d "yesterday"`
$RUNDIR/run.sh

# submit the reanalysis
RUNDIR=`date +~/work/emep/cwf/CWF_00RE-%Y%m%d -d "20 days ago"`
$RUNDIR/run.sh
```

## *nebula.nsc.liu.se* suites

Production modules are mirrored at *~sm_alvva/CAMS50/*, and tracked at [sm_alvva@nebula][].
The CAMS50 runs on the previous example can be setup and run after making the modules available with:

```bash
module use ~sm_alvva/CAMS50/modulefiles
```

[sm_alvva@nebula]: https://github.com/metno/modulefiles/tree/sm_alvva%40nebula
