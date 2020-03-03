EMEP/ZD_3DVarNN extension
=========================

3D-var data assimilation.


Enable "ajs" test code
----------------------

Extra codes enabled using "with_ajs" macro,
for example used to have per-processor log files.

Make the following changes to enable this:
- Makefile
    In "SRCS" line, add:
      tools/AJS.F90
- Makefile.conf
    Add flag to define macro:
      DFLAGS  = -D_MPI -Dwith_lapack95_mkl -Dwith_ajs
- tools/AJS.F90
    Check settings.



Configuration
-------------

Template for CAMS50 NRT analysis settings:
  ./config_ANALYSIS.nml


Changes in ZD_3DVar20
---------------------

Module names in EMEP changed form "*_ml" to "*_mod".
  Analysis/DA_3DVar_ml.f90
  Analysis/DA_Obs_ml.f90
  DA_Util_ml.f90
  My_3DVar_ml.f90
  Makefile
  
Updated arguments of "find_index".
  Analysis/DA_Obs_ml.f90
  
Share "DEBUG_DA_1STEP" with EMEP-CTM.
  DA_ml.f90


SoilNOx project - S5p simulation only
-------------------------------------

Adapt assimilation code "ZD_3Dvar20" that only observation operator remains.
This is then used to simulate S5p observations including kernels.

Introduced macro "with_assim" to enable the 3D-var specific code,
if this is not defined only the observation operator code remains.
  Analysis/DA_3DVar_ml.f90
  Analysis/DA_Obs_ml.f90
  Makefile
  Makefile.conf

Restructuredd access to S5p observation files.
  Analysis/EMIP_OMI.f90   # renamed to "EMIP.f90"
  Analysis/EMIP.f90
  Analysis/HESK_SuperObs.f90

Removed copies of GO routines from DA_Util module.
Updated GO modulde for date processing.
  DA_Util_ml.f90
  tools/go_date.F90
  tools/AJS.F90
  B_NMC/c3po_coordinates.F90
  B_NMC/c3po_datafile.F90

Support simulation of S5p super observations.  
  Analysis/DA_Obs_ml.f90
  tools/go_path.F90
  tools/go.F90

Added script to create S5p superobs listing files.
  bin/s5p-superobs-create-listing

Added template settings.
  nml/config_CAMS81_S5p_superobs.nml

