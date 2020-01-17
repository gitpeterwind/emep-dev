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


