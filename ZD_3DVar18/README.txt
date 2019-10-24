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


Commit new version
------------------

# for colors:
git config --global color.ui auto

cd emep-mscw

git status

# to replace files with unnecesary changes, 
# use "checkout" with the filename, this will discard all changes:
git checkout ZD_3DVar18/DA_Util_ml.f90

# introduce new version:
git add ZD_3DVar19

