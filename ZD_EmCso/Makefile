# ------------------------------------------------
# defaults
# ------------------------------------------------

include ../Makefile

# ------------------------------------------------
# libraries
# ------------------------------------------------

ifneq (,$(findstring $(MACHINE),stratus nebula))
  MODULES   = buildenv-intel/2018a-eb \
              netCDF-HDF5/4.3.2-1.8.12-nsc1-intel-2018.u1-bare
else
  $(error Unknown machine $(MACHINE))
endif

# ------------------------------------------------
# macro definitions
# ------------------------------------------------
#   _MPI                 # enables MPI code, used in the past for testing with serial version
# ------------------------------------------------
#~ enable MPI:
DFLAGS  = -D_MPI


