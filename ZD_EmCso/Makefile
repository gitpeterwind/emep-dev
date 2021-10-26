# ------------------------------------------------
# defaults
# ------------------------------------------------

include ../Makefile


# ------------------------------------------------
# rules
# ------------------------------------------------

## compile files with debug flags:
#cso%.o:     cso%.F90
#	$(F90) $(F90FLAGS) $(DEBUG_FLAGS) -c $< -o $@
#EMEP_CSO_mod.o: EMEP_CSO_mod.f90
#	$(F90) $(F90FLAGS) $(DEBUG_FLAGS) -c $< -o $@


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
##~ enable MPI:
#DFLAGS  = -D_MPI

# ------------------------------------------------
# targets
# ------------------------------------------------

clean-local:
	rm -f cso*.{o,mod}
	rm -f EMEP_CSO_mod.{o,mod}



