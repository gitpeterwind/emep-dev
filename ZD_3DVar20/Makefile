# ------------------------------------------------
# defaults
# ------------------------------------------------

# target:
PROG ?= Unimod_3DVar

# ------------------------------------------------
# compiler and flags
# ------------------------------------------------

# defaults:
FC      = f90
F77     = f77
FFLAGS  = 
LDFLAGS = 
INCS    = 
LIBS    = 

# machine specific settings:
include Makefile.conf

# ------------------------------------------------
# implicit rules
# ------------------------------------------------

# how to form object files from F90 source:
%.o: %.F90
	$(FC) -c -o $@ $(FFLAGS) $(INCS) $<
	@echo ' '

# how to form object files from f90 source:
%.o: %.f90
	$(FC) -c -o $@ $(FFLAGS) $(INCS) $<
	@echo ' '

# how to form object files from F source:
%.o: %.F
	$(F77) -c -o $@ $(FFLAGS) $(INCS) $<
	@echo ' '

# how to form object files from F source:
%.o: %.f
	$(F77) -c -o $@ $(FFLAGS) $(INCS) $<
	@echo ' '

../NetCDF_mod.o: FFLAGS += -fpe-all=3

# ------------------------------------------------
# rules
# ------------------------------------------------

help:
	@echo ""
	@echo "Usage:"
	@echo "    make [help]"
	@echo "    make $(PROG)"
	@echo "    make clean"
	@echo ""

clean:
	rm -f *.mod $(foreach d,$(sort $(dir $(SRCS))),$d*.o)
	rm -f emepctm

# ------------------------------------------------
# dependencies
# ------------------------------------------------

include ../Makefile.SRCS
SRCS := $(addprefix ../,$(filter-out My_3DVar_mod.f90,$(SRCS))) \
	DA_ml.f90 Analysis/DA_3DVar_ml.f90 Analysis/DA_Obs_ml.f90 Analysis/EMIP.f90 Analysis/HESK_SuperObs.f90 \
	B_NMC/c3po_coordinates.F90 B_NMC/c3po_datafile.F90 B_NMC/c3po.F90 \
        B_NMC/linalg_lapack_350.F90 B_NMC/linalg.F90 B_NMC/emep_nmc_linalg.F90 \
	B_NMC/emep_bcovarsqrt.F90 \
	DA_Util_ml.f90 \
	tools/AJS.F90 \
	tools/go_date.F90 tools/go_print.F90 tools/go_string.F90 \
	tools/go_fu.F90 tools/go_par.F90 tools/go_timer.F90 \
        tools/go_path.F90 tools/go.F90 \
	tools/mpif90.F90 tools/fftw3.F90 tools/ufftw3.F90 \
	../ZD_3DVar/m1qn3-3.3-distrib/src/m1qn3-mpi.f
#	$(filter -DUSE_GO,$(DFLAGS),$(wildcard tools/go*.F90))

#
# ORIGINAL: 
#
#Makefile.deps: $(SRCS)
#	$(MAKEDEPF90) -o '$$(PROG)' $(SRCS) $(DFLAGS) > $@
#
#  The "-o '$$(PROG)'" flag will insert following lines:
#
#    FOBJ=../AeroConstants_mod.o ../AeroFunctions.o ...
#
#    $(PROG): $(FOBJ)
#           $(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(FOBJ) $(LIBS)
#
# The "FOBJ" macro has all possible objects (based on the arguments instered with the "$(SRC)" macro),
# also the objects that are actually not used. This is not very useful when for example only the
# observation operator is needed, without all assimilation (and covariance) stuff.
#
#>>>
#
# Create dependencies for all source files, do not add rule for main program:
Makefile.deps: $(SRCS)
	$(MAKEDEPF90) $(SRCS) $(DFLAGS) > $@

# object files from sources that are no modules, the dependencies on this are not trapped by makedepf90 ...
XOBJ  = ../global2local.o
XOBJ += $(if $(findstring with_assim,$(DFLAGS)), ../ZD_3DVar/m1qn3-3.3-distrib/src/m1qn3-mpi.o )

# main object:
MAINOBJ = ../emep_Main.o

# Rule for main program, only depends on main object, also ensure that non-module objects are created;
# link rule requires all possible object files:
$(PROG): $(XOBJ) $(MAINOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) `find .. -name '*.o'` $(LIBS)
#
#<<<

# dependencies created by 'makedepf90':
include Makefile.deps
