#
#
PROG =	Unimod
###################################################

include Makefile.SRCS

###################################################

LIBS = -lnetcdf -lnetcdff
INCL = /global/apps/netcdf/4.1.1/include
LLIB = /global/apps/netcdf/4.1.1/lib

F90 = mpif90

F90FLAGS =  -shared-intel -CB -r8 -recursive -debug-parameters all -traceback -ftrapuv -g -fpe0 -O0 \
  -convert big_endian -IPF_fp_relaxed -I$(INCL)
F90FLAGS =  -shared-intel     -r8 -recursive -O2 -ftz \
  -convert big_endian -IPF_fp_relaxed -I$(INCL)
LDFLAGS =  $(F90FLAGS)  -L$(LLIB) -Wl,-rpath,$(LLIB) $(LIBS)

.SUFFIXES: $(SUFFIXES)  .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<


# Include the dependency-list created by makedepf90 below
all:  $(PROG)

include .depend

#
depend .depend: Makefile.SRCS $(SRCS)
	/home/mifapw/bin/makedepf90 $(filter %.f90,$+) \
	-o '$$(PROG)' -l '$$(F90) -o $$@ $$+ $$(LDFLAGS)' > .depend

clean: diskclean touchdepend depend

diskclean:
	rm -f $(PROG) *.o *.mod

touchdepend:
	touch .depend

# GenChem config for standard EMEP & MACC runs
.SECONDEXPANSION:
EMEP SR-EMEP: $$@-GenChem-EmChem09soa \
          ./ZD_OZONE/My_ExternalBICs_ml.f90 \
	  ./ZD_OZONE/My_Derived_ml.f90 ./ZD_OZONE/My_Outputs_ml.f90 \
	  ./ZD_OZONE/My_Aerosols_ml.f90 ./ZD_VBS/My_SOA_ml.f90
	ln -sf $(filter %.f90,$+) . && \
	$(MAKE) -j4 modules $(PROG)
MACC SR-MACC: $$@-GenChem-EmChem09soa \
	  ./ZD_OZONE/IFSMOZ_ExternalBICs_ml.f90 \
	  ./ZD_OZONE/My_Derived_ml.f90 ./ZD_OZONE/My_Outputs_ml.f90 \
	  ./ZD_OZONE/My_Aerosols_ml.f90 ./ZD_VBS/My_SOA_ml.f90
	ln -sf $(filter %.f90,$+) . && \
	ln -sf IFSMOZ_ExternalBICs_ml.f90 My_ExternalBICs_ml.f90 && \
	$(MAKE) -j4 modules $(PROG)

EMEP-GenChem-%:
	mk.GenChem -r $* -f FINNv1 -e SeaSalt,Dust,Isotopes
MACC-GenChem-%:
	mk.GenChem -r $* -f GFED   -e SeaSalt,Dust,Isotopes
SR-EMEP-GenChem-%:
	mk.GenChem -r $* -f FINNv1 -e none
SR-MACC-GenChem-%:
	mk.GenChem -r $* -f GFED   -e none

# Check if intended modules are loaded
MODULES = intel-compiler/11.1 openmpi/1.4 netcdf/4.1.1
fidmodule = $(findstring $(1),$(shell bash -c 'module list' 2>&1))
chkmodule = $(if $(call fidmodule,$(1)),$(info Found module: $(1)),\
  $(error Missing module $(1): try 'module load $(1)'))
checkmodules = $(foreach m,$(subst _,/,$(1)),$(call chkmodule,$(m)))
check-module-%:
	$(call checkmodules,$*)
modules: $(foreach m,$(MODULES),check-module-$(subst /,_,$(m)))
##########################################################
