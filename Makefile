#
#
export PROG ?= $(if $(BINDIR),$(BINDIR)/)emepctm
###################################################

include Makefile.SRCS

###################################################

F90 = mpif90
DEBUG_FLAGS = -check all -check noarg_temp_created -debug-parameters all \
              -traceback -ftrapuv -g -fpe0 -O0 -fp-stack-check
OPT_FLAGS = -O2 -ftz
F90FLAGS =  -r8  -IPF_fp_relaxed -assume noold_maxminloc
LDFLAGS =  $(F90FLAGS) $(LLIB) $(LIBS)

export MACHINE ?= stallo
export DEBUG ?= no
export ARCHIVE ?= no
ifeq ($(MACHINE),fram)
  MODULES = netCDF-Fortran/4.4.5-iimpi-2019a
  LDFLAGS +=  $(shell nc-config --flibs)
  F90FLAGS += $(shell nc-config --cflags)
  MAKEDEPF90=/cluster/projects/nn2890k/bin/makedepf90
  OPT_FLAGS = -O2 -ftz
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
  F90=mpiifort
else ifeq ($(MACHINE),stallo)
  MODULES = netCDF-Fortran/4.4.4-intel-2016b
  LDFLAGS +=  $(shell nc-config --flibs)
  F90FLAGS += $(shell nc-config --cflags)
  MAKEDEPF90=/home/mifapw/bin/makedepf90
  OPT_FLAGS = -O2 -ftz
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
  F90=mpiifort
else ifeq ($(MACHINE),gstallo)
  # Needs module swap intel gcc/4.7.2
  MODULES = gcc/4.7.2 openmpi/1.6.2 netcdf/4.2.1.1
  F90FLAGS = -fbacktrace -fdefault-real-8 -O3 -Wall \
    -ffixed-line-length-none -ffree-line-length-none \
    -fbounds-check -pedantic -fimplicit-none
  LDFLAGS += $(shell nc-config --flibs)
  F90FLAGS+= $(shell nc-config --cflags)
  MAKEDEPF90=/home/mifapw/bin/makedepf90
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
  DEBUG_FLAGS =
  OPT_FLAGS =
  FC=mpif90
  LD=mpif90
else ifeq ($(MACHINE),vilje)
  MODULES ?= intelcomp/13.0.1 mpt/2.06 netcdf/4.3.0 # fftw/3.3.3
# MODULES ?= intelcomp/14.0.1 mpt/2.09 netcdf/4.3.1 # fftw/3.3.4
# MODULES ?= intelcomp/15.0.1 mpt/2.10 netcdf/4.3.2 # fftw/3.3.4
# MODULES ?= intelcomp/16.0.1 mpt/2.13 netcdf/4.4.0 # fftw/3.3.4
  LDFLAGS += $(shell nc-config --flibs)
  F90FLAGS+= $(shell nc-config --cflags)
  MAKEDEPF90=/home/metno/mifapw/bin/makedepf90
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
else ifeq ($(MACHINE),byvind)
  MODULES = intel/12.1.0 openmpi/1.4.1-i101011 netcdf/4.1.2-i1210
  LIBS += -lnetcdf -lnetcdff
  INCL += /software/apps/netcdf/4.1.2/i1210/include
  LLIB += /software/apps/netcdf/4.1.2/i1210/lib
# MAKEDEPF90=????
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
else ifneq (,$(findstring $(MACHINE),frost alvin elvis))
  MODULES = buildenv-intel/2015-1 hdf5/1.8.14-i1501 netcdf/4.3.2-i1501-hdf5-1.8.14
  LIBS += -lnetcdf -lnetcdff
  INCL += /software/apps/netcdf/4.3.2/i1501-hdf5-1.8.14/include/
  LLIB += /software/apps/netcdf/4.3.2/i1501-hdf5-1.8.14/lib/
  MAKEDEPF90=makedepf90
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
  LD := mpif90
  LDFLAGS += -Nmpi
  F90FLAGS += -Nmpi
else ifneq (,$(findstring $(MACHINE),stratus nebula))
# MODULES = buildenv-intel/2018a-eb netCDF-HDF5/4.3.2-1.8.12-nsc1-intel-2018.u1-bare FFTW/3.3.6-nsc1
  MODULES = buildenv-intel/2018.u1-bare netCDF-HDF5/4.3.2-1.8.12-nsc1-intel-2018.u1-bare
  LDFLAGS += $(shell nf-config --flibs)
  F90FLAGS+= $(shell nf-config --fflags)
  MAKEDEPF90=makedepf90
else ifeq ($(MACHINE),abel)
  MODULES = intel/2011.10 openmpi.intel/1.6.1 netcdf.intel/4.2.1.1
  INTEL  = /cluster/software/VERSIONS/$(subst /,-,$(filter intel%,$(MODULES)))
  NETCDF = /cluster/software/VERSIONS/$(subst /,-,$(filter netcdf%,$(MODULES)))
  LIBS += -lirc -lnetcdf
  INCL += $(NETCDF)/include $(INTEL)/include/intel64
  LLIB += -L$(NETCDF)/lib -L$(INTEL)/lib/intel64
  MAKEDEPF90=/usit/$(MACHINE)/u1/mifapw/bin/makedepf90
else ifeq ($(MACHINE),xenial)  # ubuntu 16.04
  # sudo apt-get install makedepf90 libmpich-dev libnetcdf-dev libnetcdff-dev
  F90FLAGS = -fdefault-real-8 -fdefault-double-8 -ffixed-line-length-none -ffree-line-length-none -fno-range-check
  LDFLAGS += $(shell nf-config --flibs)
  F90FLAGS+= $(shell nf-config --cflags)
  MAKEDEPF90 = /usr/bin/makedepf90
  LD = gfortran
  DEBUG_FLAGS = -Wall -fbacktrace -fbounds-check -fimplicit-none -pedantic
  OPT_FLAGS = -O3
else ifeq ($(MACHINE),ppixenial) # ubuntu 16.04, ifort
  LDFLAGS +=  $(shell nc-config --flibs)
  F90FLAGS += $(shell nc-config --cflags)
  MAKEDEPF90=makedepf90
  OPT_FLAGS = -O2 -ftz
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
endif
F90FLAGS += -cpp $(DFLAGS) $(addprefix -I,$(INCL)) \
   $(if $(filter yes,$(DEBUG)),$(DEBUG_FLAGS),$(OPT_FLAGS))

.SUFFIXES: $(SUFFIXES) .f90
%.o:%.f90
	$(F90) $(F90FLAGS) -c $< -o $@

# disable div0 exeption (DEBUG=yes) on netcdf/4.3.1 .. netcdf/4.4.0
ifneq ($(LD),gfortran)
NetCDF_mod.o:NetCDF_mod.f90
	$(F90) $(F90FLAGS) -fpe-all=3 -c $< -o $@
endif

# Include the dependency-list created by makedepf90 below
all:  $(PROG)
$(PROG): .depend
include .depend

# File dependencies
.depend: depend
depend: Makefile Makefile.SRCS $(SRCS)
	test -n "$(MAKEDEPF90)" && $(MAKEDEPF90) $(SRCS) $(DFLAGS) \
	  -o '$$(PROG)' -l '$$(F90) -o $$@ $$(FOBJ) $$(LDFLAGS)' > .$@
.version: version
version: Makefile Makefile.SRCS $(SRCS)
	svnversion -n > .$@

clean: diskclean touchdepend depend

diskclean:
	rm -f $(PROG) $(foreach d,$(sort $(dir $(SRCS))),$d*.o $d*.mod)

touchdepend:
	touch .depend

###
# Model/Config specific targets
###

# Test
TEST:
	$(MAKE) -j4 PROG=ModuleTester DEBUG=yes \
	  SRCS="$(filter-out emep_Main.f90,$(SRCS)) ModuleTester.f90"

# eEMP Default Scenarios: Vents, NPPs & NUCs
Emergency: VENTS ?= DefaultVolcano
#Eyjafjoll,Vesuvius,Etna,Kr.suv.k,Katla,Askja
Emergency: NPPAS ?= 
#Olkiluoto,Loviisa,Kola,Leningrad,Ringhals,Forsmark,Oskarshamn,Torness,Sellafield
Emergency: NUCXS ?= 
#NorthKorea,Tehran
Emergency:
	ZCM_Emergency/mk.Emergency -V 7bin,$(VENTS)
#-N $(NPPAS) -X $(NUCXS)

# eEMP Default AshInversion: Vents
AshInversion: VENTS ?= Eyjafjoll
AshInversion:
	ZCM_Emergency/mk.Emergency -V 19lev,9bin,$(VENTS)

# Data assimilation: 3DVar
3DVar16 3DVar17 3DVar20: EmChem19p
	$(MAKE) -C ZD_$@/ $(PROG) DFLAGS="-D_MPI -Dwith_lapack95_mkl -Dwith_assim"

# Chemistry version
EmChem%:
	rsync --checksum ZCM_$@/CM*.{inc,f90} .

# Always re-make this targets
.PHONY: $(PHONY) all depend modules

# Check if intended modules are loaded
LOADEDMODULES ?= $(shell bash -c 'module list' 2>&1)
modulelist = $(subst :, ,$(LOADEDMODULES))
modulefind = $(filter $(1),$(modulelist))
modulecheck= $(if $(call modulefind,$(1)),$(info Found module: $(1)),\
  $(error Missing module $(1): try 'module load $(1)'))
checkmodules = $(foreach m,$(subst _,/,$(1)),$(call modulecheck,$(m)))
check-module-%:
	$(call checkmodules,$*)
modules: $(foreach m,$(MODULES),check-module-$(subst /,_,$(m)))
	@echo "Loaded modules: $(modulelist)"
##########################################################
