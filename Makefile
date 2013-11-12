#
#
PROG =	Unimod
###################################################

include Makefile.SRCS

###################################################

F90 = mpif90
DEBUG_FLAGS = -CB -debug-parameters all -traceback -ftrapuv -g -fpe0 -O0
OPT_FLAGS = -O3 -ftz
F90FLAGS = -shared-intel -r8 -recursive -convert big_endian -IPF_fp_relaxed \
           -cpp $(DFLAGS)
LDFLAGS =  $(F90FLAGS) $(LLIB) $(LIBS)

MACHINE ?= stallo
DEBUG ?= no
ifeq ($(MACHINE),stallo)
  MODULES = intel/13.0 openmpi/1.6.2 netcdf/4.2.1.1
  LIBS += -lnetcdf -lnetcdff
  INCL += $(NETCDF_ROOT)/include
  LLIB += $(NETCDF_ROOT)/lib
  MAKEDEPF90=/home/mifapw/bin/makedepf90
  OPT_FLAGS = -O2 -ftz
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
else ifeq ($(MACHINE),gstallo)
  # Needs module swap intel gcc/4.7.2
  MODULES = gcc/4.7.2 openmpi/1.6.2 netcdf/4.2.1.1
  LIBS += -lnetcdf -lnetcdff
  INCL += $(NETCDF_ROOT)/include $(filter /global/apps/openmpi%,$(subst :, ,$(CPATH)))
  LLIB += $(NETCDF_ROOT)/lib
  MAKEDEPF90=/home/mifapw/bin/makedepf90
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
  DEBUG_FLAGS =
  OPT_FLAGS =
  F90FLAGS = -fbacktrace -fdefault-real-8 -O3 -Wall -ffixed-line-length-none -ffree-line-length-none \
    -fbounds-check -pedantic -fimplicit-none
  FC=mpif90
  LD=mpif90
else ifeq ($(MACHINE),vilje)
  MODULES = intelcomp/13.0.1 mpt/2.06 netcdf/4.3.0
  LIBS += -lnetcdf -lnetcdff
  INCL += $(NETCDF_PREFIX)/include
  LLIB += $(NETCDF_PREFIX)/lib
  MAKEDEPF90=/home/metno/mifapw/bin/makedepf90
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
else ifeq ($(MACHINE),byvind)
  MODULES = intel/12.1.0 openmpi/1.4.1-i101011 netcdf/4.1.2-i1210
  LIBS += -lnetcdf -lnetcdff
  INCL += /software/apps/netcdf/4.1.2/i1210/include
  LLIB += /software/apps/netcdf/4.1.2/i1210/lib
# MAKEDEPF90=????
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
else ifeq ($(MACHINE),abel)
  MODULES = intel/2011.10 openmpi.intel/1.6.1 netcdf.intel/4.2.1.1
  INTEL  = /cluster/software/VERSIONS/$(subst /,-,$(filter intel%,$(MODULES)))
  NETCDF = /cluster/software/VERSIONS/$(subst /,-,$(filter netcdf%,$(MODULES)))
  LIBS += -lirc -lnetcdf
  INCL += $(NETCDF)/include $(INTEL)/include/intel64
  LLIB += -L$(NETCDF)/lib -L$(INTEL)/lib/intel64
  MAKEDEPF90=/usit/$(MACHINE)/u1/mifapw/bin/makedepf90
else ifeq ($(MACHINE),precise)  #ubuntu 12.04
  LIBS += -lnetcdff -lnetcdf
  INCL +=   /usr/include 
  LLIB +=   -L/usr/lib  
  MAKEDEPF90 = $(EMEPLOCAL)/bin/makedepf90
  LD = gfortran
  DEBUG_FLAGS = -Wall -fbacktrace -fbounds-check -pedantic 
  OPT_FLAGS = -O3  
  F90FLAGS = -fdefault-real-8 -ffixed-line-length-none -ffree-line-length-none -fimplicit-none
endif
F90FLAGS += $(addprefix -I,$(INCL)) \
   $(if $(filter yes,$(DEBUG)),$(DEBUG_FLAGS),$(OPT_FLAGS))

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<


# Include the dependency-list created by makedepf90 below
all:  $(PROG)
$(PROG): .depend

ifndef MAKECMDGOALS
  include .depend
else ifneq (,$(filter all $(PROG) %.o,$(MAKECMDGOALS)))
  include .depend
endif

#
.depend depend: Makefile Makefile.SRCS $(SRCS)
	test -n "$(MAKEDEPF90)" && $(MAKEDEPF90) $(SRCS) $(DFLAGS) \
	  -o '$$(PROG)' -l '$$(F90) -o $$@ $$(FOBJ) $$(LDFLAGS)' > .depend

clean: diskclean touchdepend depend

diskclean:
	rm -f $(PROG) $(foreach d,$(sort $(dir $(SRCS))),$d*.o $d*.mod)

touchdepend:
	touch .depend

# Model/Config specific targets
EMEP EMEP2010 EMEP2011 SR-EMEP SR-EMEP2010 SR-EMEP2011 EmChem09 CRI_v2_R5 \
MACC MACC-EVA2010 MACC-EVA2011 SR-MACC eEMEP eEMEP2010 eEMEP2013:
	ln -sf $(filter %.f90 %.inc,$+) . && \
	$(MAKE) MACHINE=$(MACHINE) DEBUG=$(DEBUG) -j4 $(PROG)

# My_* files pre-requisites
EMEP EMEP2010 EMEP2011 MACC MACC-EVA2010 MACC-EVA2011 eEMEP2010: \
	  ./ZD_OZONE/My_Derived_ml.f90 ./ZD_OZONE/My_Outputs_ml.f90 \
	  ./ZD_OZONE/My_Aerosols_ml.f90 ./ZD_VBS/My_SOA_ml.f90 \
	  ./ZD_3DVar/My_3DVar_ml.f90 ./ZD_Pollen/My_Pollen_ml.f90
# no SOA:
EmChem09 CRI_v2_R5 eEMEP eEMEP2013: \
	  ./ZD_OZONE/My_Derived_ml.f90 ./ZD_OZONE/My_Outputs_ml.f90 \
	  ./ZD_OZONE/My_Aerosols_ml.f90 ./ZD_OZONE/My_SOA_ml.f90 \
	  ./ZD_3DVar/My_3DVar_ml.f90 ./ZD_Pollen/My_Pollen_ml.f90
#For SR we use the small My_Derived
SR-EMEP SR-EMEP2010 SR-EMEP2011 SR-MACC: \
	  ./ZD_SR/My_Derived_ml.f90 ./ZD_OZONE/My_Outputs_ml.f90 \
	  ./ZD_OZONE/My_Aerosols_ml.f90 ./ZD_VBS/My_SOA_ml.f90 \
	  ./ZD_3DVar/My_3DVar_ml.f90 ./ZD_Pollen/My_Pollen_ml.f90
# Pollen for MACC FC runs
MACC: SRCS := $(filter-out My_Pollen_ml.f90,$(SRCS)) Pollen_ml.f90 Pollen_const_ml.f90
MACC: ./ZD_Pollen/Pollen_ml.f90 ./ZD_Pollen/Pollen_const_ml.f90 | depend

# GenChem config
.SECONDEXPANSION:
EMEP EMEP2010 EMEP2011:             GenChem-EMEP-EmChem09soa
SR-EMEP SR-EMEP2010 SR-EMEP2011:    GenChem-SR-EMEP-EmChem09soa
EmChem09 CRI_v2_R5:                 GenChem-EMEP-$$@
MACC SR-MACC:                       GenChem-$$@-EmChem09soa
eEMEP:                              GenChem-$$@-EmChem09     # GenChem-Emergency not yet ready
MACC-EVA2010 MACC-EVA2011:          GenChem-MACCEVA-EmChem09soa
eEMEP2010:                          GenChem-EMEP-EmChem09soa
eEMEP2013:                          GenChem-SR-MACC-EmChem09soa

GenChem%:
	mk.GenChem $(GenChemOptions) -q #-h
GenChem-%:          GenChemOptions += -r $(lastword $(subst -, ,$*))
GenChem-EMEP-%:     GenChemOptions += -f FINNv1 -e SeaSalt,Dust,Isotopes
GenChem-SR-EMEP-%:  GenChemOptions += -f FINNv1 -e none
GenChem-MACC-%:     GenChemOptions += -f GFASv1 -e SeaSalt,Dust,Isotopes,Pollen
GenChem-SR-MACC-%:  GenChemOptions += -f GFASv1 -e none
GenChem-MACCEVA-%:  GenChemOptions += -f GFASv1 -e SeaSalt,Dust
GenChem-eEMEP-%:    GenChemOptions += -f GFASv1 -e SeaSalt,Dust

# Emergency options
EMEP2010 SR-EMEP2010 MACC-EVA2010:  GenChemOptions += -V 2bin,Eyjafjoll
EMEP2011 SR-EMEP2011 MACC-EVA2011:  GenChemOptions += -V 2bin,Grimsvotn
eEMEP2010:                          GenChemOptions += -V 2bin,Eyjafjoll
eEMEP2013:                          GenChemOptions += -N NorthKorea

# eEMP Default Vents, NPPs & NUCs
eEMEP: VENTS ?= Vesuvius,Etna,Kr.suv.k,Katla,Askja
eEMEP: NPPAS ?= Olkiluoto,Loviisa,Kola,Leningrad,Ringhals,Forsmark,Oskarshamn,Torness,Sellafield
eEMEP: NUCXS ?= NorthKorea,Tehran
eEMEP: GenChemOptions += -V 7bin,$(VENTS) -N $(NPPAS) -X $(NUCXS)

# Data assimilation: Bnmc / 3DVar
%-Bnmc %-3DVar: $$*
	$(MAKE) MACHINE=$(MACHINE) DEBUG=$(DEBUG) -C ZD_3DVar/ $(@:$*-%=EXP_%)

# Archive: create $(PROG).tar.bz2
$(PROG) archive: $(PROG)_$(shell date +%Y%m%d).tar.bz2
%.tar.bz2: $(SRCS) Makefile Makefile.SRCS .depend $(wildcard *.inc *.pl mk.* *.nml)
	@echo "Creating archive $@"; tar --dereference -cjf $@ $+

# Check if intended modules are loaded
modulelist = $(if $(filter byvind,$(MACHINE)),/etc/cmod/bash.init,/etc/profile.d/modules.sh)
modulelist:= $(shell bash -c '. $(modulelist); module list' 2>&1)
modulefind = $(findstring $(1),$(modulelist))
modulecheck= $(if $(call modulefind,$(1)),$(info Found module: $(1)),\
  $(error Missing module $(1): try 'module load $(1)'))
checkmodules = $(foreach m,$(subst _,/,$(1)),$(call modulecheck,$(m)))
check-module-%:
	$(call checkmodules,$*)
modules: $(foreach m,$(MODULES),check-module-$(subst /,_,$(m)))
##########################################################
