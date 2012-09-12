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
           -cpp $(DFLAGS) $(addprefix -I,$(INCL))
LDFLAGS =  $(F90FLAGS) $(LLIB) $(LIBS)

MACHINE ?= stallo
DEBUG ?= no
ifeq ($(MACHINE),stallo)
# MODULES = intel-compiler/11.1   openmpi/1.4   netcdf/4.1.1
  MODULES = intel-compiler/12.1.2 openmpi/1.4.4 netcdf/4.1.3
  LIBS = -lnetcdf -lnetcdff
  INCL += $(NETCDF_ROOT)/include
  LLIB =  $(NETCDF_ROOT)/lib
  MAKEDEPF90=/home/mifapw/bin/makedepf90
  OPT_FLAGS = -O2 -ftz
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
else ifeq ($(MACHINE),vilje)
# MODULES = intelcomp/11.1.073   mpt/2.04 netcdf/4.1.3-intel.11.1.073
  MODULES = intelcomp/12.0.5.220 mpt/2.06 netcdf/4.1.3
  LIBS = -lnetcdf -lnetcdff
  INCL += $(NETCDF_PREFIX)/include
  LLIB  = $(NETCDF_PREFIX)/lib
  MAKEDEPF90=/home/metno/mifapw/bin/makedepf90
  LLIB := $(foreach L,$(LLIB),-L$(L) -Wl,-rpath,$(L))
else ifeq ($(MACHINE),titan)
  MODULES = intel/11.1u8 openmpi/1.4.3.intel
  INTEL = /site/VERSIONS/intel-11.1u8
  NETCDF = /usit/titan/u1/alfgr/netcdf-3.6.1
  LIBS = -lirc -lnetcdf
  INCL += $(NETCDF)/include $(INTEL)/include/intel64
  LLIB = -L$(NETCDF)/lib -L$(INTEL)/lib/intel64
  MAKEDEPF90=/usit/titan/u1/mifapw/bin/makedepf90
else ifeq ($(MACHINE),njord)
  LIBS =  -lnetcdf
  INCL += /home/ntnu/usrlocal/netcdf/netcdf-3.6.1/include
  LLIB = -L/home/ntnu/usrlocal/netcdf/netcdf-3.6.1/lib -L/usr/lib
  MAKEDEPF90=/home/ntnu/mifahik/local/bin/makedepf90
  F90 = mpxlf90_r
  F90FLAGS = -q64 -qrealsize=8 -O3 -qarch=pwr5 -qtune=pwr5 $(addprefix -I,$(INCL))
else ifeq ($(MACHINE),RSS)
  LIBS = -lnetcdf
  INCL += /usr/local/netcdf/include
  LLIB = -L/usr/local/netcdf/lib
  MAKEDEPF90 = /home/davids/local/bin/makedepf90
  LD = gfortran
  F90FLAGS = -fdefault-real-8 -O3 -Wall -ffixed-line-length-none -fbacktrace \
    -fbounds-check -pedantic -fimplicit-none $(addprefix -I,$(INCL))
else ifeq ($(MACHINE),EeePC)
  LIBS = -lnetcdf
  INCL += /home/davids/WRF_2009/netcdf4/include
  LLIB = -L/home/davids/WRF_2009/netcdf4/lib
  MAKEDEPF90 = /home/davids/local/bin/makedepf90
  LD = gfortran
  F90FLAGS = -fdefault-real-8 -O3 -Wall -ffixed-line-length-none -fbacktrace \
    -fbounds-check -pedantic -fimplicit-none $(addprefix -I,$(INCL))
else ifeq ($(MACHINE),TP)
  LIBS = -lnetcdf
  INCL += /home/davids/local/netcdf-4.0.1/include
  LLIB = -L/home/davids/local/netcdf-4.0.1/lib
  MAKEDEPF90 = /home/davids/local/bin/makedepf90
  LD = gfortran
  F90FLAGS = -fdefault-real-8 -O3 -Wall -ffixed-line-length-none -fbacktrace \
    -fbounds-check -pedantic -fimplicit-none $(addprefix -I,$(INCL))
else ifeq ($(MACHINE),hardy)  #ubuntu 8.04
  LIBS = -lnetcdff -lnetcdf
  INCL += /usr/include
  LLIB = -L/usr/lib
  MAKEDEPF90 = /disk1/emep_U804/local/bin/makedepf90
  LD = gfortran
  F90FLAGS = -fdefault-real-8 -O3 -Wall -ffixed-line-length-none -ffree-line-length-none \
    -fbounds-check -pedantic -fimplicit-none $(addprefix -I,$(INCL))
else ifeq ($(MACHINE),lucid)  #ubuntu 10.04
  LIBS = -lnetcdff -lnetcdf
  INCL += /usr/include
  LLIB = -L/usr/lib
  MAKEDEPF90 = /disk1/emep_U1004/local/bin/makedepf90
  LD = gfortran
  F90FLAGS = -fdefault-real-8 -O3 -Wall -ffixed-line-length-none -fbacktrace \
    -fbounds-check -pedantic -fimplicit-none $(addprefix -I,$(INCL))
endif
F90FLAGS += $(if $(filter yes,$(DEBUG)),$(DEBUG_FLAGS),$(OPT_FLAGS))

.SUFFIXES: $(SUFFIXES)  .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<


# Include the dependency-list created by makedepf90 below
all:  $(PROG)

ifndef MAKECMDGOALS
  include .depend
else ifneq (,$(filter all $(PROG) %.o,$(MAKECMDGOALS)))
  include .depend
endif

#
depend .depend: Makefile.SRCS $(SRCS)
	$(MAKEDEPF90) $(filter %.f90,$+) -o '$$(PROG)' $(DFLAGS) \
	  -l '$$(F90) -o $$@ $$+ $$(LDFLAGS)' > .depend

clean: diskclean touchdepend depend

diskclean:
	rm -f $(PROG) *.o *.mod

touchdepend:
	touch .depend

# Model/Config specific targets
EMEP EMEP2010 SR-EMEP SR-EMEP2010 eEMEP:
	ln -sf $(filter %.f90 %.inc,$+) . && \
	$(MAKE) MACHINE=$(MACHINE) -f $(firstword $(MAKEFILE_LIST)) -j4 $(PROG)
MACC MACC-EVA2010 SR-MACC eEMEP2010: ./ZD_OZONE/IFSMOZ_ExternalBICs_ml.f90
	ln -sf $(filter %.f90 %.inc,$+) . && \
	ln -sf IFSMOZ_ExternalBICs_ml.f90 My_ExternalBICs_ml.f90 && \
	$(MAKE) MACHINE=$(MACHINE) -f $(firstword $(MAKEFILE_LIST)) -j4 $(PROG)

# My_* files pre-requisites
EMEP EMEP2010 MACC MACC-EVA2010 eEMEP2010: \
	  ./ZD_OZONE/My_RunSettings.inc ./ZD_OZONE/My_Derived_ml.f90 \
	  ./ZD_OZONE/My_ExternalBICs_ml.f90 ./ZD_OZONE/My_Outputs_ml.f90 \
	  ./ZD_OZONE/My_Aerosols_ml.f90 ./ZD_VBS/My_SOA_ml.f90
#For SR we use the small My_Derived
SR-EMEP SR-EMEP2010 SR-MACC: \
	  ./ZD_SR/My_RunSettings.inc ./ZD_SR/My_Derived_ml.f90 \
	  ./ZD_OZONE/My_ExternalBICs_ml.f90 ./ZD_OZONE/My_Outputs_ml.f90 \
	  ./ZD_OZONE/My_Aerosols_ml.f90 ./ZD_VBS/My_SOA_ml.f90

eEMEP: \
	  ./ZD_OZONE/My_RunSettings.inc ./ZD_OZONE/My_Derived_ml.f90 \
	  ./ZD_OZONE/My_ExternalBICs_ml.f90 ./ZD_OZONE/My_Outputs_ml.f90 \
	  ./ZD_OZONE/My_Aerosols_ml.f90 ./ZD_OZONE/My_SOA_ml.f90
# GenChem config
.SECONDEXPANSION:
EMEP EMEP2010 SR-EMEP SR-EMEP2010 \
MACC MACC-EVA2010 SR-MACC eEMEP2010: modules $$@-GenChem-EmChem09soa
eEMEP: modules $$@-GenChem-EmChem09soa #-Emergency not yet ready
EMEP-GenChem-%:
	mk.GenChem -r $* -f FINNv1 -e SeaSalt,Dust,Isotopes
EMEP2010-GenChem-%:
	mk.GenChem -r $* -f FINNv1 -e SeaSalt,Dust,Isotopes -V 2bin,Eyjafj.ll #-h
MACC-GenChem-%:
	mk.GenChem -r $* -f GFED   -e SeaSalt,Dust,Isotopes
MACC-EVA2010-GenChem-%:
	mk.GenChem -r $* -f GFASv1 -e SeaSalt,Dust,Isotopes -V 2bin,Eyjafj.ll
SR-EMEP-GenChem-%:
	mk.GenChem -r $* -f FINNv1 -e none
SR-EMEP2010-GenChem-%:
	mk.GenChem -r $* -f FINNv1 -e none -V 2bin,Eyjafj.ll
SR-MACC-GenChem-%:
	mk.GenChem -r $* -f GFED   -e none
eEMEP-GenChem-%:
	mk.GenChem -r $* -f FINNv1 -e none -V 7bin,Vesuvius,Etna,Kr.suv.k,Katla,Askja #-h
eEMEP2010-GenChem-%:
	mk.GenChem -r $* -f FINNv1 -e SeaSalt,Dust,Isotopes -V 2bin,Eyjafj.ll #-h

# Check if intended modules are loaded
lstmodule = $(shell bash -c '. /etc/profile.d/modules.sh; module list' 2>&1)
fidmodule = $(findstring $(1),$(lstmodule))
chkmodule = $(if $(call fidmodule,$(1)),$(info Found module: $(1)),\
  $(error Missing module $(1): try 'module load $(1)'))
checkmodules = $(foreach m,$(subst _,/,$(1)),$(call chkmodule,$(m)))
check-module-%:
	$(call checkmodules,$*)
modules: $(foreach m,$(MODULES),check-module-$(subst /,_,$(m)))
##########################################################
