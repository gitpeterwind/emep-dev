#
#
###################################################
include ../Makefile
###################################################
ifeq ($(EXP),Unimod)
  override EXP :=
endif
export CLEANCPP ?= 1
export MKL ?= 1
export FFTW ?=1

#stop: info-varaiable-EXP info-varaiable-PROG
ifdef EXP
  override DFLAGS += -DOLD_FFT -DOLD_SPL $(if $(filter yes,$(DEBUG)),-DDEBUG_B)
# override DFLAGS += -DBIAS_CORRECTION
  ifneq (,$(findstring OLD_FFT,$(DFLAGS)))
    $(PROG): fftpack5/libfftpack_r8.a
    LIBS += -lfftpack_r8
    LLIB += -L./fftpack5/
  else
    SRCS.$(EXP) += fftpack5.1.f90
  endif
  ifdef FFTW
    ifeq ($(MACHINE),vilje)
      MODULES += $(FFTW:1=fftw/3.3.3)
      ifneq (,$(call modulefind,fftw/3.3.3))
        FFTW3_ROOT = /sw/sdev/Modules/fftw/fftw-3.3.3
        LIBS += -lfftw3_mpi -lfftw3
        INCL += $(FFTW3_ROOT)/include
        LLIB += -L$(FFTW3_ROOT)/lib -Wl,-rpath,$(FFTW3_ROOT)/lib
      endif
    endif
  endif
  ifneq (,$(findstring OLD_SPL,$(DFLAGS)))
    SRCS.$(EXP) += ./B_NMC/fitpack.F
  else
    SRCS.$(EXP) += spl.f90
  endif
  ifdef MKL
    SRCS.$(EXP) += ./B_NMC/My_LAPACK77_ml.f90
    LIBS += -lmkl_lapack95_lp64 -lmkl_blas95_lp64 \
            -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
            -Wl,--end-group
    override DFLAGS += -DMKL
    MKLPATH=$(firstword $(subst :, ,$(MKLROOT) $(MKL_ROOT)))
    ifeq ($(MACHINE),stallo)
      MODULES += $(MKL:1=mkl/10.3.8)
      ifneq (,$(call modulefind,mkl/10.2.0.013))
        INCL += $(MKLPATH)/include $(MKLPATH)/include/em64t/lp64
        LLIB += -L$(MKLPATH)/lib/em64t
      else ifneq (,$(call modulefind,mkl/10.3.8))
        INCL += $(MKLPATH)/include $(MKLPATH)/include/intel64/lp64
        LLIB += -L$(MKL)/lib/intel64
      else ifneq (,$(call modulefind,mkl/11.2.137))
        MKLPATH=/global/apps/intel/XE_2011/mkl
        INCL += $(MKLPATH)/include
        LLIB += -L$(MKLPATH)/lib/intel64
      else
        $(call modulecheck,$(filter mkl%,$(MODULES)))
      endif
    else ifeq ($(MACHINE),vilje)
      MODULES += $(MKL:1=intelcomp/13.0.1)
#     MODULES += $(MKL:1=intelcomp/14.0.1)
      ifneq (,$(call modulefind,intelcomp/14.0.1))     
        INCL += $(MKLPATH)/include $(MKLPATH)/include/intel64/lp64
        LLIB += -L$(MKLPATH)/lib/intel64
      else ifneq (,$(call modulefind,intelcomp/13.0.1))
        INCL += $(MKLPATH)/include $(MKLPATH)/include/intel64/lp64
        LLIB += -L$(MKLPATH)/lib/intel64
      else ifneq (,$(call modulefind,intelcomp/12.0.5.220))
        INCL += $(MKLPATH)/include $(MKLPATH)/include/intel64/lp64
        LLIB += -L$(MKLPATH)/lib/intel64
      else
        $(call modulecheck,$(filter intelcomp%,$(MODULES)))
      endif
    endif
  else
    SRCS.$(EXP) += ./B_NMC/My_LAPACK77_ml.f90
    LIBS += -llapack -lblas #-Wl,--start-group -llapack -lblas -Wl,--end-group
    ifeq ($(MACHINE),stallo)
      MODULES += netlib/1.8.0
      LLIB += /global/apps/netlib/lib
    else ifeq ($(MACHINE),vilje)
    endif
  endif
  override SRCS:=$(SRCS.$(EXP))
  export SRCS
  depend: $(SRCS)
# stop: info-varaiable-SRCS
endif

# -J/-module option sets the path for the .mod files
# e.g. $(if $(filter gfortran,$(LD)),-J,-module )$(@D)/  # same path as source
F90FLAGS += $(if $(filter gfortran,$(LD)),-J,-module ) ../
SUFFIXES += .F .f
%.o %.mod:%.f90
	$(F90) $(F90FLAGS) -c $< -o $*.o
%.o:%.f
	$(F90) $(F90FLAGS) -c $< -o $@
%.o:%.F
	$(F90) $(F90FLAGS) -c $< -o $@

fftpack5/libfftpack_r8.a:
	$(MAKE) -C $(@D) intel

PHONY+=DA_ml.o My_Pollen_ml.o My_LAPACK77_ml.o \
  $(if $(filter 1,$(CLEANCPP)),\
    $(patsubst .*,.o,$(shell grep -l -s '\#if' $(SRCS))))
.NOTPARALLEL: $(NOTPARALLEL) EXP_Unimod EXP_3DVar EXP_Bnmc

allprog: CLEANCPP=1
allprog: EXP_Unimod EXP_3DVar EXP_Bnmc
EXP_Unimod%:
	$(MAKE) -C ../ MACC$*
EXP_%:
	$(MAKE) -C ./ -j4 EXP=$* PROG=$(@:EXP%=$(PROG)%)

$(PROG): depend modules info-varaiable-MAKECMDGOALS info-varaiable-SRCS
modules: $(foreach m,$(MODULES),check-module-$(subst /,_,$(m)))
info-varaiable-%:
	@echo -e "Variable: $*='$(value $*)'"
##########################################################
