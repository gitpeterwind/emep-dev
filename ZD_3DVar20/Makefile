# ------------------------------------------------
# defaults
# ------------------------------------------------

include ../Makefile

# ------------------------------------------------
# libraries
# ------------------------------------------------

ifneq (,$(findstring $(MACHINE),stratus nebula))
  MODULES   = buildenv-intel/2018a-eb netCDF-HDF5/4.3.2-1.8.12-nsc1-intel-2018.u1-bare FFTW/3.3.6-nsc1
  ## FFTW
  FFTW_ROOT = $(shell pkg-config --variable=prefix fftw3)
  INCL     += $(FFTW_ROOT)/include 
  LDFLAGS  += -L$(FFTW_ROOT)/lib -lfftw3_mpi -lfftw3 -Wl,-rpath,$(FFTW_ROOT)/lib
  ## LAPACK/BLAS
  INCL 	   += $(MKLROOT)/include/intel64/lp64
  LDFLAGS  += -L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 \
		      -lmkl_sequential -lmkl_core -Wl,-rpath,$(MKLROOT)/lib/intel64
else ifneq (,$(findstring $(MACHINE),stratus2 nebula2))
  MODULES   = makedepf90/2.8.8-hpc1-gcc-2022a-eb buildenv-intel/2023a-eb netCDF-HDF5/4.9.2-1.12.2-hpc1 FFTW/3.3.9
  ## FFTW
  FFTW_ROOT = $(shell pkg-config --variable=prefix fftw3)
  INCL     += $(FFTW_ROOT)/include 
  LDFLAGS  += -L$(FFTW_ROOT)/lib -lfftw3_mpi -lfftw3 -Wl,-rpath,$(FFTW_ROOT)/lib
  ## LAPACK/BLAS
  INCL 	   += $(MKLROOT)/include/intel64/lp64
  LDFLAGS  += -L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 \
		      -lmkl_sequential -lmkl_core -Wl,-rpath,$(MKLROOT)/lib/intel64
else
  $(error Unknown machine $(MACHINE))
endif

# ------------------------------------------------
# macro definitions
# ------------------------------------------------
#   _MPI                 # enables MPI code, used in the past for testing with serial version
#   with_assim           # enables 3D-var assimilation, otherwise only observation operator is enabled
#   with_lapack95_mkl    # enables calls to LAPack routines in covariance calculations (only needed when "with_assim" is defined)
#   with_ajs             # special tests by AJS (output per processor, timing)
# ------------------------------------------------
##~ full 3D-var application:
DFLAGS  = -D_MPI -Dwith_lapack95_mkl -Dwith_assim
#~ observation operator only:
#DFLAGS  = -D_MPI
