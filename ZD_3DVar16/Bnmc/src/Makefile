#
# Makefile
#
#

# --- settings --------------------------------------

## define flags:
#with_omp = True
#with_mpi = True

# destination prefix:
prefix = ..

# target:
exe = emep_nmc.x

# which flags ?
flags = optim
#flags = check

# default machine:
machine = vilje

# machine specific settings; these define:
#   F2003, LD
#   FFLAGS, LDFLAGS
#   INCS, LIBS
include Makefile.$(machine)


# --- rules -------------------------------------------

%.o: %.F90
	$(F2003) -c -o $@ $(FFLAGS) $(INCS) $<
	@echo ' '


# --- help --------------------------------------------

help:
	@echo ''
	@echo 'Usage:'
	@echo '  make $(exe) \'
	@echo '        [machine=vilje] \'
	@echo '        [flags=optim|check] \'
	@echo '        [with_mpi=True] [with_omp=True]'
	@echo '  make install [prefix=$(prefix)]'
	@echo '  make clean'
	@echo ''

	
# --- targets ------------------------------------------

go_fu.o : go_fu.F90 go.inc
go_print.o : go_print.F90 go.inc go_fu.o
go_timer.o : go_timer.F90 go.inc go_print.o go_fu.o
go_string.o : go_string.F90 go.inc go_print.o
go_file.o : go_file.F90 go.inc go_print.o go_fu.o
go_rc.o : go_rc.F90 go.inc go_string.o go_print.o go_file.o
go_date.o : go_date.F90 go.inc go_string.o go_print.o
go_sampstat.o : go_sampstat.F90 go.inc go_print.o
go.o : go.F90 go.inc go_sampstat.o go_date.o go_rc.o go_file.o go_string.o go_timer.o go_print.o

GO_OBJS = go.o go_sampstat.o go_date.o go_rc.o go_file.o go_string.o go_timer.o go_print.o go_fu.o

# *

file_nc_i4.o: file_nc_i4.F90
file_nc_r4.o: file_nc_r4.F90
file_nc_r8.o: file_nc_r8.F90
file_nc.o: file_nc.F90 file_nc_i4.o file_nc_r4.o file_nc_r8.o

FILE_NC_OBJS = file_nc.o file_nc_i4.o file_nc_r4.o file_nc_r8.o

# *

c3po_datafile.o : c3po_datafile.F90 c3po.inc go.o
c3po_coordinates.o : c3po_coordinates.F90 c3po.inc c3po_datafile.o go.o
c3po.o : c3po.F90 c3po.inc c3po_coordinates.o

C3PO_OBJS = c3po.o c3po_coordinates.o c3po_datafile.o

# *

fftw3.o : fftw3.F90
ufftw3.o : ufftw3.F90 fftw3.o

FFTW3_OBJS = ufftw3.o fftw3.o

# *

linalg_lapack_350.o : linalg_lapack_350.F90
linalg.o : linalg.F90 linalg_lapack_350.o

LINALG_OBJS = linalg.o linalg_lapack_350.o

# *

emep_nmc_common.o : emep_nmc_common.F90 go.o file_nc.o
emep_nmc_modeloutput.o : emep_nmc_modeloutput.F90 go.o
emep_nmc_output.o : emep_nmc_output.F90 c3po.o go.o
emep_nmc_linalg.o : emep_nmc_linalg.F90 linalg.o go.o
emep_bcovarsqrt.o : emep_bcovarsqrt.F90 emep_nmc_output.o c3po.o go.o ufftw3.o fftw3.o
emep_nmc_driver_eps.o : emep_nmc_driver_eps.F90 emep_nmc_common.o emep_nmc_modeloutput.o emep_nmc_output.o c3po.o go.o
emep_nmc_driver_eps_eval.o : emep_nmc_driver_eps_eval.F90 emep_nmc_common.o emep_nmc_output.o c3po.o go.o
emep_nmc_driver_eta.o : emep_nmc_driver_eta.F90 emep_nmc_common.o emep_nmc_output.o c3po.o go.o
emep_nmc_driver_eta_f.o : emep_nmc_driver_eta_f.F90 emep_nmc_common.o emep_nmc_output.o c3po.o go.o ufftw3.o fftw3.o
emep_nmc_driver_c_f.o : emep_nmc_driver_c_f.F90 emep_nmc_common.o emep_nmc_output.o c3po.o go.o
emep_nmc_driver_d_f.o : emep_nmc_driver_d_f.F90 emep_nmc_output.o c3po.o go.o
emep_nmc_driver_b_f.o : emep_nmc_driver_b_f.F90 emep_nmc_output.o c3po.o go.o
emep_nmc_driver_xlx.o : emep_nmc_driver_xlx.F90 emep_nmc_linalg.o emep_nmc_output.o c3po.o go.o
emep_nmc_driver_bb.o : emep_nmc_driver_bb.F90 emep_bcovarsqrt.o emep_nmc_output.o c3po.o go.o
emep_nmc_driver_bb.o : emep_nmc_driver_bb_eval.F90 emep_bcovarsqrt.o emep_nmc_output.o c3po.o go.o
emep_nmc_driver.o : emep_nmc_driver.F90 \
                       emep_nmc_driver_eps.o emep_nmc_driver_eps_eval.o \
                       emep_nmc_driver_eta.o \
                       emep_nmc_driver_eta_f.o \
                       emep_nmc_driver_c_f.o \
                       emep_nmc_driver_d_f.o \
                       emep_nmc_driver_b_f.o \
                       emep_nmc_driver_xlx.o \
                       emep_nmc_driver_bb.o emep_nmc_driver_bb_eval.o \
                       go.o
emep_nmc.o : emep_nmc.F90 emep_bcovarsqrt.o emep_nmc_driver.o go.o

EMEP_NMC_OBJS = emep_nmc.o emep_bcovarsqrt.o \
                emep_nmc_driver.o \
                emep_nmc_driver_eps.o emep_nmc_modeloutput.o \
                emep_nmc_driver_eps_eval.o \
                emep_nmc_driver_eta.o \
                emep_nmc_driver_eta_f.o \
                emep_nmc_driver_c_f.o \
                emep_nmc_driver_d_f.o \
                emep_nmc_driver_b_f.o \
                emep_nmc_driver_xlx.o \
                emep_nmc_driver_bb.o emep_nmc_driver_bb_eval.o \
                emep_nmc_linalg.o \
                emep_nmc_output.o emep_nmc_common.o \
                $(FFTW3_OBJS) $(LINALG_OBJS) $(C3PO_OBJS) $(FILE_NC_OBJS) $(GO_OBJS)

$(exe): emep_nmc.o
	$(LD) -o $@ $(LDFLAGS) $(EMEP_NMC_OBJS) $(LIBS)

# *

install:
	@if [ ! -d $(prefix)/bin ]; then /bin/mkdir -p $(prefix)/bin ; fi
	/bin/cp $(exe) $(prefix)/bin

# *

clean:
	/bin/rm -f *.o *.mod *.x
	/bin/rm -f launcher*.out
	/bin/rm -f launcher*.err


# --- end ----------------------------------------------
