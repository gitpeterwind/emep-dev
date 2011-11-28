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

F90FLAGS =  -shared-intel -CB -r8  -recursive   -debug-parameters all -traceback  -ftrapuv -g -fpe0 -O0 -convert big_endian -IPF_fp_relaxed -I$(INCL)
#F90FLAGS =  -shared-intel -r8  -recursive   -O2 -ftz -convert big_endian  -IPF_fp_relaxed -I$(INCL)
LDFLAGS =  $(F90FLAGS)  -L$(LLIB) -Wl,-rpath,$(LLIB) $(LIBS)




.SUFFIXES: $(SUFFIXES)  .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<


# Include the dependency-list created by makedepf90 below
all:  $(PROG)

include .depend

#
depend .depend:
	/home/mifapw/bin/makedepf90 $(SRCS) \
	-o $(PROG) \
	-l "$(F90)   -o $(PROG) $(FOBJ) $(LDFLAGS) " > .depend

clean: diskclean touchdepend depend

diskclean:
	rm -f $(PROG) *.o *.mod

touchdepend:
	touch .depend
##########################################################

