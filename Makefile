#
#
PROG =	Unimod
###################################################

include Makefile.SRCS

###################################################

LIBS = -lnetcdf
INCL = -I/global/apps/netcdf/3.6.2/include
LLIB = -L/global/apps/netcdf/3.6.2/lib


F90 = mpif90

LDFLAGS =  -shared-intel -CB -r8  -recursive   -debug-parameters all -traceback  -ftrapuv -g -fpe0 -O0 -convert big_endian -IPF_fp_relaxed $(INCL)
LDFLAGS =  -shared-intel -r8  -recursive   -O2 -ftz -convert big_endian  -IPF_fp_relaxed 
F90FLAGS =  $(LDFLAGS) $(INCL)


LD = mpif90


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
	-l "$(F90) $(LDFLAGS) $(LLIB) -o $(PROG) $(FOBJ) $(INCL) $(LIBS)" > .depend

clean: diskclean touchdepend depend

diskclean:
	rm -f $(PROG) *.o *.mod

touchdepend:
	touch .depend
##########################################################

