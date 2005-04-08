#
# Problems with makemake.perl
#     My_Runmode_ml not included for GridValues 
#     gc_com stuff not needed
#     repetition of .F files in SRCS; so it gave two global2local.f
#     etc.
#     My_Chem_ml.o was included in My_Chem_ml.o definition - gives
#     predecessor cycle.
#
# Currently fix by using Makefile.insert2 and Makefile.tail
#
PROG =	Unimod
SRCS =	Aero_Rb_ml.f90 Aero_water_ml.f90 Ammonium_ml.f90 \
	Advection_ml.f90 AirEmis_ml.f90 Aqueous_ml.f90 Biogenics_ml.f90 \
	BoundaryConditions_ml.f90 Chem_ml.f90 Country_ml.f90 Dates_ml.f90 \
	DefPhotolysis_ml.f90 DepVariables_ml.f90 Derived_ml.f90 DryDep_ml.f90 \
	EmisDef_ml.f90 EmisGet_ml.f90 Emissions_ml.f90 Functions_ml.f90 \
	GlobalBCs_ml.f90 GridValues_ml.f90 Io_ml.f90 MassBudget_ml.f90 \
	Met_ml.f90 EQSAM_ml.f90 MARS_ml.f90 ModelConstants_ml.f90 My_Aerosols_ml.f90 \
	My_BoundConditions_ml.f90 My_Chem_ml.f90 My_Derived_ml.f90 \
	My_DryDep_ml.f90 My_Emis_ml.f90 My_MassBudget_ml.f90 \
	My_Outputs_ml.f90 My_WetDep_ml.f90 NetCDF_ml.f90 Nest_ml.f90 Out_restri_ml.f90 \
	Output_binary.f90 Output_hourly.f90 Par_ml.f90 \
	PhysicalConstants_ml.f90 Polinat_ml.f90 Radiation_ml.f90 \
	ReadField_ml.f90 Rsurface_ml.f90 Runchem_ml.f90 SOA_ml.f90 Setup_1d_ml.f90 \
	Setup_1dfields_ml.f90 Sites_ml.f90 SoilWater_ml.f90 Solver.f90 SeaSalt_ml.f90 \
	SubMet_ml.f90 Tabulations_ml.f90 TimeDate_ml.f90 Timefactors_ml.f90 Timing_ml.f90 \
	UK_ml.f90 UKsetup_ml.f90 Unimod.f90 Volcanos_ml.f90 Wesely_ml.f90 \
	Wrtchem.f90 global2local.f local2global.f outchem_restri.f phyche.f


EXTRA_OBJS = put_restri_i2.o putflti2.o rmycop.o \
	gc_com.o getflti2.o

F90_CONFORM_CHECK_ABORT=ON

###################################################
# Diff GRIDUR
ARCH = PGON
#_CRAY_CF77LIBS =
CF77LIBS = -lmpi
#LIBS = $(CF77LIBS)	
LIBS = -lmpi -lnetcdf
INCL = -I/home/u4/mifahik/netcdf/include
LLIB = -L/home/u4/mifahik/netcdf/lib64

CPP = cpp
#_CRAY_CPPFLAGS	= -P -N -DMPI_SRC 
CPPFLAGS	= -P -cpp -DMPI_SRC -DFLP_64B
CC = cc
#_CRAY_CFLAGS = -O3
#CFLAGS = -64 -C
CFLAGS = -64 -O 3
FC = f90
#_CRAY_FFLAGS = -O 3,fusion,aggress,bl,msgs,negmsgs,unroll2
#FFLAGS = -default64 -O3
#FFLAGS = -64 -r8 -O 3
FFLAGS = -64 -r8 -O3 -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV
#FFLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -DEBUG:conform_check=ON -DEBUG:subscript_check:verbose_runtime=ON -DEBUG:fullwarn=ON -TARG:exc_min=0ZV

F90 = f90

F90FLAGS = -64 -r8 -O3 -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV $(INCL)
#F90FLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV $(INCL)
#F90FLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -DEBUG:conform_check=ON -DEBUG:subscript_check:verbose_runtime=ON -DEBUG:fullwarn=ON -TARG:exc_min=0ZV $(INCL)
#_CRAY_F90FLAGS = -O 3,fusion,aggress,bl,unroll2,msgs,negmsgs

#_CRAY_LDFLAGS = -X8 -O 3,fusion,aggress,bl,msgs,negmsgs,unroll2 
#LDFLAGS = -default64 -O3
#QUERY LDFLAGS = -64 -r8 -O 3 -g -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV
LDFLAGS = -64 -r8 -O 3 -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV
#LDFLAGS = -64 -r8 -C -g -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV
#LDFLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -DEBUG:conform_check=ON -DEBUG:subscript_check:verbose_runtime=ON -DEBUG:fullwarn=ON -TARG:exc_min=0ZV
LD = f90

#_CRAY_$(PROG): $(OBJS)
#_CRAY_	$(F90) $(LDFLAGS) -o $@ $(OBJS)) $(LIBS)

.SUFFIXES: $(SUFFIXES) .c .f .f90 .F

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.F.o:	
	-rm -f $*.f $*.i
	$(FC) $(CPPFLAGS) $*.F
	mv $*.i $*.f
	$(F90) -c $(F90FLAGS) $*.f
#		$(CPP) $(CPPFLAGS) $*.F $*.f

.f.o:	
	$(F90) $(F90FLAGS) -c $<

.c.o:	
	$(CC) -c $(CFLAGS) $<

# Include the dependency-list created by makedepf90 below
all: $(EXTRA_OBJS) $(PROG)

include .depend

#LLIB added, ds
depend .depend:
	/home/u4/mifahik/bin/makedepf90 $(SRCS) \
	-o $(PROG) \
	-l "$$(F90) $$(LDFLAGS) $$(LLIB) -o $$(PROG) $$(FOBJ) $$(EXTRA_OBJS) $$(INCL) $$(LIBS)" > .depend

clean:
	rm -f $(PROG) *.o *.mod .depend; \
	#touch .depend
#make depend


##########################################################
#Makefile.tail  copied from Makefile of 18/10/01

gc_com.o : gc_com.F
	-rm -f gc_com.f gc_com.i
	$(FC) $(CPPFLAGS) gc_com.F
	mv gc_com.i gc_com.f
	$(F90) $(F90FLAGS) -c gc_com.f
#gc_com.o : gc_com.F
#        -rm -f gc_com.f gc_com.i
#        $(CPP) $(CPPFLAGS) gc_com.F gc_com.f
#        $(F90) $(F90FLAGS) -c gc_com.f
getflti2.o : Par_ml.o getflti2.F
	-rm -f getflti2.f getflti2.i
	$(FC) $(CPPFLAGS) getflti2.F
	mv getflti2.i getflti2.f
	$(F90) $(F90FLAGS) -c getflti2.f
#getflti2.o : Par_ml.o getflti2.F
#	-rm -f getflti2.f getflti2.i
#	$(CPP) $(CPPFLAGS) getflti2.F getflti2.f
#	$(F90) $(F90FLAGS) -c Par_ml.o getflti2.f


putflti2.o : Par_ml.o putflti2.F
	-rm -f  putflti2.f  putflti2.i
	$(FC) $(CPPFLAGS)  putflti2.F
	mv  putflti2.i  putflti2.f
	$(F90) $(F90FLAGS) -c  putflti2.f
#putflti2.o : Par_ml.o putflti2.F
#        -rm -f putflti2.f putflti2.i
#        $(CPP) $(CPPFLAGS) putflti2.F putflti2.f
#        $(F90) $(F90FLAGS) -c Par_ml.o putflti2.f


put_restri_i2.o : Par_ml.o Out_restri_ml.o put_restri_i2.F
	-rm -f put_restri_i2.f put_restri_i2.i
	$(FC) $(CPPFLAGS) put_restri_i2.F
	mv put_restri_i2.i put_restri_i2.f
	$(F90) $(F90FLAGS) -c put_restri_i2.f
#put_restri_i2.o : Par_ml.o Out_restri_ml.o put_restri_i2.F
#	-rm -f put_restri_i2.f put_restri_i2.i
#	$(CPP) $(CPPFLAGS) put_restri_i2.F put_restri_i2.f
#	$(F90) $(F90FLAGS) -c Par_ml.o Out_restri_ml.o put_restri_i2.f
