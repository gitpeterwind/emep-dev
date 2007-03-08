#
#
PROG =	Unimod
SRCS =	Aero_Rb_ml.f90 Aero_water_ml.f90 Ammonium_ml.f90 \
	Advection_ml.f90 AirEmis_ml.f90 Aqueous_ml.f90 Biogenics_ml.f90 \
	BoundaryConditions_ml.f90 Chem_ml.f90 Country_ml.f90 Dates_ml.f90 \
	DefPhotolysis_ml.f90 DepVariables_ml.f90 Derived_ml.f90 DryDep_ml.f90 \
	EmisDef_ml.f90 EmisGet_ml.f90 Emissions_ml.f90 Functions_ml.f90 \
	gc_com.f GlobalBCs_ml.f90 GridValues_ml.f90 Io_ml.f90 MassBudget_ml.f90 \
	Met_ml.f90 EQSAM_ml.f90 MARS_ml.f90 ModelConstants_ml.f90 My_Aerosols_ml.f90 \
	My_BoundConditions_ml.f90 My_Chem_ml.f90 My_Derived_ml.f90 \
	My_DryDep_ml.f90 My_Emis_ml.f90 My_MassBudget_ml.f90 \
	My_Outputs_ml.f90 My_WetDep_ml.f90 NetCDF_ml.f90 Nest_ml.f90 Out_restri_ml.f90 \
	Output_hourly.f90 Par_ml.f90 \
	PhysicalConstants_ml.f90 Trajectory_ml.f90 Radiation_ml.f90 \
	ReadField_ml.f90 Rsurface_ml.f90 Runchem_ml.f90 SOA_ml.f90 Setup_1d_ml.f90 \
	Setup_1dfields_ml.f90 Sites_ml.f90 SoilWater_ml.f90 Solver.f90 SeaSalt_ml.f90 \
	SubMet_ml.f90 Tabulations_ml.f90 TimeDate_ml.f90 Timefactors_ml.f90 Timing_ml.f90 \
	UK_ml.f90 UKsetup_ml.f90 Unimod.f90 Volcanos_ml.f90 Wesely_ml.f90 \
	Wrtchem.f90 global2local.f local2global.f PhyChem_ml.f90


F90_CONFORM_CHECK_ABORT=ON

###################################################

LIBS = -lmpi -lnetcdf
INCL = -I/home/u4/mifahik/netcdf/include
LLIB = -L/home/u4/mifahik/netcdf/lib64


F90 = f90

F90FLAGS = -64 -r8 -O3 -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV $(INCL)
#F90FLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV $(INCL)
#F90FLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -DEBUG:conform_check=ON -DEBUG:subscript_check:verbose_runtime=ON -DEBUG:fullwarn=ON -TARG:exc_min=0ZV $(INCL)

LDFLAGS = -64 -r8 -O 3 -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV
#LDFLAGS = -64 -r8 -C -g -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV
#LDFLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -DEBUG:conform_check=ON -DEBUG:subscript_check:verbose_runtime=ON -DEBUG:fullwarn=ON -TARG:exc_min=0ZV
LD = f90


.SUFFIXES: $(SUFFIXES)  .f .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f.o:	
	$(F90) $(F90FLAGS) -c $<


# Include the dependency-list created by makedepf90 below
all:  $(PROG)

include .depend

#LLIB added, ds
depend .depend:
	/home/u4/mifahik/bin/makedepf90 $(SRCS) \
	-o $(PROG) \
	-l "$$(F90) $$(LDFLAGS) $$(LLIB) -o $$(PROG) $$(FOBJ) $$(INCL) $$(LIBS)" > .depend

clean:
	rm -f $(PROG) *.o *.mod .depend; \
	#touch .depend
#make depend


##########################################################

