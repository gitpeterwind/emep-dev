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

SRCS =	Advection_ml.f90 AirEmis_ml.f90 Aqueous_ml.f90 Biogenics_ml.f90 \
	BoundaryConditions_ml.f90 Chem_ml.f90 Country_ml.f90 Dates_ml.f90 \
	DefPhotolysis_ml.f90 DepVariables_ml.f90 Derived_ml.f90 DryDep_ml.f90 \
	EmisDef_ml.f90 EmisGet_ml.f90 Emissions_ml.f90 Functions_ml.f90 \
	GlobalBCs_ml.f90 GridValues_ml.f90 Io_ml.f90 MassBudget_ml.f90 \
	Met_ml.f90 ModelConstants_ml.f90 My_Aerosols_ml.f90 \
	My_BoundConditions_ml.f90 My_Chem_ml.f90 My_Derived_ml.f90 \
	My_DryDep_ml.f90 My_Emis_ml.f90 My_MassBudget_ml.f90 \
	My_Outputs_ml.f90 My_WetDep_ml.f90 Nest_ml.f90 Out_restri_ml.f90 \
	Output_binary.f90 Output_hourly.f90 Par_ml.f90 \
	PhysicalConstants_ml.f90 Polinat_ml.f90 Radiation_ml.f90 \
	ReadField_ml.f90 Rsurface_ml.f90 Runchem_ml.f90 Setup_1d_ml.f90 \
	Setup_1dfields_ml.f90 Sites_ml.f90 SoilWater_ml.f90 Solver.f90 \
	SubMet_ml.f90 Tabulations_ml.f90 Timefactors_ml.f90 Timing_ml.f90 \
	UK_ml.f90 UKsetup_ml.f90 Unimod.f90 Volcanos_ml.f90 Wesely_ml.f90 \
	Wrtchem.f90 \
	outchem_restri.f phyche.f rmycop.c \
	global2local.f local2global.f gc_com.F \
	getflti2.F put_restri_i2.F putflti2.F

OBJS =	Advection_ml.o AirEmis_ml.o Aqueous_ml.o Biogenics_ml.o \
	BoundaryConditions_ml.o Chem_ml.o Country_ml.o Dates_ml.o \
	DefPhotolysis_ml.o DepVariables_ml.o Derived_ml.o DryDep_ml.o \
	EmisDef_ml.o EmisGet_ml.o Emissions_ml.o Functions_ml.o \
	GlobalBCs_ml.o GridValues_ml.o Io_ml.o MassBudget_ml.o Met_ml.o \
	ModelConstants_ml.o My_Aerosols_ml.o My_BoundConditions_ml.o \
	My_Chem_ml.o My_Derived_ml.o My_DryDep_ml.o My_Emis_ml.o \
	My_MassBudget_ml.o My_Outputs_ml.o My_WetDep_ml.o Nest_ml.o \
	Out_restri_ml.o Output_binary.o Output_hourly.o Par_ml.o \
	PhysicalConstants_ml.o Polinat_ml.o Radiation_ml.o ReadField_ml.o \
	Rsurface_ml.o Runchem_ml.o Setup_1d_ml.o Setup_1dfields_ml.o \
	Sites_ml.o SoilWater_ml.o Solver.o SubMet_ml.o Tabulations_ml.o \
	Timefactors_ml.o Timing_ml.o UK_ml.o UKsetup_ml.o Unimod.o \
	Volcanos_ml.o Wesely_ml.o Wrtchem.o \
	global2local.o local2global.o outchem_restri.o phyche.o \
	put_restri_i2.o putflti2.o rmycop.o \
	gc_com.o getflti2.o

###################################################
# Diff GRIDUR
ARCH = PGON
#_CRAY_CF77LIBS =
CF77LIBS = -lmpi
LIBS = $(CF77LIBS)	
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
FFLAGS = -64 -r8 -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV
FFLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV

F90 = f90
#F90FLAGS = -64 -r8 -O3 -fullwarn
F90FLAGS = -64 -r8 -O3 -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV
F90FLAGS = -64 -r8 -g -C -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV
#_CRAY_F90FLAGS = -O 3,fusion,aggress,bl,unroll2,msgs,negmsgs

#_CRAY_LDFLAGS = -X8 -O 3,fusion,aggress,bl,msgs,negmsgs,unroll2 
#LDFLAGS = -default64 -O3
LDFLAGS = -64 -r8 -O 3 -g -OPT:IEEE_arithm=3:roundoff=3 -TARG:exc_min=0ZV
LDFLAGS = -64 -r8 -C -g -DEBUG:trap_uninitialized=ON:verbose_runtime=ON -TARG:exc_min=0ZV
LD = f90

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)
#_CRAY_$(PROG): $(OBJS)
#_CRAY_	$(F90) $(LDFLAGS) -o $@ $(OBJS)) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

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

###################################################
Advection_ml.o: Chem_ml.o GridValues_ml.o MassBudget_ml.o Met_ml.o \
	ModelConstants_ml.o My_Chem_ml.o Par_ml.o Timing_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o GridValues_ml.o \
		MassBudget_ml.o Met_ml.o ModelConstants_ml.o \
		My_Chem_ml.o Par_ml.o Timing_ml.o Advection_ml.f90
AirEmis_ml.o: GridValues_ml.o Io_ml.o ModelConstants_ml.o Par_ml.o \
	PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c GridValues_ml.o Io_ml.o \
		ModelConstants_ml.o Par_ml.o PhysicalConstants_ml.o \
		AirEmis_ml.f90
Aqueous_ml.o: GridValues_ml.o Met_ml.o ModelConstants_ml.o My_Derived_ml.o \
	My_WetDep_ml.o Par_ml.o Setup_1dfields_ml.o
	$(F90) $(F90FLAGS) -c GridValues_ml.o Met_ml.o \
		ModelConstants_ml.o My_Derived_ml.o My_WetDep_ml.o \
		Par_ml.o Setup_1dfields_ml.o Aqueous_ml.f90
Biogenics_ml.o: GridValues_ml.o Io_ml.o My_Emis_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c GridValues_ml.o Io_ml.o My_Emis_ml.o \
		Par_ml.o Biogenics_ml.f90
BoundaryConditions_ml.o: Chem_ml.o GlobalBCs_ml.o GridValues_ml.o \
	ModelConstants_ml.o My_BoundConditions_ml.o My_Chem_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o GlobalBCs_ml.o \
		GridValues_ml.o ModelConstants_ml.o \
		My_BoundConditions_ml.o My_Chem_ml.o Par_ml.o \
		BoundaryConditions_ml.f90
Chem_ml.o: ModelConstants_ml.o My_Chem_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c ModelConstants_ml.o My_Chem_ml.o \
		Par_ml.o Chem_ml.f90
DefPhotolysis_ml.o: GridValues_ml.o Io_ml.o Met_ml.o ModelConstants_ml.o \
	Par_ml.o Setup_1dfields_ml.o
	$(F90) $(F90FLAGS) -c GridValues_ml.o Io_ml.o Met_ml.o \
		ModelConstants_ml.o Par_ml.o Setup_1dfields_ml.o \
		DefPhotolysis_ml.f90
Derived_ml.o: Chem_ml.o Met_ml.o ModelConstants_ml.o My_Chem_ml.o \
	My_Derived_ml.o My_Outputs_ml.o Par_ml.o PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o Met_ml.o ModelConstants_ml.o \
		My_Chem_ml.o My_Derived_ml.o My_Outputs_ml.o Par_ml.o \
		PhysicalConstants_ml.o Derived_ml.f90
DryDep_ml.o: Chem_ml.o Dates_ml.o DepVariables_ml.o GridValues_ml.o \
	MassBudget_ml.o Met_ml.o ModelConstants_ml.o My_Chem_ml.o \
	My_DryDep_ml.o My_Derived_ml.o \
	Par_ml.o PhysicalConstants_ml.o Radiation_ml.o \
	Rsurface_ml.o SoilWater_ml.o SubMet_ml.o UK_ml.o Wesely_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o Dates_ml.o DepVariables_ml.o \
		GridValues_ml.o MassBudget_ml.o Met_ml.o \
		ModelConstants_ml.o My_Chem_ml.o My_Derived_ml.o \
		My_DryDep_ml.o \
		Par_ml.o PhysicalConstants_ml.o Radiation_ml.o \
		Rsurface_ml.o SoilWater_ml.o SubMet_ml.o UK_ml.o \
		Wesely_ml.o DryDep_ml.f90
EmisGet_ml.o: Country_ml.o EmisDef_ml.o Functions_ml.o Io_ml.o My_Emis_ml.o \
	Par_ml.o Volcanos_ml.o
	$(F90) $(F90FLAGS) -c Country_ml.o EmisDef_ml.o \
		Functions_ml.o Io_ml.o My_Emis_ml.o Par_ml.o \
		Volcanos_ml.o EmisGet_ml.f90
Emissions_ml.o: Biogenics_ml.o Country_ml.o Dates_ml.o EmisDef_ml.o \
	EmisGet_ml.o GridValues_ml.o Io_ml.o Met_ml.o ModelConstants_ml.o \
	My_Emis_ml.o My_MassBudget_ml.o Par_ml.o PhysicalConstants_ml.o \
	ReadField_ml.o Timefactors_ml.o Volcanos_ml.o
	$(F90) $(F90FLAGS) -c Biogenics_ml.o Country_ml.o Dates_ml.o \
		EmisDef_ml.o EmisGet_ml.o GridValues_ml.o Io_ml.o \
		Met_ml.o ModelConstants_ml.o My_Emis_ml.o \
		My_MassBudget_ml.o Par_ml.o PhysicalConstants_ml.o \
		ReadField_ml.o Timefactors_ml.o Volcanos_ml.o \
		Emissions_ml.f90
GlobalBCs_ml.o: Dates_ml.o GridValues_ml.o Io_ml.o ModelConstants_ml.o \
	Par_ml.o
	$(F90) $(F90FLAGS) -c Dates_ml.o GridValues_ml.o Io_ml.o \
		ModelConstants_ml.o Par_ml.o GlobalBCs_ml.f90
GridValues_ml.o: Io_ml.o ModelConstants_ml.o Par_ml.o PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c Io_ml.o ModelConstants_ml.o Par_ml.o \
		PhysicalConstants_ml.o GridValues_ml.f90
MassBudget_ml.o: Chem_ml.o GridValues_ml.o Io_ml.o Met_ml.o \
	ModelConstants_ml.o My_Chem_ml.o My_Derived_ml.o My_MassBudget_ml.o \
	Par_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o GridValues_ml.o Io_ml.o \
		Met_ml.o ModelConstants_ml.o My_Chem_ml.o \
		My_Derived_ml.o My_MassBudget_ml.o Par_ml.o \
		MassBudget_ml.f90
Met_ml.o: Dates_ml.o GridValues_ml.o Io_ml.o ModelConstants_ml.o Par_ml.o \
	PhysicalConstants_ml.o ReadField_ml.o Tabulations_ml.o
	$(F90) $(F90FLAGS) -c Dates_ml.o GridValues_ml.o Io_ml.o \
		ModelConstants_ml.o Par_ml.o PhysicalConstants_ml.o \
		ReadField_ml.o Tabulations_ml.o Met_ml.f90
ModelConstants_ml.o: Dates_ml.o PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c Dates_ml.o PhysicalConstants_ml.o \
		ModelConstants_ml.f90
My_Aerosols_ml.o: ModelConstants_ml.o My_Chem_ml.o PhysicalConstants_ml.o \
	Setup_1dfields_ml.o
	$(F90) $(F90FLAGS) -c ModelConstants_ml.o My_Chem_ml.o \
		PhysicalConstants_ml.o Setup_1dfields_ml.o \
		My_Aerosols_ml.f90
My_BoundConditions_ml.o: GlobalBCs_ml.o GridValues_ml.o Met_ml.o \
	ModelConstants_ml.o My_Chem_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c GlobalBCs_ml.o GridValues_ml.o \
		Met_ml.o ModelConstants_ml.o My_Chem_ml.o Par_ml.o \
		My_BoundConditions_ml.f90
My_Chem_ml.o: Dates_ml.o Functions_ml.o ModelConstants_ml.o \
	PhysicalConstants_ml.o Radiation_ml.o
	$(F90) $(F90FLAGS) -c Dates_ml.o Functions_ml.o \
		ModelConstants_ml.o PhysicalConstants_ml.o \
		Radiation_ml.o My_Chem_ml.f90
My_Derived_ml.o: Chem_ml.o Met_ml.o ModelConstants_ml.o My_Chem_ml.o Par_ml.o \
	PhysicalConstants_ml.o Radiation_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o Met_ml.o ModelConstants_ml.o \
		My_Chem_ml.o Par_ml.o PhysicalConstants_ml.o \
		Radiation_ml.o My_Derived_ml.f90
My_DryDep_ml.o: ModelConstants_ml.o My_Chem_ml.o My_Derived_ml.o Wesely_ml.o
	$(F90) $(F90FLAGS) -c ModelConstants_ml.o My_Chem_ml.o \
		My_Derived_ml.o Wesely_ml.o My_DryDep_ml.f90
My_MassBudget_ml.o: My_Chem_ml.o My_Emis_ml.o
	$(F90) $(F90FLAGS) -c My_Chem_ml.o My_Emis_ml.o \
		My_MassBudget_ml.f90
My_Outputs_ml.o: Dates_ml.o ModelConstants_ml.o My_Chem_ml.o My_Derived_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c Dates_ml.o ModelConstants_ml.o \
		My_Chem_ml.o My_Derived_ml.o Par_ml.o My_Outputs_ml.f90
My_WetDep_ml.o: MassBudget_ml.o ModelConstants_ml.o My_Chem_ml.o \
	My_Derived_ml.o
	$(F90) $(F90FLAGS) -c MassBudget_ml.o ModelConstants_ml.o \
		My_Chem_ml.o My_Derived_ml.o My_WetDep_ml.f90
Nest_ml.o: Chem_ml.o Dates_ml.o Functions_ml.o GridValues_ml.o Io_ml.o \
	Met_ml.o ModelConstants_ml.o My_Chem_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o Dates_ml.o Functions_ml.o \
		GridValues_ml.o Io_ml.o Met_ml.o ModelConstants_ml.o \
		My_Chem_ml.o Par_ml.o Nest_ml.f90
Out_restri_ml.o: My_Outputs_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c My_Outputs_ml.o Par_ml.o \
		Out_restri_ml.f90
Output_binary.o: Derived_ml.o GridValues_ml.o Io_ml.o ModelConstants_ml.o \
	My_Derived_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c Derived_ml.o GridValues_ml.o Io_ml.o \
		ModelConstants_ml.o My_Derived_ml.o Par_ml.o \
		Output_binary.f90
Output_hourly.o: Chem_ml.o GridValues_ml.o Io_ml.o Met_ml.o \
	ModelConstants_ml.o My_Chem_ml.o My_Outputs_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o GridValues_ml.o Io_ml.o \
		Met_ml.o ModelConstants_ml.o My_Chem_ml.o \
		My_Outputs_ml.o Par_ml.o Output_hourly.f90
Polinat_ml.o: Chem_ml.o GridValues_ml.o Io_ml.o Met_ml.o ModelConstants_ml.o \
	My_Chem_ml.o My_Outputs_ml.o Par_ml.o
	$(F90) $(F90FLAGS) -c Chem_ml.o GridValues_ml.o Io_ml.o \
		Met_ml.o ModelConstants_ml.o My_Chem_ml.o \
		My_Outputs_ml.o Par_ml.o Polinat_ml.f90
Radiation_ml.o: Dates_ml.o GridValues_ml.o ModelConstants_ml.o Par_ml.o \
	PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c Dates_ml.o GridValues_ml.o \
		ModelConstants_ml.o Par_ml.o PhysicalConstants_ml.o \
		Radiation_ml.f90
ReadField_ml.o: Io_ml.o Par_ml.o GridValues_ml.o
	$(F90) $(F90FLAGS) -c Io_ml.o Par_ml.o GridValues_ml.o ReadField_ml.f90
Rsurface_ml.o: Dates_ml.o DepVariables_ml.o Functions_ml.o Io_ml.o \
	My_DryDep_ml.o PhysicalConstants_ml.o Radiation_ml.o SoilWater_ml.o \
	UKsetup_ml.o Wesely_ml.o
	$(F90) $(F90FLAGS) -c Dates_ml.o DepVariables_ml.o \
		Functions_ml.o Io_ml.o My_DryDep_ml.o \
		PhysicalConstants_ml.o Radiation_ml.o SoilWater_ml.o \
		UKsetup_ml.o Wesely_ml.o Rsurface_ml.f90
Runchem_ml.o: Aqueous_ml.o Chem_ml.o DefPhotolysis_ml.o ModelConstants_ml.o \
	My_Aerosols_ml.o My_Chem_ml.o Par_ml.o Setup_1d_ml.o \
	Setup_1dfields_ml.o Solver.o Timing_ml.o
	$(F90) $(F90FLAGS) -c Aqueous_ml.o Chem_ml.o \
		DefPhotolysis_ml.o ModelConstants_ml.o My_Aerosols_ml.o \
		My_Chem_ml.o Par_ml.o Setup_1d_ml.o \
		Setup_1dfields_ml.o Solver.o Timing_ml.o Runchem_ml.f90
Setup_1d_ml.o: AirEmis_ml.o Biogenics_ml.o Chem_ml.o Dates_ml.o \
	Emissions_ml.o GridValues_ml.o MassBudget_ml.o Met_ml.o \
	ModelConstants_ml.o My_BoundConditions_ml.o My_Chem_ml.o My_Emis_ml.o \
	My_MassBudget_ml.o Par_ml.o PhysicalConstants_ml.o Radiation_ml.o \
	Setup_1dfields_ml.o Tabulations_ml.o Volcanos_ml.o
	$(F90) $(F90FLAGS) -c AirEmis_ml.o Biogenics_ml.o Chem_ml.o \
		Dates_ml.o Emissions_ml.o GridValues_ml.o \
		MassBudget_ml.o Met_ml.o ModelConstants_ml.o \
		My_BoundConditions_ml.o My_Chem_ml.o My_Emis_ml.o \
		My_MassBudget_ml.o Par_ml.o PhysicalConstants_ml.o \
		Radiation_ml.o Setup_1dfields_ml.o Tabulations_ml.o \
		Volcanos_ml.o Setup_1d_ml.f90
Setup_1dfields_ml.o: Biogenics_ml.o ModelConstants_ml.o My_Chem_ml.o \
	My_Emis_ml.o
	$(F90) $(F90FLAGS) -c Biogenics_ml.o ModelConstants_ml.o \
		My_Chem_ml.o My_Emis_ml.o Setup_1dfields_ml.f90
Sites_ml.o: Io_ml.o Met_ml.o ModelConstants_ml.o My_Chem_ml.o My_Derived_ml.o \
	My_Outputs_ml.o \
	Par_ml.o
	$(F90) $(F90FLAGS) -c Io_ml.o Met_ml.o ModelConstants_ml.o \
		My_Chem_ml.o My_Derived_ml.o \
		My_Outputs_ml.o Par_ml.o Sites_ml.f90
Solver.o: Biogenics_ml.o DefPhotolysis_ml.o ModelConstants_ml.o \
	My_Aerosols_ml.o My_Chem_ml.o My_Emis_ml.o Par_ml.o Radiation_ml.o \
	Setup_1dfields_ml.o My_Reactions.inc
	$(F90) $(F90FLAGS) -c Biogenics_ml.o DefPhotolysis_ml.o \
		ModelConstants_ml.o My_Aerosols_ml.o My_Chem_ml.o \
		My_Emis_ml.o Par_ml.o Radiation_ml.o \
		Setup_1dfields_ml.o Solver.f90
SubMet_ml.o: Functions_ml.o Io_ml.o ModelConstants_ml.o \
	PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c Functions_ml.o Io_ml.o \
		ModelConstants_ml.o PhysicalConstants_ml.o SubMet_ml.f90
Tabulations_ml.o: ModelConstants_ml.o PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c ModelConstants_ml.o \
		PhysicalConstants_ml.o Tabulations_ml.f90
Timefactors_ml.o: Country_ml.o Dates_ml.o EmisDef_ml.o Io_ml.o My_Emis_ml.o
	$(F90) $(F90FLAGS) -c Country_ml.o Dates_ml.o EmisDef_ml.o \
		Io_ml.o My_Emis_ml.o Timefactors_ml.f90
UK_ml.o: Dates_ml.o DepVariables_ml.o Functions_ml.o GridValues_ml.o Io_ml.o \
	ModelConstants_ml.o Par_ml.o UKsetup_ml.o
	$(F90) $(F90FLAGS) -c Dates_ml.o DepVariables_ml.o \
		Functions_ml.o GridValues_ml.o Io_ml.o ModelConstants_ml.o Par_ml.o \
		UKsetup_ml.o UK_ml.f90
UKsetup_ml.o: DepVariables_ml.o Io_ml.o Met_ml.o
	$(F90) $(F90FLAGS) -c DepVariables_ml.o Io_ml.o Met_ml.o \
		UKsetup_ml.f90
Unimod.o: Advection_ml.o AirEmis_ml.o Biogenics_ml.o BoundaryConditions_ml.o \
	Chem_ml.o Dates_ml.o DefPhotolysis_ml.o Emissions_ml.o \
	GridValues_ml.o Io_ml.o MassBudget_ml.o Met_ml.o ModelConstants_ml.o \
	My_Chem_ml.o My_Emis_ml.o My_Outputs_ml.o My_WetDep_ml.o \
	Out_restri_ml.o Par_ml.o Polinat_ml.o Sites_ml.o Tabulations_ml.o \
	Timing_ml.o
	$(F90) $(F90FLAGS) -c Advection_ml.o AirEmis_ml.o \
		Biogenics_ml.o BoundaryConditions_ml.o Chem_ml.o \
		Dates_ml.o DefPhotolysis_ml.o Emissions_ml.o \
		GridValues_ml.o Io_ml.o MassBudget_ml.o Met_ml.o \
		ModelConstants_ml.o My_Chem_ml.o My_Emis_ml.o \
		My_Outputs_ml.o My_WetDep_ml.o Out_restri_ml.o \
		Par_ml.o Polinat_ml.o Sites_ml.o Tabulations_ml.o \
		Timing_ml.o Unimod.f90
Volcanos_ml.o: EmisDef_ml.o GridValues_ml.o Io_ml.o Met_ml.o \
	ModelConstants_ml.o My_Emis_ml.o Par_ml.o PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c EmisDef_ml.o GridValues_ml.o Io_ml.o \
		Met_ml.o ModelConstants_ml.o My_Emis_ml.o \
		Par_ml.o PhysicalConstants_ml.o Volcanos_ml.f90
Wesely_ml.o: PhysicalConstants_ml.o
	$(F90) $(F90FLAGS) -c PhysicalConstants_ml.o Wesely_ml.f90
Wrtchem.o: Derived_ml.o Io_ml.o ModelConstants_ml.o My_Derived_ml.o \
	My_Outputs_ml.o Out_restri_ml.o Output_binary.o Par_ml.o
	$(F90) $(F90FLAGS) -c Derived_ml.o Io_ml.o \
		ModelConstants_ml.o My_Derived_ml.o My_Outputs_ml.o \
		Out_restri_ml.o Output_binary.o Par_ml.o Wrtchem.f90
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

phyche.o :  My_Outputs_ml.o Sites_ml.o Timing_ml.o Dates_ml.o \
	DryDep_ml.o Par_ml.o My_Chem_ml.o Met_ml.o ModelConstants_ml.o \
	Emissions_ml.o Timefactors_ml.o Chem_ml.o Derived_ml.o Radiation_ml.o \
	Runchem_ml.o Advection_ml.o Polinat_ml.o  Nest_ml.o phyche.f
