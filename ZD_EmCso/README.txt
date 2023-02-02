Simulate satellite observations from EMEP using CSO tools
=========================================================


Source files
------------

The following files come with the "CSO" observation operator code:

  cso*.F90

The interface with the EMEP model is implemented in:

  EMEP_CSO_mod.f90
  
The following files are part of the standard model code,
but have extra calls to the CSO interface;
the Makefile will compile these versions instead of standard codes:

  emep_Main.f90
  PhyChem_mod.f90


Compile
-------

To compile a model version including this extension, run from here:

  make emepctm

Alternatively, this could also be called from the main directory
if the change to "Makefile" is made as described below.


Makefiles
---------

To have these modules compiled with model:

- Changes to:  ../Makefile

  - Added target to let "make EmCso" create the executable:

      # simulation of satellite retrievals
      EmCso: EmChem19rp
       	$(MAKE) -C ZD_$@ $(PROG)

  - Ensure that .F90 files are compiled 
    (some older Makefiles might not have that):
  
      %.o: %.F90
	      $(F90) $(F90FLAGS) -c $< -o $@

- Changes to: Makefile.SRCS

  - Include local files:

      NEW_SRCS = cso_comm.F90 ... EMEP_CSO_mod.F90
              
      RENEW_SRCS = PhyChem_mod.f90 emep_Main.f90


Configuration files
-------------------

The "config/" directory contains examples of:
- EMEP configuration file "config_emep.nml" with extension for CSO settings;
- example "cso-settings.rc" with CSO specific settings.



Changes
-------

Support extra types for broadcast etc.
  cso_comm.F90

Read on root, scatter to domains.
  cso_ncfile.F90
  cso_domains.F90
  cso_pixels.F90
  cso_sat.F90

Merged updates from main model code.
  emep_Main.f90
  PhyChem_mod.f90

Deleted file that is not different from main code.
  GridValues_mod.f90
	Makefile.SRCS

Added `CSO_Format` routine.
	cso_datetime.F90

Added `CSO_CheckDir` routine.
	cso_file.F90

Extended error messages.
	cso_listing.F90
	cso_string.F90

Support input and output of packed variables.
	cso_mapping.F90
	cso_ncfile.F90
	cso_pixels.F90
	cso_sat.F90

Added `GetWeightsCell` methodes to assign points to grid cells.
	cso_grid.F90

Support integer(1) and character variables.
	cso_comm.F90
	cso_domains.F90
	cso_ncfile.F90
	cso_pixels.F90
	cso_sat.F90

Trap undefined values in kernel application.
Do not pack `time` variables.
	cso_pixels.F90

Added `LinInterp` routine to perform linear interpolation.
	cso_tools.F90

Added formula `TropoPressures` to create variable with half-level pressures
at surface, top of boundary layer, and top of troposhere.
	cso_pixels.F90
	EMEP_CSO_mod.f90

Added grid mapping for point locations, this is used if corner locations are not available.
	cso_sat.F90
	EMEP_CSO_mod.f90
  
Introduced ObsCompInfo type to accumulate tracers (for example for total PM) and apply unit conversions.
	EMEP_CSO_mod.f90

Fixed problem with creation of output directories when using Intel compiler.
  cso_file.F90
