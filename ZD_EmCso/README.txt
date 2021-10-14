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

