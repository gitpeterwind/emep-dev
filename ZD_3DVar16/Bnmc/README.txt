EMEB Background covariance using NMC method
===========================================

Description
-----------

Scripts and code to create a NetCDF file with all entities
that are needed for operations with the background covariance
matrix.

Main script:
  launcher
  
This uses the settings in:
  rc/emep-nc.rc
  

History
-------

2015-08, Arjo Segers
  Asigned version number 'v1.0' to existing code.
  Development continues in 'vtrunk'.

2015-11, Arjo Segers
  Splitted driver module into sub modules.
  Support different tracer order between model and store Bsqrt files.
  Imported changes made to covariance sqrt module in EMEP assimilation code.

2016-10, Arjo Segers
  Added linalg modules as wrapper around various LAPack libraries to
  support porting to other computers.

