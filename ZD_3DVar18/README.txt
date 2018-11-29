EMEP 3D var
===========

Version 17, intended for CAMS-50 November-2017 update.
Subversion "e", support total PM observations.
Subversion "f", test per-station obs.repr. errors and scale factor for sigma.

Changes
-------

Explicitly set sigma to zero at (top) boundary cells to avoid strange updates.
(AJS, 2018-10)
  B_NMC/emep_bcovarsqrt.F90

Commented 'use MPI' statements for 'MPI_AllGather' and 'MPI_AllToAllV',
these seem not present in the MPI module on Elvis but will be linked correctly afterwards.
  tools/go_par.F90
  tools/mpif90.F90
