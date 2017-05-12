EMEP 3D var
===========

Version 17, intended for CAMS-50 November-2017 update.
Subversion "e", support total PM observations.
Subversion "f", test per-station obs.repr. errors and scale factor for sigma.

vilje:~arjos/work/emo/Unimod.GIT.3407.ajs1/ZD_3DVar17f/

 
mpif90 -D_MPI -c -o Unimod.o -extend-source -fpp -recursive -r8 -convert big_endian -shared-intel -O3 -ftz  -Dwith_lapack95_mkl -Dwith_ajs -I/sw/sdev/Modules/netcdf/netcdf-4.3.0/include -I/sw/sdev/Modules/fftw/fftw-3.3.3/include -I/sw/sdev/Modules/intelcomp/14.0.1/composer_xe_2013_sp1.1.106/mkl/include/intel64/lp64 Unimod.f90
 
mpif90 -D_MPI -o Unimod_ajs  *.o -L/sw/sdev/Modules/netcdf/netcdf-4.3.0/lib -lnetcdff -lnetcdf -Wl,-rpath -Wl,/sw/sdev/Modules/netcdf/netcdf-4.3.0/lib -L/sw/sdev/Modules/fftw/fftw-3.3.3/lib -lfftw3_mpi -lfftw3 -Wl,-rpath -Wl,/sw/sdev/Modules/fftw/fftw-3.3.3/lib -L/sw/sdev/Modules/intelcomp/14.0.1/composer_xe_2013_sp1.1.106/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath -Wl,/sw/sdev/Modules/intelcomp/14.0.1/composer_xe_2013_sp1.1.106/mkl/lib/intel64


