!
! Fortran module for FFTW,
! following instructions in documentation:
!
!   FFTW User Manual (v3.3.4)
!   7 Calling FFTW from Modern Fortran
!   7.7 Defining an FFTW module
!
! http://www.fftw.org/fftw3_doc/Defining-an-FFTW-module.html
!

module FFTW3

  use, intrinsic :: iso_c_binding

#ifdef _MPI
  include 'fftw3-mpi.f03'
#else
  include 'fftw3.f03'
#endif

end module
