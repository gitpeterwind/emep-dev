!######################################################################
!
! LinAlg - interface to various flavours of LAPack routines
!
! DESCRIPTION
!
!   Prefered interface is MKL LAPack95 style as part of 
!   the Intel compiler suite.
!   Note that that is an upgrade of the netlib LAPack95 interface
!   which seems valid for older LAPack versions (3.0 ?) only.
!   This module therefore provides wrappers following the latest
!   interface style around LAPack (F77) or other implementation.
!
! MACRO'S
!
!  The following macro's are used to switch between various
!  external libraries:
!
!    with_lapack95_mkl : Lapack95 interface from Intel MKL
!
!    with_lapack       : native lapack modules ;
!                        here interface towards version 3.5.0
!
!######################################################################

module LinAlg

  !
  ! select appropriatue module:
  !
#ifdef with_lapack95_mkl
  use Lapack95
#else
#ifdef with_lapack
  use LinAlg_Lapack_350   ! valid for lapack 3.5.0
#endif
#endif

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  SyEvR
  

end module LinAlg
