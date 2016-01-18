module MPI_Groups_ml

  USE mpi
  implicit none
  integer, public,save :: MPI_COMM_CALC, MPI_COMM_IO, MPI_COMM_SUB
  integer, dimension(MPI_STATUS_SIZE), public :: MPISTATUS
  integer :: MPIInfo
  integer, public :: request_ps_w, request_ps_e, request_xn_w, request_xn_e
  integer, public :: request_ps_s, request_ps_n, request_xn_s, request_xn_n
  integer, public :: request_s, request_n, request_w, request_e
  integer, public :: IERROR

  public :: MPI_world_init

contains
subroutine MPI_world_init(NPROC,ME)

integer,intent(out) ::NPROC,ME

  CALL MPI_INIT(IERROR)
  CALL MPI_COMM_RANK(MPI_COMM_CALC, ME, IERROR)

  CALL MPI_COMM_SIZE(MPI_COMM_CALC, NPROC, IERROR)
55 format(A,I5,A)
  if(ME==0)write(*,55)' Found ',NPROC,' MPI processes available'

end subroutine MPI_world_init
end module MPI_Groups_ml
