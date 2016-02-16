module MPI_Groups_ml

  USE mpi, only:MPI_REAL8,MPI_DOUBLE_PRECISION, MPI_SUM,MPI_INTEGER, MPI_BYTE,MPI_CHARACTER, &
                MPI_LOR,MPI_LOGICAL,MPI_MAX,MPI_MIN, MPI_COMM_WORLD,MPI_IN_PLACE,&
                MPI_ADDRESS_KIND,MPI_INFO_NULL,MPI_STATUS_SIZE
  implicit none
  integer, public,save :: MPI_COMM_CALC, MPI_COMM_IO, MPI_COMM_SUB
  integer :: MPIInfo, MPISTATUS(MPI_STATUS_SIZE)
  integer, public :: request_ps_w, request_ps_e, request_xn_w, request_xn_e
  integer, public :: request_ps_s, request_ps_n, request_xn_s, request_xn_n
  integer, public :: request_s, request_n, request_w, request_e
  integer, public :: IERROR

!dummy
  integer, public ::MPI_groups_split, MPI_COMM_TYPE_SHARED
  integer, public, save :: ME_MPI,ME_IO,ME_CALC,ME_SUB,largeLIMAX,largeLJMAX
  integer, public, save:: NPROCX_IO,NPROCY_IO,NPROC_IO,NPROCX_SUB,NPROCY_SUB,NPROC_SUB
  integer, public, save:: NPROC_MPI
  integer, public, save:: LargeSub_Ix


  public :: MPI_world_init

contains
subroutine MPI_world_init(NPROC,ME)

integer,intent(out) ::NPROC,ME

  CALL MPI_INIT(IERROR)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, IERROR)
  ME_MPI=ME
  ME_CALC=ME
  ME_SUB=ME
  ME_IO=-1
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NPROC, IERROR)
  NPROC_MPI=NPROC
  MPI_COMM_CALC=MPI_COMM_WORLD
  MPI_COMM_SUB=MPI_COMM_WORLD
55 format(A,I5,A)
  if(ME==0)write(*,55)' Found ',NPROC,' MPI processes available'

end subroutine MPI_world_init
subroutine share(shared_data,data_shape,xsize,MPI_COMM_SHARED)

!share the array shared_data
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER!fortran 2003 extensions 
  implicit none
  TYPE(C_PTR) :: baseptr,baseptr2
!  TYPE(MPI_Win) :: win
!  TYPE(MPI_Comm), intent(in) :: MPI_COMM_SHARED
  integer :: win
  integer, intent(in) :: MPI_COMM_SHARED
  real , dimension(:,:,:),pointer, intent(inout) :: shared_data
  INTEGER(KIND=MPI_ADDRESS_KIND) :: MPI_XSIZE
  INTEGER, intent(inout) :: XSIZE
  INTEGER, intent(in) :: data_shape(3)
  INTEGER DISP_UNIT, IERROR,me_shared
  integer :: i,j,n,data_size
  real :: ONE

  DISP_UNIT=sizeof(ONE)!8, number of bytes for real

  CALL MPI_COMM_RANK(MPI_COMM_SHARED, ME_shared, IERROR)
  
  nullify(shared_data)
  MPI_XSIZE=XSIZE*DISP_UNIT
  data_size=1
  do i=1,size(data_shape)
     data_size=data_size*data_shape(i)
  enddo
  if(data_size/=XSIZE)then
     write(*,*)'WARNING: incompatible dimensions in MPI_groups_ml ',data_size,XSIZE,data_shape
  endif
  if(ME_shared/=0) MPI_XSIZE = 0
  
!  CALL MPI_WIN_ALLOCATE_SHARED(MPI_XSIZE, DISP_UNIT, MPI_INFO_NULL, MPI_COMM_SHARED, BASEPTR2, WIN,IERROR)

  call MPI_Win_fence(0, win, ierror)
!  CALL MPI_Win_shared_query(win, 0, MPI_xsize, disp_unit, baseptr,IERROR)
  CALL C_F_POINTER(baseptr, shared_data, data_shape)
  call MPI_Win_fence(0, win, ierror)

!test if it works
!  shared_data(1,1,1)=0
!  shared_data(2,1,1)=0
  if(me_mpi==1000)then
     shared_data(1,1,1)=111.
  elseif(me_mpi==0)then
     shared_data(1,1,1)=11.
!  if(me_io==0.and.me_sub==0)shared_data(2,1,1)=0
!  if(me_io==1.and.me_sub==0)shared_data(2,2,1)=4
!  if(me_io==1.and.me_sub==0)shared_data(1,1,1)=0.
  else if(me_mpi==10)then
     shared_data(2,1,1)=22.
  else
     shared_data(3,1,1)=me_mpi
  endif
  CALL MPI_BARRIER(MPI_COMM_SHARED, IERROR)
  call MPI_Win_fence(0, win, ierror)
 78 format(A,5i7,12F11.2)
!  write(*,78)'data in share ',ME_MPI,me_calc,me_io,me_sub,me_shared,shared_data(1,1,1),&
!shared_data(2,1,1),shared_data(3,1,1),1.0*mpi_xsize!,shared_data(1,2,1),shared_data(2,2,1)
!  if(me_io>=0)write(*,*)' COMM',MPI_COMM_SHARED,MPI_COMM_WORLD

end subroutine share
subroutine share_logical(shared_data,data_shape,xsize,MPI_COMM_SHARED)

!share the array shared_data
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER!fortran 2003 extensions 
  implicit none
  TYPE(C_PTR) :: baseptr
!  TYPE(MPI_Win) :: win
!  TYPE(MPI_Comm), intent(in) :: MPI_COMM_SHARED
  integer :: win
  integer, intent(in) :: MPI_COMM_SHARED
  logical , pointer, intent(inout) :: shared_data
  INTEGER(KIND=MPI_ADDRESS_KIND) :: MPI_XSIZE
  INTEGER, intent(inout) :: XSIZE
  INTEGER, intent(in) :: data_shape(1)
  INTEGER DISP_UNIT, IERROR,me_shared
  integer :: i,j,n,data_size
  real :: ONE
  logical :: mybool

  DISP_UNIT=1!number of bytes for logical
  CALL MPI_COMM_RANK(MPI_COMM_SHARED, ME_shared, IERROR)

  nullify(shared_data)
  MPI_XSIZE=XSIZE*sizeof(mybool)
  if(me_mpi==0)write(*,*)'size of logical ',sizeof(mybool)
  data_size=1
  do i=1,size(data_shape)
     data_size=data_size*data_shape(i)
  enddo
  if(data_size/=XSIZE)then
     write(*,*)'WARNING: incompatible dimensions in MPI_groups_ml ',data_size,XSIZE,data_shape
  endif

!  CALL MPI_WIN_ALLOCATE_SHARED(MPI_XSIZE, DISP_UNIT, MPI_INFO_NULL, MPI_COMM_SHARED, BASEPTR, WIN, IERROR)
  call MPI_Win_fence(0, win, ierror)
!  CALL MPI_Win_shared_query(win, 0, MPI_xsize, disp_unit, baseptr, IERROR)
  CALL C_F_POINTER(baseptr, shared_data)
  call MPI_Win_fence(0, win, ierror)

!test if it works
  shared_data=.true.
  call MPI_Win_fence(0, win, ierror)
  if(me_io==1.and.me_sub==0)shared_data=.true.
  if(me_io==0.and.me_sub==0)shared_data=.false.
  call MPI_Win_fence(0, win, ierror)
 78 format(A,2i4,12F7.2)
!  write(*,*)'logical data in share ',ME_MPI,me_shared,shared_data
end subroutine share_logical
end module MPI_Groups_ml
