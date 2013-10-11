!> MODULE AllocInits allocates arrays, and initialises to a given value, 
!! Used to prevent common errors when allocations are made without being
!! initialised. (gfortran doesn't catch this.)

module AllocInits
  use CheckStops, only : CheckStop
  implicit none
  private

  public :: AllocInit
  private :: Alloc1d  ! Alloc/init for 1-D arrays
  private :: Alloc2d  ! Alloc/init for 2-D arrays
  private :: Alloc3d  ! Alloc/init for 3-D arrays

  interface AllocInit
    module procedure Alloc1d
    module procedure Alloc1d_int
    module procedure Alloc2d
    module procedure Alloc3d
  end interface AllocInit

contains

  subroutine Alloc1d( a, init, n1, txt)
    real, allocatable, dimension(:), intent(inout) :: a
    real, intent(in) :: init
    character(len=*), intent(in) :: txt
    integer :: n1, istat

    allocate( a(n1), stat = istat )
    call CheckStop ( istat /= 0, "ERROR 1d Alloc: "//trim(txt) )
    a = init
  end subroutine Alloc1d

  subroutine Alloc1d_int( a, init, n1, txt)
    integer, allocatable, dimension(:), intent(inout) :: a
    integer, intent(in) :: init
    character(len=*), intent(in) :: txt
    integer :: n1, istat

    allocate( a(n1), stat = istat )
    call CheckStop ( istat /= 0, "ERROR 1d_int Alloc: "//trim(txt) )
    a = init
  end subroutine Alloc1d_int

  subroutine Alloc2d( a, init, n1,n2, txt)
    real, allocatable, dimension(:,:), intent(inout) :: a
    real, intent(in) :: init
    character(len=*), intent(in) :: txt
    integer :: n1, n2, istat

    allocate( a(n1,n2), stat = istat )
    call CheckStop ( istat /= 0, "ERROR 2d Alloc: "//trim(txt) )
    a = init
  end subroutine Alloc2d

  subroutine Alloc3d( a, init, n1,n2,n3, txt)
    real, allocatable, dimension(:,:,:), intent(inout) :: a
    real, intent(in) :: init
    character(len=*), intent(in) :: txt
    integer :: n1, n2,n3, istat

    allocate( a(n1,n2,n3), stat = istat )
    call CheckStop ( istat /= 0, "ERROR 2d Alloc: "//trim(txt) )
    a = init
  end subroutine Alloc3d
end module AllocInits
