
!>   Definition of z-layer variables, e.g. dzmid, 

module esx_Zgrid
  !use ModelConstants, only: dp
  use esx_Variables, only:  esx

  implicit none
  private

  public ::  init_Zgrid   !> allocates arrays for  z values
  public ::  test_Zgrid   !> tests subroutines in this module

 contains

  !----------------------------------------------------------------------------
  subroutine init_Zgrid()

    integer :: nz
    real, pointer, dimension(:) :: z, dz, zbnd, dzmid

    esx%nz = maxloc(esx%zbnd,dim=1)
    nz     = esx%nz
print *, "NZ ",nz !!, " ZIN:", zin

    z => esx%z(1:nz)
    dz => esx%dz(1:nz)
    zbnd => esx%zbnd(1:nz)
    dzmid => esx%dzmid(1:nz)

!ORIG    z(:)       = zin(:)
!ORIG    zbnd(1:nz-1) = (z(1:nz-1)+z(2:nz))/2    !z_i+½ = (z_i + z_i+1)/2

    !zbnd(:)       = zin(:)
    z(1)    = 0.5*zbnd(1)                           !
    z(2:nz) = (zbnd(1:nz-1)+zbnd(2:nz))/2           !z_i+½ = (z_i + z_i+1)/2
 
    !ORIG dz(2:nz-1) = zbnd(2:nz-1)-zbnd(1:nz-2)  !dz_i = z_i+½ - z_i-½
    !ORIG dz(nz)   = 2*( z(nz)-zbnd(nz-1) )       !dz_n+1 = 2( z_n - z_n-½ )

    dz(1)    = zbnd(1)                              !dz_1 = z_1½ - 0
    dz(2:nz) = zbnd(2:nz)-zbnd(1:nz-1)              !dz_i = z_i+½ - z_i-½
 
    dzmid(1:nz-1) = z(2:nz)-z(1:nz-1)               !dz_i+½ = z_i+1 - z_i

  end subroutine init_Zgrid

  !----------------------------------------------------------------------------
  !> code to test Zgrid routines in this module

  subroutine test_Zgrid(ionum)

    integer, intent(in) :: ionum    ! Unit number for output file
    integer :: i, im
    real, dimension(14 ) :: &
      ztest = (/1.0,2.0,4.0,8.0,16.0,32.0,64.0,100.0,&
                150.0,200.0,250.0,300.0,400.0,500.0/) 
      esx%zbnd = 0.0
      esx%zbnd(1:14) = ztest(:)

    call init_Zgrid()

    write(ionum,"(3a3,4a12)") "nz","n ","n-1","z(i)","zbnd(i)", "dz(i)", "dzmid(i)"
    do i = 1,  esx%nz
      write(ionum,"(3i3,4f12.3)") esx%nz, i, im, &
                                   esx%z(i), esx%zbnd(i), esx%dz(i), esx%dzmid(i)
    end do
    write(ionum,*) "===== done test ========="
    
  end subroutine test_Zgrid
  !----------------------------------------------------------------------------


end module esx_Zgrid

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
!DSX program tester
!DSX  use esx_Zgrid,  only : test_Zgrid ! TESTING only
!DSX  implicit none
!DSX  integer :: ionum
!DSX    open(newunit=ionum,file="TestingZgrid.txt")
!DSX    call  test_Zgrid(ionum)
!DSX    close(ionum)
!DSX end program tester
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
