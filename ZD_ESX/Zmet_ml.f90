!> MODULE Zmet_ml sets "simple" array variables, e.g. temp, H2O, from the 
!! input Zmet type. Needed so that we can switch between ESX meteorology
!! and EMEP/CTM meteorology.

module Zmet_ml
  use esx_Variables, only : Zmet_t
  implicit none
  private
  
  public :: Set1Dmet

  !COMPILER PROBLEM: real, pointer, public, save,  dimension(:) :: temp, ..
  real, public, save,  allocatable, dimension(:) :: temp, rh, tinv, lt300, M, N2, O2, H2O
  real, public, save,  allocatable, dimension(:) :: log300divt, logtdiv300
  integer, public, save,  allocatable, dimension(:) :: itemp  ! inverse T

contains

  subroutine Set1Dmet(nz, Zdata)
    integer, intent(in) :: nz
    type(Zmet_t), dimension(:), target, intent(in) :: Zdata
    real, dimension(nz) :: esat            ! helper variables for H2O
    real, parameter :: RGAS_KG = 287.0     ! Molar Gas constant (J K-1 kg-1)
    real, parameter :: ATWAIR = 28.964     ! mol wt of air, g/mol
 
    !> We need to re-allocate if changes in array size
    !! Should be eonough to test just temperature.
      if ( allocated( temp) ) then
        if ( size( temp) /= nz ) deallocate(temp, rh, tinv, lt300, M, O2, H2O,&
                                            log300divt, logtdiv300, itemp )
      end if
      if ( .not. allocated( temp) ) then
          allocate(temp(nz), rh(nz), tinv(nz), lt300(nz), M(nz), O2(nz), H2O(nz) )
          allocate(log300divt(nz), logtdiv300(nz)) ! For chemistry
          allocate(itemp(nz)) ! For chemistry
      end if

      temp = Zdata(1:nz)%tzK
      itemp = nint( temp -1.0e-9 ) ! e-9 avoids some machine problems
      M    = Zdata(1:nz)%M
      rh   = Zdata(1:nz)%rh

      !> Assume we have RH in Zdata, need to calculate H2O from Teten etc.:
      !! nb 611 is sat. vapour pressure (Pa) at 273K

      esat(:) = 611.0 * exp(0.622*2.5e6*((1.0/273.15) - (1.0/temp(:)))/ RGAS_KG )
      H2O(:)  = rh(:)* esat(:)/Zdata(1:nz)%Pa * M(:)*ATWAIR/18.0


      tinv = 1.0/temp
      N2   = 0.79 * M
      O2   = 0.29 * M

      log300divt(:) = log(300.0*tinv(:))
      logtdiv300(:) = log(temp(:)/300.0)

     
  end subroutine Set1Dmet

end module Zmet_ml

!DSXprogram testz
!DSX  use esx_Variables, only : Zmet
!DSX  use Zmet_ml
!DSX  implicit none
!DSX  integer :: i
!DSX  do i = 1, size(Zmet)
!DSX    Zmet(i)%tzK = 1.0*i
!DSX  end do
!DSX  call Set1Dmet(10, Zmet)
!DSX  print *, "TEST1D", size(Zmet), size(temp), temp(4), Zmet(4)%tzK
!DSXend program testz

