!> MODULE Zmet_ml sets "simple" array variables, e.g. temp, H2O, from the 
!! input Zmet type. Needed so that we can switch between ESX meteorology
!! and EMEP/CTM meteorology.

module Zmet_ml
  use AllocInits, only : AllocInit
  use esx_Variables, only : Zmet_t, esx
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
    integer :: iz ! TMP

    call AllocInit(temp, Zdata(1:nz)%tzK, nz, "temp")
    call AllocInit(M,    Zdata(1:nz)%M,   nz, "M")
    call AllocInit(rh,   Zdata(1:nz)%rh,  nz, "rh")

    !> Assume we have RH in Zdata, need to calculate H2O from Teten etc.:
    !! nb 611 is sat. vapour pressure (Pa) at 273K

    esat(:) = 611.0 * exp(0.622*2.5e6*((1.0/273.15) - (1.0/temp(:)))/ RGAS_KG )
    call AllocInit(H2O, 0.0, nz, "H2O")
    H2O(:)  = rh(:)* esat(:)/Zdata(1:nz)%Pa * M(:)*ATWAIR/18.0

    call AllocInit(tinv, 1.0/temp(:), nz, "tinv")
    call AllocInit(N2,   0.79*M(:),   nz, "N2")
    call AllocInit(O2,   0.29*M(:),   nz, "O2")

    ! For chemistry
    call AllocInit(itemp,      nint(temp(:)-1.0e-9), nz, "itemp") ! e-9 avoids some machine problems
    call AllocInit(log300divt, log(300.0*tinv(:)),   nz, "log300divt")
    call AllocInit(logtdiv300, log(temp(:)/300.0),   nz, "logtdiv300")

    call AllocInit(lt300, 0.0, nz, "lt300")   ! TODO: remove this? it isn't used...

    if( esx%debug_Zchem > 1 ) then
      print *, "ZMET t,1/t,rh,M,H2O at z1 ", temp(1), itemp(1), rh(1), M(1), H2O(1)
      print *, "ZMET t,1/t,rh,M,H2O at nz", temp(nz), itemp(nz), rh(nz), M(nz), H2O(nz)
    end if
    !stop

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

