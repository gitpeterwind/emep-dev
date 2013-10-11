!>_________________________________________________________<

  module  ChemRates
!-----------------------------------------------------------

   
  use ChemFunctions       ! => kaero, RiemerN2O5

  use Zmet_ml        ! => tinv, h2o, m, Fgas !ESX
  use ZchemData      ! => rct
  use ChemSpecs      ! => NSPEC_TOT, PINALD, .... for FgasJ08
  use DefPhotolysis  ! => IDNO2 etc.
  implicit none
  private

  !> CM_ChemRates RCTYPE rct sets Rate-coefficients - temperature dependant

    public :: setchemrates

     integer, parameter, public :: NCHEMRATES = 9   !! No. coefficients

!> Photolysis rates
     integer, parameter, public :: NPHOTOLRATES = 3   !! No. DJ vals used
     integer, parameter, public,dimension(NPHOTOLRATES) :: photol_used= (/&
         IDBO3 &
        ,IDNO2 &
        ,IDCH3CHO &
  /)

  contains
  !------------------------------------
  subroutine setchemrates(debug_level) 
    integer, intent(in) :: debug_level

       rct(1,:) = 6.0e-34*M*O2*(TEMP/300.0)**2.6 
       rct(2,:) = 2.2e-10*H2O 
       rct(3,:) = 1.4e-12*exp(-1310.0*TINV) 
       rct(4,:) = 2.03e-17*(TEMP**2)*exp(78.0*TINV) 
       rct(5,:) = 4.4e-12*exp(365.0*TINV) 
       rct(6,:) = 3.6e-12*exp(270.0*TINV) 
       rct(7,:) = 2.54e-12*exp(360.0*TINV) 
       rct(8,:) = 7.5e-12*exp(290.0*TINV) 
       rct(9,:) = 1.95e+16*exp(-13543.0*TINV) 

  end subroutine setchemrates
end module  ChemRates
