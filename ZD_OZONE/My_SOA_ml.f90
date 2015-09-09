! <SOA_ml.f90 - A component of the EMEP MSC-W Chemical transport Model>
!*****************************************************************************! 
module OrganicAerosol_ml

   !--------------------------------------------------------------------------
   ! This module is fake - for initial 2011 public-domain ozone model only, pending
   ! decision as to which SOA scheme to release as default.
   ! Contact David.Simpson@met.no for more information if interested in SOA
   ! schemes
   !--------------------------------------------------------------------------
   use ChemSpecs,  only:  A1 => FIRST_SEMIVOL , A2 => LAST_SEMIVOL, species
   use ModelConstants_ml,    only : CHEMTMIN, CHEMTMAX, &
                                    K2 => KMAX_MID, K1 => KCHEMTOP
   use PhysicalConstants_ml, only : AVOG
   use Setup_1dfields_ml,    only : itemp, xn => xn_2d
   implicit none

   !/-- subroutines
    public  :: Init_OrganicAerosol
    public  :: OrganicAerosol


   !/-- public

    logical, public, parameter :: ORGANIC_AEROSOLS = .false.
! From Setup_1dfields now
!   real, public, dimension(A1:A2,K1:K2), save :: Fgas  ! Fraction in gas-phase

    character(len=*), public, parameter :: SOA_MODULE_FLAG="NotUsed"

   contains
   !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   !+ Driver routine for Secondary Organic Aerosol  module

   subroutine Init_OrganicAerosol(i,j,debug_flag)
     integer, intent(in) :: i,j
     logical, intent(in) :: debug_flag  ! for debugging purposes only

     ! empty 
     
   end subroutine Init_OrganicAerosol

   subroutine OrganicAerosol(i,j,debug_flag)
     integer, intent(in) :: i,j
     logical, intent(in) :: debug_flag  ! for debugging purposes only

     ! empty 

   end subroutine OrganicAerosol

end module OrganicAerosol_ml
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

