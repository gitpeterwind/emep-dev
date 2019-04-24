!>_________________________________________________________<

  module  GenOut_SOA_ml
!-----------------------------------------------------------

  use ChemSpecs_tot_ml  ! => NSPEC_TOT, species indices
  implicit none
  private
!+ Defines SOA, NONVOL and VBS params 

  type, public :: VBST
       integer     :: index    ! just for clarity
       real        :: CiStar   ! ug/m3
       !real        :: Tref    ! Assumed 300
       real        :: DeltaH   ! kJ/mole
  end type VBST


 end module GenOut_SOA_ml
 !-----------------------------------------------------------
