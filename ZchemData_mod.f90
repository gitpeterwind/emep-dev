! <ZchemData_mod.f90 - part of the EMEP MSC-W Chemical transport Model>
!_____________________________________________________________________________!
 module ZchemData_mod

  ! Arrays of meteorology and concentration for 1-D column , for input to
  ! chemical solver ........
  ! The k-dimension spans the ground (KMAX_MID) to the K-values
  ! specified by KCHEMTOP - here 2
  !
  ! - new aray added to keep o2, m, and for MADE oh, etc

  use AllocInits,    only: AllocInit
  use ChemDims_mod,  only: NCHEMRATES, nspec => NSPEC_TOT, &
      NPHOTOLRATES !A2018 QUERY new usage
  use Config_module, only: KCHEMTOP, KMAX_MID, AERO ! for %NSAREA = No types surface area
  implicit none
  private

!ESX  public :: Alloc1Dchem       !> Sets up arrays of dimension nspec x nz
!ESX  public :: SetChemMet

  !/ variables to keep track of which call

   logical, public, save :: first_call = .true.
   integer, public, save :: ncalls     =      0

  !/-- the chemistry is calculated for arrays of size:

   integer, public, save  :: CHEMSIZE  !

 !Emissions in column. We assume that these only involve advected species
  real, public, save, allocatable, dimension(:,:) :: & ! dims: nspec x nz
      rct               & !< chemical rate coeffs.  molec/cm3/s
     ,rcemis            & !< Emissions rate coeff.  molec/cm3/s
     ,rcbio             & !< Emissions rate coeff.  molec/cm3/s (BVOC, soil-NO, etc.)
!ESX     ,esxemis           & !< Emissions rate coeff.  molec/cm3/s !INTERIM -saves rcemis
!ESX     ,bcemis            & !<  rcemis provided by boundary conditions or tests
!ESX     ,DChem             & !< chemical tendency, molec/cm3/s
     ,rcphot              !< Photolysis rates

  real, public, allocatable, dimension(:,:), save :: &
       xn_2d            ! Concentrations [molecules/cm3]

! For semivolatiles we track the farction as gas and particle- used for SOA
! We use NSPEC_TOT to allow us to write Fpart for FFUEL and WOOD also -
! these may be semivol one day.
   !real, public, dimension(FIRST_SOA:LAST_SOA,KCHEMTOP:KMAX_MID), save :: &
   real, public, allocatable, dimension(:,:), save :: &
       Fgas       &! Fraction as gas-phase
      ,Fpart      ! Fraction as gas-phase

  ! We define a column array for isoprene and terpene for use in
  ! the chemical solver. All values except for k=KMAX_MID will
  ! remain zero however

   real, public, allocatable, dimension(:), save :: &
       rh                  & ! RH (fraction, 0-1)
      ,M                   & ! M - atmospheric conc. (was amk)
      ,o2, n2              & ! oxygen, nitrogen
      ,h2o                 & ! water
      ,temp                & ! temperature
      ,tinv                & ! inverse temp
      ,cN2O5               & ! mol speed, N2O5
      ,cHNO3               & ! mol speed, HNO3 
      ,cHO2                & ! mol speed, HO2  
      ,cNO2                & ! mol speed, NO2    ! kHet tests
      ,cNO3                & ! mol speed, NO2    ! kHet tests
      ,cO3                 & ! mol speed, O3   
      ,gamN2O5               ! jAero gamma values for output

   real, public, allocatable, dimension(:,:), save :: &
       DpgNw  & ! wet diameter,           dim:NSAREA,k
      ,S_m2m3   ! surface area, m2/m3     dim:NSAREA,k

   real, public, allocatable, dimension(:), save :: &
       deltaZcm             & ! layer thickness, cm
      ,aero_fom, aero_fss, aero_fdust, aero_fbc &! fractions
      ,pp                     !pressure
!         ,ugdryPM             & ! for wet radius from Gerber, etc.

   integer, public, allocatable, dimension(:), save :: &
       itemp                  ! int of temperature


!ESX  contains
!ESX  !--------------------------------------------------------------------------!
!ESX   subroutine Alloc1Dchem(nz, nrct, nrcbio, debug_level)
!ESX    integer, intent(in) :: nz, nrct, nrcbio, debug_level
!ESX
!ESX    !> We need to re-allocate if changes in array size, e.g. moving from
!ESX    !! EMEP to ESX
!ESX
!ESX      if( debug_level > 0 ) print *, "ALLOCATE xCHEM?? ", nspec, nz, size(xn_2d, 1)
!ESX
!ESX      !ESX if ( .not. allocated( xn_2d) ) then
!ESX      !ESX     if( debug_level > 0 ) print *, "ALLOCATES xCHEM ", nspec, nz
!ESX      !ESX end if
!ESX
!ESX      call AllocInit( xn_2d,  0.0, nspec, nz, "set1D:xn_2d")
!ESX      call AllocInit( rcemis, 0.0, nspec, nz, "set1D:rcemis",dbg=.true.)
!ESX      call AllocInit( Fgas,  1.0, nspec, nz, "set1D:Fgas")
!ESX      call AllocInit( Fpart, 0.0, nspec, nz, "set1D:Fpart")
!ESX
!ESX      call AllocInit( rcbio,  0.0, nrcbio, nz, "set1D:rcbio")
!ESX      call AllocInit( rct,    0.0, nrct, nz, "set1D:rct")
!ESX      call AllocInit( rcphot, 0.0, MAXRCPHOT, nz, "set1D:rcphot")
!ESX
!ESX     !A2018 Need to check bounds, if it matters...
!ESX
!ESX
 
 end module ZchemData_mod
