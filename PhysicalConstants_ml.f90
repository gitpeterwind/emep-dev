      module PhysicalConstants_ml
!----------------------------------------------------------------------------
!  Defines Physical constants 
!----------------------------------------------------------------------------
implicit none

!-- contains no subroutine:

!

  real , public, parameter ::         &
    AVOG   = 6.023e23                 & ! Avogadros number
  , ATWAIR = 28.964                   & ! mol wt of air, g/mol
  , RGAS_ATML = 0.08205               & ! Molar Gas constant (atm M-1 K-1)
  , RGAS_KG   = 287.0                 & ! Molar Gas constant (J K-1 kg-1)
  , RGAS_J    = 8.314                   ! Molar Gas constant (J mol-1 K-1)

                                        ! NB. ( J = N m2 = kg m2 s-2 )
                                        !       M = mol l-1

!  additional parameters, former set in defcon, from former eulcon.inc
!
  real, public, parameter  ::    &
       GRAV    = 9.807           &   ! Gravity, m s-2
    ,  CP      = 1004.0          &   ! Specific heat at const. pressure
    ,  R       = 287.0           &   ! Gas constants J K-1 kg-1
    ,  XKAP    = R/CP            &
    ,  KARMAN  = 0.41            &   ! Von Karman  (=0.35 elsehwere in code!)
    ,  PI      = 3.141592653589793238462643383279 & ! www.verbose.net/Pi.html
    ,  DEG2RAD = PI/180.0        &   ! COnverts degrees to radians
    ,  ROWATER = 1000.           &   ! pw density of water kg m-3
    ,  BOLTZMANN = 1.380e-23     &   ! Boltzmann'c constant[J/deg/molec]
    !ds ,  NU        = 1.46e-5       &   ! kinematic air viscosity[m2/s]
    ,  FREEPATH  = 6.5e-8        &   ! Mean Free Path of air [m]
    ,  VISCO     = 1.46e-5           ! Air viscosity [m2/s]

!=================== DEP CODE ================================Y:0

  ! CHARNOCK is used to calculate the roughness length for the 
  ! landuse category water

   real, public, parameter  :: &
       PRANDTL = 0.71,            &   ! Prandtl number (see Garratt, 1992)
       Sc_H20  = 0.6,             &   ! Schmidt number for water
    CHARNOCK = 0.032   ! Charnock's alpha:
                       ! see Nordeng (1986), p.31, 
                       ! Nordeng(1991), JGR, 96, no. C4, pp. 7167-7174.
                       ! In the second of these publications, Nordeng uses
                       ! "m" to denote Charnock's alpha whilst in the first
                       ! he specifies the value 0.032.

  ! Standard temperature :

  real, public, parameter :: T0 = 273.15   ! zero degrees Celsius in Kelvin 

!=================== DEP CODE ================================Y:0



end  module PhysicalConstants_ml

!REMOVED
!    ATWAIR = 0.79*28.0 + 0.21*32.0    & ! Atomic weight of air


