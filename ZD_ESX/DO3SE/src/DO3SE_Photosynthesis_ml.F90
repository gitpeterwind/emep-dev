module DO3SE_Photosynthesis_ml

  use DO3SE_ModelConstants_ml, only: UNDEF
  use DO3SE_GstoTypes_ml
  use DO3SE_Util_ml
#include "DO3SE_Util_ml.h"

  implicit none
  private

  type, public :: PnGstoConfig_t
    real :: g_sto_0 = UNDEF       !< Closed stomata conductance (umol/m^2/s)
    real :: m = UNDEF             !< Species-specific sensitivity to An (dimensionless)
    real :: V_cmax_25 = UNDEF     !< Maximum catalytic rate at 25 degrees (umol/m^2/s)
    real :: J_max_25 = UNDEF      !< Maximum rate of electron transport at 25 degrees (umol/m^2/s)
  end type PnGstoConfig_t

  public :: check_PnGstoConfig
  public :: gsto_pn

contains

  subroutine check_PnGstoConfig(pn_gsto)
    type(PnGstoConfig_t), intent(inout) :: pn_gsto

    ASSERT_DEFINED(pn_gsto%g_sto_0)
    ASSERT_DEFINED(pn_gsto%m)
    ASSERT_DEFINED(pn_gsto%V_cmax_25)
    ASSERT_DEFINED(pn_gsto%J_max_25)
  end subroutine check_PnGstoConfig

  !> Use Farquar photosynthesis model to calculate stomatal conductance
  !! (mmol O3 m-2 PLA s-1).
  pure real function gsto_pn(pgc, Tair_C, Tleaf_C, u, CO2, RH, PPFD, Lm)
    type(pngstoconfig_t), intent(in) :: pgc   !< Photosynthesis gsto parameters
    real, intent(in) :: Tair_C      !< Air temperature (degrees C)
    real, intent(in) :: Tleaf_C     !< Leaf temperature (degrees C)
    real, intent(in) :: u           !< Wind speed (m/s)
    real, intent(in) :: CO2         !< CO2 concentration (ppm)
    real, intent(in) :: RH          !< Relative humidity (fraction)
    real, intent(in) :: PPFD        !< PPFD (umol/m^2/s)
    real, intent(in) :: Lm          !< Leaf dimension (m)

    real :: An    ! Temporary variable for CO2 assimilation

    call farquar_photosynthesis(Tair_C, Tleaf_C, u, CO2, RH, PPFD, Lm, &
                                pgc%g_sto_0, pgc%m, pgc%V_cmax_25, pgc%J_max_25, &
                                gsto_pn, An)
  end function gsto_pn

  !> Model the net assimilation rate of C3 plants according to the model
  !! developed by Farquar (1980).  Outputs are stomatal conductance and net
  !! CO2 assimilation.
  pure subroutine farquar_photosynthesis(Tair_C, Tleaf_C, uh, c_a, h_a, Q, &
                                         d, g_sto_0, m, V_cmax_25, J_max_25, &
                                         gsto_final, pngsto_An)
    real, intent(in) :: Tair_C        !< Air temperature (degrees C)
    real, intent(in) :: Tleaf_C       !< Leaf temperature (degrees C)
    real, intent(in) :: uh            !< Wind speed (m/s)
    real, intent(in) :: c_a           !< CO2 concentration (ppm)
    real, intent(in) :: h_a           !< Relative humidity (fraction)
    real, intent(in) :: Q             !< PPFD (umol/m^2/s)
    real, intent(in) :: d             !< Leaf dimension (m)
    real, intent(in) :: g_sto_0       !< Closed stomata conductance (umol/m^2/s)
    real, intent(in) :: m             !< Species-specific sensitivity to An (dimensionless)
    real, intent(in) :: V_cmax_25     !< Maximum catalytic rate at 25 degrees (umol/m^2/s)
    real, intent(in) :: J_max_25      !< Maximum rate of electron transport at 25 degrees (umol/m^2/s)

    real, intent(out) :: gsto_final   !< Raw photosynthesis-based stomatal conductance (???)
    real, intent(out) :: pngsto_An    !< Net CO2 assimilation (???)

    ! Constants
    real, parameter :: Ts_K = 273.16    ! Offset between Celcius and Kelvin

    ! parameters considered (or defined) to be constant for all species
    real, parameter :: R = 8.314472          !universal gas constant             [J/(K*mol)]
    real, parameter :: p_O2 = 210.0          !O2 partial pressure                [mmol/mol]
    real, parameter :: E_K_C = 79430.0       !activation energy of K_C           [J/mol]            Medlyn2002
    real, parameter :: E_K_O = 36380.0       !activation energy of K_O           [J/mol]            Medlyn2002
    real, parameter :: E_R_d = 53000.0       !activation energy of R_d           [J/mol]            Leuning1995
    real, parameter :: E_Gamma_star = 37830.0 !activation energy for C-comp-point [J/mol]            Medlyn2002
    real, parameter :: K_C_25 = 404.9        !K.C at reference temperature 25    [micro mol/mol]    Medlyn2002
    real, parameter :: K_O_25 = 278.4        !K.O at reference temperature 25    [mmol/mol]         Medlyn2002
    real, parameter :: R_d_20 = 0.32         !R_d at reference temperature 20    [micro mol/(m^2*s)]Leuning1995
    real, parameter :: Gamma_star_25 = 42.75 !CO2 compensation point at T= 25    [micro mol/mol]    Medlyn2002

    ! species spedific model parameters (that don't tend to have species specific
    ! values, others are supplied as arguments)
    real, parameter :: alpha = 0.3           !efficiency light energy conversion [mol electrons/mol photons]
    real, parameter :: Teta = 0.95           !shape of J~Q determining factor    []
    real, parameter :: H_a_jmax = 50300      !activation energy for J_max        [J/mol]
    real, parameter :: H_d_jmax = 152044     !deactivation energy for J_max      [J/mol]
    real, parameter :: H_a_vcmax = 73637     !activation energy for V_cmax       [J/mol]
    real, parameter :: H_d_vcmax = 149252    !deactivation energy for V_cmax     [J/mol]
    real, parameter :: S_V_vcmax = 486       !entropy terms                      [J/(mol*K)]
    real, parameter :: S_V_jmax = 495        !entropy terms                      [J/(mol*K)

    ! Converted inputs
    real :: T_air, T_leaf, u

    ! state variables
    real :: A_n                             !netto assimilation rate            [micro mol/(m^2*s)]
    real :: A_c                             !Rub. activity. lim. ass. rate      [micro mol/(m^2*s)]
    real :: A_q                             !electr. transp. lim. ass. rate     [micro mol/(m^2*s)]
    real :: Gamma_star                      !CO2 comp. point without day resp.  [micro mol/mol]
    real :: R_d                             !day respiration rate               [micro mol/(m^2*s)]
    real :: K_C                             !Michaelis constant CO2             [micro mol/mol]
    real :: K_O                             !Michaelis constant O2              [micro mol/mol]
    real :: J                               !Rate of electron transport         [micro mol/(m^2*s)]
    real :: ratio                           !TODO: ???
    real :: Gamma                           !CO2 compensation point             [micro mol/mol]
    real :: J_max                           !Max rate of electron transport     [micro mol/(m^2*s)]
    real :: V_cmax                          !Max catalytic rate of Rubisco      [micro mol/(m^2*s)]
    real :: e_sat_i                         !internal saturation vapour pressure[Pa]
    real :: g_bl                            !two sided bound.l. conduct., vapour[micro mol/(m^2s)]
    real :: g_sto                           !two sided stomatal conduct.,vapour [micro mol/(m^2s)]
    real :: g_tot                           !TODO: ???
    real :: h_s                             !relative humidity at leaf surface  [decimal fraction]
    real :: c_s                             !CO2 concentration at leaf surface  [micromol/mol]
    real :: c_i                             !CO2 concentration inside stomata   [micromol/mol]
    real :: e_a                             !ambient vapour pressure            [Pa]


    ! iteration parameters
    integer :: iterations                   !number of the iterations bofore convergence
    real :: c_i_sup                         !CO2 concentration inside stomata possible through supply
    integer :: k                            !loop parameters

    T_air = Tair_C + Ts_K
    T_leaf = Tleaf_C + Ts_K
    u = max(0.01, uh)

    ! Calculation of the model variables which are only
    ! dependend on environmental conditions:

    Gamma_star = Gamma_star_25*exp((E_Gamma_star*(T_leaf-298))/(298*R*T_leaf))

    K_C        = K_C_25*exp((E_K_C*(T_leaf-298))/(298*R*T_leaf))

    K_O        = K_O_25*exp((E_K_O*(T_leaf-298))/(298*R*T_leaf))

    R_d        = R_d_20*exp((E_R_d*(T_leaf-293))/(293*R*T_leaf))

    J_max      = J_max_25*exp((H_a_jmax*(T_leaf-298))/(298*R*T_leaf))*&
                 (1+exp((298*S_V_jmax-H_d_jmax)/(298*R)))/&
                 (1+exp((T_leaf*S_V_jmax-H_d_jmax)/(T_leaf*R)))

    V_cmax     = V_cmax_25*exp((H_a_vcmax*(T_leaf-298))/(298*R*T_leaf))*&
                 (1+exp((298*S_V_vcmax-H_d_vcmax)/(298*R)))/&
                 (1+exp((T_leaf*S_V_vcmax-H_d_vcmax)/(T_leaf*R)))

    J          = ((alpha*Q+J_max)-sqrt((alpha*Q+J_max)**2-4*&
                 (alpha*Q*J_max*Teta)))/(2*Teta)



    e_sat_i    = 613.75*exp((17.502*(T_leaf-273.16))/(240.07+T_leaf-273.16))

    e_a        = h_a*613.75*exp((17.502*(T_air-273.15))/(240.97+T_air-273.15))

    ratio        = (9.81*abs(T_leaf-T_air))/((u**2)*(T_leaf+T_air/2))

    Gamma        = (Gamma_star+(K_C*R_d*(1+(p_O2/K_O))/V_cmax))/&
                   (1-(R_d/V_cmax))

    !aproximates the boundary layer conductance for forced convection
    !and sets the !conductance to a fixed value during free convection:

    if (ratio <= 0.1) then
      g_bl     = 2*0.147*sqrt(u/d)*1e6
    else
      g_bl     = 1.5*1e6
    end if

    !print *, Gamma_star,K_C,K_O,R_d,J_max,V_cmax
    !print *, J,e_sat_i,e_a,ratio
    !print *, Gamma,p_O2

    !The following loop guesses a start value for c_i and tests whether
    !it satisfies all the relevant restrictions. If not a new value for
    !c_i is tested:

    c_i         = 0.0

    ! gsto needs a starting point, so let's set it to g_sto_0
    g_sto = g_sto_0

    do k=1,50

      A_c        = V_cmax *((c_i-Gamma_star)/(c_i + K_C*(1+(p_O2/K_O))))

      A_q        = J*(c_i-Gamma_star)/(4*(c_i+2*Gamma_star))

      A_n        = min(A_c,A_q)-R_d

      c_s        = c_a-(A_n*(1.37/g_bl))

      h_s        = (g_sto*e_sat_i+g_bl*e_a)/(e_sat_i*(g_sto+g_bl))

      g_sto      = g_sto_0 + ( m * A_n *(h_s/(c_s - Gamma)))*1e6

      g_tot      = 1/(1.6/g_sto+1.3/g_bl)

      c_i_sup    = c_a-((A_n/g_tot)*1e6 )

      !exits the loop when c_i calculated with both ways meet the convergence
      !criterium:

      iterations = k
      if (abs(c_i - c_i_sup) < 0.001) then
          exit
      end if

      !Guesses a new c_i as the mean of the first guess and c_i resulting from
      !the supply function:

      c_i      = c_i-(c_i-c_i_sup)/2

    end do

    ! Calculate final stomatal conductances
    gsto_final = max(0.0, g_sto / 1000.0)
    pngsto_An = A_n
  end subroutine farquar_photosynthesis

end module DO3SE_Photosynthesis_ml
