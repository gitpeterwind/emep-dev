!########### OZONE model ###################################
!ToDo
! Sort out where ORG_AEROSOLS etc should be!
!>_________________________________________________________<

  module  GenSpec_bgn_ml
!-----------------------------------------------------------
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP
  implicit none
  private

!+ Defines indices and NSPEC for bgn : Background species

 ! Species which can be specified simply for each column, e.g.
 ! as function of local meteorology or zenith angle
 !   o2, m,  and for MADE-like, oh, ch3coo2

   integer, public, parameter ::  NSPEC_BGN = 0 ! No. 3D bgn species
   integer, public, parameter ::  NSPEC_COL = 0 ! total no. prescribed specs

  !/ u2 define xn_2d_bgn here.
   !u3 real, public, save, dimension(NSPEC_COL,KCHEMTOP:KMAX_MID) :: xn_2d_bgn
   !u3 assign dimesnion 1 for compilation only
   real, public, save, dimension(1,KCHEMTOP:KMAX_MID) :: xn_2d_bgn

  !/ In ACID we now specify IXBGN_xx values. Omitted here.
  !REMds real, public, parameter :: O3fix=10. !add 10 ppb to Logan data
!-----------------------------------------------------------
  end module GenSpec_bgn_ml


!>_________________________________________________________<

  module  GenSpec_adv_ml
!-----------------------------------------------------------
  implicit none
!/ ...... ..   ( from GenChem )

!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 56 
 
 ! Aerosols:

   logical, public, parameter ::  ORG_AEROSOLS = .false. 

   integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  



   integer, public, parameter ::   & 
     IXADV_O3          =   1   &
  ,  IXADV_NO          =   2   &
  ,  IXADV_NO2         =   3   &
  ,  IXADV_PAN         =   4   &
  ,  IXADV_MPAN        =   5   &
  ,  IXADV_NO3         =   6   &
  ,  IXADV_N2O5        =   7   &
  ,  IXADV_ISONO3      =   8   &
  ,  IXADV_HNO3        =   9

   integer, public, parameter ::   & 
     IXADV_CH2CCH3     =  10   &
  ,  IXADV_CH3COO2     =  11   &
  ,  IXADV_MACR        =  12   &
  ,  IXADV_ISNI        =  13   &
  ,  IXADV_ISNIR       =  14   &
  ,  IXADV_GLYOX       =  15   &
  ,  IXADV_MGLYOX      =  16   &
  ,  IXADV_MAL         =  17   &
  ,  IXADV_MEK         =  18   &
  ,  IXADV_MVK         =  19

   integer, public, parameter ::   & 
     IXADV_HCHO        =  20   &
  ,  IXADV_CH3CHO      =  21   &
  ,  IXADV_C2H6        =  22   &
  ,  IXADV_NC4H10      =  23   &
  ,  IXADV_C2H4        =  24   &
  ,  IXADV_C3H6        =  25   &
  ,  IXADV_OXYL        =  26   &
  ,  IXADV_ISOP        =  27   &
  ,  IXADV_CH3O2H      =  28   &
  ,  IXADV_C2H5OOH     =  29

   integer, public, parameter ::   & 
     IXADV_BURO2H      =  30   &
  ,  IXADV_ETRO2H      =  31   &
  ,  IXADV_PRRO2H      =  32   &
  ,  IXADV_OXYO2H      =  33   &
  ,  IXADV_MEKO2H      =  34   &
  ,  IXADV_MALO2H      =  35   &
  ,  IXADV_MVKO2H      =  36   &
  ,  IXADV_MARO2H      =  37   &
  ,  IXADV_ISRO2H      =  38   &
  ,  IXADV_H2O2        =  39

   integer, public, parameter ::   & 
     IXADV_CH3COO2H    =  40   &
  ,  IXADV_CH2CO2HCH3  =  41   &
  ,  IXADV_ISONO3H     =  42   &
  ,  IXADV_ISNIRH      =  43   &
  ,  IXADV_CH3OH       =  44   &
  ,  IXADV_C2H5OH      =  45   &
  ,  IXADV_H2          =  46   &
  ,  IXADV_CO          =  47   &
  ,  IXADV_CH4         =  48   &
  ,  IXADV_SO2         =  49

   integer, public, parameter ::   & 
     IXADV_SO4         =  50   &
  ,  IXADV_pNO3        =  51   &
  ,  IXADV_NH3         =  52   &
!hf  ,  IXADV_AMSU        =  53   &
!hf  ,  IXADV_AMNI        =  54
  ,  IXADV_aNH4        =   53   & !total NH4
  ,  IXADV_aNO3        =   54   & !total particulate nitrate (in UNI-OZONE: -NO3 unspecified)
!st
  ,  IXADV_PM25        =   55  &
  ,  IXADV_PMco        =   56    
 !-----------------------------------------------------------
  end module GenSpec_adv_ml
!>_________________________________________________________<

  module  GenSpec_shl_ml
!-----------------------------------------------------------
  implicit none
!/ ...... ..   ( from GenChem )

!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_SHL = 15 
 
 ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  



   integer, public, parameter ::   & 
     IXSHL_OD          =   1   &
  ,  IXSHL_OP          =   2   &
  ,  IXSHL_OH          =   3   &
  ,  IXSHL_HO2         =   4   &
  ,  IXSHL_CH3O2       =   5   &
  ,  IXSHL_C2H5O2      =   6   &
  ,  IXSHL_SECC4H9O2   =   7   &
  ,  IXSHL_ISRO2       =   8   &
  ,  IXSHL_ETRO2       =   9

   integer, public, parameter ::   & 
     IXSHL_PRRO2       =  10   &
  ,  IXSHL_OXYO2       =  11   &
  ,  IXSHL_MEKO2       =  12   &
  ,  IXSHL_MALO2       =  13   &
  ,  IXSHL_MVKO2       =  14   &
  ,  IXSHL_MACRO2      =  15

 !-----------------------------------------------------------
  end module GenSpec_shl_ml
!>_________________________________________________________<

  module  GenSpec_tot_ml
!-----------------------------------------------------------

  
  implicit none
!/ ...... ..   ( from GenChem )

!+ Defines indices and NSPEC for tot : All reacting species 

!   ( Output from GenChem, sub print_species ) 

   logical, public, parameter ::  ORG_AEROSOLS = .false. 

   integer, public, parameter ::  NSPEC_TOT = 71 
 
 ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  



   integer, public, parameter ::   & 
     OD          =   1   &
  ,  OP          =   2   &
  ,  OH          =   3   &
  ,  HO2         =   4   &
  ,  CH3O2       =   5   &
  ,  C2H5O2      =   6   &
  ,  SECC4H9O2   =   7   &
  ,  ISRO2       =   8   &
  ,  ETRO2       =   9

   integer, public, parameter ::   & 
     PRRO2       =  10   &
  ,  OXYO2       =  11   &
  ,  MEKO2       =  12   &
  ,  MALO2       =  13   &
  ,  MVKO2       =  14   &
  ,  MACRO2      =  15   &
  ,  O3          =  16   &
  ,  NO          =  17   &
  ,  NO2         =  18   &
  ,  PAN         =  19

   integer, public, parameter ::   & 
     MPAN        =  20   &
  ,  NO3         =  21   &
  ,  N2O5        =  22   &
  ,  ISONO3      =  23   &
  ,  HNO3        =  24   &
  ,  CH2CCH3     =  25   &
  ,  CH3COO2     =  26   &
  ,  MACR        =  27   &
  ,  ISNI        =  28   &
  ,  ISNIR       =  29

   integer, public, parameter ::   & 
     GLYOX       =  30   &
  ,  MGLYOX      =  31   &
  ,  MAL         =  32   &
  ,  MEK         =  33   &
  ,  MVK         =  34   &
  ,  HCHO        =  35   &
  ,  CH3CHO      =  36   &
  ,  C2H6        =  37   &
  ,  NC4H10      =  38   &
  ,  C2H4        =  39

   integer, public, parameter ::   & 
     C3H6        =  40   &
  ,  OXYL        =  41   &
  ,  ISOP        =  42   &
  ,  CH3O2H      =  43   &
  ,  C2H5OOH     =  44   &
  ,  BURO2H      =  45   &
  ,  ETRO2H      =  46   &
  ,  PRRO2H      =  47   &
  ,  OXYO2H      =  48   &
  ,  MEKO2H      =  49

   integer, public, parameter ::   & 
     MALO2H      =  50   &
  ,  MVKO2H      =  51   &
  ,  MARO2H      =  52   &
  ,  ISRO2H      =  53   &
  ,  H2O2        =  54   &
  ,  CH3COO2H    =  55   &
  ,  CH2CO2HCH3  =  56   &
  ,  ISONO3H     =  57   &
  ,  ISNIRH      =  58   &
  ,  CH3OH       =  59

   integer, public, parameter ::   & 
     C2H5OH      =  60   &
  ,  H2          =  61   &
  ,  CO          =  62   &
  ,  CH4         =  63   &
  ,  SO2         =  64   &
  ,  SO4         =  65   &
  ,  pNO3        =  66   &
  ,  NH3         =  67   &
!hf  ,  AMSU        =  68   &
!hf  ,  AMNI        =  69
  ,  aNH4        =   68   &
  ,  aNO3        =   69  &
  ,  PM25        =   70  &
  ,  PMco        =   71   
 !-----------------------------------------------------------
  end module GenSpec_tot_ml
!>_________________________________________________________<
!>_________________________________________________________<

  module  GenChemicals_ml
!-----------------------------------------------------------

   use GenSpec_tot_ml, only: NSPEC_TOT    ! Total number of species for chemistry
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !/--   Characteristics of species: 
  !/--   Number, name, molwt, carbon num, nmhc (1) or not(0)
 
  public :: define_chemicals    ! Sets names, molwts, carbon num, advec, nmhc

  type, public :: Chemical 
       character(len=12) :: name
       integer           :: molwt
       integer           :: nmhc      ! nmhc (1) or not(0)
       integer           :: carbons   ! Carbon-number
       real              :: nitrogens ! Nitrogen-number
       integer           :: sulphurs  ! Sulphur-number
  end type Chemical
  type(Chemical), public, dimension(NSPEC_TOT) :: species

  contains
    subroutine define_chemicals()
    !+
    ! Assigns names, mol wts, carbon numbers, advec,  nmhc to user-defined Chemical
    ! array, using indices from total list of species (advected + short-lived).
    !
    ! ds1 - Carbon's added in 5/12/2003 - for VOC calculation
    !                                           MW NMHC  C    N   S
       species(  1) = Chemical("OD          ",  16,  0,  0,   0,  0 ) 
       species(  2) = Chemical("OP          ",  16,  0,  0,   0,  0 ) 
       species(  3) = Chemical("OH          ",  17,  0,  0,   0,  0 ) 
       species(  4) = Chemical("HO2         ",  33,  0,  0,   0,  0 ) 
       species(  5) = Chemical("CH3O2       ",  47,  0,  1,   0,  0 ) 
       species(  6) = Chemical("C2H5O2      ",  61,  0,  2,   0,  0 ) 
       species(  7) = Chemical("SECC4H9O2   ",  89,  0,  4,   0,  0 ) 
       species(  8) = Chemical("ISRO2       ", 101,  0,  5,   0,  0 ) 
       species(  9) = Chemical("ETRO2       ",  77,  0,  2,   0,  0 ) 
       species( 10) = Chemical("PRRO2       ",  91,  0,  3,   0,  0 ) 
       species( 11) = Chemical("OXYO2       ",   1,  0,  8,   0,  0 ) !ds1
       species( 12) = Chemical("MEKO2       ", 103,  0,  4,   0,  0 ) 
       species( 13) = Chemical("MALO2       ", 147,  0,  5,   0,  0 ) 
       species( 14) = Chemical("MVKO2       ", 119,  0,  4,   0,  0 ) 
       species( 15) = Chemical("MACRO2      ", 102,  0,  4,   0,  0 ) 
       species( 16) = Chemical("O3          ",  48,  0,  0,   0,  0 ) 
       species( 17) = Chemical("NO          ",  30,  0,  0,   1,  0 ) 
       species( 18) = Chemical("NO2         ",  46,  0,  0,   1,  0 ) 
       species( 19) = Chemical("PAN         ", 121,  0,  2,   1,  0 ) 
       species( 20) = Chemical("MPAN        ", 132,  0,  4,   1,  0 ) 
       species( 21) = Chemical("NO3         ",  62,  0,  0,   1,  0 ) 
       species( 22) = Chemical("N2O5        ", 108,  0,  0,   2,  0 ) 
       species( 23) = Chemical("ISONO3      ", 110,  0,  5,   1,  1 ) !ds1
       species( 24) = Chemical("HNO3        ",  63,  0,  0,   1,  0 ) 
       species( 25) = Chemical("CH2CCH3     ",  73,  0,  3,   0,  0 ) 
       species( 26) = Chemical("CH3COO2     ",  75,  0,  2,   0,  0 ) 
       species( 27) = Chemical("MACR        ",  70,  0,  4,   0,  0 ) 
       species( 28) = Chemical("ISNI        ",  46,  0,  4,   1,  1 ) !ds1
       species( 29) = Chemical("ISNIR       ",  46,  0,  4,   1,  1 ) !ds1 
       species( 30) = Chemical("GLYOX       ",  58,  0,  2,   0,  0 ) 
       species( 31) = Chemical("MGLYOX      ",  72,  0,  3,   0,  0 ) 
       species( 32) = Chemical("MAL         ",  98,  0,  5,   0,  0 ) 
       species( 33) = Chemical("MEK         ",  72,  0,  4,   0,  0 ) 
       species( 34) = Chemical("MVK         ",  70,  0,  4,   0,  0 ) 
       species( 35) = Chemical("HCHO        ",  30,  0,  1,   0,  0 ) 
       species( 36) = Chemical("CH3CHO      ",  44,  0,  2,   0,  0 ) 
       species( 37) = Chemical("C2H6        ",  30,  1,  2,   0,  0 ) 
       species( 38) = Chemical("NC4H10      ",  58,  1,  4,   0,  0 ) 
       species( 39) = Chemical("C2H4        ",  28,  1,  2,   0,  0 ) 
       species( 40) = Chemical("C3H6        ",  42,  1,  3,   0,  0 ) 
       species( 41) = Chemical("OXYL        ", 106,  1,  8,   0,  0 ) 
       species( 42) = Chemical("ISOP        ",  68,  1,  5,   0,  0 ) 
       species( 43) = Chemical("CH3O2H      ",  48,  0,  1,   0,  0 ) 
       species( 44) = Chemical("C2H5OOH     ",  62,  0,  2,   0,  0 ) 
       species( 45) = Chemical("BURO2H      ",  90,  0,  4,   0,  0 ) 
       species( 46) = Chemical("ETRO2H      ",  78,  0,  2,   0,  0 ) 
       species( 47) = Chemical("PRRO2H      ",  92,  0,  3,   0,  0 ) 
       species( 48) = Chemical("OXYO2H      ",   1,  0,  8,   0,  0 ) !ds1 
       species( 49) = Chemical("MEKO2H      ", 104,  0,  4,   0,  0 ) 
       species( 50) = Chemical("MALO2H      ", 147,  0,  5,   0,  0 ) 
       species( 51) = Chemical("MVKO2H      ",   1,  0,  4,   0,  0 ) !ds1
       species( 52) = Chemical("MARO2H      ",   1,  0,  5,   0,  0 ) !ds1 
       species( 53) = Chemical("ISRO2H      ",   1,  0,  5,   0,  0 ) !ds1
       species( 54) = Chemical("H2O2        ",  34,  0,  0,   0,  0 ) 
       species( 55) = Chemical("CH3COO2H    ",  76,  0,  2,   0,  0 ) 
       species( 56) = Chemical("CH2CO2HCH3  ",  74,  0,  3,   0,  0 ) 
       species( 57) = Chemical("ISONO3H     ",   1,  0,  5,   0,  0 ) !ds1
       species( 58) = Chemical("ISNIRH      ",   1,  0,  5,   0,  0 ) !ds1
       species( 59) = Chemical("CH3OH       ",  32,  0,  1,   0,  0 ) 
       species( 60) = Chemical("C2H5OH      ",  46,  0,  2,   0,  0 ) 
       species( 61) = Chemical("H2          ",   2,  0,  0,   0,  0 ) 
       species( 62) = Chemical("CO          ",  28,  0,  1,   0,  0 ) 
       species( 63) = Chemical("CH4         ",  16,  0,  1,   0,  0 ) 
       species( 64) = Chemical("SO2         ",  64,  0,  0,   0,  1 ) 
       species( 65) = Chemical("SO4         ",  96,  0,  0,   0,  1 ) 
       species( 66) = Chemical("pNO3        ",  62,  0,  0,   1,  0 ) 
       species( 67) = Chemical("NH3         ",  17,  0,  0,   1,  0 ) 
!hf       species( 68) = Chemical("AMSU        ", 114,  0,  0,   1,  1 ) 
!hf       species( 69) = Chemical("AMNI        ",  80,  0,  0,   2,  0 ) 
       species( 68) = Chemical("aNH4        ", 18,   0,  0,   1,  0 ) 
       species( 69) = Chemical("aNO3        ", 62,   0,  0,   1,  0 ) 
       species( 70) = Chemical("PM25        ", 100,  0,  0,   0,  0 ) 
       species( 71) = Chemical("PMCO        ", 100,  0,  0,   0,  0 ) 
   end subroutine define_chemicals
 end module GenChemicals_ml
 !-----------------------------------------------------------
!>_________________________________________________________<

  module  GenRates_rcmisc_ml
!-----------------------------------------------------------

  
  use PhysicalConstants_ml,  only : PI, RGAS_J
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP
!                                      VOLFAC        ! for N2O5-> NO3-
  use Dates_ml,              only : daynumber !u7.4vg  for so2ox
  use Functions_ml,          only : troe
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !+ Tabulates Rate-coefficients - complex dependancies 

    public :: set_rcmisc_rates

    !u1 integer, parameter, public :: NRCMISC = 7   !! No. coefficients
    integer, parameter, public :: NRCMISC = 18   !! No. coefficients

    real, save, public, dimension(NRCMISC) :: rcvmisc 
    real, save, public, dimension(366)  :: &
      tab_so2ox   ! Tabulated so2->so4 rate for 366 days (leap-year safe!)

  contains
  !------------------------------------
  subroutine set_rcmisc_rates(itemp,tinv,m,o2,h2o,rh,rcmisc) 
  integer, intent(in), dimension(KCHEMTOP:KMAX_MID) :: itemp
  real, intent(in), dimension(KCHEMTOP:KMAX_MID) :: tinv,m,o2,h2o,rh
  real, intent(out),dimension(NRCMISC,KCHEMTOP:KMAX_MID) :: rcmisc
  integer :: k ! local   
  real,  dimension(KCHEMTOP:KMAX_MID) ::  n2   ! nitrogen
  real :: lt3(KCHEMTOP:KMAX_MID)   ! u1 - for Troe
  n2 = m - o2
 
       rcmisc(1,:) = 6.0e-34*m*o2*(300.0*tinv)**2.3 
       rcmisc(2,:) = 1.8e-11*n2*exp(107.0*tinv) 
       rcmisc(3,:) = 3.2e-11*o2*exp(67.0*tinv) 
       rcmisc(4,:) = 2.2e-10*h2o 
       rcmisc(5,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*2.3e-13*exp(600.0*tinv) 
       rcmisc(6,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*1.7e-33*exp(1000.0*tinv)*m 
       rcmisc(7,:) = 1.3e-13*(1+0.6*m/2.55e19) 

!v1 from odin.x3

       do k = KCHEMTOP, KMAX_MID
          if ( rh(k) > 0.4) then
!6sj     !       rcmisc(8,k) = VOLFAC * tab_vav_n2o5( itemp(k) ) * rh(k)
!            rcmisc(8,k) = VOLFAC * rh(k) &
!                * sqrt(3.0 * RGAS_J * itemp(k) / 0.108)  ! m/s !
            rcmisc(8,k) = (2.5 - rh(k)*1.25) & !density, corrected for rh (moderate approx.)
                                               !VOLFAC now in My_Reactions
                * sqrt(3.0 * RGAS_J * itemp(k) / 0.108)  ! m/s !
          else
            rcmisc(8,k) = 0.0
          endif

     if (rh(k) > 0.9 ) then
         rcmisc(10,k) = 1.0e-4
     else
          rcmisc(10,k) = 5.0e-6
     end if
       end do ! k

    !ux - new SO2 -> SO4 method from old MADE/hf code

    !ds rcmisc(9,:) = tab_so2ox(daynumber)


    !u1 - troe stuff put here to simplify .....

  lt3(:) = log(300.0*tinv(:))


  rcmisc(11,:) = troe(1.0e-31*exp(1.6*lt3(:)),3.0e-11*exp(-0.3*lt3(:)), -0.1625,m(:))
  rcmisc(12,:) = troe(2.7e-30*exp(3.4*lt3(:)),2.0e-12*exp(-0.2*lt3(:)),  -1.109,m(:))
  rcmisc(13,:) = troe(1.0e-3*exp(3.5*lt3(:))*exp(-11000*tinv(:)),9.70e14*exp(-0.1*lt3(:))*exp(-11080*tinv(:)),  -1.109,m(:)) 
    rcmisc(14,:) = troe(2.6e-30*exp(2.9*lt3(:)),6.7e-11*exp(0.6*lt3(:)), -0.844,m(:))
    rcmisc(15,:) = troe(2.7e-28*exp(7.1*lt3(:)),1.2e-11*exp(0.1*lt3(:)), -1.204,m(:)) 
    rcmisc(16,:) = 1*troe(4.9e-3*exp(-12100*tinv(:)),5.4e16*exp(-13830*tinv(:)),  -1.204,m(:)) 
    rcmisc(17,:) = troe(7.0e-29*exp(3.1*lt3(:)),9.0e-12, -0.3567,m(:)) 
    rcmisc(18,:) = troe(8.0e-17*exp(3.5*lt3(:)),3.0e-11,-0.6931,m(:)) 

  end subroutine set_rcmisc_rates
end module  GenRates_rcmisc_ml
!>_________________________________________________________<

  module  GenRates_rct_ml
!-----------------------------------------------------------
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP &
                                    , CHEMTMIN, CHEMTMAX !u3
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !+ Tabulates Rate-coefficients - temperature dependant 

!u3    public :: set_rct_rates
    public :: set_rct_rates, set_night_rct

    integer, parameter, public :: NRCT = 37   !! No. coefficients

    real, save, public, dimension(NRCT) :: rcvt 

!u3
!/ Output gas-phase chemical rates:   !u3 - from Tabulations

  real, save, public, &
         dimension(NRCT,CHEMTMIN:CHEMTMAX) :: rcit  ! rate-coefficients

!hf u2!u3 - added for ozone model also
   logical, public, parameter ::  ONLY_NIGHT = .false. ! Only nighttime NO2->NO3
       

  contains
  !------------------------------------
  subroutine set_rct_rates(tinv) 
  real, intent(in) :: tinv
       rcvt(1) = 1.8e-12*exp(-1370.0*tinv) 
       rcvt(2) = 1.2e-13*exp(-2450.0*tinv) 
       rcvt(3) = 1.9e-12*exp(-1000.0*tinv) 
       rcvt(4) = 1.4e-14*exp(-600.0*tinv) 
       rcvt(5) = 1.8e-11*exp(110.0*tinv) 
       rcvt(6) = 3.7e-12*exp(240.0*tinv) 
       rcvt(7) = 7.2e-14*exp(-1414.0*tinv) 
       rcvt(8) = 4.8e-11*exp(250.0*tinv) 
       rcvt(9) = 2.9e-12*exp(-160.0*tinv) 
       rcvt(10) = 7.7e-12*exp(-2100.0*tinv) 
       rcvt(11) = 1.05e-14*exp(785.0*tinv) 
       rcvt(12) = 3.9e-12*exp(-1765.0*tinv) 
       rcvt(13) = 4.2e-12*exp(180.0*tinv) 
       rcvt(14) = 5.9e-14*exp(509.0*tinv) 
       rcvt(15) = 7.04e-14*exp(365.0*tinv) 
       rcvt(16) = 3.1e-12*exp(-360.0*tinv) 
       rcvt(17) = 3.8e-13*exp(780.0*tinv) 
       rcvt(18) = 1e-12*exp(190.0*tinv) 
       rcvt(19) = 1.9e-12*exp(190.0*tinv) 
       rcvt(20) = 8.6e-12*exp(20.0*tinv) 
       rcvt(21) = 7.9e-12*exp(-1030.0*tinv) 
       rcvt(22) = 2.7e-13*exp(1000.0*tinv) 
       rcvt(23) = 5.8e-12*exp(190.0*tinv) 
       rcvt(24) = 5.6e-12*exp(310.0*tinv) 
       rcvt(25) = 2.8e-12*exp(530*tinv) 
       rcvt(26) = 1.3e-13*exp(1040.0*tinv) 
       rcvt(27) = 3e-13*exp(1040.0*tinv) 
       rcvt(28) = 3.69e-12*exp(-70*tinv) 
       rcvt(29) = 1.64e-11*exp(-559.0*tinv) 
       rcvt(30) = 1.2e-14*exp(-2630.0*tinv) 
       rcvt(31) = 6.5e-15*exp(-1880.0*tinv) 
       rcvt(32) = 1.23e-14*exp(-2013*tinv) 
       rcvt(33) = 2.54e-11*exp(410.0*tinv) 
       rcvt(34) = 4.13e-12*exp(452.0*tinv) 
       rcvt(35) = 1.86e-11*exp(175.0*tinv) 
       rcvt(36) = 1.34e+16*exp(-13330.0*tinv) 
       rcvt(37) = 4.32e-15*exp(-2016.0*tinv) 

  end subroutine set_rct_rates
  !------------------------------------------------------
  subroutine set_night_rct(rct,rh,i,j)
  implicit none
  integer,intent(in) :: i,j
  real,intent(in) :: rct(NRCT,KCHEMTOP:KMAX_MID)
  real,intent(in)    :: rh(KCHEMTOP:KMAX_MID)

   ! Dummy for OZONE 

  end subroutine set_night_rct
  !------------------------------------------------------
end module GenRates_rct_ml

!u3 - MOVED HERE  ********
!>_________________________________________________________<

  module  MyChem_ml
!-----------------------------------------------------------
!+u2-u3
! Module containijng initial setup routine calls (Init_mychem)
! and intended to allow the user to specify miscelanneaous
! bits of extra code as needed. Here we have so far included
! Set_2dBgnd in orer to get xn_2d:bgnd for MADE.
! u3 - NEW ********
! We have a new subroutine Init_mychem for all model versions
! which now does tabulations previously done in Tabulations_ml

!u4 - NEW
  use Functions_ml,   only : Daily_sine  ! to specify so2ox
  use GenSpec_bgn_ml,        only : NSPEC_COL ! u3 - nothing more needed
                                   ! for OZONE , xn_2d_bgn, IXBGN_OH, 

 use GenRates_rct_ml, only : &  ! u3
                 NRCT, &         ! No. temperature dependant coefficients
                 rcvt, &         ! Temperature dependant coefficients
                 rcit, &         ! Rate coeffs as rc(n, temp(k) )
                 set_rct_rates   ! Gives RCT as function of temp t

 use GenRates_rcmisc_ml, only : tab_so2ox   !u7.1

!u3   use My_BoundConditions_ml, only : set_daily, h2o2conc
!u3   use Setup_1dfields_ml, only        : amk
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP, KCLOUDTOP &
                                     ,CHEMTMIN, CHEMTMAX  !u3 temp. range
  !u7.4vg                                   ,daynumber  !u3
  use Dates_ml,              only : daynumber !u7.4vg  for so2ox
  use PhysicalConstants_ml,  only : PI, DEG2RAD
  implicit none
  private
                                    !depending on clouds

  public :: Init_mychem          ! Calls model-specific routines
  public :: Set_2dBgnd   ! Sets model-specific background concs.
                         ! (dummy for OZONE so far)


  contains
    !------------------------------------------------------------------

    subroutine  Init_mychem()

    !u3 - moved from Tabulations_ml....
    !+1) Temperature-dependant rates (rct). Only needs to be called once
    !    at beginning of simulations to set up table

      integer :: it          !   Local loop variable
      real    ::  tinv       !   temperature in K

      do it = CHEMTMIN, CHEMTMAX
        tinv = 1.0/real(it)
        call set_rct_rates(tinv)
        rcit(:,it) = rcvt(:)
      end do

     !ACID !+2)
     !ACID ! Tabulate H2O2 values for the full year, with
     !ACID ! a safe 366 value for ndays
     !ACID
     !ACID  tab_h2o2 = Daily_halfsine(0.35e-9,0.3e-9,366)

     !u7.1
     !+2) 
     ! Tabulate SO2 oxidation rates with a safe 366 value for ndays
     ! Coefficients taken from Eliassen+Saltbones (1983) (also in
     ! Berge and Jakobsen, 1998

        !ds tab_so2ox = Daily_sine(3.0e-6,2.0e-6,80,366)
        !ds Changed to use real dmax instead of day of mean in argument

       !rv1.2.1 reset
        !rv1.4.1: tab_so2ox = Daily_sine(3.0e-6,2.5e-6,80+91,366)
        tab_so2ox = Daily_sine(4.0e-6,2.5e-6,80+91,366)
      

    end subroutine  Init_mychem
    !------------------------------------------------------------------

    !u3 - added m as argument:
    subroutine Set_2dBgnd(izen,cloud,m)
      integer, intent(in) :: izen
      real,dimension(KMAX_MID), intent(in) :: cloud ! cloud-cover fraction
      real, intent(in), dimension(KCHEMTOP:KMAX_MID) :: m ! air density

       ! Dummy for OZONE 
    end subroutine Set_2dBgnd

 end module MyChem_ml
 !-----------------------------------------------------------
!>_________________________________________________________<
