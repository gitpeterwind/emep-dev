!######### ACID model ######################################
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
 
   integer, public, parameter ::  NSPEC_BGN = 2 ! No. 3D bgn species
   integer, public, parameter ::  NSPEC_COL = 4 ! total no. prescribed specs


!**/ First species from global CTM
   integer, public, parameter ::   & 
     IXBGN_O3            =   1     

!**/ THEN (!!) species prescribed by zenith angle 
   integer, public, parameter ::   & 
     IXBGN_OH            =   2& !3     &
    ,IXBGN_CH3COO2       =   3       
!**/ THEN (!!) simple function reset daily
!u7.1    ,IXBGN_H2O2           =   4   

   real, public, save, dimension(NSPEC_COL,KCHEMTOP:KMAX_MID) :: xn_2d_bgn


!-----------------------------------------------------------
  end module GenSpec_bgn_ml

!>_________________________________________________________<

  module  GenSpec_adv_ml
!-----------------------------------------------------------
  implicit none
!/ ...... ..   ( from GenChem )

!+ Defines indices and NSPEC for adv : Advected species

!   ( Output from GenChem, sub print_species ) 

   integer, public, parameter ::  NSPEC_ADV = 9
 
 ! Aerosols:

   logical, public, parameter ::  ORG_AEROSOLS = .false. 

   integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  



   integer, public, parameter ::   & 
     IXADV_NO          =   1   &
  ,  IXADV_NO2         =   2   &
  ,  IXADV_PAN         =   3   &
  ,  IXADV_HNO3        =   4   &
  ,  IXADV_SO2         =   5   &
  ,  IXADV_SO4         =   6   &
  ,  IXADV_NH3         =   7   &
  ,  IXADV_AMSU        =   8   &
  ,  IXADV_AMNI        =   9

 !-----------------------------------------------------------
  end module GenSpec_adv_ml
!>_________________________________________________________<

  module  GenSpec_shl_ml
!-----------------------------------------------------------
  implicit none
!/ ...... ..   ( from GenChem )

!+ Defines indices and NSPEC for shl : Short-lived (non-advected) species 

!   ( Output from GenChem, sub print_species ) 

   logical, public, parameter ::  ORG_AEROSOLS = .false. 

   integer, public, parameter ::  NSPEC_SHL = 0 
 
 ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  


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

   integer, public, parameter ::  NSPEC_TOT = 9 
 
 ! Aerosols:
           integer, public, parameter :: &
                NAEROSOL  = 0,         &!   Number of aerosol species
                FIRST_SOA = -99,       &!   First aerosol species
                LAST_SOA  = -99         !   Last  aerosol species  



   integer, public, parameter ::   & 
     NO          =   1   &
  ,  NO2         =   2   &
  ,  PAN         =   3   &
  ,  HNO3        =   4   &
  ,  SO2         =   5   &
  ,  SO4         =   6   &
  ,  NH3         =   7   &
  ,  AMSU        =   8   &
  ,  AMNI        =   9

 !-----------------------------------------------------------
  end module GenSpec_tot_ml
!>_________________________________________________________<

  module  GenChemicals_ml
!-----------------------------------------------------------

   use GenSpec_tot_ml, only: NSPEC_TOT ! Total number of species for chemistry
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
       species( 1) = Chemical("NO          ",  30,  0,  0,   1,  0 ) 
       species( 2) = Chemical("NO2         ",  46,  0,  0,   1,  0 ) 
       species( 3) = Chemical("PAN         ", 121,  0,  2,   1,  0 ) 
       species( 4) = Chemical("HNO3        ",  63,  0,  0,   1,  0 ) 
       species( 5) = Chemical("SO2         ",  64,  0,  0,   0,  1 ) 
       species( 6) = Chemical("SO4         ",  96,  0,  0,   0,  1 ) 
       species( 7) = Chemical("NH3         ",  17,  0,  0,   1,  0 ) 
       species( 8) = Chemical("AMSU        ", 114,  0,  0,   1,  1 ) 
       species( 9) = Chemical("AMNI        ",  80,  0,  0,   2,  0 ) 
   end subroutine define_chemicals
 end module GenChemicals_ml
 !-----------------------------------------------------------
!>_________________________________________________________<

  module  GenRates_rcmisc_ml
!-----------------------------------------------------------

  
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP,&
                                      daynumber, & ! for so2ox
                                      VOLFAC ! for N2O5-> NO3-
  use Functions_ml,          only : troe
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !+ Tabulates Rate-coefficients - complex dependancies 

    public :: set_rcmisc_rates

    integer, parameter, public :: NRCMISC = 18   !! No. coefficients
!hf made For now I keep them all
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

    !ux - new SO2 -> SO4 method from old MADE/hf code

      rcmisc(9,:) = tab_so2ox(daynumber)

    !u1 - troe stuff put here to simplify .....

  lt3(:) = log(300.0*tinv(:))


  rcmisc(11,:) = troe(1.0e-31*exp(1.6*lt3(:)),3.0e-11*exp(-0.3*lt3(:)), -0.1625,m(:))
  rcmisc(12,:) = troe(2.7e-30*exp(3.4*lt3(:)),2.0e-12*exp(-0.2*lt3(:)),  -1.109,m(:))
  rcmisc(13,:) = troe(1.0e-3*exp(3.5*lt3(:))*exp(-11000*tinv(:)),9.70e14*exp(-0.1*lt3(:))*exp(-11080*tinv(:)),  -1.109,m(:)) 
    rcmisc(14,:) = troe(2.6e-30*exp(2.9*lt3(:)),6.7e-11*exp(0.6*lt3(:)), -0.844,m(:))
    rcmisc(15,:) = troe(2.7e-28*exp(7.1*lt3(:)),1.2e-11*exp(0.1*lt3(:)), -1.204,m(:)) 
    rcmisc(16,:) = troe(4.9e-3*exp(-12100*tinv(:)),5.4e16*exp(-13830*tinv(:)),  -1.204,m(:)) 
    rcmisc(17,:) = troe(7.0e-29*exp(3.1*lt3(:)),9.0e-12, -0.3567,m(:)) 
    rcmisc(18,:) = troe(8.0e-17*exp(3.5*lt3(:)),3.0e-11,-0.6931,m(:)) 

  end subroutine set_rcmisc_rates
end module  GenRates_rcmisc_ml
!<
!>_________________________________________________________<

  module  GenRates_rct_ml
!-----------------------------------------------------------
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP &
                                   ,CHEMTMIN, CHEMTMAX  ! u3
!hf u2:
  use Radiation_ml,           only : zen
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !+ Tabulates Rate-coefficients - temperature dependant 

!hf u2    public :: set_rct_rates
          public :: set_rct_rates, set_night_rct

    integer, parameter, public :: NRCT = 2   !! No. coefficients

    real, save, public, dimension(NRCT) :: rcvt 

!/ Output gas-phase chemical rates:   !u3 - from Tabulations

  real, save, public, &
         dimension(NRCT,CHEMTMIN:CHEMTMAX) :: rcit  ! rate-coefficients

!hf u2
   logical, public, parameter ::  ONLY_NIGHT = .true. ! Only nighttime NO2->NO3
                                                      ! reaction
  contains
  !------------------------------------
  subroutine set_rct_rates(tinv) 
  real, intent(in) :: tinv
!hf Need only 1 and 2 for MADE like version
       rcvt(1) = 1.8e-12*exp(-1370.0*tinv) 
       rcvt(2) = 1.2e-13*exp(-2450.0*tinv) 
  end subroutine set_rct_rates
  !----------------------------------------
  subroutine set_night_rct(rct,rh,i,j)
  implicit none
  real,intent(inout) :: rct(NRCT,KCHEMTOP:KMAX_MID)
  real,intent(in)    :: rh(KCHEMTOP:KMAX_MID)
  integer,intent(in) :: i,j
  integer            ::  k

  !**/ Only NO2+O3->H+ +NO3- at night time and in the
  !    8 lowest layers and if rh>0.5
  ! rct(2,13:KMAX_MID) already defined

    if ( zen(i,j)<90. ) rct(2,:)=0. !Day 
    if ( zen(i,j)>90. ) then !Night
       rct(2,KCHEMTOP:12) =0.!Too high up
       do k=KCHEMTOP,KMAX_MID
          if ( rh(k)<0.5 ) rct(2,k)= 0.!Too low rel. hum.
       enddo
    endif !night
  end subroutine set_night_rct
  !----------------------------------------
end module  GenRates_rct_ml
!>_________________________________________________________<
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

!u3 - NEW
  use Functions_ml,   only : Daily_sine  ! to specify so2ox
  use GenSpec_bgn_ml,        only : NSPEC_COL, xn_2d_bgn, &  
                                    IXBGN_OH, IXBGN_CH3COO2  !u7.1  & !u3
                                   !u7.1 ,IXBGN_H2O2
 use GenRates_rct_ml, only : &  ! u3
                 NRCT, &         ! No. temperature dependant coefficients
                 rcvt, &         ! Temperature dependant coefficients
                 rcit, &         ! Rate coeffs as rc(n, temp(k) )
                 set_rct_rates   ! Gives RCT as function of temp t

 use GenRates_rcmisc_ml, only : tab_so2ox    ! u7.1

!hf h2o2
!u3   use My_BoundConditions_ml, only : set_daily, h2o2conc
!u3   use Setup_1dfields_ml, only        : amk
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP, KCLOUDTOP &
                                     ,CHEMTMIN, CHEMTMAX & !u3 temp. range
                                     ,daynumber  !u3
  use PhysicalConstants_ml,  only : PI, DEG2RAD
  implicit none
  private
                                    !depending on clouds

  public :: Init_mychem          ! Calls model-specific routines
  public :: Set_2dBgnd   ! Sets model-specific background concs.

  !u7.1 real, public, dimension(366), save :: &
  !u7.1     tab_h2o2   ! Tabulated h2o2 (ppb) for 366 days (leap-year safe!)


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

     !+2)
     ! Tabulate H2O2 values for the full year, with
     ! a safe 366 value for ndays

     !u7.1   tab_h2o2 = Daily_halfsine(0.35e-9,0.3e-9,366)

     !u7.1
     !+2) 
     ! Tabulate SO2 oxidation rates with a safe 366 value for ndays
     ! Coefficients taken from Eliassen+Saltbones (1983) (also in
     ! Berge and Jakobsen, 1998

        tab_so2ox = Daily_sine(4.0e-6,2.5e-6,80,366)

    end subroutine  Init_mychem
    !------------------------------------------------------------------

    !u3 - added m as argument:
    subroutine Set_2dBgnd(izen,cloud,m)
      integer, intent(in) :: izen
      real,dimension(KMAX_MID), intent(in) :: cloud ! cloud-cover fraction
      real, intent(in), dimension(KCHEMTOP:KMAX_MID) :: m ! air density

!hf MADE
  real    :: diur       ! zenith angle dependency of bgn conc.
  real    :: radzen     ! zenith angle in radians
  real, dimension(KCHEMTOP:KMAX_MID) &
          :: reduce_fac ! Factor to reduce oh and ch3coo
     
 

   !  Background species prescribed from solar zenith angle
   !  or simple function

     !u2 if (k >  8) then
     !/ reduce conc below layer 8

     reduce_fac(KCHEMTOP:KCLOUDTOP)   = 1.0
     reduce_fac(KCLOUDTOP+1:KMAX_MID) = 1.0-0.5*cloud(KCLOUDTOP+1:KMAX_MID)

     radzen= izen * DEG2RAD        !convert to radians

     if( izen < 90  ) then !daytime
        diur=exp( -0.25/cos(radzen) )
        xn_2d_bgn(IXBGN_OH,:)      = ( 1.e4 + 4.e6*diur)*reduce_fac
        xn_2d_bgn(IXBGN_CH3COO2,:) = ( 0.5*(1.e6 + 5.e6*diur) )*reduce_fac
     else
        xn_2d_bgn(IXBGN_OH,:)      = 1.e4
        xn_2d_bgn(IXBGN_CH3COO2,:) = 1.e6
     endif
 
     !u3 call set_daily
     !u7.1  xn_2d_bgn(IXBGN_H2O2,:) = tab_h2o2(daynumber) * m(:) 
    end subroutine Set_2dBgnd

 end module MyChem_ml
 !-----------------------------------------------------------
!>_________________________________________________________<
