!_____________________________________________________________________________!
 module Setup_1dfields_ml

  ! Arrays of meteorology and concentration for 1-D column , for input to 
  ! chemical solver ........
  ! The k-dimension spans the ground (KMAX_MID) to the K-values
  ! specified by KCHEMTOP - 2 in MACHO.  
  !
  !u2 - new aray added to keep o2, m, and for MADE oh, etc

  use ModelConstants_ml,     only :  KMAX_MID, KCHEMTOP, KUPPER
  use My_Emis_ml,            only :  NRCEMIS, NSS  !SeaS
  use Biogenics_ml,          only :  NBIO     ! Not used yet...
  use GenSpec_tot_ml,        only :  NSPEC_TOT
  use GenSpec_bgn_ml,        only :  NSPEC_COL
  use GenRates_rct_ml,       only :  NRCT
  use GenRates_rcmisc_ml,    only :  NRCMISC
  implicit none
  private


  !/ variables to keep track of which call

   logical, public, save :: first_call = .true.
   integer, public, save :: ncalls     =      0

  !/-- the chemistry is calculated for arrays of size: 

   integer, public, parameter  :: CHEMSIZE = KMAX_MID-KCHEMTOP+1 !  
 
   real, public, dimension(NSPEC_TOT,KCHEMTOP:KMAX_MID), save :: &
                   xn_2d            ! Concentrations [molecules/cm3]  

   real, public, dimension(NRCEMIS,KCHEMTOP:KMAX_MID), save :: rcemis   !emissions
   real, public, dimension(NRCT   ,KCHEMTOP:KMAX_MID), save :: rct    ! T-dependant
   real, public, dimension(NRCMISC,KCHEMTOP:KMAX_MID), save :: rcmisc ! T,M,H2O-dependant
   real, public, dimension(NBIO ,KCHEMTOP:KMAX_MID), save   :: rcbio  !  Biogenic emissions
   real, public, dimension(KCHEMTOP:KMAX_MID), save   :: rc_Rn222  ! 210Pb emissions, ds Pb210
   real, public, dimension(NSS,KCHEMTOP:KMAX_MID),     save :: rcss   ! Sea salt emissions

   real, public, dimension(KCHEMTOP:KMAX_MID), save :: &
          rh                  & ! RH (fraction, 0-1)
         ,amk                 & ! M - atmospheric conc.
         ,temp                & ! temperature
         ,f_Riemer            &  ! weighting factor for N2O5 hydrolysis
         ,pp                     !pressure
   integer, public, dimension(KCHEMTOP:KMAX_MID), save :: &
          itemp                  ! int of temperature



   integer, public, save :: izen           ! integer of zenith angle
                                           !ds rv1.6.3 added here....
   real,    public, save :: Idrctt, Idfuse ! Direct-total and diffuse radiation

 end module Setup_1dfields_ml








