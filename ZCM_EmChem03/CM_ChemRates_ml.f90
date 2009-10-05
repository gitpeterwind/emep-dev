!>_________________________________________________________<

  module  ChemRates_rcmisc_ml
!-----------------------------------------------------------

  
  use ChemFunctions_ml       ! => kaero, RiemerN2O5
  use Setup_1dfields_ml      ! => tinv, h2o, m, Fgas
  use Setup_1dfields_ml, m=> amk
  use ChemSpecs_tot_ml         ! => PINALD, .... for FgasJ08
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP,DebugCell,DEBUG_RUNCHEM
  implicit none
  private

  !+ Tabulates Rate-coefficients - complex dependancies 

    public :: set_rcmisc_rates

    integer, parameter, public :: NRCMISC = 18   !! No. coefficients

    real, save, public, dimension(NRCMISC,KCHEMTOP:KMAX_MID) :: rcmisc 

  contains
  !------------------------------------
  subroutine set_rcmisc_rates() 
     real, dimension(KCHEMTOP:KMAX_MID) :: lt300
!     real, dimension(KCHEMTOP:KMAX_MID) :: BranchingNO
!     real, dimension(KCHEMTOP:KMAX_MID) :: BranchingHO2
       lt300(:) = log(300.0*tinv(:))
!       BranchingNO(:) = 1.0e-11*xn_2d(NO,:)/ &
                !( 1.0e-11*xn_2d(NO,:) + 4.2e-12*exp(180*TINV(:))*xn_2d(HO2,:) )
!       BranchingHO2(:) = 1.0 - BranchingNO(:)

       rcmisc(1,:) = 6.0e-34*m*o2*(300.0*tinv)**2.3 
       rcmisc(2,:) = 1.8e-11*n2*exp(107.0*tinv) 
       rcmisc(3,:) = 3.2e-11*o2*exp(67.0*tinv) 
       rcmisc(4,:) = 2.2e-10*h2o 
       rcmisc(5,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*2.3e-13*exp(600.0*tinv) 
       rcmisc(6,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*1.7e-33*exp(1000.0*tinv)*m 
       rcmisc(7,:) = 1.3e-13*(1+0.6*m/2.55e19) 
       rcmisc(8,:) = riemern2o5() 
       rcmisc(9,:) = 1e-12*h2o 
       rcmisc(10,:) = kaero() 
       rcmisc(11,:) = troe(1.0e-31*exp(1.6*lt300),3.0e-11*exp(-0.3*lt300),0.85,m)!a97,j 
       rcmisc(12,:) = troe(2.7e-30*exp(3.4*lt300),2.0e-12*exp(-0.2*lt300),0.33,m)!;a97,j 
       rcmisc(13,:) = troe(1.0e-3*exp(3.5*lt300)*exp(-11000*tinv),9.70e14*exp(-0.1*lt300)*exp(-11080*tinv),0.33,m) 
       rcmisc(14,:) = troe(2.6e-30*exp(2.9*lt300),6.7e-11*exp(0.6*lt300),0.43,m)!;a97,j 
       rcmisc(15,:) = troe(2.7e-28*exp(7.1*lt300),1.2e-11*exp(0.1*lt300),0.3,m) 
       rcmisc(16,:) = troe(4.9e-3*exp(-12100*tinv),5.4e16*exp(-13830*tinv),0.3,m) 
       rcmisc(17,:) = troe(7.0e-29*exp(3.1*lt300),9.0e-12,0.7,m) 
       rcmisc(18,:) = troe(8.0e-17*exp(3.5*lt300),3.0e-11,0.5,m) 

  end subroutine set_rcmisc_rates
end module  ChemRates_rcmisc_ml
!>_________________________________________________________<

  module  ChemRates_rct_ml
!-----------------------------------------------------------

  
  use Setup_1dfields_ml      ! => tinv, h2o, m, Fgas
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP
  implicit none
  private

  !+ Tabulates Rate-coefficients - temperature dependant 

    public :: set_rct_rates

    integer, parameter, public :: NRCT = 37   !! No. coefficients

    real, save, public, dimension(NRCT,KCHEMTOP:KMAX_MID) :: rct 

  contains
  !------------------------------------
  subroutine set_rct_rates() 
       rct(1,:) = 1.8e-12*exp(-1370.0*TINV) 
       rct(2,:) = 1.2e-13*exp(-2450.0*TINV) 
       rct(3,:) = 1.9e-12*exp(-1000.0*TINV) 
       rct(4,:) = 1.4e-14*exp(-600.0*TINV) 
       rct(5,:) = 1.8e-11*exp(110.0*TINV) 
       rct(6,:) = 3.7e-12*exp(240.0*TINV) 
       rct(7,:) = 7.2e-14*exp(-1414.0*TINV) 
       rct(8,:) = 4.8e-11*exp(250.0*TINV) 
       rct(9,:) = 2.9e-12*exp(-160.0*TINV) 
       rct(10,:) = 7.7e-12*exp(-2100.0*TINV) 
       rct(11,:) = 1.05e-14*exp(785.0*TINV) 
       rct(12,:) = 3.9e-12*exp(-1765.0*TINV) 
       rct(13,:) = 4.2e-12*exp(180.0*TINV) 
       rct(14,:) = 5.9e-14*exp(509.0*TINV) 
       rct(15,:) = 7.04e-14*exp(365.0*TINV) 
       rct(16,:) = 3.1e-12*exp(-360.0*TINV) 
       rct(17,:) = 3.8e-13*exp(780.0*TINV) 
       rct(18,:) = 1e-12*exp(190.0*TINV) 
       rct(19,:) = 1.9e-12*exp(190.0*TINV) 
       rct(20,:) = 8.6e-12*exp(20.0*TINV) 
       rct(21,:) = 7.9e-12*exp(-1030.0*TINV) 
       rct(22,:) = 2.7e-13*exp(1000.0*TINV) 
       rct(23,:) = 5.8e-12*exp(190.0*TINV) 
       rct(24,:) = 5.6e-12*exp(310.0*TINV) 
       rct(25,:) = 2.8e-12*exp(530*TINV) 
       rct(26,:) = 1.3e-13*exp(1040.0*TINV) 
       rct(27,:) = 3e-13*exp(1040.0*TINV) 
       rct(28,:) = 3.69e-12*exp(-70*TINV) 
       rct(29,:) = 1.64e-11*exp(-559.0*TINV) 
       rct(30,:) = 1.2e-14*exp(-2630.0*TINV) 
       rct(31,:) = 6.5e-15*exp(-1880.0*TINV) 
       rct(32,:) = 1.23e-14*exp(-2013*TINV) 
       rct(33,:) = 2.54e-11*exp(410.0*TINV) 
       rct(34,:) = 4.13e-12*exp(452.0*TINV) 
       rct(35,:) = 1.86e-11*exp(175.0*TINV) 
       rct(36,:) = 1.34e+16*exp(-13330.0*TINV) 
       rct(37,:) = 4.32e-15*exp(-2016.0*TINV) 

  end subroutine set_rct_rates
end module  ChemRates_rct_ml
