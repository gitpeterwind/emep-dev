! <My_ChemRates_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
!* 
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*  
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!* 
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!* 
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************! 
!>_________________________________________________________<

!>_________________________________________________________<

  module  GenRates_rcmisc_ml
!-----------------------------------------------------------

  
  use ChemFunctions_ml, only :  troeInLog, kaero, RiemerN2O5
  use Setup_1dfields_ml      ! => tinv, h2o, m, Fgas
  use Setup_1dfields_ml, m=> amk
  use GenSpec_tot_ml         ! => PINALD, .... for Fgas/SOA
  use PhysicalConstants_ml,  only : PI, RGAS_J
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP,DebugCell,DEBUG_RUNCHEM
!                                      VOLFAC        ! for N2O5-> NO3-

  implicit none
  private
!/ ...... ..   ( from GenChem )

  !+ Tabulates Rate-coefficients - complex dependancies 

    public :: set_rcmisc_rates

    integer, parameter, public :: NRCMISC = 18   !! No. coefficients

    real, save, public, dimension(NRCMISC,KCHEMTOP:KMAX_MID) :: rcmisc 
!DSGC    real, save, public, dimension(NRCMISC) :: rcvmisc 
!dsx   real, save, public, dimension(366)  :: &
!dsx  tab_so2ox   ! Tabulated so2->so4 rate for 366 days (leap-year safe!)

  contains
  !------------------------------------
  subroutine set_rcmisc_rates() 
!DSGC  subroutine set_rcmisc_rates(itemp,tinv,m,o2,h2o,rh,rcmisc) 
!DSGC  integer, intent(in), dimension(KCHEMTOP:KMAX_MID) :: itemp
!DSGC  real, intent(in), dimension(KCHEMTOP:KMAX_MID) :: tinv,m,o2,h2o,rh
!DSGC  real, intent(out),dimension(NRCMISC,KCHEMTOP:KMAX_MID) :: rcmisc
     real, dimension(KCHEMTOP:KMAX_MID) :: lt300
     real, dimension(KCHEMTOP:KMAX_MID) :: BranchingNO
     real, dimension(KCHEMTOP:KMAX_MID) :: BranchingHO2
       lt300(:) = log(300.0*tinv(:))
!DSGC NOT YET:
!       BranchingNO(:) = 1.0e-11*xn_2d(NO,:)/ &
!                ( 1.0e-11*xn_2d(NO,:) + 4.2e-12*exp(180*TINV(:))*xn_2d(HO2,:) )
!       BranchingHO2(:) = 1.0 - BranchingNO(:)
!

       rcmisc(1,:) = 6.0e-34*m*o2*(300.0*tinv)**2.3 
       rcmisc(2,:) = 1.8e-11*n2*exp(107.0*tinv) 
       rcmisc(3,:) = 3.2e-11*o2*exp(67.0*tinv) 
       rcmisc(4,:) = 2.2e-10*h2o 
       rcmisc(5,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*2.3e-13*exp(600.0*tinv) 
       rcmisc(6,:) = (1.0+1.4e-21*h2o*exp(2200.0*tinv))*1.7e-33*exp(1000.0*tinv)*m 
       rcmisc(7,:) = 1.3e-13*(1+0.6*m/2.55e19) 
       rcmisc(8,:) = RiemerN2O5()
       rcmisc(10,:) = kaero()

!DSGC       do k = KCHEMTOP, KMAX_MID
!DSGC          if ( rh(k) > 0.4) then
!DSGC            rcmisc(8,k) = &
!DSGC                sqrt(3.0 * RGAS_J * itemp(k) / 0.108) & ! mean molecular speed,m/s !
!DSGC            /(4*(2.5 - rh(k)*1.25))!density, corrected for rh (moderate approx.)
!DSGC                                   !VOLFAC now in My_Reactions
!DSGC          else
!DSGC            rcmisc(8,k) = 0.0
!DSGC          endif
!DSGC
!DSGC     if (rh(k) > 0.9 ) then
!DSGC         rcmisc(10,k) = 1.0e-4
!DSGC     else
!DSGC          rcmisc(10,k) = 5.0e-6
!DSGC     end if
!DSGC       end do ! k


  rcmisc(11,:) = troeInLog(1.0e-31*exp(1.6*lt300(:)),3.0e-11*exp(-0.3*lt300(:)), -0.1625,m(:))
  rcmisc(12,:) = troeInLog(2.7e-30*exp(3.4*lt300(:)),2.0e-12*exp(-0.2*lt300(:)),  -1.109,m(:))
  rcmisc(13,:) = troeInLog(1.0e-3*exp(3.5*lt300(:))*exp(-11000*tinv(:)),&
        9.70e14*exp(-0.1*lt300(:))*exp(-11080*tinv(:)),  -1.109,m(:)) 
    rcmisc(14,:) = troeInLog(2.6e-30*exp(2.9*lt300(:)),6.7e-11*exp(0.6*lt300(:)), -0.844,m(:))
    rcmisc(15,:) = troeInLog(2.7e-28*exp(7.1*lt300(:)),1.2e-11*exp(0.1*lt300(:)), -1.204,m(:)) 
    rcmisc(16,:) = troeInLog(4.9e-3*exp(-12100*tinv(:)),5.4e16*exp(-13830*tinv(:)),  -1.204,m(:)) 
    rcmisc(17,:) = troeInLog(7.0e-29*exp(3.1*lt300(:)),9.0e-12, -0.3567,m(:)) 
    rcmisc(18,:) = troeInLog(8.0e-17*exp(3.5*lt300(:)),3.0e-11,-0.6931,m(:)) 

  end subroutine set_rcmisc_rates
end module  GenRates_rcmisc_ml
!>_________________________________________________________<

  module  GenRates_rct_ml
!-----------------------------------------------------------

  use Setup_1dfields_ml      ! => tinv, h2o, m, Fgas
  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP

!DSGC  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP &
!DSGC                                    , CHEMTMIN, CHEMTMAX !u3
  implicit none
  private
!/ ...... ..   ( from GenChem )


  !+ Tabulates Rate-coefficients - temperature dependant 

    public :: set_rct_rates   !DSGC, set_night_rct

    integer, parameter, public :: NRCT = 37   !! No. coefficients

!DSGC    real, save, public, dimension(NRCT) :: rcvt 
    real, save, public, dimension(NRCT,KCHEMTOP:KMAX_MID) :: rct

!/ Output gas-phase chemical rates:   ! - from Tabulations

!DSGC  real, save, public, &
!DSGC         dimension(NRCT,CHEMTMIN:CHEMTMAX) :: rcit  ! rate-coefficients

!- added for ozone model also (consistency with ACID)
! Only nighttime NO2->NO3
!DSGC   logical, public, parameter ::  ONLY_NIGHT = .false. 
       

  contains
  !------------------------------------
  subroutine set_rct_rates() 
!DSGC  real, intent(in) :: tinv
       rct(1,:) = 1.8e-12*exp(-1370.0*tinv) 
       rct(2,:) = 1.2e-13*exp(-2450.0*tinv) 
       rct(3,:) = 1.9e-12*exp(-1000.0*tinv) 
       rct(4,:) = 1.4e-14*exp(-600.0*tinv) 
       rct(5,:) = 1.8e-11*exp(110.0*tinv) 
       rct(6,:) = 3.7e-12*exp(240.0*tinv) 
       rct(7,:) = 7.2e-14*exp(-1414.0*tinv) 
       rct(8,:) = 4.8e-11*exp(250.0*tinv) 
       rct(9,:) = 2.9e-12*exp(-160.0*tinv) 
       rct(10,:) = 7.7e-12*exp(-2100.0*tinv) 
       rct(11,:) = 1.05e-14*exp(785.0*tinv) 
       rct(12,:) = 3.9e-12*exp(-1765.0*tinv) 
       rct(13,:) = 4.2e-12*exp(180.0*tinv) 
       rct(14,:) = 5.9e-14*exp(509.0*tinv) 
       rct(15,:) = 7.04e-14*exp(365.0*tinv) 
       rct(16,:) = 3.1e-12*exp(-360.0*tinv) 
       rct(17,:) = 3.8e-13*exp(780.0*tinv) 
       rct(18,:) = 1e-12*exp(190.0*tinv) 
       rct(19,:) = 1.9e-12*exp(190.0*tinv) 
       rct(20,:) = 8.6e-12*exp(20.0*tinv) 
       rct(21,:) = 7.9e-12*exp(-1030.0*tinv) 
       rct(22,:) = 2.7e-13*exp(1000.0*tinv) 
       rct(23,:) = 5.8e-12*exp(190.0*tinv) 
       rct(24,:) = 5.6e-12*exp(310.0*tinv) 
       rct(25,:) = 2.8e-12*exp(530*tinv) 
       rct(26,:) = 1.3e-13*exp(1040.0*tinv) 
       rct(27,:) = 3e-13*exp(1040.0*tinv) 
       rct(28,:) = 3.69e-12*exp(-70*tinv) 
       rct(29,:) = 1.64e-11*exp(-559.0*tinv) 
       rct(30,:) = 1.2e-14*exp(-2630.0*tinv) 
       rct(31,:) = 6.5e-15*exp(-1880.0*tinv) 
       rct(32,:) = 1.23e-14*exp(-2013*tinv) 
       rct(33,:) = 2.54e-11*exp(410.0*tinv) 
       rct(34,:) = 4.13e-12*exp(452.0*tinv) 
       rct(35,:) = 1.86e-11*exp(175.0*tinv) 
       rct(36,:) = 1.34e+16*exp(-13330.0*tinv) 
       rct(37,:) = 4.32e-15*exp(-2016.0*tinv) 

  end subroutine set_rct_rates
  !------------------------------------------------------
!DSGC  subroutine set_night_rct(rct,rh,i,j)
!DSGC  implicit none
!DSGC  integer,intent(in) :: i,j
!DSGC  real,intent(in) :: rct(NRCT,KCHEMTOP:KMAX_MID)
!DSGC  real,intent(in)    :: rh(KCHEMTOP:KMAX_MID)
!DSGC
!DSGC   ! Dummy for OZONE 
!DSGC
!DSGC  end subroutine set_night_rct
  !------------------------------------------------------
end module GenRates_rct_ml

!>_________________________________________________________<

!DSGC  module  MyChem_ml
!DSGC!-----------------------------------------------------------
!DSGC! Module containijng initial setup routine calls (Init_mychem)
!DSGC! and intended to allow the user to specify miscelanneaous
!DSGC! bits of extra code as needed. Here we have so far included
!DSGC! Set_2dBgnd in orer to get xn_2d:bgnd for MADE.
!DSGC! 
!DSGC! We have a new subroutine Init_mychem for all model versions
!DSGC! which now does tabulations previously done in Tabulations_ml
!DSGC
!DSGC!dsx  use Functions_ml,   only : Daily_sine  ! to specify so2ox
!DSGC  use GenSpec_bgn_ml,        only : NSPEC_COL ! - nothing more needed
!DSGC                                   ! for OZONE , xn_2d_bgn, IXBGN_OH, 
!DSGC
!DSGC use GenRates_rct_ml, only : &  !
!DSGC                 NRCT, &         ! No. temperature dependant coefficients
!DSGC                 rcvt, &         ! Temperature dependant coefficients
!DSGC                 rcit, &         ! Rate coeffs as rc(n, temp(k) )
!DSGC                 set_rct_rates   ! Gives RCT as function of temp t
!DSGC
!DSGC!dsx use GenRates_rcmisc_ml, only : tab_so2ox
!DSGC
!DSGC  use ModelConstants_ml,     only : KMAX_MID,KCHEMTOP, KCLOUDTOP &
!DSGC                                     ,CHEMTMIN, CHEMTMAX  !u3 temp. range
!DSGC  use PhysicalConstants_ml,  only : PI, DEG2RAD
!DSGC  implicit none
!DSGC  private
!DSGC                                    !depending on clouds
!DSGC
!DSGC  public :: Init_mychem          ! Calls model-specific routines
!DSGC  public :: Set_2dBgnd   ! Sets model-specific background concs.
!DSGC                         ! (dummy for OZONE so far)
!DSGC
!DSGC
!DSGC  contains
!DSGC    !------------------------------------------------------------------
!DSGC
!DSGC    subroutine  Init_mychem()
!DSGC
!DSGC    !+1) Temperature-dependant rates (rct). Only needs to be called once
!DSGC    !    at beginning of simulations to set up table
!DSGC
!DSGC      integer :: it          !   Local loop variable
!DSGC      real    ::  tinv       !   temperature in K
!DSGC
!DSGC      do it = CHEMTMIN, CHEMTMAX
!DSGC        tinv = 1.0/real(it)
!DSGC        call set_rct_rates(tinv)
!DSGC        rcit(:,it) = rcvt(:)
!DSGC      end do
!DSGC
!DSGC     !+2) 
!DSGC     ! Tabulate SO2 oxidation rates with a safe 366 value for ndays
!DSGC     ! Coefficients taken from Eliassen+Saltbones (1983) (also in
!DSGC     ! Berge and Jakobsen, 1998
!DSGC
!DSGC!dsx        tab_so2ox = Daily_sine(4.0e-6,2.5e-6,80+91,366)
!DSGC      
!DSGC
!DSGC    end subroutine  Init_mychem
!DSGC    !------------------------------------------------------------------
!DSGC
!DSGC    subroutine Set_2dBgnd(izen,cloud,m)
!DSGC      integer, intent(in) :: izen
!DSGC      real,dimension(KMAX_MID), intent(in) :: cloud ! cloud-cover fraction
!DSGC      real, intent(in), dimension(KCHEMTOP:KMAX_MID) :: m ! air density
!DSGC
!DSGC       ! Dummy for OZONE 
!DSGC    end subroutine Set_2dBgnd
!DSGC
!DSGC end module MyChem_ml
!DSGC !-----------------------------------------------------------
!>_________________________________________________________<
