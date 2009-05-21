! <ChemFunctions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module ChemFunctions_ml
!____________________________________________________________________
! Miscellaneous collection of "standard" (or guessed ) functions
! Including Troe, sine and cosine curves, 
! bilinear-interpolation routines, 
! and Standard Atmosphere p -> H conversion
!____________________________________________________________________
!
!** includes
!   troe - standrad chemical function
!
!   Depends on: none - self-contained.
!   Language: F
!   History:
!   ds - 2000-Jan. 2008
!____________________________________________________________________
 use ModelConstants_ml,     only : K1  => KCHEMTOP, K2 => KMAX_MID
 use PhysicalConstants_ml,  only : AVOG, RGAS_J
 use Setup_1dfields_ml,     only : itemp, rh, x=> xn_2d
 use GenSpec_tot_ml,        only : SO4,aNO3,aNH4
  implicit none
  private

  public :: troe
  public :: troeInLog  ! When log(Fc) provided
  public ::  kaero
  public ::  RiemerN2O5


! weighting factor for N2O5 hydrolysis
! Mass of sulfate relative to sulfate+nitrate
! according to  Riemer N, Vogel H, Vogel B, 
! Schell B, Ackermann I, Kessler C, Hass H
! JGR 108 (D4): FEB 27 2003 

  real, parameter, private  :: VOLFACSO4 = 96.0/(AVOG) * 0.90236 *0.02/0.034e-6 
  real, parameter, private  :: VOLFACNO3 = 62.0/(AVOG) * 0.90236 *0.02/0.034e-6 
  real, parameter, private  :: VOLFACNH4 = 18.0/(AVOG) * 0.90236 *0.02/0.034e-6 


  !========================================
  contains
  !========================================

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  elemental function troe(k0,kinf,Fc,M) result (rctroe)

  !+ Calculates Troe expression
  ! -----------------------------------------------------------
  ! ds note - this isn't checked or optimised yet. Taken from
  ! Seinfeld+Pandis, 1998, pp 283, eqn. 5.98. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.

     real, intent(in)  :: k0,kinf,Fc,M
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacament, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf	! used also below
     x    = log10(y)
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

!	give Fc already as log(Fc)
!DSGC - no! Keep as proper Fc. Slower but less confusing

!     rctroe = k0M / ( 1.0 + k0M/kinf) * exp(x*log(Fc))
     rctroe = k0M / ( 1.0 + y) * exp(x*log(Fc))

  end function troe
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  elemental function troeInLog(k0,kinf,LogFc,M) result (rctroe)

  !+ Calculates Troe expression
  ! -----------------------------------------------------------
  ! ds note - this isn't checked or optimised yet. Taken from
  ! Seinfeld+Pandis, 1998, pp 283, eqn. 5.98. 

  ! Input arguments are intended to represent:
  !   M may be O2+N2 or just N2 or just O2.

     real, intent(in)  :: k0,kinf,LogFc,M
     real              :: rctroe

     !-- local
     real :: x,y, K0M               ! temp variable

     k0M = k0 * M
     

     !- use the power function replacament, m**n == exp(n*log m) 
     !-k0M   = a*(T/300.0)**(-2.3) * M
     !-kinf = p*(T/300.0)**(-1.4)

     ! k0M   = a * exp( b*log(t/300.0) ) * M
     ! kinf = p * exp( q*log(t/300.0) )

     ! factors for Fc:
     y    = k0M/kinf	! used also below
     x    = log10(y)
     x    = 1.0/( 1.0 + x*x )

     !- F**x == exp(x*logF)

!	give Fc already as log(Fc)

!     rctroe = k0M / ( 1.0 + k0M/kinf) * exp(x*log(Fc))
     rctroe = k0M / ( 1.0 + y) * exp(x*LogFc)

  end function troeInLog
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


! N2O5 -> nitrate calculation
!===========================================================================
! N2O5 -> nitrate calculation. Some constants for
! calculation of volume fraction of sulphate aerosol, and rate of uptake
! 
!
! The first order reaction coefficient K (corrected for gas phase diffusion, 
! Schwartz, 1986) is given by
!
! K= A* alpha* v/4
!    alpha=sticking coeff. for N2O5 =0.02
!    v=mean molecular speed for N2O5
!    A=aerosol surfac
!
! The surface area of the aerosols can be calculated as
! 
! A = V * surface/volume of aerosols
!     V=volume fraction of sulphate (cm3 aerosol/cm3 air)
!     (similar for nitrate and ammonium):
!
!     e.g.
!     V = (so4 in moleculescm-3) x atw sulphate
!         ---------------------------------------------------------
!        AVOG X specific density of aerosols (assumed 2g/cm3*rh correction)
!
!    Or, shorter, V = S x M0/(AVOG*rho)
!
!    where S is conc. e.g. sulphate (molecule/cm3), M0 is molwt. 
!
!
!    We do not want to include  concentrations  or rho yet, so:
!
!     Let VOL =  M0/AVOG
!   
! The surface/volume ratio is calculated using Whitby particle distribution
! with number mean radius 0.034  and standars deviation (Sigma)=2. 
! Then surface/volume=3/r *  exp( -5/2 *(lnSigma)^2)=26.54 
! 3* exp( -5/2 *(lnSigma)^2)=0.90236
! (monodisperse aerosols; 4*pi*r^2/(4/3 pi*r^3)= 3/r =88.2)
!
! Then 
!      A = VOL * S * 0.90236 /(0.034e-6*rho) 
! and
!      K = VOL * S * 0.90236 /(0.034e-6*rho)    * alpha* v/4
! Set
!      VOLFAC= VOL*0.90236/0.034e-6 *alpha    
! Then
!      K = VOLFAC *S *v/(4*rho)
!
! rcmisc(8,k) (in My_Chem)=v/(4*rho) 
!
!      K = VOLFAC *rcmisc(8,k) *S
! According to Riemer et al, 2003, we weight the reaction probability
! according to the composition of the aerosol
!
! alpha(N2O5)=f*alpha1 +(1-f)alpha2
!   alpha1=0.02
!   alpha2=0.002
!   f= Mso4/(Mso4+Mno3), M=aerosol mass concentration
 
! N2O5 -> aerosol based upon  based on Riemer 2003
! J08 - changed from use of more accurate xnew to xn_2d, since
! surface area won't change so much, and anyway the uncertainties
! are large. (and xn_2d leads to fewer dependencies)

  function RiemerN2O5() result(rate) 
     real, dimension(K1:K2) :: rate
     real    :: rc   ! was rcmisc(8,)
     real    :: f   ! Was f_Riemer
     integer :: k
     
! old Setup_ml had:
! setup weighting factor for hydrolysis  
!     f_Riemer(k)=96.*xn_2d(SO4,k)/( (96.*xn_2d(SO4,k))+(62.*xn_2d(aNO3,k)) )
!Query - no aNH4 above?
! then FastReactions had for L(N2O5):
!     + (0.9*f_Riemer(k)+0.1) * rcmisc(8,k)* &
!                 ( VOLFACSO4*xnew(SO4)      & !Total sulpate aerosol surface
!                 + VOLFACNO3*xnew(aNO3)     & !Total sulpate aerosol surface
!                 + VOLFACNH4*xnew(aNH4)  )  & !Total sulpate aerosol surface
   

     do k = K1, K2
       if ( rh(k)  > 0.4) then

          rc = sqrt(3.0 * RGAS_J * itemp(k) / 0.108) & ! mean mol. speed,m/s
             /(4*(2.5 - rh(k)*1.25)) !density, corrected for rh (moderate approx.)

          f = 96.0*x(SO4,k)/( 96.*x(SO4,k) + 62.0*x(aNO3,k) )


          rate(k) =  (0.9*f + 0.1) * rc *  &
             ( VOLFACSO4 * x(SO4,k) + VOLFACNO3 * x(aNO3,k) &
              + VOLFACNH4 * x(aNH4,k) )    !Total aerosol surface
        else
          rate(k) = 0.0
        endif
    end do ! k

  end function RiemerN2O5
  !---------------------------------------------------------------------
  function kaero() result(rate) 
     real, dimension(K1:K2) :: rate
     integer :: k
     
    do k = K1, K2
      if ( rh(k)  > 0.9) then
         rate(k) = 1.0e-4
      else
         rate(k) = 5.0e-6
      end if
    end do !k


  end function kaero
end module ChemFunctions_ml
