! <Functions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Functions_ml
!____________________________________________________________________
! Miscellaneous collection of "standard" (or guessed ) functions
! Including Troe, sine and cosine curves, 
! and Standard Atmosphere p -> H conversion
!____________________________________________________________________
!
!** includes
!
!   Depends on: none - self-contained.
!   Language: F
!   History:
!   ds - 2000-Jan. 2001
!____________________________________________________________________
  use PhysicalConstants_ml, only : KAPPA, PI, DEG2RAD
  implicit none
  private

  public ::  Daily_cosine   ! Generates daily values of a variable
                            ! specified as a cosine curve over the year.
  public ::  Daily_sine     ! Generates daily values of a variable
                            ! specified as a sine curve over the year.
  public ::  Daily_halfsine ! Similar, but only half-sine curve (0..pi)
                            ! used. (E.g. for H2O2 in ACID versions)

  public :: StandardAtmos_kPa_2_km   ! US Standard Atmosphere conversion
  public :: StandardAtmos_km_2_kPa   ! US Standard Atmosphere conversion

!  public :: inside_1234 !test wether a point is inside the quadrilateral 1234 

  public :: great_circle_distance!distance between two points following the surface on a unit sphere

 !/- Exner subroutines: ------------------------------------------------------

 public :: Exner_nd        ! (p/P0)**KAPPA
 public :: Tpot_2_T        ! Same as Exner_nd - but easier to remember
 public :: T_2_Tpot        ! Inverse as Exner_nd
 public :: Exner_tab       ! Tabulation. Must be called first


 !/- Interpolation constants

 real, private, parameter  ::    &
       PINC=1000.0              &  
      ,P0  =1.0e5               &    ! Standard pressure
      ,PBAS=-PINC

 real, save, private, dimension(131) ::  tab_exf  ! Tabulated Exner



  !========================================
  contains
  !========================================

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function Daily_cosine(mean, amp, dmax, ndays) result (daily)
  !+
  !   Specifies cosine curve for a variable over a year

     real,    intent(in)  :: mean, amp    ! Annual mean and amplitude of sine
     integer, intent(in)  :: dmax         ! Day where maximum occurs
     integer, intent(in)  :: ndays        ! No. days per year   (365/366)

     real, dimension(ndays) :: daily
     integer    :: d           
     real, save :: twopi                  ! Could use PhysiclConstants_ml
     twopi = 8.0 * atan(1.0)              ! but I prefer to keep Functions_ml
                                          ! standalone

     do d = 1, ndays
      daily(d) = mean + amp * cos ( twopi * (d - dmax)/ ndays )
     end do

  end function Daily_cosine
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function Daily_sine(mean, amp, dmax, ndays) result (daily)
  !+
  !   Specifies sine curve for a variable over a year
  !   25/9/2002, ds, dmax redifined to be true dmax. Before it was
  !   80 and actually the day when the mean ocrrurred (spotted by hf)

     real,    intent(in)  :: mean, amp    ! Annual mean and amplitude of sine
     integer, intent(in)  :: dmax         ! Day where maximum occurs
     integer, intent(in)  :: ndays        ! No. days per year   (365/366)

     real, dimension(ndays) :: daily
     integer    :: d           
     real, save :: shift                  ! Shifts sine curve to give max 
                                          ! when d = dmax
     real, save :: twopi                  ! Could use PhysiclConstants_ml
     twopi = 8.0 * atan(1.0)              ! but I prefer to keep Functions_ml
                                          ! standalone
     shift = ndays/4.0

     do d = 1, ndays
      daily(d) = mean + amp * sin ( twopi * (d + shift - dmax)/ ndays )
     end do

  end function Daily_sine
  
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  function Daily_halfsine(base, amp, ndays) result (daily)
  !+
  !   Specifies half-sine curve for a variable over a year, with
  !   values 1.0 at start and end, and max in mid-summer.

     real,    intent(in)  :: base, amp    ! Annual base and amplitude of sine
     integer, intent(in)  :: ndays        ! No. days per year   (365/366)

     real, dimension(ndays) :: daily
     integer    :: d           
     real, save :: pi                     ! Could use PhysiclConstants_ml
     pi = 4.0 * atan(1.0)                 ! but I prefer to keep Functions_ml
                                          ! standalone

     do d = 1, ndays
      daily(d) = base + amp * sin ( pi * (ndays - d )/ ndays )
     end do

  end function Daily_halfsine

 !=======================================================================

 elemental function StandardAtmos_km_2_kPa(h_km) result (p_kPa)
 !=======================================================================
   implicit none

  !+ Converts height (km)  to pressure (kPa) for a US standard Atmosphere
  !  Valid up to 20 km
  !
  ! pw 07/4/2010

   real :: p_kPa
   real   , intent(in)          :: h_km
   real :: t    ! Temperature (K)

   if( h_km < 11.0 ) then   ! = p_kPa > 22.632
      ! t = 288.15/(p_kPa/101.325)**(-1.0/5.255876)
      !- use the power function replacament, m**n == exp(n*log m)
      p_kPa = 101.325*exp(-5.255876*log(288.15/(288.15-6.5*h_km)))
   else
      p_kPa =  22.632*exp(-0.1576884*(h_km - 11.0)  )
      
   end if
   
 end function StandardAtmos_km_2_kPa

 !=======================================================================
 elemental function StandardAtmos_kPa_2_km(p_kPa) result (h_km)
 !=======================================================================
   implicit none

  !+ Converts pressure (kPa)  to height (km) for a US standard Atmosphere
  !  Valid up to 20 km
  !
  ! ds 27/7/2003

   real, intent(in) :: p_kPa
   real             :: h_km
   real :: t    ! Temperature (K)

   if( p_kPa > 22.632 ) then   ! = 11 km height
         ! t = 288.15/(p_kPa/101.325)**(-1.0/5.255876)
         !- use the power function replacament, m**n == exp(n*log m)

         t = 288.15/exp(-1.0/5.255876*log(p_kPa/101.325))
         h_km = (288.15-t)/6.5
   else
         h_km = 11.0 + log( p_kPa/22.632)/(-0.1576884)
   end if

 end function StandardAtmos_kPa_2_km

 !=======================================================================

 !+
 !  Exner functions
 !
 !  Defined here as (p/P0)**KAPPA

 ! Where KAPPA = R/CP = 0.286
 ! P0 = 1.0e5 Pa

 ! CAREFUL:  The term Exner function  can also be used for CP * (p/P0)**KAPPA
 ! Hence notation Exner_nd  - non dimensional version
 !
 ! Tabulate  :
    ! defines the exner-function for every 1000 pa from zero to 1.3e+5 pa
    ! in a table for efficient interpolation (same procedure as used in
    ! the nwp-model, see mb1e.f)
 !
 ! Exner_nd returns the non-dimesnional excner function (p/p0)**R/CP
 !
 ! Added 7/4/2005, Dave, based upon tpi code from tiphys
 ! Test prog at end along with results.
 !----------------------------------------------------------------------------



  !-------------------------------------------------------------------
  subroutine Exner_tab()
  !
    real    :: p
    integer :: i

    do i = 1,131
      p = PBAS + i*PINC
      ! tpi(i) = CP*(p/1.0e+5)**KAPPA ! With CP!!!!
      tab_exf(i) = (p/1.0e+5)**KAPPA  ! Without CP
    enddo

  end subroutine Exner_tab
  !-------------------------------------------------------------------

  elemental function Exner_nd(p) result(exf)

     real, intent(in) :: p    ! Pressure, p
     real :: exf, x1
     integer :: ix1

        x1 = (p-PBAS)/PINC
        ix1 = floor(x1)
        exf =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))

  end function Exner_nd
  !-------------------------------------------------------------------

  elemental function Tpot_2_T(p) result(fTpot)
     ! Identical to Exner_nd
     ! Usage:   T = Tpot * Tpot_2_T(p)

     real, intent(in) :: p    ! Pressure, p
     real :: fTpot, x1
     integer :: ix1

        x1 = (p-PBAS)/PINC
        ix1 = x1
        fTpot =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))

  end function Tpot_2_T
  !-------------------------------------------------------------------
  elemental function T_2_Tpot(p) result(fT)
     ! Iinvese of Exner_nd
     ! Usage:   Tpot = T * T_2_Tpot(p)

     real, intent(in) :: p    ! Pressure, p
     real :: fT, exf, x1
     integer :: ix1

        x1 = (p-PBAS)/PINC
        ix1 = x1
        exf =  tab_exf(ix1) + (x1-ix1)*(tab_exf(ix1+1) - tab_exf(ix1))
        fT = 1.0/exf

  end function T_2_Tpot
  !-------------------------------------------------------------------


  function great_circle_distance(fi1,lambda1,fi2,lambda2) result(dist)

    !compute the great circle distance between to points given in 
    !spherical coordinates. Sphere has radius 1.
    real, intent(in) ::fi1,lambda1,fi2,lambda2 !NB: in DEGREES here
    real :: dist


    !ds sind not allowed in gfortran, so replaced. Also, 360 removed
    !dist=2*asin(sqrt(sind(0.5*(lambda1-lambda2+360.0))**2+&
    !     cosd(lambda1+360.0)*cosd(lambda2+360.0)*sind(0.5*(fi1-fi2+360.0))**2))

    dist=2*asin(sqrt(sin(DEG2RAD*0.5*(lambda1-lambda2))**2+&
         cos(DEG2RAD*lambda1)*cos(DEG2RAD*lambda2)*&
           sin(DEG2RAD*0.5*(fi1-fi2))**2))

  end function great_circle_distance


!program Test_exn
!  use Exner_ml
!  use PhysicalConstants_ml, only : KAPPA
!  implicit none
!
!  real :: p, exf1, exf2
!  integer :: i
!
!  call Exner_tab()
!
!   do i = 1, 20
!     p = 0.05 * i * 1.0e5
!     exf1 = Exner_nd(p)
!     exf2 = (p*1.0e-5)**KAPPA
!     print "(f8.3,4f12.5)", 1.0e-2*p, exf1, exf2, Tpot_2_T(p), T_2_Tpot(p)
!  end do
!end program Test_exn
    
! Results:
! p(mb)        exf1         exf2      Tpot_2_T    T_2_Tpot
!  50.000     0.42471     0.42471     0.42471     2.35455
! 100.000     0.51778     0.51778     0.51778     1.93133
! 150.000     0.58141     0.58141     0.58141     1.71997
! 200.000     0.63124     0.63124     0.63124     1.58418
! 250.000     0.67282     0.67282     0.67282     1.48629
! 300.000     0.70881     0.70881     0.70881     1.41081
! 350.000     0.74075     0.74075     0.74075     1.34999
! 400.000     0.76957     0.76957     0.76957     1.29943
! 450.000     0.79592     0.79592     0.79592     1.25641
! 500.000     0.82025     0.82025     0.82025     1.21913
! 550.000     0.84291     0.84291     0.84291     1.18637
! 600.000     0.86414     0.86414     0.86414     1.15722
! 650.000     0.88414     0.88414     0.88414     1.13105
! 700.000     0.90307     0.90307     0.90307     1.10734
! 750.000     0.92105     0.92105     0.92105     1.08571
! 800.000     0.93820     0.93820     0.93820     1.06587
! 850.000     0.95461     0.95461     0.95461     1.04755
! 900.000     0.97033     0.97033     0.97033     1.03058
! 950.000     0.98544     0.98544     0.98544     1.01477
!1000.000     1.00000     1.00000     1.00000     1.00000
end module Functions_ml
