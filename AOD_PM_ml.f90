! <AOD_PM_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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


 module AOD_PM_ml

  !-----------------------------------------------------------------------! 
  ! Calculates Aerosol Optical Depth (AOD) for 0.5 um radiation based on 
  ! Mie calculated average scattering efficiency and particle size 
  ! distribution. Complex refractive indices are from GADS
  !-----------------------------------------------------------------------!
 use Chemfields_ml,        only : AOD
 use ChemChemicals_ml,     only : species  
 use ChemGroups_ml,        only : AOD_GROUP
 use ChemSpecs_tot_ml
 use GridValues_ml,        only : i_fdom, j_fdom
 use MetFields_ml,         only : z_bnd
 use ModelConstants_ml,    only : KMAX_MID, KMAX_BND, KCHEMTOP,   &
                                  MasterProc, NPROC, DEBUG_i, DEBUG_j
 use Par_ml,               only : MAXLIMAX,MAXLJMAX   ! => x, y dimensions
 use PhysicalConstants_ml, only : AVOG
 use Setup_1dfields_ml,    only : xn_2d

  implicit none
  private
  !-----------------------------------------------------------------------!
 !// subroutines
  public ::   AOD_calc 

  real, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: kext 
!.. Specific extinction cross-sections are in My_Chem_ml
!  real, parameter, dimension(7) ::    &
!.. Specific extinction cross-sections [m2/g]
!        EXTSP  = (/8.5,  8.5,  8.5,  5.7,  9.0,  1.0, 3.0 /),  &
       !test EXTSP  = (/8.5,  8.5,  8.5,  5.7,  9.0,  0.9, 0.4 /),  &

!.. Specific extinction cross-sections for dry particles
!        EXTdry = (/5.0,  5.0,  5.0,  5.0,  9.0,  0.9, 0.2 /),  &
!.. Growth factor
!        GF     = (/1.7,  1.7,  1.7,  1.3,  1.0,  1.0, 2.0 /)

  contains

! <---------------------------------------------------------->

    subroutine AOD_calc (i,j,debug)

 !------------------------------------------------
 ! Calculates AOD
 !-------------------------------------------------

 implicit none


   integer, intent(in) :: i,j    ! coordinates of column
   logical, intent(in) :: debug

   integer :: k, n, itot
   real    :: kext_help
   real, parameter ::  lambda = 0.55e-6

!-----------------------------------------------------------------
!   AOD_GROUP = (/ SO4, NO3_F, NH4_F, EC_F_NEW, EC_F_AGE, POC_F, &
!       EXTC  = (/ 8.5, 8.5,   8.5,   7.5,      11.0,     5.7,   &
!                  SEASALT_F, SEASALT_C, DUST_NAT_F, DUST_NAT_C /)
!                  3.0,        0.4       1.0,         0.3,      /)
!__________________________________________________________________


  AOD(i,j)     = 0.0
  kext(i,j,:)  = 0.0

  do k =  KCHEMTOP, KMAX_MID   !_______________ vertical layer loop

!   kext(i,j,k)  = 0.0
 
!.. ===========================================================================
!..  Extinction coefficients: Kext [1/m] = SpecExtCross [m2/g] * mass [g/m3]
!..                           summed up for all components  
!..    xn_2d(ispec,k)*1.e15 *species(ispec)%molwt/AVOG   [molec/m3] -> [ng/m3]
!..                                                   [ng/m3 ] * 1e-9 -> [g/m3]
!..=>  xn_2d(ispec,k) * species(ispec)%molwt * 1.e6 / AVOG  [g/m3]
!.. ===========================================================================

    do n = 1, size(AOD_GROUP)
      itot = AOD_GROUP(n)

      kext(i,j,k) = kext(i,j,k) +   &
                    xn_2d(itot,k) * species(itot)%molwt * species(itot)%ExtC        
    enddo

     kext(i,j,k) = kext(i,j,k) * 1.0e6 / AVOG 

     if(debug .and. (k == 18 .or. k == KMAX_MID) )  &
            write(6,'(a17,i4,es15.3)') '> Ext. coeff', k, kext(i,j,k)

!.. Aerosol extinction optical depth : integral over all vertical layers
!.. [1/m} * [m]

      AOD(i,j) = AOD(i,j) + kext(i,j,k) * (z_bnd(i,j,k)-z_bnd(i,j,k+1))

!      if(debug .and. (k == 18 .or. k == KMAX_MID) )  & 
!      write(6,'(a25,i4,2es15.4,2f8.1)') '>> Kext AOD for layer', k,  &
!                kext(i,j,k), AOD(i,j), z_bnd(i,j,k), z_bnd(i,j,k+1)

  enddo                        !_______________ vertical layer loop

!  if(debug )  write(6,'(a30,2i4,es15.3)') '>>>  AOD  <<<',   &
!              i_fdom(i), j_fdom(j), AOD(i,j)


   end subroutine AOD_calc

 end module AOD_PM_ml

