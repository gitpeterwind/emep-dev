! <SoilWater_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module SoilWater_ml
 use GridValues_ml,     only : debug_proc, debug_li, debug_lj
 use MetFields_ml,      only : SoilWater, SoilWater_deep
 use ModelConstants_ml, only : DEBUG_SOILWATER
    implicit none

    public :: Set_SoilWater
    real, dimension(366), public, save :: SWP = 0.0  ! daily soil water potential 
                                              ! in  MPa

contains
   subroutine Set_SoilWater()
      integer :: i, j
      if ( DEBUG_SOILWATER .and. debug_proc ) then
         i = debug_li
         j = debug_lj
         write(*,*) "DEBUG_SW: ", SoilWater(i,j,1), SoilWater_deep(i,j,1)
      end if
   end subroutine Set_SoilWater
end module SoilWater_ml
