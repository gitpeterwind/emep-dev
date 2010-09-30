! <My_DryDep_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module My_DryDep_ml    ! DryDep_ml
!+
! Module to define the dry deposition components and rates. We
! define the min (vd_min) and max dep. velocity (Vg) first and then derive the
! daytime addition (vd_day).
!/**************************************************************************
!  Specifies which of the possible species (from Wesely's list)
!  are required in the current air pollution model   
!/**************************************************************************

 use CheckStop_ml,        only : CheckStop, StopAll
 use ChemSpecs_tot_ml,    only : O3
 use ChemSpecs_adv_ml     ! Needed for all IXADV_ 
 use OwnDataTypes_ml,     only : print_Deriv_type, depmap
 use Wesely_ml
 implicit none
 private

 ! WE NEED A FLUX_CDDEP, FLUX_ADV FOR OZONE;
 ! (set to one for non-ozone models)

  logical, public, parameter :: STO_FLUXES = .true.
  integer, public, parameter :: FLUX_CDDEP  = CDDEP_O3
  integer, public, parameter :: FLUX_ADV   = IXADV_O3
  integer, public, parameter :: FLUX_TOT   = O3


  logical, public, parameter :: COMPENSATION_PT = .false. 

!SKIP   real, public, save, dimension(NSPEC_ADV) :: DepLoss   ! Amount lost

  ! Maps from adv index to one of calc indices
  !SKI integer, public, save, dimension(NSPEC_ADV) :: DepAdv2Calc

   logical, private, save :: first_call = .true.

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   ! .... Define the mapping between the advected species and
   !      the specied for which the calculation needs to be done.
   !  We also define the number of species which will be deposited in
   ! total, NDRYDEP_ADV. This number should be >= NDRYDEP_GASES
   ! The actual species used and their relation to the CDDEP_ indices
   ! above will be defined in Init_DepMap

       include 'CM_DryDep.inc'

   !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


end module My_DryDep_ml

