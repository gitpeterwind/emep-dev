! <Pollen_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2013 met.no
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
module Pollen_const_ml
!-----------------------------------------------------------------------!
! Birch pollen emission calculation based on
! M. Sofiev et al. 2006, doi:10.1007/s00484-006-0027-x
!
! Pollen emission based upon meteorology paparameters, and heatsum.
! Pollen particles are assumed of 22 um diameter and 800 kg/m3 density. 
!-----------------------------------------------------------------------!
  use PhysicalConstants_ml, only: PI
  use ModelConstants_ml,    only: USE_POLLEN,DEBUG=>DEBUG_POLLEN
  use CheckStop_ml,         only: CheckStop
  use ChemChemicals_ml,     only: species_adv
  implicit none
  public

  real, parameter  :: &
    T_cutoff = 273.2+3.5, & ! Cut-off temperature [K]
    dH       = 50*24*3600,& ! Flowering period in degree seconds
    PREC_MIN = 0.0,       & ! Min cut-off precipitation [mm/h]
    PREC_MAX = 0.5,       & ! Max cut-off precipitation [mm/h]
    N_TOT    = 3.7e8,     & ! Available pollen [grains/m2] (Could be a fild)
    RH_LOW   = 0.50,      & ! Min cut-off relative humidity [ ]
    RH_HIGH  = 0.80,      & ! Max cut-off relative humidity [ ]
    PROB_IN  = 0.2,       & ! Probability for flowering to start
    PROB_OUT = 0.2,       & ! Probability for flowering to end 
                            ! (could be assumed to be larger than PROB_IN)
    D_POLL   = 22.0,      & ! Pollen grain diameter [um]
    POLL_DENS= 800e3        ! Pollen density [g/m3]

  real, parameter  :: &
    grain_wt = POLL_DENS*PI*(D_POLL*1e-6)**3/6.0, &! 1 grain weight [g]
    ug2grains= 1e-6/grain_wt

contains
subroutine pollen_check()
  logical,save :: fist_call=.true.
  if(.not.fist_call)return
  fist_call=.false.
  call CheckStop(USE_POLLEN.and.all(species_adv(:)%name/="POLLEN_B"),&
    "USE_POLLEN on model compiled without pollen")
  call CheckStop(DEBUG.and..not.USE_POLLEN,&
    "DEBUG_POLLEN on run without USE_POLLEN")
endsubroutine pollen_check
endmodule Pollen_const_ml
