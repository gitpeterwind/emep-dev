! <EmisDef_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                         module EmisDef_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!_____________________________________________________________________________
implicit none

  !! No. emission files to be read for non-speciated

   integer, public, parameter :: NEMIS_FILES =      7 

   !/** The names used below must have length of 6 characters and must
   !    belong to the full list given in Emissions_ml, as they will
   !    be used together with information in that module. An error
   !    message (fatal)  will be produced if this is not the case.
   !-----------------------------------------------------------------

    character(len=6), public, save, dimension(NEMIS_FILES) :: &
      EMIS_NAME  = &
      (/ "sox   ", "co    "   &   ! =non-split first
       , "nh3   ", "pm25  ", "pmco  "   &
       , "nox   ", "voc   "  /)                       ! =to be split

    !----------------- basic emissions file definitions --------------------!
    !  Here we define the parameters *not* likely to change often           !
    !  between different  model versions - e.g. the size and characteristics!
    !  of emission files and sector splits.                                 !
    !-----------------------------------------------------------------------!
    !  What is "Flat emissions":                                            !
    !  Most emission sources will have a seasonal, weekly, daily cycle. For !
    !  some sources there is no cycle, or the emission cycle is not known.  !
    !  For these source emissions will be constant (or flat) throughout the !
    !  year.                                                                !
    !-----------------------------------------------------------------------!


    ! Note on SNAP sectors:
    ! ----------------------
    ! SNAP1  =  public power stations,150m nat gas
    ! SNAP2  =  Comm./inst. combustion
    ! SNAP3  =  Industrial combustion !60m nat gas
    ! SNAP4  =  Production processes
    ! SNAP5  =  Extracton fossil fuels
    ! SNAP6  =  Solvents
    ! SNAP7  =  Road traffic
    ! SNAP8  =  Other mobile (trains+planes, ...)
    ! SNAP9  =  Waste! .. ds + some ground level
    ! SNAP10 =  Agriculture
    ! SNAP11 =  Nature



  !/.. First, define here the characteristics of the EMEP/CORINAIR
  !    SNAP sector emissions data-bases:



   integer, public, parameter :: NCMAX  =  11  ! Max. No. countries per grid
   integer, public, parameter :: FNCMAX =  20  ! Max. No. countries (with
                                               ! flat emissions) per grid

   !/.. Sector specific information

   integer, public, parameter :: &
          NSECTORS  = 11       ! Number of SNAP-sectors in emissions
!
! SCENARIO
   integer, public, parameter :: &
          ANTROP_SECTORS=10, &   ! Non-natural sectors
          ISNAP_NAT  = 11,&      ! SNAP index for volcanoe emissions
          ISNAP_SHIP = 8         ! SNAP index for flat emissions,e.g ship
                                 ! Note that flat emissions do NOT necessarily
                                 ! belong to the same SNAP sector

!  New vertical allocation from SNAP sectors.
!       - allocations are guesswork  - should be investigated

   integer, public, parameter :: NEMISLAYERS = 7
   real, public, parameter, &
     dimension(NEMISLAYERS,NSECTORS) ::    &

   VERTFAC =                    &  ! Vertical distribution of SNAP emissions
             reshape (          &  ! Vertical distribution of SNAP emissions
!       Ground , ...       High
    (/  0.0    , 0.00, 0.08, 0.46, 0.29, 0.17, 0.00,  & ! SNAP1
        0.5    , 0.50, 0.00, 0.00, 0.00, 0.00, 0.0,   & ! SNAP2
        0.0    , 0.04, 0.19, 0.41, 0.30, 0.06, 0.0,   & ! SNAP3
        0.9    , 0.10, 0.00, 0.00, 0.00, 0.00, 0.0,   & ! SNAP4
        0.9    , 0.10, 0.00, 0.00, 0.00, 0.00, 0.0,   & ! SNAP5
        1.0    , 0.00, 0.00, 0.00, 0.00, 0.00, 0.0,   & ! SNAP6
        1.0    , 0.00, 0.00, 0.00, 0.00, 0.00, 0.0,   & ! SNAP7
        1.0    , 0.00, 0.00, 0.00, 0.00, 0.00, 0.0,   & ! SNAP8
        0.1    , 0.15, 0.40, 0.35, 0.00, 0.00, 0.0,   & ! SNAP9
        1.0    , 0.00, 0.00, 0.00, 0.00, 0.00, 0.0,   & ! SNAP10
        1.0    , 0.00, 0.00, 0.00, 0.00, 0.00, 0.0    & ! SNAP11
        /), &
        (/NEMISLAYERS,NSECTORS /) )!hf stakheigth


  ! Biogenics
   integer, public, parameter ::   NBVOC = 2   
   character(len=8),public, save, dimension(NBVOC) :: &
                                   BVOC_USED = (/ "isoprene","terpene "/)   

   !SeaSalt
   integer, public, parameter ::  NSS   = 2 &   ! number of sea salt size modes
                                 ,QSSFI = 1 &   ! production of fine SS
                                 ,QSSCO = 2     ! production of coarse SS  


   !/** Lightning and aircraft NOx. 
    logical, public, parameter :: AIRNOX   = .true.   ! Gives NOx emission

   !/** Volcanos. 
    logical, public, parameter :: VOLCANOES  = .true.  ! Gives Volcanos

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
                     end module EmisDef_ml
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!_____________________________________________________________________________
