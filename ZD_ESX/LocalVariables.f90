!> MODULE  <LocalVariables.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2007-2011 met.no
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

module LocalVariables
! -----------------------------------------------------------------------
! Near-surface meteorology and other variables for local area,
! e.g. for a measurement site or for a specific landuse within a grid square
! -----------------------------------------------------------------------

 !use ModelConstants, only: dp !! NLANDUSEMAX
!ESX use Wesely_ml,         only: NDRYDEP_CALC

implicit none
private

real,    private, parameter ::  NOT_SET = -999.     ! Fake value to check for
integer, private, parameter :: INOT_SET = -999      ! variables being set

integer, public, save :: iL = INOT_SET              ! Landuse index

! -----------------------------------------------------------------------
! 1) Grid data which should be okay for local
! -----------------------------------------------------------------------
!ESX SKIPPED  GridDat

! -----------------------------------------------------------------------
! 2) Near-surface & Sub-grid Veg/landcover and meteo data
! -----------------------------------------------------------------------
type, public :: LocDat !> Location-specific Data, was SubDat
  integer ::              &
   iL  = INOT_SET         & ! Landcover index
  ,SGS = INOT_SET         & ! Start, growing seasons (day num)
  ,EGS = INOT_SET           ! End, growing seasons (day num)
  logical :: &
    is_forest, is_water , is_veg, is_ice, is_crop

  real :: & 
     latitude  = NOT_SET  &  ! deg N   !ESX-added
    ,longitude = NOT_SET  &  ! deg E   !ESX-added
    ,t2C       = NOT_SET  & ! Surface (2m) temperature in degrees C
    ,t2        = NOT_SET  & ! Surface (2m) temperature in degrees K
    ,rh        = NOT_SET  & ! Relative humidity, fraction (0-1)
    ,rho_s     = NOT_SET  & ! Air density (kg/m3) at surface, here 2m
    ,psurf     = NOT_SET  & ! ESX added. Surface pressure (Pa)
    ,vpd       = NOT_SET  & ! Vapour pressure deficit  (kPa) ! CHECK UNITS
    ,EvapTransp= NOT_SET  & ! Evapotranspiration             ! CHECK UNITS
    ,precip    = NOT_SET  & ! Precipitation (mm/hour) 
    ,SWP       = NOT_SET  & ! SWP  ! CHECK UNITS
    ,fSW       = NOT_SET  & ! function for fSWP or fSMD or...
    ,ustar     = NOT_SET  & ! friction velocity, m/s
    ,wstar     = NOT_SET  & ! convective velocity scale, m/s
    ,invL      = NOT_SET  & ! 1/L, where L is Obukhiov length (1/m)
    ,Hd        = NOT_SET  & !  Sensible Heat flux, *away from* surface
    ,LE        = NOT_SET  & !  Latent Heat flux, *away from* surface
    ,hSL       = NOT_SET  & ! Height of surface layer (m), ESX added
    ,Hmix      = NOT_SET  & ! Height of surface layer (m), ESX added
    ,Ra_ref    = NOT_SET  &
    ,Ra_2m     = NOT_SET  &
    ,Ra_3m     = NOT_SET  &
    ,RgsO      = NOT_SET  & ! ground-surface resistances - set in DO3SE
    ,RgsS      = NOT_SET  & ! ground-surface resistances - set in DO3SE
    ,coverage  = NOT_SET  & ! Area covered (fraction)
    ,LAI       = NOT_SET  & ! Leaf area index (m2/m2)
    ,SAI       = NOT_SET  & ! Surface area index (m2/m2)
    ,hveg      = NOT_SET  & ! Height of veg.      (m)
    ,d         = NOT_SET  & ! displacement height (m)
    ,z_refd    = NOT_SET  & ! z_ref - d (m)
    ,z0        = NOT_SET  & ! roughness length    (m)
! Canopy-Associated Radiation
    ,PARsun    = NOT_SET  & ! photosynthetic active radn. for sun-leaves
    ,PARshade  = NOT_SET  & !  " " for shade leaves
    ,LAIsunfrac= NOT_SET  & ! fraction of LAI in sun
    ,alpha      = NOT_SET  & ! ESX added, light response
    ,gMax      = NOT_SET  & ! ESX added,  max gleaf, m/s
! outputs from Rsurface will include:
    ,g_sto     = NOT_SET  & ! stomatal conductance (m/s)
    ,g_sun     = NOT_SET  & ! g_sto for sunlit upper-canopy (flag) leaves
    ,f_sun     = NOT_SET  & ! f_env for SPOD?
    ,f_shade   = NOT_SET  & ! f_env for SPOD?
    ,f_env     = NOT_SET  & ! f_env for SPOD?
    ,f_vpd     = NOT_SET  & ! f_env for SPOD?
    ,f_light     = NOT_SET  & ! f_env for SPOD?
    ,f_temp      = NOT_SET  & ! f_env for SPOD?
    ,f_phen      = NOT_SET  & ! f_env for SPOD?
    ,f_min       = NOT_SET  & ! f_env for SPOD?
! and enable concentrations at canopy height:
    ,cano3_ppb  = 0.0     & ! Use 0.0 to make d_2d behave better
    ,cano3_nmole= 0.0     & ! Use 0.0 to make d_2d behave better
    ,FstO3      = 0.0       ! leaf O3 flux, nmole/m2/s
  !ESX real, dimension(NDRYDEP_CALC) :: & ! for species subject to dry depostion
!ESX8  real, allocatable,  dimension(:) :: & ! for species subject to dry depostion
!ESX8     Vg_ref   &  ! Grid average of Vg at ref ht. (effective Vg for cell)
!ESX8    ,Vg_3m    &  ! and at 3m
!ESX8    ,Gsur,Gns
end type LocDat

!ESX type(SubDat), public, dimension(0:NLANDUSEMAX), save :: Sub
type(LocDat), public, allocatable, dimension(:), save :: Sub ! NOTE- was ZERO:..
type(LocDat), public, save :: L         ! For just one land-class
type(LocDat), public, save :: ResetSub  ! Keeps NOT_SET values
end module LocalVariables
