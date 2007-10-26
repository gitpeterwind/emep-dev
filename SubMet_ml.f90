! <SubMet_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module SubMet_ml
!=============================================================================
!+
! Description
!  Module for setting some local grid-cell data (mainly from NWP)
!  and for calculating sub-grid  meteorology for each land-use. 
!  The sub-grid part of this module is also undergoing constant change!!
!=============================================================================


use CheckStop_ml, only : CheckStop
use LandDefs_ml,   only: LandType
use Landuse_ml,    only: LandCover
use LocalVariables_ml, only: Grid, Sub
use MicroMet_ml, only :  PsiM, AerRes    !functions
use PhysicalConstants_ml, only : PI, RGAS_KG, CP, GRAV, KARMAN, CHARNOCK, T0

implicit none
private


  public :: Get_Submet    ! calculates met. data for sub-grid areas


  logical, private, parameter ::  DEBUG_SUB = .false.  ! for extra tests/printouts

contains
!=======================================================================

  subroutine Get_Submet(iL, debug_flag )

!---------------------------------------------------------------
!  Sub-grid calculations of  stability, Ra and ustar for this landuse
!---------------------------------------------------------------
!
!..The profile manipulation is introduced to calculate different Ra and Rb
!..values for different landuse types (e.g. from SEI). Thus, different z0,
!..are assumed within each EMEP square, from which one averaged value is 
!..available as the basic input from the NWP-model. Therefore, we first have a
!..grid square average for u*, T*, L and z0 (i.e. from the NWP-model).  Then,
!..we calculate the wind speed u from the wind profile (according to the
!..M-O similarity) and the NWP-model data, and from this u we calculate
!..new u* values for each z0, i.e. landuse, within each grid square. Here,
!..we have to use the old L, even although L=f(u*), since the new u* is
!..unknown before it is calculated. (An iteration process is possible.) 
!..Anyway, by using the new u*'s, we obtain
!..new T*'s, L's for each landuse type and finally Ra's and ustars for each
!..subgrid area.
!
!.. For further details, see EMEP report 1/2003 (and before that 3/95 ...)
!
! The subroutine Get_Submet also generates output for 
! two terms utilized in Rsurface, namely rh (the relative humidity
! term required for the evaluation of the stomatal compensation point 
! point for NH_3) and vpd (the vapour pressure deficit term required
! for the evaluation of stomatal conductance).

!----------------------------------------------------------------- 
!..In

   integer, intent(in) :: iL      ! lu index
   logical, intent(in) :: debug_flag   ! set true for wanted grid square

   ! IMPORTANT - ASSUMES INITIAL VALUES SET FOR USTAR, INVL, ....

!.. Local
    real :: rho_surf               ! Density at surface (2 m), kg/m3
    real :: z_1m                   !  1m above vegetation
    real :: z_3m                   !  3m above ground, or top of trees
    real :: z_3md                  !  minus displacemt ht.
    
    
    logical, save ::  my_first_call = .true.
    integer, parameter ::  NITER = 1           ! no. iterations to be performed

    integer :: iter                ! iteration variable

   ! For vapour pressure calculations

    real, parameter :: ESAT0=611.0   ! saturation vapour pressure at 
                                     ! T=0 deg. C (Pa)
      
    real :: qw                       ! specific humidity (kg/kg) corrected  
                                     ! down to z_0+d metres above the ground 
    real :: esat    ! saturation vapour pressure  (Pa)
    real :: e       ! vapour pressure at surface
    real :: Ra_2m   ! to get 2m qw


    ! initial guesses for u*, t*, 1/L
        Sub(iL)%ustar  = Grid%ustar      ! First guess = NWP value
        Sub(iL)%invL   = Grid%invL       ! First guess = NWP value
        Sub(iL)%Hd     = Grid%Hd         ! First guess = NWP value
        Sub(iL)%LE     = Grid%LE         ! First guess = NWP value
        Sub(iL)%t2     = Grid%t2         ! First guess = NWP value
        Sub(iL)%t2C    = Grid%t2C        ! First guess = NWP value

        Sub(iL)%is_water  = LandType(iL)%is_water
        Sub(iL)%is_forest = LandType(iL)%is_forest

     ! If NWP thinks this is a sea-square, but we anyway have land,
     ! the surface temps will be wrong and so will stability gradients.
     ! As a simple substitute, we assume neutral conditions for these
     ! situations.
        if( Grid%is_NWPsea  .and. (.not. Sub(iL)%is_water) ) then
             Sub(iL)%invL = 0.0
             Sub(iL)%Hd   = 0.0
        end if


!    The zero-plane displacement (d) is the height that
!     needs to be added in order to make the profile theories
!     work over a tall vegetation (see e.g. Stull (1988), p.381). 
!     This corresponds to moving our co-ordinate system upwards 
!     by a distance d.
! 
!     by definition: u(z0+d) = 0
!
!     it has been observed that d is approximately 0.7 times
!     the mean height of the vegetation (h) and z0=h/10
!     (see e.g. Stull, 1998, Garratt, 1992.).
!     The reference height for u* transformation is then
!     taken arbitrarily at about 45m, the height of the
!     centre of the EMEP grid cell.
          

!.. For water, we introduce a new zero plane 
!   displacement and use the Charnock relation to calculate the z0-values
! nb - addded max 1cm limit for z0 over sea, because of problems
! caused by z0>1m. Garratt (section 4.1, Fig 4.2) suggested that Charnock's
! relation is only valid for  u* < 1 m/s, which gives z0 < 1 cm/s.


        if ( Sub(iL)%is_water ) then ! water
             Sub(iL)%d  = 0.0
             Sub(iL)%z0 = CHARNOCK * Sub(iL)%ustar * Sub(iL)%ustar/GRAV
             Sub(iL)%z0 = max( Sub(iL)%z0 ,0.01)
             z_1m   = 1.0       ! 1m above sea surface
             z_3m   = 3.0       ! 3m above sea surface
        else
             Sub(iL)%d  =  0.7 * Sub(iL)%hveg
             Sub(iL)%z0 = max( 0.1 * Sub(iL)%hveg, 0.001) !  Fix for deserts, 
               ! ice, snow (where, for bare ground, h=0 and hence z0=0)

           !Heights relative to displacement height, d:

             z_1m   = (Sub(iL)%hveg + 1.0) - Sub(iL)%d
             z_3m   = max(3.0,Sub(iL)%hveg)
!CHECK!!!! z_3m z_3md....

        end if
          
        Sub(iL)%z_refd = Grid%z_ref - Sub(iL)%d  !  minus displacement height
        z_3md  = z_3m  - Sub(iL)%d               !  minus displacement height


    do iter = 1, NITER 

        ! ****
        !   PsiM calculates the stability functions for momentum
        !   at heights z_ref (about 45m) & z0
        ! ****               
        !..calculate friction velocity based first on NWP-model PsiM-values 
        !..and u_ref. The NWP-model PsiM-values are used despite the fact that
        !..L=F(u*), since we do not know the EMEP subgrid averaged 
        !..z0-values ...

        if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
            write(6,"(a12,i2,5f8.3,f12.3)") "UKDEP SUBI", iter, &
                       Sub(iL)%hveg, Sub(iL)%z0, Sub(iL)%d, &
                         Sub(iL)%z_refd, z_3md, Sub(iL)%invL
        end if

       Sub(iL)%ustar = Grid%u_ref * KARMAN/ &
        (log( Sub(iL)%z_refd/Sub(iL)%z0 ) - PsiM( Sub(iL)%z_refd*Sub(iL)%invL)&
          + PsiM( Sub(iL)%z0*Sub(iL)%invL ) )

       Sub(iL)%ustar = max( Sub(iL)%ustar, 1.0e-2)

    !  We must use L (the Monin-Obukhov length) to calculate deposition,
    ! Thus, we calculate T* and then L, based on sub-grid data. 


        rho_surf = Grid%psurf/(RGAS_KG * Sub(iL)%t2 )


    ! New 1/L value ....


        Sub(iL)%invL =  -KARMAN * GRAV * Sub(iL)%Hd / &
      ( CP*rho_surf*Sub(iL)%ustar*Sub(iL)%ustar*Sub(iL)%ustar * Sub(iL)%t2)

      !.. we limit the range of 1/L to prevent numerical and printout problems
      !   This range is very wide anyway.

        Sub(iL)%invL  = max( -1.0, Sub(iL)%invL ) !! limit very unstable
        Sub(iL)%invL  = min(  1.0, Sub(iL)%invL ) !! limit very stable

    end do ! iter


    if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
        write(6,"(a12,10f9.3)") "UKDEP SUBL", Sub(iL)%z0, Sub(iL)%d, &
          Sub(iL)%z_refd, 0.001*Grid%psurf, Sub(iL)%t2, rho_surf, Sub(iL)%Hd,&
              Sub(iL)%ustar, Sub(iL)%t2, Sub(iL)%invL
    end if

    if ( DEBUG_SUB .and. debug_flag ) then 
         if ( my_first_call ) then ! title line

                write(unit=*, fmt="(a6,3a3, a6, 3a8,2a7, 2a6)") &
                 "STAB ", "mm", "dd", "hh", "t2_C", "Hd", &
                 "L_nwp", "L  ", "z/L_nwp", "z/L ", "u*_nwp", "u*"
                my_first_call = .false.
         end if

            write(unit=*, &
              fmt="(a6,4i3, f6.1, 3f8.2, 2f7.2, 2f6.2)") "SUBB", iL, &
              999, & !SUBcurrent_date%month, &
              999, & !SUBcurrent_date%day, &
              999, & !SUBcurrent_date%hour, &
              Sub(iL)%t2C, Sub(iL)%Hd, Grid%invL, Sub(iL)%invL, &
              Sub(iL)%z_refd*Grid%invL, Sub(iL)%z_refd*Sub(iL)%invL, &
                 Grid%ustar, Sub(iL)%ustar
    end if



!     *** Aerodynamic resistances for each landuse type k ***
!      Ra_ref is used to estimate the aerodynamic resistance to latent
!      heat transfer from height z_ref to z0+d and from height  h+3 to 
!      z0+d, respectively.
!      Only Ra_ref and Ra_3m are used in main code.
      
        Sub(iL)%Ra_ref = AerRes(Sub(iL)%z0,Sub(iL)%z_refd,Sub(iL)%ustar,&
            Sub(iL)%invL,KARMAN)
        Sub(iL)%Ra_3m  = AerRes(Sub(iL)%z0,z_3md,Sub(iL)%ustar,Sub(iL)%invL,KARMAN)
        Ra_2m  = AerRes(Sub(iL)%z0,1.0+z_1m,Sub(iL)%ustar,Sub(iL)%invL,KARMAN)

    if ( DEBUG_SUB ) then
       if ( Sub(iL)%Ra_ref < 0 .or. Sub(iL)%Ra_3m < 0 &
           .or. Ra_2m < 0  ) call CheckStop("RAREF NEG ")
      if ( Sub(iL)%Ra_3m > Sub(iL)%Ra_ref ) &
           call CheckStop("ERROR!!! Ra_ref<Ra_3")
    end if


!  *****  Calculate rh and vpd  *********

!  
!....The model has the specific humidity qw_ref for the lowest model layer as 
!    an input obtained from the NWP model.  We now correct this down to 
!    z_0+d metres above the ground using q(z_0+d) = q(z2) + Q*Ra_ref,
!    where Ra_ref is the same as calculated above. Q is latent heat flux 
!    (W/m^2) divided by latent heat of vap. (2.5e6  J/kg). From input data, 
!    Hl is directed downwards (!), so  we used Q = -Hl/2.5e6
!    but above ds reversed it, so we go back to the normal micromet. 
!    formula Q = Hl/2.5e6

                       
! 2m qw:
        qw = Grid%qw_ref  + Sub(iL)%LE/2.5e6 * ( Sub(iL)%Ra_ref - Ra_2m)  

      !..   qw is in kg/kg  so  e = qw*psurf/epsilon
      !..   to get e in Pascal.

       e = qw * Grid%psurf/0.622

!The equation below relies on the following assumptions: epsilon=0.622
!and L=2.5x10^6 Joules/Kg, where "epsilon" and "L" denote the ratio of 
!the gas constants for dry air against water vapour and the latent
!heat of vaporization, respectively.


     esat = ESAT0 * exp(0.622*2.5e6*((1.0/T0) - (1.0/Sub(iL)%t2))/RGAS_KG )

    if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
        print "(a15,2f12.6,2f12.3)", "UKDEP SUB water", Grid%qw_ref, qw, &
              Sub(iL)%LE, 100.0*e/esat
    end if

   ! Straighforward calculation sometimes gives rh<0 or rh>1.0 -
   ! probably due to mismatches between the assumptions used for the stability
   ! profile here and in HIRLAM. Here we set crude limits on e to prevent
   ! impossible rh values at least:


     e = max(0.001*esat,e)    ! keeps rh >= 0.1%
     e = min(esat,e)          ! keeps rh <= 1
     Sub(iL)%rh = e/esat

!  ****  leaf sat. vapour pressure

      Sub(iL)%vpd    =  0.001*(esat-e)     ! gives vpd in kPa !
      Sub(iL)%vpd    =  max(Sub(iL)%vpd,0.0) 


    if ( DEBUG_SUB .and. debug_flag ) then !!  .and. &
        write(6,"(a22,2f12.4)") "UKDEP SUB7 e/esat, rh", e/esat, Sub(iL)%rh
    end if

  end subroutine Get_Submet
! =====================================================================

end module SubMet_ml
