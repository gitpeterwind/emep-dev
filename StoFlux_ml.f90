! <StoFlux_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module StoFlux_ml
  use My_outputs_ml, only : 
  use DO3SE_ml, only : do3se
  use LandDefs_ml, only : LandType, STUBBLE,NLanduse_DEF, iLC_grass
  use LocalVariables_ml, only : L, Grid, Sub
  use MicroMet_ml, only : AerRes, Wind_at_h
  use ModelConstants_ml, only : NLANDUSEMAX, dt_advec, DEBUG_STOFLUX
  use Par_ml, only : MAXLIMAX, MAXLJMAX
  use PhysicalConstants_ml, only : AVOG, KARMAN
  use TimeDate_ml, only : current_date
  use Wesely_ml, only : WES_O3, Rb_cor
  implicit none
  private

!  public :: Init_StoFlux
  public :: Setup_StoFlux
  public :: Calc_StoFlux

  logical, public, parameter  :: STO_FLUXES = .true.
  logical, public,  dimension(NLANDUSEMAX), save :: luflux_wanted

  real, private :: c_hveg, u_hveg ! Values at canopy top, for fluxes
  real, public,  dimension(NLANDUSEMAX), save :: &
      lai_flux,      & ! Fluxes to total LAI
      unit_flux        ! Fluxes per m2 of leaf area (flag-leaf)

 real,   public,save,dimension(MAXLIMAX,MAXLJMAX) :: &
       SumVPD ,   &   ! For critical VPD calcs, reset each day  
       old_gsun       !

  real, private, save :: nmole_o3, ppb_o3    ! O3 in nmole/m3, ppb
  real, private, save :: gext_leaf = 1.0/2500.0
  real, private :: rc_leaf, rb_leaf, Fst

! Converts from mol/cm3 to nmole/m3
  real, private, parameter :: NMOLE_M3 = 1.0e6*1.0e9/AVOG  




contains
!  subroutine Init_StoFlux()
!    integer :: iL
!
!     luflux_wanted(:) = .false.
!     do iL = 1, NLanduse_DEF
!        if ( LandType(iL)%is_iam ) luflux_wanted(iL)  = .true.
!        if ( LandType(iL)%is_seminat ) luflux_wanted(iL)  = .true.
!     end do

!  end subroutine Init_StoFlux

  subroutine Setup_StoFlux(jd,xo3,xm)
 ! xo3 is in #/cm3. xo3/xm gives mixing ratio, and 1.0e9*xo3/xm = ppb
    real, intent(in) :: xo3   ! O3 mixing ratio
    real, intent(in) :: xm    ! Air mixing ratio
    integer, intent(in) :: jd ! daynumber
    integer, save :: old_daynumber

       ppb_o3   =  1.0e9* xo3/xm
       nmole_o3 =  xo3 * NMOLE_M3

       Sub(:)%FstO3 = 0.0


       ! resets whole grid on 1st day change
       if ( jd /= old_daynumber ) then
           SumVPD   = 0.0    ! For Critical VPD stuff, wheat
           old_gsun = 1.0e99 ! "     "
           old_daynumber = jd
       end if


  end subroutine Setup_StoFlux


  subroutine Calc_StoFlux(iL,Vg_ref,debug_flag)
    integer, intent(in) :: iL
    real, intent(in) :: Vg_ref
    logical, intent(in) :: debug_flag

    real :: loss   !SKIP ,sto_frac
    real :: Ra_diff,tmp_gsun
    integer :: i,j
    ! Evapotranspiration needs:
    real :: gv, gvcms  ! conductace for water vapour, mmole/ms/s and cm/s

    i = Grid%i
    j = Grid%j

! take care of  temperate crops, outside growing season
    if ( L%hveg < 1.1 * L%z0 ) then 

          Sub(iL)%FstO3 = 0.0
          Sub(iL)%cano3_ppb   = 0.0  !! Can't do better?
          Sub(iL)%EvapTransp = 0.0   ! No evapo-transpiration ?
          if ( DEBUG_STOFLUX .and. debug_flag ) &
            write(6,"(a15,f8.3,a8,f8.3)")  "FST - too low ", &
                L%hveg, " <  1.1*",  L%z0

    else 
   !=======================
   ! The fraction going to the stomata = g_sto/g_sur = g_sto * R_sur.
   ! Vg*nmole_o3 is the total deposition flux of ozone, but
   ! we calculate the actual flux later (once we know DepLoss(O3)).
   ! For now we just calculate the g_sto*R_sur bit:
   ! (Caution - g_sto is for O3 only)


          Ra_diff = AerRes(max( L%hveg-L%d, STUBBLE) , Grid%z_ref-L%d,&
                       L%ustar,L%invL,KARMAN)

          c_hveg         = nmole_o3 * ( 1.0 - Ra_diff * Vg_ref )

         ! Need to be careful with scope. L is within iL loop, whereas Sub
         ! will be kept throughout i,j calculations:
          Sub(iL)%cano3_ppb   = ppb_o3   * ( 1.0 - Ra_diff * Vg_ref )


         !Could be coded faster with Ra....
          u_hveg  = Wind_at_h( Grid%u_ref, Grid%z_ref, L%hveg, &
                      L%d, L%z0, L%invL )

          rc_leaf = 1.0/(L%g_sto+ gext_leaf)

          !McNaughton + van den Hurk:
 
          if ( do3se(iL)%Lw > 0 )  then
             rb_leaf = 1.3 * 150.0 * sqrt(do3se(iL)%Lw/u_hveg)
          else ! default (CAREFUL!)
             rb_leaf = 1.3 * 150.0 * sqrt(0.05/u_hveg)
          end if

          ! VPD limitation for wheat

         if ( LandType(iL)%is_iam .and. LandType(iL)%is_crop ) then
              if( L%g_sun > 0.0 ) SumVPD(i,j) = SumVPD(i,j) + L%vpd*dt_advec/3600.0
              tmp_gsun = L%g_sun
              if ( SumVPD(i,j) > 8.0 ) L%g_sun = min( L%g_sun, old_gsun(i,j) )
              !if( DEBUG_STOFLUX .and. debug_flag ) then
              !   write(6,"(a8,3i3,4f8.3,4es10.2)") "SUMVPD ", &
              !    current_date%month, current_date%day, current_date%hour, &
              !      L%rh, L%t2C,  L%vpd, SumVPD(i,j), &
              !       old_gsun(i,j), tmp_gsun, L%g_sun , L%g_sto
              !endif

              old_gsun(i,j) = L%g_sun
          end if

         ! Flux in nmole/m2/s:
          !ds leaf_flux(iL) = c_hveg * rc_leaf/(rb_leaf+rc_leaf) * L%g_sun 
          Sub(iL)%FstO3 = c_hveg * rc_leaf/(rb_leaf+rc_leaf) * L%g_sun 

! ======   CLOVER  ===========================================================
      ! For Clover we have a very special procedure, using O3 from grassland
      ! to scale the fluxes. As grassland is entered in Inputs.Landuse before
      ! clover we can assume that Fst and O3 for grassland are available and 
      ! correct
          if( LandType(iL)%is_clover) then

             Sub(iL)%FstO3 = Sub(iLC_grass)%cano3_ppb/Sub(iL)%cano3_ppb * Sub(iL)%FstO3

             if ( DEBUG_STOFLUX  .and. debug_flag ) then
                write(6,"(a,es12.3,2f12.3)") "CLOVER ",  &
                 Sub(iL)%FstO3, Sub(iLC_grass)%cano3_ppb, &
                       Sub(iLC_grass)%cano3_ppb/Sub(iL)%cano3_ppb
             end if
          end if ! clover
! ======   CLOVER  =========================================================

! ======   Evapo-transpiration =============================================
!
         ! E = g * (Cleaf - Cair) = g * D(kPa)/101(KPa), cf. Cambell+Norman,p82
         ! If E in mole/m2/s

         ! and use 1/1.51 instead of 1/1.6. Checl
         ! Lisa had R in s/mm, hence 1000.0
         !
         !gv = 1/(Rb+Rc) = 1(Rb+1/gsto)
         !  Rb(O3) = 2.0 * Rb_cor(WES_O3)/(KARMAN*ustar) 
         !
         ! Step 1: Get gv in m/s units:

            gv = 1.0/ (2.0 * Rb_cor(WES_O3) /(KARMAN*L%ustar) + 1.0/ ( L%g_sto * L%LAI ) )
            gvcms = gv

         ! Step 2: Convert to mole/m2/s and for H2O
         ! mmol2sm = 8.3144e-8*L%t2  for O3, plus factor 1.6 for H2O
         ! * 0.001 -> mole/m2/s
         ! ms2molm2s = 1.6*1.0e-3/(8.3144e-8*L%t2)
         
            gv = gv * 1.6e-3/(8.3144e-8*L%t2)
         
         ! Step 3: 
         ! Mass flux density is E x 0.018 kg/mole -> kg/m2/s
         ! 1 kg/m2 = 1mm
         ! ie ms2kgm2s  = 1.6*1.0e-3/(8.3144e-8*L%t2)
         ! ie ms2mm     = 1.6*1.0e-3/(8.3144e-8*L%t2)

            Sub(iL)%EvapTransp = 0.018 *  L%vpd/101.0  * gv

         !!L%g_sto * L%LAI  * 1.6e-3/(8.3144e-8*L%t2)   ! Evapo-transpiration
         !                 (rb_leaf/1.6 + 0.0224*(L%t2/273.0) 


          if ( DEBUG_STOFLUX .and. debug_flag ) then 
            write(6,*)  "STOFLUX? ", iL, LandType(iL)%flux_wanted
            write(6,"(a8,3i3,i4,f6.2,2f8.1,2es10.2,f6.2,3es12.3)")  "STOFLUX ",&
            iL, current_date%month, current_date%day, current_date%hour, &
            L%LAI, nmole_o3, c_hveg, L%g_sto, L%g_sun, u_hveg,&
            Sub(iL)%cano3_ppb, Sub(iL)%FstO3, Sub(iL)%EvapTransp
            write(6,"(a8,2es12.3)")  "STOFLUXG ", gv, gvcms
          end if

    end if

  end subroutine Calc_StoFlux
!..............................................................................
end module StoFlux_ml
