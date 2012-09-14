! <StoFlux_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module StoFlux_ml
  use CheckStop_ml
  use DO3SE_ml, only : do3se, nSumVPD, SumVPD_LC
  use Io_Progs_ml, only : current_date, datewrite
  use LandDefs_ml, only : LandType, STUBBLE,NLanduse_DEF, iLC_grass
  use LocalVariables_ml, only : L, Grid, Sub
  use MicroMet_ml, only : AerRes, Wind_at_h
  use ModelConstants_ml, only : NLANDUSEMAX, dt_advec, DEBUG_STOFLUX
  use Par_ml, only : MAXLIMAX, MAXLJMAX
  use PhysicalConstants_ml, only : AVOG, KARMAN
  use SmallUtils_ml, only : find_index
  use Wesely_ml, only : WES_O3, Rb_cor
  implicit none
  private

  public :: Setup_StoFlux
  public :: Calc_StoFlux

  logical, public,  dimension(NLANDUSEMAX), save :: luflux_wanted

  real, private :: u_hveg ! Values at canopy top, for fluxes
  real, public,  dimension(NLANDUSEMAX), save :: &
      lai_flux,      & ! Fluxes to total LAI
      unit_flux        ! Fluxes per m2 of leaf area (flag-leaf)

 real,   public,save, allocatable,dimension(:,:,:) :: &
       SumVPD ,   &   ! For critical VPD calcs, reset each day  
       old_gsun       !
  integer, private, save, dimension(NLANDUSEMAX) :: mapSumVPD

  real, private, save :: gext_leaf = 1.0/2500.0
  real, private :: rc_leaf, rb_leaf, Fst



contains

  subroutine Setup_StoFlux(jd )

    integer, intent(in) :: jd ! daynumber
    integer, save :: old_daynumber
    logical, save :: my_first_call = .true.
    integer ::  istat, iL

     if ( my_first_call ) then
       allocate(SumVPD(MAXLIMAX,MAXLJMAX,nSumVPD),stat=istat)
       allocate(old_gsun(MAXLIMAX,MAXLJMAX,nSumVPD),stat=istat)
       do iL = 1, NLANDUSEMAX
         if ( do3se(iL)%VPDcrit > 0.0  ) then
           mapSumVPD(iL) = find_index( iL, SumVPD_LC )
         end if
       end do
     end if


       Sub(:)%FstO3 = 0.0


       ! resets whole grid on 1st day change
       if ( jd /= old_daynumber ) then
           SumVPD   = 0.0    ! For Critical VPD stuff, wheat
           old_gsun = 1.0e99 ! "     "
           old_daynumber = jd
       end if


  end subroutine Setup_StoFlux


  subroutine Calc_StoFlux(nLC,iL_used,debug_flag)
    integer, intent(in) :: nLC
    integer, dimension(nLC), intent(in) :: iL_used
    logical, intent(in) :: debug_flag

    real :: tmp_gsun
    integer :: i,j, iiL, iL, ivpd
    ! Evapotranspiration needs:
    real :: gv, gvcms  ! conductace for water vapour, mmole/ms/s and cm/s

    i = Grid%i
    j = Grid%j

    LC_LOOP: do iiL = 1, nLC
        iL = iL_used(iiL) 
        L = Sub(iL)

       ! take care of  temperate crops, outside growing season
        if ( L%hveg < 1.1 * L%z0 ) then 

          Sub(iL)%FstO3 = 0.0
          Sub(iL)%cano3_ppb   = 0.0  !! Can't do better?
          Sub(iL)%EvapTransp = 0.0   ! No evapo-transpiration ?
          if ( DEBUG_STOFLUX .and. debug_flag ) &
            call datewrite("FST - hveg < z0 ", iL, (/ L%hveg, L%z0 /) )

        else 
       !=======================
       ! The fraction going to the stomata = g_sto/g_sur = g_sto * R_sur.
       ! Vg*nmole_o3 is the total deposition flux of ozone, but
       ! we calculate the actual flux later (once we know DepLoss(O3)).
       ! For now we just calculate the g_sto*R_sur bit:
       ! (Caution - g_sto is for O3 only)


         !Could be coded faster with Ra....

          u_hveg = Wind_at_h( Grid%u_ref, Grid%z_ref, L%hveg, &
                      L%d, L%z0, L%invL )

          rc_leaf = 1.0/(L%g_sto+ gext_leaf)

          !McNaughton + van den Hurk:
 
          if ( do3se(iL)%Lw > 0 )  then
             rb_leaf = 1.3 * 150.0 * sqrt(do3se(iL)%Lw/u_hveg)
          else ! default (CAREFUL!)
             rb_leaf = 1.3 * 150.0 * sqrt(0.05/u_hveg)
          end if

          ! VPD limitation for wheat

         if ( do3se(iL)%VPDcrit > 0.0  ) then
              ivpd = mapSumVPD(iL)
              if( L%g_sun > 0.0 ) SumVPD(i,j,ivpd) = SumVPD(i,j,ivpd) + L%vpd*dt_advec/3600.0
              tmp_gsun = L%g_sun
              if ( SumVPD(i,j,ivpd) > 8.0 ) L%g_sun = min( L%g_sun, old_gsun(i,j,ivpd) )
              if( DEBUG_STOFLUX .and. debug_flag ) then
                 call datewrite("StoFlux SUMVPD ", iL, &
                    (/ real(ivpd), L%rh, L%t2C,  L%vpd, SumVPD(i,j,ivpd) /) ) 
              !       old_gsun(i,j), tmp_gsun, L%g_sun , L%g_sto
              endif

              old_gsun(i,j,ivpd) = L%g_sun
          end if

         ! Flux in nmole/m2/s:

          Sub(iL)%FstO3 = L%cano3_nmole * rc_leaf/(rb_leaf+rc_leaf) * L%g_sun 

          if( DEBUG_STOFLUX .and. debug_flag ) then
            call datewrite("STOFLUX ", iL, (/ L%cano3_nmole, L%cano3_ppb /) )
          end if

! ======   CLOVER  ===========================================================
      ! For Clover we have a very special procedure, using O3 from grassland
      ! to scale the fluxes. As grassland is entered in Inputs.Landuse before
      ! clover we can assume that Fst and O3 for grassland are available and 
      ! correct
          if( LandType(iL)%is_clover) then

             Sub(iL)%FstO3 = Sub(iLC_grass)%cano3_ppb/Sub(iL)%cano3_ppb * Sub(iL)%FstO3

             if ( DEBUG_STOFLUX  .and. debug_flag ) then
                call datewrite("CLOVER ", iL, &
                   (/ Sub(iL)%FstO3, Sub(iLC_grass)%cano3_ppb, &
                       Sub(iLC_grass)%cano3_ppb/Sub(iL)%cano3_ppb /) )
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

            gv = 0.0
            gvcms = 0.0
            if( L%g_sto > 1.0e-10 ) then
              gv = 1.0/ (2.0 * Rb_cor(WES_O3) /(KARMAN*L%ustar) + 1.0/ ( L%g_sto * L%LAI ) )
              gvcms = gv

         ! Step 2: Convert to mole/m2/s and for H2O
         ! mmol2sm = 8.3144e-8*L%t2  for O3, plus factor 1.6 for H2O
         ! * 0.001 -> mole/m2/s
         ! ms2molm2s = 1.6*1.0e-3/(8.3144e-8*L%t2)
         
              gv = gv * 1.6e-3/(8.3144e-8*L%t2)
            end if
         
         ! Step 3: 
         ! Mass flux density is E x 0.018 kg/mole -> kg/m2/s
         ! 1 kg/m2 = 1mm
         ! ie ms2kgm2s  = 1.6*1.0e-3/(8.3144e-8*L%t2)
         ! ie ms2mm     = 1.6*1.0e-3/(8.3144e-8*L%t2)

            Sub(iL)%EvapTransp = 0.018 *  L%vpd/101.0  * gv

         !!L%g_sto * L%LAI  * 1.6e-3/(8.3144e-8*L%t2)   ! Evapo-transpiration
         !                 (rb_leaf/1.6 + 0.0224*(L%t2/273.0) 


          if ( DEBUG_STOFLUX .and. debug_flag ) then 
            call datewrite("STO ", iL, (/ L%LAI, L%g_sto, L%g_sun, u_hveg,&
                             Sub(iL)%cano3_ppb, Sub(iL)%FstO3, gvcms /) )
          end if

      end if

    end do LC_LOOP
  end subroutine Calc_StoFlux

!..............................................................................
end module StoFlux_ml
