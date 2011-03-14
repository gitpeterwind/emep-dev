! <SoilWater_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2010-2011 met.no
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
 use GridValues_ml,     only : debug_proc, debug_li, debug_lj, i_fdom, j_fdom,&
                             longitude => glon
 use Landuse_ml,        only : water_cover
 use LocalVariables_ml, only: Grid
 use Met_ml,            only : extendarea
 use MetFields_ml,      only : SoilWater_deep, nwp_sea, SoilWaterSource
 use ModelConstants_ml, only : USE_SOILWATER, DEBUG_SOILWATER
 use Par_ml,            only : limax, ljmax, MAXLIMAX, MAXLJMAX, me
 use TimeDate_ml,       only : current_date, daynumber


    implicit none

    public :: Set_SoilWater


!
!                        :-------------------
!                       /:
!                      / :
!                     /  :
!                    /   :
!                   /    :
!                  /     :
!  ----------------------:-----------------x
!                 0      D                MAM
!
!   MAM = 0.02 in HIRLAM - max possible water
!    or ca. 0.7 in IFS (organic soils, wetlands have v. high values)
!   DAM = D/MAM, min 0, max 1.0

    !real, private, parameter ::   SoilMAM = 0.02  ! HIRLAM
    real, private, save ::   SoilDAM  ! assumed fraction of
                            ! MAXM when decline begins, D on figure
                            ! DO NOT SET TO ZERO!
!    real, dimension(366), public, save :: SWP = 0.0  ! daily soil water potential
                                              ! in  MPa
   real,public, save, dimension(MAXLIMAX,MAXLJMAX) :: &
    fSW = 1.0    ! fSW= f(relative extractable water) =  (sw-swmin)/(swFC-swmin)


contains
   subroutine Set_SoilWater()
      integer :: i, j, ii, jj, ii2, jj2
      logical :: my_first_call = .true.
      logical :: mydebug
      real :: land, sumland,  REW ! For SW tests
      integer :: old_day=-999, old_hour = -999,  hourloc

     ! We will average SW over a 3x3 area, excluding nwp sea squares. This
     ! helps to avoid coastal effects, and in any case we only aim to capture
     ! large scale dryness. Soil water modelling is too tricky, and within a 
     ! grid square there may be many soil textures. We can never capture this
     ! properly, so try to get SW 'about-right'

      real, dimension(MAXLIMAX+4,MAXLJMAX+4) :: & ! extended areas for averaging
         xsw, xsea   ! soil water and sea-fraction values

      if( DEBUG_SOILWATER .and. debug_proc ) write(*,*) "DEBUG_SW START: ", &
        current_date%day, current_date%hour, current_date%seconds, old_day

      if (my_first_call ) then
        ! Need to assume some values when soil drying affects gsto
         SoilDAM = 0.0
         if( SoilWaterSource=="IFS")    SoilDAM = 0.3
         if( SoilWaterSource=="PARLAM") SoilDAM = 0.5
      end if
         
      ! We rest once per day, but need to loop through the cells to find
      ! the 3am reset point
      !if ( current_date%day      == old_day .and. .not. my_first_call ) return

      if ( current_date%seconds  /= 0       .and. .not. my_first_call ) return

     ! If NWP thinks this is a sea-square, but we anyway have land,
     ! the soil moisture might be very strange.  We search neighbouring
     ! grids and make a land-weighted mean SW

        if ( .not. USE_SOILWATER  ) return ! and fSW has been set to 1. at start

        call extendarea( SoilWater_deep(:,:,1), xsw, debug_flag=.true.)
        call extendarea( water_cover(:,:),   xsea, debug_flag=.true.)

        do j = 1, ljmax
           do i = 1, limax

             hourloc= mod(nint(current_date%hour+24*(1+longitude(i,j)/360.0)),24)

             mydebug = ( DEBUG_SOILWATER .and. debug_proc.and. i==debug_li.and.j==debug_lj ) 
             if ( mydebug ) write(*,*) "CHECK_SWF", hourloc, " date ", current_date

!TMP          !if ( nwp_sea(i,j) .and. water_cover(i,j) < 0.9 ) then
!TMP          if ( water_cover(i,j) < 0.9 ) then

             if ( my_first_call ) hourloc = 3 ! fake to get started
             if ( hourloc /= 3  ) cycle  ! Only set one per day, at 3am


     ! Take a 3x3 average of the land-weighted values for SW. Seems best not
     ! to "believe" NWP models too much for this param, and the variation in
     ! a grid is so big anyway. We aim at the broad effect.
     ! (Alternative might be to find max values?)

             REW = 0.0
             sumland  = 0.0
             do ii = -1, 1
             do jj = -1, 1
               ii2=i+ii+2  ! coord in extended array of thickenss 2
               jj2=j+jj+2
               if( xsw(ii2,jj2) > 1.0e-10 ) then ! have some SW to work with
                  land    = 1.0 - xsea(ii2,jj2)
                  sumland = sumland + land
                  REW = REW + land*xsw(ii2,jj2)
                 if ( mydebug ) then
                  write(*,"(a,2i5,8f12.4)") "SUBSWF: ", ii2, jj2,&
                    water_cover(i,j), xsea(ii2,jj2), SoilWater_deep(i,j,1),&
                      xsw(ii2,jj2), REW, land, sumland
                 end if
                end if
             end do
             end do

             if( sumland > 0.1 ) then
                REW = REW/sumland
             else
                REW = 1.0
             end if

             !REW      = SoilWater_deep(i,j,1) !!!!/ SoilMAM ! Now done in Met_ml

             if ( REW < SoilDAM ) then
               fSW(i,j) = REW/SoilDAM
             else
               fSW(i,j) = 1.0
             end if

             if ( mydebug ) then
               write(*,"(a,i4,f7.4,i4,2f12.4)") "RESET_SWF: ", &
                 daynumber, SoilWater_deep(i,j,1), hourloc, REW, fSW(i,j)
             end if

        end do
        end do

      if ( DEBUG_SOILWATER .and. debug_proc ) then
         i = debug_li
         j = debug_lj
         hourloc= mod(nint(current_date%hour+24*(1+longitude(i,j)/360.0)),24)
         REW   = SoilWater_deep(i,j,1) !done:/ SoilMAM

         write(*,"(a,f7.4,i4,f7.4,i4,2f12.4,L8,f12.4)") "DEBUG_SWF: ", &
           water_cover(i,j), daynumber, SoilWater_deep(i,j,1), hourloc,&
            REW, fSW(i,j), nwp_sea(i,j)
             
      end if

      my_first_call = .false.

   end subroutine Set_SoilWater
end module SoilWater_ml
