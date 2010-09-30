! <SoilWater_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************!
!*
!*  Copyright (C) 2010 met.no
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
                             longitude => gl
 use Landuse_ml,        only : water_fraction
 use LocalVariables_ml, only: Grid
 use MetFields_ml,      only : SoilWater_deep, nwp_sea
 use ModelConstants_ml, only : USE_SOILWATER, DEBUG_SOILWATER
 use Par_ml,            only : limax, ljmax, MAXLIMAX, MAXLJMAX
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
!   DAM = D/MAM, min 0, max 1.0

    real, private, parameter ::   SoilMAM = 0.02  ! HIRLAM
    real, private, parameter ::   SoilDAM = 1.0 ! assumed fraction of
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
      real :: land, sumland, newsw, REW ! For SW tests
      integer :: old_day=-999, old_hour = -999,  hourloc

      if( DEBUG_SOILWATER .and. debug_proc ) write(*,*) "DEBUG_SW START: ", &
        current_date%day, current_date%hour, current_date%seconds, old_day
      ! We rest once per day, but need to loop through the cells to find
      ! the 3am reset point
      !if ( current_date%day      == old_day .and. .not. my_first_call ) return
      if ( current_date%seconds  /= 0       .and. .not. my_first_call ) return
      !old_day = current_date%day

     ! If NWP thinks this is a sea-square, but we anyway have land,
     ! the soil moisture might be very strange.  We search neighbouring
     ! grids and make a land-weighted mean SW

     !   fSW =  1.0

       !TEST if ( USE_SOILWATER == .true. ) return  ! FAKER
       !TEST fSW = 0.1
       !  return
        if ( USE_SOILWATER == .false. ) return

        do j = 1, ljmax
           do i = 1, limax
!TMP
              hourloc= mod(nint(current_date%hour+24*(1+longitude(i,j)/360.0)),24)
             if ( DEBUG_SOILWATER .and. debug_proc.and. i==debug_li.and.j==debug_lj ) then
               write(*,*) "CHECK_SWF", hourloc, " date ", current_date
             end if
!TMP          !if ( nwp_sea(i,j) .and. water_fraction(i,j) < 0.9 ) then
!TMP          if ( water_fraction(i,j) < 0.9 ) then
             if ( my_first_call ) hourloc = 3 ! fake to get started
             if ( hourloc /= 3  ) cycle  ! Only set one per day, at 3am
             newsw = 0.0
             sumland  = 0.0
             REW      = SoilWater_deep(i,j,1)/ SoilMAM ! HIRLAM specific!

             if ( REW < SoilDAM ) then
               fSW(i,j) = REW/SoilDAM  ! HIRLAM specific!
             else
               fSW(i,j) = 1.0
             end if
             if ( DEBUG_SOILWATER .and. debug_proc.and. i==debug_li.and.j==debug_lj ) then
               write(*,"(a,i4,f7.4,i4,2f12.4)") "RESET_SWF: ", &
                 daynumber, SoilWater_deep(i,j,1), hourloc, REW, fSW(i,j)
             end if

  !!!! WILL RE-VISIT LATER - with Peter's readneighbour routine
!TMP             do ii = -1, 1
!TMP             do jj = -1, 1
!TMP               ii2=max(1,i+ii)
!TMP               ii2=min(ii2, limax)
!TMP               jj2=max(1,j+jj)
!TMP               jj2=min(jj2, ljmax)
!TMP               land = (1-water_fraction(ii2,jj2) )
!TMP               sumland = sumland + land
!TMP               newsw = newsw + land*SoilWater_deep(ii2,jj2,1)
!TMP               write(*,"(a,2i5,5f12.4)") "SUBSWF: ", ii2, jj2,&
!TMP                 water_fraction(ii2,jj2), SoilWater_deep(ii2,jj2,1), newsw, land, sumland
!TMP             end do
!TMP             end do
!TMP             newsw = newsw/sumland

        end do
        end do
newsw = -999.9

      if ( DEBUG_SOILWATER .and. debug_proc ) then
         i = debug_li
         j = debug_lj
         hourloc= mod(nint(current_date%hour+24*(1+longitude(i,j)/360.0)),24)
         REW      = SoilWater_deep(i,j,1)/ SoilMAM ! HIRLAM specific!

             write(*,"(a,f7.4,i4,f7.4,i4,2f12.4,L8,f12.4)") "DEBUG_SWF: ", &
                 water_fraction(i,j), daynumber, SoilWater_deep(i,j,1), hourloc, REW, fSW(i,j), nwp_sea(i,j) !!, newsw
             
      end if

      my_first_call = .false.

   end subroutine Set_SoilWater
end module SoilWater_ml
