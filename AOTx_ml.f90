! <AOTx_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module AOTx_ml
  use Chemfields_ml, only : xn_adv, cfac
  use ChemSpecs_adv_ml, only : IXADV_O3
  use GridValues_ml, only : debug_li, debug_lj
  use Io_Progs_ml,   only : datewrite
  use LandDefs_ml,   only : LandType
  use Landuse_ml,    only : WheatGrowingSeason
  use LocalVariables_ml, only : L, Grid, Sub
  use MetFields_ml, only: zen
  use ModelConstants_ml, only : dt_advec, KMAX_MID &
     ,DEBUG_AOT  & 
     ,PPBINV ! 1.0e9, for conversion from mixing ratio to ppb
  use Par_ml, only : MAXLIMAX, MAXLJMAX, limax, ljmax
  use TimeDate_ml, only : current_date, effectivdaynumber
  implicit none
  private

  public :: Calc_AOTx          ! called from AddMosaicOutput, My_DryDep
  public :: Calc_GridAOTx      ! called from Derived_ml
  public :: setaccumulate_2dyear

! Limit of daylight zenith angle for AOTs
  integer, private,  parameter :: AOT_HORIZON  = 89         
  integer, public, parameter:: STARTMONTH_FOREST=4,ENDMONTH_FOREST=9&
                                ,STARTMONTH_CROPS=5,ENDMONTH_CROPS=7 ! EU only!




contains
 !=========================================================================
  subroutine Calc_AOTx(defn,iLC, X,aot, debug_flag, debug_txt ) !,accumulate_2dyear)
    character(len=*), intent(in) :: defn
    logical, intent(in) :: debug_flag
    character(len=*), intent(in), optional :: debug_txt
    integer, intent(in) :: iLC
    real, intent(in) :: X
    real, intent(out)    :: aot
    ! logical, intent(inout) :: accumulate_2dyear

    real    :: o3
    integer :: i,j



    i = Grid%i
    j = Grid%j

    aot = 0.0


    if (   LandType(iLC)%is_forest .and. &
         (   current_date%month<STARTMONTH_FOREST&
         .or.current_date%month>ENDMONTH_FOREST)   )  then    ! both EU and MM-EMEP
             ! accumulate_2dyear = .false.
             return
    end if

   !If night, or outside growing season, we simply exit with aot=0
    if ( defn == "EU" .and.  (current_date%hour < 9 .or. &
         current_date%hour > 21 )) then ! 8-20 CET, assuming summertime
          ! accumulate_2dyear = .false.
          return

    else if ( Grid%izen >= AOT_HORIZON ) then  !UN or MM use daylight
          ! accumulate_2dyear = .false.
          return
    end if

  ! the wheat growing season is based upon the latitude function
    if ( defn == "MM" ) then
      if (   LandType(iLC)%is_crop .and.  LandType(iLC)%is_iam .and. &
             WheatGrowingSeason(i,j) ==  0 ) then
             ! accumulate_2dyear = .false.
             return
      end if

      o3 = Sub(iLC)%cano3_ppb   ! canopy top O3 for MM

    else if ( defn == "EU" ) then
      if (   LandType(iLC)%is_crop .and. &
            ( current_date%month<STARTMONTH_CROPS&
              .or.current_date%month>ENDMONTH_CROPS) ) then
           return
      end if

      o3 = Grid%surf_o3_ppb 

    end if

  ! Temporary setup - accumulation period for NS Clover for Mills et al paper
!    if ( defn == "CV" ) then
!      if ( current_date%month < 5 .or.  &
!           current_date%month == 5 .and.  current_date%day < 15) .or.
!           current_date%month == 5 .and.  current_date%day < 15) .or.
!
!           ) then
!             ! accumulate_2dyear = .false.
!             return
!      end if
!    end if


    !========== Calculate ========================

    if ( o3>X ) then  ! Add AOT, scaling for time-fraction

        aot = (o3-X)  ! BUG-fix. dt_advec takes care of:* dt_advec/3600.0

    end if
    if ( DEBUG_AOT .and. debug_flag .and. present( debug_txt )) then
       call datewrite("XXX"//trim(debug_txt) // "defn" // trim(defn), iLC, (/ real(Grid%izen), X,  o3, aot /) )
    end if

  end subroutine Calc_AOTx

 !=========================================================================
 !
  function Calc_GridAOTx( X, debug_flag, debug_txt ) result (aot)

    !/-- Calcuates AOT values for input threshold. Daylight values calculated
    !    only, for zenith < AOT_HORIZON ( e.g. 89 )
    !    Only relevant in ozone models, so far....

    real, intent(in) :: X           ! Threshold in AOTX (ppb)
    logical, intent(in) :: debug_flag
    character(len=*), intent(in), optional :: debug_txt

    real, dimension(limax,ljmax) :: aot

    integer :: izen                    ! integer of zenith angle
    real :: o3, o3_ref                 ! Ozone (ppb) - needed if AOTs
    integer :: i, j 

    aot = 0.0

    !--------- ONLY April-Sept -------------------------
    if(   current_date%month<STARTMONTH_FOREST&
      .or.current_date%month>ENDMONTH_FOREST) return
    !--------- ONLY April-Sept -------------------------

    do i=1, limax
        do j=1, ljmax

           izen = max(1,int( zen(i,j) + 0.5))

           if ( izen < AOT_HORIZON ) then
                o3 = xn_adv(IXADV_O3,i,j,KMAX_MID) &
                     * cfac(IXADV_O3,i,j) * PPBINV

                aot(i,j) = max( o3 - X , 0.0 )   ! Definition of AOTs

           end if
        end do
    end do
    if ( DEBUG_AOT .and. debug_flag .and. present( debug_txt )) then
       o3_ref = xn_adv(IXADV_O3,debug_li,debug_lj,KMAX_MID) * PPBINV
       o3   = o3_ref * cfac(IXADV_O3,debug_li,debug_lj)
       write(*,"(a, 4i5,f7.1,3f10.3)") "CalcGridAOT:"//debug_txt,&
          current_date%month,current_date%day, &
          current_date%hour,current_date%seconds, &
               zen(debug_li,debug_lj), o3_ref, o3, aot(debug_li,debug_lj)
    end if

  end function Calc_GridAOTx

!=========================================================================
  subroutine setaccumulate_2dyear(name,accumulate_2dyear)

! We don't want the yearly output to accumulate over the whole year
!OLD text:
!AOTXXc: accumulated into yearly_output from May to July
!AOTXXf: accumulated into yearly_output from April to September
!D2_EUAOTXXWH: accumulated into yearly_output from May to July
!D2_EUAOTXXDF: accumulated into yearly_output from April to September
!D2_UNAOTXXWH: accumulated into yearly_output from May to July
!D2_UNAOTXXDF: accumulated into yearly_output from April to September
!D2_MMAOTXXWH: accumulated into yearly_output over growing season

     character(len=*), intent(in) :: name
     logical, intent(inout) :: accumulate_2dyear !flag to know when to
                                                  !accumulate d_2d (case "EXT")

     if( name=="D2_EUAOT30DF".or.&
          name=="D2_EUAOT40DF".or.&
          name=="D2_UNAOT30DF".or.&
          name=="D2_UNAOT40DF"    &
          )then
         if(   current_date%month<STARTMONTH_FOREST&
              .or.current_date%month>ENDMONTH_FOREST)then
            accumulate_2dyear=.false.
         endif
     endif



     if(name=="D2_EUAOT30WH".or.&
        name=="D2_EUAOT40WH".or.&
        name=="D2_UNAOT30WH".or.&
        name=="D2_UNAOT40WH"    &
          )then
         if(   current_date%month<STARTMONTH_CROPS&
              .or.current_date%month>ENDMONTH_CROPS)then
            accumulate_2dyear=.false.

         endif
     endif

  end subroutine setaccumulate_2dyear

end module AOTx_ml
