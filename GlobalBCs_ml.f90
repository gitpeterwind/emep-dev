! <GlobalBCs_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module GlobalBCs_ml
!
!+ DATA/SUBROUTINES FOR USING Logan climatology for BOUNDARY
!  CONDITIONS  (bcs)
!  From Hilde's uni.2_Logan_O3 but:
!  module name changed to GlobalBCs_ml, so that both UiO_ml and
!  public subroutine names changed back:
!        GetLoganData -> GetGlobalData   : UiO and Logan modules
!        (no L suffix now).
!
! Time-of-year, height, and latitude functions applied to non-oozne BCs
! Comments:
!   In using these functions I have had to accept some known problems,
!   for example that the seasonal variation is often weaker at higher
!   altitude. The "vmin" minimum concentration can correct for some
!   of this, but for now (this week) this simplification might be ok.
!   It is after all (I hope) the near-surface BCs which matter most.
!   In principal we could specify al concentrations as complex 3-D fieldds
!   here, but that can wait for the new setup!
!
! (the private routines, e.g. emeplat2Logan are left as they
!  are since these names are not seen outside this module.
!   -- for use with BoundConditions_ml and My_BoundConditions_ml --
!____________________________________________________________________________
  use My_Derived_ml,only: model

  use CheckStop_ml,      only: CheckStop
  use GridValues_ml, only: gbacmax,gbacmin,glacmax,glacmin&
                          ,gl,gb_glob,GlobalPosition &
                          ,i_fdom,j_fdom &
                          ,sigma_mid             ! for use in Hz scaling
  use Functions_ml, only: StandardAtmos_kPa_2_km ! for use in Hz scaling
  use GridValues_ml, only: lb2ij, xp, yp, AN, fi, gb_glob, gl_glob
  use Io_ml,             only : IO_GLOBBC, ios, open_file
  use ModelConstants_ml, only:  PPB, PPT,PPBINV &
                               ,KMAX_MID, PT &
                               ,MasterProc  &
                               ,IIFULLDOM,JJFULLDOM &
                               ,FORECAST     ! FORECAST mode

  use NetCDF_ml,only :GetCDF,Read_Inter_CDF
  use Par_ml, only:GIMAX,GJMAX,IRUNBEG,JRUNBEG, me
  use PhysicalConstants_ml, only: PI
  use TimeDate_ml,   only : daynumber


  implicit none
  private
  !/-- subroutines

  public :: GetGlobalData         ! Opens, reads bc_data, closes global data
  public :: setgl_actarray

   INCLUDE 'mpif.h'
   INTEGER STATUS(MPI_STATUS_SIZE),INFO

  logical, parameter, private :: DEBUG_GLOBBC= .false.
  logical, parameter, private :: DEBUG_Logan = .false.
  logical, parameter, private :: DEBUG_HZ    = .false.

  ! A. Define parameters and indices of global-model species
  ! ==========================================================================
  !-- definitions in Jostein's grid. Generally, these will be from
  !   a Txx model, where xx is currently 21.

  integer, parameter, public  :: &
  ! 15/3-2005 We assume BC defined in large domain:
         IGLOB = IIFULLDOM, &  ! number of large domain grids cells, longitude
         JGLOB = JJFULLDOM     ! number of large domain grids cells, latitude

  ! Chemical species:
  ! -- IBC indices text generated by perl script mkp.jost - ds
  ! ** usually only changed when global-model output changes **

  integer, public, parameter ::  NGLOB_BC  = 18 ! No. species setup here
  integer, public, parameter :: &
    IBC_O3       =  1   &
   ,IBC_NO       =  2   &
   ,IBC_NO2      =  3   &
   ,IBC_PAN      =  4   &
   ,IBC_HNO3     =  5   &  ! used for nitrate too
   ,IBC_SO2      =  6   &
   ,IBC_SO4      =  7   &
   ,IBC_CO       =  8   &
   ,IBC_C2H6     =  9   &
   ,IBC_C4H10    = 10   &
   ,IBC_HCHO     = 11   &
   ,IBC_CH3CHO   = 12   &
   ,IBC_H2O2     = 13   &
   ,IBC_aNH4     = 14   &
   ,IBC_aNO3     = 15   &
   ,IBC_pNO3     = 16   &
   ,IBC_CH3COO2  = 17   &
   ,IBC_OH       = 18


  ! we define some concentrations in terms of sine curves and other
  ! simple data:

   type, private :: sineconc
      real :: surf       ! Mean surface conc. (ppb)
      integer :: dmax    ! Day when concentrations peak
      real :: amp        ! amplitude of surface conc. (ppb)
      real :: hz         ! Scale-height (km) - height to drop 1/e concentration
      real :: vmin       ! background , minimum conc., in vertical direction
      real :: hmin       ! background , minimum conc., in horiz direction
      real :: conv_fac   ! factor to convert input data to mixing ratio
   end type sineconc
   type(sineconc), private, save, dimension(NGLOB_BC) :: SpecBC


! the actual values - do not use IGLOB,JGLOB, but the actual one's
  integer, save, private  :: iglbeg,iglend
  integer, save, private  :: jglbeg,jglend

  ! ==========================================================================

contains


  subroutine setgl_actarray(iglobact,jglobact)
 !set actual domain in 150*150 emep coord

  real hel1,hel2
  integer,intent(out)   :: iglobact,jglobact
    hel1 = IRUNBEG
    hel2 = IRUNBEG+GIMAX-1
    iglbeg = nint(hel1)  ! global i emep coord of start of domain
    iglend = nint(hel2)
    iglobact = GIMAX

    hel1 = JRUNBEG
    hel2 = JRUNBEG+GJMAX-1
    jglbeg = nint(hel1)
    jglend = nint(hel2)
    jglobact = GJMAX


 end subroutine setgl_actarray
 !-------
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine GetGlobalData(year,iyr_trend,month,ibc,used        &
                ,iglobact,jglobact,bc_data,io_num,errcode)

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   !== HANDLES READ_IN OF GLOBAL DATA. We read in the raw data from the
   !   global model, and do the vertical interpolation to EMEP k values
   !   here if the species is to be used.
   integer,             intent(in) :: year   !for Mace Head correction
   integer,             intent(in) :: iyr_trend !Allows future/past years
   integer,             intent(in) :: month
   integer,             intent(in) :: ibc    ! Index of BC, u3
   integer,             intent(in) :: used   ! set to 1 if species wanted
   integer,             intent(in) :: iglobact,jglobact !u3 intent added
   real, dimension(iglobact,jglobact,KMAX_MID), &
                       intent(out) :: bc_data   ! Data  from Logan model
   integer,            intent(out) :: io_num    !  i/o number
   integer,          intent(inout) :: errcode   !  i/o number
   logical, save :: my_first_call = .true.      ! u3

   real, dimension(IGLOB,JGLOB,KMAX_MID) :: bc_rawdata   ! Data (was rtcdmp)


! Now we want fake data on EMEP 50*50 grid....
   real:: USso2trend=1.0,USnoxtrend=1.0,USnh4trend=1.0

   integer, dimension(IGLOB,JGLOB), save :: lat5     !u3  for latfunc below

   real, dimension(NGLOB_BC,6:14), save  :: latfunc  !u3  lat. function
   real, save ::  twopi_yr, cosfac   ! for time-variations !u3
   real, dimension(12) :: macehead_O3
   real :: O3fix
   integer :: i, j , k,i1,j1,icount
   character(len=30) :: fname    ! input filename
   character(len=99) :: errmsg   ! For error messages
   character(len=30) :: BCpoll   ! pollutant name
   real :: trend_o3, trend_co, trend_voc
  !Use of standard atmosphere
   real, dimension(KMAX_MID), save :: p_kPa, h_km
   real :: scale_old, scale_new,iMH,jMH
   logical :: notfound !set true if NetCDF BIC are not found
   real, parameter :: macehead_lat = 53.3 !latitude of Macehead station
   real, parameter :: macehead_lon = -9.9 !longitude of Macehead station
   logical,parameter :: MACEHEADFIX=.true.

   io_num = IO_GLOBBC          ! for closure in BoundCOnditions_ml

!==================================================================
!Trends 1980-2003 derived from EPA emissions of so2,nox. nh4 derived from
!2/3so3+1/3nox
!Support for SO2 can be found in Hicks, Artz, Meyer and Hosker, 2002
!Figure 7 (Eastern US) which show 'close' correspondance between national
!emissions and concentration trend

!1920-1970 BCs derived from:
!NH4: nh3 emissions
!SOx: winter ice cores, Col du dome
!NOx: winter ice cores
!1890-1920: trender fra utslipp for SOx,NOx, NH3, Aardenne USA

if (iyr_trend == 1890) then
   USso2trend=0.12
   USnoxtrend=0.15
   USnh4trend=0.44
elseif (iyr_trend == 1900) then
   USso2trend=0.18
   USnoxtrend=0.20
   USnh4trend=0.48
elseif (iyr_trend == 1910) then
   USso2trend=0.27
   USnoxtrend=0.27
   USnh4trend=0.52
elseif (iyr_trend == 1920) then
   USso2trend=0.32
   USnoxtrend=0.33
   USnh4trend=0.59
elseif(iyr_trend == 1930)then
   USso2trend=0.35
   USnoxtrend=0.33
   USnh4trend=0.55
elseif(iyr_trend == 1940)then
   USso2trend=0.46
   USnoxtrend=0.25
   USnh4trend=0.59
elseif(iyr_trend == 1950)then
   USso2trend=0.59
   USnoxtrend=0.33
   USnh4trend=0.69
elseif(iyr_trend == 1960)then
   USso2trend=0.76
   USnoxtrend=0.5
   USnh4trend=0.76
elseif(iyr_trend == 1970)then
   USso2trend=0.95
   USnoxtrend=0.75
   USnh4trend=0.90
elseif(iyr_trend == 1980)then
   USso2trend=1.
   USnoxtrend=1.
   USnh4trend=1.
else if( iyr_trend == 1985) then
   USso2trend=0.91
   USnoxtrend=0.95
   USnh4trend=0.94
else if( iyr_trend == 1990) then
   USso2trend=0.89
   USnoxtrend=0.94
   USnh4trend=0.93
else if( iyr_trend == 1995) then
   USso2trend=0.72
   USnoxtrend=0.92
   USnh4trend=0.88
else if( iyr_trend == 1996) then
   USso2trend=0.71
   USnoxtrend=0.92
   USnh4trend=0.88
else if( iyr_trend == 1997) then
   USso2trend=0.73
   USnoxtrend=0.91
   USnh4trend=0.88
else if( iyr_trend == 1998) then
   USso2trend=0.73
   USnoxtrend=0.90
   USnh4trend=0.87
else if( iyr_trend == 1999) then
   USso2trend=0.68
   USnoxtrend=0.84
   USnh4trend=0.81
else if( iyr_trend == 2000) then
   USso2trend=0.63
   USnoxtrend=0.83
   USnh4trend=0.80
else if( iyr_trend == 2001) then
   USso2trend=0.62
   USnoxtrend=0.80
   USnh4trend=0.76
else if( iyr_trend == 2002) then
   USso2trend=0.59
   USnoxtrend=0.78
   USnh4trend=0.74
else if( iyr_trend >= 2003) then
   USso2trend=0.62
   USnoxtrend=0.77
   USnh4trend=0.74
else
   write(unit=errmsg,fmt=*) "Unspecified trend BCs for this year:", ibc, year
   call CheckStop(errmsg)

endif


!==================================================================
! Trends - derived from EMEP report 3/97
! adjustment for years outside the range 1990-2000.


  if ( iyr_trend >=  1990 ) then
         trend_o3 = 1.0
         trend_co = 1.0
         trend_voc= 1.0
  else
         trend_o3 = exp(-0.01*1.0 *(1990-iyr_trend))
         trend_co = exp(-0.01*0.85*(1990-iyr_trend)) ! Zander:CO
         trend_voc= exp(-0.01*0.85*(1990-iyr_trend)) ! Zander,1975-1990
  end if
  if (MasterProc .and. my_first_call) then
     write(6,"(a20,i5)") "GLOBAL TREND YEAR ",  iyr_trend
     write(6,"(a20,3f8.3)") "TRENDS O3,CO,VOC ", trend_o3, trend_co, trend_voc
  end if

!==================================================================
!=========== BCs Generated from Mace Head Data =======================
!
! Mace Head ozone concentrations for backgroudn sectors
! from Fig 5.,  Derwent et al., 1998, AE Vol. 32, No. 2, pp 145-157
!
! Here we use the meteorology year to get a reaslistic O3.
!      Later we use iyr_trend to adjust for other years, say for 2050.

! For 2010, 2020 "trend" runs  - use 10 yr average as base-O3
! then later scale by trend_o3:

 if ( iyr_trend /= year ) then

! For trends, use defaults from 1990-2000 average
!   macehead_O3 = (/  37.6, 40.0, 42.9, 43.2, 41.9, 33.9, &
!                     29.4, 30.1, 33.3, 36.5, 35.1, 37.8 /)

! For trends, use defaults from 1998-2007 average (Agnes, 01.06.2009)
   macehead_O3 = (/  40.2, 41.8, 45.1, 46.5, 43.5, 35.8, &
                     30.4, 30.0, 33.9, 37.3, 39.6, 38.9 /)

! 01.06.2009, Joffen and Dave suggested to drop the +3 (4.5) ppb scaling
! for future "trends". The 1998-2007 average will be used for these.
!    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    if ( iyr_trend == 2010   ) then
!       macehead_O3 = macehead_O3  + 3.0    !ASSUMPTION FOR IIASA SR runs
!    elseif ( iyr_trend > 2010   )then
!       macehead_O3 = macehead_O3  + 4.5    !ASSUMPTION FOR JUN06 SR runs
!    endif
!    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 else if( year == 1990) then
   macehead_O3 = (/    35.3,    36.3,    38.4,    43.0,    41.2,    33.4 &
    ,    35.1,    27.8,    33.7,    36.2,    28.4,    37.7/)
 else if( year == 1991) then
   macehead_O3 = (/    36.1,    38.7,    37.7,    45.8,    38.8,    36.3 &
    ,    29.6,    33.1,    33.4,    35.7,    37.3,    36.7/)
 else if( year == 1992) then
   macehead_O3 = (/    36.1,    37.3,    41.8,    39.6,    41.2,    31.5 &
    ,    28.3,    30.3,    31.3,    34.2,    36.1,    34.9/)
 else if( year == 1993) then
   macehead_O3 = (/    37.6,    40.4,    44.4,    42.6,    43.4,    29.2 &
    ,    28.5,    29.6,    32.2,     0.0,    37.3,    38.3/)
 else if( year == 1994) then
   macehead_O3 = (/    38.6,    37.3,    45.7,    43.8,    42.9,    35.1 &
    ,    30.8,    30.5,    33.8,    36.5,    34.0,    37.3/)
 else if( year == 1995) then
   macehead_O3 = (/    37.5,    37.1,    41.6,    42.4,    41.1,    33.1 &
    ,    29.1,    28.7,    33.7,    34.8,    35.0,    36.0/)
 else if( year == 1996) then
   macehead_O3 = (/    37.0,    40.1,    42.9,    44.6,    41.3,    38.3 &
    ,    29.3,    29.4,    35.6,    38.4,    37.8,    38.4/)
 else if( year == 1997) then
   macehead_O3 = (/    36.2,    41.9,    41.8,    40.4,    40.6,    34.4 &
    ,    26.2,    29.3,    31.3,    35.2,    25.7,    39.5/)
 else if( year == 1998) then
   macehead_O3 = (/    38.6,    42.0,    44.6,    45.1,    44.2,    33.0 &
    ,    29.7,    32.9,    35.7,    38.8,    39.7,    40.4/)
 else if( year == 1999) then
   macehead_O3 = (/    39.9,    44.5,    49.4,    45.0,    42.8,    34.3 &
    ,    29.0,    30.0,    31.8,    36.9,    39.6,    39.2/)
 else if( year == 2000) then
   macehead_O3 = (/    39.5,    42.1,    41.8,    43.8,    43.4,    34.5 &
    ,    28.0,    27.3,    33.6,    37.4,    35.6,    35.8/)
 else if( year == 2001) then
   macehead_O3 = (/    37.3,    38.0,    42.2,    44.8,    42.6,    34.9 &
        ,    28.9,    29.4,    29.9,    35.3,    37.3,    37.5/)
!---------------------------------------------------------------------------
! Preliminary BCs generated using Mace Head CFC and other greenhouse gases
! data to define clean air masses. Data cover all of 2002 and 9 months
! of 2003. What to do for Oct-Dec 2003?
! Could use (1) 2002 data or (2) 10-year average?
! Simmonds paper would support (1), simplicity (2).
! After seeing earlier 2003 plots, chose (2).
!
elseif ( year == 2002 ) then
   macehead_O3 = (/  42.4 ,   44.4 ,   45.5 ,     45.0 ,     45.9 ,     39.8 &
               ,     32.5 ,   28.7 ,   37.7 ,     39.3 ,     40.5 ,     42.3 /)

elseif ( year == 2003 ) then
  macehead_O3 = (/  39.8 ,   40.1 ,    44.7 ,     45.4 ,     45.7 ,     41.7 &
               ,     33.3 ,  31.0 ,    35.7,      37.9 ,     40.9 ,     38.1 /)

elseif( year == 2004) then
   macehead_O3 = (/ 40.8,    42.0,    48.3,    46.6,    39.9,    31.9 &
               ,    32.4,    32.1,    33.9,    36.7,    40.2,    39.8/)

elseif( year == 2005) then
   macehead_O3 = (/ 40.9,    41.4,    44.1,    45.6,    42.7,    32.9 &
               ,    26.7,    30.0,    33.2,    37.7,    39.5,    38.0/)

!  elseif ( year == 2006 )  then ! Derwent et al., (2007)
!    macehead_O3 = (/  42.1, 44.6, 48.7, 49.9, 44.5, 39.3, &
!                      32.9, 31.2, 36.7, 39.9, 42.8, 43.1 /)

! 2006 and 2007 are now calculated with 1998-2007 averages
! Modified by Agnes, 01.06.2009
 elseif ( year == 2006 )  then
   macehead_O3 = (/  39.8, 42.4, 44.2, 48.3, 41.3, 39.0, &
                     31.9, 29.5, 34.8, 37.4, 41.9, 39.9 /)

 elseif ( year == 2007 )  then
   macehead_O3 = (/  40.7, 38.2, 46.1, 46.4, 40.9, 34.5, &
                     31.2, 28.8, 33.3, 36.1, 40.6, 41.7 /)
!---------------------------------------------------------------------------
 elseif ( year > 2007 .and. &   ! use defoult Mace-Head correction
          .not.FORECAST ) then  ! (10-yr average) on FORECAST mode
   write(unit=errmsg,fmt=*) &
       "ERROR: No Mace Head correction for this year yet! ", year
   call CheckStop(errmsg)
!  else  ! Defaults, from 1990-2000 average !
!    macehead_O3 = (/  37.6, 40.0, 42.9, 43.2, 41.9, 33.9, &
!                      29.4, 30.1, 33.3, 36.5, 35.1, 37.8 /)
 else  ! Defaults, from 1998-2007 average (Agnes, 01.06.2009)
   macehead_O3 = (/  40.2, 41.8, 45.1, 46.5, 43.5, 35.8, &
                     30.4, 30.0, 33.9, 37.3, 39.6, 38.9 /)
 end if
!=========== Generated from Mace Head Data =======================

   errcode = 0
   errmsg = "ok"

 if ( DEBUG_Logan ) write(*,*) "DEBUG_LOgan ibc, mm", ibc, month

 ! ========= first call =========================================
   if ( my_first_call ) then
      ! Set up arrays to contain Logan's grid as lat/long
      !/ COnversions derived from emeplat2Logan etc.:

     twopi_yr = 2.0 * PI / 365.25

   call GlobalPosition  !get gb for global domaib
   do i = 1, IGLOB
     do j = 1, JGLOB    ! Don't bother with south pole complications

        lat5(i,j) = gb_glob(i,j)/5      ! lat/5 used in latfunc below
        lat5(i,j) = max(lat5(i,j),6)   ! Min value in latfunc
        lat5(i,j) = min(lat5(i,j),14)  ! Max value in latfunc
     end do
   end do

   ! Define concs where a simple  specification based on lat/mm
   ! etc. will be given
   !                           surf   dmax   amp   hz    vmin  hmin conv_fac!ref
   !                            ppb          ppb   km   hmin,vmin:same units as input data=conv_fac

  SpecBC(IBC_SO2  ) = sineconc( 0.15 , 15.0, 0.05, 999.9, 0.03 , 0.03,PPB) !W99, bcKz vmin
  SpecBC(IBC_SO4  ) = sineconc( 0.15 ,180.0, 0.00, 1.6,  0.05 , 0.03,PPB) !W99
  SpecBC(IBC_NO   ) = sineconc( 0.1  , 15.0, 0.03, 4.0  , 0.03, 0.02,PPB)
  SpecBC(IBC_NO2  ) = sineconc( 0.1  , 15.0, 0.03, 4.0  , 0.05, 0.04,PPB)
  SpecBC(IBC_PAN  ) = sineconc( 0.20 ,120.0, 0.15, 999.9, 0.20, 0.1 ,PPB)  !bcKz change vmin
  SpecBC(IBC_CO   ) = sineconc( 125.0, 75.0, 35.0,25.0  , 70.0, 30.0,PPB )!JEJ-W
  SpecBC(IBC_C2H6 ) = sineconc( 2.0  , 75.0, 1.0 , 10.0 , 0.05, 0.05,PPB )
  SpecBC(IBC_C4H10) = sineconc( 2.0  , 45.0, 1.0 , 6.0  , 0.05, 0.05,PPB )
  SpecBC(IBC_HCHO ) = sineconc( 0.7  ,180.0, 0.3 , 6.0  , 0.05, 0.05,PPB )
  SpecBC(IBC_CH3CHO) = sineconc(2.0 ,180.0, 0.5 , 6.0  , 0.05, 0.05,PPB )
!  SpecBC(IBC_HNO3 ) = sineconc( 0.1  , 15.0, 0.03, 999.9, 0.05, 0.05,PPB )!M, bcKz vmin
  SpecBC(IBC_HNO3 ) = sineconc( 0.07 , 180.0, 0.03, 999.9, 0.025, 0.03 ,PPB)!changed to be approx eq to NO3,
                                                                        !but with opposite seasonal var.
  SpecBC(IBC_aNO3  )= sineconc( 0.07 , 15.0, 0.03, 1.6,  0.025 , 0.02,PPB) !ACE-2
  SpecBC(IBC_pNO3  )= sineconc( 0.07 , 15.0, 0.00, 1.6,  0.025 , 0.02,PPB) !ACE-2
  SpecBC(IBC_aNH4  ) = sineconc( 0.15 , 180.0, 0.00, 1.6,  0.05 , 0.03,PPB) !ACE-2(SO4/NH4=1)

 ! all BCs read in are in mix. ratio, thus hmin,vmin needs to be in mix. ratio for those
  SpecBC(IBC_O3   ) = sineconc(-99.9 ,-99.9,-99.9,-99.9 ,-99.9,10.0*PPB,1.)!N1
  SpecBC(IBC_H2O2 ) = sineconc(-99.9 ,-99.9,-99.9,-99.9 ,-99.9,0.01*PPB,1.)
  SpecBC(IBC_OH   ) = sineconc(-99.9 ,-99.9,-99.9,-99.9 ,-99.9,1.0e-7*PPB,1.)
  SpecBC(IBC_CH3COO2)=sineconc(-99.9 ,-99.9,-99.9,-99.9 ,-99.9,1.0e-7*PPB,1.)

 ! Consistency check:

   if (DEBUG_GLOBBC) print *, "SPECBC NGLB ", NGLOB_BC
   do i = 1, NGLOB_BC
      if (DEBUG_GLOBBC) print *, "SPECBC i, hmin ",  i, SpecBC(i)%surf, SpecBC(i)%hmin

      if( SpecBC(i)%hmin*SpecBC(i)%conv_fac < 1.0e-17) then

         write(unit=errmsg,fmt=*) "PECBC: Error: No SpecBC set for species", i
         call CheckStop(errmsg)
      end if
   end do


 !refs:
 ! N1 - for ozone we read Logan's data, so the only paramater specified
 !      is a min value of 10 ppb. I hope this doesn't come into effect in
 !      Europe as presumably any such min values are in the S. hemisphere.
 !      Still, giving O3 such a value let's us use the same code for
 !      all species.
 ! W99: Warneck, Chemistry of the Natural Atmosphere, 2nd edition, 1999
 !    Academic Press. Fig 10-6 for SO2, SO4.
 ! JEJ - Joffen's suggestions from Mace/Head, UiO and other data..
 !      with scale height estimated large from W99, Isaksen+Hov (1987)
 ! M -Mozart-obs comparison
 ! ACE-2 Lots of conflicting measurements exist, from NH4/SO4=2 to NH4/SO4=0.5
 ! A 'mean' value of NH4/SO4=1 is therefore selected. Otherwise NH4 is assumed to
 ! act as SO4
 ! aNO3 is assumed to act as SO4, but with 1/2 concentrations and seasonal var.
 ! pNO3 is assumed to act like seasalt, with decreasing conc with height, with
 ! approx same conc. as fine nitrate.
 ! The seasonal var of HNO3 is now assumed to be opposite of aNO3.

   ! Latitude functions taken from Lagrangian model, see Simpson (1992)
   latfunc(:,6:14) = 1.0   ! default
                            !  30       40        50       60         70 degN
   latfunc(IBC_SO2,6:14) = (/ 0.05,0.15,0.3 ,0.8 ,1.0 ,0.6 ,0.2 ,0.12,0.05/)
   latfunc(IBC_HNO3,6:14)= (/ 1.00,1.00,1.00,0.85,0.7 ,0.55,0.4 ,0.3 ,0.2 /)
   latfunc(IBC_PAN,6:14) = (/ 0.15,0.33,0.5 ,0.8 ,1.0 ,0.75,0.5 ,0.3 ,0.1 /)
   latfunc(IBC_CO ,6:14) = (/ 0.6 ,0.7 ,0.8 ,0.9 ,1.0 ,1.0 ,0.95,0.85,0.8 /)

   latfunc(IBC_SO4,:) = latfunc(IBC_SO2,:)
   latfunc(IBC_NO ,:) = latfunc(IBC_SO2,:)
   latfunc(IBC_NO2,:) = latfunc(IBC_SO2,:)
   latfunc(IBC_HCHO,:) = latfunc(IBC_HNO3,:)
   latfunc(IBC_CH3CHO,:) = latfunc(IBC_HNO3,:)
   latfunc(IBC_aNH4,:) = latfunc(IBC_SO2,:)
   latfunc(IBC_aNO3,:) = latfunc(IBC_SO2,:)
   latfunc(IBC_pNO3,:) = latfunc(IBC_SO2,:)



  !Use Standard Atmosphere to get average heights of layers

    p_kPa(:) = 0.001*( PT + sigma_mid(:)*(101325.0-PT) ) ! Pressure in kPa
    h_km     = StandardAtmos_kPa_2_km(p_kPa)

    my_first_call = .false.

   end if ! my_first_call
 ! ========= end of first call ===================================


   !+
   !  Specifies concentrations for a fake set of Logan data.

   !/-- tmp and crude - we associate 1 km of scaleht with 1 Logan level
   !    - okayish for 1st 5-6 km in any case......

      fname = "none"           ! dummy for printout


      select case ( ibc )

           case  ( IBC_OH)
                write(unit=fname,fmt="(a6,i2.2)") "D3_OH.",month

             BCpoll='D3_OH'
             call ReadBC_CDF(BCpoll,month,bc_rawdata,IGLOB,JGLOB,KMAX_MID,notfound)
             if(notfound)then
               call open_file(IO_GLOBBC,"r",fname,needed=.true.,skip=1)
               if ( ios /= 0 ) errmsg = "BC Error H2O2"
               read(IO_GLOBBC,*) bc_rawdata
               close(IO_GLOBBC)
             endif

           case  ( IBC_CH3COO2)
              write(unit=fname,fmt="(a11,i2.2)") "D3_CH3COO2.",month

             BCpoll='D3_CH3COO2'
             call ReadBC_CDF(BCpoll,month,bc_rawdata,IGLOB,JGLOB,KMAX_MID,notfound)
             if(notfound)then
                call open_file(IO_GLOBBC,"r",fname,needed=.true.,skip=1)
                if ( ios /= 0 ) errmsg = "BC Error H2O2"
                read(IO_GLOBBC,*) bc_rawdata
                close(IO_GLOBBC)
             endif


            case  ( IBC_H2O2 )

              write(unit=fname,fmt="(a8,i2.2)") "D3_H2O2.",month

             BCpoll='D3_H2O2'
             call ReadBC_CDF(BCpoll,month,bc_rawdata,IGLOB,JGLOB,KMAX_MID,notfound)
             if(notfound)then
                call open_file(IO_GLOBBC,"r",fname,needed=.true.,skip=1)
                if ( ios /= 0 ) errmsg = "BC Error H2O2"
                read(IO_GLOBBC,*) bc_rawdata
                close(IO_GLOBBC)
             endif


            case  ( IBC_O3 )

              write(unit=fname,fmt="(a6,i2.2)") "D3_O3.",month

             BCpoll='D3_O3_Logan'
             call ReadBC_CDF(BCpoll,month,bc_rawdata,IGLOB,JGLOB,KMAX_MID,notfound)

             if(notfound)then
                call open_file(IO_GLOBBC,"r",fname,needed=.true.,skip=1)
                if ( ios /= 0 ) errmsg = "BC Error O3"
                read(IO_GLOBBC,*) bc_rawdata
                close(IO_GLOBBC)
             endif

              write(*,*) "dsOH READ OZONE3 ", fname, ": ", bc_rawdata(IGLOB/2,JGLOB/2,20)

       ! Mace Head adjustment: get mean ozone from Eastern sector

              O3fix=0.0
              icount=0
              if(MACEHEADFIX)then
              do j=1,JGLOB
                 do i=1,IGLOB
                    if(       gb_glob(i,j)<macehead_lat+20.0 &
                         .and.gb_glob(i,j)>macehead_lat-25.0 &
                         .and.gl_glob(i,j)<macehead_lon      &
                         .and.gl_glob(i,j)>macehead_lon-40.0)then
                       O3fix=O3fix+bc_rawdata(i,j,20)
                       icount=icount+1
                    endif
                 enddo
              enddo

! grid coordinates of Mace Head

              call lb2ij(macehead_lon,macehead_lat, iMH,jMH)

              if(icount>=1) O3fix = O3fix/(icount)  - macehead_O3(month)*PPB

              write(6,"(a10,2f7.2,i4,i6,3f8.3)") "O3FIXes ",iMH,jMH, &
                   month,icount,bc_rawdata(nint(iMH),nint(jMH),20)/PPB,&
                   macehead_O3(month),O3fix/PPB
              endif

              if (model=='ZD_ACID') then

                 bc_rawdata=max(15.*PPB,bc_rawdata)

              elseif (model=='ZD_OZONE') then
                 bc_rawdata=max(15.0*PPB,bc_rawdata-O3fix)
              else
                 call CheckStop("Problem with Mace Head Correction in GlobalBCs_ml")
              endif

              bc_rawdata=trend_o3 * bc_rawdata


            case ( IBC_NO, IBC_NO2, IBC_HNO3, IBC_CO, &
                     IBC_C2H6, IBC_C4H10, IBC_PAN, IBC_pNO3 )

             ! NB since we only call once per month we add 15 days to
             ! day-number to get a mid-month value

               cosfac = cos( twopi_yr * (daynumber+15.0-SpecBC(ibc)%dmax))

               bc_rawdata(:,:,KMAX_MID) =    SpecBC(ibc)%surf + &
                                       ( SpecBC(ibc)%amp * cosfac)

             !/ - correct for other heights
               do k = 1, KMAX_MID-1

                 scale_new = exp( -h_km(k)/SpecBC(ibc)%hz )

                 bc_rawdata(:,:,k) = bc_rawdata(:,:,KMAX_MID)*scale_new

                 if (DEBUG_HZ) then
                   scale_old = exp( -(KMAX_MID-k)/SpecBC(ibc)%hz )
                   write(6,"(a8,2i3,2f8.3,i4,f8.2,f8.3,2f8.3)") "SCALE-HZ ",&
                            month, ibc, SpecBC(ibc)%surf, SpecBC(ibc)%hz, k,&
                            h_km(k), p_kPa(k), scale_old, scale_new
                 end if ! DEBUG_HZ
               end do

               bc_rawdata = max( bc_rawdata, SpecBC(ibc)%vmin )


             !/ - correct for latitude functions --------------------------

               do i = 1, IGLOB
                  do j = 1, JGLOB
                     bc_rawdata(i,j,:) = bc_rawdata(i,j,:) * latfunc(ibc,lat5(i,j))
                  enddo
               end do
             !/ - end of correction for latitude functions ---------------

             !/ trend adjustments:
               if( ibc == IBC_C4H10 .or. ibc == IBC_C2H6 )then
                 bc_rawdata = trend_voc * bc_rawdata
               else if( ibc == IBC_CO )then
                 bc_rawdata = trend_co * bc_rawdata
               end if

          case ( IBC_SO2, IBC_SO4, IBC_HCHO, IBC_CH3CHO, IBC_aNH4, IBC_aNO3  )

              ! (No vertical variation for S in marine atmosphere, see W99)
              ! aNO3 and NH4 assumed to act as SO4
              !  and PAN is just temporary, with some guessing that
              !  since sources decrease with altitude, but lifetime
              !  increases the concs don't change much.

               cosfac = cos( twopi_yr * (daynumber+15.0-SpecBC(ibc)%dmax))

               bc_rawdata(:,:,:) =    SpecBC(ibc)%surf  + &
                                       ( SpecBC(ibc)%amp * cosfac)

             !/ - correct for latitude functions --------------------------
               if ( DEBUG_Logan ) write(*,*) "LOGAN HORIZ", &
                       ibc, SpecBC(ibc)%surf, cosfac
               do i = 1, IGLOB
                 do j = 1, JGLOB

                  bc_rawdata(i,j,:) = bc_rawdata(i,j,:) * latfunc(ibc,lat5(i,j))

                  if ( DEBUG_Logan ) then
                     write(6,"(2i4,f8.2,f15.4)")j, lat5(i,j), &
                         latfunc(ibc,lat5(i,j)), bc_rawdata(36,j,1)
                  end if
                 enddo
               end do

             !/ - end of correction for latitude functions ----------------

         case  default
             print *,"Error with specified BCs:", ibc
             errmsg = "BC Error UNSPEC"
          end select
          !================== end select ==================================
         if (DEBUG_GLOBBC) write(*,*) "dsOH FACTOR ", ibc, fname

         call CheckStop(errmsg)
   if( DEBUG_Logan )then

      write(*,"(a15,3i4,f8.3)") "DEBUG:LOGAN: ",ibc, used, month, cosfac
      write(*,*) "LOGAN BC MAX ", maxval ( bc_rawdata ), &
                         " MIN ", minval ( bc_rawdata )
      do k = KMAX_MID, 1, -1        ! print out for equator, mid-lat
         write(*, "(i4,f12.3)") k,  bc_rawdata(15,15,k)
      end do
   end if ! DEBUG


  !/ - correction for latitude functions -----------------------
   bc_rawdata = max( bc_rawdata, SpecBC(ibc)%hmin )


!  =========================================
   !trend data
   select case ( ibc )
   case(IBC_SO2, IBC_SO4)
      bc_rawdata=bc_rawdata*USso2trend
   case(IBC_aNH4)
      bc_rawdata=bc_rawdata*USnh4trend
   case(IBC_aNO3, IBC_pNO3, IBC_HNO3, IBC_NO2, IBC_NO, IBC_PAN)
      bc_rawdata=bc_rawdata*USnoxtrend
   end select
!  =========================================


   if( DEBUG_Logan ) write(*,*) "LOGAN NEW MIN", minval ( bc_rawdata )


   if ( used == 1 ) then

       bc_rawdata(:,:,:)=bc_rawdata(:,:,:) * SpecBC(ibc)%conv_fac !Convert to mixing ratio

        j1 = 1
        do j = jglbeg,jglend
          i1 = 1
          do i = iglbeg,iglend

          bc_data(i1,j1,:) =  bc_rawdata(i,j,:)
          i1 = i1+1
          enddo
          j1 = j1+1
        enddo


   end if

 end subroutine GetGlobalData
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   subroutine ReadBC_CDF(varname,month,bc_rawdata,varGIMAX,varGJMAX,varKMAX,notfound)

     logical, intent(out) :: notfound
     integer, intent(in) :: month
     integer, intent(in) :: varGIMAX,varGJMAX,varKMAX!dimensions of bc_rawdata
     real, dimension(*), intent(inout) :: bc_rawdata !NB written as one-dimensional
     character (len = *),intent(in) ::varname

     integer :: nstart,nfetch
     character (len = 100)::filename
     nstart=month
     nfetch=1
     fileName='Boundary_and_Initial_Conditions.nc'
     call GetCDF(varname,fileName,bc_rawdata,varGIMAX,varGJMAX,varKMAX,nstart,nfetch,needed=.false.)
     notfound=.false.
     if(nfetch==0)notfound=.true.
     if(notfound)then
        !interpolate from global file
        nstart=month
        nfetch=1
        fileName='GLOBAL_Boundary_and_Initial_Conditions.nc'
        call Read_Inter_CDF(fileName,varname,bc_rawdata,varGIMAX,varGJMAX,varKMAX,nstart,nfetch,interpol='bilinear',needed=.false.)
        notfound=.false.
        if(nfetch==0)notfound=.true.
     endif
   end subroutine ReadBC_CDF
end module GlobalBCs_ml

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
