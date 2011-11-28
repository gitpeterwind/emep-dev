! <Timefactors_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                    module Timefactors_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!.......................................................................
!  DESCRIPTION:
!  Calculates emission temporal variation.
!  Reads monthly and daily (GENEMIS) factors for all emissions from files
!  Monthly.sox, Daily.sox, Monthly.nox, etc., -> in fac_emm, fac_edd arrays 
!  For every day, calculates emission factor "timefac" per country, emission 
!  sector, emission component
!  
!  Sets the day/night emissions variation in day_factor
!
!  D. Simpson,    3/2/99-11 0 H. Fagerli, 2011
!_____________________________________________________________________________

  use ModelConstants_ml, only : iyr_trend
  use PhysicalConstants_ml, only : PI
  use CheckStop_ml, only : CheckStop
  use Country_ml,   only : NLAND
  use EmisDef_ml,   only : NSECTORS, NEMIS_FILE, EMIS_FILE, ISNAP_DOM
  use GridValues_ml    , only : i_fdom,j_fdom, debug_proc,debug_li,debug_lj
  use Met_ml,       only : Getmeteofield
  use ModelConstants_ml, only : MasterProc, DEBUG => DEBUG_EMISTIMEFACS
  use Par_ml,       only : MAXLIMAX,MAXLJMAX, me, li0, lj0, li1, lj1
  use TimeDate_ml,  only:            &  ! subroutine, sets:
                     date,           &  ! date-type definition
                     nmdays, nydays, &  ! days per month (12), days per year
                     day_of_week,&      ! weekday
                     day_of_year        ! day count in year
  use Io_ml,        only :            &
                     open_file,       & ! subroutine
                     check_file,       & ! subroutine
                     PrintLog,        &
                     ios,  IO_TIMEFACS  ! i/o error number, i/o label

  implicit none
  private

  !-- subroutines:

  public :: NewDayFactors
  public :: timefactors
  public :: DegreeDayFactors

  !-- time factor stuff: 

  real, public, save, &
     dimension(NLAND,NSECTORS,NEMIS_FILE) :: timefac ! overall emission 
                                                      ! timefactor 
                                                      ! calculated daily
  real, public, save,  &
     dimension(NLAND,12,NSECTORS,NEMIS_FILE) :: fac_emm  ! Monthly factors
 ! We keep track of min value for degree-day work
 !ds
  real, public, save,  &
     dimension(NLAND,NSECTORS,NEMIS_FILE) :: fac_min ! Min of Monthly factors
!hf
  real, public, save,  &
     dimension(12) :: fac_cemm  ! Change in monthly factors over the years

  real, public, save,  &
     dimension(NLAND, 7,NSECTORS,NEMIS_FILE) :: fac_edd  ! Daily factors

  real, public, save, dimension(NSECTORS,0:1):: day_factor  ! Day/night factor 

  ! Heating-degree day factor for SNAP-2. Independent of country:
  logical, public, save :: Gridded_SNAP2_Factors = .false.
  real, public, dimension (MAXLIMAX,MAXLJMAX), save :: gridfac_HDD
  real, private, dimension (MAXLIMAX,MAXLJMAX), save :: tmpt2

  ! Used for general file calls and mpi routines below

  character(len=30), private :: fname2   ! input filename - do not change 

contains


  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine timefactors(year)

   !.......................................................................
   !  DESCRIPTION:
   !  Read in monthly and daily factors, -> fac_emm, fac_edd arrays
   !  The input files are Monthly.sox, Daily.sox, Monthly.nox, etc.
   !  Sets the day/night variation in day_factor
   !
   !  D. Simpson,    3/2/99
   !.......................................................................

  !-- Input
  integer, intent(in) :: year

  !-- Outputs -  module's fac_emm, fac_edd, day_factor, etc.

  !-- local
  integer ::  inland, insec     ! Country and sector value read from femis
  integer ::  i, ic, isec, n, idd, iday, mm, mm2 ! Loop and count variables
  integer ::  iemis             ! emission count variables

  integer :: weekday            ! 1=monday, 2=tuesday etc.
  real    :: xday, sumfac       ! used in interpolation, testing
  character(len=100) :: errmsg
!hf
  real :: fracchange
  real, dimension(NLAND,NEMIS_FILE):: sumfacc !factor to normalize monthly changes                                                         


! Factor giving nighttime  emission ratio. 
! Note this is hard-coded with NSECTORS=11. 

   real, parameter, dimension(NSECTORS) ::  & 
        DAY_NIGHT = (/      & 
                       1.0  &! 1.  Power production
                     , 0.8  &! 2.  Comm/res. combustion
                     , 0.8  &! 3.  Industrial combustion
                     , 1.0  &! 4.  Non-industrial combustion
                     , 1.0  &! 5.  Processes
                     , 0.5  &! 6.  Solvent use
                     , 0.5  &! 7.  Road transport
                     , 0.8  &! 8.  Other transport
                     , 1.0  &! 9.  Waste
                     , 0.6  &! 10. Agriculture
                     , 1.0  &! 11. Nature
                     /)

  if (DEBUG) write(unit=6,fmt=*) "into timefactors "

   call CheckStop( nydays < 365, &
      "Timefactors: ERR:Call set_nmdays before timefactors?")

   call CheckStop(  NSECTORS /= 11 , &
      "Timefactors: ERR:Day-Night dimension wrong!")


!  #################################
!  1) Read in Monthly factors, and determine min value (for baseload)

   fac_emm(:,:,:,:) = 1.0
   fac_min(:,:,:) = 1.0
!hf
   fac_cemm(:) = 1.0
   fracchange=0.005*(iyr_trend -1990)
   fracchange=max(0.0,fracchange) !do not change before 1990
   fracchange=min(0.1,fracchange) !stop change after 2010 
                                  !equal 1.1/0.9=1.22 summer/winter change
   write(unit=6,fmt=*) "Change summer/winter ratio in SNAP1 by ", fracchange

   do mm=1,12
      !Assume max change for august and february
      fac_cemm(mm)  = 1.0 + fracchange * cos ( 2 * PI * (mm - 8)/ 12 )
      write(unit=6,fmt=*) "Change in emis mm, fac_cemm ", mm,fac_cemm(mm)
   enddo
!hf end


   do iemis = 1, NEMIS_FILE

       fname2 = "MonthlyFac." // trim ( EMIS_FILE(iemis) )
       call open_file(IO_TIMEFACS,"r",fname2,needed=.true.)

       call CheckStop( ios, &
        "Timefactors: IOS error in Monthlyfac")

       n = 0
       do 
           read(IO_TIMEFACS,fmt=*,iostat=ios) inland, insec, &
             (fac_emm(inland,mm,insec,iemis),mm=1,12)
           if ( ios <  0 ) exit     ! End of file

           fac_min(inland,insec,iemis) = minval( fac_emm(inland,:,insec,iemis) )
           if( DEBUG ) write(*,"(a,3i4,f8.3)") "MIN tfac ", &
                   inland,insec,iemis, fac_min(inland,insec,iemis)

           call CheckStop( ios, "Timefactors: Read error in Monthlyfac")

           n = n + 1
       enddo

       close(IO_TIMEFACS)

!hf Apply change in monthly factors for SNAP 1
          sumfacc(:,:)=0.0
          do ic = 1, NLAND
             do mm=1,12
                fac_emm(ic,mm,1,iemis)=fac_emm(ic,mm,1,iemis)*fac_cemm(mm)
                sumfacc(ic,iemis)=sumfacc(ic,iemis)+fac_emm(ic,mm,1,iemis)
             enddo
          enddo
!hf normalize
          do ic = 1, NLAND
             do mm=1,12
                fac_emm(ic,mm,1,iemis)=fac_emm(ic,mm,1,iemis)*12./sumfacc(ic,iemis)
             enddo
          enddo
!hf end
       if (DEBUG) write(unit=6,fmt=*) "Read ", n, " records from ", fname2 
   enddo  ! iemis


! #################################
! 2) Read in Daily factors

  fac_edd(:,:,:,:) = 1.0

  do iemis = 1, NEMIS_FILE

       fname2 = "DailyFac." // trim ( EMIS_FILE(iemis) )
       call open_file(IO_TIMEFACS,"r",fname2,needed=.true.)

       call CheckStop( ios, &
        "Timefactors: Opening error in Dailyfac")

       n = 0
       do
         read(IO_TIMEFACS,fmt=*,iostat=ios) inland, insec, &
             (fac_edd(inland,i,insec,iemis),i=1,7)
           if ( ios <  0 ) exit   ! End of file

           call CheckStop( ios, "Timefactors: Read error in Dailyfac")

           n = n + 1

           !-- Sum over days 1-7
           xday =  sum( fac_edd(inland,1:7,insec,iemis) ) / 7.0

           call CheckStop( xday > 1.001 .or. xday < 0.999, &
                "Timefactors: ERROR: Dailyfac - not normalised")

       enddo

       close(IO_TIMEFACS)
       if (DEBUG) write(unit=6,fmt=*) "Read ", n, " records from ", fname2

  enddo  ! NEMIS_FILE

! #######################################################################
! 3) Normalise the monthly-daily factors. This is needed in order to
!    account for leap years (nydays=366) and for the fact that different
!    years have different numbers of e.g. Saturdays/Sundays. 
!    Here we execute the same interpolations which are later done
!    in "NewDayFactors", and scale efac_mm if necessary.


  write(unit=6,fmt="(a,I6,a,I5)")" Time factors normalisation: ",nydays,' days in ',year 

  do iemis = 1, NEMIS_FILE
       n = 0
       do isec = 1, NSECTORS
           do ic = 1, NLAND
             iday = 0
             sumfac = 0.0

             do mm = 1, 12     ! Jan - Dec
                do idd = 1, nmdays(mm)

                   weekday=day_of_week (year,mm,idd)

                   if ( weekday == 0 ) weekday = 7  ! restores Sunday to 7

                   mm2 = mm + 1 
                   if( mm2  > 12 ) mm2 = 1          ! December+1 => January

                   xday = real(idd-1) /real(nmdays(mm))

                   sumfac = sumfac +                            &  ! timefac 
                      ( fac_emm(ic,mm,isec,iemis) +             &
                        ( fac_emm(ic,mm2,isec,iemis)            &
                         - fac_emm(ic,mm,isec,iemis) ) * xday ) &
                      * fac_edd(ic,weekday,isec,iemis)   

                end do ! idd
             end do ! mm

             sumfac = real(nydays)/sumfac    


              if ( sumfac < 0.97 .or. sumfac > 1.03 ) then
                 write(unit=errmsg,fmt=*) &
                   "Time-factor error! for ",iemis, isec, ic," sumfac ",sumfac
                 call CheckStop(errmsg)
              end if

             if ( sumfac < 0.999 .or. sumfac > 1.001 ) then
                 n = n+1
             ! Slight adjustment of monthly factors
                  do mm = 1, 12
                    fac_emm(ic,mm,isec,iemis)  =  &
                       fac_emm(ic,mm,isec,iemis) * sumfac
                  end do ! mm
             end if

          end do ! ic
       enddo ! isec

       if ( n ==  0 ) &
            write(unit=6,fmt=*)  &
           "Correction not needed for iemis, sumfac = " ,iemis, sumfac

      enddo ! iemis


!#########################################################################
!
! Day/night factors are set from parameter DAY_NIGHT in emisdef_ml
! daytime = 2 - nightime :

  day_factor(:,0)  =  DAY_NIGHT(:)             ! Night
  day_factor(:,1) = 2.0 - day_factor(:,0)      ! Day

!#################################

    if (DEBUG) write(unit=6,fmt=*) "End of subroutine timefactors"

    if (DEBUG ) then 
       print *, " test of time factors, UK: "
       do mm = 1, 12
           print "(i2,i6,f8.3,3f8.4)", mm, nydays, sumfac,  &
            fac_emm(27,mm,2,1), fac_edd(27,1,2,1), fac_edd(27,7,2,1)
       end do ! mm
       print *, " day factors traffic are", day_factor(7,0), day_factor(7,1)
    end if ! DEBUG

 end subroutine timefactors
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine NewDayFactors(newdate)
  
  ! Calculates the monthly and daily factors for emission temporal variation
  ! for each country, emission, and sector.  Called at midnight every day.
  !
  ! Uses arays: 
  !     fac_emm(NLAND,NM,NSECTORS,NEMIS_FILES)    ! Jan - Dec.
  !     fac_edd(NLAND,7,NSECTORS,NEMIS_FILES)     ! Monday=1, Sunday=7
  !
  ! Outputs:
  !    real timefac(NLAND,NSECTORS,NEMIS_FILES)
  !
  !...........................................................................
  ! nyear(1) - year of simulation 
  !...........................................................................

  type(date), intent(in) :: newdate
  integer :: isec           ! index over emission sectors
  integer :: iemis          ! index over emissions (so2,nox,..)
  integer :: iland          ! index over countries 
  integer :: nmnd, nmnd2    ! this month, next month.
  integer :: weekday        ! 1=Monday, 2=Tuesday etc.
  real    :: xday           ! used in interpolation
  integer :: yyyy,dd 

 !-----------------------------

   yyyy=newdate%year
   nmnd=newdate%month
   dd=newdate%day

   weekday = day_of_week(yyyy,nmnd,dd)
   if ( weekday == 0 ) weekday = 7  ! restores Sunday to 7

!   Parameters for time interpolation

    nmnd2 = nmnd + 1                   ! Next month
    if( nmnd2 > 12 ) nmnd2 = 1         ! December+1 => January

    xday = real( newdate%day - 1 ) / real( nmdays(nmnd) ) 

!   Calculate monthly and daily factors for emissions 

    do iemis = 1, NEMIS_FILE
      do isec = 1, NSECTORS 
         do iland = 1, NLAND

             timefac(iland,isec,iemis) =                           &
                ( fac_emm(iland,nmnd,isec,iemis)  +                &
                   ( fac_emm(iland,nmnd2,isec,iemis) -             &
                      fac_emm(iland,nmnd,isec,iemis ) ) * xday )   &
               *  fac_edd(iland,weekday,isec,iemis) 

         enddo ! iland  
      enddo ! isec   
   enddo ! iemis 

 end subroutine NewDayFactors
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 subroutine DegreeDayFactors(daynumber)

!.....................................................................
!**    DESCRIPTION:
!   Generally called with daynumber, and then reads the gridded degree-day
!   based factors for emissions.
!   If called with daynumber = 0, just checks existance of file. If not
!   found, can use default country-based (GENEMIS) factors.

    integer, intent(in) :: daynumber   ! daynumber (1,..... ,365)
   

    integer :: nstart, alloc_err
    logical :: my_first_call = .true.   ! DSFF
    integer :: dd_old = -1,  n
    integer,dimension(2)  :: ijloc   ! debug only 
    integer :: iii, jjj              ! debug only 
    real :: checkmax
    character(len=80) :: errmsg, units

!      Gridded_SNAP2_Factors = .false.
!      return

   !/ See if we have a file to work with....
    if ( daynumber == 0 ) then
      call check_file("DegreeDayFac.nc", Gridded_SNAP2_Factors,&
        needed=.false., errmsg=errmsg )
      if ( Gridded_SNAP2_Factors ) then
         call PrintLog("Found DEGREE-day factors", MasterProc)
      else
         call PrintLog("Not-found: DEGREE-day factors", MasterProc)
      end if
      return
    end if

    !===============================================
    if ( .not. Gridded_SNAP2_Factors )  return !
    !===============================================

   !/ We have a file, calculate every day ... .

    if (dd_old == daynumber) return   ! Only calculate once per day max
    dd_old= daynumber
!
!    call ReadField_CDF('DegreeDayFac.nc',"DegreeDayFac",&
!              gridfac_HDD,daynumber,interpol='zero_order',needed=.true.,debug_flag=DEBUG)

   ! DegreeDays have the same domain/grid as the met data, so we can use:

    call Getmeteofield('DegreeDayFac.nc',"DegreeDayFac1000",nrec=daynumber,ndim=2,&
         unit=units,validity=errmsg, field=gridfac_HDD(:,:))
       call CheckStop(errmsg=="field_not_found", "DegreeDay field not found:")
    call Getmeteofield('DegreeDayFac.nc',"temperature_24h",nrec=daynumber,ndim=2,&
         unit=units,validity=errmsg, field=tmpt2(:,:))
       call CheckStop(errmsg=="field_not_found", "DegreeDayT2 field not found:")

    if ( DEBUG ) then
       ijloc = maxloc( gridfac_HDD(li0:li1,lj0:lj1))
       iii = ijloc(1)+li0-1
       jjj = ijloc(2)+lj0-1
       checkmax = maxval( gridfac_HDD(li0:li1,lj0:lj1))

       !print "(a,i3, 4x, 3i4,4x,3i4)", "DEBUG GRIDFAC IJS", me,  li0, li1, MAXLIMAX, lj0, lj1, MAXLJMAX
       print "(a,2i4,3f10.2,20i4)", "DEBUG GRIDFAC MAx", me, daynumber, &
           checkmax, gridfac_HDD(iii,jjj), tmpt2(iii,jjj),  ijloc(1), ijloc(2), i_fdom(iii), j_fdom(jjj)
  
       if( me==3) then
           write(*,"(a,3i4,2f12.3)") "GRIDFACDAY ", daynumber, &
             i_fdom(iii), j_fdom(jjj),  tmpt2(iii,jjj), 0.001*gridfac_HDD(iii,jjj)
       end if
    end if
    gridfac_HDD = 0.001 * gridfac_HDD ! CRUDE and TMP


    if ( DEBUG .and. debug_proc ) then
       iii = debug_li
       jjj = debug_lj
       print *, "DEBUG GRIDFAC", me, daynumber, iii, jjj, tmpt2(iii,jjj),  gridfac_HDD(iii, jjj)
    end if


   end subroutine DegreeDayFactors

end module Timefactors_ml

