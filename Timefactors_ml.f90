!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                    module Timefactors_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!.......................................................................
!**    DESCRIPTION:
!  Read in monthly and daily factors, -> fac_emm, fac_edd arrays
!  The input files are Monthly.sox, Daily.sox, Monthly.nox, etc.
!  Sets the day/night variation in day_factor
!
!  Now, time-zone stuff is set in Country_ml, not here.
!
!**    REVISION HISTORY:
!      24/1/01   ds: removed timezone to Countries_ml
!      23/1/01   ds: recoded as independant module from Emissions_ml
!      25/4/00   ds: recoded as module
!      8/2/99    ds: added timezone
!      4/12/98   su: minor change, only on bcast less
!      D. Simpson,    3/2/99
!_____________________________________________________________________________

  use Country_ml,   only : NLAND
  use My_Emis_ml,   only : NEMIS, EMIS_NAME
  use EmisDef_ml,   only : NSECTORS

  use Dates_ml, only:       & ! subroutine, sets:
                   date,           & ! date-type definition
                   nmdays, nydays, & ! days per month (12), days per year
                   STARTDAY,       & ! specifies start day of year (MON, etc.)
                   SUN               ! =7, index for Sunday
  use Io_ml, only : &
                   open_file       &   ! subroutine
                 , ios,  IO_TIMEFACS   ! i/o error number, i/o label
  implicit none
  private

  !-- subroutines:

  public :: NewDayFactors
  public :: timefactors

  !-- time factor stuff: 

  real, public, save, &
      dimension(NLAND,NSECTORS,NEMIS) :: timefac ! overall time-factor, 
                                                 ! calculated daily

  !/*   The main code should not need to know about the following, except
  !     we need to do a gc_send  and I prefer to keep the parallel stuff
  !     out of this module ... ds */

  real, public, save,  &
           dimension(NLAND,12,NSECTORS,NEMIS) :: fac_emm    ! Monthly
  real, public, save,  &
           dimension(NLAND, 7,NSECTORS,NEMIS) :: fac_edd    ! Daily

  real, public, save, dimension(NSECTORS,0:1):: day_factor  ! Day/night factor 

  real, public, save, dimension(NLAND)  :: timezone  ! time-diff from GMT

!hf u2:
  logical, private, parameter :: DEBUG = .false.

  !/** used for general file calls and gc routines below **/
  character(len=30), private :: fname   ! input filename

contains


  ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine timefactors(year)

   !.......................................................................
   !**    DESCRIPTION:
   !  Read in monthly and daily factors, -> fac_emm, fac_edd arrays
   !  The input files are Monthly.sox, Daily.sox, Monthly.nox, etc.
   !  Sets the day/night variation in day_factor
   !  Set the local time correction (diff. from GMT) for each country
   !
   !**    REVISION HISTORY:
   !      23/1/01   ds: recoded as independant module from Emissions_ml
   !      25/4/00   ds: recoded as module
   !      8/2/99    ds: added timezone
   !      4/12/98   su: minor change, only on bcast less
   !      D. Simpson,    3/2/99
   !.......................................................................

  !--arguments
  integer, intent(in) :: year

  !-- Outputs -  module's fac_emm, fac_edd, day_factor, etc.

  !-- local
  integer ::  inland, insec     ! Country and sector value read from femis
  integer ::  i, ic, isec, n, idd, iday, mm, mm2 ! Loop and count variables
  integer ::  iemis                        ! emission count variables

  integer :: weekday         ! 1=monday, 2=
  real    :: xday, sumfac    ! used in interpolation, testing

!/** Factor giving nighttime  emission ratio. 
! ** note this is hard-coded with NSECTORS=11. Checked in code

   real, parameter, dimension(NSECTORS) ::  & 
        DAY_NIGHT = (/      & 
                       1.0  &! 1. comm/res. combustion
                     , 0.8  &! 2. comm/res. combustion
                     , 0.8  &! 3. Industrial combustion
                     , 1.0  &! 4.
                     , 1.0  &! 5. 
                     , 0.5  &! 6. Solvent use
                     , 0.5  &! 7. Road transport
                     , 0.8  &! 8. Other transport
                     , 1.0  &! 9. 
                     , 1.0  &! 10. Agriculture
                     , 1.0  &! 11. Nature
                     /)
! MADE look-alike ...............................................................
!                      1.0 & ! High source
!                    , 1.0 & ! Low  source

  write(unit=6,fmt=*) "into timefactors.f "

  ios = 0
  if ( year  <  1990 .or. year > 2002  ) then
	print *,"ERROR:fix year< 1990 or > 2002"
	ios = 1
  endif
  if ( nydays < 365 )  then
	print *,"ERR:Call set_nmdays before timefactors?"
	ios = 1
  endif
  if ( NSECTORS /= 11 ) then
	print *,"Day-Night dimension wrong!"
	ios = 1
  endif
  if(ios>0) return

!  #################################
!  *** 1) Read in Monthly factors

   fac_emm(:,:,:,:) = 1.0

   do iemis = 1, NEMIS

       fname = "MonthlyFac." // trim ( EMIS_NAME(iemis) )
       call open_file(IO_TIMEFACS,"r",fname,needed=.true.)
       if ( ios /= 0 ) then
          print *,"ios error: Monthlyfac"
          return
       endif

       n = 0
       do 
           read(IO_TIMEFACS,fmt=*,iostat=ios) inland, insec, &
             (fac_emm(inland,mm,insec,iemis),mm=1,12)
           if ( ios <  0 ) exit    ! End of file
           if ( ios >  0 ) then     ! A different problem..
               print *,"Read error on Monthlyfac"
               return
           end if
           n = n + 1
       enddo
       close(IO_TIMEFACS)
       write(unit=6,fmt=*) "Read ", n, " records from ", fname 
  enddo  ! iemis


! #################################
!CCC*** 2) Read in Daily factors

  fac_edd(:,:,:,:) = 1.0
	ios = 0

  do iemis = 1, NEMIS

       fname = "DailyFac." // trim ( EMIS_NAME(iemis) )
       call open_file(IO_TIMEFACS,"r",fname,needed=.true.)
       if ( ios /= 0 ) then
	  print *,"ios error: Dailyfac"
	  return
	endif

       n = 0
	ios = 0
       do !gv while(.true.)
           read(IO_TIMEFACS,fmt=*,iostat=ios) inland, insec, &
            (fac_edd(inland,i,insec,iemis),i=1,7)
           if ( ios <  0 ) exit   ! End of file
           if ( ios >  0 ) then
		print *,"Read error on Dailyfac"
		return
	  endif

           n = n + 1
           !-- Sum over days 1-7
           xday =  sum( fac_edd(inland,1:7,insec,iemis) ) / 7.0
           if ( xday > 1.001 .or. xday < 0.999 ) then
             print *,"ERROR: Dailyfac - not normalised"
	     ios = 7
	     return
           endif
       enddo
       close(IO_TIMEFACS)
       write(unit=6,fmt=*) "Read ", n, " records from ", fname
  enddo  ! NEMIS

	ios = 0
! #######################################################################
!cccc  3) Normalise the monthly-daily factors. This is needed in order to
!         account for leap years (nydays=366) and for the fact that different
!         years have different numbers of e.g. saturdays/sundays. 
!            Here we execute the same interpolations which are later done
!         in "NewDayFactors", and scale efac_mm if necessary


  write(unit=6,fmt=*) "Time factor interpolation "
  write(unit=6,fmt=*) "for nmdays(2) = ", nmdays(2), " gives nydays= ", nydays

  do iemis = 1, NEMIS
       n = 0
       do isec = 1, NSECTORS
           do ic = 1, NLAND
             iday = 0
             sumfac = 0.0

             do mm = 1, 12     ! Jan - Dec
                do idd = 1, nmdays(mm)

                   iday = iday + 1
                   weekday = STARTDAY( year ) + iday -1
                   weekday = modulo(weekday,7)      ! restores day to 1--6
                   if ( weekday == 0 ) weekday = SUN  ! restores sunday

                   mm2 = mm + 1 
                   if( mm2  > 12 ) mm2 = 1      ! December+1 => January
                   xday = real(idd-1) /real(nmdays(mm))

                   sumfac = sumfac +                         &  ! timefac :
                      ( fac_emm(ic,mm,isec,iemis) +          &
                      xday * (fac_emm(ic,mm2,isec,iemis)     &
                              -fac_emm(ic,mm,isec,iemis) ) ) &
                   *  fac_edd(ic,weekday,isec,iemis)   

                end do ! idd
             end do ! mm

             sumfac = real(nydays)/sumfac    

             if ( sumfac  <  0.97 .or. sumfac  >  1.03 ) then
                 print "(a8,3i3,a6,f12.7)",  &
                 "Error?  ",iemis, isec, ic," with ",sumfac 
                  print *,"sumfac problem"
		  ios = 2
		  return
             end if 

             if ((sumfac /= 1.0) .and.  &   ! ds -why??
                ( sumfac  <  0.999 .or. sumfac  > 1.001 )) then
                 n = n+1
                 do mm = 1, 12
                    fac_emm(ic,mm,isec,iemis)  =  &
                       fac_emm(ic,mm,isec,iemis)*sumfac
                 end do ! mm
             end if

          end do ! ic
       enddo ! isec
       if ( n ==  0 ) &
            write(unit=6,fmt=*)  &
           "Correction not needed for iemis, sumfac = " ,iemis, sumfac
      enddo ! iemis

	ios = 0


!#########################################################################
!
! Day/night factors are set from parameter DAY_NIGHT in emisdef_ml
!     daytime = 2 - nightime :

  day_factor(:,0)  =  DAY_NIGHT(:)             ! Night
  day_factor(:,1) = 2.0 - day_factor(:,0)      ! Day

!     #################################

    write(unit=6,fmt=*) "End of subroutine timefactors"

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
  !
  !  Calculates the monthly and daily correction factors for each
  !  country, emission, and sector.  Called once per day.
  !
  !  Uses arays from eulemis.inc:
  ! 
  !     fac_emm(NLAND,NM,NSECTORS,NEMIS)    ! Jan - Dec.
  !     fac_edd(NLAND,7,NSECTORS,NEMIS)     ! Monday=1, Sunday=7
  !
  ! Outputs:
  !    real timefac(NLAND,NSECTORS,NEMIS)
  !
  !...........................................................................
  !nyear(1) - year of simulation 
  !...........................................................................

  type(date), intent(in) :: newdate
  integer :: isec          ! index over emission sectors
  integer :: iemis         ! index over emissions (so2,nox,..)
  integer :: iland         ! index over countries 
  integer :: nmnd, nmnd2   ! this month, next month.
  integer :: weekday,nday,n       ! 1=monday, 2=
  real    :: xday          ! used in interpolation

  nmnd  = newdate%month              ! gv: nmonth(1)
  nday=0
  do n=1,nmnd-1
     nday=nday+nmdays(n)
  enddo
  nday=nday+newdate%day

!pw  weekday = STARTDAY( newdate%year ) + nydays - 1 
  weekday = STARTDAY( newdate%year ) + nday - 1 
  weekday = modulo(weekday,7)        ! restores day to 1--6
  if ( weekday == 0 ) weekday = SUN  ! restores sunday 
    
!   Parameters for time interpolation

  nmnd  = newdate%month              ! gv: nmonth(1)
  nmnd2 = nmnd + 1                   ! ds Next month
  if( nmnd2 > 12 ) nmnd2 = 1         ! December+1 => January

  xday = real( newdate%day - 1 ) / real( nmdays(nmnd) )
 

!   Calculate monthly and daily factors for emissions 

    do iemis = 1, NEMIS
      do isec = 1, NSECTORS 
         do iland = 1, NLAND
 
             timefac(iland,isec,iemis) =   &
                ( fac_emm(iland,nmnd,isec,iemis) + xday &
                      * (fac_emm(iland,nmnd2,isec,iemis) &
                        -fac_emm(iland,nmnd,isec,iemis) ) ) &
               *    fac_edd(iland,weekday,isec,iemis) 
 
         enddo ! iland  
      enddo ! isec   
   enddo ! iemis 
 end subroutine NewDayFactors
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module Timefactors_ml
