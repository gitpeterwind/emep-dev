 subroutine Wrtchem(numt)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   ! Writes out binary data fields, and ascii fields for AOT if
   ! present. 
   !
   ! HISTORY
   ! Wrtchem.f90 version created by ds to make use of new Derived
   ! field methods and Output_binary_ml.
   ! wrtchem routines originally from MACHO, modified by su to use IOU_
   ! type outputs.
   !-----------------------------------------------------------------------
   use My_Derived_ml, only: IOU_INST, IOU_YEAR, IOU_MON, IOU_DAY &
               ,NDERIV_2D, f_2d, d_2d

   use Derived_ml, only: ResetDerived  ! 6c d_2d
             !ODIN6  ,valin_day,valin_month
   use Io_ml   , only : IO_AOT
   use ModelConstants_ml , only : nprint,current_date
   use My_Outputs_ml, only: NBDATES, wanted_dates_bi
   use out_restri_ml, only: to_out_restri    ! su - allows 3-h output 
   use Output_binary_ml, only : Output_binary   ! ODIN6 replaces out_typ1_4
   use Par_ml   , only : MAXLIMAX,MAXLJMAX,GIMAX,GJMAX ,limax,ljmax,me
   implicit none

   integer, intent(in) ::  numt

   real, dimension(MAXLIMAX, MAXLJMAX)  ::  local_2d  ! copy of local array
   real, dimension(GIMAX, GJMAX)        ::  glob_2d   ! array for whole domain
   integer msnr1, msnr2,nmonpr
   integer i,j,n
   character*30 outfilename
   integer nyear,nmonth,nday,nhour,nmonpr

   nyear  = current_date%year
   nmonth = current_date%month
   nday   = current_date%day
   nhour  = current_date%hour
!
!   actualize identi:
!   identi(12) = nyear
!   identi(13) = nday + 100*nmonth
!   identi(14) = 100 * nhour
!
   !su - allow 3h output of all concentrations for restricted domain

    if ( to_out_restri ) then

       write(outfilename,fmt &
               ='(''outrestri'',i4.4,i2.2,i2.2,i2.2,''.dat'')') &
               nyear,nmonth,nday,nhour
       call outchem_restri(outfilename)
     endif

!su   possible actual array output for specified days and hours - outchem
!ds   now defined in wanted_dates_bi array in My_Outputs

     do n = 1, NBDATES
         if ( wanted_dates_bi(n)%month == nmonth .and. &
              wanted_dates_bi(n)%day   == nday   .and. &
              wanted_dates_bi(n)%hour  == nhour ) then

               write(outfilename,fmt &
                  ='(''outchem'',i4.4,i2.2,i2.2,i2.2,''.dat'')') &
                  nyear,nmonth,nday,nhour
               call Output_binary(IOU_INST,outfilename)

         end if
     end do


!su   daily output - outday

   if (nhour ==  6) then

     write(outfilename,fmt='(''out_c'',i2.2,i2.2''.dat'')') nmonth,nday

     call Output_binary(IOU_DAY,outfilename)
     call ResetDerived(IOU_DAY)

   end if

!su   end of run:

   if (mod(numt,nprint) == 0) then

      !su   write the remaining part of outday for the hours 6-0 at end of run

       write(outfilename,fmt='(''out_c'',i2.2,i2.2''.dat'')') nmonth,nday
       call Output_binary(IOU_DAY,outfilename)

       if(nmonth == 1 .and.nday == 1 .and.nhour == 0)then
            write(outfilename,fmt='(''outyear'',i4.4,''.dat'')') nyear-1

       else
            write(outfilename,fmt='(''outterm'',i4.4,''.dat'')') numt
       endif
       call Output_binary(IOU_YEAR,outfilename)

   end if    ! last numt check

!su   aot each month or end of run
!ds ODIN6 - more flexible version?: Here we search all the derived 2d
!   fields. If the class is AOT we print-out. However, it is straightfoward
!   to test for other classes here also, e.g. we could add ACCSU also.
!
! The output file name is built up from the class and index of the 2d field,
!  e.g. AOT40.0589, AOT60.0590, or ACCSU-1.0590
!   ( and hope the index stays below 100 for the near future...)

! note - I introduced an extra array local_2d before passing to local2global.
! This is hopefully easier to understand, although not so efficient, as
! passing a portion of the 4d array
 
  if ( (mod(numt,nprint) == 0) .or. (nday == 1.and.nhour == 0) )then
    do n = 1, NDERIV_2D
      select case ( f_2d(n)%class )
        case ( "AOT" )

            msnr1 = 2000
            local_2d(:,:) = d_2d(n,:,:,IOU_MON)   ! make 2-d copy
            call local2global(local_2d,glob_2d,msnr1)
            if (me  ==  0) then

               if(nday == 1.and.nhour == 0)then
                   nmonpr = nmonth-1
                   if( nmonpr == 0) nmonpr=12
                   write(outfilename,fmt='(a,i2.2,i2.2)')  &  ! Build-up 
                                                              ! output name
                               trim( f_2d(n)%class ),  &
                                     f_2d(n)%index ,  &
                                     nmonpr
               else if ( mod(numt,nprint) == 0)then
                   write(outfilename,fmt='(a,i2.2,i4.4,i2.2,i2.2,i2.2)')  &  
                               trim( f_2d(n)%class ),  &
                                     f_2d(n)%index ,  &
                                     nyear,nmonth,nday,nhour
               endif

               open(IO_AOT,file=outfilename)
               do i = 1,GIMAX
                 do j = 1,GJMAX
                     write(IO_AOT,'(1x,i4,1x,i4,f17.6)') i,j,glob_2d(i,j)
                 end do
               end do
               close(IO_AOT)

             end if  ! me loop
            
       end select
     end do 
  end if    ! last numt check


!su   move the monthly output from the nn.ne.nold to this place
!caverage depositions and surface conc. calculated for each month. Hence, 
!cmonthly values are printed out and
!cnav and depositions are set to zero for a new month.

   if(nday == 1.and.nhour == 0)then

     nmonpr = nmonth-1
     if(nmonpr.eq.0)nmonpr=12
     write(outfilename,fmt='(''outmonth'',i4.4,''.dat'')') nmonpr
     call Output_binary(IOU_MON,outfilename)
     call ResetDerived(IOU_MON)
   endif

   end subroutine Wrtchem

