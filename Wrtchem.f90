 subroutine Wrtchem(numt)
   !-----------------------------------------------------------------------
   ! DESCRIPTION
   ! Writes out data fields as NetCDF and (if out_binary defined)
   ! as binary. The latter option will be removed soon (ds 18/3/2005).
   !
   ! HISTORY
   ! Wrtchem.f90 version created by ds to make use of new Derived
   ! field methods and Output_binary_ml.
   ! wrtchem routines originally from MACHO, modified by su to use IOU_
   ! type outputs.
   !
   ! 23/10/2002 - ds changes
   !    Introduced END_OF_EMEPDAY into ModelConstants_ml.
   !      -- for "EMEP" days which end between 0 and 6am we need to use the previous 
   !         day as the output name. We assume that this is wanted for days
   !         where the END_OF_EMEPDAY is less than or equal to 7am
   !      -- Jan_1st logical introduced to deal with this  special case
   !         since first output should occur just as Jan 2nd starts (e.g.
   !         at 6am on 2nd Jan), and Jan_1st also marks end of a year run.
   !         (For runs starting in other months, one partial write-out will
   !          occur at 6am of the 1st day, but this should be over-written
   !          as soon as a full day of data is available).
   !      -- End_of_Run logical introduced to help readability.
   ! 21/3/2005 - ds changes
   !    AOT outputs removed.
   !-----------------------------------------------------------------------
!ds New deriv system, 21/12/2003:
   use My_Derived_ml, only: NDERIV_2D, & !dsAscii3D
                            NDERIV_3D

   use Dates_ml,           only: nmdays             !ds-out
!ds New Deriv:
   use Derived_ml, only: IOU_INST, IOU_YEAR, IOU_MON, IOU_DAY, f_2d, d_2d &
                                ,f_3d, d_3d, nav_3d  & !dsAscii3D
                                ,ResetDerived
   use Io_ml   ,           only: IO_WRTCHEM
   use ModelConstants_ml , only: nprint,current_date, END_OF_EMEPDAY , & !dsAscii3d
                                  KMAX_MID
   use My_Outputs_ml,      only: NBDATES, wanted_dates_bi, & !dsAcii3D
                                 Ascii3D_WANTED
   use out_restri_ml,      only: to_out_restri    ! su - allows 3-h output 
   use Output_binary_ml,   only: Output_binary
   use Par_ml,             only: MAXLIMAX,MAXLJMAX,GIMAX,GJMAX ,limax,ljmax,me,&
                                  ISMBEG,JSMBEG   !dsOH
   implicit none

   integer, intent(in) ::  numt

   real, dimension(MAXLIMAX, MAXLJMAX)  ::  local_2d  ! copy of local array
   real, dimension(GIMAX, GJMAX)        ::  glob_2d   ! array for whole domain
   !!real, allocatable, dimension(:,:,:)  ::  glob_3d   ! 3-D array for whole domain
   integer msnr1, msnr2
   integer i,j,n,k
   character*30 outfilename
   integer nyear,nmonth,nday,nhour,nmonpr
   integer :: yy_out, mm_out, dd_out   !ds - after allowance for END_OF_EMEPDAY
   logical :: Jan_1st, End_of_Run

   nyear  = current_date%year
   nmonth = current_date%month
   nday   = current_date%day
   nhour  = current_date%hour

   dd_out = nday
   mm_out = nmonth
   Jan_1st    = ( nmonth == 1 .and. nday == 1 )
   End_of_Run = ( mod(numt,nprint) == 0       )

   if(me==0) then
          write(6,"(a12,i5,5i4)") "DAILY DD_OUT ", numt, nmonth, mm_out, nday, dd_out, nhour
   end if
   if ( END_OF_EMEPDAY  <= 7 ) then

         dd_out = nday - 1     ! only used for daily outputs
         if(me==0)write(6,"(a12,i5,5i4)") "DAILY SET ", numt, nmonth, mm_out, nday, dd_out, nhour 
         if ( dd_out == 0 ) then
             mm_out = nmonth - 1
             if ( nmonth == 1 ) mm_out = 12
             dd_out = nmdays( mm_out )  !  Last day of month
             if(me==0)write(6,"(a12,i5,4i4)") "DAILY FIX ", numt, nmonth, mm_out, nday, dd_out 
         end if

   end if
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

   !ds-out if (nhour ==  6) then
   if (nhour ==  END_OF_EMEPDAY ) then

     !ds write(outfilename,fmt='(''out_c'',i2.2,i2.2''.dat'')') nmonth,nday

     if ( numt > 1 .and. .not. Jan_1st ) then   !ds don't write out
       write(outfilename,fmt='(''out_c'',i2.2,i2.2''.dat'')') mm_out, dd_out
       if(me==0)write(6,"(a12,i5,4i4,i5,2x,a30)") "DAILY FILE ", numt, nmonth, &
                   mm_out,nday,dd_out, nhour, outfilename

       call Output_binary(IOU_DAY,outfilename)

     end if

     call ResetDerived(IOU_DAY)    ! ds reset even on 1st jan.

   end if

!su   end of run:

   !ds if (mod(numt,nprint) == 0) then

   if ( End_of_Run ) then

      !su   write the remaining part of outday for the hours 6-0 at end of run

       !ds-out write(outfilename,fmt='(''out_c'',i2.2,i2.2''.dat'')') nmonth,nday

       write(outfilename,fmt='(''out_c'',i2.2,i2.2''.dat'')') mm_out,dd_out
       call Output_binary(IOU_DAY,outfilename)

       !ds if(nmonth == 1 .and.nday == 1 .and.nhour == 0)then

       if( Jan_1st .and.nhour == 0)then  ! End of year !
            write(outfilename,fmt='(''outyear'',i4.4,''.dat'')') nyear-1

       else
            write(outfilename,fmt='(''outterm'',i4.4,''.dat'')') numt
       endif
       call Output_binary(IOU_YEAR,outfilename)

   end if    ! last numt check


!su   move the monthly output from the nn.ne.nold to this place
!caverage depositions and surface conc. calculated for each month. Hence, 
!cmonthly values are printed out and
!cnav and depositions are set to zero for a new month.

   if(nday == 1.and.nhour == 0)then

     nmonpr = nmonth-1
     if(nmonpr.eq.0)nmonpr=12

    !dsAcii3D =============================================================
    if( NDERIV_3D > 0) then
          !!!if (me  ==  0) allocate(glob_3d(GIMAX, GJMAX,KMAX_MID), stat=msnr2)
          msnr1 = 2000

          do n = 1, NDERIV_3D


            if( me == 0 ) then
!change output name to be different from input
!hf               write(outfilename,fmt='(a,a1,i2.2)')  trim( f_3d(n)%name ), ".",  nmonpr
               write(outfilename,fmt='(a,a5,i2.2)')  trim( f_3d(n)%name ), ".out.",  nmonpr
               open(IO_WRTCHEM,file=outfilename)
               write(IO_WRTCHEM,fmt="(4i4)")  ISMBEG, GIMAX+ISMBEG-1,&
                                              JSMBEG, GJMAX+JSMBEG-1 ! domain
            end if

            if( nav_3d(n,IOU_MON) == 0 ) then
                    write(IO_WRTCHEM,*) "ERRROROR"!
            else 

              do k = 1, KMAX_MID
                   local_2d(:,:) = d_3d(n,:,:,k,IOU_MON) / nav_3d(n,IOU_MON)  ! make 2-d copy
                   call local2global(local_2d,glob_2d,msnr1)

                   if (me  ==  0) then

                       write(*,*) "XF FINISH ", n, f_3d(n)%name, nav_3d(n,IOU_MON) 
        
                       do j=1,GJMAX
                          do i=1,GIMAX
                             write(IO_WRTCHEM,"(es10.3)") glob_2d(i,j)
                          end do
                       end do

                    end if  ! me loop
               end do ! k
            end if ! nav == 0

            if( me == 0 ) close(IO_WRTCHEM)
            
        end do 
     end if
    !dsAcii3D =============================================================

     write(outfilename,fmt='(''outmonth'',i4.4,''.dat'')') nmonpr
     call Output_binary(IOU_MON,outfilename)
     call ResetDerived(IOU_MON)
   endif

   end subroutine Wrtchem

