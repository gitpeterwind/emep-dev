!***********************************************************************
  subroutine hourly_out() !!  spec,ofmt,ix1,ix2,iy1,iy2,unitfac)
!***********************************************************************
!**    DESCRIPTION:
!       Calculates and 
!       Outputs hourly concentration (or met) values for a sub-set of the grid.
!
!**    REVISION HISTORY:
!      Extended to produce new file, Hourly.mmyy, every month, 10/5/01 ds
!      stop_test used instead of stop_all, su, 05/01
!      Extended for variable format, met, xn_adv or xn_shl, ds, and to use
!       Asc2D type 19/4/01
!      Corrected for ISMBEG, etc., su, 4/01
!      New, ds, 5/3/99
!
!*************************************************************************
!
!hf u2   use My_Runmode_ml,    only : stop_test, DEBUG
   use My_Derived_ml,    only : d_2d, D2_HMIX, IOU_INST  !u7.4vg 
!rv1.2 TMP     D2_VG_REF,D2_VG_1M, D2_VG_STO,D2_FX_REF,D2_FX_STO
   use My_Outputs_ml,    only : NHOURLY_OUT, &      ! No. outputs
                                 Asc2D, hr_out      ! Required outputs

   use Par_ml ,          only : MAXLIMAX,MAXLJMAX,GIMAX,GJMAX    &
                                ,li0,li1,lj0,lj1 &  ! u7.5vg FIX
                                ,me,ISMBEG,JSMBEG,limax,ljmax,NPROC
   use ModelConstants_ml,only : current_date,KMAX_MID
   use Chemfields_ml ,   only : xn_adv,xn_shl, cfac
   ! tmp use Dates_ml,         only : date ,add_dates
   use Met_ml,           only : t2,th, roa   !u7.4vg temp2m, th
   use GenSpec_shl_ml , only : NSPEC_SHL  ! Maps indices
   !u1 use GenSpec_maps_ml , only : MAP_ADV2TOT, MAP_SHL2TOT   ! Maps indices
   use GenChemicals_ml , only : species                    ! Gives names
   use GridValues_ml,    only :  i_glob, j_glob   ! Gives emep coordinates
   use Io_ml,            only : IO_HOURLY
   implicit none

   !*.. Components of  hr_out
   !*  character(len=3) :: type   ! "ADVp" or "ADVu" or "SHL" or "T2 "
   !*  integer          :: spec   ! Species number in xn_adv or xn_shl array
   !* character(len=12) :: ofmt   ! Output format (e.g. es12.4)
   !*  integer          :: ix1    ! bottom-left x
   !*  integer          :: iy1    ! bottom-left y
   !*  integer          :: ix2    ! upper-right x
   !*  integer          :: iy2    ! upper-right y
   !*  real             :: unitconv   !  conv. factor
   !*  real             :: max    ! max allowed value

   ! local variables
   logical, save     :: my_first_call = .true. ! Set false after file opened
   integer msnr                        ! Message number for gc_rsend
   real hourly(MAXLIMAX,MAXLJMAX)      ! Local hourly value  (e.g. ppb)
   real ghourly(GIMAX,GJMAX)           ! Global hourly value (e.g. ppb)
   real :: arrmax                      ! Maximum value from array
   real :: unit_conv                   ! Unit conversion (ppb ug etc.)
   integer, dimension(2) :: maxpos     ! Location of max value 
   integer i,j,ih,ispec,itot           ! indices
   integer ist,ien,jst,jen             ! start and end coords
   character(len=20) :: errmsg = "ok"  ! For  consistecny check
   character(len=20) :: name           ! For output file, species names
   character(len=4)  :: suffix         ! For date "mmyy"
   integer, save :: prev_month = -99   ! Initialise with non-possible month
!hf u2:
   logical, parameter :: DEBUG = .false.

   !tmp type(date)    :: nextop_time        ! Time for next output

   ! if ( my_first_call ) then
   !     hourly(:,:) = 0.0      ! Initialise (ljmax+1:MAXLJMAX, limax+1:LIMAX
   !                            !  would have done,  but this is simpler)
   ! else                        
   ! Mask the edges of the hourly array, so that we can use maxval later
   ! This makes the code a bit neater below, but costs some CPU time here,
   ! and in evaluating maxval over the whole MAXLIMAX*MAXLJMAX dimension.

        !u7.5vg FIX hourly(limax+1:MAXLIMAX,:) = 0.0
        !u7.5vg FIX hourly(1:limax,ljmax+1:MAXLJMAX) = 0.0
        hourly(:,:) = 0.0

   ! end if

   if(me == 0 .and. current_date%month /= prev_month ) then

        if ( prev_month > 0 ) close(IO_HOURLY)      ! Close last-months file

       !/.. Open new file for write-out

        write(suffix,fmt="(2i2.2)") current_date%month, &
                           modulo ( current_date%year, 100 )
        name = "Hourly" // "." // suffix
        open(file=name,unit=IO_HOURLY,action="write")
        prev_month = current_date%month
   end if


!......... Uses concentration/met arrays from Chem_ml or Met_ml ..................
!
!        real xn_adv(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID)
!        real cfac(NSPEC_ADV,MAXLIMAX,MAXLJMAX)
! or...
!        real xn_shl(NSPEC_ADV,MAXLIMAX,MAXLJMAX,KMAX_MID)
! or...
!        real temp2m(MAXLIMAX,MAXLJMAX)
!
!..........................................................................


   HLOOP: do ih = 1, NHOURLY_OUT

      msnr  = 3475 + ih
      ispec = hr_out(ih)%spec 

       OPTIONS: select case ( hr_out(ih)%type ) 
         case ( "ADVppbv" )
            !u1 itot = MAP_ADV2TOT( ispec )
            !u1 advected species follow short-lived in totals array
            itot = NSPEC_SHL + ispec 
            name = species(itot)%name
            unit_conv =  hr_out(ih)%unitconv
            do j = 1, ljmax
               do i = 1, limax
                  hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                                 * cfac(ispec,i,j) &    ! 50m->1m conversion
                                 * unit_conv            ! Units conv.
		if(DEBUG .and. i.eq.3.and.j.eq.3) then
                  print *,"HOURLY",me,ih,ispec,xn_adv(ispec,i,j,KMAX_MID),&
                     cfac(ispec,i,j), hr_out(ih)%unitconv
		end if
               enddo
            enddo

         case ( "ADVugm3" )
            !u1 itot = MAP_ADV2TOT( ispec )
            !u1 advected species follow short-lived in totals array
            itot = NSPEC_SHL + ispec 
            name = species(itot)%name
            unit_conv =  hr_out(ih)%unitconv * species(itot)%molwt
            do j = 1, ljmax
               do i = 1, limax
                  hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                                 * cfac(ispec,i,j) &     ! 50m->1m conversion
                                 * unit_conv       &     ! Units conv.
                                 * roa(i,j,KMAX_MID,1)   ! density.
		if(DEBUG .and. i.eq.3.and.j.eq.3) then
                  print *,"HOURLYug",me,ih,ispec,itot,xn_adv(ispec,i,j,KMAX_MID),&
                     cfac(ispec,i,j), hr_out(ih)%unitconv,species(itot)%molwt, &
                     roa(i,j,KMAX_MID,1)
		end if
               enddo
            enddo

          case ( "SHLmcm3" )        ! No cfac for short-lived species
            !u1 itot = MAP_SHL2TOT( ispec )
            itot = ispec 
            name = species(itot)%name
            do j = 1, ljmax
                do i = 1, limax
                     hourly(i,j) = xn_shl(ispec,i,j,KMAX_MID) &
                                    * hr_out(ih)%unitconv  ! Units conv.
		if(DEBUG .and. i.eq.3.and.j.eq.3) then
                  print *,"HOURLY",me,ih,ispec,xn_shl(ispec,i,j,KMAX_MID),&
                      hr_out(ih)%unitconv
		end if
                enddo
            enddo

          case ( "T2 " )        ! No cfac for short-lived species
            name = "T2"
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = t2(i,j) - 273.15     ! Skip Units conv.
            end forall

          case ( "th " )        ! No cfac for short-lived species
            name = "theta"      ! Potential tempeature
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = th(i,j,KMAX_MID,1)  ! Skip Units conv.
            end forall

          case ( "D2D" )        ! No cfac for short-lived species
            name = "D_2D"       ! Default

            !rv1.2 if( ispec == D2_HMIX) name = "Hmix"
            !rv1.2 if( ispec == D2_VG_REF) name = "Vg_Ref"
            !rv1.2 if( ispec == D2_VG_1M ) name = "Vg_1m"
            !rv1.2 if( ispec == D2_VG_STO) name = "Vg_Sto"
            !rv1.2 if( ispec == D2_FX_REF) name = "Flux_ref"
            !rv1.2 if( ispec == D2_FX_STO) name = "Flux_sto"

            !if ( me == 17 .and. ispec == D2_VG_REF ) then
            !   print "(a10,2es12.3)", "HRFLUXES", d_2d(ispec,2,6,IOU_INST), hr_out(ih)%unitconv
            !end if
            !if ( me == 17 .and. ispec == D2_FX_REF ) then
            !   print "(a10,2es12.3)", "HXFLUXES", d_2d(ispec,2,6,IOU_INST), hr_out(ih)%unitconv
            !end if
            !u7.5vg FIX forall ( i=1:limax, j=1:ljmax)
            !u7.5vgc forall ( i=li0:li1, j=lj0:lj1)

            forall ( i=1:limax, j=1:ljmax)
               ! hourly(i,j) = pzpbl(i,j)
               hourly(i,j) = d_2d(ispec,i,j,IOU_INST) * hr_out(ih)%unitconv
            end forall

          case DEFAULT 
             errmsg = "ERROR! Hourly_out: " // hr_out(ih)%type 

       end select OPTIONS 



      !/ Get maximum value of hourly array

       arrmax = maxval(hourly)
       if ( arrmax  >   hr_out(ih)%max ) then
            write(6,*) "Hourly value too big!: ", ih, hr_out(ih)%type, arrmax
            write(6,*) "Species : ", name," : ",  " ispec ", ispec
            write(6,*) "max allowed is : ",  hr_out(ih)%max
            write(6,*) "unitconv was   : ", hr_out(ih)%unitconv
            write(6,*) " me, limax, ljmax : ",   me, limax, ljmax 
            maxpos = maxloc(hourly)
            write(6,*) "Location is i=", maxpos(1), " j=", maxpos(2)
            write(6,*) "EMEP coords ix=", i_glob(maxpos(1)), " iy=", j_glob(maxpos(2))
            if ( hr_out(ih)%type == "ADV" ) then
              write(6,*) "hourly is ", hourly(maxpos(1),maxpos(2))
              write(6,*) "xn_ADV is ", xn_adv(ispec,maxpos(1),maxpos(2),KMAX_MID)
              write(6,*) "cfac   is ",   cfac(ispec,maxpos(1),maxpos(2))
            end if
            !call stop_all('hourly too big')
!f u2            errmsg = "hourly too big"
            call gc_abort(me,NPROC,"hourly too big")
       endif

      !/ Consistency check on first_call

!hf u2       if ( errmsg .ne. "ok" ) then
!hf u2             call stop_test(.false.,me,NPROC,99,errmsg)
!hf u2       end if


      !/ Send to ghourly

       call local2global(hourly,ghourly,msnr)

       if (me ==  0) then

            !....   write out for a sub-section of the grid:

            !6d: write(IO_HOURLY,"('Spec ',i3,' = ',a12,' Date:',i5,7i4,es12.5)")  &
            !6d:       ispec,  name                                                &
            !6d:      ,current_date%year,current_date%month,current_date%day       &
            !6d:      ,current_date%hour                                           &
            !6d:      ,hr_out(ih)%ix1, hr_out(ih)%ix2                              &
            !6d:      ,hr_out(ih)%iy1, hr_out(ih)%iy2, hr_out(ih)%unitconv
  
            !/** We need to correct for small run-domains and the asked-for
            !    output domain. We can only print out the intersection of
            !    these two rectangles.

            !/ In emep coordinates we have:

            ist = max(ISMBEG,hr_out(ih)%ix1)
            jst = max(JSMBEG,hr_out(ih)%iy1)
            ien = min(GIMAX+ISMBEG-1,hr_out(ih)%ix2)
            jen = min(GJMAX+JSMBEG-1,hr_out(ih)%iy2)

            write(IO_HOURLY,"('Spec ',i3,' = ',a12,' Date:',i5,7i4,es12.5)")  &
                  ispec,  name                                                &
                 ,current_date%year,current_date%month,current_date%day       &
                 ,current_date%hour                                           &
                 ,ist, ien, jst, jen,         &
                  unit_conv

            if ( DEBUG ) print *, "TTTHOUR ISTS", me, ist, ien, jst, jen 

            !/ Then in model coordinates we have:

            ist = max(1,hr_out(ih)%ix1-ISMBEG+1)
            jst = max(1,hr_out(ih)%iy1-JSMBEG+1)
            ien = min(GIMAX,hr_out(ih)%ix2-ISMBEG+1)
            jen = min(GJMAX,hr_out(ih)%iy2-JSMBEG+1)

            do i = ist,ien
              do j = jst,jen

                write(IO_HOURLY, fmt=hr_out(ih)%ofmt ) ghourly(i,j)

              end do ! j
            end do   ! i

       end if  ! me loop

      end do HLOOP

     if ( my_first_call ) then
          my_first_call = .false.
     end if

  end subroutine hourly_out
