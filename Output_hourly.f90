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
   use My_Derived_ml,    only : d_2d, D2_HMIX, IOU_INST,IOU_HOUR,Deriv  !u7.4vg 
   use My_Outputs_ml,    only : NHOURLY_OUT, &      ! No. outputs
                                 NLEVELS_HOURLY, &  ! ds rv1_8_2 
                                 FREQ_HOURLY, &     ! No. hours between outputs
                                 Asc2D, hr_out      ! Required outputs

   use Par_ml ,          only : MAXLIMAX,MAXLJMAX,GIMAX,GJMAX    &
                                ,li0,li1,lj0,lj1 &  ! u7.5vg FIX
                                ,me,ISMBEG,JSMBEG,limax,ljmax,NPROC
   use ModelConstants_ml,only : current_date,KMAX_MID,DEBUG_i,DEBUG_j,identi
   use Chemfields_ml ,   only : xn_adv,xn_shl, cfac
   use Met_ml,           only : t2,th, roa, surface_precip   !u7.4vg temp2m, th
   use GenSpec_shl_ml ,  only : NSPEC_SHL  ! Maps indices
   use GenChemicals_ml , only : species                    ! Gives names
   use GridValues_ml,    only : i_glob, j_glob   ! Gives emep coordinates
   use Io_ml,            only : IO_HOURLY
   use Radiation_ml,     only : Idirectt, Idiffuse
   use NetCDF_ml,        only : Out_netCDF,Init_new_netCDF &
                                ,Int1,Int2,Int4,Real4,Real8 !Output data type to choose

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
   logical, save     :: debug_flag = .false.
   integer, save     :: i_debug, j_debug       ! Coords matching i,j
   integer msnr                        ! Message number for gc_rsend
   real hourly(MAXLIMAX,MAXLJMAX)      ! Local hourly value  (e.g. ppb)
   real ghourly(GIMAX,GJMAX)           ! Global hourly value (e.g. ppb)
   real :: arrmax                      ! Maximum value from array
   real :: unit_conv                   ! Unit conversion (ppb ug etc.)
   real :: g                           ! tmp - saves value of ghourly(i,j)
   integer, dimension(2) :: maxpos     ! Location of max value 
   integer i,j,ih,ispec,itot           ! indices
   integer :: ik                       ! Index for vertical level
   integer ist,ien,jst,jen             ! start and end coords
   character(len=20) :: errmsg = "ok"  ! For  consistecny check
   character(len=20) :: name,netCDFName ! For output file, species names
   character(len=4)  :: suffix         ! For date "mmyy"
   integer, save :: prev_month = -99   ! Initialise with non-possible month
   logical, parameter :: DEBUG = .false.
   integer :: NLEVELS_HOURLYih
   type(Deriv) :: def1 !for NetCDF
   real :: scale !for NetCDF
   integer ::CDFtype,nk,klevel,ist,jst,ien,jen !for NetCDF

    if ( my_first_call ) then

      !/ds rv1.6.2
      !/ Ensure that domain limits specified in My_Outputs lie within
      !  model domain. In emep coordinates we have:

        do ih = 1, NHOURLY_OUT

           hr_out(ih)%ix1 = max(ISMBEG,hr_out(ih)%ix1)
           hr_out(ih)%iy1 = max(JSMBEG,hr_out(ih)%iy1)
           hr_out(ih)%ix2 = min(GIMAX+ISMBEG-1,hr_out(ih)%ix2)
           hr_out(ih)%iy2 = min(GJMAX+JSMBEG-1,hr_out(ih)%iy2)
           hr_out(ih)%nk = min(KMAX_MID,hr_out(ih)%nk)

        end do ! ih
        if ( DEBUG ) then
           do j = 1, ljmax
              do i = 1, limax
                  if ( i_glob(i)==DEBUG_i .and. j_glob(j)==DEBUG_j) then
                       debug_flag = .true.
                       i_debug = i
                       j_debug = j
                       !print *, "DEBUG FOUNDIJ me ", me, " IJ ", i, j
                  end if
              end do
            end do
        end if ! DEBUG
        my_first_call = .false.
    end if  ! first_call

   !     hourly(:,:) = 0.0      ! Initialise (ljmax+1:MAXLJMAX, limax+1:LIMAX
   !                            !  would have done,  but this is simpler)
   ! else                        
   ! Mask the edges of the hourly array, so that we can use maxval later
   ! This makes the code a bit neater below, but costs some CPU time here,
   ! and in evaluating maxval over the whole MAXLIMAX*MAXLJMAX dimension.

        !u7.5vg FIX hourly(limax+1:MAXLIMAX,:) = 0.0
        !u7.5vg FIX hourly(1:limax,ljmax+1:MAXLJMAX) = 0.0

     hourly(:,:) = 0.0

    !end if  ! first_call

   if(me == 0 .and. current_date%month /= prev_month ) then

        if ( prev_month > 0 ) close(IO_HOURLY)      ! Close last-months file

       !/.. Open new file for write-out

        write(suffix,fmt="(2i2.2)") current_date%month, &
                           modulo ( current_date%year, 100 )
        name = "Hourly" // "." // suffix
        open(file=name,unit=IO_HOURLY,action="write")
        prev_month = current_date%month

        netCDFName ="out_hour" // "."// suffix // ".nc"
        call Init_new_netCDF(netCDFName,IOU_HOUR)

       !ds rv1.6.2: Write summary of outputs to top of Hourly file
       !  - remember - with corrected domain limits here

        write(IO_HOURLY,*) NHOURLY_OUT, " Outputs"
        write(IO_HOURLY,*) FREQ_HOURLY, " Hours betwen outputs"
        write(IO_HOURLY,*) NLEVELS_HOURLY, "Max Level(s)"    !ds rv1_8_2

        do ih = 1, NHOURLY_OUT
           write(IO_HOURLY,fmt="(a12,a8,a10,i4,5i4,a13,es12.5,es10.3)") hr_out(ih)
        end do

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
      NLEVELS_HOURLYih=hr_out(ih)%nk
   KVLOOP: do ik = KMAX_MID, KMAX_MID-NLEVELS_HOURLYih+1, -1

      msnr  = 3475 + ih
      ispec = hr_out(ih)%spec 
      name  = hr_out(ih)%name   !ds rv1.6.1 
      if ( debug_flag ) print *, "DEBUG OH ", me, ispec, name,  hr_out(ih)%type

       if(DEBUG ) print *, "INTO HOUR TYPE ",hr_out(ih)%type

   !----------------------------------------------------------------
   !ds rv1_8_2: Added possibility of multi-layer output. Specify
   ! NLEVELS_HOURLY here, and in hr_out defs use either:
   !
   !      ADVppbv to get surface concentrations (onyl relevant for
   !              layer k=20 of course - gives meaningless  number f
   !               or higher levels.
   !Or,
   !      BCVppbv to get grid-centre concentrations (relevant for
   !      all layers.
   !----------------------------------------------------------------


       OPTIONS: select case ( hr_out(ih)%type ) 
         case ( "ADVppbv" )
            itot = NSPEC_SHL + ispec 
            name = species(itot)%name
            unit_conv =  hr_out(ih)%unitconv
            forall ( i=1:limax, j=1:ljmax)
                  hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                                 * cfac(ispec,i,j) &    ! 50m->1m conversion
                                 * unit_conv            ! Units conv.
            end forall

         case ( "BCVppbv" )            !ds rv1_8_2 
            itot = NSPEC_SHL + ispec 
            name = species(itot)%name
            unit_conv =  hr_out(ih)%unitconv
            forall ( i=1:limax, j=1:ljmax)
                  hourly(i,j) = xn_adv(ispec,i,j,ik) & !BCV:KMAX_MID) &
                                 !BCV * cfac(ispec,i,j) &    ! 50m->1m conversion
                                 * unit_conv            ! Units conv.
            end forall
            if ( DEBUG ) print *, "K-level", ik, name, itot

         case ( "ADVugm3" )
            itot = NSPEC_SHL + ispec 
            name = species(itot)%name
            unit_conv =  hr_out(ih)%unitconv * species(itot)%molwt
            forall ( i=1:limax, j=1:ljmax)
                  hourly(i,j) = xn_adv(ispec,i,j,KMAX_MID) &
                                 * cfac(ispec,i,j) &     ! 50m->1m conversion
                                 * unit_conv       &     ! Units conv.
                                 * roa(i,j,KMAX_MID,1)   ! density.
            end forall

          case ( "SHLmcm3" )        ! No cfac for short-lived species
            itot = ispec 
            name = species(itot)%name
            forall ( i=1:limax, j=1:ljmax)
                     hourly(i,j) = xn_shl(ispec,i,j,KMAX_MID) &
                                    * hr_out(ih)%unitconv  ! Units conv.
            end forall

          case ( "T2_C   " )        ! No cfac for short-lived species
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = t2(i,j) - 273.15     ! Skip Units conv.
            end forall

          case ( "theta  " )        ! No cfac for short-lived species
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = th(i,j,KMAX_MID,1)  ! Skip Units conv.
            end forall

          case ( "PRECIP " )        ! No cfac for short-lived species
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = surface_precip(i,j)     ! Skip Units conv.
            end forall

          case ( "Idirect" )        !  Direct radiation (W/m2)
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = Idirectt(i,j)    ! Skip Units conv.
            end forall

          case ( "Idiffus" )        !  Diffuse radiation (W/m2)
            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = Idiffuse(i,j)    ! Skip Units conv.
            end forall

          case ( "D2D" )        ! No cfac for short-lived species

            forall ( i=1:limax, j=1:ljmax)
               hourly(i,j) = d_2d(ispec,i,j,IOU_INST) * hr_out(ih)%unitconv
            end forall

          case DEFAULT 
             errmsg = "ERROR-DEF! Hourly_out: " // hr_out(ih)%type 
             call gc_abort(me,NPROC,"ABORT! hourly type not found")

       end select OPTIONS 

	if(DEBUG .and. debug_flag ) then
             i = i_debug
             j = j_debug
             print *,"DEBUG-HOURLY-TH ",me,ih,ispec,hourly(i,j),&
                      hr_out(ih)%unitconv
        end if


     !ds rv1.6.2 ---- why needed?
        hourly(limax+1:MAXLIMAX,:) = 0.0
        hourly(1:limax,ljmax+1:MAXLJMAX) = 0.0

      !/ Get maximum value of hourly array

       arrmax = maxval(hourly)
       if ( arrmax  >   hr_out(ih)%max ) then
            write(6,*) "Hourly value too big!: ", ih, hr_out(ih)%type, arrmax
            write(6,*) "Species : ", name," : ",  " ispec ", ispec
            write(6,*) "max allowed is : ",  hr_out(ih)%max
            write(6,*) "unitconv was   : ", hr_out(ih)%unitconv
            write(6,*) " me, limax, ljmax, MAXLIMAX,MAXLJMAX : ",  me, &
                             limax, ljmax ,MAXLIMAX,MAXLJMAX
            maxpos = maxloc(hourly)
            write(6,*) "Location is i=", maxpos(1), " j=", maxpos(2)
            write(6,*) "EMEP coords ix=", i_glob(maxpos(1)), " iy=", j_glob(maxpos(2))
            write(6,*) "hourly is ", hourly(maxpos(1),maxpos(2))
            if ( hr_out(ih)%type == "ADV" ) then
              write(6,*) "xn_ADV is ", xn_adv(ispec,maxpos(1),maxpos(2),KMAX_MID)
              write(6,*) "cfac   is ",   cfac(ispec,maxpos(1),maxpos(2))
            end if

            call gc_abort(me,NPROC,"hourly too big")
       endif

!NetCDF hourly output
       def1%name=hr_out(ih)%name
       def1%unit=hr_out(ih)%unit
       def1%class=hr_out(ih)%type 
       ist = max(1,hr_out(ih)%ix1-ISMBEG+1)
       jst = max(1,hr_out(ih)%iy1-JSMBEG+1)
       ien = min(GIMAX,hr_out(ih)%ix2-ISMBEG+1)
       jen = min(GJMAX,hr_out(ih)%iy2-JSMBEG+1)
       nk = min(KMAX_MID,hr_out(ih)%nk)
       CDFtype=Real4 ! can be choosen as Int1,Int2,Int4,Real4 or Real8
       scale=1.
       if(nk==1)then !write as 2D
       call Out_netCDF(IOU_HOUR,def1,2,identi &
            ,1,1,hourly(:,:),1,scale,CDFtype,ist,jst,ien,jen)
       elseif(nk>1)then   !write as 3D
          klevel=KMAX_MID-ik+1
       call Out_netCDF(IOU_HOUR,def1,3,identi &
            ,1,1,hourly(:,:),1,scale,CDFtype,ist,jst,ien,jen,klevel)
       else
          !nk<1 : no output
       endif

      !/ Send to ghourly

       call local2global(hourly,ghourly,msnr)

       if (me ==  0) then

            !....   write out for a sub-section of the grid:

            !/** We need to correct for small run-domains and the asked-for
            !    output domain. We can only print out the intersection of
            !    these two rectangles.

            !ds!/ In emep coordinates we have:

            !dsist = max(ISMBEG,hr_out(ih)%ix1)
            !dsjst = max(JSMBEG,hr_out(ih)%iy1)
            !dsien = min(GIMAX+ISMBEG-1,hr_out(ih)%ix2)
            !dsjen = min(GJMAX+JSMBEG-1,hr_out(ih)%iy2)

            !ds rv1_8_2 Extra info:
            write(IO_HOURLY,"('Spec ',i3,' Var ',i2,' = ',2a12,'Lev: ',i2,' Date:',i5,3i3)")  &
                  ispec,  ih, name, hr_out(ih)%name                          &
                 ,ik  &  ! ds rv1_8_2 
                 ,current_date%year,current_date%month,current_date%day       &
                 ,current_date%hour !ds                                           &
                 !ds ,ist, ien, jst, jen,         &
                 !ds  unit_conv

            if ( DEBUG ) print *, "TTTHOUR ISTS", me, ist, ien, jst, jen 

            !/ In model coordinates we have:

            ist = max(1,hr_out(ih)%ix1-ISMBEG+1)
            jst = max(1,hr_out(ih)%iy1-JSMBEG+1)
            ien = min(GIMAX,hr_out(ih)%ix2-ISMBEG+1)
            jen = min(GJMAX,hr_out(ih)%iy2-JSMBEG+1)

            do i = ist,ien
              do j = jst,jen

                g = ghourly(i,j)
                if ( g /= 0.0 ) then
                    write(IO_HOURLY, fmt=hr_out(ih)%ofmt ) g
                else ! Save disc-space used by thousands of  0.00000
                    write(IO_HOURLY, fmt="(i1)" ) 0
                end if

              end do ! j
            end do   ! i

       end if  ! me loop

      end do KVLOOP
      end do HLOOP


  end subroutine hourly_out
