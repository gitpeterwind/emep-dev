c
cds - call to setemis added
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     (3) physical and chemical calculations begin
c
      subroutine phyche(numt)
      use My_Outputs_ml , only : NHOURLY_OUT
     &                   ,FREQ_SITE, FREQ_SONDE, FREQ_HOURLY
      use Sites_ml, only: siteswrt_surf, siteswrt_sondes
      use My_Timing_ml, only : Code_timer, Add_2timing, 
     &                         tim_before, tim_after  
      use Dates_ml,     only : date, add_dates,dayno, daynumber
      use Derived_ml, only :wdep,ddep,IOU_INST , 
     &                       DerivedProds,Derived
      use DryDep_ml,    only : drydep,init_drydep
      use Par_ml   ,    only : me, MAXLIMAX, MAXLJMAX
      use Met_ml ,      only : roa,z_bnd,z_mid,metint, psurf, cc3dmax,
     &                            zen,coszen,Idirect,Idiffuse
      use ModelConstants_ml , only : KMAX_MID
     &			, dt_advec    ! time-step for phyche/advection
     &			,nmax,nstep
     &                  ,current_date, END_OF_EMEPDAY   ! date-type and end of "EMEP" day (usually 6am)
      use Nest_ml,         only : readxn, wrtxn
      use Emissions_ml,    only : EmisSet  
      use Timefactors_ml,  only : NewDayFactors  
      use Chemfields_ml,   only : xn_adv,cfac,xn_shl
      use Radiation_ml,  only : SolarSetup,        !ds mar2005 - sets up radn params
     &                          ZenithAngle,       ! gets zenith angle
     &                          ClearSkyRadn,      ! Idirect, Idiffuse
     &                          CloudAtten         ! Idirect, Idiffuse
!ds  &                          SolBio2D       ! Idirect, Idiffuse
!ds  &                          ,zen, coszen    ! TMP test of RAD
      use Runchem_ml  , only : runchem   ! Calls setup subs and runs chemistry
      use Advection_ml, only: advecdiff,adv_int
      use Polinat_ml, only : polinat_out
      use GridValues_ml, only : i_glob, j_glob, gl, gb ! DEBUG RAD !TMP

      implicit none
      integer i,j,k,n  ! DEBUG
      logical, parameter :: DEBUG = .false.
      logical, save :: End_of_Day = .false.    ! ds rv1_9_23

c

cds
      integer  numt
	integer ndays
	real :: thour
c
c
        if ( DEBUG ) 
     &    write(6,*)"enters phyche proc:", me,  
     &        "  current_date:", current_date
c
c     start of inner time loop which calls the physical and
c     chemical routines.
c
c
	DO_OUTER: do nstep = 1,nmax

c	  istep = istep+1
c
c     Hours since midnight at any time-step
c	using current_date we have already nstep taken into account
c
	  thour = real(current_date%hour) + current_date%seconds/3600. 
     &			+ 0.5*dt_advec/60./60.
c	  thour = real(current_date%hour) + (nstep-0.5)*dt_advec/60./60.

          if ( DEBUG ) then
              write(6,*) "me, thour-CURRENT_DATE CHECK: ", me, thour,  
     &                   current_date%hour, current_date%seconds
c	since we check now for current_date we can omit the nstep-check
	      if (me .eq. 0.and.current_date%hour .eq. 12) then
		call dayno(current_date%month,current_date%day,ndays)
	        write(6,*) 'thour,ndays,nstep,dt',
     &                      thour,ndays,nstep,dt_advec
	      endif

          endif

	  if (me .eq. 0) write(6,"(a15,i6,f8.3,2i5)") 
     >        'timestep nr.',nstep,thour
c
          call readxn(current_date)

!        ==================
	call Code_timer(tim_before)

    	call EmisSet(current_date)
        call Add_2timing(15,tim_after,tim_before,"phyche:EmisSet")
!jej
          !ds rv1_9_16 call sumDerived(dt_advec)    !< =====  Should add these lines

          wdep(:,:,:,IOU_INST) = 0.
          ddep(:,:,:,IOU_INST) = 0.
cc!        ==================
c
c
c
c
c   date needed for the calculation of solar zenith angle
c
!     ds mar2005 - new system using elemental subroutines:

        call SolarSetup(current_date%year,current_date%month,
     &                     current_date%day,thour)

        call ZenithAngle(thour, gb, gl, zen, coszen )

        if( DEBUG .and. me == 0 ) then
          write(*,*) "ZenRad ", current_date%day, current_date%hour, 
     &        thour, gl(3,3),gb(3,3),zen(3,3),coszen(3,3)
        end if
        !ds mar2005 call SolBio2D(daynumber,cc3dmax(:,:,KMAX_MID),psurf)

        call ClearSkyRadn(psurf,coszen,Idirect,Idiffuse)

        call CloudAtten(cc3dmax(:,:,KMAX_MID),Idirect,Idiffuse)

        call Add_2timing(16,tim_after,tim_before,"phyche:ZenAng")

c
        !================
	call advecdiff
        call Add_2timing(17,tim_after,tim_before,"phyche:advecdiff")
        !================
c
        call Code_timer(tim_before)
c

         !/ See if we are calculating any before-after chemistry productions:

          !=============================
          if ( nstep == nmax ) call DerivedProds("Before",dt_advec)
          !=============================

          call Add_2timing(26,tim_after,tim_before,"phyche:MACHO-prod")

         !hf===================================
           call init_drydep()
         !===================================

         !=========================================================

          call runchem(numt)   !  calls setup subs and runs chemistry

          call Add_2timing(28,tim_after,tim_before,"Runchem")

          !=========================================================
                           

         !/ See if we are calculating any before-after chemistry productions:

          !=============================
          if ( nstep == nmax ) call DerivedProds("After",dt_advec)
          !=============================
c
	  call Code_timer(tim_before)
          !=============================
!hf	  call drydep()
          !=============================
          call Add_2timing(34,tim_after,tim_before,"phyche:drydep")

c
c
c	this output needs the 'old' current_date_hour
	call polinat_out
c
c	the following partly relates to end of time step - hourly output
c	partly not depends on current_date
c	=> add dt_advec to current_date already here
c

          !====================================
          current_date = add_dates(current_date,nint(dt_advec))
          !====================================

        !u7.4vg - move here to allow some d_2d variables to be set
        !         and then used in ascii printouts.

!jej          call SumDerived(dt_advec)

          !ds rv1_9_28: to be consisten with reset of IOU_DAY

          End_of_Day = (current_date%seconds==0.and.
     &                  current_date%hour==END_OF_EMEPDAY)
          if( End_of_Day .and. me == 0 ) then
              write(*,"(a20,2i4,i6)") "END_OF_EMEPDAY, Hour,seconds=",
     &          END_OF_EMEPDAY, current_date%hour,current_date%seconds
          endif


          call Derived(dt_advec,End_of_Day)  !rv1_9_28 change

	  if ( current_date%seconds == 0 ) then

              if ( modulo(current_date%hour, FREQ_SITE) == 0 ) 
     &                           call siteswrt_surf(xn_adv,cfac,xn_shl)

              if ( modulo(current_date%hour, FREQ_SONDE) == 0 ) 
     &            call siteswrt_sondes(xn_adv,xn_shl)

              if ( NHOURLY_OUT > 0 .and. 
     &             modulo(current_date%hour, FREQ_HOURLY) == 0 ) 
     &             call hourly_out()

	  end if

          call Add_2timing(35,tim_after,tim_before,"phyche:outs")


        !u7.4vg   call SumDerived(dt_advec)
	  call metint
	  call adv_int

c	  if (mod(nstep,nhelp).eq.0) call massbud

          call wrtxn(current_date)

          call Add_2timing(36,tim_after,tim_before,"phyche:ints")
c
	enddo DO_OUTER

	end
