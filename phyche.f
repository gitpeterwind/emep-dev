c
cds - call to setemis added
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     (3) physical and chemical calculations begin
c
      subroutine phyche(numt)
      !ds New Deriv use My_Derived_ml, only :wdep,ddep,IOU_INST 
      use My_Outputs_ml , only : NHOURLY_OUT
     &                   ,FREQ_SITE, FREQ_SONDE, FREQ_HOURLY
      use Sites_ml, only: siteswrt_surf, siteswrt_sondes
      use My_Timing_ml, only : Code_timer, Add_2timing, 
     &                         tim_before, tim_after  
      use Dates_ml,     only : date, add_dates,dayno, daynumber
!ds New Deriv
      use Derived_ml, only :wdep,ddep,IOU_INST , 
     &                       DerivedProds,Derived
!ds  &                       SumDerived, DerivedProds,Derived
      use DryDep_ml,    only : drydep,init_drydep
      use Par_ml   ,    only : me, MAXLIMAX, MAXLJMAX
      use Met_ml ,      only : roa,z_bnd,z_mid,metint, psurf, cc3dmax  ! dsrv1_6_x for SolBio2D
      use ModelConstants_ml , only : KMAX_MID
     &			, dt_advec    ! time-step for phyche/advection
     &			,nmax,nstep
     &                  ,current_date   ! date-type
      use Nest_ml,         only : readxn, wrtxn
      use Emissions_ml,    only : EmisSet  
      use Timefactors_ml,  only : NewDayFactors  
      use Chemfields_ml,   only : xn_adv,cfac,xn_shl
      use Radiation_ml,  only : ZenAng,    ! gets zenith angle
     &                          SolBio2D   ! Idirect, Idiffuse
      use Runchem_ml  , only : runchem   ! Calls setup subs and runs chemistry
      use Advection_ml, only: advecdiff,adv_int
      use Polinat_ml, only : polinat_out

      implicit none
      integer i,j,k,n  ! DEBUG
      logical, parameter :: DEBUG = .false.
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
        call ZenAng(thour)
        call SolBio2D(daynumber,cc3dmax(:,:,KMAX_MID),psurf)

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
          call Derived(dt_advec)

	  if ( current_date%seconds == 0 ) then

              if ( modulo(current_date%hour, FREQ_SITE) == 0 ) 
     &                           call siteswrt_surf(xn_adv,cfac,xn_shl)

              if ( modulo(current_date%hour, FREQ_SONDE) == 0 ) 
     &            call siteswrt_sondes(xn_adv,xn_shl)
!!odin.3      &            call siteswrt_sondes(z_bnd,z_mid,roa,xn_adv,xn_shl)

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
