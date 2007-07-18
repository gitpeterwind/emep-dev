
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
module RunChem_ml
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
   use DryDep_ml, only : drydep
   use My_Aerosols_ml, only: My_MARS,My_EQSAM,AERO_DYNAMICS,INORGANIC_AEROSOLS&
                            ,RUN_MARS,RUN_EQSAM,ORGANIC_AEROSOLS, Aero_water
   use Par_ml,           only : lj0,lj1,li0,li1  &
                        ,gi0, gj0, me & !! for testing
                        ,IRUNBEG, JRUNBEG    !! for testing
   use ModelConstants_ml, only :  PPB, KMAX_MID, dt_advec, &
                                  nprint, END_OF_EMEPDAY, &
                            DEBUG_i, DEBUG_j,nstep, NPROC

   use Setup_1d_ml,       only: setup_1d, &
                                setup_bio, setup_rcemis, reset_3d

   use Setup_1dfields_ml, only: first_call  &
      ,amk , rcemis, rcbio, xn_2d  ! DEBUG for testing
   use Aqueous_ml,        only: Setup_Clouds, prclouds_present, WetDeposition
   use Ammonium_ml,       only: Ammonium!hf, INORGANIC_AEROSOLS
   use OrganicAerosol_ml, only: OrganicAerosol!hf, ORGANIC_AEROSOLS
   use Chemsolver_ml,     only: chemistry
   use My_Timing_ml,      only: Code_timer, Add_2timing, tim_before, tim_after
   use DefPhotolysis_ml,  only: setup_phot
   use TimeDate_ml,       only: current_date
!DEBUG ONLY:
   use GenSpec_tot_ml
   use GenSpec_adv_ml
   use Chemfields_ml,     only: xn_adv  ! DEBUG XXXXX
   use SeaSalt_ml,        only : SeaSalt_flux
   use My_Aerosols_ml,    only : SEASALT
   use CheckStop_ml,      only: CheckStop

!
   implicit none
   private

    public :: runchem


  logical, private, save :: MYDEBUG         = .false.

contains

subroutine runchem(numt)

!/ Definitions 
!/
   integer, intent(in) :: numt       !water

!/ local

   integer :: i, j
   integer :: i_emep, j_emep ! EMEP coordinates - for testing
   integer :: dt_chem,nchem  !  Time-step for chemistry solver
   integer :: Niter          !  No. iterations used in 2step
   logical , save :: reset_chem = .true.  ! =>  more iterns. at start of chem.
   logical ::  debug_flag    ! =>   Set true for selected i,j

   integer ::  iday, ihour, ispec, n
   integer ::  advec,errcode
   integer ::  nmonth,nday,nhour      !water
   logical ::  Jan_1st, End_of_Run    !water

! ========================= STARTS HERE ================================
      

	Niter = 2

!water
   nmonth = current_date%month
   nday   = current_date%day
   nhour  = current_date%hour

   Jan_1st    = ( nmonth == 1 .and. nday == 1 )
   End_of_Run = ( mod(numt,nprint) == 0       )
!water

 !......................................................................
 !..  Get emission rates (molec/cm2/s) - man-made and biogenic

 !                call setemis()  ! done in phyche

 !......................................................................
 !     Deposition velocities

!TMP !         call DepVel

!TMP !......................................................................

!              reset_chem = .true.  ! Gives more iterations at start of hour

              !.... ************************************************
              !.... ****   chemistry calls *************************

                errcode = 0
                do j = lj0, lj1
                   do i = li0, li1


                     call Code_Timer(tim_before)

                     !****** debug cell set here *******
                     debug_flag =  .false.   !ds rv1_9_23 for safety
                     if ( MYDEBUG ) then

                       i_emep = i + IRUNBEG + gi0 - 2  ! EMEP coordinates
                       j_emep = j + JRUNBEG + gj0 - 2  ! EMEP coordinates

                       debug_flag = ( i_emep == DEBUG_i .and. j_emep == DEBUG_j )
                     end if
                     !****** debug cell set here *******

                     call setup_1d(i,j)   ! g12 - now calculates izen

                     call Add_2timing(27,tim_after,tim_before,&
                                              "Runchem:setup_1d")

                     call Setup_Clouds(i,j)

                     call setup_bio(i,j)    ! new

                     call Add_2timing(28,tim_after,tim_before,&
                                         "Runchem:setup_cl/bio")

                     call setup_phot(i,j,errcode)

                     call CheckStop(errcode,"setup_photerror in Runchem") 
                     call Add_2timing(29,tim_after,tim_before,&
                                           "Runchem:1st setups")

                     call setup_rcemis(i,j)
!SeaS
                     if ( SEASALT ) &
                      call SeaSalt_flux(i,j)

                     if ( ORGANIC_AEROSOLS ) &
                      call OrganicAerosol(debug_flag)

                     call Add_2timing(30,tim_after,tim_before,"Runchem:2nd setups")


                     !==================================================
                     call chemistry(i,j)
                     !==================================================

                     call Add_2timing(31,tim_after,tim_before,&
                                               "Runchem:chemistry")

                     !????????????????????????????????????????????????????
                     !ALTERNATING DRYDEP AND EQ
!hf dec-2002 Add check that one and only one eq is chosen
                     if(mod(nstep,2) /= 0 )then !do eq first, then drydep

                        if ( INORGANIC_AEROSOLS ) call ammonium() 
                        if ( RUN_MARS )           call My_MARS(debug_flag)
                        if ( RUN_EQSAM )          call My_EQSAM(debug_flag) 

                        call DryDep(i,j)
                     else !do drydep first, then eq

                        call DryDep(i,j)
                        if ( INORGANIC_AEROSOLS ) call ammonium() 
                        if ( RUN_MARS )           call My_MARS(debug_flag)
                        if ( RUN_EQSAM )          call My_EQSAM(debug_flag) 
                     endif
                    !????????????????????????????????????????????????????

                     call Add_2timing(32,tim_after,tim_before,&
                                                 "Runchem:ammonium+Drydep")

                      if ( MYDEBUG .and. debug_flag  ) then
                         !write(6,*) "DEBUG_RUN me pre WetDep", me, prclouds_present
                           write(6,"(a20,2i3,i5,2es12.3)") "DEBUG_RUN me OH", &
                            current_date%day, current_date%hour,&
                             current_date%seconds, &
                             xn_2d(OH,20), xn_2d(PHNO3,20)

                      end if

                     if ( prclouds_present)  &
                        call WetDeposition(i,j)

            !water
                     if ( nhour == END_OF_EMEPDAY .or.  End_of_Run )    &
                        call Aero_water(i,j)
            !water                        

                     call reset_3d(i,j)

                     call Add_2timing(33,tim_after,tim_before,&
                                            "Runchem:post stuff")

                     first_call = .false.   !** end of first call **** !

                   end do ! j
                end do ! i

                 !.... ************************************************
                reset_chem = .false.  

   end subroutine runchem

!------------------------------------------------------------------------------

end module RunChem_ml

