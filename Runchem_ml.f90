
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
module RunChem_ml
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------

!
   implicit none
   private

    public :: runchem


!    logical, private, save :: my_first_call = .true.
    logical, private, save :: MYDEBUG         = .false.

contains

subroutine runchem()

!/ Definitions 
!hf
   use DryDep_ml, only : drydep
   use My_Aerosols_ml, only: My_MARS,My_EQSAM,AERO_DYNAMICS,INORGANIC_AEROSOLS&
                            ,RUN_MARS,RUN_EQSAM,ORGANIC_AEROSOLS  
   use Par_ml,           only : lj0,lj1,li0,li1  &
                        ,gi0, gj0, me,NPROC & !! for testing
                        ,ISMBEG, JSMBEG    !! for testing
   use ModelConstants_ml, only :  PPB, KMAX_MID, dt_advec, &
                            DEBUG_i, DEBUG_j,nstep    ! rv1.2

   use Setup_1d_ml,       only: setup_1d, &
                                setup_bio, setup_rcemis, reset_3d

   use Setup_1dfields_ml, only: first_call  &  !DEBUG, ncalls & ! IS here best?
      ,amk , rcemis, rcbio, xn_2d  ! DEBUG for testing
   use Aqueous_ml,        only: Setup_Clouds, prclouds_present, WetDeposition
   use Ammonium_ml,       only: Ammonium!hf, INORGANIC_AEROSOLS
   use OrganicAerosol_ml, only: OrganicAerosol!hf, ORGANIC_AEROSOLS
   use Chemsolver_ml,     only: chemistry
   use My_Timing_ml,      only: Code_timer, Add_2timing, tim_before, tim_after
   use DefPhotolysis_ml,  only: setup_phot
!DEBUG ONLY:
   use GenSpec_tot_ml
   use GenSpec_adv_ml
  use Chemfields_ml,     only: xn_adv  ! DEBUG XXXXX


!/ local

   integer :: i, j
   integer :: i_emep, j_emep ! EMEP coordinates - for testing
   integer :: dt_chem,nchem  !  Time-step for chemistry solver
   integer :: Niter          !  No. iterations used in 2step
   logical , save :: reset_chem = .true.  ! =>  more iterns. at start of chem.
   logical ::  debug_flag    ! =>   Set true for selected i,j

   integer ::  iday, ihour, ispec, n
   integer ::  advec,errcode

! ========================= STARTS HERE ================================
      
!      call check_dts(dt_advec,dt_chem)
 
	dt_chem = 60   !DAVETMP 60

!pw u3 check dt_chem 
        nchem = max(2, nint(dt_advec)/dt_chem )
        dt_chem = nint(dt_advec)/nchem
        if(dt_chem*nchem.ne.nint(dt_advec))&
           write(*,*)'WARNING: using wrong dt_chem!', dt_chem

	Niter = 2


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

                       i_emep = i + ISMBEG + gi0 - 2  ! EMEP coordinates
                       j_emep = j + JSMBEG + gj0 - 2  ! EMEP coordinates

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

                    if(errcode /= 0)call gc_abort(me,NPROC,"setup_phot error")
                     call Add_2timing(29,tim_after,tim_before,&
                                           "Runchem:1st setups")

                     call setup_rcemis(i,j)

                     if ( ORGANIC_AEROSOLS ) &
                      call OrganicAerosol(debug_flag)

                     call Add_2timing(30,tim_after,tim_before,"Runchem:2nd setups")


                     !==================================================
                     call chemistry(i,j,dt_chem,Niter,reset_chem)
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
                         write(6,*) "DEBUG_RUN me pre WetDep", me, prclouds_present
                      end if

                     if ( prclouds_present)  &
                        call WetDeposition(i,j)

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

