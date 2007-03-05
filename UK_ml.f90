module UKdep_ml

use Dates_ml,       only: daynumber, nydays
use DepVariables_ml,only: NLANDUSE         &  ! No. UK land-classes
                      ,luname              &
                      ,crops, bulk, water  & ! logical variables
                      ,IAM_MEDOAK,IAM_WHEAT& !
                      ,LU_WATER, LU_ICE    & ! Pb210
                      ,forest,conif_forest & !    "      "
                      ,vegetation, urban   & !
                      ,luflux_wanted       & ! 
                      ,STUBBLE             & ! Ht. of stubble (m)
                      ,SAIadd              & ! surface-area index, rv1.4.6
                      ,hveg_max, b_inc, albedo, NH4_pl, SGS50, DSGS   &
                      ,EGS50, DEGS, LAImin, LAImax, SLAIlen, ELAIlen  &
                      ,f_phen_Slen , f_phen_Elen     &
                      ,f_phen_a, f_phen_b, f_phen_c, f_phen_d &
                      ,g_max     , f_min     , f_lightfac    &
                      ,f_temp_min, f_temp_opt, f_temp_max  &
                      ,RgsS      , RgsO        &
                      ,VPD_max   , VPD_min     &
                      ,SWP_max   , PWP       , rootdepth 
use Functions_ml,   only: GridAllocate, Polygon
use GridValues_ml,  only: gb_glob, gb, i_glob, j_glob, & ! latitude, coordinates
                            debug_proc, debug_li, debug_lj   !JUN06 
use Io_ml,          only: open_file, ios, IO_FORES
use ModelConstants_ml,  only : current_date, DEBUG_i, DEBUG_j, NNLANDUSE
use UKsetup_ml,     only: ukdep_init, get_growing_season, fPhenology
use Par_ml,         only: GIMAX, GJMAX, ISMBEG, JSMBEG, &
                          li0, lj0, IILARDOM, JJLARDOM, &
                          li1, lj1, MAXLIMAX, MAXLJMAX, &
                          me, NPROC, MSG_READ1, MSG_READ2, MSG_READ3
implicit none
private


!/- subroutines:

  public :: Init_ukdep
  public :: ReadLanduse
  public :: SetLanduse
  private :: Conif_fphen
 INCLUDE 'mpif.h'
 INTEGER STATUS(MPI_STATUS_SIZE),INFO

 integer, public, parameter :: NLUMAX = 17 ! max no. landuse per grid

 integer,public,save,dimension(MAXLIMAX,MAXLJMAX)        :: landuse_ncodes 
 integer,public,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: landuse_codes  
 real,   public,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: landuse_data 

 integer,public,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: &
          landuse_SGS   &    ! Start of growing season (days)
         ,landuse_EGS        ! End of growing season (days)

 integer,public,save,dimension(MAXLIMAX,MAXLJMAX) :: &
          InGrowingSeason   ! Growing season (days), IAM_WHEAT =1 for true
 
 integer,private,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: &
          Astart        &    ! JUN06 
         ,Aend               ! JUN06

 real,   public,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: &
          landuse_LAI   &    ! Leaf-area-index (m2/m2)
         ,landuse_hveg  &    ! Max. height of veg.
         ,landuse_fphen      ! Potential (age) factor for Jarvis-calc

 !ds 2006/06/25 - from Unimod.rv2_3_mmcV (2005/12/28)

 real,   public,save,dimension(MAXLIMAX,MAXLJMAX) :: &
             SumVPD ,   &   ! For critical VPD calcs, reset each day
             old_gsun       !

 !ds Pb210: Emissions from water (v.small) and from ice zero.
 real,public,save,dimension(MAXLIMAX,MAXLJMAX) :: water_fraction, ice_fraction 

 logical, private, parameter :: DEBUG_DEP = .false.
 character(len=30), private :: errmsg


contains

 !--------------------------------------------------------------------------
  subroutine Init_ukdep()
    character(len=20) :: errmsg
    integer :: i,j,lu    ! indices


  !/ 1./ -- Read in basic data for UK model

       !=====================================
        call ukdep_init(errmsg,me)
       !=====================================
        if ( errmsg /= "ok" ) then
           errmsg = "ukdep_init: " // errmsg
             WRITE(*,*) 'MPI_ABORT: ', errmsg 
             call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
        end if

end subroutine Init_ukdep
 !--------------------------------------------------------------------------
subroutine ReadLanduse()

! arrays for whole EMEP area ( allocated for me=0 only)

  integer, allocatable, dimension(:,:)      :: g_ncodes 
  integer, allocatable, dimension(:,:,:)    :: g_codes  
  real,    allocatable, dimension(:,:,:)    :: g_data 

   integer :: i,j,n,lu, index_lu, maxlufound
   integer :: err1, err2, err3
   integer,parameter :: BIG = IILARDOM*JJLARDOM*NNLANDUSE
   real, dimension(NNLANDUSE) :: tmp

   !JUN06 integer, dimension(NNLANDUSE) :: rivm2uk  ! maps RIVM landuse to
   !JUN06       ! uk. Need to set perm crops to arable, other to moorland, etc.
   integer :: uklu    ! UK (SEI) landuse code
   integer :: i_in, j_in ! for debug
    
   integer :: ncodes, kk ! debug
   logical :: debug_flag

   if ( DEBUG_DEP ) write(*,*) "UKDEP Starting ReadLandUse, me ",me
   if ( NLANDUSE /= NNLANDUSE )   WRITE(*,*) 'MPI_ABORT: ', "NNLANDuseerror"  ! Why need both? 
   if ( NLANDUSE /= NNLANDUSE ) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 

   maxlufound = 0   
   if ( me == 0 ) then


       allocate(g_ncodes(GIMAX,GJMAX),stat=err1)
       allocate(g_codes (GIMAX,GJMAX,NLUMAX),stat=err2)
       allocate(g_data  (GIMAX,GJMAX,NLUMAX),stat=err3)
       if ( err1 /= 0 .or. err2/=0 .or. err3/=0 )then
          WRITE(*,*) 'MPI_ABORT: ', "ioserror: landuse"
          call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
       endif
       g_ncodes(:,:)   = 0       !/**  initialise  **/
       g_data  (:,:,:) = 0.0     !/**  initialise  **/

      call open_file(IO_FORES,"r","landuse.JUN06",needed=.true.,skip=1)
      if ( ios /= 0 )   WRITE(*,*) 'MPI_ABORT: ', "ioserror: landuse"  
      if ( ios /= 0 ) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 


      do n = 1, BIG
         read(IO_FORES,*,iostat=ios) i,j, ( tmp(lu), lu=1,NNLANDUSE)
         if ( ios /= 0 ) exit   ! likely end of file
         if ( DEBUG_DEP ) debug_flag = ( i == DEBUG_i  .and. j == DEBUG_j  )

         i_in   = i ! for debug
         j_in   = j

         i = i - ISMBEG + 1
         j = j - JSMBEG + 1
         if ( i >= 1 .and. i <= GIMAX .and. &
              j >= 1 .and. j <= GJMAX  ) then
          if ( debug_flag ) then
            write(*,*) "MEDOAK", gb_glob(i,j), IAM_MEDOAK, tmp(IAM_MEDOAK) 
          endif
             if ( gb_glob(i,j) > 50.0 ) tmp(IAM_MEDOAK) = 0.0 !JUN06 
          if ( debug_flag ) then
            write(*,*) "MEDOAK", gb_glob(i,j), IAM_MEDOAK, tmp(IAM_MEDOAK) 
          endif
             do lu = 1, NNLANDUSE
                 if ( tmp(lu) > 0.0 ) then

                    uklu = lu

                    call GridAllocate("LANDUSE",i,j,uklu,NLUMAX, &
                      index_lu, maxlufound, g_codes, g_ncodes,errmsg)

                    if (errmsg /= "ok" )   WRITE(*,*) 'MPI_ABORT: ', errmsg 
                      if (errmsg /= "ok" ) call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
   
                    g_data(i,j,index_lu) = &
                       g_data(i,j,index_lu) + 0.01 * tmp(lu)
                 end if
             end do ! lu

         end if
      end do

      close(IO_FORES)
      write(6,*) "UK_ml: maxlufound = ", maxlufound
      if ( DEBUG_DEP ) write(*,*)  "DEBUG_DEP FILE CLOSED nrecords=", n
       
   end if  !  ( me == 0 )


   call global2local_int(g_ncodes,landuse_ncodes,MSG_READ1,   &
                                      GIMAX, GJMAX, 1,1,1)
   call global2local_int(g_codes,landuse_codes,MSG_READ2,   &
                                      GIMAX, GJMAX, NLUMAX, 1,1)
   call global2local(g_data,landuse_data,MSG_READ3,   &
                                    1,GIMAX, GJMAX, NLUMAX, 1,1)

   if ( me == 0 ) then
       deallocate(g_ncodes,stat=err1)
       deallocate(g_codes ,stat=err2)
       deallocate(g_data  ,stat=err3)

       if ( err1 /= 0 .or. err2 /= 0 .or. err3 /= 0 ) then
            WRITE(*,*) 'MPI_ABORT: ', "De-Allocerror - g_landuse" 
            call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
       end if ! errors
   end if ! me==0

  end subroutine  ReadLanduse
 
  !-------------------------------------------------------------------------
  subroutine  SetLandUse()
    integer :: i,j,ilu,lu, nlu, n ! indices
    logical, save :: my_first_call = .true.
    logical :: debug_flag = .true.
    real :: hveg!
    integer :: effectivdaynumber !6 monthes shift in Southern hemisphere.
!pw Treatment of growing seasons in the southern hemisphere:
!   all the static definitions (SGS,EGS...) refer to northern hemisphere, but the actual 
!   simulation dates are shifted by 6 monthes in the southern hemisphere by using
!   uses effectivdaynumber and mod(current_date%month+5,12)+1 in southern hemis


    if ( DEBUG_DEP .and. debug_proc ) write(*,*) "UKDEP SetLandUse, me, day ", me, daynumber

    if ( my_first_call ) then
        if ( DEBUG_DEP .and. debug_proc ) write(*,*) "UKDEP FIrst Start SetLandUse, me ", me

!pw effectiv daynumber to shift 6 month when in southern hemisphere
    effectivdaynumber=daynumber

       ! crops are distinguished since their height varies
       ! over the growing season. Note changed rule below:

        SAIadd = 0.0
        luflux_wanted = .false.
        do lu = 1, NLANDUSE
         crops(lu) = ( hveg_max(lu) < 5 .and. SGS50(lu) > 1 )
         bulk (lu) = ( LAImax(lu)   < 0.0 )   ! Set neg. in ukdep_biomass.dat
         water(lu) = ( hveg_max(lu) < 0.0 )   ! Set neg. in ukdep_biomass.dat
         forest(lu) = ( hveg_max(lu) > 4.99 .and. LAImax(lu) > 0.1 )
         conif_forest(lu) = ( forest(lu) .and. SGS50(lu) <=1 .and. &
              index(luname(lu),"broadleaf") == 0 )  !JUN06 fixed Bug?
         urban (lu) = ( hveg_max(lu) > 5.0 .and. LAImax(lu) < 0.0 )
         vegetation (lu) = ( hveg_max(lu) > 0.0 .and. .not.urban(lu) )

        !/ Set input negative values to physical ones, since we have
        !  now set bulk, water, etc.

         hveg_max(lu) = max( hveg_max(lu), 0.0)
         LAImax(lu)   = max( LAImax(lu),   0.0)
         if( forest(lu) ) SAIadd(lu) = 1.0   ! Addition to LAI to get surface area

         if( lu==IAM_MEDOAK)  conif_forest(lu) = .false.  !JUN06 
         if ( index(luname(lu),"IAM_") > 0  ) luflux_wanted(lu)=.true.   ! JUN06
         if ( DEBUG_DEP .and. debug_proc ) write(*,*) "UKDEP_VEG", lu, &
           crops(lu), bulk(lu), forest(lu), conif_forest(lu), luflux_wanted(lu)
        end do


     ! ICP - define whuich landuse we might be interested in stomtal
     ! fluxes for:   (inly used if STO_FLUXES set in My_DryDep)
      !where ( forest .and. .not. conif_forest ) !ds_sep27  .or. crops )
      !JUN06 luflux_wanted(WHEAT) = .true.

      !/ 2./ -- Calculate growing seasons where needed and water_fraction
      !          (for Rn emissions)

        water_fraction(:,:) = 0.0         !ds Pb210 
        ice_fraction(:,:)   = 0.0         !ds Pb210 

        do i = li0, li1
          do j = lj0, lj1

             do ilu= 1, landuse_ncodes(i,j)
                lu      = landuse_codes(i,j,ilu)

                if ( SGS50(lu) > 0 ) then  ! need to set growing seasons 

                    call get_growing_season( lu,abs(gb(i,j)),&  
                            landuse_SGS(i,j,ilu),landuse_EGS(i,j,ilu) )
                else
                   landuse_SGS(i,j,ilu) =  SGS50(lu)
                   landuse_EGS(i,j,ilu) =  EGS50(lu)
                end if

             !JUN06
             if ( index(luname(lu),"IAM_MEDOAK") > 0 ) then
                Astart(i,j,ilu) = 80
                Aend(i,j,ilu)   = 320
             else if ( forest(lu) .and. .not. conif_forest(lu) ) then
                Astart(i,j,ilu) = landuse_SGS(i,j,ilu)
                Aend(i,j,ilu)   = landuse_EGS(i,j,ilu) - 10.0   !! End of discolaration  
             else ! conif_forest(lu)
                Astart(i,j,ilu) = landuse_SGS(i,j,ilu)
                Aend(i,j,ilu)   = landuse_EGS(i,j,ilu)  
             end if

               !/ for landuse classes with bulk-resistances, we only
               !  need to specify height once. Dummy values are assigned
               !  to LAI and gpot:

                if ( bulk(lu) ) then
                    landuse_hveg(i,j,ilu) =  hveg_max(lu)
                    landuse_LAI(i,j,ilu)  =  0.0          
                    landuse_fphen(i,j,ilu) =  0.0          
                 end if

               !ds Pb210 : water fraction

                if ( lu == LU_WATER ) water_fraction(i,j) = landuse_data(i,j,ilu)
                if ( lu == LU_ICE   )   ice_fraction(i,j) = landuse_data(i,j,ilu)

                if ( DEBUG_DEP .and. debug_proc .and. &
                    i == debug_li .and. j == debug_lj ) then
                      write(*,*) "DEBUG WATER ", ilu, lu, &
                          water_fraction(i,j), ice_fraction(i,j)
                end if

            end do ! ilu
          end do ! j
        end do ! i

        my_first_call = .false.
    end if ! my_first_call
    
     do i = li0, li1
       do j = lj0, lj1

          !pw effectiv daynumber to shift 6 month when in southern hemisphere
          if(gb(i,j)<0.0)effectivdaynumber=mod(daynumber+182,nydays)+1 

          debug_flag = ( debug_proc .and. i == debug_li .and. j == debug_lj ) 
          if ( DEBUG_DEP .and. debug_flag ) then
                 write(*,"(a12,i3,i4)") "LANDUSE N Day? ", landuse_ncodes(i,j), daynumber
          end if
          do ilu= 1, landuse_ncodes(i,j)
             lu      = landuse_codes(i,j,ilu)
             if ( lu <= 0 .or. lu > NLANDUSE ) then
                WRITE(*,*) 'MPI_ABORT: ', "SetLandUselu<0" 
                call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
             endif

             if ( bulk(lu) ) cycle    !else Growing veg present:

             landuse_LAI(i,j,ilu) = Polygon(effectivdaynumber, &
                                      0.0, LAImin(lu), LAImax(lu),&
                                      landuse_SGS(i,j,ilu), SLAIlen(lu), &
                                      landuse_EGS(i,j,ilu), ELAIlen(lu))

             landuse_fphen(i,j,ilu) = fPhenology( debug_flag, effectivdaynumber, &
                          f_phen_a(lu), f_phen_b(lu), f_phen_c(lu),&
                            f_phen_d(lu), f_phen_Slen(lu), f_phen_Elen(lu), &
                              landuse_SGS(i,j,ilu), landuse_EGS(i,j,ilu) , Astart(i,j,ilu), Aend(i,j,ilu))


            !if ( DEBUG_DEP .and. debug_flag ) then
            !       write(*,"(a12,i3,i4,f7.2,2f8.3,4i4)") "LANDPhen0 ", lu,  &
            !         daynumber, -99.9,  &
            !          landuse_LAI(i,j,ilu), landuse_fphen(i,j,ilu), & 
            !          landuse_SGS(i,j,ilu), landuse_EGS(i,j,ilu), Astart(i,j,ilu), Aend(i,j,ilu)
            ! end if

          ! For coniferous forest we need to correct for old needles.

             if ( conif_forest(lu) )then
                if(gb(i,j)<0.0)then
                   !southern hemisphere
                   call Conif_fphen( mod(current_date%month+5,12)+1 ,landuse_fphen(i,j,ilu))
                else
                   call Conif_fphen(current_date%month,landuse_fphen(i,j,ilu))
                endif
             endif


             hveg = hveg_max(lu)   ! default

      !ds SAIadd code moved from here to DryDep_ml to fix bug spotted by Peter,
      !   19/8/2003
      !SAIadd for other vegetation still defined in UK_ml

             if (  crops(lu) ) then

                if ( lu == IAM_WHEAT ) then ! for NEWAOT
                    if  ( effectivdaynumber >= landuse_SGS(i,j,ilu) .and. &
                          effectivdaynumber <= landuse_EGS(i,j,ilu)  ) then
                            InGrowingSeason(i,j) =  1
                    else
                            InGrowingSeason(i,j) =  0
                    end if
                end if

                if ( effectivdaynumber < landuse_SGS(i,j,ilu) .or. &
                     effectivdaynumber > landuse_EGS(i,j,ilu)  ) then
                   hveg = STUBBLE
                else if ( effectivdaynumber < &
                     (landuse_SGS(i,j,ilu) + SLAIlen(lu)) ) then
                   hveg=  hveg_max(lu) * landuse_LAI(i,j,ilu)/LAImax(lu)
                else if ( effectivdaynumber < landuse_EGS(i,j,lu) ) then
                   hveg = hveg_max(lu)                  ! not needed?
                end if
             end if ! crops

             landuse_hveg(i,j,ilu) =  hveg
            if ( DEBUG_DEP .and. debug_flag ) then
                   if(lu==IAM_WHEAT) write(*,*) "GROWSEASON ", effectivdaynumber, InGrowingSeason(i,j)
                   write(*,"(a12,i3,i4,f7.2,2f8.3,4i4)") "LANDPhen ", lu, daynumber, &
                     hveg, landuse_LAI(i,j,ilu), landuse_fphen(i,j,ilu), &
                     landuse_SGS(i,j,ilu), landuse_EGS(i,j,ilu), &
                     Astart(i,j,ilu), Aend(i,j,ilu)
            end if
                   

         end do ! lu
       end do ! j
    end do ! i
    if ( DEBUG_DEP .and. me==0 ) write(*,*)"UKDEP Finishing SetLandUse "
    if(debug_proc ) write(*,*) "LAST GROWSEASON ", effectivdaynumber, InGrowingSeason(debug_li,debug_lj)

  end subroutine  SetLandUse
  !-------------------------------------------------------------------------
! =====================================================================
    subroutine Conif_fphen(imm,f_phen)
! =====================================================================
!   modifies g_pot (g_age) for effect of older needles, with the simple
!   assumption that g_age(old) = 0.5.
!
   !/ arguments

    integer, intent(in) :: imm    ! month
    real,   intent(inout) :: f_phen  ! Requires initial input of f_phen 
                                     ! (once obtained as output from g_stomatal)

   !/ Some parameters:
   !  Proportion of needles which are from current year:
    real, parameter, dimension(12) :: Pc = (/  &
                   0.53, 0.535, 0.54, 0.545, 0.1, 0.15,  &
                   0.27,  0.36, 0.42,  0.48, 0.5,  0.5  /)

    real, parameter :: F_OLD = 0.5  ! value of f_phen for old needles



!needles from current year assumed to have g_pot as evaluated above;
!needles from previous years assumed to have f_phen of 0.5
!The sum of the f_phen's for the current year is added to the sum of the
!f_phen's for previous years to obtain the overall f_phen for the landuse 
!category temp. conif. forests.

    f_phen = Pc(imm)*f_phen + (1.0-Pc(imm))*F_OLD

  end subroutine Conif_fphen
! =====================================================================

                        
end module UKdep_ml


