module UKdep_ml

use Dates_ml,       only: daynumber, nydays
use DepVariables_ml,only: NLANDUSE             &  ! No. UK land-classes
                      ,luname              &
                      ,crops, bulk, water  & ! logical variables
                      ,forest,conif_forest & !    "      "
                      ,vegetation, urban   & !
                      ,luflux_wanted       & ! 
                      ,STUBBLE             & ! Ht. of stubble (m)
                      ,SAIadd              & ! surface-area index, rv1.4.6
                      ,hveg_max, b_inc, albedo, NH4_pl, SGS50, DSGS   &
                      ,EGS50, DEGS, LAImin, LAImax, SLAIlen, ELAIlen  &
                      ,g_pot_min , Sg_potlen , Eg_potlen     &
                      ,g_max     , g_min     , g_lightfac    &
                      ,g_temp_min, g_temp_opt, g_temp_max  &
                      ,RgsS      , RgsO        &
                      ,VPD_max   , VPD_min     &
                      ,SWP_max   , PWP       , rootdepth 
use Functions_ml,   only: GridAllocate, Polygon
use GridValues_ml,  only: &!u7.4vg GlobalPosition, &
                           gb_glob, gb, i_glob, j_glob  ! latitude
use Io_ml,          only: open_file, ios, IO_FORES
use ModelConstants_ml,  only : current_date, debug_i, debug_j, NNLANDUSE
use UKsetup_ml,     only: ukdep_init, get_growing_season
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
  private :: Conif_Gpot

 integer, public, parameter :: NLUMAX = 15 ! max no. landuse per grid

 ! Keep crops and bulk variables here until I figure our
 ! how to use logicals in gc_send

!rv1.2 logical, public, dimension(NLANDUSE), save :: &
!rv1.2            crops   &! true for veg which grows...
!rv1.2           ,bulk     ! true for land-classes without LAI

!rv1.4.6 real, public, dimension(NLANDUSE), save :: &
!rv1.4.6       SAIadd                  ! Additional surface area for bark, twigs

!.. very crude for now:
!rv1.2 integer, public, parameter :: NDEP_FLUX = 5,  &  ! 
!rv1.2     DEP_VG_REF = 1, DEP_VG_1M = 2, DEP_VG_STO = 3, &
!rv1.2     DEP_FL_REF = 4, DEP_FL_STO = 5
!rv1.2 real, public, save, dimension(NDEP_FLUX,MAXLIMAX,MAXLJMAX,2) :: dep_flux

 integer,public,save,dimension(MAXLIMAX,MAXLJMAX)        :: landuse_ncodes 
 integer,public,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: landuse_codes  
 real,   public,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: landuse_data 

 integer,public,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: &
          landuse_SGS   &    ! Start of growing season (days)
         ,landuse_EGS        ! End of growing season (days)

 real,   public,save,dimension(MAXLIMAX,MAXLJMAX,NLUMAX) :: &
          landuse_LAI   &    ! Leaf-area-index (m2/m2)
         ,landuse_hveg  &    ! Max. height of veg.
         ,landuse_gpot        ! Potential (age) factor for Jarvis-calc

 logical, private, parameter :: DEBUG_DEP = .true.
 character(len=30), private :: errmsg


contains

 !--------------------------------------------------------------------------
  subroutine Init_ukdep()
    character(len=20) :: errmsg
    integer :: gc_info   ! for error messages from gc
    integer :: i,j,lu    ! indices


  !/ 1./ -- Read in basic data for UK model

    if ( me == 0 ) then 
       !=====================================
        call ukdep_init(errmsg)
       !=====================================
        if ( errmsg /= "ok" ) then
           errmsg = "ukdep_init: " // errmsg
           call gc_abort(me,NPROC,errmsg)
        end if
    end if

    call gc_rbcast(101,NLANDUSE,0,NPROC,gc_info,hveg_max)
    call gc_rbcast(102,NLANDUSE,0,NPROC,gc_info,b_inc)
    call gc_rbcast(103,NLANDUSE,0,NPROC,gc_info,albedo)
    call gc_rbcast(104,NLANDUSE,0,NPROC,gc_info,NH4_pl)
    call gc_rbcast(105,NLANDUSE,0,NPROC,gc_info,SGS50)
    call gc_rbcast(106,NLANDUSE,0,NPROC,gc_info,DSGS)
    call gc_rbcast(107,NLANDUSE,0,NPROC,gc_info,EGS50)
    call gc_rbcast(108,NLANDUSE,0,NPROC,gc_info,DEGS)
    call gc_rbcast(109,NLANDUSE,0,NPROC,gc_info,LAImin)
    call gc_rbcast(110,NLANDUSE,0,NPROC,gc_info,LAImax)
    call gc_rbcast(111,NLANDUSE,0,NPROC,gc_info,SLAIlen)
    call gc_rbcast(112,NLANDUSE,0,NPROC,gc_info,ELAIlen)

                      
    call gc_rbcast(120,NLANDUSE,0,NPROC,gc_info,g_pot_min)
    call gc_rbcast(121,NLANDUSE,0,NPROC,gc_info, Sg_potlen) 
    call gc_rbcast(122,NLANDUSE,0,NPROC,gc_info, Eg_potlen)
    call gc_rbcast(123,NLANDUSE,0,NPROC,gc_info,g_max)    
    call gc_rbcast(124,NLANDUSE,0,NPROC,gc_info, g_min)    
    call gc_rbcast(125,NLANDUSE,0,NPROC,gc_info, g_lightfac)
    call gc_rbcast(126,NLANDUSE,0,NPROC,gc_info,g_temp_min)
    call gc_rbcast(127,NLANDUSE,0,NPROC,gc_info, g_temp_opt)
    call gc_rbcast(128,NLANDUSE,0,NPROC,gc_info, g_temp_max)
                      
    call gc_rbcast(131,NLANDUSE,0,NPROC,gc_info,RgsS)
    call gc_rbcast(132,NLANDUSE,0,NPROC,gc_info, RgsO)
    call gc_rbcast(133,NLANDUSE,0,NPROC,gc_info,VPD_max)
    call gc_rbcast(134,NLANDUSE,0,NPROC,gc_info, VPD_min)
    call gc_rbcast(135,NLANDUSE,0,NPROC,gc_info,SWP_max)   
    call gc_rbcast(136,NLANDUSE,0,NPROC,gc_info, PWP)       
    call gc_rbcast(137,NLANDUSE,0,NPROC,gc_info, rootdepth) 

end subroutine Init_ukdep
 !--------------------------------------------------------------------------
subroutine ReadLanduse()

! arrays for whole EMEP area ( allocated for me=0 only)

  integer, allocatable, dimension(:,:)      :: g_ncodes 
  integer, allocatable, dimension(:,:,:)    :: g_codes  
  real,    allocatable, dimension(:,:,:)    :: g_data 

   integer :: i,j,n,lu, index_lu, maxlufound
   integer :: err1, err2, err3
!TFMM   integer,parameter :: NNLANDUSE = 17  !temp.
   integer,parameter :: BIG = IILARDOM*JJLARDOM*NNLANDUSE
   real, dimension(NNLANDUSE) :: tmp

   integer, dimension(NNLANDUSE) :: rivm2uk  ! maps RIVM landuse to
         ! uk. Need to set perm crops to arable, other to moorland, etc.
   integer :: uklu    ! UK (SEI) landuse code
   integer :: i_in, j_in ! for debug
    
   integer :: ncodes, kk ! debug
   logical :: debug_flag

   if ( DEBUG_DEP ) print *, "UKDEP Starting ReadLandUse, me ", me

   maxlufound = 0   
   if ( me == 0 ) then

       allocate(g_ncodes(GIMAX,GJMAX),stat=err1)
       allocate(g_codes (GIMAX,GJMAX,NLUMAX),stat=err2)
       allocate(g_data  (GIMAX,GJMAX,NLUMAX),stat=err3)
       if ( err1 /= 0 .or. err2/=0 .or. err3/=0 ) &
            call gc_abort(me,NPROC,"ios error: landuse") 

       g_ncodes(:,:)   = 0       !/**  initialise  **/
       g_data  (:,:,:) = 0.0     !/**  initialise  **/

      !TFMM call open_file(IO_FORES,"r","landuse.170",needed=.true.)
      call open_file(IO_FORES,"r","landuse.tf2",needed=.true.,skip=1)
      if ( ios /= 0 ) call gc_abort(me,NPROC,"ios error: landuse") 

!c   rivm2uk
!c       based on rivm/lbg land use data, 8 classes
!c       plus extra classes for special use:
!c
!c       1 = grass
!c       2 = arable
!c       3 = permanent crops
!c       4 = coniferous forest
!c       5 = deciduous forest
!c       6 = water
!c       7 = urban
!c       8 = other i.e. short grassy area
!c       9 = desert      !used for n-africa region emep-grid
!c       10= ice
!c

      do n = 1, BIG
         read(IO_FORES,*,iostat=ios) i,j, ( tmp(lu), lu=1,NNLANDUSE)
         if ( ios /= 0 ) exit   ! likely end of file
         if ( DEBUG_DEP ) debug_flag = ( i == debug_i .and. j == debug_j )

       !TFIAM changes: treat permanent crops as grass. Assign all
       ! arable crops to TC - wheat etc.
       !TFMM  if ( gb_glob(i,j) > 60.0 ) then   ! use Tundra for 8, 9
       !TFMM       !TFIAM rivm2uk = (/ 10, 5, 5, 1, 2, 15, 17,13, 13, 16 /)
       !TFMM       rivm2uk = (/ 10, 5, 10, 1, 2, 15, 17,13, 13, 16 /)
       !TFMM  else if ( gb_glob(i,j) > 45.0 ) then   !use moorland for 8
       !TFMM       !TFIAM rivm2uk = (/ 10, 5, 5, 1, 2, 15, 17, 8, 14, 16 /)
       !TFMM       rivm2uk = (/ 10, 5, 10, 1, 2, 15, 17, 8, 14, 16 /)
       !TFMM  else  ! Mediterranean 
       !TFMM       !TFIAM rivm2uk = (/ 10, 6, 6, 3, 4, 15, 17,11, 14, 16 /)
       !TFMM       !TFMM rivm2uk = (/ 10, 5, 10, 3, 4, 15, 17,11, 14, 16 /)
       !TFMM       !TFMM - use temp forests overall
       !TFMM       rivm2uk = (/ 10, 5, 10, 1, 2, 15, 17,11, 14, 16 /)
       !TFMM  end if
         i_in   = i ! for debug
         j_in   = j

         if ( DEBUG_DEP .and. debug_flag ) then
             print *, "UKDEP-ukluS", gb_glob(i,j), rivm2uk(2)
         endif
      
         i = i - ISMBEG + 1
         j = j - JSMBEG + 1
         if ( i >= 1 .and. i <= GIMAX .and. &
              j >= 1 .and. j <= GJMAX  ) then
             do lu = 1, NNLANDUSE
                 if ( tmp(lu) > 0.0 ) then

                    !TFMMuklu = rivm2uk(lu)  ! uk code
                    uklu = lu

                    call GridAllocate("LANDUSE",i,j,uklu,NLUMAX, &
                      index_lu, maxlufound, g_codes, g_ncodes,errmsg)

                    if ( DEBUG_DEP .and. debug_flag ) then
                       print "(a12,2i3,2i4,4i3,f12.2)", "UKDEP-Grid", lu, uklu, i_in,j_in,i,j,&
                              index_lu,g_ncodes(i,j), tmp(lu)
                       ncodes = g_ncodes(i,j)
                       do kk = 1, ncodes
                         print *, "g_codes ", kk, g_codes(i,j,kk)
                       end do
                    end if ! DEBUG

                    if (errmsg /= "ok" ) call gc_abort(me,NPROC,errmsg)
   
                    g_data(i,j,index_lu) = &
                       g_data(i,j,index_lu) + 0.01 * tmp(lu)
                 end if
             end do ! lu
             !============ WHEAT FIX ===============================!
             !TFMM if( tmp(2) > 0.001 ) then
             !TFMM    uklu = 9  ! Wheat in dummy files - add 1% of arable here
             !TFMM       call GridAllocate("LANDUSE",i,j,uklu,NLUMAX, &
             !TFMM         index_lu, maxlufound, g_codes, g_ncodes,errmsg)
             !TFMM       if (errmsg /= "ok" ) call gc_abort(me,NPROC,errmsg)
             !TFMM       g_data(i,j,index_lu) = &
             !TFMM          g_data(i,j,index_lu) + 0.0001 * tmp(lu)
             !TFMM end if
             !============ WHEAT FIX ===============================!

         end if
      end do

      close(IO_FORES)
      write(6,*) "UK_ml: maxlufound = ", maxlufound
      if ( DEBUG_DEP ) print *,  "DEBUG_DEP FILE CLOSED nrecords=", n
       
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
          call gc_abort(me,NPROC,"De-Alloc error - g_landuse")
       end if ! errors
   end if ! me==0

   !rv1.2 if ( DEBUG_DEP .and. local_i > 0 ) then
   !rv1.2     print "(a20,i4)", "UKDEP FinReadLU, me ", me
   !rv1.2     print "(a20,i4,i6,i3,f8.3)", "UKDEP FinReadLU, me ", me, &
   !rv1.2       landuse_ncodes(local_i,local_j), landuse_codes(local_i,local_j,1), &
   !rv1.2          landuse_data(local_i,local_j,1)
   !rv1.2 end if ! DEBUG

  end subroutine  ReadLanduse
 
  !-------------------------------------------------------------------------
  subroutine  SetLandUse()
    integer :: i,j,ilu,lu, nlu, n ! indices
    logical, save :: my_first_call = .true.
    real :: hveg

    if ( DEBUG_DEP ) print *, "UKDEP SetLandUse, me, day ", me, daynumber

    if ( my_first_call ) then
        if ( DEBUG_DEP ) print *, "UKDEP FIrst Start SetLandUse, me ", me


       !u7.4vf crops are distinguished since their height varies
       !       over the growing season. Note changed rule below:

        do lu = 1, NLANDUSE
         crops(lu) = ( hveg_max(lu) < 5 .and. SGS50(lu) > 1 )
         bulk (lu) = ( LAImax(lu)   < 0.0 )   ! Set neg. in ukdep_biomass.dat
         water(lu) = ( hveg_max(lu) < 0.0 )   ! Set neg. in ukdep_biomass.dat
         forest(lu) = ( hveg_max(lu) > 4.99 .and. LAImax(lu) > 0.1 )
         conif_forest(lu) = ( forest(lu) .and. SGS50(lu) <=1 )  !Bug?
                                                    ! Includes Med. broadleaf!
         urban (lu) = ( hveg_max(lu) > 5.0 .and. LAImax(lu) < 0.0 )
         vegetation (lu) = ( hveg_max(lu) > 0.0 .and. .not.urban(lu) )

        !/ Set input negative values to physical ones, since we have
        !  now set bulk, water, etc.

         hveg_max(lu) = max( hveg_max(lu), 0.0)
         LAImax(lu)   = max( LAImax(lu),   0.0)

!Tried 4l3, but ....
           if ( DEBUG_DEP .and. me == 0 ) print *, "UKDEP_VEG", lu, &
                           crops(lu), bulk(lu), forest(lu), conif_forest(lu)
        end do

      !/ 2./ -- Calculate additional surface area for trees

      SAIadd = 0.0
      where ( forest )
        SAIadd = 1.0     ! Addition to LAI to get surface area
      end where

     ! ICP - define whuich landuse we might be interested in stomtal
     ! fluxes for:   (inly used if STO_FLUXES set in My_DryDep)

      luflux_wanted = .false.
      where ( forest .or. crops )
        luflux_wanted = .true.
      end where
    ! CRUDE AND TMP !!!!!
      luflux_wanted(10) = .true.

      !/ 3./ -- Calculate growing seasons where needed

        do i = li0, li1
          do j = lj0, lj1

             do ilu= 1, landuse_ncodes(i,j)
                lu      = landuse_codes(i,j,ilu)

                if ( SGS50(lu) > 0 ) then  ! need to set growing seasons 

                    call get_growing_season( lu,gb(i,j),&
                            landuse_SGS(i,j,ilu),landuse_EGS(i,j,ilu) )
                else
                   landuse_SGS(i,j,ilu) =  SGS50(lu)
                   landuse_EGS(i,j,ilu) =  EGS50(lu)
                end if

               !/ for landuse classes with bulk-resistances, we only
               !  need to specify height once. Dummy values are assigned
               !  to LAI and gpot:

                if ( bulk(lu) ) then
                    landuse_hveg(i,j,ilu) =  hveg_max(lu)
                    landuse_LAI(i,j,ilu)  =  0.0          
                    landuse_gpot(i,j,ilu) =  0.0          
                 end if

            end do ! ilu
          end do ! j
        end do ! i

        my_first_call = .false.
    end if ! my_first_call
    
     do i = li0, li1
       do j = lj0, lj1

          do ilu= 1, landuse_ncodes(i,j)
             lu      = landuse_codes(i,j,ilu)
             if ( lu <= 0 ) call gc_abort(me,NPROC,"SetLandUse lu<0")

             !rv1.2 if ( bulk(lu) ) then
             !rv1.2    landuse_LAI(i,j,ilu) =  -99.99
             !rv1.2    landuse_hveg(i,j,ilu) =  hveg_max(lu)
             !rv1.2    hveg =  hveg_max(lu)

             if ( bulk(lu) ) cycle 
                !rv1.2 else ! Growing veg present

             landuse_LAI(i,j,ilu) = Polygon(daynumber, &
                                      0.0, LAImin(lu), LAImax(lu),&
                                      landuse_SGS(i,j,ilu), SLAIlen(lu), &
                                      landuse_EGS(i,j,ilu), ELAIlen(lu))

             landuse_gpot(i,j,ilu) =  Polygon(daynumber,&
                                       0.0,g_pot_min(lu),1.0,&
                                       landuse_SGS(i,j,ilu), Sg_potlen(lu), &
                                       landuse_EGS(i,j,ilu), Eg_potlen(lu))



          ! For coniferous forest we need to correct for old needles.

             if ( conif_forest(lu) )  &
                   call Conif_gpot(current_date%month,landuse_gpot(i,j,ilu))



             hveg = hveg_max(lu)   ! default

             if (  crops(lu) ) then

                if ( daynumber < landuse_SGS(i,j,ilu) .or. &
                     daynumber > landuse_EGS(i,j,ilu)  ) then
                     hveg = STUBBLE
                else if ( daynumber < &
                          (landuse_SGS(i,j,ilu) + SLAIlen(lu)) ) then
                     hveg=  hveg_max(lu) * landuse_LAI(i,j,ilu)/LAImax(lu)
                  !ICP - crude....
                     SAIadd(lu) = ( 5.0/3.5 - 1.0) * landuse_LAI(i,j,ilu)
                else if ( daynumber < landuse_EGS(i,j,lu) ) then
                     hveg = hveg_max(lu)                  ! not needed?
                     SAIadd(lu) = 1.5   ! Sensescent
                end if
             end if ! crops
             !rv1.2 end if ! bulk

             landuse_hveg(i,j,ilu) =  hveg

         end do ! lu
       end do ! j
    end do ! i
    if ( DEBUG_DEP ) print *, "UKDEP Finishing SetLandUse, me ", me

  end subroutine  SetLandUse
  !-------------------------------------------------------------------------
! =====================================================================
    subroutine Conif_gpot(imm,g_pot)
! =====================================================================
!   modifies g_pot (g_age) for effect of older needles, with the simple
!   assumption that g_age(old) = 0.5.
!
   !/ arguments

    integer, intent(in) :: imm    ! month
    real,   intent(inout) :: g_pot   ! Requires initial input of g_pot 
                                     ! (once obtained as output from g_stomatal)

   !/ Some parameters:
   !  Proportion of needles which are from current year:
    real, parameter, dimension(12) :: Pc = (/  &
                   0.53, 0.535, 0.54, 0.545, 0.1, 0.15,  &
                   0.27,  0.36, 0.42,  0.48, 0.5,  0.5  /)

    real, parameter :: G_POTOLD = 0.5  ! value of g_pot for old needles



!needles from current year assumed to have g_pot as evaluated above;
!needles from previous years assumed to have g_pot of 0.5
!The sum of the g_pot's for the current year is added to the sum of the
!g_pot's for previous years to obtain the overall g_pot for the landuse 
!category temp. conif. forests.

    g_pot = Pc(imm)*g_pot + (1.0-Pc(imm))*G_POTOLD

  end subroutine Conif_gpot
! =====================================================================

                        
end module UKdep_ml


