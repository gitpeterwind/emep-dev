
!==============================================================================
module Derived_ml

  !---------------------------------------------------------------------------
  ! DESCRIPTION
  ! This module performs the calculations associated with "derived" 2D and 3D,
  ! such as accumulated precipitation or sulphate, daily, monthly or yearly 
  ! averages, depositions. These fields are all typically output as binary 
  ! fields.
  !
  ! The definitions of the derived fields shoulöd have been specified in the
  ! user-defined My_Derived_ml.
  !
  ! These modules combine Steffen's rewriting of the output routines
  ! (IOU notation), elements of chemint_mach, Bud_ml, etc.
  ! Re-coded to use only 2 types of data (d_2 and d_3)
  ! and F90 types for Deriv by ds, Sept. 2001. 

  ! User-defined routines and treatments are often needed here. So far I have 
  ! added stuff for VOC, AOTs, accsu, taken from the MACHO chemint code. In 
  ! general
  ! such code should be added in such a way that it isn't activated if not
  ! needed. It then doesn't need to be commented out if not used.
  ! For "accsu" I haven't found a way to do this though.
  ! See the examples given.
  ! ds, 28/9/01
  !---------------------------------------------------------------------------

use My_Derived_ml  ! Definitions of derived fields, NWDEP, etc., f_wdep, etc.
use Chemfields_ml, only : xn_adv, xn_shl, cfac,xn_bgn
use GenSpec_adv_ml         ! Use NSPEC_ADV amd any of IXADV_ indices
use GenSpec_shl_ml,  only : NSPEC_SHL  ! Use NSPEC_SHL to map adv to tot.indices
!u1 use GenSpec_maps_ml, only : MAP_ADV2TOT
use GenChemicals_ml, only : species
use ModelConstants_ml, &
                   only: KMAX_MID &   ! =>  z dimension
                         , MFAC   &   ! for conversion of units
                         , PPBINV &   ! for conversion of units for DEBUG
                         , current_date
use Par_ml,    only: MAXLIMAX,MAXLJMAX, &   ! => max. x, y dimensions
                     me, NPROC,         &   ! for gc_abort checks
                     limax, ljmax           ! => used x, y area 
use PhysicalConstants_ml,  only : PI
!hf hmix
use Met_ml, only :   roa,pzpbl,xksig

implicit none
private

 public  :: SumDerived           ! old chemint
 public  :: ResetDerived         ! Resets values to zero
 public  :: DerivedProds         ! Calculates any production terms
 private :: Setups 
 private :: Consistency_checks   !ds  - checks index numbers from My_Derived
 private :: Consistency_count    !ds  - checks index numbers from My_Derived
 private :: Setup_VOC            ! Defines VOC group
!jej private :: Derived              ! Calculations of sums, avgs etc.
 public :: Derived              ! Calculations of sums, avgs etc.
!6c private :: acc_sulphate         ! Sums sulphate column
!6c private :: aot_calc             ! Calculates daylight AOTs
 private :: voc_2dcalc           ! Calculates sum of VOC for 2d fields
 private :: voc_3dcalc           ! Calculates sum of VOC for 3d fields



  ! Counters to keep track of averaging  (old nav, nav8to16, nav_8to16, etc
  ! Initialise here to zero.

    integer, public, dimension(NWDEP,LENOUT2D),     save :: nav_wdep = 0
    integer, public, dimension(NDDEP,LENOUT2D),     save :: nav_ddep = 0
    integer, public, dimension(NDERIV_2D,LENOUT2D), save :: nav_2d   = 0
    integer, public, dimension(NDERIV_3D,LENOUT3D), save :: nav_3d   = 0

   ! Note - previous versions did not have the LENOUT2D dimension
   ! for wet and dry deposition. Why not?  Are annual or daily
   ! depositions never printed? Since I prefer to keep all 2d
   ! fields as similar as posisble, I have kept this dimension
   ! for now - ds


!!    real, private :: timefrac     ! timestep (dt) as fraction of hour 
                                    ! (3600/dt) -- Decleared in individual 
                                    ! subroutines

   !-- some variables for the VOC sum done for ozone models
   !   (have no effect in non-ozone models - leave in code)

   integer, private, save :: nvoc   ! No. VOCs 
   integer, private, dimension(NSPEC_ADV), save :: &
             voc_index, &     ! Index of VOC in xn_adv
             voc_carbon       ! Number of C atoms

   logical, private, parameter :: MY_DEBUG = .false.
   integer, private :: i,j,k,n, ivoc, index    ! Local loop variables

   contains

    !=========================================================================
    subroutine SumDerived(dt)
      logical, save :: my_first_call = .true.
      real, intent(in) :: dt  !  time-step used in intergrations

      if ( my_first_call ) then
          print *, "INITIALISE My DERIVED STUFF"
          call Set_My_Derived()
          call Consistency_checks()  !ds - checks index numbers from My_Derived
          call Setups()
          my_first_call = .false. 
!jej      else
!jej          call Derived(dt)
      end if

    end subroutine SumDerived

    !=========================================================================
     subroutine Consistency_checks()  

    !/** ds - checks index numbers from My_Derived to look for duplicates

       integer, parameter :: MAX_INDEX = 2000
       integer, dimension(MAX_INDEX) :: index_used
       integer :: i

       index_used = 0
       call  Consistency_count(MAX_INDEX, index_used, NWDEP, f_wdep) 
       call  Consistency_count(MAX_INDEX, index_used, NDDEP, f_ddep)
       call  Consistency_count(MAX_INDEX, index_used, NDERIV_2D, f_2d) 
       call  Consistency_count(MAX_INDEX, index_used, NDERIV_3D, f_3d) 

       if ( any( index_used > 1 ) ) then
           do i = 1, MAX_INDEX
             if( index_used(i) > 1 ) print *,  &
                   "Derived code problem, index: ",i, index_used(i)
           end do
           call gc_abort(me,NPROC,"Derived code problem!!")
       end if
     end subroutine

    !=========================================================================
     subroutine Consistency_count(max,index_used, n,data)

    !/** ds adds up number of times each code from derived data array
    !    is used.

      integer, intent(in) :: max, n
      integer, dimension(:), intent(inout)  :: index_used
      type(Deriv), dimension(n), intent(in) :: data

      integer :: code, i

        do i = 1, n
           code = data(i)%code
           if ( code > max ) call gc_abort(me,NPROC,"My_Derived code >max!")
           index_used( code )  = index_used( code )  + 1
        end do
     end subroutine Consistency_count
    
    !=========================================================================
     subroutine Setups()

    !/** flexibility note. By making use of character-based tests such
    !    as for "VOC" below, we achieve code which can stay for both MADE and
    !    MACHO without having to define non-used indices. 
    !    Similarly, we avoid the previous "if NUM_ACCSU eq 1" type test,
    !    since the appropriate code will now only activate 

    !/ ** if voc wanted, set up voc_array. Works for all ozone chemistries
    !     (and code not called for MADE-type).

      if ( any(  f_2d(:)%class == "VOC" ) .or. &
           any(  f_3d(:)%class == "VOC" )  ) then
            call Setup_VOC()
            print *, "Derived VOC setup returns ", nvoc, "vocs"
            print "(a12,30i3)",  "indices ", voc_index(1:nvoc)
            print "(a12,30i3)",  "carbons ", voc_carbon(1:nvoc)
      end if


    end subroutine Setups
    !=========================================================================

    subroutine Derived(dt)

    !/** DESCRIPTION
    !  Integration and averaging of chemical fields. Intended to be
    !  a more flexible version of the old chemint routine.
    !  Includes AOT40, AOT60 if present

      real, intent(in) :: dt  !  time-step used in intergrations

      !u1 character(len=size(f_2d%class)) :: typ  !  See defs of f_2d
      character(len=len(f_2d%class)) :: typ  !  See defs of f_2d
      real :: thour                          ! Time of day (GMT)
      real :: timefrac                       ! dt as fraction of hour (3600/dt)
      real, dimension(MAXLIMAX,MAXLJMAX) :: density !  roa (kgair m-3 when 
                                                    ! scale in ug,  else 1

      !bug: timefrac = 3600.0/dt
      timefrac = dt/3600.0


     !/***** 2-D fields **************************

     do n = 1, NDERIV_2D

        
        typ = f_2d(n)%class

        !u4 density(:,:) = 1.  ! u4 - only do if needed:

        if ( f_2d(n)%rho ) then
            forall ( i=1:limax, j=1:ljmax )
                density(i,j) = roa(i,j,KMAX_MID,1)
            end forall
        else !u4
            density(:,:) = 1.0
        end if

        !/** user-defined time-averaging. Here we have defined TADV and TVOC
        !    so that 8-hour daytime averages will be calculated. 
        !    Just comment out if not wanted, or (better!) don't define any
        !    f_2d as TADV or TVOC

        if ( typ == "TADV" .or. typ == "TVOC" ) then
             thour = current_date%hour+current_date%seconds/3600.0
             if(thour <= 8.0 .or. thour > 16.0 ) cycle  ! Start next species
        end if
!hf hmix average at 00 and 12

        if ( typ == "HMIX00" .or. typ == "XKSIG00" ) then
             thour = current_date%hour+current_date%seconds/3600.0
             if(thour /= 0.0 ) cycle  ! Start next species
        end if

        if ( typ == "HMIX12" .or. typ == "XKSIG12" ) then
             thour = current_date%hour+current_date%seconds/3600.0
             if(thour /= 12.0 ) cycle  ! Start next species
        end if

        index = f_2d(n)%index
        select case ( typ )

!hf hmix
          case ( "HMIX", "HMIX00", "HMIX12" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = pzpbl(i,j)  
            end forall

  

            if ( MY_DEBUG .and. me == 0 ) then
             print "(a12,2i4,4f12.3)", "HMIX" , n &
                      , d_2d(n,3,3,IOU_INST)       
            end if

         ! Simple advected species:
          case ( "ADV", "TADV" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_INST) = xn_adv(index,i,j,KMAX_MID)  &
                                     * cfac(index,i,j) * density(i,j)  
            end forall

            if ( MY_DEBUG .and. me == 0 ) then
             print "(a12,2i4,4f12.3)", "JUST ADV" , n, index  &
                      , d_2d(n,3,3,IOU_INST) * PPBINV      &
                      ,  xn_adv(index,3,3,KMAX_MID)* PPBINV  &
                      ,  density(3,3), cfac(index,3,3)
            end if

         !u4 Daily maxima - advected
         !u4 ds - derived from jej's MAX_O3 stuff, but added density
         !u4       Is this right ...??!!!

          case ( "MAXADV" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_DAY) = max( d_2d( n, i,j,IOU_DAY), &
                                xn_adv(index,i,j,KMAX_MID)  &
                                     * cfac(index,i,j) * density(i,j) )
            end forall
            if ( MY_DEBUG .and. me == 0 ) then
             print "(a12,2i4,4f12.3)", "ADV MAX. ", n, index  &
                      , d_2d(n,3,3,IOU_DAY) * PPBINV      &
                      ,  xn_adv(index,3,3,KMAX_MID)* PPBINV  &
                      ,  density(3,3), cfac(index,3,3)

            end if


         !u4 Daily maxima - short-lived
         !u4 ds - derived from jej's MAX_OH stuff
          case ( "MAXSHL" )

            forall ( i=1:limax, j=1:ljmax )
              d_2d( n, i,j,IOU_DAY) = max( d_2d( n, i,j,IOU_DAY), &
                                xn_shl(index,i,j,KMAX_MID)  &
                                    / (density(i,j)*MFAC) )
                                   !u4  / (roa(:,:,KMAX_MID,1)*MFAC) )
            end forall

            if ( MY_DEBUG .and. me == 0 ) then
               print *, "SHL:MAX. ", n, index  
               print *, "SHL:MFAC ",  MFAC
               print "(a12,2i4,4es12.3)", "SHL MAX. ", n, index  &
                      , d_2d(n,3,3,IOU_DAY)       &
                      ,  xn_shl(index,3,3,KMAX_MID)  &
                      ,  density(3,3), MFAC

            end if


          case ( "VOC", "TVOC" )

            call voc_2dcalc()

          case ( "EXT" )

          ! Externally set for IOU_INST (in other routines); so no new work needed
            if ( MY_DEBUG .and. me == 0 ) then
                 print *, "EXTDer:Externally set d_2d should already have values"
                 print *, "EXTDer:", typ
                 print *, "EXTDer: d_2d(2,2) is ", d_2d(n,2,2,IOU_INST)
            end if

          case  default

            if ( MY_DEBUG .and. me == 0 ) then
                 print *, "My_Deriv called for n=", n, "Type ",typ
                 print *, "My_Deriv index", index
                 print *, "My_Deriv avg? ", f_2d(n)%avg,  &
                                   " for nav ", nav_2d(n,IOU_INST)
                 print *, "Deriv: Length of f_2d is ", len(f_2d%class)
                 print *, "Deriv: f_2d class is ", f_2d(n)%class
             end if 

             call My_DerivFunc( n, typ, timefrac, density ) 

        end select


        !/** add to daily, monthly and yearly average, and increment counters
        !    /wdep, ddep not done here ???)
        !u4 - note that the MAXADV and MAXSHL needn't be summed here, but since
        !     the INST values are zero it doesn't harm, and the code is shorter

        d_2d(n,:,:,IOU_DAY )  = d_2d(n,:,:,IOU_DAY )  + d_2d(n,:,:,IOU_INST) 
        d_2d(n,:,:,IOU_MON )  = d_2d(n,:,:,IOU_MON )  + d_2d(n,:,:,IOU_INST) 
        d_2d(n,:,:,IOU_YEAR ) = d_2d(n,:,:,IOU_YEAR ) + d_2d(n,:,:,IOU_INST) 

        if ( f_2d(n)%avg ) nav_2d(n,:) = nav_2d(n,:) + 1

        !DEBUG print *, n, d_2d(n,2,2,1)
     end do   ! NDERIV_2D
     !/***** WET DEPOSITION **************************

!u7oozne - rsinstated
       do n = 1, NWDEP
!hf     do n=2,NWDEP !precip in Met_ml
        wdep(n,:,:,IOU_DAY )  = wdep(n,:,:,IOU_DAY )  + wdep(n,:,:,IOU_INST) 
        wdep(n,:,:,IOU_MON )  = wdep(n,:,:,IOU_MON )  + wdep(n,:,:,IOU_INST) 
        wdep(n,:,:,IOU_YEAR ) = wdep(n,:,:,IOU_YEAR ) + wdep(n,:,:,IOU_INST) 

        if ( MY_DEBUG .and. me == 0 ) then
          print *, "wet deposition  ",wdep(n,3,3,IOU_DAY)
        end if 

     end do  ! WET DEP.
     !/***** WET DEPOSITION **************************

     do n = 1, NDDEP

        ddep(n,:,:,IOU_DAY )  = ddep(n,:,:,IOU_DAY )  + ddep(n,:,:,IOU_INST) 
        ddep(n,:,:,IOU_MON )  = ddep(n,:,:,IOU_MON )  + ddep(n,:,:,IOU_INST) 
        ddep(n,:,:,IOU_YEAR ) = ddep(n,:,:,IOU_YEAR ) + ddep(n,:,:,IOU_INST) 

         if ( MY_DEBUG .and. me == 0 ) then
             print *, "dry deposition  ",ddep(n,3,3,IOU_DAY)
         end if 

     end do  ! DRY DEP.

     !/***** 3-D fields **************************


     do n = 1, NDERIV_3D

        index = f_3d(n)%index

        select case ( f_3d(n)%class )

         ! Simple advected species:
          case ( "ADV" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
            end forall

         case ( "BGN" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xn_bgn(index,i,j,k)
            end forall
!hf hmix xksig
         case ("XKSIG00", "XKSIG12" )

            forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
              d_3d( n, i,j,k,IOU_INST) = xksig(i,j,k)
            end forall

          case ( "VOC" )

            call voc_3dcalc()

          case  default

            print *, "Derived 3D class NOT  FOUND", n

        end select
     

       !/** add to monthly and yearly average, and increment counters
       !    ( no daily averaging done for 3-D fields so far).

        d_3d(n,:,:,:,IOU_MON ) = d_3d(n,:,:,:,IOU_MON ) &
                               + d_3d(n,:,:,:,IOU_INST)
        d_3d(n,:,:,:,IOU_YEAR) = d_3d(n,:,:,:,IOU_YEAR) &
                               + d_3d(n,:,:,:,IOU_INST)

        if ( f_3d(n)%avg )  nav_3d(n,:) = nav_3d(n,:) + 1
 
      end do
    end subroutine Derived
    !=========================================================================

    subroutine DerivedProds(text,dt)

    !/** DESCRIPTION
    !  Calculates chemical changes by comparing values before and  after 
    !  chemistry subroutine. Intended to be a more flexible version of the old 
    !  PRODO3  calculation

      character(len=*), intent(in) :: text  ! "Before" or "After"
      real,             intent(in) :: dt    ! timestep (s)

      real :: timefrac                      ! dt as fraction of hour (3600/dt)



      if (.not. any( f_3d%class == "PROD" ) ) return

      !bug: timefrac = 3600.0/dt
      timefrac = dt/3600.0
     !/***** 3-D fields **************************

     do n = 1, NDERIV_3D

        if ( f_3d(n)%class  == "PROD " ) then
           index = f_3d(n)%index

           select case ( text )

               case ( "Before" )   !! Initialise to xn_adv

                 forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                   d_3d( n, i,j,k,IOU_INST) = xn_adv(index,i,j,k)
                 end forall

               case ( "After" )    !! Calculate change

                 forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
                   d_3d( n, i,j,k,IOU_INST) = &
                      d_3d( n, i,j,k,IOU_INST) - xn_adv(index,i,j,k)
                 end forall

           end select
        end if
      end do
     
    end subroutine DerivedProds
    !=========================================================================

    subroutine ResetDerived(period)
      integer, intent(in) :: period   ! Either IOU_DAY or IOU_MON

       if ( period <= LENOUT2D ) then
           nav_wdep(:,period) = 0.0
           nav_ddep(:,period) = 0.0
           nav_2d  (:,period) = 0.0

           wdep(:,:,:,period) = 0.0
           ddep(:,:,:,period) = 0.0
           d_2d(:,:,:,period) = 0.0
       end if 


       if ( period <= LENOUT3D ) then
           nav_3d    (:,period) = 0.0
           d_3d(:,:,:,:,period) = 0.0
       end if

    end subroutine ResetDerived
 !=========================================================================

  subroutine Setup_VOC()
      !--------------------------------------------------------
      ! Searches through the advected species and colects the
      ! index and carbon content of nmhc species, as they were
      ! defined in GenOut_ml
      !
      ! Works for jej and ds chem
      !--------------------------------------------------------
       integer :: n
   
      do n = 1, NSPEC_ADV
        !u1 if ( species( MAP_ADV2TOT(n) )%nmhc == 1 ) then
        if ( species( NSPEC_SHL+n )%nmhc == 1 ) then
             nvoc = nvoc + 1
             voc_index(nvoc) = n
             voc_carbon(nvoc) = species( NSPEC_SHL+n )%carbons
        end if
      end do
  end subroutine Setup_VOC
 !=========================================================================

   subroutine voc_2dcalc()

    !/-- Sums up voc species using the indices defined earlier in Setup_VOCs

     ! We initialise d_2d first, the use a simple loop
     ! over voc. Some CPU could be saved by initialising
     ! with the 1st voc, then looping over 2, nvoc, but who cares...

      
      d_2d( n, 1:limax,1:ljmax,IOU_INST) =  0.0

      do ivoc = 1, nvoc

         index = voc_index(ivoc)           ! Gives which IXADV_ to use.
         forall ( i=1:limax, j=1:ljmax )
             d_2d( n, i,j,IOU_INST) = d_2d( n, i,j,IOU_INST)      &
                                    + xn_adv(index,i,j,KMAX_MID)  &
                                    * voc_carbon(ivoc) * cfac(index,i,j)
                               ! multiplied by nr. of C and "reduced to surface"
         end forall
      end do ! ivoc
   end subroutine voc_2dcalc

 !=========================================================================
   subroutine voc_3dcalc()

    !/-- as for voc_2dcalc

      d_3d( n, 1:limax,1:ljmax,1:KMAX_MID,IOU_INST) =  0.0

      do ivoc = 1, nvoc

         index = voc_index(ivoc)
         forall ( i=1:limax, j=1:ljmax, k=1:KMAX_MID )
             d_3d( n, i,j,k,IOU_INST) = d_3d( n, i,j,k,IOU_INST) + &
                     xn_adv(index,i,j,k)*voc_carbon(ivoc)
         end forall
      end do ! ivoc

   end subroutine voc_3dcalc
 !=========================================================================

end module Derived_ml










