! <BoundaryConditions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
!          Chemical transport Model>
!*****************************************************************************! 
!* 
!*  Copyright (C) 2007 met.no
!* 
!*  Contact information:
!*  Norwegian Meteorological Institute
!*  Box 43 Blindern
!*  0313 OSLO
!*  NORWAY
!*  email: emep.mscw@met.no
!*  http://www.emep.int
!*  
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!* 
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!* 
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*****************************************************************************! 
! ToDo: ????
! Check all allocates ?
! Make one subroutine to do all stuff rtepeated for adv and bgn??
! xn_Adv_changed doesn't include msic. Does this matter?
! Is 2nd call to Set_Misc... needed???
! Answer:
! Since bc_bgn and bc_adv are initiated each month,
! the second call to Set_Misc is needed. This can be re written
! but it is safer as it is....and do not cost too much
!____________________________________________________________________________
module BoundaryConditions_ml
!____________________________________________________________________________
! This module is the main driver module for defining and setting boundary
! conditions (bcs) for the chemical species. 
!
! The main code calls up these routines with just:
!
!           call BoundaryConditions(month)  !once per month, after say newmonth
!
!  On first call, this routine runs some intialisation routines in related
!  modules, reads the global data, and sets full 3-D concentration fields
!  of the advected and background concentration fields (xn_adv, xn_bgn).
!  On subsequent calls (my_first_call=.false.), the routine reads new
!  global input data, and rsets the concentrations at the top level and
!  lateral boundaries for advected species. For background species it reset 
!  full 3-D concentration fields. 
!
!  A bilinear interpolation routine is used to extrapolate from the coarser
!  global data to the EMEP model arrays, from the module Interpolations_ml.
!  Vertical interpolation is done on reading in the data.
!
!  The main related bc modules are:
!
!   1. UiO_ml.f90  - sets indices of global model data, e.g. IBC_O3, as well
!    as the number of global model fdata (NGLOB_BC).
!
!   2. My_bcmap_ml.f90 - assigns mappings, telling which EMEP species the
!      global model fdata are assigned to (bc2xn_adv, bc2xn_bgn arrays).
!  
! Language : F
! History  :
!
!
!ds Summer  2003: Added Mace Head corrections. Corrected bugs in twopi_year
!                 and vertical scaling.
!ds January 2002: modified for case with no BCs to be set (num_changed). 
!                 Re-formatted, replaced stop_test by mpi_abort   
!hf september-01: 2-dim mask changed into a 3-dim mask, and new restrictions 
!                 are made
!                 All adv species are only reset in top and edges, bgn species 
!                 in 3d.
!hf october-01:   Set_MiscBoundaryConditions renamed to MiscBoundaryConditions
! ds - December 2000-January 2001: added interpolations, removed txxlat, etc.
!   (See !REM at end for removals). Reorganised to put UiO stuff  in UiO_ml, 
!   "my" stuff (model-dependant) in My_BoundConditions_ml.
! jej - summer 2000 - original code called globinit.f
!
! ToDo:
! More flexible update times (update date array maybe?) ! (should make it 
! Dates_ml....) 
! Could restrict calculations of wt_00 etc. to limax, ljmax, but this was too 
! complicated in F90 (due to problems with dummary arguments and non-fixed 
! bounds)? To finish in January, I gave up!
!____________________________________________________________________________
! IMPORTANT NOTES:
! 1. The routines given here are constructed around the global model
! fields from the University of Oslo (T21) global model. In order
! to use other models as bcs  then usually these routines will have to be
! replaced by model-specific routines. The important thing is
! that the inputs and outputs from the routine are independant of the
! global module ufor one bc speciessed.
! 2. The routines make use of a "feature" of the model: that the concentration
! (xn) values  along boundaries are not changed due to advection or chemistry.
! Thus, bc values only need to be set once per month, firstly for the whole 3-D
! domain, then monthly for the sides and top.
! Background species must be reset in 3-D each month.
!_____________________________________________________________________________

  use CheckStop_ml,      only: CheckStop
  use Chemfields_ml,     only: xn_adv, xn_bgn, NSPEC_BGN  ! emep model concs.
!ds  use ChemSpecs_bgn_ml,        only :NSPEC_BGN
  use ChemSpecs_adv_ml         ! Lots, including NSPEC_ADV and IXADV_
  use GlobalBCs_ml,                only: &
          NGLOB_BC                 &  ! Number of species from global-model
          ,GetGlobalData           &  ! Sub., reads global data+vert interp.
          ,IBC_SO2, IBC_SO4, IBC_HCHO, IBC_CH3CHO &
          ,IBC_O3,IBC_HNO3,IBC_PAN,IBC_CO,IBC_C2H6   &
          ,IBC_C4H10, IBC_NO ,IBC_NO2,IBC_aNH4,IBC_aNO3,IBC_pNO3&
          ,IBC_H2O2,IBC_CH3COO2 &
          ,setgl_actarray
  use GridValues_ml,         only: gl, gb    &! lat, long
                                  ,sigma_mid &  !sigma layer midpoint
                                  ,debug_proc & !DSGC
                                  ,i_fdom, j_fdom  !u1 for testing
  use MetFields_ml               ,only : z_mid       ! height of half layers
  use ModelConstants_ml ,    only: KMAX_MID  &  ! Number of levels in vertical
                    ,NPROC   &    ! Number of processors
                    ,DEBUG_BCS, DEBUG_i, DEBUG_j, MasterProc
  use Par_ml,                only : &
           MAXLIMAX, MAXLJMAX, limax, ljmax, me &
          ,neighbor, NORTH, SOUTH, EAST, WEST   &  ! domain neighbours
          ,NOPROC&
          ,IRUNBEG,JRUNBEG,li1,li0,lj0,lj1

  implicit none
  private


  ! -- subroutines in this module:

  public  :: BoundaryConditions         ! call every month

  private :: Set_bcmap                  ! sets xn2adv_changed, etc.
  private :: MiscBoundaryConditions     ! misc bcs, not from global model.
  private :: Set_BoundaryConditions     ! assigns concentrations (xn) from bcs
  private :: My_bcmap                   ! sets bc2xn_adv, bc2xn_bc, and  misc_bc


  !/-- Allow different behaviour on 1st call - full 3-D asimilation done
   INCLUDE 'mpif.h'
   INTEGER STATUS(MPI_STATUS_SIZE),INFO

  logical, private, save :: my_first_call  = .true.  
  logical, public, parameter  :: BGN_2D = .false. !No 2d bgn species
  logical, private, parameter :: DEBUG_MYBC = .false.

  !/ - for debugging
  real    :: ppb = 1.0e-9

  !From Old My_BoundaryConditions_ml
 ! A. Set indices
 ! ===========================================================================
 ! For species which have constant mixing ratios:

 integer, public, parameter ::  NMISC_BC  =  2      ! H2, CH4   ! OC

 integer, public, parameter :: IBC_H2  = NGLOB_BC + 1   &
                              ,IBC_CH4 = NGLOB_BC + 2
                              !6s ,IBC_OC  = NGLOB_BC + 2

 integer, public, parameter ::  NTOT_BC  = NGLOB_BC + NMISC_BC

 ! We also need the array misc_bc to specify concentrations of these species:

 real, public, save, dimension(NGLOB_BC+1:NTOT_BC,KMAX_MID) :: misc_bc
!real, public, save, dimension(NGLOB_BC+1:NTOT_BC) :: misc_bc 

 ! B. Define mapping arrays
 ! ===========================================================================
 ! The mapping is done through the arrays bc2xn_adv and bc2xn_bgn, such that 
 ! the emep species are given along the x-dimension and the bc species along 
 ! the y.  e.g., the statement
 !
 !     bc2xn_adv(IBC_NOX,IXADV_NO2) = 0.55 
 !
 ! would assign the BC concentration of NOX  to the EMEP model concentration
 ! of NO2 after multiplication with a factor 0.55.
 !(The CTM2 concentration of NOx used as BC after a multilication 
 ! with a factor 0.55.)
 !_______________________________________________________________________

    real, public, save, dimension(NTOT_BC,NSPEC_ADV) :: bc2xn_adv  ! see above
    real, public, save, dimension(NTOT_BC,NSPEC_BGN) :: bc2xn_bgn ! see above


  ! Arrays for mapping from global bc to emep xn concentrations:
  ! ---------------------------------------------------------------------------

    integer, private,save, dimension(NTOT_BC) :: &
             bc_used         &! set to 1 if bc used
            ,bc_used_adv     &! set to 1 if bc used
            ,bc_used_bgn      ! set to 1 if bc used

    integer, private, save  :: &
        num_adv_changed        &       !Num. adv. species that have bc's
       ,num_used_adv                   !Max times a conc. from 
                                       !e.g CTM2 is used as bc 
    integer, private, save  :: &
        num_bgn_changed,num_used_bgn, &
        num_changed   !u1  - sum of adv and bgn 

   !In general "changed" means "bc used for this species"  

    logical, private, save, dimension(NSPEC_ADV)  :: &
        xn_adv_changed                 ! true if emep xn_adv changed by bcs
    logical, private, save, dimension(NSPEC_BGN )  :: &
        xn_bgn_changed                 ! true if emep xn_adv changed by bcs

    integer, private, save, dimension(NSPEC_ADV)  :: &
        spc_adv2changed            !index of advected specie is converted to 
                                   !index in the row of advected species 
                                   !that have bc 
    integer, private, save, dimension(NSPEC_BGN)  :: &
        spc_bgn2changed                 ! 

    integer, allocatable, dimension(:),save ::  &
        spc_changed2adv            !index of adv. specie that have bc is
                                   !converted to index in the row of
                                   !advected species 
    integer, allocatable, dimension(:),save :: &
        spc_changed2bgn

    integer, allocatable, dimension(:,:),save :: &
        spc_used_adv               !1.dimension(ibc) runs through bc adv.species
                                   !2.dim.(i) through those who get same conc.
                                   !from bc (eg i=1,2 for ibc=HNO3 when HNO3 
                                   !from CTM2 is used as bc both for HNO3 and SO4)
                                   !spc_used_adv gives the index in the row of 
                                   !advected species
    integer, allocatable, dimension(:,:),save ::  &
        spc_used_bgn


contains

 !-------
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine BoundaryConditions(year,iyr_trend,month)

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !** DESCRIPTION
 !   read in monthly-average global mixing ratios, and if found, collect the 
 !   data in  bc_adv, bc_bgn arrays for later interpolations
 !   (NB!!  if mixing ratio by mass the scale by molcular weight)
 !   ds- comment - so far no scaling is done, but this could be done
 !   in Set_bcmap with atomic weights.... for the future..
 !
 ! On the first call, we also run the setup-subroutines
 !
 !ds rv1.6.11 change: year is now obtained from the iyr_trend set in grun.pl
 ! Allows say runs with BCs for 2100 and met of 1990.
 !____________________________________________________________________________
  integer, intent(in) :: year        ! "meteorology" year        ds  rv1.6.11
  integer, intent(in) :: iyr_trend   ! "trend" year              ds  rv1.6.11 
  integer, intent(in) :: month
  integer :: ibc, iem, k,iem1,i,j     ! loop variables
  integer :: info              !  used in rsend  
  integer :: io_num            !  i/o number used for reading global data

  !/ data arrays for boundary data (bcs) - quite large, so NOT saved
  integer alloc_err1,alloc_err2,alloc_err3

  real, allocatable,dimension(:,:,:)   :: bc_data   ! for one bc species
  real, allocatable,dimension(:,:,:,:) :: bc_adv
  real, allocatable,dimension(:,:,:,:) :: bc_bgn

! (  nb dimensions correspond to:
!          IGLOB,JGLOB,KMAX_MID            :: bc_data
!          NSPEC_ADV,IGLOB,JGLOB,KMAX_MID  :: bc_adv
!          NSPEC_BGN,IGLOB,JGLOB,KMAX_MID  :: bc_bgn   )

  integer  :: iglobact,jglobact
  integer  ::  errcode
  integer, save :: idebug = 0, itest=1, i_test=0, j_test=0
  character(len=30) :: fname 

  if ( my_first_call ) then

    if (DEBUG_BCS) write(*,*) "FIRST CALL TO BOUNDARY CONDITIONS, me :", me
    if (DEBUG_BCS) write(*,*) "TREND YR, me ", iyr_trend, me

    call My_bcmap(iyr_trend)      ! assigns bc2xn_adv and bc2xn_bgn mappings
    call Set_bcmap()              ! assigns xn2adv_changed, etc.

    num_changed = num_adv_changed + num_bgn_changed   !u1

    if (DEBUG_BCS) write(*,*) "BCs: num_adv_changed: ", num_adv_changed
    if (DEBUG_BCS) write(*,*) "BCs: num_bgn_changed: ", num_bgn_changed
    if (DEBUG_BCS) write(*,*) "BCs: num     changed: ", num_changed

  end if ! first call
   if (DEBUG_BCS) write(*,*) "CALL TO BOUNDARY CONDITIONS, me, month :", me, month
   if (DEBUG_BCS) write(*,*) "TREND2 YR, me ", iyr_trend, me
  
  if ( num_changed == 0 ) then
      write(*,*) "BCs: No species requested"
      return
  end if

!MUST CONTAIN DECIDED DIMENSION FOR READ-IN DATA
! iglobac and jglobac are no the actual domains (the chosen domain)
! given in the same coord as the data we read - now 50*50
  call setgl_actarray(iglobact,jglobact)

  allocate(bc_data(iglobact,jglobact,KMAX_MID),stat=alloc_err1)
  call CheckStop(alloc_err1, "alloc1 failed in BoundaryConditions_ml") 

  ! - check if anything has changed before allocating:

      allocate(bc_adv(num_adv_changed,iglobact,jglobact,KMAX_MID), &
         stat=alloc_err2)
   call CheckStop(alloc_err2, "alloc2 failed in BoundaryConditions_ml") 
     bc_adv(:,:,:,:) = 0.0

      allocate(bc_bgn(num_bgn_changed,iglobact,jglobact,KMAX_MID), &
        stat=alloc_err3)
   call CheckStop(alloc_err3, "alloc3 failed in BoundaryConditions_ml") 
      bc_bgn(:,:,:,:) = 0.0


  errcode = 0
  if ( DEBUG_BCS .and. debug_proc ) then
       do i = 1, limax
            do j = 1, ljmax
                if ( i_fdom(i) == DEBUG_i .and. j_fdom(j) == DEBUG_j ) then
                    i_test = i
                    j_test = j
                end if
            end do
       end do
  end if



 !== BEGIN READ_IN OF GLOBAL DATA 

  do ibc = 1, NGLOB_BC

                   !=================================================
      if( MasterProc ) call GetGlobalData(year,iyr_trend,month,ibc,bc_used(ibc) &
                 ,iglobact,jglobact,bc_data,io_num,errcode)
                   !=================================================
      if ( DEBUG_BCS .and. MasterProc ) then
          write(*,*)'Calls GetGlobalData: year,iyr_trend,ibc,month,bc_used=', &
                    year,iyr_trend,ibc,month,bc_used(ibc)

      end if
      call CheckStop(ibc == 1 .and. errcode /= 0 ,&
             "ERRORBCs: GetGlobalData, failed in BoundaryConditions_ml") 

    !-- If the read-in bcs are required, we broadcast and use:

      if ( bc_used(ibc) > 0 ) then

            CALL MPI_BCAST(                              bc_data  ,8*iglobact*jglobact*KMAX_MID,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
          
          ! - set bc_adv, - advected species
          do i = 1, bc_used_adv(ibc)
             iem = spc_used_adv(ibc,i)
             iem1 = spc_adv2changed(iem)
             bc_adv (iem1,:,:,:) = bc_adv(iem1,:,:,:) + &
                                  bc2xn_adv(ibc,iem) * bc_data(:,:,:)
          end do

          ! set bc_bgn, - background (prescribed) species
          do i = 1, bc_used_bgn(ibc) 
             iem = spc_used_bgn(ibc,i)
             iem1 = spc_bgn2changed(iem)
             bc_bgn(iem1,:,:,:) = bc_bgn(iem1,:,:,:) + &
                                  bc2xn_bgn(ibc,iem) * bc_data(:,:,:)
          end do

        endif    ! bc_used
   end do  ! ibc

   if( MasterProc ) close(io_num)

   if ( my_first_call ) then

      idebug = 1
      if (DEBUG_BCS) print *, "RESET 3D BOUNDARY CONDITIONS", me
      !===================================
      !  -> 3-D arrays of new BCs
      call MiscBoundaryConditions(iglobact,jglobact,bc_adv,bc_bgn) 
      call Set_BoundaryConditions("3d",iglobact,jglobact        &
                ,bc_adv,bc_bgn)
      !===================================

      !===================================
      my_first_call = .false.
      !===================================

   else
      idebug = idebug + 1
      !===================================

      ! Call to set misc conditions added here on 7/08/01. This uses
      !     more CPU than needed since bc_adv is reset for the whole
      !     domain, so in future a "3d"/"lateral" mask could be used
      !     as for the normal boun. conds. 

      !  -> lateral (edge and top) arrays of new BCs
      call MiscBoundaryConditions(iglobact,jglobact,bc_adv,bc_bgn) 
      call Set_BoundaryConditions("lateral",iglobact,jglobact        &
                       ,bc_adv,bc_bgn)
      !===================================
   endif


   if ( DEBUG_BCS .and. debug_proc .and. i_test > 0  ) then
         write(6,"(a20,3i4,2f8.2)") "DEBUG BCS Rorvik", me, i,j,gl(i,j),gb(i,j)
         write(6,"(a20,3i4)")       "DEBUG BCS Rorvik DIMS", & 
                    num_adv_changed,iglobact,jglobact
         do k = 1, KMAX_MID
              write(6,"(a20,i4,f8.2)") "DEBUG O3  Debug-site ", k, &
                          xn_adv(IXADV_O3,i_test,j_test,k)/ppb
         end do
    end if ! DEBUG

   if ( DEBUG_BCS .and.  debug_proc ) then
       itest  = 1

       write(6,*) "BoundaryConditions: No CALLS TO BOUND Cs", &
                                        my_first_call, idebug
        !/** the following uses hard-coded  IXADV_ values for testing. 
        !    Remove later **/
      info = 1   ! index for ozone in bcs
       write(6,*) "BCs: bc2xn(info,itest) : ", bc2xn_adv(info,itest)
       do k = KMAX_MID, 1, -1
         write(6,"(a,2i3,f12.4)") "BCs: After Set_3d BOUND: me, itest: " , &
           me, itest, bc_adv(spc_adv2changed(itest),i_test,i_test,k)/ppb
       end do

      info = 43   ! index for NO in bcs
       write(6,*) "BCs: NSPECS: BC, ADV, BG, ", NTOT_BC, NSPEC_ADV, NSPEC_BGN
       write(6,*) "BCs: Number  of bc_used: ", sum(bc_used)
       write(6,*) "BCs: limax, ljmax",  limax, ljmax

       if ( NSPEC_BGN > 0 ) then
 
          do k = KMAX_MID, 1, -1
             write(6,"(a23,i3,e14.4)") "BCs NO :",   k  & 
                 ,xn_bgn(itest,i_test,j_test,k)/PPB 
          end do
       else 
             write(6,"(a)") "No SET BACKGROUND BCs"
      end if
    end if !  DEBUG

    deallocate(bc_data,stat=alloc_err1)
    call CheckStop(alloc_err1,"de-alloc1 failed in BoundaryConditions_ml") 

    if ( num_adv_changed > 0 ) then
       deallocate(bc_adv,stat=alloc_err2)
       call CheckStop(alloc_err2,"de-alloc2 failed in BoundaryConditions_ml") 
    end if

    if ( num_bgn_changed > 0 ) then
       deallocate(bc_bgn,stat=alloc_err3)
       call CheckStop(alloc_err3,"de-alloc3 failed in BoundaryConditions_ml") 
    end if

 end subroutine BoundaryConditions

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 ! From Old My_BoundaryConditions
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine My_bcmap(iyr_trend)    ! sets bc2xn_adv, bc2xn_bc, and  misc_bc
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    integer, intent(in) :: iyr_trend !ds Year for which BCs are wanted 
    real    :: ppb = 1.0e-9
    real :: trend_ch4  !ds rv1.6.11
    integer :: ii,i,j,k 
    real :: decrease_factor(NGLOB_BC+1:NTOT_BC) ! Decrease factor for misc bc's
                                      ! Gives the factor for how much of 
                                      ! the top-layer conc. that is left 
                                      ! at bottom layer

    real :: top_misc_bc(NGLOB_BC+1:NTOT_BC) ! Conc. at top of misc bc
!    real :: ratio_length(KMAX_MID)    ! Vertical length of the actual layer
                                      ! divided by length from midpoint of 
                                      ! layer 1 to layer KMAX_MID

    misc_bc            = 0.0  ! Initialise
    bc2xn_adv          = 0.0  ! Initialise
    bc2xn_bgn          = 0.0  ! Initialise

    ! Own (constant mixing ratio) boundary conditions **********

    ! NOTE - these species have to have the bc2xn_ indices set to 1.0 for either
    ! the advected or the background concentrations, in order that the 
    ! concentrations specified in misc_bc are transferred correctly into the 
    ! boundary conditions.
    !
    ! 18.09.01 -hf- misc bc'c as function of sigma possible 

        !ds top_misc_bc(IBC_CH4) = 1760.0 * ppb

        !ds set values of 1625 in 1980, 1780 in 1990, and 1820 in 2000. Interpolate
        ! between these for other years. Values from EMEP Rep 3/97, Table 6.2 for
        ! 1980, 1990, and from CDIAC (Mace Head) data for 2000.
 
        if ( iyr_trend >= 1990 ) then

            top_misc_bc(IBC_CH4) = 1780.0 + &
                        (iyr_trend-1990) * 0.1*(1820-1780.0)

        else 

            top_misc_bc(IBC_CH4) = 1780.0 * &
                        exp(-0.01*0.91*(1990-iyr_trend)) ! Zander,1975-1990
                     !exp(-0.01*0.6633*(1975-iyr_trend)) !Zander,1951-1975, check
        end if
        trend_ch4 = top_misc_bc(IBC_CH4)/1780.0  ! Crude for now.
        if ( MasterProc ) write(6,"(a20,i5,2f12.3)") "TREND CH4", &
               iyr_trend, trend_ch4, top_misc_bc(IBC_CH4)

        top_misc_bc(IBC_CH4)  =  top_misc_bc(IBC_CH4) * ppb

        top_misc_bc(IBC_H2)  =  600.0 * ppb

        !! top_misc_bc(IBC_OC) =   0.0 !!! 1.0 * ppb

        decrease_factor(IBC_H2) =0.0 !No increase/decrease with height
        decrease_factor(IBC_CH4)=0.0 !No increase/decrease with height


!a) Function of height(not included yet!):
!ratio_length(i,j,k)=(z_mid(i,j,1)-z_mid(i,j,k))/ &
!                      (z_mid(i,j,1)-z_mid(i,j,KMAX_MID))
!Replace sigma_mid with ratio_length and make misc_bc 4 dimentional
!
!b)Function of sigma_mid (I assume that top_misc_bc is given for
!top(sigma_bnd(1)=0.), and that at ground (sigma_bnd(KMAX_BND)=1.) the conc.
!is top_misc_bc -decrease_factor*top_misc_bc. Since I choose to set the 
! concentration as a factor of sigma_mid, the concentration in the lowest 
! grid cell will not be 
! excactly  top_misc_bc -decrease_factor*top_misc_bc, but close.

        do ii=NGLOB_BC+1,NTOT_BC
           do k=1,KMAX_MID
              misc_bc(ii,k) = top_misc_bc(ii) - &
                  top_misc_bc(ii)*decrease_factor(ii)*sigma_mid(k) 
              if ( MasterProc .and. DEBUG_MYBC ) then
                    write(*,"(a20,2es12.4,i4)")"height,misc_vert,k", &
                                               sigma_mid(k),misc_bc(ii,k),k
              endif
           enddo
        enddo

!        misc_bc(IBC_H2)     = 600.0 * ppb
!        misc_bc(IBC_OC)     =   1.0*ppb !!! 0.0

        bc2xn_adv(IBC_H2,  IXADV_H2)    = 1.0
        bc2xn_adv(IBC_CH4, IXADV_CH4)   = 1.0

        !/-- check, just in case we forgot something...!

        if ( DEBUG_MYBC  ) then
             print *, "In My_bcmap, NGLOB_BC, NTOT_BC is", NGLOB_BC, NTOT_BC
             do i = NGLOB_BC+1 , NTOT_BC
                print *, "In My_bcmap, sum-adv", i, " is", sum(bc2xn_adv(i,:))
                print *, "In My_bcmap, sum-bgn", i, " is", sum(bc2xn_bgn(i,:))
             end do
        end if ! DEBUG

        do i = NGLOB_BC+1 , NTOT_BC
           if ( sum(bc2xn_adv(i,:)) + sum(bc2xn_bgn(i,:)) /= 1.0 )then 
              WRITE(*,*) 'MPI_ABORT: ', "BCproblem - my" 
              call  MPI_ABORT(MPI_COMM_WORLD,9,INFO) 
           endif
        end do


     !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      include 'CM_BoundaryConditions.inc'
     !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!DSGC ! mappings for species from LOgan + obs model ***********
!DSGC
!DSGC  bc2xn_adv(IBC_O3      ,IXADV_O3      )   =   1.0
!DSGC  bc2xn_adv(IBC_HNO3    ,IXADV_HNO3    )   =   1.0
!.....
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 end subroutine My_bcmap
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   subroutine Set_bcmap()

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !_______________________________________________________________________
    ! - returns some 1-D arrays which say if a bc is used or if an 
    ! emep xn_adv or xn_bgn is affected.  This information is derived from
    ! the mapping arrays bc2xn_adv and bc2xn_bgn, such that the emep species 
    ! are given along the x-dimension and the bc species along the y.   
    ! e.g., the statement
    !
    !     bc2xn_adv(IBC_NOX,IXADV_NO2) = 0.55 
    !
    ! would assign the BC concentration of NOX  to the EMEP model concentration
    ! of NO2 after multiplication with a factor 0.55.
    !
    ! These arrays have been set in the My_BoundConditions_ml.f90 file.
    !_______________________________________________________________________

     integer :: ibc, iem   ! local loop variables
     integer :: i

    !-- bc_used set to one where a bc species is to be used.

    bc_used      = 0    ! Initialise
    bc_used_adv  = 0    ! Initialise
    bc_used_bgn  = 0    ! Initialise

    do ibc = 1, NTOT_BC
        if ( any( bc2xn_adv(ibc,:) > 0)  .or.  &
             any( bc2xn_bgn (ibc,:) > 0) ) then
                bc_used(ibc) = 1
        end if
      do iem = 1, NSPEC_ADV
         if(bc2xn_adv(ibc,iem) > 0)bc_used_adv(ibc) = bc_used_adv(ibc)+1
      enddo
      do iem = 1, NSPEC_BGN
         if(bc2xn_bgn(ibc,iem) > 0)bc_used_bgn(ibc) = bc_used_bgn(ibc)+1
      enddo
    end do ! ibc
    num_used_adv = maxval(bc_used_adv)
    num_used_bgn = maxval(bc_used_bgn)

    xn_adv_changed    = .false.    ! Initialise
    xn_bgn_changed    = .false.    ! Initialise
    num_adv_changed = 0
    num_bgn_changed = 0

    do iem = 1, NSPEC_ADV
         if ( any( bc2xn_adv(:,iem) > 0) ) then
             xn_adv_changed(iem) = .true.  
             num_adv_changed = num_adv_changed + 1
         end if
    end do ! iem

    do iem = 1, NSPEC_BGN
         if ( any( bc2xn_bgn(:,iem) > 0) ) then
              xn_bgn_changed(iem) = .true.  
              num_bgn_changed = num_bgn_changed + 1
         end if
    end do ! iem

    if ( DEBUG_BCS )  then 
       write(6,*) "TEST SET_BCMAP bc_used: "
       write(6,"(10i5)")  (bc_used(ibc),ibc=1, NTOT_BC)
    end if
    if (me==0) write(unit=6,fmt=*) "Finished Set_bcmap: Nbcused is ", sum(bc_used)

    allocate(spc_changed2adv(num_adv_changed))
    allocate(spc_changed2bgn(num_bgn_changed))
    i = 0
    spc_adv2changed = 0
    do iem = 1, NSPEC_ADV
      if(xn_adv_changed(iem))then
        i = i+1
        spc_changed2adv(i) = iem
        spc_adv2changed(iem) = i
      endif
    enddo
    i = 0
    spc_bgn2changed = 0
    do iem = 1, NSPEC_BGN
      if(xn_bgn_changed(iem))then
        i = i+1
        spc_changed2bgn(i) = iem
        spc_bgn2changed(iem) = i
      endif
    enddo
    allocate(spc_used_adv(NTOT_BC,num_used_adv))
    allocate(spc_used_bgn(NTOT_BC,num_used_bgn))

    spc_used_adv = 0
    spc_used_bgn = 0

    do ibc = 1, NTOT_BC

       if ( bc_used(ibc) > 0 ) then

          
          ! - set bc_adv, - advected species
          i = 0
          do iem = 1, NSPEC_ADV
             if ( bc2xn_adv(ibc,iem) > 0.0 ) then
                i = i+1
                spc_used_adv(ibc,i) = iem
             end if
          end do

          ! set bc_bgn, - background (prescribed) species
          i = 0
         do iem = 1, NSPEC_BGN 
             if ( bc2xn_bgn(ibc,iem) > 0.0 ) then
                i = i+1
                spc_used_bgn(ibc,i) = iem
             end if
          end do

        endif    ! bc_used
   end do  ! ibc

 end subroutine Set_bcmap
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine MiscBoundaryConditions(iglobact,jglobact,bc_adv,bc_bgn) 

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 ! - set bc_adv, bc_bgn,::::::::::::::::::::::::::::::::::::
 ! Note - this subroutine is a first draft and is only used for species
 ! which have constant mixing ratios to start with - here CH4, H2.
 ! I guess it should also work to set more complex variations (e.g.
 ! vertical gradients) though - this would then be better off in
 ! the module My_BoundaryConditions, but then we wouldn't have the
 ! bc_adv, bc_bgn arrays available   :-(


  integer, intent(in) ::   iglobact,jglobact 
  real, intent(inout), &
    dimension(num_adv_changed,iglobact,jglobact,KMAX_MID) :: bc_adv
  real, intent(inout), &
    dimension(num_bgn_changed,iglobact,jglobact,KMAX_MID) :: bc_bgn

  integer ::     itest  ! Used to specify species index
  integer :: ibc, iem,i,iem1,k   ! local loop variables


  if (NTOT_BC > NGLOB_BC) then
  do k=1,KMAX_MID
     do ibc = NGLOB_BC+1, NTOT_BC
        do i = 1,bc_used_adv(ibc)
          iem = spc_used_adv(ibc,i)
          iem1 = spc_adv2changed(iem)
          bc_adv(iem1,:,:,k) = misc_bc(ibc,k)
        enddo
        do i = 1,bc_used_bgn(ibc)
          iem = spc_used_bgn(ibc,i)
          iem1 = spc_bgn2changed(iem)
          bc_bgn(iem1,:,:,k) = misc_bc(ibc,k)
        enddo
     enddo
  enddo
endif

  itest = 1
  if (DEBUG_BCS .and. debug_proc) write(6,"(a50,i4,/,(5es12.4))") &
      "From MiscBoundaryConditions: ITEST (ppb): ",&
          itest, ((bc_adv(spc_adv2changed(itest),1,1,k)/1.0e-9),k=1,20) 

 end subroutine MiscBoundaryConditions

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine Set_BoundaryConditions(mode,iglobact,jglobact        &
                ,bc_adv,bc_bgn)  
 

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 ! Assign the global values to the interior model domain only first time
 ! (mode=3d) or only the edge and top boundaries are reset for other
 ! calls (mode=lateral). The emep model concentrations (xn) are only
 ! changed for those species where bcs are available, as given in the
 ! xn_adv_changed and xn_bgn_changed arrays.
 !
 ! Programming note:  we use F95 "forall" and logical masks to say which grid 
 ! squares are to be assigned.  This is probably slower than the more-explicit
 ! do-loops used previously, but is neater and shouldn't make too much 
 ! difference while this reset is done only once per month.

  character(len=*), intent(in) :: mode            ! "3d" or "lateral"
  integer,          intent(in) :: iglobact,jglobact
  real, intent(in), &
    dimension(num_adv_changed,iglobact,jglobact,KMAX_MID) :: bc_adv
  real, intent(in), &
    dimension(num_bgn_changed,iglobact,jglobact,KMAX_MID) :: bc_bgn
  integer, dimension(0:MAXLIMAX+1) :: &
         i150 !EMEP 150*150 coord of emep point (i)
  integer, dimension(0:MAXLJMAX+1) :: &
         j150 !EMEP 150*150 coord of emep point (j)
   logical, dimension(MAXLIMAX,MAXLJMAX,KMAX_MID) :: mask  
   integer :: i, j, k, n, ix, iy, ix1, iy1, itest


   if ( mode == "3d" ) then
      mask(:,:,:)       = .true.    ! We set everything 

   else if ( mode == "lateral" ) then
      mask(:,:,:)  = .false.        ! Initial


! Set edges (except on the top) 
! there may be no neighbor, but no external boundary (Poles in lat lon)
      if(neighbor(SOUTH) == NOPROC)   mask(:,1:(lj0-1),2:KMAX_MID)     = .true.
      if(neighbor(NORTH) == NOPROC)   mask(:,(lj1+1):ljmax,2:KMAX_MID) = .true.
      if(neighbor(EAST)  == NOPROC)   mask((li1+1):limax,:,2:KMAX_MID) = .true.
      if(neighbor(WEST)  == NOPROC)   mask(1:(li0-1),:,2:KMAX_MID)     = .true.

      mask(:,:,1) = .true.        !Set top layer
   else
      call CheckStop("BCs:Illegal option failed in BoundaryConditions_ml") 
   endif

   !** Set concentrations (xn) from boundary conditions (bcs) 

   ! Note on domains: although the geographical stuff has been specified
   ! for the whole MAXLIMAX,MAXLJMAX grid, the interpolations take time
   ! and are only needed for the sub-domain actually used, i.e. for
   ! limax, ljmax.


   !a) Advected species


   do k = 1, KMAX_MID
      do j = 1, ljmax
         do i = 1,limax
            if (mask(i,j,k)) then                 
               do n = 1, num_adv_changed
                  xn_adv(spc_changed2adv(n),i,j,k) =   &
                       bc_adv(n,(i_fdom(i)-IRUNBEG+1),(j_fdom(j)-JRUNBEG+1),k)
               end do
            endif
         end do
      end do
   end do
   



   !b) Non-advected background species

   forall(i=1:limax, j=1:ljmax, k=1:KMAX_MID, n=1:num_bgn_changed)

          xn_bgn(spc_changed2bgn(n),i,j,k) = & 
             bc_bgn(n,(i_fdom(i)-IRUNBEG+1),(j_fdom(j)-JRUNBEG+1),k)

   end forall    


 end subroutine Set_BoundaryConditions    ! call every 3-hours
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module BoundaryConditions_ml
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
