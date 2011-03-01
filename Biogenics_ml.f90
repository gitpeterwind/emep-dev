! <Biogenics_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
module Biogenics_ml

  !/-- Reads in BVOC emisions factors 
  !
  !     1) From LPJ or defaults globally
  !
  !     2) from local file if available (e.g. Europe, used by default)
  !
  !   Terminology:
  !
  !    LCC = land cover class, e.g. DF = decid forest
  !
  !    EF = Emission factor at 30 deg C, full sunlight (1000 uE)
  !         ug/g/hr
  !
  !    Em = Emissions, = EF * LAI * fraction of grid
  !         ug/m2(ground)/hr
  !
  !    The code will assign EFs from the local data if available, otherwise
  !    use defaults which have been readfrom Inputs_Landuse, Eiso, Emtp, Emtl. 
  !    Note that we use emissions per m2 of vegetation, not per 
  !    m2 of grid. This lets the model use the landcover from any veg-map
  !    without having to make this consistent with the EF maps - the latter
  !    are regarded as smoothly varying fields which can be interpolated 
  !    by the ReadField_CDF interpolation routines. No need to worry about
  !    conserving these very imperfect numbers accurately ;-)
  !
  !  ( MEGAN assumes rates for 5 m2/m2 LAI.... think about later)
  !
  !    Dave Simpson, 2010
  !---------------------------------------------------------------------------


  use CheckStop_ml,      only: CheckStop
  use GridValues_ml    , only : xm2, gb, &
          i_fdom,j_fdom,debug_proc,debug_li,debug_lj
  use Io_ml            , only : IO_FORES, open_file, ios, Read2DN, PrintLog, datewrite
  use KeyValue_ml,       only : KeyVal,KeyValue
  use LandDefs_ml,       only: LandType, LandDefs
  use LandPFT_ml,        only: MapPFT_LAI, pft_lai
  use Landuse_ml,        only : LandCover
  use LocalVariables_ml, only : Grid  ! -> izen, DeltaZ
  use ModelConstants_ml, only : NPROC, MasterProc, TINY, &
                           USE_PFT_MAPS, NLANDUSEMAX, IOU_INST, & 
                           KT => KCHEMTOP, KG => KMAX_MID, & ! DSBIO
                           DEBUG_BIO, BVOC_USED, USE_BVOC_2010, MasterProc
  use NetCDF_ml,         only : ReadField_CDF, Out_netCDF,  Real4
  use OwnDataTypes_ml,  only : Deriv, TXTLEN_SHORT
  use Par_ml   , only :  MAXLIMAX,MAXLJMAX,MSG_READ1,li0,li1,lj0,lj1,me
  use PhysicalConstants_ml,  only :  AVOG, GRAV ! , PI
  use Radiation_ml,          only : PARfrac, Wm2_uE
  use Setup_1dfields_ml,     only : rcbio  
  use SmallUtils_ml, only : find_index
  use TimeDate_ml,       only :  current_date
  implicit none
  private

  !/-- subroutines
  public ::  Init_BVOC
  private :: Get_LCinfo
  public ::  GetEuroBVOC
  private :: MergedBVOC
  public ::  setup_bio
  public ::  SetDailyBVOC
  private :: TabulateECF
  private :: Export_Bio

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO
  integer, public, parameter :: N_ECF=2, ECF_ISOP=1, ECF_TERP=2
  integer, public, parameter :: NBIO_DEF=3, BIO_ISOP=1, BIO_MTP=2, BIO_MTL=3
  integer, public, parameter :: BIO_TERP=2 ! Used for final emis, sum of MTP+MTL
  integer, public, save ::  last_bvoc_LC   !max index land-cover with BVOC (min 4)

 ! Set true if LCC read from e.g. Euro_.nc:
 ! (Currently for 1st four LCC)
  logical, private, dimension(NLANDUSEMAX), save :: HaveLocalEF 

  !TEST
  real, public, save, dimension(MAXLIMAX,MAXLJMAX,size(BVOC_USED)) :: &
  !TEST real, public, save, dimension(MAXLIMAX,MAXLJMAX,2) :: &
      emforest    & !  Gridded standard (30deg. C, full light) emissions
     ,EmisNat       !  will be transferred to d_2d emis sums


  !standard emission factors per LC  for LAI=5 m2/m2
  !Need to dimension later for Emtp, Emtl, last_bvoc_LC
  real, private, save, allocatable, dimension(:,:,:,:) :: &
     bvocEF       !  Gridded std. emissions per PFT

  !standard emission factors per LC  for daily LAI
  real, private, save, dimension(MAXLIMAX,MAXLJMAX,size(BVOC_USED)) :: &
     day_embvoc   !  emissions scaled by daily LAI

  logical, private, dimension(MAXLIMAX,MAXLJMAX) :: EuroMask

  !/-- Canopy environmental correction factors-----------------------------
  !
  !    - to correct for temperature and light of the canopy 
  !    - from Guenther's papers. (Limit to 0<T<40 deg C.)

  real, public, save, dimension(N_ECF,40) :: canopy_ecf  ! Canopy env. factors
                                                        
  integer, public, parameter :: IQ_DMS = 35  ! code for DMS emissions
  logical, public, save :: first_dms_read

  contains
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    subroutine Init_BVOC()

!    Read natural BVOC emission potentials
!-----------------------------------------------------------------------------
!   Emissions now read from 50x50 landuse file, forests.dat, derived from 
!   landuse.mar2004. Emission rates now based upon Simpson et al., 1999, JGR,
!   Vol 104, D7, 8113-8152.


    integer i, j, n ,d, info, alloc_err
    real sumland

    real ::  bvocsum, bvocsum1

    ! Specify the assumed coords and units - Read2DN will check that the data
    ! conform to these.
    type(keyval), dimension(2) :: CheckValues = &
        (/ keyval("Units","ug/m2/h"), keyval("Coords","ModelCoords") /)

      if ( size(BVOC_USED) == 0 ) then
        call PrintLog("No Biogenic Emissions ", MasterProc)
        return
      end if
      if( USE_BVOC_2010 ) then
         call PrintLog("BVOC Emission Style: 2010", MasterProc)
      else
         call PrintLog("BVOC Emission Style: 2003", MasterProc)
      end if

   !====================================
   ! Initialise chemical emission rates to zero. In rest of  code
   ! only k=KMAX_MID will be set

     rcbio(:,:) = 0.0  
   !====================================
 
    call TabulateECF()   ! Tabulates temp functions
   !====================================

    call Get_LCinfo() ! Gets landcover info, last_bvoc_LC

    allocate(  bvocEF(MAXLIMAX,MAXLJMAX,last_bvoc_LC,size(BVOC_USED)),&
        stat=alloc_err )
    call CheckStop( alloc_err , "bvocEF alloc failed"  )

    bvocEF(:,:,:,:) = 0.0

   !========= Read in Standard (30 deg C, full sunlight emissions factors = !
   ! Remember factor used in Emissions_ml to get from ug/m2/s
   ! to molecules/cm2/s  (needs more documentation!)

   !====================================

 if ( USE_BVOC_2010  .or. DEBUG_BIO ) then
    call GetEuroBVOC()
   !====================================

   !====================================
   ! Merges Local and global/defaults, and scales with land-cover area
   ! Emissions factors shoudl now by ug/m2(grid)/h

    call MergedBVOC() 
 end if ! BVOC_2010
   !====================================

    emforest = 0.0
    if ( USE_BVOC_2010 .eqv. .false. ) then
      call Read2DN("Inputs.BVOC",2,emforest,CheckValues)
      call CheckStop( minval(emforest) < 0.0, "Negative BVOC emis!")
    end if
   !========================================================================!

      if( debug_proc .and. DEBUG_BIO ) then
          write(*,"(a8,i3,2i4,9f9.2)") "BIONEW ", size(BVOC_USED), &
              i_fdom(debug_li), j_fdom(debug_lj), &
            ( emforest(debug_li,debug_lj,i), i=1,2) !!!!,&
!BIO           ( emforest(debug_li,debug_lj,i), i=1,size(BVOC_USED)) !!!!,&
!Z               (euro_bvoc(debug_li,debug_lj,i,BIO_ISOP), i=1,size(VegName))
      end if

     !output sums. Remember that "shadow" area not included here.
      do i = 1,  2 !! size(BVOC_USED) 
         bvocsum   = sum ( emforest(li0:li1,lj0:lj1,i) )
         CALL MPI_ALLREDUCE(bvocsum,bvocsum1, 1, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
         if ( MasterProc  ) write(6,"(a20,i4,2es12.4)") &
              'Biogenics_ml, ibio, sum1',i, bvocsum, bvocsum1
      end do

! No, we would need the day:Embvoc every ady to get the equivalent annual
! sum
!      if( DEBUG_BIO .and. USE_BVOC_2010 ) then
!            do n = 1, size(BVOC_USED)
!               bvocsum   = sum ( bvocEF(li0:li1,lj0:lj1,n) )
!               CALL MPI_ALLREDUCE(bvocsum,bvocsum1, 1, &
!                 MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, INFO) 
!               if ( MasterProc  ) &
!                write(6,"(a20,i4,2es12.4)") 'Biogenics_ml, Merge',n, &
!                  bvocsum, bvocsum1
!            end do !n
!      end if

   end subroutine Init_BVOC
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine Get_LCinfo() ! Gets landcover info, last_bvoc_LC
      integer :: iL

        do iL= 1, size(LandType(:)%pft )
            
            if( LandDefs(iL)%Eiso > 0 ) last_bvoc_LC = iL
            if( LandDefs(iL)%Emtp > 0 ) last_bvoc_LC = iL
            if( LandDefs(iL)%Emtl > 0 ) last_bvoc_LC = iL

         end do
         if( MasterProc ) write(*,*) " LAST BVOC LC (pre 4):", last_bvoc_LC

       ! We need at least 4 for CF, DF, NF, BF in Euro file
         last_bvoc_LC =  max(last_bvoc_LC, 4 ) 

    end subroutine Get_LCinfo
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   subroutine  GetEuroBVOC()

!.....................................................................
!**    DESCRIPTION:

!    Reads the processed LPJ-based LAIv and BVOC emission potentials.
!    The LPJ data have been merged into 4 EMEP forest classes and two
!    other veg, for either C3 or C4 vegetation.
!    Normed_LAIv is relative LAI, with max value 1.0


    real    :: loc(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
    logical :: my_first_call = .true.
    integer ::  n, pft, ivar, iVeg, iEmis, ibvoc
    character(len=1000) :: varname
    character(len=2), dimension(4) :: VegName = (/ "CF", "DF", "NF", "BF" /)


       varname = "Fake"
       HaveLocalEF(:) = .false.

       call ReadField_CDF('LOCAL_BVOC.nc',varname,&
           loc,1,interpol='zero_order',needed=.true.,debug_flag=.true.)

       if( debug_proc ) write(*,*)  "LOCAL_BVOC i,j fake ", &
          loc(debug_li, debug_lj)

     do iVeg = 1, size(VegName)
       ibvoc = find_index( VegName(iveg), LandDefs(:)%code )
       HaveLocalEF(ibvoc) = .true.
       do iEmis = 1, size(BVOC_USED)
          varname = trim(BVOC_USED(iEmis)) // "_" // trim(VegName(iVeg))
          call ReadField_CDF('LOCAL_BVOC.nc',varname,&
             loc,1,interpol='zero_order',needed=.true.,debug_flag=.false.)
         if( debug_proc ) write(*, "(2a,f12.3,3i2)") "EURO-BVOC:E ", &
             trim(varname), loc(debug_li, debug_lj), iVeg, ibvoc, iEmis
         bvocEF(:,:,ibvoc,iEmis) = loc(:,:)
       end do
     end do

      ! Make a mask where we can use the local bvoc. Should be the same from
      ! all EFs, since only non-def areas set to -999, otherwise zero or +

      where(loc>-1.0)
          EuroMask = .true.
      else where 
          EuroMask = .false.
      end where

  end subroutine GetEuroBVOC
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   subroutine MergedBVOC()

      integer :: i, j, n, nlu, iL, iiL
      integer :: pft
      real :: biso, bmt    !  Just for printout
      logical :: use_local, debug_flag

      if( MasterProc ) write(*,*) "Into MergedBVOC"

      do i = 1, MAXLIMAX
      do j = 1, MAXLJMAX

        nlu = LandCover(i,j)%ncodes

        use_local = EuroMask(i,j) 
        debug_flag = ( debug_proc .and. debug_li == i .and. debug_lj == j )

        LULOOP: do iiL= 1, nlu

            iL      = LandCover(i,j)%codes(iiL)
            pft     = LandType(iL)%pft

           if( debug_flag ) then
               write(*,"(a,2i7,2L)") &
                   "TryMergeBVOC" //trim(LandDefs(iL)%name), iL, pft, &
                     use_local, HaveLocalEF(iL)
               write(*,*) "TryMergeBVOC last_bvoc_LC = ", last_bvoc_LC
           end if

           if( use_local .and. HaveLocalEF(iL) ) then 

                ! Keep EFs from EuroBVOC
                if( debug_flag ) write(*,*) "MergeBVOC: Inside local"

           else if ( iL <= last_bvoc_LC ) then ! otherwise use defaults
               bvocEF(i,j,iL,BIO_ISOP) = LandDefs(iL)%Eiso * LandDefs(iL)%BiomassD 
               bvocEF(i,j,iL,BIO_MTP)  = LandDefs(iL)%Emtp * LandDefs(iL)%BiomassD
               bvocEF(i,j,iL,BIO_MTL)  = LandDefs(iL)%Emtl * LandDefs(iL)%BiomassD
                if( debug_flag ) write(*,*) "MergeBVOC: Outside local", iL
           else
                if( debug_flag ) write(*,*) "MergeBVOC: Outside LCC", iL
           end if


           if( debug_flag ) then

              biso = 0.0
              bmt  = 0.0
              if ( iL <= last_bvoc_LC ) then
                biso   = bvocEF(i, j,iL, BIO_ISOP) 
                bmt    = bvocEF(i, j,iL, BIO_TERP) 
              end if
              write(*,"(a,2i4,2L,f9.4,9f10.3)") "MergeBVOC", &
                 iL, pft,  use_local, HaveLocalEF(iL),  &
                   LandCover(i,j)%fraction(iiL), biso, bmt, LandDefs(iL)%Eiso, &
                     LandDefs(iL)%Emtp, LandDefs(iL)%Emtl
           end if 
        end do LULOOP
      end do
      end do

   end subroutine MergedBVOC
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine SetDailyBVOC(daynumber)

      integer, intent(in) :: daynumber
      integer, save :: last_daynumber = -999, alloc_err
      integer :: i, j, n, nlu, iL, iiL, ibvoc
      integer :: pft  ! pft associated with LAI data
      real :: LAIfac  ! Multiplies by land-fraction and then divides by 5
      real :: b    !  Just for printout
      logical :: use_local, debug
      logical :: my_first_call = .true.
      real, allocatable, dimension(:,:) ::  workarray

      if( MasterProc ) write(*,"(a,3i5)") "Into SetDailyBVOC", &
            daynumber, last_daynumber, last_bvoc_LC

      if ( daynumber == last_daynumber ) return
      last_daynumber = daynumber

      do i = 1, MAXLIMAX
      do j = 1, MAXLJMAX

        nlu = LandCover(i,j)%ncodes

        day_embvoc(i,j,:) = 0.0
        debug = ( DEBUG_BIO .and. debug_proc .and.  &
                   debug_li == i .and. debug_lj == j )

        LULOOP: do iiL= 1, nlu

            iL      = LandCover(i,j)%codes(iiL)

            if ( iL >  last_bvoc_LC ) cycle
            if ( LandCover(i,j)%LAI(iiL)<1.0e-5 ) cycle ! no veg:
             
              LAIfac = LandCover(i,j)%LAI(iiL)/LandDefs(IL)%LAImax
              LAIfac= min(LAIfac, 1.0)
              LAIfac = LAIfac * LandCover(i,j)%fraction(iiL)

              do ibvoc = 1, size(BVOC_USED) 
                day_embvoc(i,j,ibvoc) = day_embvoc(i,j,ibvoc) + &
                   LAIfac * max(1.0e-10,bvocEF(i,j,iL,ibvoc))
                   !done above LandCover(i,j)%fraction(iiL) *  &
              end do

              if ( debug ) then
                 b = 0.0
                 if ( iL <= last_bvoc_LC ) b = bvocEF(i, j,iL, BIO_ISOP)
                 write(*,"(a,a10,2i5,f9.5,2f7.3,8f10.3)") "SetBVOC", &
                  trim(LandDefs(iL)%name), daynumber, iL, &
                   LandCover(i,j)%fraction(iiL), &
                   LandCover(i,j)%LAI(iiL), LandDefs(iL)%LAImax, b,&
                     ( day_embvoc(i, j, ibvoc), ibvoc = 1, size(BVOC_USED) ) 
                  
              end if
        end do LULOOP
      end do
      end do

      if ( DEBUG_BIO ) then
         if ( my_first_call  ) then ! print out 1st day
              allocate(  workarray(MAXLIMAX,MAXLJMAX),&
                        stat=alloc_err )
              write(*,"(a,i3)") "workarray success????", alloc_err
              call CheckStop( alloc_err , "workarray alloc failed"  )
                  workarray(:,:) = day_embvoc(:,:,1)
                  call Export_Bio("Eiso", workarray )
                  workarray(:,:) = day_embvoc(:,:,2)
                  call Export_Bio("Emt", workarray )
              deallocate(  workarray )
         end if
       end if
       my_first_call = .false.

   end subroutine SetDailyBVOC
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    subroutine TabulateECF()
    !/-- Tabulate canopy environmental correction factors
    !-----------------------------------------------------
    ! ag suffix  = Alex Guenther's parameter values

    real    :: agts, agr, agtm, agct1, agct2, agct ,itk, fac
    integer :: it, i

    agts = 303.
    agr = 8.314
    agtm = 314.
    agct1 = 95000.
    agct2 = 230000.

    do it = 1,40
      itk = it + 273.15
      agct = exp(agct1*(itk - agts)/(agr*agts*itk)) / &
                (1. + exp(agct2*(itk - agtm)/(agr*agts*itk)))

      canopy_ecf(ECF_ISOP,it) = agct

      ! Terpenes
      agct = exp( 0.09*(itk-agts) )
      canopy_ecf(ECF_TERP,it) = agct

      !?? for terpene fac = 0.5*fac(iso): as mass terpene = 2xmass isoprene

      if(DEBUG_BIO  .and.  MasterProc ) &
             write(6,"(A12,i4,5g12.3)") 'Biogenic ecfs: ', &
                  it, ( canopy_ecf(i,it), i=1, N_ECF )
    end do
    end subroutine TabulateECF
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  subroutine setup_bio(i,j)
  !
  !---- assign isoprene rates  ------------------------------------------------
  !
  !  So far, assigns isoprene using surface (2m) temperature, and for all
  !  zenith angles <90. Should include light dependance at some stage
  !
  !  Output : rcbio added to rcemis - isoprene emissions for 1d column
  !
  !  Called from setup_ml, every  advection step.
  !----------------------------------------------------------------------------

  integer, intent(in) ::  i,j

  integer la,it2m,n,k,base,top,iclcat
  real :: E2003, E2010 , MT2003, MTP, MTL

! To get from ug/m2/h to molec/cm3/s
! ug -> g  1.0e-9; g -> mole / MW; x AVOG
  real, save :: biofac_ISOP = 1.0e-12*AVOG/64.0/3600.0   ! needs /Grid%DeltaZ
  real, save :: biofac_TERP = 1.0e-12*AVOG/136.0/3600.0  ! needs /Grid%DeltaZ

 ! Light effects added for isoprene emissions

  real            :: par   ! Photosynthetically active radiation
  real            :: cL    ! Factor for light effects
  real, parameter :: &
      CL1 = 1.066  , &    ! Guenther et al's params
      ALPHA = 0.0027      ! Guenther et al's params
  real :: tmpdz ! dzTEST

  if ( size(BVOC_USED) == 0  ) return   ! e.g. for ACID only


  it2m = nint( Grid%t2C - TINY )
  it2m = max(it2m,1)
  it2m = min(it2m,40)

  rcbio(:,KG) = 0.0
  EmisNat(i,j,:) = 0.0
  MTP = 0.0
  MTL = 0.0


  if ( Grid%izen <= 90) then ! Isoprene in daytime only:

     ! Light effects from Guenther G93

      par = (Grid%Idirect + Grid%Idiffuse) * PARfrac * Wm2_uE

      cL = ALPHA * CL1 * par/ sqrt( 1 + ALPHA*ALPHA * par*par)

     ! E in ug/m2/h

      ! both just for print out
      E2003 = emforest(i,j,BIO_ISOP)  *canopy_ecf(BIO_ISOP,it2m) * cL 
      E2010 = day_embvoc(i,j,BIO_ISOP)*canopy_ecf(BIO_ISOP,it2m) * cL 
      ! Add light-dependent terpenes to pool-only
      if(BIO_TERP > 0) MTL = &
             day_embvoc(i,j,BIO_MTL)*canopy_ecf(ECF_TERP,it2m)*cL

     !  molecules/cm3/s
     ! And we scale EmisNat to get units kg/m2 consistent with
     ! Emissions_ml (snapemis).  ug/m2/h -> kg/m2/s needs 1.0-9/3600.0. 
     if ( USE_BVOC_2010 ) then

        E2010 = day_embvoc(i,j,BIO_ISOP)*canopy_ecf(BIO_ISOP,it2m) * cL 
        rcbio(BIO_ISOP,KG) = E2010 * biofac_ISOP/Grid%DeltaZ
        EmisNat(i,j,BIO_ISOP)= E2010 * 1.0e-9/3600.0
     else 
        E2003 = emforest(i,j,BIO_ISOP)  *canopy_ecf(BIO_ISOP,it2m) * cL 
        rcbio(BIO_ISOP,KG) = E2003     * biofac_ISOP/Grid%DeltaZ
        EmisNat(i,j,BIO_ISOP)= E2003 * 1.0e-9/3600.0
     end if

  endif ! daytime

  if ( BIO_TERP > 0 ) then
    if ( USE_BVOC_2010 ) then
     ! add pool-only terpenes rate;
        MTP = day_embvoc(i,j,BIO_MTP)*canopy_ecf(ECF_TERP,it2m)
        rcbio(BIO_TERP,KG)    = (MTL+MTP) * biofac_TERP/Grid%DeltaZ
        EmisNat(i,j,BIO_TERP) = (MTL+MTP) * 1.0e-9/3600.0
     else 

       MT2003 = emforest(i,j,BIO_TERP)*canopy_ecf(BIO_TERP,it2m) 
       rcbio(BIO_TERP,KG) =  MT2003 * biofac_TERP/Grid%DeltaZ
       EmisNat(i,j,BIO_TERP)=  MT2003 * 1.0e-9/3600.0  ! Nov4th
     end if
  end if
 
     if ( DEBUG_BIO .and. debug_proc .and.  &
            i==debug_li .and. j==debug_lj .and. &
         current_date%seconds == 0 ) then

      ! just for print out
      ! both just for print out
      E2003 = emforest(i,j,BIO_ISOP)  *canopy_ecf(BIO_ISOP,it2m) * max(cL,0.0)
      E2010 = day_embvoc(i,j,BIO_ISOP)*canopy_ecf(BIO_ISOP,it2m) * max(cL,0.0) 
       MTP = day_embvoc(i,j,BIO_MTP)*canopy_ecf(ECF_TERP,it2m)
       MT2003 = emforest(i,j,BIO_TERP)*canopy_ecf(BIO_TERP,it2m) 

        call datewrite("DBIO NatEmis env ", it2m,  &
          (/ max(par,0.0), max(cL,0.0), &
            canopy_ecf(BIO_ISOP,it2m),canopy_ecf(BIO_TERP,it2m) /) )
        call datewrite("DBIO NatEmis for ", &
          (/ emforest(i,j,1),  emforest(i,j,2) /) )
        call datewrite("DBIO NatEmis EIcmp ", (/  E2003, E2010 /) ) 
        call datewrite("DBIO NatEmis EMTcmp ", (/  MT2003, MTP /) ) 
        call datewrite("DBIO NatEmis rc  ", &
          (/ rcbio(BIO_ISOP,KG), rcbio(BIO_TERP,KG) /) )
        call datewrite("DBIO NatEmis EmisNat  ", EmisNat(i,j,:) )

!        write(*,"(a,3i3,f6.1,3f7.3,9es9.2)") "DBIO NatEmis-2003 ", &
!        write(*,"(a,3i3,f6.1,3f7.3,9es9.2)") "DBIO NatEmis ", &
!          current_date%month, current_date%day, current_date%hour, &
!          par, cL, canopy_ecf(BIO_ISOP,it2m),canopy_ecf(BIO_TERP,it2m), E2003, E2010, MTP, MTL, &
!          rcbio(BIO_ISOP,KG), rcbio(BIO_TERP,KG),  EmisNat(i,j,BIO_ISOP), EmisNat(i,j,BIO_TERP)

     end if


  end subroutine setup_bio

  !----------------------------------------------------------------------------
 subroutine Export_Bio(name, array)
    real, dimension(:,:), intent(in) :: array
    character(len=*), intent(in) :: name
    character(len=60) :: fname
    type(Deriv) :: def1 ! definition of fields
    
    def1%class='Biogenics' !written
    def1%avg=.false.      !not used
    def1%index=0          !not used
    def1%scale=1.0      !not used
!FEB2011    def1%inst=.true.      !not used
!FEB2011    def1%year=.false.     !not used
!FEB2011    def1%month=.false.    !not used
!FEB2011    def1%day=.false.      !not used
    def1%name=trim(name)   ! eg 'EmisPot'        !written
    def1%unit='ug/m2/h'       !written
    
    fname = "OUTBIO_" // trim(name) // ".nc"

    !Out_netCDF(iotyp,def1,ndim,kmax,dat,scale,CDFtype,ist,jst,ien,jen,ik,fileName_given)

    if(MasterProc) write(*,*) "TEST EXPORT BIO:"//trim(fname),  maxval(array)
    call Out_netCDF(IOU_INST,def1,2,1, array,1.0,&
           CDFtype=Real4,fileName_given=fname)
  end subroutine Export_Bio


end module Biogenics_ml
