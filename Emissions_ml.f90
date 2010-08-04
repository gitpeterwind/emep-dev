! <Emissions_ml.f90 - A component of the EMEP MSC-W Unified Eulerian
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
!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                    module Emissions_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!+ Calls up emission read/set routines
!  This routine interfaces the stand-alone emission-file reading routines
!  with the 3D model.
!_____________________________________________________________________________

  use Biogenics_ml, only: first_dms_read,IQ_DMS,emnat,emforest
  use CheckStop_ml,only : CheckStop
  use ChemSpecs_shl_ml, only: NSPEC_SHL
  use ChemSpecs_tot_ml, only: NSPEC_TOT,NO2
  use ChemChemicals_ml, only: species
  use Country_ml,    only : NLAND,Country_Init,Country, IC_NAT
  use EmisDef_ml, only : NSECTORS & ! No. sectors
                     ,NEMISLAYERS & ! No. vertical layers for emission
                     ,NCMAX &       ! Max. No. countries per grid
                     ,FNCMAX &      ! Max. No. countries (with flat ! emissions) per grid
                     ,NEMIS_FILES & ! No. emission files
                     ,EMIS_NAME   & ! Names of species ("sox  ",...)
                     ,NBVOC       & ! > 0 if forest voc wanted
                     ,VOLCANOES   & ! 
                     ,ISNAP_SHIP  & ! snap index for ship emissions
                     ,ISNAP_NAT   & ! snap index for nat. (dms) emissions
                     ,VERTFAC     & ! vertical emission split
                     ,AIRNOX        !aircraft emissions or not
  use EmisGet_ml, only : EmisGet, EmisSplit, &
         nrcemis, nrcsplit, emisfrac &  ! speciation routines and array
        ,iqrc2itot                   &  !maps from split index to total index
        ,emis_masscorr               &  ! 1/molwt for most species
        ,emis_nsplit                    ! No. species per emis file
  use GridValues_ml, only:  GRIDWIDTH_M    & !  size of grid (m)
                           ,xm2            & ! map factor squared
                           ,debug_proc,debug_li,debug_lj & 
                           ,sigma_bnd, xmd, gl,dA,dB
  use Io_Nums_ml,      only : IO_LOG, IO_DMS, IO_EMIS
  use Io_Progs_ml,     only : ios, open_file
  use MetFields_ml,    only : roa, ps, z_bnd   ! ps in Pa, roa in kg/m3
  use ModelConstants_ml, only : KMAX_MID, KMAX_BND, PT ,dt_advec, &
                              IS_GLOBAL, & 
                              DEBUG => DEBUG_EMISSIONS,  MasterProc, & 
                              NPROC, IIFULLDOM,JJFULLDOM 
  use Par_ml,     only : MAXLIMAX,MAXLJMAX,me,gi0,gi1,gj0,gj1, &
                             GIMAX, GJMAX, IRUNBEG, JRUNBEG,  &   
                             limax,ljmax,li0,lj0,li1,lj1, &
                             MSG_READ1,MSG_READ7
  use PhysicalConstants_ml,  only :  GRAV,  AVOG
  use ReadField_ml, only : ReadField    ! Reads ascii fields
  use TimeDate_ml,only : nydays, nmdays, date, current_date   ! No. days per year, date-type 
  use Timefactors_ml, only : &
               NewDayFactors   &         ! subroutines
              ,timefac, day_factor  ! time-factors
  use Timefactors_ml, only : timefactors   &                 ! subroutine
                              ,fac_emm, fac_edd, day_factor  ! time-factors
  use Volcanos_ml


  implicit none
  private


 !/* subroutines:

  public :: Emissions                ! Main emissions module 
  public :: newmonth
  public :: EmisSet                  ! Sets emission rates every hour/time-step

  !/*   The main code does not need to know about the following  */
  private :: consistency_check       ! Safety-checks

  INCLUDE 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO


 !** land-code information in each grid square - needed to know which country
 !   is emitting.                        
 ! nlandcode = No. countries in grid square
 ! landcode  = Country codes for that grid square
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX)       :: nlandcode
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX,NCMAX) :: landcode
! for `flat emissions, i.e. no vertical extent:
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX)       :: flat_nlandcode
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX,FNCMAX):: flat_landcode

 !
 !  The output emission matrix for the 11-SNAP data is snapemis:
 !

  real, private, dimension(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES) &
            , save ::  snapemis !/* main emission arrays, in kg/m2/s

  real, private, dimension(MAXLIMAX,MAXLJMAX,FNCMAX,NEMIS_FILES) &
            , save ::  snapemis_flat !/* main emission arrays, in kg/m2/s  

 !ds May 2010 - we store the emissions for output to d_2d files and netcdf
 ! in kg/m2/s

  real, public, dimension(MAXLIMAX,MAXLJMAX,NEMIS_FILES), save :: SumSnapEmis

  !/-- emissions for input to chemistry routines

   ! KEMISTOP added to avoid hard-coded KMAX_MID-3:

   integer, public, parameter :: KEMISTOP = KMAX_MID - NEMISLAYERS + 1
   real, public, allocatable, save, dimension(:,:,:,:) :: &
        gridrcemis      & ! varies every time-step (as ps changes)
       ,gridrcemis0       ! varies every hour

  !/-- and for budgets (not yet used - not changed dimension)

   !DSRC real, public,  save, dimension(NRCEMIS) ::  totemadd
   real, public,  save, dimension(NSPEC_SHL+1:NSPEC_TOT) ::  totemadd



contains
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine Emissions(year)


 !+ calls main emission reading routines
 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !***********************************************************************
 !**    DESCRIPTION:
 !   0) call set_molwts and set_emisconv_and_iq, followed by
 !      consistency check
 !   1) Calls some setups:
 !         Country_Init
 !         timefactors: monthly and daily factors, +time zone
 !                            -> fac_emm, fac_edd arrays, timezone
 !   2) Read in emission correction file femis
 !   3) call emis_split for speciations
 !   4) Reads the annual emission totals in each grid-square, converts
 !      units and corrects using femis data. 
 !
 !  The output emission matrix for the 11-SNAP data is snapemis:
 !
 !       real    snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES)
 !  
 !**********************************************************************

  !--arguments
  integer, intent(in)   :: year        !  Year ( 4-digit)

  !-- local variables
  !DSRC integer, dimension(NEMIS_FILES) :: eindex   ! Index of emissions in EmisDef
  real    :: conv              ! Conversion factor
  integer :: iqrc, k, kused    ! index over emitted species, QRCSO2.. 
  integer :: i, j, n           ! Loop variables
  integer :: i_l,j_l           ! Local i,j
  real   :: tonne_to_kgm2s    ! Converts tonnes/grid to kg/m2/s
  real   :: ccsum             ! Sum of emissions for one country !ds, rv1_9_3
 
  ! arrays for whole EMEP area:
  !--    additional arrays on host only for landcode,nlandcode
  !..  BIG arrays ... will be needed only on me=0. Make allocatable
  ! to reduce static memory requirements.

  real,    allocatable, dimension(:,:,:,:)  :: globemis 
  integer, allocatable, dimension(:,:)      :: globnland 
  integer, allocatable, dimension(:,:,:)    :: globland  
  real,    allocatable, dimension(:,:,:)    :: globemis_flat
  integer, allocatable, dimension(:,:)      :: flat_globnland 
  integer, allocatable, dimension(:,:,:)    :: flat_globland 
  integer :: err1, err2, err3, err4, err5, err6! Error messages
  integer :: fic 
  integer :: ic        ! country codes 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over pollutants (1..NEMIS_FILES)
  integer :: eindex_vol       ! for volcanos


  !/** emission sums (after e_fact adjustments):
  real, dimension(NEMIS_FILES)       :: emsum    ! Sum emis over all countries
  real, dimension(NLAND,NEMIS_FILES) :: sumemis  ! Sum of emissions per country

  if (MasterProc) write(6,*) "Emissions called with me, year", me, year

  ! ** 0) set molwts, conversion factors (e.g. tonne NO2 -> tonne N), and
  !      emission indices (IQSO2=.., )

  !=========================
  !DSRC call EmisDef_Init()                      ! In EmisDef_ml
  call Country_Init()    ! In Country_ml, => NLAND, country codes and 
                         !                   names, timezone
  !=========================

  !DSRC do i = 1, NEMIS_FILES
  !DSRC    eindex(i) = EmisDef_Index( EMIS_NAME(i) )
  !DSRC end do

  !=========================
  !  Check that all is well!
    call consistency_check()               ! Below
    !DSRC call consistency_check(eindex)               ! Below
  !=========================
  ios = 0

  if( MasterProc) then   !::::::: ALL READ-INS DONE IN HOST PROCESSOR ::::

    ! ** 1)
    !=========================
     call timefactors(year)               ! => fac_emm, fac_edd, day_factor
    !=========================


    !** 2) 
    !=========================
    !DSRC call EmisSplit()    ! In EmisGet_ml, => emisfrac
    !=========================


  endif

  !DSRC - do read ins of splits on all procs
    !=========================
     call EmisSplit()    ! In EmisGet_ml, => emisfrac
    !=========================
  call CheckStop(ios, "ioserror: EmisSplit")


  ! #################################
  !    *** Broadcast  monthly and Daily factors ****
!DSRC    CALL MPI_BCAST( emisfrac ,8*NRCSPLIT*NSECTORS*NLAND,MPI_BYTE,  0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST( fac_emm ,8*NLAND*12*NSECTORS*NEMIS_FILES,MPI_BYTE,  0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST( fac_edd ,8*NLAND*7*NSECTORS*NEMIS_FILES,MPI_BYTE,   0,MPI_COMM_WORLD,INFO) 
    CALL MPI_BCAST( day_factor ,8*2*NSECTORS,MPI_BYTE,               0,MPI_COMM_WORLD,INFO) 

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !c4b) Set up DMS factors here - to be used in newmonth
  !    Taken from IQ_DMS=35 for SO2 nature (sector 11)
  !    first_dms_read is true until first call to newmonth finished.


  first_dms_read = .true. 
    
  !** 4) Read emission files ***************
  ! ******************************************************************

   ! allocate for me=0 only:
   err1 = 0
   if ( MasterProc ) then

       if (DEBUG) write(unit=6,fmt=*) "TTT me ", me , "pre-allocate" 
       allocate(globnland(GIMAX,GJMAX),stat=err1)
       allocate(globland(GIMAX,GJMAX,NCMAX),stat=err2)
       allocate(globemis(NSECTORS,GIMAX,GJMAX,NCMAX),stat=err3)
!
       allocate(flat_globnland(GIMAX,GJMAX),stat=err4)
       allocate(flat_globland(GIMAX,GJMAX,FNCMAX),stat=err5)
       allocate(globemis_flat(GIMAX,GJMAX,FNCMAX),stat=err6)

       call CheckStop(err1, "Allocation error 1 - globland")
       call CheckStop(err2, "Allocation error 2 - globland")
       call CheckStop(err3, "Allocation error 3 - globland")
       call CheckStop(err4, "Allocation error 4 - globland")
       call CheckStop(err5, "Allocation error 5 - globland")
       call CheckStop(err6, "Allocation error 6 - globland")


       !/**  and  initialise  **/
       globnland(:,:) = 0     ! csu  initialise globnland with 0
       flat_globnland(:,:)=0  !hf
       globland(:,:,:) = 0    !hk
       globemis(:,:,:,:) = 0  !hk
       flat_globland(:,:,:)=0 !hk
       globemis_flat(:,:,:) =0!hk

   end if

   do iem = 1, NEMIS_FILES
      ! now again test for me=0
      if ( MasterProc ) then

           ! read in global emissions for one pollutant
           ! *****************
             call EmisGet(iem,EMIS_NAME(iem),IRUNBEG,JRUNBEG,GIMAX,GJMAX, &
                          globemis,globnland,globland,sumemis,&
                          globemis_flat,flat_globnland,flat_globland)
           ! *****************


           emsum(iem) = sum( globemis(:,:,:,:) ) + &
                        sum( globemis_flat(:,:,:) )    ! hf
      endif  ! MasterProc

      call CheckStop(ios, "ios error: EmisGet")

      !**  Send data to processors ........
      !
      !     as  e.g. snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,iem)
      !
      !.................................
      !      send to nodes

       call global2local(globemis,snapemis(1,1,1,1,iem),MSG_READ1,   &
               NSECTORS,GIMAX,GJMAX,NCMAX,1,1)
!
       call global2local(globemis_flat,snapemis_flat(1,1,1,iem),MSG_READ1,   &
               1,GIMAX,GJMAX,FNCMAX,1,1)

    end do ! iem = 1, NEMIS_FILES-loop


    if ( MasterProc ) then
        write(unit=6,fmt=*) "Country totals"
        write(unit=IO_LOG,fmt=*) "Country totals"
        write(unit=6,fmt="(2a4,3x,30a12)")  "  N "," CC ",(EMIS_NAME(iem),iem=1,NEMIS_FILES)
        write(unit=IO_LOG,fmt="(2a4,3x,30a12)") "  N "," CC ",(EMIS_NAME(iem),iem=1,NEMIS_FILES)

        do ic = 1, NLAND
           ccsum = sum( sumemis(ic,:) )
           if ( ccsum > 0.0 ) then
                    write(unit=6,fmt="(i3,1x,a4,3x,30f12.2)") &
                        ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS_FILES)
                    write(unit=IO_LOG,fmt="(i3,1x,a4,3x,30f12.2)")& 
                        ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS_FILES)
           end if
        end do
    end if

    !    now all values are read, snapemis is distributed, globnland and 
    !    globland are ready for distribution
    !    print *, "calling glob2local_int for iem", iem, " me ", me

      call global2local_int(globnland,nlandcode,326, GIMAX,GJMAX,1,1,1)
      call global2local_int(globland, landcode ,326, GIMAX,GJMAX,NCMAX,1,1)
!
      call global2local_int(flat_globnland,flat_nlandcode,326,&
                            GIMAX,GJMAX,1,1,1)!extra array
      call global2local_int(flat_globland,flat_landcode,326,&
                            GIMAX,GJMAX,FNCMAX,1,1)

     !/**  broadcast volcanoe info derived in EmisGet 

        CALL MPI_BCAST(nvolc,4*1,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
        CALL MPI_BCAST(i_volc,4*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
        CALL MPI_BCAST(j_volc,4*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
        CALL MPI_BCAST(emis_volc,8*nvolc,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 



!    Conversions --
!
!     The emission-data file are so far in units of 
!     tonnes per grid-square. The conversion factor from tonnes per 50*50km2
!     annual emission values to surface flux (kg/m2/s) is found by division
!     with (nydays*24*60*60)s and (h*h)m2 and multiply by 1.e+3.
!     the conversion factor then equals 1.27e-14
!

    tonne_to_kgm2s  = 1.0e3 / (nydays * 24.0 * 3600.0 * GRIDWIDTH_M * GRIDWIDTH_M)

    if ( DEBUG .and.  MasterProc ) then
        write(unit=6,fmt=*) "CONV:me, nydays, gridwidth = ",me,nydays,GRIDWIDTH_M
        write(unit=6,fmt=*) "No. days in Emissions: ", nydays
        write(unit=6,fmt=*) "tonne_to_kgm2s in Emissions: ", tonne_to_kgm2s
        write(unit=6,fmt=*) "Emissions sums:"
        do iem = 1, NEMIS_FILES
           write(unit=6,fmt="(a15,f12.2)") EMIS_NAME(iem),emsum(iem)
        end do
    endif


    do iem = 1, NEMIS_FILES
       conv = tonne_to_kgm2s !DSRC * EmisDef( eindex(iem) )%conv
 
       forall (ic=1:NCMAX, j=lj0:lj1, i=li0:li1, isec=1:NSECTORS)
          snapemis (isec,i,j,ic,iem) = &
                 snapemis (isec,i,j,ic,iem) * conv * xm2(i,j)
       end forall

       forall (fic=1:FNCMAX, j=lj0:lj1, i=li0:li1)
          snapemis_flat(i,j,fic,iem) = &
                 snapemis_flat(i,j,fic,iem) * conv * xm2(i,j)
       end forall
    enddo !iem

    if ( VOLCANOES ) then

       !DSRC eindex_vol = EmisDef_Index( "sox" )
       conv = tonne_to_kgm2s !DSRC * EmisDef(eindex_vol)%conv

       do volc_no=1,nvolc 
          i=i_volc(volc_no)
          j=j_volc(volc_no)
          !Find global<->local coordinates for xm2
             if ((i >= gi0).and.(i<=gi1).and.(j>= gj0).and.&
                 (j<= gj1))then !on the correct processor
                if ( DEBUG ) write(*,*)'i,j for volcanoe is',i,j
                if ( DEBUG ) write(*,*)'EMIS_VOLC is',emis_volc(volc_no)
                i_l = i -gi0 +1
                j_l = j- gj0 +1
                 if ( DEBUG ) write(*,*)'Local coord is',i_l,j_l,gi0,gj0
                emis_volc(volc_no) = emis_volc(volc_no)* conv * xm2(i_l,j_l)
            endif
        enddo !volc_no

        !/** Read  Volcano.dat to get volcano height
        if (MasterProc)then
            call VolcGet(height_volc)
             if (DEBUG) write(*,*)'Volcano heights',height_volc
        endif

        !/** broadcast volcano heights
          CALL MPI_BCAST(height_volc,4*NMAX_VOLC,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
     endif ! VOLCANOES

    err1 = 0
    if ( MasterProc ) then
       deallocate(globnland,stat=err1)
       deallocate(globland ,stat=err2)
       deallocate(globemis ,stat=err3)

       deallocate(flat_globnland,stat=err4)
       deallocate(flat_globland,stat=err5)
       deallocate(globemis_flat,stat=err6)

       call CheckStop(err1, "De-Allocation error 1 - globland") 
       call CheckStop(err2, "De-Allocation error 2 - globland")
       call CheckStop(err3, "De-Allocation error 3 - globland")
       call CheckStop(err4, "De-Allocation error 4 - globland")
       call CheckStop(err5, "De-Allocation error 5 - globland")
       call CheckStop(err6, "De-Allocation error 6 - globland")

    end if

    !DSRC now we now nrecmis and can allocate for gridrcemis:
    !print *, "ALLOCATING GRIDRC", me, NRCEMIS
   allocate(gridrcemis(NRCEMIS,KEMISTOP:KMAX_MID,MAXLIMAX,MAXLJMAX),stat=err1)
   allocate(gridrcemis0(NRCEMIS,KEMISTOP:KMAX_MID,MAXLIMAX,MAXLJMAX),stat=err2)
   call CheckStop(err1, "Allocation error 1 - gridrcemis") 
   call CheckStop(err2, "Allocation error 2 - gridrcemis0")

  end subroutine Emissions

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine consistency_check() !DSRC eindex)
  !------------------------------------------------------------------!
  !    checks that all the values given so far are consistent        !
  !------------------------------------------------------------------!
!DSRC  integer, dimension(NEMIS_FILES), intent(in) :: eindex
! Should add more checks
  character(len=30) :: errormsg
  integer :: i

  errormsg = "ok"
  if ( size(EMIS_NAME) /= NEMIS_FILES    ) errormsg = " size EMISNAME wrong "

  call CheckStop(errormsg,"Failed consistency check")

 end subroutine consistency_check
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  !*****************************************************************

  subroutine EmisSet(indate)   !  emission re-set every time-step/hour

  !
  !***********************************************************************
  !**    DESCRIPTION:
  !     calculates the emmision-tendencies and the local (instantaneous) dry 
  !     deposition in the emission squares.
  !       emis set once per hour to allow for day/night variation (and voc 
  !       speciation) (based on local time)  for each snap sector.
  !     gridrcemis0 calculated every time-step to allow for ps changes.
  !     inputs from Emissions in EMISSIONS_ML:
  !      country and snap-specific array : 
  !          snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES) 
  !  
  !*** Units:
  !     snapemis has units of kg/m2/s, SO2 as S, NO2 as N, NH3 as N. 
  !     Map factor (xm2) already accounted for. 
  !  
  !    Data on how many countries contribute per grid square is stored in
  !    nlandcode(MAXLIMAX,MAXLJMAX) , with the country-codes given by
  !    landcode(MAXLIMAX,MAXLJMAX,NCMAX).
  !     
  !    Monthly and weekday factors are pre-multiplied and stored in:
  !       real timefac(NLAND,NSECTORS,NEMIS_FILES)
  !    And day-night factors are applied here:
  !       day_factor(11,0:1)                  ! 0=night, 1=day
  ! ..........................................................................
  !
  !**    REVISION HISTORY:
  !       Revised , 30/5/01, jej/st found problem on gridur - split NEMIS_FILES loop 
  !       into separate NEMIS_PLAIN and NEMIS_SPLIT loops.
  !       Revised, ds,  Feb. 2001 for unified model. Use of date%seconds replaces
  !       thourloc.
  !       Revised  : d. simpson 4/2/98 to act as common subrouinte
  !       between 3-D models, and to avoid hard-coded emissions
  !       !uni - revised to F90 and for more flexible handling of emissions
  !       fraction through NEMIS_FRAC
  !
  !      25/3-2002, pw changed test for hour and day change (Now the first day
  !      does not need to start at 0 hours)
  !
  !      11/2005,pw if timezone=-100 use timezone based on longitude
  !
  !*************************************************************************

  implicit none
  type(date), intent(in) :: indate           ! Gives year..seconds
  integer :: i, j, n, k, f                   ! cooridnates, loop variables
  integer :: icc, ncc                        ! No. of countries in grid.
!
  integer :: ficc,fncc                       ! No. of countries with
  integer :: iqrc                            ! emis indices 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILES
  integer :: itot             ! index in xn()
!
  !uni - save daytime value between calls, intiialise to zero
  integer, save, dimension(NLAND) ::  daytime = 0  !  0=night, 1=day
  integer                         ::  hourloc      !  local hour 
  logical                         ::  hourchange   !      "     "           
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions

  real ::  ehlpcom,ehlpcom0(KEMISTOP:KMAX_MID)
  real ::  tfac, dtgrid    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s               ! source term (emis) before splitting
  integer :: iland, iland_timefac  ! country codes, and codes for timefac 

  real ::  ftfac           ! time-factor for flat emissions
  real ::  sf              ! source term (emis) before splitting (for flat emissions)
  integer :: flat_iland    ! country codes (countries with flat emissions)

  integer, save :: oldday = -1, oldhour = -1

! If timezone=-100, calculate daytime based on longitude rather than timezone
  integer :: daytime_longitude,daytime_iland
 
! Initialize
    ehlpcom0(:)=0.0

   do k=KEMISTOP,KMAX_MID
      ehlpcom0(k) = GRAV* 0.001*AVOG/ (sigma_bnd(k+1) - sigma_bnd(k))
   enddo

   !/** scaling for totemadd:
     dtgrid = dt_advec * GRIDWIDTH_M * GRIDWIDTH_M 


   !/** The emis array only needs to be updated either every full hour. The 
   !    time-factor calculation needs to know if a local-time causes a shift 
   !    from day to night.  In addition, we reset an overall day's Time-factors
   !    at midnight every day. 

     hourchange = .false.
     if ( indate%hour /= oldhour .or. indate%day /= oldday ) then

           hourchange = .true.
           oldhour = indate%hour

           if ( indate%day /= oldday  )then

              !==========================
               call NewDayFactors(indate)
              !==========================

               oldday = indate%day
           endif
     end if


    !..........................................
    !  Look for day-night changes, after local time correction
    !  (daytime(iland) not used if  LONGITUDE_TIME=true)

    do iland = 1, NLAND

       daytime(iland) = 0
       hourloc        = indate%hour + Country(iland)%timezone

       if ( hourloc  >=   7 .and.  hourloc <= 18 ) daytime(iland)=1

    end do ! iland


    if ( hourchange ) then 

         totemadd(:)  = 0.
         gridrcemis0(:,:,:,:) = 0.0 
         SumSnapEmis(:,:,:) = 0.0


        !..........................................
        ! Process each grid:

        do j = lj0,lj1
            do i = li0,li1

               ncc = nlandcode(i,j)            ! No. of countries in grid

                ! find the approximate local time:
                  hourloc= mod(nint(indate%hour+24*(1+gl(i,j)/360.0)),24)
                  daytime_longitude=0
                  if( hourloc>=7.and.hourloc<= 18) daytime_longitude=1
    
         
              !*************************************************
              ! First loop over non-flat(one sector) emissions
              !*************************************************

              tmpemis(:)=0.
              do icc = 1, ncc
                  iland = landcode(i,j,icc)     ! 1=Albania, etc.
                  iland_timefac = Country(iland)%timefac_index

                if(Country(iland)%timezone==-100)then
                   daytime_iland=daytime_longitude
                else
                   daytime_iland=daytime(iland)
                endif

                 !  As each emission sector has a different diurnal profile
                 !  and possibly speciation, we loop over each sector, adding
                 !  the found emission rates to gridrcemis as we go.
                 !  ==================================================


                do isec = 1, NSECTORS       ! Loop over snap codes

                   !--  Calculate emission rates from snapemis, time-factors, 
                   !    and if appropriate any speciation  fraction (NEMIS_FRAC)
                   iqrc = 0   ! index over emisfrac

                   do iem = 1, NEMIS_FILES 

                      tfac = timefac(iland_timefac,isec,iem) * &
                                 day_factor(isec,daytime_iland)

                      s =  tfac * snapemis(isec,i,j,icc,iem)

            !prelim emis sum kg/m2/s
                       SumSnapEmis(i,j,iem) = SumSnapEmis(i,j,iem) + s

                      do f = 1, emis_nsplit( iem )
                           iqrc = iqrc + 1
                           itot = iqrc2itot(iqrc)
                           tmpemis(iqrc) = s * emisfrac(iqrc,isec,iland)
                        !--  Add up emissions in ktonne ......
                           totemadd(itot) = totemadd(itot) + &
                                     tmpemis(iqrc) * dtgrid * xmd(i,j)
                      end do ! f

                   end do ! iem


                   !..   Assign to height levels 1-KEMISTOP

                   do k=KEMISTOP,KMAX_MID
                      do iqrc =1, nrcemis
                         gridrcemis0(iqrc,k,i,j) =   &
                            gridrcemis0(iqrc,k,i,j) + tmpemis(iqrc)*   &
                            ehlpcom0(k)*VERTFAC(KMAX_BND-k,isec) &
                            * emis_masscorr(iqrc)   !DSRC was /molwt
                      end do ! iem
                   end do   ! k

                enddo  ! isec
!      ==================================================

       end do ! icc  
 
       !************************************
       ! Then loop over flat emissions
       !************************************
       tmpemis(:)=0.
       fncc = flat_nlandcode(i,j) ! No. of countries with flat 
                                          ! emissions in grid

       do ficc = 1, fncc
          flat_iland = flat_landcode(i,j,ficc) ! 30=BAS etc.

          if ( Country(flat_iland)%is_sea ) then   ! - saves if statements below
               isec = ISNAP_SHIP 
          else
               isec = ISNAP_NAT
          end if

          !  As each emission sector has a different diurnal profile
          !  and possibly speciation, we loop over each sector, adding
          !  the found emission rates to gridrcemis as we go.
          !  ==================================================


          !--  Calculate emission rates from snapemis, time-factors, 
          !    and if appropriate any speciation  fraction (NEMIS_FRAC)

            iqrc  = 0   ! index over emis

            do iem = 1, NEMIS_FILES 

                sf =  snapemis_flat(i,j,ficc,iem)    

            !prelim emis sum kg/m2/s
                SumSnapEmis(i,j,iem) = SumSnapEmis(i,j,iem) + sf

                do f = 1, emis_nsplit( iem )
                   iqrc = iqrc + 1
                   itot = iqrc2itot(iqrc)
                   tmpemis(iqrc) = sf * emisfrac(iqrc,isec,flat_iland)

                  !--   Add flat emissions in ktonne ......
                   totemadd(itot) = totemadd(itot) + &
                             tmpemis(iqrc) * dtgrid * xmd(i,j)

                end do ! f

            end do ! iem

         !..   Assign flat emissions to height levels 1-4
         !..   Note, no VERTFAC

             do iqrc =1, nrcemis

                gridrcemis0(iqrc,KMAX_MID,i,j) =   &
                  gridrcemis0(iqrc,KMAX_MID,i,j) + tmpemis(iqrc)*&
                    ehlpcom0(KMAX_MID) * emis_masscorr(iqrc)   !DSRC was /molwt
             end do ! iem

!      ==================================================
       end do !ficc 
   end do ! i
 end do ! j

    if ( VOLCANOES ) call Set_Volc !set hourly volcano emission(rcemis_volc0)

  end if ! hourchange 


  !/** we now scale gridrcemis to get emissions in molecules/cm3/s

   do k= KEMISTOP, KMAX_MID
     do j = lj0,lj1
        do i = li0,li1

              ehlpcom= roa(i,j,k,1)/(ps(i,j,1)-PT)
              do iqrc =1, NRCEMIS
                 gridrcemis(iqrc,k,i,j) =  &
                    gridrcemis0(iqrc,k,i,j)* ehlpcom
              enddo  ! iqrc
        end do ! i
      end do ! j
   end do ! k

 !/** Scale volc emissions to get emissions in molecules/cm3/s (rcemis_volc)
   if ( VOLCANOES ) call Scale_Volc

 if( NBVOC > 0  )then

    do j = lj0,lj1
      do i = li0,li1

        ehlpcom = ehlpcom0(KMAX_MID) * roa(i,j,KMAX_MID,1)/(ps(i,j,1)-PT)

        emnat(i,j,1:NBVOC) = emforest(i,j,1:NBVOC)*ehlpcom

      end do ! i
    end do ! j

    if ( DEBUG .and. debug_proc ) then 
       !print "(a12,2i4,/,(g12.3))",  "bio-setemis",li0,lj0, &
        !  (gridrcemis(i,KMAX_MID,2,2),i=1, NRCEMIS), 
       write(*,"(a12,2g12.3,3x,2g12.3)")  "bio-setemis", &
         ( emforest(debug_li,debug_li,i), i = 1, NBVOC),&
         ( emnat(debug_li,debug_li,i), i = 1, NBVOC) 

      ! dz = dp/(rho.g) = dsigma*pstar/(rho.g)
      write(*,"(a,2f10.4)") "DEBUG BIODZ ", &
             (sigma_bnd(21) - sigma_bnd(20)) * (ps(i,j,1)-PT)/ &
                ( roa(i,j,KMAX_MID,1) * GRAV ) , z_bnd(i,j,20)
    end if

  endif ! NBVOC 

 end subroutine EmisSet
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    subroutine newmonth

!.....................................................................
!**    DESCRIPTION:
!     Reads in natural DMS emissions at start of each month. Update
!     landcode and nlandcode arrays as needed.

!     Reads in snow cover at start of each month. 
!
!     April 2010: read monthly aircraft NOx emissions
!
!...........................................................................
      use AirEmis_ml , only : airn
      use ModelConstants_ml    , only : KCHEMTOP, KMAX_MID
      use NetCDF_ml, only : ReadField_CDF

    integer i, j,k
    integer ijin(2) 
    integer n, flat_ncmaxfound      ! Max. no. countries w/flat emissions
    real :: rdemis(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
    character*20 fname
    real ktonne_to_kgm2s, tonnemonth_to_kgm2s          ! Units conversion
    integer :: IQSO2             ! Index of sox in  EMIS_NAME
    integer errcode
    real,    allocatable, dimension(:,:,:,:)  :: globemis 
!DSRC    integer, dimension(NEMIS_FILES) :: eindex   ! Index of emissions in EmisDef
    integer:: month,iem,ic,iic,isec, err3,ic1,icc
    real ::duml,dumh,tmpsec(NSECTORS),conv
        logical ,save ::first_call=.true.
    real, dimension(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS_FILES) &
            ::  snapemis_month !/*monthly emissions tonne/month

   !dsx For now, only the global runs use the Monthly files.
   !   Will need to reconsider later.
  !dsx logical, parameter ::MONTHLY_GRIDEMIS=.true. ! GLOBAL

        logical, parameter ::MONTHLY_GRIDEMIS= IS_GLOBAL      
        integer ::kstart,kend

if(AIRNOX)then
!AIRCRAFT
kstart=KCHEMTOP
kend=KMAX_MID
do k=KEMISTOP,KMAX_MID
do j=1,ljmax
do i=1,limax

airn(k,i,j)=0.0

enddo
enddo
enddo
call ReadField_CDF('AircraftEmis_FL.nc','NOx',airn,nstart=current_date%month,kstart=kstart,kend=kend,interpol='mass_conservative', &
     needed=.true.,debug_flag=.true.)

!convert from kg(NO2)/month into molecules/cm3/s
!from kg to molecules: 0.001*AVOG/species(NO2)%molwt
!use roa to find dz for consistency with other emissions (otherwise could have used z_bnd directly)
!dz=dP/(roa*GRAV)  dP=dA(k) + dB(k)*ps(i,j,1)
!dV=dz*dx*dy=dz*gridwidth**2/xm**2 *1e6 (1e6 for m3->cm3)
!from month to seconds: ndaysmonth*24*3600

conv=0.001*AVOG/species(NO2)%molwt*GRAV/gridwidth_m**2*1.0e-6/(nmdays(current_date%month)*24*3600)
do k=KEMISTOP,KMAX_MID
do j=1,ljmax
do i=1,limax

airn(k,i,j)=airn(k,i,j)*conv*(roa(i,j,k,1))/(dA(k) + dB(k)*ps(i,j,1))*xm2(i,j)

enddo
enddo
enddo

endif


!DMS
!*** Units:
!    Input files seem to be in ktonne PER YEAR. We convert here to kg/m2/s
!     to save CPU in setemis.f.
!      The conversion factor from 50*50km2
!      annual emission values to surface flux (kg/m2/s) is found by division
!      with (nydays*24*60*60)s and (h*h)m2 and multiply by 1.e+6.
!      the conversion factor (ktonne_to_kgm2s) then equals 1.27e-8 
!      NB: a new file is read every month; this means that total emissions 
!          are NOT the sum of the 12 files emissions (but about 12 times less than the sum). 
!          More precisely: year_emis=sum_months(emis_month*nmdays/nydays)

        ktonne_to_kgm2s  = 1.0e6 /        &
        (nydays*24.*60.*60.*GRIDWIDTH_M*GRIDWIDTH_M)

    if ( MasterProc .and. DEBUG) then
      write(6,*) 'Enters newmonth, mm, ktonne_to_kgm2s = ',    &
          current_date%month,ktonne_to_kgm2s
      write(6,*) ' first_dms_read = ', first_dms_read
    end if ! me 
!...........................................................................

!...........................................................................
!        DMS Input - land 35 - SNAP sector 11
!...........................................................................
     flat_ncmaxfound = 0 ! Max. no. countries(w/flat emissions) per grid
!    natural so2 emissions

          IQSO2 = 0
          do i = 1, NEMIS_FILES
            if ( trim( EMIS_NAME(i) ) == "sox" ) IQSO2 = i
          end do

          if ( IQSO2 < 1 ) then
              write(*,*) " No SO2 emissions - need to skip DMS also"
              return    ! No need to read DMS fields ...

          else    
            !/--- we have so2 emission so need DMS also

             if ( MasterProc ) then

              write(fname,fmt='(''natso2'',i2.2,''.dat'')')     &
                current_date%month
              write(6,*) 'filename for nat-so2',fname
              endif

              call ReadField(IO_DMS,fname,rdemis)

          errcode = 0
          do j=1,ljmax
              do i=1,limax

!            Add DMS for country code IQ_DMS=35  to snap sector 11=Nature.
                ! First time we read we must add DMS to the "countries" 
                ! contributing within the grid square. 
 
                  ! - for flat emissions:

                  if ( first_dms_read ) then 
                     flat_nlandcode(i,j) = flat_nlandcode(i,j) + 1 
                     n = flat_nlandcode(i,j)
                     flat_landcode(i,j,n) = IQ_DMS   ! country code 35 
                     if ( n > flat_ncmaxfound ) then
                          flat_ncmaxfound = n 
                          if (DEBUG) write(6,*)'DMS Increased flat_ncmaxfound to ',n 
                          call CheckStop( n > FNCMAX, "IncreaseFNCMAX for dms")
                     endif 
                  else  ! We know that DMS lies last in the array, so:
                      n = flat_nlandcode(i,j)
                      call CheckStop(flat_landcode(i,j,n), IQ_DMS, &
                          "Newmonth:DMS not last!")
                  endif 

                  snapemis_flat(i,j,n,IQSO2) = rdemis(i,j) * ktonne_to_kgm2s &
                     * xm2(i,j)
            enddo ! i
          enddo ! j


              if ( first_dms_read ) then
                 if (DEBUG) write(6,*)'me ',me, ' Increased flat_ncmaxfound to ' &
                  ,flat_ncmaxfound 
              first_dms_read = .false.
          end if

           end if  ! IQSO2>0


           if(MONTHLY_GRIDEMIS)then

   !Read monthly emission files

   if(first_call)then
      do j=lj0,lj1
         do i=li0,li1
            nlandcode(i,j)=nlandcode(i,j)+1
            icc=nlandcode(i,j)
            landcode(i,j,icc)=67 
         enddo
      enddo
      first_call=.false.
   endif
    month = current_date%month
    tonnemonth_to_kgm2s= 1.0e3 /         &
         (nmdays(month)*24.*60.*60.*GRIDWIDTH_M*GRIDWIDTH_M)

    if ( MasterProc ) then
       allocate(globemis(NSECTORS,GIMAX,GJMAX,NCMAX),stat=err3)
       call CheckStop(err3, "Allocation error err3 - globland")
    end if
!DSRC    do i = 1, NEMIS_FILES
!DSRC       eindex(i) = EmisDef_Index( EMIS_NAME(i) )
!DSRC    end do
    do iem = 1, NEMIS_FILES
!      if (trim(EMIS_NAME(iem)).ne.'nox' .and. trim(EMIS_NAME(iem)).ne.'co'.and.&
!           trim(EMIS_NAME(iem)).ne.'pm25'.and.&
!           trim(EMIS_NAME(iem)).ne.'voc'.and.trim(EMIS_NAME(iem)).ne.'nh3'.and.trim(EMIS_NAME(iem)).ne.'sox')cycle !
!      snapemis (:,:,:,:,iem) = 0.0 !NB all previous (snap)emis set to 0
      if ( MasterProc ) then

         globemis = 0.0

         write(fname,fmt='(''grid'',A,i2.2)')     &
              trim(EMIS_NAME(iem))//'.',month
         write(6,*) 'filename for GLOBAL emission',fname
         call open_file(IO_EMIS,"r",fname,needed=.true.)
         call CheckStop( ios , "ios error: emislist" // fname )

         READEMIS: do   ! ************* Loop over emislist files **********************

            read(unit=IO_EMIS,fmt=*,iostat=ios) iic,i,j, duml,dumh,  &
                 (tmpsec(isec),isec=1,NSECTORS)
            if ( ios <  0 ) exit READEMIS            ! End of file
            call CheckStop( ios , "GetEmis ios: error on " // fname ) ! exits if ios>0

            ic=1 !NBNB default country



            i = i-IRUNBEG+1     ! for RESTRICTED domain
            j = j-JRUNBEG+1     ! for RESTRICTED domain

            if ( i  <=  0 .or. i  >  GIMAX .or.   & 
                 j  <=  0 .or. j  >  GJMAX .or.   &
                 ic <=  0 .or. ic >  NLAND .or.   &
                 ic == IC_NAT )                   &  ! Excludes DMS
                 cycle READEMIS
            globemis(:,i,j,ic) = globemis(:,i,j,ic) & !NB: put everything in land "1"
                 !NB no femis                  + e_fact(:,ic,iemis) *  tmpsec(:)!e_fact=femis
                 + tmpsec(:)

         end do READEMIS
         !
         close(IO_EMIS)
         ios = 0
!         globemis(:,i,j,ic1) = 0.0

     endif  ! MasterProc
      call global2local(globemis,snapemis_month(1,1,1,1,iem),MSG_READ1,   &
           NSECTORS,GIMAX,GJMAX,NCMAX,1,1)
   end do ! iem = 1, NEMIS-loop
!   nlandcode=1
!   landcode=67 !default land
   ic=1 !BUG corrected in this directory 3/11-2008 (ic=1 only for ME=0 otherwise)
   do iem = 1, NEMIS_FILES
!         write(*,*)'iem=',iem
!      if (trim(EMIS_NAME(iem)).ne.'nox' .and. trim(EMIS_NAME(iem)).ne.'co'.and.&
!           trim(EMIS_NAME(iem)).ne.'pm25'.and.&
!           trim(EMIS_NAME(iem)).ne.'voc'.and.trim(EMIS_NAME(iem)).ne.'nh3'.and.trim(EMIS_NAME(iem)).ne.'sox')cycle !
      conv = tonnemonth_to_kgm2s !DSRC * EmisDef( eindex(iem) )%conv
      do j=lj0,lj1
         do i=li0,li1
            icc=nlandcode(i,j)!67
           do isec=1,NSECTORS
              snapemis (isec,i,j,icc,iem) = &
                    snapemis_month (isec,i,j,ic,iem) * conv * xm2(i,j)
            enddo
         enddo
     enddo
   enddo !iem
   if ( MasterProc ) then
      deallocate(globemis)
   end if
   endif

    end subroutine newmonth

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module Emissions_ml
