!_____________________________________________________________________________
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD

                    module Emissions_ml

! MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD MOD  MOD MOD MOD MOD MOD MOD MOD
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! October 2001 hf added call to ReadField
! 6b: Nov. 2001 - dead code removed, some renaming (rsnd -> rdemis), ds
! Modified for Joffen's new mass budget method - Bud_ml removed here.
! HF Mars-02 Flat emission traetment and volcanos added
! ToDo: remove DMS treatment, add point sources, volcanoes
! 6d: correction to tfac spotted by hf, 29/11/01
!_____________________________________________________________________________
!+ Calls up emission read/set routines
!  This routine interfaces the stand-alone emission-file reading routines
!  with the 3D model.
!_____________________________________________________________________________
  use My_Emis_ml,  only : &
          NEMIS,        &   ! No. emission files
          NEMIS_PLAIN,  &   ! No. emission files for non-speciated emissions
          NEMIS_SPLIT,  &   ! No. emission files for speciated emissions
          NRCEMIS,&         ! Total No. emission species after speciation
          VOLCANOES         ! u3 - changed from VOL
  use My_MassBudget_ml, only : set_mass_eqvs   ! Some equivalences bewteen 
                                               ! indices

  use EmisDef_ml, only : NSECTORS,  &  ! No. sectors
                         NEMISLAYERS,& !rv1.3b No. vertical layers for emission
                         NCMAX,&       ! Max. No. countries per grid
                         FNCMAX,&      ! Max. No. countries (with flat 
                                       ! emissions) per grid
                         ISNAP_SHIP,&  ! snap index for ship emissions
                         ISNAP_NAT     ! snap index for nat. (dms) emissions
  use Io_ml,      only : IO_LOG
  use Par_ml,     only : MAXLIMAX,MAXLJMAX,NPROC,me,gi0,gi1,gj0,gj1

  use ModelConstants_ml,     only :   KMAX_MID
  use Volcanos_ml                   !hf

  implicit none
  private


 !/* subroutines:

  public :: Emissions                ! Main emissions module 
  public :: newmonth
  public :: EmisSet                  ! Sets emission rates every hour/time-step

  !/*   The main code does not need to know about the following  */
  private :: consistency_check       ! Safety-checks

  logical, private, parameter :: DEBUG = .false.

 !** land-code information in each grid square - needed to know which country
 !   is emitting.                        
 ! nlandcode = No. countries in grid square
 ! landcode  = Country codes for that grid square
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX)       :: nlandcode
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX,NCMAX) :: landcode
!hf
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX)       :: flat_nlandcode
   integer, private, save, dimension(MAXLIMAX,MAXLJMAX,FNCMAX):: flat_landcode
 !
 !  The output emission matrix for the 11-SNAP data is snapemis:
 !

  real, private, dimension(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS) &
            , save ::  snapemis !/* main emission arrays, in kg/m2/s

!hf
  real, private, dimension(MAXLIMAX,MAXLJMAX,FNCMAX,NEMIS) &
            , save ::  snapemis_flat !/* main emission arrays, in kg/m2/s  

  !/-- emissions for input to chemistry routines

   ! KEMISTOP added to avoid hard-coded KMAX_MID-3:
   !rv1.2 now calculated from new NEMISLAYERS param

   integer, public, parameter :: KEMISTOP = KMAX_MID - NEMISLAYERS + 1
   real, public, save, dimension(NRCEMIS,KEMISTOP:KMAX_MID,MAXLIMAX,MAXLJMAX) :: &
         gridrcemis         ! varies every time-step (as ps changes)
   real, private, save, dimension(NRCEMIS,KEMISTOP:KMAX_MID,MAXLIMAX,MAXLJMAX) :: &
          gridrcemis0       ! varies every hour

  !/-- and for budgets (not yet used - not changed dimension)

   real, public,  save, dimension(NRCEMIS) ::  totemadd



contains


 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine Emissions(year)

  use My_Emis_ml,  only : &
          EMIS_NAME,    &   ! Names of species ("sox  ",...)
          NEMIS_SPLIT,  &   ! No. emission files to be speciated
          NRCSPLIT,     &   ! No. emis species from split emissions species
          QRCVOL,       &   ! For volcanoes
          set_molwts      ! subroutine to set molwt
  use Biogenics_ml, only:               &
          Forests_Init                  &
          , first_dms_read 
  use Country_ml,    only : NLAND,Country_Init,Country
  use Dates_ml,      only : nydays     ! No. days per year, date-type 
  use EmisDef_ml,    only : &  
                     EmisDef_Init    &! Sub to define conversion factors
                     ,EmisDef_Index   &! Sub to get index of emis name 
                     , EmisDef        ! Superset of names/factors
  use EmisGet_ml, only : EmisGet, EmisSplit, emisfrac  ! speciation routines and array

  use Timefactors_ml, only : timefactors   &                 ! subroutine
                              ,fac_emm, fac_edd, day_factor  ! time-factors

  use Par_ml,         only : GIMAX, GJMAX, ISMBEG, JSMBEG,  &   
                              li0,lj0,li1,lj1, MSG_READ1
  use GridValues_ml,  only : xm2,GRIDWIDTH_M
  use Io_ml,          only : ios  

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
 !       real    snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS)
 !  
 !**    REVISION HISTORY:
 !         Original from MADE, adapted to MACHO j.e.jonson
 !         Call to femis added, j.e. jonson
 !         1-2/99 Rearranged d. simpson, for 11-sector input. Condensed. 
 !         Monthly and daily input read in in subroutine timefactors:
 !         efac_mm (monthly) and efac_dd (daily)
 !         added . Note  that emissions are no longer multiplied by
 !         monthly factors here - this is done in subroutine emission where
 !         also daily and hourly factors are applied.
 !         4/2/99 - checked/corrected and re-arranged by s. unger.
 !         8/2/99 - timezone and eulxxxx.inc method added by ds
 !
 !**********************************************************************

  !--arguments
  integer, intent(in)   :: year        !  Year ( 4-digit)

  !-- local variables
  integer, dimension(NEMIS) :: eindex   ! Index of emissions in EmisDef
  real    :: conv              ! Conversion factor
  integer :: iqrc, k, kused    ! index over emitted species, QRCSO2.. 
  integer :: i, j, n           ! Loop variables
!hf
  integer :: i_l,j_l           ! Local i,j
  real   :: tonne_to_kgm2s    ! Converts tonnes/grid to kg/m2/s
  real   :: ccsum             ! Sum of emissions for one country !ds, rv1_9_3
 
  !/ emissions for loal area and 1 pollutant. Temporary array to
  !  help tranfer through global2local subroutine:
!su  real, dimension(NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX) :: localemis

  ! arrays for whole EMEP area:
  !su--    additional arrays on host only for landcode,nlandcode
  !.. ds BIG arrays ... will be needed only on me=0. Make allocatable
  ! to reduce static memory requirements.

  real,    allocatable, dimension(:,:,:,:)  :: globemis 
  integer, allocatable, dimension(:,:)      :: globnland 
  integer, allocatable, dimension(:,:,:)    :: globland  
!hf
  real,    allocatable, dimension(:,:,:)    :: globemis_flat
  integer, allocatable, dimension(:,:)      :: flat_globnland 
  integer, allocatable, dimension(:,:,:)    :: flat_globland 
  integer :: err1, err2, err3, err4, err5, err6! Error messages
  integer :: fic 
  integer :: ic        ! country codes 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over pollutants (1..NEMIS)
  integer :: info             ! Used in gc_rsend


  !/** emission sums (after e_fact adjustments):
  real, dimension(NEMIS)       :: emsum    ! Sum emis over all countries
  real, dimension(NLAND,NEMIS) :: sumemis  ! Sum of emissions per country

  write(6,*) "Emissions called with me, year", me, year

  ! ** 0) set molwts, conversion factors (e.g. tonne NO2 -> tonne N), and
  !      emission indices (IQSO2=.., )

  !=========================
  call EmisDef_Init()                      ! In EmisDef_ml
  call set_molwts()                        ! In My_Emis_ml
  call set_mass_eqvs()                     ! In My_MassBudget_ml
  call Country_Init()    ! In Country_ml, => NLAND, country codes and 
                         !                   names, timezone
  !=========================

  do i = 1, NEMIS
     eindex(i) = EmisDef_Index( EMIS_NAME(i) )
  end do
  !u4 - error checks moved to consistency_check

  !=========================
  !  Check that all is well!
    call consistency_check(eindex)               ! Below
  !=========================
  ios = 0

  if( me ==  0) then   !::::::: ALL READ-INS DONE IN HOST PROCESSOR ::::

    ! ** 1)
    !=========================
     call timefactors(year)               ! => fac_emm, fac_edd, day_factor
    !=========================

      if ( ios /= 0) call gc_abort(me,NPROC,"ios error: timefactors")

    !** 2) 
    !=========================
     if ( NEMIS_SPLIT > 0 ) call EmisSplit()    ! In EmisGet_ml, => emisfrac
    !=========================


  endif !(me=0)
  if ( ios /= 0 ) call gc_abort(me,NPROC,"ios error: EmisSplit")

  ! #################################
  !uni    *** Broadcast  monthly and Daily factors ****
  call gc_rbcast(601, NRCSPLIT*NSECTORS*NLAND,  0, NPROC, info, emisfrac)
  call gc_rbcast(701, NLAND*12*NSECTORS*NEMIS,  0, NPROC, info, fac_emm)
  call gc_rbcast(702, NLAND*7*NSECTORS*NEMIS,   0, NPROC, info, fac_edd)
  call gc_rbcast(703, 2*NSECTORS,               0, NPROC, info, day_factor)

  !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !c4b) Set up DMS factors here - to be used in newmonth
  !    Taken from IQ_DMS=35 for SO2 nature (sector 11)
  !    first_dms_read is true until first call to newmonth finished.

  !gv print *, "TTT, after emisconv_and_iq, IQSO2 = ", IQSO2, emisconv(IQSO2)
  !gv dms_fact = e_fact(11,IQ_DMS,IQSO2)

  !u4 dms_fact = 1.0            !! ** ds - temporary solution !!!!!
  first_dms_read = .true. 
    
  !** 4) Read emission files ***************
  ! ******************************************************************

   !uni allocate for me=0 only:
   err1 = 0
   if ( me == 0 ) then

       write(unit=6,fmt=*) "TTT me ", me , "pre-allocate" 
       allocate(globnland(GIMAX,GJMAX),stat=err1)
       allocate(globland(GIMAX,GJMAX,NCMAX),stat=err2)
       allocate(globemis(NSECTORS,GIMAX,GJMAX,NCMAX),stat=err3)
!hf
       allocate(flat_globnland(GIMAX,GJMAX),stat=err4)
       allocate(flat_globland(GIMAX,GJMAX,FNCMAX),stat=err5)
       allocate(globemis_flat(GIMAX,GJMAX,FNCMAX),stat=err6)
       if ( err1 /= 0 .or. err2 /= 0 .or. err3 /= 0 .or. &
            err4 /= 0 .or. err5 /= 0 .or. err6 /= 0) then
            call gc_abort(me,NPROC,"Alloc error - globland")
       end if


       !/**  and  initialise  **/
       globnland(:,:) = 0     ! csu  initialise globnland with 0
       flat_globnland(:,:)=0  !hf
   end if

   do iem = 1, NEMIS
      !su    now again test for me=0
      if ( me == 0 ) then

           ! read in global emissions for one pollutant
           ! *****************
!hf             call EmisGet(iem,EMIS_NAME(iem),ISMBEG,JSMBEG,GIMAX,GJMAX, &
!hf                          globemis,globnland,globland,sumemis)
             call EmisGet(iem,EMIS_NAME(iem),ISMBEG,JSMBEG,GIMAX,GJMAX, &
                          globemis,globnland,globland,sumemis,&
                          globemis_flat,flat_globnland,flat_globland)
           ! *****************


           emsum(iem) = sum( globemis(:,:,:,:) ) + &
                        sum( globemis_flat(:,:,:) )    ! hf
      endif  ! me==0
      if ( ios /= 0) call gc_abort(me,NPROC,"ios error: EmisGet")

      !CC**  Send data to processors ........
      !
      !     as  e.g. snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,iem)
      !
      !.................................
      !      send to nodes

       call global2local(globemis,snapemis(1,1,1,1,iem),MSG_READ1,   &
               NSECTORS,GIMAX,GJMAX,NCMAX,1,1)
!hf F
       call global2local(globemis_flat,snapemis_flat(1,1,1,iem),MSG_READ1,   &
               1,GIMAX,GJMAX,FNCMAX,1,1)

    end do ! iem = 1, NEMIS-loop


    if ( me == 0 ) then
        write(unit=6,fmt=*) "Country totals"
        write(unit=IO_LOG,fmt=*) "Country totals"
        write(unit=6,fmt="(2a4,3x,10a12)")  " N"," CC ",(EMIS_NAME(iem),iem=1,NEMIS)
        write(unit=IO_LOG,fmt="(2a4,3x,10a12)") " N"," CC ",(EMIS_NAME(iem),iem=1,NEMIS)

        do ic = 1, NLAND
           ccsum = sum( sumemis(ic,:) )
           if ( ccsum > 0.0 ) then
                    write(unit=6,fmt="(i3,1x,a4,3x,8f12.2)") &
                        ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS)
                    write(unit=IO_LOG,fmt="(i3,1x,a4,3x,8f12.2)")& 
                        ic, Country(ic)%code, (sumemis(ic,i),i=1,NEMIS)
           end if
        end do
    end if

    !su    now all values are read, snapemis is distributed, globnland and 
    !su    globland are ready for distribution
    !gv    print *, "calling glob2local_int for iem", iem, " me ", me

      call global2local_int(globnland,nlandcode,326, GIMAX,GJMAX,1,1,1)
      call global2local_int(globland, landcode ,326, GIMAX,GJMAX,NCMAX,1,1)
!hf
      call global2local_int(flat_globnland,flat_nlandcode,326,&
                            GIMAX,GJMAX,1,1,1)!extra array
      call global2local_int(flat_globland,flat_landcode,326,&
                            GIMAX,GJMAX,FNCMAX,1,1)

     !/**hf  broadcast volcanoe info derived in EmisGet 

      call gc_ibcast(327,1,0,NPROC,info,nvolc)
      call gc_ibcast(327,nvolc,0,NPROC,info,i_volc)
      call gc_ibcast(327,nvolc,0,NPROC,info,j_volc)
      call gc_rbcast(327,nvolc,0,NPROC,info,emis_volc)



!    Conversions --
!
!     The emission-data file are so far in units of 
!     tonnes per grid-square. The conversion factor from tonnes per 50*50km2
!     annual emission values to surface flux (kg/m2/s) is found by division
!     with (nydays*24*60*60)s and (h*h)m2 and multiply by 1.e+3.
!     the conversion factor then equals 1.27e-14
!
    if ( DEBUG) print *,  "CONV:me, nydays, gridwidth = ",me,nydays,GRIDWIDTH_M

    tonne_to_kgm2s  = 1.0e3 / (nydays * 24.0 * 3600.0 * GRIDWIDTH_M * GRIDWIDTH_M)

    if ( DEBUG .and.  me ==  0 ) then
        write(unit=6,fmt=*) "No. days in Emissions: ", nydays
        write(unit=6,fmt=*) "tonne_to_kgm2s in Emissions: ", tonne_to_kgm2s
        write(unit=6,fmt=*) "Emissions sums:"
        write(unit=6,fmt=*) (EMIS_NAME(iem),iem=1,NEMIS)
        write(unit=6,fmt=*) (emsum(iem),iem=1,NEMIS)
    endif


    do iem = 1, NEMIS
       conv = tonne_to_kgm2s * EmisDef( eindex(iem) )%conv
 
       forall (ic=1:NCMAX, j=lj0:lj1, i=li0:li1, isec=1:NSECTORS)
          snapemis (isec,i,j,ic,iem) = &
                 snapemis (isec,i,j,ic,iem) * conv * xm2(i,j)
       end forall

       forall (fic=1:FNCMAX, j=lj0:lj1, i=li0:li1)
          snapemis_flat(i,j,fic,iem) = &
                 snapemis_flat(i,j,fic,iem) * conv * xm2(i,j)
       end forall
    enddo !iem

!hf Volc
    if ( VOLCANOES ) then    !u3

       conv = tonne_to_kgm2s * EmisDef( eindex(QRCVOL) )%conv   !st - bug found

       do volc_no=1,nvolc 
          i=i_volc(volc_no)
          j=j_volc(volc_no)
          !Find global<->local coordinates for xm2
             if ((i >= gi0).and.(i<=gi1).and.(j>= gj0).and.&
                 (j<= gj1))then !on the correct processor
                write(*,*)'i,j for volcanoe is',i,j
                write(*,*)'EMIS_VOLC is',emis_volc(volc_no)
                i_l = i -gi0 +1
                j_l = j- gj0 +1
                write(*,*)'Local coord is',i_l,j_l,gi0,gj0
                emis_volc(volc_no) = emis_volc(volc_no)* conv * xm2(i_l,j_l)
            endif
        enddo !volc_no

        !/** Read  Volcano.dat to get volcano height
        if (me==0)then
            call VolcGet(height_volc)
            write(*,*)'Volcano heights',height_volc
        endif

        !/** broadcast volcano heights
        call gc_ibcast(327,NMAX_VOLC,0,NPROC,info,height_volc)
     endif ! VOLCANOES

    err1 = 0
    if ( me == 0 ) then
       deallocate(globnland,stat=err1)
       deallocate(globland ,stat=err2)
       deallocate(globemis ,stat=err3)
!hf
       deallocate(flat_globnland,stat=err4)
       deallocate(flat_globland,stat=err5)
       deallocate(globemis_flat,stat=err6)
       if ( err1 /= 0 .or. err2 /= 0 .or. err3 /= 0 .or. &
            err4 /= 0 .or. err5 /= 0 .or. err6 /= 0        ) then
          call gc_abort(me,NPROC,"De-Alloc error - globland")
       end if
    end if

  end subroutine Emissions

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine consistency_check(eindex)
  use My_Emis_ml,  only : &
          EMIS_NAME,    &   ! Names of species ("sox  ",...)
          NEMIS_PLAIN,  &   ! No. emission files which are not speciated
          EMIS_NSPLIT,  &   ! No. emission files to be speciated
          molwt             ! Mol. wts
  !------------------------------------------------------------------!
  !    checks that all the values given so far are consistent        !
  !------------------------------------------------------------------!
  integer, dimension(NEMIS), intent(in) :: eindex
  character(len=30) :: errormsg
  integer :: i

  errormsg = "ok"
  do i = 1, NEMIS
     if ( eindex(i) < 0 ) then
          print *, "EmisIndex: Mis-match for ", i, eindex(i)
          errormsg = "EmisIndex: Mismatch"
     end if
  end do
  if ( NRCEMIS < NEMIS             ) errormsg = " NRCEMIS < NEMIS"
  if ( size(EMIS_NAME) /= NEMIS    ) errormsg = " size EMISNAME wrong "
  if ( NEMIS_PLAIN+sum(EMIS_NSPLIT) /= NRCEMIS   ) errormsg = "sum ne NRCEMIS"
  if ( any( molwt < 1.0 )          ) errormsg = " Mol. wt not assigned "

  if ( errormsg /= "ok" ) then
     write(6,*) "Failed consistency check:", errormsg
     call gc_abort(me,NPROC,errormsg)
  endif

 end subroutine consistency_check
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  !*****************************************************************

  subroutine EmisSet(indate)   !  emission re-set every time-step/hour

  use My_Emis_ml,  only : &
          NEMIS_PLAIN,  &   ! No. emission files which are not speciated
          EMIS_NSPLIT,  &   ! No. emission files to be speciated
          NFORESTVOC,   &   ! > 0 if forest voc wanted
          molwt             ! Mol. wts
  use Biogenics_ml, only:               &
          emnat,emforest ,&               ! Emission arrays
          IQ_DMS
  use Country_ml,    only : NLAND,Country
  use Dates_ml,      only : date     ! No. days per year, date-type 
  use EmisDef_ml,    only : &  
                     VERTFAC          ! vertical emission split

  use EmisGet_ml, only : emisfrac  ! speciation routines and array

  use Timefactors_ml, only : &
               NewDayFactors   &         ! subroutines
              ,timefac, day_factor  ! time-factors

  use Par_ml,                only : li0,lj0,li1,lj1
  use PhysicalConstants_ml,  only :  GRAV,  AVOG
  use ModelConstants_ml,     only :  KMAX_BND, PT ,dt_advec
  use GridValues_ml  ,       only :  sigma_bnd, xmd,GRIDWIDTH_M
  use Met_ml,                only :  ps, roa   ! ps in Pa, roa in kg/m3
  !
  !***********************************************************************
  !**    DESCRIPTION:
  !     calculates the emmision-tendencies and the local (instantaneous) dry 
  !     deposition in the emission squares.
  !       emis set once per hour to allow for day/night variation (and voc 
  !       speciation) (based on local time)  for each snap sector.
  !     gridrcemis0 calculated every time-step to allow for ps changes.
  !uni inputs from Emissions in EMISSIONS_ML:
  !      country and snap-specific array : 
  !          snapemis (NSECTORS,MAXLIMAX,MAXLJMAX,NCMAX,NEMIS) 
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
  !       real timefac(NLAND,NSECTORS,NEMIS)
  !    And day-night factors are applied here:
  !       day_factor(11,0:1)                  ! 0=night, 1=day
  !ds ..........................................................................
  !
  !**    REVISION HISTORY:
  !       Revised , 30/5/01, jej/st found problem on gridur - split NEMIS loop 
  !       into separate NEMIS_PLAIN and NEMIS_SPLIT loops.
  !       Revised, ds,  Feb. 2001 for unified model. Use of date%seconds replaces
  !       thourloc.
  !       Revised  : d. simpson 4/2/98 to act as common subrouinte
  !       between 3-D models, and to avoid hard-coded emissions
  !       !uni - revised to F90 and for more flexible handling of emissions
  !       fraction through NEMIS_FRAC
  !       Originally emission.f from MADE and MACHO.
  !
  !      25/3-2002, pw changed test for hour and day change (Now the first day
  !      does not need to start at 0 hours)
  !
  !*************************************************************************

  implicit none
  type(date), intent(in) :: indate           ! Gives year..seconds
  integer :: i, j, n, k, f                   ! cooridnates, loop variables
  integer :: icc, ncc                        ! No. of countries in grid.
!hf F
  integer :: ficc,fncc                       ! No. of countries with
  integer :: iqrc, ifrac                     ! emis indices 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over pollutants (1..NEMIS)
!hf
  integer :: i_l,j_l           ! Local i,j
  !uni - save daytime value between calls, intiialise to zero
  integer, save, dimension(NLAND) ::  daytime = 0  !  0=night, 1=day
  integer                         ::  hourloc      !  local hour 
!u4   integer                         ::  olddaytime   !      "     " 
  logical                         ::  hourchange   !      "     "           
  real, dimension(NRCEMIS)        ::  emis         !  local array for emissions

  real ::  deploc,ehlpcom,ehlpcom0(KEMISTOP:KMAX_MID)
  real ::  tfac, dtgrid    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s               ! source term (emis) before splitting
  integer :: iland        ! country codes 

!hf
  real ::  ftfac           ! time-factor for flat emissions
  real ::  sf              ! source term (emis) before splitting (for flat emissions)
  integer :: flat_iland    ! country codes (countries with flat emissions)

  integer, save :: oldday = -1, oldhour = -1
 
!hf initialize
    ehlpcom0(:)=0.0
!TIL KMAX_BND???(endre dim)
!hf u2   do k=1,4
!u4   do k=1,KMAX_MID !due to volcanos
!u4 - ds - ehlpcom0 isn't used in Volvano_ml - it has its own version.
   !u4 do k=1,KMAX_MID !due to volcanos
   !u4 - and now we make use of the new KEMISTOP parameter to define the
   !u4   vertical range. Note that ehlpcom0 nopw has value KMAX_MID=20 at ground level.

   do k=KEMISTOP,KMAX_MID !due to volcanos
      ehlpcom0(k) = GRAV* 0.001*AVOG/ (sigma_bnd(k+1) - sigma_bnd(k))
      !u4 ehlpcom0(k) = GRAV* 0.001*AVOG/ (sigma_bnd(KMAX_BND-k+1) - sigma_bnd(KMAX_BND-k))
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
!pw u3    if ( indate%seconds == 0 ) then
!
!           hourchange = .true.
!                                     !==========================
!           if ( indate%hour == 0  )  call NewDayFactors(indate)
!                                     !==========================
!pw u3     end if


    !..........................................
    !  Look for day-night changes, after local time correction

    do iland = 1, NLAND

       !u4 olddaytime     = daytime(iland)
       daytime(iland) = 0
       hourloc        = indate%hour + Country(iland)%timezone

       if ( hourloc  >=   7 .and.  hourloc <= 18 ) daytime(iland)=1

    end do ! iland


    if ( hourchange ) then 

         totemadd(:)  = 0.
         gridrcemis0(:,:,:,:) = 0.0 


        !..........................................
        ! Process each grid:

        do j = lj0,lj1
            do i = li0,li1

               ncc = nlandcode(i,j)            ! No. of countries in grid
               !u4 ncc = nlandcode(i,j)       ! No. of countries in grid
              !*************************************************
              ! First loop over non-flat(one sector) emissions
              !*************************************************

              emis(:)=0.
              do icc = 1, ncc
                  iland = landcode(i,j,icc)     ! 1=Albania, etc.

                 !  As each emission sector has a different diurnal profile
                 !  and possibly speciation, we loop over each sector, adding
                 !  the found emission rates to gridrcemis as we go.
                 !  ==================================================
                do isec = 1, NSECTORS       ! Loop over snap codes


                   !--  Calculate emission rates from snapemis, time-factors, 
                   !    and if appropriate any speciation  fraction (NEMIS_FRAC)

                   iqrc  = 0   ! index over emis
                   ifrac = 0   ! index over emisfrac

		    !/.. First, the simple emissions
                   do iem = 1, NEMIS_PLAIN

                      tfac = timefac(iland,isec,iem) * &
                                 day_factor(isec,daytime(iland))

                      iqrc = iqrc + 1
                      emis(iqrc) = snapemis(isec,i,j,icc,iem) * tfac 

                   end do ! iem=1,NEMIS_PLAIN

		    !/.. Then , the split (speciated) emissions if NEMIS_SPLIT>0

                   do iem = 1, NEMIS_SPLIT

                      tfac = timefac(iland,isec,iem+NEMIS_PLAIN ) * &
                                 day_factor(isec,daytime(iland))

                      s =  tfac * snapemis(isec,i,j,icc,iem+NEMIS_PLAIN)

                      do f = 1, EMIS_NSPLIT( iem )
                           ifrac = ifrac + 1
                           iqrc  = iqrc  + 1
                           emis(iqrc) = s * emisfrac(ifrac,isec,iland)
                      end do ! f

                   end do ! iem=1,NEMIS_SPLIT

                !--  Add up emissions in ktonne ......

                   do iqrc = 1, NRCEMIS
                      totemadd(iqrc) = totemadd(iqrc) + &
                                          emis(iqrc) * dtgrid * xmd(i,j)
                   end do

                   !..   Assign to height levels 1-4

                   !u4 do k=1,4
                   !u4 - rewrite with k ranging over range defined in ModelConstants_ml

                   do k=KEMISTOP,KMAX_MID
                      do iqrc =1, NRCEMIS
                         gridrcemis0(iqrc,k,i,j) =   &
                            gridrcemis0(iqrc,k,i,j) + emis(iqrc)*   &
                            ehlpcom0(k)*VERTFAC(KMAX_BND-k,isec) /molwt(iqrc)  
                      end do ! iem
                   end do   ! k

                enddo  ! isec
!      ==================================================

       end do ! icc  
 
       !************************************
       ! Then loop over flat emissions
       !************************************
       emis(:)=0.
       fncc = flat_nlandcode(i,j) ! No. of countries with flat 
                                          ! emissions in grid

       do ficc = 1, fncc
          flat_iland = flat_landcode(i,j,ficc) ! 30=BAS etc.

          if ( Country(flat_iland)%is_sea ) then   ! u4 - saves if statements below
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
            ifrac = 0   ! index over emisfrac

          !/.. First, plain emissions

             do iem = 1, NEMIS_PLAIN
                iqrc = iqrc + 1

                !u4 - ds we don't have time-factors for ships, just for
                !u4 "other mobile sources", which is likely based upon
                !u4 tractors and other land-based vehicles. Therefore
                !u4 for now we can just skip the time-factors.

                !u4 if ( Country(flat_iland)%is_sea ) then ! add ship emissions
                !u4    ftfac = timefac(flat_iland,ISNAP_SHIP,iem) * &
                !u4          day_factor(ISNAP_SHIP,daytime(flat_iland))
                !u4    emis(iqrc) = snapemis_flat(i,j,ficc,iem) * ftfac
                !u4elseif (flat_iland == IQ_DMS)then ! DMS
                !No timevariation for DMS
!Bor man teste paa sox?                          

                    emis(iqrc) =  snapemis_flat(i,j,ficc,iem) 

                !u4 endif !ship emissions

             end do ! iem=1,NEMIS_PLAIN

            !/.. Then , the split (speciated) emissions if NEMIS_SPLIT>0

             do iem = 1, NEMIS_SPLIT

                !u4 ftfac = timefac(flat_iland,ISNAP_SHIP,iem+NEMIS_PLAIN) * &
                !u4              day_factor(ISNAP_SHIP,daytime(flat_iland))
                !u4 sf =  ftfac * snapemis_flat(i,j,ficc,iem+NEMIS_PLAIN)    

                sf =  snapemis_flat(i,j,ficc,iem+NEMIS_PLAIN)    


                do f = 1, EMIS_NSPLIT( iem )
                   ifrac = ifrac + 1
                   iqrc  = iqrc  + 1

                   !u4 if ( Country(flat_iland)%is_sea )then ! ship emissions
                   !u4      emis(iqrc) = sf * emisfrac(ifrac,ISNAP_SHIP,flat_iland)
                   !u4 else ! other flat emissions than ships, &
                   !u4      ! does not have split emissions
                   !u4      emis(iqrc) = 0.
                   !u4 endif
!kunne satt alle=sf*emisfrac siden snapemis 0 for alle andre
!u4
!u4 - ds - yes, and since EMIS_NSPLIT would be zero there is no need for an 
!     extra if-test within the loop. It is enough to say simply:

                     emis(iqrc) = sf * emisfrac(ifrac,isec,flat_iland)

                end do ! f

             end do ! iem=1,NEMIS_SPLIT

         !--   Add flat emissions in ktonne  (to non-flat emissions)......

             do iqrc = 1, NRCEMIS

                totemadd(iqrc) = totemadd(iqrc) + &
                                   emis(iqrc) * dtgrid * xmd(i,j)
             end do


         !..   Assign flat emissions to height levels 1-4
         !..   Note, no VERTFAC


             do iqrc =1, NRCEMIS
                !u4 gridrcemis0(iqrc,KMAX_BND-1,i,j) =   &
                !u4 gridrcemis0(iqrc,KMAX_BND-1,i,j) + emis(iqrc)*&
                !u4 ehlpcom0(1)/molwt(iqrc)
                !u4 - simpler use of KMAX_MId and new ehlpcom0:

                gridrcemis0(iqrc,KMAX_MID,i,j) =   &
                  gridrcemis0(iqrc,KMAX_MID,i,j) + emis(iqrc)*&
                    ehlpcom0(KMAX_MID)/molwt(iqrc)
             end do ! iem

!      ==================================================
       end do !ficc 
   end do ! i
 end do ! j

!hf
    if ( VOLCANOES ) call Set_Volc !set hourly volcano emission(rcemis_volc0)

  end if ! hourchange 


  !/** we now scale gridrcemis to get emissions in molecules/cm3/s
  !u4   do k=1,4
   do k= KEMISTOP, KMAX_MID
     do j = lj0,lj1
        do i = li0,li1
              !u4 ehlpcom= roa(i,j,KMAX_BND-k,1)/(ps(i,j,1)-PT)
              ehlpcom= roa(i,j,k,1)/(ps(i,j,1)-PT)
              do iqrc =1, NRCEMIS
                 !u4 gridrcemis(iqrc,KMAX_BND-k,i,j) =  &
                 !u4    gridrcemis0(iqrc,KMAX_BND-k,i,j)* ehlpcom
                 gridrcemis(iqrc,k,i,j) =  &
                    gridrcemis0(iqrc,k,i,j)* ehlpcom
              enddo  ! iqrc
        end do ! i
      end do ! j
   end do ! k

!hf u2
 !/** Scale volc emissions to get emissions in molecules/cm3/s (rcemis_volc)
   if ( VOLCANOES ) call Scale_Volc

 if( NFORESTVOC > 0  )then

    do j = lj0,lj1
      do i = li0,li1

        !u4 ehlpcom = ehlpcom0(1) * roa(i,j,KMAX_MID,1)/(ps(i,j,1)-PT)

        ehlpcom = ehlpcom0(KMAX_MID) * roa(i,j,KMAX_MID,1)/(ps(i,j,1)-PT)

        emnat(1:NFORESTVOC,i,j) = emforest(1:NFORESTVOC,i,j)*ehlpcom

      end do ! i
    end do ! j
!   if ( DEBUG ) then 
!       print "(a12,2i4,/,(g12.3))",  "bio-setemis",li0,lj0, &
!          (gridrcemis(i,KMAX_MID,2,2),i=1, NRCEMIS), emnat(1,2,2),emnat(2,2,2)
!   end if

  endif ! NFORESTVOC 

 end subroutine EmisSet
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	subroutine newmonth              ! (exclu)

!ds.....................................................................
!**    DESCRIPTION:
!     Reads in natural DMS emissions at start of each month. Update
!     landcode and nlandcode arrays as needed.

!     Reads in snow cover at start of each month. (not yet used in
!	MACHO, but read-in code implemented here).

!**    REVISION HISTORY:
!     Original from MADE

!     ds, Feb 99:
!     Previous code to adjust the time-series factors has been removed
!     this work is now done in setemis.f. Thus newmonth has substantially
!     fewer jobs than in previous code.
!      ds Feb 01
!      MADE/MACHO differences removed. Snow and DMS read-in for both.

!ds...........................................................................

	use Par_ml   , only : IILARDOM,JJLARDOM		&
		,MSG_READ7,ISMBEG,JSMBEG	&
		,limax,ljmax
	use Dates_ml, only:		&
		nydays          ! = 365.0 or 366.0  (days per year)
        use ModelConstants_ml, only : current_date
        use GridValues_ml, only:			&
		GRIDWIDTH_M    & !  size of grid (m)
		,xm2             ! map factor squared

	use Io_ml, only :   IO_DMS,ios, open_file

	use My_Emis_ml,    only : EMIS_NAME

        use Biogenics_ml, only:		&
		IQ_DMS              & ! 
		!u4 ,dms_fact            & !  ... document/remove ...??
		,first_dms_read       !  ... document/remove ...??

        use ReadField_ml, only : ReadField    ! Reads ascii fields
	implicit none

	integer i, j
	integer ijin(2) 
!hf	integer n, ncmaxfound		! Max. no. countries found in grid
        integer n, flat_ncmaxfound      ! Max. no. countries w/flat emissions
!6b	real rtmp(IILARDOM,JJLARDOM), rsnd(MAXLIMAX,MAXLJMAX)
	real :: rdemis(MAXLIMAX,MAXLJMAX)  ! Emissions read from file
	character*20 fname
	real ktonne_to_kgm2s         ! Units conversion
        integer :: IQSO2             ! Index of sox in  EMIS_NAME
	integer errcode

!*** Units:
!	Input files seem to be in ktonne. We convert here to kg/m2/s
! 	to save CPU in setemis.f.
!      The conversion factor from 50*50km2
!      annual emission values to surface flux (kg/m2/s) is found by division
!      with (nydays*24*60*60)s and (h*h)m2 and multiply by 1.e+6.
!      the conversion factor (ktonne_to_kgm2s) then equals 1.27e-8 

        ktonne_to_kgm2s  = 1.0e6 / 				&
		(nydays*24.*60.*60.*GRIDWIDTH_M*GRIDWIDTH_M)

	if ( me .eq. 0 ) then
	  write(6,*) '    Newmonth           me ', me
	  write(6,*) 'Enters newmonth, ktonne_to_kgm2s = '	&
			,ktonne_to_kgm2s
	  write(6,*) ' first_dms_read = ', first_dms_read
	end if ! me 
!ds...........................................................................

!ds...........................................................................
!		DMS Input - land 35 - SNAP sector 11
!ds...........................................................................
     flat_ncmaxfound = 0 ! Max. no. countries(w/flat emissions) per grid
!    natural so2 emissions

          IQSO2 = 0
          do i = 1, NEMIS
            if ( trim( EMIS_NAME(i) ) == "sox" ) IQSO2 = i
          end do

          if ( IQSO2 < 1 ) then
              write(*,*) " No SO2 emissions - need to skip DMS also"
              return    ! No need to read DMS fields ...
              !u1 call stop_test(.true.,me,NPROC,99 ,"IQ_SO2  not assigned")

          else    
            !/--- we have so2 emission so need DMS also

             if ( me == 0 ) then

	          write(fname,fmt='(''natso2'',i2.2,''.dat'')') 	&
			current_date%month
	          write(6,*) 'filename for nat-so2',fname
              endif

              call ReadField(IO_DMS,fname,rdemis)

	      errcode = 0
	      do j=1,ljmax
	          do i=1,limax

!ds            Add DMS for country code IQ_DMS=35  to snap sector 11=Nature.
                ! First time we read we must add DMS to the "countries" 
                ! contributing within the grid square. 
 
                  !hf - for flat emissions:

                  if ( first_dms_read ) then 
                     flat_nlandcode(i,j) = flat_nlandcode(i,j) + 1 
                     n = flat_nlandcode(i,j)
                     flat_landcode(i,j,n) = IQ_DMS   ! country code 35 
                     if ( n > flat_ncmaxfound ) then
                          flat_ncmaxfound = n 
                          write(6,*)'DMS Increased flat_ncmaxfound to ',n 
                          if ( n > FNCMAX ) then
		  	    call gc_abort(me,NPROC,'Increase FNCMAX for dms')
			  endif
                     endif 
                  else  ! We know that DMS lies last in the array, so:
                      n = flat_nlandcode(i,j) 
                      if ( flat_landcode(i,j,n) /= IQ_DMS )then
                         call gc_abort(me,NPROC,'newmonth: DMS not last!')
		      endif 
                  endif 

                  snapemis_flat(i,j,n,IQSO2) = rdemis(i,j) * ktonne_to_kgm2s &
			* xm2(i,j)
	        enddo ! i
	      enddo ! j


              if ( first_dms_read ) then
   	          write(6,*)'me ',me, ' Increased flat_ncmaxfound to ' &
			,flat_ncmaxfound 
	          first_dms_read = .false.
	      end if

           end if  ! IQSO2>0

	end subroutine newmonth

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module Emissions_ml
