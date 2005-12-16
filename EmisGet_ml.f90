module EmisGet_ml
  use My_Emis_ml,   only : NEMIS  &
                          , NRCSPLIT & 
                          , EMIS_NAME, SPLIT_NAME  &
                          , NEMIS_PLAIN, NEMIS_SPLIT, EMIS_NSPLIT
  use Country_ml,   only : NLAND    & !u4 movde other stuff here
                            ,IC_NAT,IC_VUL, Country
  use EmisDef_ml,   only : NSECTORS, ANTROP_SECTORS, NCMAX, FNCMAX & 
                            ,ISNAP_SHIP, ISNAP_NAT  !u3 for NAT
  use Functions_ml, only : GridAllocate       !ds u7.2
  use Io_ml,        only : open_file,  wordsplit &      ! subs
                          ,NO_FILE, ios, IO_EMIS         ! variables
  use Par_ml,       only : me, NPROC
  use Volcanos_ml                      !hf

  implicit none
  private


  public :: EmisGet            ! Collects emissions of 1 pollutant
  public :: EmisSplit          ! => emisfrac, Speciation of voc, pm25, etc.
  private :: femis             ! Sets emissions control factors 


  logical, private, save :: my_first_call = .true.

  logical, private, parameter :: DEBUG = .false.

  ! e_fact is read in from the femis file and applied within EmisGet
  real, private, save, &
         dimension(NSECTORS,NLAND,NEMIS)  :: e_fact ! emission factor

  !/* emisfrac is used at each time-step of the model run to split
  !    emissions such as VOC; PM into species. 

   real, public, dimension(NRCSPLIT,NSECTORS,NLAND), save :: emisfrac


  !/ some common variables
  character(len=40), private :: fname             ! File name
  character(len=40), private :: errmsg            ! Error message!

contains
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  subroutine EmisGet(iemis,emisname,ISMBEG,JSMBEG,GIMAX,GJMAX, &
                      globemis,globnland,globland,sumemis,     &
                      globemis_flat,flat_globnland,flat_globland)

!.......................................................................
!**    DESCRIPTION:
!  Read in emissions from one file, specified by iemis. The arrays read
!  in here are the global arrays (allocatable)
!
!**    REVISION HISTORY:
!      1/5/02 u7.2 Code moved to Functions/GridAllocate, ds
!      .../2002  Flat emissions added, hf
!      25/1/01   Re-coded as separate module, and for F
!      30/5/00   Created from earlier readem, ds
!
!** Language: F
!.......................................................................

  !--arguments
  integer, intent(in)                      :: iemis     ! emis index 
  character(len=*), intent(in)             :: emisname  ! emission name
  integer, intent(in) :: ISMBEG,JSMBEG,GIMAX,GJMAX        ! domain limits
  real,    intent(out), dimension(:,:,:,:) :: globemis  ! emission values
  integer, intent(inout), dimension(:,:,:) :: globland  ! C'try-code -emitters
  integer, intent(inout), dimension(:,:)   :: globnland ! No. emit's in grid
  real,    intent(inout), dimension(:,:)   :: sumemis   ! Sums (after e_fact) 
!hf V
  integer, intent(inout), dimension(:,:,:) :: flat_globland  ! C'try-code -emitters
  integer, intent(inout), dimension(:,:)   :: flat_globnland ! No. flat emit's in grid
  real,    intent(out), dimension(:,:,:)   ::globemis_flat !Flat emissions

  !--local
!hf V
  integer :: flat_iland, flat_nc
  integer :: i, j, isec, iland, k, nc        ! loop variables
  integer :: iic,ic                          ! country code (read from file)
  real    :: duml,dumh                   ! dummy variables, low/high emis.
  real, dimension(NSECTORS)  :: tmpsec   ! array for reading emission files
  integer, save :: ncmaxfound = 0        ! Max. no. countries found in grid
!hf V
  integer, save :: flat_ncmaxfound = 0    ! Max. no. countries(with flat emissions) found in grid
     !>============================
       if ( my_first_call ) then
            sumemis(:,:) =  0.0       ! Initialise sums
	    ios = 0
            call femis()              ! emission factors (femis.dat file).
	    if ( ios /= 0 )return
            my_first_call = .false.
       endif
     !>============================

        globemis   (:,:,:,:) = 0.0
        globemis_flat(:,:,:) = 0.0 !hf F

        write(unit=6,fmt=*) "Called EmisGet with index, name", iemis, emisname
        fname = "emislist." // emisname
        call open_file(IO_EMIS,"r",fname,needed=.true.)
           if ( ios /= 0 ) call gc_abort(me,NPROC,"ios error: emislist") !u4

READEMIS: do   ! ************* Loop over emislist files **********************

              read(unit=IO_EMIS,fmt=*,iostat=ios) iic,i,j, duml,dumh,  &
                                 (tmpsec(isec),isec=1,NSECTORS)

              if ( ios <  0 ) exit READEMIS            ! End of file
              if ( ios >  0  ) then                     ! A different problem..
                  errmsg = "GetEmis ios: error on " // fname
                  call gc_abort(me,NPROC,errmsg)
              end if

!find country number (ic) corresponding to index as written in emisfile (iic)
              do ic=1,NLAND
                 if(Country(ic)%index==iic)goto 543
              enddo 
              if(me==0)write(*,*)'COUNTRY CODE NOT RECOGNIZED',iic
              ic=0
543           continue

              i = i-ISMBEG+1     ! for RESTRICTED domain
              j = j-JSMBEG+1     ! for RESTRICTED domain

              if ( i  <=  0 .or. i  >  GIMAX .or.   & 
                   j  <=  0 .or. j  >  GJMAX .or.   &
                   ic <=  0 .or. ic >  NLAND .or.   &
                   ic == IC_NAT )                   &  ! Excludes DMS
               cycle READEMIS
             !/**hf ship emissions

             if ( Country(ic)%is_sea ) then    ! ship emissions

                ! ..........................................................
                ! generate new land allocation in 50 km grid for FLAT 
                ! EMISSIONS(ships). First, we check if
                ! country "ic" has already  been found within that grid. If not,
                ! then ic is added to flat_landcode and flat_nlandcode 
                ! increased by one.
   
                !/** Test that ship emissions are only in sector ISNAP_SHIP
                do isec=1,(ISNAP_SHIP-1) 
                   if ( tmpsec(isec) /= 0.) &
                      call gc_abort(me,NPROC,"EMISGET1, NOT FLAT EMISSION") 
                enddo
                do isec=ISNAP_SHIP+1,NSECTORS
                   if ( tmpsec(isec) /= 0. ) &
                      call gc_abort(me,NPROC,"EMISGET2, NOT FLAT EMISSION") 
                enddo
                !/** end test

                call GridAllocate("FLat",i,j,ic,FNCMAX, flat_iland, &
                    flat_ncmaxfound,flat_globland,flat_globnland,errmsg)
                if ( errmsg /= "ok" ) call gc_abort(me,NPROC,errmsg)

              ! ...................................................
              ! Assign e_fact corrected emissions to global FLAT 
              ! emission matrices.
              ! ...................................................

                 globemis_flat(i,j,flat_iland) = globemis_flat(i,j,flat_iland) &
                     + e_fact(ISNAP_SHIP,ic,iemis) *  tmpsec(ISNAP_SHIP)

              !......................................................
              !..        Sum over all sectors, store as ktonne:
              !......................................................

                 sumemis(ic,iemis) = sumemis(ic,iemis)   &
                      + 0.001 * globemis_flat(i,j,flat_iland)  

                cycle READEMIS
             endif !ship emissions            


              !.......................................................
              !/**hf  Volcanos
              !.......................................................

              if ( trim ( emisname ) == "sox" )then
                if (ic == IC_VUL) then
                   volc_no=volc_no+1
                   write(*,*)'Volcano no. ',volc_no
                   i_volc(volc_no)=i
                   j_volc(volc_no)=j
                   !write(*,*)'Volcano i,j ',i_volc(volc_no),j_volc(volc_no)
                   emis_volc(volc_no)=tmpsec(ISNAP_NAT) * &
                                         e_fact(ISNAP_NAT,IC_VUL,iemis)
                   nvolc=volc_no
                   if (nvolc>NMAX_VOLC)then
                      call gc_abort (me,NPROC,"EMISGET, NMAX_VULC: STOP")
                   endif
                   write(*,*)'Found ',nvolc,' volcanoes on sox file'


                   sumemis(IC_VUL,iemis) = sumemis(IC_VUL,iemis)   &
                    + 0.001 * emis_volc(volc_no)
                   cycle READEMIS ! do not want to count volcano "landcode"
                endif ! ic
             endif ! so2

             !..............................................................
             !/**hf end volcanoes
             !..............................................................


             !  For VOC natural and agricultur emissions (managed forests) 
             !  set to  zero

!rv2_4_1         if ( trim ( emisname ) == "voc" ) tmpsec(10:11) = 0.0
             if ( trim ( emisname ) == "voc" ) tmpsec(11:11) = 0.0

             ! ..........................................................
             ! generate new land allocation in 50 km grid. First, we check if
             ! country "ic" has already  been found within that grid. If not,
             ! then ic is added to landcode and nlandcode increased by one.

              call GridAllocate("SNAP",i,j,ic,NCMAX, &
                                 iland,ncmaxfound,globland,globnland,errmsg)
                if ( errmsg /= "ok" ) call gc_abort(me,NPROC,errmsg)

              ! ...................................................
              ! ...................................................
              !su        corrected:

                 globemis(:,i,j,iland) = globemis(:,i,j,iland) &
                        + e_fact(:,ic,iemis) *  tmpsec(:)


              !..        Sum over all sectors, store as ktonne:

               sumemis(ic,iemis) = sumemis(ic,iemis)   &
                 + 0.001 * sum ( globemis (:,i,j,iland) )
              
        end do READEMIS 
        !
        close(IO_EMIS)
	ios = 0
  end subroutine EmisGet
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!u7.2 ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!u7.2   subroutine LandAllocate(label,i,j,ic,ncmax,iland,ncmaxfound,land,nland)
!u7.2 
!u7.2     !-- Checks if a country "ic" whose data has just been read in has
!u7.2     !   already been found within the given emissions square.
!u7.2     !   If not, the array "nland" is incremented by one and the
!u7.2     !   country index added to "land".
!u7.2  
!u7.2      character(len=*), intent(in) :: label   ! Type of data
!u7.2      integer, intent(in) :: i,j
!u7.2      integer, intent(in) :: ic        ! Index of country just read in
!u7.2      integer, intent(in) :: ncmax     ! Max. no countries allowed
!u7.2 
!u7.2      integer, intent(out)   :: iland         ! Index of country in that grid
!u7.2      integer, intent(inout) :: ncmaxfound    ! No. countries found so far
!u7.2      integer, dimension(:,:,:), intent(inout) :: land   ! Land-codes
!u7.2      integer, dimension(:,:),   intent(inout) ::nland   ! No. countries
!u7.2 
!u7.2      integer :: nc, k, iland      ! local variables
!u7.2 
!u7.2        nc=nland(i,j)       ! nc = no. countries known so far
!u7.2 
!u7.2        do k = 1,nc
!u7.2           if( land(i,j,k) == ic) then
!u7.2               iland = k        ! country is already in the list
!u7.2               goto 100
!u7.2           endif
!u7.2        enddo
!u7.2 
!u7.2        nland(i,j) = nland(i,j) + 1    ! country is a new one
!u7.2        land(i,j,nc+1) = ic
!u7.2        iland=nc+1
!u7.2 
!u7.2        if( iland >  ncmaxfound) then
!u7.2            ncmaxfound = iland
!u7.2            write(*,*) "LandAlloc ", label, "increased ncmaxfound:",i,j,iland
!u7.2            write(*,*) "LandAlloc ", label," now have:", &
!u7.2                            (land(i,j,k),k=1,ncmaxfound)
!u7.2            if ( ncmaxfound >  ncmax ) then
!u7.2 		  print *,"LandAlloc ncmax ERROR ", label, "NCMAX"
!u7.2                   call gc_abort (me,NPROC,"EMISGET, NCMAX")
!u7.2            endif
!u7.2         endif
!u7.2  100    continue
!u7.2 
!u7.2   end subroutine LandAllocate
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine femis()
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

! Reads adjustment factors for the emissions from (optional) file femis. 
! Emission factors are applied to specified country and/or emission sector. 
!  e.g. the femis file with input:
!      Code  5  co  sox  nox    nh3   voc
!      27    6  1.0  0.5  1.0   1.0   1.0
!   would reduce SO2 emissions from sector 6 by a factor 0.5 for the UK 
!   (country 27). A zero for country or sector means apply factors to all 
!   countries and/or sectors. The number following the first text on line 1 
!   (number 5 above) gives the number of pollutants treated in the file (ncols
!   below). Note that ncols can be greater than the number of emitted species
!   we actually have, and the species can be specified in any order, as given
!   on the top line. (Thus, we can use the same femis file for MADE,MACHO
!   or AERO.

  !su    (for compatibility NEMIS as last index)

  !/** local variables **/
  integer            :: ie, iq, ic, iland1, iland2 ! loop variables
  integer            :: inland                     ! Country read from femis
  integer            :: isec, isec1 , isec2        ! loop vars: emis sectors
  character(len=80)  :: txt                        ! For read-in 
  integer, parameter :: NCOLS_MAX = 12           ! Max. no. cols. in "femis"
  integer            :: ncols, n, oldn             ! No. cols. in "femis", 
  integer, dimension(NEMIS) :: qc    ! index for sorting femis columns
  character(len=5), dimension(NCOLS_MAX)::  polltxt  ! to read line 1
  real, dimension(NCOLS_MAX)   :: e_f              ! factors read from femis


  e_fact(:,:,:) = 1.0            !/*** default value = 1 ***/


  call open_file(IO_EMIS,"r","femis.dat",needed=.false.)

  !**************************
  if ( ios == NO_FILE ) then
	ios = 0
	return !/** if no femis file, e_fact=1 as default **/ 
  endif
  if ( ios < 0  ) then
       errmsg = "GetEmis ios: femis.dat"
       call gc_abort(me,NPROC,errmsg)
  endif
  !**************************


  !/** read in header line, e.g. name sec sox nox voc. Pollutant
  !    names wil be searched against those  defined in moddef_ml   **/

  read(unit=IO_EMIS,fmt="(a80)") txt

  call wordsplit(txt,NCOLS_MAX,polltxt,ncols,ios)
	if(ios>0)return
  write(unit=6,fmt=*) "In femis, header is: ",  txt
  write(unit=6,fmt=*) "In femis, file has ", ncols, " columns (-2)"

  !/** we allow the femis file to give factors in any order, and
  !    for pollutants not needed, so we need to work out the indices
  !    for each column. Remember also that ncols includes the 1st
  !    2 columns (country_code and sector), which are not e_factors

  ncols = ncols - 2
  if ( ncols > NCOLS_MAX .or. ncols < 1) then
       call gc_abort(me,NPROC,"femis ncols")
  endif

  n = 0
  COLS: do ic=1,ncols
      oldn = n
      EMLOOP: do ie=1, NEMIS
                if ( polltxt(ic+2) == trim ( EMIS_NAME(ie) ) ) then
                     qc(ie) = ic
                     n = n + 1
                     write(unit=6,fmt=*) "In femis: ", polltxt(ic+2),  &
                              " assigned to ", ie, EMIS_NAME(ie)
                  exit EMLOOP
              end if
       end do EMLOOP ! ie
       if (oldn == n) write(unit=6,fmt=*) "femis: ",polltxt(ic+2)," NOT assigned"
  end do COLS   ! ic

  if ( n < NEMIS ) call gc_abort(me,NPROC,"too few femis items")
  
  n = 0

  READFILE: do  ! ************ read lines of femis ****************************** 

      read(unit=IO_EMIS,fmt=*,iostat=ios) inland, isec, (e_f(ic),ic=1,ncols)

      if ( ios <  0 ) exit READFILE                   ! End of file
      if ( ios >  0 ) then                             ! A different problem..
          print *, "ios: ", ios, "n  : ", n
          call gc_abort(me,NPROC,"Read error on femis")
      end if

      n = n + 1
      write(unit=6,fmt=*) "FEMIS READ", inland, isec, (e_f(ic),ic=1,ncols)

      if (inland == 0 ) then    ! Process all countries
          iland1 = 1 
          iland2 = NLAND
      else                       ! one country: inland
!find country number  corresponding to index as written in emisfile
         do iland1=1,NLAND
            if(Country(iland1)%index==inland)goto 544
         enddo
         if(me==0)write(*,*)'COUNTRY CODE NOT RECOGNIZED',inland
         iland1 = 0
         iland2 =-1
544      continue
         if(iland1/=0) iland2 = iland1
      end if

      if (isec == 0 ) then    ! All sectors
          isec1 = 1
          isec2 = NSECTORS
      elseif (isec==100) then    !hf scenario
          isec1 = 1
          isec2 = ANTROP_SECTORS
      else                       ! one sector: isec
          isec1 = isec
          isec2 = isec
      end if

      do ie = 1,NEMIS

          do iq = iland1,iland2
              do isec = isec1, isec2
                e_fact(isec,iq,ie) = e_fact(isec,iq,ie) * e_f( qc(ie) )
              end do !isec
          end do !iq

          if (DEBUG ) then
              write(unit=6,fmt=*) "IN NEMIS LOOP WE HAVE : ", ie, &
                                       qc(ie), e_f( qc(ie) )
              write(unit=6,fmt=*) "loops over ", isec1, isec2, iland1, iland2
          end if ! DEBUG
      end do !ie
          
  enddo READFILE ! Loop over femis

  close(IO_EMIS)

  write(unit=6,fmt=*) "In femis, read ", n, "records from femis."
  if ( DEBUG ) then
     write(unit=6,fmt=*) " For UK this gives: "
     write(unit=6,fmt="(6x, 10a12)") (EMIS_NAME(ie), ie=1,NEMIS)
     do isec = 1, 11
       write(unit=6,fmt="(i6, 10f12.6)") isec, (e_fact(isec,27,ie), ie=1,NEMIS)
     end do
  end if ! DEBUG
	ios = 0
  end subroutine femis
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        
!***********************************************************************
    subroutine EmisSplit()
!***********************************************************************
!**    DESCRIPTION:
!* Sets voc or pm speciation per source category for each country from input
!  files. Done for species where NEMISFRAC > 1
!
! Input files:
!    xxxsplit.defaults   (e.g. vocsplit.defaults, pm25split.defaults)
!    xxxsplit.special
!        where xxx can be voc or pm or whatever has NEMISFRAC > 1
!
! Output: (module public variable)
!    real emisfrac(NRCSPLIT,NSECTORS,NLAND)
!
! Language : F
! REVISION HISTORY:
!    Use My_Emis_ml splitting,            ds 26/1/01
!    Re-written for unified model,        ds 5/6/00
!    Adapted from Lagramgian ozone model, ds 9/2/99
!
!*******************************************************************

  !-- local
  integer ::  ie                  ! emission index in EMIS_NAME   (1..NEMIS)
  integer ::  isp                 ! emission index in 1..NEMIS_SPLIT
  integer ::  ifr0, ifr           ! index of split compound in emisfrac   

  integer, parameter :: NONREACTIVE = 1   ! No.columns of non-reactive species
                                          ! enforced for read-ins.

  ! for read-ins, dimension for max possible number of columns: 
  character(len=12), dimension(0:1, NRCSPLIT + NONREACTIVE ) :: intext
  real             ,      dimension(NRCSPLIT + NONREACTIVE ) :: tmp 
  integer           :: nsplit   ! No.columns data to be read
  real              :: sumtmp
  integer           :: iland,isec,i,n,nn

  logical  :: defaults    ! Set to true for defaults, false for specials
  integer  :: idef        ! Set to   0  for defaults, 1 for specials
  integer  :: iland1, iland2    ! loop variables over countries


  if ( NEMIS_SPLIT == 0 ) return  !/** for safety **/


  ifr0 =    1                 ! Starting index in  emisfrac array

  do isp = 1, NEMIS_SPLIT
    ie =   NEMIS_PLAIN + isp   ! Split species, index in 1..NEMIS

    nsplit =    EMIS_NSPLIT(isp) + NONREACTIVE 

    if ( isp > 1 ) ifr0 = ifr0 + EMIS_NSPLIT(isp-1) !start index of next species

    !/ Just in case ....
    if ( EMIS_NAME(ie) /= SPLIT_NAME(isp) )  &
          call gc_abort(me,NPROC,"Mis-match SPLIT")

    IDEF_LOOP: &
    do idef = 0, 1

       defaults = ( idef == 0 )

       if ( defaults ) then

          fname = trim( EMIS_NAME(ie) ) // "split.defaults"
          call open_file(IO_EMIS,"r",fname,needed=.true.)
          if (ios /= 0 ) &
             call gc_abort(me,NPROC,"ios error:split.defaults")

       else  ! Specials - overwrite the defaults if they exist.

          fname = trim( EMIS_NAME(ie) ) // "split.special"
          call open_file(IO_EMIS,"r",fname,needed=.false.)

          if ( ios == NO_FILE ) then   ! g12:.not. fexist ) then
              write(unit=6,fmt=*) "emis_split: no specials for",EMIS_NAME(ie)
		ios = 0
              exit IDEF_LOOP
          endif
       end if
       if ( DEBUG ) write(unit=6,fmt=*) "TTT split defaults=", defaults, fname
 
       !/ Read text line and speciation:
       !  the following lines expect one line of header text with the
       !  species names, followed by lines of the format:
       !    iland, isec, tmp1, tmp2.... tmpN+1, where the N+1'th
       !  column is for non-reactive species. These non-reactives are not used 
       !  in the rest of the program, but are required to check mass-balance.

       read(unit=IO_EMIS,fmt=*,iostat=ios) iland, isec ,(intext(idef,i), i=1, nsplit)
 
       if ( ios /= 0 ) then
           print *, "ERROR ios: ", ios, "values  : ", idef, EMIS_NAME(ie)
           call gc_abort(me,NPROC,"Read error on hearer, emis_split")
       end if

       write(unit=6,fmt="(a25,i3,/,(12a7))") "SPLIT species for idef=", &
                      idef, (intext(idef,i), i=1, nsplit)
       write(unit=6,fmt=*) "Will try to split ", EMIS_NSPLIT(isp) , " times"
           
       n = 0
       READ_DATA: &
       do 
          read(unit=IO_EMIS,fmt=*,iostat=ios) iland, isec, (tmp(i),i=1, nsplit)

           if ( ios <  0 ) exit READ_DATA     ! End of file
           if ( ios >  0 ) then               ! A different problem..
               print *, "ios: ", ios, "values  : ", idef, n,  EMIS_NAME(ie)
               call gc_abort(me,NPROC,"Read error on emis_split")
           end if
           n = n + 1

           !/... some checks:
           sumtmp = sum( tmp(1:nsplit) )
!hkpw           if ( ( sumtmp  >    100.01 .or. sumtmp   <   99.99   )  .or. &
!                ( defaults .and. iland /= 0                     )  .or. &
!                ( defaults .and. isec  /= n                     )  .or. &
!                ( .not. defaults .and. &
!                 any( intext(1,1:nsplit) /= intext(0,1:nsplit) ))      &
!                                                                  ) then
!               print *, "ERROR: emisfrac:", idef, iland, isec, sumtmp 
!               call gc_abort(me,NPROC,"emfracisec")
!           end if
          if ( ( sumtmp  >    100.01 .or. sumtmp   <   99.99   )  .or. &
                ( defaults .and. iland /= 0                     )  .or. &
                ( defaults .and. isec  /= n                     )  &
                                                                  ) then
               print *, "ERROR: emisfrac:", idef, iland, isec, sumtmp 
               call gc_abort(me,NPROC,"emfracisec")
           end if
           if (  .not. defaults ) then
              do nn=1,nsplit
                 if(intext(1,nn) /= intext(0,nn))then
                    write(*,*)"ERROR: emisfrac:", intext(1,nn), intext(0,nn)
                    call gc_abort(me,NPROC,"emfracisec")
                 endif
              enddo
           end if

           if ( defaults .or. iland == 0 ) then
             iland1 = 1
             iland2 = NLAND
           else  ! specials for one country
             iland1 = iland
             iland2 = iland
           end if
               
           do iland = iland1, iland2
             do i = 1, EMIS_NSPLIT(isp)
                ifr = ifr0 + i - 1    ! => index in emisfrac array

                !/** assign and convert from percent to fractions: **/

                emisfrac(ifr,isec,iland) = 0.01 * tmp(i)

                if ( DEBUG .and. iland == 27 ) then 
                    print "(a15,3i3,f10.4)", "TTT splitdef UK", isec, &
                                   ifr0, ifr, emisfrac(ifr,isec,iland)
                endif
             enddo ! i
           enddo ! iland

       enddo READ_DATA 
       close(IO_EMIS)

       if ( defaults .and. n  /=  NSECTORS ) &
               call gc_abort(me,NPROC,"emfracread1")

       write(unit=6,fmt=*) "Read ", n, " records from ",fname

     end do IDEF_LOOP 
   end do ! ie
  ios = 0

 end subroutine EmisSplit
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

end module EmisGet_ml
