module Io_ml

!_____________________________________________________________________________
! -- sets global file io numbers and provides routines to check/open files
!_____________________________________________________________________________
implicit none

! -- subroutines in this module:

public :: check_file   !  checks that file exists and stops if required
public :: open_file    !  checks that file exists and opens if required
public :: wordsplit    !  Splits input text into words

logical, public :: fexist                      ! true if file exists 
integer, public, parameter :: NO_FILE = 777    ! code for non-existing file
integer, public, save :: ios                   ! i/o error status number

! The idea is to keep all unit numbers used for input and output
! stored here. So that new programmers can quickly see which numbers are
! in use or not.
!
! Assign unit number, e.g. io_xxx, here, and use read(io_xxx,*) in
! main program.
!
! (c) - file opened AND closed in subroutine
! (o) - file remains open in subroutine
! ...........................................................................

  integer, parameter, public  :: &
    IO_LOG      = 7   &! General output log (o)
   ,IO_SITES    = 8   &! sites module, first for input(c)
   ,IO_COMPA    = 9   &! o3mod.f(c)-input comp_array
   ,IO_MYTIM    = 20  &! o3mod.f(c)-output mytim.out 
   ,IO_RES      = 25  &! o3mod,massbud(o) - ! out eulmod.res
   ,IO_CHECK    = 26  &! readpar(o?) - check output
   ,IO_TMP      = 27   ! General IO number (files *must* be type (c))


! tmp units 30-39 reserved for sondes. Will be set of
!     one unit number io_sonde later

  integer, parameter, public  :: &
    IO_SONDES   = 30  ! siteswrt_ml(o)  for output of sonde data


  integer, parameter, public  :: &
    IO_FORES    = 49  &! rforest.f(c)-read land use %
   ,IO_AIRN     = 49  &! airnox.f(c) - read aircr. em.
   ,IO_LIGHT    = 49  &! lightning.f(c) - read lightning. emiss.
   ,IO_JOST     = 49  &! newjostinit(c) - read  global mixing ratios
   ,IO_GLOBBC   = 49  &!u3 - read  global mixing ratios e.g. Logan
   ,IO_GLOBBC2  = 91  &!u3 - read  global mixing ratios e.g. h2o2
   ,IO_INFIELD  = 50  &! infield.F(c) -reads fil000xx  
   ,IO_ROUGH    = 52  &! inpar.f -reads roughn. class  
   ,IO_SNOW     = 53  &! newmonth(c): for snow
   ,IO_DJ       = 55  &! readdiss.f(c) - inp. solar r.
   ,IO_AIRCR    = 66  &! phyche.f(c) - write aircraft conc.
   ,IO_OUT      = 80  &! (c)write outday etc.
   ,IO_UKDEP    = 81  &! (o)write fluxes, etc.
   ,IO_STAB     = 82  &! (o)write fluxes, etc.
   ,IO_EMIS     = 84  &! Used for femis , emis_split(c)
   ,IO_TIMEFACS = 85   ! Used for monthly

  integer, parameter, public  :: &
    IO_DMS      = 90  &! newmonth(c): for DMS 
   ,IO_AOT      = 118 &! o3mod(c): for 6-monthly values of AOT, mean
   ,IO_HOURLY   = 119 &! hourly_out(o)
!hf
   ,IO_VOLC     = 54  &
   ,IO_NEST     = 88   


contains

!===========================================================>

subroutine check_file(fname,fexist,needed)
  !**************************************************************
  ! Checks for the existence of a file. If the file is
  ! specified as "needed", and missing, an error message is
  ! printed and the run is stopped.

  character (len=*),   intent(in)  :: fname    ! file name
  logical,             intent(in)  :: needed   ! see below
  logical,             intent(out) :: fexist   ! file exists

  inquire(file=fname,exist=fexist)

  write(unit=6,fmt=*) "check_file::: ", fname
  if ( .not. fexist .and. .not. needed ) then
     write(unit=6,fmt=*) "not needed, skipping....."
     ios = 0

  else if ( .not. fexist .and. needed ) then
     ios = -1
     print *, "ERROR: Missing!!! in check-file"

  else
     write(unit=6,fmt=*) "ok. File exists"
  end if
end subroutine check_file

!===========================================================>

subroutine open_file(io_num,mode,fname,needed,skip)
  !**************************************************************
  ! Checks for the existence of a file and opens if present. If the 
  ! file is specified as "needed", and missing, an error message is
  ! printed and the run is stopped.

  integer,             intent(in)  :: io_num   ! i/o number
  character (len=*),   intent(in)  :: mode     ! "r" for read, "w" for write
  character (len=*),   intent(in)  :: fname    ! file name
  logical,             intent(in)  :: needed   ! see below
  integer, optional,   intent(in)  :: skip     ! No. text lines to be skipped


  integer :: i  ! local loop counter

  ios = 0   ! Start with  assumed ok status
  inquire(file=fname,exist=fexist)

  if ( .not. fexist ) then

     if ( mode == "r" ) then
         if ( needed ) then
             print *, "OPEN_FILE ::: missing file is :: ", fname
             ios = -1
         else
            ios = NO_FILE
         end if
     else
        write(unit=6,fmt=*) "File created: ", io_num, fname
        open(io_num,file=fname,status="new",action="write",iostat=ios)
     end if
  else
     write(unit=6,fmt=*) "File exists: ", fname
     select case (mode) 
     case ("r")
        open(io_num,file=fname,status="old",action="read",iostat=ios)
        ! *** skip header lines if requested ****
        if ( present( skip ) ) then ! Read (skip) some lines at start of file
           do i = 1, skip
               read(unit=io_num,fmt=*)
           end do
        end if ! skip
        ! ***  
     case ("w")
        write(unit=6,fmt=*) "File created: ", fname
        open(io_num,file=fname,status="new",action="write",iostat=ios)
     case default
        print *, "OPEN FILE: Incorrect mode: ", mode
        ios = -1
     end select

  end if
end subroutine open_file

!===========================================================>

  subroutine wordsplit(text,nword_max,wordarray,nwords,errcode)
  !**************************************************************
  !   Subroutine takes in a character string and splits it into
  !   a word-array, of length nwords        
  !**************************************************************

 !-- arguments
  character(len=*), intent(in) ::  text       ! to be split
  integer,          intent(in) ::  nword_max  ! Max. no. words expected

  character(len=*), dimension(:), intent(out) :: wordarray
  integer,          intent(out) :: nwords      ! No. words found
  integer,          intent(out) :: errcode      ! No. words found

  !-- local
  logical   :: wasinword   ! true if we are in or have just left a word
  integer   :: i, is, iw
  character(len=1) ::  c

  errcode = 0
  wasinword = .false.   !To be safe, with spaces at start of line (spotted-hf)
  is = 0 ! string index
  iw = 1 ! Word index
  wordarray(1) = ""

  do i = 1, len_trim(text)
      c = text(i:i)
      if( c /= " " ) then
          is = is + 1
          wordarray(iw)(is:is) = c
          wasinword = .true.
      else 
         if ( wasinword ) then
             iw = iw + 1
             wordarray(iw) = ""
             wasinword = .false.
             is = 0
         endif
      endif
  enddo
  nwords = iw
  if (  nwords >= nword_max ) then
	errcode = 2
      print *, "ERROR in WORDSPLIT : Problem at ", text
      print *,"Too many words"
!      call stop_all("Too many words")
  endif

 end subroutine wordsplit
end module Io_ml
