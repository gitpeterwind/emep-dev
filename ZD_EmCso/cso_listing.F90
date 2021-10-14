!#################################################################
!
! Listing with CSO files.
!
! Example content:
!
!   filename;start_time;end_time;orbit
!   2018/06/S5p_RPRO_NO2_03272.nc;2018-06-01T01:32:46.673000000;2018-06-01T01:36:12.948000000;03272
!   2018/06/S5p_RPRO_NO2_03273.nc;2018-06-01T03:12:53.649000000;2018-06-01T03:17:43.082000000;03273
!   2018/06/S5p_RPRO_NO2_03274.nc;2018-06-01T04:52:43.586000000;2018-06-01T04:59:12.377000000;03274

!  
!#################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!#################################################################

module CSO_Listing

  use CSO_Logging  , only : csol, csoPr, csoErr
  use CSO_DateTimes, only : T_CSO_DateTime

  implicit none

  ! --- in/out ---------------------

  private
  
  public  ::  T_CSO_Listing


  ! --- const ---------------------------------
  
  character(len=*), parameter  ::  mname = 'CSO_Listing'


  ! --- types ---------------------------------
  
  type T_CSO_Record
    ! filename:
    character(len=256)      ::  filename
    ! timerange:
    type(T_CSO_DateTime)    ::  t1, t2
    ! average:
    type(T_CSO_DateTime)    ::  taver
  contains
    procedure ::  Init        =>  Record_Init
    procedure ::  Done        =>  Record_Done
    procedure ::  Match       =>  Record_Match
  end type T_CSO_Record
  
  ! *
  
  type T_CSO_Listing
    ! filename:
    character(len=1024)               ::  filename
    ! directory part:
    character(len=1024)               ::  dirname
    ! formatting:
    character(len=1)                  ::  sep
    ! number of records:
    integer                           ::  nrec
    ! storage:
    type(T_CSO_Record), allocatable   ::  rec(:)  ! (n)
  contains
    procedure ::  Init        =>  Listing_Init
    procedure ::  Done        =>  Listing_Done
    procedure ::  Show        =>  Listing_Show
    procedure ::  SearchFile  =>  Listing_SearchFile
  end type T_CSO_Listing



contains

 

  ! ==============================================================
  ! ===
  ! === listing record
  ! ===
  ! ==============================================================


  subroutine Record_Init( self, filename, t1, t2, status )
  
    use CSO_DateTimes, only : operator(+), operator(-), operator(*)

    ! --- in/out ------------------------

    class(T_CSO_Record), intent(out)          ::  self
    character(len=*), intent(in)              ::  filename
    type(T_CSO_DateTime), intent(in)          ::  t1, t2
    integer, intent(out)                      ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Record_Init'
    
    ! --- local ----------------------------

    ! --- begin ----------------------------
    
    ! store:
    self%filename = trim(filename)
    self%t1 = t1
    self%t2 = t2
    ! fill average:
    self%taver = self%t1 + 0.5*( self%t2 - self%t1 )

    ! ok
    status = 0

  end subroutine Record_Init


  ! ***


  subroutine Record_Done( self, status )

    ! --- in/out -----------------

    class(T_CSO_Record), intent(inout)   ::  self
    integer, intent(out)                 ::  status

    ! --- const ----------------------
    
    character(len=*), parameter  ::  rname = mname//'/Record_Done'
    
    ! --- begin ------------------------

    ! ok
    status = 0

  end subroutine Record_Done


  ! ***

  !
  ! Try if record can be assigned to target interval [t1,t2].
  !
  ! How the record is matched depends on the timekey:
  !   'aver'  : record time average in (t1,t2]
  !
  ! Return value:
  !   match :  .true. if matched
  !   status: 0 if ok, >0 error
  !

  subroutine Record_Match( self, t1, t2, timekey, match, status )

    use CSO_DateTimes, only : Pretty, operator(<), operator(<=)

    ! --- in/out -----------------

    class(T_CSO_Record), intent(in)       ::  self
    type(T_CSO_Datetime), intent(in)      ::  t1, t2
    character(len=*), intent(in)          ::  timekey
    logical, intent(out)                  ::  match
    integer, intent(out)                  ::  status

    ! --- const ----------------------
    
    character(len=*), parameter  ::  rname = mname//'/Record_Match'
    
    ! --- local ------------------------
    
    ! --- begin ------------------------
    
    ! switch:
    select case ( trim(timekey) )
    
      !~ average?
      case ( 'aver' )
        ! average should be in interval:
        match = (t1 < self%taver) .and. (self%taver <= t2)
        
      !~ unknown ...
      case default
        write (csol,'("unsupported timekey `",a,"`")') trim(timekey); call csoErr
        TRACEBACK; status=1; return
        
    end select

    ! ok
    status = 0

  end subroutine Record_Match





  ! ==============================================================
  ! ===
  ! === listing file
  ! ===
  ! ==============================================================


  subroutine Listing_Init( self, filename, status )

    use CSO_File     , only : CSO_GetFU
    use CSO_File     , only : CSO_GetDirname
    use CSO_String   , only : CSO_ReadFromLine
    use CSO_DateTimes, only : CSO_ReadFromLine
    
    ! --- in/out ------------------------

    class(T_CSO_Listing), intent(out)         ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status

    ! --- const ---------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Listing_Init'
    
    ! --- local ----------------------------

    logical               ::  exist
    integer               ::  fu
    character(len=1024)   ::  line
    integer               ::  irec
    character(len=1024)   ::  rec_filename
    type(T_CSO_DateTime)  ::  rec_t1, rec_t2

    ! --- begin ----------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! directory part:
    call CSO_GetDirname( self%filename, self%dirname, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! formatting:
    self%sep = ';'

    ! file exist ?
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (csol,'("listing file not found:")'); call csoErr
      write (csol,'("  file name : ",a)') trim(filename); call csoErr
      write (csol,'("in ",a)') rname; call csoErr; status=1; return
    end if

    ! select free file unit:
    Call CSO_GetFU( fu, status )
    IF_NOT_OK_RETURN(status=1)
    
    !~ count records

    ! open file:
    open( unit=fu, file=trim(filename), iostat=status, form='formatted' )
    if ( status /= 0 ) then
      write (csol,'("from opening : ",a)') trim(filename); call csoErr
      TRACEBACK; status=1; return
    end if

    ! init counter:
    self%nrec = 0
    ! loop over records:
    do
      ! read next line:
      read (fu,*,iostat=status) line
      if ( status < 0 ) then
        exit ! eof
      else if ( status > 0 ) then
        write (csol,'("from reading line ",i0," from file: ",a)') self%nrec+1,trim(self%filename); call csoErr
        TRACEBACK; status=1; return
      end if
      ! increase counter:
      self%nrec = self%nrec + 1
    end do  ! records

    ! close file:
    close( unit=fu, iostat=status )
    if ( status /= 0 ) then
      write (csol,'("from closing file:")'); call csoErr
      write (csol,'("  ",a)') trim(self%filename); call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! check ..
    if ( self%nrec == 0 ) then
      write (csol,'("found empty listing file: ",a)') trim(self%filename); call csoErr
      TRACEBACK; status=1; return
    end if

    ! one less for header:
    self%nrec = self%nrec - 1
    
    ! ~ records
    
    ! any?
    if ( self%nrec > 0 ) then
    
      ! storage:
      allocate( self%rec(self%nrec), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! open file:
      open( unit=fu, file=trim(filename), iostat=status, form='formatted' )
      if ( status /= 0 ) then
        write (csol,'("from file open :")'); call csoErr
        write (csol,'("  file name : ",a)') trim(filename); call csoErr
        TRACEBACK; status=1; return
      end if

      ! loop over records, including header:
      do irec = 0, self%nrec
        ! read next line:
        read (fu,'(a)',iostat=status) line
        if ( status < 0 ) then
          exit ! eof
        else if ( status > 0 ) then
          write (csol,'("from reading line ",i0," from file: ",a)') self%nrec+1,trim(self%filename); call csoErr
          TRACEBACK; status=1; return
        end if
        
        ! header line?
        if ( irec == 0 ) then
          ! check ...
          if ( line(1:8) /= 'filename' ) then
            write (csol,'("unexpected header:")'); call csoErr
            write (csol,'("  ",a)') trim(line); call csoErr
            write (csol,'("  file name : ",a)') trim(filename); call csoErr
            TRACEBACK; status=1; return
          end if
          ! skip:
          cycle
        end if
        
        ! extract elements:
        call CSO_ReadFromLine( line, rec_filename, status, sep=self%sep )
        IF_NOT_OK_RETURN(status=1)
        call CSO_ReadFromLine( line, rec_t1, status, sep=self%sep )
        IF_NOT_OK_RETURN(status=1)
        call CSO_ReadFromLine( line, rec_t2, status, sep=self%sep )
        IF_NOT_OK_RETURN(status=1)

        ! store:
        call self%rec(irec)%Init( rec_filename, rec_t1, rec_t2, status )
        IF_NOT_OK_RETURN(status=1)
        
      end do  ! records

      ! close file:
      close( unit=fu, iostat=status )
      if ( status /= 0 ) then
        write (csol,'("from closing file:")'); call csoErr
        write (csol,'("  ",a)') trim(self%filename); call csoErr
        TRACEBACK; status=1; return
      end if
      
    end if ! n>0

    ! ok
    status = 0

  end subroutine Listing_Init


  ! ***


  subroutine Listing_Done( self, status )

    ! --- in/out -----------------

    class(T_CSO_Listing), intent(inout)   ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------
    
    character(len=*), parameter  ::  rname = mname//'/Listing_Done'
    
    ! --- local ------------------------
    
    integer     ::  irec
    
    ! --- begin ------------------------
    
    ! any?
    if ( self%nrec > 0 ) then
      ! loop:
      do irec = 1, self%nrec
        ! clear:
        call self%rec(irec)%Done( status )
        IF_NOT_OK_RETURN(status=1)
      end do
      ! clear:
      deallocate( self%rec, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! empty:
      self%nrec = 0
    end if ! n>0
    
    ! ok
    status = 0

  end subroutine Listing_Done


  ! ***


  subroutine Listing_Show( self, status )

    use CSO_DateTimes, only : Pretty

    ! --- in/out -----------------

    class(T_CSO_Listing), intent(in)      ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------
    
    character(len=*), parameter  ::  rname = mname//'/Listing_Show'
    
    ! --- local ------------------------
    
    integer     ::  irec
    
    ! --- begin ------------------------
    
    ! info ...
    write (csol,'("---------------------------------------------------------------------------------------")'); call csoPr
    write (csol,'("listing file      : ",a)') trim(self%filename); call csoPr
    write (csol,'("number of records : ",i0)') self%nrec; call csoPr
    do irec = 1, self%nrec
      write (csol,'(i6,"  [ ",a,", ",a," ]  ",a)') irec, &
                trim(Pretty(self%rec(irec)%t1)), trim(Pretty(self%rec(irec)%t2)), &
                trim(self%rec(irec)%filename); call csoPr
    end do
    write (csol,'("---------------------------------------------------------------------------------------")'); call csoPr
    
    ! ok
    status = 0

  end subroutine Listing_Show


  ! ***

  !
  ! Search record that assigned to time range [t1,t2] .
  !
  ! How the record is matched depends on the timekey:
  !   'aver'  : record time average in [t1,t2]
  !
  ! Return value:
  !   filename: full path to file, empty if not file found
  !   status: 0 if ok, >0 error
  !

  subroutine Listing_SearchFile( self, t1, t2, timekey, filename, status )

    use CSO_DateTimes, only : Pretty

    ! --- in/out -----------------

    class(T_CSO_Listing), intent(in)      ::  self
    type(T_CSO_Datetime), intent(in)      ::  t1, t2
    character(len=*), intent(in)          ::  timekey
    character(len=*), intent(out)         ::  filename
    integer, intent(out)                  ::  status

    ! --- const ----------------------
    
    character(len=*), parameter  ::  rname = mname//'/Listing_SearchFile'
    
    ! --- local ------------------------
    
    integer     ::  irec
    logical     ::  match
    
    ! --- begin ------------------------
    
    ! init result:
    filename = ''
    ! loop until first match:
    do irec = 1, self%nrec
      ! match?
      call self%rec(irec)%Match( t1,t2, timekey, match, status )
      IF_NOT_OK_RETURN(status=1)
      if ( match ) then
        ! match, set result:
        write (filename,'(a,"/",a)') trim(self%dirname), trim(self%rec(irec)%filename)
        ! leave:
        exit
      end if ! match?
    end do ! records
    
    ! ok
    status = 0

  end subroutine Listing_SearchFile



end module CSO_Listing
