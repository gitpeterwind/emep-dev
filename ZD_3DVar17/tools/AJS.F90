!#######################################################################
!
! AJS tools
!
!#######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!#######################################################################

module AJS

  use GO             , only : gol, goPr, goErr

  implicit none
  
  ! --- in/out -----------------------------------------
  
  private
  
  public  ::  AJS_Init, AJS_Done

  
  ! --- const -----------------------------------------

  character(len=*), parameter  ::  mname = 'AJS'

contains


  !-----------------------------------------------------------------------
  ! module init/done
  !-----------------------------------------------------------------------

  
  subroutine AJS_Init( status )

    use GO               , only : GO_Init
    use GO               , only : TrcFile, Init, Done, ReadRc
    use GO               , only : goGetFU, goStdErr
    use GO               , only : GO_Par_Setup, me
    use GO               , only : GO_Print_Set, GO_Print_Logfile
    use ModelConstants_ml, only : masterProc
    use MPI_Groups_ml    , only : MPI_COMM_CALC
    use DA_ml            , only : DEBUG_DA_1STEP

    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/AJS_Init'
    
    ! settings used by 3D-var:
    character(len=*), parameter   ::  rcfile = 'emo.rc'

    ! --- local -----------------------------
    
    character(len=1024)     ::  fname
    character(len=32)       ::  key
    type(TRcFile)           ::  rcF
    logical                 ::  flag
    character(len=1024)     ::  mefile
    
    ! --- begin -----------------------------
    
    ! check ...
#ifndef _MPI
    write (gol,'("MPI code should be enabled, define _MPI macro!")'); call goErr
    TRACEBACK; status=1; return
#endif
  
    ! initialize GO tools:
    call GO_Init( status )
    if ( status /= 0 ) then
      write (*,'("in ",a," (line",i5,")")') __FILE__, __LINE__
      stop
    end if
    ! from now on, the gol/goPr/goPr logging could be used ...
  
    ! setup parallel tools in GO modules:
    call GO_Par_Setup( MPI_COMM_CALC, status )
    IF_NOT_OK_RETURN(status=1)
  
    ! extra settings:
    call Init( rcF, rcfile, status )
    IF_NOT_OK_RETURN(status=1)

    ! log all processes ?
    call ReadRc( rcF, 'emo.log.all', flag, status )
    IF_NOT_OK_RETURN(status=1)
    ! set flag if processes should be logged:
    call GO_Print_Set( status, apply=(MasterProc .or. flag) )
    IF_NOT_OK_RETURN(status=1)

    ! per-process log files:
    write (mefile,'("go-pe",i3.3,".log")') me
    ! write gol messages to log file, also echo to stdout from root:
    call GO_Print_LogFile( status, file=trim(mefile), echo=MasterProc )
    IF_NOT_OK_RETURN(status=1)
    
    ! testing ...
    write (gol,'(a,": test output message ...")') rname; call goPr
    write (gol,'(a,": test error  message ...")') rname; call goErr

    ! single analysis step only?
    call ReadRc( rcF, 'test.da_1step', DEBUG_DA_1STEP, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! done with settings:
    call Done( rcF, status )
    IF_NOT_OK_RETURN(status=1)
  
    ! ok
    status = 0
    
  end subroutine AJS_Init
  
  
  ! ***
  
  
  subroutine AJS_Done( status )

    use ModelConstants_ml, only : masterProc
    use GO               , only : GO_Timer_Post
    use GO               , only : GO_Done

    ! --- in/out ----------------------------
    
    integer, intent(out)           ::  status
    
    ! --- const ----------------------------
    
    character(len=*), parameter  ::  rname = mname//'/AJS_Done'
    
    ! --- begin -----------------------------
  
    ! show and write time profile, sufficient on root:
    if ( MasterProc ) then
      ! done with timing, write profile, and echo to log file:
      call GO_Timer_Post( status, file='timing.prf', verbose=.true. )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! done with GO modules:
    call GO_Done( status )
    if ( status /= 0 ) then
      write (*,'("in ",a," (line",i5,")")') __FILE__, __LINE__
      stop
    end if

    ! ok
    status = 0
    
  end subroutine AJS_Done


end module AJS
