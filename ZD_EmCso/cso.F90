!#######################################################################
!
! CSO - CAMS Satellite Operator
!
! Part of the CSO toolbox.
! See documentation included in toolbox for info.
!
! For usage, see "tutorial.F90".
!
!
!### macro's ###########################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!#######################################################################

module CSO

  use CSO_Comm
  use CSO_Logging
  use CSO_Domains
  use CSO_DateTimes
  use CSO_String
  use CSO_Rc
  use CSO_Listing
  use CSO_Sat
  use CSO_Tools
  use CSO_Grid
  use CSO_Profile
  use CSO_NcFile
  
  
  ! --- in/out ---------------------------------
  
  ! all public ...
  public
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter, private   ::  mname = 'CSO'


  ! --- types --------------------------------
  
  ! main object
  type :: T_CSO
    ! dummy ...
    logical               ::  enabled
  contains
    procedure :: Init           =>  CSO_Init
    procedure :: Done           =>  CSO_Done
    procedure :: SetLogging     =>  CSO_SetLogging
  end type T_CSO
  

  
contains


  ! ====================================================================
  ! ===
  ! === module init/done
  ! ===
  ! ====================================================================


#ifdef _MPI
  subroutine CSO_Init( self, status, comm, icomm )
#else
  subroutine CSO_Init( self, status )
#endif

#ifdef _MPI
    use MPI_F08, only : MPI_Comm
#endif
    
    use CSO_Comm   , only : csoc
    !use CSO_Logging, only : lgr
  
    ! --- in/out ---------------------------------
    
    class(T_CSO), intent(out)         ::  self
    integer, intent(out)              ::  status
    
#ifdef _MPI
    type(MPI_Comm), intent(in), optional    ::  comm
    integer, intent(in), optional           ::  icomm
#endif

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_Init'
    
    ! --- local ----------------------------------
    
#ifdef _MPI
    type(MPI_Comm)   ::  xcomm
#endif
    
    ! --- begin ----------------------------------
    
    ! dummy:
    self%enabled = .true.

#ifdef _MPI
    !~ F08 communicator provided?
    if ( present(comm) ) then
      ! init communication:
      call csoc%Init( status, comm=comm )
      IF_NOT_OK_RETURN(status=1)
    !~ F77 integer communicator provided?
    else if ( present(icomm) ) then
      ! copy integer value into F08 type:
      xcomm%MPI_VAL = icomm
      ! init communication:
      call csoc%Init( status, comm=xcomm )
      IF_NOT_OK_RETURN(status=1)
    !~ serial:
    else
      ! init for serial run:
      call csoc%Init( status )
      IF_NOT_OK_RETURN(status=1)
    end if ! with_comm
#else
    ! init for serial run:
    call csoc%Init( status )
    IF_NOT_OK_RETURN(status=1)
#endif

    ! init logger:
    call lgr%Init( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine CSO_Init


  ! ***


  subroutine CSO_Done( self, status )
  
    use CSO_Comm   , only : csoc
    !use CSO_Logging, only : lgr
    
    ! --- in/out ---------------------------------
    
    class(T_CSO), intent(inout)         ::  self
    integer, intent(out)                ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/CSO_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! done with logger:
    call lgr%Done( status )
    if (status/=0) stop
  
    ! done with communication:
    call csoc%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine CSO_Done
  
  
  ! **


  subroutine CSO_SetLogging( self, status, &
                               file, unit, root_only )
    
    use CSO_Comm   , only : csoc
    !use CSO_Logging, only : lgr
  
    ! --- in/out ---------------------------------
    
    class(T_CSO), intent(in)                ::  self
    integer, intent(out)                    ::  status
    
    character(len=*), optional              ::  file
    integer, optional                       ::  unit
    logical, optional                       ::  root_only

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/CSO_SetLogging'
    
    ! --- local ----------------------------------
    
    logical                 ::  lgr_apply
    character(len=1024)     ::  lgr_file
    
    ! --- begin ----------------------------------
    
    ! set unit?
    if ( present(unit) ) then
      ! set unit:
      call lgr%Set( status, unit=unit )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! message from root only?
    if ( present(root_only) ) then
      ! set flag:
      if ( root_only ) then
        lgr_apply = csoc%id == 0
      else
        lgr_apply = .true.
      end if
      ! set flag:
      call lgr%Set( status, apply=lgr_apply )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! write messages to file?
    if ( present(file) ) then
      ! redirect:
      call lgr%LogFile( status, file=file )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine CSO_SetLogging


end module CSO
