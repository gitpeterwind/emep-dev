!######################################################################
!
! EMEP DA NMC - driver routines
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!######################################################################

module EMEP_NMC_Common

  use GO                  , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  nrun

  public  ::  maxvar
  public  ::  T_ModelVars

  public  ::  maxhours

  public  ::  maxpoint
  public  ::  T_Points
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Common'
  
  ! number of model runs, probably 2 for NMC test:
  integer, parameter  ::  nrun = 2

  ! assumed maximum number of variables:
  integer, parameter  ::  maxvar = 10

  ! assumed maximum number of times within days:
  integer, parameter  ::  maxhours = 24

  ! assumed maximum number of point:
  integer, parameter  ::  maxpoint = 10


  ! --- types ----------------------------------------
  
  type T_ModelVar
    character(len=32)           ::  name
    character(len=32)           ::  units
  end type T_ModelVar
  
  type T_ModelVars
    integer                           ::  n
    type(T_ModelVar), allocatable     ::  value(:)
  contains
    procedure   ::  Init        => ModelVars_Init
    procedure   ::  Done        => ModelVars_Done
  end type T_ModelVars
  
  ! *
  
  type T_Point
    character(len=32)           ::  name
    real                        ::  lon, lat   ! deg
    integer                     ::  ilon, ilat, ilev
  end type T_Point
  
  type T_Points
    integer                     ::  n
    type(T_Point), allocatable  ::  value(:)
  contains
    procedure   ::  Init        => Points_Init
    procedure   ::  Done        => Points_Done
  end type T_Points
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** ModelVars
  ! ***
  ! ********************************************************************


  subroutine ModelVars_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine, goSplitString
  
    ! --- in/out ---------------------------------
    
    class(T_ModelVars), intent(out)           ::  self
    type(TrcFile), intent(in)              ::  rcF
    integer, intent(out)                   ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelVars_Init'
    
    ! --- local ----------------------------------
    
    character(len=1024)   ::  line
    character(len=32)     ::  vnames(maxvar)
    integer               ::  ivar
    
    ! --- begin ----------------------------------
    
    ! variable names:
    call ReadRc( rcF, 'nmc.variables', line, status )
    IF_NOT_OK_RETURN(status=1)
    ! split:
    call goSplitString( trim(line), self%n, vnames, status )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%value(self%n) )
    ! fill names:
    do ivar = 1, self%n
      self%value(ivar)%name = trim(vnames(ivar))
    end do
    
    ! ok
    status = 0
  
  end subroutine ModelVars_Init


  ! ***
  

  subroutine ModelVars_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ModelVars), intent(inout)        ::  self
    integer, intent(out)                     ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelVars_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    deallocate( self%value )
    
    ! ok
    status = 0
  
  end subroutine ModelVars_Done


  ! ********************************************************************
  ! ***
  ! *** Points
  ! ***
  ! ********************************************************************


  subroutine Points_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine, goSplitString
  
    ! --- in/out ---------------------------------
    
    class(T_Points), intent(out)           ::  self
    type(TrcFile), intent(in)              ::  rcF
    integer, intent(out)                   ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Points_Init'
    
    ! --- local ----------------------------------
    
    character(len=1024)   ::  line
    character(len=32)     ::  pnames(maxpoint)
    integer               ::  ipoint
    
    ! --- begin ----------------------------------
    
    ! point names:
    call ReadRc( rcF, 'nmc.points', line, status )
    IF_NOT_OK_RETURN(status=1)
    ! split:
    call goSplitString( trim(line), self%n, pnames, status )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%value(self%n) )
    ! loop:
    do ipoint = 1, self%n
      ! store:
      self%value(ipoint)%name = trim(pnames(ipoint))
      ! location:
      call ReadRc( rcF, 'nmc.point.'//trim(pnames(ipoint))//'.lon', self%value(ipoint)%lon, status )
      IF_NOT_OK_RETURN(status=1)
      call ReadRc( rcF, 'nmc.point.'//trim(pnames(ipoint))//'.lat', self%value(ipoint)%lat, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! ok
    status = 0
  
  end subroutine Points_Init


  ! ***
  

  subroutine Points_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Points), intent(inout)           ::  self
    integer, intent(out)                     ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Points_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear point info:
    deallocate( self%value )
    
    ! ok
    status = 0
  
  end subroutine Points_Done
  

end module EMEP_NMC_Common

