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

module EMEP_NMC_Driver

  use GO                      , only : gol, goPr, goErr
  use EMEP_NMC_Driver_Eps     , only : T_Driver_Eps
  use EMEP_NMC_Driver_Eps_Eval, only : T_Driver_Eps_Eval
  use EMEP_NMC_Driver_Eta     , only : T_Driver_Eta
  use EMEP_NMC_Driver_Eta_f   , only : T_Driver_Eta_f
  use EMEP_NMC_Driver_C_f     , only : T_Driver_C_f
  use EMEP_NMC_Driver_D_f     , only : T_Driver_D_f
  use EMEP_NMC_Driver_B_f     , only : T_Driver_B_f
  use EMEP_NMC_Driver_XLX     , only : T_Driver_XLX
  use EMEP_NMC_Driver_BB      , only : T_Driver_BB
  use EMEP_NMC_Driver_BB_Eval , only : T_Driver_BB_Eval

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  ParseArguments
  public  ::  T_Driver_Eps
  public  ::  T_Driver_Eps_Eval
  public  ::  T_Driver_Eta
  public  ::  T_Driver_Eta_f
  public  ::  T_Driver_C_f
  public  ::  T_Driver_D_f
  public  ::  T_Driver_B_f
  public  ::  T_Driver_XLX
  public  ::  T_Driver_BB
  public  ::  T_Driver_BB_Eval
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver'
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** arguments
  ! ***
  ! ********************************************************************


  !
  ! Return status:
  !     -1 : usage displayed
  !      0 : ok
  !  other : error
  !
  
  subroutine ParseArguments( rcfile, status )
  
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    character(len=*), intent(out)    ::  rcfile
    integer, intent(out)             ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/ParseArguments'
    
    ! assumed target name of executable fromed:
    character(len=*), parameter  ::  xname = 'emep_nmc.x'
    
    ! --- local ----------------------------------
    
    integer                         ::  narg, iarg
    integer                         ::  arglength
    character(len=1024)             ::  argvalue
    integer                         ::  k
    
    ! --- begin ----------------------------------
    
    ! dummy values at start:
    rcfile = ''
    
    ! count:
    narg = Command_Argument_Count()

    ! loop:
    do iarg = 1, narg
      ! get length:
      call Get_Command_Argument( iarg, length=arglength, status=status )
      if ( status /= 0 ) then
        write (gol,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength; call goErr
        TRACEBACK; stop 1
      end if
      ! get value: 
      call Get_Command_Argument( iarg, value=argvalue, status=status )
      if ( status /= 0 ) then
        write (gol,'("non-zero status ",i6," from getting lenth of argument ",i6)') status, arglength; call goErr
        TRACEBACK; stop 1
      end if
      ! switch:
      select case ( trim(argvalue) )
        !~ help:
        case ( '-h', '--help' )
          write (gol,'("")'); call goPr
          write (gol,'("Usage:")'); call goPr
          write (gol,'("  ",a," rcfile")') xname; call goPr
          write (gol,'("  ",a," -h|--help")') xname; call goPr
          write (gol,'("")'); call goPr
          status = -1; return
        !~ values:
        case default
          ! switch:
          if ( len_trim(rcfile) == 0 ) then
            rcfile = trim(argvalue)
          else
            write (*,'("unsupported argument `",a,"`")') trim(argvalue)
            TRACEBACK; stop 1
          end if
      end select
    end do   ! arguments

    ! check ...
    if ( len_trim(rcfile) == 0 ) then
      write (*,'("no input rcfile argument specified")')
      TRACEBACK; stop 1
    end if
    
    ! ok
    status = 0
  
  end subroutine ParseArguments


end module EMEP_NMC_Driver

