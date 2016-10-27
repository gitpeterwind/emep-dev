!#######################################################################
!
! Postprocessor for EMEP Background error covariance using NMC method.
!
! DESCRIPTION
!
!   Input is formed by output of two runs driven by different meteo.
!   This programs:
!    - reads the output files
!    - computes the difference
!    - computes statistics:
!      - sample mean
!      - sample variance
!      - selected covariance patterns
!    - perform FFT
!    - etc
!
! MACRO'S
!
!   The code uses the following fpp macro's:
!
!   _MPI   :   Enable MPI-specific code. This variable is named similar
!              to the OpenMP macro '_OPENMP' that is defined automatically
!              if OpenMP is enabled; for MPI such is not present by default,
!              but the user needs to define it as compile flag, probably:
!                mpif90 ... -D_MPI ..
!
!#######################################################################
!
#define TRACEBACK write (gol,'("ERROR in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; stop 1; end if
!
#define IF_MPI_NOT_OK_RETURN(action) if (status/=MPI_SUCCESS) then; TRACEBACK; action; stop 1; end if
!
!#######################################################################

program EMEP_NMC

#ifdef _MPI
  use MPI, only : MPI_COMM_WORLD
  use MPI, only : MPI_SUCCESS
  use MPI, only : MPI_DOUBLE
  use MPI, only : MPI_Init, MPI_Finalize
  use MPI, only : MPI_Comm_Size, MPI_Comm_Rank
#endif

  use GO, only : gol, goPr, goErr
  use GO, only : GO_Print_Init, GO_Print_Done
  use GO, only : GO_Timer_Init, GO_Timer_Done, &
                 GO_Timer_Def, GO_Timer_Start, GO_Timer_End
  use GO, only : TrcFile, Init, Done, ReadRc

  use EMEP_NMC_Driver, only : ParseArguments
  use EMEP_NMC_Driver, only : T_Driver_Eps
  use EMEP_NMC_Driver, only : T_Driver_Eps_Eval
  use EMEP_NMC_Driver, only : T_Driver_Eta
  use EMEP_NMC_Driver, only : T_Driver_Eta_f
  use EMEP_NMC_Driver, only : T_Driver_C_f
  use EMEP_NMC_Driver, only : T_Driver_D_f
  use EMEP_NMC_Driver, only : T_Driver_B_f
  use EMEP_NMC_Driver, only : T_Driver_XLX
  use EMEP_NMC_Driver, only : T_Driver_BB
  use EMEP_NMC_Driver, only : T_Driver_BB_Eval

  implicit none


  ! --- const ------------------------------------
  
  character(len=*), parameter   ::  rname = 'EMEP_NMC'


  ! --- external -------------------------------

#ifdef _OPENMP
  integer, external  ::  omp_get_thread_num
  integer, external  ::  omp_get_num_threads
#endif


  ! --- local ------------------------------------
  
  integer                   ::  status
  character(len=1024)       ::  rcfile
  type(TrcFile)             ::  rcF
  logical                   ::  apply
  type(T_Driver_Eps)        ::  Driver_Eps
  type(T_Driver_Eps_Eval)   ::  Driver_Eps_Eval
  type(T_Driver_Eta)        ::  Driver_Eta
  type(T_Driver_Eta_f)      ::  Driver_Eta_f
  type(T_Driver_C_f)        ::  Driver_C_f
  type(T_Driver_D_f)        ::  Driver_D_f
  type(T_Driver_B_f)        ::  Driver_B_f
  type(T_Driver_XLX)        ::  Driver_XLX
  type(T_Driver_BB)         ::  Driver_BB
  type(T_Driver_BB_Eval)    ::  Driver_BB_Eval
#ifdef _MPI
  integer                   ::  comm
  integer                   ::  npes, mype
#endif
  
  ! --- begin ------------------------------------
  
#ifdef _MPI
  ! start:
  call MPI_Init( status )
  IF_MPI_NOT_OK_RETURN(status=1)
  ! communicator:
  comm = MPI_COMM_WORLD
  ! count:
  call MPI_Comm_Size( comm, npes, status )
  IF_MPI_NOT_OK_RETURN(status=1)
  ! who am i ?
  call MPI_Comm_Rank( comm, mype, status )
  IF_MPI_NOT_OK_RETURN(status=1)
#endif
  
  ! init logging:
  call GO_Print_Init( status )
  IF_NOT_OK_RETURN(status=1)
  
  ! arguments:
  call ParseArguments( rcfile, status )
  ! usage displayed ?
  if ( status == -1 ) call Exit(0)
  ! check:
  IF_NOT_OK_RETURN(status=1)  
  
  ! info ...
  write (gol,'("")'); call goPr
  write (gol,'("***********************************************")'); call goPr
  write (gol,'("***")'); call goPr
  write (gol,'("*** EMEP DA NMC PostProcessing")'); call goPr
  write (gol,'("***")'); call goPr
  write (gol,'("***********************************************")'); call goPr
  write (gol,'("")'); call goPr
#ifdef _OPENMP
  write (gol,'("OpenMP: enabled        : yes")'); call goPr
  !$OMP PARALLEL
  write (gol,'("OpenMP: thread number  : ",i4," / ",i4)') omp_get_thread_num(), omp_get_num_threads(); call goPr
  !$OMP END PARALLEL
#else
  write (gol,'("OpenMP: enabled        : no")'); call goPr
#endif

  
  ! info ...
  write (gol,'("settings from : ",a)') trim(rcfile); call goPr
  
  ! settings:
  call Init( rcF, trim(rcfile), status )
  IF_NOT_OK_RETURN(status=1)
  
  ! collect eps samples?
  call ReadRc( rcF, 'nmc.eps.compute', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_Eps%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! compute:
    call Driver_Eps%Compute( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_Eps%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if
  
  ! evaluate sample covariance at selected points ?
  call ReadRc( rcF, 'nmc.sample-B', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_Eps_Eval%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! compute:
    call Driver_Eps_Eval%Evaluate( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_Eps_Eval%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if
  
  ! collect eps samples?
  call ReadRc( rcF, 'nmc.eta.compute', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_Eta%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! compute:
    call Driver_Eta%Compute( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_Eta%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if
  
  ! compute spectral transforms of eta samples?
  call ReadRc( rcF, 'nmc.eta_f.compute', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_Eta_f%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! compute:
    call Driver_Eta_f%Compute( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_Eta_f%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if

  ! compute level/comp covariances of spectral indices?
  call ReadRc( rcF, 'nmc.C_f.compute', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_C_f%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! compute:
    call Driver_C_f%Compute( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_C_f%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if
  
  ! angular average level/comp covariances?
  call ReadRc( rcF, 'nmc.D_f.compute', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_D_f%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! compute:
    call Driver_D_f%Compute( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_D_f%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if

  ! gamma and B_f :
  call ReadRc( rcF, 'nmc.B_f.compute', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_B_f%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! compute:
    call Driver_B_f%Compute( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_B_f%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if
  
  ! eigenvalue decompo of B_f:
  call ReadRc( rcF, 'nmc.XLX.compute', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_XLX%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! compute:
    call Driver_XLX%Compute( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_XLX%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if

  ! collection:
  call ReadRc( rcF, 'nmc.BB.collect', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_BB%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! collect:
    call Driver_BB%Collect( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_BB%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if
  
  ! approx parts of B matrix ?
  call ReadRc( rcF, 'nmc.approx-B', apply, status )
  IF_NOT_OK_RETURN(status=1)
  ! apply ?
  if ( apply ) then
    ! init:
    call Driver_BB_Eval%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    ! collect:
    call Driver_BB_Eval%Evaluate( status )
    IF_NOT_OK_RETURN(status=1)
    ! done:
    call Driver_BB_Eval%Done( status )
    IF_NOT_OK_RETURN(status=1)
  end if
  
  ! done with settings:
  call Done( rcF, status )
  IF_NOT_OK_RETURN(status=1)
  
  ! info ...
  write (gol,'("")'); call goPr
  write (gol,'("*** End")'); call goPr
  write (gol,'("")'); call goPr

  ! done with logging:
  call GO_Print_Done( status )
  IF_NOT_OK_RETURN(status=1)
  
#ifdef _MPI
  ! start:
  call MPI_Finalize( status )
  IF_MPI_NOT_OK_RETURN(status=1)
#endif

end program EMEP_NMC

