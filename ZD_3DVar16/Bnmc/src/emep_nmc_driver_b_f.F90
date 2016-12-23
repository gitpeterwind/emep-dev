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

module EMEP_NMC_Driver_B_f

  use GO                  , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_B_f
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_B_f'
  
  ! --- types ----------------------------------------
  
  type T_Driver_B_f
    ! statistics files created or read:
    character(len=1024)               ::  D_f_filename
    character(len=1024)               ::  gamma_filename
    character(len=1024)               ::  B_f_filename
  contains
    procedure   ::  Init            => Driver_B_f_Init
    procedure   ::  Done            => Driver_B_f_Done
    procedure   ::  Compute         => Driver_B_f_Compute
  end type T_Driver_B_f
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_B_f_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_B_f), intent(out)          ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_B_f_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! input/output:
    call ReadRc( rcF, 'nmc.D_f.filename', self%D_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.gamma.filename', self%gamma_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.B_f.filename', self%B_f_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_B_f_Init


  ! ***
  

  subroutine Driver_B_f_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_B_f), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_B_f_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
  
  end subroutine Driver_B_f_Done


  ! ***
  
  
  !
  ! Compute B_f = D_f/sqrt(gamma*gamma)
  ! 

  subroutine Driver_B_f_Compute( self, status )
  
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_Driver_B_f), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_B_f_Compute'
  
    ! gonio:
    real, parameter     ::  pi = 4.0 * atan(1.0)
    
    ! --- local ----------------------------------

    type(T_NMC_Output)               ::  D_f_file
    type(RealCoordinate)             ::  akstar_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(LabelCoordinate)            ::  tracer_coor
    type(TimeCoordinate)             ::  ctime_coor
    integer                          ::  nakstar
    integer                          ::  nlev
    integer                          ::  ntracer
    integer                          ::  ntime
    real                             ::  p0
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    integer, allocatable             ::  fixed(:,:)  ! (nlev,nvar)
    character(len=32)                ::  fixed__units
    integer                          ::  fixed__varid
    real, allocatable                ::  akstar(:)  ! (nakstar)
    real, allocatable                ::  D_f(:,:,:,:,:)  ! (nakstar,nlev,nlev,ntr,ntr)
    character(len=32)                ::  D_f__units
    integer                          ::  itime

    character(len=1024)              ::  msg
    real, allocatable                ::  ps(:,:)   ! (nakstar,ntime)
    integer                          ::  ps__varid

    type(T_NMC_Output)               ::  gamma_file
    real, allocatable                ::  gamma(:,:,:)  ! (nakstar,nlev,ntr)
    integer                          ::  gamma__varid
    
    type(T_NMC_Output)               ::  B_f_file
    real, allocatable                ::  B_f(:,:,:,:,:)  ! (nakstar,nlev,nlev,ntr,ntr)
    integer                          ::  B_f__varid

    integer                          ::  ilev, ilev1, ilev2
    integer                          ::  itracer, itracer1, itracer2
    integer                          ::  ik
    real                             ::  D_f_diagelem_sum
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Compute gamma and B_f **")'); call goPr
    write (gol,'("")'); call goPr
    
    ! ~ input samples
          
    ! open file with D_f:
    call D_f_file%Open( trim(self%D_f_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create average k* coordinate:
    call akstar_coor%Init( 'akstar', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call akstar_coor%Read( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call akstar_coor%Get_Dim( status, n=nakstar )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( akstar(nakstar), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! extract:
    call akstar_coor%Get_Values( status, values=akstar )
    IF_NOT_OK_RETURN(status=1)

    ! create level coordinate:
    call lev_coor%Init( 'lev', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lev_coor%Read( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lev_coor%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=ntracer )
    IF_NOT_OK_RETURN(status=1)

    ! create (climatological) time coordinate:
    call ctime_coor%Init( 'time', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call ctime_coor%Read( D_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call ctime_coor%Get_Dim( status, n=ntime )
    IF_NOT_OK_RETURN(status=1)
    
    ! get name of pressure variable:
    call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p, p0=p0 )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage:
    allocate( D_f(nakstar,nlev,nlev,ntracer,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( fixed(nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! read:
    call D_f_file%Get_IField2D( 'fixed', fixed, fixed__units, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ output general
    
    ! storage:
    allocate( ps(nakstar,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill with reference pressure:
    ps = p0

    ! ~ output: gamma

    ! create new result file:
    call gamma_file%Create( trim(self%gamma_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("normalization of angular averages of spectral level/tracer covariances from `",a,"`")') trim(self%D_f_filename)
    ! add:
    call gamma_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define ak* coordinate in output file:
    call akstar_coor%Def( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define level coordinate in output file:
    call lev_coor%Def( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define tracer coordinate in output file:
    call tracer_coor%Def( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call gamma_file%Def_Field2D( trim(vname_ps), trim(units_p), &
                                 akstar_coor, ctime_coor, &
                                 ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call gamma_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call gamma_file%Def_IField2D( 'fixed', trim(fixed__units), &
                                        lev_coor, tracer_coor, &
                                        fixed__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call gamma_file%Put_Att( fixed__varid, 'long_name', 'fixed layer (1=fixed, 0=dynamic)', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call gamma_file%Def_AVar( 'gamma', '1', &
                              akstar_coor, lev_coor, tracer_coor, ctime_coor, &
                              gamma__varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call gamma_file%Put_Att( gamma__varid, 'long_name', 'scaling for angular average of covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)
    
    ! end:
    call gamma_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call akstar_coor%Write( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Write( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Write( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Write( gamma_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write:
    call gamma_file%Put_Field2D( ps__varid, ps, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write flags:
    call gamma_file%Put_IField2D( fixed__varid, fixed, status )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( gamma(nakstar,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ output: B_f

    ! create new result file:
    call B_f_file%Create( trim(self%B_f_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("angular averages of spectral level/tracer covariances from `",a,"`")') trim(self%D_f_filename)
    write (msg,'(a," scalled with `",a,"`")') trim(msg), trim(self%gamma_filename)
    ! add:
    call B_f_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define ak* coordinate in output file:
    call akstar_coor%Def( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! define level coordinate in output file:
    call lev_coor%Def( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define tracer coordinate in output file:
    call tracer_coor%Def( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call B_f_file%Def_Field2D( trim(vname_ps), trim(units_p), &
                                 akstar_coor, ctime_coor, &
                                 ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call B_f_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call B_f_file%Def_IField2D( 'fixed', trim(fixed__units), &
                                   lev_coor, tracer_coor, &
                                   fixed__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call B_f_file%Put_Att( fixed__varid, 'long_name', 'fixed layer (1=fixed, 0=dynamic)', status )
    IF_NOT_OK_RETURN(status=1)

    ! define, obtain fill value for no-data:
    call B_f_file%Def_ACovar( 'B_f', '1', &
                              akstar_coor, lev_coor, tracer_coor, ctime_coor, &
                              B_f__varid, status )
    IF_NOT_OK_RETURN(status=1)
    ! annote:
    call B_f_file%Put_Att( B_f__varid, 'long_name', 'scaled angular average of covariance of spectral values', status )
    IF_NOT_OK_RETURN(status=1)
    
    ! end:
    call B_f_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call akstar_coor%Write( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Write( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Write( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Write( B_f_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write:
    call B_f_file%Put_Field2D( ps__varid, ps, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! write flags:
    call B_f_file%Put_IField2D( fixed__varid, fixed, status )
    IF_NOT_OK_RETURN(status=1)
      
    ! storage:
    allocate( B_f(nakstar,nlev,nlev,ntracer,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ time loop

    ! info ...
    write (gol,'("dimensions:")'); call goPr
    write (gol,'("  number of k* values   : ",i0)') nakstar; call goPr
    write (gol,'("  number of time values : ",i0)') ntime; call goPr
    write (gol,'("  number of tracers     : ",i0)') ntracer; call goPr
    
    ! loop over hours:
    do itime = 1, ntime
    
      ! info ...
      write (gol,'("  itime ",i0)') itime; call goPr

      ! read angular averaged spectral covariances:
      call D_f_file%Get_ACovar( 'D_f', itime, D_f, D_f__units, status )
      IF_NOT_OK_RETURN(status=1)
      
      !
      !                        D_f(k,l;k,l;k*)
      ! gamma(k,l;k*) = ---------------------------------
      !                 sum_k*' 2 pi k*' D_f(k,l;k,l;k*')
      !
      ! loop over diagonal elements:
      do itracer = 1, ntracer
        do ilev = 1, nlev
        
          ! sum:
          D_f_diagelem_sum = sum( 2*pi*akstar * D_f(:,ilev,ilev,itracer,itracer) )
        
          ! non-zero variance?
          if ( D_f_diagelem_sum > 0.0 ) then
            ! fill scale factor:
            gamma(:,ilev,itracer) = D_f(:,ilev,ilev,itracer,itracer) / D_f_diagelem_sum
          else
            ! no scaling:
            gamma(:,ilev,itracer) = 0.0
          end if
          
        end do ! levels
      end do  ! tracers
      
      !
      !                                D_f(k,l;k',l';k*)
      !  B_f(k,l;k',l';k*) = -------------------------------------
      !                      sqrt( gamma(k,l;k*) gamma(k',l';k*) )
      !
      ! loop over covariance elements:
      do itracer1 = 1, ntracer
        do itracer2 = 1, ntracer
          do ilev1 = 1, nlev
            do ilev2 = 1, nlev
            
              ! loop over wave numbers:
              do ik = 1, nakstar
              
                ! trap zero variances:
                if ( (gamma(ik,ilev1,itracer1) > 0.0) .and. &
                     (gamma(ik,ilev2,itracer2) > 0.0)       ) then
                  ! fill scalled version:
                  B_f(ik,ilev1,ilev2,itracer1,itracer2) = &
                    D_f(ik,ilev1,ilev2,itracer1,itracer2) / &
                    sqrt( gamma(ik,ilev1,itracer1) * gamma(ik,ilev2,itracer2) )
                else
                  ! no variance (top level?), so remain zero:
                  B_f(ik,ilev1,ilev2,itracer1,itracer2) = 0.0
                end if  ! zero variance
                
              end do  ! k*
              
            end do  ! ilev2
          end do ! ilev1
        end do ! itracer2
      end do ! itracer1

      ! write scaling:
      call gamma_file%Put_AVar( gamma__varid, itime, gamma, status )
      IF_NOT_OK_RETURN(status=1)
      ! write covar:
      call B_f_file%Put_ACovar( B_f__varid, itime, B_f, status )
      IF_NOT_OK_RETURN(status=1)

    end do  ! times

    ! ~ done with output
    
    ! clear:
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( gamma, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! close output:
    call gamma_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( B_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! close output:
    call B_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ done with input
    
    ! clear:
    deallocate( akstar, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( fixed, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( D_f, stat=status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    call akstar_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call tracer_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! close input:
    call D_f_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~ 

    ! ok
    status = 0
  
  end subroutine Driver_B_f_Compute
  
  
end module EMEP_NMC_Driver_B_f
