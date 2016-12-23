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

module EMEP_NMC_Driver_Eta

  use GO                  , only : gol, goPr, goErr
  use EMEP_NMC_Common     , only : T_ModelVars
  use GO                  , only : TDate

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_Eta
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_Eta'
  
  ! --- types ----------------------------------------
  
  type T_Driver_Eta
    ! time range:
    type(TDate)                       ::  tday1, tday2
    ! variables:
    type(T_ModelVars)                 ::  mvars
    ! statistics files created or read:
    character(len=1024)               ::  eps_samples_filename
    character(len=1024)               ::  eps_stats_filename
    character(len=1024)               ::  eta_samples_filename
  contains
    procedure   ::  Init            => Driver_Eta_Init
    procedure   ::  Done            => Driver_Eta_Done
    procedure   ::  Compute         => Driver_Eta_Compute
  end type T_Driver_Eta
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_Eta_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine, goSplitString
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eta), intent(out)          ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eta_Init'
    
    ! --- local ----------------------------------
    
    character(len=32)     ::  tvalue
    
    ! --- begin ----------------------------------
    
    ! start time:
    call ReadRc( rcF, 'nmc.timerange.t1', tvalue, status )
    IF_NOT_OK_RETURN(status=1)
    ! extract:
    call goReadFromLine( tvalue, self%tday1, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! end time:
    call ReadRc( rcF, 'nmc.timerange.t2', tvalue, status )
    IF_NOT_OK_RETURN(status=1)
    ! extract:
    call goReadFromLine( tvalue, self%tday2, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! variable names:
    call self%mvars%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! input/output:
    call ReadRc( rcF, 'nmc.eps.samples.filename', self%eps_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.eps.stats.filename', self%eps_stats_filename, status )
    IF_NOT_OK_RETURN(status=1)
    call ReadRc( rcF, 'nmc.eta.samples.filename', self%eta_samples_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_Eta_Init


  ! ***
  

  subroutine Driver_Eta_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eta), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eta_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%mvars%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_Eta_Done


  ! ***
  
  
  !
  ! Compute eta samples:
  ! 
  !    eta_comp(sample,hour,lev,lat,lon)  =  eps_comp(sample,hour,lev,lat,lon) / eps_comp_stdv(hour,lev,lat,lon)
  !

  subroutine Driver_Eta_Compute( self, status )
  
    use GO             , only : TDate, IncrDate, operator(+), operator(>), wrtgol
    use GO             , only : goReplace
    use EMEP_NMC_Output, only : T_NMC_Output
    use C3PO           , only : RealCoordinate, HybrideLevelCoordinate, TimeCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_Eta), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_Eta_Compute'
    
    ! --- local ----------------------------------

    type(TDate)                      ::  tday, t
    character(len=1024)              ::  filename
    character(len=1024)              ::  eps_filename
    character(len=1024)              ::  eta_filename
    integer                          ::  eps_isample
    integer                          ::  eta_isample
    logical                          ::  setup
    type(T_NMC_Output)               ::  eps_stats_file
    type(T_NMC_Output)               ::  eps_file
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(TimeCoordinate)             ::  ctime_coor
    type(TimeCoordinate)             ::  eps_sample_coor
    type(TimeCoordinate)             ::  eta_sample_coor
    integer                          ::  nlon, nlat, nlev
    integer                          ::  ntime
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    integer, allocatable             ::  eps__varids(:)   ! (nvar)
    real, allocatable                ::  ps(:,:)
    real, allocatable                ::  eps(:,:,:)
    real, allocatable                ::  eta(:,:,:)
    character(len=32)                ::  vname, units
    real, allocatable                ::  eps_mean(:,:,:,:,:)
    real, allocatable                ::  eps_stdv(:,:,:,:,:)
    integer, allocatable             ::  fixed(:,:)  ! (nlev,nvar)
    character(len=32)                ::  fixed__units
    integer, allocatable             ::  fixed__varids(:)
    integer                          ::  ivar
    integer                          ::  itime
    
    type(T_NMC_Output)               ::  eta_file
    character(len=1024)              ::  msg
    integer                          ::  ps__varid
    integer, allocatable             ::  eta__varids(:)   ! (nvar)
    
    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Compute eta samples **")'); call goPr
    write (gol,'("")'); call goPr

    ! start time:
    tday = self%tday1
    ! set flag:
    setup = .false.
    ! init current filenames:
    eps_filename = ''
    eta_filename = ''
    ! time loop:
    do
      ! info ...
      call wrtgol( 'day ', tday ); call goPr
      
      ! ~ input samples

      ! input file:
      filename = trim(self%eps_samples_filename)
      ! replace keys:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{hh}'  , '(i2.2)', tday%hour , status )
      IF_NOT_OK_RETURN(status=1)
      
      ! new ?
      if ( trim(filename) /= trim(eps_filename) ) then
      
        ! open?
        if ( len_trim(eps_filename) > 0 ) then
          ! close:
          call eps_file%Close( status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! copy:
        eps_filename = trim(filename)

        ! open file with eps samples:
        call eps_file%Open( trim(filename), status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init counter:
        eps_isample = 1
        
      else
      
        ! next sample for already open file:
        eps_isample = eps_isample + 1
        
      end if
      
      ! ~ setup

      ! initial setup?
      if ( .not. setup ) then

        ! create longitude coordinate:
        call lon_coor%Init( 'longitude', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lon_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lon_coor%Get_Dim( status, n=nlon )
        IF_NOT_OK_RETURN(status=1)

        ! create latitude coordinate:
        call lat_coor%Init( 'latitude', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lat_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lat_coor%Get_Dim( status, n=nlat )
        IF_NOT_OK_RETURN(status=1)

        ! create level coordinate:
        call lev_coor%Init( 'lev', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call lev_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call lev_coor%Get_Dim( status, n=nlev )
        IF_NOT_OK_RETURN(status=1)

        ! create hour coordinate:
        call ctime_coor%Init( 'time', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file, here includes climatology bounds:
        call ctime_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)
        ! size:
        call ctime_coor%Get_Dim( status, n=ntime )
        IF_NOT_OK_RETURN(status=1)

        ! create sample coordinate:
        call eps_sample_coor%Init( 'sample', status )
        IF_NOT_OK_RETURN(status=1)
        ! read from file:
        call eps_sample_coor%Read( eps_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! get name of pressure variable:
        call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( eps__varids(self%mvars%n) )
        ! loop over variables:
        do ivar = 1, self%mvars%n
          ! get id:
          call eps_file%Get_VarID( trim(self%mvars%value(ivar)%name), eps__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do

        ! storage:
        allocate( ps(nlon,nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( eps(nlon,nlat,nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ eps statistics:

        ! open file with eps samples:
        call eps_stats_file%Open( trim(self%eps_stats_filename), status )
        IF_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( eps_mean(nlon,nlat,nlev,ntime,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( eps_stdv(nlon,nlat,nlev,ntime,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( fixed(nlev,self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%mvars%n
          !
          ! variable:
          write (vname,'(a,"_mean")') trim(self%mvars%value(ivar)%name)
          ! loop over hours:
          do itime = 1, ntime
            ! read:
            call eps_stats_file%Get_Field3D_Series( trim(vname), itime, eps_mean(:,:,:,itime,ivar), units, status )
            IF_NOT_OK_RETURN(status=1)
          end do ! time
          !
          ! variable:
          write (vname,'(a,"_stdv")') trim(self%mvars%value(ivar)%name)
          ! loop over times:
          do itime = 1, ntime
            ! read:
            call eps_stats_file%Get_Field3D_Series( trim(vname), itime, eps_stdv(:,:,:,itime,ivar), units, status )
            IF_NOT_OK_RETURN(status=1)
          end do ! time
          !
          ! variable:
          write (vname,'(a,"_fixed")') trim(self%mvars%value(ivar)%name)
          call eps_stats_file%Get_IField1D( trim(vname), fixed(:,ivar), fixed__units, status )
          IF_NOT_OK_RETURN(status=1)
          !
        end do ! var

        ! close:
        call eps_stats_file%Close( status )
        IF_NOT_OK_RETURN(status=1)

        ! ~ output storage

        ! storage:
        allocate( eta(nlon,nlat,nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! access:
        allocate( eta__varids(self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( fixed__varids(self%mvars%n), stat=status )
        IF_NOT_OK_RETURN(status=1)

        ! ~

        ! reset flag:
        setup = .true.

      end if  ! not setup yet

      ! ~ output
        
      ! input file:
      filename = trim(self%eta_samples_filename)
      ! replace keys:
      call goReplace( filename, '%{yyyy}', '(i4.4)', tday%year , status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{mm}'  , '(i2.2)', tday%month, status )
      IF_NOT_OK_RETURN(status=1)
      call goReplace( filename, '%{dd}'  , '(i2.2)', tday%day  , status )
      IF_NOT_OK_RETURN(status=1)
      
      ! new ?
      if ( trim(filename) /= trim(eta_filename) ) then
      
        ! open?
        if ( len_trim(eta_filename) > 0 ) then
          ! write times of samples:
          call eta_sample_coor%Write( eta_file, status )
          IF_NOT_OK_RETURN(status=1)
          ! close:
          call eta_file%Close( status )
          IF_NOT_OK_RETURN(status=1)
        end if
        
        ! copy:
        eta_filename = trim(filename)

        ! create new result file:
        call eta_file%Create( trim(filename), status )
        IF_NOT_OK_RETURN(status=1)

        ! description:
        write (msg,'("eta samples from `",a,"` and `",a,"`")') &
                 trim(self%eps_samples_filename), trim(self%eps_stats_filename)
        ! add:
        call eta_file%Extend_History( trim(msg), status )
        IF_NOT_OK_RETURN(status=1)

        ! define coordinates in output file:
        call lon_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lev_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call ctime_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! output sample coor:
        call eta_sample_coor%Init( 'sample', status )
        IF_NOT_OK_RETURN(status=1)
        ! set unlimited size:
        call eta_sample_coor%Set_Dim( status, unlimited=.true. )
        IF_NOT_OK_RETURN(status=1)
        ! fill attributes:
        call eta_sample_coor%Set_Values( status, units_step='days', units_t0=self%tday1 )
        IF_NOT_OK_RETURN(status=1)
        ! define sample coordinate in output file:
        call eta_sample_coor%Def( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! init counter:
        eta_isample = 1

        ! define:
        call eta_file%Def_Sample2D( trim(vname_ps), trim(units_p), &
                                      lon_coor, lat_coor, ctime_coor, eta_sample_coor, &
                                      ps__varid, status )
        IF_NOT_OK_RETURN(status=1)
        call eta_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%mvars%n
          ! define, result will be unit-less:
          call eta_file%Def_Sample3D( trim(self%mvars%value(ivar)%name), '1', &
                                        lon_coor, lat_coor, lev_coor, ctime_coor, eta_sample_coor, &
                                        eta__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
          call eta_file%Put_Att( eta__varids(ivar), 'long_name', trim(self%mvars%value(ivar)%name)//' (eps-eps_mean)/eps_stdv', status )
          IF_NOT_OK_RETURN(status=1)

          ! define:
          call eta_file%Def_IField1D( trim(self%mvars%value(ivar)%name)//'_fixed', trim(fixed__units), &
                                        lev_coor, fixed__varids(ivar), status )
          IF_NOT_OK_RETURN(status=1)
          call eta_file%Put_Att( fixed__varids(ivar), 'long_name', trim(self%mvars%value(ivar)%name)//' fixed layer (1=fixed, 0=dynamic)', status )
          IF_NOT_OK_RETURN(status=1)

        end do  ! var

        ! end:
        call eta_file%EndDef( status )
        IF_NOT_OK_RETURN(status=1)

        ! write coordinate variables:
        call lon_coor%Write( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lat_coor%Write( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call lev_coor%Write( eta_file, status )
        IF_NOT_OK_RETURN(status=1)
        call ctime_coor%Write( eta_file, status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%mvars%n
          ! write:
          call eta_file%Put_IField1D( fixed__varids(ivar), fixed(:,ivar), status )
          IF_NOT_OK_RETURN(status=1)
        end do
        
      else
      
        ! next record in already opened file:
        eta_isample = eta_isample + 1
        
      end if  ! new output file

      ! ~
      
      ! get sample info:
      call eps_sample_coor%Get_Value( eps_isample, status, t=t )
      IF_NOT_OK_RETURN(status=1)
      ! copy to output coordinate:
      call eta_sample_coor%Set_Value( eta_isample, status, t=t )
      IF_NOT_OK_RETURN(status=1)

      ! loop over hours:
      do itime = 1, ntime

        ! info ...
        write (gol,'("  time ",i4," ...")') itime; call goPr
    
        ! read sample:
        call eps_file%Get_Sample2D( vname_ps, itime, eps_isample, ps, units, status )
        IF_NOT_OK_RETURN(status=1)
        
        ! write sample:
        call eta_file%Put_Sample2D( ps__varid, itime, eta_isample, ps, status )
        IF_NOT_OK_RETURN(status=1)

        ! loop over variables:
        do ivar = 1, self%mvars%n
    
          ! info ...
          write (gol,'("    variable ",i0," ",a)') itime, trim(self%mvars%value(ivar)%name); call goPr
        
          ! read sample:
          call eps_file%Get_Sample3D( trim(self%mvars%value(ivar)%name), itime, eps_isample, eps, units, status )
          IF_NOT_OK_RETURN(status=1)
          
          ! convert:
          where ( eps_stdv(:,:,:,itime,ivar) > 0.0 )
            eta = ( eps - eps_mean(:,:,:,itime,ivar) ) / eps_stdv(:,:,:,itime,ivar)
          elsewhere
            eta = 0.0
          end where
          
          ! write sample:
          call eta_file%Put_Sample3D( eta__varids(ivar), itime, eta_isample, eta, status )
          IF_NOT_OK_RETURN(status=1)
          
        end do ! var
        
      end do  ! hour

      ! next:
      tday = tday + IncrDate(day=1)
      ! end ?
      if ( tday > self%tday2 ) exit
      
    end do  ! samples
    
    ! info ...
    write (gol,'("ok")'); call goPr

    ! ~ close output
    
    ! write times of samples:
    call eta_sample_coor%Write( eta_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! close ouptut:
    call eta_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! done:
    call eta_sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( eta__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( fixed__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( eta, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ clear stats

    ! clear stats:
    deallocate( eps_mean, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( eps_stdv, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( fixed )
    IF_NOT_OK_RETURN(status=1)

    ! ~ close input
    
    ! clear:
    deallocate( ps )
    deallocate( eps )

    ! clear:
    deallocate( eps__varids )
      
    ! clear:
    call lon_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call eps_sample_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! close input:
    call eps_file%Close( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ~

    ! ok
    status = 0
  
  end subroutine Driver_Eta_Compute


end module EMEP_NMC_Driver_Eta
  
  
