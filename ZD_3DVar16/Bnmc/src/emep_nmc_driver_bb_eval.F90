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

module EMEP_NMC_Driver_BB_Eval

  use GO                  , only : gol, goPr, goErr
  use EMEP_NMC_Common     , only : T_Points
  
  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_Driver_BB_Eval
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Driver_BB_Eval'
  
  ! --- types ----------------------------------------
  
  type T_Driver_BB_Eval
    ! points:
    type(T_Points)                    ::  points
    ! flags:
    logical                           ::  correlation
    logical                           ::  levtr
    ! statistics files created or read:
    character(len=1024)               ::  BB_filename
    character(len=1024)               ::  B_point_filename
  contains
    procedure   ::  Init            => Driver_BB_Eval_Init
    procedure   ::  Done            => Driver_BB_Eval_Done
    procedure   ::  Evaluate        => Driver_BB_Eval_Evaluate
  end type T_Driver_BB_Eval
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** model runs
  ! ***
  ! ********************************************************************


  subroutine Driver_BB_Eval_Init( self, rcF, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReadFromLine
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_BB_Eval), intent(out)          ::  self
    type(TrcFile), intent(in)                 ::  rcF
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_BB_Eval_Init'
    
    ! --- local ----------------------------------
    
    character(len=32)     ::  tvalue
    
    ! --- begin ----------------------------------
    
    ! point locations:
    call self%points%Init( rcF, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! input:
    call ReadRc( rcF, 'nmc.BB.filename', self%BB_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! flag:
    call ReadRc( rcF, 'nmc.approx-B.correlation', self%correlation, status )
    IF_NOT_OK_RETURN(status=1)
    ! flag:
    call ReadRc( rcF, 'nmc.approx-B.levtr', self%levtr, status )
    IF_NOT_OK_RETURN(status=1)

    ! target file:
    call ReadRc( rcF, 'nmc.approx-B.filename', self%B_point_filename, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_BB_Eval_Init


  ! ***
  

  subroutine Driver_BB_Eval_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Driver_BB_Eval), intent(inout)           ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_BB_Eval_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  irun
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%points%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine Driver_BB_Eval_Done


  ! ***
  
  
  !
  ! Evaluate parameterized covariance matrix at selected locations.
  !

  subroutine Driver_BB_Eval_Evaluate( self, status )
  
    use EMEP_NMC_Output, only : T_NMC_Output
    use EMEP_BCovarSqrt, only : T_BCovarSqrt
    use C3PO           , only : Dimension
    use C3PO           , only : HybrideLevelCoordinate, TimeCoordinate
    use C3PO           , only : LabelCoordinate
    use C3PO           , only : RealCoordinate
    
    ! --- in/out ---------------------------------
    
    class(T_Driver_BB_Eval), intent(inout)      ::  self
    integer, intent(out)                        ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/Driver_BB_Eval_Evaluate'
  
    ! --- local ----------------------------------
    
    type(T_BCovarSqrt)               ::  BCovarSqrt
    character(len=32)                ::  cv_label
    character(len=1)                 ::  cv_id
    character(len=32)                ::  cv_units

    type(T_NMC_Output)               ::  inp_file
    type(RealCoordinate)             ::  lon_coor, lat_coor
    type(HybrideLevelCoordinate)     ::  lev_coor
    type(LabelCoordinate)            ::  tracer_coor
    type(TimeCoordinate)             ::  ctime_coor
    integer                          ::  nlon, nlat, nlev
    integer                          ::  ntime
    integer                          ::  ntracer
    character(len=32)                ::  tracer_name
    character(len=32)                ::  tracer_units
    character(len=32)                ::  units_p
    character(len=32)                ::  vname_ps
    real, allocatable                ::  ps(:,:,:)   ! (nlon,nlat,ntime)
    character(len=32)                ::  vname, units
    real, allocatable                ::  S(:,:,:,:)  ! (nlon,nlat,nlev,ntracer)
    integer, allocatable             ::  S__varids(:)
    integer                          ::  itime
    integer                          ::  ipoint
    integer                          ::  ilon, ilat, ilev
    integer                          ::  itracer, itracer2
    
    type(T_NMC_Output)               ::  out_file
    character(len=1024)              ::  msg
    integer                          ::  ps__varid
    real, allocatable                ::  x(:,:,:,:)  ! (nlon,nlat,nlev,ntracer)
    real, allocatable                ::  Bx(:,:,:,:)  ! (nlon,nlat,nlev,ntracer)
    integer, allocatable             ::  Bx__varids(:,:)  ! (npoint,ntracer)

    integer                          ::  nw
    complex, allocatable             ::  w(:)
    
    real, allocatable                ::  Bzs(:,:,:,:)  ! (nlev,nlev,ntracer,ntracer)
    integer, allocatable             ::  Bzs__varids(:)  ! (npoint)
    integer                          ::  k
    
    ! --- begin ----------------------------------

    ! info ...
    write (gol,'("")'); call goPr
    write (gol,'("** Evaluate B e_i **")'); call goPr
    write (gol,'("")'); call goPr
    
    ! ~ covar sqrt
    
    ! info ...
    write (gol,'("read B sqrt ...")'); call goPr

    ! init covariance structure:
    call BCovarSqrt%Init( trim(self%BB_filename), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! set label:
    if ( self%correlation ) then
      cv_label = 'correlation'
      cv_id    = 'C'
      cv_units = '1'
    else
      cv_label = 'covariance'
      cv_id    = 'B'
      cv_units = 'tracer*tracer'
    end if
    

    ! ~ coordiantes etc:
    
    ! info ...
    write (gol,'("read coordinates ...")'); call goPr
          
    ! open file with covar sqrt:
    call inp_file%Open( trim(self%BB_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! create longitude coordinate:
    call lon_coor%Init( 'longitude', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lon_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lon_coor%Get_Dim( status, n=nlon )
    IF_NOT_OK_RETURN(status=1)

    ! create latitude coordinate:
    call lat_coor%Init( 'latitude', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lat_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lat_coor%Get_Dim( status, n=nlat )
    IF_NOT_OK_RETURN(status=1)

    ! create level coordinate:
    call lev_coor%Init( 'lev', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call lev_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call lev_coor%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)
    ! get name of pressure variable:
    call lev_coor%Get_Values( status, vname_ps=vname_ps, units_p=units_p )
    IF_NOT_OK_RETURN(status=1)
    
    ! create hour coordinate:
    call ctime_coor%Init( 'time', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file, here includes climatology bounds:
    call ctime_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call ctime_coor%Get_Dim( status, n=ntime )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( ps(nlon,nlat,ntime), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! loop over times:
    do itime = 1, ntime
      ! read surface pressure:
      call inp_file%Get_Field2D_Series( trim(vname_ps), itime, ps(:,:,itime), units, status )
      IF_NOT_OK_RETURN(status=1)
    end do

    ! create tracer coordinate:
    call tracer_coor%Init( 'tracer', status )
    IF_NOT_OK_RETURN(status=1)
    ! read from file:
    call tracer_coor%Read( inp_file, status )
    IF_NOT_OK_RETURN(status=1)
    ! size:
    call tracer_coor%Get_Dim( status, n=ntracer )
    IF_NOT_OK_RETURN(status=1)
      
    ! close:
    call inp_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ points
    
    ! loop over points:
    do ipoint = 1, self%points%n
      ! nearby index:
      call lon_coor%Get_Index( self%points%value(ipoint)%lon, self%points%value(ipoint)%ilon, status )
      IF_NOT_OK_RETURN(status=1)
      ! nearby index:
      call lat_coor%Get_Index( self%points%value(ipoint)%lat, self%points%value(ipoint)%ilat, status )
      IF_NOT_OK_RETURN(status=1)
      ! surface (top-down order!)
      self%points%value(ipoint)%ilev = nlev
    end do
    
    ! ~ output
    
    ! info ...
    write (gol,'("create B_point file ...")'); call goPr
    
    ! create new result file:
    call out_file%Create( trim(self%B_point_filename), status )
    IF_NOT_OK_RETURN(status=1)

    ! description:
    write (msg,'("evaluate ",a," from `",a,"`")') &
            trim(cv_label), trim(self%BB_filename)
    ! add:
    call out_file%Extend_History( trim(msg), status )
    IF_NOT_OK_RETURN(status=1)

    ! define coordinates in output file:
    call lon_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define level coordinate in output file:
    call lev_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define tracer coordinate in output file:
    call tracer_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define hour coordinate in output file:
    call ctime_coor%Def( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call out_file%Def_Field2D_Series( trim(vname_ps), trim(units_p), &
                                  lon_coor, lat_coor, ctime_coor, &
                                  ps__varid, status )
    IF_NOT_OK_RETURN(status=1)
    call out_file%Put_Att( ps__varid, 'standard_name', 'surface_air_pressure', status )
    IF_NOT_OK_RETURN(status=1)

    ! storage:
    allocate( S__varids(ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bx__varids(self%points%n,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over tracers:
    do itracer = 1, ntracer

      ! name and units:
      call BCovarSqrt%TracerGet( itracer, status, &
                                 name=tracer_name, units=tracer_units )
      IF_NOT_OK_RETURN(status=1)
      
      ! variable:
      write (vname,'("S_",a)') trim(tracer_name)
      ! define:
      call out_file%Def_Field3D_Series( trim(vname), trim(tracer_units), &
                              lon_coor, lat_coor, lev_coor, ctime_coor, &
                              S__varids(itracer), status )
      IF_NOT_OK_RETURN(status=1)
      call out_file%Put_Att( S__varids(itracer), 'long_name', &
                     trim(tracer_name)//' grid point standard deviation', status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over points:
      do ipoint = 1, self%points%n

        ! variable:
        write (vname,'(a,"_",a,"_",a)') trim(cv_id), trim(tracer_name), trim(self%points%value(ipoint)%name)
        ! define:
        call out_file%Def_Sample3D( trim(vname), trim(cv_units), &
                                lon_coor, lat_coor, lev_coor, ctime_coor, tracer_coor, &
                                Bx__varids(ipoint,itracer), status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'long_name', &
                       trim(cv_label)//' with '//trim(tracer_name)//' in '//trim(self%points%value(ipoint)%name), status )
        IF_NOT_OK_RETURN(status=1)
        ! write as attributes:
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'lon', self%points%value(ipoint)%lon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'lat', self%points%value(ipoint)%lat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'ilon', self%points%value(ipoint)%ilon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'ilat', self%points%value(ipoint)%ilat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bx__varids(ipoint,itracer), 'ilev', self%points%value(ipoint)%ilev, status )
        IF_NOT_OK_RETURN(status=1)

      end do ! points

    end do ! tracers
    
    ! evaluate level/tracer covariance?
    if ( self%levtr ) then

      ! storage for result:
      allocate( Bzs(nlev,nlev,ntracer,ntracer), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! output:
      allocate( Bzs__varids(self%points%n), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over points:
      do ipoint = 1, self%points%n

        ! variable:
        write (vname,'(a,"zs_",a)') trim(cv_id), trim(self%points%value(ipoint)%name)
        ! define:
        call out_file%Def_ACovar1( trim(vname), trim(cv_units), &
                                    lev_coor, tracer_coor, ctime_coor, &
                                    Bzs__varids(ipoint), status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'long_name', &
                       'vertical '//trim(cv_label)//' in '//trim(self%points%value(ipoint)%name), status )
        IF_NOT_OK_RETURN(status=1)
        ! write as attributes:
        call out_file%Put_Att( Bzs__varids(ipoint), 'lon', self%points%value(ipoint)%lon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'lat', self%points%value(ipoint)%lat, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'ilon', self%points%value(ipoint)%ilon, status )
        IF_NOT_OK_RETURN(status=1)
        call out_file%Put_Att( Bzs__varids(ipoint), 'ilat', self%points%value(ipoint)%ilat, status )
        IF_NOT_OK_RETURN(status=1)
        
      end do  ! points
      
    end if  ! level/tracer covar

    ! end:
    call out_file%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lon_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call lev_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call ctime_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! write coordinate variables:
    call tracer_coor%Write( out_file, status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over hours:
    do itime = 1, ntime
      ! write:
      call out_file%Put_Field2D_Series( ps__varid, itime, ps(:,:,itime), status )
      IF_NOT_OK_RETURN(status=1)
    end do ! times

    ! storage for std.dev.
    allocate( S(nlon,nlat,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for real state:
    allocate( x(nlon,nlat,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Bx(nlon,nlat,nlev,ntracer), stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over times:
    do itime = 1, ntime

      ! info ...
      write (gol,'("  time ",i0," / ",i0," ...")') itime, ntime; call goPr
      
      ! extract:
      call BCovarSqrt%TimeGet( itime, status, S=S )
      IF_NOT_OK_RETURN(status=1)
      ! loop over tracers:
      do itracer = 1, ntracer
        ! write sample:
        call out_file%Put_Field3D_Series( S__varids(itracer), itime, &
                                           S(:,:,:,itracer), status )
        IF_NOT_OK_RETURN(status=1)
      end do

      ! count:
      call BCovarSqrt%TimeGet( itime, status, nw=nw )
      IF_NOT_OK_RETURN(status=1)

      ! storage:
      if ( allocated(w) .and. (size(w) /= nw) ) then
        deallocate( w, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      if ( .not. allocated(w) ) then
        allocate( w(nw), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! loop over selected points:
      do ipoint = 1, self%points%n

        ! info ...
        write (gol,'("    point ",i0," / ",i0," ...")') ipoint, self%points%n; call goPr
      
        ! grid cell index:
        ilon = self%points%value(ipoint)%ilon
        ilat = self%points%value(ipoint)%ilat
        ilev = self%points%value(ipoint)%ilev
        
        ! loop over tracers (observed in point):
        do itracer = 1, ntracer

          ! info ...
          write (gol,'("      tracer ",i0," / ",i0," ...")') itracer, ntracer; call goPr

          ! unity vector:
          x = 0.0
          x(ilon,ilat,ilev,itracer) = 1.0

          ! evaluate:  B x =  B^{1/2} B^{H/2} x
          call BCovarSqrt%Forward( itime, x, w, status, correlation=self%correlation )
          IF_NOT_OK_RETURN(status=1)
          call BCovarSqrt%Reverse( itime, w, Bx, status, correlation=self%correlation )
          IF_NOT_OK_RETURN(status=1)
          
          ! info ...
          write (gol,'("        value range : ",f8.2," - ",f8.2)') minval(Bx), maxval(Bx); call goPr
          
          ! loop over correlated tracers:
          do itracer2 = 1, ntracer
            ! write sample:
            call out_file%Put_Sample3D( Bx__varids(ipoint,itracer), itime, itracer2, &
                                               Bx(:,:,:,itracer2), status )
            IF_NOT_OK_RETURN(status=1)
          end do  ! correlated tracers
          
        end do ! tracers
        
        ! ~
        
        ! evaluate full level/tracer covar?
        if ( self%levtr ) then

          ! info ...
          write (gol,'("      evaluate level/tracer covariance ...")'); call goPr

          ! loop over tracers (observed in point):
          do itracer = 1, ntracer
            ! loop:
            do k = 1, nlev
              ! info ...
              write (gol,'("        level ",i0," / ",i0)') k, nlev; call goPr
              ! unity vector:
              x = 0.0
              x(ilon,ilat,k,itracer) = 1.0
              ! evaluate:  B x =  B^{1/2} B^{H/2} x
              call BCovarSqrt%Forward( itime, x, w, status, correlation=self%correlation )
              IF_NOT_OK_RETURN(status=1)
              call BCovarSqrt%Reverse( itime, w, Bx, status, correlation=self%correlation )
              IF_NOT_OK_RETURN(status=1)
              ! info ...
              write (gol,'("          value range : ",f8.2," - ",f8.2)') minval(Bx), maxval(Bx); call goPr
              ! copy:
              Bzs(k,:,itracer,:) = Bx(ilon,ilat,:,:)
              Bzs(:,k,:,itracer) = Bx(ilon,ilat,:,:)
            end do ! levels
          end do ! tracers

          ! write sample:
          call out_file%Put_ACovar1( Bzs__varids(ipoint), itime, Bzs, status  )
          IF_NOT_OK_RETURN(status=1)

        end if ! level/tracer covariance
        
        ! ~
          
      end do  ! points

    end do ! times
    
    ! clear:
    deallocate( w, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( Bx, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( S, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! evaluate level/tracer covariance?
    if ( self%levtr ) then
      ! clear:
      deallocate( Bzs, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! clear:
      deallocate( Bzs__varids, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if ! level/tracer covar

    ! ~ done with input
      
    ! clear:
    call lon_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lat_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call lev_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
    call ctime_coor%Done( status )
    IF_NOT_OK_RETURN(status=1)
      
    ! clear:
    deallocate( ps, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ done with output

    ! close output:
    call out_file%Close( status )
    IF_NOT_OK_RETURN(status=1)

    deallocate( Bx__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( S__varids, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ done with covar
    
    ! clear:
    call BCovarSqrt%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ~ 

    ! ok
    status = 0
  
  end subroutine Driver_BB_Eval_Evaluate
  


end module EMEP_NMC_Driver_BB_Eval
