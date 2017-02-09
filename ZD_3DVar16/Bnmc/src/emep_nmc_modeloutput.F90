!######################################################################
!
! EMEP DA NMC - model output routines
!
! Usage:
!
!    ! external:
!    use GO, only : TrcFile
!    use GO, only : TDate
!
!    ! storage type:
!    type(T_ModelOutputSeries)   ::   mos
!    ! settings:
!    type(TrcFile)               ::   rcF
!    ! time value:
!    type(TDate)                 ::   t
!    ! dims:
!    integer                     ::  nlon, nlat, nlev
!    ! data:
!    real, allocatable           ::  data(:,:,:)
!
!    ! init settings:
!    ... rcF ...
!
!    ! initialize from rcfile settings, 
!    ! specify base key and run id:
!    call mos%Init( rcF, 'myruns', '00', status )
!    if (status/=0) stop
!
!    ! target time:
!    ... t ...
!
!    ! fill buffer for requested time,
!    ! probably all records from a single file:
!    call mos%Init( t, status )
!    if (status/=0) stop
!
!    ! get dimensions:
!    call mos%Get( status, nlon=nlon, nlat=nlat, nlev=nlev )
!    if (status/=0) stop
!
!    ! storage:
!    allocate( data(nlon,nlat,nlev) )
!
!    ! read record from buffer:
!    call mos%Get( t, data, status )
!    if (status/=0) stop
!
!    ! done:
!    call mos%Done( status )
!    if (status/=0) stop
!
! Rcfile settings:
!
!    ! template for filename, keys are replaced:
!    !    %{id}       :  id
!    !    %{yyyy}     :  4-digit year
!    !    %{mm}       :  2-digit month
!    myruns.00.file       :  /data/output/NMC_%{id}_%{yyyy}%{mm}.nc
!
!    ! variable name:
!    myruns.00.variable   :  NO2
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
!######################################################################

module EMEP_NMC_ModelOutput

  use GO    , only : gol, goPr, goErr
  use GO    , only : TDate
  use NetCDF, only : NF90_NOERR, NF90_StrError

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_ModelOutputSeries
  

  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_ModelOutput'
  

  ! --- types ----------------------------------------
  
  type T_ModelOutput
    character(len=32)           ::  name
    real, allocatable           ::  data(:,:,:,:)  ! (nlon,nlat,nlev,ntime)
    character(len=32)           ::  units
  end type T_ModelOutput
  
  ! *
  
  type T_ModelOutputSeries
    ! filename template:
    character(len=1024)         ::  filename_template
    ! current filename:
    character(len=1024)         ::  filename_curr
    ! annote:
    character(len=1024)         ::  label
    ! flags:
    logical                     ::  gridsetup
    logical                     ::  filled
    ! grid:
    integer                     ::  nlon, nlat
    real, allocatable           ::  lons(:)
    real, allocatable           ::  lats(:)
    ! levels:
    integer                     ::  nlev
    !real, allocatable           ::  sigma(:)   ! (nlev)
    !real                        ::  ptop
    !real, allocatable           ::  ps(:,:,:)  ! (nlon,nlat,ntime)
    !character(len=32)           ::  vname_ps
    !character(len=32)           ::  vname_ptop
    !character(len=32)           ::  units_p
    character(len=32)           ::  vname_ap, vname_b
    character(len=32)           ::  vname_ps, vname_p0
    character(len=32)           ::  units_p
    real                        ::  p0
    real, allocatable           ::  sigma(:)   ! (nlev)
    real, allocatable           ::  hyam(:), hybm(:) ! (nlev)
    real, allocatable           ::  isigma(:)   ! (nlev+1)
    real, allocatable           ::  hyai(:), hybi(:) ! (nlev+1)
    real, allocatable           ::  ps(:,:,:)  ! (nlon,nlat,ntime)
    ! time values:
    integer                     ::  ntime
    type(TDate), allocatable    ::  tt(:)
    ! data:
    integer                             ::  nvar
    type(T_ModelOutput), allocatable    ::  var(:)
    !
  contains
    procedure   ::  Init        => ModelOutputSeries_Init
    procedure   ::  Done        => ModelOutputSeries_Done
    procedure   ::  Get         => ModelOutputSeries_Get
    procedure   ::  ReadBuffer  => ModelOutputSeries_ReadBuffer
    procedure   ::  GetRecord   => ModelOutputSeries_GetRecord
  end type T_ModelOutputSeries
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** output series
  ! ***
  ! ********************************************************************


  subroutine ModelOutputSeries_Init( self, rcF, rcbase, id, mvars, status )
  
    use GO, only : TrcFile, ReadRc
    use GO, only : goReplace, goSplitString
    use EMEP_NMC_Common, only : T_ModelVars
  
    ! --- in/out ---------------------------------
    
    class(T_ModelOutputSeries), intent(out)   ::  self
    type(TrcFile), intent(in)                 ::  rcF
    character(len=*), intent(in)              ::  rcbase
    character(len=*), intent(in)              ::  id
    type(T_ModelVars), intent(in)             ::  mvars
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelOutputSeries_Init'
    
    ! --- local ----------------------------------
    
    character(len=32)   ::  rckey
    integer             ::  ivar
    
    ! --- begin ----------------------------------
  
    ! combine keys:
    write (rckey,'(a,".",a)') trim(rcbase), trim(id)
      
    ! read filename template:
    call ReadRc( rcF, trim(rckey)//'.file', self%filename_template, status )
    IF_NOT_OK_RETURN(status=1)
    ! replace id values now:
    call goReplace( self%filename_template, '%{id}', trim(id), status )
    IF_NOT_OK_RETURN(status=1)
    
    ! info ...
    write (gol,'("filename template : ",a)') trim(self%filename_template); call goPr
    
    ! nothing read yet:
    self%filename_curr = ''
    self%gridsetup = .false.
    self%filled = .false.
    
    ! number of variables:
    self%nvar = mvars%n
    ! storage:
    allocate( self%var(self%nvar) )
    ! copy names:
    do ivar = 1, mvars%n
      ! name of input variable:
      call ReadRc( rcF, trim(rckey)//'.var.'//trim(mvars%value(ivar)%name), self%var(ivar)%name, status )
      IF_NOT_OK_RETURN(status=1)
    end do
    
    ! ok
    status = 0
  
  end subroutine ModelOutputSeries_Init


  ! ***
  

  subroutine ModelOutputSeries_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_ModelOutputSeries), intent(inout)   ::  self
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelOutputSeries_Done'
    
    ! --- local ----------------------------------
    
    integer             ::  ivar
    
    ! --- begin ----------------------------------
    
    ! grid setup ?
    if ( self%gridsetup ) then
      ! clear:
      deallocate( self%lons )
      deallocate( self%lats )
      deallocate( self%sigma )
      deallocate( self%hyam, self%hybm )
      deallocate( self%isigma )
      deallocate( self%hyai, self%hybi )
      ! reset:
      self%gridsetup = .false.
    end if
    
    ! buffer filled ?
    if ( self%filled ) then
      ! clear coordinate variables:
      deallocate( self%tt )
      deallocate( self%ps )
      ! variables:
      do ivar = 1, self%nvar
        deallocate( self%var(ivar)%data )
      end do
      ! reset:
      self%filled = .false.
    end if
    
    ! clear variable description:
    deallocate( self%var )
    
    ! ok
    status = 0
  
  end subroutine ModelOutputSeries_Done


  ! ***
  

  subroutine ModelOutputSeries_Get( self, status, &
                                     nlon, nlat, nlev, &
                                     lons, lats, &
                                     sigma, isigma, hyam, hybm, hyai, hybi, p0, &
                                     vname_ps, units_p, &
                                     label )
  
    ! --- in/out ---------------------------------
    
    class(T_ModelOutputSeries), intent(inout)   ::  self
    integer, intent(out)                        ::  status
    integer, intent(out), optional              ::  nlon, nlat, nlev
    real, intent(out), optional                 ::  lons(:), lats(:)
    real, intent(out), optional                 ::  sigma(:), isigma(:)
    real, intent(out), optional                 ::  hyam(:), hybm(:), hyai(:), hybi(:)
    real, intent(out), optional                 ::  p0
    character(len=*), intent(out), optional     ::  vname_ps
    character(len=*), intent(out), optional     ::  units_p
    character(len=*), intent(out), optional     ::  label

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelOutputSeries_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. self%filled ) then
      write (gol,'("buffer not filled yet ...")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! return values:
    if ( present(nlon    ) ) nlon    = self%nlon
    if ( present(nlat    ) ) nlat    = self%nlat
    if ( present(nlev    ) ) nlev    = self%nlev
    if ( present(lons    ) ) lons    = self%lons
    if ( present(lats    ) ) lats    = self%lats
    if ( present(sigma   ) ) sigma   = self%sigma
    if ( present(isigma  ) ) isigma  = self%isigma
    if ( present(hyam    ) ) hyam    = self%hyam
    if ( present(hybm    ) ) hybm    = self%hybm
    if ( present(hyai    ) ) hyai    = self%hyai
    if ( present(hybi    ) ) hybi    = self%hybi
    !if ( present(ptop    ) ) ptop    = self%ptop
    if ( present(p0      ) ) p0      = self%p0
    if ( present(vname_ps) ) vname_ps = self%vname_ps
    if ( present(units_p ) ) units_p  = self%units_p
    if ( present(label   ) ) label    = self%label
    
    ! ok
    status = 0
  
  end subroutine ModelOutputSeries_Get


  ! ***
  
  
  !
  ! Read file content for requested time.
  ! Here probably a month file, read complete month at once.
  !

  subroutine ModelOutputSeries_ReadBuffer( self, t, status )
  
    use GO, only : TDate, TIncrDate, Get, Extract_Ref_and_Step, operator(+), operator(*)
    use GO, only : goReplace, goReadFromLine

    use NetCDF, only : NF90_NOWRITE
    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_Open, NF90_Close
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_INQ_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_ModelOutputSeries), intent(inout)   ::  self
    type(TDate), intent(in)                     ::  t
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelOutputSeries_ReadBuffer'
    
    ! --- local ----------------------------------
    
    integer                 ::  yyyy, mm, dd, hh
    character(len=1024)     ::  filename
    logical                 ::  exist
    integer                 ::  ncid, dimid, varid
    character(len=32)       ::  vname
    character(len=64)       ::  standard_name
    character(len=64)       ::  formula_terms
    character(len=64)       ::  line
    character(len=32)       ::  key, value
    real(8), allocatable    ::  time(:)
    character(len=64)       ::  units
    type(TDate)             ::  tref
    type(TIncrDate)         ::  dtstep
    integer                 ::  itime
    integer                 ::  ivar
    
    ! --- begin ----------------------------------
    
    ! get time values:
    call Get( t, year=yyyy, month=mm, day=dd, hour=hh )
    
    ! filename, replace keys:
    filename = trim(self%filename_template)
    call goReplace( filename, '%{yyyy}', '(i4.4)', yyyy, status )
    IF_NOT_OK_RETURN(status=1)
    call goReplace( filename, '%{mm}'  , '(i2.2)', mm  , status )
    IF_NOT_OK_RETURN(status=1)
    call goReplace( filename, '%{dd}'  , '(i2.2)', dd  , status )
    IF_NOT_OK_RETURN(status=1)
    call goReplace( filename, '%{hh}'  , '(i2.2)', hh  , status )
    IF_NOT_OK_RETURN(status=1)
    
    ! different from curren ?
    if ( trim(filename) /= trim(self%filename_curr) ) then
    
      ! info ..
      write (gol,'("      read ",a," ...")') trim(filename); call goPr
    
      ! check ...
      inquire( file=trim(filename), exist=exist )
      if ( .not. exist ) then
        write (gol,'("file not found :")'); call goErr
        write (gol,'("  ",a)') trim(filename); call goErr
        TRACEBACK; status=1; return
      end if

      ! open file for reading:
      status = NF90_Open( trim(filename), NF90_NOWRITE, ncid )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! setup grid ?
      if ( .not. self%gridsetup ) then
    
        ! read label of this run:
        status = NF90_Get_Att( ncid, NF90_GLOBAL, 'run_label', self%label )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! info ..
        write (gol,'("        run label: ",a)') trim(self%label); call goPr
    
        ! info ..
        write (gol,'("        read grid ...")'); call goPr

        ! longitudes:
        vname = 'lon'
        ! access dimension:
        status = NF90_Inq_DimID( ncid, vname, dimid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! size:
        status = NF90_Inquire_Dimension( ncid, dimid, len=self%nlon )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! storage:
        allocate( self%lons(self%nlon), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! get variable id:
        status = NF90_INQ_VarID( ncid, vname, varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! get standard name:
        status = NF90_Get_Att( ncid, varid, 'standard_name', standard_name )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! check ...
        if ( trim(standard_name) /= 'longitude' ) then
          write (gol,'("unexpected standard name for longitude coordinate: ",a)') trim(standard_name); call goErr
          TRACEBACK; status=1; return
        end if
        ! read:
        status = NF90_Get_Var( ncid, varid, self%lons )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! latgitudes:
        vname = 'lat'
        ! access dimension:
        status = NF90_Inq_DimID( ncid, vname, dimid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! size:
        status = NF90_Inquire_Dimension( ncid, dimid, len=self%nlat )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! storage:
        allocate( self%lats(self%nlat), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! get variable id:
        status = NF90_INQ_VarID( ncid, vname, varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! get standard name:
        status = NF90_Get_Att( ncid, varid, 'standard_name', standard_name )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! check ...
        if ( trim(standard_name) /= 'latitude' ) then
          write (gol,'("unexpected standard name for latitude coordinate: ",a)') trim(standard_name); call goErr
          TRACEBACK; status=1; return
        end if
        ! read:
        status = NF90_Get_Var( ncid, varid, self%lats )
        IF_NF90_NOT_OK_RETURN(status=1)
        
        ! ~

        !! sigma values:
        !vname = 'k'
        !! access dimension:
        !status = NF90_Inq_DimID( ncid, vname, dimid )
        !IF_NF90_NOT_OK_RETURN(status=1)
        !! size:
        !status = NF90_Inquire_Dimension( ncid, dimid, len=self%nlev )
        !IF_NF90_NOT_OK_RETURN(status=1)
        !! storage:
        !allocate( self%sigma(self%nlev), stat=status )
        !IF_NOT_OK_RETURN(status=1)
        !! get variable id:
        !status = NF90_INQ_VarID( ncid, vname, varid )
        !IF_NF90_NOT_OK_RETURN(status=1)
        !! read values:
        !status = NF90_Get_Var( ncid, varid, self%sigma )
        !IF_NF90_NOT_OK_RETURN(status=1)
        !! get standard name:
        !status = NF90_Get_Att( ncid, varid, 'standard_name', standard_name )
        !IF_NF90_NOT_OK_RETURN(status=1)
        !! check ...
        !if ( trim(standard_name) /= 'atmosphere_sigma_coordinate' ) then
        !  write (gol,'("unexpected standard name for vertical coordinate: ",a)') trim(standard_name); call goErr
        !  TRACEBACK; status=1; return
        !end if
        !
        !! formula: 
        !!   p(n,k,j,i) = ptop + sigma(k)*(ps(n,j,i)-ptop)
        !! get formula terms description:
        !!   sigma: k ps: PS ptop: PT
        !status = NF90_Get_Att( ncid, varid, 'formula_terms', formula_terms )
        !IF_NF90_NOT_OK_RETURN(status=1)
        !! init dummy values:
        !self%vname_ps = ''
        !self%vname_ptop = ''
        !! loop:
        !do
        !  ! keyword:
        !  call goReadFromLine( formula_terms, key, status, sep=' ' )
        !  IF_NOT_OK_RETURN(status=1)
        !  ! value:
        !  call goReadFromLine( formula_terms, value, status, sep=' ' )
        !  IF_NOT_OK_RETURN(status=1)
        !  ! assign:
        !  select case ( trim(key) )
        !    ! weights:
        !    case ( 'sigma:' )
        !      ! should be this one ...
        !      if ( trim(value) /= trim(vname) ) then
        !        write (gol,'("sigma variable is named `",a,"` but current variable is : ",a)') trim(key), trim(vname); call goErr
        !        TRACEBACK; status=1; return
        !      end if
        !    ! pressures:
        !    case ( 'ps:' )
        !      self%vname_ps = trim(value)
        !    case ( 'ptop:' )
        !      self%vname_ptop = trim(value)
        !    ! unkown ...
        !    case default
        !      write (gol,'("unexpected formula_terms keyword : ",a)') trim(key); call goErr
        !      TRACEBACK; status=1; return
        !  end select
        !  ! end ?
        !  if ( len_trim(formula_terms) == 0 ) exit
        !end do ! formula terms
        !! check ...
        !if ( len_trim(self%vname_ps) == 0 ) then
        !  write (gol,'("no ps defined in formula terms : ",a)') trim(formula_terms); call goErr
        !  TRACEBACK; status=1; return
        !end if
        !! check ...
        !if ( len_trim(self%vname_ptop) == 0 ) then
        !  write (gol,'("no ptop defined in formula terms : ",a)') trim(formula_terms); call goErr
        !  TRACEBACK; status=1; return
        !end if
        !
        !! get variable id:
        !status = NF90_INQ_VarID( ncid, self%vname_ptop, varid )
        !IF_NF90_NOT_OK_RETURN(status=1)
        !! read values:
        !status = NF90_Get_Var( ncid, varid, self%ptop )
        !IF_NF90_NOT_OK_RETURN(status=1)
        !! units:
        !status = NF90_Get_Att( ncid, varid, 'units', self%units_p )
        !IF_NF90_NOT_OK_RETURN(status=1)
        
        ! ~

        ! hybride values:
        vname = 'lev'
        ! access dimension:
        status = NF90_Inq_DimID( ncid, vname, dimid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! size:
        status = NF90_Inquire_Dimension( ncid, dimid, len=self%nlev )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! storage:
        allocate( self%sigma(self%nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! variable:
        status = NF90_INQ_VarID( ncid, vname, varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! read values:
        status = NF90_Get_Var( ncid, varid, self%sigma )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! get standard name:
        status = NF90_Get_Att( ncid, varid, 'standard_name', standard_name )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! check ...
        if ( trim(standard_name) /= 'atmosphere_hybrid_sigma_pressure_coordinate' ) then
          write (gol,'("unexpected standard name for vertical coordinate: ",a)') trim(standard_name); call goErr
          TRACEBACK; status=1; return
        end if
        ! formula: 
        !   p(n,k,j,i) = ptop + sigma(k)*(ps(n,j,i)-ptop)
        ! get formula terms description:
		    !   ap: hyai b: hybi ps: PS p0: P0
        status = NF90_Get_Att( ncid, varid, 'formula_terms', formula_terms )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! init dummy values:
        self%vname_ap = ''
        self%vname_b  = ''
        self%vname_ps = ''
        self%vname_p0 = ''
        ! loop:
        line = trim(formula_terms)
        do
          ! keyword:
          call goReadFromLine( line, key, status, sep=' ' )
          IF_NOT_OK_RETURN(status=1)
          ! value:
          call goReadFromLine( line, value, status, sep=' ' )
          IF_NOT_OK_RETURN(status=1)
          ! assign:
          select case ( trim(key) )
            ! weights:
            case ( 'ap:' )
              self%vname_ap = trim(value)
            case ( 'b:' )
              self%vname_b = trim(value)
            ! pressures:
            case ( 'ps:' )
              self%vname_ps = trim(value)
            case ( 'p0:' )
              self%vname_p0 = trim(value)
            ! unkown ...
            case default
              write (gol,'("unexpected formula_terms keyword : ",a)') trim(key); call goErr
              TRACEBACK; status=1; return
          end select
          ! end ?
          if ( len_trim(line) == 0 ) exit
        end do ! formula terms
        ! check ...
        if ( len_trim(self%vname_ap) == 0 ) then
          write (gol,'("no ap defined in formula terms : ",a)') trim(formula_terms); call goErr
          TRACEBACK; status=1; return
        end if
        ! check ...
        if ( len_trim(self%vname_b) == 0 ) then
          write (gol,'("no b defined in formula terms : ",a)') trim(formula_terms); call goErr
          TRACEBACK; status=1; return
        end if
        ! check ...
        if ( len_trim(self%vname_ps) == 0 ) then
          write (gol,'("no ps defined in formula terms : ",a)') trim(formula_terms); call goErr
          TRACEBACK; status=1; return
        end if
        ! check ...
        if ( len_trim(self%vname_p0) == 0 ) then
          write (gol,'("no ptop defined in formula terms : ",a)') trim(formula_terms); call goErr
          TRACEBACK; status=1; return
        end if

        ! hybride values:
        vname = 'ilev'
        ! storage:
        allocate( self%isigma(self%nlev+1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! variable:
        status = NF90_INQ_VarID( ncid, vname, varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! read values:
        status = NF90_Get_Var( ncid, varid, self%isigma )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! get variable id:
        status = NF90_INQ_VarID( ncid, self%vname_p0, varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! read values:
        status = NF90_Get_Var( ncid, varid, self%p0 )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! units:
        status = NF90_Get_Att( ncid, varid, 'units', self%units_p )
        IF_NF90_NOT_OK_RETURN(status=1)

        ! storage:
        allocate( self%hyam(self%nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%hybm(self%nlev), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%hyai(self%nlev+1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%hybi(self%nlev+1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! read:
        status = NF90_INQ_VarID( ncid, 'hyam', varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Get_Var( ncid, varid, self%hyam )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! read:
        status = NF90_INQ_VarID( ncid, 'hybm', varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Get_Var( ncid, varid, self%hybm )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! read:
        status = NF90_INQ_VarID( ncid, 'hyai', varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Get_Var( ncid, varid, self%hyai )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! read:
        status = NF90_INQ_VarID( ncid, 'hybi', varid )
        IF_NF90_NOT_OK_RETURN(status=1)
        status = NF90_Get_Var( ncid, varid, self%hybi )
        IF_NF90_NOT_OK_RETURN(status=1)
        
        ! adhoc fixes ...
        self%hyai(self%nlev+1) = 0.0
        self%hybi(self%nlev+1) = 1.0
        
        ! ~

        ! reset:
        self%gridsetup = .true.
      end if
    
      ! info ..
      write (gol,'("        read times ...")'); call goPr
      ! time values:
      vname = 'time'
      ! access dimension:
      status = NF90_Inq_DimID( ncid, vname, dimid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! size:
      status = NF90_Inquire_Dimension( ncid, dimid, len=self%ntime )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! storage:
      allocate( time(self%ntime), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! get variable id:
      status = NF90_INQ_VarID( ncid, vname, varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! read values:
      status = NF90_Get_Var( ncid, varid, time )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! units:
      status = NF90_Get_Att( ncid, varid, 'units', units )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! extract content:
      !   days since 1900-1-1 0:0:0
      call Extract_Ref_and_Step( units, tref, dtstep, status )
      IF_NOT_OK_RETURN(status=1)

      ! clear buffer if necessary:
      if ( self%filled ) then
        ! info ..
        write (gol,'("        clear existing buffer ...")'); call goPr
        deallocate( self%tt )
        deallocate( self%ps )
      end if
      
      ! new storage:
      allocate( self%tt(self%ntime), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( self%ps(self%nlon,self%nlat,self%ntime), stat=status )
      IF_NOT_OK_RETURN(status=1)

      ! loop over time records:
      do itime = 1, self%ntime
        ! fill:
        self%tt(itime) = tref + dtstep * time(itime)
      end do
      ! clear:
      deallocate( time )
      
      ! get variable id:
      status = NF90_INQ_VarID( ncid, trim(self%vname_ps), varid )
      if (status/=NF90_NOERR) then
        ! info ..
        write (gol,'("        fill dummy surface pressure ...")'); call goPr
        !! dummy, not in files yet:
        !select case ( trim(self%units_p) )
        !  case ( 'hPa' )
        !    self%ps = 1000.0 ! hPa
        !  case ( 'Pa' )
        !    self%ps = 1000.0e2 ! hPa
        !  case default
        !    write (gol,'("unsupported pressure units : ",a)') trim(self%units_p); call goErr
        !    TRACEBACK; status=1; return
        !end select
        ! P0 is a reference surface pressure, use for entire field:
        self%ps = self%p0
      else
        ! read:
        status = NF90_Get_Var( ncid, varid, self%ps )
        IF_NF90_NOT_OK_RETURN(status=1)
      end if
      
      ! info ..
      write (gol,'("        read data ...")'); call goPr
      ! loop:
      do ivar = 1, self%nvar
        ! info ..
        write (gol,'("        var ",i0," (",a,")")') ivar, trim(self%var(ivar)%name); call goPr
        ! clear buffer if necessary:
        if ( self%filled ) then
          deallocate( self%var(ivar)%data )
        end if
        ! new storage:
        allocate( self%var(ivar)%data(self%nlon,self%nlat,self%nlev,self%ntime), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! data field:
        status = NF90_INQ_VarID( ncid, trim(self%var(ivar)%name), varid )
        if ( status /= NF90_NOERR ) then
          write (gol,'("could not find variable `",a,"` in file:")') trim(self%var(ivar)%name); call goErr
          write (gol,'("  ",a)') trim(filename); call goErr
          TRACEBACK; status=1; return
        end if
        ! read values:
        status = NF90_Get_Var( ncid, varid, self%var(ivar)%data )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! units:
        status = NF90_Get_Att( ncid, varid, 'units', units )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! convert to CF units if necessary:
        if ( trim(units) == 'mix_ratio' ) units = 'mol mol^-1'
        ! store:
        self%var(ivar)%units = trim(units)
      end do  ! var
      
      ! close file:
      status = NF90_Close( ncid )
      IF_NF90_NOT_OK_RETURN(status=1)   
      
      ! store current name:
      self%filename_curr = trim(filename)
      ! buffer is filled now:
      self%filled = .true.

    end if  ! new file ?

    ! ok
    status = 0
  
  end subroutine ModelOutputSeries_ReadBuffer


  ! ***
  

  subroutine ModelOutputSeries_GetRecord( self, t, ps, ivar, data, units, status )
  
    use GO, only : TDate, operator(==), wrtgol
  
    ! --- in/out ---------------------------------
    
    class(T_ModelOutputSeries), intent(inout)   ::  self
    type(TDate), intent(in)                     ::  t
    real, intent(out)                           ::  ps(:,:)
    integer, intent(in)                         ::  ivar
    real, intent(out)                           ::  data(:,:,:)
    character(len=*), intent(out)               ::  units
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/ModelOutputSeries_GetRecord'
    
    ! --- local ----------------------------------
    
    logical     ::  found
    integer     ::  itime
    
    ! --- begin ----------------------------------
    
    ! init flag:
    found = .false.
    ! search:
    do itime = 1, self%ntime
      ! match ?
      found = t == self%tt(itime)
      ! found ?
      if ( found ) exit
    end do
    ! check ..
    if ( .not. found ) then
      call wrtgol( 'time record not found: ', t ); call goErr
      write (gol,'("file        : ",a)') trim(self%filename_curr); call goErr
      write (gol,'("time records: ")'); call goErr
      do itime = 1, self%ntime
        call wrtgol( '  ', self%tt(itime) ); call goErr
      end do
      TRACEBACK; status=1; return
    end if
    
    ! extract:
    ps    = self%ps  (:,:,itime)
    
    ! check ...
    if ( (ivar < 1) .or. (ivar > self%nvar) ) then
      write (gol,'("variable index ",i4," outside expected range 1 .. ",i4)') ivar, self%nvar; call goErr
      TRACEBACK; status=1; return
    end if
    data  = self%var(ivar)%data(:,:,:,itime)
    units = trim(self%var(ivar)%units)
    
    ! ok
    status = 0
  
  end subroutine ModelOutputSeries_GetRecord


end module EMEP_NMC_ModelOutput


