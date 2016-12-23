!######################################################################
!
! EMEP DA NMC - output routines
!
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
!######################################################################

module EMEP_NMC_Output

  use GO    , only : gol, goPr, goErr
  use GO    , only : TDate
  use C3PO  , only : Datafile
  
  use NetCDF, only : NF90_NOERR, NF90_StrError

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_NMC_Output
  

  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_NMC_Output'
  
  ! fill values:
  integer(4), parameter  ::  fill_value_int    = huge(1_4)
  real(4), parameter     ::  fill_value_float  = huge(1.0_4)
  real(8), parameter     ::  fill_value_double = huge(1.0_8)
  

  ! --- types ----------------------------------------
  
  type, extends(Datafile) :: T_NMC_Output
    ! ref time:
    type(TDate)                 ::  t0
    ! access:
!    integer                     ::  dimid_nv
!    integer                     ::  dimid_ri
!    character(len=32)           ::  vname_tracer
!    integer                     ::  dimid_tracer, dimid_tracer_len
!    integer                     ::  varid_tracer
  contains
    procedure   ::  Init              => NMC_Output_Init
    procedure   ::  Done              => NMC_Output_Done
    procedure   ::  Create            => NMC_Output_Create
    procedure   ::  EndDef            => NMC_Output_EndDef
    procedure   ::  Extend_History    => NMC_Output_Extend_History
    procedure   ::  Get               => NMC_Output_Get
    procedure   ::  Get_VarID         => NMC_Output_Get_VarID
    !
    ! integer sclar
    procedure   ::  Def_IValue         => NMC_Output_Def_IValue
    procedure   ::  Put_IValue         => NMC_Output_Put_IValue
    procedure   ::  Get_IValue         => NMC_Output_Get_IValue
    !
    ! real (x,y)
    procedure   ::  Def_Field2D        => NMC_Output_Def_Field2D
    procedure   ::  Put_Field2D        => NMC_Output_Put_Field2D
    procedure   ::  Get_Field2D        => NMC_Output_Get_Field2D
    !
    ! integer (x)
    procedure   ::  Def_IField1D       => NMC_Output_Def_IField1D
    procedure   ::  Put_IField1D       => NMC_Output_Put_IField1D
    procedure   ::  Get_IField1D       => NMC_Output_Get_IField1D
    !
    ! integer (x,y)
    procedure   ::  Def_IField2D       => NMC_Output_Def_IField2D
    procedure   ::  Put_IField2D       => NMC_Output_Put_IField2D
    procedure   ::  Get_IField2D       => NMC_Output_Get_IField2D
    !
    ! integer (x,y,t), write per (x,y) slab
    procedure   ::  Def_IField2D_Series => NMC_Output_Def_IField2D_Series
    procedure   ::  Put_IField2D_Series => NMC_Output_Put_IField2D_Series
    procedure   ::  Get_IField2D_Series => NMC_Output_Get_IField2D_Series
    !
    ! real (x,y,t), write per (x,y) slab
    procedure   ::  Def_Field2D_Series => NMC_Output_Def_Field2D_Series
    procedure   ::  Put_Field2D_Series => NMC_Output_Put_Field2D_Series
    procedure   ::  Get_Field2D_Series => NMC_Output_Get_Field2D_Series
    !
    ! real (x,y,z,t), write per (x,y,z) slab
    procedure   ::  Def_Field3D_Series => NMC_Output_Def_Field3D_Series
    procedure   ::  Put_Field3D_Series => NMC_Output_Put_Field3D_Series
    procedure   ::  Get_Field3D_Series => NMC_Output_Get_Field3D_Series
    !
    ! real (x,y,z,z,label,label,t)
    procedure   ::  Def_Covar         => NMC_Output_Def_Covar
    procedure   ::  Put_Covar         => NMC_Output_Put_Covar
    procedure   ::  Get_Covar         => NMC_Output_Get_Covar
    !
    ! real (k,z,z,label,label,t)
    procedure   ::  Def_ACovar        => NMC_Output_Def_ACovar
    procedure   ::  Put_ACovar        => NMC_Output_Put_ACovar
    procedure   ::  Get_ACovar        => NMC_Output_Get_ACovar
    !
    ! real (z,z,label,label,t)
    procedure   ::  Def_ACovar1       => NMC_Output_Def_ACovar1
    procedure   ::  Put_ACovar1       => NMC_Output_Put_ACovar1
    procedure   ::  Get_ACovar1       => NMC_Output_Get_ACovar1
    !
    ! real (k,z,label,t)
    procedure   ::  Def_AVar          => NMC_Output_Def_AVar
    procedure   ::  Put_AVar          => NMC_Output_Put_AVar
    procedure   ::  Get_AVar          => NMC_Output_Get_AVar
    !
    ! real (k,t)
    procedure   ::  Def_Nev           => NMC_Output_Def_Nev
    procedure   ::  Put_Nev           => NMC_Output_Put_Nev
    procedure   ::  Get_Nev           => NMC_Output_Get_Nev
    !
    ! real (k,v,t)
    procedure   ::  Def_EVal          => NMC_Output_Def_EVal
    procedure   ::  Put_EVal          => NMC_Output_Put_EVal
    procedure   ::  Get_EVal          => NMC_Output_Get_EVal
    !
    ! real (k,z,label,v,t)
    procedure   ::  Def_EVec          => NMC_Output_Def_EVec
    procedure   ::  Put_EVec          => NMC_Output_Put_EVec
    procedure   ::  Get_EVec          => NMC_Output_Get_EVec
    !
    ! real (x,y,t,n), write per (x,y) slab
    procedure   ::  Def_Sample2D      => NMC_Output_Def_Sample_2d_r
    procedure   ::  Put_Sample2D      => NMC_Output_Put_Sample_2d_r
    procedure   ::  Get_Sample2D      => NMC_Output_Get_Sample_2d_r
    !
    ! real (x,y,z,t,n), write per (x,y,z) slab
    procedure   ::  Def_Sample3D      => NMC_Output_Def_Sample_3d_r
    procedure   ::  Put_Sample3D      => NMC_Output_Put_Sample_3d_r
    procedure   ::  Get_Sample3D      => NMC_Output_Get_Sample_3d_r
    !
    ! complex (x,y,z,t,n), write per (x,y,z) slab
    procedure   ::  Def_CSample3D     => NMC_Output_Def_Sample_3d_c
    procedure   ::  Put_CSample3D     => NMC_Output_Put_Sample_3d_c
    procedure   ::  Get_CSample3D     => NMC_Output_Get_Sample_3d_c
    !
    procedure   ::  Close             => NMC_Output_Close
    !
    procedure   ::                       NMC_Output_Put_Att_s
    procedure   ::                       NMC_Output_Put_Att_i
    procedure   ::                       NMC_Output_Put_Att_r
    generic     ::  Put_Att           => NMC_Output_Put_Att_s, &
                                         NMC_Output_Put_Att_i, &
                                         NMC_Output_Put_Att_r
    !
    procedure   ::  Get_Att           => NMC_Output_Get_Att_s
  end type T_NMC_Output
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** output file
  ! ***
  ! ********************************************************************


  subroutine NMC_Output_Init( self, filename, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(out)          ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Init
  
  
  ! ***


  subroutine NMC_Output_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Done
  
  
  ! ***

  
  ! 
  ! Arguments:
  !   fmt  :  netcdf format: 
  !             'netcdf' | 'classic'  : default
  !             'netcdf4'             : supports parallel i/o
  !

  subroutine NMC_Output_Create( self, filename, status, fmt )
  
    use NetCDF, only : NF90_Create
    use NetCDF, only : NF90_CLOBBER
    use NetCDF, only : NF90_NETCDF4, NF90_CLASSIC_MODEL
    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(out)          ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status
    
    character(len=*), intent(in), optional    ::  fmt

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Create'
    
    ! --- local ----------------------------------
    
    integer                 ::  cmode
    integer                 ::  cfmt
    
    ! --- begin ----------------------------------
    
    ! init main:
    call self%Init( filename, status )
    IF_NOT_OK_RETURN(status=1)

    ! init creation mode:
    cmode = 0
   
    ! overwite existing files:
    cmode = cmode + NF90_CLOBBER        ! overwrite existing
    
    ! default format:
    cfmt = NF90_CLASSIC_MODEL
    ! overwrite?
    if ( present(fmt) ) then
      ! switch:
      select case ( trim(fmt) )
        case ( 'netcdf', 'classic' )
          cfmt = NF90_CLASSIC_MODEL
        case ( 'netcdf4' )
          cfmt = NF90_NETCDF4
        case default
          write (gol,'("unsupported fmt `",a,"`")') trim(fmt); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    ! add:
    cmode = cmode + cfmt
    
    ! open new file:
    status = NF90_Create( trim(self%filename), cmode, self%ncid )
    if (status/=NF90_NOERR) then
      write (gol,'("creating: ",a)') trim(self%filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! some global attributes:
    status = NF90_Put_Att( self%ncid, NF90_GLOBAL, 'Conventions', 'CF-1.6' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, NF90_GLOBAL, 'title', 'Statistics of NMC output.' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, NF90_GLOBAL, 'institution', 'MetNorway, Oslo' )
    IF_NF90_NOT_OK_RETURN(status=1)
    status = NF90_Put_Att( self%ncid, NF90_GLOBAL, 'source', 'EMEP_MSC-W' )
    IF_NF90_NOT_OK_RETURN(status=1)
    
!    ! no dims defined yet ...
!    self%dimid_nv = -999
!    self%dimid_ri = -999

    ! ok
    status = 0
  
  end subroutine NMC_Output_Create


  ! ***
  

  subroutine NMC_Output_Close( self, status )

    use NetCDF, only : NF90_Close
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)          ::  self
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Close'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    status = NF90_Close( self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1) 
    
    ! end main:
    call self%Done( status )
    IF_NOT_OK_RETURN(status=1)

    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Close


  ! ***


  subroutine NMC_Output_EndDef( self, status )
  
    use NetCDF, only : NF90_EndDef
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_EndDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! end of definition phase:
    status = NF90_EndDef( self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1) 

    ! ok
    status = 0
  
  end subroutine NMC_Output_EndDef
  
  
  ! ***
  
  
  !
  ! Prepend line to history:
  !   Wed Apr 23 16:32:03 2014: Message\n
  !
  ! Only allowed in definition mode.
  !
  
  subroutine NMC_Output_Extend_History( self, message, status )

    use NetCDF, only : NF90_GLOBAL
    use NetCDF, only : NF90_Get_Att, NF90_Put_Att
    
    use GO, only : TDate, SystemDate
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)  ::  self
    character(len=*), intent(in)        ::  message
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Extend_History'
    
    ! month names:
    character(len=3),parameter :: monthname(12) = &
                         (/'Jan','Feb','Mrc','Apr','May','Jun',&
                           'Jul','Aug','Sep','Oct','Nov','Dec'/)
   
    ! newline character:
    character(len=1), parameter  ::  newline = char(10)
     
    ! --- local ----------------------------------
    
    character(len=64)     ::  attname
    character(len=64)     ::  tstamp
    character(len=1024)   ::  history
    type(TDate)           ::  t
    
    ! --- begin ----------------------------------
    
    ! target attribute:
    attname = 'history'
    
    ! current:
    status = NF90_Get_Att( self%ncid, NF90_GLOBAL, attname, history )
    if ( status /= NF90_NOERR ) history = ''
    
    ! current time:
    t = SystemDate()
    ! write time stamp:
    write (tstamp,'("Day ",a," ",i2," ",i2,2(":",i2.2)," ",i4)') &
             monthname(t%month), t%day, t%hour, t%min, t%sec, t%year
    ! fill history:
    if ( len_trim(history)== 0 ) then
      write (history,'(a,": ",a,a,a)') trim(tstamp), trim(message)
    else
      write (history,'(a,": ",a,a,a)') trim(tstamp), trim(message), &
                                               newline, trim(history)
    end if
    
    ! write:
    status = NF90_Put_Att( self%ncid, NF90_GLOBAL, attname, history )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine NMC_Output_Extend_History


  ! ***


  subroutine NMC_Output_Get( self, status )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Get


  ! ***


  subroutine NMC_Output_Get_VarID( self, varname, varid, status )
  
    use NetCDF, only : NF90_Inq_VarID
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  varname
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_VarID'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(varname), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_VarID


  ! ***


!  subroutine NMC_Output_Def_Tracer( self, ntracer, tracer_len, status )
!  
!    use GO, only : TDate, NewDate, Get
!  
!    use NetCDF, only : NF90_CHAR
!    use NetCDF, only : NF90_Def_Dim
!    use NetCDF, only : NF90_Def_Var
!    use NetCDF, only : NF90_Put_Att
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_NMC_Output), intent(inout)        ::  self
!    integer, intent(in)                       ::  ntracer, tracer_len
!    integer, intent(out)                      ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Tracer'
!    
!    ! --- local ----------------------------------
!    
!    integer                 ::  varid
!    
!    ! --- begin ----------------------------------
!    
!    ! size:
!    status = NF90_Def_Dim( self%ncid, trim(self%vname_tracer), ntracer, self%dimid_tracer )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! size:
!    status = NF90_Def_Dim( self%ncid, trim(self%vname_tracer)//'_len', tracer_len, self%dimid_tracer_len )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    
!    ! coordinate:
!    status = NF90_Def_Var( self%ncid, trim(self%vname_tracer), NF90_CHAR, &
!                           (/self%dimid_tracer_len,self%dimid_tracer/), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! annote:
!    status = NF90_Put_Att( self%ncid, varid, 'long_name', 'tracer' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    status = NF90_Put_Att( self%ncid, varid, 'units', '1' )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! store:
!    self%varid_tracer = varid
!
!    ! ok
!    status = 0
!  
!  end subroutine NMC_Output_Def_Tracer
!
!
!  ! ***
!
!
!  subroutine NMC_Output_Put_Tracer( self, itracer, name, status )
!  
!    use GO    , only : TDate, operator(-), rTotal
!    use NetCDF, only : NF90_Put_Var
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_NMC_Output), intent(inout)        ::  self
!    integer, intent(in)                       ::  itracer
!    character(len=*), intent(in)              ::  name
!    integer, intent(out)                      ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Tracer'
!    
!    ! --- local ----------------------------------
!    
!    ! --- begin ----------------------------------
!    
!    ! write:
!    status = NF90_Put_Var( self%ncid, self%varid_tracer, name, &
!                             start=(/1,itracer/), count=(/len(name),1/) )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! ok
!    status = 0
!  
!  end subroutine NMC_Output_Put_Tracer
!
!
!  ! ***
!
!
!  subroutine NMC_Output_Get_Tracer( self, names, status )
!  
!    use GO    , only : TDate, Num_to_Date
!    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var, NF90_Get_Att
!  
!    ! --- in/out ---------------------------------
!    
!    class(T_NMC_Output), intent(inout)        ::  self
!    character(len=*), intent(out)             ::  names(:)
!    integer, intent(out)                      ::  status
!
!    ! --- const ----------------------------------
!
!    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Tracer'
!    
!    ! --- local ----------------------------------
!    
!    integer             ::  varid
!    
!    ! --- begin ----------------------------------
!    
!    ! access variable:
!    status = NF90_Inq_VarID( self%ncid, trim(self%vname_tracer), varid )
!    IF_NF90_NOT_OK_RETURN(status=1)
!    ! read:
!    status = NF90_Get_Var( self%ncid, varid, names )
!    IF_NF90_NOT_OK_RETURN(status=1)
!
!    ! ok
!    status = 0
!  
!  end subroutine NMC_Output_Get_Tracer


  ! ***


  subroutine NMC_Output_Def_Field2D( self, name, units, lons, lats, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  lons, lats
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Field2D'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lon, dimid_lat
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call lons%Get_Dim( status, dimid=dimid_lon )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats%Get_Dim( status, dimid=dimid_lat )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
                             (/dimid_lon,dimid_lat/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_Field2D


  ! ***


  subroutine NMC_Output_Put_Field2D( self, varid, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    real, intent(in)                          ::  field(:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Field2D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Field2D


  ! ***


  subroutine NMC_Output_Get_Field2D( self, name, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    real, intent(out)                         ::  field(:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Field2D'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field  )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Field2D


  ! ***


  subroutine NMC_Output_Def_IValue( self, name, units, &
                                     varid, status )
  
    use NetCDF, only : NF90_INT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_IValue'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! variable, no dimensions:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_INT, varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_IValue


  ! ***


  subroutine NMC_Output_Put_IValue( self, varid, value, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  value
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_IValue'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, value )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_IValue


  ! ***


  subroutine NMC_Output_Get_IValue( self, name, value, units, status )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  value
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_IValue'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, value  )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_IValue


  ! ***


  subroutine NMC_Output_Def_IField1D( self, name, units, dims, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_INT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  dims
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    integer, intent(out), optional            ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_IField1D'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call dims%Get_Dim( status, dimid=dimid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_INT, &
                             (/dimid/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_int
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_int )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_IField1D


  ! ***


  subroutine NMC_Output_Put_IField1D( self, varid, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  field(:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_IField1D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_IField1D


  ! ***


  subroutine NMC_Output_Get_IField1D( self, name, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  field(:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    integer, intent(out), optional            ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_IField1D'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field  )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_IField1D


  ! ***


  subroutine NMC_Output_Def_IField2D( self, name, units, lons, lats, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_INT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  lons, lats
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    integer, intent(out), optional            ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_IField2D'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lon, dimid_lat
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call lons%Get_Dim( status, dimid=dimid_lon )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats%Get_Dim( status, dimid=dimid_lat )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_INT, &
                             (/dimid_lon,dimid_lat/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_int
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_int )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_IField2D


  ! ***


  subroutine NMC_Output_Put_IField2D( self, varid, field, status, irec )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  field(:,:)
    integer, intent(out)                      ::  status

    integer, intent(in), optional             ::  irec  ! index of 2nd dim

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_IField2D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! record or all ?
    if ( present(irec) ) then
      ! write full field:
      status = NF90_Put_Var( self%ncid, varid, field, &
                         start=(/1,irec/), count=(/size(field,1),1/) )
      IF_NF90_NOT_OK_RETURN(status=1)
    else
      ! write full field:
      status = NF90_Put_Var( self%ncid, varid, field )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_IField2D


  ! ***


  subroutine NMC_Output_Get_IField2D( self, name, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  field(:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    integer, intent(out), optional            ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_IField2D'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field  )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_IField2D


  ! ***


  subroutine NMC_Output_Def_IField2D_Series( self, name, units, lons, lats, times, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_INT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  lons, lats, times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    integer, intent(out), optional            ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_IField2D_Series'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lon, dimid_lat, dimid_time
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call lons%Get_Dim( status, dimid=dimid_lon )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats%Get_Dim( status, dimid=dimid_lat )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_INT, &
                             (/dimid_lon,dimid_lat,dimid_time/), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_int
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_int )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_IField2D_Series


  ! ***


  subroutine NMC_Output_Put_IField2D_Series( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    integer, intent(in)                       ::  field(:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_IField2D_Series'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write record:
    status = NF90_Put_Var( self%ncid, varid, field, &
                       start=(/1,1,itime/), count=(/size(field,1),size(field,2),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_IField2D_Series


  ! ***


  subroutine NMC_Output_Get_IField2D_Series( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    integer, intent(out)                      ::  field(:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    integer, intent(out), optional            ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_IField2D_Series'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
                            start=(/1,1,itime/), count=(/size(field,1),size(field,2),1/)  )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_IField2D_Series


  ! ***


  subroutine NMC_Output_Def_Field2D_Series( self, name, units, &
                                              lons, lats, times, &
                                              varid, status, &
                                              fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  lons, lats, times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Field2D_Series'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lon, dimid_lat
    integer                 ::  dimid_time
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call lons%Get_Dim( status, dimid=dimid_lon )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats%Get_Dim( status, dimid=dimid_lat )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
               (/dimid_lon,dimid_lat,dimid_time/), &
               varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_Field2D_Series


  ! ***


  subroutine NMC_Output_Put_Field2D_Series( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    real, intent(in)                          ::  field(:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Field2D_Series'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
                             start=(/1,1,itime/), &
                             count=(/size(field,1),size(field,2),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Field2D_Series


  ! ***


  subroutine NMC_Output_Get_Field2D_Series( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    real, intent(out)                         ::  field(:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Field2D_Series'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
                             start=(/1,1,itime/), &
                             count=(/size(field,1),size(field,2),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Field2D_Series


  ! ***


  subroutine NMC_Output_Def_Field3D_Series( self, name, units, &
                                       lons, lats, levs, times, &
                                       varid, status, &
                                       fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  lons, lats, levs, times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Field3D_Series'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lon, dimid_lat, dimid_lev, dimid_time
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call lons%Get_Dim( status, dimid=dimid_lon )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats%Get_Dim( status, dimid=dimid_lat )
    IF_NF90_NOT_OK_RETURN(status=1)
    call levs%Get_Dim( status, dimid=dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
               (/dimid_lon,dimid_lat,dimid_lev,dimid_time/), &
               varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_Field3D_Series


  ! ***


  subroutine NMC_Output_Put_Field3D_Series( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    real, intent(in)                          ::  field(:,:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Field3D_Series'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
                             start=(/1,1,1,itime/), &
                             count=(/size(field,1),size(field,2),size(field,3),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Field3D_Series


  ! ***


  subroutine NMC_Output_Get_Field3D_Series( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    real, intent(out)                         ::  field(:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Field3D_Series'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
                             start=(/1,1,1,itime/), &
                             count=(/size(field,1),size(field,2),size(field,3),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Field3D_Series


  ! ***


  subroutine NMC_Output_Def_Covar( self, name, units, &
                                     lons_f, lats_f, levs, tracers, times, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension, LabelCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  lons_f, lats_f, levs
    class(LabelCoordinate), intent(in)        ::  tracers
    class(Dimension), intent(in)              ::  times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Covar'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lon_f, dimid_lat_f
    integer                 ::  dimid_lev, dimid_time
    integer                 ::  dimid_tracer
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call lons_f%Get_Dim( status, dimid=dimid_lon_f )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats_f%Get_Dim( status, dimid=dimid_lat_f )
    IF_NF90_NOT_OK_RETURN(status=1)
    call levs%Get_Dim( status, dimid=dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call tracers%Get_Dim( status, dimid=dimid_tracer )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
               (/dimid_lon_f,dimid_lat_f,&
                   dimid_lev,dimid_lev,dimid_tracer,dimid_tracer, &
                   dimid_time/), &
               varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'coordinates', trim(tracers%name)//'_name' )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_Covar


  ! ***


  subroutine NMC_Output_Put_Covar( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    real, intent(in)                          ::  field(:,:,:,:,:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Covar'
    
    ! --- local ----------------------------------
    
    real(4), allocatable   ::  field_r4(:,:,:,:,:,:)
    
    ! --- begin ----------------------------------

    ! convert kind, otherwise errors from netcdf routine:
    allocate( field_r4(size(field,1),size(field,2),size(field,3),size(field,4),size(field,5),size(field,6)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! convert:
    field_r4 = real( field, 4 )
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field_r4, &
               start=(/1,1,1,1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),size(field,4),size(field,5),size(field,6),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( field_r4, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Covar


  ! ***


  subroutine NMC_Output_Get_Covar( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    real, intent(out)                         ::  field(:,:,:,:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Covar'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
               start=(/1,1,1,1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),size(field,4),size(field,5),size(field,6),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Covar


  ! ***


  subroutine NMC_Output_Def_ACovar( self, name, units, &
                                     kstars, levs, tracers, times, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension, LabelCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  kstars
    class(Dimension), intent(in)              ::  levs
    class(LabelCoordinate), intent(in)        ::  tracers
    class(Dimension), intent(in)              ::  times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_ACovar'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_kstar
    integer                 ::  dimid_lev, dimid_time
    integer                 ::  dimid_tracer
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call kstars%Get_Dim( status, dimid=dimid_kstar )
    IF_NF90_NOT_OK_RETURN(status=1)
    call levs%Get_Dim( status, dimid=dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call tracers%Get_Dim( status, dimid=dimid_tracer )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
                 (/dimid_kstar,&
                   dimid_lev,dimid_lev,dimid_tracer,dimid_tracer, &
                   dimid_time/), &
               varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'coordinates', trim(tracers%name)//'_name' )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_ACovar


  ! ***


  subroutine NMC_Output_Put_ACovar( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    real, intent(in)                          ::  field(:,:,:,:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_ACovar'
    
    ! --- local ----------------------------------
    
    real(4), allocatable   ::  field_r4(:,:,:,:,:)
    
    ! --- begin ----------------------------------

    ! convert kind, otherwise errors from netcdf routine:
    allocate( field_r4(size(field,1),size(field,2),size(field,3),size(field,4),size(field,5)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! convert:
    field_r4 = real( field, 4 )
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field_r4, &
               start=(/1,1,1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),size(field,4),size(field,5),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! clear:
    deallocate( field_r4, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_ACovar


  ! ***


  subroutine NMC_Output_Get_ACovar( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    real, intent(out)                         ::  field(:,:,:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_ACovar'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
               start=(/1,1,1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),size(field,4),size(field,5),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_ACovar


  ! ***


  subroutine NMC_Output_Def_ACovar1( self, name, units, &
                                     levs, tracers, times, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension, LabelCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  levs
    class(LabelCoordinate), intent(in)        ::  tracers
    class(Dimension), intent(in)              ::  times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_ACovar1'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lev, dimid_time
    integer                 ::  dimid_tracer
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call levs%Get_Dim( status, dimid=dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call tracers%Get_Dim( status, dimid=dimid_tracer )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
                 (/dimid_lev,dimid_lev,dimid_tracer,dimid_tracer, &
                   dimid_time/), &
               varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'coordinates', trim(tracers%name)//'_name' )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_ACovar1


  ! ***


  subroutine NMC_Output_Put_ACovar1( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    real, intent(in)                          ::  field(:,:,:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_ACovar1'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
               start=(/1,1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),size(field,4),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_ACovar1


  ! ***


  subroutine NMC_Output_Get_ACovar1( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    real, intent(out)                         ::  field(:,:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_ACovar1'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
               start=(/1,1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),size(field,4),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_ACovar1


  ! ***


  subroutine NMC_Output_Def_AVar( self, name, units, &
                                     kstars, levs, tracers, times, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension, LabelCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  kstars
    class(Dimension), intent(in)              ::  levs
    class(LabelCoordinate), intent(in)        ::  tracers
    class(Dimension), intent(in)              ::  times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_AVar'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_kstar
    integer                 ::  dimid_lev, dimid_time
    integer                 ::  dimid_tracer
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call kstars%Get_Dim( status, dimid=dimid_kstar )
    IF_NF90_NOT_OK_RETURN(status=1)
    call levs%Get_Dim( status, dimid=dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call tracers%Get_Dim( status, dimid=dimid_tracer )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
                 (/dimid_kstar,dimid_lev,dimid_tracer, dimid_time/), &
                  varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'coordinates', trim(tracers%name)//'_name' )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_AVar


  ! ***


  subroutine NMC_Output_Put_AVar( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    real, intent(in)                          ::  field(:,:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_AVar'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
               start=(/1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_AVar


  ! ***


  subroutine NMC_Output_Get_AVar( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    real, intent(out)                         ::  field(:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_AVar'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
               start=(/1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_AVar


  ! ***


  subroutine NMC_Output_Def_Nev( self, name, units, &
                                     kstars, times, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_INT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  kstars
    class(Dimension), intent(in)              ::  times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Nev'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_kstar
    integer                 ::  dimid_time
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call kstars%Get_Dim( status, dimid=dimid_kstar )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_INT, &
                 (/dimid_kstar,dimid_time/), &
                  varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_int
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_int )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_Nev


  ! ***


  subroutine NMC_Output_Put_Nev( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    integer, intent(in)                       ::  field(:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Nev'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
               start=(/1,itime/), &
               count=(/size(field),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Nev


  ! ***


  subroutine NMC_Output_Get_Nev( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    integer, intent(out)                      ::  field(:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Nev'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
               start=(/1,itime/), &
               count=(/size(field),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Nev


  ! ***


  subroutine NMC_Output_Def_EVal( self, name, units, &
                                     kstars, evals, times, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  kstars
    class(Dimension), intent(in)              ::  evals
    class(Dimension), intent(in)              ::  times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_EVal'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_kstar
    integer                 ::  dimid_time
    integer                 ::  dimid_eval
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call kstars%Get_Dim( status, dimid=dimid_kstar )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call evals%Get_Dim( status, dimid=dimid_eval )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
                 (/dimid_kstar,dimid_eval,dimid_time/), &
                  varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_EVal


  ! ***


  subroutine NMC_Output_Put_EVal( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    real, intent(in)                          ::  field(:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_EVal'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
               start=(/1,1,itime/), &
               count=(/size(field,1),size(field,2),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_EVal


  ! ***


  subroutine NMC_Output_Get_EVal( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    real, intent(out)                         ::  field(:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_EVal'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
               start=(/1,1,itime/), &
               count=(/size(field,1),size(field,2),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_EVal


  ! ***


  subroutine NMC_Output_Def_EVec( self, name, units, &
                                     kstars, levs, tracers, evals, times, &
                                     varid, status, &
                                     fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension, LabelCoordinate
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  kstars
    class(Dimension), intent(in)              ::  levs
    class(LabelCoordinate), intent(in)        ::  tracers
    class(Dimension), intent(in)              ::  evals
    class(Dimension), intent(in)              ::  times
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_EVec'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_kstar
    integer                 ::  dimid_lev, dimid_time
    integer                 ::  dimid_tracer
    integer                 ::  dimid_eval
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call kstars%Get_Dim( status, dimid=dimid_kstar )
    IF_NF90_NOT_OK_RETURN(status=1)
    call levs%Get_Dim( status, dimid=dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call tracers%Get_Dim( status, dimid=dimid_tracer )
    IF_NF90_NOT_OK_RETURN(status=1)
    call evals%Get_Dim( status, dimid=dimid_eval )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
                 (/dimid_kstar,dimid_lev,dimid_tracer,dimid_eval,dimid_time/), &
                  varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'coordinates', trim(tracers%name)//'_name' )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_EVec


  ! ***


  subroutine NMC_Output_Put_EVec( self, varid, itime, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    real, intent(in)                          ::  field(:,:,:,:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_EVec'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
               start=(/1,1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),size(field,4),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_EVec


  ! ***


  subroutine NMC_Output_Get_EVec( self, name, itime, field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    real, intent(out)                         ::  field(:,:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_EVec'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
               start=(/1,1,1,1,itime/), &
               count=(/size(field,1),size(field,2),size(field,3),size(field,4),1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_EVec


  ! ***


  subroutine NMC_Output_Def_Sample_2d_r( self, name, units, &
                                          lons, lats, times, samples, &
                                          varid, status )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  lons, lats, times, samples
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Sample_2d_r'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lon, dimid_lat
    integer                 ::  dimid_time
    integer                 ::  dimid_sample
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call lons%Get_Dim( status, dimid=dimid_lon )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats%Get_Dim( status, dimid=dimid_lat )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call samples%Get_Dim( status, dimid=dimid_sample )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
               (/dimid_lon,dimid_lat,dimid_time,dimid_sample/), &
               varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_Sample_2d_r


  ! ***


  subroutine NMC_Output_Put_Sample_2d_r( self, varid, itime, isample, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    integer, intent(in)                       ::  isample
    real, intent(in)                          ::  field(:,:)  ! (nlon,nlat)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Sample_2d_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
                             start=(/1,1,itime,isample/), &
                             count=(/size(field,1),size(field,2),1,1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Sample_2d_r


  ! ***


  subroutine NMC_Output_Get_Sample_2d_r( self, name, itime, isample, field, units, status )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    integer, intent(in)                       ::  isample
    real, intent(out)                         ::  field(:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Sample_2d_r'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
                             start=(/1,1,itime,isample/), &
                             count=(/size(field,1),size(field,2),1,1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Sample_2d_r


  ! ***


  subroutine NMC_Output_Def_Sample_3d_r( self, name, units, &
                                          lons, lats, levs, times, samples, &
                                          varid, status )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  lons, lats, levs, times, samples
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Sample_3d_r'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_lon, dimid_lat, dimid_lev
    integer                 ::  dimid_time
    integer                 ::  dimid_sample
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call lons%Get_Dim( status, dimid=dimid_lon )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats%Get_Dim( status, dimid=dimid_lat )
    IF_NF90_NOT_OK_RETURN(status=1)
    call levs%Get_Dim( status, dimid=dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call samples%Get_Dim( status, dimid=dimid_sample )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
               (/dimid_lon,dimid_lat,dimid_lev,dimid_time,dimid_sample/), &
               varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_Sample_3d_r


  ! ***


  subroutine NMC_Output_Put_Sample_3d_r( self, varid, itime, isample, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    integer, intent(in)                       ::  isample
    real, intent(in)                          ::  field(:,:,:)  ! (nlon,nlat,nlev)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Sample_3d_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! write:
    status = NF90_Put_Var( self%ncid, varid, field, &
                             start=(/1,1,1,itime,isample/), &
                             count=(/size(field,1),size(field,2),size(field,3),1,1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Sample_3d_r


  ! ***


  subroutine NMC_Output_Get_Sample_3d_r( self, name, itime, isample, field, units, status )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    integer, intent(in)                       ::  isample
    real, intent(out)                         ::  field(:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Sample_3d_r'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, &
                             start=(/1,1,1,itime,isample/), &
                             count=(/size(field,1),size(field,2),size(field,3),1,1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Sample_3d_r


  ! ***


  subroutine NMC_Output_Def_Sample_3d_c( self, name, units, &
                                           cplx, lons_f, lats_f, levs, times, samples, &
                                           varid, status, &
                                           fill_value )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO, only : Dimension
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  units
    class(Dimension), intent(in)              ::  cplx, lons_f, lats_f, levs, times, samples
    integer, intent(out)                      ::  varid
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Def_Sample_3d_c'
    
    ! --- local ----------------------------------
    
    integer                 ::  dimid_cplx
    integer                 ::  dimid_lon_f, dimid_lat_f
    integer                 ::  dimid_lev
    integer                 ::  dimid_time
    integer                 ::  dimid_sample
    
    ! --- begin ----------------------------------
    
    ! dimids:
    call cplx%Get_Dim( status, dimid=dimid_cplx )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lons_f%Get_Dim( status, dimid=dimid_lon_f )
    IF_NF90_NOT_OK_RETURN(status=1)
    call lats_f%Get_Dim( status, dimid=dimid_lat_f )
    IF_NF90_NOT_OK_RETURN(status=1)
    call levs%Get_Dim( status, dimid=dimid_lev )
    IF_NF90_NOT_OK_RETURN(status=1)
    call times%Get_Dim( status, dimid=dimid_time )
    IF_NF90_NOT_OK_RETURN(status=1)
    call samples%Get_Dim( status, dimid=dimid_sample )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! coordinate:
    status = NF90_Def_Var( self%ncid, trim(name), NF90_FLOAT, &
               (/dimid_cplx,dimid_lon_f,dimid_lat_f,dimid_lev,dimid_time,dimid_sample/), &
               varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! annote:
    status = NF90_Put_Att( self%ncid, varid, 'units', trim(units) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! no-data value:
    if ( present(fill_value) ) then
      ! set suitable value for datatype of output variable:
      fill_value = fill_value_float
      ! write as attribute:
      status = NF90_Put_Att( self%ncid, varid, '_FillValue', fill_value_float )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
  
  end subroutine NMC_Output_Def_Sample_3d_c


  ! ***


  subroutine NMC_Output_Put_Sample_3d_c( self, varid, itime, isample, field, status )
  
    use NetCDF, only : NF90_Put_Var
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    integer, intent(in)                       ::  itime
    integer, intent(in)                       ::  isample
    complex, intent(in)                       ::  field(:,:,:)  ! (nlon,nlat,nlev)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Sample_3d_c'
    
    ! output kind for reals:
    integer, parameter         ::  rknd = 4
    
    ! --- local ----------------------------------
    
    integer                    ::  nlon, nlat, nlev
    real(rknd), allocatable    ::  values(:,:,:)
    
    ! --- begin ----------------------------------
    
    ! dims:
    nlon = size(field,1)
    nlat = size(field,2)
    nlev = size(field,3)
    
    ! storage for real/imag parts:
    allocate( values(nlon,nlat,nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! extract real part, convert kind:
    values = real(real(field),kind=rknd)
    ! write real part:
    status = NF90_Put_Var( self%ncid, varid, values, &
                             start=(/1,1,1,1,itime,isample/), &
                             count=(/1,nlon,nlat,nlev,1,1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! extract imag part, convert kind:
    values = real(aimag(field),kind=rknd)
    ! write imag part:
    status = NF90_Put_Var( self%ncid, varid, values, &
                             start=(/2,1,1,1,itime,isample/), &
                             count=(/1,nlon,nlat,nlev,1,1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( values, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Sample_3d_c


  ! ***


  subroutine NMC_Output_Get_Sample_3d_c( self, name, itime, isample, &
                                          field, units, status, &
                                          fill_value )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  itime
    integer, intent(in)                       ::  isample
    complex, intent(out)                      ::  field(:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Sample_3d_c'
    
    ! --- local ----------------------------------
    
    integer             ::  varid
    integer             ::  nlon, nlat, nlev
    real, allocatable   ::  values(:,:,:,:)
    
    ! --- begin ----------------------------------
    
    ! dims:
    nlon = size(field,1)
    nlat = size(field,2)
    nlev = size(field,3)
    
    ! storage:
    allocate( values(2,nlon,nlat,nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( self%ncid, varid, values, &
                             start=(/1,1,1,1,itime,isample/), &
                             count=(/2,nlon,nlat,nlev,1,1/) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! store:
    field = cmplx( values(1,:,:,:), values(2,:,:,:) )

    ! annote:
    status = NF90_Get_Att( self%ncid, varid, 'units', units )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! no-data value:
    if ( present(fill_value) ) then
      status = NF90_Get_Att( self%ncid, varid, '_FillValue', fill_value )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! clear:
    deallocate( values, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Sample_3d_c


  ! ***


  subroutine NMC_Output_Put_Att_s( self, varid, name, value, status )
  
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  value
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Att_s'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! put attribute:
    status = NF90_Put_Att( self%ncid, varid, trim(name), trim(value) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Att_s


  ! ***


  subroutine NMC_Output_Get_Att_s( self, varid, name, value, status )
  
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    character(len=*), intent(in)              ::  name
    character(len=*), intent(out)             ::  value
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Get_Att_s'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! put attribute:
    status = NF90_Get_Att( self%ncid, varid, trim(name), value )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Get_Att_s


  ! ***


  subroutine NMC_Output_Put_Att_i( self, varid, name, value, status )
  
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    character(len=*), intent(in)              ::  name
    integer, intent(in)                       ::  value
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Att_i'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! put attribute:
    status = NF90_Put_Att( self%ncid, varid, trim(name), value )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Att_i


  ! ***


  subroutine NMC_Output_Put_Att_r( self, varid, name, value, status )
  
    use NetCDF, only : NF90_Put_Att
  
    ! --- in/out ---------------------------------
    
    class(T_NMC_Output), intent(inout)        ::  self
    integer, intent(in)                       ::  varid
    character(len=*), intent(in)              ::  name
    real, intent(in)                          ::  value
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/NMC_Output_Put_Att_r'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! put attribute:
    status = NF90_Put_Att( self%ncid, varid, trim(name), value )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine NMC_Output_Put_Att_r


  ! ***
  

end module EMEP_NMC_Output


