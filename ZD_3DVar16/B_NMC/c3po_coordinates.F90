!######################################################################
!
! C3PO - CF Convention Compliance Program Objects
!
!
! CLASS HIERARCHY
!
!   The classes defined here have the following hierarchy:
!
!     Dimension
!       Coordinate
!         RealCoordinate
!         IntegerCoordinate
!         LabelCoordinate
!         HybrideLevelCoordinate
!         TimeCoordinate
!
!
! DIMENSIONS
!
!   A dimension consists of a name and associated length.
!
!   Usage:
!
!     ! import:
!     use C3PO_Coordinates, only : Dimension
!
!     ! variables:
!     type(Dimension)    ::  dim
!     integer              ::  n
!     integer              ::  ncid
!     integer              ::  status
!
!     ! ~ create new dimension
!
!     ! init:
!     call dim%Init( 'value', status )
!     if (status/=0) stop
!
!     ! fill:
!     n = 10
!
!     ! set length:
!     call dim%Set_Dim( status, n=n )
!     if (status/=0) stop
!
!     ! create new netcdf file, obtain file id:
!     ncid = ...
!
!     ! define in file:
!     call dim%Def( ncid, status )
!     if (status/=0) stop
!
!     ! end of definition phase:
!     status = NF90_EndDef( ncid )
!     if (status/=NF90_NOERR) stop
!
!     ! write in file:
!     call dim%Write( ncid, status )
!     if (status/=0) stop
!
!     ! done:
!     call dim%Done( status )
!     if (status/=0) stop
!
!     ! ~ read existing dimension:
!
!     ! init:
!     call dim%Init( 'value', status )
!     if (status/=0) stop
!
!     ! open existing netcdf file, obtain file id:
!     ncid = ...
!
!     ! read definition from file for the initialized name:
!     call dim%Read( file%ncid, status )
!     if (status/=0) stop
!
!     ! extract value:
!     call dim%Get_Dim( status, n=n )
!     if (status/=0) stop
!
!     ! done:
!     call dim%Done( status )
!     if (status/=0) stop
!
!
! COORDINATES
!
!   A coordinate is a dimension (a name and associated length)
!   with a corresponding set of data values.
!
!   Following the CF conventions, each coordinate variable
!   should have at least two attributes:
!     units                        : character string with the units
!     standard_name or long_name   : describes the content
!
!   Depending on the type of data, different coordinates are available.
!
!   * RealCoordinate
!       This type has a 1D array of real values.
!
!   * IntegerCoordinate
!       This type has a 1D array of integer values.
!
!   * HybrideLevelCoordinate
!       This type has 1D arrays with real values 'a' (or 'ap') and 'b',
!       and a scalar reference pressure 'p0' ; in addition it stores
!       the name of a surface pressure variable.
!       Eventually this type also has interface values 'ai' (or 'api') and 'bi'.
!
!   * TimeCoordinate
!       This type has a 1D array with time values.
!       Units are in the form "seconds since 2000-01-01 00:00:00" with appropriate
!       steps ("seconds", "hours", "days", etc) and time offset.
!
!   * LabelCoordinate
!       List of character labels.
! 
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; gol=NF90_StrError(status); call goErr; TRACEBACK; action; return; end if
!
#include "c3po.inc"
!
!######################################################################

module C3PO_Coordinates

  use GO    , only : gol, goPr, goErr
  use GO    , only : TDate
  use NetCDF, only : NF90_NOERR, NF90_StrError

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  Dimension
  public  ::  Coordinate
  public  ::  RealCoordinate
  public  ::  IntegerCoordinate
  public  ::  LabelCoordinate
  public  ::  HybrideLevelCoordinate
  public  ::  TimeCoordinate
  

  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_RealCoordinates'

  ! lengths:
  integer, parameter    ::  LEN_NAME  = 64
  integer, parameter    ::  LEN_UNITS = 64
  integer, parameter    ::  LEN_LINE  = 1024

  ! number of boundary values:
  integer, parameter           ::  nv = 2
  ! dimension name:
  character(len=*), parameter  ::  dname_nv = 'nv'

    
  ! --- types ----------------------------------------
  
  type Dimension
    ! default attributes:
    character(len=LEN_NAME)   ::  name
    integer                   ::  n
    logical                   ::  unlimited
    ! variable access:
    integer                   ::  dimid
  contains
    procedure   ::  Init    => Dimension_Init
    procedure   ::  Done    => Dimension_Done
    procedure   ::  Set_Dim => Dimension_Set_Dim
    procedure   ::  Get_Dim => Dimension_Get_Dim
    procedure   ::  Def     => Dimension_Def
    procedure   ::  Write   => Dimension_Write
    procedure   ::  Read    => Dimension_Read
  end type Dimension
  
  ! *
  
  type, extends(Dimension) ::  Coordinate
    ! default attributes:
    character(len=LEN_UNITS)    ::  units
    character(len=LEN_LINE)     ::  long_name
    character(len=LEN_LINE)     ::  standard_name
    character(len=LEN_LINE)     ::  formula
    ! variable access:
    integer                     ::  varid_in   ! defined by 'Coordinate_Read', 
                                               ! eventually used by other 'Read'
    integer                     ::  varid_out  ! defined by 'Coordinate_Def',
                                               ! used by '*Coordinate_Def' and 
                                               !   '*Coordinate_Write'
    ! boundaries:
    character(len=LEN_NAME)     ::  bounds
    integer                     ::  dimid_nv
    integer                     ::  varid_in_bnds
    integer                     ::  varid_out_bnds
  contains
    procedure   ::  Init    =>   Coordinate_Init
    procedure   ::  Done    =>   Coordinate_Done
    procedure   ::  Set_Attrs => Coordinate_Set_Attrs
    procedure   ::               Coordinate_Get_Attrs
    procedure   ::               Coordinate_Def
    procedure   ::               Coordinate_Write
    procedure   ::               Coordinate_Read
  end type Coordinate
  
  ! *
  
  type, extends(Coordinate) :: RealCoordinate
    ! data values:
    real, allocatable       ::  values(:)
    ! boundary values:
    real, allocatable       ::  value_bnds(:,:)
  contains
    procedure   ::  Done       => RealCoordinate_Done
    procedure   ::  Set_Values => RealCoordinate_Set_Values
    procedure   ::  Get_Values => RealCoordinate_Get_Values
    procedure   ::  Get_Index  => RealCoordinate_Get_Index
    procedure   ::  Def        => RealCoordinate_Def
    procedure   ::  Write      => RealCoordinate_Write
    procedure   ::  Read       => RealCoordinate_Read
  end type RealCoordinate
  
  ! *
  
  type, extends(Coordinate) :: IntegerCoordinate
    ! data values:
    integer, allocatable       ::  values(:)
  contains
    procedure   ::  Done       => IntegerCoordinate_Done
    procedure   ::  Set_Values => IntegerCoordinate_Set_Values
    procedure   ::  Get_Values => IntegerCoordinate_Get_Values
    procedure   ::  Def        => IntegerCoordinate_Def
    procedure   ::  Write      => IntegerCoordinate_Write
    procedure   ::  Read       => IntegerCoordinate_Read
  end type IntegerCoordinate
  
  ! *
  
  type, extends(Coordinate) :: LabelCoordinate
    ! data values:
    character(len=LEN_NAME), allocatable       ::  values(:)
    integer                                    ::  maxlen
  contains
    procedure   ::  Done       => LabelCoordinate_Done
    procedure   ::  Set_Value  => LabelCoordinate_Set_Value
    procedure   ::  Get_Value  => LabelCoordinate_Get_Value
    procedure   ::  Def        => LabelCoordinate_Def
    procedure   ::  Write      => LabelCoordinate_Write
    procedure   ::  Read       => LabelCoordinate_Read
  end type LabelCoordinate
  
  ! *
  
  type, extends(Coordinate) :: HybrideLevelCoordinate
    ! data values:
    logical                   ::  with_ap
    !logical                   ::  with_bnds
    logical                   ::  with_levi
    real, allocatable         ::  lev(:), a(:), ap(:), b(:)
    real, allocatable         ::  levi(:), ai(:), api(:), bi(:)
    !real, allocatable         ::  lev_bnds(:,:), a_bnds(:,:), ap_bnds(:,:), b_bnds(:,:)
    real                      ::  p0
    character(len=LEN_NAME)   ::  units_p
    character(len=LEN_NAME)   ::  positive
    character(len=LEN_NAME)   ::  vname_a
    character(len=LEN_NAME)   ::  vname_ap
    character(len=LEN_NAME)   ::  vname_b
    character(len=LEN_NAME)   ::  vname_p0
    character(len=LEN_NAME)   ::  vname_ps
    !character(len=LEN_LINE)   ::  long_name_bnds
    character(len=LEN_LINE)   ::  long_name_i
    character(len=LEN_LINE)   ::  formula_terms
    !character(len=LEN_LINE)   ::  formula_terms_bnds
    character(len=LEN_LINE)   ::  formula_terms_i
    ! variable access:
    integer                   ::  varid_out_a
    integer                   ::  varid_out_ap
    integer                   ::  varid_out_b
    integer                   ::  varid_out_i
    integer                   ::  varid_out_ai
    integer                   ::  varid_out_api
    integer                   ::  varid_out_bi
    !integer                   ::  varid_out_bnds
    !integer                   ::  varid_out_a_bnds
    !integer                   ::  varid_out_ap_bnds
    !integer                   ::  varid_out_b_bnds
    integer                   ::  varid_out_p0
  contains
    procedure   ::  Init       => HybrideLevelCoordinate_Init
    procedure   ::  Done       => HybrideLevelCoordinate_Done
    procedure   ::  Set_Values => HybrideLevelCoordinate_Set_Values
    procedure   ::  Get_Values => HybrideLevelCoordinate_Get_Values
    procedure   ::  Def        => HybrideLevelCoordinate_Def
    procedure   ::  Write      => HybrideLevelCoordinate_Write
    procedure   ::  Read       => HybrideLevelCoordinate_Read
  end type HybrideLevelCoordinate
  
  ! *
  
  type, extends(Coordinate) :: TimeCoordinate
    ! data values:
    real, allocatable         ::  values(:)
    character(len=LEN_NAME)   ::  units_step
    type(TDate)               ::  units_t0
    character(len=LEN_NAME)   ::  calendar
    logical                   ::  climatology
    real, allocatable         ::  climat(:,:)
    ! variable access:
    character(len=LEN_LINE)   ::  vname_climat
    integer                   ::  varid_out_climat
  contains
    procedure   ::  Init        => TimeCoordinate_Init
    procedure   ::  Done        => TimeCoordinate_Done
    procedure   ::  Set_Values  => TimeCoordinate_Set_Values
    procedure   ::  Set_Value   => TimeCoordinate_Set_Value
    procedure   ::  Get_Value   => TimeCoordinate_Get_Value
    procedure   ::  Def         => TimeCoordinate_Def
    procedure   ::  Write       => TimeCoordinate_Write
    procedure   ::  Write_Value => TimeCoordinate_Write_Value
    procedure   ::  Read        => TimeCoordinate_Read
  end type TimeCoordinate

! adhoc ...
#ifdef without_f2003
#define XTYPE type
#else
#define XTYPE class
#endif
  

contains


  ! ********************************************************************
  ! ***
  ! *** dimension
  ! ***
  ! ********************************************************************


  subroutine Dimension_Init( self, name, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(Dimension), intent(out)     ::  self
    character(len=*), intent(in)      ::  name
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Dimension_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%name = trim(name)
    ! zero size:
    self%n = 0
    ! undefined:
    self%unlimited = .false.
    self%dimid      = -999
    
    ! ok
    status = 0
    
  end subroutine Dimension_Init


  ! ***
  

  subroutine Dimension_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(Dimension), intent(inout)   ::  self
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Dimension_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! undefined:
    self%name       = ''
    self%n          = -999
    self%unlimited = .false.
    self%dimid      = -999
    
    ! ok
    status = 0
    
  end subroutine Dimension_Done


  ! ***
  

  subroutine Dimension_Set_Dim( self, status, &
                                  n, unlimited )
  
    ! --- in/out ---------------------------------
    
    XTYPE(Dimension), intent(inout)           ::  self
    integer, intent(out)                      ::  status

    integer, intent(in), optional             ::  n
    logical, intent(in), optional             ::  unlimited

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Dimension_Set_Dim'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! length:
    if ( present(n) ) then
      ! check ...
      if ( n < 0 ) then
        write (gol,'("RealCoordinate length should be zero or positive, not ",i0)') n
        TRACEBACK; status=1; return
      end if
      ! store:
      self%n = n
    end if
    
    ! flag:
    if ( present(unlimited) ) self%unlimited = unlimited
    
    ! ok
    status = 0
    
  end subroutine Dimension_Set_Dim


  ! ***
  

  subroutine Dimension_Get_Dim( self, status, &
                              name, n, unlimited, dimid )
  
    ! --- in/out ---------------------------------
    
    XTYPE(Dimension), intent(in)              ::  self
    integer, intent(out)                      ::  status
    
    character(len=*), intent(out), optional   ::  name
    integer, intent(out), optional            ::  n
    logical, intent(out), optional            ::  unlimited
    integer, intent(out), optional            ::  dimid

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Dimension_Get_Dim'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! name:
    if ( present(name) ) name = trim(self%name)
    
    ! flag:
    if ( present(unlimited) ) unlimited = self%unlimited
    
    ! length:
    if ( present(n) ) then
      ! check ...
      if ( self%n < 0 ) then
        write (gol,'("dimension length for `",a,"` not defined yet")') trim(self%name); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      n = self%n
    end if
    
    ! length:
    if ( present(dimid) ) then
      ! check ...
      if ( self%dimid <= 0 ) then
        write (gol,'("dimension id for `",a,"` not defined yet")') trim(self%name); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      dimid = self%dimid
    end if
    
    ! ok
    status = 0
    
  end subroutine Dimension_Get_Dim


  ! ***
  

  subroutine Dimension_Def( self, file, status )
  
    use NetCDF, only : NF90_Def_Dim, NF90_UNLIMITED
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(Dimension), intent(inout)     ::  self
    class(Datafile), intent(in)         ::  file
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Dimension_Def'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! switch:
    if ( self%unlimited ) then
      ! create dimension in file:
      status = NF90_Def_Dim( file%ncid, trim(self%name), NF90_UNLIMITED, self%dimid )
      IF_NF90_NOT_OK_RETURN(status=1)
    else
      ! create dimension in file:
      status = NF90_Def_Dim( file%ncid, trim(self%name), self%n, self%dimid )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine Dimension_Def


  ! ***
  

  subroutine Dimension_Write( self, file, status )
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(Dimension), intent(in)      ::  self
    class(Datafile), intent(in)       ::  file
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Dimension_Write'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! dimension already written by Dimension_Def,
    ! so nothing to be done ...
    
    ! ok
    status = 0
    
  end subroutine Dimension_Write


  ! ***
  

  subroutine Dimension_Read( self, file, status )
  
    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(Dimension), intent(inout)     ::  self
    class(Datafile), intent(in)         ::  file
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Dimension_Read'
    
    ! --- local ----------------------------------
    
    integer     ::  dimid
    
    ! --- begin ----------------------------------
    
    ! obtain dimension id for requested name:
    status = NF90_Inq_DimID( file%ncid, trim(self%name), dimid )
    if ( status /= NF90_NOERR) then
      gol=NF90_StrError(status); call goErr
      write (gol,'("could not inquire id for dimension `",a,"` in file `",a,"`")') &
              trim(self%name), trim(file%filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! obtain info:
    status = NF90_Inquire_Dimension( file%ncid, dimid, len=self%n )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Dimension_Read


  ! ********************************************************************
  ! ***
  ! *** Coordinate
  ! ***
  ! ********************************************************************


  subroutine Coordinate_Init( self, name, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(Coordinate), intent(out)    ::  self
    character(len=*), intent(in)      ::  name
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Coordinate_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init dimension with same name:
    call self%Dimension%Init( name, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! undefined:
    self%units          = ''
    self%long_name      = ''
    self%standard_name  = ''
    self%formula        = ''
    self%bounds         = ''
    self%varid_in       = -999
    self%varid_out      = -999
    self%varid_out_bnds = -999
    self%dimid_nv       = -999
    
    ! ok
    status = 0
    
  end subroutine Coordinate_Init


  ! ***
  

  subroutine Coordinate_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(Coordinate), intent(inout)      ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Coordinate_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! clear:
    call self%Dimension%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! undefined:
    self%units          = ''
    self%long_name      = ''
    self%standard_name  = ''
    self%formula        = ''
    self%bounds         = ''
    self%varid_in       = -999
    self%varid_out      = -999
    self%varid_out_bnds = -999
    self%dimid_nv       = -999
    
    ! ok
    status = 0
    
  end subroutine Coordinate_Done


  ! ***
  

  subroutine Coordinate_Set_Attrs( self, status, &
                                    units, long_name, standard_name, &
                                    bounds, formula )
  
    ! --- in/out ---------------------------------
    
    XTYPE(Coordinate), intent(inout)        ::  self
    integer, intent(out)                    ::  status
    
    character(len=*), intent(in), optional  ::  units
    character(len=*), intent(in), optional  ::  long_name
    character(len=*), intent(in), optional  ::  standard_name
    character(len=*), intent(in), optional  ::  bounds
    character(len=*), intent(in), optional  ::  formula

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Coordinate_Set_Attrs'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    if ( present(units        ) ) self%units         = trim(units)
    if ( present(long_name    ) ) self%long_name     = trim(long_name)
    if ( present(standard_name) ) self%standard_name = trim(standard_name)
    if ( present(formula      ) ) self%formula       = trim(formula)
    if ( present(bounds       ) ) self%bounds        = trim(bounds)
    
    ! ok
    status = 0
    
  end subroutine Coordinate_Set_Attrs


  ! ***
  

  subroutine Coordinate_Get_Attrs( self, status, &
                                    units, long_name, standard_name, &
                                    bounds, formula )
  
    ! --- in/out ---------------------------------
    
    XTYPE(Coordinate), intent(in)         ::  self
    integer, intent(out)                      ::  status
    
    character(len=*), intent(out), optional   ::  units
    character(len=*), intent(out), optional   ::  long_name
    character(len=*), intent(out), optional   ::  standard_name
    character(len=*), intent(out), optional   ::  formula
    character(len=*), intent(out), optional   ::  bounds

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Coordinate_Get_Attrs'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! attribute:
    if ( present(units) ) then
      ! check ...
      if ( len_trim(self%units) == 0 ) then
        write (gol,'("coordinate units for `",a,"` not defined yet")') trim(self%name); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      units = trim(self%units)
    end if
    
    ! attribute:
    if ( present(long_name) ) then
      ! copy:
      long_name = trim(self%long_name)
    end if
    
    ! attribute:
    if ( present(standard_name) ) then
      ! copy:
      standard_name = trim(self%standard_name)
    end if
    
    ! attribute:
    if ( present(formula) ) then
      ! copy:
      formula = trim(self%formula)
    end if
    
    ! attribute:
    if ( present(bounds) ) then
      ! copy:
      bounds = trim(self%bounds)
    end if
    
    ! ok
    status = 0
    
  end subroutine Coordinate_Get_Attrs


  ! ***
  

  subroutine Coordinate_Def( self, file, varid, status )
  
    use NetCDF, only : NF90_INQ_DimID, NF90_Def_Dim
    use NetCDF, only : NF90_Put_Att
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(Coordinate), intent(inout)  ::  self
    class(Datafile), intent(in)       ::  file
    integer, intent(in)               ::  varid
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Coordinate_Def'
    
    ! --- local ----------------------------------

    ! --- begin ----------------------------------
    
    ! units defined ?
    if ( len_trim(self%units) > 0 ) then
      ! define:
      status = NF90_Put_Att( file%ncid, varid, 'units', trim(self%units) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! check ...
    if ( (len_trim(self%long_name) == 0) .and. (len_trim(self%standard_name) == 0) ) then
      write (gol,'("at least one of `long_name` or `standard_name` should be defined")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! write:
    if ( len_trim(self%long_name) > 0 ) then
      status = NF90_Put_Att( file%ncid, varid, 'long_name', trim(self%long_name) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    ! write:
    if ( len_trim(self%standard_name) > 0 ) then
      status = NF90_Put_Att( file%ncid, varid, 'standard_name', trim(self%standard_name) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! write:
    if ( len_trim(self%formula) > 0 ) then
      status = NF90_Put_Att( file%ncid, varid, 'formula', trim(self%formula) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! boundaries defined ?
    if ( len_trim(self%bounds) > 0 ) then
      ! inquire dimension, error if not present:
      status = NF90_INQ_DimID( file%ncid, dname_nv, self%dimid_nv )
      if (status/=NF90_NOERR) then
        ! create new:
        status = NF90_Def_Dim( file%ncid, dname_nv, nv, self%dimid_nv )
        IF_NF90_NOT_OK_RETURN(status=1)
      end if
    end if  ! bounds

    ! ok
    status = 0
    
  end subroutine Coordinate_Def


  ! ***
  

  subroutine Coordinate_Write( self, file, status )

    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(Coordinate), intent(in)     ::  self
    class(Datafile), intent(in)       ::  file
    integer, intent(out)              ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Coordinate_Write'
    
    ! --- local ----------------------------------

    integer     ::  dimid
    
    ! --- begin ----------------------------------

    ! write dimension:
    call self%Dimension%Write( file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( self%varid_out < 0 ) then
      write (gol,'("no varid defined for coordinate values yet")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! boundaries defined?
    if ( len_trim(self%bounds) > 0 ) then
      ! check ...
      if ( self%varid_out_bnds < 0 ) then
        write (gol,'("no varid defined for coordinate bounds yet")'); call goErr
        TRACEBACK; status=1; return
      end if
    end if ! bounds

    ! ok
    status = 0
    
  end subroutine Coordinate_Write


  ! ***
  

  subroutine Coordinate_Read( self, file, status, varname )

    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable
    use NetCDF, only : NF90_Inq_AttName, NF90_Get_Att
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(Coordinate), intent(inout)        ::  self
    class(Datafile), intent(in)             ::  file
    integer, intent(out)                    ::  status
    
    ! alternative variable name,needed for labelcoordinate:
    character(len=*), intent(in), optional  ::  varname

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Coordinate_Read'
    
    ! --- local ----------------------------------

    character(len=LEN_NAME)   ::  name
    integer                   ::  natt, iatt
    character(len=64)         ::  attname
    character(len=LEN_NAME)   ::  vname
    
    ! --- begin ----------------------------------
    
    ! read dimension parameters (length):
    call self%Dimension%Read( file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! obtain name:
    call self%Dimension%Get_Dim( status, name=name )
    IF_NOT_OK_RETURN(status=1)

    ! variable name, default or explicit:
    vname = trim(self%name)
    if ( present(varname) ) vname = trim(varname)
    ! get variable id:
    status = NF90_INQ_VarID( file%ncid, trim(vname), self%varid_in )
    if ( status /= NF90_NOERR) then
      gol=NF90_StrError(status); call goErr
      write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
              trim(vname), trim(file%filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! count number of attributes:
    status = NF90_Inquire_Variable( file%ncid, self%varid_in, nAtts=natt )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! loop over attributes:
    do iatt = 1, natt
      ! name:
      status = NF90_Inq_AttName( file%ncid, self%varid_in, iatt, attname )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! switch:
      select case ( trim(attname) )
        !
        case ( 'units' )
          ! read:
          status = NF90_Get_Att( file%ncid, self%varid_in, trim(attname), self%units )
          IF_NF90_NOT_OK_RETURN(status=1)
        !
        case ( 'long_name' )
          ! read:
          status = NF90_Get_Att( file%ncid, self%varid_in, trim(attname), self%long_name )
          IF_NF90_NOT_OK_RETURN(status=1)
        !
        case ( 'standard_name' )
          ! read:
          status = NF90_Get_Att( file%ncid, self%varid_in, trim(attname), self%standard_name )
          IF_NF90_NOT_OK_RETURN(status=1)
        !
        case ( 'formula' )
          ! read:
          status = NF90_Get_Att( file%ncid, self%varid_in, trim(attname), self%formula )
          IF_NF90_NOT_OK_RETURN(status=1)
        !
        case ( 'bounds' )
          ! read:
          status = NF90_Get_Att( file%ncid, self%varid_in, trim(attname), self%bounds )
          IF_NF90_NOT_OK_RETURN(status=1)
        !
        case default
          write (gol,'("unsupported attribute `",a,"`")') trim(attname); call goErr
          TRACEBACK; status=1; return
      end select
    end do ! attributes
    
    ! boundaries defined ?
    if ( len_trim(self%bounds) > 0 ) then
      ! get variable id:
      status = NF90_INQ_VarID( file%ncid, trim(self%bounds), self%varid_in_bnds )
      if ( status /= NF90_NOERR) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
                trim(self%bounds), trim(file%filename); call goErr
        TRACEBACK; status=1; return
      end if
    end if ! bounds
    
    ! ok
    status = 0
    
  end subroutine Coordinate_Read


  ! ********************************************************************
  ! ***
  ! *** real valued coordinate
  ! ***
  ! ********************************************************************


  subroutine RealCoordinate_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(RealCoordinate), intent(inout)  ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/RealCoordinate_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! general coordinate:
    call self%Coordinate%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    if ( allocated(self%values) ) then
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine RealCoordinate_Done


  ! ***
  

  subroutine RealCoordinate_Set_Values( self, status, values, value_bnds )
  
    ! --- in/out ---------------------------------
    
    XTYPE(RealCoordinate), intent(inout)    ::  self
    integer, intent(out)                    ::  status
    real, intent(in), optional              ::  values(:)
    real, intent(in), optional              ::  value_bnds(:,:)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/RealCoordinate_Set_Values'
    
    ! --- local ----------------------------------
    
    integer   ::  nvalue
    
    ! --- begin ----------------------------------
    
    ! get size:
    call self%Coordinate%Get_Dim( status, n=nvalue )
    IF_NOT_OK_RETURN(status=1)
    
    ! boundary values provided ?
    if ( present(value_bnds) ) then
      ! set variable name if not present yet:
      if ( len_trim(self%bounds) == 0 ) then
        ! extend name:
        self%bounds = trim(self%name)//'_bnds'
      end if
      ! check ...
      if ( any( shape(value_bnds) /= (/nv,nvalue/) ) ) then
        write (gol,'("values array has shape (",i0,",",i0,") while expected (",i0,",",i0,")")') &
                shape(value_bnds), nv, nvalue; call goErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      if ( .not. allocated(self%value_bnds) ) then
        allocate( self%value_bnds(nv,nvalue), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! copy:
      self%value_bnds = value_bnds
    end if
    
    ! values provided ?
    if ( present(values) ) then
      ! check ...
      if ( size(values) /= nvalue ) then
        write (gol,'("values array has size ",i0," while dimension is ",i0)') size(values), nvalue; call goErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      if ( .not. allocated(self%values) ) then
        allocate( self%values(nvalue), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! copy:
      self%values = values
    end if
    
    ! no values yet, but boundary values provided ? then set to average:
    if ( (.not. allocated(self%values)) .and. present(value_bnds) ) then
      ! storage:
      if ( .not. allocated(self%values) ) then
        allocate( self%values(nvalue), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! average:
      self%values = sum(self%value_bnds,1) / real(nv)
    end if
    
    ! ok
    status = 0
    
  end subroutine RealCoordinate_Set_Values


  ! ***
  

  subroutine RealCoordinate_Get_Values( self, status, values, value_bnds )
  
    ! --- in/out ---------------------------------
    
    XTYPE(RealCoordinate), intent(in)         ::  self
    integer, intent(out)                      ::  status
    real, intent(out), optional               ::  values(:)
    real, intent(out), optional               ::  value_bnds(:,:)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/RealCoordinate_Get_Values'
    
    ! --- local ----------------------------------
    
    integer     ::  nvalue
    
    ! --- begin ----------------------------------
    
    ! get size:
    call self%Coordinate%Get_Dim( status, n=nvalue )
    IF_NOT_OK_RETURN(status=1)
    
    ! return values?
    if ( present(values) ) then
      ! check ...
      if ( .not. allocated(self%values) ) then
        write (gol,'("no values stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( size(values) /= nvalue ) then
        write (gol,'("values array has size ",i0," while dimension is ",i0)') size(values), nvalue; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      values = self%values
    end if
    
    ! return values?
    if ( present(value_bnds) ) then
      ! check ...
      if ( .not. allocated(self%value_bnds) ) then
        write (gol,'("no value_bnds stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( any( shape(value_bnds) /= (/2,nvalue/) ) ) then
        write (gol,'("values array has shape (",i0,",",i0,") while dimension is (",i0,",",i0,")")') &
                    shape(value_bnds), 2,nvalue; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      value_bnds = self%value_bnds
    end if

    ! ok
    status = 0
    
  end subroutine RealCoordinate_Get_Values


  ! ***
  

  ! return index with smallest distance to value
  
  subroutine RealCoordinate_Get_Index( self, value, ind, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(RealCoordinate), intent(in)         ::  self
    real, intent(in)                          ::  value
    integer, intent(out)                      ::  ind
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/RealCoordinate_Get_Index'
    
    ! --- local ----------------------------------
    
    integer     ::  i
    real        ::  d, dmin
    
    ! --- begin ----------------------------------
    
    ! loop over values:
    do i = 1, self%n
      ! distance:
      d = abs( value - self%values(i) )
      ! init or update ?
      if ( i == 1 ) then
        ind = i
        dmin = d
      else if ( d < dmin ) then
        ind = i
        dmin = d
      end if
    end do  ! values
    
    ! ok
    status = 0
    
  end subroutine RealCoordinate_Get_Index


  ! ***


  subroutine RealCoordinate_Def( self, file, status )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(RealCoordinate), intent(inout)    ::  self
    class(Datafile), intent(in)             ::  file
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/RealCoordinate_Def'
    
    ! --- local ----------------------------------

    character(len=LEN_NAME)   ::  name
    integer                   ::  dimid
    
    ! --- begin ----------------------------------
    
    ! define dimension:
    call self%Dimension%Def( file, status )
    IF_NOT_OK_RETURN(status=1)
    ! get access id:
    call self%Dimension%Get_Dim( status, name=name, dimid=dimid )
    IF_NOT_OK_RETURN(status=1)

    ! define new variable:
    status = NF90_Def_Var( file%ncid, trim(name), NF90_FLOAT, (/dimid/), self%varid_out )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! add default attributes:
    call self%Coordinate_Def( file, self%varid_out, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! boundaries defined?
    if ( len_trim(self%bounds) > 0 ) then
      ! extra attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out, 'bounds', trim(self%bounds) )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! define new variable:
      status = NF90_Def_Var( file%ncid, trim(self%bounds), NF90_FLOAT, &
                                 (/self%dimid_nv,dimid/), self%varid_out_bnds )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! add default attributes:
      call self%Coordinate_Def( file, self%varid_out_bnds, status )
      IF_NOT_OK_RETURN(status=1)
    end if  ! bounds
    
    ! ok
    status = 0
    
  end subroutine RealCoordinate_Def


  ! ***
  

  subroutine RealCoordinate_Write( self, file, status )

    use NetCDF, only : NF90_Put_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(RealCoordinate), intent(in)   ::  self
    class(Datafile), intent(in)         ::  file
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/RealCoordinate_Write'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! write general coordinate, perform checks:
    call self%Coordinate_Write( file, status )
    IF_NOT_OK_RETURN(status=1)

    ! check ...
    if ( .not. allocated(self%values) ) then
      write (gol,'("no values stored yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! write:
    status = NF90_Put_Var( file%ncid, self%varid_out, self%values )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! boundaries defined?
    if ( len_trim(self%bounds) > 0 ) then
      ! check ...
      if ( .not. allocated(self%value_bnds) ) then
        write (gol,'("no value_bnds stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! write:
      status = NF90_Put_Var( file%ncid, self%varid_out_bnds, self%value_bnds )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if ! bounds
    
    ! ok
    status = 0
    
  end subroutine RealCoordinate_Write


  ! ***
  

  subroutine RealCoordinate_Read( self, file, status )

    use NetCDF, only : NF90_Get_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(RealCoordinate), intent(inout)    ::  self
    class(Datafile), intent(in)             ::  file
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/RealCoordinate_Read'
    
    ! --- local ----------------------------------

    integer                   ::  nvalue
    
    ! --- begin ----------------------------------

    ! read general coordinate parameters:
    !   length, units, varid_in[_bnds]
    call self%Coordinate_Read( file, status )
    IF_NOT_OK_RETURN(status=1)

    ! get size:
    call self%Coordinate%Get_Dim( status, n=nvalue )
    IF_NOT_OK_RETURN(status=1)
    
    ! storage:
    if ( .not. allocated(self%values) ) then
      allocate( self%values(nvalue), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! read:
    status = NF90_Get_Var( file%ncid, self%varid_in, self%values )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! boundaries ?
    if ( len_trim(self%bounds) > 0 ) then
      ! storage:
      if ( .not. allocated(self%value_bnds) ) then
        allocate( self%value_bnds(nv,nvalue), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! read:
      status = NF90_Get_Var( file%ncid, self%varid_in_bnds, self%value_bnds )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! ok
    status = 0
    
  end subroutine RealCoordinate_Read


  ! ********************************************************************
  ! ***
  ! *** integer valued coordinate
  ! ***
  ! ********************************************************************


  subroutine IntegerCoordinate_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(IntegerCoordinate), intent(inout)  ::  self
    integer, intent(out)                  ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntegerCoordinate_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! general coordinate:
    call self%Coordinate%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    if ( allocated(self%values) ) then
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine IntegerCoordinate_Done


  ! ***
  

  subroutine IntegerCoordinate_Set_Values( self, status, values )
  
    ! --- in/out ---------------------------------
    
    XTYPE(IntegerCoordinate), intent(inout)   ::  self
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  values(:)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntegerCoordinate_Set_Values'
    
    ! --- local ----------------------------------
    
    integer   ::  nvalue
    
    ! --- begin ----------------------------------
    
    ! get size:
    call self%Coordinate%Get_Dim( status, n=nvalue )
    IF_NOT_OK_RETURN(status=1)
    
    ! values provided ?
    if ( present(values) ) then
      ! check ...
      if ( size(values) /= nvalue ) then
        write (gol,'("values array has size ",i0," while dimension is ",i0)') size(values), nvalue; call goErr
        TRACEBACK; status=1; return
      end if
      ! storage:
      if ( .not. allocated(self%values) ) then
        allocate( self%values(nvalue), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! copy:
      self%values = values
    end if
    
    ! ok
    status = 0
    
  end subroutine IntegerCoordinate_Set_Values


  ! ***
  

  subroutine IntegerCoordinate_Get_Values( self, values, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(IntegerCoordinate), intent(in)      ::  self
    integer, intent(out)                      ::  values(:)
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntegerCoordinate_Get_Values'
    
    ! --- local ----------------------------------
    
    integer     ::  nvalue
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( .not. allocated(self%values) ) then
      write (gol,'("no values stored yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! get size:
    call self%Coordinate%Get_Dim( status, n=nvalue )
    IF_NOT_OK_RETURN(status=1)
    ! check ...
    if ( size(values) /= nvalue ) then
      write (gol,'("values array has size ",i0," while dimension is ",i0)') size(values), nvalue; call goErr
      TRACEBACK; status=1; return
    end if
    ! copy:
    values = self%values
    
    ! ok
    status = 0
    
  end subroutine IntegerCoordinate_Get_Values


  ! ***


  subroutine IntegerCoordinate_Def( self, file, status )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(IntegerCoordinate), intent(inout)   ::  self
    class(Datafile), intent(in)               ::  file
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntegerCoordinate_Def'
    
    ! --- local ----------------------------------

    character(len=LEN_NAME)   ::  name
    integer                   ::  dimid
    
    ! --- begin ----------------------------------
    
    ! define dimension:
    call self%Dimension%Def( file, status )
    IF_NOT_OK_RETURN(status=1)
    ! get access id:
    call self%Dimension%Get_Dim( status, name=name, dimid=dimid )
    IF_NOT_OK_RETURN(status=1)

    ! define new variable:
    status = NF90_Def_Var( file%ncid, trim(name), NF90_FLOAT, (/dimid/), self%varid_out )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! define default attributes:
    call self%Coordinate_Def( file, self%varid_out, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine IntegerCoordinate_Def


  ! ***
  

  subroutine IntegerCoordinate_Write( self, file, status )

    use NetCDF, only : NF90_Put_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(IntegerCoordinate), intent(in)    ::  self
    class(Datafile), intent(in)             ::  file
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntegerCoordinate_Write'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! write general coordinate:
    call self%Coordinate_Write( file, status )
    IF_NOT_OK_RETURN(status=1)

    ! check ...
    if ( .not. allocated(self%values) ) then
      write (gol,'("no values stored yet")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! write:
    status = NF90_Put_Var( file%ncid, self%varid_out, self%values )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine IntegerCoordinate_Write


  ! ***
  

  subroutine IntegerCoordinate_Read( self, file, status )

    use NetCDF, only : NF90_Get_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(IntegerCoordinate), intent(inout)     ::  self
    class(Datafile), intent(in)                 ::  file
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/IntegerCoordinate_Read'
    
    ! --- local ----------------------------------

    integer     ::  nvalue
    
    ! --- begin ----------------------------------

    ! read general coordinate parameters:
    !   length, units, varid_in
    call self%Coordinate_Read( file, status )
    IF_NOT_OK_RETURN(status=1)

    ! get size:
    call self%Coordinate%Get_Dim( status, n=nvalue )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    if ( .not. allocated(self%values) ) then
      allocate( self%values(nvalue), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! read:
    status = NF90_Get_Var( file%ncid, self%varid_in, self%values )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine IntegerCoordinate_Read


  ! ********************************************************************
  ! ***
  ! *** character valued coordinate
  ! ***
  ! ********************************************************************


  subroutine LabelCoordinate_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(LabelCoordinate), intent(inout)  ::  self
    integer, intent(out)                   ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/LabelCoordinate_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! general coordinate:
    call self%Coordinate%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    if ( allocated(self%values) ) then
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
    
  end subroutine LabelCoordinate_Done


  ! ***
  

  subroutine LabelCoordinate_Set_Value( self, irec, status, &
                                          value )
  
    ! --- in/out ---------------------------------
    
    XTYPE(LabelCoordinate), intent(inout)     ::  self
    integer, intent(in)                       ::  irec
    integer, intent(out)                      ::  status
    
    character(len=*), intent(in), optional    ::  value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/LabelCoordinate_Set_Value'
    
    ! --- local ----------------------------------
    
    character(len=LEN_NAME), allocatable   ::  tmp_values(:)
    integer                                ::  k
    
    ! --- begin ----------------------------------
    
    ! * check record index:
    
    ! switch:
    if ( self%unlimited ) then
      ! check ...
      if ( irec < 1 ) then
        write (gol,'("argument `irec` out of range 1 .. ")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! reset:
      self%n = max( self%n, irec )
    else
      ! check ...
      if ( self%n <= 0 ) then
        write (gol,'("dimension length for `",a,"` not defined yet")') trim(self%name); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( (irec < 1) .or. (irec > self%n) ) then
        write (gol,'("argument `irec` out of range 1 .. ",i0)') self%n; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( self%n < 1 ) then
      write (gol,'("number of time records still undefined")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! * value
    
    ! provided as argument ?
    if ( present(value) ) then
      ! initialize storage if necessary:
      if ( .not. allocated(self%values) ) then
        ! new:
        allocate( self%values(self%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else if ( size(self%values) < self%n ) then
        ! storage for copy:
        allocate( tmp_values(size(self%values)), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! copy:
        do k = 1, size(self%values)
          tmp_values(k) = trim(self%values(k))
        end do
        ! clear existing:
        deallocate( self%values, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! new storage:
        allocate( self%values(self%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! copy:
        do k = 1, size(tmp_values)
          self%values(k) = trim(tmp_values(k))
        end do
        ! clear copy:
        deallocate( tmp_values, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! fill:
      self%values(irec) = trim(value)
    end if
    
    ! ok
    status = 0
    
  end subroutine LabelCoordinate_Set_Value


  ! ***
  

  subroutine LabelCoordinate_Get_Value( self, irec, status, value )
  
    ! --- in/out ---------------------------------
    
    XTYPE(LabelCoordinate), intent(in)        ::  self
    integer, intent(in)                       ::  irec
    integer, intent(out)                      ::  status
    
    character(len=*), intent(out), optional   ::  value

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/LabelCoordinate_Get_Value'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( self%n <= 0 ) then
      write (gol,'("dimension length for `",a,"` not defined yet")') trim(self%name); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( .not. allocated(self%values) ) then
      write (gol,'("no values stored yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( (irec < 1) .or. (irec > self%n) ) then
      write (gol,'("argument `irec` out of range 1 .. ",i0)') self%n; call goErr
      TRACEBACK; status=1; return
    end if
    ! copy:
    value = self%values(irec)
    
    ! ok
    status = 0
    
  end subroutine LabelCoordinate_Get_Value


  ! ***


  subroutine LabelCoordinate_Def( self, file, status )
  
    use NetCDF, only : NF90_Def_Dim
    use NetCDF, only : NF90_CHAR
    use NetCDF, only : NF90_Def_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(LabelCoordinate), intent(inout)     ::  self
    class(Datafile), intent(in)               ::  file
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/LabelCoordinate_Def'
    
    ! --- local ----------------------------------

    character(len=LEN_NAME)   ::  name
    integer                   ::  dimid_value
    integer                   ::  dimid_maxlen
    integer                   ::  k
    
    ! --- begin ----------------------------------
    
    ! define dimension:
    call self%Dimension%Def( file, status )
    IF_NOT_OK_RETURN(status=1)
    ! get access id:
    call self%Dimension%Get_Dim( status, name=name, dimid=dimid_value )
    IF_NOT_OK_RETURN(status=1)

    ! check ...
    if ( self%n <= 0 ) then
      write (gol,'("dimension length for `",a,"` not defined yet")') trim(self%name); call goErr
      TRACEBACK; status=1; return
    end if

    ! max length of tracer names:
    self%maxlen = -999
    do k = 1, self%n
      self%maxlen = max( self%maxlen, len_trim(self%values(k)) )
    end do
    ! define dimension:
    status = NF90_Def_Dim( file%ncid, trim(name)//'_len', self%maxlen, dimid_maxlen )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! define new variable:
    status = NF90_Def_Var( file%ncid, trim(name)//'_name', NF90_CHAR, &
                             (/dimid_maxlen,dimid_value/), self%varid_out )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! define default attributes:
    call self%Coordinate_Def( file, self%varid_out, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine LabelCoordinate_Def


  ! ***
  

  subroutine LabelCoordinate_Write( self, file, status )

    use NetCDF, only : NF90_Put_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(LabelCoordinate), intent(in)      ::  self
    class(Datafile), intent(in)             ::  file
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/LabelCoordinate_Write'
    
    ! --- local ----------------------------------
    
    character(len=LEN_NAME)   ::  value
    integer                   ::  k
    integer                   ::  vlen
    
    ! --- begin ----------------------------------

    ! write general coordinate:
    call self%Coordinate_Write( file, status )
    IF_NOT_OK_RETURN(status=1)

    ! check ...
    if ( .not. allocated(self%values) ) then
      write (gol,'("no values stored yet")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! loop over labels:
    do k = 1, self%n
      ! current lenth:
      vlen = len_trim(self%values(k))
      ! white space:
      value = repeat( ' ', self%maxlen )
      ! replace first part:
      value(1:vlen) = trim(self%values(k))
      ! write padded value:
      status = NF90_Put_Var( file%ncid, self%varid_out, value, &
                               start=(/1,k/), count=(/self%maxlen,1/) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end do
    
    ! ok
    status = 0
    
  end subroutine LabelCoordinate_Write


  ! ***
  

  subroutine LabelCoordinate_Read( self, file, status )

    use NetCDF, only : NF90_Inq_DimID, NF90_Inquire_Dimension
    use NetCDF, only : NF90_Get_Var

    use GO           , only : goCharToString
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(LabelCoordinate), intent(inout)       ::  self
    class(Datafile), intent(in)                 ::  file
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/LabelCoordinate_Read'
    
    ! --- local ----------------------------------

    integer                         ::  nvalue
    integer                         ::  dimid_maxlen
    character(len=1), allocatable   ::  chars(:)
    integer                         ::  k
    
    ! --- begin ----------------------------------

    ! read general coordinate parameters:
    !   length, units, varid_in
    ! use alternative name:
    call self%Coordinate_Read( file, status, varname=trim(self%name)//'_name' )
    IF_NOT_OK_RETURN(status=1)

    ! get size:
    call self%Coordinate%Get_Dim( status, n=nvalue )
    IF_NOT_OK_RETURN(status=1)
    ! storage:
    if ( .not. allocated(self%values) ) then
      allocate( self%values(nvalue), stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! get dimension for maximum length:
    status = NF90_Inq_DimID( file%ncid, trim(self%name)//'_len', dimid_maxlen )
    if ( status /= NF90_NOERR) then
      gol=NF90_StrError(status); call goErr
      write (gol,'("could not inquire id for dimension `",a,"` in file `",a,"`")') &
              trim(self%name)//'_len', trim(file%filename); call goErr
      TRACEBACK; status=1; return
    end if
    ! get size:
    status = NF90_Inquire_Dimension( file%ncid, dimid_maxlen, len=self%maxlen )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! storage:
    allocate( chars(self%maxlen), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! loop:
    do k = 1, self%n
      ! read padded label:
      status = NF90_Get_Var( file%ncid, self%varid_in, chars, &
                                start=(/1,k/), count=(/self%maxlen,1/) )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! copy, take care of null characters:
      call goCharToString( chars, self%values(k), status )
      IF_NOT_OK_RETURN(status=1)
    end do
    ! clear:
    deallocate( chars, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine LabelCoordinate_Read


  ! ********************************************************************
  ! ***
  ! *** hybride level coordinate
  ! ***
  ! ********************************************************************


  subroutine HybrideLevelCoordinate_Init( self, name, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(HybrideLevelCoordinate), intent(out)    ::  self
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/HybrideLevelCoordinate_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init general coordinate:
    call self%Coordinate%Init( name, status )
    
    ! undefined:
    self%with_ap           = .false.
    !self%with_bnds         = .false.
    self%with_levi         = .false.
    self%p0                = -999.9
    self%units_p           = ''
    self%vname_a           = 'a'
    self%vname_ap          = 'ap'
    self%vname_b           = 'b'
    self%vname_p0          = 'p0'
    self%vname_ps          = 'ps'
    self%varid_out_a       = -999
    self%varid_out_ap      = -999
    self%varid_out_b       = -999
    self%varid_out_i       = -999
    self%varid_out_ai      = -999
    self%varid_out_api     = -999
    self%varid_out_bi      = -999
    !self%varid_out_bnds    = -999
    !self%varid_out_a_bnds  = -999
    !self%varid_out_ap_bnds = -999
    !self%varid_out_b_bnds  = -999
    self%varid_out_p0      = -999
    
    ! ok
    status = 0
    
  end subroutine HybrideLevelCoordinate_Init


  ! ***


  subroutine HybrideLevelCoordinate_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    XTYPE(HybrideLevelCoordinate), intent(inout)  ::  self
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/HybrideLevelCoordinate_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! general coordinate:
    call self%Coordinate%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    if ( allocated(self%lev) ) then
      deallocate( self%lev, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%a) ) then
      deallocate( self%a, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%ap) ) then
      deallocate( self%ap, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%b) ) then
      deallocate( self%b, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%levi) ) then
      deallocate( self%levi, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%ai) ) then
      deallocate( self%ai, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%api) ) then
      deallocate( self%api, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%bi) ) then
      deallocate( self%bi, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    !if ( allocated(self%lev_bnds) ) then
    !  deallocate( self%lev_bnds, stat=status )
    !  IF_NOT_OK_RETURN(status=1)
    !end if
    !if ( allocated(self%a_bnds) ) then
    !  deallocate( self%a_bnds, stat=status )
    !  IF_NOT_OK_RETURN(status=1)
    !end if
    !if ( allocated(self%ap_bnds) ) then
    !  deallocate( self%ap_bnds, stat=status )
    !  IF_NOT_OK_RETURN(status=1)
    !end if
    !if ( allocated(self%b_bnds) ) then
    !  deallocate( self%b_bnds, stat=status )
    !  IF_NOT_OK_RETURN(status=1)
    !end if
    
    ! undefined:
    self%with_ap   = .false.
    self%with_levi = .false.
    !self%with_bnds = .false.
    self%p0        = -999.9
    self%units_p   = ''
    self%vname_a   = ''
    self%vname_ap  = ''
    self%vname_b   = ''
    self%vname_p0  = ''
    self%vname_ps  = ''
    self%varid_out_a       = -999
    self%varid_out_ap      = -999
    self%varid_out_b       = -999
    self%varid_out_i       = -999
    self%varid_out_ai      = -999
    self%varid_out_api     = -999
    self%varid_out_bi      = -999
    !self%varid_out_bnds    = -999
    !self%varid_out_a_bnds  = -999
    !self%varid_out_ap_bnds = -999
    !self%varid_out_b_bnds  = -999
    self%varid_out_p0      = -999
    
    ! ok
    status = 0
    
  end subroutine HybrideLevelCoordinate_Done


  ! ***
  
  
  !
  ! Define hybrid level values:
  !    p = ap     + b * ps  , lev = ap/p0 + b
  !    p = a * p0 + b * ps  , lev = a     + b
  !
  ! Required:
  !   p0       :   reference surface pressure
  !   units_p  :   units of p0 and ps
  !    
  ! Computed from arguments:
  !    a|ap, b  :  computed from ai|api and bi if only these are present
  !    lev      :  comptued from: a or ap/p0, and b
  !
  

  subroutine HybrideLevelCoordinate_Set_Values( self, status, &
                                                  a, ap, b, &
                                                  ai, api, bi, &
                                                  !a_bnds, ap_bnds, b_bnds, &
                                                  p0, units_p, vname_ps )
  
    ! --- in/out ---------------------------------
    
    XTYPE(HybrideLevelCoordinate), intent(inout)   ::  self
    integer, intent(out)                           ::  status

    real, intent(in), optional                     ::  a(:), ap(:), b(:)  ! (nlev)
    real, intent(in), optional                     ::  ai(:), api(:), bi(:)  ! (nlev+1)
    !real, intent(in), optional                     ::  a_bnds(:,:), ap_bnds(:,:), b_bnds(:,:) ! (nv,nlev)
    real, intent(in), optional                     ::  p0
    character(len=*), intent(in), optional         ::  units_p
    character(len=*), intent(in), optional         ::  vname_ps

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/HybrideLevelCoordinate_Set_Values'
    
    ! --- local ----------------------------------
    
    integer               ::  nlev, ilev
    character(len=1)      ::  mark
    
    ! --- begin ----------------------------------
    
    ! pre-defined:
    call self%Coordinate%Set_Attrs( status, units='1', &
                   standard_name='atmosphere_hybrid_sigma_pressure_coordinate' )
    IF_NOT_OK_RETURN(status=1)
    
    ! get size:
    call self%Coordinate%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)
    
    ! * check argument shapes
    
    ! check ...
    if ( present(a) ) then
      if ( any( shape(a) /= (/nlev/) ) ) then
        write (gol,'("argument `a` has shape (",i0,") while nlev is ",i0)') shape(a), nlev; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( present(ap) ) then
      if ( any( shape(ap) /= (/nlev/) ) ) then
        write (gol,'("argument `ap` has shape (",i0,") while nlev is ",i0)') shape(ap), nlev; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( present(b) ) then
      if ( any( shape(b) /= (/nlev/) ) ) then
        write (gol,'("argument `b` has shape (",i0,") while nlev is ",i0)') shape(b), nlev; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( present(ai) ) then
      if ( any( shape(ai) /= (/nlev+1/) ) ) then
        write (gol,'("argument `ai` has shape (",i0,") while nlev+1 is ",i0)') shape(ai), nlev+1; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( present(api) ) then
      if ( any( shape(api) /= (/nlev+1/) ) ) then
        write (gol,'("argument `api` has shape (",i0,") while nlev+1 is ",i0)') shape(api), nlev+1; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( present(bi) ) then
      if ( any( shape(bi) /= (/nlev+1/) ) ) then
        write (gol,'("argument `bi` has shape (",i0,") while nlev+1 is ",i0)') shape(bi), nlev+1; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    !! check ...
    !if ( present(a_bnds) ) then
    !  if ( any( shape(a_bnds) /= (/nv,nlev/) ) ) then
    !    write (gol,'("argument `a_bnds` has shape (",i0,",",i0,") while nlev is ",i0)') shape(a_bnds), nlev; call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !end if
    !! check ...
    !if ( present(ap_bnds) ) then
    !  if ( any( shape(ap_bnds) /= (/nv,nlev/) ) ) then
    !    write (gol,'("argument `ap_bnds` has shape (",i0,",",i0,") while nlev is ",i0)') shape(ap_bnds), nlev; call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !end if
    !! check ...
    !if ( present(b_bnds) ) then
    !  if ( any( shape(b_bnds) /= (/nv,nlev/) ) ) then
    !    write (gol,'("argument `b_bnds` has shape (",i0,",",i0,") while nlev is ",i0)') shape(b_bnds), nlev; call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !end if
    
    ! * check on ap or a*p0 mode
    
    ! ap or a*p0 ?
    self%with_ap = present(ap) .or. present(api) !.or. present(ap_bnds)
    
    ! check ...
    if ( self%with_ap ) then
      ! check ...
      if ( present(a) ) then
        !write (gol,'("provided `a` while one of `ap/api/ap_bnds` was provided")'); call goErr
        write (gol,'("provided `a` while one of `ap/api` was provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( present(ai) ) then
        !write (gol,'("provided `ai` while one of `ap/api/ap_bnds` was provided")'); call goErr
        write (gol,'("provided `ai` while one of `ap/api` was provided")'); call goErr
        TRACEBACK; status=1; return
      end if
      !! check ...
      !if ( present(a_bnds) ) then
      !  !write (gol,'("provided `a_bnds` while one of `ap/api/ap_bnds` was provided")'); call goErr
      !  write (gol,'("provided `a_bnds` while one of `ap/api` was provided")'); call goErr
      !  TRACEBACK; status=1; return
      !end if
    end if
    
    ! * check if boundaries are requested
    
    ! levi provided ?
    self%with_levi = present(ai) .or. present(api) .or. present(bi)
    
    !! bnds provided ?
    !self%with_bnds = present(ap_bnds) .or. present(a_bnds) .or. present(b_bnds) &
    !                   .or. present(ai) .or. present(api) .or. present(bi)
    
    ! * fill bounds if necessary
    
    !! fill bounds?
    !if ( self%with_bnds ) then
    !  ! switch ..
    !  if ( self%with_ap ) then
    !    ! storage:
    !    allocate( self%ap_bnds(nv,nlev), stat=status )
    !    IF_NOT_OK_RETURN(status=1)
    !    ! select:
    !    if ( present(api) ) then
    !      ! copy:
    !      self%ap_bnds(1,:) = api(1:nlev  )
    !      self%ap_bnds(2,:) = api(2:nlev+1)
    !    else if ( present(ap_bnds) ) then
    !      ! copy:
    !      self%ap_bnds = ap_bnds
    !    else
    !      write (gol,'("could not fill `ap_bnds`")'); call goErr
    !      TRACEBACK; status=1; return
    !    end if
    !  else
    !    ! storage:
    !    allocate( self%a_bnds(nv,nlev), stat=status )
    !    IF_NOT_OK_RETURN(status=1)
    !    ! select:
    !    if ( present(ai) ) then
    !      ! copy:
    !      self%a_bnds(1,:) = ai(1:nlev  )
    !      self%a_bnds(2,:) = ai(2:nlev+1)
    !    else if ( present(a_bnds) ) then
    !      ! copy:
    !      self%a_bnds = a_bnds
    !    else
    !      write (gol,'("could not fill `a_bnds`")'); call goErr
    !      TRACEBACK; status=1; return
    !    end if
    !  end if
    !  ! storage:
    !  allocate( self%b_bnds(nv,nlev), stat=status )
    !  IF_NOT_OK_RETURN(status=1)
    !  ! select:
    !  if ( present(bi) ) then
    !    ! copy:
    !    self%b_bnds(1,:) = bi(1:nlev  )
    !    self%b_bnds(2,:) = bi(2:nlev+1)
    !  else if ( present(b_bnds) ) then
    !    ! copy:
    !    self%b_bnds = b_bnds
    !  else
    !    write (gol,'("could not fill `b_bnds`")'); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !end if
    
    ! fill interfaces ?
    if ( self%with_levi ) then
      ! switch ..
      if ( self%with_ap ) then
        ! storage:
        allocate( self%api(nlev+1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! select:
        if ( present(api) ) then
          ! copy:
          self%api = api
        else
          write (gol,'("could not fill `api`")'); call goErr
          TRACEBACK; status=1; return
        end if
      else
        ! storage:
        allocate( self%ai(nlev+1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! select:
        if ( present(ai) ) then
          ! copy:
          self%ai = ai
        else
          write (gol,'("could not fill `ai`")'); call goErr
          TRACEBACK; status=1; return
        end if
      end if
      ! storage:
      allocate( self%bi(nlev+1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! select:
      if ( present(bi) ) then
        ! copy:
        self%bi = bi
      else
        write (gol,'("could not fill `bi`")'); call goErr
        TRACEBACK; status=1; return
      end if
    end if  ! with bounds
    
    ! * fill mid values, eventually from bounds
    
    ! switch:
    if ( self%with_ap ) then
      ! storage:
      allocate( self%ap(nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy or compute ?
      if ( present(ap) ) then
        ! copy:
        self%ap = ap
      else if ( self%with_levi ) then
        ! average:
        self%ap = 0.5 * ( self%api(1:nlev) + self%api(2:nlev+1) )
      !else if ( self%with_bnds ) then
      !  ! average:
      !  self%ap = 0.5 * ( self%ap_bnds(1,:) + self%ap_bnds(2,:) )
      else
        write (gol,'("could not fill `ap`")'); call goErr
        TRACEBACK; status=1; return
      end if
    else
      ! storage:
      allocate( self%a(nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy or compute ?
      if ( present(a) ) then
        ! copy:
        self%a = a
      else if ( self%with_levi ) then
        ! average:
        self%a = 0.5 * ( self%ai(1:nlev) + self%ai(2:nlev+1) )
      !else if ( self%with_bnds ) then
      !  ! average:
      !  self%a = 0.5 * ( self%a_bnds(1,:) + self%a_bnds(2,:) )
      else
        write (gol,'("could not fill `a`")'); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! storage:
    allocate( self%b(nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! copy or compute ?
    if ( present(b) ) then
      ! copy:
      self%b = b
    else if ( self%with_levi ) then
      ! average:
      self%b = 0.5 * ( self%bi(1:nlev) + self%bi(2:nlev) )
    !else if ( self%with_bnds ) then
    !  ! average:
    !  self%b = 0.5 * ( self%b_bnds(1,:) + self%b_bnds(2,:) )
    else
      write (gol,'("could not fill `b`")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! * fill p0
    
    ! check ...
    if ( .not. present(p0) ) then
      write (gol,'("argument `p0` not provided")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! copy:
    self%p0 = p0
    
    ! * surface pressure variable:
    
    ! reset from default if provided as argument:
    if ( present(vname_ps) ) self%vname_ps = trim(vname_ps)
    
    ! * fill pressure units
    
    ! check ...
    if ( .not. present(units_p) ) then
      write (gol,'("argument `units_p` not provided")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! copy:
    self%units_p = units_p
    
    ! * fill reference sigma levels

    ! storage:
    allocate( self%lev(nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    if ( self%with_ap ) then
      ! fill:
      self%lev = self%ap / self%p0 + self%b
      ! attribute:
      write (self%long_name,'("hybride level at layer midpoints (",a,"/",a,"+",a,")")') &
              trim(self%vname_ap), trim(self%vname_p0), trim(self%vname_b)
    else
      ! fill:
      self%lev = self%a            + self%b
      ! attribute:
      write (self%long_name,'("hybride level at layer midpoints (",a,"+",a,")")') &
              trim(self%vname_a), trim(self%vname_b)
    end if
    
    !! bounds?
    !if ( self%with_bnds ) then
    !  ! storage:
    !  allocate( self%lev_bnds(nv,nlev), stat=status )
    !  IF_NOT_OK_RETURN(status=1)
    !  ! fill:
    !  if ( self%with_ap ) then
    !    ! fill:
    !    self%lev_bnds = self%ap_bnds / self%p0 + self%b_bnds
    !    ! attribute:
    !    write (self%long_name,'("hybride level at layer boundaries (",a,"_bnds/",a,"+",a,"_bnds)")') &
    !            trim(self%vname_ap), trim(self%vname_p0), trim(self%vname_b)
    !  else
    !    ! fill:
    !    self%lev_bnds = self%a_bnds            + self%b_bnds
    !    ! attribute:
    !    write (self%long_name_bnds,'("hybride level at layer boundaries (",a,"_bnds+",a,"_bnds)")') &
    !          trim(self%vname_a), trim(self%vname_b)
    !  end if
    !end if
    
    ! interfaces?
    if ( self%with_levi ) then
      ! storage:
      allocate( self%levi(nlev+1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! fill:
      if ( self%with_ap ) then
        ! fill:
        self%levi = self%api / self%p0 + self%bi
        ! attribute:
        write (self%long_name_i,'("hybride level at layer interfaces (",a,"i/",a,"+",a,"i)")') &
              trim(self%vname_ap), trim(self%vname_p0), trim(self%vname_b)
      else
        ! fill:
        self%levi = self%ai + self%bi
        ! attribute:
        write (self%long_name_i,'("hybride level at layer interfaces (",a,"i+",a,"i)")') &
              trim(self%vname_a), trim(self%vname_b)
      end if
    end if
    
    ! * fill direction attribute
    
    ! direction:
    if ( self%lev(1) > self%lev(nlev) ) then
      self%positive = 'up'
    else
      self%positive = 'down'
    end if
    
    ! * check consistency
    
    ! check ...
    select case ( trim(self%positive) )
      ! upward:
      case ( 'up' )
        ! mids should decrease:
        if ( any( self%lev(2:nlev) >= self%lev(1:nlev-1) ) ) then
          write (gol,'("upward direction requires strictly decreasing mid levels:")'); call goErr
          do ilev = 1, nlev
            mark = ' '
            if ( ilev > 1 ) then
              if ( self%lev(ilev) >= self%lev(ilev-1) ) then
                mark = '*'
              end if
            end if
            write (gol,'("  level ",i4,f12.4," ",a)') ilev, self%lev(ilev), mark; call goErr
          end do
          TRACEBACK; status=1; return
        end if
        !! bounds?
        !if ( self%with_bnds ) then
        !  ! bounds should decrease:
        !  if ( any( self%lev_bnds(2,:) >= self%lev_bnds(1,:) ) ) then
        !    write (gol,'("upward direction requires strictly decreasing bound levels:")'); call goErr
        !    do ilev = 1, nlev
        !      mark = ' '
        !      if ( self%lev_bnds(2,ilev) >= self%lev_bnds(1,ilev) ) mark = '*'
        !      write (gol,'("  level ",i4,2f12.4," ",a)') ilev, self%lev_bnds(:,ilev), mark; call goErr
        !    end do
        !    TRACEBACK; status=1; return
        !  end if
        !end if
        ! interfaces?
        if ( self%with_levi ) then
          ! bounds should decrease:
          if ( any( self%levi(2:nlev+1) >= self%levi(1:nlev) ) ) then
            write (gol,'("upward direction requires strictly decreasing bound levels:")'); call goErr
            do ilev = 1, nlev
              mark = ' '
              if ( self%levi(ilev+1) >= self%levi(ilev) ) mark = '*'
              write (gol,'("  level ",i4,f12.4," ",a)') ilev, self%levi(ilev), mark; call goErr
            end do
            TRACEBACK; status=1; return
          end if
        end if
      ! downward:
      case ( 'down' )
        ! mids should increase:
        if ( any( self%lev(2:nlev) <= self%lev(1:nlev-1) ) ) then
          write (gol,'("downward direction requires strictly increasing mid levels:")'); call goErr
          do ilev = 1, nlev
            mark = ' '
            if ( ilev > 1 ) then
              if ( self%lev(ilev) <= self%lev(ilev-1) ) mark = '*'
            end if
            write (gol,'("  level ",i4,f12.4," ",a)') ilev, self%lev(ilev), mark; call goErr
          end do
          TRACEBACK; status=1; return
        end if
        !! bounds?
        !if ( self%with_bnds ) then
        !  ! bounds should increase:
        !  if ( any( self%lev_bnds(2,:) <= self%lev_bnds(1,:) ) ) then
        !    write (gol,'("downward direction requires strictly increasing bound levels:")'); call goErr
        !    do ilev = 1, nlev
        !      mark = ' '
        !      if ( self%lev_bnds(2,ilev) <= self%lev_bnds(1,ilev) ) mark = '*'
        !      write (gol,'("  level ",i4,2f12.4," ",a)') ilev, self%lev_bnds(:,ilev), mark; call goErr
        !    end do
        !    TRACEBACK; status=1; return
        !  end if
        !end if
        ! interfaces?
        if ( self%with_levi ) then
          ! bounds should increase:
          if ( any( self%levi(2:nlev+1) <= self%levi(1:nlev) ) ) then
            write (gol,'("downward direction requires strictly increasing bound levels:")'); call goErr
            do ilev = 1, nlev
              mark = ' '
              if ( self%levi(ilev+1) <= self%levi(ilev) ) mark = '*'
              write (gol,'("  level ",i4,f12.4," ",a)') ilev, self%levi(ilev), mark; call goErr
            end do
            TRACEBACK; status=1; return
          end if
        end if
      ! unknown ...
      case default
        write (gol,'("unsupported positive direction `",a,"`")') trim(self%positive); call goErr
        TRACEBACK; status=1; return
    end select
    
    !! bounds?
    !if ( self%with_bnds ) then
    !  ! bounds should follow:
    !  if ( any( self%lev_bnds(1,2:nlev) /= self%lev_bnds(2,1:nlev-1) ) ) then
    !    write (gol,'("boundary values should be equal at level interfaces:")'); call goErr
    !    do ilev = 1, nlev
    !      mark = ' '
    !      if ( ilev > 1 ) then
    !        if ( self%lev_bnds(1,ilev) /= self%lev_bnds(2,ilev-1) ) mark = '*'
    !      end if
    !      write (gol,'("  level ",i4,2f12.4," ",a)') ilev, self%lev_bnds(:,ilev), mark; call goErr
    !    end do
    !    TRACEBACK; status=1; return
    !  end if
    !  ! not in between ?
    !  if ( any( self%lev < minval(self%lev_bnds,dim=1) ) .or. any( self%lev > maxval(self%lev_bnds,dim=1) ) ) then
    !    write (gol,'("lev values not between bounds:")'); call goErr
    !    write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !    write (gol,'("    level         lev                  lev_bnds")'); call goErr
    !    write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !    do ilev = 1, nlev
    !      mark = ' '
    !      if ( (self%lev(ilev) < minval(self%lev_bnds(:,ilev))) .or. (self%lev(ilev) > maxval(self%lev_bnds(:,ilev))) ) mark = '*'
    !      write (gol,'("    ",i5,3f12.4," ",a)') ilev, self%lev(ilev), self%lev_bnds(:,ilev); call goErr
    !    end do
    !    TRACEBACK; status=1; return
    !    write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !  end if
    !  ! switch:
    !  if ( self%with_ap ) then
    !    ! not in between ?
    !    if ( any( self%ap < minval(self%ap_bnds,dim=1) ) .or. any( self%ap > maxval(self%ap_bnds,dim=1) ) ) then
    !      write (gol,'("ap values not between bounds:")'); call goErr
    !      write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !      write (gol,'("    level          ap                   ap_bnds")'); call goErr
    !      write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !      do ilev = 1, nlev
    !        mark = ' '
    !        if ( (self%ap(ilev) < minval(self%ap_bnds(:,ilev))) .or. (self%ap(ilev) > maxval(self%ap_bnds(:,ilev))) ) mark = '*'
    !        write (gol,'("    ",i5,3f12.4," ",a)') ilev, self%ap(ilev), self%ap_bnds(:,ilev); call goErr
    !      end do
    !      TRACEBACK; status=1; return
    !      write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !    end if
    !  else
    !    ! not in between ?
    !    if ( any( self%a < minval(self%a_bnds,dim=1) ) .or. any( self%a > maxval(self%a_bnds,dim=1) ) ) then
    !      write (gol,'("lev values not between bounds:")'); call goErr
    !      write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !      write (gol,'("    level           a                    a_bnds")'); call goErr
    !      write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !      do ilev = 1, nlev
    !        mark = ' '
    !        if ( (self%a(ilev) < minval(self%a_bnds(:,ilev))) .or. (self%a(ilev) > maxval(self%a_bnds(:,ilev))) ) mark = '*'
    !        write (gol,'("    ",i5,3f12.4," ",a)') ilev, self%a(ilev), self%a_bnds(:,ilev); call goErr
    !      end do
    !      TRACEBACK; status=1; return
    !      write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !    end if
    !  end if
    !  ! not in between ?
    !  if ( any( self%b < minval(self%b_bnds,dim=1) ) .or. any( self%b > maxval(self%b_bnds,dim=1) ) ) then
    !    write (gol,'("b values not between bounds:")'); call goErr
    !    write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !    write (gol,'("    level           b                    b_bnds")'); call goErr
    !    write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !    do ilev = 1, nlev
    !      mark = ' '
    !      if ( (self%b(ilev) < minval(self%b_bnds(:,ilev))) .or. (self%b(ilev) > maxval(self%b_bnds(:,ilev))) ) mark = '*'
    !      write (gol,'("    ",i5,3f12.4," ",a)') ilev, self%b(ilev), self%b_bnds(:,ilev); call goErr
    !    end do
    !    TRACEBACK; status=1; return
    !    write (gol,'("    -----  ----------  ------------------------")'); call goErr
    !  end if
    !end if
    
    ! * fill formula attributes
    
    ! set formula terms:
    if ( self%with_ap ) then
      write (self%formula,'("p(n,k,j,i) = ap(k) + b(k) * ps(n,j,i)")')
      write (self%formula_terms,'("ap: ",a," b: ",a," ps: ",a," p0: ",a)') &
               trim(self%vname_ap), trim(self%vname_b), trim(self%vname_ps), trim(self%vname_p0)
      !write (self%formula_terms_bnds,'("ap: ",a,"_bnds b: ",a,"_bnds ps: ",a," p0: ",a)') &
      !         trim(self%vname_ap), trim(self%vname_b), trim(self%vname_ps), trim(self%vname_p0)
      write (self%formula_terms_i,'("ap: ",a,"i b: ",a,"i ps: ",a," p0: ",a)') &
               trim(self%vname_ap), trim(self%vname_b), trim(self%vname_ps), trim(self%vname_p0)
    else
      write (self%formula,'("p(n,k,j,i) = a(k) * p0 + b(k) * ps(n,j,i)")')
      write (self%formula_terms,'("a: ",a," b: ",a," ps: ",a," p0: ",a)') &
               trim(self%vname_a), trim(self%vname_b), trim(self%vname_ps), trim(self%vname_p0)
      !write (self%formula_terms_bnds,'("a: ",a,"_bnds b: ",a,"_bnds ps: ",a," p0: ",a)') &
      !         trim(self%vname_a), trim(self%vname_b), trim(self%vname_ps), trim(self%vname_p0)
      write (self%formula_terms_i,'("a: ",a,"i b: ",a,"i ps: ",a," p0: ",a)') &
               trim(self%vname_a), trim(self%vname_b), trim(self%vname_ps), trim(self%vname_p0)
    end if
    
    ! *
    
    ! ok
    status = 0
    
  end subroutine HybrideLevelCoordinate_Set_Values


  ! ***
  

  subroutine HybrideLevelCoordinate_Get_Values( self, status, &
                                                  a, ap, b, &
                                                  ai, api, bi, &
                                                  !a_bnds, ap_bnds, b_bnds, &
                                                  with_ap, with_levi, & !with_bnds, &
                                                  p0, units_p, vname_ps )
  
    ! --- in/out ---------------------------------
    
    XTYPE(HybrideLevelCoordinate), intent(in)     ::  self
    integer, intent(out)                          ::  status

    real, intent(out), optional                   ::  a(:), ap(:), b(:)
    real, intent(out), optional                   ::  ai(:), api(:), bi(:)
    !real, intent(out), optional                   ::  a_bnds(:,:), ap_bnds(:,:), b_bnds(:,:)
    logical, intent(out), optional                ::  with_ap
    logical, intent(out), optional                ::  with_levi
    !logical, intent(out), optional                ::  with_bnds
    real, intent(out), optional                   ::  p0
    character(len=*), intent(out), optional       ::  units_p
    character(len=*), intent(out), optional       ::  vname_ps

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/HybrideLevelCoordinate_Get_Values'
    
    ! --- local ----------------------------------
    
    integer     ::  nlev
    
    ! --- begin ----------------------------------
    
    ! flags:
    if ( present(with_ap  ) ) with_ap   = self%with_ap
    if ( present(with_levi) ) with_levi = self%with_levi
    !if ( present(with_bnds) ) with_bnds = self%with_bnds

    ! get size:
    call self%Coordinate%Get_Dim( status, n=nlev )
    IF_NOT_OK_RETURN(status=1)
    
    ! value ?
    if ( present(a) ) then
      ! check ...
      if ( self%with_ap ) then
        write (gol,'("requested `a` while in `ap` mode")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( .not. allocated(self%a) ) then
        write (gol,'("no `a` stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( any( shape(a) /= (/nlev/) ) ) then
        write (gol,'("argument `a` has shape (",i0,") while nlev is ",i0)') shape(a), nlev; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      a = self%a
    end if
    
    ! value ?
    if ( present(ap) ) then
      ! check ...
      if ( .not. self%with_ap ) then
        write (gol,'("requested `ap` while not in `ap` mode")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( .not. allocated(self%ap) ) then
        write (gol,'("no `ap` stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( any( shape(ap) /= (/nlev/) ) ) then
        write (gol,'("argument `a` has shape (",i0,") while nlev is ",i0)') shape(ap), nlev; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      ap = self%ap
    end if
    
    ! value ?
    if ( present(b) ) then
      ! check ...
      if ( .not. allocated(self%b) ) then
        write (gol,'("no `b` stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( any( shape(b) /= (/nlev/) ) ) then
        write (gol,'("argument `b` has shape (",i0,") while nlev is ",i0)') shape(a), nlev; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      b = self%b
    end if
    
    ! value ?
    if ( present(ai) ) then
      ! check ...
      if ( self%with_ap ) then
        write (gol,'("requested `ai` while in `ap` mode")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( .not. allocated(self%ai) ) then
        write (gol,'("no `ai` stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( any( shape(ai) /= (/nlev+1/) ) ) then
        write (gol,'("argument `ai` has shape (",i0,") while nlev+1 is ",i0)') shape(ai), nlev+1; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      ai = self%ai
    end if
    
    ! value ?
    if ( present(api) ) then
      ! check ...
      if ( .not. self%with_ap ) then
        write (gol,'("requested `api` while not in `ap` mode")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( .not. allocated(self%api) ) then
        write (gol,'("no `api` stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( any( shape(api) /= (/nlev+1/) ) ) then
        write (gol,'("argument `a` has shape (",i0,") while nlev+1 is ",i0)') shape(api), nlev+1; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      api = self%api
    end if
    
    ! value ?
    if ( present(bi) ) then
      ! check ...
      if ( .not. allocated(self%bi) ) then
        write (gol,'("no `bi` stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( any( shape(bi) /= (/nlev+1/) ) ) then
        write (gol,'("argument `bi` has shape (",i0,") while nlev+1 is ",i0)') shape(ai), nlev+1; call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      bi = self%bi
    end if
    
    !! value ?
    !if ( present(a_bnds) ) then
    !  ! check ...
    !  if ( self%with_ap ) then
    !    write (gol,'("requested `a_bnds` while in `ap` mode")'); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! check ...
    !  if ( .not. allocated(self%a_bnds) ) then
    !    write (gol,'("no `a_bnds` stored yet")'); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! check ...
    !  if ( any( shape(a_bnds) /= (/nv,nlev/) ) ) then
    !    write (gol,'("argument `a_bnds` has shape (",i0,",",i0,") while nlev is ",i0)') shape(a_bnds), nlev; call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! copy:
    !  a_bnds = self%a_bnds
    !end if
    
    !! value ?
    !if ( present(ap_bnds) ) then
    !  ! check ...
    !  if ( .not. self%with_ap ) then
    !    write (gol,'("requested `ap_bnds` while not in `ap` mode")'); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! check ...
    !  if ( .not. allocated(self%ap_bnds) ) then
    !    write (gol,'("no `ap_bnds` stored yet")'); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! check ...
    !  if ( any( shape(ap_bnds) /= (/nv,nlev/) ) ) then
    !    write (gol,'("argument `a_bnds` has shape (",i0,",",i0,") while nlev is ",i0)') shape(ap_bnds), nlev; call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! copy:
    !  ap_bnds = self%ap_bnds
    !end if
    
    !! value ?
    !if ( present(b_bnds) ) then
    !  ! check ...
    !  if ( .not. allocated(self%b_bnds) ) then
    !    write (gol,'("no `b_bnds` stored yet")'); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! check ...
    !  if ( any( shape(b_bnds) /= (/nv,nlev/) ) ) then
    !    write (gol,'("argument `b_bnds` has shape (",i0,",",i0,") while nlev is ",i0)') shape(b_bnds), nlev; call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! copy:
    !  b_bnds = self%b_bnds
    !end if
    
    ! value?
    if ( present(p0) ) p0 = self%p0
    
    ! value?
    if ( present(units_p) ) then
      ! check ...
      if ( len_trim(self%units_p) == 0 ) then
        write (gol,'("no `units_p` defined yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      units_p = trim(self%units_p)
    end if
    
    ! value?
    if ( present(vname_ps) ) vname_ps = trim(self%vname_ps)
    
    ! ok
    status = 0
    
  end subroutine HybrideLevelCoordinate_Get_Values


  ! ***


  subroutine HybrideLevelCoordinate_Def( self, file, status )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Dim, NF90_Inq_DimID
    use NetCDF, only : NF90_Def_Var  
    use NetCDF, only : NF90_Put_Att
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(HybrideLevelCoordinate), intent(inout)    ::  self
    class(Datafile), intent(in)                     ::  file
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/HybrideLevelCoordinate_Def'
    
    ! --- local ----------------------------------

    character(len=LEN_NAME)   ::  name
    integer                   ::  dimid_lev
    integer                   ::  dimid_levi
    integer                   ::  dimid_nv
    integer                   ::  nlev
    
    ! --- begin ----------------------------------
    
    ! define dimension:
    call self%Dimension%Def( file, status )
    IF_NOT_OK_RETURN(status=1)
    ! get access id:
    call self%Dimension%Get_Dim( status, name=name, n=nlev, dimid=dimid_lev )
    IF_NOT_OK_RETURN(status=1)
    
    ! define new variable:
    status = NF90_Def_Var( file%ncid, trim(name), NF90_FLOAT, (/dimid_lev/), self%varid_out )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! define default attributes:
    call self%Coordinate_Def( file, self%varid_out, status )
    IF_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Put_Att( file%ncid, self%varid_out, 'positive', trim(self%positive) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Put_Att( file%ncid, self%varid_out, 'formula_terms', trim(self%formula_terms) )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Put_Att( file%ncid, self%varid_out, 'formula', trim(self%formula) )
    IF_NF90_NOT_OK_RETURN(status=1)
    !! bounds ?
    !if ( self%with_bnds ) then
    !  ! attribute:
    !  status = NF90_Put_Att( file%ncid, self%varid_out, 'bounds', trim(name)//'_bnds' )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !end if

    ! mode?
    if ( self%with_ap ) then
      ! define new variable:
      status = NF90_Def_Var( file%ncid, trim(self%vname_ap), NF90_FLOAT, (/dimid_lev/), self%varid_out_ap )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_ap, 'long_name', 'hybrid ap coefficient at layer midpoints' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_ap, 'units', trim(self%units_p) )
      IF_NF90_NOT_OK_RETURN(status=1)
      !! bounds ?
      !if ( self%with_bnds ) then
      !  ! attribute:
      !  status = NF90_Put_Att( file%ncid, self%varid_out_ap, 'bounds', trim(self%vname_ap)//'_bnds' )
      !  IF_NF90_NOT_OK_RETURN(status=1)
      !end if
    else
      ! define new variable:
      status = NF90_Def_Var( file%ncid, trim(self%vname_a), NF90_FLOAT, (/dimid_lev/), self%varid_out_a )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_a, 'long_name', 'hybrid a coefficient at layer midpoints' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_a, 'units', '1' )
      IF_NF90_NOT_OK_RETURN(status=1)
      !! bounds ?
      !if ( self%with_bnds ) then
      !  ! attribute:
      !  status = NF90_Put_Att( file%ncid, self%varid_out_a, 'bounds', trim(self%vname_a)//'_bnds' )
      !  IF_NF90_NOT_OK_RETURN(status=1)
      !end if
    end if

    ! define new variable:
    status = NF90_Def_Var( file%ncid, trim(self%vname_b), NF90_FLOAT, (/dimid_lev/), self%varid_out_b )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Put_Att( file%ncid, self%varid_out_b, 'long_name', 'hybrid b coefficient at layer midpoints' )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Put_Att( file%ncid, self%varid_out_b, 'units', '1' )
    IF_NF90_NOT_OK_RETURN(status=1)
    !! bounds ?
    !if ( self%with_bnds ) then
    !  ! attribute:
    !  status = NF90_Put_Att( file%ncid, self%varid_out_b, 'bounds', trim(self%vname_b)//'_bnds' )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !end if

    !! bounds ?
    !if ( self%with_bnds ) then
    !
    !  ! create dimension if not done yet:
    !  status = NF90_INQ_DimID( file%ncid, dname_nv, dimid_nv )
    !  if (status/=NF90_NOERR) then
    !    status = NF90_Def_Dim( file%ncid, dname_nv, nv, dimid_nv )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !  end if
    !
    !  ! define new variable:
    !  status = NF90_Def_Var( file%ncid, trim(name)//'_bnds', NF90_FLOAT, (/dimid_nv,dimid_lev/), self%varid_out_bnds )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !  ! attribute:
    !  status = NF90_Put_Att( file%ncid, self%varid_out_bnds, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate' )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !  ! attribute:
    !  status = NF90_Put_Att( file%ncid, self%varid_out_bnds, 'long_name', trim(self%long_name_bnds) )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !  ! attribute:
    !  status = NF90_Put_Att( file%ncid, self%varid_out_bnds, 'units', '1' )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !  !! attribute:
    !  !status = NF90_Put_Att( file%ncid, self%varid_out_bnds, 'positive', trim(self%positive) )
    !  !IF_NF90_NOT_OK_RETURN(status=1)
    !  !! attribute:
    !  !status = NF90_Put_Att( file%ncid, self%varid_out_bnds, 'formula_terms', trim(self%formula_terms_bnds) )
    !  !IF_NF90_NOT_OK_RETURN(status=1)
    !  !! attribute:
    !  !status = NF90_Put_Att( file%ncid, self%varid_out_bnds, 'formula', trim(self%formula) )
    !  !IF_NF90_NOT_OK_RETURN(status=1)
    !
    !  ! mode?
    !  if ( self%with_ap ) then
    !    ! define new variable:
    !    status = NF90_Def_Var( file%ncid, trim(self%vname_ap)//'_bnds', NF90_FLOAT, (/dimid_nv,dimid_lev/), self%varid_out_ap_bnds )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !    ! attribute:
    !    status = NF90_Put_Att( file%ncid, self%varid_out_ap_bnds, 'long_name', 'hybrid a coefficient at layer boundaries' )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !    ! attribute:
    !    status = NF90_Put_Att( file%ncid, self%varid_out_ap_bnds, 'units', trim(self%units_p) )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !  else
    !    ! define new variable:
    !    status = NF90_Def_Var( file%ncid, trim(self%vname_a)//'_bnds', NF90_FLOAT, (/dimid_nv,dimid_lev/), self%varid_out_a_bnds )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !    ! attribute:
    !    status = NF90_Put_Att( file%ncid, self%varid_out_a_bnds, 'long_name', 'hybrid a coefficient at layer boundaries' )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !    ! attribute:
    !    status = NF90_Put_Att( file%ncid, self%varid_out_a_bnds, 'units', '1' )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !  end if
    !
    !  ! define new variable:
    !  status = NF90_Def_Var( file%ncid, trim(self%vname_b)//'_bnds', NF90_FLOAT, (/dimid_nv,dimid_lev/), self%varid_out_b_bnds )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !  ! attribute:
    !  status = NF90_Put_Att( file%ncid, self%varid_out_b_bnds, 'long_name', 'hybrid b coefficient at layer boundaries' )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !  ! attribute:
    !  status = NF90_Put_Att( file%ncid, self%varid_out_b_bnds, 'units', '1' )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !
    !end if ! bounds
    
    ! interfaces?
    if ( self%with_levi ) then
    
      ! create dimension if not done yet:
      status = NF90_INQ_DimID( file%ncid, trim(name)//'i', dimid_levi )
      if (status/=NF90_NOERR) then
        status = NF90_Def_Dim( file%ncid, trim(name)//'i', nlev+1, dimid_levi )
        IF_NF90_NOT_OK_RETURN(status=1)
      end if
    
      ! define new variable:
      status = NF90_Def_Var( file%ncid, trim(name)//'i', NF90_FLOAT, (/dimid_levi/), self%varid_out_i )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_i, 'standard_name', 'atmosphere_hybrid_sigma_pressure_coordinate' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_i, 'long_name', trim(self%long_name_i) )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_i, 'units', '1' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_i, 'positive', trim(self%positive) )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_i, 'formula_terms', trim(self%formula_terms_i) )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_i, 'formula', trim(self%formula) )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! mode?
      if ( self%with_ap ) then
        ! define new variable:
        status = NF90_Def_Var( file%ncid, trim(self%vname_ap)//'i', NF90_FLOAT, (/dimid_levi/), self%varid_out_api )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! attribute:
        status = NF90_Put_Att( file%ncid, self%varid_out_api, 'long_name', 'hybrid a coefficient at layer interfaces' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! attribute:
        status = NF90_Put_Att( file%ncid, self%varid_out_api, 'units', trim(self%units_p) )
        IF_NF90_NOT_OK_RETURN(status=1)
      else
        ! define new variable:
        status = NF90_Def_Var( file%ncid, trim(self%vname_a)//'i', NF90_FLOAT, (/dimid_levi/), self%varid_out_ai )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! attribute:
        status = NF90_Put_Att( file%ncid, self%varid_out_ai, 'long_name', 'hybrid a coefficient at layer interfaces' )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! attribute:
        status = NF90_Put_Att( file%ncid, self%varid_out_ai, 'units', '1' )
        IF_NF90_NOT_OK_RETURN(status=1)
      end if

      ! define new variable:
      status = NF90_Def_Var( file%ncid, trim(self%vname_b)//'i', NF90_FLOAT, (/dimid_levi/), self%varid_out_bi )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_bi, 'long_name', 'hybrid b coefficient at layer interfaces' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out_bi, 'units', '1' )
      IF_NF90_NOT_OK_RETURN(status=1)
      
    end if  ! bounds

    ! define scalar variable, no 'dimids' array is passed:
    status = NF90_Def_Var( file%ncid, trim(self%vname_p0), NF90_FLOAT, self%varid_out_p0 )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Put_Att( file%ncid, self%varid_out_p0, 'long_name', 'reference surface pressure' )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Put_Att( file%ncid, self%varid_out_p0, 'units', trim(self%units_p) )
    IF_NF90_NOT_OK_RETURN(status=1)

    
    ! ok
    status = 0
    
  end subroutine HybrideLevelCoordinate_Def


  ! ***
  

  subroutine HybrideLevelCoordinate_Write( self, file, status )

    use NetCDF, only : NF90_Put_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(HybrideLevelCoordinate), intent(in)   ::  self
    class(Datafile), intent(in)         ::  file
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/HybrideLevelCoordinate_Write'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! write general coordinate:
    call self%Coordinate_Write( file, status )
    IF_NOT_OK_RETURN(status=1)

    ! check ...
    if ( .not. allocated(self%lev) ) then
      write (gol,'("no `lev` values stored yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! write:
    status = NF90_Put_Var( file%ncid, self%varid_out, self%lev )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! switch:
    if ( self%with_ap ) then
      ! check ...
      if ( .not. allocated(self%ap) ) then
        write (gol,'("no `ap` values stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! write:
      status = NF90_Put_Var( file%ncid, self%varid_out_ap, self%ap )
      IF_NF90_NOT_OK_RETURN(status=1)
    else
      ! check ...
      if ( .not. allocated(self%a) ) then
        write (gol,'("no `a` values stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! write:
      status = NF90_Put_Var( file%ncid, self%varid_out_a, self%a )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! check ...
    if ( .not. allocated(self%b) ) then
      write (gol,'("no `b` values stored yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! write:
    status = NF90_Put_Var( file%ncid, self%varid_out_b, self%b )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    !! boundaries?
    !if ( self%with_bnds ) then
    !
    !  ! check ...
    !  if ( .not. allocated(self%lev_bnds) ) then
    !    write (gol,'("no `lev_bnds` values stored yet")'); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! write:
    !  status = NF90_Put_Var( file%ncid, self%varid_out_bnds, self%lev_bnds )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !
    !  ! switch:
    !  if ( self%with_ap ) then
    !    ! check ...
    !    if ( .not. allocated(self%ap_bnds) ) then
    !      write (gol,'("no `ap_bnds` values stored yet")'); call goErr
    !      TRACEBACK; status=1; return
    !    end if
    !    ! write:
    !    status = NF90_Put_Var( file%ncid, self%varid_out_ap_bnds, self%ap_bnds )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !  else
    !    ! check ...
    !    if ( .not. allocated(self%a_bnds) ) then
    !      write (gol,'("no `a_bnds` values stored yet")'); call goErr
    !      TRACEBACK; status=1; return
    !    end if
    !    ! write:
    !    status = NF90_Put_Var( file%ncid, self%varid_out_a_bnds, self%a_bnds )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !  end if
    !
    !  ! check ...
    !  if ( .not. allocated(self%b_bnds) ) then
    !    write (gol,'("no `b_bnds` values stored yet")'); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! write:
    !  status = NF90_Put_Var( file%ncid, self%varid_out_b_bnds, self%b_bnds )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !
    !end if ! bounds
    
    ! interfaces
    if ( self%with_levi ) then

      ! check ...
      if ( .not. allocated(self%levi) ) then
        write (gol,'("no `levi` values stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! write:
      status = NF90_Put_Var( file%ncid, self%varid_out_i, self%levi )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! switch:
      if ( self%with_ap ) then
        ! check ...
        if ( .not. allocated(self%api) ) then
          write (gol,'("no `api` values stored yet")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! write:
        status = NF90_Put_Var( file%ncid, self%varid_out_api, self%api )
        IF_NF90_NOT_OK_RETURN(status=1)
      else
        ! check ...
        if ( .not. allocated(self%ai) ) then
          write (gol,'("no `ai` values stored yet")'); call goErr
          TRACEBACK; status=1; return
        end if
        ! write:
        status = NF90_Put_Var( file%ncid, self%varid_out_ai, self%ai )
        IF_NF90_NOT_OK_RETURN(status=1)
      end if

      ! check ...
      if ( .not. allocated(self%bi) ) then
        write (gol,'("no `bi` values stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! write:
      status = NF90_Put_Var( file%ncid, self%varid_out_bi, self%bi )
      IF_NF90_NOT_OK_RETURN(status=1)
      
    end if  ! bounds
    
    ! write:
    status = NF90_Put_Var( file%ncid, self%varid_out_p0, self%p0 )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine HybrideLevelCoordinate_Write


  ! ***
  

  subroutine HybrideLevelCoordinate_Read( self, file, status )

    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Inq_AttName, NF90_Get_Att
  
    use GO           , only : goReadFromLine
    use C3PO_Datafile, only : Datafile

    ! --- in/out ---------------------------------
    
    XTYPE(HybrideLevelCoordinate), intent(inout)    ::  self
    class(Datafile), intent(in)                     ::  file
    integer, intent(out)                            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/HybrideLevelCoordinate_Read'
    
    ! --- local ----------------------------------

    integer                   ::  nlev
    integer                   ::  varid
    real, allocatable         ::  a(:), ap(:), b(:)    ! (nlev)
    real, allocatable         ::  ai(:), api(:), bi(:)    ! (nlev+1)
    !real, allocatable         ::  a_bnds(:,:), ap_bnds(:,:), b_bnds(:,:)  ! (nv,nlev)
    real                      ::  p0
    character(len=LEN_UNITS)  ::  units_p
    character(len=LEN_LINE)   ::  formula_terms
    character(len=LEN_LINE)   ::  line
    character(len=LEN_NAME)   ::  name
    character(len=LEN_NAME)   ::  key, value
    !integer                   ::  natt, iatt
    !character(len=LEN_NAME)   ::  attname
    logical                   ::  with_ap
    logical                   ::  with_levi
    !logical                   ::  with_bnds
    
    ! --- begin ----------------------------------
    
    ! read dimension parameters (length):
    call self%Dimension%Read( file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! obtain name and size:
    call self%Dimension%Get_Dim( status, name=name, n=nlev )
    IF_NOT_OK_RETURN(status=1)

    ! get variable id:
    status = NF90_INQ_VarID( file%ncid, trim(name), varid )
    if ( status /= NF90_NOERR) then
      gol=NF90_StrError(status); call goErr
      write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
              trim(name), trim(file%filename); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! attribute:
    status = NF90_Get_Att( file%ncid, varid, 'units', self%units )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Get_Att( file%ncid, varid, 'standard_name', self%standard_name )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! attribute:
    status = NF90_Get_Att( file%ncid, varid, 'formula_terms', formula_terms )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! reset names:
    self%vname_ap = ''
    self%vname_a  = ''
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
        case ( 'a:' )
          self%vname_a = trim(value)
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
    if ( (len_trim(self%vname_ap) == 0) .and. (len_trim(self%vname_a) == 0) ) then
      write (gol,'("neither ap or a defined in formula terms : ",a)') trim(formula_terms); call goErr
      TRACEBACK; status=1; return
    end if
    ! check ...
    if ( (len_trim(self%vname_ap) > 0) .and. (len_trim(self%vname_a) > 0) ) then
      write (gol,'("both ap and a defined in formula terms : ",a)') trim(formula_terms); call goErr
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
      write (gol,'("no p0 defined in formula terms : ",a)') trim(formula_terms); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! set flag:
    with_ap = len_trim(self%vname_ap) > 0

    !! check if bounds are present ...
    !with_bnds = .false.
    !! count number of attributes:
    !status = NF90_Inquire_Variable( file%ncid, varid, nAtts=natt )
    !IF_NF90_NOT_OK_RETURN(status=1)
    !! loop over attributes:
    !do iatt = 1, natt
    !  ! name:
    !  status = NF90_Inq_AttName( file%ncid, varid, iatt, attname )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !  ! compare:
    !  with_bnds = trim(attname) == 'bounds'
    !  ! leave?
    !  if ( with_bnds ) exit
    !end do
    
    ! check if interface levels are present:
    status = NF90_INQ_VarID( file%ncid, trim(name)//'i', varid )
    with_levi = status == NF90_NOERR

    ! mode ?
    if ( with_ap ) then
      ! storage:
      allocate( ap(nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! get variable id:
      status = NF90_INQ_VarID( file%ncid, trim(self%vname_ap), varid )
      if ( status /= NF90_NOERR) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
                trim(self%vname_ap), trim(file%filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! read:
      status = NF90_Get_Var( file%ncid, varid, ap )
      IF_NF90_NOT_OK_RETURN(status=1)
    else
      ! storage:
      allocate( a(nlev), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! get variable id:
      status = NF90_INQ_VarID( file%ncid, trim(self%vname_a), varid )
      if ( status /= NF90_NOERR) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
                trim(self%vname_a), trim(file%filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! read:
      status = NF90_Get_Var( file%ncid, varid, a )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if

    ! storage:
    allocate( b(nlev), stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! get variable id:
    status = NF90_INQ_VarID( file%ncid, trim(self%vname_b), varid )
    if ( status /= NF90_NOERR) then
      gol=NF90_StrError(status); call goErr
      write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
              trim(self%vname_b), trim(file%filename); call goErr
      TRACEBACK; status=1; return
    end if
    ! read:
    status = NF90_Get_Var( file%ncid, varid, b )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    !! bounds ?
    !if ( with_bnds ) then
    !
    !  ! mode ?
    !  if ( with_ap ) then
    !    ! storage:
    !    allocate( ap_bnds(nv,nlev), stat=status )
    !    IF_NOT_OK_RETURN(status=1)
    !    ! get variable id:
    !    status = NF90_INQ_VarID( file%ncid, trim(self%vname_ap)//'_bnds', varid )
    !    if ( status /= NF90_NOERR) then
    !      gol=NF90_StrError(status); call goErr
    !      write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
    !              trim(self%vname_ap)//'_bnds', trim(file%filename); call goErr
    !      TRACEBACK; status=1; return
    !    end if
    !    ! read:
    !    status = NF90_Get_Var( file%ncid, varid, ap_bnds )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !  else
    !    ! storage:
    !    allocate( a_bnds(nv,nlev), stat=status )
    !    IF_NOT_OK_RETURN(status=1)
    !    ! get variable id:
    !    status = NF90_INQ_VarID( file%ncid, trim(self%vname_a)//'_bnds', varid )
    !    if ( status /= NF90_NOERR) then
    !      gol=NF90_StrError(status); call goErr
    !      write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
    !              trim(self%vname_a)//'_bnds', trim(file%filename); call goErr
    !      TRACEBACK; status=1; return
    !    end if
    !    ! read:
    !    status = NF90_Get_Var( file%ncid, varid, a_bnds )
    !    IF_NF90_NOT_OK_RETURN(status=1)
    !  end if
    !
    !  ! storage:
    !  allocate( b_bnds(nv,nlev), stat=status )
    !  IF_NOT_OK_RETURN(status=1)
    !  ! get variable id:
    !  status = NF90_INQ_VarID( file%ncid, trim(self%vname_b)//'_bnds', varid )
    !  if ( status /= NF90_NOERR) then
    !    gol=NF90_StrError(status); call goErr
    !    write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
    !            trim(self%vname_b)//'_bnds', trim(file%filename); call goErr
    !    TRACEBACK; status=1; return
    !  end if
    !  ! read:
    !  status = NF90_Get_Var( file%ncid, varid, b_bnds )
    !  IF_NF90_NOT_OK_RETURN(status=1)
    !
    !end if  ! bounds

    ! interfaces ?
    if ( with_levi ) then

      ! mode ?
      if ( with_ap ) then
        ! storage:
        allocate( api(nlev+1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! get variable id:
        status = NF90_INQ_VarID( file%ncid, trim(self%vname_ap)//'i', varid )
        if ( status /= NF90_NOERR) then
          gol=NF90_StrError(status); call goErr
          write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
                  trim(self%vname_ap)//'i', trim(file%filename); call goErr
          TRACEBACK; status=1; return
        end if
        ! read:
        status = NF90_Get_Var( file%ncid, varid, api )
        IF_NF90_NOT_OK_RETURN(status=1)
      else
        ! storage:
        allocate( ai(nlev+1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! get variable id:
        status = NF90_INQ_VarID( file%ncid, trim(self%vname_a)//'i', varid )
        if ( status /= NF90_NOERR) then
          gol=NF90_StrError(status); call goErr
          write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
                  trim(self%vname_b)//'i', trim(file%filename); call goErr
          TRACEBACK; status=1; return
        end if
        ! read:
        status = NF90_Get_Var( file%ncid, varid, ai )
        IF_NF90_NOT_OK_RETURN(status=1)
      end if

      ! storage:
      allocate( bi(nlev+1), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! get variable id:
      status = NF90_INQ_VarID( file%ncid, trim(self%vname_b)//'i', varid )
      if ( status /= NF90_NOERR) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
                trim(self%vname_b)//'i', trim(file%filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! read:
      status = NF90_Get_Var( file%ncid, varid, bi )
      IF_NF90_NOT_OK_RETURN(status=1)

    end if  ! interfaces

    ! get variable id:
    status = NF90_INQ_VarID( file%ncid, trim(self%vname_p0), varid )
    if ( status /= NF90_NOERR) then
      gol=NF90_StrError(status); call goErr
      write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
              trim(self%vname_p0), trim(file%filename); call goErr
      TRACEBACK; status=1; return
    end if
    ! read:
    status = NF90_Get_Var( file%ncid, varid, p0 )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Get_Att( file%ncid, varid, 'units', units_p )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! setup:
    !if ( with_bnds ) then
    !  if ( with_ap ) then
    !    ! setup:
    !    call self%Set_Values( status, ap=ap, b=b, p0=p0, units_p=trim(units_p), &
    !                             ap_bnds=ap_bnds, b_bnds=b_bnds )
    !    IF_NOT_OK_RETURN(status=1)
    !  else
    !    ! setup:
    !    call self%Set_Values( status, a=a, b=b, p0=p0, units_p=trim(units_p), &
    !                             ap_bnds=ap_bnds, b_bnds=b_bnds )
    !    IF_NOT_OK_RETURN(status=1)
    !  end if
    if ( with_levi ) then
      if ( with_ap ) then
        ! setup:
        call self%Set_Values( status, ap=ap, b=b, p0=p0, units_p=trim(units_p), &
                                 api=api, bi=bi )
        IF_NOT_OK_RETURN(status=1)
      else
        ! setup:
        call self%Set_Values( status, a=a, b=b, p0=p0, units_p=trim(units_p), &
                                 ai=ai, bi=bi )
        IF_NOT_OK_RETURN(status=1)
      end if
    else
      if ( with_ap ) then
        ! setup:
        call self%Set_Values( status, ap=ap, b=b, p0=p0, units_p=trim(units_p) )
        IF_NOT_OK_RETURN(status=1)
      else
        ! setup:
        call self%Set_Values( status, a=a, b=b, p0=p0, units_p=trim(units_p) )
        IF_NOT_OK_RETURN(status=1)
      end if
    end if
    
    ! clear:
    if ( with_ap ) then
      deallocate( ap, stat=status )
      IF_NOT_OK_RETURN(status=1)
    else
      deallocate( a, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    deallocate( b, stat=status )
    IF_NOT_OK_RETURN(status=1)
    if ( with_levi ) then
      if ( with_ap ) then
        deallocate( api, stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        deallocate( ai, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      deallocate( bi, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    !if ( with_bnds ) then
    !  if ( with_ap ) then
    !    deallocate( ap_bnds, stat=status )
    !    IF_NOT_OK_RETURN(status=1)
    !  else
    !    deallocate( a_bnds, stat=status )
    !    IF_NOT_OK_RETURN(status=1)
    !  end if
    !  deallocate( b_bnds, stat=status )
    !  IF_NOT_OK_RETURN(status=1)
    !end if

    ! ok
    status = 0
    
  end subroutine HybrideLevelCoordinate_Read


  ! ********************************************************************
  ! ***
  ! *** time coordinate
  ! ***
  ! ********************************************************************


  subroutine TimeCoordinate_Init( self, name, status )
  
    use GO, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(out)    ::  self
    character(len=*), intent(in)                  ::  name
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! init general coordinate:
    call self%Coordinate%Init( name, status )
    
    ! undefined:
    self%climatology       = .false.
    self%calendar          = ''
    self%units_step        = ''
    self%units_t0          = AnyDate()
    self%vname_climat      = ''
    
    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Init


  ! ***


  subroutine TimeCoordinate_Done( self, status )
  
    use GO, only : AnyDate
  
    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(inout)  ::  self
    integer, intent(out)                          ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! general coordinate:
    call self%Coordinate%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    if ( allocated(self%values) ) then
      deallocate( self%values, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( allocated(self%climat) ) then
      deallocate( self%climat, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! undefined:
    self%climatology       = .false.
    self%calendar          = ''
    self%units_step        = ''
    self%units_t0          = AnyDate()
    self%vname_climat      = ''
    
    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Done


  ! ***
  
  
  subroutine TimeCoordinate_Set_Values( self, status, &
                                          units, units_step, units_t0, &
                                          calendar, climatology )
  
    use GO, only : TDate, IsAnyDate, TIncrDate, Extract_Ref_and_Step
    
    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Set_Values'
    
    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(inout)           ::  self
    integer, intent(out)                           ::  status

    character(len=*), intent(in), optional         ::  units      ! 'days since 2000-01-01'
    character(len=*), intent(in), optional         ::  units_step ! 'days'
    type(TDate), intent(in), optional              ::  units_t0   !  TDate(2000,1,1)
    character(len=*), intent(in), optional         ::  calendar
    logical, intent(in), optional                  ::  climatology

    ! --- local ----------------------------------
    
    type(TIncrDate)     ::  tstep
    
    ! --- begin ----------------------------------
    
    ! pre-defined:
    call self%Coordinate%Set_Attrs( status, standard_name='time' )
    IF_NOT_OK_RETURN(status=1)

    ! * store units
    
    ! store arguments:
    if ( present(units     ) ) self%units      = trim(units)
    if ( present(units_step) ) self%units_step = trim(units_step)
    if ( present(units_t0  ) ) self%units_t0   = units_t0
    if ( present(calendar  ) ) self%calendar   = trim(calendar)
    
    ! * check and set units
    
    ! either full 'units' should be defined or the individual elements:
    if ( len_trim(self%units) > 0 ) then
      ! others both undefined ? then fill from units:
      if ( (len_trim(self%units_step) == 0) .and. IsAnyDate(self%units_t0) ) then
        ! extract elements:
        call Extract_Ref_and_Step( self%units, self%units_t0, tstep, status, &
                                     stepunits=self%units_step )
        IF_NOT_OK_RETURN(status=1)
      ! at least one undefined? both undefined already trapped:
      else if ( (len_trim(self%units_step) == 0) .or. IsAnyDate(self%units_t0) ) then
        write (gol,'("`units` defined while either `units_step` or `units_t0` undefined")'); call goErr
        TRACEBACK; status=1; return
      end if
      !
    ! individual elements defined ?
    else if ( (len_trim(self%units_step) > 0) .and. (.not. IsAnyDate(self%units_t0)) ) then
      ! fill units as either:
      !   days since 2000-01-01
      !   hours|minutes|seconds since 2000-01-01 00:00:00
      if ( (trim(self%units_step) == 'days') .and. &
             all( (/self%units_t0%hour,self%units_t0%min,self%units_t0%sec/) == 0 ) ) then
        write (self%units,'(a," since ",i4,2("-",i2.2))') &
                 trim(self%units_step), &
                 self%units_t0%year, self%units_t0%month, self%units_t0%day
      else
        write (self%units,'(a," since ",i4,2("-",i2.2)," ",i2.2,2(":",i2.2))') &
                 trim(self%units_step), &
                 self%units_t0%year, self%units_t0%month, self%units_t0%day, &
                 self%units_t0%sec, self%units_t0%min, self%units_t0%sec
      end if
    else
      write (gol,'("neither `units` or (`units_step`,`units_t0`) defined")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! * set calender:
    
    ! not defined yet ? then default:
    if ( len_trim(self%calendar) == 0 ) self%calendar = 'standard'
    
    ! * climatology
    
    ! (re)set from argument ?
    if ( present(climatology) ) self%climatology = self%climatology .or. climatology
    
    ! apply ?
    if ( self%climatology ) then
      ! name not defined yet ?
      if ( len_trim(self%vname_climat) == 0 ) then
        ! default:
        self%vname_climat = trim(self%name)//'_climatology_bounds'
      end if
    end if
    
    ! *
    
    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Set_Values


  ! ***
  
  
  subroutine TimeCoordinate_Set_Value( self, irec, status, &
                                         t, climat )
  
    use GO, only : TDate, IsAnyDate, operator(-), rTotal
    use GO, only : TIncrDate
    use GO, only : Extract_Ref_and_Step
    
    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(inout)           ::  self
    integer, intent(in)                            ::  irec
    integer, intent(out)                           ::  status

    type(TDate), intent(in), optional              ::  t
    type(TDate), intent(in), optional              ::  climat(nv)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Set_Value'
    
    ! --- local ----------------------------------
    
    integer             ::  iv
    real, allocatable   ::  tmp_values(:)
    real, allocatable   ::  tmp_climat(:,:)
    
    ! --- begin ----------------------------------
    
    ! * check record index:
    
    ! switch:
    if ( self%unlimited ) then
      ! check ...
      if ( irec < 1 ) then
        write (gol,'("argument `irec` out of range 1 .. ")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! reset:
      self%n = max( self%n, irec )
    else
      ! check ...
      if ( (irec < 1) .or. (irec > self%n) ) then
        write (gol,'("argument `irec` out of range 1 .. ",i0)') self%n; call goErr
        TRACEBACK; status=1; return
      end if
    end if
    ! check ...
    if ( self%n < 1 ) then
      write (gol,'("number of time records still undefined")'); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! * time
    
    ! provided as argument ?
    if ( present(t) ) then
      ! initialize storage if necessary:
      if ( .not. allocated(self%values) ) then
        ! new:
        allocate( self%values(self%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else if ( size(self%values) < self%n ) then
        ! storage for copy:
        allocate( tmp_values(size(self%values)), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! copy:
        tmp_values = self%values
        ! clear existing:
        deallocate( self%values, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! new storage:
        allocate( self%values(self%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! copy:
        self%values(1:size(tmp_values)) = tmp_values
        ! clear copy:
        deallocate( tmp_values, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! check ...
      if ( IsAnyDate(self%units_t0) ) then
        write (gol,'("units time offset undefined")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( len_trim(self%units_step) == 0 ) then
        write (gol,'("units step undefined")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! convert from structure to real value:
      self%values(irec) = rTotal( t - self%units_t0, trim(self%units_step) )
    end if
    
    ! * climatology
    
    ! (re)set flag:
    self%climatology = self%climatology .or. present(climat)
    
    ! set name if not done yet:
    if ( len_trim(self%vname_climat) == 0 ) then
      ! add extension:
      self%vname_climat = trim(self%name)//'_climatology_bounds'
    end if
    
    ! provided as argument ?
    if ( present(climat) ) then
      ! initialize storage if necessary:
      if ( .not. allocated(self%climat) ) then
        ! new storage:
        allocate( self%climat(nv,self%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else if ( size(self%climat,2) < self%n ) then
        ! storage for copy:
        allocate( tmp_climat(nv,size(self%climat,2)), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! copy:
        tmp_climat = self%climat
        ! clear existing:
        deallocate( self%climat, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! new storage:
        allocate( self%climat(nv,self%n), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! copy:
        self%climat(:,1:size(tmp_climat,2)) = tmp_climat
        ! clear copy:
        deallocate( tmp_climat, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      ! check ...
      if ( IsAnyDate(self%units_t0) ) then
        write (gol,'("units time offset undefined")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! check ...
      if ( len_trim(self%units_step) == 0 ) then
        write (gol,'("units step undefined")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! loop over boundaries:
      do iv = 1, nv
        ! convert from structure to real value:
        self%climat(iv,irec) = rTotal( climat(iv) - self%units_t0, trim(self%units_step) )
      end do
    end if
    
    ! *
    
    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Set_Value


  ! ***
  

  subroutine TimeCoordinate_Get_Value( self, irec, status, &
                                          t, climat )
  
    use GO, only : TDate, Num_to_Date

    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(in)             ::  self
    integer, intent(in)                           ::  irec
    integer, intent(out)                          ::  status

    type(TDate), intent(out), optional            ::  t
    type(TDate), intent(out), optional            ::  climat(nv)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Get_Value'
    
    ! --- local ----------------------------------
    
    integer     ::  nt
    integer     ::  iv
    
    ! --- begin ----------------------------------
    
    ! get size:
    call self%Coordinate%Get_Dim( status, n=nt )
    IF_NOT_OK_RETURN(status=1)
    
    ! check ...
    if ( (irec < 1) .or. (irec > nt) ) then
      write (gol,'("argument `irec` out of range 1 .. ",i0)') nt; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! arguments that need units ?
    if ( present(t) .or. present(climat) ) then
      if ( len_trim(self%units) == 0 ) then
        write (gol,'("units not defined yet")'); call goErr
        TRACEBACK; status=1; return
      end if
    end if
    
    ! return time value ?
    if ( present(t) ) then
      ! check ...
      if ( .not. allocated(self%values) ) then
        write (gol,'("time values not stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! convert:
      call Num_to_Date( self%values(irec), self%units, t, status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! return climatology bounds ?
    if ( present(climat) ) then
      ! check ...
      if ( .not. allocated(self%climat) ) then
        write (gol,'("climat values not stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! loop over bounds:
      do iv = 1, nv
        ! convert:
        call Num_to_Date( self%climat(iv,irec), self%units, climat(iv), status )
        IF_NOT_OK_RETURN(status=1)
      end do
    end if
    
    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Get_Value


  ! ***


  subroutine TimeCoordinate_Def( self, file, status )
  
    use NetCDF, only : NF90_FLOAT
    use NetCDF, only : NF90_Def_Dim, NF90_Inq_DimID
    use NetCDF, only : NF90_Def_Var  
    use NetCDF, only : NF90_Put_Att
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(inout)    ::  self
    class(Datafile), intent(in)         ::  file
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Def'
    
    ! --- local ----------------------------------

    character(len=LEN_NAME)   ::  name
    integer                   ::  dimid_time
    integer                   ::  dimid_nv
    integer                   ::  nt
    
    ! --- begin ----------------------------------
    
    ! define dimension:
    call self%Dimension%Def( file, status )
    IF_NOT_OK_RETURN(status=1)
    ! get access id:
    call self%Dimension%Get_Dim( status, name=name, n=nt, dimid=dimid_time )
    IF_NOT_OK_RETURN(status=1)
    
    ! define new variable:
    status = NF90_Def_Var( file%ncid, trim(name), NF90_FLOAT, (/dimid_time/), self%varid_out )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! define default attributes:
    call self%Coordinate_Def( file, self%varid_out, status )
    IF_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Put_Att( file%ncid, self%varid_out, 'calendar', trim(self%calendar) )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! climatology bounds defined ?
    if ( self%climatology ) then

      ! check ...
      if ( len_trim(self%vname_climat) == 0 ) then
        write (gol,'("vname_climat not defind")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! attribute:
      status = NF90_Put_Att( file%ncid, self%varid_out, 'climatology', trim(self%vname_climat) )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! create dimension if not done yet:
      status = NF90_INQ_DimID( file%ncid, dname_nv, dimid_nv )
      if (status/=NF90_NOERR) then
        status = NF90_Def_Dim( file%ncid, dname_nv, nv, dimid_nv )
        IF_NF90_NOT_OK_RETURN(status=1)
      end if

      ! define new variable:
      status = NF90_Def_Var( file%ncid, trim(self%vname_climat), NF90_FLOAT, (/dimid_nv,dimid_time/), self%varid_out_climat )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! define same default attributes:
      call self%Coordinate_Def( file, self%varid_out_climat, status )
      IF_NOT_OK_RETURN(status=1)
      
    end if

    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Def


  ! ***
  

  subroutine TimeCoordinate_Write( self, file, status )

    use NetCDF, only : NF90_Put_Var
    
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(in)   ::  self
    class(Datafile), intent(in)         ::  file
    integer, intent(out)                ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Write'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! write general coordinate:
    call self%Coordinate_Write( file, status )
    IF_NOT_OK_RETURN(status=1)

    ! check ...
    if ( .not. allocated(self%values) ) then
      write (gol,'("no values stored yet")'); call goErr
      TRACEBACK; status=1; return
    end if
    ! write:
    status = NF90_Put_Var( file%ncid, self%varid_out, self%values )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! climatology boundaries?
    if ( self%climatology ) then
    
      ! check ...
      if ( .not. allocated(self%climat) ) then
        write (gol,'("no `climat` values stored yet")'); call goErr
        TRACEBACK; status=1; return
      end if
      ! write:
      status = NF90_Put_Var( file%ncid, self%varid_out_climat, self%climat )
      IF_NF90_NOT_OK_RETURN(status=1)
    
    end if  ! bounds

    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Write


  ! ***
  

  subroutine TimeCoordinate_Write_Value( self, file, irec, status, &
                                           t, climat )

    use NetCDF, only : NF90_Put_Var
    
    use GO           , only : TDate
    use C3PO_Datafile, only : Datafile
  
    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(inout)  ::  self
    class(Datafile), intent(in)           ::  file
    integer, intent(in)                   ::  irec
    integer, intent(out)                  ::  status
    
    type(TDate), intent(in), optional     ::  t
    type(TDate), intent(in), optional     ::  climat(nv)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Write_Value'

    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    call self%Set_Value( irec, status, t=t, climat=climat )
    IF_NOT_OK_RETURN(status=1)

    ! write element:
    status = NF90_Put_Var( file%ncid, self%varid_out, self%values(irec:irec), &
                             start=(/irec/), count=(/1/) )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! climatology boundaries?
    if ( self%climatology ) then
      ! write elements:
      status = NF90_Put_Var( file%ncid, self%varid_out_climat, self%climat(:,irec), &
                               start=(/1,irec/), count=(/nv,1/) )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if  ! bounds

    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Write_Value


  ! ***
  

  subroutine TimeCoordinate_Read( self, file, status )

    use NetCDF, only : NF90_INQ_VarID, NF90_Inquire_Variable, NF90_Get_Var
    use NetCDF, only : NF90_Inq_AttName, NF90_Get_Att
  
    use GO           , only : goReadFromLine   
    use C3PO_Datafile, only : Datafile

    ! --- in/out ---------------------------------
    
    XTYPE(TimeCoordinate), intent(inout)    ::  self
    class(Datafile), intent(in)             ::  file
    integer, intent(out)                    ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/TimeCoordinate_Read'
    
    ! --- local ----------------------------------

    character(len=LEN_NAME)   ::  name
    integer                   ::  nt
    integer                   ::  varid
    integer                   ::  natt, iatt
    character(len=LEN_NAME)   ::  attname
    
    ! --- begin ----------------------------------
    
    ! read dimension parameters (length):
    call self%Dimension%Read( file, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! obtain name:
    call self%Dimension%Get_Dim( status, name=name, n=nt )
    IF_NOT_OK_RETURN(status=1)

    ! get variable id:
    status = NF90_INQ_VarID( file%ncid, trim(name), varid )
    if ( status /= NF90_NOERR) then
      gol=NF90_StrError(status); call goErr
      write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
              trim(name), trim(file%filename); call goErr
      TRACEBACK; status=1; return
    end if
    ! attribute:
    status = NF90_Get_Att( file%ncid, varid, 'units', self%units )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Get_Att( file%ncid, varid, 'standard_name', self%standard_name )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Get_Att( file%ncid, varid, 'calendar', self%calendar )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! storage:
    allocate( self%values(nt), stat=status )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! read:
    status = NF90_Get_Var( file%ncid, varid, self%values )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! check if climatology bounds are present ...
    self%climatology = .false.
    ! count number of attributes:
    status = NF90_Inquire_Variable( file%ncid, varid, nAtts=natt )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! loop over attributes:
    do iatt = 1, natt
      ! name:
      status = NF90_Inq_AttName( file%ncid, varid, iatt, attname )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! compare:
      self%climatology = trim(attname) == 'climatology'
      ! found?
      if ( self%climatology ) then
        ! get name of variable with climatology bounds:
        status = NF90_Get_Att( file%ncid, varid, trim(attname), self%vname_climat )
        IF_NF90_NOT_OK_RETURN(status=1)
        ! found, leave now:
        exit
      end if
    end do
    
    ! climatology bounds ?
    if ( self%climatology ) then

      ! storage:
      allocate( self%climat(nv,nt), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! get variable id:
      status = NF90_INQ_VarID( file%ncid, trim(self%vname_climat), varid )
      if ( status /= NF90_NOERR) then
        gol=NF90_StrError(status); call goErr
        write (gol,'("could not inquire id for variable `",a,"` in file `",a,"`")') &
                trim(self%vname_climat), trim(file%filename); call goErr
        TRACEBACK; status=1; return
      end if
      ! read:
      status = NF90_Get_Var( file%ncid, varid, self%climat )
      IF_NF90_NOT_OK_RETURN(status=1)
    
    end if  ! bounds

    ! ok
    status = 0
    
  end subroutine TimeCoordinate_Read


end module C3PO_Coordinates

