!######################################################################
!
! C3PO - CF Convention Compliance Program Objects
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

module C3PO_Datafile

  use DA_Util_ml, only : gol, goPr, goErr
  use NetCDF, only : NF90_NOERR, NF90_StrError

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  Datafile
  

  ! --- const ------------------------------------

  character(len=*), parameter  ::  mname = 'C3PO_Datafile'

  ! lengths:
  integer, parameter    ::  LEN_FILENAME = 1024


  ! --- types ----------------------------------------
  
  type Datafile
    ! filename:
    character(len=LEN_FILENAME)   ::  filename
    ! handle:
    integer                       ::  ncid
    ! flag:
    logical                       ::  parallel
  contains
    procedure   ::  Open          => C3PO_Datafile_Open
    procedure   ::  Create        => C3PO_Datafile_Create
    procedure   ::  Close         => C3PO_Datafile_Close
    procedure   ::  EndDef        => C3PO_Datafile_EndDef
    !
    procedure   ::                   C3PO_Datafile_Get_Var_i_2d
    procedure   ::                   C3PO_Datafile_Get_Var_r_3d
    procedure   ::                   C3PO_Datafile_Get_Var_r_4d
    procedure   ::                   C3PO_Datafile_Get_Var_r_5d
    generic     ::  Get_Var       => C3PO_Datafile_Get_Var_i_2d, &
                                     C3PO_Datafile_Get_Var_r_3d, &
                                     C3PO_Datafile_Get_Var_r_4d, &
                                     C3PO_Datafile_Get_Var_r_5d
    procedure   ::  Get_Var_Attr  => C3PO_Datafile_Get_Var_Attr_s
  end type Datafile


contains


  ! ********************************************************************
  ! ***
  ! *** data file
  ! ***
  ! ********************************************************************


  subroutine C3PO_Datafile_Open( self, filename, status, &
                                   writable, comm )
  
#ifdef _MPI
    use MPI   , only : MPI_INFO_NULL
#endif
    use NetCDF, only : NF90_Open
    use NetCDF, only : NF90_WRITE, NF90_NOWRITE
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(out)              ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status
    
    logical, intent(in), optional             ::  writable
    integer, intent(in), optional             ::  comm

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_Open'
    
    ! --- local ----------------------------------
    
    integer                 ::  cmode
    logical                 ::  exist
    
    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! check ...
    inquire( file=trim(filename), exist=exist )
    if ( .not. exist ) then
      write (gol,'("file not found : ",a)') trim(self%filename); call goErr
      TRACEBACK; status=1; return
    end if

    ! creation mode:
    cmode = 0
    if ( present(writable) ) then
      if ( writable ) then
        cmode = cmode + NF90_WRITE         ! read and write
      else
        cmode = cmode + NF90_NOWRITE       ! read-only
      end if
    else
      cmode = cmode + NF90_NOWRITE       ! default read-only
    end if
    
    ! set flag:
    self%parallel = present(comm)

    ! parallel access?
    if ( self%parallel ) then
      ! open new file for parallel access:
#ifdef _MPI
      status = NF90_Open( trim(self%filename), cmode, self%ncid, &
                            comm=comm, info= MPI_INFO_NULL )
      IF_NF90_NOT_OK_RETURN(status=1)
#else
      write (gol,'("Not compiled with _MPI defined, could not open in parallel.")'); call goErr
      TRACEBACK; status=1; return
#endif
    else
      ! open new file:
      status = NF90_Open( trim(self%filename), cmode, self%ncid )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
  
  end subroutine C3PO_Datafile_Open


  ! ***


  subroutine C3PO_Datafile_Create( self, filename, status, comm )
  
#ifdef _MPI
    use MPI   , only : MPI_INFO_NULL
#endif
    use NetCDF, only : NF90_Create
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(out)              ::  self
    character(len=*), intent(in)              ::  filename
    integer, intent(out)                      ::  status
    
    integer, intent(in), optional             ::  comm

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_Create'
    
    ! --- local ----------------------------------
    
    integer                 ::  cmode
    
    ! --- begin ----------------------------------
    
    ! store:
    self%filename = trim(filename)
    
    ! creation mode:
    cmode = 0
    
    ! set flag:
    self%parallel = present(comm)

    ! parallel access?
    if ( self%parallel ) then
      ! create new file for parallel access:
#ifdef _MPI
      status = NF90_Create( trim(self%filename), cmode, self%ncid, &
                            comm=comm, info=MPI_INFO_NULL )
      IF_NF90_NOT_OK_RETURN(status=1)
#else
      write (gol,'("Not compiled with _MPI defined, could not create in parallel.")'); call goErr
      TRACEBACK; status=1; return
#endif
    else
      ! create new file:
      status = NF90_Create( trim(self%filename), cmode, self%ncid )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    
    ! ok
    status = 0
  
  end subroutine C3PO_Datafile_Create


  ! ***
  

  subroutine C3PO_Datafile_Close( self, status )

    use NetCDF, only : NF90_Close
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(inout)              ::  self
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_Close'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! close file:
    status = NF90_Close( self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine C3PO_Datafile_Close


  ! ***
  

  subroutine C3PO_Datafile_EndDef( self, status )

    use NetCDF, only : NF90_EndDef
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(inout)              ::  self
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_EndDef'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! end of definition phase:
    status = NF90_EndDef( self%ncid )
    IF_NF90_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine C3PO_Datafile_EndDef


  ! ********************************************************************
  ! ***
  ! *** get var
  ! ***
  ! ********************************************************************


  subroutine C3PO_Datafile_Get_Var_i_2d( self, name, field, units, status, &
                                          fill_value, start, count )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_Var_Par_Access, NF90_COLLECTIVE
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(inout)            ::  self
    character(len=*), intent(in)              ::  name
    integer, intent(out)                      ::  field(:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    integer, intent(out), optional            ::  fill_value
    integer, intent(in), optional             ::  start(2)
    integer, intent(in), optional             ::  count(2)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_Get_Var_i_2d'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! file opened for parallel acces?
    if ( self%parallel ) then
      ! enable collective read:
      status = NF90_Var_Par_Access( self%ncid, varid, NF90_COLLECTIVE )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, start=start, count=count  )
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
  
  end subroutine C3PO_Datafile_Get_Var_i_2d


  ! ***


  subroutine C3PO_Datafile_Get_Var_r_3d( self, name, field, units, status, &
                                          fill_value, start, count )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_Var_Par_Access, NF90_COLLECTIVE
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(inout)            ::  self
    character(len=*), intent(in)              ::  name
    real, intent(out)                         ::  field(:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value
    integer, intent(in), optional             ::  start(3)
    integer, intent(in), optional             ::  count(3)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_Get_Var_r_3d'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! file opened for parallel acces?
    if ( self%parallel ) then
      ! enable collective read:
      status = NF90_Var_Par_Access( self%ncid, varid, NF90_COLLECTIVE )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, start=start, count=count )
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
  
  end subroutine C3PO_Datafile_Get_Var_r_3d


  ! ***


  subroutine C3PO_Datafile_Get_Var_r_4d( self, name, field, units, status, &
                                          fill_value, start, count )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_Var_Par_Access, NF90_COLLECTIVE
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(inout)            ::  self
    character(len=*), intent(in)              ::  name
    real, intent(out)                         ::  field(:,:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value
    integer, intent(in), optional             ::  start(4)
    integer, intent(in), optional             ::  count(4)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_Get_Var_r_4d'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! file opened for parallel acces?
    if ( self%parallel ) then
      ! enable collective read:
      status = NF90_Var_Par_Access( self%ncid, varid, NF90_COLLECTIVE )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, start=start, count=count )
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
  
  end subroutine C3PO_Datafile_Get_Var_r_4d


  ! ***


  subroutine C3PO_Datafile_Get_Var_r_5d( self, name, field, units, status, &
                                          fill_value, start, count )
  
    use NetCDF, only : NF90_Inq_VarID, NF90_Get_Var
    use NetCDF, only : NF90_Get_Att
    use NetCDF, only : NF90_Var_Par_Access, NF90_COLLECTIVE
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(inout)            ::  self
    character(len=*), intent(in)              ::  name
    real, intent(out)                         ::  field(:,:,:,:,:)
    character(len=*), intent(out)             ::  units
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  fill_value
    integer, intent(in), optional             ::  start(5)
    integer, intent(in), optional             ::  count(5)

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_Get_Var_r_5d'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! file opened for parallel acces?
    if ( self%parallel ) then
      ! enable collective read:
      status = NF90_Var_Par_Access( self%ncid, varid, NF90_COLLECTIVE )
      IF_NF90_NOT_OK_RETURN(status=1)
    end if
    ! read:
    status = NF90_Get_Var( self%ncid, varid, field, start=start, count=count )
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
  
  end subroutine C3PO_Datafile_Get_Var_r_5d


  ! ********************************************************************
  ! ***
  ! *** get attr
  ! ***
  ! ********************************************************************


  subroutine C3PO_Datafile_Get_Var_Attr_s( self, name, attr, values, status )
  
    use NetCDF, only : NF90_Inq_VarID
    use NetCDF, only : NF90_Get_Att
  
    ! --- in/out ---------------------------------
    
    class(Datafile), intent(inout)            ::  self
    character(len=*), intent(in)              ::  name
    character(len=*), intent(in)              ::  attr
    character(len=*), intent(out)             ::  values
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/C3PO_Datafile_Get_Var_Attr_s'
    
    ! --- local ----------------------------------
    
    integer       ::  varid
    
    ! --- begin ----------------------------------
    
    ! access variable:
    status = NF90_Inq_VarID( self%ncid, trim(name), varid )
    IF_NF90_NOT_OK_RETURN(status=1)
    ! attribute:
    status = NF90_Get_Att( self%ncid, varid, attr, values )
    IF_NF90_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine C3PO_Datafile_Get_Var_Attr_s

  
end module C3PO_Datafile

