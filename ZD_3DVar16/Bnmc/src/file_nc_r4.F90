!### macro's #####################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,i6,")")') rname, __FILE__, __LINE__ ; call goErr
!
#define IF_NOTOK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOTOK_RETURN(action) if (status/=NF90_NOERR) then; gol=nf90_strerror(status); call goErr; TRACEBACK; action; return; end if
!
!#################################################################

module file_nc_r4

  use GO, only : gol, goPr, goErr
  
  use NetCDF

  implicit none

 
  ! --- in/out ------------------------------------
  
  private
  
  public  ::  nc_dump


  ! --- const -------------------------------------

  character(len=*), parameter  ::  mname = 'file_nc_r4'
  

  ! --- interfaces --------------------------------
  
  interface nc_dump
    module procedure nc_dump_r4_1d
    module procedure nc_dump_r4_2d
    module procedure nc_dump_r4_3d
    module procedure nc_dump_r4_4d
    module procedure nc_dump_r4_5d
    module procedure nc_dump_r4_6d
    module procedure nc_dump_r4_7d
  end interface
  

contains


  ! =============================================================
  
  
  subroutine nc_dump_r4_1d( file_name, x, var_name, dim_names, status, typ )

    use NetCDF
      
    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/nc_dump_r4_1d'
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  file_name
    real(4), intent(in)    ::  x(:)
    character(len=*), intent(in)      ::  var_name
    character(len=*), intent(in)      ::  dim_names(1)
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  typ
    
    ! --- local ----------------------------
    
    integer             ::  var_type
    character(len=4)    ::  transfer_func
    integer(4)          ::  transfer_kind
    integer             ::  ncid
    integer             ::  dimids(1)
    integer             ::  varid
    integer             ::  idim
    
    ! --- begin ----------------------------
    
    ! set default target type and transfer function:
    var_type = NF90_REAL  
    transfer_func = 'none'
    
    ! overwrite for special requests:
    if ( present(typ) ) then
      select case ( typ )
        case ( 'int', 'integer', 'integer(4)' )
          var_type = NF90_INT
          transfer_func = 'int '
          transfer_kind = 4
        case ( 'real', 'float', 'real(4)' )
          var_type = NF90_REAL
          transfer_func = 'real'
          transfer_kind = 4
        case ( 'double', 'real(8)' )
          var_type = NF90_DOUBLE
          transfer_func = 'real'
          transfer_kind = 8
        case default
          write (gol,'("Unsupported type : ",a)') trim(typ); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    
    ! open file, overwrite if necessary:
    status = NF90_Create( trim(file_name), NF90_CLOBBER, ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, 1
      ! create dimension:
      status = NF90_Def_Dim( ncid, dim_names(idim), size(x,idim), dimids(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do
    
    ! create variable:
    status = NF90_Def_Var( ncid, var_name, var_type, dimids, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! end definition phase:
    status = NF90_EndDef( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    select case ( transfer_func )
      case ( 'none' )
        status = NF90_Put_Var( ncid, varid, x )
        IF_NF90_NOTOK_RETURN(status=1)
      case ( 'int ' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, int(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case ( 'real' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 8 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=8) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case default
        write (gol,'("unsupported transfer_func : ",a)') trim(transfer_func); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! close:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine nc_dump_r4_1d
  
  
  ! ***
  
  
  subroutine nc_dump_r4_2d( file_name, x, var_name, dim_names, status, typ )

    use NetCDF
      
    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/nc_dump_r4_2d'
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  file_name
    real(4), intent(in)    ::  x(:,:)
    character(len=*), intent(in)      ::  var_name
    character(len=*), intent(in)      ::  dim_names(2)
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  typ
    
    ! --- local ----------------------------
    
    integer             ::  var_type
    character(len=4)    ::  transfer_func
    integer(4)          ::  transfer_kind
    integer             ::  ncid
    integer             ::  dimids(2)
    integer             ::  varid
    integer             ::  idim
    
    ! --- begin ----------------------------
    
    ! set default target type and transfer function:
    var_type = NF90_REAL  
    transfer_func = 'none'
    
    ! overwrite for special requests:
    if ( present(typ) ) then
      select case ( typ )
        case ( 'int', 'integer', 'integer(4)' )
          var_type = NF90_INT
          transfer_func = 'int '
          transfer_kind = 4
        case ( 'real', 'float', 'real(4)' )
          var_type = NF90_REAL
          transfer_func = 'real'
          transfer_kind = 4
        case ( 'double', 'real(8)' )
          var_type = NF90_DOUBLE
          transfer_func = 'real'
          transfer_kind = 8
        case default
          write (gol,'("Unsupported type : ",a)') trim(typ); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    
    ! open file, overwrite if necessary:
    status = NF90_Create( trim(file_name), NF90_CLOBBER, ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, 2
      ! create dimension:
      status = NF90_Def_Dim( ncid, dim_names(idim), size(x,idim), dimids(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do
    
    ! create variable:
    status = NF90_Def_Var( ncid, var_name, var_type, dimids, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! end definition phase:
    status = NF90_EndDef( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    select case ( transfer_func )
      case ( 'none' )
        status = NF90_Put_Var( ncid, varid, x )
        IF_NF90_NOTOK_RETURN(status=1)
      case ( 'int ' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, int(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case ( 'real' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 8 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=8) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case default
        write (gol,'("unsupported transfer_func : ",a)') trim(transfer_func); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! close:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine nc_dump_r4_2d
  
  
  ! ***
  
  
  subroutine nc_dump_r4_3d( file_name, x, var_name, dim_names, status, typ )

    use NetCDF
      
    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/nc_dump_r4_3d'
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  file_name
    real(4), intent(in)    ::  x(:,:,:)
    character(len=*), intent(in)      ::  var_name
    character(len=*), intent(in)      ::  dim_names(3)
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  typ
    
    ! --- local ----------------------------
    
    integer             ::  var_type
    character(len=4)    ::  transfer_func
    integer(4)          ::  transfer_kind
    integer             ::  ncid
    integer             ::  dimids(3)
    integer             ::  varid
    integer             ::  idim
    
    ! --- begin ----------------------------
    
    ! set default target type and transfer function:
    var_type = NF90_REAL  
    transfer_func = 'none'
    
    ! overwrite for special requests:
    if ( present(typ) ) then
      select case ( typ )
        case ( 'int', 'integer', 'integer(4)' )
          var_type = NF90_INT
          transfer_func = 'int '
          transfer_kind = 4
        case ( 'real', 'float', 'real(4)' )
          var_type = NF90_REAL
          transfer_func = 'real'
          transfer_kind = 4
        case ( 'double', 'real(8)' )
          var_type = NF90_DOUBLE
          transfer_func = 'real'
          transfer_kind = 8
        case default
          write (gol,'("Unsupported type : ",a)') trim(typ); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    
    ! open file, overwrite if necessary:
    status = NF90_Create( trim(file_name), NF90_CLOBBER, ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, 3
      ! create dimension:
      status = NF90_Def_Dim( ncid, dim_names(idim), size(x,idim), dimids(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do
    
    ! create variable:
    status = NF90_Def_Var( ncid, var_name, var_type, dimids, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! end definition phase:
    status = NF90_EndDef( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    select case ( transfer_func )
      case ( 'none' )
        status = NF90_Put_Var( ncid, varid, x )
        IF_NF90_NOTOK_RETURN(status=1)
      case ( 'int ' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, int(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case ( 'real' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 8 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=8) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case default
        write (gol,'("unsupported transfer_func : ",a)') trim(transfer_func); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! close:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine nc_dump_r4_3d
  
  
  ! ***
  
  
  subroutine nc_dump_r4_4d( file_name, x, var_name, dim_names, status, typ )

    use NetCDF
      
    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/nc_dump_r4_4d'
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  file_name
    real(4), intent(in)    ::  x(:,:,:,:)
    character(len=*), intent(in)      ::  var_name
    character(len=*), intent(in)      ::  dim_names(4)
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  typ
    
    ! --- local ----------------------------
    
    integer             ::  var_type
    character(len=4)    ::  transfer_func
    integer(4)          ::  transfer_kind
    integer             ::  ncid
    integer             ::  dimids(4)
    integer             ::  varid
    integer             ::  idim
    
    ! --- begin ----------------------------
    
    ! set default target type and transfer function:
    var_type = NF90_REAL  
    transfer_func = 'none'
    
    ! overwrite for special requests:
    if ( present(typ) ) then
      select case ( typ )
        case ( 'int', 'integer', 'integer(4)' )
          var_type = NF90_INT
          transfer_func = 'int '
          transfer_kind = 4
        case ( 'real', 'float', 'real(4)' )
          var_type = NF90_REAL
          transfer_func = 'real'
          transfer_kind = 4
        case ( 'double', 'real(8)' )
          var_type = NF90_DOUBLE
          transfer_func = 'real'
          transfer_kind = 8
        case default
          write (gol,'("Unsupported type : ",a)') trim(typ); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    
    ! open file, overwrite if necessary:
    status = NF90_Create( trim(file_name), NF90_CLOBBER, ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, 4
      ! create dimension:
      status = NF90_Def_Dim( ncid, dim_names(idim), size(x,idim), dimids(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do
    
    ! create variable:
    status = NF90_Def_Var( ncid, var_name, var_type, dimids, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! end definition phase:
    status = NF90_EndDef( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    select case ( transfer_func )
      case ( 'none' )
        status = NF90_Put_Var( ncid, varid, x )
        IF_NF90_NOTOK_RETURN(status=1)
      case ( 'int ' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, int(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case ( 'real' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 8 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=8) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case default
        write (gol,'("unsupported transfer_func : ",a)') trim(transfer_func); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! close:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine nc_dump_r4_4d
  
  
  ! ***
  
  
  subroutine nc_dump_r4_5d( file_name, x, var_name, dim_names, status, typ )

    use NetCDF
      
    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/nc_dump_r4_5d'
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  file_name
    real(4), intent(in)    ::  x(:,:,:,:,:)
    character(len=*), intent(in)      ::  var_name
    character(len=*), intent(in)      ::  dim_names(5)
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  typ
    
    ! --- local ----------------------------
    
    integer             ::  var_type
    character(len=4)    ::  transfer_func
    integer(4)          ::  transfer_kind
    integer             ::  ncid
    integer             ::  dimids(5)
    integer             ::  varid
    integer             ::  idim
    
    ! --- begin ----------------------------
    
    ! set default target type and transfer function:
    var_type = NF90_REAL  
    transfer_func = 'none'
    
    ! overwrite for special requests:
    if ( present(typ) ) then
      select case ( typ )
        case ( 'int', 'integer', 'integer(4)' )
          var_type = NF90_INT
          transfer_func = 'int '
          transfer_kind = 4
        case ( 'real', 'float', 'real(4)' )
          var_type = NF90_REAL
          transfer_func = 'real'
          transfer_kind = 4
        case ( 'double', 'real(8)' )
          var_type = NF90_DOUBLE
          transfer_func = 'real'
          transfer_kind = 8
        case default
          write (gol,'("Unsupported type : ",a)') trim(typ); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    
    ! open file, overwrite if necessary:
    status = NF90_Create( trim(file_name), NF90_CLOBBER, ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, 5
      ! create dimension:
      status = NF90_Def_Dim( ncid, dim_names(idim), size(x,idim), dimids(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do
    
    ! create variable:
    status = NF90_Def_Var( ncid, var_name, var_type, dimids, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! end definition phase:
    status = NF90_EndDef( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    select case ( transfer_func )
      case ( 'none' )
        status = NF90_Put_Var( ncid, varid, x )
        IF_NF90_NOTOK_RETURN(status=1)
      case ( 'int ' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, int(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case ( 'real' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 8 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=8) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case default
        write (gol,'("unsupported transfer_func : ",a)') trim(transfer_func); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! close:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine nc_dump_r4_5d
  
  
  ! ***
  
  
  subroutine nc_dump_r4_6d( file_name, x, var_name, dim_names, status, typ )

    use NetCDF
      
    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/nc_dump_r4_6d'
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  file_name
    real(4), intent(in)    ::  x(:,:,:,:,:,:)
    character(len=*), intent(in)      ::  var_name
    character(len=*), intent(in)      ::  dim_names(6)
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  typ
    
    ! --- local ----------------------------
    
    integer             ::  var_type
    character(len=4)    ::  transfer_func
    integer(4)          ::  transfer_kind
    integer             ::  ncid
    integer             ::  dimids(6)
    integer             ::  varid
    integer             ::  idim
    
    ! --- begin ----------------------------
    
    ! set default target type and transfer function:
    var_type = NF90_REAL  
    transfer_func = 'none'
    
    ! overwrite for special requests:
    if ( present(typ) ) then
      select case ( typ )
        case ( 'int', 'integer', 'integer(4)' )
          var_type = NF90_INT
          transfer_func = 'int '
          transfer_kind = 4
        case ( 'real', 'float', 'real(4)' )
          var_type = NF90_REAL
          transfer_func = 'real'
          transfer_kind = 4
        case ( 'double', 'real(8)' )
          var_type = NF90_DOUBLE
          transfer_func = 'real'
          transfer_kind = 8
        case default
          write (gol,'("Unsupported type : ",a)') trim(typ); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    
    ! open file, overwrite if necessary:
    status = NF90_Create( trim(file_name), NF90_CLOBBER, ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, 6
      ! create dimension:
      status = NF90_Def_Dim( ncid, dim_names(idim), size(x,idim), dimids(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do
    
    ! create variable:
    status = NF90_Def_Var( ncid, var_name, var_type, dimids, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! end definition phase:
    status = NF90_EndDef( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    select case ( transfer_func )
      case ( 'none' )
        status = NF90_Put_Var( ncid, varid, x )
        IF_NF90_NOTOK_RETURN(status=1)
      case ( 'int ' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, int(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case ( 'real' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 8 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=8) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case default
        write (gol,'("unsupported transfer_func : ",a)') trim(transfer_func); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! close:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine nc_dump_r4_6d
  
  
  ! ***
  
  
  subroutine nc_dump_r4_7d( file_name, x, var_name, dim_names, status, typ )

    use NetCDF
      
    ! --- const ---------------------------
    
    character(len=*), parameter   ::  rname = mname//'/nc_dump_r4_7d'
    
    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)      ::  file_name
    real(4), intent(in)    ::  x(:,:,:,:,:,:,:)
    character(len=*), intent(in)      ::  var_name
    character(len=*), intent(in)      ::  dim_names(7)
    integer, intent(out)              ::  status
    
    character(len=*), intent(in), optional    ::  typ
    
    ! --- local ----------------------------
    
    integer             ::  var_type
    character(len=4)    ::  transfer_func
    integer(4)          ::  transfer_kind
    integer             ::  ncid
    integer             ::  dimids(7)
    integer             ::  varid
    integer             ::  idim
    
    ! --- begin ----------------------------
    
    ! set default target type and transfer function:
    var_type = NF90_REAL  
    transfer_func = 'none'
    
    ! overwrite for special requests:
    if ( present(typ) ) then
      select case ( typ )
        case ( 'int', 'integer', 'integer(4)' )
          var_type = NF90_INT
          transfer_func = 'int '
          transfer_kind = 4
        case ( 'real', 'float', 'real(4)' )
          var_type = NF90_REAL
          transfer_func = 'real'
          transfer_kind = 4
        case ( 'double', 'real(8)' )
          var_type = NF90_DOUBLE
          transfer_func = 'real'
          transfer_kind = 8
        case default
          write (gol,'("Unsupported type : ",a)') trim(typ); call goErr
          TRACEBACK; status=1; return
      end select
    end if
    
    ! open file, overwrite if necessary:
    status = NF90_Create( trim(file_name), NF90_CLOBBER, ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! loop over dimensions:
    do idim = 1, 7
      ! create dimension:
      status = NF90_Def_Dim( ncid, dim_names(idim), size(x,idim), dimids(idim) )
      IF_NF90_NOTOK_RETURN(status=1)
    end do
    
    ! create variable:
    status = NF90_Def_Var( ncid, var_name, var_type, dimids, varid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! end definition phase:
    status = NF90_EndDef( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! write data:
    select case ( transfer_func )
      case ( 'none' )
        status = NF90_Put_Var( ncid, varid, x )
        IF_NF90_NOTOK_RETURN(status=1)
      case ( 'int ' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, int(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case ( 'real' )
        select case ( transfer_kind )
          case ( 4 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=4) )
            IF_NF90_NOTOK_RETURN(status=1)
          case ( 8 )
            status = NF90_Put_Var( ncid, varid, real(x,kind=8) )
            IF_NF90_NOTOK_RETURN(status=1)
          case default
            write (gol,'("unsupported kind ",i2," for ",a)') transfer_kind, trim(transfer_func); call goErr
            TRACEBACK; status=1; return
        end select
      case default
        write (gol,'("unsupported transfer_func : ",a)') trim(transfer_func); call goErr
        TRACEBACK; status=1; return
    end select
    
    ! close:
    status = NF90_Close( ncid )
    IF_NF90_NOTOK_RETURN(status=1)
    
    ! ok
    status = 0
  
  end subroutine nc_dump_r4_7d
  
  
  ! ***
  
  
end module file_nc_r4
