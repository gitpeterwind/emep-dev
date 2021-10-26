!*****************************************************************************!
!
! EMEP interface to CSO (CAMS Satellite Operator) modules.
! Used to simulate satellite retrievals from model.
!
! HISTORY
!   2021-09, Arjo Segers
!     Extracted content from 3D-var analysis codes to allow
!     simulation from standard model.
!
!*****************************************************************************!
!
#define TRACELINE rname//' (__FILE__, line __LINE__)'
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
!
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action)  if (status> 0) then; TRACEBACK; action; return; end if
!
!*****************************************************************************!

module EMEP_CSO_mod

  use CSO             , only : T_CSO
  use CSO             , only : T_CSO_RcFile
  use CSO             , only : T_CSO_Listing
  use CSO             , only : T_CSO_Sat_Data
  use CSO             , only : T_CSO_Sat_State

  implicit none


  ! --- in/out -----------------------------------------

  private

  public  ::  NTIMING_CSO
  public  ::  TIMING_CSO

  public  ::  EMEP_CSO_Init, EMEP_CSO_Done
  public  ::  EMEP_CSO_Simulate
  public  ::  EMEP_CSO_Clear


  ! --- const -----------------------------------------

  character(len=*), parameter  ::  mname = 'EMEP_CSO_mod'

  ! timing parameters:
  !~ number of timing parameters:
  integer, parameter  ::  NTIMING_CSO = 0
  !~ timer index in 1,..,NTIMING_CSO  
  !  (outside to not collect timing)
  integer, parameter  ::  TIMING_CSO  = 0

  ! maximum number of sets:
  integer, parameter                    ::  maxcso = 10


  ! --- types -----------------------------------------
  
  ! storage for entities used for mapping
  ! from EMEP grid to pixel footprints:
  type :: T_EMEP_GridMapping
    ! grid description:
    integer                   ::  nlon, nlat
    ! mapping type:
    integer                   ::  levels
    ! sampling points:
    real, allocatable         ::  xxp(:), yyp(:)
    real, allocatable         ::  wwp(:)
    ! source cells, weights:
    integer, pointer          ::  ii(:), jj(:)
    real, pointer             ::  ww(:)
    ! work array:
    real, allocatable         ::  wmap(:,:)  ! (nlon,nlat)
    !
    ! idem for pixel arrays:
    integer                   ::  nall
    integer                   ::  mxall
    integer, pointer          ::  all_iw0(:), all_nw(:)   ! (npix)
    integer, pointer          ::  all_ii(:), all_jj(:)   ! (mxall)
    real, pointer             ::  all_ww(:)    ! (mxall)
    real, pointer             ::  all_area(:)  ! (npix)
    !
  contains
    procedure ::  Init         =>  EMEP_GridMapping_Init
    procedure ::  Done         =>  EMEP_GridMapping_Done
    procedure ::                   EMEP_GridMapping_GetWeights_0d
    procedure ::                   EMEP_GridMapping_GetWeights_1d
    generic   ::  GetWeights   =>  EMEP_GridMapping_GetWeights_0d, &
                                   EMEP_GridMapping_GetWeights_1d
  end type T_EMEP_GridMapping


  ! --- var -------------------------------------------
  
  ! message line:
  character(len=1024)                   ::  gol

  ! main object for CSO data:  
  type(T_CSO)                           ::  csod
  ! settings:
  type(T_CSO_RcFile)                    ::  cso_rcf
  ! actual number of datasets:
  integer                               ::  ncso
  ! datasets:
  character(len=32)                     ::  cso_keys(maxcso)
  ! listing of satellite fiels:
  type(T_CSO_Listing)                   ::  cso_listing(maxcso)
  ! current open file, or empty:
  character(len=1024)                   ::  cso_orbit_filename(maxcso)
  ! data and state:
  type(T_CSO_Sat_Data)                  ::  cso_sdata(maxcso)
  type(T_CSO_Sat_State)                 ::  cso_sstate(maxcso)
  ! observed tracer:
  character(len=32)                     ::  cso_tracer(maxcso)
  ! tool for mapping regular grid to footprint:
  type(T_EMEP_GridMapping)              ::  emep_grid_mapping
  ! how to decide if pixel overlaps (partly) with domain:
  !  'center'     :  based on center only
  !  'footprint'  :  based on corners
  character(len=16)                     ::  cso_mapping_select_by
  
  
contains


  !=====================================================================
  !===
  !=== logging
  !===
  !=====================================================================
  
  ! write message stored in "gol" to model logging system,
  ! only on root:
  subroutine goPr()
    use Io_RunLog_mod, only : PrintLog
    use Config_module, only : MasterProc
    call PrintLog( trim(gol), OutputProc=MasterProc )
  end subroutine goPr
  
  ! write error message stored in "gol" to model logging system,
  ! all processors:
  subroutine goErr()
    use Io_RunLog_mod, only : PrintLog
    call PrintLog( 'ERROR - '//trim(gol) )
  end subroutine goErr


  !=====================================================================
  !===
  !=== module init/done
  !===
  !=====================================================================


  subroutine EMEP_CSO_Init( status )

    use CSO           , only : CSO_SplitString

    use Io_mod        , only : IO_NML               ! file unit for namelist file
    use MPI_Groups_mod, only : MPI_COMM_CALC
    use Par_Mod       , only : me                   ! processor id: 0,..,nproc-1
    use Par_mod       , only : limax, ljmax         ! local  cell range:  (1:limax,1:ljmax)
    use Par_mod       , only : gi0, gi1, gj0, gj1   ! global cell range:  (gi0:gi1,gj0:gj1)

    ! --- in/out ---------------------------------

    integer, intent(out)           ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/EMEP_CSO_Init'

    ! --- local ----------------------------------

    ! logfile:
    character(len=1024)                   ::  logfile
    ! settings:
    character(len=1024)                   ::  cso_rcfile
    ! storage for keyword values:
    character(len=1024)                   ::  cso_line
    character(len=1024)                   ::  cso_listing_filename
    integer                               ::  cso_mapping_levels
    ! loop variables:
    integer                               ::  icso

    ! --- begin ----------------------------------
    
    ! info ...
    write (gol,'(a,": initialize interface to CSO tools ...")') rname; call goPr
    
    !
    ! Namelist section:
    !
    !    &CSO_CONFIG
    !      ! CSO settings
    !      cso_rcfile                = 'cso-settings.rc',
    !    &end
    !
    ! name list that defines the actual values: 
    namelist /CSO_CONFIG/ cso_rcfile
    ! back to start:
    rewind(IO_NML)
    ! read using namelist:
    read( unit=IO_NML, nml=CSO_CONFIG, iostat=status )
    if ( status /= 0 ) then
      write (gol,'("reading namelist `CSO_CONFIG` from `config_emep.nml`")'); call goErr
      TRACEBACK; status=1; return
    end if

    ! init CSO module including MPI communication:
    call csod%Init( status, icomm=MPI_COMM_CALC )
    IF_NOT_OK_RETURN(status=1)

    ! processor specific logging:
    write (logfile,'("cso-pe",i3.3,".log")') me
    ! redirect CSO messages:
    call csod%SetLogging( status, file=trim(logfile) )
    IF_NOT_OK_RETURN(status=1)
    !! root only?
    !call csod%SetLogging( status, root_only=.true. )
    !IF_NOT_OK_RETURN(status=1)

    ! info ...
    write (gol,'(a,": CSO settings file: ",a)') rname, trim(cso_rcfile); call goPr
    ! read:
    call cso_rcf%Init( cso_rcfile, status )
    IF_NOT_OK_RETURN(status=1)

    ! list with rckey prefixes:
    call cso_rcf%Get( 'cso.keys', cso_line, status )
    IF_NOT_OK_RETURN(status=1)
    ! split:
    call CSO_SplitString( cso_line, ncso, cso_keys, status )
    IF_NOT_OK_RETURN(status=1)
    ! loop over keys:
    do icso = 1, ncso
    
      ! observed tracer:
      call cso_rcf%Get( trim(cso_keys(icso))//'.tracer', cso_tracer(icso), status )
      IF_NOT_OK_RETURN(status=1)

      ! listing file:
      call cso_rcf%Get( trim(cso_keys(icso))//'.listing', cso_listing_filename, status )
      IF_NOT_OK_RETURN(status=1)
      ! info ...
      write (gol,'(a,": read CSO listing file: ",a)') rname, trim(cso_listing_filename); call goPr
      ! read listing file:
      call cso_listing(icso)%Init( cso_listing_filename, status )
      IF_NOT_OK_RETURN(status=1)
      
      ! no data read yet ...
      cso_orbit_filename(icso) = ''
      
    end do ! cso sets

    ! how to decide if pixel overlaps with domain:
    call cso_rcf%Get( 'cso.mapping.select_by', cso_mapping_select_by, status )
    IF_NOT_OK_RETURN(status=1)

    ! mapping level:
    call cso_rcf%Get( 'cso.mapping.levels', cso_mapping_levels, status )
    IF_NOT_OK_RETURN(status=1)

    ! info ...
    write (gol,'(a,": define CSO grid mapping:")') rname; call goPr
    write (gol,'(a,":   global cell range : [",i3,",",i3,"] x [",i3,",",i3,"]")') rname, gi0, gi1, gj0, gj1; call goPr
    write (gol,'(a,":   local domain shape: ",i4," x ",i4)') rname, limax, ljmax; call goPr
    ! init grid mapping using:
    ! - local grid shape
    ! - recursion level
    call emep_grid_mapping%Init( limax, ljmax, cso_mapping_levels, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine EMEP_CSO_Init


  ! ***


  subroutine EMEP_CSO_Done( status )

    ! --- in/out ----------------------------

    integer, intent(out)           ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/EMEP_CSO_Done'

    ! --- local ----------------------------

    integer                               ::  icso

    ! --- begin -----------------------------
    
    ! info ...
    write (gol,'(a,": done with interface to CSO tools ...")') rname; call goPr

    ! done with grid mapping:
    call emep_grid_mapping%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! loop over datasets:
    do icso = 1, ncso
      ! done with listing:
      call cso_listing(icso)%Done( status )
      IF_NOT_OK_RETURN(status=1)
    end do ! icso

    ! done with CSO settings:
    call cso_rcf%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! done with CSO:
    call csod%Done( status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine EMEP_CSO_Done


  ! ***
  
  !
  ! Simulate pixels from model arrays.
  !
  ! Arguments:
  !
  ! * outkey = 'fg' | 'an'
  !     Selection keyword for subset of simulated variables;
  !     here used by assimilation run to distinguish analyzed simulation from first-guess
  ! * putout = .true. | .false.
  !     Flag to enable output. All variables will be saved, thus:
  !     - in standard model run set this to .true.
  !     - in assimilation, set this to .false. for first-guess output and .true. after analysis
  !


  subroutine EMEP_CSO_Simulate( outkey, putout, status )

    use CSO                   , only : T_CSO_DateTime, T_CSO_TimeDelta
    use CSO                   , only : CSO_DateTime, CSO_TimeDelta
    use CSO                   , only : Pretty
    use CSO                   , only : Precisely
    use CSO                   , only : operator(+), operator(-), operator(*), operator(>=)

    use Config_module         , only : KMAX_MID
    use TimeDate_mod          , only : current_date
    use Par_mod               , only : gi0, gj0
    use GridValues_mod        , only : coord_in_domain
    use ChemFields_mod        , only : xn_adv           ! (NSPEC_ADV,LIMAX,LJMAX,KMAX_MID)
    use ChemSpecs_mod         , only : species_adv
    use SmallUtils_mod        , only : find_index
    use Units_mod             , only : Units_Scale
    use GridValues_mod        , only : A_bnd, B_bnd
    use MetFields_mod         , only : ps
    use MetFields_mod         , only : roa

    ! --- in/out ----------------------------

    character(len=*), intent(in)    ::  outkey    ! output selection key: 'fg', 'an'
    logical, intent(in)             ::  putout
    integer, intent(out)            ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/EMEP_CSO_Simulate'

    ! --- local ----------------------------

    ! time variables:
    type(T_CSO_DateTime)              ::  cso_t1, cso_t2
    type(T_CSO_DateTime)              ::  cso_t
    type(T_CSO_TimeDelta)             ::  cso_dt
    
    ! loop variables:
    integer                           ::  icso
  
    ! number of global pixels:
    integer                           ::  nglb
    ! global pixel arrays:
    real, pointer                     ::  glb_lon(:)       ! (nglb)
    real, pointer                     ::  glb_lat(:)       ! (nglb)
    real, pointer                     ::  glb_clons(:,:)   ! (nglb,4)
    real, pointer                     ::  glb_clats(:,:)   ! (nglb,4)
    logical, pointer                  ::  glb_select(:)    ! (nglb)
    ! loop variable:
    integer                           ::  iglb
    ! corner variables:
    integer                           ::  ncrnr, icrnr

    ! user defined dimensions:
    integer                           ::  nudim
    integer                           ::  iudim
    character(len=32)                 ::  udimname

    ! pixel counter and index:
    integer                           ::  npix
    integer                           ::  ipix 
    ! pixel footprints:
    real, pointer                     ::  lons(:), lats(:)         ! (npix)
    real, pointer                     ::  clons(:,:), clats(:,:)   ! (ncorner,npix)
    ! pixel arrays:
    real, pointer                     ::  data0(:)          ! (npix)
    real, pointer                     ::  data1(:,:)        ! (:,npix)
  
    ! mapping arrays:
    real, pointer                     ::  areas(:)
    integer, pointer                  ::  ii(:), jj(:)
    real, pointer                     ::  ww(:)
    integer, pointer                  ::  iw0(:), nw(:)
    ! mapping index:
    integer                           ::  iw

    ! user defined variables:
    integer                           ::  nuvar
    integer                           ::  iuvar
    character(len=64), allocatable    ::  uvarnames(:)  ! (nuvar)
    character(len=64), allocatable    ::  uvarunits(:)  ! (nuvar)
    
    ! conversion:
    integer                           ::  ispec
    real                              ::  fscale
    logical                           ::  needroa

    ! target file:
    character(len=1024)               ::  cso_output_filename

    ! --- begin -----------------------------
    
    ! make orbit filenames empty, since "Clear" checks on length;
    ! loop over data sets:
    do icso = 1, ncso
      ! no file read:
      cso_orbit_filename(icso) = ''
    end do ! icso
    
    ! current time ;
    ! "current_date" is type "date" with fields: year,hour,day, hour,seconds
    ! (thus no "minutes"!)
    cso_t  = CSO_DateTime( year=current_date%year, month=current_date%month, day=current_date%day, &
                               hour=current_date%hour, sec=current_date%seconds )
    !! testing ...
    !write (gol,'(a,": current time : ",a)') rname, trim(Pretty(cso_t)); call goPr

    ! only at whole hour:
    if ( Precisely( cso_t, 1.0, 'hour' ) ) then

      ! info ...
      write (gol,'(a,": CSO target time : ",a)') rname, trim(Pretty(cso_t)); call goPr

      ! hourly time step:
      cso_dt = CSO_TimeDelta( hour=1 )
      ! time range half-hour before/after current time:
      cso_t1 = cso_t - 0.5 * cso_dt
      cso_t2 = cso_t + 0.5 * cso_dt

      ! loop ..
      do icso = 1, ncso
        ! info ...
        write (gol,'(a,": CSO set `",a,"` ..")') rname, trim(cso_keys(icso)); call goPr

        ! get name of CSO file for orbit with time average
        ! in interval around current time:
        call cso_listing(icso)%SearchFile( cso_t1, cso_t2, 'aver', cso_orbit_filename(icso), status )
        IF_NOT_OK_RETURN(status=1)

        ! orbit file found for this time?
        if ( len_trim(cso_orbit_filename(icso)) > 0 ) then

          ! info ..
          write (gol,'(a,":   orbit file : ",a)') rname, trim(cso_orbit_filename(icso)); call goPr

          ! ~ orbit data

          ! inititalize orbit data,
          ! read all footprints available in file:
          call cso_sdata(icso)%Init( cso_rcF, cso_keys(icso), cso_orbit_filename(icso), status )
          IF_NOT_OK_RETURN(status=1)

          ! obtain from the orbit:
          ! - number of footprints (global, all pixels in file)
          ! - pointer to to selection flags, should be used to select pixels overlapping with local domain
          call cso_sdata(icso)%Get( status, nglb=nglb, glb_select=glb_select )
          IF_NOT_OK_RETURN(status=1)

          ! select pixels that overlap with local domain:
          select case ( trim(cso_mapping_select_by) )

            !~ just on centers:
            case ( 'center' )

              ! obtain from the orbit:
              ! - pointers to arrays with footprint centers
              call cso_sdata(icso)%Get( status, glb_lon=glb_lon, glb_lat=glb_lat )
              IF_NOT_OK_RETURN(status=1)

              ! loop over global pixels:
              do iglb = 1, nglb
                !~ decide if pixel center is part of local domain;
                !  do not adjust to prefered domain since this gives problems at date-line:
                glb_select(iglb) = coord_in_domain( 'local', glb_lon(iglb), glb_lat(iglb), fix=.false. )
              end do ! iglb

            !~ based on entire footprint
            case ( 'footprint' )

              ! obtain from the orbit:
              ! - pointers to arrays with footprint corners
              call cso_sdata(icso)%Get( status, glb_clons=glb_clons, glb_clats=glb_clats )
              IF_NOT_OK_RETURN(status=1)

              ! number of corners:
              ncrnr = size(glb_clons,1)

              ! loop over global pixels:
              do iglb = 1, nglb
                ! by default no overlap:
                glb_select(iglb) = .false.
                !~ check if footprint overlaps (partly) with domain;
                !  do not adjust to prefered domain since this gives problems at date-line:
                do icrnr = 1, ncrnr
                  glb_select(iglb) = glb_select(iglb) .or. &
                      coord_in_domain( 'local' , glb_clons(icrnr,iglb), glb_clats(icrnr,iglb), fix=.false. )
                end do ! icrnr
                ! ... and if it is entirely in global domain;
                !  do not adjust to prefered domain:
                do icrnr = 1, ncrnr
                  glb_select(iglb) = glb_select(iglb) .and. &
                      coord_in_domain( 'global', glb_clons(icrnr,iglb), glb_clats(icrnr,iglb), fix=.false. )
                end do ! icrnr
              end do ! iglb
              
            !~
            case default
              write (gol,'("unsupported selection method `",a,"`")') trim(cso_mapping_select_by); call goPr
              TRACEBACK; status=1; return
          end select

          ! read orbit, locally store only pixels that are flagged in 'glb_select'
          ! as having overlap with this domain:
          call cso_sdata(icso)%Read( cso_rcF, cso_keys(icso), status )
          IF_NOT_OK_RETURN(status=1)

          ! obtain info on track:
          !  - number of local pixels
          call cso_sdata(icso)%Get( status, npix=npix )
          IF_NOT_OK_RETURN(status=1)
          ! info ..
          write (gol,'(a,":   number of local pixels: ",i0)') rname, npix; call goPr

          ! initialize simulation state;
          ! optional arguments:
          !   description='long name'    : used for output attributes
          call cso_sstate(icso)%Init( cso_sdata(icso), cso_rcF, cso_keys(icso), status, &
                                        description='simulated retrievals' )
          IF_NOT_OK_RETURN(status=1)

          ! number of user defined dimensions:
          call cso_sstate(icso)%Get( status, nudim=nudim )
          IF_NOT_OK_RETURN(status=1)
          ! loop:
          do iudim = 1, nudim
            ! get dimension name to be defined:
            call cso_sstate(icso)%GetDim( iudim, status, name=udimname )
            IF_NOT_OK_RETURN(status=1)
            ! switch:
            select case ( trim(udimname) )
              !~ model layers:
              case ( 'model_layer' )
                ! define number of model layers, one extra for strato:
                call cso_sstate(icso)%SetDim( udimname, KMAX_MID+1, status )
                IF_NOT_OK_RETURN(status=1)
              !~ model layer interfaces:
              case ( 'model_layeri' )
                ! define number of layer interfaces (model layers + 1 extra strato layer):
                call cso_sstate(icso)%SetDim( udimname, KMAX_MID+2, status )
                IF_NOT_OK_RETURN(status=1)
              !~ unknown ...
              case default
                write (gol,'("unsupported dimension `",a,"`")') trim(udimname); call goErr
                TRACEBACK; status=1; return
            end select
          end do ! idim

          ! dimension defined, allocate storage:
          call cso_sstate(icso)%EndDef( status )
          IF_NOT_OK_RETURN(status=1)

          ! any pixels?
          if ( npix > 0 ) then

            ! pointers to (*,pixel) arrays:
            ! - footprint centers  (not used, for inspiration ..)
            ! - footprint corners
            ! - half level pressure profiles
            call cso_sdata(icso)%Get( status, lons=lons, lats=lats, clons=clons, clats=clats )
            IF_NOT_OK_RETURN(status=1)

            ! info ...
            write (gol,'(a,": compute mapping weights ...")') rname; call goPr
            ! get pointers to mapping arrays for all pixels:
            ! - areas(1:npix)            : pixel area [m2]
            ! - iw0(1:npix), nw(1:npix)  : offset and number of elements in ii/jj/ww
            ! - ii(:), jj(:), ww(:)      : cell and weight arrays for mapping to footprint,
            call emep_grid_mapping%GetWeights( clons, clats, &
                                                areas, iw0, nw, ii, jj, ww, status )
            IF_NOT_OK_RETURN(status=1)

            ! info ...
            write (gol,'(a,":   store ...")') rname; call goPr
            ! store mapping weights, might be saved to check,
            ! or used for swap operation;
            ! cell indices ii/jj need to be the global index numbers;
            ! here use that (gi0,gj0) is 1-based index of the lower-left cell 
            ! in the global index space:
            call cso_sdata(icso)%SetMapping( areas, nw, &
                                              gi0-1+ii, gj0-1+jj, ww, &
                                              status )
            IF_NOT_OK_RETURN(status=1)

            ! number of user defined variables:
            call cso_sstate(icso)%Get( status, nuvar=nuvar, outkey=outkey )
            IF_NOT_OK_RETURN(status=1)
            ! any user defined?
            if ( nuvar > 0 ) then
  
              ! storage for variable names and units:
              allocate( uvarnames(nuvar), stat=status )
              IF_NOT_OK_RETURN(status=1)
              allocate( uvarunits(nuvar), stat=status )
              IF_NOT_OK_RETURN(status=1)
              ! fill:
              call cso_sstate(icso)%Get( status, uvarnames=uvarnames, uvarunits=uvarunits, outkey=outkey )
              IF_NOT_OK_RETURN(status=1)
  
              ! loop over variables:
              do iuvar = 1, nuvar
                ! info ..
                write (gol,'(a,": user defined variable: ",a)') rname, trim(uvarnames(iuvar)); call goPr
                ! switch:
                select case ( trim(uvarnames(iuvar)) )
  
                  !~ model concentrations
                  case ( 'mod_conc', 'mod_conc_an' )
                  
                    ! info ...
                    write (gol,'(a,":   simulate concentrations [",a,"] from tracer `",a,"` ...")') &
                            rname, trim(uvarunits(iuvar)), trim(cso_tracer(icso)); call goPr
  
                    ! find index of tracer:
                    ispec = find_index( cso_tracer(icso), species_adv(:)%name )
  
                    ! conversion factor, multipication with air density?
                    call Units_Scale( uvarunits(iuvar), ispec, fscale, &
                                             needroa=needroa, debug_msg=TRACELINE )
  
                    ! get pointer to target array with shape (nlev+1,npix):
                    call cso_sstate(icso)%GetData( status, name=uvarnames(iuvar), data1=data1 )
                    IF_NOT_OK_RETURN(status=1)
                    ! loop over pixels:
                    do ipix = 1, npix
                      ! any source contributions?
                      if ( nw(ipix) > 0 ) then
                        ! init sum, strato in 1 (top down order!) will remain zero:
                        data1(:,ipix) = 0.0
                        ! loop over source contributions:
                        do iw = iw0(ipix)+1, iw0(ipix)+nw(ipix)
                          ! convert concentrations using air density?
                          if ( needroa ) then
                            data1(2:KMAX_MID+1,ipix) = data1(2:KMAX_MID+1,ipix) &
                                                     + xn_adv(ispec,ii(iw),jj(iw),:) &
                                                         * roa(ii(iw),jj(iw),:,1) &
                                                         * fscale * ww(iw)/areas(ipix)
                          else
                            data1(2:KMAX_MID+1,ipix) = data1(2:KMAX_MID+1,ipix) &
                                                     + xn_adv(ispec,ii(iw),jj(iw),:) &
                                                         * fscale * ww(iw)/areas(ipix)
                          end if ! need roa
                        end do ! iw
                      end if ! nw > 0
                    end do ! ipix
  
                  !~ model half-level pressures:
                  case ( 'mod_hp' )
                  
                    ! info ...
                    write (gol,'(a,":   simulate pressure using hybride coeffs ...")') rname; call goPr

                    ! check units:
                    if ( trim(uvarunits(iuvar)) /= 'Pa' ) then
                      write (gol,'(a,": variable `",a,"` requires conversion to `",a,"` from `",a,"`")') &
                              rname, trim(uvarnames(iuvar)), trim(uvarunits(iuvar)), 'Pa'; call goErr
                    end if
                    ! get pointer to target array with shape (nlev+1,npix):
                    call cso_sstate(icso)%GetData( status, name=uvarnames(iuvar), data1=data1 )
                    IF_NOT_OK_RETURN(status=1)
                    ! loop over pixels:
                    do ipix = 1, npix
                      ! any source contributions?
                      if ( nw(ipix) > 0 ) then
                        ! init sum, top of atmosphere in 1 (top down order!) will remain zero:
                        data1(:,ipix) = 0.0
                        ! loop over source contributions:
                        do iw = iw0(ipix)+1, iw0(ipix)+nw(ipix)
                          ! add contribution,
                          ! comppute half level pressures using hybride coeff;
                          ! first record of surface pressure is said to be the correct one ..
                          data1(2:KMAX_MID+2,ipix) = data1(2:KMAX_MID+2,ipix) + (A_bnd + B_bnd * ps(ii(iw),jj(iw),1)) * ww(iw)/areas(ipix)
                        end do ! iw
                      end if ! nw > 0
                    end do ! ipix
  
                  !~
                  case default
                    write (gol,'("unsupported variable name `",a,"`")') trim(uvarnames(iuvar)); call goErr
                    TRACEBACK; status=1; return
                end select
  
              end do ! user defined variables
  
              ! info ...
              write (gol,'("end ...")'); call goPr
  
              ! clear:
              deallocate( uvarnames, stat=status )
              IF_NOT_OK_RETURN(status=1)
              deallocate( uvarunits, stat=status )
              IF_NOT_OK_RETURN(status=1)
  
            end if ! nuvar > 0

          end if ! npix > 0

          ! ~ formula

          ! inquire info on pixels covering multiple domains,
          ! setup exchange parameters:
          call cso_sdata(icso)%SetupExchange( status )
          IF_NOT_OK_RETURN(status=1)
          ! exchange simulations:
!... selection
          call cso_sstate(icso)%Exchange( cso_sdata(icso), status )
          IF_NOT_OK_RETURN(status=1)
  
          ! apply formula: kernel convolution etc:
!... selection
          call cso_sstate(icso)%ApplyFormulas( cso_sdata(icso), status )
          IF_NOT_OK_RETURN(status=1)
  
          ! ~ put out
          
          ! enabled?
          if ( putout ) then
  
            ! target file for selected data:
            write (cso_output_filename,'("CSO_output_",i4.4,2i2.2,"_",2i2.2,"_",a,"_data.nc")') &
                       cso_t%year, cso_t%month, cso_t%day, cso_t%hour, cso_t%minute, &
                       trim(cso_keys(icso))
            ! setup output arrays and output weights:
            call cso_sdata(icso)%PutOut( cso_output_filename, status )
            IF_NOT_OK_RETURN(status=1)

            ! target file for simulation state:
            write (cso_output_filename,'("CSO_output_",i4.4,2i2.2,"_",2i2.2,"_",a,"_state.nc")') &
                       cso_t%year, cso_t%month, cso_t%day, cso_t%hour, cso_t%minute, &
                       trim(cso_keys(icso))
            ! write output selection:
            call cso_sstate(icso)%PutOut( cso_sdata(icso), cso_output_filename, status )
            IF_NOT_OK_RETURN(status=1)
            
          end if ! put out
          
          !! testing ..        
          !write (gol,'("break after first cso simulations ...")'); call goErr
          !TRACEBACK; status=1; return
  
          ! ~ clear
  
          ! clear:
          !nullify( glb_lon )
          !nullify( glb_lat )
          nullify( glb_clons )
          nullify( glb_clats )
          nullify( glb_select )

        else

          ! info ..
          write (gol,'(a,":   no orbit file found ...")') rname; call goPr
          !! no pixels:
          !npix = 0

        end if  ! orbit found in listing

      end do ! cso sets
      
    end if  ! time at whole hour
    
    ! ok
    status = 0

  end subroutine EMEP_CSO_Simulate


  ! ***


  subroutine EMEP_CSO_Clear( status )

    ! --- in/out ----------------------------

    integer, intent(out)           ::  status

    ! --- const ----------------------------

    character(len=*), parameter  ::  rname = mname//'/EMEP_CSO_Clear'

    ! --- local ----------------------------

    integer                               ::  icso

    ! --- begin -----------------------------
    
    ! info ...
    write (gol,'(a,": clear CSO data ...")') rname; call goPr

    ! loop over datasets:
    do icso = 1, ncso

      ! orbit file found for this time?
      if ( len_trim(cso_orbit_filename(icso)) > 0 ) then

        ! done with state:
        call cso_sstate(icso)%Done( cso_sdata(icso), status )
        IF_NOT_OK_RETURN(status=1)

        ! done with data:
        call cso_sdata(icso)%Done( status )
        IF_NOT_OK_RETURN(status=1)

      end if  ! orbit found in listing

    end do  ! icso

    ! ok
    status = 0

  end subroutine EMEP_CSO_Clear
  
  
  ! ====================================================================
  ! ===
  ! === GridMapping
  ! ===
  ! ====================================================================

  !
  ! Original file: cso_tools.F90
  ! This version uses the grid cel search routine from EMEP model.
  !

  subroutine EMEP_GridMapping_Init( self, nlon, nlat, levels, status )

    use CSO_PArray, only : CSO_PArray_Init
  
    ! --- in/out ---------------------------------
    
    class(T_EMEP_GridMapping), intent(out)    ::  self
    integer, intent(in)                       ::  nlon, nlat
    integer, intent(in)                       ::  levels
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/EMEP_GridMapping_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%nlon  = nlon
    self%nlat  = nlat
    
    ! store:
    self%levels = levels
    
    ! maximum storage:
    allocate( self%ii(self%nlon*self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%jj(self%nlon*self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( self%ww(self%nlon*self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
              
    ! allocate 2D array with weight sum:
    allocate( self%wmap(self%nlon,self%nlat), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! no arrays for "all" pixels yet:
    self%nall = 0
    self%mxall = 0
    call CSO_PArray_Init( self%all_ii, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%all_jj, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%all_ww, status )
    IF_NOT_OK_RETURN(status=1)
    ! no area per pixel yet:
    call CSO_PArray_Init( self%all_iw0, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%all_nw, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%all_area, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine EMEP_GridMapping_Init
  
  
  ! *
  

  subroutine EMEP_GridMapping_Done( self, status )

    use CSO_PArray, only : CSO_PArray_Done
  
    ! --- in/out ---------------------------------
    
    class(T_EMEP_GridMapping), intent(inout)  ::  self
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/EMEP_GridMapping_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! sampling points?
    if ( allocated(self%xxp) ) then
      deallocate( self%xxp, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%yyp, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%wwp, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    
    ! source info:
    deallocate( self%ii, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%jj, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%ww, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( self%wmap, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear arrays per pixel:
    call CSO_PArray_Done( self%all_iw0, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Done( self%all_nw, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Done( self%all_area, status )
    IF_NOT_OK_RETURN(status=1)
    ! clear mapping arrays for all pixels:
    call CSO_PArray_Done( self%all_ii, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Done( self%all_jj, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Done( self%all_ww, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine EMEP_GridMapping_Done
  
  ! *
  
  !
  ! Devide polygon defined by corners (xx(:),yy(:)) in triangles or quadrangles,
  ! and store centroids (xxp(:),yyp(:)) such that these can be used as sampling points.
  ! Assign centroids (xxp(:),yyp(:)) to grid cells (ii(1:n),jj(1:n)),
  ! the weights ww(1:n) have approximately the part of the polygon area
  ! that covers a grid cell; the total polygon area is returned in "area"
  ! Note that the polygon might only partly overlap with the domain, thus:
  !   sum(ww) <= area
  ! On input the arrays ii/jj/ww should have sufficient size,
  ! number of elements filled on exit is n.
  !
  
  subroutine EMEP_GridMapping_GetWeights_0d( self, xx, yy, area, n, ii, jj, ww, status )
  
    use CSO                 , only : GetPolygonPoints
    use GridValues_mod      , only : coord_in_domain

    ! --- in/out ---------------------------------
    
    class(T_EMEP_GridMapping), intent(inout)  ::  self
    real, intent(in)                          ::  xx(:), yy(:)  ! [degree]
    real, intent(out)                         ::  area          ! [m2]
    integer, intent(out)                      ::  n
    integer, pointer                          ::  ii(:), jj(:)  ! (n)
    real, pointer                             ::  ww(:)         ! (n) [m2]
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/EMEP_GridMapping_GetGridWeights_0d'
    
    ! --- local ----------------------------------
    
    integer                 ::  i1, i2, j1, j2
    integer                 ::  ip
    logical                 ::  inside
    integer                 ::  i, j
    integer                 ::  k
    integer                 ::  np

    ! --- begin ----------------------------------
        
    ! devide footprint in triangles or quadrangles,
    ! return centroids that can be used as sampling points,
    ! arrays xxp and yyp are allocated on output to maximum size needed:
    call GetPolygonPoints( xx, yy, self%xxp, self%yyp, self%wwp, status, levels=self%levels )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of sampling points:
    np = size(self%xxp)
    
    ! total area of polygon:
    area = sum(self%wwp)
 
    !! index box:
    !i1 = max(         1, ceiling( (minval(self%xxp) - self%west )/self%dlon ) - self%ilon0 )
    !i2 = min( self%nlon, ceiling( (maxval(self%xxp) - self%west )/self%dlon ) - self%ilon0 )
    !j1 = max(         1, ceiling( (minval(self%yyp) - self%south)/self%dlat ) - self%ilat0 )
    !j2 = min( self%nlat, ceiling( (maxval(self%yyp) - self%south)/self%dlat ) - self%ilat0 )
    
    ! init range:
    i1 = self%nlon + 1000
    i2 = -999
    j1 = self%nlat + 1000
    j2 = -999
    ! reset sum:
    self%wmap = 0    
    ! loop over points:
    do ip = 1, np
      !! target cell:
      !i = ceiling( (self%xxp(ip)-self%west )/self%dlon ) - self%ilon0
      !j = ceiling( (self%yyp(ip)-self%south)/self%dlat ) - self%ilat0
      ! test if location is in local domain, also return indices (if defined),
      ! do not adjust to prefered domain since this gives problems at date-line:
      inside = coord_in_domain( 'local', self%xxp(ip), self%yyp(ip), &
                                  iloc=i, jloc=j, fix=.false. )
      ! in local domain?
      if ( inside ) then
        ! increase weights on map:
        self%wmap(i,j) = self%wmap(i,j) + self%wwp(ip)
        ! update range:
        i1 = min( i, i1 )
        i2 = max( i, i2 )
        j1 = min( j, j1 )
        j2 = max( j, j2 )
      end if
    end do

    ! number of cells with contributions:
    n = count( self%wmap > 0.0 )
    ! loop over grid cells:
    k = 0
    do i = i1, i2
      do j = j1, j2
        ! any contribution?
        if ( self%wmap(i,j) > 0.0 ) then
          ! increase counter:
          k = k + 1
          ! store location:
          self%ii(k) = i
          self%jj(k) = j
          ! copy total area over this cell:
          self%ww(k) = self%wmap(i,j)
        end if
      end do
    end do
    
    ! assign pointers:
    ii => self%ii
    jj => self%jj
    ww => self%ww
    
    ! ok
    status = 0
    
  end subroutine EMEP_GridMapping_GetWeights_0d

  ! *
  
  !
  ! Idem for arrays with pixels.
  ! Input:
  !   xx, yy  : pixel footprints
  ! Ouptut:
  !   area       : pixel area
  !   iw0, nw    : per pixel the offset and number of elements in ii/jj/ww arrays
  !   ii, jj     : source cell indices
  !   ww         : source cell weights
  !
  
  subroutine EMEP_GridMapping_GetWeights_1d( self, xx, yy, &
                                               area, iw0, nw, ii, jj, ww, status )

    use CSO_PArray, only : CSO_PArray_Reshape

    ! --- in/out ---------------------------------
    
    class(T_EMEP_GridMapping), intent(inout)  ::  self
    real, intent(in)                          ::  xx(:,:), yy(:,:)  ! (ncorner,npix) [degree]
    real, pointer                             ::  area(:)           ! (npix) [m2]
    integer, pointer                          ::  iw0(:)            ! (npix)
    integer, pointer                          ::  nw(:)             ! (npix)
    integer, pointer                          ::  ii(:), jj(:)      ! (nw)
    real, pointer                             ::  ww(:)             ! (nw) [m2]
    integer, intent(out)                      ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/EMEP_GridMapping_GetWeights_1d'
    
    ! --- local ----------------------------------
    
    integer                   ::  npix
    integer                   ::  ipix
    real                      ::  pix_area
    integer, pointer          ::  pix_ii(:), pix_jj(:)
    real, pointer             ::  pix_ww(:)
    integer                   ::  pix_nw
    integer                   ::  nnew

    ! --- begin ----------------------------------

    ! number of pixels:
    npix = size(xx,2)
    ! check ...
    if ( size(yy,2) /= npix ) then
      write (gol,'("arrays xx (",i0,",",i0,") and yy (",i0,",",i0,") should have same shape")') &
                size(xx,1), size(xx,2), size(yy,1), size(yy,2); call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage for pixel area; check current storage:
    call CSO_PArray_Reshape( self%all_area, npix, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%all_iw0 , npix, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%all_nw  , npix, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! reset counter:
    self%nall = 0
    self%all_nw = 0
    
    ! loop over pixels:
    do ipix = 1, npix
      ! pixel mapping weights:
      call self%GetWeights( xx(:,ipix), yy(:,ipix), &
                             pix_area, pix_nw, pix_ii, pix_jj, pix_ww, status )
      IF_NOT_OK_RETURN(status=1)

      ! store pixel area:
      self%all_area(ipix) = pix_area

      ! offset and number of overlapping cells (might be 0 ...):
      self%all_iw0 (ipix) = self%nall
      self%all_nw  (ipix) = pix_nw

      ! any overlap? some pixels might not overlap with domain ...
      if ( pix_nw > 0 ) then
        ! exceeds maximum storage?
        if ( self%nall + pix_nw > self%mxall ) then
          ! new size, extend with 1 value extra per cell until it fits ...
          do
            self%mxall = self%mxall + self%nlon*self%nlat
            if ( self%nall + pix_nw <= self%mxall ) exit
          end do
          ! extend arrays, copy current:
          call CSO_PArray_Reshape( self%all_ii , self%mxall, status )
          IF_NOT_OK_RETURN(status=1)
          call CSO_PArray_Reshape( self%all_jj , self%mxall, status )
          IF_NOT_OK_RETURN(status=1)
          call CSO_PArray_Reshape( self%all_ww , self%mxall, status )
          IF_NOT_OK_RETURN(status=1)
        end if
        ! store pixel mapping:
        self%all_ii  (self%nall+1:self%nall+pix_nw) = pix_ii(1:pix_nw)
        self%all_jj  (self%nall+1:self%nall+pix_nw) = pix_jj(1:pix_nw)
        self%all_ww  (self%nall+1:self%nall+pix_nw) = pix_ww(1:pix_nw)
        ! increase counter:
        self%nall = self%nall + pix_nw
      end if ! nw > 0
      
    end do ! ipix
    
    ! truncate to length that is actually used:
    call CSO_PArray_Reshape( self%all_ii , self%nall, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%all_jj , self%nall, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%all_ww , self%nall, status )
    IF_NOT_OK_RETURN(status=1)
    ! reset maximum size:
    self%mxall = self%nall

    ! set pointers:
    area => self%all_area
    iw0  => self%all_iw0
    nw   => self%all_nw
    ii   => self%all_ii
    jj   => self%all_jj
    ww   => self%all_ww
    
    ! ok
    status = 0
    
  end subroutine EMEP_GridMapping_GetWeights_1d


end module EMEP_CSO_mod


