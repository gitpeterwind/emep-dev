!###############################################################################
!
! CSO_Grid - tools for regular lon/lat grid
!
! HISTORY
!
!   2022-09, Arjo Segers
!     Updated call to initialize file for writing.
!     Fixed wrong gather of lon/lat arrays.
!     Support input and output of packed variables.
!
!   2023-01, Arjo Segers
!     Added "GetWeightsCell" methodes to assign points to grid cells.
!
!  
!###############################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!###############################################################################


module CSO_Grid

  use CSO_Logging     , only : csol, csoPr, csoErr

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private
  
  public    ::  T_GridMapping
  public    ::  GridSampleFile
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Grid'
  
  ! --- types ------------------------------
  
  type :: T_GridMapping
    ! grid description:
    real                      ::  west, south
    real                      ::  dlon, dlat
    integer                   ::  ilon0, ilat0
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
    ! idem for arrays for pixel arrays:
    integer                   ::  nall
    integer                   ::  mxall
    integer, pointer          ::  all_iw0(:), all_nw(:)   ! (npix)
    integer, pointer          ::  all_ii(:), all_jj(:)   ! (mxall)
    real, pointer             ::  all_ww(:)    ! (mxall)
    real, pointer             ::  all_area(:)  ! (npix)
    !
  contains
    procedure ::  Init            =>  GridMapping_Init
    procedure ::  Done            =>  GridMapping_Done
    procedure ::                      GridMapping_GetWeights_0d
    procedure ::                      GridMapping_GetWeights_1d
    generic   ::  GetWeights      =>  GridMapping_GetWeights_0d, &
                                      GridMapping_GetWeights_1d
    procedure ::                      GridMapping_GetWeightsCell_0d
    procedure ::                      GridMapping_GetWeightsCell_1d
    generic   ::  GetWeightsCell  =>  GridMapping_GetWeightsCell_0d, &
                                      GridMapping_GetWeightsCell_1d
  end type T_GridMapping

  
contains


  ! ====================================================================
  ! ===
  ! === GridMapping
  ! ===
  ! ====================================================================


  !
  ! Regular grid defined by:
  ! - south west corner (corner of first grid cell) in global domain
  ! - grid resolution
  ! - grid shape
  ! - cell offset in global domain
  !

  subroutine GridMapping_Init( self, west , dlon, ilon0, nlon, &
                                     south, dlat, ilat0, nlat, &
                                     levels, status )

    use CSO_PArray, only : CSO_PArray_Init
  
    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(out)    ::  self
    real, intent(in)                     ::  west, south
    real, intent(in)                     ::  dlon, dlat
    integer, intent(in)                  ::  ilon0, ilat0
    integer, intent(in)                  ::  nlon, nlat
    integer, intent(in)                  ::  levels
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! store:
    self%west  = west
    self%dlon  = dlon
    self%ilon0 = ilon0
    self%nlon  = nlon
    self%south = south
    self%dlat  = dlat
    self%ilat0 = ilat0
    self%nlat  = nlat
    
    ! store:
    self%levels = levels
    
    ! initialize with maximum storage for single mapping: all cells needed
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
    
  end subroutine GridMapping_Init
  
  
  ! *
  

  subroutine GridMapping_Done( self, status )

    use CSO_PArray, only : CSO_PArray_Done
  
    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! sampling points?
    if ( allocated(self%xxp) ) then
      deallocate( self%xxp, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%yyp, stat=status )
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
    
  end subroutine GridMapping_Done
  
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
  
  subroutine GridMapping_GetWeights_0d( self, xx, yy, area, n, ii, jj, ww, status )
  
    use CSO_Tools, only : GetPolygonPoints

    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    real, intent(in)                     ::  xx(:), yy(:)  ! [degree]
    real, intent(out)                    ::  area          ! [m2]
    integer, intent(out)                 ::  n
    integer, pointer                     ::  ii(:), jj(:)  ! (n)
    real, pointer                        ::  ww(:)         ! (n)  [m2]
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_GetWeights_0d'
    
    ! --- local ----------------------------------
    
    integer                 ::  i1, i2, j1, j2
    integer                 ::  ip
    integer                 ::  i, j
    integer                 ::  k
    integer                 ::  np

    ! --- begin ----------------------------------
        
    ! devide polygon in triangles or quadrangles,
    ! return centroids that can be used as sampling points,
    ! arrays xxp and yyp are allocated on output to maximum size needed:
    call GetPolygonPoints( xx, yy, self%xxp, self%yyp, self%wwp, status, levels=self%levels )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of sampling points:
    np = size(self%xxp)
    
    ! total area of polygon:
    area = sum(self%wwp)
 
    ! index box:
    i1 = max(         1, ceiling( (minval(self%xxp) - self%west )/self%dlon ) - self%ilon0 )
    i2 = min( self%nlon, ceiling( (maxval(self%xxp) - self%west )/self%dlon ) - self%ilon0 )
    j1 = max(         1, ceiling( (minval(self%yyp) - self%south)/self%dlat ) - self%ilat0 )
    j2 = min( self%nlat, ceiling( (maxval(self%yyp) - self%south)/self%dlat ) - self%ilat0 )
    
    ! reset sum:
    self%wmap = 0.0
    ! loop over points:
    do ip = 1, np
      ! target cell:
      i = ceiling( (self%xxp(ip)-self%west )/self%dlon ) - self%ilon0
      j = ceiling( (self%yyp(ip)-self%south)/self%dlat ) - self%ilat0
      ! in range?
      if ( (i >= 1) .and. (i <= self%nlon) .and. (j >= 1) .and. (j <= self%nlat) ) then
        ! increase area on map:
        self%wmap(i,j) = self%wmap(i,j) + self%wwp(ip)
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
    
  end subroutine GridMapping_GetWeights_0d

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
  
  subroutine GridMapping_GetWeights_1d( self, xx, yy, &
                                         area, iw0, nw, ii, jj, ww, status )

    use CSO_PArray, only : CSO_PArray_Reshape

    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    real, intent(in)                     ::  xx(:,:), yy(:,:)  ! (ncorner,npix) [degree]
    real, pointer                        ::  area(:)           ! (npix) [m2]
    integer, pointer                     ::  iw0(:)            ! (npix)
    integer, pointer                     ::  nw(:)             ! (npix)
    integer, pointer                     ::  ii(:), jj(:)      ! (nw)
    real, pointer                        ::  ww(:)             ! (nw) [m2]
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_GetWeights_1d'
    
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
      write (csol,'("arrays xx (",i0,",",i0,") and yy (",i0,",",i0,") should have same shape")') &
                size(xx,1), size(xx,2), size(yy,1), size(yy,2); call csoErr
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
      if ( status /= 0 ) then
        write (csol,'("could not compute mapping weights for pixel ",i0)') ipix; call csoErr
        write (csol,*) '  xx = ', xx(:,ipix); call csoErr
        write (csol,*) '  yy = ', yy(:,ipix); call csoErr
        TRACEBACK; status=ipix; return
      end if

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
    
  end subroutine GridMapping_GetWeights_1d
  
  !
  ! ***
  !
  
  !
  ! Assign points (x,y) to grid cell (i,j).
  ! The interpolation weight w is set to 1.0.
  ! The area is for footprints used to normalize,
  ! here set to 1.0 to avoid division by zero.
  !
  
  subroutine GridMapping_GetWeightsCell_0d( self, x, y, area, i, j, w, status )
  
    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    real, intent(in)                     ::  x, y          ! [degree]
    real, intent(out)                    ::  area
    integer, intent(out)                 ::  i, j
    real , intent(out)                   ::  w
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_GetWeightsCell_0d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
        
    ! target cell, first cell if exactly on first edge:
    if ( x == self%west ) then
      i =               1                     - self%ilon0
    else
      i = ceiling( (x-self%west )/self%dlon ) - self%ilon0
    end if
    if ( y == self%south ) then
      j =               1                     - self%ilat0
    else
      j = ceiling( (y-self%south)/self%dlat ) - self%ilat0
    end if
    ! in range?
    if ( (i >= 1) .and. (i <= self%nlon) .and. (j >= 1) .and. (j <= self%nlat) ) then
      ! full weight:
      w = 1.0
    else
      ! outside map, no weight:
      w = 0.0
    end if
    
    ! dummy "pixel" area; used for normization with weights, so set to unity:
    area = 1.0

    ! ok
    status = 0
    
  end subroutine GridMapping_GetWeightsCell_0d

  ! *
  
  !
  ! Array of pixel centers, select corresponding cell.
  ! Input:
  !   x, y       : pixel centers
  ! Ouptut:
  !   area       : pixel area (here 1.0, corners are unknown)
  !   iw0, nw    : per pixel the offset and number of elements in ii/jj/ww arrays ;
  !                here nw==1 for all pixels
  !   ii, jj     : source cell indices
  !   ww         : source cell weights (here always 1.0)
  !
  
  subroutine GridMapping_GetWeightsCell_1d( self, x, y, &
                                              area, iw0, nw, ii, jj, ww, status )

    use CSO_PArray, only : CSO_PArray_Reshape

    ! --- in/out ---------------------------------
    
    class(T_GridMapping), intent(inout)  ::  self
    real, intent(in)                     ::  x(:), y(:)        ! (npix) [degree]
    real, pointer                        ::  area(:)           ! (npix) [m2]
    integer, pointer                     ::  iw0(:)            ! (npix)
    integer, pointer                     ::  nw(:)             ! (npix)
    integer, pointer                     ::  ii(:), jj(:)      ! (nw)
    real, pointer                        ::  ww(:)             ! (nw) [m2]
    integer, intent(out)                 ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridMapping_GetWeightsCell_1d'
    
    ! --- local ----------------------------------
    
    integer                   ::  npix
    integer                   ::  ipix
    real                      ::  pix_area
    integer                   ::  pix_i, pix_j
    real                      ::  pix_w
    integer                   ::  pix_nw
    integer                   ::  nnew

    ! --- begin ----------------------------------

    ! number of pixels:
    npix = size(x)
    ! check ...
    if ( size(y) /= npix ) then
      write (csol,'("arrays x (",i0,") and y (",i0,") should have same shape")') &
                size(x), size(y); call csoErr
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
      call self%GetWeightsCell( x(ipix), y(ipix), &
                                   pix_area, pix_i, pix_j, pix_w, status )
      if ( status /= 0 ) then
        write (csol,'("could not compute mapping weights for pixel ",i0)') ipix; call csoErr
        write (csol,*) '  xx = ', x(ipix); call csoErr
        write (csol,*) '  yy = ', y(ipix); call csoErr
        TRACEBACK; status=ipix; return
      end if
      ! single source cell:
      pix_nw = 1

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
        self%all_ii  (self%nall+1:self%nall+pix_nw) = pix_i
        self%all_jj  (self%nall+1:self%nall+pix_nw) = pix_j
        self%all_ww  (self%nall+1:self%nall+pix_nw) = pix_w
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
    
  end subroutine GridMapping_GetWeightsCell_1d
  
  
  ! ====================================================================
  ! ===
  ! === grid defintion file
  ! ===
  ! ====================================================================
  
  !
  ! Create netcdf file with grid definition.
  ! Used by postprocessing for averaging pixels over grid.
  !
  
  subroutine GridSampleFile( filename, lons, lats, area, status, &
                                 ilon0, ilat0 )
  
    use CSO_NcFile, only : T_NcFile

    ! --- in/out ---------------------------------
    
    character(len=*), intent(in)          ::  filename
    real, intent(in)                      ::  lons(:)
    real, intent(in)                      ::  lats(:)
    real, intent(in)                      ::  area(:,:)
    integer, intent(out)                  ::  status
    
    ! offsets in global grid:
    integer, intent(in)                   ::  ilon0, ilat0

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/GridSampleFile'
    
    ! --- local ----------------------------------
  
    ! grid definition file:
    type(T_NcFile)                    ::  GridFile
    integer                           ::  nlon, nlat
    integer                           ::  ivar_lon, ivar_lat, ivar_area

    ! --- begin ----------------------------------
  
    ! init file for creation:
    call GridFile%Init( trim(filename), 'w', status )
    IF_NOT_OK_RETURN(status=1)

    ! global attributes:
    call GridFile%Set_Attr( 0, 'conventions', 'CF-1.7', status )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Set_Attr( 0, 'title', 'CSO grid description', status )
    IF_NOT_OK_RETURN(status=1)
    
    ! shape:
    nlon = size(lons)
    nlat = size(lats)

    ! add dimensions, provide local size and offset in global grid:
    call GridFile%Def_Dim( 'longitude', nlon, status, offset=ilon0 )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Def_Dim( 'latitude', nlat, status, offset=ilat0 )
    IF_NOT_OK_RETURN(status=1)

    ! variable, get internal index back:
    call GridFile%Def_Var( 'longitude', (/'longitude'/), status, ivar=ivar_lon )
    IF_NOT_OK_RETURN(status=1)
    !~ add attributes:
    call GridFile%Set_Attr( ivar_lon, 'standard_name', 'longitude', status )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Set_Attr( ivar_lon, 'units', 'degrees_east', status )
    IF_NOT_OK_RETURN(status=1)

    ! variable, get internal index back:
    call GridFile%Def_Var( 'latitude', (/'latitude'/), status, ivar=ivar_lat )
    IF_NOT_OK_RETURN(status=1)
    !~ add attributes:
    call GridFile%Set_Attr( ivar_lat, 'standard_name', 'latitude', status )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Set_Attr( ivar_lat, 'units', 'degrees_north', status )
    IF_NOT_OK_RETURN(status=1)

    ! variable, get internal index back:
    call GridFile%Def_Var( 'cell_area', (/'longitude','latitude '/), status, ivar=ivar_area )
    IF_NOT_OK_RETURN(status=1)
    !~ add attributes:
    call GridFile%Set_Attr( ivar_area, 'standard_name', 'area', status )
    IF_NOT_OK_RETURN(status=1)
    call GridFile%Set_Attr( ivar_area, 'units', 'm2', status )
    IF_NOT_OK_RETURN(status=1)

    ! define:
    call GridFile%EndDef( status )
    IF_NOT_OK_RETURN(status=1)

    ! write 1D array, gathered on root from first processor row:
    call GridFile%Put_Var1D( ivar_lon, lons, status, empty=(ilat0 > 0) )
    IF_NOT_OK_RETURN(status=1)
    ! write 1D array, gathered on root from first processor column:
    call GridFile%Put_Var1D( ivar_lat, lats, status, empty=(ilon0 > 0) )
    IF_NOT_OK_RETURN(status=1)
    ! write 2D array, gathered on root:
    call GridFile%Put_Var2D( ivar_area, area, status )
    IF_NOT_OK_RETURN(status=1)

    ! close:
    call GridFile%Done( status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine GridSampleFile


end module CSO_Grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! End
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
