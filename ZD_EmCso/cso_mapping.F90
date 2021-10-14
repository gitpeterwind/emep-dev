!###############################################################################
!
! Mapping from model cells to pixel footprings
!
!###############################################################################
!
#define TRACEBACK write (csol,'("in ",a," (",a,", line",i5,")")') rname, __FILE__, __LINE__; call csoErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
#define IF_ERROR_RETURN(action) if (status> 0) then; TRACEBACK; action; return; end if
!
#define IF_NF90_NOT_OK_RETURN(action) if (status/=NF90_NOERR) then; csol=nf90_strerror(status); call csoErr; TRACEBACK; action; return; end if
!
#include "cso.inc"
!
!###############################################################################

module CSO_Mapping

  use CSO_Logging, only : csol, csoPr, csoErr
  use NetCDF     , only : NF90_StrError, NF90_NOERR

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private

  public  ::  T_Mapping
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Mapping'


  ! --- types --------------------------------
  
  !
  ! mapping from grid cells to pixels
  !
  type T_Mapping
    ! count:
    integer                           ::  npix
    ! number of grid cell indices used for mapping for each pixel:
    integer, pointer                  ::  map_n(:)      ! (npix)
    integer, pointer                  ::  map_i0(:)     ! (npix)
    ! total number of mappings:
    integer                           ::  nmap          ! mappings from local domain
    integer                           ::  nmap_all      ! mappings from all domains
    ! grid cell indices and weight:
    integer, pointer                  ::  map_iglb(:)    ! (nmap+)
    integer, pointer                  ::  map_ii(:)      ! (nmap+)
    integer, pointer                  ::  map_jj(:)      ! (nmap+)
    real, pointer                     ::  map_ww(:)      ! (nmap+)
    ! output data, on root only:
    integer, allocatable              ::  map_n_out(:)    ! (nout)
    integer, allocatable              ::  map_ii_out(:)   ! (nmap_all)
    integer, allocatable              ::  map_jj_out(:)   ! (nmap_all)
    real, allocatable                 ::  map_ww_out(:)   ! (nmap_all)
    ! netcdf variables:
    integer                           ::  varid_n_out
    integer                           ::  varid_ii_out
    integer                           ::  varid_jj_out
    integer                           ::  varid_ww_out
    !
  contains
    procedure :: Init            => Mapping_Init
    procedure :: InitSwap        => Mapping_InitSwap
    procedure :: Done            => Mapping_Done
    procedure :: GetPixel        => Mapping_GetPixel
    procedure :: SetPixel        => Mapping_SetPixel
    procedure :: SetData         => Mapping_SetData
    procedure :: Collect         => Mapping_Collect
    procedure :: NcDef           => Mapping_NcDef
    procedure :: NcPut           => Mapping_NcPut
  end type T_Mapping
  

  
contains


  ! ====================================================================
  ! ===
  ! === mapping info
  ! ===
  ! ====================================================================


  subroutine Mapping_Init( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(out)               ::  self
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Mapping_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! no pixels yet:
    self%npix = -999
    ! no mapping info:
    self%nmap     = 0
    nullify( self%map_n )
    nullify( self%map_i0 )
    nullify( self%map_iglb )
    nullify( self%map_ii )
    nullify( self%map_jj )
    nullify( self%map_ww )
    
    ! ok
    status = 0
    
  end subroutine Mapping_Init


  ! ***


  !
  ! Init copy by swap.
  ! - swp  : swapping info for pixel arrays (npix)
  ! - mswp : swapping info for mapping arrays (nmap)
  !

  subroutine Mapping_InitSwap( self, mapping, swp, mswp, status )
  
    use CSO_Comm    , only : csoc
    use CSO_Swapping, only : T_Swapping
  
    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(out)               ::  self
    class(T_Mapping), intent(in)                ::  mapping
    type(T_Swapping), intent(in)                ::  swp
    type(T_Swapping), intent(in)                ::  mswp
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Mapping_InitSwap'
    
    ! --- local ----------------------------------
    
    real, allocatable      ::  sendbuf(:,:)  ! (4,nsend)
    real, allocatable      ::  recvbuf(:,:)  ! (4,nrecv)
    integer                ::  k
    integer                ::  imap
    
    ! --- begin ----------------------------------
    
    ! init base:
    call self%Init( status )
    IF_NOT_OK_RETURN(status=1)

    !~ pixel arrays
    
    ! logical size:
    self%npix = swp%nrecv
    
    ! storage for array with number of overlapping cells per pixel:
    allocate( self%map_n(max(1,swp%nrecv)), source=0, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill "swapped" values, not from "mapping%map_n",
    ! but from the special array 'mswp%nwsend' with number of weights per swapped pixel;
    ! do not use the 'swp%Swap' routine, since input is already "sorted" ..
    call csoc%AllToAllV( mswp%nwsend, swp%sendcounts, &
                         self%map_n , swp%recvcounts, status )
    IF_NOT_OK_RETURN(status=1)
    ! count:
    self%nmap = sum(self%map_n)
    
    ! offset is cumulative sum:
    allocate( self%map_i0(max(1,swp%nrecv)), source=0, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do k = 2, max(1,swp%nrecv)
      self%map_i0(k) = self%map_i0(k-1) + self%map_n(k-1)
    end do

    !~ mapping arrays
    
    ! check ...
    if ( self%nmap /= mswp%nrecv ) then
      write (csol,'("nmap ",i0," while receiving ",i0," mapping elements")') self%nmap, mswp%nrecv; call csoErr
      TRACEBACK; status=1; return
    end if
    
    ! storage for swapping mapping arrays:
    allocate( sendbuf(4,max(1,mswp%nsend)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( recvbuf(4,max(1,mswp%nrecv)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! any to be send?
    if ( mswp%nsend > 0 ) then
      ! loop over elements, in order of sending (grouped per pe):
      do k = 1, mswp%nsend
        ! source element:
        imap = mswp%iimap(k)
        ! fill:
        sendbuf(1,k) = mapping%map_iglb(imap)
        sendbuf(2,k) = mapping%map_ii  (imap)
        sendbuf(3,k) = mapping%map_jj  (imap)
        sendbuf(4,k) = mapping%map_ww  (imap)
      end do ! sending elements
    end if
    
    ! swap values:
    call csoc%AllToAllV( sendbuf, mswp%sendcounts, &
                         recvbuf, mswp%recvcounts, status )
    IF_NOT_OK_RETURN(status=1)

    ! any received?
    if ( mswp%nrecv > 0 ) then
      ! storage:
      allocate( self%map_iglb(mswp%nrecv), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( self%map_ii(mswp%nrecv), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( self%map_jj(mswp%nrecv), stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( self%map_ww(mswp%nrecv), stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! copy:
      self%map_iglb = nint(recvbuf(1,:))
      self%map_ii   = nint(recvbuf(2,:))
      self%map_jj   = nint(recvbuf(3,:))
      self%map_ww   = recvbuf(4,:)
    end if  ! nrecv > 0
    
    ! clear:
    deallocate( sendbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( recvbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Mapping_InitSwap


  ! ***


  subroutine Mapping_Done( self, status )

    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(inout)               ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Mapping_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! footprint mapping info?
    if ( associated(self%map_n) ) then
      deallocate( self%map_n, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%map_i0, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    if ( associated(self%map_iglb) ) then
      deallocate( self%map_iglb, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%map_ii, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%map_jj, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%map_ww, stat=status )
      IF_NOT_OK_RETURN(status=1)
      ! reset counter:
      self%nmap = 0
    end if

    ! ok
    status = 0
    
  end subroutine Mapping_Done
  
  
  ! ***


  subroutine Mapping_GetPixel( self, ipix, status, &
                                            nw, ii, jj, ww )
    
    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(in)                  ::  self
    integer, intent(in)                           ::  ipix
    integer, intent(out)                          ::  status
    
    integer, intent(out), optional                ::  nw
    integer, pointer, optional                    ::  ii(:)   ! (nmap)
    integer, pointer, optional                    ::  jj(:)   ! (nmap)
    real, pointer, optional                       ::  ww(:)   ! (nmap)
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Mapping_GetPixel'
    
    ! --- local ----------------------------------
    
    integer           ::  iw1, iw2
    
    ! --- begin ----------------------------------
    
    ! number of weights:
    if ( present(nw) ) then
      ! check ...
      if ( .not. associated(self%map_n) ) then
        !! testing ...
        !write (csol,'("requested nw but no mapping weights stored yet")'); call csoErr
        !TRACEBACK; status=1; return
        ! no mapping, pixel footprint probably outside (but bounding box partly overlaps ...)
        nw = 0
      else
        ! copy:
        nw = self%map_n(ipix)
      end if
    end if ! nw
    
    ! mapping element:
    if ( any((/present(ii),present(jj),present(ww)/)) ) then
      ! check ...
      if ( .not. associated(self%map_n) ) then
        write (csol,'("requested mapping info but no mapping weights stored yet")'); call csoErr
        TRACEBACK; status=1; return
      end if
      ! index range:
      if ( ipix == 1 ) then
        iw1 = 1
      else
        iw1 = sum(self%map_n(1:ipix-1)) + 1
      end if
      iw2 = iw1-1 + self%map_n(ipix)
      ! copy:
      if ( present(ii) ) ii => self%map_ii(iw1:iw2)
      if ( present(jj) ) jj => self%map_jj(iw1:iw2)
      if ( present(ww) ) ww => self%map_ww(iw1:iw2)
    end if ! ii
    
    ! ok
    status = 0
    
  end subroutine Mapping_GetPixel
  
  
  ! ***


  subroutine Mapping_SetPixel( self, ipix, npix, iglb, ii, jj, ww, status )
    
    use CSO_PArray, only : CSO_PArray_Reshape
    
    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(inout)               ::  self
    integer, intent(in)                           ::  ipix
    integer, intent(in)                           ::  npix
    integer, intent(in)                           ::  iglb
    integer, intent(in)                           ::  ii(:)   ! (nmap)
    integer, intent(in)                           ::  jj(:)   ! (nmap)
    real, intent(in)                              ::  ww(:)   ! (nmap)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Mapping_SetPixel'
    
    ! --- local ----------------------------------
    
    integer             ::  nmap
    integer             ::  ninc
    integer             ::  nnew
    
    ! --- begin ----------------------------------
    
    ! new storage?
    if ( self%nmap == 0 ) then
      ! storage per pixel:
      allocate( self%map_n(npix), source=0, stat=status )
      IF_NOT_OK_RETURN(status=1)
      allocate( self%map_i0(npix), source=0, stat=status )
      IF_NOT_OK_RETURN(status=1)
    end if
    ! check ...
    if ( self%map_n(ipix) > 0 ) then
      write (csol,'("already ",i0," mappings for pixel ",i0," defined")') self%map_n(ipix), ipix; call csoErr
      TRACEBACK; status=1; return
    end if
    ! store number of elements:
    self%map_n(ipix) = size(ww)
    ! cumulative sum:
    if ( ipix > 1 ) self%map_i0(ipix) = self%map_i0(ipix-1) + self%map_n(ipix-1)

    ! new size:
    nmap = self%nmap + size(ww)
    ! need to extend?
    if ( (self%nmap == 0) .or. (nmap > size(self%map_ww)) ) then
      ! size increment, not too modest ...
      ninc = maxval( (/ size(ww), npix, 100 /) )
      ! new size:
      if ( self%nmap == 0 ) then
        nnew = ninc
      else
        nnew = size(self%map_ww) + ninc
      end if
      ! reshape:
      call CSO_PArray_Reshape( self%map_iglb, nnew, status, source=0 )
      IF_NOT_OK_RETURN(status=1)
      call CSO_PArray_Reshape( self%map_ii  , nnew, status, source=0 )
      IF_NOT_OK_RETURN(status=1)
      call CSO_PArray_Reshape( self%map_jj  , nnew, status, source=0 )
      IF_NOT_OK_RETURN(status=1)
      call CSO_PArray_Reshape( self%map_ww  , nnew, status, source=0.0 )
      IF_NOT_OK_RETURN(status=1)
    end if

    ! extension:
    self%map_iglb(self%nmap+1:nmap) = iglb
    self%map_ii  (self%nmap+1:nmap) = ii
    self%map_jj  (self%nmap+1:nmap) = jj
    self%map_ww  (self%nmap+1:nmap) = ww
    ! reset size:
    self%nmap = nmap
    
    ! ok
    status = 0
    
  end subroutine Mapping_SetPixel
  
  
  ! ***


  subroutine Mapping_SetData( self, npix, iglb, nw, ii, jj, ww, status )
    
    use CSO_PArray, only : CSO_PArray_Reshape
    
    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(inout)               ::  self
    integer, intent(in)                           ::  npix
    integer, intent(in)                           ::  iglb(:)  ! (npix)
    integer, intent(in)                           ::  nw(:)    ! (npix)
    integer, intent(in)                           ::  ii(:)    ! (nmap)
    integer, intent(in)                           ::  jj(:)    ! (nmap)
    real, intent(in)                              ::  ww(:)    ! (nmap)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Mapping_SetData'
    
    ! --- local ----------------------------------
    
    integer             ::  ipix
    integer             ::  iw
    integer             ::  imap
    
    ! --- begin ----------------------------------
    
    ! store:
    self%npix = npix
    
    ! storage per pixel:
    allocate( self%map_n(npix), source=0, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! copy:
    self%map_n = nw(1:npix)
    
    ! offset is cumulative sum:
    allocate( self%map_i0(npix), source=0, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! fill:
    do ipix = 2, npix
      self%map_i0(ipix) = self%map_i0(ipix-1) + self%map_n(ipix-1)
    end do
    
    ! total size:
    self%nmap = sum(nw)
    ! reshape arrays if necessary:
    call CSO_PArray_Reshape( self%map_iglb, self%nmap, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%map_ii  , self%nmap, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%map_jj  , self%nmap, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%map_ww  , self%nmap, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! copy global index per pixel;
    ! init index in mapping arrays:
    imap = 0
    ! loop over pixels:
    do ipix = 1, npix
      ! loop over weights:
      do iw = 1, nw(ipix)
        ! increase counter:
        imap = imap + 1
        ! copy global index:
        self%map_iglb(imap) = iglb(ipix)
      end do ! iw
    end do ! ipix
    
    ! copy arrays:
    self%map_ii(1:self%nmap) = ii(1:self%nmap)
    self%map_jj(1:self%nmap) = jj(1:self%nmap)
    self%map_ww(1:self%nmap) = ww(1:self%nmap)
    
    ! ok
    status = 0
    
  end subroutine Mapping_SetData


  ! ***


  subroutine Mapping_Collect( self, npix, iout_all, nout, iout_glb, status )

    use CSO_Comm        , only : csoc
    
    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(inout)               ::  self
    integer, intent(in)                           ::  npix
    ! npix_all is sum(npix) over all domains,
    ! mapping from 1:npix_all to 1:nout:
    integer, intent(in)                           ::  iout_all(:)  ! (npix_all)
    ! number of unique pixels used at some domain,
    ! equal or less than nglb in case some pixels are outside domain:
    integer, intent(in)                           ::  nout
    ! mapping from 1:nglb to 1:nout,
    ! for example used to write glb_lon/etc arrays
    integer, intent(in)                           ::  iout_glb(:)  ! (nglb)

    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Mapping_Collect'
    
    ! --- local ----------------------------------
    
    integer, allocatable      ::  map_n_all(:)     ! (npix_all)
    integer, allocatable      ::  map_iglb_all(:)  ! (nmap_all)
    integer, allocatable      ::  map_ii_all(:)    ! (nmap_all)
    integer, allocatable      ::  map_jj_all(:)    ! (nmap_all)
    real, allocatable         ::  map_ww_all(:)    ! (nmap_all)
    integer, allocatable      ::  i0(:)     ! (nout)
    integer, allocatable      ::  nk(:)     ! (nout)
    integer                   ::  npix_all
    integer                   ::  ipix_all
    integer                   ::  iglb
    integer                   ::  iout
    integer                   ::  k
    integer                   ::  i
    
    ! --- begin ----------------------------------
    
    !! testing ..
    !write (csol,'(a,": start")') rname; call csoPr

    ! total number of mapping elements:
    call csoc%ParInfo( self%nmap, status, ntot=self%nmap_all )
    IF_NOT_OK_RETURN(status=1)
    
    ! any mappings at all?
    if ( self%nmap_all > 0 ) then
    
      ! storage for collection arrays:
      if ( csoc%root ) then
        ! size:
        npix_all = size(iout_all)
        ! number of mapping elements per pixel on local domain:
        allocate( map_n_all(npix_all), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! all mapping arrays:
        allocate( map_iglb_all(self%nmap_all), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( map_ii_all(self%nmap_all), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( map_jj_all(self%nmap_all), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( map_ww_all(self%nmap_all), stat=status )
        IF_NOT_OK_RETURN(status=1)
      else
        ! dummy ...
        allocate( map_n_all(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! dummy ...
        allocate( map_iglb_all(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( map_ii_all(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( map_jj_all(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( map_ww_all(1), stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if
      
      ! no mappings? probably pixel footprints outside domain (but selected anyway ..)
      if ( .not. associated(self%map_n) ) then
        ! dummy storage ...
        allocate( self%map_n(max(1,npix)), source=0, stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%map_i0(max(1,npix)), source=0, stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%map_iglb(1), source=-999, stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%map_ii(1), source=-999, stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%map_jj(1), source=-999, stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%map_ww(1), source=0.0, stat=status )
        IF_NOT_OK_RETURN(status=1)
      end if ! no mappings on this domain

      ! gather on root:
      call csoc%GatherV( self%map_n, map_n_all, status, nloc=npix )
      IF_NOT_OK_RETURN(status=1)
      call csoc%GatherV( self%map_iglb, map_iglb_all, status, nloc=self%nmap )
      IF_NOT_OK_RETURN(status=1)
      call csoc%GatherV( self%map_ii, map_ii_all, status, nloc=self%nmap )
      IF_NOT_OK_RETURN(status=1)
      call csoc%GatherV( self%map_jj, map_jj_all, status, nloc=self%nmap )
      IF_NOT_OK_RETURN(status=1)
      call csoc%GatherV( self%map_ww, map_ww_all, status, nloc=self%nmap )
      IF_NOT_OK_RETURN(status=1)

      ! collect into output arrays on root:
      if ( csoc%root ) then
        ! output arrays:
        allocate( self%map_n_out(nout), source=0, stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%map_ii_out(self%nmap_all), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%map_jj_out(self%nmap_all), stat=status )
        IF_NOT_OK_RETURN(status=1)
        allocate( self%map_ww_out(self%nmap_all), stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over pixels:
        do ipix_all = 1, npix_all
          ! output entry:
          iout = iout_all(ipix_all)
          ! increase:
          self%map_n_out(iout) = self%map_n_out(iout) + map_n_all(ipix_all)
        end do
        ! offsets:
        allocate( i0(nout), stat=status )
        IF_NOT_OK_RETURN(status=1)
        i0(1) = 0
        do i = 2, nout
          i0(i) = i0(i-1) + self%map_n_out(i-1)
        end do
        ! counters for mappings already stored, init to zero:
        allocate( nk(nout), source=0, stat=status )
        IF_NOT_OK_RETURN(status=1)
        ! loop over mapping elements:
        do k = 1, self%nmap_all
          ! corresponding global pixel:
          iglb = map_iglb_all(k)
          ! corresponding output entry:
          iout = iout_glb(iglb)
          ! check ...
          if ( iout < 1 ) then
            write (csol,'("invalid output entry ",i0," for iglb ",i0)') iout, iglb; call csoErr
            TRACEBACK; status=1; return
          end if
          ! increase counter:
          nk(iout) = nk(iout) + 1
          ! target index:
          i = i0(iout) + nk(iout)
          ! re-order:
          self%map_ii_out(i) = map_ii_all(k)
          self%map_jj_out(i) = map_jj_all(k)
          self%map_ww_out(i) = map_ww_all(k)
        end do ! mapping elements
        ! clear:
        deallocate( i0, stat=status )
        IF_NOT_OK_RETURN(status=1)
        deallocate( nk, stat=status )
        IF_NOT_OK_RETURN(status=1)

      end if ! root

      ! clear
      deallocate( map_n_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( map_iglb_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( map_ii_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( map_jj_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( map_ww_all, stat=status )
      IF_NOT_OK_RETURN(status=1)
      
    end if  ! nmap_all > 0
    
    !! testing ..
    !write (csol,'(a,": end")') rname; call csoPr

    ! ok
    status = 0
    
  end subroutine Mapping_Collect


  ! ***


  subroutine Mapping_NcDef( self, ncid, dimid_pixel, status )

    use NetCDF, only : NF90_INT, NF90_FLOAT
    use NetCDF, only : NF90_Def_Dim
    use NetCDF, only : NF90_Def_Var
    use NetCDF, only : NF90_Put_Att

    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(inout)               ::  self
    integer, intent(in)                           ::  ncid
    integer, intent(in)                           ::  dimid_pixel
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Mapping_NcDef'
    
    ! --- local ----------------------------------

    integer                   ::  dimid_mapping
    integer                   ::  varid
    
    ! --- begin ----------------------------------
    
    ! any mappings defined?
    if ( self%nmap_all > 0 ) then
      ! define dimension:
      status = NF90_Def_Dim( ncid, 'mapping', self%nmap_all, dimid_mapping )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! define variable:
      status = NF90_Def_Var( ncid, 'mapping_n', NF90_INT, (/dimid_pixel/), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attributes:
      status = NF90_Put_Att( ncid, varid, 'long_name', 'number of mapping elements' )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( ncid, varid, 'units', '1' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! store:
      self%varid_n_out = varid
      
      ! define variable:
      status = NF90_Def_Var( ncid, 'mapping_i', NF90_INT, (/dimid_mapping/), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attributes:
      status = NF90_Put_Att( ncid, varid, 'long_name', 'mapping source cell longitude index' )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( ncid, varid, 'units', '1' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! store:
      self%varid_ii_out = varid

      ! define variable:
      status = NF90_Def_Var( ncid, 'mapping_j', NF90_INT, (/dimid_mapping/), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attributes:
      status = NF90_Put_Att( ncid, varid, 'long_name', 'mapping source cell latitude index' )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( ncid, varid, 'units', '1' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! store:
      self%varid_jj_out = varid

      ! define variable:
      status = NF90_Def_Var( ncid, 'mapping_w', NF90_FLOAT, (/dimid_mapping/), varid )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! attributes:
      status = NF90_Put_Att( ncid, varid, 'long_name', 'mapping source cell overlap' )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Att( ncid, varid, 'units', 'm2' )
      IF_NF90_NOT_OK_RETURN(status=1)
      ! store:
      self%varid_ww_out = varid

    end if ! mapping
    
    ! ok
    status = 0
    
  end subroutine Mapping_NcDef


  ! ***


  subroutine Mapping_NcPut( self, ncid, status )

    use NetCDF, only : NF90_Put_Var

    ! --- in/out ---------------------------------
    
    class(T_Mapping), intent(inout)               ::  self
    integer, intent(in)                           ::  ncid
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Mapping_NcPut'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
        
    ! any mappings?
    if ( self%nmap_all > 0 ) then

      ! write variables:
      status = NF90_Put_Var( ncid, self%varid_n_out, self%map_n_out )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Var( ncid, self%varid_ii_out, self%map_ii_out )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Var( ncid, self%varid_jj_out, self%map_jj_out )
      IF_NF90_NOT_OK_RETURN(status=1)
      status = NF90_Put_Var( ncid, self%varid_ww_out, self%map_ww_out )
      IF_NF90_NOT_OK_RETURN(status=1)

      ! clear output arrays:
      deallocate( self%map_n_out, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%map_ii_out, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%map_jj_out, stat=status )
      IF_NOT_OK_RETURN(status=1)
      deallocate( self%map_ww_out, stat=status )
      IF_NOT_OK_RETURN(status=1)

    end if ! mappings defined
    
    ! ok
    status = 0
    
  end subroutine Mapping_NcPut



end module CSO_Mapping

