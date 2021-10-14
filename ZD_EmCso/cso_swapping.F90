!###############################################################################
!
! Exchange pixels between domains
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

module CSO_Swapping

  use CSO_Logging, only : csol, csoPr, csoErr

  implicit none
  
  
  ! --- in/out ----------------------------
  
  private

  public  ::  T_Swapping
  
  
  ! --- const ------------------------------
  
  character(len=*), parameter   ::  mname = 'CSO_Swapping'


  ! --- types --------------------------------
  
  ! Info on swapping between domains
  type T_Swapping
    ! type: 'pix', 'map'
    character(len=3)          ::  atype
    ! total number of values to send/recv:
    integer                   ::  nsend
    integer                   ::  nrecv
    ! per targer/source:
    integer, allocatable      ::  sendcounts(:)  ! (0:npes-1)
    integer, allocatable      ::  recvcounts(:)  ! (0:npes-1)
    ! pixel indices in order of sending:
    integer, pointer          ::  iipix(:)       ! (nsend)
    integer, pointer          ::  iimap(:)       ! (nsend)
    ! extra for atype 'pix': number of selected weights:
    integer, pointer          ::  nwsend(:)       ! (nsend)
    !
  contains
    procedure :: Init            => Swapping_Init
    procedure :: Done            => Swapping_Done
    procedure ::                    Swapping_Swap_i_1d
    procedure ::                    Swapping_Swap_r_2d
    generic   :: Swap            => Swapping_Swap_i_1d, &
                                    Swapping_Swap_r_2d
  end type T_Swapping
  

  
contains


  ! ====================================================================
  ! ===
  ! === exchange info
  ! ===
  ! ====================================================================

  ! Init swapping:
  ! - atype  : defines array type to be swapped:
  !     'pix'   : pixel arrays, local size (npix)
  !     'map'   : mapping weights, local size (nmap)

  subroutine Swapping_Init( self, atype, nw, ii, jj, doms, doms_f, status )
  
    use CSO_Comm   , only : csoc
    use CSO_Domains, only : T_CSO_Domains
    use CSO_PArray , only : CSO_PArray_Init, CSO_PArray_Done, CSO_PArray_Reshape

    ! --- in/out ---------------------------------
    
    class(T_Swapping), intent(out)              ::  self
    character(len=*), intent(in)                ::  atype
    integer, intent(in)                         ::  nw(:)  ! (npix)
    integer, intent(in)                         ::  ii(:)  ! (nmap)
    integer, intent(in)                         ::  jj(:)  ! (nmap)
    type(T_CSO_Domains), intent(in)             ::  doms
    type(T_CSO_Domains), intent(in)             ::  doms_f
    integer, intent(out)                        ::  status

    ! --- const ----------------------------------
    
    character(len=*), parameter   :: rname = mname//'/Swapping_Init'
    
    ! --- local ----------------------------------

    integer                   ::  npix
    integer                   ::  nsend_max
    integer                   ::  displ
    integer                   ::  ipix
    integer                   ::  iw
    integer                   ::  imap
    integer                   ::  id
    integer                   ::  id_f
    logical                   ::  assigned_to_id
    integer                   ::  iwsend
    
    ! --- begin ----------------------------------
    
    ! store:
    self%atype = trim(atype)
    ! per pixel or per mapping?
    select case ( self%atype )
      case ( 'pix', 'map' )
        ! ok
      case default
        write (csol,'("unsupported atype `",a,"`")') trim(self%atype); call csoErr
        TRACEBACK; status=1; return
    end select
    
    !! testing ...
    !if ( self%atype == 'pix' ) then
    !  print *, '---- define swap of pixel data ...'
    !else
    !  print *, '---- define swap of mapping data ...'
    !end if
    
    ! storage for sending counts and displacements:
    allocate( self%sendcounts(0:csoc%npes-1), source=0, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! storage for receiving countns and displacements:
    allocate( self%recvcounts(0:csoc%npes-1), source=0, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of pixels:
    npix = size(nw)
    
    ! init send buffer and pixel-to-send mapping,
    ! actual size is dynamically extended when needed:
    call CSO_PArray_Init( self%iipix, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%iimap, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Init( self%nwsend, status )
    IF_NOT_OK_RETURN(status=1)
    
    ! number of pixels send to other domains;
    ! some pixels might be send to multiple domains,
    ! so more than 'npix' needed ...
    nsend_max = max( 100, int( 1.1 * npix ) )
    ! initial storage:
    call CSO_PArray_Reshape( self%iipix, nsend_max, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%iimap, nsend_max, status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%nwsend, nsend_max, status )
    IF_NOT_OK_RETURN(status=1)

    ! init counter for sorted obs on this proc:
    self%nsend = 0
    iwsend = 0
    ! loop over target procs:
    do id = 0, csoc%npes-1
      !! testing ..
      !print *, '  target pe ', id
      ! init counter:
      self%sendcounts(id) = 0
      ! offset in sorted array:
      displ = self%nsend
      ! init mapping index:
      imap = 0
      ! loop over pixels:
      do ipix = 1, npix
        ! not assigned yet:
        assigned_to_id = .false.
        ! any overlap with grid cells?
        if ( nw(ipix) == 0 ) then
          !! testing ...
          !print *, '  -- no w ipix ', ipix
        else
          ! loop over overlapping grid cells:
          do iw = 1, nw(ipix)
            ! increase index:
            imap = imap + 1
            ! find pe of target domain:
            call doms_f%Find( (/ii(imap),jj(imap)/), id_f, status )
            IF_NOT_OK_RETURN(status=1)
            ! current target ?
            if ( id_f == id ) then
              ! one pixel per processor, and already assigned? then skip:
              if ( (self%atype == 'pix') .and. assigned_to_id ) cycle
              ! increase counters:
              self%nsend          = self%nsend + 1
              self%sendcounts(id) = self%sendcounts(id) + 1
              ! extend storage?
              if ( self%nsend > nsend_max ) then
                ! extend:
                nsend_max = nsend_max + 100
                ! reshape:
                call CSO_PArray_Reshape( self%iipix, nsend_max, status )
                IF_NOT_OK_RETURN(status=1)
                call CSO_PArray_Reshape( self%iimap, nsend_max, status )
                IF_NOT_OK_RETURN(status=1)
                call CSO_PArray_Reshape( self%nwsend, nsend_max, status )
                IF_NOT_OK_RETURN(status=1)
              end if
              ! store source indices:
              self%iipix(self%nsend) = ipix
              self%iimap(self%nsend) = imap
              ! init number of weights, only used for atype 'pix' actually ...
              if ( .not. assigned_to_id ) then
                iwsend = iwsend + 1
                self%nwsend(iwsend) = 1
              else
                self%nwsend(iwsend) = self%nwsend(iwsend) + 1
              end if
              !! testing ...
              !print *, '  -- swap ipix ', ipix, ' to domain ', id, ' ipix ', self%nsend
              !if ( id == 0 ) then
              !  print *, '  -- swap ipix ', ipix, ' iw ', iw, ' to domain ', id, ' ipix ', self%nsend, &
              !                ';', self%nwsend(iwsend), sum(self%nwsend(1:iwsend))
              !end if
              !! testing ..
              !to_npe(ipix) = to_npe(ipix) + 1
              ! reset flag:
              assigned_to_id = .true.
            end if  ! current target
          end do ! iw
        end if ! nw > 0
      end do  ! ipix
    end do  ! target pe
    
    ! total number to be send:
    self%nsend = sum(self%sendcounts)
    ! reset to exact shape:
    call CSO_PArray_Reshape( self%iipix, max(1,self%nsend), status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%iimap, max(1,self%nsend), status )
    IF_NOT_OK_RETURN(status=1)
    call CSO_PArray_Reshape( self%nwsend, max(1,iwsend), status )
    IF_NOT_OK_RETURN(status=1)
  
    ! breakup "sendcounts" in single number per target pe,
    ! broadcast to all pe, and store the number received from each source pe:
    call csoc%AllToAll( self%sendcounts, 1, &
                        self%recvcounts, 1, status )
    IF_NOT_OK_RETURN(status=1)   

    ! total number received:
    self%nrecv = sum(self%recvcounts)
    
    !! testing ...
    !print *, '-- iwsend     = ', iwsend
    !print *, '-- sendcounts = ', self%sendcounts, ';', self%nsend
    !print *, '-- recvcounts = ', self%recvcounts, ';', self%nrecv
    ! testing ...
    do id = 0, csoc%npes-1
      write (csol,'(a,": send ",i0," pixels to pe ",i0,", receive ",i0)') &
                            rname, self%sendcounts(id), id, self%recvcounts(id); call csoPr
    end do
  
    ! ok
    status = 0
    
  end subroutine Swapping_Init


  ! ***


  subroutine Swapping_Done( self, status )

    ! --- in/out ---------------------------------
    
    class(T_Swapping), intent(inout)              ::  self
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Swapping_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------

    ! clear:
    deallocate( self%sendcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%recvcounts, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( self%iipix, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( self%iimap, stat=status )
    IF_NOT_OK_RETURN(status=1)
    ! clear:
    deallocate( self%nwsend, stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! ok
    status = 0
    
  end subroutine Swapping_Done


  ! ***


  subroutine Swapping_Swap_i_1d( self, values, recvbuf, status )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------
    
    class(T_Swapping), intent(in)                 ::  self
    integer, intent(in)                           ::  values(:)    ! (npix)
    integer, intent(out)                          ::  recvbuf(:)   ! (npix_f)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Swapping_Swap_i_1d'
    
    ! --- local ----------------------------------
    
    integer, allocatable      ::  sendbuf(:)
    integer                   ::  k
    
    ! --- begin ----------------------------------

    ! storage for data to be send and received:
    allocate( sendbuf(max(1,self%nsend)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! any data to be send?
    if ( self%nsend > 0 ) then
      select case ( self%atype )
        !~ pixel arrays
        case ( 'pix' )
          ! copy values into send buffer in order of targe pe:
          do k = 1, self%nsend
            sendbuf(k) = values(self%iipix(k))
          end do
        !~ mapping arrays
        case ( 'map' )
          ! copy values into send buffer in order of targe pe:
          do k = 1, self%nsend
            sendbuf(k) = values(self%iimap(k))
          end do
        !~
        case default
          write (csol,'("unsupported atype `",a,"`")') trim(self%atype); call csoErr
          TRACEBACK; status=1; return
      end select
    end if ! nsend > 0
  
    ! swap values:
    call csoc%AllToAllV( sendbuf, self%sendcounts, &
                         recvbuf, self%recvcounts, &
                         status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Swapping_Swap_i_1d


  ! ***


  subroutine Swapping_Swap_r_2d( self, values, recvbuf, status )
  
    use CSO_Comm, only : csoc

    ! --- in/out ---------------------------------
    
    class(T_Swapping), intent(in)                 ::  self
    real, intent(in)                              ::  values(:,:)    ! (n,npix)
    real, intent(out)                             ::  recvbuf(:,:)   ! (n,npix_f)
    integer, intent(out)                          ::  status
  
    ! --- const ----------------------------------
    
    character(len=*), parameter  ::  rname = mname//'/Swapping_Swap_i_2d'
    
    ! --- local ----------------------------------
    
    integer                   ::  n
    real, allocatable         ::  sendbuf(:,:)
    integer                   ::  k
    
    ! --- begin ----------------------------------
    
    ! number of elements per pixel:
    n = size(values,1)

    ! storage for data to be send and received:
    allocate( sendbuf(n,max(1,self%nsend)), stat=status )
    IF_NOT_OK_RETURN(status=1)
    
    ! any data to be send?
    if ( self%nsend > 0 ) then
      select case ( self%atype )
        !~ pixel arrays
        case ( 'pix' )
          ! copy values into send buffer in order of targe pe:
          do k = 1, self%nsend
            sendbuf(:,k) = values(:,self%iipix(k))
          end do
        !~ mapping arrays
        case ( 'map' )
          ! copy values into send buffer in order of targe pe:
          do k = 1, self%nsend
            sendbuf(:,k) = values(:,self%iimap(k))
          end do
        !~
        case default
          write (csol,'("unsupported atype `",a,"`")') trim(self%atype); call csoErr
          TRACEBACK; status=1; return
      end select
    end if ! nsend > 0
  
    ! swap values:
    call csoc%AllToAllV( sendbuf, self%sendcounts, &
                         recvbuf, self%recvcounts, &
                         status )
    IF_NOT_OK_RETURN(status=1)
    
    ! clear:
    deallocate( sendbuf, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
    
  end subroutine Swapping_Swap_r_2d


end module CSO_Swapping

