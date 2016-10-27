!######################################################################
!
! GO_SampStat  -  incremental sample statistics
!
!
!
! USAGE
!
! BACKGROUND
!
!          1  n       
!   m[n] = - sum x[i] 
!          n i=1      
!
!          1 
!        = - { (n-1) m[n-1] + x[n] }
!          n 
!
!        = m[n-1] + { x[n] - m[n-1] } /n
!                   \__________________/
!                           d[n]
!
!           1  n       
!   ms[n] = - sum x[i]^2
!           n i=1      
!
!          1 
!        = - { (n-1) ms[n-1] + x[n]^2 }
!          n 
!
!        = ms[n-1] + { x[n]^2 - ms[n-1] } /n
!                     \____________________/
!                             ds[n]
!
!
!           1   n
!   P[n] = --- sum ( x[i] - m[n] ) ( x[i] - m[n] )'
!          n-1 i=1
!
!           1   n
!        = --- sum ( x[i] - m[n-1] - d[n] ) ( x[i] - m[n-1] - d[n] )'
!          n-1 i=1
!
!           1   n
!        = --- sum (x[i]-m[n-1])(x[i]-m[n-1])'
!          n-1 i=1
!
!             1   n
!          - --- sum { (x[i]-m[n-1])d[n]' + d[n](x[i]-m[n-1])' }
!            n-1 i=1
!
!             1   n
!          + --- sum d[n]d[n]'
!            n-1 i=1
!
!           1                 
!        = --- { (n-2)P[n-1]  +  n^2 d[n]d[n]' }
!          n-1               
!
!             n  
!          - --- { (m[n] - m[n-1])d[n]' + d[n](m[n] - m[n-1])' }
!            n-1 
!
!             n  
!          + --- d[n]d[n]'
!            n-1 
!
!          n-2            n^2+n
!        = --- P[n-1]  +  ----- d[n]d[n]'
!          n-1             n-1
!
!             n  
!          - --- { (m[n] - m[n-1])d[n]' + d[n](m[n] - m[n-1])' }
!            n-1 
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!######################################################################

module GO_SampStat

  use GO_Print, only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  T_SampStat, T_SampStat_c
  public  ::  T_SampCovr_3D_1D
  public  ::  T_SampCovr_2D_2D
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'SampStat'
  
  ! max rank:
  integer, parameter    ::  max_rank = 3

  ! --- types ----------------------------------------
  
  ! * statistics per element
  
  type T_SampStat
    ! rank:
    integer                       ::  rank
    ! problem size:
    integer                       ::  shp(max_rank)
    ! sample counters:
    ! (per element to prepare for masked samples)
    integer, allocatable          ::  n(:,:,:)
    ! data:
    real, allocatable             ::  m_n (:,:,:)  ! mean
    real, allocatable             ::  ms_n(:,:,:)  ! mean square
    real, allocatable             ::  v_n (:,:,:)  ! variance
    ! temporary:
    real, allocatable             ::  m_nm1 (:,:,:)
    real, allocatable             ::  ms_nm1(:,:,:)
    real, allocatable             ::  v_nm1 (:,:,:)
    real, allocatable             ::  x     (:,:,:)
    real, allocatable             ::  d     (:,:,:)
    real, allocatable             ::  dm    (:,:,:)
    !
  contains
    !
    procedure   ::  Init      => SampStat_Init
    procedure   ::  Done      => SampStat_Done
    !
    procedure   ::  Get2d     => SampStat_Get2d
    procedure   ::  Get3d     => SampStat_Get3d
    !
    procedure   ::               SampStat_AddSample_2d
    procedure   ::               SampStat_AddSample_3d
    generic     ::  AddSample => SampStat_AddSample_2d, &
                                 SampStat_AddSample_3d
  end type T_SampStat
  
  ! ~
  
  type T_SampStat_c
    ! rank:
    integer                       ::  rank
    ! problem size:
    integer                       ::  shp(max_rank)
    ! sample counters:
    ! (per element to prepare for masked samples)
    integer, allocatable          ::  n(:,:,:)
    ! data:
    complex, allocatable          ::  m_n (:,:,:)  ! mean
    complex, allocatable          ::  ms_n(:,:,:)  ! mean square
    complex, allocatable          ::  v_n (:,:,:)  ! variance
    ! temporary:
    complex, allocatable          ::  m_nm1 (:,:,:)
    complex, allocatable          ::  ms_nm1(:,:,:)
    complex, allocatable          ::  v_nm1 (:,:,:)
    complex, allocatable          ::  x     (:,:,:)
    complex, allocatable          ::  d     (:,:,:)
    complex, allocatable          ::  dm    (:,:,:)
    !
  contains
    !
    procedure   ::  Init      => SampStat_c_Init
    procedure   ::  Done      => SampStat_c_Done
    !
    procedure   ::  Get2d     => SampStat_c_Get2d
    procedure   ::  Get3d     => SampStat_c_Get3d
    !
    procedure   ::               SampStat_c_AddSample_2d
    procedure   ::               SampStat_c_AddSample_3d
    generic     ::  AddSample => SampStat_c_AddSample_2d, &
                                 SampStat_c_AddSample_3d
  end type T_SampStat_c
  
  ! * covariance of field 'x(:,:,:)' with elements 'y(:)'
  
  type T_SampCovr_3D_1D
    ! rank:
    integer                       ::  xrank
    ! problem size:
    integer                       ::  xshp(max_rank)
    integer                       ::  yshp
    ! sample counter, do not prepare for masked:
    integer                       ::  n
    ! data for x:
    real, allocatable             ::  xm_n (:,:,:)    ! x mean
    real, allocatable             ::  xv_n (:,:,:)    ! x variance
    ! data for y:
    real, allocatable             ::  ym_n (:)        ! y mean
    real, allocatable             ::  yv_n (:)        ! y variance
    ! data for x and y:
    real, allocatable             ::  cv_n (:,:,:,:)  ! xy covariance
    ! temporary:
    real, allocatable             ::  xm_nm1 (:,:,:)
    real, allocatable             ::  xv_nm1 (:,:,:)
    real, allocatable             ::  ym_nm1 (:)
    real, allocatable             ::  yv_nm1 (:)
    real, allocatable             ::  cv_nm1 (:,:,:,:)
    real, allocatable             ::  xd     (:,:,:)
    real, allocatable             ::  xdm    (:,:,:)
    real                          ::  yd 
    real                          ::  ydm
    !
  contains
    !
    procedure   ::  Init      => SampCovr_3D_1D_Init
    procedure   ::  Done      => SampCovr_3D_1D_Done
    !
    procedure   ::  GetPoint  => SampCovr_3D_1D_GetPoint
    procedure   ::  Get3d     => SampCovr_3D_1D_Get3D
    !
    procedure   ::               SampCovr_3D_1D_AddSample_3d
    generic     ::  AddSample => SampCovr_3D_1D_AddSample_3d
  end type T_SampCovr_3D_1D
  
  ! *
  
  type T_SampCovr_2D_2D
    ! ranks:
    integer                       ::  xrank
    integer                       ::  yrank
    ! problem size:
    integer                       ::  xshp(max_rank)
    integer                       ::  yshp(max_rank)
    ! sample counter, do not prepare for masked:
    integer                       ::  n
    ! data for x:
    real, allocatable             ::  xm_n (:,:)      ! x mean
    real, allocatable             ::  xv_n (:,:)      ! x variance
    ! data for y:
    real, allocatable             ::  ym_n (:,:)      ! y mean
    real, allocatable             ::  yv_n (:,:)      ! y variance
    ! data for x and y:
    real, allocatable             ::  cv_n (:,:,:,:)  ! xy covariance
    ! temporary:
    real, allocatable             ::  xm_nm1 (:,:)
    real, allocatable             ::  xv_nm1 (:,:)
    real, allocatable             ::  ym_nm1 (:,:)
    real, allocatable             ::  yv_nm1 (:,:)
    real, allocatable             ::  cv_nm1 (:,:,:,:)
    real, allocatable             ::  xd     (:,:)
    real, allocatable             ::  xdm    (:,:)
    real                          ::  yd 
    real                          ::  ydm
    !
  contains
    !
    procedure   ::  Init      => SampCovr_2D_2D_Init
    procedure   ::  Done      => SampCovr_2D_2D_Done
    !
    procedure   ::  GetXY     => SampCovr_2D_2D_GetXY
    !
    procedure   ::               SampCovr_2D_2D_AddSample_2d
    generic     ::  AddSample => SampCovr_2D_2D_AddSample_2d
  end type T_SampCovr_2D_2D
  
  

contains


  ! ********************************************************************
  ! ***
  ! *** statistics per element
  ! ***
  ! ********************************************************************


  subroutine SampStat_Init( self, shp, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat), intent(out)      ::  self
    integer, intent(in)                 ::  shp(:)
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! count:
    self%rank = size(shp)
    ! check ...
    if ( self%rank > max_rank ) then
      write (*,'("ERROR - size of shp ",i6," exceeds max rank ",i6)') self%rank, max_rank
      TRACEBACK; status=1; return
    end if
    
    ! store:
    self%shp = 1
    self%shp(1:self%rank) = shp
    
    ! storage for end results:
    allocate( self%n     (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%m_n   (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%ms_n  (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%v_n   (self%shp(1),self%shp(2),self%shp(3)) )
    ! work arrays:
    allocate( self%m_nm1 (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%ms_nm1(self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%v_nm1 (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%x     (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%d     (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%dm    (self%shp(1),self%shp(2),self%shp(3)) )
    
    ! no samples yet:
    self%n    = 0
    self%m_n  = 0.0
    self%ms_n = 0.0
    self%v_n  = 0.0

    ! ok
    status = 0

  end subroutine SampStat_Init


  ! ***
  
  
  subroutine SampStat_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat), intent(inout)   ::  self
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! dummy:
    self%rank = -999
    self%shp  = -999
    
    ! clear:
    deallocate( self%n      )
    deallocate( self%m_n    )
    deallocate( self%ms_n   )
    deallocate( self%v_n    )
    ! clear work arrays:
    deallocate( self%m_nm1  )
    deallocate( self%ms_nm1 )
    deallocate( self%v_nm1  )
    deallocate( self%x      )
    deallocate( self%d      )
    deallocate( self%dm     )

    ! ok
    status = 0

  end subroutine SampStat_Done


  ! ***
  
  
  subroutine SampStat_Get2d( self, status, n, mean, var, stdv )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat), intent(inout)   ::  self
    integer, intent(out)               ::  status
    
    real, intent(out), optional        ::  n   (:,:)
    real, intent(out), optional        ::  mean(:,:)
    real, intent(out), optional        ::  var (:,:)
    real, intent(out), optional        ::  stdv(:,:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_Get2d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! copy:
    if ( present(n   ) ) n    = self%n  (:,:,1)
    if ( present(mean) ) mean = self%m_n(:,:,1)
    if ( present(var ) ) var  = self%v_n(:,:,1)
    
    ! process:
    if ( present(stdv) ) then
      ! compute from variance if possible:
      where ( self%n(:,:,1) > 1 )
        stdv = sqrt( self%v_n(:,:,1) )
      else where
        stdv = 0.0
      end where
    end if

    ! ok
    status = 0

  end subroutine SampStat_Get2d


  ! ***
  
  
  subroutine SampStat_Get3d( self, status, n, mean, var, stdv )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat), intent(inout)   ::  self
    integer, intent(out)               ::  status
    
    real, intent(out), optional        ::  n   (:,:,:)
    real, intent(out), optional        ::  mean(:,:,:)
    real, intent(out), optional        ::  var (:,:,:)
    real, intent(out), optional        ::  stdv(:,:,:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_Get3d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! copy:
    if ( present(n   ) ) n    = self%n  (:,:,:)
    if ( present(mean) ) mean = self%m_n(:,:,:)
    if ( present(var ) ) var  = self%v_n(:,:,:)
    
    ! process:
    if ( present(stdv) ) then
      ! compute from variance if possible:
      where ( self%n(:,:,:) > 1 )
        stdv = sqrt( self%v_n(:,:,:) )
      else where
        stdv = 0.0
      end where
    end if

    ! ok
    status = 0

  end subroutine SampStat_Get3d
  
  
  ! ***
  
  
  subroutine SampStat_AddSample_2d( self, x_n, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat), intent(inout)        ::  self
    real, intent(in)                        ::  x_n(:,:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_AddSample_2d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! copy:
    self%x(:,:,1) = x_n
    ! add:
    call SampStat_AddSample_3d( self, self%x, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SampStat_AddSample_2d
  
  
  ! ***
  
  
  subroutine SampStat_AddSample_3d( self, x_n, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat), intent(inout)        ::  self
    real, intent(in)                        ::  x_n(:,:,:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_AddSample_3d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    !$OMP parallel
    !$OMP   workshare
    
    ! copy current values:
    self%m_nm1  = self%m_n
    self%ms_nm1 = self%ms_n
    self%v_nm1  = self%v_n
    
    ! increase sample counters :
    self%n = self%n + 1

    ! update means and mean squares:
    self%m_n  = ( real(self%n-1)*self%m_nm1  + x_n    ) / real(self%n)
    self%ms_n = ( real(self%n-1)*self%ms_nm1 + x_n**2 ) / real(self%n)

    ! difference fields:
    self%d  = (      x_n - self%m_nm1 ) / real(self%n)
    self%dm =   self%m_n - self%m_nm1
    ! update variance, only if at least 2 samples present:
    where ( self%n > 1 )
      self%v_n =              real(self%n-2)/real(self%n-1) * self%v_nm1 + &
                 real(self%n)*real(self%n+1)/real(self%n-1) * self%d * self%d  - &
                              real(self%n*1)/real(self%n-1) * ( self%dm * self%d  + self%d * self%dm )
    end where
    
    !$OMP   end workshare
    !$OMP end parallel

    ! ok
    status = 0

  end subroutine SampStat_AddSample_3d


  ! ********************************************************************
  ! ***
  ! *** statistics per element
  ! ***
  ! ********************************************************************


  subroutine SampStat_c_Init( self, shp, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat_c), intent(out)    ::  self
    integer, intent(in)                 ::  shp(:)
    integer, intent(out)                ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_c_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! count:
    self%rank = size(shp)
    ! check ...
    if ( self%rank > max_rank ) then
      write (*,'("ERROR - size of shp ",i6," exceeds max rank ",i6)') self%rank, max_rank
      TRACEBACK; status=1; return
    end if
    
    ! store:
    self%shp = 1
    self%shp(1:self%rank) = shp
    
    ! storage for end results:
    allocate( self%n     (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%m_n   (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%ms_n  (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%v_n   (self%shp(1),self%shp(2),self%shp(3)) )
    ! work arrays:
    allocate( self%m_nm1 (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%ms_nm1(self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%v_nm1 (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%x     (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%d     (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%dm    (self%shp(1),self%shp(2),self%shp(3)) )
    
    ! no samples yet:
    self%n    = 0
    self%m_n  = 0.0
    self%ms_n = 0.0
    self%v_n  = 0.0

    ! ok
    status = 0

  end subroutine SampStat_c_Init


  ! ***
  
  
  subroutine SampStat_c_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat_c), intent(inout)      ::  self
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_c_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! dummy:
    self%rank = -999
    self%shp  = -999
    
    ! clear:
    deallocate( self%n      )
    deallocate( self%m_n    )
    deallocate( self%ms_n   )
    deallocate( self%v_n    )
    ! clear work arrays:
    deallocate( self%m_nm1  )
    deallocate( self%ms_nm1 )
    deallocate( self%v_nm1  )
    deallocate( self%x      )
    deallocate( self%d      )
    deallocate( self%dm     )

    ! ok
    status = 0

  end subroutine SampStat_c_Done


  ! ***
  
  
  subroutine SampStat_c_Get2d( self, status, n, mean, var, stdv )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat_c), intent(inout)    ::  self
    integer, intent(out)                  ::  status
    
    complex, intent(out), optional        ::  n   (:,:)
    complex, intent(out), optional        ::  mean(:,:)
    complex, intent(out), optional        ::  var (:,:)
    complex, intent(out), optional        ::  stdv(:,:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_c_Get2d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! copy:
    if ( present(n   ) ) n    = self%n  (:,:,1)
    if ( present(mean) ) mean = self%m_n(:,:,1)
    if ( present(var ) ) var  = self%v_n(:,:,1)
    
    ! process:
    if ( present(stdv) ) then
      ! compute from variance if possible:
      where ( self%n(:,:,1) > 1 )
        stdv = sqrt( self%v_n(:,:,1) )
      else where
        stdv = 0.0
      end where
    end if

    ! ok
    status = 0

  end subroutine SampStat_c_Get2d


  ! ***
  
  
  subroutine SampStat_c_Get3d( self, status, n, mean, var, stdv )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat_c), intent(inout)    ::  self
    integer, intent(out)                  ::  status
    
    complex, intent(out), optional        ::  n   (:,:,:)
    complex, intent(out), optional        ::  mean(:,:,:)
    complex, intent(out), optional        ::  var (:,:,:)
    complex, intent(out), optional        ::  stdv(:,:,:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_c_Get3d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! copy:
    if ( present(n   ) ) n    = self%n  (:,:,:)
    if ( present(mean) ) mean = self%m_n(:,:,:)
    if ( present(var ) ) var  = self%v_n(:,:,:)
    
    ! process:
    if ( present(stdv) ) then
      ! compute from variance if possible:
      where ( self%n(:,:,:) > 1 )
        stdv = sqrt( self%v_n(:,:,:) )
      else where
        stdv = 0.0
      end where
    end if

    ! ok
    status = 0

  end subroutine SampStat_c_Get3d
  
  
  ! ***
  
  
  subroutine SampStat_c_AddSample_2d( self, x_n, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat_c), intent(inout)      ::  self
    complex, intent(in)                     ::  x_n(:,:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_c_AddSample_2d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! copy:
    self%x(:,:,1) = x_n
    ! add:
    call SampStat_c_AddSample_3d( self, self%x, status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0

  end subroutine SampStat_c_AddSample_2d
  
  
  ! ***
  
  
  subroutine SampStat_c_AddSample_3d( self, x_n, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampStat_c), intent(inout)      ::  self
    complex, intent(in)                     ::  x_n(:,:,:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampStat_c_AddSample_3d'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    !$OMP parallel
    !$OMP   workshare
    
    ! copy current values:
    self%m_nm1  = self%m_n
    self%ms_nm1 = self%ms_n
    self%v_nm1  = self%v_n
    
    ! increase sample counters :
    self%n = self%n + 1

    ! update means and mean squares:
    self%m_n  = ( real(self%n-1)*self%m_nm1  + x_n    ) / real(self%n)
    self%ms_n = ( real(self%n-1)*self%ms_nm1 + x_n**2 ) / real(self%n)

    ! difference fields:
    self%d  = (      x_n - self%m_nm1 ) / real(self%n)
    self%dm =   self%m_n - self%m_nm1
    ! update variance, only if at least 2 samples present:
    where ( self%n > 1 )
      self%v_n =              real(self%n-2)/real(self%n-1) * self%v_nm1 + &
                 real(self%n)*real(self%n+1)/real(self%n-1) * self%d * conjg(self%d)  - &
                              real(self%n*1)/real(self%n-1) * ( self%dm * conjg(self%d)  + self%d * conjg(self%dm) )
    end where
    
    !$OMP   end workshare
    !$OMP end parallel

    ! ok
    status = 0

  end subroutine SampStat_c_AddSample_3d


  ! ********************************************************************
  ! ***
  ! *** covariance between field and points
  ! ***
  ! ********************************************************************


  subroutine SampCovr_3D_1D_Init( self, xshp, yshp, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_3D_1D), intent(out)  ::  self
    integer, intent(in)                   ::  xshp(:)
    integer, intent(in)                   ::  yshp
    integer, intent(out)                  ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_3D_1D_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! count:
    self%xrank = size(xshp)
    ! check ...
    if ( self%xrank > max_rank ) then
      write (*,'("ERROR - size of shp ",i6," exceeds max rank ",i6)') self%xrank, max_rank
      TRACEBACK; status=1; return
    end if
    ! store:
    self%xshp = 1
    self%xshp(1:self%xrank) = xshp
    
    ! count:
    self%yshp = yshp
    
    ! storage for end results:
    allocate( self%xm_n  (self%xshp(1),self%xshp(2),self%xshp(3)) )
    allocate( self%xv_n  (self%xshp(1),self%xshp(2),self%xshp(3)) )
    allocate( self%ym_n  (self%yshp) )
    allocate( self%yv_n  (self%yshp) )
    allocate( self%cv_n  (self%xshp(1),self%xshp(2),self%xshp(3),self%yshp) )
    ! work arrays:
    allocate( self%xm_nm1 (self%xshp(1),self%xshp(2),self%xshp(3)) )
    allocate( self%xv_nm1 (self%xshp(1),self%xshp(2),self%xshp(3)) )
    allocate( self%ym_nm1 (self%yshp) )
    allocate( self%yv_nm1 (self%yshp) )
!    allocate( self%ms_nm1(self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%cv_nm1 (self%xshp(1),self%xshp(2),self%xshp(3),self%yshp) )
!    allocate( self%x     (self%shp(1),self%shp(2),self%shp(3)) )
    allocate( self%xd     (self%xshp(1),self%xshp(2),self%xshp(3)) )
    allocate( self%xdm    (self%xshp(1),self%xshp(2),self%xshp(3)) )
    
    ! no samples yet:
    self%n    = 0
    self%xm_n = 0.0
    self%ym_n = 0.0
    self%cv_n = 0.0

    ! ok
    status = 0

  end subroutine SampCovr_3D_1D_Init


  ! ***
  
  
  subroutine SampCovr_3D_1D_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_3D_1D), intent(inout)  ::  self
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_3D_1D_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! dummy:
    self%xrank = -999
    self%xshp  = -999
    self%yshp  = -999
    
    ! clear:
    deallocate( self%xm_n   )
    deallocate( self%xv_n   )
    deallocate( self%ym_n   )
    deallocate( self%yv_n   )
    deallocate( self%cv_n   )
    ! clear work arrays:
    deallocate( self%xm_nm1  )
    deallocate( self%xv_nm1  )
    deallocate( self%ym_nm1  )
    deallocate( self%yv_nm1  )
!    deallocate( self%ms_nm1 )
    deallocate( self%cv_nm1  )
!    deallocate( self%x      )
    deallocate( self%xd      )
    deallocate( self%xdm     )

    ! ok
    status = 0

  end subroutine SampCovr_3D_1D_Done


  ! ***
  
  
  subroutine SampCovr_3D_1D_GetPoint( self, iy, status, stdv )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_3D_1D), intent(inout)    ::  self
    integer, intent(in)                       ::  iy
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  stdv

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_3D_1D_GetPoint'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (iy < 1) .or. (iy > self%yshp) ) then
      write (gol,'("y index ",i4," outside expected range 1 .. ",i4)') iy, self%yshp; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! process:
    if ( present(stdv) ) then
      ! compute from variance if possible:
      if ( self%n > 1 ) then
        stdv = sqrt( self%yv_n(iy) )
      else
        stdv = 0.0
      end if
    end if

    ! ok
    status = 0

  end subroutine SampCovr_3D_1D_GetPoint


  ! ***
  
  
  subroutine SampCovr_3D_1D_Get3D( self, iy, status, covr, corr )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_3D_1D), intent(inout)  ::  self
    integer, intent(in)                     ::  iy
    integer, intent(out)                    ::  status
    
    real, intent(out), optional             ::  covr(:,:,:)
    real, intent(out), optional             ::  corr(:,:,:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_3D_1D_Get3D'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! check ...
    if ( (iy < 1) .or. (iy > self%yshp) ) then
      write (gol,'("y index ",i4," outside expected range 1 .. ",i4)') iy, self%yshp; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! extract covariance?
    if ( present(covr) )  then
      ! check:
      if ( any( shape(covr) /= self%xshp(1:3) ) ) then
        write (gol,'("shape of covr (",i0,2(",",i0),") does not match with x (",i0,2(",",i0),")")') &
                shape(covr), self%xshp(1:3); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      covr = self%cv_n(:,:,:,iy)
    end if
    
    ! extract correlation?
    if ( present(corr) )  then
      ! check:
      if ( any( shape(corr) /= self%xshp(1:3) ) ) then
        write (gol,'("shape of covr (",i0,2(",",i0),") does not match with x (",i0,2(",",i0),")")') &
                shape(corr), self%xshp(1:3); call goErr
        TRACEBACK; status=1; return
      end if
      ! init with covariance:
      corr = self%cv_n(:,:,:,iy)
      ! scale with x-stdv (sqrt of x-variance):
      where ( self%xv_n > 0.0 )
        corr = corr / sqrt(self%xv_n)
      endwhere
      ! scale with y-stdv (sqrt of y-variance):
      if ( self%yv_n(iy) > 0.0 ) then
        corr = corr / sqrt(self%yv_n(iy))
      end if
    end if

    ! ok
    status = 0

  end subroutine SampCovr_3D_1D_Get3D
  
  
  ! ***
  
  
  subroutine SampCovr_3D_1D_AddSample_3d( self, x_n, y_n, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_3D_1D), intent(inout)  ::  self
    real, intent(in)                        ::  x_n(:,:,:)
    real, intent(in)                        ::  y_n(:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_3D_1D_AddSample_3d'
    
    ! --- local ----------------------------------
    
    integer         ::  i
    
    ! --- begin ----------------------------------
    
    !xOMP parallel
    !xOMP   workshare
    
    ! copy current values:
    self%xm_nm1  = self%xm_n
    self%xv_nm1  = self%xv_n
    self%ym_nm1  = self%ym_n
    self%yv_nm1  = self%yv_n
    self%cv_nm1  = self%cv_n
    
    ! increase sample counter :
    self%n = self%n + 1

    ! update mean field, same for every point:
    self%xm_n  = ( real(self%n-1)*self%xm_nm1 + x_n ) / real(self%n)
    !self%ms_n = ( real(self%n-1)*self%ms_nm1 + x_n**2 ) / real(self%n)

    ! difference fields:
    self%xd  = (       x_n - self%xm_nm1 ) / real(self%n)
    self%xdm =   self%xm_n - self%xm_nm1

    ! update (co)variance, only if at least 2 samples present:
    if ( self%n > 1 ) then
      ! variance of x values:
      self%xv_n    =              real(self%n-2)/real(self%n-1) *      self%xv_nm1  + &
                     real(self%n)*real(self%n+1)/real(self%n-1) *   self%xd  * self%xd  - &
                                  real(self%n*1)/real(self%n-1) * ( self%xdm * self%xd  + self%xd * self%xdm )
    end if

    ! loop over point values:
    do i = 1, self%yshp

      ! update mean point values, same over every cell:
      self%ym_n(i) = ( real(self%n-1)*self%ym_nm1(i) + y_n(i)    ) / real(self%n)
      !self%ms_n(i) = ( real(self%n-1)*self%ms_nm1(i) + y_n(i)**2 ) / real(self%n)

      ! difference points:
      self%yd  = (       y_n(i) - self%ym_nm1(i) ) / real(self%n)
      self%ydm =   self%ym_n(i) - self%ym_nm1(i)

      ! update (co)variance, only if at least 2 samples present:
      if ( self%n > 1 ) then
        ! variance of y values:
        self%yv_n(i) =              real(self%n-2)/real(self%n-1) *      self%yv_nm1(i) + &
                       real(self%n)*real(self%n+1)/real(self%n-1) *   self%yd  * self%yd  - &
                                    real(self%n*1)/real(self%n-1) * ( self%ydm * self%yd  + self%yd * self%ydm )
        ! covar:
        self%cv_n(:,:,:,i) =              real(self%n-2)/real(self%n-1) *      self%cv_nm1(:,:,:,i) + &
                             real(self%n)*real(self%n+1)/real(self%n-1) *   self%xd  * self%yd  - &
                                          real(self%n*1)/real(self%n-1) * ( self%xdm * self%yd  + self%xd * self%ydm )
      end if

    end do

    !xOMP   end workshare
    !xOMP end parallel

    ! ok
    status = 0

  end subroutine SampCovr_3D_1D_AddSample_3d


  ! ********************************************************************
  ! ***
  ! *** covariance between pair of 2D fields
  ! ***
  ! ********************************************************************


  subroutine SampCovr_2D_2D_Init( self, xshp, yshp, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_2D_2D), intent(out)      ::  self
    integer, intent(in)                       ::  xshp(:)
    integer, intent(in)                       ::  yshp(:)
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_2D_2D_Init'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! count:
    self%xrank = size(xshp)
    ! check ...
    if ( self%xrank > max_rank ) then
      write (*,'("ERROR - size of shp ",i6," exceeds max rank ",i6)') self%xrank, max_rank
      TRACEBACK; status=1; return
    end if
    ! store:
    self%xshp = 1
    self%xshp(1:self%xrank) = xshp
    
    ! count:
    self%yrank = size(yshp)
    ! check ...
    if ( self%yrank > max_rank ) then
      write (*,'("ERROR - size of shp ",i6," eyceeds max rank ",i6)') self%yrank, max_rank
      TRACEBACK; status=1; return
    end if
    ! store:
    self%yshp = 1
    self%yshp(1:self%yrank) = yshp
    
    ! storage for end results:
    allocate( self%xm_n  (self%xshp(1),self%xshp(2)) )
    allocate( self%xv_n  (self%xshp(1),self%xshp(2)) )
    allocate( self%ym_n  (self%yshp(1),self%yshp(2)) )
    allocate( self%yv_n  (self%yshp(1),self%yshp(2)) )
    allocate( self%cv_n  (self%xshp(1),self%xshp(2),self%yshp(1),self%yshp(2)) )
    ! work arrays:
    allocate( self%xm_nm1 (self%xshp(1),self%xshp(2)) )
    allocate( self%xv_nm1 (self%xshp(1),self%xshp(2)) )
    allocate( self%ym_nm1 (self%yshp(1),self%yshp(2)) )
    allocate( self%yv_nm1 (self%yshp(1),self%yshp(2)) )
    allocate( self%cv_nm1 (self%xshp(1),self%xshp(2),self%yshp(1),self%yshp(2)) )
    allocate( self%xd     (self%xshp(1),self%xshp(2)) )
    allocate( self%xdm    (self%xshp(1),self%xshp(2)) )
    
    ! no samples yet:
    self%n    = 0
    self%xm_n = 0.0
    self%ym_n = 0.0
    self%cv_n = 0.0

    ! ok
    status = 0

  end subroutine SampCovr_2D_2D_Init


  ! ***
  
  
  subroutine SampCovr_2D_2D_Done( self, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_2D_2D), intent(inout)    ::  self
    integer, intent(out)                      ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_2D_2D_Done'
    
    ! --- local ----------------------------------
    
    ! --- begin ----------------------------------
    
    ! dummy:
    self%xrank = -999
    self%xshp  = -999
    self%yshp  = -999
    
    ! clear:
    deallocate( self%xm_n   )
    deallocate( self%xv_n   )
    deallocate( self%ym_n   )
    deallocate( self%yv_n   )
    deallocate( self%cv_n   )
    ! clear work arrays:
    deallocate( self%xm_nm1  )
    deallocate( self%xv_nm1  )
    deallocate( self%ym_nm1  )
    deallocate( self%yv_nm1  )
    deallocate( self%cv_nm1  )
    deallocate( self%xd      )
    deallocate( self%xdm     )

    ! ok
    status = 0

  end subroutine SampCovr_2D_2D_Done


  ! ***
  
  
  subroutine SampCovr_2D_2D_GetXY( self, status, covr, corr )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_2D_2D), intent(inout)    ::  self
    integer, intent(out)                      ::  status
    
    real, intent(out), optional               ::  covr(:,:,:,:)
    real, intent(out), optional               ::  corr(:,:,:,:)

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_2D_2D_GetXY'
    
    ! --- local ----------------------------------
    
    integer     ::  ix1, ix2
    integer     ::  iy1, iy2
    
    ! --- begin ----------------------------------
    
    ! extract covariance?
    if ( present(covr) ) then
      ! check:
      if ( any( shape(covr) /= (/self%xshp(1),self%xshp(2),self%yshp(1),self%yshp(2)/) ) ) then
        write (gol,'("shape of covr (",i0,3(",",i0),") does not match with x (",i0,",",i0,") and y (",i0,",",i0,")")') &
                shape(covr), self%xshp(1:2), self%yshp(1:2); call goErr
        TRACEBACK; status=1; return
      end if
      ! copy:
      covr = self%cv_n
    end if
    
    ! extract correlation?
    if ( present(corr) ) then
      ! check:
      if ( any( shape(corr) /= (/self%xshp(1),self%xshp(2),self%yshp(1),self%yshp(2)/) ) ) then
        write (gol,'("shape of covr (",i0,3(",",i0),") does not match with x (",i0,",",i0,") and y (",i0,",",i0,")")') &
                shape(corr), self%xshp(1:2), self%yshp(1:2); call goErr
        TRACEBACK; status=1; return
      end if
      ! init with covariance:
      corr = self%cv_n
      ! scale with x-stdv (sqrt of x-variance):
      do ix2 = 1, self%xshp(2)
        do ix1 = 1, self%xshp(1)
          if ( self%xv_n(ix1,ix2) > 0.0 ) then
            corr(ix1,ix2,:,:) = corr(ix1,ix2,:,:) / sqrt(self%xv_n(ix1,ix2))
          end if
        end do ! ix1
      end do ! ix2
      ! scale with y-stdv (sqrt of y-variance):
      do iy2 = 1, self%yshp(2)
        do iy1 = 1, self%yshp(1)
          if ( self%yv_n(iy1,iy2) > 0.0 ) then
            corr(:,:,iy1,iy2) = corr(:,:,iy1,iy2) / sqrt(self%yv_n(iy1,iy2))
          end if
        end do ! iy1
      end do ! iy2
    end if ! corr

    ! ok
    status = 0

  end subroutine SampCovr_2D_2D_GetXY
  
  
  ! ***
  
  
  subroutine SampCovr_2D_2D_AddSample_2d( self, x_n, y_n, status )
  
    ! --- in/out ---------------------------------
    
    class(T_SampCovr_2D_2D), intent(inout)  ::  self
    real, intent(in)                        ::  x_n(:,:)
    real, intent(in)                        ::  y_n(:,:)
    integer, intent(out)                    ::  status

    ! --- const --------------------------------------

    character(len=*), parameter  ::  rname = mname//'/SampCovr_2D_2D_AddSample_2d'
    
    ! --- local ----------------------------------
    
    integer         ::  i, j
    
    ! --- begin ----------------------------------
    
    !xOMP parallel
    !xOMP   workshare
    
    ! copy current values:
    self%xm_nm1  = self%xm_n
    self%xv_nm1  = self%xv_n
    self%ym_nm1  = self%ym_n
    self%yv_nm1  = self%yv_n
    self%cv_nm1  = self%cv_n
    
    ! increase sample counter :
    self%n = self%n + 1

    ! update mean field, same for every point:
    self%xm_n = ( real(self%n-1)*self%xm_nm1 + x_n    ) / real(self%n)
    !self%ms_n = ( real(self%n-1)*self%ms_nm1 + x_n**2 ) / real(self%n)

    ! difference fields:
    self%xd  = (       x_n - self%xm_nm1 ) / real(self%n)
    self%xdm =   self%xm_n - self%xm_nm1

    ! update covariance, only if at least 2 samples present:
    if ( self%n > 1 ) then
      ! variance of x values:
      self%xv_n    =              real(self%n-2)/real(self%n-1) *      self%xv_nm1  + &
                     real(self%n)*real(self%n+1)/real(self%n-1) *   self%xd  * self%xd  - &
                                  real(self%n*1)/real(self%n-1) * ( self%xdm * self%xd  + self%xd * self%xdm )
    end if

    ! loop over y values:
    do j = 1, self%yshp(2)
      do i = 1, self%yshp(1)

        ! update mean point values, same over every cell:
        self%ym_n(i,j) = ( real(self%n-1)*self%ym_nm1(i,j) + y_n(i,j)    ) / real(self%n)
        !self%ms_n(i,j) = ( real(self%n-1)*self%ms_nm1(i,j) + y_n(i,j)**2 ) / real(self%n)

        ! difference points:
        self%yd  = (       y_n(i,j) - self%ym_nm1(i,j) ) / real(self%n)
        self%ydm =   self%ym_n(i,j) - self%ym_nm1(i,j)

        ! update (co)variance, only if at least 2 samples present:
        if ( self%n > 1 ) then
          ! variance of y values:
          self%yv_n(i,j) =              real(self%n-2)/real(self%n-1) *      self%yv_nm1(i,j) + &
                         real(self%n)*real(self%n+1)/real(self%n-1) *   self%yd  * self%yd  - &
                                      real(self%n*1)/real(self%n-1) * ( self%ydm * self%yd  + self%yd * self%ydm )
          ! covar:
          self%cv_n(:,:,i,j) =              real(self%n-2)/real(self%n-1) *      self%cv_nm1(:,:,i,j) + &
                               real(self%n)*real(self%n+1)/real(self%n-1) *   self%xd  * self%yd  - &
                                            real(self%n*1)/real(self%n-1) * ( self%xdm * self%yd  + self%xd * self%ydm )
        end if

      end do ! i
    end do ! j
    
!    ! testing diagonal ...
!    if ( self%n > 1 ) then
!      write (*,'(" ")')
!      write (*,'("---> samples and variance")')
!      do j = 1, self%yshp(2)
!        do i = 1, self%yshp(1)
!          write (*,'("samples ",2i3,4es16.6)') i, j, x_n(i,j), self%xv_n(i,j), y_n(i,j), self%yv_n(i,j)
!        end do
!      end do
!      write (*,'(" ")')
!      write (*,'("---> covr diag and variances")')
!      do j = 1, self%yshp(2)
!        do i = 1, self%yshp(1)
!          if ( (self%xv_n(i,j) > 0.0) .and. (self%yv_n(i,j) > 0.0) ) then
!            write (*,'("diag ",2i3,4es16.6)') i, j, self%cv_n(i,j,i,j), self%xv_n(i,j), self%yv_n(i,j), self%cv_n(i,j,i,j)/sqrt(self%xv_n(i,j))/sqrt(self%yv_n(i,j))
!          else
!            write (*,'("diag ",2i3,4es16.6)') i, j, self%cv_n(i,j,i,j), self%xv_n(i,j), self%yv_n(i,j), 0.0
!          end if
!        end do
!      end do
!      write (*,'(" ")')
!      stop 'break after first diag'
!    end if
  
    !xOMP   end workshare
    !xOMP end parallel

    ! ok
    status = 0

  end subroutine SampCovr_2D_2D_AddSample_2d


end module GO_SampStat



