!######################################################################
!
! LinAlg_LaPack_350 - interface to selected lapack 3.5.0 routines
!
!######################################################################
!
#define TRACEBACK write (*,'("ERROR in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__
!
!######################################################################

module LinAlg_LaPack_350

#ifdef with_lapack

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  SyEvR
  

  ! --- interface -----------------------------------
  
  interface SyEvR
    module procedure lapack_350_dSyEvR
  end interface
  
  
contains


  !===================================================================
  ! SyEvR
  !
  ! Follow MKL implementation:
  !  https://software.intel.com/en-us/node/469188#DE6F1082-B885-4C07-BAC2-4E844ACC14E3
  !
  ! A       real(wp), dimension(n,n), intent(inout)
  !         On input:
  !           Array containing either upper or lower triangiuar part 
  !           of the symmetric matrix A (size n) as specified by uplo.
  !           The second dimension of a must be at least max(1, n).
  !         On exit:
  !           The lower triangle (if uplo = 'L') or the upper triangle (if uplo = 'U') 
  !           of A, including the diagonal, is overwritten.
  !
  ! w       Selected eigenvalues in ascending order, stored in w(1) to w(m) .
  !
  ! il,iu   integer, intent(in), optional
  !         Only compute eigen values from number il onwards (defaiut 1)
  !         or up to iu (defaiut n).
  !       
  ! vl,vu   integer, intent(in), optional
  !         Only compute eigen values from value vl onwares (defaiut -huge(vl))
  !         or up to vu (defaiut +huge(vu)).
  !
  ! uplo    character(len=1), intent(in), optional
  !         Must be 'U' (defaiut) or 'L'.
  !         If uplo = 'U', a stores the upper triangiuar part of A.
  !         If uplo = 'L', a stores the lower triangiuar part of A.
  !
  ! m       integer, intent(out), optional
  !         The total number of eigenvalues found, 0 <= m <= n.
  !
  ! isuppz  integer, dimension(:), optional
  !         Size at least 2 *max(1, m).
  !         The support of the eigenvectors in Z, i.e., 
  !         the indices indicating the nonzero elements in z. 
  !         The i-th eigenvector is nonzero only in elements isuppz( 2i-1) 
  !         through isuppz( 2i ). 
  !         Referenced only if eigenvectors are needed 
  !         and all eigenvalues are needed
  !
  ! abstol  real(wp), intent(in), optional
  !         The absolute error tolerance to which each eigenvalue/eigenvector 
  !         is required.
  !
  ! info    integer, intent(out), optional
  !         If info =  0, the execution is successfiu.
  !         If info = -i, the i-th parameter had an illegal value.
  !         If info =  i, an internal error has occurred.
  !
  !===================================================================
  
  subroutine lapack_350_dSyEvR( A, w, &
                                uplo, Z, vl, vu, il, iu, m, &
                                isuppz, abstol, info )
                                
    ! --- const ----------------------------------
    
    ! routine name:
    character(len=*), parameter   ::  rname = 'lapack_350_dSyEvR'
    
    ! working precision:
    integer, parameter    ::  wp = 8
                                
    ! --- in/out ---------------------------------
    
    real(wp), intent(inout)                   ::  A(:,:)
    real(wp), intent(out)                     ::  w(:)
    character(len=1), intent(in), optional    ::  uplo
    real(wp), intent(out), optional           ::  Z(:,:)
    integer, intent(in), optional             ::  il, iu
    real(wp), intent(in), optional            ::  vl, vu
    integer, intent(out), optional            ::  m
    integer, intent(out), optional            ::  isuppz(:)
    real(wp), intent(in), optional            ::  abstol
    integer, intent(out), optional            ::  info
    
    ! --- local ----------------------------------
    
    character(len=1)        ::  jobz
    character(len=1)        ::  range
    character(len=1)        ::  uplo_
    integer                 ::  n
    integer                 ::  lda
    integer                 ::  il_, iu_
    real(wp)                ::  vl_, vu_
    real(wp)                ::  abstol_
    integer                 ::  m_
    integer                 ::  ldz, sdz
    real(wp), allocatable   ::  Z_(:,:)
    integer                 ::  lisuppz
    integer, allocatable    ::  isuppz_(:)
    integer                 ::  lwork
    real(wp), allocatable   ::  work(:)
    integer                 ::  liwork
    integer, allocatable    ::  iwork(:)
    integer                 ::  status
    
    ! --- local ----------------------------------
    
    ! what to compute?
    if ( present(Z) ) then
      ! eigenvalues and vectors:
      jobz = 'V'
    else
      ! only eigenvalues:
      jobz = 'N'
    end if
    
    ! range of computed eigenvalues:
    if ( present(vl) .or. present(vu) ) then
      ! value range:
      range = 'V'
    else if ( present(il) .or. present(iu) ) then
      ! index range:
      range = 'I'
    else
      ! all:
      range = 'A'
    end if
    
    ! storage:
    if ( present(uplo) ) then
      uplo_ = uplo
    else
      uplo_ = 'U' ! upper triangle
    end if
    
    ! logical shape:
    n = size(A,1)
    ! actual shape:
    lda = size(A,1)
    
    ! value range:
    if ( present(vl) ) then
      vl_ = vl
    else
      vl_ = -huge(vl_)
    end if
    if ( present(vu) ) then
      vu_ = vu
    else
      vu_ = huge(vu_)
    end if
    
    ! index range:
    if ( present(il) ) then
      il_ = il
    else
      il_ = 1
    end if
    if ( present(iu) ) then
      iu_ = iu
    else
      iu_ = n
    end if
    
    ! storage for eigenvectors:
    if ( present(Z) ) then
      ! shape:
      ldz = size(Z,1)
      sdz = size(Z,2)
    else
      ! shape of input:
      ldz = n
      sdz = n
    end if
    ! storage:
    allocate( Z_(ldz,sdz), stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - allocating Z_ with shape (",i0,",",i0,")")') ldz, sdz
      if ( present(info) ) info = 1
      TRACEBACK; return
    end if
    
    ! storage support values:
    if ( present(isuppz) ) then
      lisuppz = size(isuppz)
    else
      lisuppz = 2 * n
    end if
    ! storage:
    allocate( isuppz_(lisuppz), stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - allocating isuppz_ with size (",i0,")")') lisuppz
      if ( present(info) ) info = 1
      TRACEBACK; return
    end if

    ! length of work array:
    lwork = 26*n
    ! storage:
    allocate( work(lwork), stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - allocating work with size (",i0,")")') lwork
      if ( present(info) ) info = 1
      TRACEBACK; return
    end if

    ! length of work array:
    liwork = 10*n
    ! storage:
    allocate( iwork(liwork), stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - allocating iwork with size (",i0,")")') liwork
      if ( present(info) ) info = 1
      TRACEBACK; return
    end if
    
    ! call f77 routine:
    call dSyEvR( jobz, range, uplo_, n, A, lda, &
                  vl_, vu_, il_, iu_, abstol_, &
                  m_, w, Z_, ldz, isuppz_, &
                  work, lwork, iwork, liwork, status )
    ! check ...
    if ( status /= 0 ) then
      write (*,'("ERROR - from dSyEvR ; info = ",i0)') status
      if ( present(info) ) info = status
      TRACEBACK; return
    end if
  
    ! copy output:
    if ( present(m) ) m = m_
    if ( present(Z) ) Z = Z_
    if ( present(info) ) info = status
    if ( present(isuppz) ) isuppz = isuppz_

    ! clear:
    deallocate( Z_, stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - deallocating Z_")')
      if ( present(info) ) info = 1
      TRACEBACK; return
    end if
    deallocate( isuppz_, stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - deallocating isuppz_")')
      if ( present(info) ) info = 1
      TRACEBACK; return
    end if
    deallocate( work, stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - deallocating work")')
      if ( present(info) ) info = 1
      TRACEBACK; return
    end if
    deallocate( iwork, stat=status )
    if ( status /= 0 ) then
      write (*,'("ERROR - deallocating iwork")')
      if ( present(info) ) info = 1
      TRACEBACK; return
    end if
    
  end subroutine lapack_350_dSyEvR

#endif
    
end module LinAlg_LaPack_350
