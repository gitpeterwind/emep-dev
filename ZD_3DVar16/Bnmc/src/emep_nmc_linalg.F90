!######################################################################
!
! EMEP NMC background covariance.
! Linear Algebra tools for matrices with (tracer,level) pair as dimension.
!
!######################################################################
!
#define TRACEBACK write (gol,'("in ",a," (",a,", line ",i5,")")') rname, __FILE__, __LINE__; call goErr
#define IF_NOT_OK_RETURN(action) if (status/=0) then; TRACEBACK; action; return; end if
!
!######################################################################

module EMEP_NMC_LinAlg

  use GO    , only : gol, goPr, goErr

  implicit none
  
  
  ! --- in/out -----------------------------------
  
  private
  
  public  ::  Decompose_Bhat
  

  ! --- const ----------------------------------------

  character(len=*), parameter  ::  mname = 'MEP_BCovarSqrt_LinAlg'
  

contains


  ! ********************************************************************
  ! ***            ^
  ! *** decommpose B = X Lambda X^T
  ! ***
  ! ********************************************************************

  
  subroutine Decompose_Bhat( Bhat, X, Lambda, nev, status )
  
    use LinAlg, only : SyEvR

    ! --- in/out ---------------------------------
    
    real, intent(in)                ::  Bhat(:,:,:,:)   ! (nlev,nlev,ntr,ntr)
    real, intent(out)               ::  X(:,:,:)        ! (nlev,ntr,nv)
    real, intent(out)               ::  Lambda(:)       ! (nv)
    integer, intent(out)            ::  nev
    integer, intent(out)            ::  status

    ! --- const ----------------------------------

    character(len=*), parameter  ::  rname = mname//'/Decompose_Bhat'
  
    ! --- local ---------------------------------- 
    
    integer                          ::  nlev
    integer                          ::  ntracer
    integer                          ::  ik
    integer                          ::  nv, iv, iv1, iv2
    integer                          ::  ilev, ilev1, ilev2
    integer                          ::  itracer, itracer1, itracer2
    real, allocatable                ::  A(:,:)  ! (nv,nv)
    real, allocatable                ::  Z(:,:)  ! (nv,nv)
    real, allocatable                ::  w(:)  ! (nv)
    integer                          ::  m, m1, m2
    integer                          ::  nzero

    ! --- local ----------------------------------
    
    ! sizes:
    nlev    = size(Bhat,1)
    ntracer = size(Bhat,3)
    
    ! maximum number of eigenvalues:
    nv = nlev * ntracer
    
    ! check ..
    if ( any( shape(Bhat) /= (/nlev,nlev,ntracer,ntracer/) ) ) then
      write (gol,'("shape of Bhat is (",i0,3(",",i0),") while nlev=",i0," and ntracer=",i0)') &
              shape(Bhat), nlev, ntracer; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    if ( any( shape(X) /= (/nlev,ntracer,nv/) ) ) then
      write (gol,'("shape of X is (",i0,2(",",i0),") while nlev=",i0," and ntracer=",i0)') &
              shape(X), nlev, ntracer; call goErr
      TRACEBACK; status=1; return
    end if
    ! check ..
    if ( size(Lambda) /= nv ) then
      write (gol,'("shape of Lambda is (",i0,") while nlev=",i0," and ntracer=",i0)') &
              shape(Lambda), nlev, ntracer; call goErr
      TRACEBACK; status=1; return
    end if
    
    ! storage for eigenvalue decompo:
    allocate( A(nv,nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( w(nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
    allocate( Z(nv,nv), stat=status )
    IF_NOT_OK_RETURN(status=1)
      
    ! pack covar in matrix array;
    A = -999.9
    ! first dim:
    iv1 = 0
    do itracer1 = 1, ntracer
      do ilev1 = 1, nlev
        iv1 = iv1 + 1
        ! second dim:
        iv2 = 0
        do itracer2 = 1, ntracer
          do ilev2 = 1, nlev
            iv2 = iv2 + 1
            ! copy:
            A(iv1,iv2) = Bhat(ilev1,ilev2,itracer1,itracer2)
          end do ! ilev2
        end do ! itracer2
      end do ! ilev1
    end do ! itracer1

    ! count number of zero rows:
    nzero = 0
    do iv = 1, nv
      if ( all( A(iv,:) == 0.0 ) ) nzero = nzero + 1
    end do

    ! eigenvalue decomposition,
    ! only positive eigenvalues:
    call SyEvr( A, w, Z=Z, m=m, vl=0.0, abstol=1.0e-4, info=status )
    if ( status < 0 ) then
      write (gol,'("SyEvr: illegal value for parameter ",i0)') abs(status); call goErr
      TRACEBACK; status=1; return
    else if ( status > 0 ) then
      write (gol,'("SyEvr: internal error, info=",i0)') status; call goErr
      TRACEBACK; status=1; return
    end if

    ! skip nearly zero values:
    if ( m > nv-nzero ) then
      m1 = m - (nv-nzero) + 1
      m2 = m
    else
      m1 = 1
      m2 = m
    end if

    ! reset counter:
    m = m2 - m1 + 1

    ! number of valid eigenvalues:
    nev = m

    ! loop over valid eigenvalues
    do iv = 1, m
      ! lev/tracer dim:
      iv2 = 0
      do itracer = 1, ntracer
        do ilev = 1, nlev
          iv2 = iv2 + 1
          ! copy eigenvector, reverse order:
          X(ilev,itracer,iv) = Z(iv2,m2+1-iv)
        end do
      end do
      ! copy eigenvalue, reverse order:
      Lambda(iv) = w(m2+1-iv)
    end do  ! eigenvalues
    
    ! clear:
    deallocate( A, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( w, stat=status )
    IF_NOT_OK_RETURN(status=1)
    deallocate( Z, stat=status )
    IF_NOT_OK_RETURN(status=1)

    ! ok
    status = 0
  
  end subroutine Decompose_Bhat


end module EMEP_NMC_LinAlg


