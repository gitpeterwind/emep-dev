module DA_Obs_ml
use MetFields_ml,     only: z_bnd,z_mid,roa,th,q,ps
use TimeDate_ml,      only: current_date
use CheckStop_ml,     only: CheckStop
use Functions_ml,     only: great_circle_distance
use GridValues_ml,    only: glon,glat,coord_in_processor
use Par_ml,           only: limax,ljmax,me
use Chemfields_ml,    only: xn_adv,cfac,PM25_water,PM25_water_rh50
!se ChemSpecs_shl_ml, only: NSPEC_SHL        ! Maps indices
use ChemSpecs_adv_ml, only: NSPEC_ADV
use ChemChemicals_ml, only: species
use ChemGroups_ml,    only: chemgroups
use Io_ml,            only: IO_TMP,PrintLog
use InterpolationRoutines_ml,  only : point2grid_coeff
use ModelConstants_ml,only: KMAX_MID,MasterProc,runlabel1
use SmallUtils_ml,    only: find_index
use TimeDate_ExtraUtil_ml, only: date2string
use Units_ml,         only: Units_Scale,Group_Units
use DA_ml,            only: debug=>DA_DEBUG,dafmt=>da_fmt_msg,damsg=>da_msg
use Util_ml,          only: norm
use spectralcov,      only: nx,ny,nlev,nchem,nxex,nyex,nv1,FGSCALE,&
                            nChemObs,nChemNoObs,iChemInv,iChemObs,iChemNoObs
use chitox_ml,       only:  chitox_adj
implicit none
integer, save :: nobs!,matidx
! integer, dimension(:), allocatable, save :: pidx
real, dimension(:), allocatable, save :: innov,obsstddev
real, dimension(:,:,:,:), allocatable, save :: H_jac
type H_jac_idx
  integer :: i(4),j(4),l(0:1)
end type
type(H_jac_idx), dimension(:), allocatable, save :: pidx
!-----------------------------------------------------------------------
type obs_data
  logical            :: set=.false.,found=.false.,unitroa=.false.
  integer            :: ispec=-1,ichem=-1,ichemObs=-1,iobs(0:1)=-1,nobs
  real               :: min=0.0,max=1e3,error_rel=-1e0,error_rep=-1e0,unitconv=-1e0
  character(len=064) :: name='',unit='',deriv=''
  character(len=256) :: file='',tag=''
endtype
integer            :: nobsData=0
integer, parameter :: nobsDataMax=6
type(obs_data)     :: obsData(nobsDataMax)
character(len=032) :: interpolation_method='distance-weighted'
namelist /OBSERVATIONS/nobsData,obsData,interpolation_method
!-----------------------------------------------------------------------
real, dimension(:,:,:,:), allocatable :: xn_adv_ex
integer               :: varSpec(NSPEC_ADV),varSpecInv(NSPEC_ADV)
character(len=016)    :: varName(NSPEC_ADV)='',obsVarName(NSPEC_ADV)=''
logical               :: observedVar(NSPEC_ADV)=.false.
contains
subroutine allocate_obs()
!-----------------------------------------------------------------------
! @description
! Allocate the arrays innov and obsstddev. This step is
! required, since the number of observations can vary from
! one time step to the next
! @author AMVB
!-----------------------------------------------------------------------
  integer :: ierr
! no4=nobs*4
! no8=nobs*8
! no12=nobs*12
! matidx=0
  if(.not.allocated(innov))then
    allocate(innov(nobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: INNOV.')
    innov=0.0
  endif
  if(.not.allocated(obsstddev))then
    allocate(obsstddev(nobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: OBSSTDDEV.')
    obsstddev=0.0
  endif
!-----------------------------------------------------------------------
! H_jac is a matrix with components H_{n;i,j,l,k}, where
!   n=1,...,nobs     (observation index)
!   i=1,...,nx       (longitude index)
!   j=1,...,ny       (latitude index)
!   l=1,...,nlev     (altitude index)
!   k=1,...,nchemobs (observed chemical component)
! Due to the 4-point interpolation in the horizontal direction
! of the model results to the observation point, the Jacobian is
! highly sparse wrt i,j,l. For this reason, the matrix is stored using
! the following compact storage scheme:
! H_{n;i,j,l,k}=H_jac(n,p,l,k), where n=1:nobs,p=1:4,l=l(0):l(1)
!   pidx(n)%i(1)=i_n,   pidx(n)%j(1)=j_n,   pidx(n)%l(0)=l0_n,
!   pidx(n)%i(2)=i_n+1, pidx(n)%j(2)=j_n,   pidx(n)%l(1)=l1_n,
!   pidx(n)%i(3)=i_n,   pidx(n)%j(3)=j_n+1,
!   pidx(n)%i(4)=i_n+1, pidx(n)%j(4)=j_n+1,
!-----------------------------------------------------------------------
  if(.not.allocated(H_jac))then
    allocate(H_jac(nobs,4,nlev,nchemobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: H_JAC.')
    H_jac=0.0
  endif
  if(.not.allocated(pidx))then
    allocate(pidx(nobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: PIDX.')
    pidx=H_jac_idx((/0,0,0,0/),(/0,0,0,0/),(/0,0/))
  endif
end subroutine allocate_obs
subroutine deallocate_obs()
!-----------------------------------------------------------------------
! @description
! Deallocate the arrays innov and obsstddev. This step is
! required, since the number of observations can vary from
! one time step to the next
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none
  if(allocated(innov))    deallocate(innov)
  if(allocated(obsstddev))deallocate(obsstddev)
  if(allocated(H_jac))    deallocate(H_jac)
  if(allocated(pidx))     deallocate(pidx)
endsubroutine deallocate_obs
subroutine read_obs(maxobs,flat,flon,falt,y,stddev,ipar)
!-----------------------------------------------------------------------
! @description
! Read in observations y in standard format
! (Station-No flat flon falt obs-value stddev).
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none

  integer :: maxobs,ipar(maxobs)
  real :: y(maxobs),stddev(maxobs),flat(maxobs),flon(maxobs),falt(maxobs)

  integer            :: no,iobsData,ierr
  character(len=256) :: file
  logical            :: exist
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  nobs=0
  obsData%iobs(0)=0
  obsData%iobs(1)=0
  obsData%nobs=0
!-----------------------------------------------------------------------
! open observational data file in standard format, if applicable
!-----------------------------------------------------------------------
  do iobsData=1,nobsData
    if(obsData(iobsData)%file=="")then
      write(damsg,"(A,I0,1X,A)")'WARNING: obsData#',iobsData,"file not specified"
      if(MasterProc) print dafmt,damsg
      cycle
    endif
    file=date2string(obsData(iobsData)%file,current_date)
    inquire(file=file,exist=exist)
    if(.not.exist)then
      write(damsg,"(A,I0,1X,A)")'WARNING: obsData#',iobsData,"file not found"
      if(MasterProc) print dafmt,damsg
      cycle
    endif
    open(IO_TMP,file=file,form='formatted',status='old',iostat=ierr)
    if(ierr/=0)then
      write(damsg,"(A,I0,1X,A)")'WARNING: obsData#',iobsData,"file not opened"
      if(MasterProc) print dafmt,damsg
      close(IO_TMP)
      cycle
    endif
    print dafmt,"Obsdata opened "//trim(file)
!-----------------------------------------------------------------------
! read data
!-----------------------------------------------------------------------
    obsData(iobsData)%iobs(0)=nobs+1
    do while(ierr==0)
      if(nobs>=maxobs)then
        print dafmt,'ERROR reading observations, maxobs too small'
        ierr=-2
        cycle
      endif
      nobs=nobs+1
      ipar(nobs)=iobsData
      read(IO_TMP,*,iostat=ierr)no,flat(nobs),flon(nobs),falt(nobs),&
                                   y(nobs),stddev(nobs)
      select case(ierr)
      case(0)
        if((.not.coord_in_processor(flon(nobs),flat(nobs))).or.&
           (y(nobs)<obsData(iobsData)%min).or.&
           (y(nobs)>obsData(iobsData)%max))then
         nobs=nobs-1
         cycle
        endif
        stddev(nobs)=max(stddev(nobs),&
          y(nobs)*obsData(iobsData)%error_rel,&
          obsData(iobsData)%error_rep)
        if(stddev(nobs)<=0e0)then
          print dafmt,'WARNING obs stddev <= 0'
          stddev(nobs)=1e-9
        endif
      case(:-1)     !reached EOF: nothing readed
        nobs=nobs-1
      case(1:)
        print dafmt,'ERROR reading observations'
        nobs=nobs-1
        ierr=-1
      endselect
    enddo
    close(IO_TMP)
    obsData(iobsData)%iobs(1)=nobs
  enddo
!-----------------------------------------------------------------------
! Setup obsData &  Write #obs to log file
!-----------------------------------------------------------------------
  do iobsData=1,nobsData
    if(.not.obsData(iobsData)%set)then
      no=find_index(obsData(iobsData)%name,obsVarName(:))
      obsData(iobsData)%found=(no>0)
      if(obsData(iobsData)%found)then
        obsData(iobsData)%ichem=no                    ! index observed
        obsData(iobsData)%ichemObs=ichemobs(no)       ! index observed/unobserved
        obsData(iobsData)%ispec=varSpec(ichemobs(no)) ! index species
        obsData(iobsData)%unitconv=Units_Scale(obsData(iobsData)%unit,&
          obsData(iobsData)%ispec,needroa=obsData(iobsData)%unitroa)
      endif
      call CheckStop(obsData(iobsData)%deriv=='',&
        "read_obs: Inconsistent obsData%deriv")
      write(obsData(iobsData)%tag,"(2(A,1X),'(',A,')')")&
        trim(obsData(iobsData)%name),trim(obsData(iobsData)%deriv),&
        trim(obsData(iobsData)%unit)
      obsData(iobsData)%set=.true.
    endif
    if(any(obsData(iobsData)%iobs/=0)) obsData(iobsData)%nobs=&
      DIM(obsData(iobsData)%iobs(1)+1,obsData(iobsData)%iobs(0))
    call CheckStop(obsData(iobsData)%nobs,count(ipar(:nobs)==iobsData),&
      "read_obs: Inconsistent obsData%nobs")
    if(MasterProc)then
      file=date2string(obsData(iobsData)%file,current_date)
      write(damsg,"('obsData(',I0,') contains ',I0,1X,A,' observations')")&
        iobsData,obsData(iobsData)%nobs,trim(obsData(iobsData)%tag)
      write(damsg,dafmt)trim(damsg)
      call PrintLog(damsg)
      call PrintLog(file)
    endif
  enddo
endsubroutine read_obs
subroutine H_op(lat,lon,alt,n,yn,rho0,ipar)
!-----------------------------------------------------------------------
! @description
! H operator mapping from model to observation space;
! Further: interpolate air density to observation point
! @author AMVB
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: n,ipar
  real, intent(in)    :: lat,lon,alt
  real, intent(inout) :: yn
  real, intent(out)   :: rho0

  integer i(4),j(4),p,l,ll,k,kk,igrp,g,ispec,nn
  real xn,Hj,H_int(4),T1,QML1,PS1,Hchem(nchem),unitconv
  integer, pointer, dimension(:) :: gspec=>null()      ! group array of indexes
  real,    pointer, dimension(:) :: gunit_conv=>null() ! group array of unit conv. factors
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! no4=4*nobs
! no8=8*nobs
!-----------------------------------------------------------------------
! operator for horizontal four-point interpolation between nearest
! neighbouring grid points:
!-----------------------------------------------------------------------
  call get_interpol_const(lon,lat,i,j,H_int,interpolation_method)
!-----------------------------------------------------------------------
! z grid coordinate of observations:
!-----------------------------------------------------------------------
  ll=KMAX_MID
  do while(alt>z_bnd(i(1),j(1),ll).and.ll>=1)
    ll=ll-1
  enddo
  l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
  if(debug.and.MasterProc) then
    print dafmt,'Observation/Model Location:Weight[%],lon[deg],lat[deg],alt[km]'
    print "('#',I0,5(1X,A3,':',F0.1,3(',',F0.2)))",&
    n,'Observation',100.0,lon,lat,alt,&
     ('Model',H_int(p)*100,glon(i(p),j(p)),glat(i(p),j(p)),z_mid(i(p),j(p),ll)*1e-3,p=1,4)
  endif
! interpolated air density:
  rho0=dot_product(H_int(:),(/(roa(i(p),j(p),ll,1),p=1,4)/))
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  nn=n-obsData(ipar)%iobs(0)+1
  if(n==obsData(ipar)%iobs(0)) &
    print dafmt,"Analysis of "//trim(obsData(ipar)%tag)
  if(debug.or.nn==obsData(ipar)%nobs) print '(5(A,I0))',&
    'Obs(',me,'):',n,'/',nobs,',',nn,'/',obsData(ipar)%nobs
!-----------------------------------------------------------------------
! analysis of direct observations: model mid-levels
!-----------------------------------------------------------------------
  if(obsData(ipar)%found.and.obsData(ipar)%deriv=='mod-lev')then
    unitconv=obsData(ipar)%unitconv
    if(obsData(ipar)%unitroa)&
      H_int(:)=H_int(:)*(/(roa(i(p),j(p),ll,1),p=1,4)/)
    k =obsData(ipar)%ichem
    kk=obsData(ipar)%ichemObs!ichemobs(k)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    do p=1,4
      H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
        * unitconv                                    ! unit conversion
      xn=xn_adv_ex(i(p),j(p),l,kk)
      yn=yn+H_jac(n,p,l,k)*xn
    enddo
!-----------------------------------------------------------------------
! analysis of direct observations: surface level
!-----------------------------------------------------------------------
  elseif(obsData(ipar)%found.and.obsData(ipar)%deriv=='surface')then
    unitconv=obsData(ipar)%unitconv
    if(obsData(ipar)%unitroa)&
      H_int(:)=H_int(:)*(/(roa(i(p),j(p),ll,1),p=1,4)/)
    ispec=obsData(ipar)%ispec
    k =obsData(ipar)%ichem
    kk=obsData(ipar)%ichemObs!ichemobs(k)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    do p=1,4
      H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
        * unitconv                                  & ! unit conversion
        * cfac(ispec,i(p),j(p))                       ! 50m->3m conversion
      xn=xn_adv_ex(i(p),j(p),l,kk)
      yn=yn+H_jac(n,p,l,k)*xn
    enddo
!-----------------------------------------------------------------------
! analysis of indirect observations: COLUMN (eg NO2tc)
!-----------------------------------------------------------------------
  elseif(obsData(ipar)%found.and.obsData(ipar)%deriv=='Trop.Col.')then
    unitconv=obsData(ipar)%unitconv
    ispec=obsData(ipar)%ispec
    k=obsData(ipar)%ichem
    kk=ichemobs(k)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=(/1,nlev/)
    do ll=1,KMAX_MID-nlev   ! B might have a reduced number of levels...
      do p=1,4
        Hj=H_int(p)                                   & ! interpolation
          * unitconv                                  & ! unit conversion
          * roa(i(p),j(p),ll,1)                       & ! density.
          * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
        xn=xn_adv(ispec,i(p),j(p),ll)*FGSCALE
        yn=yn+Hj*xn
      enddo
    enddo
    do l=1,nlev
      ll=KMAX_MID-nlev+l
      do p=1,4
        H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
          * unitconv                                  & ! unit conversion
          * roa(i(p),j(p),ll,1)                       & ! density.
          * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
        xn=xn_adv_ex(i(p),j(p),l,kk)
        yn=yn+H_jac(n,p,l,k)*xn
      enddo
    enddo
!-----------------------------------------------------------------------
! analysis of indirect observations: model mid-levels
!  PM* (eg PM10), BSC and EXT GROUPs (wavelength is in the %unit)
!-----------------------------------------------------------------------
  elseif((index(obsData(ipar)%name,'PM')>0.or.&
          any(obsData(ipar)%name==(/'BSC','EXT'/)))&
         .and.obsData(ipar)%deriv=='mod-lev')then
!   call CheckStop(.not.obsData(ipar)%unitroa.or.&
!     (index(obsData(ipar)%name,'PM')>0.and.obsData(ipar)%unit(1:2)/='ug'),&
!     "Unsupported obsData unit for: "//trim(obsData(ipar)%tag))
    igrp=find_index(obsData(ipar)%name,chemgroups(:)%name)
    if(igrp<1)then
      if(MasterProc) print dafmt,"ERROR obsData not yet supported"
      return
    endif
    call Group_Units(igrp,obsData(ipar)%unit,gspec,gunit_conv,debug)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    do g=1,size(gspec)
      kk=varSpecInv(gspec(g))
      if(kk>0.and.observedVar(kk))then
        k=ichemInv(kk)
        do p=1,4
          H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
            * gunit_conv(g)                             & ! unit conversion
            * roa(i(p),j(p),ll,1)                         ! density.
          xn=xn_adv_ex(i(p),j(p),l,kk)
          yn=yn+H_jac(n,p,l,k)*xn
        enddo
      else
        do p=1,4
          Hj=H_int(p)                                   & ! interpolation
            * gunit_conv(g)                             & ! unit conversion
            * roa(i(p),j(p),ll,1)                         ! density.
          xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
          yn=yn+Hj*xn
        enddo
      endif
    enddo
    if(index(obsData(ipar)%name,'PM')>0)then
      do p=1,4
        yn=yn+PM25_water(i(p),j(p),ll)*FGSCALE ! PM water content in ug/m3 at model mid-levels
      enddo
    endif
!-----------------------------------------------------------------------
! analysis of indirect observations: model mid-levels
!  PM* (eg PM10), BSC and EXT GROUPs (wavelength is in the %unit)
!-----------------------------------------------------------------------
  elseif((index(obsData(ipar)%name,'PM')>0.or.&
          any(obsData(ipar)%name==(/'BSC','EXT'/)))&
         .and.obsData(ipar)%deriv=='surface')then
!   call CheckStop(.not.obsData(ipar)%unitroa.or.&
!     (index(obsData(ipar)%name,'PM')>0.and.obsData(ipar)%unit(1:2)/='ug'),&
!     "Unsupported obsData unit for: "//trim(obsData(ipar)%tag))
    igrp=find_index(obsData(ipar)%name,chemgroups(:)%name)
    if(igrp<1)then
      if(MasterProc) print dafmt,"ERROR obsData not yet supported"
      return
    endif
    print *,igrp,chemgroups(igrp)%name
    call Group_Units(igrp,obsData(ipar)%unit,gspec,gunit_conv,debug)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    do g=1,size(gspec)
      kk=varSpecInv(gspec(g))
      if(kk>0.and.observedVar(kk))then
        k=ichemInv(kk)
        do p=1,4
          H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
            * gunit_conv(g)                             & ! unit conversion
            * roa(i(p),j(p),ll,1)                       & ! density.
            * cfac(gspec(g),i(p),j(p))                    ! 50m->3m conversion
          xn=xn_adv_ex(i(p),j(p),l,kk)
          yn=yn+H_jac(n,p,l,k)*xn
        enddo
      else
        do p=1,4
          Hj=H_int(p)                                   & ! interpolation
            * gunit_conv(g)                             & ! unit conversion
            * roa(i(p),j(p),ll,1)                       & ! density.
            * cfac(gspec(g),i(p),j(p))                    ! 50m->3m conversion
          xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
          yn=yn+Hj*xn
        enddo
      endif
    enddo
    if(index(obsData(ipar)%name,'PM')>0)then
      do p=1,4
        yn=yn+PM25_water_rh50(i(p),j(p))*FGSCALE ! PM water content in ug/m3 at surface level
      enddo
    endif
!-----------------------------------------------------------------------
! analysis of indirect observations: column
!  AOD GROUP (wavelength is in the %unit)
!-----------------------------------------------------------------------
  elseif((obsData(ipar)%name=='AOD').or.&
         (obsData(ipar)%name=='EXT'.and.obsData(ipar)%deriv=='Trop.Col.'))then
    igrp=find_index(obsData(ipar)%name,chemgroups(:)%name)
    if(igrp<1)then
      if(MasterProc) print dafmt,"ERROR obsData not yet supported"
      return
    endif
    call Group_Units(igrp,obsData(ipar)%unit,gspec,gunit_conv,debug)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    pidx(n)%l(:)=(/1,nlev/)
    do g=1,size(gspec)
      kk=varSpecInv(gspec(g))
      if(kk>0.and.observedVar(kk))then
        k=ichemInv(kk)
        do ll=1,KMAX_MID-nlev
          do p=1,4
            Hj=H_int(p)                                   & ! interpolation
              * gunit_conv(g)                             & ! unit conversion
              * roa(i(p),j(p),ll,1)                       & ! density.
              * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
            xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
            yn=yn+Hj*xn
          enddo
        enddo
        do l=1,nlev
          ll=KMAX_MID-nlev+l
          do p=1,4
            H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
              * gunit_conv(g)                             & ! unit conversion
              * roa(i(p),j(p),ll,1)                       & ! density.
              * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
            xn=xn_adv_ex(i(p),j(p),l,kk)
            yn=yn+H_jac(n,p,l,k)*xn
          enddo
        enddo
      else
        do ll=1,KMAX_MID
          do p=1,4
            Hj=H_int(p)                                   & ! interpolation
              * gunit_conv(g)                             & ! unit conversion
              * roa(i(p),j(p),ll,1)                       & ! density.
              * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
            xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
            yn=yn+Hj*xn
          enddo
        enddo
      endif
    enddo
  else
    if(MasterProc) print dafmt,"ERROR obsData not yet supported"
    return
  endif
end subroutine H_op
!-----------------------------------------------------------------------
! @description
! Compute CHI^2/nobs on an observation by observation
! basis by taken a subset of the total number of observations.
! This should be of order 0.5.
! Calculated separately for each observation type.
! @author M.Kahnert (SMHI) and L.Robertson (SMHI)
!-----------------------------------------------------------------------
subroutine chisq_over_nobs2(nobs,dep)
  integer :: nobs                    ! number of observations,...
  real    :: dep(nobs)               ! Analysis departure from background
  integer :: i,j,n,p,frac,ipar
  real :: hbh
  integer, allocatable :: nm(:)      ! number of used observations
  real, allocatable :: chisq(:),&    ! chi-square
                       h(:,:,:,:)    ! Work array for observation operator
  complex, allocatable :: uh(:,:,:)  ! .. on spectral space
  character(len=128) :: file
!-----------------------------------------------------------------------
  if((nobsData<1).or.(nobs<1))then
    if(MasterProc) print dafmt,'No Obs for CHI^2 calculations'
    return
  endif
!----------------------------------------------------------------------- 
  if(MasterProc) print dafmt,'CHI^2 calculations' 
  allocate(chisq(nobsData),nm(nobsData),&
           h(nxex,nyex,nlev,nchemobs),uh(nv1,nxex,nyex))
!-----------------------------------------------------------------------
! Determine (H^T)_n; compute U^{-\dagger}*(H^T)_n, n: nth component of H^T, n=1..nobs
!-----------------------------------------------------------------------
  chisq(:)=0e0
  nm(:)=0
! frac=nobs/max(nobs/5,min(200,nobs)) ! 1 (up to 200 obs) .. 5 (from 1000 obs)
  do ipar=1,nobsData
    n=obsData(ipar)%nobs
    if(n==0)cycle
    frac=n/max(n/5,min(200,n))        ! 1 (up to 200 obs) .. 5 (from 1000 obs)
    do n=obsData(ipar)%iobs(0)+frac,& ! only use every frac'th observation, i.e.
         obsData(ipar)%iobs(1),frac   ! between all obs (up to 200) and 1/5th
      h(:,:,:,:)=0e0 ! observation operator for this observations
      do p=1,4
        i = pidx(n)%i(p)
        j = pidx(n)%j(p)
        h(i,j,:,:)=H_jac(n,p,:,ichemobs(:))
      enddo
      call chitox_adj(uh,h)                 ! to spectral space
      hbh=norm(uh,squared=.true.)+obsstddev(n)**2  ! HBH+R
      chisq(ipar)=chisq(ipar)+dep(n)**2/hbh ! chi-square for this observation & sum
      nm(ipar)=nm(ipar)+1                   ! number of used observations
    enddo
  enddo
  deallocate(h,uh)
! CALL MPI_ALLREDUCE(MPI_IN_PLACE,chisq,ipar,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,INFO)
! CALL MPI_ALLREDUCE(MPI_IN_PLACE,nm   ,ipar,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,INFO)
  where(nm(:)>0) chisq(:)=chisq(:)/nm(:)     ! Normalize by number of used observations

  if(MasterProc) print '((X,I3,A,X,G12.4,X,3(A,I0)))',&
    (ipar,'CHI^2/nobs=',chisq(ipar),'nobs=',nm(ipar),&
      '/',obsData(ipar)%nobs,'/',nobs,ipar=1,nobsData)

  file=trim(runlabel1)//"_chisqr.dat"
  open(IO_TMP,file=file,form="FORMATTED",position="APPEND")
  write(IO_TMP,"(A,';',I3,50(';',G12.4,';',I7,:))"),&
    date2string("YYYY-MM-DD hh:mm",current_date),nobsData,&
    (chisq(ipar),nm(ipar),ipar=1,nobsData)
  close(IO_TMP)
  deallocate(chisq,nm)
endsubroutine chisq_over_nobs2
subroutine get_interpol_const(lon,lat,IIij,JJij,Weight,method)
! find the four closest points & calulcate distance weights
  real,    intent(in)         :: lon,lat
  integer, intent(out)        :: IIij(4),JJij(4)
  real,    intent(out)        :: Weight(4)
  character(len=*),intent(in) :: method
!!logical, parameter :: debug=.true.
  integer :: n

  select case(method)
  case('old','local')
    call get_coeff()
  case('dis','nearest4','distance-weighted') ! 4-point distance-weighted
    call point2grid_coeff(lon,lat,IIij,JJij,Weight,&
      glon(:limax,:ljmax),glat(:limax,:ljmax),limax,ljmax,debug)
  case('nnn','nearest1','nearest-neighbor')  ! 1-point nearest-neighbor
    call point2grid_coeff(lon,lat,IIij,JJij,Weight,&
      glon(:limax,:ljmax),glat(:limax,:ljmax),limax,ljmax,debug)
    Weight(:)=(/1.0,0.0,0.0,0.0/) ! Only use nearest grid cell
    IIij(2:4)=IIij(1)             ! Reset unused indices, 
    JJij(2:4)=JJij(1)             ! these elements have weight zero.
  case default  
    call CheckStop("get_interpol_const: "//&
      "Unknown interpolation method '"//trim(method)//"'.")
  endselect
  if(debug) print "(A,1X,A,4(1X,I0,1X,I0,1X,F5.3,','))",&
    'get_interpol_const',trim(method),(IIij(n),JJij(n),Weight(n),n=1,4)
contains
subroutine get_coeff()
  integer :: i,j,n
  real :: d,dist(4)
  dist=1.0E40
  do j=1,ljmax
    do i=1,limax
      !distance between (lon(i,j),lat(i,j) and (lon,lat)
      d=great_circle_distance(lon,lat,glon(i,j),glat(i,j))
      if(d>=dist(4))cycle
      n=MINVAL([1,2,3,4],MASK=(d<dist))
      dist(n:4)=EOSHIFT(dist(n:4),-1,BOUNDARY=d)
      IIij(n:4)=EOSHIFT(IIij(n:4),-1,BOUNDARY=i)
      JJij(n:4)=EOSHIFT(JJij(n:4),-1,BOUNDARY=j)
    enddo
  enddo
  Weight(1)=1.0-3.0*dist(1)/sum(dist(1:4))
  Weight(2)=(1.0-Weight(1))*(1.0-2.0*dist(2)/sum(dist(2:4)))
  Weight(3)=(1.0-Weight(1)-Weight(2))*(1.0-dist(3)/sum(dist(3:4)))
  Weight(4)=1.0-Weight(1)-Weight(2)-Weight(3)
endsubroutine get_coeff
endsubroutine get_interpol_const
endmodule DA_Obs_ml
