#define STRING2(x) #x
#define STRING(x) STRING2(x)
#define HERE(MSG) MSG//" ("//__FILE__//":"//STRING(__LINE__)//")."
module DA_Obs_ml
use MetFields_ml,     only: z_bnd,z_mid,roa,th,q,ps
use TimeDate_ml,      only: current_date
use CheckStop_ml,     only: CheckStop
use Functions_ml,     only: great_circle_distance
use GridValues_ml,    only: glon,glat,glon_fdom,glat_fdom,&
                            i_local,j_local,coord_in_domain
use Par_ml,           only: limax,ljmax,me
use Chemfields_ml,    only: xn_adv,cfac,PM25_water,PM25_water_rh50
!se ChemSpecs_shl_ml, only: NSPEC_SHL        ! Maps indices
use ChemSpecs_adv_ml, only: NSPEC_ADV
use ChemChemicals_ml, only: species
use ChemGroups_ml,    only: chemgroups
use Io_ml,            only: IO_TMP,PrintLog
use InterpolationRoutines_ml,  only : point2grid_coeff
use ModelConstants_ml,only: KMAX_MID,IIFULLDOM,JJFULLDOM,MasterProc,runlabel1
use SmallUtils_ml,    only: find_index
use Functions_ml,     only: norm
use TimeDate_ExtraUtil_ml, only: date2string
use Units_ml,         only: Units_Scale,Group_Units
use DA_ml,            only: debug=>DA_DEBUG,debug_obs=>DA_DEBUG_OBS,&
                            dafmt=>da_fmt_msg,damsg=>da_msg
use spectralcov,      only: nx,ny,nlev,nchem,nxex,nyex,nv1,FGSCALE,&
                            nChemObs,nChemNoObs,iChemInv,iChemObs,iChemNoObs
use chitox_ml,        only: matched_domain,iyf,chitox_adj
use mpi,              only: MPI_COMM_WORLD,MPI_ALLREDUCE,MPI_SUM,&
                            MPI_IN_PLACE,MPI_INTEGER
implicit none
integer, save :: nobs!,matidx
real, dimension(:), allocatable, save :: innov,obsstddev
real, dimension(:,:,:,:), allocatable, save :: H_jac
type H_jac_idx
  integer :: i(4),j(4),l(0:1)
  logical :: local(4)
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
  logical, parameter :: F=.false.
  integer :: ierr
  if(.not.allocated(innov))then
    allocate(innov(nobs),stat=ierr)
    call CheckStop(ierr,HERE('Allocate INNOV'))
    innov=0.0
  endif
  if(.not.allocated(obsstddev))then
    allocate(obsstddev(nobs),stat=ierr)
    call CheckStop(ierr,HERE('Allocation OBSSTDDEV'))
    obsstddev=0.0
  endif
!-----------------------------------------------------------------------
! H_jac is a matrix with components H_{n;i,j,l,k}, where
!   n=1,...,nobs     (observation index)
!   i=1,...,nx       (longitude index in full/global domain)
!   j=1,...,ny       (latitude  index in full/global domain)
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
    call CheckStop(ierr,HERE('Allocate H_JAC'))
    H_jac=0.0
  endif
  if(.not.allocated(pidx))then
    allocate(pidx(nobs),stat=ierr)
    call CheckStop(ierr,HERE('Allocate PIDX'))
    pidx=H_jac_idx([0,0,0,0],[0,0,0,0],[0,0],[F,F,F,F])
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
subroutine read_obs(domain,maxobs,flat,flon,falt,y,stddev,ipar)
!-----------------------------------------------------------------------
! @description
! Read in observations y in standard format
! (Station-No flat flon falt obs-value stddev).
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none

  character(len=*), intent(in) :: domain
  integer :: maxobs,ipar(maxobs)
  real :: y(maxobs),stddev(maxobs),flat(maxobs),flon(maxobs),falt(maxobs)

  integer            :: no,nd,ierr
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
  do nd=1,nobsData
    if(obsData(nd)%file=="")then
      write(damsg,"(A,I0,1X,A)")'WARNING: obsData#',nd,"file not specified"
      if(MasterProc) print dafmt,damsg
      cycle
    endif
    file=date2string(obsData(nd)%file,current_date)
    inquire(file=file,exist=exist)
    if(.not.exist)then
      write(damsg,"(A,I0,1X,A)")'WARNING: obsData#',nd,"file not found"
      if(MasterProc) print dafmt,damsg
      cycle
    endif
    open(IO_TMP,file=file,form='formatted',status='old',iostat=ierr)
    if(ierr/=0)then
      write(damsg,"(A,I0,1X,A)")'WARNING: obsData#',nd,"file not opened"
      if(MasterProc) print dafmt,damsg
      close(IO_TMP)
      cycle
    endif
    print dafmt,"Obsdata opened "//trim(file)
!-----------------------------------------------------------------------
! read data
!-----------------------------------------------------------------------
    obsData(nd)%iobs(0)=nobs+1
    do while(ierr==0)
      if(nobs>=maxobs)then
        print dafmt,'ERROR reading observations, maxobs too small'
        ierr=-2
        cycle
      endif
      nobs=nobs+1
      ipar(nobs)=nd
      read(IO_TMP,*,iostat=ierr)no,flat(nobs),flon(nobs),falt(nobs),&
                                   y(nobs),stddev(nobs)
      select case(ierr)
      case(0)
        if((.not.coord_in_domain(domain,flon(nobs),flat(nobs))).or.&
           (y(nobs)<obsData(nd)%min).or.&
           (y(nobs)>obsData(nd)%max))then
         nobs=nobs-1
         cycle
        endif
        stddev(nobs)=max(stddev(nobs),&
          y(nobs)*obsData(nd)%error_rel,&
          obsData(nd)%error_rep)
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
    obsData(nd)%iobs(1)=nobs
  enddo
!-----------------------------------------------------------------------
! Setup obsData &  Write #obs to log file
!-----------------------------------------------------------------------
  do nd=1,nobsData
    if(.not.obsData(nd)%set)then
      no=find_index(obsData(nd)%name,obsVarName(:))
      obsData(nd)%found=(no>0)
      if(obsData(nd)%found)then
        obsData(nd)%ichem=no                    ! index observed
        obsData(nd)%ichemObs=ichemobs(no)       ! index observed/unobserved
        obsData(nd)%ispec=varSpec(ichemobs(no)) ! index species
        obsData(nd)%unitconv=Units_Scale(obsData(nd)%unit,&
          obsData(nd)%ispec,needroa=obsData(nd)%unitroa)
      else
        obsData(nd)%ichem=find_index(obsData(nd)%name,chemgroups(:)%name)
        obsData(nd)%unitconv=Units_Scale(obsData(nd)%unit,&
          -1,needroa=obsData(nd)%unitroa)
      endif
      call CheckStop(obsData(nd)%deriv=='',&
        HERE("Inconsistent obsData%deriv"))
      write(obsData(nd)%tag,"(2(A,1X),'(',A,')')")&
        trim(obsData(nd)%name),trim(obsData(nd)%deriv),&
        trim(obsData(nd)%unit)
      call CheckStop(obsData(nd)%ichem<1,&
        HERE("Unsupported obsData "//trim(obsData(nd)%tag)))
      obsData(nd)%set=.true.
    endif
    if(any(obsData(nd)%iobs/=0)) obsData(nd)%nobs=&
      DIM(obsData(nd)%iobs(1)+1,obsData(nd)%iobs(0))
    call CheckStop(obsData(nd)%nobs,count(ipar(:nobs)==nd),&
      HERE("Inconsistent obsData%nobs"))
    if(MasterProc)then
      file=date2string(obsData(nd)%file,current_date)
      write(damsg,"('obsData(',I0,') contains ',I0,1X,A,' observations')")&
        nd,obsData(nd)%nobs,trim(obsData(nd)%tag)
      write(damsg,dafmt)trim(damsg)
      call PrintLog(damsg)
      call PrintLog(file)
    endif
  enddo
endsubroutine read_obs
subroutine H_op(lat,lon,alt,n,yn,ipar)
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

  integer :: i(4),j(4),p,l,ll,k,kk,igrp,g,ispec,nn,INFO
  real :: xn,Hj,H_int(4),T1,QML1,PS1,Hchem(nchem),unitconv
  integer, pointer, dimension(:) :: gspec=>null()      ! group array of indexes
  real,    pointer, dimension(:) :: gunit_conv=>null() ! group array of unit conv. factors
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  nn=n-obsData(ipar)%iobs(0)+1
  if(n==obsData(ipar)%iobs(0)) &
    print dafmt,"Analysis of "//trim(obsData(ipar)%tag)
  if(debug_obs.or.(nn==obsData(ipar)%nobs.and.MasterProc)) print '(5(A,I0))',&
    'Obs(',me,'):',n,'/',nobs,',',nn,'/',obsData(ipar)%nobs
!-----------------------------------------------------------------------
! operator for horizontal four-point interpolation between nearest
! neighbouring grid points:
!-----------------------------------------------------------------------
  call get_interpol_const(lon,lat,i,j,H_int,interpolation_method)
  pidx(n)%i=i;i=i_local(i)  ! save global indexing
  pidx(n)%j=j;j=j_local(j)  ! and use local indexing
  pidx(n)%local=(H_int>0.0).and.(i>=1).and.(i<=limax).and.(j>=1).and.(j<=ljmax)
!-----------------------------------------------------------------------
! z grid coordinate of observations:
!-----------------------------------------------------------------------
  if(pidx(n)%local(1))then  ! p:=1 has H_int maxval
    ll=KMAX_MID             ! calcualte level on p:=1
    do while(alt>z_bnd(i(1),j(1),ll).and.ll>=1)
      ll=ll-1
    enddo
  else
    ll=0
  endif ! share ll from p:=1 host CPU to all CPUs
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,ll,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,INFO)
  l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
  if((l<1).or.(debug_obs.and.any(pidx(n)%local))) then
    print dafmt,'Observation/Model Location:Weight[%],lon[deg],lat[deg],alt[km]'
    print "('#',I0,5(1X,A3,':',F0.1,3(',',F0.2)))",&
      n,'Observation',100.0,lon,lat,alt,&
       ('Model',H_int(p)*100,glon(i(p),j(p)),glat(i(p),j(p)),&
                z_mid(i(p),j(p),ll)*1e-3,p=1,4)
    call CheckStop(l<1,HERE("Out of bounds level"))
  endif
!-----------------------------------------------------------------------
! analysis of direct observations: NO2, O3, SO2, etc
!-----------------------------------------------------------------------
  if(obsData(ipar)%found)then
    unitconv=obsData(ipar)%unitconv
    ispec=obsData(ipar)%ispec
    k =obsData(ipar)%ichem
    yn=0e0
    select case(obsData(ipar)%deriv)
!-----------------------------------------------------------------------
! direct observations: model mid-levels
!-----------------------------------------------------------------------
    case('mod-lev')
      pidx(n)%l(:)=l
      if(obsData(ipar)%unitroa)&
        forall(p=1:4,pidx(n)%local(p)) &
          H_int(p)=H_int(p)*roa(i(p),j(p),ll,1)
      do p=1,4
        if(.not.pidx(n)%local(p))cycle
        Hj=H_int(p)*unitconv                       ! interpolation,unit conversion        
        xn=xn_adv(ispec,i(p),j(p),ll)*FGSCALE
        yn=yn+Hj*xn
        H_jac(n,p,l,k)=Hj
      enddo
!-----------------------------------------------------------------------
! direct observations: surface level
!-----------------------------------------------------------------------
    case('surface')
      pidx(n)%l(:)=l
      if(obsData(ipar)%unitroa)&
        forall(p=1:4,pidx(n)%local(p)) &
          H_int(p)=H_int(p)*roa(i(p),j(p),ll,1)
      do p=1,4
        if(.not.pidx(n)%local(p))cycle
        Hj=H_int(p)*unitconv*cfac(ispec,i(p),j(p))  ! interpolation,unit conversion,50m->3m
        xn=xn_adv(ispec,i(p),j(p),ll)*FGSCALE
        yn=yn+Hj*xn
        H_jac(n,p,l,k)=Hj
      enddo
!-----------------------------------------------------------------------
! direct observations: COLUMN (eg NO2tc)
!-----------------------------------------------------------------------
    case('Trop.Col.')
      pidx(n)%l(:)=(/1,nlev/)
      call CheckStop(.not.obsData(ipar)%unitroa,&
        "Unsupported obsData unit for: "//trim(obsData(ipar)%tag))     
      do p=1,4
        if(.not.pidx(n)%local(p))cycle
        do ll=1,KMAX_MID
          Hj=H_int(p)*unitconv            & ! interpolation,unit conversion
            * roa(i(p),j(p),ll,1)                       & ! density.
            * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
          xn=xn_adv(ispec,i(p),j(p),ll)*FGSCALE
          yn=yn+Hj*xn
          l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
          if(l>0)H_jac(n,p,l,k)=Hj
      enddo
    enddo
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    case default
      call CheckStop(HERE("obsData not yet supported "//trim(obsData(ipar)%tag)))
    endselect
!-----------------------------------------------------------------------
! analysis of indirect observations:
!   PM* (eg PM10), BSC and EXT GROUPs (wavelength is in %unit)
!-----------------------------------------------------------------------
  elseif(index(obsData(ipar)%name,'PM')>0.or.&
         any(obsData(ipar)%name==['BSC','EXT','AOD']))then
    call CheckStop(any(obsData(ipar)%name==['BSC','EXT','AOD']),&
      HERE("obsData not yet supported "//trim(obsData(ipar)%tag)))
!   call CheckStop(.not.obsData(ipar)%unitroa.or.&
!     (index(obsData(ipar)%name,'PM')>0.and.obsData(ipar)%unit(1:2)/='ug'),&
!     "Unsupported obsData unit for: "//trim(obsData(ipar)%tag))
    call CheckStop((obsData(ipar)%name=='AOD').and.(obsData(ipar)%deriv/='Trop.Col.'),&
      HERE("Unsupported obsData "//trim(obsData(ipar)%tag)))
    igrp=obsData(ipar)%ispec
!!  print *,igrp,chemgroups(igrp)%name
    call Group_Units(igrp,obsData(ipar)%unit,gspec,gunit_conv,debug)
    yn=0e0
    select case(obsData(ipar)%deriv)
!-----------------------------------------------------------------------
! indirect observations: model mid-levels
!-----------------------------------------------------------------------
    case('mod-lev')
      if(index(obsData(ipar)%name,'PM')>0)then
        do p=1,4 ! PM water content in ug/m3 at model mid-levels
          if(.not.pidx(n)%local(p))cycle
          yn=yn+H_int(p)*PM25_water(i(p),j(p),ll)*FGSCALE
        enddo
      endif
      pidx(n)%l(:)=l
      if(obsData(ipar)%unitroa)&
        forall(p=1:4,pidx(n)%local(p)) &
          H_int(p)=H_int(p)*roa(i(p),j(p),ll,1)
      do p=1,4
        if(.not.pidx(n)%local(p))cycle
        do g=1,size(gspec)
          Hj=H_int(p)*gunit_conv(g)       ! interpolation,unit conversion
          xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
          yn=yn+Hj*xn
          kk=varSpecInv(gspec(g))
          k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
          if(k>0)H_jac(n,p,l,k)=Hj
        enddo
      enddo
!-----------------------------------------------------------------------
! indirect observations: surface level
!-----------------------------------------------------------------------
    case('surface')
      if(index(obsData(ipar)%name,'PM')>0)then
        do p=1,4 ! PM water content in ug/m3 at surface level
          if(.not.pidx(n)%local(p))cycle
          yn=yn+H_int(p)*PM25_water_rh50(i(p),j(p))*FGSCALE
        enddo
      endif
      pidx(n)%l(:)=l
      if(obsData(ipar)%unitroa)&
        forall(p=1:4,pidx(n)%local(p)) &
          H_int(p)=H_int(p)*roa(i(p),j(p),ll,1)
      do p=1,4
        if(.not.pidx(n)%local(p))cycle
        do g=1,size(gspec)
          Hj=H_int(p)*gunit_conv(g)     & ! interpolation,unit conversion
            * cfac(gspec(g),i(p),j(p))    ! 50m->3m conversion
          xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
          yn=yn+Hj*xn
          kk=varSpecInv(gspec(g))
          k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
          if(k>0)H_jac(n,p,l,k)=Hj
        enddo
      enddo
!-----------------------------------------------------------------------
! direct observations: COLUMN (eg AOD)
!-----------------------------------------------------------------------
    case('Trop.Col.')
      if(index(obsData(ipar)%name,'PM')>0)then
        do p=1,4 ! PM water content in ug/m3 at model mid-levels
          if(.not.pidx(n)%local(p))cycle
          xn=0.0
          do ll=1,KMAX_MID
            xn=xn+PM25_water(i(p),j(p),ll) &
              * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
          enddo
          yn=yn+H_int(p)*xn*FGSCALE
        enddo
      endif
      pidx(n)%l(:)=(/1,nlev/)
      do p=1,4
        if(.not.pidx(n)%local(p))cycle
        do ll=1,KMAX_MID
          Hj=H_int(p)*roa(i(p),j(p),ll,1)           & ! interpolation,density
              * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
          l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
          do g=1,size(gspec)
            xn=xn_adv(gspec(g),i(p),j(p),ll)*gunit_conv(g)*FGSCALE
            yn=yn+Hj*xn
            kk=varSpecInv(gspec(g))
            k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
            if(k>0.and.l>0)H_jac(n,p,l,k)=Hj*gunit_conv(g)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    case default
      call CheckStop(HERE("obsData not yet supported "//trim(obsData(ipar)%tag)))
    endselect
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  else
    call CheckStop(HERE("obsData not yet supported "//trim(obsData(ipar)%tag)))
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
  integer :: i,j,n,p,frac,ipar,iyf0,iyf1
  real :: hbh
  integer, allocatable :: nm(:)      ! number of used observations
  real, allocatable :: chisq(:),&    ! chi-square
                   h_loc(:,:,:,:)    ! Work array for observation operator
  complex, allocatable :: uh(:,:,:)  ! .. on spectral space
  character(len=128) :: file
!-----------------------------------------------------------------------
  if((nobsData<1).or.(nobs<1))then
    if(MasterProc) print dafmt,'No Obs for CHI^2 calculations'
    return
  endif
!----------------------------------------------------------------------- 
  if(MasterProc) print dafmt,'CHI^2 calculations' 
  iyf0=iyf(0,me);iyf1=iyf(1,me)       ! local y-slab
  allocate(chisq(nobsData),nm(nobsData),h_loc(nxex,iyf0:iyf1,nlev,nchemobs),&
           uh(nv1,nxex,nyex))
  chisq(:)=0e0
  nm(:)=0
  h_loc(:,:,:,:)=0e0
!-----------------------------------------------------------------------
! Determine (H^T)_n; compute U^{-\dagger}*(H^T)_n, n: nth component of H^T, n=1..nobs
!-----------------------------------------------------------------------
  do ipar=1,nobsData
    n=obsData(ipar)%nobs
    if(n==0)cycle
    frac=n/max(n/5,min(200,n))        ! 1 (up to 200 obs) .. 5 (from 1000 obs)
    do n=obsData(ipar)%iobs(0)+frac,& ! only use every frac'th observation, i.e.
         obsData(ipar)%iobs(1),frac   ! between all obs (up to 200) and 1/5th
      do p=1,4                        ! observation operator for this observation
        i = pidx(n)%i(p)
        j = pidx(n)%j(p)
        if((j>=iyf0).and.(j<=iyf1)) & ! local y-slab
          h_loc(i,j,:,:)=H_jac(n,p,:,ichemobs(:))
      enddo
      call chitox_adj(uh,h_loc)             ! to spectral space
      hbh=norm(uh,squared=.true.)+obsstddev(n)**2  ! HBH+R
      chisq(ipar)=chisq(ipar)+dep(n)**2/hbh ! chi-square for this observation & sum
      nm(ipar)=nm(ipar)+1                   ! number of used observations
     h_loc(pidx(n)%i,pidx(n)%j,:,:)=0e0     ! clean up for next iteration
    enddo
  enddo
  deallocate(h_loc,uh)
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
  logical, parameter :: debug=.false.
  integer :: n

  select case(method)
  case('old','local')
    call get_coeff()
  case('dis','nearest4','distance-weighted') ! 4-point distance-weighted
    call point2grid_coeff(lon,lat,IIij,JJij,Weight,&
      glon_fdom,glat_fdom,IIFULLDOM,JJFULLDOM,debug)
  case('nnn','nearest1','nearest-neighbor')  ! 1-point nearest-neighbor
    call point2grid_coeff(lon,lat,IIij,JJij,Weight,&
      glon_fdom,glat_fdom,IIFULLDOM,JJFULLDOM,debug)
    Weight(:)=(/1.0,0.0,0.0,0.0/) ! Only use nearest grid cell
    IIij(2:4)=IIij(1)             ! Reset unused indices, 
    JJij(2:4)=JJij(1)             ! these elements have weight zero.
  case default  
    call CheckStop("get_interpol_const: "//&
      HERE("Unknown interpolation method '"//trim(method)))
  endselect
  if(debug) print "(A,1X,A,4(1X,I0,1X,I0,1X,F5.3,','))",&
    'get_interpol_const',trim(method),(IIij(n),JJij(n),Weight(n),n=1,4)
contains
subroutine get_coeff()
  integer :: i,j,n
  real :: d,dist(4)
  dist=1.0E40
  do j=1,JJFULLDOM
    do i=1,IIFULLDOM
      !distance between (lon(i,j),lat(i,j) and (lon,lat)
      d=great_circle_distance(lon,lat,glon_fdom(i,j),glat_fdom(i,j))
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
