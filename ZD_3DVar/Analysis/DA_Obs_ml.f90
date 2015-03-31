#define STRING2(x) #x
#define STRING(x) STRING2(x)
#define HERE(MSG) MSG//" ("//__FILE__//":"//STRING(__LINE__)//")."
module DA_Obs_ml
use MetFields_ml,     only: z_bnd,z_mid,roa,th,q,ps
use TimeDate_ml,      only: current_date
use CheckStop_ml,     only: CheckStop
use Functions_ml,     only: great_circle_distance
use GridValues_ml,    only: glon,glat,coord_in_domain
use Par_ml,           only: me
use Chemfields_ml,    only: xn_adv,cfac,PM25_water,PM25_water_rh50
!se ChemSpecs_shl_ml, only: NSPEC_SHL        ! Maps indices
use ChemSpecs_adv_ml, only: NSPEC_ADV
use ChemChemicals_ml, only: species
use ChemGroups_ml,    only: chemgroups
use Io_ml,            only: IO_TMP,PrintLog
use ModelConstants_ml,only: KMAX_MID,IIFULLDOM,JJFULLDOM,MasterProc,runlabel1
use SmallUtils_ml,    only: find_index
use Functions_ml,     only: norm
use TimeDate_ExtraUtil_ml, only: date2string
use Units_ml,         only: Units_Scale,Group_Units
use DA_ml,            only: debug=>DA_DEBUG,dafmt=>da_fmt_msg,damsg=>da_msg,&
                            debug_obs=>DA_DEBUG_OBS
use spectralcov,      only: nx,ny,nlev,nchem,nxex,nyex,nv1,FGSCALE,&
                            nChemObs,nChemNoObs,iChemInv,iChemObs,iChemNoObs
use chitox_ml,        only: matched_domain,iyf,chitox_adj
use mpi,              only: MPI_COMM_WORLD,MPI_ALLREDUCE,MPI_SUM,&
                            MPI_IN_PLACE,MPI_INTEGER
implicit none
integer, save :: nobs!,matidx
real, dimension(:), allocatable, save :: innov,obsstddev
real, dimension(:,:,:), allocatable, save :: H_jac
type H_jac_idx
  integer :: i,j,l0,l1
  logical :: in_mgrid,in_yslab
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
namelist /OBSERVATIONS/nobsData,obsData
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
! The Jacobian is highly sparse wrt i,j,l, so the matrix is stored using
! the following compact storage scheme:
! H_{n;i,j,l,k}=H_jac(n,l,k), where n=1:nobs,l=l0:l1
!   pidx(n)%i=i_n,   pidx(n)%j=j_n,   pidx(n)%l0=l0_n,  pidx(n)%l1=l1_n,
!-----------------------------------------------------------------------
  if(.not.allocated(H_jac))then
    allocate(H_jac(nobs,nlev,nchemobs),stat=ierr)
    call CheckStop(ierr,HERE('Allocate H_JAC'))
    H_jac=0.0
  endif
  if(.not.allocated(pidx))then
    allocate(pidx(nobs),stat=ierr)
    call CheckStop(ierr,HERE('Allocate PIDX'))
    pidx=H_jac_idx(0,0,0,0,F,F)
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
      if(MasterProc.or.debug_obs) print dafmt,trim(damsg)
      cycle
    endif
    file=date2string(obsData(nd)%file,current_date)
    inquire(file=file,exist=exist)
    if(.not.exist)then
      write(damsg,"(A,I0,1X,A)")'WARNING: obsData#',nd,"file not found"
      if(MasterProc.or.debug_obs) print dafmt,trim(damsg)
      cycle
    endif
    open(IO_TMP,file=file,form='formatted',status='old',iostat=ierr)
    if(ierr/=0)then
      write(damsg,"(A,I0,1X,A)")'WARNING: obsData#',nd,"file not opened"
      if(MasterProc.or.debug_obs) print dafmt,trim(damsg)
      close(IO_TMP)
      cycle
    endif
    if(MasterProc.or.debug_obs)print dafmt,"Obsdata opened "//trim(file)
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
  real, intent(inout) :: lat,lon,alt
  real, intent(inout) :: yn

  integer :: i,j,l,ll,k,kk,igrp,g,ispec,nn,INFO
  real :: xn,Hj,unitconv
  integer, pointer, dimension(:) :: gspec=>null()      ! group array of indexes
  real,    pointer, dimension(:) :: gunit_conv=>null() ! group array of unit conv. factors
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  nn=n-obsData(ipar)%iobs(0)+1
  if((debug_obs.or.MasterProc).and.(n==obsData(ipar)%iobs(0))) &
    print dafmt,"Analysis of "//trim(obsData(ipar)%tag)
  if(debug_obs.or.(nn==obsData(ipar)%nobs.and.MasterProc)) print '(5(A,I0))',&
    'Obs(',me,'):',n,'/',nobs,',',nn,'/',obsData(ipar)%nobs
!-----------------------------------------------------------------------
! Find nearest grid cell
!-----------------------------------------------------------------------
  pidx(n)%in_mgrid=coord_in_domain("processor",lon,lat,iloc=i,jloc=j,&
    iglob=pidx(n)%i,jglob=pidx(n)%j) ! save global indexing for FFT domain
  pidx(n)%in_yslab=(pidx(n)%j>=iyf(0,me)).and.(pidx(n)%j<=iyf(1,me))
!-----------------------------------------------------------------------
! z grid coordinate of observations:
!-----------------------------------------------------------------------
  if(pidx(n)%in_mgrid)then
    ll=KMAX_MID
    do while(alt>z_bnd(i,j,ll).and.ll>=1)
      ll=ll-1
    enddo
  else
    ll=0
  endif ! share ll from in_mgrid-host CPU to all CPUs
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,ll,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,INFO)
  l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
  if(((l<1).or.debug_obs).and.pidx(n)%in_mgrid) then
    print dafmt,'Observation/Model Location:Weight[%],lon[deg],lat[deg],alt[km]'
    print "('#',I0,2(1X,A3,':',F0.2,',',F0.2,',',F0.2))",&
      n,'Observation',lon,lat,alt,&
        'Model',glon(i,j),glat(i,j),z_mid(i,j,ll)*1e-3
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
      pidx(n)%l0=l;pidx(n)%l1=l
      if(.not.pidx(n)%in_mgrid)return
      if(obsData(ipar)%unitroa)then   ! unit conversion        
        Hj=unitconv*roa(i,j,ll,1)
      else
        Hj=unitconv
      endif
      xn=xn_adv(ispec,i,j,ll)*FGSCALE
      yn=yn+Hj*xn
      H_jac(n,l,k)=Hj
!-----------------------------------------------------------------------
! direct observations: surface level
!-----------------------------------------------------------------------
    case('surface')
      pidx(n)%l0=l;pidx(n)%l1=l
      if(.not.pidx(n)%in_mgrid)return
      if(obsData(ipar)%unitroa)then  ! unit conversion,50m->3m
        Hj=unitconv*cfac(ispec,i,j)*roa(i,j,ll,1)
      else
        Hj=unitconv*cfac(ispec,i,j)
      endif
      xn=xn_adv(ispec,i,j,ll)*FGSCALE
      yn=yn+Hj*xn
      H_jac(n,l,k)=Hj
!-----------------------------------------------------------------------
! direct observations: COLUMN (eg NO2tc)
!-----------------------------------------------------------------------
    case('Trop.Col.')
      pidx(n)%l0=1;pidx(n)%l1=nlev
      if(.not.pidx(n)%in_mgrid)return
      call CheckStop(.not.obsData(ipar)%unitroa,&
        "Unsupported obsData unit for: "//trim(obsData(ipar)%tag))     
      do ll=1,KMAX_MID
        Hj=unitconv                       & ! unit conversion
          * roa(i,j,ll,1)                 & ! density.
          * (z_bnd(i,j,ll)-z_bnd(i,j,ll+1)) ! level thickness
        xn=xn_adv(ispec,i,j,ll)*FGSCALE
        yn=yn+Hj*xn
        l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
        if(l>0)H_jac(n,l,k)=Hj
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
      pidx(n)%l0=l;pidx(n)%l1=l
      if(.not.pidx(n)%in_mgrid)return
      if(obsData(ipar)%unitroa)&    ! unit conversion
        gunit_conv(:)=gunit_conv(:)*roa(i,j,ll,1)
      do g=1,size(gspec)
        Hj=gunit_conv(g)
        xn=xn_adv(gspec(g),i,j,ll)*FGSCALE
        yn=yn+Hj*xn
        kk=varSpecInv(gspec(g))
        k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
        if(k>0)H_jac(n,l,k)=Hj
      enddo
      ! PM water content in ug/m3 at model mid-levels
      if(index(obsData(ipar)%name,'PM')>0)&
        yn=yn+PM25_water(i,j,ll)*FGSCALE
!-----------------------------------------------------------------------
! indirect observations: surface level
!-----------------------------------------------------------------------
    case('surface')
      pidx(n)%l0=l;pidx(n)%l1=l
      if(.not.pidx(n)%in_mgrid)return
      if(obsData(ipar)%unitroa)&    ! unit conversion
        gunit_conv(:)=gunit_conv(:)*roa(i,j,ll,1)
      do g=1,size(gspec)
        Hj=gunit_conv(g)          & ! unit conversion
          *cfac(gspec(g),i,j)       ! 50m->3m conversion
        xn=xn_adv(gspec(g),i,j,ll)*FGSCALE
        yn=yn+Hj*xn
        kk=varSpecInv(gspec(g))
        k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
        if(k>0)H_jac(n,l,k)=Hj
      enddo
      ! PM water content in ug/m3 at surface level
      if(index(obsData(ipar)%name,'PM')>0)&
        yn=yn+PM25_water_rh50(i,j)*FGSCALE
!-----------------------------------------------------------------------
! direct observations: COLUMN (eg AOD)
!-----------------------------------------------------------------------
    case('Trop.Col.')
      pidx(n)%l0=l;pidx(n)%l1=nlev
      if(.not.pidx(n)%in_mgrid)return
      do ll=1,KMAX_MID
        Hj=roa(i,j,ll,1)                 & ! interpolation,density
          *(z_bnd(i,j,ll)-z_bnd(i,j,ll+1)) ! level thickness
        l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
        do g=1,size(gspec)
          xn=xn_adv(gspec(g),i,j,ll)*gunit_conv(g)*FGSCALE
          yn=yn+Hj*xn
          kk=varSpecInv(gspec(g))
          k=0;if(kk>0.and.observedVar(kk))k=ichemInv(kk)
          if(k>0.and.l>0)H_jac(n,l,k)=Hj*gunit_conv(g)
        enddo
      enddo
      ! PM water content in ug/m3 at model mid-levels
      if(index(obsData(ipar)%name,'PM')>0)then
        xn=0.0
        do ll=1,KMAX_MID
          xn=xn+PM25_water(i,j,ll)&
            *(z_bnd(i,j,ll)-z_bnd(i,j,ll+1)) ! level thickness
        enddo
        yn=yn+xn*FGSCALE
      endif
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
  integer :: i,j,n,frac,ipar,iyf0,iyf1
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
      i=pidx(n)%i;j=pidx(n)%j         ! observation operator for this observation
      if(pidx(n)%in_yslab) &          ! local y-slab
        h_loc(i,j,:,:)=H_jac(n,:,ichemobs(:))
      call chitox_adj(uh,h_loc)             ! to spectral space
      hbh=norm(uh,squared=.true.)+obsstddev(n)**2  ! HBH+R
      chisq(ipar)=chisq(ipar)+dep(n)**2/hbh ! chi-square for this observation & sum
      nm(ipar)=nm(ipar)+1                   ! number of used observations
      h_loc(pidx(n)%i,pidx(n)%j,:,:)=0e0    ! clean up for next iteration
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
endmodule DA_Obs_ml
