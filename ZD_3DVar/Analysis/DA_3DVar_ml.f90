#define STRING2(x) #x
#define STRING(x)  STRING2(x)
#define HERE(MSG)  MSG//" ("//__FILE__//":"//STRING(__LINE__)//")."
module DA_3DVar_ml
use ModelConstants_ml,only: KMAX_MID,KCHEMTOP,RUNDOMAIN,PPB,PPBINV,ANALYSIS,&
                            MasterProc,NPROC
use ChemSpecs,        only: NSPEC_SHL, species  ! Maps indices, SPC names
use ChemGroups_ml,    only: chemgroups          ! group names
use Chemfields_ml,    only: xn_adv
use GridValues_ml,    only: coord_in_domain
use TimeDate_ml,      only: date,current_date
use TimeDate_ExtraUtil_ml,only: date2string, compare_date
use Io_ml,            only: IO_TMP
use My_Timing_ml,     only: Code_timer,Add_2timing,NTIMING_UNIMOD
use Par_ml,           only: me,gi0,gi1,gj0,gj1,li0,li1,lj0,lj1,limax,ljmax
use CheckStop_ml,     only: CheckStop
use SmallUtils_ml,    only: find_index
use Util_ml,          only: norm
use DA_ml,            only: debug=>DA_DEBUG,DAFMT_DEF=>DA_FMT_DEF,&
                            debug_3dv=>DA_DEBUG_3DV,debug_obs=>DA_DEBUG_OBS,&
                            danml=>da_namelist,dafmt=>da_fmt_msg,damsg=>da_msg,&
                            tim_before=>datim_before,tim_after=>datim_after
use DA_Obs_ml,        only: varName,obsVarName,observedVar,varSpec,varSpecInv,&
                            OBSERVATIONS,nobs,nobsData,obsData,&
                            innov,obsstddev,H_jac,pidx,H_op,chisq_over_nobs2,&
                            allocate_obs,read_obs,deallocate_obs
#ifdef BIAS_CORRECTION
use DA_Bias_ml
#endif
use covmat_ml,        only: set_chemobs_idx,read_speccov
use spectralcov,      only: nx,ny,nlev,nchem,nxex,nyex,nex,iChemInv,&
                            nchemObs,nv1,FGSCALE,FGSCALE_INV,DAPREC
use chitox_ml,        only: matched_domain, iyf, nhcrTot, nhcrLoc, &
                            chitox, chitou, chitox_adj, chi_Init, chi_Done, &
                            chi_vec2hc, chi_hc2vec, l2wR, l2wC
use MPI_Groups_ml,    only: MPI_COMM_WORLD,MPI_BARRIER,MPI_SUM,MPI_LAND,&
                            MPI_IN_PLACE,MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION
!                           MPI_ALLREDUCE
implicit none
private
public :: NTIMING_3DVAR,T_3DVAR,main_3dvar
integer, parameter :: ANALYSIS_NDATE = 4 ! Number of analysis to perform
real,    parameter :: ANALYSIS_RELINC_MAX=5.0 ! Limit dx,du to 500%
logical :: use_chisq=.true.,use_unobserved=.true.
integer :: maxiter=500,maxsim=500
type(date) :: analysis_date(ANALYSIS_NDATE)=(/& ! when to perform the analysys.
  date(-1,-1,-1,00,0),&              ! By default an analysis
  date(-1,-1,-1,06,0),&              ! it will be done every
  date(-1,-1,-1,12,0),&              ! 00,06,12,18 UTC
  date(-1,-1,-1,18,0)/)
integer, parameter :: debug_n=1,debug_p=1,debug_k=1
integer, allocatable, dimension(:) :: ipar
real, pointer, dimension(:,:,:,:) :: dx=>null(),du=>null()
integer, parameter ::       &
  T_DOMEX=NTIMING_UNIMOD+01,& ! 40: Domain extension
  T_INNOV=NTIMING_UNIMOD+02,& ! 41: Get innovations from observations
  T_MPIOB=NTIMING_UNIMOD+03,& ! 42: MPI observations
  T_OPTIM=NTIMING_UNIMOD+04,& ! 43: Optimization
  T_COSTF=NTIMING_UNIMOD+05,& ! 44: costFunction
  T_MPICF=NTIMING_UNIMOD+06,& ! 45: MPI costFunction
  T_CHI2X=NTIMING_UNIMOD+07,& ! 46: FFT trasformations
  T_MPIFT=NTIMING_UNIMOD+08,& ! 47: MPI FFT
  T_OBSPC=NTIMING_UNIMOD+09,& ! 48: Update observed species
  T_UOSPC=NTIMING_UNIMOD+10,& ! 49: Update unobserved species
  T_CHISQ=NTIMING_UNIMOD+11,& ! 50: CHI^2 evaluation
  T_3DVAR=NTIMING_UNIMOD+12,& ! 51: Total
  NTIMING_3DVAR=12
character(len=*), parameter, dimension(NTIMING_3DVAR) ::    &
  TIMING_3DVAR=['3DVar: Domain extension',                 &
                '3DVar: Get innovations from observations',&
                '3DVar:   MPI',                            &
                '3DVar: Optimization',                     &
                '3DVar: costFunction',                     &
                '3DVar:   MPI',                            &
                '3DVar: FFT trasformations',               &
                '3DVar:   MPI',                            &
                '3DVar: Update observed species',          &
                '3DVar: Update unobserved species',        &
                '3DVar: CHI^2 evaluation',                 &
                '3DVar: Total']
contains
!-----------------------------------------------------------------------
subroutine main_3dvar()
!-----------------------------------------------------------------------
! @description
! Main module for starting 3dvar program.
! @author M.Kahnert & AMVB
!-----------------------------------------------------------------------
implicit none
logical, save :: first_call=.true.
integer :: ierr
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  if(.not.ANALYSIS)return
  dafmt=date2string(DAFMT_DEF,current_date)
  if(first_call) then
    if(debug.and.MasterProc)print dafmt,'Initialisation'
    call init_3dvar()
    first_call=.false.
  endif

 !if(debug.and.MasterProc)print dafmt,'Test analysis'
  if(.not.compare_date(ANALYSIS_NDATE,current_date,analysis_date,wildcard=-1))return
  if(debug.and.MasterProc)print dafmt,'Start analysis'
!-----------------------------------------------------------------------
! Spectral decomposition
!-----------------------------------------------------------------------
  ! init new fft, will be applied on extended domain
  call chi_Init()
  call generic3dvar(ierr)
  call chi_Done()
endsubroutine main_3dvar
subroutine init_3dvar()
!-----------------------------------------------------------------------
! @description
! Initialise 3dvar variables.
! @author AMVB
!-----------------------------------------------------------------------
implicit none
integer(4) :: inNml=172
integer    :: nvar,k,ierr
namelist /DA_CONFIG/ analysis_date, nChem, nChemObs,&
                     varName, obsVarName, observedVar,&
                     use_chisq,use_unobserved,maxiter,maxsim
!-----------------------------------------------------------------------
! Set timming messages
!-----------------------------------------------------------------------
  call Code_timer(tim_before)
  do k=1,NTIMING_3DVAR
    call Add_2timing(NTIMING_UNIMOD+k,tim_after,tim_before,TIMING_3DVAR(k))
  enddo
!-----------------------------------------------------------------------
! Read config: variable names (observed & unobserved)
!-----------------------------------------------------------------------
  open(unit=inNml,file=danml,status='OLD',action='READ',&
       form='FORMATTED',iostat=ierr)
  call CheckStop(ierr,HERE('open namelist'))
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  k=find_index("DAOBS",chemgroups(:)%name)
  call CheckStop(k<1,HERE('DA group not found: "DAOBS"'))
  nChemObs=size(chemgroups(k)%ptr)
  obsVarName(:nChemObs)=species(chemgroups(k)%ptr)%name
  varName(:nChemObs)=obsVarName(:nChemObs)
  k=find_index("DAUNOBS",chemgroups(:)%name)
  if(k>0)then
  nChem=nChemObs+size(chemgroups(k)%ptr(:))
    varName(nChemObs+1:nChem)=species(chemgroups(k)%ptr)%name
  else
    use_unobserved=.false.
    nChem=nChemObs
    if(MasterProc)print dafmt,'WARNING: DA group not found: "DAUNOBS".'
  endif
  observedVar(:)=.false.
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  read(unit=inNml,nml=DA_CONFIG,iostat=ierr)
  call CheckStop(ierr,HERE('read namelist: DA_CONFIG'))
  call CheckStop(nChemObs>0.eqv.any(observedVar),'Incomplete/Redundant &
    &definition of nChemObs and observedVar on DA_CONFIG namelist'//HERE(''))
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  do nvar=1,nChemObs
    k=find_index(obsVarName(nvar),varName(:nChem))
    call CheckStop(k<1,HERE('Unknown observed variable: '//trim(obsVarName(nvar))))
    observedVar(k)=.true.
  enddo
  call CheckStop(.not.any(observedVar),&
    HERE('No observed variables found (observedVar)'))
  ! sort obsVarName following varName order
  nChemObs=0 ! count(observedVar)
  do nvar=1,nChem
    if(observedVar(nvar))then
      nChemObs=nChemObs+1
      obsVarName(nChemObs)=varName(nvar)
    endif
  enddo
  call CheckStop(nChemObs<1,HERE('No observed variables found (nChemObs)'))
  if(debug.and.MasterProc) then
    print dafmt,'B matrix description'
    print "(2(A,:,'(',I0,')'))",&
      'Variable: Observed',nChemObs,'/Unobserved',nChem-nChemObs
    do nvar=1,nChem
      if(observedVar(nvar))then
        print "(I4,': ',A10,'=O',I3.3,$)",nvar,trim(varName(nvar)),&
          count(observedVar(:nvar))
      else
        print "(I4,': ',A10,'=U',I3.3,$)",nvar,trim(varName(nvar)),&
          count(.not.observedVar(:nvar))
      endif
      if(mod(nvar,5)==0)print *,''
    enddo
    if(mod(nChem,5)/=0)print *,''
    write(*,nml=DA_CONFIG,delim='QUOTE')
  endif
!-----------------------------------------------------------------------
! From observed/unobserved variable names (varName) to model species
!-----------------------------------------------------------------------
  varSpec=-1
  varSpecInv=-1
  do nvar=1,nChem
    call CheckStop(count(varName(:nvar)==varName(nvar))>1,&
         HERE('Multiple definitions of variable name: '//trim(varName(nvar))))
    k=find_index(varName(nvar),species(NSPEC_SHL+1:)%name)
    call CheckStop(k<1,HERE('Wrong variable name: '//trim(varName(nvar))))
    varSpec(nvar)=k
    varSpecInv(k)=nvar
  enddo
!-----------------------------------------------------------------------
! Read observation parameters
!-----------------------------------------------------------------------
  read(unit=inNml,nml=OBSERVATIONS,iostat=ierr)
  call CheckStop(ierr,HERE('read namelist: OBSERVATIONS'))
!-----------------------------------------------------------------------
! Read bias parameters
!-----------------------------------------------------------------------
#if BIAS_CORRECTION
  read(unit=inNml,nml=BIAS_PREDICTOR,iostat=ierr)
  call CheckStop(ierr,HERE('read namelist: BIAS_PREDICTOR'))
  call CheckStop(nbiasData,nobsData,HERE('init_3dvar: nbiasData'))
#endif
!-----------------------------------------------------------------------
! Read covariance matrix
!-----------------------------------------------------------------------
  call set_chemobs_idx(nchem,observedVar(1:nchem))
  call read_speccov() ! & set extended domain dimensions (nxex,nyex)
  call CheckStop(nX,RUNDOMAIN(2)-RUNDOMAIN(1)+1,HERE("Inconsistent NX"))
  call CheckStop(nY,RUNDOMAIN(4)-RUNDOMAIN(3)+1,HERE("Inconsistent NY"))
  call CheckStop(nLev,KMAX_MID-KCHEMTOP+1      ,HERE("Inconsistent NLEV"))
endsubroutine init_3dvar
subroutine generic3dvar(ierr)
!-----------------------------------------------------------------------
! @description
! Generic version of 3d variational analysis.
! @author M.Kahnert
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
! Formal parameters
!-----------------------------------------------------------------------
integer ierr
!-----------------------------------------------------------------------
! Work arrrays
!-----------------------------------------------------------------------
real(kind=8) :: Jcost=0.0!huge(0.0)
real(kind=8), allocatable :: chi(:)
real(kind=8) :: dzs=0.0
!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
integer :: i,j,k,kLev0,kLev1,kk,ntot,maxobs,nv,nv2,nvTot,nvLoc,nvar,nv2np,INFO,gIJ(4)
logical :: laux
real, allocatable, dimension(:) :: flat,flon,falt,obs,obsstddev1
real, dimension(:,:,:), pointer :: an_bgd=>null(),an_inc=>null()
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  ntot=nxex*nyex*nlev*nchemobs
  nv=min(nv1,ntot)
  nvTot=nhcrTot*nv1  ! total number of half-complex-to-real values
  nvLoc=nhcrLoc*nv1  ! local number of half-complex-to-real values
#ifdef BIAS_CORRECTION
  nv2=nv*nxex*nyex
  nv2np=nv2+(NBIAS_PREDICTORS+1)*nobsData
#endif
  maxobs=nx*ny
!-----------------------------------------------------------------------
! read observations
!-----------------------------------------------------------------------
  if(debug.and.MasterProc)print dafmt,'Read observations'
  allocate(ipar(maxobs),flat(maxobs),flon(maxobs),falt(maxobs),&
           obs(maxobs),obsstddev1(maxobs),stat=ierr)
  call CheckStop(ierr,HERE('Allocate IPAR'))
  call read_obs("global",maxobs,flat,flon,falt,obs,obsstddev1,ipar)
  if(nobs==0)then
    call my_deallocate(MasterProc,'WARNING: No obserations found')
    return
  endif
  if(debug)print "(2(1X,A,1X,I0))","observations @",me,'read',nobs
!-----------------------------------------------------------------------
! scale bias data
!-----------------------------------------------------------------------
#ifdef BIAS_CORRECTION
  call allocate_bias()
  if(FGSCALE/=1e0)then
    bias(:)=bias(:)*FGSCALE
    biasStdDev(:)=biasStdDev(:)*FGSCALE
  endif
#endif
!-----------------------------------------------------------------------
! extend & scale background field
!-----------------------------------------------------------------------
  call Code_timer(tim_before)
!-----------------------------------------------------------------------
! compute innovations H(xb)-y, save in file
!-----------------------------------------------------------------------
  call get_innovations(maxobs,flat,flon,falt,obs,obsstddev1)
  if(all(innov==0.0))then
    call my_deallocate(MasterProc,'WARNING: No innovations found')
    return
  endif
  if(all(H_jac==0.0))print dafmt,'WARNING: H_jac==0.0'
!-----------------------------------------------------------------------
! allocate work arrrays
!-----------------------------------------------------------------------
  if(.not.associated(dx).and.nchemobs>0)then
    allocate(dx(nxex,nyex,nlev,nchemobs),stat=ierr)
    call CheckStop(ierr,HERE('Allocate DX'))
    dx=0e0
  endif
  if(.not.associated(du).and.nchem>nchemobs)then
    allocate(du(nxex,nyex,nlev,nchem-nchemobs),stat=ierr)
    call CheckStop(ierr,HERE('Allocate DU'))
    du=0e0
  endif
  if(.not.allocated(chi))then
#ifdef BIAS_CORRECTION
    allocate(chi(nv2np),stat=ierr)
#else
    if(nvLoc>0)then               ! defined local y-slab
      allocate(chi(nvLoc),stat=ierr)
    else
      allocate(chi(1),stat=ierr)  ! dummy
    endif
#endif
    call CheckStop(ierr,HERE('Allocate CHI'))
    chi=0e0
  endif
!-----------------------------------------------------------------------
! perform variational analysis
!-----------------------------------------------------------------------
  if(debug.and.MasterProc) print dafmt,'call var3d'
#ifdef BIAS_CORRECTION
  call var3d(nv2np,Jcost,chi,ntot,dzs)
#else
  call var3d(nvTot,nvLoc,Jcost,chi,ntot,dzs)
#endif
!-----------------------------------------------------------------------
! update (local) background field:
!   Only update whithin the non-extended zone (since there are no obs.
!   in that zone => zero increment). Leave BCs (lateral & top) unchanged.
!   Process only one SPC at the time.
!-----------------------------------------------------------------------
  kLev0=max(KCHEMTOP,KMAX_MID-nLev+1)
  kLev1=min(KMAX_MID,kLev0+nLev-1)
  call CheckStop(kLev1-kLev0+1,nLev,HERE('Size AN_BGD'))
! x/y indexes, outer bnd excluded:
!   global(gIJ(1):gIJ(2),gIJ(3):gIJ(4)) corresponds to local(li0:li1,lj0:lj1)
  gIJ=[max(gi0,2),min(gi1,nX-1),max(gj0,2),min(gj1,nY-1)]
  allocate(an_bgd(gIJ(1):gIJ(2),gIJ(3):gIJ(4),kLev0:kLev1),stat=ierr)
  call CheckStop(ierr,HERE('Allocation AN_BGD'))
!-----------------------------------------------------------------------
  do nvar=1,nChem
!-----------------------------------------------------------------------
! read result for \delta x and \delta
!-----------------------------------------------------------------------
    kk=ichemInv(nvar)                 ! 1..nchemObs||1..nchemNoObs
    if(observedVar(nvar))then
      if(.not.associated(dx))cycle
      an_inc=>dx(gIJ(1):gIJ(2),gIJ(3):gIJ(4),:,kk)
    else
      if(.not.associated(du))cycle
      an_inc=>du(gIJ(1):gIJ(2),gIJ(3):gIJ(4),:,kk)
    endif
!-----------------------------------------------------------------------
! scale background field
!-----------------------------------------------------------------------
    kk=varSpec(nvar)
    an_bgd=xn_adv(kk,li0:li1,lj0:lj1,kLev0:kLev1)*FGSCALE
!-----------------------------------------------------------------------
! add \delta x and \delta to background field:
!   ensure result between 0 and ANALYSIS_RELINC_MAX times background
!-----------------------------------------------------------------------
    an_bgd=min(max(an_bgd+an_inc,0.0),an_bgd*ANALYSIS_RELINC_MAX)
!-----------------------------------------------------------------------
! descale background field
!-----------------------------------------------------------------------
    xn_adv(kk,li0:li1,lj0:lj1,kLev0:kLev1)=an_bgd*FGSCALE_INV
  enddo
  call Add_2timing(T_DOMEX,tim_after,tim_before)
#ifdef BIAS_CORRECTION
  if(FGSCALE/=1e0)then
    bias(:)=(bias(:)+dbias(:))*FGSCALE_INV
    biasStdDev(:)=biasStdDev(:)*FGSCALE_INV
  endif
  call write_bias("bias_YYYYMMDD_hh.nml")
#endif
!-----------------------------------------------------------------------
! deallocate observation arrays
!-----------------------------------------------------------------------
  call my_deallocate(debug.and.MasterProc,'end generic3dvar')
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
subroutine my_deallocate(verb,msg)
  logical, intent(in) :: verb
  character(len=*), intent(in) :: msg
  integer :: ierr=0
  if(verb) print dafmt,msg
  if(allocated(ipar))       deallocate(ipar)
  if(allocated(flat))       deallocate(flat)
  if(allocated(flon))       deallocate(flon)
  if(allocated(falt))       deallocate(falt)
  if(allocated(obs))        deallocate(obs)
  if(allocated(obsstddev1)) deallocate(obsstddev1)
  call deallocate_obs()
  if(associated(an_bgd))    deallocate(an_bgd)
!!if(associated(an_inc))    nullify(an_inc)
  if(associated(dx))        deallocate(dx)
  if(associated(du))        deallocate(du)
  if(allocated(chi))        deallocate(chi)
#ifdef BIAS_CORRECTION
  call deallocate_bias()
#endif
endsubroutine  
endsubroutine generic3dvar
subroutine var3d(nvTot,nvLoc,Jcost,chi,ntot,dzs)
!-----------------------------------------------------------------------
! @description
! 3-D variational analysis
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
! Formal parameters
!-----------------------------------------------------------------------
  integer,      intent(in)    :: nvTot,nvLoc,ntot
  real(kind=8), intent(inout) :: Jcost,chi(nvLoc),dzs
  real :: rzs,Jcost0
!-----------------------------------------------------------------------
! Work arrays
!-----------------------------------------------------------------------
  integer, parameter :: nupdates=20
  real(kind=8), allocatable :: gradJcost(:),rz(:)
!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
  character(len=*),parameter :: LOGFILE="YYYYMMDDhh-PPP_m1qn3.log"
  character(len=*),parameter :: omode_str(-1:7)=(/&
  "('ERROR on m1qn3. omode==',I0,': ',A)",&
    "ERROR on m1qn3. omode==0: The simulator asks to stop by returning the value (indic=0)",&
    "Optimum found.  omode==1: Normal m1qn3 exit: successfull gradient test",&
    "ERROR on m1qn3. omode==2: One of the input arguments is not well initialized",&
    "ERROR on m1qn3. omode==3: Line-search blocked on tmax = 10**20",&
    "ERROR on m1qn3. omode==4: Reached maximal number of iterations (maxiter)",&
    "ERROR on m1qn3. omode==5: Reached maximal number of simulations (maxsim)",&
    "ERROR on m1qn3. omode==6: Stop on dxmin during the line-search",&
    "ERROR on m1qn3. omode==7: Either <g,d> is nonnegative or <y,s> is nonpositive"/)
  character(len=3) :: normtype='dfn'
  integer, parameter :: &
    OKmode(4)=[1,4,5,6],& ! omode accepted codes
    impres=3,&            ! verbosity level: 0=No print .. 5=Full verbosity
    io=IO_TMP,&           ! write unit for m1qn3 printouts
    imode(3)=[1,&         ! run in SIS mode
              0,&         ! cold start
              0]          ! calculate Jcost,gradJcost every interation
  real(kind=8), parameter :: dxmin=DAPREC ! norm resolution

  logical      :: first_call
  integer      :: iz(5),nrzTot,nrzLoc,indic,omode,niter,nsim,reverse,ierr
  real(kind=8) :: epsg,df1
  integer      :: n,p,l,k
  external simul_rc,euclid_weight_mpi,ctonbe,ctcabe
!-----------------------------------------------------------------------
! Make one first call to the objective function
!-----------------------------------------------------------------------
  allocate(gradJcost(nvLoc),stat=ierr)
  call CheckStop(ierr,HERE('Allocate gradJcost'))
!-----------------------------------------------------------------------
! Call minimization routine (m1qn3-3.3;2009)
! subroutine m1qn3 (simul, prosca, ctonb, ctcab, n, x, f, g,
!                   dxmin, df1, epsg, normtype, impres, io,
!                   imode, omode, niter, nsim, iz, dz, ndz,
!                   reverse, indic, izs, rzs, dzs)
!-----------------------------------------------------------------------
  open(io,file=date2string(LOGFILE,current_date),iostat=ierr)
  call CheckStop(ierr,HERE('open m1qn3 logfile'))

  call Code_timer(tim_before)
  omode=-1          ! trigger .not.OKmode, unless m1qn3 is called
  reverse=1         ! reverse comunication mode
  indic=4           ! calculate costFunction
  first_call=.true. ! 1st call to costFunction
  epsg=dxmin**0.9   ! precision of the stopping criterion, based on ||gradient||
  niter=maxiter     ! maximal number of iterations accepted
  nsim =maxsim      ! maximal number of simulations accepted
  nrzLoc=4*nvLoc+nupdates*(2*nvLoc+1) ! rz dimension
  nrzTot=4*nvTot+nupdates*(2*nvTot+1) ! for m1qn3 to estimate nupdates 
  allocate(rz(nrzLoc),stat=ierr)      ! working array for m1qn3
  call CheckStop(ierr,HERE('Allocate RZ.'))

  write(damsg,"(A,'#',I3.3,': ',A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
    'Start  m1qn3',me,'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
  if(MasterProc.or.debug_3dv)print dafmt,trim(damsg)

  do while(reverse>=0)
    call costFunction(indic,nvLoc,chi,Jcost,gradJcost)
    if(first_call)then
      first_call=.false.
      Jcost0=Jcost  ! initial cost function
      df1=Jcost*0.5 ! expected decrease estimate in f during the first iteration
      if(norm(gradJcost,l2wR,reduce=.true.)<=1e-20)then
        if(MasterProc)&
          print dafmt,'WARNING Starting point is almost optimal.&
            & No optimization needed..'
        reverse=-9  ! skipp m1qn3 call
        omode=1     ! fake OKmode
      endif
    endif
    if(reverse>=0)&
      call m1qn3(simul_rc,euclid_weight_mpi,ctonbe,ctcabe,&
                nvTot,nvLoc,chi,Jcost,gradJcost,&
                dxmin,df1,epsg,normtype,impres,io,&
                imode,omode,niter,nsim,iz,rz,nrzTot,nrzLoc,&
                reverse,indic,l2wR,rzs,dzs)
    if(mod(niter,50)==0)FLUSH(io)
  enddo

  write(damsg,"(A,'#',I3.3,': ',A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
    'Finish m1qn3',me,'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
  if(MasterProc.or.debug_3dv)print dafmt,trim(damsg)

  select case(omode)
  case(1:7)
    damsg=trim(omode_str(omode))
  case default
    write(damsg,omode_str(-1))omode,"Unknown omode"
  endselect
  write(io,*)trim(damsg)
  if(debug)then
    write(damsg,"('m1qn3#',I3.3,'. ',A)")me,trim(damsg)
    print dafmt,trim(damsg)
  endif
  call CheckStop(all(omode/=OKmode),HERE(trim(damsg)))

  close(io)
!-----------------------------------------------------------------------
! Make a final call to costFunction for converting \chi to \delta x.
!-----------------------------------------------------------------------
  indic=-1
  call costFunction(indic,nvLoc,chi,Jcost,gradJcost)
  call my_deallocate(.false.,"NoMessage")
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  write(damsg,*)Jcost0,'-->',Jcost,'=',(1.0-Jcost/Jcost0)*100
  if(MasterProc.or.debug_3dv)&
    print dafmt,'Cost function '//trim(ADJUSTL(damsg))//'% Reduction'
  call Add_2timing(T_OPTIM,tim_after,tim_before)
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
subroutine my_deallocate(verb,msg)
  logical, intent(in) :: verb
  character(len=*), intent(in) :: msg
  if(verb) print dafmt,msg
  if(allocated(gradJcost)) deallocate(gradJcost)
  if(allocated(rz))        deallocate(rz)
endsubroutine my_deallocate
endsubroutine var3d
subroutine get_innovations(maxobs,flat,flon,falt,obs,obs_stddev)
!-----------------------------------------------------------------------
! @description
! Compute innovations H(xb)-y, where xb is
! the background field, y denotes the observations, and
! H is an operator mapping from model to observation space.
! On input, innov contains the observations y. On output,
! innov contains the innovations.
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none
  integer :: maxobs
  real, dimension(maxobs) :: flat,flon,falt,obs,obs_stddev
  integer :: n,INFO
  real :: alt(nx,ny,nlev)
  real, allocatable :: yn(:)
!-----------------------------------------------------------------------
! Allocate & initialise dynamic arrays:
!-----------------------------------------------------------------------
  call allocate_obs()
  allocate(yn(nobs))
!-----------------------------------------------------------------------
! x,y grid coordinates of observations:
!-----------------------------------------------------------------------
  do n=1,nobs
    call CheckStop(.not.coord_in_domain("global",flon(n),flat(n)),&
      HERE("Observation outside geographical domain"))
    call CheckStop(obs(n),[obsData(ipar(n))%min,obsData(ipar(n))%max],&
      HERE("Observation outside accepted value range"))
!-----------------------------------------------------------------------
! mapping from model to obs-space:
!-----------------------------------------------------------------------
! Units Conversion:
!   Model [mol/mol] to Obs [obsData(ipar)%unit] conversion is handeled 
!   by H_op (DA_Obs_ml) via Units_Scale (Units_ml). Conversion info
!   is stored on obsData(ipar)%unitconv and obsData(ipar)%unitroa.
!-----------------------------------------------------------------------
    call H_op(flat(n),flon(n),falt(n),n,yn(n),ipar(n))
    if(debug_obs.and.pidx(n)%in_mgrid)&
      print "(I3,'#',I0,2(1X,A3,':',ES12.3))",me,&
        n,'Observation',obs(n),'Model',yn(n)*FGSCALE_INV
  enddo
  call Add_2timing(T_INNOV,tim_after,tim_before)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,yn,nobs,MPI_DOUBLE_PRECISION,MPI_SUM,&
    MPI_COMM_WORLD,INFO)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,H_jac,size(H_jac),MPI_DOUBLE_PRECISION,MPI_SUM,&
    MPI_COMM_WORLD,INFO)
  call Add_2timing(T_MPIOB,tim_after,tim_before)
!-----------------------------------------------------------------------
! Innovation = yn-y
!-----------------------------------------------------------------------
  innov=yn-obs*FGSCALE
  obsstddev=obs_stddev*FGSCALE
  deallocate(yn)
  if(debug)& ! local obs
    print "(2(1X,A,1X,I0))","observations @",me,'%in_mgrid',count(pidx%in_mgrid)
  call Add_2timing(T_INNOV,tim_after,tim_before)
end subroutine get_innovations
#ifdef BIAS_CORRECTION
subroutine costFunction(ind,nvLoc,npLoc,chi,Jcost,gradJcost)
#else
subroutine costFunction(ind,nvLoc,chi,Jcost,gradJcost)
#endif
!-----------------------------------------------------------------------
! @description
! Compute the costfunction and its gradient
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
! Formal parameters
!-----------------------------------------------------------------------
  integer, intent(inout) :: ind
#ifdef BIAS_CORRECTION
  integer, intent(in)    :: nvLoc,npLoc!!=(NBIAS_PREDICTORS+1)*nobsData in NPROC-1
  real, intent(in)       :: chi(nvLoc+npLoc)
  real, intent(inout)    :: Jcost,gradJcost(nvLoc+npLoc)
#else
  integer, intent(in)    :: nvLoc
  real, intent(in)       :: chi(nvLoc)
  real, intent(inout)    :: Jcost,gradJcost(nvLoc)
#endif
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
  complex(kind=8),pointer, dimension(:,:,:) :: chi_hc=>null()
  integer lenwork,i,ik,j,j0,j1,k,kk,l,l0,l1,m,n,m1,n1,mneg,nneg
  integer :: p, INFO, iyf0, iyf1
  real, pointer, dimension(:,:,:,:) :: dx_loc=>null(),du_loc=>null()
  real, allocatable, dimension(:) :: yn,dep,Oinvdep,bn,Binvdep
  real :: jobs,jb,jbias
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  call Add_2timing(T_OPTIM,tim_after,tim_before)
  select case(ind)
  case (1)
    if(debug.and.MasterProc) print dafmt,"Free call to costFunction"
    return
  case(-1)
    if(debug.and.MasterProc) print dafmt,"Update observed species"
  case default
    if(debug.and.MasterProc) print dafmt,"Optimization costFunction"
  endselect
!-----------------------------------------------------------------------
! conversion from spectral to model space
!-----------------------------------------------------------------------
#ifdef BIAS_CORRECTION
  call CheckStop((me<(NPROC-1)).and.(npLoc>0),&
    HERE("Only NPROC-1 allocates bias indexes on CHI"))
  call chi_vec2hc(nvLoc,chi(:nvLoc),chi_hc)
  ! only NPROC-1 allocates bias indexes on CHI
  if(me==(NPROC-1))dbias=chi(nvLoc+1:)
  CALL MPI_BCAST(dbias,1,MPI_DOUBLE_PRECISION,NPROC-1,MPI_COMM_WORLD,INFO)
#else
  call chi_vec2hc(nvLoc,chi,chi_hc)
#endif
  iyf0=iyf(0,me);iyf1=iyf(1,me) ! local y-slab
  dx_loc=>dx(:,iyf0:iyf1,:,:)
  call chitox(chi_hc,dx_loc)
  if(ind==-1)then
! Brodcast dx: y-slab --> xy-domain
    if(.not.matched_domain)then
      do n=0,NPROC-1
        if(iyf(0,n)==iyf(1,n))cycle
        CALL MPI_BCAST(dx(:,iyf(0,n):iyf(1,n),:,:),&
                  size(dx(:,iyf(0,n):iyf(1,n),:,:)),&
          MPI_DOUBLE_PRECISION,n,MPI_COMM_WORLD,INFO)
      enddo
    endif
    if(debug)&
      write(damsg,"(I3,3(1X,A,'=',ES10.3))")me,&
        '||chi||',norm(chi,l2wR,reduce=.true.),&
        '||chi_arr||',norm(chi_hc,l2wC,reduce=.true.),&
        '||dx||',norm(dx)
    if(MasterProc.and.debug) print *,trim(damsg)
    call Add_2timing(T_DOMEX,tim_after,tim_before)
 else
    if(debug_3dv)&
      write(damsg,"(I3,3(1X,A,'=',ES10.3))")me,&
        '||chi||',norm(chi,l2wR,reduce=.true.),&
        '||chi_arr||',norm(chi_hc,l2wC,reduce=.true.),&
        '||dx||',norm(dx_loc,reduce=.true.)
    if(MasterProc.and.debug_3dv) print *,trim(damsg)
    call Add_2timing(T_CHI2X,tim_after,tim_before)
  endif
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  if(ind==-1)then
! update unobserved species:
    if(use_unobserved)then
      if(debug.and.MasterProc) print dafmt,"Update unobserved species"
      du_loc=>du(:,iyf0:iyf1,:,:)
      call chitou(chi_hc,du_loc)
! Brodcast du: y-slab --> xy-domain
      if(.not.matched_domain)then
        do n=0,NPROC-1
          if(iyf(0,n)==iyf(1,n))cycle
          CALL MPI_BCAST(du(:,iyf(0,n):iyf(1,n),:,:),&
                    size(du(:,iyf(0,n):iyf(1,n),:,:)),&
            MPI_DOUBLE_PRECISION,n,MPI_COMM_WORLD,INFO)
        enddo
        if(debug) &
          print "(I3,1(1X,A,'=',E10.3))",me,'||du||',norm(du)
        call Add_2timing(T_DOMEX,tim_after,tim_before)
      endif
    call Add_2timing(T_UOSPC,tim_after,tim_before)
    endif
! chi^2 eval and return:
    if(use_chisq)then
      ! departures: dep=H(xb)+H_jac*dx-y=innov+yn
      ! dx|yn defined over full y-range
      allocate(dep(nobs))
      do n=1,nobs
        i=pidx(n)%i;j=pidx(n)%j;l0=pidx(n)%l0;l1=pidx(n)%l1
#ifdef BIAS_CORRECTION
        dep(n)=innov(n)+sum(H_jac(n,l0:l1,:)*dx(i,j,l0:l1,:))+bn(ipar(n))
#else
        dep(n)=innov(n)+sum(H_jac(n,l0:l1,:)*dx(i,j,l0:l1,:))
#endif
      enddo
      call chisq_over_nobs2(nobs,dep)
    endif
    call Add_2timing(T_CHISQ,tim_after,tim_before)
    call my_deallocate(.false.,"NoMessage")
    return
  endif
!-----------------------------------------------------------------------
! Jbias = 0.5 * [ beta_b-beta ]^{T} * Bb^{-1} * [ beta_b-beta ]
!-----------------------------------------------------------------------
  jbias=0.0
#ifdef BIAS_CORRECTION
  allocate(bn(nobsData),Binvdep(nobsData))
  bn=bias+dbias
! do n=1,nobsData
!   if(obsData(n)%nobs<=0)cycle
!   i=obsData(n)%iobs(0);j=obsData(n)%iobs(1)
!   call CheckStop(any(ipar(i:j)/=n),"Inconsistent ipar/obsData")
!   print "('#',3(1X,I0,1X,ES12.3))",&
!     n,bn(n)*FGSCALE_INV,i,dep(i)*FGSCALE_INV,j,dep(j)*FGSCALE_INV
! enddo
  Binvdep=0.0   ! only NPROC-1 has jbias>0
  if(me==(NPROC-1))Binvdep=dbias*biasWeight/(biasStdDev**2)
  jbias=0.5*dot_product(dbias,Binvdep)
! print "(('#',2(1X,I0),2(1X,ES12.3)))",&
!   (me,n,dbias(n)*FGSCALE_INV,Binvdep(n)*FGSCALE_INV,n=1,nobsData)
#endif
!-----------------------------------------------------------------------
! Jb = 0.5 * \chi^{\dagger}*\chi
!-----------------------------------------------------------------------
  jb=0.0
  if(nvLoc>0)&    ! defined local y-slab
    jb=0.5*norm(chi_hc,l2wC,squared=.true.,reduce=.false.)
!-----------------------------------------------------------------------
! conversion from model to observation space
!-----------------------------------------------------------------------
  allocate(yn(nobs),dep(nobs),Oinvdep(nobs))
! dx|yn|dep defined only over y-slab (iyf0:iyf1)
  yn(:)=0e0
  dep(:)=0e0
  do n=1,nobs
    if(.not.pidx(n)%in_yslab)cycle
    i=pidx(n)%i;j=pidx(n)%j;l0=pidx(n)%l0;l1=pidx(n)%l1
    yn(n)=sum(H_jac(n,l0:l1,:)*dx(i,j,l0:l1,:))
! departures: h(xb)+bias(x,beta)-y=H(xb)-y+H_jac*dx+bias(x,beta)
#ifdef BIAS_CORRECTION
    dep(n)=innov(n)+yn(n)+bn(ipar(n))
#else
    dep(n)=innov(n)+yn(n)
#endif
  enddo
!-----------------------------------------------------------------------
! Jobs = 0.5 * [ H(xb)+H_jac*dx-y ]^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  Oinvdep=dep/(obsstddev**2)  ! O^{-1} * [ H(xb)+H_jac*dx-y ]
  jobs=0.5*dot_product(dep,Oinvdep)
! if(debug)print "(I3,5(1X,A,'=',ES10.3))",me,&
!   '||H(xb)-y||',norm(innov),&
!   '||H_jac*dx||',norm(yn),&
!   '||H(xb)+H_jac*dx-y||',norm(dep),&
!   '||O^-1*[H(xb)+H_jac*dx-y]||',norm(Oinvdep),&
!   'Jo=1/2 * [H(xb)+H_jac*dx-y]^T * O^-1 * [H(xb)+H_jac*dx-y]',Jobs
!-----------------------------------------------------------------------
! J = Jb + Jobs + Jbias
!-----------------------------------------------------------------------
  Jcost = jb + jobs + jbias
  call Add_2timing(T_COSTF,tim_after,tim_before)
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,Jcost,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
    MPI_COMM_WORLD,INFO)
  call Add_2timing(T_MPICF,tim_after,tim_before) ! MPI costFunction
!-----------------------------------------------------------------------
! grad J = grad Jb + grad Jobs + grad Jbias
!-----------------------------------------------------------------------
! grad Jb = independent elements of \chi
!-----------------------------------------------------------------------
  gradJcost(:nvLoc)=chi(:nvLoc)
!-----------------------------------------------------------------------
!  Bb^{-1} * [ beta_b-beta ] + bias_jac^{T} * O^{-1} * [ h(x)+bias(x,beta)-y ]
!-----------------------------------------------------------------------
#ifdef BIAS_CORRECTION 
  gradJbias(:)=0.0
  do n=1,nobsData
    if(obsData(n)%nobs<=0)cycle
    i=obsData(n)%iobs(0);j=obsData(n)%iobs(1)
    gradJbias(n)=Binvdep(n)+sum(Oinvdep(i:j))
  enddo
  ! only NPROC-1 allocates bias indexes on CHI
  CALL MPI_REDUCE(gradJbias,gradJcost(nvLoc+1:),nobsData,&
       MPI_DOUBLE_PRECISION,MPI_SUM,NPROC-1,MPI_COMM_WORLD,INFO)
#endif
!-----------------------------------------------------------------------
! grad Jobs = U^{-+}*H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
! dx := H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  dx=0e0
  do n=1,nobs
    if(.not.pidx(n)%in_yslab)cycle
    i=pidx(n)%i;j=pidx(n)%j;l0=pidx(n)%l0;l1=pidx(n)%l1
    dx(i,j,l0:l1,:)=dx(i,j,l0:l1,:)+H_jac(n,l0:l1,:)*Oinvdep(n)
  enddo
  call Add_2timing(T_COSTF,tim_after,tim_before)
!-----------------------------------------------------------------------
! apply fft etc on dx field:
!-----------------------------------------------------------------------
  call chitox_adj(chi_hc,dx_loc)   ! half-complex from (y-slab)
!-----------------------------------------------------------------------
! add independent elements of grad Jobs to array gradJcost:
! grad Jobs = U^{-+}*H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
! and store L2 weights for half-complex-real packaging (l2w)
!-----------------------------------------------------------------------
  call chi_hc2vec(nvLoc,gradJcost,chi_hc)
  call Add_2timing(T_CHI2X,tim_after,tim_before)
!-----------------------------------------------------------------------
  if(debug)then
    write(damsg,"(5(1X,A,'=',ES10.3))"),&
      'J',Jcost,'Jb',Jb,'Jo',Jobs,'Jbias',Jbias,'nobs/2',nobs*0.5
#ifdef BIAS_CORRECTION
    Jb=norm(chi(:nvLoc),l2wR,reduce=.true.)
    ! only NPROC-1 allocates bias indexes on CHI
    if(me==(NPROC-1))Jbias=norm(chi(nvLoc+1:),reduce=.true.)
    CALL MPI_BCAST(Jbias,1,MPI_DOUBLE_PRECISION,NPROC-1,MPI_COMM_WORLD,INFO)
#else
    Jb=norm(chi,l2wR,reduce=.true.)
    Jbias=0.0
#endif
    Jobs=DIM(norm(gradJcost,l2wR,reduce=.true.),Jb+Jbias)
    print "(I3,A,4(1X,A,'=',ES10.3))",me,trim(damsg),&
      '||DJ||',Jb+Jobs+Jbias,'||DJb||',Jb,'||DJo||',Jobs,'||DJbias||',Jbias
  endif
  call Add_2timing(T_COSTF,tim_after,tim_before)
  call my_deallocate(.false.,"NoMessage")
! CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)
! call CheckStop(MasterProc,HERE("AMVB-TEST"))
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
subroutine my_deallocate(verb,msg)
  logical, intent(in) :: verb
  character(len=*), intent(in) :: msg
  integer :: ierr=0
  if(verb) print dafmt,msg
!!if(associated(chi_arr)) deallocate(chi_arr)
  if(associated(chi_hc))  deallocate(chi_hc)
  if(allocated(yn))       deallocate(yn)
  if(allocated(dep))      deallocate(dep)
  if(allocated(Oinvdep))  deallocate(Oinvdep)
  if(allocated(Binvdep))  deallocate(Binvdep)
endsubroutine my_deallocate
endsubroutine costFunction
endmodule DA_3DVar_ml
