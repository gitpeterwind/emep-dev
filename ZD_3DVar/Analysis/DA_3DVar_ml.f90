#define STRING2(x) #x
#define STRING(x)  STRING2(x)
#define HERE(MSG)  MSG//" ("//__FILE__//":"//STRING(__LINE__)//")."
module DA_3DVar_ml
use ModelConstants_ml,only: KMAX_MID,KCHEMTOP,RUNDOMAIN,PPB,PPBINV,ANALYSIS,&
                            MasterProc,NPROC
use ChemSpecs_shl_ml, only: NSPEC_SHL        ! Maps indices
use ChemChemicals_ml, only: species          ! Gives names
use ChemGroups_ml,    only: chemgroups       ! group  names
use Chemfields_ml,    only: xn_adv
use GridValues_ml,    only: coord_in_domain
use TimeDate_ml,      only: date,current_date
use TimeDate_ExtraUtil_ml,only: date2string, compare_date
use Io_ml,            only: IO_TMP
use My_Timing_ml,     only: Code_timer,Add_2timing,NTIMING_UNIMOD
use Par_ml,           only: me,gi0,gi1,gj0,gj1,li0,li1,lj0,lj1,limax,ljmax
use CheckStop_ml,     only: CheckStop
use SmallUtils_ml,    only: find_index
use Functions_ml,     only: norm
use DA_ml,            only: debug=>DA_DEBUG,DAFMT_DEF=>DA_FMT_DEF,&
                            danml=>da_namelist,dafmt=>da_fmt_msg,damsg=>da_msg,&
                            tim_before=>datim_before,tim_after=>datim_after
use DA_Obs_ml
#if BIAS_CORRECTION
use DA_Bias_ml
#endif
use covmat_ml,        only: set_chemobs_idx,read_speccov
use spectralcov,      only: nx,ny,nlev,nchem,nxex,nyex,nex,&
                            nchemObs,ichemObs,nchemNoObs,ichemNoObs,&
                            nv1,mm,nn,nkstar,ikstar,vt,&
                            ucovmat,sqrt_lambda,sqrt_gamma,stddev,&
                            kx,kxmin,ky,kymin,lensav,wsave,&
                            FGSCALE,FGSCALE_INV,DAPREC
use chitox_ml,        only: matched_domain,iyf, &
                            chitox, chitox_adj, chi_Init, chi_Done, &
                            chi_vec2fc, chi_fc2vec, chi_fc2hc, chi_hc2fc
use mpi,              only: MPI_COMM_WORLD,MPI_BARRIER,MPI_ALLREDUCE,MPI_SUM,MPI_LAND,&
                            MPI_IN_PLACE,MPI_LOGICAL,MPI_INTEGER,MPI_DOUBLE_PRECISION
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
integer :: i,j,k,kLev0,kLev1,kk,ntot,maxobs,nv,nv2,nvar,nv2np,INFO,gIJ(4)
logical :: laux
real, allocatable, dimension(:) :: flat,flon,falt,obs,obsstddev1
real, dimension(:,:,:), pointer :: an_bgd=>null(),an_inc=>null()
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  ntot=nxex*nyex*nlev*nchemobs
  nv=min(nv1,ntot)
  nv2=nv*nxex*nyex
#ifdef BIAS_CORRECTION
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
  call Add_2timing(T_INNOV,tim_after,tim_before)
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
    allocate(chi(nv2),stat=ierr)
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
  call var3d(nv2,Jcost,chi,ntot,dzs)
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
subroutine var3d(nv2,Jcost,chi,ntot,dzs)
!-----------------------------------------------------------------------
! @description
! 3-D variational analysis
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
! Formal parameters
!-----------------------------------------------------------------------
  integer,      intent(in)    :: nv2,ntot
  real(kind=8), intent(inout) :: Jcost,chi(nv2),dzs
  real :: rzs,Jcost0
!-----------------------------------------------------------------------
! Work arrays
!-----------------------------------------------------------------------
  integer, parameter ::nupdates=4
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
  integer      :: iz(5),nrz,indic,omode,niter,nsim,reverse,ierr
  real(kind=8) :: epsg,df1
  integer      :: n,p,l,k
  external simul_rc,euclid,ctonbe,ctcabe
!-----------------------------------------------------------------------
! Make one first call to the objective function
!-----------------------------------------------------------------------
  allocate(gradJcost(nv2),stat=ierr)
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
  omode=-1                    ! trigger .not.OKmode, unless m1qn3 is called
  reverse=1                   ! reverse comunication mode
  indic=4                     ! calculate costFunction
  first_call=.true.           ! 1st call to costFunction
  epsg=dxmin**0.9             ! precision of the stopping criterion, based on ||gradient||
  niter=maxiter               ! maximal number of iterations accepted
  nsim =maxsim                ! maximal number of simulations accepted
  nrz=4*nv2+nupdates*(2*nv2+1)! rz dimension
  allocate(rz(nrz),stat=ierr) ! working array for m1qn3
  call CheckStop(ierr,HERE('Allocate RZ.'))

  write(damsg,"(A,'#',I3.3,': ',A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
    'Start  m1qn3',me,'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
  print dafmt,trim(damsg)

  do while(reverse>=0)
    call costFunction(indic,nv2,chi,Jcost,gradJcost,ntot,rzs,dzs)
    if(first_call)then
      first_call=.false.
      Jcost0=Jcost  ! initial cost function
      df1=Jcost*0.5 ! expected decrease estimate in f during the first iteration
      if(norm(gradJcost)<=1e-20)then
        print dafmt,'WARNING Starting point is almost optimal.&
          & No optimization needed..'
        reverse=-9  ! skipp m1qn3 call
        omode=1     ! fake OKmode
      endif
    endif
    if(reverse>=0)&
      call m1qn3(simul_rc,euclid,ctonbe,ctcabe,nv2,chi,Jcost,gradJcost,&
                dxmin,df1,epsg,normtype,impres,io,&
                imode,omode,niter,nsim,iz,rz,nrz,&
                reverse,indic,ntot,rzs,dzs)
    if(mod(niter,50)==0)FLUSH(io)
  enddo

  write(damsg,"(A,'#',I3.3,': ',A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
    'Finish m1qn3',me,'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
  print dafmt,trim(damsg)

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
  call Add_2timing(T_OPTIM,tim_after,tim_before)
!-----------------------------------------------------------------------
! Make a final call to costFunction for converting \chi to \delta x.
!-----------------------------------------------------------------------
  indic=-1
  call costFunction(indic,nv2,chi,Jcost,gradJcost,ntot,rzs,dzs)
  call my_deallocate(.false.,"NoMessage")
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  write(damsg,*)Jcost0,'-->',Jcost,'=',(1.0-Jcost/Jcost0)*100
  print dafmt,'Cost function '//trim(ADJUSTL(damsg))//'% Reduction'
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

  integer maxobs!,ipar(maxobs)
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
    if(debug_obs.and.any(pidx(n)%local))&
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
  if(debug)then
    info=sum((/(count(pidx%local(n)),n=1,4)/))
    print "(2(1X,A,1X,I0))","observations @",me,&
    '%local',nint(info*25.0/nobs) ! 0..4 pidx%local per obs
  endif
end subroutine get_innovations
#ifdef BIAS_CORRECTION
subroutine costFunction(ind,nv2np,chi,Jcost,gradJcost,ntot,rzs,dzs)
#else
subroutine costFunction(ind,nv2,chi,Jcost,gradJcost,ntot,rzs,dzs)
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
  integer, intent(in)    :: nv2np
  integer                :: nv2!!=nv2np-(NBIAS_PREDICTORS+1)*nobsData
  real, intent(in)       :: chi(nv2np)
  real, intent(inout)    :: Jcost,gradJcost(nv2np)
#else
  integer, intent(in)    :: nv2
  real, intent(in)       :: chi(nv2)
  real, intent(inout)    :: Jcost,gradJcost(nv2)
#endif
!-----------------------------------------------------------------------
! Dummy parameters
!-----------------------------------------------------------------------
  integer, intent(in)      :: ntot
  real, intent(in)         :: rzs
  real(kind=8), intent(in) :: dzs
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
 !complex chi_arr(nv1,nxex,nyex)
  complex(kind=8), pointer, dimension(:,:,:) :: chi_arr=>null(),chi_hc=>null()
  integer lenwork,i,ik,j,j0,j1,k,kk,l,l0,l1,m,n,m1,n1,mneg,nneg
 !integer maxobs,p,ij
 !real yn(nobs),dep(nobs),Oinvdep(nobs)
 !real dx(nxex,nyex,nlev,nchemobs)
 !real work(2*nxex*nyex),jobs,jb
  integer :: p, INFO, iyf0, iyf1
  real, pointer, dimension(:,:,:,:) :: dx_loc=>null()
  real, allocatable, dimension(:) :: yn,dep,Oinvdep,bn,Binvdep
  real :: jobs,jb,jbias
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  select case(ind)
  case (1)
    if(debug.and.MasterProc) print dafmt,"Free call to costFunction"
    return
  case(-1)
    if(debug.and.MasterProc) print dafmt,"Update observed species"
    call Add_2timing(T_OPTIM,tim_after,tim_before)
  case default
    if(debug.and.MasterProc) print dafmt,"Optimization costFunction"
  endselect
  call Add_2timing(T_OPTIM,tim_after,tim_before)
#ifdef BIAS_CORRECTION
  nv2=nv2np-(NBIAS_PREDICTORS+1)*nobsData
#endif
!-----------------------------------------------------------------------
! Copy chi (in compact storage format) into chi_arr (in matrix format),
! and select non reduntant part for FFTW3 (chi_hc)
!-----------------------------------------------------------------------
#ifdef BIAS_CORRECTION
  call chi_vec2fc(nv2np,chi,chi_arr)
#else
  call chi_vec2fc(nv2,chi,chi_arr)
#endif
!-----------------------------------------------------------------------
! conversion from spectral to model space
!-----------------------------------------------------------------------
  call chi_fc2hc(chi_arr,chi_hc)
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
    print "(I3,3(1X,A,'=',ES9.3))",me,&
      '||chi||',norm(chi),'||chi_arr||',norm(chi_arr),'||dx||',norm(dx)
    call Add_2timing(T_DOMEX,tim_after,tim_before)
  else
    print "(I3,3(1X,A,'=',ES9.3))",me,&
      '||chi||',norm(chi),'||chi_arr||',norm(chi_arr),'||dx_loc||',norm(dx_loc)
    call Add_2timing(T_CHI2X,tim_after,tim_before)
  endif
#ifdef BIAS_CORRECTION
  dbias=chi(nv2+1:)
#endif
!-----------------------------------------------------------------------
! update unobserved species:
!-----------------------------------------------------------------------
  if(ind==-1)then
    if(use_unobserved)then
    if(debug.and.MasterProc) print dafmt,"Update unobserved species"
      call update_unobserved(chi_arr)
! Brodcast du: y-slab --> xy-domain
      if(.not.matched_domain)then
        do n=0,NPROC-1
          if(iyf(0,n)==iyf(1,n))cycle
          CALL MPI_BCAST(du(:,iyf(0,n):iyf(1,n),:,:),&
                    size(du(:,iyf(0,n):iyf(1,n),:,:)),&
            MPI_DOUBLE_PRECISION,n,MPI_COMM_WORLD,INFO)
        enddo
        if(debug) print "(I3,1(1X,A,'=',E10.3))",me,'||du||',norm(du)
        call Add_2timing(T_DOMEX,tim_after,tim_before)
      endif
    call Add_2timing(T_UOSPC,tim_after,tim_before)
    endif
  endif ! return after chi^2 eval
!-----------------------------------------------------------------------
! chi^2 eval and return:
!-----------------------------------------------------------------------
  if(ind==-1)then
    if(use_chisq)then
      ! departures: dep=H(xb)+H_jac*dx-y=innov+yn
      ! dx|yn defined over full y-range
  do n=1,nobs
    do p=1,4
      i=pidx(n)%i(p)
      j=pidx(n)%j(p)
      l0=pidx(n)%l(0)
      l1=pidx(n)%l(1)
          innov(n)=innov(n)+sum(H_jac(n,p,l0:l1,:)*dx(i,j,l0:l1,:))
    enddo
  enddo
      call chisq_over_nobs2(nobs,innov)
    end if
    call Add_2timing(T_CHISQ,tim_after,tim_before)
    call my_deallocate(.false.,"NoMessage")
    return
  endif
!-----------------------------------------------------------------------
! Jb = 0.5 * \chi^{\dagger}*\chi
!-----------------------------------------------------------------------
  jb=0.5*norm(chi_arr,squared=.true.)
  call Add_2timing(T_OPTIM,tim_after,tim_before)
!-----------------------------------------------------------------------
! conversion from model to observation space
!-----------------------------------------------------------------------
  allocate(yn(nobs),dep(nobs),Oinvdep(nobs))
#ifdef BIAS_CORRECTION
  allocate(bn(nobsData),Binvdep(nobsData))
#endif
! dx|yn defined only over y-slab (iyf0:iyf1)
  do n=1,nobs
    yn(n)=0e0
    do p=1,4
      i=pidx(n)%i(p)
      j=pidx(n)%j(p)
      l0=pidx(n)%l(0)
      l1=pidx(n)%l(1)
      if(j>=iyf0.and.j<=iyf1)&
        yn(n)=yn(n)+sum(H_jac(n,p,l0:l1,:)*dx(i,j,l0:l1,:))
    enddo
  enddo
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,yn,nobs,MPI_DOUBLE_PRECISION,MPI_SUM,&
    MPI_COMM_WORLD,INFO)
  call Add_2timing(T_MPICF,tim_after,tim_before) ! MPI costFunction
!-----------------------------------------------------------------------
! Jobs = 0.5 * [ H(xb)+H_jac*dx-y ]^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  dep=innov+yn                  ! departures: h(xb)-y=H(xb)-y+H_jac*dx
#ifdef BIAS_CORRECTION
  bn=bias+dbias
! do n=1,nobsData
!   if(all(obsData(n)%iobs==0))cycle
!   i=obsData(n)%iobs(0);j=obsData(n)%iobs(1)
!   call CheckStop(any(ipar(i:j)/=n),"Inconsistent ipar/obsData")
!   print "('#',3(1X,I0,1X,ES12.3))",&
!     n,bn(n)*FGSCALE_INV,i,dep(i)*FGSCALE_INV,j,dep(j)*FGSCALE_INV
! enddo
  dep=dep+bn(ipar(:nobs))       ! departures: h(xb)+bias(x,beta)-y
#endif
  Oinvdep=dep/(obsstddev**2)    ! O^{-1} * [ H(xb)+H_jac*dx-y ]
  jobs=0.5*dot_product(dep,Oinvdep)
! if(debug)print "(I3,5(1X,A,'=',ES9.3))",me,&
!   '||H(xb)-y||',norm(innov),&
!   '||H_jac*dx||',norm(yn),&
!   '||H(xb)+H_jac*dx-y||',norm(dep),&
!   '||O^-1*[H(xb)+H_jac*dx-y]||',norm(Oinvdep),&
!   'Jo=1/2 * [H(xb)+H_jac*dx-y]^T * O^-1 * [H(xb)+H_jac*dx-y]',Jobs
!-----------------------------------------------------------------------
! Jbias = 0.5 * [ beta_b-beta ]^{T} * Bb^{-1} * [ beta_b-beta ]
!-----------------------------------------------------------------------
#ifdef BIAS_CORRECTION
  Binvdep=dbias*biasWeight/(biasStdDev**2)
! print "(('#',2(1X,I0),2(1X,ES12.3)))",&
!   (me,n,dbias(n)*FGSCALE_INV,Binvdep(n)*FGSCALE_INV,n=1,nobsData)
  jbias=0.5*dot_product(dbias,Binvdep)
#else
  jbias=0.0
#endif
!-----------------------------------------------------------------------
! J = Jb + Jobs + Jbias
!-----------------------------------------------------------------------
  Jcost = jb + jobs + jbias
!-----------------------------------------------------------------------
! grad J = grad Jb + grad Jobs + grad Jbias
!-----------------------------------------------------------------------
! grad Jb = independent elements of \chi
!-----------------------------------------------------------------------
  gradJcost(:nv2)=chi(:nv2)
!-----------------------------------------------------------------------
!  Bb^{-1} * [ beta_b-beta ] + bias_jac^{T} * O^{-1} * [ h(x)+bias(x,beta)-y ]
!-----------------------------------------------------------------------
#ifdef BIAS_CORRECTION 
  gradJcost(nv2+1:)=0.0
  do n=1,nobsData
    if(obsData(n)%nobs<=0)CYCLE
    i=obsData(n)%iobs(0);j=obsData(n)%iobs(1)
    gradJcost(nv2+n)=Binvdep(n)+sum(Oinvdep(i:j))
  enddo
#endif
!-----------------------------------------------------------------------
! grad Jobs = U^{-+}*H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
! dx := H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  dx=0e0
  do n=1,nobs
    do p=1,4
      i=pidx(n)%i(p)
      j=pidx(n)%j(p)
      l0=pidx(n)%l(0)
      l1=pidx(n)%l(1)
      dx(i,j,l0:l1,:)=dx(i,j,l0:l1,:)+H_jac(n,p,l0:l1,:)*Oinvdep(n)
    enddo
  enddo
  call Add_2timing(T_OPTIM,tim_after,tim_before)
!-----------------------------------------------------------------------
! apply fft etc on dx field:
!-----------------------------------------------------------------------
  call chitox_adj(chi_hc,dx_loc)   ! half-complex from (y-slab)
  call chi_hc2fc(chi_arr,chi_hc)   ! full complex form (full y)
! if(debug) print "(I3,3(1X,A,'=',ES9.3))",me,&
!  '||chi||',norm(chi),'||chi_arr||',norm(chi_arr),'||dx||',norm(dx)
!-----------------------------------------------------------------------
! add independent elements of grad Jobs to array gradJcost:
! grad Jobs = U^{-+}*H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  call chi_fc2vec(nv2,gradJcost,chi_arr)
  call Add_2timing(T_CHI2X,tim_after,tim_before)
!-----------------------------------------------------------------------
  if(debug)then
    write(damsg,"(5(1X,A,'=',ES9.3))"),&
      'J',Jcost,'Jb',Jb,'Jo',Jobs,'Jbias',Jbias,'nobs/2',nobs*0.5
#ifdef BIAS_CORRECTION
    Jb=norm(chi(:nv2))
    Jbias=norm(chi(nv2+1:))
#else
    Jb=norm(chi)
    Jbias=0.0
#endif
    Jobs=DIM(norm(gradJcost),Jb+Jbias)
    print "(I3,A,4(1X,A,'=',ES9.3))",me,trim(damsg),&
      '||DJ||',Jb+Jobs+Jbias,'||DJb||',Jb,'||DJo||',Jobs,'||DJbias||',Jbias
  endif
 !ind=4
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
  if(associated(chi_arr)) deallocate(chi_arr)
  if(associated(chi_hc))  deallocate(chi_hc)
  if(allocated(yn))       deallocate(yn)
  if(allocated(dep))      deallocate(dep)
  if(allocated(Oinvdep))  deallocate(Oinvdep)
  if(allocated(Binvdep))  deallocate(Binvdep)
endsubroutine my_deallocate
endsubroutine costFunction
subroutine update_unobserved(chi)
!-----------------------------------------------------------------------
! @description
! Update unobserved species once the variational analysis has converged
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none

  complex chi(nv1,nxex,nyex)

 !real du(nxex,nyex,nlev,nchem-nchemobs)
  integer i,j,k,kk,l,ik,m,n,lenwork,ierr,ibox,jbox,ndim
  real work(2*nxex*nyex)
  complex :: c(nxex,nyex)
  complex, allocatable :: temp1(:), temp2(:), temp3(:,:,:)

  if(nchem==nchemObs)then
    print dafmt,'WARNING no unobserved to update'
    return
  endif
  lenwork=2*nxex*nyex
  nchemNoObs=nchem-nchemObs
  ndim=nchemObs*nlev
  du=0e0
!-----------------------------------------------------------------------
! Compute c = X * \Lambda^{-1/2} * \chi, where \Lambda is the
! diagonal matrix containing the reduced set of eigenvalues of the
! spectral covariance matrix, and X is the matrix containing the
! corresponding eigenvectors:
!-----------------------------------------------------------------------
  allocate(temp1(nv1),temp2(ndim),temp3(ndim,nxex,nyex))
  do n=1,nyex
    do m=1,nxex
      ik=ikstar(m,n)
      if(ik<=nkstar)then
        temp1(:)=chi(:,m,n)/sqrt_lambda(:,ik)
        temp2(:)=matmul(vt(:,:,ik),temp1)
        temp3(:,m,n)=temp2(:)/sqrt_gamma(:,ik)
      endif
    enddo
  enddo
  do k=1,nchemnoobs
    kk=ichemnoobs(k)
    do l=1,nlev
      ibox=(k-1)*nlev+l
! nested double loop over 2D spectral space coordinates:
!-----------------------------------------------------------------------
! Multiply from the left with L, where L is a diagonal matrix
! with elements 1/sqrt(\gamma_i) along the diagonal, where \gamma_i
! denotes the bi-Fourier spectral density of the variance for
! level/component ibox:
!-----------------------------------------------------------------------
      c=cmplx(0e0,0e0)
      forall(n=1:nyex,m=1:nxex,ikstar(m,n)<=nkstar) &
        c(m,n)=sum(ucovmat(ibox,:,ikstar(m,n))*temp3(:,m,n))
!-----------------------------------------------------------------------
! Perform inverse 2d-FFT:
!-----------------------------------------------------------------------
      call cfft2b(nxex,nxex,nyex,c,wsave,lensav,work,lenwork,ierr)
! nested double loop over horizontal physical space coordinates:
!-----------------------------------------------------------------------
! Only update whithin the non-extended zone (since there are no
! observations in that zone => zero increment):
!-----------------------------------------------------------------------
! Multiply from the left with S, where S is a diagonal matrix
! containing the standard deviations sqrt(sigma).
!-----------------------------------------------------------------------
      forall(i=1:nx,j=1:ny) du(i,j,l,k)=real(c(i,j))*stddev(i,j,l,kk)
    enddo ! level l
  enddo   ! component k
  call my_deallocate(.false.,"NoMessage")
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
subroutine my_deallocate(verb,msg)
  logical, intent(in) :: verb
  character(len=*), intent(in) :: msg
  if(verb) print dafmt,msg
  if(allocated(temp1))  deallocate(temp1)
  if(allocated(temp2))  deallocate(temp2)
  if(allocated(temp3))  deallocate(temp3)
endsubroutine my_deallocate 
endsubroutine update_unobserved
endmodule DA_3DVar_ml
