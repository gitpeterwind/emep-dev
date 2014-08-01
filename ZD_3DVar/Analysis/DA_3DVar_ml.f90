module DA_3DVar_ml
use TimeDate_ml,      only: date,current_date
use My_Timing_ml,     only: Code_timer,Add_2timing
use ModelConstants_ml,only: KMAX_MID,KCHEMTOP,RUNDOMAIN,PPB,PPBINV,ANALYSIS,MasterProc
use ChemSpecs_shl_ml, only: NSPEC_SHL        ! Maps indices
!se ChemSpecs_tot_ml, only: NSPEC_TOT
use ChemChemicals_ml, only: species          ! Gives names
use ChemGroups_ml,    only: chemgroups       ! group  names
use Chemfields_ml,    only: xn_adv
use GridValues_ml,    only: coord_in_processor
use CheckStop_ml,     only: CheckStop
use SmallUtils_ml,    only: find_index
use TimeDate_ExtraUtil_ml,only: date2string
use DA_ml,            only: debug=>DA_DEBUG,DAFMT_DEF=>DA_FMT_DEF,&
                            dafmt=>da_fmt_msg,damsg=>da_msg,&
                            tim_before=>datim_before,tim_after=>datim_after
use DA_Obs_ml
use DA_Bias_ml
use Util_ml,          only: io_check,compare_date,norm
use exd_domain_ml,    only: EXT_DOMAIN,EXT_DOMAIN_INV
use covmat_ml,        only: set_chemobs_idx,read_speccov
use spectralcov,      only: nx,ny,nlev,nchem,nxex,nyex,nex,&
                            nchemObs,ichemObs,nchemNoObs,ichemNoObs,&
                            nv1,mm,nn,nkstar,ikstar,vt,&
                            ucovmat,sqrt_lambda,sqrt_gamma,stddev,&
                            kx,kxmin,ky,kymin,lensav,wsave,&
                            FGSCALE,FGSCALE_INV,DAPREC
implicit none
integer,private,parameter :: ANALYSIS_NDATE = 4 ! Number of analysis to perform
real, parameter, private  :: ANALYSIS_RELINC_MAX=5.0 ! Limit dx,du to 500%
logical,   private :: calculate_chisq = .true.
type(date),private :: analysis_date(ANALYSIS_NDATE)=(/& ! when to perform the analysys.
                 date(-1,-1,-1,00,0),&              ! By default an analysis
                 date(-1,-1,-1,06,0),&              ! it will be done every
                 date(-1,-1,-1,12,0),&              ! 00,06,12,18 UTC
                 date(-1,-1,-1,18,0)/)
integer, parameter, private :: maxiter=100*5,maxsim=100*5
integer, parameter, private :: debug_n=1,debug_p=1,debug_k=1
integer, allocatable, dimension(:), private :: ipar
real, allocatable, dimension(:,:,:,:), private :: dx,du
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
  if(.not.compare_date(ANALYSIS_NDATE,current_date,analysis_date,wildcard=-1))return
  if(debug.and.MasterProc)print dafmt,'Start analysis'
  call generic3dvar(ierr)
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
                     calculate_chisq
!-----------------------------------------------------------------------
! Read config: variable names (observed & unobserved)
!-----------------------------------------------------------------------
  open(unit=inNml,file='namelist.nml',status='OLD',action='READ',&
       form='FORMATTED',iostat=ierr)
  call io_check(ierr,'open namelist')
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
! nChemObs=0;observedVar(:)=.false.
! call define_chemicals()
! call Init_ChemGroups()
! dafmt=DAFMT_DEF
  k=find_index("DAOBS",chemgroups(:)%name)
  call CheckStop(k<1,'DA group not found: "DAOBS".')
  nChemObs=size(chemgroups(k)%ptr)
  obsVarName(:nChemObs)=species(chemgroups(k)%ptr)%name
  varName(:nChemObs)=obsVarName(:nChemObs)
  k=find_index("DAUNOBS",chemgroups(:)%name)
  if(k>0)then
  nChem=nChemObs+size(chemgroups(k)%ptr(:))
    varName(nChemObs+1:nChem)=species(chemgroups(k)%ptr)%name
  else
    nChem=nChemObs
    if(MasterProc)print dafmt,'WARNING: DA group not found: "DAUNOBS".'
  endif
  observedVar(:)=.false.
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  read(unit=inNml,nml=DA_CONFIG,iostat=ierr)
  call io_check(ierr,'read namelist: DA_CONFIG')
  call CheckStop(nChemObs>0.eqv.any(observedVar),&
    'Incomplete/Redundant definition of nChemObs and observedVar on DA_CONFIG namelist.')
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  do nvar=1,nChemObs
    k=find_index(obsVarName(nvar),varName(:nChem))
    call CheckStop(k<1,'Unknown observed variable: '//trim(obsVarName(nvar)))
    observedVar(k)=.true.
  enddo
  call CheckStop(.not.any(observedVar),'No observed variables found (observedVar).')
  ! sort obsVarName following varName order
  nChemObs=0 ! count(observedVar)
  do nvar=1,nChem
    if(observedVar(nvar))then
      nChemObs=nChemObs+1
      obsVarName(nChemObs)=varName(nvar)
    endif
  enddo
  call CheckStop(nChemObs<1,'No observed variables found (nChemObs).')
  if(debug.and.MasterProc) then
    print dafmt,'B matrix description'
    print "(2(A,:,'(',I0,')'))",&
      'Variable: Observed',nChemObs,'/Unobserved',nChem-nChemObs
    do nvar=1,nChem
      if(observedVar(nvar))then
        print "(I4,': ',A10,'=O',I3.3,$)",nvar,trim(varName(nvar)),count(observedVar(:nvar))
      else
        print "(I4,': ',A10,'=U',I3.3,$)",nvar,trim(varName(nvar)),count(.not.observedVar(:nvar))
      endif
      if(mod(nvar,5)==0)print *,''
    enddo
    if(mod(nChem,5)/=0)print *,''
  endif
  if(debug) write(*,nml=DA_CONFIG,delim='APOSTROPHE')
!-----------------------------------------------------------------------
! From observed/unobserved variable names (varName) to model species
!-----------------------------------------------------------------------
  varSpec=-1
  varSpecInv=-1
  do nvar=1,nChem
    call CheckStop(count(varName(:nvar)==varName(nvar))>1,&
         'Multiple definitions of variable name: '//trim(varName(nvar)))
    k=find_index(varName(nvar),species(NSPEC_SHL+1:)%name)
    call CheckStop(k<1,'Wrong variable name: '//trim(varName(nvar)))
    varSpec(nvar)=k
    varSpecInv(k)=nvar
  enddo
!-----------------------------------------------------------------------
! Read observation parameters
!-----------------------------------------------------------------------
  read(unit=inNml,nml=OBSERVATIONS,iostat=ierr)
  call io_check(ierr,'read namelist: OBSERVATIONS')
!-----------------------------------------------------------------------
! Read bias parameters
!-----------------------------------------------------------------------
#if BIAS_CORRECTION
  read(unit=inNml,nml=BIAS_PREDICTOR,iostat=ierr)
  call io_check(ierr,'read namelist: BIAS_PREDICTOR')
  call CheckStop(nbiasData,nobsData,'init_3dvar: nbiasData')
#endif
!-----------------------------------------------------------------------
! Read the rest of the variables
!-----------------------------------------------------------------------
  call set_chemobs_idx(nchem,observedVar(1:nchem))
  call read_speccov
  call CheckStop(nX,RUNDOMAIN(2)-RUNDOMAIN(1)+1,"init_3dvar: Inconsistent NX")
  call CheckStop(nY,RUNDOMAIN(4)-RUNDOMAIN(3)+1,"init_3dvar: Inconsistent NY")
  call CheckStop(nLev,KMAX_MID-KCHEMTOP+1      ,"init_3dvar: Inconsistent NLEV")
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
integer i,j,ilev,ilev0,k,kk,ntot,ii,maxobs,nv,nv2,nvar,nv2np
real, allocatable, dimension(:) :: flat,flon,falt,obs,obsstddev1
!external costFunction
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  ntot=nxex*nyex*nlev*nchemobs
  nv=min(nv1,ntot)
  nv2=nv*nxex*nyex
  nv2np=nv2+(NBIAS_PREDICTORS+1)*nobsData
  maxobs=nx*ny
!-----------------------------------------------------------------------
! read observations
!-----------------------------------------------------------------------
  if(debug.and.MasterProc)print dafmt,'Read observations'
  allocate(ipar(maxobs),flat(maxobs),flon(maxobs),falt(maxobs),&
           obs(maxobs),obsstddev1(maxobs),stat=ierr)
  call CheckStop(ierr,'Allocation error: IPAR.')
  call read_obs(maxobs,flat,flon,falt,obs,obsstddev1,ipar,ierr)
  if(nobs==0)then
    call my_deallocate(MasterProc,'WARNING: No obserations found')
    return
  endif
  if(debug.and.MasterProc)print "(1X,I0,1X,A)",nobs,'observations read.'
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
  if(.not.allocated(xn_adv_ex))then
    allocate(xn_adv_ex(NXEX,NYEX,NLEV,NCHEM),stat=ierr)
    call CheckStop(ierr,'Allocation error: XN_ADV_EX.')
    xn_adv_ex=0.0
  endif
  if(debug.and.MasterProc)write(damsg,dafmt)'Extending..'
  do nvar=1,nChem
    if(debug.and.MasterProc)then
      if(mod(nvar,10)==1)print "(/A,$)",trim(damsg)
      print "(1X,I0,':',A,$)",nvar,trim(varName(nvar))
    endif
    CALL EXT_DOMAIN(NX,NY,NLEV,NXEX,NYEX,&
           xn_adv(varSpec(nvar),:,:,KMAX_MID-nlev+1:KMAX_MID),&
           xn_adv_ex(:,:,:,nvar))
  enddo
  if(debug.and.MasterProc)print *,''
  if(FGSCALE/=1e0) xn_adv_ex(:,:,:,:)=xn_adv_ex(:,:,:,:)*FGSCALE
  call Add_2timing(40,tim_after,tim_before,'3DVar: Domain extension.')
!-----------------------------------------------------------------------
! compute innovations H(xb)-y, save in file
!-----------------------------------------------------------------------
  call get_innovations(maxobs,flat,flon,falt,obs,obsstddev1)
  call Add_2timing(41,tim_after,tim_before,'3DVar: Get innovations from observations.')
  if(all(innov==0.0))then
    call my_deallocate(MasterProc,'WARNING: No innovations found')
    return
  endif
  if(all(H_jac==0.0))&
    print dafmt,'WARNING: H_jac==0.0'
!-----------------------------------------------------------------------
! allocate work arrrays
!-----------------------------------------------------------------------
  if(.not.allocated(dx).and.nchemobs>0)then
    allocate(dx(nxex,nyex,nlev,nchemobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: dx.')
    dx=0e0
  endif
  if(.not.allocated(du).and.nchem>nchemobs)then
!   allocate(du(nxex,nyex,nlev,nchem-nchemobs),stat=ierr)
    allocate(du(nx,ny,nlev,nchem-nchemobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: du.')
    du=0e0
  endif
  if(.not.allocated(chi))then
#ifdef BIAS_CORRECTION
    allocate(chi(nv2np),stat=ierr)
#else
    allocate(chi(nv2),stat=ierr)
#endif
    call CheckStop(ierr,'Allocation error: CHI.')
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
! read result for \delta x and \delta u and add to background field;
!-----------------------------------------------------------------------
! Only update whithin the non-extended zone (since there are no obs. in
! that zone => zero increment). Leave BCs (lateral & top) unchanged.
!-----------------------------------------------------------------------
! forall(i=1:nxex,j=1:nyex,ilev=1:nlev)
! forall(i=2:nx-1,j=2:ny-1,ilev=2:nlev)
  ilev0=max(KCHEMTOP,KMAX_MID-nlev+1)
  if(allocated(dx))forall(i=2:nx-1,j=2:ny-1,ilev=ilev0:nlev)&
    xn_adv_ex(i,j,ilev,ichemObs  )=&
      min(max(xn_adv_ex(i,j,ilev,ichemObs  )+dx(i,j,ilev,:),0.0),&
              xn_adv_ex(i,j,ilev,ichemObs  )*ANALYSIS_RELINC_MAX)
  if(allocated(du))forall(i=2:nx-1,j=2:ny-1,ilev=ilev0:nlev)&
    xn_adv_ex(i,j,ilev,ichemNoObs)=&
      min(max(xn_adv_ex(i,j,ilev,ichemNoObs)+du(i,j,ilev,:),0.0),&
              xn_adv_ex(i,j,ilev,ichemNoObs)*ANALYSIS_RELINC_MAX)
! if(debug.and.MasterProc) then
!   k=ichemobs(1)
!   kk=varSpec(k)
!   do i=1,nx
!     do j=1,ny
!       if(xn_adv_ex(i,j,nlev,k)==0.0) &
!        print *,i,j,dx(i,j,nlev,1),xn_adv(kk,i,j,KMAX_MID)
!     enddo
!   enddo
! endif
!-----------------------------------------------------------------------
! descale & crop background field:
!-----------------------------------------------------------------------
  if(FGSCALE_INV/=1e0) xn_adv_ex(:,:,:,:)=xn_adv_ex(:,:,:,:)*FGSCALE_INV
  do nvar=1,nChem
    CALL EXT_DOMAIN_INV(NX,NY,NLEV,NXEX,NYEX,&
           xn_adv_ex(:,:,:,nvar),&
           xn_adv(varSpec(nvar),:,:,KMAX_MID-nlev+1:KMAX_MID))
  enddo
!-----------------------------------------------------------------------
! update & descale bias data
!-----------------------------------------------------------------------
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
  if(verb) print dafmt,msg
  if(allocated(ipar))       deallocate(ipar)
  if(allocated(flat))       deallocate(flat)
  if(allocated(flon))       deallocate(flon)
  if(allocated(falt))       deallocate(falt)
  if(allocated(obs))        deallocate(obs)
  if(allocated(obsstddev1)) deallocate(obsstddev1)
  if(allocated(xn_adv_ex))  deallocate(xn_adv_ex)
  call deallocate_obs()
  call deallocate_bias()
  if(allocated(dx))         deallocate(dx)
  if(allocated(du))         deallocate(du)
  if(allocated(chi))        deallocate(chi)
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
! logical lqc
! external simul
!-----------------------------------------------------------------------
! Work arrays
!-----------------------------------------------------------------------
  integer, parameter ::nupdates=4
 !real(kind=8) :: gradJcost(nv2),rz(4*nv2+nupdates*(2*nv2+1))
  real(kind=8), allocatable :: gradJcost(:),rz(:)
!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
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
    impres=5,&        ! verbosity level: 0=No print .. 5=Full verbosity
    io=6,&            ! write unit for m1qn3 printouts
    imode(3)=(/ 1,&   ! run in SIS mode
                0,&   ! cold start
                0/)   ! calculate Jcost,gradJcost every interation
  real(kind=8), parameter :: dxmin=DAPREC ! norm resolution
  integer      :: iz(5),nrz,indic,omode,niter,nsim,reverse,ierr
  real(kind=8) :: epsg,df1
  integer      :: n,p,l,k
  external simul_rc,euclid,ctonbe,ctcabe
!-----------------------------------------------------------------------
! Make one first call to the objective function
!-----------------------------------------------------------------------
  allocate(gradJcost(nv2),stat=ierr)
  call CheckStop(ierr,'Allocation error: gradJcost.')
  indic=4
  if(debug.and.MasterProc) print dafmt,'call costFunction'
  call costFunction(indic,nv2,chi,Jcost,gradJcost,ntot,rzs,dzs)
  Jcost0=Jcost  ! Initial cost function
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  if(norm(gradJcost)<=1e-20)then
    print dafmt,'WARNING Starting point is almost optimal.&
      & No optimization needed..'
!   if(debug.and.MasterProc) then
!     do n=1,nobs
!       do l=pidx(n)%l(0),pidx(n)%l(1)
!         do p=1,4
!           print "(2(A,I0,2(',',I0)),1X,A,30(1X,A,'=',ES12.3,:))",&
!                 '#',n,p,l,"-->",pidx(n)%i(p),pidx(n)%j(p),KMAX_MID-nlev+l,&
!                 'H_jac',(trim(obsVarName(k)),H_jac(n,p,l,k),k=1,nChemObs)
!         enddo
!       enddo
!     enddo
!     print "(2(1X,A,'=',ES12.3))",'||DJ||',norm(gradJcost),'||DH||',norm(H_jac)
!     call CheckStop('Debug Jcost')
!   endif
  else
!-----------------------------------------------------------------------
! Call minimization routine (m1qn3-3.3;2009)
! subroutine m1qn3 (simul, prosca, ctonb, ctcab, n, x, f, g,
!                   dxmin, df1, epsg, normtype, impres, io,
!                   imode, omode, niter, nsim, iz, dz, ndz,
!                   reverse, indic, izs, rzs, dzs)
! Very old call minimization routine (m1qn3-2.0d;1996)
! subroutine m1qn3 (simul, prosca, ctonb, ctcab, n, x, f, g,
!                   dxmin, df1, epsg,           impres, io,
!                   mode,         niter, nsim, iz, rz, nrz)
!-----------------------------------------------------------------------
    epsg=dxmin**0.9
    df1=Jcost*0.5
    reverse=0!;if(debug)reverse=1
    niter=maxiter
    nsim =maxsim
   !nrz=size(rz)
    nrz=4*nv2+nupdates*(2*nv2+1)
    allocate(rz(nrz),stat=ierr)
    call CheckStop(ierr,'Allocation error: RZ.')
    if(debug.and.MasterProc.and.impres>0)then
      write(damsg,"(A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
        'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
      print dafmt,'Calling m1qn3 '//trim(damsg)
    endif
    call Code_timer(tim_before)
    if(reverse==0)then
      call m1qn3(costFunction,euclid,ctonbe,ctcabe,nv2,chi,Jcost,gradJcost,&
                dxmin,df1,epsg,normtype,impres,io,&
                imode,omode,niter,nsim,iz,rz,nrz,&
                reverse,indic,ntot,rzs,dzs)
    else
      m1qn3_call: do
        call m1qn3(simul_rc,euclid,ctonbe,ctcabe,nv2,chi,Jcost,gradJcost,&
                  dxmin,df1,epsg,normtype,impres,io,&
                  imode,omode,niter,nsim,iz,rz,nrz,&
                  reverse,indic,ntot,rzs,dzs)
        if(reverse<0) exit m1qn3_call
        call costFunction(indic,nv2,chi,Jcost,gradJcost,ntot,rzs,dzs)
      enddo m1qn3_call
    endif

    if(debug.and.MasterProc.and.impres>0)then
      write(damsg,"(A,'=',3(I0,:,','),3(1X,A,'=',I0,:,','))"),&
        'mode',imode,'reverse',reverse,'niter',niter,'nsim',nsim
      print dafmt,'Finish  m1qn3 '//trim(damsg)
    endif
    if(omode<0.or.omode>7)then
      write(damsg,omode_str(-1))omode,"Unknown omode"
    else
      damsg=trim(omode_str(omode))
    endif
    if(debug.and.MasterProc)print *,trim(damsg)
    call CheckStop(.not.any(omode==(/1,4,5,6/)),trim(damsg))
    call Add_2timing(42,tim_after,tim_before,'3DVar: Optimization.')
  endif
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

  integer i,j,l,k,n,ilev,err,ierr
  real yn,alt(nx,ny,nlev),rho0
! logical inside
!-----------------------------------------------------------------------
! Allocate & initialise dynamic arrays:
!-----------------------------------------------------------------------
  call allocate_obs()
! print *,size(flat),size(flon),size(falt),size(ipar)
!-----------------------------------------------------------------------
! x,y grid coordinates of observations:
!-----------------------------------------------------------------------
  DO_OBS: do n=1,nobs
!   if((.not.coord_in_processor(flon(n),flat(n))).or.&
!      (obs(n)<obsData(ipar(n))%min).or.&
!      (obs(n)>obsData(ipar(n))%max))then
!     innov(n)=0e0
!     cycle DO_OBS
!   endif
    call CheckStop(.not.coord_in_processor(flon(n),flat(n)),&
      "Observation outside geographical domain")
    call CheckStop((obs(n)<obsData(ipar(n))%min).or.&
                   (obs(n)>obsData(ipar(n))%max),&
      "Observation outside accepted value range")
!-----------------------------------------------------------------------
! mapping from model to obs-space:
!-----------------------------------------------------------------------
    call H_op(flat(n),flon(n),falt(n),n,yn,rho0,ipar(n))
!   if(obsData(ipar(n))%error_rel>0)&
!     obs_stddev(n)=max(obs_stddev(n),obs(n)*obsData(ipar(n))%error_rel)
!   if(obsData(ipar(n))%error_rep>0)&
!     obs_stddev(n)=max(obs_stddev(n),obsData(ipar(n))%error_rep)
!-----------------------------------------------------------------------
!**************** check  Unimod units **********************************
! Innovation = yn-y
!  for chemical observations: on input, obs(n)=y=observation [g/m3];
!    on output innov(n)=innovation [g/kg].
!    Convert standard deviation from g/m3 to g/kg
!  for NO2 Trop.Col. observation are converted to [g/kg] uppon reading
!    with the use of sfalling factors
!  for lidar observations: on input, obs(n)=y [1/m/sr]
!    likewise on output
!-----------------------------------------------------------------------
!   if(trim(obspararr(ipar(n)))=='BSCA')then
!     innov(n)=yn-obs(n)
!     obsstddev(n)=obs_stddev(n)
!   else
!     innov(n)=yn-obs(n)/rho0
!     obsstddev(n)=obs_stddev(n)/rho0
!   endif
    innov(n)=yn-obs(n)*FGSCALE
    obsstddev(n)=obs_stddev(n)*FGSCALE
    if(debug.and.MasterProc) print "('#',I0,2(1X,A3,':',ES12.3))",&
      n,'Observation',obs(n),'Model',yn*FGSCALE_INV
  enddo DO_OBS
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
  complex(kind=8), allocatable, dimension(:,:,:) :: chi_arr
  integer lenwork,i,ik,j,k,kk,l,l0,l1,m,n,m1,n1,mneg,nneg
 !integer maxobs,p,ij
 !real yn(nobs),dep(nobs),Oinvdep(nobs)
 !real dx(nxex,nyex,nlev,nchemobs)
 !real work(2*nxex*nyex),jobs,jb
  integer :: p
  real, allocatable, dimension(:) :: yn,dep,Oinvdep,bn,Binvdep 
  real :: jobs,jb,jbias
  integer :: kx1,ky1,kxm1,kym1,kx2
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  select case(ind)
  case (1)
    if(debug.and.MasterProc) print dafmt,"Free call to costFunction"
    return
  case(-1)
    if(debug.and.MasterProc) print dafmt,"Update observed species"
  case default
    if(debug.and.MasterProc) print dafmt,"Optimization costFunction"
  endselect
  call Add_2timing(42,tim_after,tim_before,'3DVar: Optimization.')
#ifdef BIAS_CORRECTION
  nv2=nv2np-(NBIAS_PREDICTORS+1)*nobsData
#endif
! maxobs=nx*ny
! no4=4*nobs
! no8=8*nobs
  kx1=kx+1
  ky1=ky+1
  kxm1=kxmin+1
  kym1=kymin+1
  kx2=nxex-kxmin+1
!   if(debug.and.MasterProc) print "(4(1X,A,'=',I0))",&
!     'kx',kx,'ky',ky,'kxmin',kxmin,'kymin',kymin
!-----------------------------------------------------------------------
! Copy chi (in compact storage format) into chi_arr (in matrix format),
! taking into account the relation chi_arr(i,m,n)=conjg(chi_arr(i,-m,-n))
!-----------------------------------------------------------------------
  allocate(chi_arr(nv1,nxex,nyex))
  chi_arr=cmplx(0e0,0e0)
!-----------------------------------------------------------------------
! independent elements:
!-----------------------------------------------------------------------
  k=0
  do i=1,nv1
    do m=1,kx1
      k=k+1
      chi_arr(i,m,1)=cmplx(chi(k),0e0)
    enddo
    do m=1,kx1
      do n=2,ky1
        k=k+1
        chi_arr(i,m,n)=cmplx(chi(k),0e0)
      enddo
    enddo
    do m=kx2,nxex
      do n=2,kym1
        k=k+1
        chi_arr(i,m,n)=cmplx(chi(k),0e0)
      enddo
    enddo
    do m=2,kxm1
      k=k+1
      chi_arr(i,m,1)=chi_arr(i,m,1)+cmplx(0e0,chi(k))
    enddo
    do n=2,kym1
      k=k+1
      chi_arr(i,1,n)=chi_arr(i,1,n)+cmplx(0e0,chi(k))
    enddo
    do m=2,kx1
      do n=2,ky1
        if(mod(nxex,2)==0.and.mod(nyex,2)==0.and.m==kx1.and.n==ky1)then
          continue
        else
          k=k+1
          chi_arr(i,m,n)=chi_arr(i,m,n)+cmplx(0e0,chi(k))
        endif
      enddo
    enddo
    do m=kx2,nxex
      do n=2,kym1
        k=k+1
        chi_arr(i,m,n)=chi_arr(i,m,n)+cmplx(0e0,chi(k))
      enddo
    enddo
  enddo
!-----------------------------------------------------------------------
! dependent elements:
!-----------------------------------------------------------------------
  do i=1,nv1
    do m1=-kxmin,kx
      m=mm(m1)
      mneg=mm(-m1)
      do n1=-kymin,-1
        n=nn(n1)
        nneg=nn(-n1)
        chi_arr(i,m,n)=conjg(chi_arr(i,mneg,nneg))
      enddo
    enddo
    do m1=-kxmin,-1
      m=mm(m1)
      mneg=mm(-m1)
      n=nn(0)
      chi_arr(i,m,n)=conjg(chi_arr(i,mneg,n))
    enddo
    if(mod(nyex,2)==0)then
      do m1=-kxmin,-1
        m=mm(m1)
        mneg=mm(-m1)
        n=nn(ky)
        chi_arr(i,m,n)=conjg(chi_arr(i,mneg,n))
      enddo
    endif
!     zero elements:
!AMVB 2009-12-11: Bug found my MK 2009-12-09
    m=mm(kx);n=nn(ky)
!   m=mm(kx);n=mm(ky)
    chi_arr(i,1,1)=cmplx(real(chi_arr(i,1,1)),0e0)
    if(mod(nxex,2)==0)&
      chi_arr(i,m,1)=cmplx(real(chi_arr(i,m,1)),0e0)
    if(mod(nyex,2)==0)&
      chi_arr(i,1,n)=cmplx(real(chi_arr(i,1,n)),0e0)
    if(mod(nxex,2)==0.and.mod(nyex,2)==0)&
      chi_arr(i,m,n)=cmplx(real(chi_arr(i,m,n)),0e0)
   enddo
!-----------------------------------------------------------------------
! conversion from spectral to model space
!-----------------------------------------------------------------------
  call chitox(chi_arr,dx)
#ifdef BIAS_CORRECTION
  dbias=chi(nv2+1:)
#endif
!-----------------------------------------------------------------------
! update unobserved species:
!-----------------------------------------------------------------------
  if(ind==-1)then
    call Add_2timing(44,tim_after,tim_before,'3DVar: Update observed species.')
#ifndef NO_UPDATE_UNOBSERVED
    if(debug.and.MasterProc) print dafmt,"Update unobserved species"
    call update_unobserved(chi_arr)
    call Add_2timing(45,tim_after,tim_before,'3DVar: Update unobserved species.')
#endif
    if(.not.calculate_chisq)return
  endif ! return after chi^2 eval
!-----------------------------------------------------------------------
! conversion from model to observation space
!-----------------------------------------------------------------------
  allocate(yn(nobs),dep(nobs),Oinvdep(nobs))
#ifdef BIAS_CORRECTION
  allocate(bn(nobsData),Binvdep(nobsData))
#endif

! do n=1,nobs
!   yn(n)=0e0
!   do p=1,no4
!     i=pidx(p)
!     j=pidx(p+no4)
!     l=pidx(p+no8)
!     do k=1,nchemobs
!       kk=ichemobs(k)
!       yn(n)=yn(n)+H_jac(n,p,kk)*dx(i,j,l,k)
!     enddo
!   enddo
! enddo
  do n=1,nobs
    yn(n)=0e0
    do p=1,4
      i=pidx(n)%i(p)
      j=pidx(n)%j(p)
#ifdef OLD_MATCH_CODE
      do l=pidx(n)%l(0),pidx(n)%l(1)
        yn(n)=yn(n)+dot_product(H_jac(n,p,l,:),dx(i,j,l,:))
      enddo
#else
     l0=pidx(n)%l(0)
     l1=pidx(n)%l(1)
     yn(n)=yn(n)+sum(H_jac(n,p,l0:l1,:)*dx(i,j,l0:l1,:))
#endif
    enddo
  enddo
  if(debug.and.MasterProc)then
    n=min(debug_n,nobs);p=min(debug_p,4)
    k=min(debug_k,nchemobs);if(obsData(ipar(n))%found) k=obsData(ipar(n))%ichem
    kk=ichemObs(k)
    write(damsg,"(A3,'#',I0,2(1X,A,'=',I0,:,':',A))")&
      'Observation',n,'varName',kk,trim(varName(kk)),'obsVarName',k,trim(obsVarName(k))
    print dafmt,trim(damsg)
    print "(1X,A,3(1X,A1,'=',I0,3(:,',',I0) ))",&
              'xxx0:','i',pidx(n)%i,'j',pidx(n)%j,'l',pidx(n)%l
    i=pidx(n)%i(p);j=pidx(n)%j(p);l0=pidx(n)%l(0);l1=pidx(n)%l(1)
    if(p==0)then
      print *,'xxx1:',yn(n),norm(H_jac(n,:,l0:l1,k)),&
              norm(dx(pidx(n)%i(:),pidx(n)%j(:),l0:l1,k)) ! <-- NOT what inteded print
    else
      print *,'xxx1:',yn(n),H_jac(n,p,l1,k),dx(i,j,l1,k)
    end if
  endif
!-----------------------------------------------------------------------
! chi^2 eval and return:
!----------------------------------------------------------------------- 
  if(ind==-1)then
    if(calculate_chisq) &
      call chisq_over_nobs2(nobs,innov+yn) ! departures: dep=H(xb)+H_jac*dx-y
    call Add_2timing(46,tim_after,tim_before,'3DVar: CHI^2 evaluation.')
    call my_deallocate(.false.,"NoMessage")
    return
  endif
#ifdef OLD_MATCH_CODE
!-----------------------------------------------------------------------
! departures H(xb)+H_jac*dx-y
!-----------------------------------------------------------------------
  do n=1,nobs
    dep(n)=innov(n)+yn(n)
  enddo
!-----------------------------------------------------------------------
! O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  do n=1,nobs
    Oinvdep(n)=dep(n)/(obsstddev(n)**2)
  enddo
!-----------------------------------------------------------------------
! Jobs = 0.5 * [ H(xb)+H_jac*dx-y ]^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  jobs=0e0
  do n=1,nobs
    jobs=jobs+dep(n)*Oinvdep(n)
  enddo
  jobs=0.5e0*jobs
!-----------------------------------------------------------------------
! Jb = 0.5 * \chi^{\dagger}*\chi
!-----------------------------------------------------------------------
  jb=0e0
  do i=1,nv1
    do m=1,nxex
      do n=1,nyex
        jb=jb+real(chi_arr(i,m,n))**2+aimag(chi_arr(i,m,n))**2
      enddo
    enddo
  enddo
  jb=0.5e0*jb
#else
!-----------------------------------------------------------------------
! Jb = 0.5 * \chi^{\dagger}*\chi
!-----------------------------------------------------------------------
 !jb=0.5*sum(conjg(chi_arr)*chi_arr)
  jb=0.5*norm(chi_arr,squared=.true.)
!-----------------------------------------------------------------------
! Jobs = 0.5 * [ H(xb)+H_jac*dx-y ]^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  dep=innov+yn                  ! departures: h(xb)-y=H(xb)-y+H_jac*dx
# ifdef BIAS_CORRECTION
  bn=bias+dbias
do n=1,nobsData
  if(all(obsData(n)%iobs==0))cycle
  i=obsData(n)%iobs(0);j=obsData(n)%iobs(1)
  call CheckStop(any(ipar(i:j)/=n),"Inconsistent ipar/obsData")
  print "('#',3(1X,I0,1X,ES12.3))",&
    n,bn(n)*FGSCALE_INV,i,dep(i)*FGSCALE_INV,j,dep(j)*FGSCALE_INV
enddo
  dep=dep+bn(ipar(:nobs))       ! departures: h(xb)+bias(x,beta)-y
# endif
  Oinvdep=dep/(obsstddev**2)    ! O^{-1} * [ H(xb)+H_jac*dx-y ]
  jobs=0.5*dot_product(dep,Oinvdep)
! if(debug.and.MasterProc)print "(5(1X,A,'=',ES12.3))",&
!   '||H(xb)-y||',norm(innov),&
!   '||H_jac*dx||',norm(yn),&
!   '||H(xb)+H_jac*dx-y||',norm(dep),&
!   '||O^-1*[H(xb)+H_jac*dx-y]||',norm(Oinvdep),&
!   'Jo=1/2 * [H(xb)+H_jac*dx-y]^T * O^-1 * [H(xb)+H_jac*dx-y]',Jobs
!-----------------------------------------------------------------------
! Jbias = 0.5 * [ beta_b-beta ]^{T} * Bb^{-1} * [ beta_b-beta ]
!-----------------------------------------------------------------------
# ifdef BIAS_CORRECTION
  Binvdep=dbias*biasWeight/(biasStdDev**2)
  print "('#',1X,I0,2(1X,ES12.3))",&
    (n,dbias(n)*FGSCALE_INV,Binvdep(n)*FGSCALE_INV,n=1,nobsData)
  jbias=0.5*dot_product(dbias,Binvdep)
# else
  jbias=0.0
# endif
!-----------------------------------------------------------------------
! J = Jb + Jobs + Jbias
!-----------------------------------------------------------------------
  Jcost = jb + jobs + jbias
#endif
!-----------------------------------------------------------------------
! grad Jb = independent elements of \chi
!-----------------------------------------------------------------------
#ifdef OLD_MATCH_CODE
  gradJcost=0e0
  do i=1,nv2
    gradJcost(i)=chi(i)
  enddo
#else
  gradJcost(:)=chi(:)
!-----------------------------------------------------------------------
!  Bb^{-1} * [ beta_b-beta ] + bias_jac^{T} * O^{-1} * [ h(x)+bias(x,beta)-y ]
!-----------------------------------------------------------------------
# ifdef BIAS_CORRECTION 
  do n=1,nobsData
    i=obsData(n)%iobs(0);j=obsData(n)%iobs(1)
    gradJcost(nv2+n)=Binvdep(n)+sum(Oinvdep(i:j))
  enddo
# endif
#endif
!-----------------------------------------------------------------------
!  H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ]
!-----------------------------------------------------------------------
  dx=0e0
! do p=1,no4
!   i=pidx(p)
!   j=pidx(p+no4)
!   l=pidx(p+no8)
!   do k=1,nchemobs
!     kk=ichemobs(k)
!     do n=1,nobs
!       dx(i,j,l,k)=dx(i,j,l,k)+H_jac(n,p,kk)*Oinvdep(n)
!     enddo
!   enddo
! enddo
  do n=1,nobs
    do p=1,4
      i=pidx(n)%i(p)
      j=pidx(n)%j(p)
#ifdef OLD_MATCH_CODE
      do l=pidx(n)%l(0),pidx(n)%l(1)
        dx(i,j,l,:)=dx(i,j,l,:)+H_jac(n,p,l,:)*Oinvdep(n)
      enddo
#else
     l0=pidx(n)%l(0)
     l1=pidx(n)%l(1)
     dx(i,j,l0:l1,:)=dx(i,j,l0:l1,:)+H_jac(n,p,l0:l1,:)*Oinvdep(n)
#endif
    enddo
  enddo
  if(debug.and.MasterProc)then
    n=min(debug_n,nobs);p=min(debug_p,4)
    k=min(debug_k,nchemobs);if(obsData(ipar(n))%found)k=obsData(ipar(n))%ichem
    if(p==0)then
      print *,'xxx2:',norm(&!(/4,4,pidx(n)%l(1)-pidx(n)%l(0)+1/),&
                      dx(pidx(n)%i(:),pidx(n)%j(:),pidx(n)%l(0):pidx(n)%l(1),k)),&
                      Oinvdep(n)
    else
      print *,'xxx2:',dx(pidx(n)%i(p),pidx(n)%j(p),pidx(n)%l(1),k),Oinvdep(n)
    endif
  endif
!-----------------------------------------------------------------------
! grad Jobs = U^{-+}*H_jac^{T} * O^{-1} * [ H(xb)+H_jac*dx-y ],
! grad J = grad Jb + grad Jobs
!-----------------------------------------------------------------------
  call chitox_adj(chi_arr,dx)
! if(debug.and.MasterProc) print "(2(1X,A,'=',ES12.3))",&
!   '||dx||',norm(dx),'||chi_arr||',norm(chi_arr)
!-----------------------------------------------------------------------
! only add independent elements of grad Jobs to array gradJcost:
!-----------------------------------------------------------------------
  k=0
  do i=1,nv1
    do m=1,kx1
      k=k+1
      gradJcost(k)=gradJcost(k)+real(chi_arr(i,m,1))
    enddo
    do m=1,kx1
      do n=2,ky1
        k=k+1
        gradJcost(k)=gradJcost(k)+real(chi_arr(i,m,n))
      enddo
    enddo
    do m=kx2,nxex
      do n=2,kym1
        k=k+1
        gradJcost(k)=gradJcost(k)+real(chi_arr(i,m,n))
      enddo
    enddo
    do m=2,kxm1
      k=k+1
      gradJcost(k)=gradJcost(k)+aimag(chi_arr(i,m,1))
    enddo
    do n=2,kym1
      k=k+1
      gradJcost(k)=gradJcost(k)+aimag(chi_arr(i,1,n))
    enddo
    do m=2,kx1
      do n=2,ky1
        if(mod(nxex,2)==0.and.mod(nyex,2)==0.and.m==kx1.and.n==ky1)then
          continue
        else
          k=k+1
          gradJcost(k)=gradJcost(k)+aimag(chi_arr(i,m,n))
        endif
      enddo
    enddo
    do m=kx2,nxex
      do n=2,kym1
        k=k+1
        gradJcost(k)=gradJcost(k)+aimag(chi_arr(i,m,n))
      enddo
    enddo
  enddo
!-----------------------------------------------------------------------
  if(debug.and.MasterProc)then
    print "(5(1X,A,'=',ES12.3),$)",'J',Jcost,&
      'Jb',Jb,'Jo',Jobs,'Jbias',Jbias,'nobs/2',nobs*0.5
#ifdef BIAS_CORRECTION
    Jb=norm(chi(:nv2))
    Jbias=norm(chi(nv2+1:))
#else
    Jb=norm(chi)
    Jbias=0.0
#endif
    Jobs=DIM(norm(gradJcost),Jb+Jbias)
    print "(4(1X,A,'=',ES12.3))",'||DJb||',Jb,&
      '||DJo||',Jobs,'||DJbias||',Jbias,'||DJ||',Jb+Jobs+Jbias
  endif
 !ind=4
  call Add_2timing(43,tim_after,tim_before,'3DVar: costFunction.')
  call my_deallocate(.false.,"NoMessage")
!-----------------------------------------------------------------------
contains
!-----------------------------------------------------------------------
subroutine my_deallocate(verb,msg)
  logical, intent(in) :: verb
  character(len=*), intent(in) :: msg
  if(verb) print dafmt,msg
  if(allocated(chi_arr))  deallocate(chi_arr)
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
