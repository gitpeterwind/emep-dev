PROGRAM unimod_B_nmc
  use DA_ml,             only: debug=>DA_DEBUG,DAFMT_DEF=>NMC_FMT_DEF,&
                               dafmt=>da_fmt_msg,damsg=>da_msg
  use ChemChemicals_ml,     only: define_chemicals,species    ! specie names
  use ChemGroups_ml,        only: Init_ChemGroups,chemgroups  ! group  names
  use CheckStop_ml,         only: CheckStop
  use SmallUtils_ml,        only: find_index
  use TimeDate_ml,          only: current_date
  use TimeDate_ExtraUtil_ml,only: date2nctime, nctime2date, date2string, nctime2string
  use Util_ml
  use exd_domain_ml
  use stddev_ml
  use covmat_ml,            only: dlon=>dxdim,dlat=>dydim, &
                                  set_chemobs_idx,allocate_covmat,&
                                  update_covmat,update_unobs_covmat,&
                                  normalise_covmat,diagonalise_covmat
  use spectralcov
  implicit none
!+------------------------------------------------------------------
  integer, parameter :: metnoID=88, maxLev=20, maxVar=100
  integer(4) :: inNml=172!, inFile=70, outFile=68
  logical, parameter :: debug_extdomain=.false.
  logical, parameter :: debug_ij=.false.
  integer, parameter :: dOut=10, ij_deb=1000
!+------------------------------------------------------------------
  integer, dimension(5) :: ncSDate, ncEDate
  integer(4) :: secSTime, secETime, secTime, secTime24, secTime48, &
                deltaHour(2), numNMC, ncFileID
  character(len=112) :: inFileName(2)=''!, ncFileName=''
  character(len=016) :: varName(maxVar)='',obsVarName(maxVar)=''
  logical            :: observedVar(maxVar)=.false.
! real(8) :: dlon=0.25, dlat=0.125
  namelist /NMC_VAR/ nChem, nChemObs, varName, obsVarName, observedVar
  namelist /NMC_CONFIG/ numNMC, nX,  nY, nLev,&
                    nXex, nYex, nex, nCorr, dlon, dlat
  namelist /NMC/ ncSDATE, ncEDATE, inFileName, deltaHour
  integer(4) :: h, k, rec, nnmc, nvar, ierr!, i, j, ij,
 !integer(4),    dimension(:), pointer :: time
 !real(8),       dimension(:), pointer :: lon, lat, lev
!+------------------------------------------------------------------
! integer(4) :: year,month,day,hour,minute,second,forecastStep,&
!               firstRec,lastRec
 !real(8), dimension(:,:,:),   pointer :: psrf,pres
  real(8), dimension(:,:,:,:,:), pointer :: var_ex=>NULL()
  real(8), dimension(:,:,:,:),   pointer :: var_nc=>NULL()
  real(8), dimension(:,:,:),     pointer :: domain=>NULL(),dom_ex=>NULL()
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
#ifdef gFortran
  open(unit=inNml,file='namelist.nml',status='OLD',action='READ',&
       form='FORMATTED',iostat=ierr)
#else
  open(unit=inNml,file='namelist.nml',status='OLD',action='READ',&
       form='FORMATTED',delim='APOSTROPHE',iostat=ierr)
#endif
  call io_check(ierr,'open namelist')
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  call define_chemicals()
  call Init_ChemGroups()
  dafmt=DAFMT_DEF
  k=find_index("DAOBS",chemgroups(:)%name)
  call CheckStop(k<1,'Unknown DAOBS group.')
  nChemObs=size(chemgroups(k)%ptr)
  obsVarName(:nChemObs)=species(chemgroups(k)%ptr)%name
  varName(:nChem)=(/obsVarName(:nChemObs),species(chemgroups(k)%ptr)%name/)
  k=find_index("DAUNOBS",chemgroups(:)%name)
 !CheckStop(k<1,'Unknown DAUNOBS group.')
  if(k>0)then
    nChem=nChemObs+size(chemgroups(k)%ptr(:))
    varName(:nChem)=(/obsVarName(:nChemObs),species(chemgroups(k)%ptr)%name/)
  else
    print dafmt,'WARNING Unknown DAUNOBS group'
    nChem=nChemObs
    varName(:nChem)=obsVarName(:nChemObs)
  endif
  observedVar(:)=.false.
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  deltaHour(1:2)=(/48,24/)
  read(unit=inNml,nml=NMC_VAR,iostat=ierr)
  call io_check(ierr,'read namelist: NMC_VAR')
  call CheckStop(nChemObs>0.eqv.any(observedVar),&
    'Incomplete/Redundant definition of nChemObs and observedVar on NMC_VAR namelist.')
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  do nvar=1,nChemObs
    k=find_index(obsVarName(nvar),varName(:nChem))
    call CheckStop(k<1,'Unknown observed variable: '//trim(obsVarName(nvar)))
    observedVar(k)=.true.
  enddo
  ! sort obsVarName following varName order
  nChemObs=0 ! count(observedVar)
  do nvar=1,nChem
    if(observedVar(nvar))then
      nChemObs=nChemObs+1
      obsVarName(nChemObs)=varName(nvar)
    endif
  enddo
  if(debug) then
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
#ifdef gFortran
  if(debug) write(*,nml=NMC_VAR)
#else
  if(debug) write(*,nml=NMC_VAR,delim='QUOTE')
#endif
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  nx=0;ny=0;nLev=0;nXex=0;nYex=0;dLon=0.0;dLat=0.0;
  read(unit=inNml,nml=NMC_CONFIG,iostat=ierr)
  call io_check(ierr,'read namelist: NMC_CONFIG')
  call CheckStop(any((/nX,nY,nLev/)<1),'Dimensions error: nX<1.or.nY<1.or.nLev<1.')
!+------------------------------------------------------------------
! EXTENDED DOMAIN FOR FFT: default size
!+------------------------------------------------------------------
!   NX=numLon
!   NY=numLat
!   NLEV=numLev
!   NCHEM=numVar
  if(nXex==0 .or. nYex==0)then
    print*,'WARNING: using extended domain default size'
    if(dLon==0.0) dLon=0.25
    if(dLat==0.0) dLat=0.25
    if(nXex==0) NXEX=NX+INT(12.0/dLon)
    if(nYex==0) NYEX=NY+INT(12.0/dLat)
  endif
  NEX=MAX(NXEX,NYEX)
  NCORR=((NLEV*NCHEM)*(NLEV*NCHEM+1))/2
#ifdef gFortran
  if(debug) write(*,nml=NMC_CONFIG)
#else
  if(debug) write(*,nml=NMC_CONFIG,delim='QUOTE')
#endif
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  !initialise FFT:
  call initfft()
  !initialise spectral space indices/variables:
  call initspec()
  nttot=0
!+------------------------------------------------------------------
!
!+------------------------------------------------------------------
  allocate(var_ex(nXex,nYex,nLev,2,nChem),stat=ierr)
  call CheckStop(ierr,'Allocation error: VAR_EX.')
!-----------------------------------------------
! FIRST TIME-LOOP FOR COMPUTING BIASES AND STD DEV
!-----------------------------------------------
  do nnmc=1,numNMC
    read(unit=inNml,nml=NMC,iostat=ierr)
    call io_check(ierr,'read namelist: NMC')
#ifdef gFortran
    if(debug) write(*,nml=NMC)
#else
    if(debug) write(*,nml=NMC,delim='QUOTE')
#endif
    call date2nctime(ncSDate,secSTime)
    call date2nctime(ncEDate,secETime)
    do secTime=secSTime,secETime,3600*24
      call nctime2date(current_date,secTime)
      dafmt=date2string(DAFMT_DEF,current_date)
      secTime24=secTime-deltaHour(1)*3600
      secTime48=secTime-deltaHour(2)*3600
!+------------------------------------------------------------------
      call nc_check(nf90_open(nctime2string(inFileName(1),secTime24,debug=debug),&
              nf90_nowrite,ncFileID))
      call GetNCDim(ncFileID,'lon' ,numLon)
      call GetNCDim(ncFileID,'lat' ,numLat)
      call GetNCDim(ncFileID,'k'   ,numLev)
      call CheckStop(any((/numLon,numLat/)/=(/NX,NY/)),&
                'Dimensions error: numLon/=NX.or.numLat/=NY for '//&
                 trim(inFileName(1)))
      call GetNCRec(ncFileID,current_date,rec,exact=debug)
      if(debug.and.nnmc==1.and.secTime==secSTime)call PrintNCDim(ncFileID,dOut)
      do nvar=1,nChem
        !24h forecast
        call GetNCVar(ncFileID,varName(nvar),var3D=var_nc,rec=rec,unitconv=FGSCALE)
        !extended domain for periodic boundary conditions:
        CALL EXT_DOMAIN(NX,NY,NLEV,NXEX,NYEX,var_nc(:,:,numLev-nlev+1:numLev,rec),var_ex(:,:,:,1,nvar))
        if(debug.and.debug_extdomain) then
          write(damsg,"(A,':',2(1X,A,':',I0))")'Debug 24h Values',&
            trim(varName(nvar)),nvar,'record',rec
          print dafmt,trim(damsg)
          print "((1X,I3,1X,A))",&
            (k,infovar(var_nc(:,:,numLev-nlev+k,rec),var_ex(:,:,k,1,nvar),&
                       trimdomain=.true.),k=1,nLev,nLev-1)
        endif
      enddo
      call nc_check(nf90_close(ncFileID))
      call nc_check(nf90_open(nctime2string(inFileName(2),secTime48,debug=debug),&
              nf90_nowrite,ncFileID))
      call GetNCRec(ncFileID,current_date,rec,exact=debug)
      if(debug.and.nnmc==1.and.secTime==secSTime) call PrintNCDim(ncFileID,dOut)
      do nvar=1,nChem
        !48h forecast
        call GetNCVar(ncFileID,varName(nvar),var3D=var_nc,rec=rec,unitconv=FGSCALE)
        !extended domain for periodic boundary conditions:
        CALL EXT_DOMAIN(NX,NY,NLEV,NXEX,NYEX,var_nc(:,:,numLev-nlev+1:numLev,rec),var_ex(:,:,:,2,nvar))
        if(debug.and.debug_extdomain) then
          write(damsg,"(A,':',2(1X,A,':',I0))")'Debug 48h Values',&
            trim(varName(nvar)),nvar,'record',rec
          print dafmt,trim(damsg)
          print "((1X,I3,1X,A))",&
            (k,infovar(var_nc(:,:,numLev-nlev+k,rec),var_ex(:,:,k,2,nvar),&
                       trimdomain=.true.),k=1,nLev,nLev-1)
        endif
      enddo
      call nc_check(nf90_close(ncFileID))
!+------------------------------------------------------------------
      !compute standard deviations:
      do nvar=1,nChem
        call update_stddev(nxex,nyex,nlev,2,nchem,nvar,current_date%hour,var_ex(:,:,:,:,nvar))
      enddo
!+------------------------------------------------------------------
    enddo
  enddo
  !normalise standard deviations:
  call normalise_stddev(nxex,nyex,nlev,nchem)
  if(debug) then
    h=current_date%hour
    do nvar=1,nChem
      write(damsg,"(A,':',1(1X,A,':',I0))"),'Debug NMC Statistics',&
        trim(varName(nvar)),nvar
      print dafmt,trim(damsg)
      print "((1X,I3,1X,A))",&
        (k,infovar(bias    (:,:,k,nvar,h),&
                   variance(:,:,k,nvar,h),&
                   stddev  (:,:,k,nvar)),k=1,nLev,nLev-1)
    enddo
  endif
  deallocate(variance)
!-----------------------------------------------
! SECOND TIME-LOOP FOR COMPUTING SPECTRAL COVARIANCES
!-----------------------------------------------
  call set_chemobs_idx(nchem,observedVar(1:nchem))
  call allocate_covmat(nex,nlev,nkstar)

  nvar=1
  rewind(unit=inNml,iostat=ierr)
  call io_check(ierr,'rewind namelist')
  do nnmc=1,numNMC
    read(unit=inNml,nml=NMC,iostat=ierr)
    call io_check(ierr,'read namelist: NMC')
#ifdef gFortran
    if(debug) write(*,nml=NMC)
#else
    if(debug) write(*,nml=NMC,delim='QUOTE')
#endif
    call date2nctime(ncSDate,secSTime)
    call date2nctime(ncEDate,secETime)
    do secTime=secSTime,secETime,3600*24
      call nctime2date(current_date,secTime)
      dafmt=date2string(DAFMT_DEF,current_date)
      secTime24=secTime-deltaHour(1)*3600
      secTime48=secTime-deltaHour(2)*3600
!+------------------------------------------------------------------
      call nc_check(nf90_open(nctime2string(inFileName(1),secTime24,debug=debug),&
              nf90_nowrite,ncFileID))
      call GetNCRec(ncFileID,current_date,rec,exact=debug)
      do nvar=1,nChem
        !24h forecast
        call GetNCVar(ncFileID,varName(nvar),var3D=var_nc,rec=rec,unitconv=FGSCALE)
        !extended domain for periodic boundary conditions:
        CALL EXT_DOMAIN(NX,NY,NLEV,NXEX,NYEX,var_nc(:,:,numLev-nlev+1:numLev,rec),var_ex(:,:,:,1,nvar))
      enddo
      call nc_check(nf90_close(ncFileID))
      call nc_check(nf90_open(nctime2string(inFileName(2),secTime48,debug=debug),&
              nf90_nowrite,ncFileID))
      if(debug.and.nnmc==1.and.secTime==secSTime) call PrintNCDim(ncFileID,dOut)
      call GetNCRec(ncFileID,current_date,rec,exact=debug)
      do nvar=1,nChem
        !48h forecast
        call GetNCVar(ncFileID,varName(nvar),var3D=var_nc,rec=rec,unitconv=FGSCALE)
        !extended domain for periodic boundary conditions:
        CALL EXT_DOMAIN(NX,NY,NLEV,NXEX,NYEX,var_nc(:,:,numLev-nlev+1:numLev,rec),var_ex(:,:,:,2,nvar))
      enddo
      call nc_check(nf90_close(ncFileID))
      !compute normalised model errors:
      call get_moderr(nxex,nyex,nlev,2,nchem,1,current_date%hour,var_ex)
      !Fourier-transform model errors from physical to spectral space:
      do nvar=1,nChem
        do k=1,nlev
          call cplx_fft_2d(nxex,nyex,nmin,nmax,var_ex(:,:,k,1,nvar),var_ex(:,:,k,2,nvar),wsave,lensav)
          if(debug) then
            if(k==1)then
              write(damsg,"(A,':',2(1X,A,':',I0))"),'Debug cplx_fft_2d',&
                trim(varName(nvar)),nvar,'record',rec
              print dafmt,trim(damsg)
            endif
            print "((1X,I3,1X,A))",&
              k,infovar(var_ex(:,:,k,1,nvar),var_ex(:,:,k,2,nvar))
          endif
        enddo
      enddo
!     call update_covmat(nxex,nyex,nex,nlev,2,nchem,nttot,var_ex,&
!                        nmin,nmax,ikstar,nkstar,nstar,intweight,mt,nt)
      call update_covmat(nxex,nyex,    nlev,2,nchem,nttot,var_ex)
!     call update_unobs_covmat(nxex,nyex,nex,nlev,2,nchem,var_ex,&
!                       nmin,nmax,nkstar,nstar,intweight,mt,nt)
      call update_unobs_covmat(nxex,nyex,    nlev,2,nchem,var_ex)
    enddo
!+------------------------------------------------------------------
  enddo
  deallocate(bias)

  !normalise covariances
! call normalise_covmat(nex,nxex,nyex,nlev,nttot,nkstar,nstar,kstar,ikstar)
  call normalise_covmat(    nxex,nyex,nlev,nttot)
  if(allocated(ucovmat))deallocate(ucovmat)
!-----------------------------------------------
! DIAGONALISE COVARIANCE MATRICES
! (one cov-matrix for each wavenumber kstar):
!-----------------------------------------------
! call diagonalise_covmat(nex,nxex,nyex,nx,ny,nlev,nchem,&
!       nkstar,kstar,ikstar,nstarmax)
  call diagonalise_covmat(nex,nxex,nyex,nx,ny,nlev,nchem,&
        nkstar,             nstarmax)
END PROGRAM unimod_B_nmc
