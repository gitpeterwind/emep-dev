module ColumnSource_ml
!-----------------------------------------------------------------------!
! Emissions for volcanic (pasive) degassing and emergency (eEEMP) scenarios
!   Volcanic eruption: ASH and SO2 (optional)
!   NPP accident: Radiaton exposure tracers
! Any scenario can set emissions dates for other species as well.
! Varaible naming follows the Volcanic Eruption scenario, e.g.
!   loc: location of the source
!   ems: emission parameters
!-----------------------------------------------------------------------!
use CheckStop_ml,         only: CheckStop
use ChemSpecs,            only: NSPEC_TOT, NSPEC_SHL, species
use ChemGroups_ml,        only: chemgroups
use EmisDef_ml,           only: VOLCANOES_LL
use GridValues_ml,        only: xm2,sigma_bnd,GridArea_m2,&
                                GRIDWIDTH_M,&
                                coord_in_processor,coord_in_gridbox
use Io_ml,                only: open_file,read_line,IO_TMP,PrintLog
use MetFields_ml,         only: roa, z_bnd, u_xmj, v_xmi
use ModelConstants_ml,    only: KCHEMTOP,KMAX_MID,MasterProc,NPROC, &
                                USE_ASH,DEBUG,USE_PREADV,&
                                TXTLEN_NAME,dt_advec,dt_advec_inv,&
                                startdate,enddate
use NetCDF_ml,            only: GetCDF_modelgrid
use MPI_Groups_ml
use Par_ml,               only: LIMAX, LJMAX, me
use PhysicalConstants_ml, only: AVOG
use SmallUtils_ml,        only: wordsplit,find_index
use TimeDate_ml,          only: nydays,&           ! No. days per year
                                current_date,tdif_secs
use TimeDate_ExtraUtil_ml,only: date2string,string2date,to_stamp

implicit none
private

!** subroutines:
public :: ColumnRate      ! Emission rate
public :: getWinds   ! Wind speeds at locations

logical, save ::      &
  source_found=.true.,&   ! Are sources found on this processor/subdomain?
  topo_found=.false.      ! topo_nc file found? (vent elevation-model surface height)

integer, parameter :: &
  NMAX_LOC = 24, &  ! Max number of locations on processor/subdomain
  NMAX_EMS =6000   ! Max number of events def per location

integer, save ::   &        ! No. of ... found on processor/subdomain
  nloc             = -1,&   ! Source locations
  nems(0:NMAX_LOC) = -1     ! Events per location

character(len=*), parameter :: SDATE_FMT="YYYY-MM-DD hh:mm:ss"
integer, parameter :: SLEN=max(len(SDATE_FMT),TXTLEN_NAME)+1

type :: loc
  character(len=9)    :: id  =''      ! e.g. V1702A02B
  character(len=SLEN) :: name=''      ! e.g. Eyjafjöll
  real :: lat=-1.0,lon=-1.0,elev=-1.0 ! vent coords and elevation
  character(len=SLEN) :: etype=''     ! e.g. S0, 10kt
  integer :: grp=-1,iloc=-1,jloc=-1   ! Which (ash)goup,local i,j indes
end type loc
type(loc), save, dimension(NMAX_LOC):: &
  locdef=loc('UNDEF','UNKNOWN',-999.0,-999.0,-999.0,"??",-99,-99,-99)

type :: ems
  character(len=9)    :: id  =''      ! e.g. V1702A02B
  character(len=SLEN) :: name='',&    ! e.g. SO2
                         htype=''     ! Height reference/type for base/top:
                                      !   'VENT'=height from vent
                                      !   'SURF'=height from model surface
                                      !   'MLEV'=explicit model level
  real :: base=-1.0,top=-1.0,&        ! Column Base & Height [m asl]
          rate=-1.0                   ! Source strenght: Total release for period [kg/s]
  character(len=SLEN) :: &
    sbeg=SDATE_FMT,send=SDATE_FMT     ! event begin,end
  integer :: loc=-1,spc=-1,grp=-1     ! Which loc,(adv)spc,(ash)goup
  logical :: edef=.true.,&            ! default setings?
             dsec=.false.             ! correct rate by 1/secs(send-sbeg)
end type ems
type(ems), save,  allocatable ,dimension(:,:):: emsdef

character(len=*),parameter :: &
  mname = "ColumnSource",&
  topo_nc="topography.nc",&
  flocdef="columnsource_location.csv",   &  ! see locdef
  femsdef="columnsource_emission.csv",   &  ! see emsdef
  ERR_LOC_CSV=mname//" LOC def.file "//flocdef,          &
  ERR_EMS_CSV=mname//" EMS def.file "//femsdef,          &
  ERR_TOPO_NC=mname//" EMS def.file "//topo_nc,          &
  ERR_LOC_MAX=mname//" NMAX_LOC exceeded in "//flocdef,  &
  ERR_EMS_MAX=mname//" NMAX_EMS exceeded in "//femsdef,  &
  MSG_FMT="('"//mname//":',:,1X,A,5(:,1X,I0,':',A),3(:,1X,ES10.3,':',A))"

character(len=3), parameter :: &  ! Expand variable name for multy sceario runs
! EXPAND_SCENARIO_NAME(4)=["ASH","NUC","###","***"] ! e.g. ASH_F --> V1702A02B_F
  EXPAND_SCENARIO_NAME(1)=""                        ! do not expand

real,save,allocatable,dimension(:,:) :: surf_height ! [m], read from topo_nc
integer, public, save :: PROC_LOC(NMAX_LOC)=-1
real, allocatable, public, save  :: Winds(:,:,:)

contains
!-----------------------------------------------------------------------!
! Emergency scenarios:
!  Volcanic Eruption: Ash & SO2.
!  NPP Accident: Radiaton exposure tracers.
!  Any scenario can set emissions for ther species as well.
!-----------------------------------------------------------------------!
function ColumnRate(i,j,REDUCE_VOLCANO) result(emiss)
!-----------------------------------------------------------------------!
! Emissions from Emergency senarios.
!-----------------------------------------------------------------------!
  integer, intent(in) :: i,j
  real, intent(in), optional :: REDUCE_VOLCANO
  real, dimension(NSPEC_SHL+1:NSPEC_TOT,KCHEMTOP:KMAX_MID) :: emiss
  logical, save :: first_call=.true.
  character(len=SLEN)::&! Time strings in SDATE_FMT format
    sbeg=SDATE_FMT,&    ! Begin
    snow=SDATE_FMT,&    ! Now (current date)
    send=SDATE_FMT      ! End
  integer :: v,e,itot,k1,k0,k
  real    :: uconv
  integer, save          :: iSO2=-1
  integer, pointer, save :: iASH(:)=>null()
  real :: Ncells, frac
!----------------------------!
!
!----------------------------!
  if(first_call)then
    allocate(surf_height(LIMAX,LJMAX))
    call GetCDF_modelgrid("topography",topo_nc,surf_height,&
                          1,1,1,1,needed=.false.,found=topo_found)
    if(.not.topo_found)surf_height=0.0
    call setRate()
    deallocate(surf_height)
    first_call=.false.
    itot=find_index("SO2",species(:)%name)
    if(itot<1)then
      call PrintLog("WARNING: "//mname//" SO2 not found",MasterProc)
    else
      iSO2=itot
    end if
    itot=find_index("ASH",chemgroups(:)%name)
    if(itot<1)then
      call PrintLog("WARNING: "//mname//" ASH not found",MasterProc)
    else
      iASH=>chemgroups(itot)%specs
    end if

  end if
!----------------------------!
!
!----------------------------!
  emiss(:,:)=0.0
  if(.not.source_found)return
  snow=date2string(SDATE_FMT,current_date)

  doLOC: do v=1,nloc
     if(USE_PREADV)then
!cannot use formula below directly, because location may be in another subdomain
!           Winds(:,1,l)=u_xmj(locdef(l)%iloc,locdef(l)%jloc,:,1)*xm2(locdef(l)%iloc,locdef(l)%jloc)*dt_advec/gridwidth_m
!           Winds(:,2,l)=v_xmi(locdef(l)%iloc,locdef(l)%jloc,:,1)*xm2(locdef(l)%iloc,locdef(l)%jloc)*dt_advec/gridwidth_m
        !spread emissions in case of strong winds
        
        if(nems(v)<1) cycle doLOC                 ! Not erupting

    else
        if((i/=locdef(v)%iloc).or.(j/=locdef(v)%jloc) & ! Wrong gridbox
             .or.(nems(v)<1)) cycle doLOC                 ! Not erupting
     endif
    if(DEBUG%COLSRC .and. .not. USE_PREADV) &
      write(*,MSG_FMT)snow//' Vent',me,'me',v,trim(locdef(v)%id),i,"i",j,"j"
    doEMS: do e=1,nems(v)
      sbeg=date2string(emsdef(v,e)%sbeg,current_date)
      send=date2string(emsdef(v,e)%send,current_date)
      if(snow<sbeg.or.send<=snow)& ! Outside time window
        cycle doEMS
      itot=emsdef(v,e)%spc
      if(emsdef(v,e)%htype=="MLEV")then
         k0=emsdef(v,e)%base
         k1=emsdef(v,e)%top
      else
         k0=getModLev(i,j,emsdef(v,e)%base)
         k1=getModLev(i,j,emsdef(v,e)%top)
      end if
      uconv=1e-3                                          ! Kg/s --> ton/s=1e6 g/s
      if(emsdef(v,e)%dsec)uconv=1e6/max(dt_advec,&        ! Tg   --> ton/s=1e6 g/s
        tdif_secs(to_stamp(sbeg,SDATE_FMT),to_stamp(send,SDATE_FMT)))
      uconv=uconv/(GridArea_m2(i,j)*DIM(z_bnd(i,j,k1),z_bnd(i,j,k0+1))) ! --> g/s/cm3=1e-6 g/s/m3
      uconv=uconv*AVOG/species(itot)%molwt                              ! --> molecules/s/cm3
      if(USE_PREADV)then

         do k=k1,k0

            !only a fraction of the emission is used in each reachable gridcell
            !test if in range. NB: Winds have sign
            if(Winds(k,1,v)>=0.0 .and. ((i-locdef(v)%iloc)>floor(Winds(k,1,v)) .or. (i-locdef(v)%iloc)<0))cycle
            if(Winds(k,1,v)<=0.0 .and. ((i-locdef(v)%iloc)<floor(Winds(k,1,v)) .or. (i-locdef(v)%iloc)>0))cycle
            if(Winds(k,2,v)>=0.0 .and. ((j-locdef(v)%jloc)>floor(Winds(k,2,v)) .or. (j-locdef(v)%jloc)<0))cycle
            if(Winds(k,2,v)<=0.0 .and. ((j-locdef(v)%jloc)<floor(Winds(k,2,v)) .or. (j-locdef(v)%jloc)>0))cycle

            !test if along the line of wind, i.e. i/j = Winds(:,1,v)/Winds(:,2,v)
            if(nint((i-locdef(v)%iloc)*Winds(k,2,v))/=nint((j-locdef(v)%jloc)*Winds(k,1,v))) cycle
            Ncells=max(abs(Winds(k,1,v)),abs(Winds(k,2,v)))!NB: not an integer
            if(DEBUG%COLSRC) &
                 write(*,MSG_FMT)snow//' Vent',me,'me',v,trim(locdef(v)%id),i,"i",j,"j"
            if(Ncells<=1.0)then
               emiss(itot,k)=emiss(itot,k)+emsdef(v,e)%rate*uconv
            else
               frac=1.0/Ncells
               if(abs(i-locdef(v)%iloc)-floor(abs(Winds(k,1,v))) ==0 .and. &
                    abs(j-locdef(v)%jloc)-floor(abs(Winds(k,2,v)))==0)frac=frac*mod(Ncells,1.0)!last cell take only what is left
               if(DEBUG%COLSRC)write(*,*)'including fraction ',frac,' for ',i,j,k,v,locdef(v)%iloc,locdef(v)%jloc 
            endif
         enddo
      else
         emiss(itot,k1:k0)=emiss(itot,k1:k0)+emsdef(v,e)%rate*uconv
      endif
      if(DEBUG%COLSRC) &
        write(*,MSG_FMT)snow//' Erup.',me,'me',e,emsdef(v,e)%sbeg,&
          itot,trim(species(itot)%name),k1,'k1',k0,'k0',&
          emiss(itot,k1),'emiss',emsdef(v,e)%rate,'rate',uconv,'uconv'
    end do doEMS
  end do doLOC
!----------------------------!
! Volcanic emission reduction for SR run
!----------------------------!
  if(present(REDUCE_VOLCANO))then
    if(iSO2>0)&
      emiss(iSO2,:)=emiss(iSO2,:)*REDUCE_VOLCANO
    if(associated(iASH))&
      emiss(iASH,:)=emiss(iASH,:)*REDUCE_VOLCANO
  end if
!----------------------------!
! Disable Volcanic emissions
!----------------------------!
  if(.not.VOLCANOES_LL.and.iSO2>0)&       ! read from emislist.sox instead
    emiss(iSO2,:)=0.0
  if(.not.USE_ASH.and.associated(iASH))&  ! do not use ASH emissions
    emiss(iASH,:)=0.0
!----------------------------!
 contains
!----------------------------!
! Model level for a given height
!----------------------------!
function getModLev(i,j,height) result(k)
  integer, intent(in) :: i,j
  real, intent(in) :: height
  integer :: k
  k=KMAX_MID
  if(height<=0.0)return
  do while(k>0.and.height>z_bnd(i,j,k))
    k=k-1
  end do
end function getModLev
!----------------------------!
! Set Volcanic Eruption Param.
!----------------------------!
subroutine setRate()
  logical, save :: second_call=.true.
  character(len=SLEN*10):: txtline       ! Long enough for a full line
  type(loc) :: dloc
  type(ems) :: dems
  integer :: stat,l,v,e,g
! Particles classes & default split as London VAAC setup for NAME
! from Witham&al:2011 Table 1
!   Evolution of volcanic ash modelling at the London VAAC April 2010 --- April 2011
!   UK Met Office. Technical Summary (v1.0). May 2011.
!   C. Witham, M. Hort, D. Thomson, S. Leadbetter, B. Devenish and H. Webster.
  real, target :: & !0.0<0.1<0.3<1.0<3.0<10.0<30.0<100.0
    VAAC_7BIN_SPLIT(7)=[0.0,0.1,0.5,5.0,20.0,70.0,4.4]*1e-2,&
    VAAC_2BIN_SPLIT(2)=[            5.6,20.0         ]*1e-2,&
                    !0.0<4.0<6.0<8.0<10.0<12.0<14.0<16.0<18.0<25.0
    NILU_9BIN_SPLIT(9)=[0.157333,0.180197,0.151462,0.128799,0.101539,&
                        0.079647,0.060689,0.07127 ,0.069065],&
    NILU_1BIN_SPLIT(1)=[1.0]
  real, pointer, dimension(:) :: binsplit => NULL()

!----------------------------!
!
!----------------------------!
  if(first_call .and. USE_PREADV)allocate(Winds(KMAX_MID,2,NMAX_LOC))
  if(.not.first_call)then
    if(MasterProc.and.DEBUG%COLSRC.and.second_call) &
      write(*,MSG_FMT)'No need for reset volc.def...'
    second_call=.false.
    return
  end if
  first_call=.false.
  if(.not.allocated(emsdef)) then
     allocate(emsdef(0:NMAX_LOC,NMAX_EMS))
     emsdef(:,:)=ems('UNDEF','UNKNOWN','??',-999.0,-999.0,-999.0,&
                     '??','??',-1,-1,-1,.true.,.false.)
  end if
!----------------------------!
! Read Vent CVS
!----------------------------!
  if(DEBUG%COLSRC) CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
  if(MasterProc)then
    call open_file(IO_TMP,"r",flocdef,needed=.true.,iostat=stat)
    call CheckStop(stat,ERR_LOC_CSV//' not found')
  end if
  nloc=0
  l = 1
  doLOC: do while (l<=NMAX_LOC)
    call read_line(IO_TMP,txtline,stat)
    if(stat/=0) exit doLOC            ! End of file
    txtline=ADJUSTL(txtline)          ! Remove leading spaces
    if(txtline(1:1)=='#')cycle doLOC  ! Comment line
    dloc=getVent(txtline)
    if(coord_in_processor(dloc%lon,dloc%lat,iloc=dloc%iloc,jloc=dloc%jloc))then
      PROC_LOC(l) = ME!The source is located on this proc 
      nloc=nloc+1
      call CheckStop(nloc>NMAX_LOC,ERR_LOC_MAX//" read")
      ! remove model surface height from vent elevation
      dloc%elev=dloc%elev-surf_height(dloc%iloc,dloc%jloc)
      locdef(nloc)=dloc
      if(DEBUG%COLSRC) &
        write(*,MSG_FMT)'Vent',me,'in',nloc,trim(dloc%id),&
          dloc%grp,trim(dloc%name),dloc%iloc,"i",dloc%jloc,"j",&
          dloc%lon,"lon",dloc%lat,"lat"
    elseif(MasterProc)then
      PROC_LOC(l) = -1!The source is not located on this proc 
                      !we use -1, so that we know that if the sum all  PROC_LOC(l)<0, it is not in the rundomain
      if(DEBUG%COLSRC) &
        write(*,MSG_FMT)'Vent',me,'out',-1,trim(dloc%id),&
          dloc%grp,trim(dloc%name),dloc%iloc,"i",dloc%jloc,"j",&
          dloc%lon,"lon",dloc%lat,"lat"
    end if
    l = l+1
  end do doLOC
  if(MasterProc) close(IO_TMP)
  source_found=(nloc>0).or.(MasterProc.and.DEBUG%COLSRC)
  if(USE_PREADV)then
     !broadcast the PROC_LOC
     CALL MPI_ALLREDUCE(MPI_IN_PLACE,PROC_LOC,NMAX_LOC,MPI_INTEGER, &
          MPI_SUM,MPI_COMM_CALC,IERROR)
     do l=1,NMAX_LOC
        PROC_LOC(l) = PROC_LOC(l)+NPROC-1
        if(MasterProc .and. PROC_LOC(l)>=0)write(*,*)'source ',l,' on PROC ',PROC_LOC(l)
     enddo
     !Now PROC_LOC=ME, because
     !(ME) + (NPROC-1)*(-1) + NPROC-1= ME
     !if not in the rundomain, PROC_LOC = -1
     call getWinds

  endif
!----------------------------!
! Read Eruption CVS
!----------------------------!
  if(DEBUG%COLSRC) CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
  if(MasterProc)then
    call open_file(IO_TMP,"r",femsdef,needed=.true.,iostat=stat)
    call CheckStop(stat,ERR_EMS_CSV//' not found')
  end if
  nems(:)=0
  l = 1
  sbeg=date2string(SDATE_FMT,startdate)
  send=date2string(SDATE_FMT,enddate)
  doEMS: do while(l<=NMAX_EMS)
    call read_line(IO_TMP,txtline,stat)
    if(stat/=0) exit doEMS            ! End of file
    if(.not.source_found)cycle doEMS  ! There is no vents on sub-domain
    txtline=ADJUSTL(txtline)          ! Remove leading spaces
    if(txtline(1:1)=='#')cycle doEMS  ! Comment line
    dems=getErup(txtline)
    if(sbeg>date2string(dems%send,enddate  ).or.& ! starts after end of run
       send<date2string(dems%sbeg,startdate).or.& ! ends before start of run
       dems%rate<=0.0)then                        ! nothing to emit
      if(DEBUG%COLSRC.and.dems%loc>0) &
        write(*,MSG_FMT)'Erup.skip',me,'in',dems%loc,trim(dems%id),&
          0,trim(dems%sbeg),1,trim(dems%send)
      cycle doEMS
    elseif(dems%edef)then                                 ! Default
      nems(0)=nems(0)+1
      call CheckStop(nems(0)>NMAX_EMS,ERR_EMS_MAX//" read")
      emsdef(0,nems(dems%loc))=dems
      if(DEBUG%COLSRC) &
        write(*,MSG_FMT)'Erup.Default',me,'in',nems(0),trim(dems%id)
    elseif(dems%loc>0.and.(dems%spc>0.or.dems%grp>0))then ! Specific
      nems(dems%loc)=nems(dems%loc)+1
      call CheckStop(nems(dems%loc)>NMAX_EMS,ERR_EMS_MAX//" read")
      emsdef(dems%loc,nems(dems%loc))=dems
      if(DEBUG%COLSRC) &
        write(*,MSG_FMT)'Erup.Specific',me,'in',nems(dems%loc),trim(dems%id),&
          dems%spc,trim(dems%name),dems%grp,trim(dems%name)//"_GROUP"
      call CheckStop(dems%spc<1.and.dems%grp>0,"Erup.Specific unsupported group")
    elseif(MasterProc)then                         ! or Unknown Vent/SPC/GROUP
      if(DEBUG%COLSRC) &
        write(*,MSG_FMT)'Erup.Specific',me,'out',-1,trim(dems%id),&
          dems%spc,trim(dems%name),dems%grp,trim(dems%name)//"_GROUP"
    end if
    l = l+1
  end do doEMS
  if(MasterProc) close(IO_TMP)
  source_found=any(nems(1:nloc)>0)
!----------------------------!
! Expand Eruption Defaults
!----------------------------!
  if(DEBUG%COLSRC) CALL MPI_BARRIER(MPI_COMM_CALC, IERROR)
  if(nems(0)<1)then
  ! if(DEBUG%COLSRC) write(*,MSG_FMT)'Erup.Default',me,'not found'
    return
  end if
  doLOCe: do v=1,nloc
    if(nems(v)>0)cycle doLOCe   ! Specific found --> no need for Default
    e=nems(0)+1       ! A single default can have multiple lines, e.g.
    do                ! each line with a definition for a different specie
      e=find_index(locdef(v)%etype,emsdef(0,:e-1)%id)
      if(e<1)    cycle doLOCe   ! No Default found
      if(DEBUG%COLSRC) &
        write(*,MSG_FMT)'Erup.Default',me,'Expand',&
          v,trim(locdef(v)%id),e,trim(emsdef(0,e)%id)
      dems=emsdef(0,e)
      if(dems%htype=='VENT')then
!!      call CheckStop(.not.topo_found,ERR_TOPO_NC//' not found')
        dems%base=dems%base+locdef(v)%elev
        dems%top =dems%top +locdef(v)%elev
        if(DEBUG%COLSRC)&
          write(*,MSG_FMT)'Erup.Default',me,'Add loc%elev',&
            nint(dems%base),'ems%base',nint(dems%top),'ems%top'
      end if
      if(dems%spc<1.and.any(dems%name(1:3)==EXPAND_SCENARIO_NAME))then ! Expand variable name
        dems%name=trim(locdef(v)%id)//trim(dems%name(4:)) ! e.g. ASH_F --> V1702A02B_F
        dems%spc=find_index(dems%name,species(:)%name)    ! Specie (total)
        if(DEBUG%COLSRC)&
          write(*,MSG_FMT)'Erup.Default',me,'Expand',&
            dems%spc,trim(dems%name)
      end if
      if(dems%spc>0)then        ! Expand single SPC
        nems(v)=nems(v)+1
        call CheckStop(nems(v)>NMAX_EMS,ERR_EMS_MAX//" expand")
        emsdef(v,nems(v))=dems
        if(DEBUG%COLSRC) &
          write(*,MSG_FMT)'Erup.Default',me,'Expand SPC',nems(v),trim(dems%id)
      elseif(dems%grp>0)then    ! Expand GROUP of SPC
        select case (size(chemgroups(dems%grp)%specs))
        case(2);binsplit=>VAAC_2BIN_SPLIT
        case(7);binsplit=>VAAC_7BIN_SPLIT
        case(9);binsplit=>NILU_9BIN_SPLIT
        case(1);binsplit=>NILU_1BIN_SPLIT
        case default
          call CheckStop(ERR_EMS_CSV//' can not expand '//trim(locdef(v)%id))
        end select
        do g=1,size(chemgroups(dems%grp)%specs)
          dems%spc=chemgroups(dems%grp)%specs(g)  ! Specie (total)
          dems%name=species(dems%spc)%name
          dems%rate=emsdef(0,e)%rate*binsplit(g)
          nems(v)=nems(v)+1
          call CheckStop(nems(v)>NMAX_EMS,ERR_EMS_MAX//" expand")
          emsdef(v,nems(v))=dems
          if(DEBUG%COLSRC) &
          write(*,MSG_FMT)'Erup.Default',me,'Expand GRP',nems(v),trim(dems%id),&
            dems%spc,trim(dems%name)
        end do
      else
        if(DEBUG%COLSRC) &
          write(*,MSG_FMT)'Erup.Default',me,'not found'
      end if
    end do
  end do doLOCe
  source_found=any(nems(1:nloc)>0)
end subroutine setRate
!----------------------------!
! Extract Vent info from CVS line
!----------------------------!
function getVent(line) result(def)
  character(len=*)    :: line
  type(loc)           :: def
  character(len=SLEN) :: words(10)='' ! Array of paramaters
  real    :: lat,lon,elev             ! vent coord and elevation
  integer :: stat,nwords,igrp
  call wordsplit(line,size(words),words,nwords,stat,strict_separator=",",empty_words=.true.)
  call CheckStop(stat,"EMERGENCY: Wrong/Unknown line format "//trim(line))
  call CheckStop(nwords,size(words),"EMERGENCY: Missing data in line "//trim(line))
!#1:NUMBER,2:NAME,3:LOCATION,4:LATITUDE,5:NS,6:LONGITUDE,7:EW,8:ELEV,9:TYPE,10:ERUPTION TYPE
!V1702A02B,Eyjafjöll,Iceland-S,63.63,N,19.62,W,1666,Stratovolcano,S0
  read(words(4),*)lat
  select case (words(5)) ! NS
  case("N","n","degN")            ! degN
  case("S","s","degS");lat=-lat   ! degS
  case default
    call CheckStop("EMERGENCY: Unknown degN/S "//trim(words(5)))
  end select
  read(words(6),*)lon
  select case (words(7)) ! EW
  case("E","e","degE")            ! degE
  case("W","w","degW");lon=-lon   ! degW
  case default
    call CheckStop("EMERGENCY: Unknown degE/W "//trim(words(7)))
  end select
  read(words(8),*)elev            ! [m]
  igrp=find_index(words(1),chemgroups(:)%name)
  def=loc(trim(words(1)),trim(words(2)),lat,lon,elev,trim(words(10)),igrp)
end function getVent
!----------------------------!
! Extract Erup. info from CVS line
!----------------------------!
function getErup(line) result(def)
  character(len=*)    :: line
  type(ems)           :: def
  character(len=SLEN) :: words(10)=''   ! Array of paramaters
  logical :: edef=.true.,&              ! default setings?
             dsec=.false.               ! correct rate by 1/secs(send-sbeg)
  integer :: stat,nwords,iloc,ispc=0,igrp=0
  real    :: base,top,rate,frac,dhh
  call wordsplit(line,size(words),words,nwords,stat,strict_separator=",",empty_words=.true.)
  call CheckStop(stat,"EMERGENCY: Wrong/Unknown line format "//trim(line))
  call CheckStop(nwords,size(words),"EMERGENCY: Missing data in line "//trim(line))
!#1:TYPE/VOLCANO,2:VARIABLE,3:BASE[km],4:H[km above BASE],5:D[h],6:dM/dt[kg/s],7:m63[-],8:START[code/date],9:END[code/date],10:DESCRIPTION
!S0       ,     ,  , 11.000,   3.00, 4e6, 0.40,SR                 ,SR+D,Silicic standard
!V1702A02B,SO2  , 0,  8.000,  24.00,  15,     ,2010-04-14 00:00:00,SE+D,Eyja 20100414 SO2
!V1702A02B,ASH_F, 0,  2.000,  24.00,   0,     ,2010-05-23 00:00:00,SE+D,Eyja 20100523 PM fine
  iloc=find_index(words(1),locdef(:nloc)%id)                ! Vent Specific
  edef=(iloc<1).and.any(locdef(:nloc)%etype==words(1))      ! Vent Default
  if(iloc>0.and.any(words(2)(1:3)==EXPAND_SCENARIO_NAME))&  ! Expand variable name
    words(2)=trim(words(1))//trim(words(2)(4:)) ! e.g. ASH_F --> V1702A02B_F
  ispc=find_index(words(2),species(:)%name)     ! Specie (total)
  igrp=find_index(words(2),chemgroups(:)%name)  ! Group  (total)
  select case (words(3))        ! base
  case("MLEV","model")          ! Explicit model level
    words(3)="MLEV"
    read(words(4),*)top         ! [model level]
    base=top
  case("VENT"," ")              ! From the vent
    words(3)="VENT"
! vent specific: base/top from vent%elev
! emiss default: vent%elev is added on expansion (doLOCe: in setRate)
    base=0.0
    if(iloc>0)then
!!    call CheckStop(.not.topo_found,ERR_TOPO_NC//' not found')
      base=locdef(iloc)%elev    ! [m]
    end if
    read(words(4),*)top         ! [km]
    top=top*1e3                 ! [m]
    top=top+base                ! [m]
  case("SURF","0")              ! From the model surface
    words(3)="SURF"
    base=0.0
    read(words(4),*)top         ! [km]
    top=top*1e3                 ! [m]
  case default
    read(words(3),*)base        ! [km]
    base=base*1e3               ! [m]
    read(words(4),*)top         ! [km]
    top=top*1e3                 ! [m]
  end select
  read(words(6),*)rate
  select case (words(7))        ! m63 or effect.fraction
    case(" ")   ;frac=1.0
    case default;read(words(7),*)frac
  end select
  select case (words(5))        ! dt[h]
  case("1dt","1DT","1adv","1ADV")
    dhh=dt_advec/3600           ! only one time step
    frac=frac*dt_advec_inv      ! assume rate=total emission in [Kg]
    dsec=.false.
  case("total","TOTAL","event","EVENT")
    dsec=.true.
  case default
    read(words(5),*)dhh         ! assume rate in [Kg/s]
    dhh=max(dhh,dt_advec/3600)  ! at least 1 time step
    dsec=.false.
  end select
  words(8)=getDate(words(8),words(8),words(9),dhh,debug=DEBUG%COLSRC) ! Start [date/code]
  words(9)=getDate(words(9),words(8),words(9),dhh,debug=DEBUG%COLSRC) ! End   [date/code]
  def=ems(trim(words(1)),trim(words(2)),trim(words(3)),base,top,rate*frac,&
    trim(words(8)),trim(words(9)),max(iloc,0),max(ispc,0),max(igrp,0),edef,dsec)
end function getErup
!----------------------------!
! Time/Date CODE--> YYYY-MM-DD hh:mm:ss
!----------------------------!
function getDate(code,se,ee,dh,debug) result(str)
  character(len=*), intent(in)  :: code,se,ee
  real, intent(in)              :: dh  ! [hours]
  logical, intent(in), optional :: debug
  character(len=SLEN)           :: str
  logical :: dbg=.false.
  integer :: hh
  dbg=.false.;if(present(debug))dbg=debug.and..false.
  select case(code(1:4))
  case("SR")    ! Start of the simulation (model actually starts on 2nd time step)
    str=date2string(SDATE_FMT,startdate,addsecs=dt_advec,debug=dbg)
  case("SR+D")  ! Start of the simulation + dh
    str=date2string(SDATE_FMT,startdate,addsecs=dh*36e2+dt_advec,debug=dbg)
  case("SR+H")  ! Start of the simulation + Hhh hours
    read(code(5:6),*)hh
    str=date2string(SDATE_FMT,startdate,addsecs=hh*36e2+dt_advec,debug=dbg)
  case("SE+D")  ! Start eruption + dh; no wildcards allowed in SE
    str=date2string(SDATE_FMT,string2date(se,SDATE_FMT,debug=dbg),&
                    addsecs=dh*36e2,debug=dbg)
! case("SE+H")  ! Start eruption + Hhh; no wildcards allowed in SE
!   read(code(5:6),*)hh
!   str=date2string(SDATE_FMT,string2date(se,SDATE_FMT,debug=dbg),&
!                   addsecs=hh*36e2,debug=dbg)
  case("EE-D")  ! End eruption   - dh; no wildcards allowed in EE
    str=date2string(SDATE_FMT,string2date(ee,SDATE_FMT,debug=dbg),&
                    addsecs=-dh*36e2,debug=dbg)
  case("ER-D")  ! End of the simulation - dh
    str=date2string(SDATE_FMT,enddate,addsecs=-dh*36e2,debug=dbg)
  case("ER")    ! End of the simulation
    str=date2string(SDATE_FMT,enddate,debug=dbg)
  case default
    str=code
  end select
end function getDate
end function ColumnRate

subroutine getWinds
  integer ::l
  !broadcast wind speeds at emis locations
  do l=1,NMAX_LOC
     if(PROC_LOC(l)>=0)then
        !exchange wind speed. Converted to Courant number
        if(PROC_LOC(l) == ME)then
           Winds(:,1,l)=u_xmj(locdef(l)%iloc,locdef(l)%jloc,:,1)*xm2(locdef(l)%iloc,locdef(l)%jloc)*dt_advec/gridwidth_m
           Winds(:,2,l)=v_xmi(locdef(l)%iloc,locdef(l)%jloc,:,1)*xm2(locdef(l)%iloc,locdef(l)%jloc)*dt_advec/gridwidth_m
        endif
        CALL MPI_BCAST(Winds(1,1,l),2*8*KMAX_MID,MPI_BYTE,PROC_LOC(l),MPI_COMM_CALC,IERROR)
     endif
  enddo
end subroutine getWinds

endmodule ColumnSource_ml
