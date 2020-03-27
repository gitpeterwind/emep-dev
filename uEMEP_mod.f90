module uEMEP_mod
!
! all subroutines for uEMEP
!
use CheckStop_mod,     only: CheckStop,StopAll
use Chemfields_mod,    only: xn_adv
use ChemDims_mod,      only: NSPEC_ADV, NSPEC_SHL,NEMIS_File
use ChemSpecs_mod,     only: species_adv,species
use Country_mod,       only: MAXNLAND,NLAND,Country
use EmisDef_mod,       only: lf, lf_emis, lf_emis_tot, loc_frac_src_1d,&
                            lf_src_acc,lf_src_tot,lf_src_full,loc_tot_full, NSECTORS,EMIS_FILE, &
                            nlandcode,landcode,sec2tfac_map,sec2hfac_map, sec2split_map,&
                            ISNAP_DOM,secemis, roaddust_emis_pot,KEMISTOP,&
                            NEmis_sources, Emis_source_2D, Emis_source
use EmisGet_mod,       only: nrcemis, iqrc2itot, emis_nsplit,nemis_kprofile, emis_kprofile,&
                             make_iland_for_time,itot2iqrc,iqrc2iem, emisfrac
use GridValues_mod,    only: dA,dB,xm2, dhs1i, glat, glon, projection, extendarea_N,i_fdom,j_fdom
use MetFields_mod,     only: ps,roa,EtaKz
use Config_module,     only: KMAX_MID, KMAX_BND,USES, uEMEP, IOU_HOUR&
                             , IOU_HOUR_INST,IOU_INST,IOU_YEAR,IOU_MON,IOU_DAY&
                             ,IOU_HOUR,IOU_HOUR_INST, IOU_MAX_MAX, KMAX_MID &
                             ,MasterProc,dt_advec, RUNDOMAIN, runlabel1 &
                             ,HOURLYFILE_ending
use MPI_Groups_mod
use NetCDF_mod,        only: Real4,Out_netCDF
use OwnDataTypes_mod,  only: Deriv, Npoll_uemep_max, Nsector_uemep_max, TXTLEN_FILE
use Par_mod,           only: me,LIMAX,LJMAX,MAXLIMAX,MAXLJMAX,gi0,gj0,li0,li1,lj0,lj1,GIMAX,GJMAX
use PhysicalConstants_mod, only : GRAV, ATWAIR 
use SmallUtils_mod,    only: find_index
use TimeDate_mod,      only: date, current_date,day_of_week
use TimeDate_ExtraUtil_mod,only: date2string
use Timefactors_mod,   only: &
    DegreeDayFactors       & ! degree-days used for SNAP-2
    ,Gridded_SNAP2_Factors, gridfac_HDD & 
    ,GridTfac &!array with monthly gridded time factors
    ,fac_min,timefactors   &                  ! subroutine
    ,fac_ehh24x7 ,fac_emm, fac_edd, timefac  ! time-factors
use My_Timing_mod,     only: Add_2timing, Code_timer, NTIMING

!(dx,dy,i,j) shows contribution of pollutants from (i+dx,j+dy) to (i,j)

implicit none
!external advection_mod_mp_vertdiffn_k

private

public  :: init_uEMEP
public  :: out_uEMEP
public  :: av_uEMEP
public  :: uemep_adv_x
public  :: uemep_adv_y
public  :: uemep_adv_k
public  :: uemep_diff
public  :: uemep_emis
public  :: add_lf_emis

integer, public, save :: uEMEP_Size1=0 !total size of the first 3 dimensions of loc_frac
!uEMEP_Size1=uEMEP%Nsec_poll*(2*uEMEP%dist+1)**2
integer, public, save :: uEMEP_Sizedxdy=0 !total size of the first 3 dimensions of loc_frac
!uEMEP_Sizedxdy=(2*uEMEP%dist+1)**2

real, private, save ::av_fac_hour,av_fac_day,av_fac_month,av_fac_full
real, allocatable, save ::loc_poll_to(:,:,:,:,:)

logical, public, save :: COMPUTE_LOCAL_TRANSPORT=.false.
integer , private, save :: uEMEPNvertout = 1!number of vertical levels to save in output
integer, public, save :: NTIMING_uEMEP=7
real, private :: tim_after,tim_before
integer, public, save :: Ndiv_coarse=0, Ndiv_rel=0, Ndiv2_coarse=0
integer, public, save :: Nsources=0

type, public :: lf_sources
  character(len=4) :: emis = 'NONE' !pollutants to include (NB: 4 char)
  character(len=20) :: type = 'relative' !Qualitatively different type of sources: "coarse", "relative", "country"
  integer :: dist = -1 ! window dimension, if defined 
  integer :: Nvert = -1 ! vertical extend of the tracking/local rwindow
  integer :: sector=-1    ! sector for this source. Zero is sum of all sectors
  integer :: poll = 1 !index of pollutant in loc_tot (set by model)
  integer :: start = 1 ! first position index in lf_src (set by model)
  integer :: end = 1 ! last position index in lf_src (set by model)
  integer :: iem = 0 ! index of emitted pollutant, emis (set by model)
  integer :: Npos = 0 ! number of position indices in lf_src (set by model)
  integer :: Nsplit = 0 ! into how many species the emitted pollutant is split into (set by model)
  integer :: species_ix = -1 !species index, if single pollutant (for example NO or NO2, instead of nox) 
  integer :: iqrc = -1 !index for emissplits, if single pollutant (for example NO or NO2, instead of nox) 
  integer, dimension(4) :: DOMAIN = -1 ! DOMAIN which will be outputted
  integer, dimension(15) :: ix = -1 ! internal index of the  (splitted) species (set by model)
  integer, dimension(15) :: mw=0.0  ! molecualr weight of the (splitted) species (set by model)
end type lf_sources
integer, public, parameter :: MAXSRC=100
type(lf_sources), public, save :: lf_src(MAXSRC)
integer, public, save :: LF_SRC_TOTSIZE
integer,  public, save :: iotyp2ix(IOU_MAX_MAX)
integer,  public, save :: Niou_ix = 0 ! number of time periods to consider (hourly, monthly, full ...)
integer,  public, save :: Npoll = 0 !Number of different pollutants to consider
integer,  public, save :: iem2ipoll(NEMIS_File)

contains
subroutine init_uEMEP
  integer :: i, ix, itot, iqrc, iem, iemis, isec, ipoll, ixnh3, ixnh4, size, IOU_ix, isrc


  call Code_timer(tim_before)
  ix=0
  if(lf_src(1)%emis == 'NONE')then
     !Temporary: we keep compatibilty with uEMEP input
     lf_src(:)%dist = uEMEP%dist !Temporary
     lf_src(:)%Nvert = uEMEP%Nvert !Temporary
     lf_src(:)%DOMAIN(1) = uEMEP%DOMAIN(1) !Temporary
     lf_src(:)%DOMAIN(2) = uEMEP%DOMAIN(2) !Temporary
     lf_src(:)%DOMAIN(3) = uEMEP%DOMAIN(3) !Temporary
     lf_src(:)%DOMAIN(4) = uEMEP%DOMAIN(4) !Temporary
     Ndiv_rel = 2*lf_src(1)%dist+1
     Ndiv_coarse = 2*lf_src(1)%dist+1
     Ndiv2_coarse = Ndiv_coarse*Ndiv_coarse
     uEMEP%Nsec_poll = 0
     uEMEP%Npoll = 0
     do isrc=1,Npoll_uemep_max
        if(uEMEP%poll(isrc)%emis=='none')then
           call CheckStop(isrc==1,"init_uEMEP: no pollutant specified")
           exit
        else
           uEMEP%Npoll = uEMEP%Npoll + 1
           uEMEP%poll(isrc)%Nsectors = 0
           uEMEP%poll(isrc)%sec_poll_ishift=uEMEP%Nsec_poll
           do isec=1,Nsector_uemep_max
              if(uEMEP%poll(isrc)%sector(isec)<0)then
                 call CheckStop(isec==0,"init_uEMEP: nosector specified for "//uEMEP%poll(isrc)%emis)
                 exit
              else
                 uEMEP%Nsec_poll = uEMEP%Nsec_poll + 1
                 uEMEP%poll(isrc)%Nsectors = uEMEP%poll(isrc)%Nsectors +1
                 ix=ix+1
                 lf_src(ix)%emis = uEMEP%poll(isrc)%emis
                 lf_src(ix)%sector = uEMEP%poll(isrc)%sector(isec)
                 lf_src(ix)%Npos = (2*lf_src(ix)%dist+1)* (2*lf_src(ix)%dist+1)     
              endif
           enddo
        endif
     enddo
  endif
  
  Nsources = 0
  do i = 1, Npoll_uemep_max
     if(lf_src(i)%emis /= 'NONE') Nsources = Nsources + 1
  enddo

  ipoll=0
  iem2ipoll = -1
  do isrc = 1, Nsources
     do i=1,4
        if(lf_src(i)%DOMAIN(i)<0)lf_src(i)%DOMAIN(i) = RUNDOMAIN(i)
     enddo
     iem=find_index(lf_src(isrc)%emis ,EMIS_FILE(1:NEMIS_FILE))
     if(iem<1)then
        !defined as single species
        ix=find_index(lf_src(isrc)%emis ,species(:)%name)
        if(ix<0)then
           ix=find_index(lf_src(isrc)%emis ,species(:)%name, any_case=.true.)
           if(me==0 .and. ix>0)then
              write(*,*)'WARNING: '//trim(lf_src(isrc)%emis)//' not found, replacing with '//trim(species(ix)%name)
              lf_src(isrc)%emis=trim(species(ix)%name)
           endif
        endif
        call CheckStop( ix<1, "uEMEP did not find corresponding pollutant: "//trim(lf_src(isrc)%emis) )
        iem=-1
        lf_src(isrc)%species_ix = ix
        lf_src(isrc)%iqrc = itot2iqrc(ix)        
     endif

     lf_src(isrc)%iem = iem

     if(iem2ipoll(iem)<0)then
        ipoll = ipoll + 1
        iem2ipoll(iem) = ipoll     
        Npoll = ipoll
     endif
     lf_src(isrc)%poll = iem2ipoll(iem)
     
     lf_src(isrc)%Nsplit=emis_nsplit(iem)
     do i=1,lf_src(isrc)%Nsplit
        iqrc=sum(emis_nsplit(1:iem-1)) + i
        itot=iqrc2itot(iqrc)
        ix=itot-NSPEC_SHL
        lf_src(isrc)%ix(i)=ix
        lf_src(isrc)%mw(i)=species_adv(ix)%molwt
        if(lf_src(isrc)%emis=="nox ")then
           ix=find_index("NO2",species_adv(:)%name)
           call CheckStop(ix<0,'Index for NO2 not found')
           lf_src(isrc)%mw(i)=species_adv(ix)%molwt
        endif
        if(lf_src(isrc)%emis=="sox ")then
           ix=find_index("SO2",species_adv(:)%name)
           call CheckStop(ix<0,'Index for SO2 not found')
           lf_src(isrc)%mw(i)=species_adv(ix)%molwt
        endif

        if(lf_src(isrc)%emis=="nh3 ")then
           lf_src(isrc)%Nsplit = 0
           ixnh4=find_index("NH4_F",species_adv(:)%name , any_case=.true.)
           ixnh3=find_index("NH3",species_adv(:)%name)
           do ix=1,NSPEC_ADV
              if(ix/=ixnh4.and.ix/=ixnh3)cycle!not reduced nitrogen
              if(species_adv(ix)%nitrogens>0)then
                 lf_src(isrc)%Nsplit = lf_src(isrc)%Nsplit + 1
                 lf_src(isrc)%ix(lf_src(isrc)%Nsplit) = ix
                 lf_src(isrc)%mw(lf_src(isrc)%Nsplit) = species_adv(ixnh3)%molwt !use NH3 mw also for NH4
              endif
           enddo
        endif
     end do

     if(MasterProc)then
        write(*,*)'uEMEP pollutant : ',lf_src(isrc)%emis
        write(*,*)'uEMEP number of species in '//trim(lf_src(isrc)%emis)//' group: ',lf_src(isrc)%Nsplit
!        write(*,"(A,30(A,F6.2))")'including:',('; '//trim(species_adv(lf_src(isrc)%ix(i))%name)//', mw ',lf_src(isrc)%mw(i),i=1,lf_src(isrc)%Nsplit)
!        write(*,"(A,30I4)")'sectors:',(uEMEP%poll(ipoll)%sector(i),i=1,uEMEP%poll(ipoll)%Nsectors)
!        write(*,"(A,30I4)")'ix:',(uEMEP%poll(ipoll)%ix(i),i=1,uEMEP%poll(ipoll)%Nix)
     end if
  end do

!find indices in EMIS_File
  do ipoll=1,uEMEP%Npoll        
     do iem=1,NEMIS_FILE 
        if(trim(EMIS_File(iem))==trim(uEMEP%poll(ipoll)%emis))then
           uEMEP%poll(ipoll)%EMIS_File_ix = iem
           exit
        endif
     enddo
     call CheckStop(iem>NEMIS_FILE,"uemep pollutant not found: "//trim(uEMEP%poll(ipoll)%emis))
  enddo
 
  do i=1,4
     if(uEMEP%DOMAIN(i)<0)uEMEP%DOMAIN(i) = RUNDOMAIN(i)
  enddo

  uEMEP_Sizedxdy = (2*lf_src(1)%dist+1)**2
  uEMEP_Size1 = uEMEP%Nsec_poll*uEMEP_Sizedxdy

  if(MasterProc)then
     write(*,*)'uEMEP pollutants : ',uEMEP%Npoll
     write(*,*)'total uEMEP pollutants and sectors : ',uEMEP%Nsec_poll
  end if
  do ipoll=1,uEMEP%Npoll        
     iem=find_index(uEMEP%poll(ipoll)%emis ,EMIS_FILE(1:NEMIS_FILE))
     call CheckStop( iem<1, "uEMEP did not find corresponding emission file: "//trim(uEMEP%poll(ipoll)%emis) )
     call CheckStop( iem/=uEMEP%poll(ipoll)%EMIS_File_ix, "uEMEP wrong emis file index for: "//trim(uEMEP%poll(ipoll)%emis) )
     uEMEP%poll(ipoll)%Nix=emis_nsplit(iem)
     do i=1,uEMEP%poll(ipoll)%Nix
        iqrc=sum(emis_nsplit(1:iem-1)) + i
        itot=iqrc2itot(iqrc)
        ix=itot-NSPEC_SHL
        uEMEP%poll(ipoll)%ix(i)=ix
        uEMEP%poll(ipoll)%mw(i)=species_adv(ix)%molwt
        if(uEMEP%poll(ipoll)%emis=="nox ")then
           ix=find_index("NO2",species_adv(:)%name)
           call CheckStop(ix<0,'Index for NO2 not found')
           uEMEP%poll(ipoll)%mw(i)=species_adv(ix)%molwt
        endif
        if(uEMEP%poll(ipoll)%emis=="sox ")then
           ix=find_index("SO2",species_adv(:)%name)
           call CheckStop(ix<0,'Index for SO2 not found')
           uEMEP%poll(ipoll)%mw(i)=species_adv(ix)%molwt
        endif

     if(uEMEP%poll(ipoll)%emis=="nh3 ")then
        uEMEP%poll(ipoll)%Nix=0
        ixnh4=find_index("NH4_F",species_adv(:)%name , any_case=.true.)
        ixnh3=find_index("NH3",species_adv(:)%name)
        do ix=1,NSPEC_ADV
           if(ix/=ixnh4.and.ix/=ixnh3)cycle!not reduced nitrogen
           if(species_adv(ix)%nitrogens>0)then
              uEMEP%poll(ipoll)%Nix =  uEMEP%poll(ipoll)%Nix + 1
              uEMEP%poll(ipoll)%ix(uEMEP%poll(ipoll)%Nix)=ix
              uEMEP%poll(ipoll)%mw(uEMEP%poll(ipoll)%Nix)=species_adv(ixnh3)%molwt!use NH3 mw also for NH4
           endif
        enddo
     endif
     end do
     if(MasterProc)then
        write(*,*)'uEMEP pollutant : ',uEMEP%poll(ipoll)%emis
        write(*,*)'uEMEP number of species in '//trim(uEMEP%poll(ipoll)%emis)//' group: ',uEMEP%poll(ipoll)%Nix
        write(*,"(A,30(A,F6.2))")'including:',('; '//trim(species_adv(uEMEP%poll(ipoll)%ix(i))%name)//', mw ',uEMEP%poll(ipoll)%mw(i),i=1,uEMEP%poll(ipoll)%Nix)
        write(*,"(A,30I4)")'sectors:',(uEMEP%poll(ipoll)%sector(i),i=1,uEMEP%poll(ipoll)%Nsectors)
        write(*,"(A,30I4)")'ix:',(uEMEP%poll(ipoll)%ix(i),i=1,uEMEP%poll(ipoll)%Nix)
     end if
  end do

  COMPUTE_LOCAL_TRANSPORT = uEMEP%COMPUTE_LOCAL_TRANSPORT

  av_fac_hour=0.0
  av_fac_day=0.0
  av_fac_month=0.0
  av_fac_full=0.0
  
  LF_SRC_TOTSIZE = 0
  do isrc = 1, Nsources
     lf_src(isrc)%start = LF_SRC_TOTSIZE + 1
     lf_src(isrc)%end = LF_SRC_TOTSIZE + lf_src(isrc)%Npos
     LF_SRC_TOTSIZE = LF_SRC_TOTSIZE + lf_src(isrc)%Npos
     if(me==0)then
        write(*,*)isrc,' ',trim(lf_src(isrc)%emis)," start ",lf_src(isrc)%start," end ",lf_src(isrc)%end,LF_SRC_TOTSIZE
     endif
  enddo
  allocate(lf(LF_SRC_TOTSIZE,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
  lf=0.0

  if(uEMEP%HOUR)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_HOUR)=Niou_ix
  endif
  if(uEMEP%HOUR_INST)then
     !Niou_ix = Niou_ix + 1 !should not be accumulated
     iotyp2ix(IOU_HOUR_inst) = -1; !should not be accumulated
  endif  
  if(uEMEP%DAY)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_DAY)=Niou_ix
  endif
  if(uEMEP%MONTH)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_MON)=Niou_ix
  endif
  if(uEMEP%YEAR)then
     Niou_ix = Niou_ix + 1
     iotyp2ix(IOU_YEAR)=Niou_ix
  endif
  allocate(lf_src_acc(LF_SRC_TOTSIZE,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,Niou_ix))
  lf_src_acc = 0.0
  allocate(lf_src_tot(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,Npoll,Niou_ix))
  lf_src_tot = 0.0
   allocate(loc_frac_src_1d(LF_SRC_TOTSIZE,0:max(LIMAX,LJMAX)+1))
  loc_frac_src_1d=0.0
  allocate(lf_emis(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,Nsources))
  lf_emis = 0.0
  allocate(lf_emis_tot(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID,Npoll))
  lf_emis_tot = 0.0
 
  call Add_2timing(NTIMING-9,tim_after,tim_before,"uEMEP: init")
end subroutine init_uEMEP


subroutine out_uEMEP(iotyp)
  integer, intent(in) :: iotyp
  character(len=200) ::filename, varname
  real :: xtot,scale,invtot,t1,t2
  integer ::i,j,k,n,n1,dx,dy,ix,iix,isec,iisec,isec_poll,ipoll,isec_poll1,isrc,iou_ix,iter
  integer ::ndim,kmax,CDFtype,dimSizes(10),chunksizes(10)
  integer ::ndim_tot,dimSizes_tot(10),chunksizes_tot(10)
  character (len=20) ::dimNames(10),dimNames_tot(10)
  type(Deriv) :: def1 ! definition of fields
  type(Deriv) :: def2 ! definition of fields
  logical ::overwrite, create_var_only
  logical,save :: first_call(10)=.true.
  real,allocatable ::tmp_out(:,:,:)!allocate since it may be heavy for the stack TEMPORARY
  type(date) :: onesecond = date(0,0,0,0,1)
  character(len=TXTLEN_FILE),save :: oldhourlyname = 'NOTSET'
  character(len=TXTLEN_FILE),save :: oldhourlyInstname = 'NOTSET'
  character(len=TXTLEN_FILE),save :: oldmonthlyname
  real :: fracsum(LIMAX,LJMAX)
  logical :: pollwritten(Npoll_uemep_max)
  
  call Code_timer(tim_before)

  if(iotyp==IOU_HOUR_INST .and. uEMEP%HOUR_INST)then
     fileName = trim(runlabel1)//'_uEMEP_hourInst'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
     if(oldhourlyInstname/=fileName)then
        first_call(iotyp) = .true.
        oldhourlyInstname = fileName
     endif
  else if(iotyp==IOU_HOUR .and. uEMEP%HOUR)then
     fileName = trim(runlabel1)//'_uEMEP_hour'//date2string(trim(HOURLYFILE_ending),current_date,-1.0)
     if(oldhourlyname/=fileName)then
        first_call(iotyp) = .true.
        oldhourlyname = fileName
     endif
  else if(iotyp==IOU_DAY .and. uEMEP%DAY)then
     fileName=trim(runlabel1)//'_uEMEP_day.nc'
  else if(iotyp==IOU_MON .and. uEMEP%MONTH)then
     if(uEMEP%MONTH_ENDING /= "NOTSET")then
        fileName=trim(runlabel1)//'_uEMEP_month'//date2string(trim(uEMEP%MONTH_ENDING),current_date,-1.0)
        if(oldmonthlyname/=fileName)then
           first_call(iotyp) = .true.
           oldmonthlyname = fileName
        endif
     else
        fileName=trim(runlabel1)//'_uEMEP_month.nc'
     endif
  else if(iotyp==IOU_YEAR .and. uEMEP%YEAR)then
     fileName=trim(runlabel1)//'_uEMEP_full.nc'
  else
     return
  endif
  ndim=5
  ndim_tot=3
  kmax=uEMEPNvertout
  scale=1.0
  CDFtype=Real4
  !  dimSizes(1)=uEMEP%Nsec_poll
  !  dimNames(1)='sector'
  dimSizes(1)=2*lf_src(1)%dist+1
  dimNames(1)='x_dist'
  dimSizes(2)=2*lf_src(1)%dist+1
  dimNames(2)='y_dist'

  uEMEP%DOMAIN(1) = max(RUNDOMAIN(1),uEMEP%DOMAIN(1))
  uEMEP%DOMAIN(2) = min(RUNDOMAIN(2),uEMEP%DOMAIN(2))
  uEMEP%DOMAIN(3) = max(RUNDOMAIN(3),uEMEP%DOMAIN(3))
  uEMEP%DOMAIN(4) = min(RUNDOMAIN(4),uEMEP%DOMAIN(4))

  dimSizes(3)=min(GIMAX,uEMEP%DOMAIN(2)-uEMEP%DOMAIN(1)+1)
  dimSizes(4)=min(GJMAX,uEMEP%DOMAIN(4)-uEMEP%DOMAIN(3)+1)

  dimSizes_tot(1)=min(GIMAX,uEMEP%DOMAIN(2)-uEMEP%DOMAIN(1)+1)
  dimSizes_tot(2)=min(GJMAX,uEMEP%DOMAIN(4)-uEMEP%DOMAIN(3)+1)

  select case(projection)
  case('Stereographic')
     dimNames(3)='i'
     dimNames(4)='j'      
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'      
  case('lon lat')
     dimNames(3)='lon'
     dimNames(4)='lat'
     dimNames_tot(1)='lon'
     dimNames_tot(2)='lat'      
  case('Rotated_Spherical')
     dimNames(3)='i'
     dimNames(4)='j'      
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'      
  case('lambert')
     dimNames(3)='i'
     dimNames(4)='j'      
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'      
  case default
     dimNames(3)='i'
     dimNames(4)='j'      
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'      
  end select

  dimSizes(5)=kmax
  dimNames(5)='klevel'
  dimSizes_tot(3)=kmax
  dimNames_tot(3)='klevel'
  def1%class='uEMEP' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=1.0      !not used
  def1%name=trim(varName)
  def1%unit=''
  def2=def1
  def2%unit='ug/m3'
  chunksizes=1
  !chunksizes(1)=dimSizes(1) !slower!!
  !chunksizes(2)=dimSizes(2) !slower!!
  chunksizes(3)=MAXLIMAX
  chunksizes(4)=MAXLJMAX
  chunksizes(5)=dimSizes(5)
  chunksizes_tot=1
  chunksizes_tot(1)=MAXLIMAX
  chunksizes_tot(2)=MAXLJMAX
  chunksizes_tot(3)=dimSizes_tot(3)

  allocate(tmp_out(Ndiv2_coarse,LIMAX,LJMAX)) !NB; assumes KMAX=1 TEMPORARY
 
  iou_ix = iotyp2ix(iotyp)

  !first loop only create all variables before writing into them (faster for NetCDF)
  do iter=1,2
     if(iter==1 .and. .not. first_call(iotyp))cycle
     if(iter==1)then
        !only create all variables before writing into them
        overwrite=.true.
        create_var_only=.true.
     else
        overwrite=.false.
        create_var_only=.false.
     endif
     pollwritten = .false.
     do isrc = 1, Nsources
        isec=lf_src(isrc)%sector
        ipoll=lf_src(isrc)%poll
        if(.not. pollwritten(ipoll))then !one pollutant may be used for several sources
           def2%name=trim(lf_src(isrc)%emis)
           scale=1.0/av_fac_full
           call Out_netCDF(iotyp,def2,ndim_tot,kmax,lf_src_tot(1,1,KMAX_MID-uEMEPNvertout+1,ipoll,iou_ix),scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
                fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes_tot)                
           pollwritten(ipoll) = .true.
        endif
        if(iter==2)then
           do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
              do j=1,ljmax
                 do i=1,limax
                    invtot=1.0/(lf_src_tot(i,j,k,ipoll,iou_ix)+1.E-20)
                    n1=0
                    do n=lf_src(isrc)%start, lf_src(isrc)%end
                       n1=n1+1
                       tmp_out(n1,i,j) = lf_src_acc(n,i,j,k,iou_ix)*invtot
                    enddo
                 enddo
              enddo
           enddo
        endif
        write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%emis)//'_sec',isec,'_fraction_'//trim(lf_src(isrc)%type)
        if(isec==0) write(def1%name,"(A,I2.2,A)")trim(lf_src(isrc)%emis)//'_fraction_'//trim(lf_src(isrc)%type)
        scale=1.0
        call Out_netCDF(iotyp,def1,ndim,kmax,tmp_out,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=overwrite,create_var_only=create_var_only,chunksizes=chunksizes)
     enddo
  enddo
  deallocate(tmp_out)
  
  do ipoll=1,Npoll
     do k = KMAX_MID-uEMEPNvertout+1,KMAX_MID
        do j=1,ljmax
           do i=1,limax
              lf_src_tot(i,j,k,ipoll,iou_ix) = 0.0
           enddo
        enddo
     enddo
  enddo
  
  !reset the cumulative counters
  if(iotyp==IOU_HOUR)then
     av_fac_hour=0
  else  if(iotyp==IOU_DAY)then
     av_fac_day=0.0
  else  if(iotyp==IOU_MON)then
     av_fac_month=0.0
  else  if(iotyp==IOU_YEAR)then
     av_fac_full=0.0
  endif

  first_call(iotyp)=.false.

  call Add_2timing(NTIMING-2,tim_after,tim_before,"uEMEP: output")

! CALL MPI_BARRIER(MPI_COMM_CALC, I)
 
!stop
end subroutine out_uEMEP

subroutine av_uEMEP(dt,End_of_Day)
  real, intent(in)    :: dt                   ! time-step used in integrations
  logical, intent(in) :: End_of_Day           ! e.g. 6am for EMEP sites
  real :: xtot
  integer ::i,j,k,n,dx,dy,ix,iix,ipoll,isec_poll1, iou_ix, isrc
  integer ::isec_poll
  logical :: pollwritten(Npoll_uemep_max)
  
  call Code_timer(tim_before)
  if(.not. uEMEP%HOUR.and.&
     .not. uEMEP%DAY .and.&
     .not. uEMEP%MONTH .and.&
     .not. uEMEP%YEAR       )return

  !do the averaging
  pollwritten = .false.
  do iou_ix = 1, Niou_ix
     do isrc=1,Nsources
        ipoll = lf_src(isrc)%poll
        do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
           do j=1,ljmax
              do i=1,limax
                 xtot=0.0
                 do iix=1,lf_src(isrc)%Nsplit
                    ix=lf_src(isrc)%ix(iix)
                    xtot=xtot+(xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix))/ATWAIR&
                         *roa(i,j,k,1)*1.E9 !for ug/m3
                    !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
                 end do
                 if(.not. pollwritten(ipoll))then !one pollutant may be used for several sources
                    lf_src_tot(i,j,k,ipoll,iou_ix) = lf_src_tot(i,j,k,ipoll,iou_ix) + xtot
                 endif
                 do n=lf_src(isrc)%start, lf_src(isrc)%end
                    lf_src_acc(n,i,j,k,iou_ix)=lf_src_acc(n,i,j,k,iou_ix)+xtot*lf(n,i,j,k)
                 end do
              enddo
           enddo
        enddo
        pollwritten(ipoll) = .true.
      enddo
  enddo
  
  av_fac_hour=av_fac_hour+1.0
  av_fac_day=av_fac_day+1.0
  av_fac_month=av_fac_month+1.0
  av_fac_full=av_fac_full+1.0

  call Add_2timing(NTIMING-8,tim_after,tim_before,"uEMEP: averaging")

end subroutine av_uEMEP

subroutine uemep_adv_x(fluxx,i,j,k)
  real, intent(in)::fluxx(NSPEC_ADV,-1:LIMAX+1)
  integer, intent(in)::i,j,k
  real ::x,xn,xx,f_in,inv_tot
  integer ::n,iix,ix,dx,dy,isrc

  call Code_timer(tim_before)
  do isrc=1,Nsources
     xn=0.0
     x=0.0
     xx=0.0
     !positive x or xx means incoming, negative means outgoing
     do iix=1,lf_src(isrc)%Nsplit
        ix=lf_src(isrc)%ix(iix)
        xn=xn+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
        x=x-xm2(i,j)*fluxx(ix,i)*lf_src(isrc)%mw(iix)!flux through "East" face (Right)
        xx=xx+xm2(i,j)*fluxx(ix,i-1)*lf_src(isrc)%mw(iix)!flux through "West" face (Left)
     end do
     !NB: here xn already includes the fluxes. Remove them!
     xn=xn-xx-x
     xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux 
     f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
     inv_tot=1.0/(xn+f_in+1.e-40)!incoming dilutes

     x =max(0.0,x)*inv_tot!factor due to flux through "East" face (Right)
     xx=max(0.0,xx)*inv_tot!factor due to flux through "West" face (Left)
     xn = xn * inv_tot
     !often either x or xx is zero
     if(lf_src(isrc)%type=='coarse')then
        if(x>1.E-20)then
           do n = lf_src(isrc)%start, lf_src(isrc)%end
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,i+1)*x
           enddo
           if(xx>1.E-20)then
              do n = lf_src(isrc)%start, lf_src(isrc)%end
                 lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n,i-1)*xx
              enddo
           endif
        else if (xx>1.E-20)then
           do n = lf_src(isrc)%start, lf_src(isrc)%end
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,i-1)*xx
           enddo
        endif
     else if(lf_src(isrc)%type=='relative')then
        if(x>1.E-20)then
           n = lf_src(isrc)%start
           do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
              lf(n,i,j,k) = lf(n,i,j,k)*xn ! when dx=-uEMEP%dist there are no local fractions to transport
              n=n+1
              do dx=-lf_src(isrc)%dist+1,lf_src(isrc)%dist
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n-1,i+1)*x
                 n=n+1
              enddo
           enddo

           if(xx>1.E-20)then
              n = lf_src(isrc)%start
              do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
                 do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                    lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+1,i-1)*xx   
                    n=n+1
                 enddo                
                 n=n+1! when dx=lf_src(isrc)%dist there are no local fractions to transport
              enddo
           endif
        else if (xx>1.E-20)then
           n = lf_src(isrc)%start
           do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist
              do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n+1,i-1)*xx 
                 n=n+1
              enddo
              lf(n,i,j,k) = lf(n,i,j,k)*xn! when dx=lf_src(isrc)%dist there are no local fractions to transport
              n=n+1
           enddo
        else
          !nothing to do if no incoming fluxes
        endif
      else
        if(me==0)write(*,*)'LF type not recognized)'
        stop
     endif
  enddo

  call Add_2timing(NTIMING-7,tim_after,tim_before,"uEMEP: adv_x")

end subroutine uemep_adv_x

subroutine uemep_adv_y(fluxy,i,j,k)
  real, intent(in)::fluxy(NSPEC_ADV,-1:LJMAX+1)
  integer, intent(in)::i,j,k
  real ::x,xn,xx,f_in,inv_tot
  integer ::n,iix,ix,dx,dy,isrc

  call Code_timer(tim_before)
  do isrc=1,Nsources
     xn=0.0
     x=0.0
     xx=0.0
     !positive x or xx means incoming, negative means outgoing
     do iix=1,lf_src(isrc)%Nsplit
        ix=lf_src(isrc)%ix(iix)
        xn=xn+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
        x=x-xm2(i,j)*fluxy(ix,j)*lf_src(isrc)%mw(iix)!flux through "North" face (Up)
        xx=xx+xm2(i,j)*fluxy(ix,j-1)*lf_src(isrc)%mw(iix)!flux through "South" face (Bottom)
     end do
     !NB: here xn already includes the fluxes. Remove them!
     xn=xn-xx-x
     xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux 
     f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
     inv_tot=1.0/(xn+f_in+1.e-40)!incoming dilutes

     x =max(0.0,x)*inv_tot!factor due to flux through "East" face (Right)
     xx=max(0.0,xx)*inv_tot!factor due to flux through "West" face (Left)
     xn = xn * inv_tot
     if(lf_src(isrc)%type=='coarse')then
        !often either x or xx is zero
        if(x>1.E-20)then
           do n = lf_src(isrc)%start, lf_src(isrc)%end
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,j+1)*x
           enddo
           if(xx>1.E-20)then
              do n = lf_src(isrc)%start, lf_src(isrc)%end
                 lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n,j-1)*xx
              enddo
           endif
        else if (xx>1.E-20)then
           do n = lf_src(isrc)%start, lf_src(isrc)%end
              lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n,j-1)*xx
           enddo
        endif
     else if(lf_src(isrc)%type=='relative')then
        if(x>1.E-20)then
           n = lf_src(isrc)%start
           dy = -lf_src(isrc)%dist
           do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
              lf(n,i,j,k) = lf(n,i,j,k)*xn
              n=n+1
           enddo
           do dy=-lf_src(isrc)%dist+1,lf_src(isrc)%dist
              do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n-Ndiv_rel,j+1)*x
                 n=n+1
              enddo
           enddo
           if(xx>1.E-20)then
              n = lf_src(isrc)%start
              do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
                 do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                    lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_1d(n+Ndiv_rel,j-1)*xx
                    n=n+1
                 enddo
              enddo
           endif
        else if (xx>1.E-20)then
           n = lf_src(isrc)%start
           do dy=-lf_src(isrc)%dist,lf_src(isrc)%dist-1
              do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
                 lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_1d(n+Ndiv_coarse,j-1)*xx 
                 n=n+1
              enddo
           enddo
           dy=lf_src(isrc)%dist
           do dx=-lf_src(isrc)%dist,lf_src(isrc)%dist
              lf(n,i,j,k) = lf(n,i,j,k)*xn
              n=n+1
           enddo
        else
           !nothing to do if no incoming fluxes          
        endif

     else
        if(me==0)write(*,*)'LF type not recognized)'
        stop
     endif

  enddo
  call Add_2timing(NTIMING-6,tim_after,tim_before,"uEMEP: adv_y")

end subroutine uemep_adv_y

subroutine uemep_adv_k(fluxk,i,j)
    real, intent(in)::fluxk(NSPEC_ADV,KMAX_MID)
    integer, intent(in)::i,j
    real ::x,xn,xx,f_in,inv_tot
    integer ::n,k,iix,ix,dx,dy,isrc
    real loc_frac_src_km1(LF_SRC_TOTSIZE,KMAX_MID-uEMEP%Nvert:KMAX_MID-1)

    call Code_timer(tim_before)
    !need to be careful to always use non-updated values on the RHS
    do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID-1
       loc_frac_src_km1(:,k)=lf(:,i,j,k)
    enddo
    loc_frac_src_km1(:,KMAX_MID-uEMEP%Nvert)=0.0!Assume zero local fractions coming from above

    do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID!k is increasing-> can use k+1 to access non-updated value
       do isrc=1,Nsources
          xn=0.0
          x=0.0
          xx=0.0
          !positive x or xx means incoming, negative means outgoing
          do iix=1,lf_src(isrc)%Nsplit
             ix=lf_src(isrc)%ix(iix)
             xn=xn+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
             if(k<KMAX_MID)x=x-dhs1i(k+1)*fluxk(ix,k+1)*lf_src(isrc)%mw(iix)
             xx=xx+dhs1i(k+1)*fluxk(ix,k)*lf_src(isrc)%mw(iix)
          end do
          !NB: here xn already includes the fluxes. Remove them!
          xn=xn-xx-x
          xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux 
          f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
          inv_tot=1.0/(xn+f_in+1.e-40)!incoming dilutes

          x =max(0.0,x)*inv_tot!factor due to flux bottom facethrough
          xx=max(0.0,xx)*inv_tot!factor due to flux through top face
          xn = xn * inv_tot
          !often either x or xx is zero
          if(lf_src(isrc)%type=='coarse' .or. lf_src(isrc)%type=='relative')then
             if(x>1.E-20 .and. k<KMAX_MID)then
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   lf(n,i,j,k) = lf(n,i,j,k)*xn +lf(n,i,j,k+1)*x! loc_frac(isec_poll,dx,dy,i,j,k+1)*x
                enddo
                if(xx>1.E-20 .and. k>KMAX_MID-uEMEP%Nvert+1)then
                   do n = lf_src(isrc)%start, lf_src(isrc)%end
                      lf(n,i,j,k) = lf(n,i,j,k) + loc_frac_src_km1(n,k-1)*xx
                   enddo
                endif
             else if (xx>1.E-20 .and. k>KMAX_MID-uEMEP%Nvert+1)then
                do n = lf_src(isrc)%start, lf_src(isrc)%end
                   lf(n,i,j,k) = lf(n,i,j,k)*xn + loc_frac_src_km1(n,k-1)*xx
                enddo
             else
                !nothing to do if no incoming fluxes
             endif
          else
             if(me==0)write(*,*)'LF type not recognized)'
             stop
          endif
       enddo
       
    end do
    
    call Add_2timing(NTIMING-5,tim_after,tim_before,"uEMEP: adv_k")
  end subroutine uemep_adv_k
 
  subroutine uemep_diff(i,j,ds3,ds4,ndiff)
    
    implicit none
    interface 
       subroutine vertdiffn(xn_k,NSPEC,Nij,KMIN_in,SigmaKz,ds3,ds4,ndiff)
         real,intent(inout) :: xn_k(NSPEC,0:*)!dummy
         real,intent(in)::  SigmaKz(*)!dummy
         real,intent(in)::  ds3(*),ds4(*)!dummy
         integer,intent(in)::  NSPEC,ndiff,Nij,KMIN_in
       end subroutine vertdiffn
    end interface

    real, intent(in) :: ds3(2:KMAX_MID),ds4(2:KMAX_MID)
    integer, intent(in) :: i,j,ndiff
    real :: xn_k(LF_SRC_TOTSIZE + Npoll,KMAX_MID),x
    integer ::isec_poll1,isrc
    integer ::k,n,ix,iix,dx,dy
    !how far diffusion should take place above uEMEP%Nvert. 
    ! KUP = 2 gives less than 0.001 differences in locfrac, except sometimes over sea, because
    !ship emission are higher up and need to come down to diminish locfrac
    integer, parameter :: KUP = 2

    call Code_timer(tim_before)
    xn_k = 0.0
    do k = 1,KMAX_MID
       do isrc=1,Nsources
          x=0.0
          do iix=1,lf_src(isrc)%Nsplit
             ix=lf_src(isrc)%ix(iix)
             !assumes mixing ratios units, but weight by mass
             x=x+xn_adv(ix,i,j,k)*lf_src(isrc)%mw(iix)
          end do
          if(k>KMAX_MID-uEMEP%Nvert)then ! lf zero above
             do n=lf_src(isrc)%start, lf_src(isrc)%end
                xn_k(n,k)=x*lf(n,i,j,k)
             enddo
          endif
          xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k) = x
      enddo
    enddo

    call vertdiffn(xn_k,LF_SRC_TOTSIZE+Npoll,1,KMAX_MID-uEMEP%Nvert-KUP,EtaKz(i,j,1,1),ds3,ds4,ndiff)
 
    do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
     do isrc=1,Nsources
          x =  1.0/(xn_k(LF_SRC_TOTSIZE+lf_src(isrc)%poll,k)+1.E-30)
          do n=lf_src(isrc)%start, lf_src(isrc)%end
             lf(n,i,j,k) = xn_k(n,k)*x                               
          enddo
       enddo
    end do
    call Add_2timing(NTIMING-4,tim_after,tim_before,"uEMEP: diffusion")

end subroutine uemep_diff

subroutine uEMEP_emis(indate)
!include emission contributions to local fractions

!NB: should replace most of the stuff and use gridrcemis instead!

  implicit none
  type(date), intent(in) :: indate  ! Gives year..seconds
  integer :: i, j, k, n       ! coordinates, loop variables
  integer :: icc, ncc         ! No. of countries in grid.
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILE

  ! Save daytime value between calls, initialise to zero
  integer, save, dimension(MAXNLAND) ::  daytime(1:MAXNLAND) = 0  !  0=night, 1=day
  integer, save, dimension(MAXNLAND) ::  localhour(1:MAXNLAND) = 1  ! 1-24 local hour in the different countries
  integer                         ::  hourloc      !  local hour 
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions
  real ::  tfac    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s       ! source term (emis) before splitting
  integer :: iland, iland_timefac, iland_timefac_hour  ! country codes, and codes for timefac 
  integer :: hour_iland
  integer ::icc_uemep, iqrc, itot
  integer, save :: wday , wday_loc ! wday = day of the week 1-7
  integer ::ix,iy,iix, iiix,dx, dy, isec_poll, iisec_poll, isec_poll1, isrc
  real::dt_uemep, xtot, emis_uemep(KMAX_MID,NEMIS_FILE,NSECTORS),emis_tot(KMAX_MID,NEMIS_FILE)
  real :: lon
  integer :: jmin,jmax,imin,imax,n0
  
  call Code_timer(tim_before)
 
  do j = lj0,lj1
    do i = li0,li1
      ix=(gi0+i-2)/((GIMAX+Ndiv_coarse-1)/Ndiv_coarse)+1 !i coordinate in coarse domain
      iy=(gj0+j-2)/((GJMAX+Ndiv_coarse-1)/Ndiv_coarse)+1 !j coordinate in coarse domain
      do isrc=1,Nsources   
         iem = lf_src(isrc)%iem
         isec = lf_src(isrc)%sector
         do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
            if(lf_emis_tot(i,j,k,lf_src(isrc)%poll)<1.E-20)cycle
            xtot=0.0
            do iix=1,lf_src(isrc)%Nsplit
               iiix=lf_src(isrc)%ix(iix)
               xtot=xtot+(xn_adv(iiix,i,j,k)*lf_src(isrc)%mw(iix))*(dA(k)+dB(k)*ps(i,j,1))/ATWAIR/GRAV
            end do
            if(lf_src(isrc)%type=='coarse')then
               n0 = lf_src(isrc)%start+ix-1+(iy-1)*Ndiv_coarse !emissions included as local           
               lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot+lf_emis(i,j,k,isrc))/(xtot+lf_emis_tot(i,j,k,lf_src(isrc)%poll)+1.e-20)
            else if(lf_src(isrc)%type=='relative')then
               n0 = lf_src(isrc)%start + (lf_src(isrc)%Npos - 1)/2 !"middle" point is dx=0 dy=0
               lf(n0,i,j,k)=(lf(n0,i,j,k)*xtot+lf_emis(i,j,k,isrc))/(xtot+lf_emis_tot(i,j,k,lf_src(isrc)%poll)+1.e-20)
            else
               if(me==0)write(*,*)'LF type not recognized)'
               stop
            endif
            do n = lf_src(isrc)%start, lf_src(isrc)%end
               if(n==n0)cycle  !counted above               
               lf(n,i,j,k)=(lf(n,i,j,k)*xtot)/(xtot+lf_emis_tot(i,j,k,lf_src(isrc)%poll)+1.e-20)!fractions are diluted
             enddo
         enddo
      enddo
    end do ! i
 end do ! j

  call Add_2timing(NTIMING-3,tim_after,tim_before,"uEMEP: emissions")

end subroutine uEMEP_emis

subroutine add_lf_emis(s,i,j,iem,isec,iland)
  real, intent(in) :: s
  integer, intent(in) :: i,j,iem,isec,iland
  integer :: isrc, k, ipoll

  ipoll = iem2ipoll(iem)
  if(ipoll>0)then
     do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
        lf_emis_tot(i,j,k,ipoll) = lf_emis_tot(i,j,k,ipoll) + s * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) * dt_advec!total for each pollutant
     enddo
   endif

  do isrc = 1, Nsources
     if(lf_src(isrc)%iem /= iem) cycle
     if(lf_src(isrc)%sector /= isec .and. lf_src(isrc)%sector /= 0) cycle
     if(lf_src(isrc)%iqrc>0)then
        !single pollutant, part of emitted group of pollutant        
        lf_emis(i,j,k,isrc) = lf_emis(i,j,k,isrc) + s * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) &
             * emisfrac(lf_src(isrc)%iqrc,sec2split_map(isec),iland) *dt_advec       
     else
        do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
           lf_emis(i,j,k,isrc) = lf_emis(i,j,k,isrc) + s * emis_kprofile(KMAX_BND-k,sec2hfac_map(isec)) * dt_advec
        enddo
     endif
  enddo
  
end subroutine add_lf_emis

end module uEMEP_mod
