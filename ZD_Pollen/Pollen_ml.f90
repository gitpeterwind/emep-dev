module Pollen_ml
!-----------------------------------------------------------------------!
! Birch pollen emission calculation based on
! M. Sofiev et al. 2006, doi:10.1007/s00484-006-0027-x
!
! Pollen emission based upon meteorology paparameters, and heatsum.
! Pollen particles are assumed of 22 um diameter and 800 kg/m3 density.
!-----------------------------------------------------------------------!
use Pollen_const_ml
use PhysicalConstants_ml, only: AVOG
use Biogenics_ml,         only: EMIS_BioNat, EmisNat
use CheckStop_ml,         only: CheckStop,CheckNC
use ChemSpecs,            only: NSPEC_SHL, species_adv
use Chemfields_ml,        only: xn_adv    ! emep model concs.
use DerivedFields_ml,     only: f_2d,d_2d ! D2D houtly (debug) output
use GasParticleCoeffs_ml, only: AERO_SIZE,CDDEP_BIRCH,CDDEP_OLIVE,CDDEP_GRASS,CDDEP_RWEED
use Functions_ml,         only: heaviside
use GridValues_ml ,       only: glon, glat, debug_proc, debug_li, debug_lj
use Landuse_ml,           only: LandCover
use LocalVariables_ml,    only: Grid
use MetFields_ml,         only: surface_precip, ws_10m ,rh2m,t2_nwp,&
                                foundws10_met,foundprecip,pr,u_ref,z_bnd,z_mid
use MicroMet_ml,          only: Wind_at_h
use ModelConstants_ml,    only: AERO, KMAX_MID, nstep, FORECAST, &
                                METSTEP, MasterProc, IOU_INST, RUNDOMAIN, &
                                dt=>dt_advec, DEBUG=>DEBUG_POLLEN
use MPI_Groups_ml,        only: MPI_INTEGER,MPI_LOGICAL,MPI_COMM_CALC,&
                                MasterPE,IERROR
use Nest_ml,              only: outdate,FORECAST_NDUMP,out_DOMAIN,&
                                template_read_IC=>template_read_3D,&
                                template_write_IC=>template_write
use NetCDF_ml,            only: ReadField_CDF,Out_netCDF,GetCDF_modelgrid,&
                                ReadTimeCDF,CDFtype=>Real4
use netcdf,               only: nf90_close
use Par_ml,               only: limax, ljmax, LIMAX, LJMAX, me
use PhysicalConstants_ml, only: PI
use OwnDataTypes_ml,      only: Deriv
use Setup_1dfields_ml,    only: rcemis
use SmallUtils_ml,        only: find_index
use SubMet_ml,            only: Sub
use TimeDate_ml,          only: current_date,daynumber,date,day_of_year
use TimeDate_ExtraUtil_ml,only: date2string,compare_date,date2nctime
use Io_ml,                only: IO_NML, PrintLog
use mpi,                  only: MPI_COMM_WORLD,MPI_DOUBLE_PRECISION,&
                                MPI_IN_PLACE,MPI_MIN,MPI_MAX!,MPI_ALLREDUCE
! openMPI has no explicit interface for MPI_ALLREDUCE
!-------------------------------------
implicit none
private
public:: pollen_flux,pollen_dump,pollen_read

!** 1) Public (saved) Variables from module:

real, public,save , allocatable, dimension(:,:,:)::&
  heatsum,      & ! heatsum, needs to be remembered for forecast
  AreaPOLL,     & ! emission of pollen
  R               ! pollen released so far

! pollen arrays indexing, order must match with POLLEN_GROUP: birch,olive,grass
integer, save, dimension(POLLEN_NUM) :: &
  inat=-1,iadv=-1,itot=-1

integer, parameter :: &
  max_string_length=200 ! large enough for paths to be set on Pollen_config namelist
!-------------------------------------
! Variables read from NetCDF Files
! pollen_flux   NetCDF file    NetCDF var
!   pollen_frac birch_frac_nc  birch_frac
!   birch_corr  birch_corr_nc  scale_factor (year specific version)
!   birch_corr  birch_data_nc  cross
!   birch_h_c   birch_data_nc  h_c
!   pollen_frac olive_data_nc  olive_frac
!   olive_h_c   olive_data_nc  olive_th
!   pollen_frac grass_field_nc grass_frac
!   grass_start grass_time_nc  grass_start
!   grass_end   grass_time_nc  grass_end
!   grass_end   grass_time_nc  grass_length+grass_start if grass_end not found
! Variables write/read from dump/restart Files
! pollen_dump/pollen_read   pollen ype (SPC)    NetCDF var
!   xn_adv(i,:,:,:)         BIRCH,OLIVE,GRASS   SPC
!   N_TOT(i)-R(:,:,i)       BIRCH,OLIVE,GRASS   SPC//'_rest'
!   heatsum(:,:,i)          BIRCH,OLIVE         SPC//'_heatsum'
character(len=max_string_length), save :: &
  birch_frac_nc ='birch_frac.nc',         &
  birch_data_nc ='pollen_data.nc',        &
  birch_corr_nc ='birch_factor_YYYY.nc',  &
  olive_data_nc ='olive_YYYY.nc',    &
  grass_field_nc='grass_frac.nc',         &
  grass_time_nc ='grass_time.nc',         &
  grass_mode    ='linear',      & ! 'linear' (old) | 'gamma' (new)
  template_read ='POLLEN_IN.nc',& ! dump/restart input
  template_write='POLLEN_OUT.nc'  ! dump/restart output

type(date), parameter :: &
  date_first_birch=date(-1,3,1,0,0),date_last_birch=date(-1,8,1,0,0),&
  date_first_olive=date(-1,1,1,0,0),date_last_olive=date_last_birch
type(date), save :: & ! will be updated when grass_time.nc is read
  date_first_grass=date(-1,2,7,0,0),date_last_grass=date(-1,9,23,0,0)

integer, save :: day_first=1,day_last=366

real, save, allocatable, dimension(:,:) :: &
  grass_start,grass_end ! Stard/End day of grass, read in

logical, parameter :: DEBUG_NC=.false.

!-------------------------------------------------------------------------!
contains
!-------------------------------------------------------------------------!
subroutine Config_Pollen()
  integer :: ios,g,n
  logical, save :: first_call = .true.
  NAMELIST /Pollen_config/&
    birch_frac_nc,birch_data_nc,birch_corr_nc,&
    olive_data_nc,grass_field_nc,grass_time_nc,grass_mode,&
    template_read,template_write

  if(.not.first_call)return
  first_call = .false.

  ! check consistency between Pollen_const_ml and species_adv definitions
  call pollen_check()

  template_read =template_read_IC   ! by default read/write
  template_write=template_write_IC  ! to Next IC/restart file
  rewind(IO_NML)
  read(IO_NML,NML=Pollen_config,iostat=ios)
  call CheckStop(ios,"NML=Pollen_config")
  if(debug.and.MasterProc)then
    write(*,*) "NAMELIST IS "
    write(*,NML=Pollen_config)
  end if

  do g=1,POLLEN_NUM
    inat(g) = find_index(POLLEN_GROUP(g),EMIS_BioNat(:))
    iadv(g) = find_index(POLLEN_GROUP(g),species_adv(:)%name)
    itot(g) = iadv(g)+NSPEC_SHL
    call CheckStop(inat(g)<0,"EMIS_BioNat misses: "//POLLEN_GROUP(g))
    call CheckStop(iadv(g)<0,"species_adv misses: "//POLLEN_GROUP(g))

    ! override gravitational setling params
    select case(POLLEN_GROUP(g))
      case(BIRCH);n=AERO_SIZE(CDDEP_BIRCH)
      case(OLIVE);n=AERO_SIZE(CDDEP_OLIVE)
      case(RWEED);n=AERO_SIZE(CDDEP_RWEED)
      case(GRASS);n=AERO_SIZE(CDDEP_GRASS)
    end select
    AERO%DpgV(n)=D_POLL(g)*1e-6  ! um to m
  ! AERO%sigma(n)=0.01
    AERO%PMdens(n)=POLL_DENS*1e-3 ! g/m3 to kg/m3
  end do
  if(MasterProc)write(*,"(A,10(' adv#',I3,'=',A,1X,es10.3,:))") &
    "Pollen: ",(iadv(g),POLLEN_GROUP(g),grain_wt(g),g=1,POLLEN_NUM)

  allocate(heatsum(LIMAX,LJMAX,POLLEN_NUM-1),&     ! Grass does not need heatsum
           R(LIMAX,LJMAX,POLLEN_NUM))
  heatsum(:,:,:) = 0.00
  R(:,:,:) = 0.0
end subroutine Config_Pollen
!-------------------------------------------------------------------------!
function checkdates(nday,spc,update) result(ok)
  integer, intent(in)           :: nday
  character(len=*), intent(in)  :: spc
  logical, intent(in), optional :: update
  logical :: ok
  integer, save :: day_first_b=1,day_last_b=366,&
                   day_first_o=1,day_last_o=366,&
                   day_first_g=1,day_last_g=366
  logical, save :: first_call = .true.
  if(present(update))first_call=first_call.and.update
  if(first_call)then
    day_first_b=day_of_year(current_date%year,date_first_birch%month,date_first_birch%day)
    day_first_o=day_of_year(current_date%year,date_first_olive%month,date_first_olive%day)
    day_first_g=day_of_year(current_date%year,date_first_grass%month,date_first_grass%day)
    day_first_g=day_first_g-uncert_grass_day
    day_first=min(day_first_b,day_first_o,day_first_g)

    day_last_b=day_of_year(current_date%year,date_last_birch%month,date_last_birch%day)
    day_last_o=day_of_year(current_date%year,date_last_olive%month,date_last_olive%day)
    day_last_g=day_of_year(current_date%year,date_last_grass%month,date_last_grass%day)
    day_last_g=day_last_g+uncert_grass_day
    day_last=max(day_last_b,day_last_o,day_last_g)
    first_call = .false.
  end if
  select case(spc)
  case("POLLEN","P","p","pollen")
    ok=(day_first  <=nday).and.(nday<=day_last  )
  case(BIRCH,"B","b","birch")
    ok=(day_first_b<=nday).and.(nday<=day_last_b)
  case(OLIVE,"O","o","olive")
    ok=(day_first_o<=nday).and.(nday<=day_last_o)
  case(GRASS,"G","g","grass")
    ok=(day_first_g<=nday).and.(nday<=day_last_g)
  case default
    call CheckStop("Unknown pollen type: "//trim(spc))
  end select
end function checkdates
!-------------------------------------------------------------------------!
subroutine pollen_flux(i,j,debug_flag)
  implicit none
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug_flag

  real, parameter :: &
    kgm2h(POLLEN_NUM)=grain_wt(:)*1e-6*3600, & ! EmisNat: grains/m2/s --> kg/m2/h
    Z10  = 10.0,              & ! 10m height
    UnDef=  0.0

  logical, save :: first_call=.true.
  logical :: debug_ij=.false.,found=.false.,pollen_out(POLLEN_NUM)=.false.

  integer :: ii, jj, nlu, ilu, lu, info, g, n
  real :: scale,lim_birch,lim_olive,rcpoll,&
       n2m(POLLEN_NUM),u10,prec,relhum,dfirst_g,dlast_g

  real, save, allocatable, dimension(:,:,:) :: &
    pollen_frac   ! fraction of pollen (birch/olive/grass), read in
  real, save, allocatable, dimension(:,:) :: &
    birch_h_c,  & ! temperature treshold birch, read in
    birch_corr, & ! correction field for p.emission, read in
    olive_h_c,  & ! temperature treshold olive, read in
    olive_dH,   & ! flowering period [degree days]
    h_day         ! temperature summed over a day

  integer,save,allocatable, dimension(:,:) :: &
    p_day

  ! Read in the different fields
  if(first_call) then
    if(.not.checkdates(daynumber,"pollen"))return
    call Config_Pollen()
    if(FORECAST) call pollen_read()

    allocate(p_day(LIMAX,LJMAX),AreaPOLL(LIMAX,LJMAX,POLLEN_NUM),h_day(LIMAX,LJMAX))
    p_day(:,:) =current_date%day
    AreaPOLL(:,:,:)=0.0
    h_day(:,:)=0.0

    allocate(pollen_frac(LIMAX,LJMAX,POLLEN_NUM))
    allocate(birch_h_c(LIMAX,LJMAX),birch_corr(LIMAX,LJMAX))
    allocate(olive_h_c(LIMAX,LJMAX),olive_dH(LIMAX,LJMAX))
    allocate(grass_start(LIMAX,LJMAX),grass_end(LIMAX,LJMAX))
    call ReadField_CDF(birch_frac_nc,'birch_frac',pollen_frac(:,:,1),1, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(birch_data_nc,'h_c',birch_h_c,1, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
! birch: cross_corr for specific year
    birch_corr_nc=date2string(birch_corr_nc,current_date,debug=DEBUG_NC.and.MasterProc)
    call ReadField_CDF(birch_corr_nc,'scale_factor',birch_corr,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
! birch: cross_corr default
    if(.not.found)&
      call ReadField_CDF(birch_data_nc,'cross',birch_corr,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
! olive
    olive_data_nc=date2string(olive_data_nc,current_date,debug=DEBUG_NC.and.MasterProc)
    call ReadField_CDF(olive_data_nc,'olive_frac',pollen_frac(:,:,2),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(olive_data_nc,'olive_th',olive_h_c,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(olive_data_nc,'olive_len',olive_dH,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
    if(.not.found) olive_dH(:,:) = dH_d_olive
! grass
    call ReadField_CDF(grass_field_nc,'grass_frac',pollen_frac(:,:,4),1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(grass_time_nc,'grass_start',grass_start,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(grass_time_nc,'grass_end',grass_end,1, &
        interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef,found=found)
    if(.not.found)then
      call ReadField_CDF(grass_time_nc,'grass_length',grass_end,1, &
        interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
      where(grass_end/=UnDef.and.grass_start/=UnDef)&
        grass_end=grass_end+grass_start
    end if

    ! reduce birch from 60 degrees north:
    forall(ii=1:limax,jj=1:ljmax,glat(ii,jj)>=60.0) &
      pollen_frac(ii,jj,1)=pollen_frac(ii,jj,1)&
                       *max(0.3,1.0-0.005*(glat(ii,jj)-60.0))

    ! olive fraction [%] --> [1/1]
    forall(ii=1:limax,jj=1:ljmax,pollen_frac(ii,jj,2)/=UnDef) &
      pollen_frac(ii,jj,2)=pollen_frac(ii,jj,2)/100.0

    ! start/end of grass pollen season
    dfirst_g=minval(grass_start,MASK=(grass_start/=UnDef))
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,dfirst_g,1,MPI_DOUBLE_PRECISION,&
                       MPI_MIN,MPI_COMM_WORLD,INFO)
    date_first_grass=date(-1,1,floor(dfirst_g),0,0)

    dlast_g=maxval(grass_end   ,MASK=(grass_start/=UnDef))
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,dlast_g ,1,MPI_DOUBLE_PRECISION,&
                       MPI_MAX,MPI_COMM_WORLD,INFO)
    date_last_grass=date(-1,1,ceiling(dlast_g),0,0)

    ! Set D2D/USET output
    call write_uset()

    ! start of pollen forecast
    if(checkdates(daynumber,"pollen",update=.true.).and.MasterProc)&
      write(*,*) "POLLEN setup ",date2string("YYYY-MM-DD hh:mm",current_date)
    first_call = .false.
  end if !first_call
  debug_ij=all([DEBUG.or.debug_flag,debug_proc,i==debug_li,j==debug_lj])

  pollen_out(1)=.not.checkdates(daynumber,BIRCH)&
    .or.any([pollen_frac(i,j,1),birch_h_c(i,j)]==UnDef)

  pollen_out(2)=.not.checkdates(daynumber,OLIVE)&
    .or.any([pollen_frac(i,j,2),olive_h_c(i,j)]==UnDef)

  pollen_out(4)=.not.checkdates(daynumber,GRASS)&
    .or.(daynumber<(grass_start(i,j)-uncert_grass_day))&
    .or.(daynumber>(grass_end  (i,j)+uncert_grass_day))&
    .or.any([pollen_frac(i,j,3),grass_start(i,j),grass_end(i,j)]==UnDef)


  EmisNat(inat(:),i,j)      = 0.0
  rcemis (itot(:),KMAX_MID) = 0.0
  AreaPOLL(i,j,:)           = 0.0
  if(all(pollen_out))then
    call write_uset()
    return
  end if

 !  Heatsum calculations
 !  Sums up the temperatures that day for each timestep (20 min).
 !  The heatsum are then divided by timesteps, and added using the
 !  function heatsum_calc.
  if(nstep==1) then ! Heatsum calculation only each meteorological timestep
    if(p_day(i,j)==current_date%day) then
      h_day(i,j) = h_day(i,j) + t2_nwp(i,j,1)
      if(current_date%hour==21)then
        if(.not.pollen_out(1)) &
          heatsum(i,j,1)=heatsum(i,j,1)&
            +heatsum_calc((h_day(i,j)/(24/METSTEP)),T_cutoff_birch)
        if(.not.pollen_out(2)) &
          heatsum(i,j,2)=heatsum(i,j,2)&
            +heatsum_calc((h_day(i,j)/(24/METSTEP)),T_cutoff_olive)
      end if
    else
      p_day(i,j) = current_date%day
      h_day(i,j) = t2_nwp(i,j,1)
    end if
    ! End of heatsum calculations
  end if ! heatsum calc each timestep

  ! if heatsum is over heatsum threhold for the grid cell, the pollen
  ! emission can start calculating
  lim_birch = (1-prob_in_birch)*birch_h_c(i,j)
  lim_olive = (1-prob_in_olive)*olive_h_c(i,j)

  !Need to set the meteorological fields if not defined in nwp.
  relhum = 0.0
  u10    = 0.0
  prec   = 0.0

  if(foundws10_met) u10=ws_10m(i,j,1)
  nlu = LandCover(i,j)%ncodes
  do ilu=1,nlu
    lu = LandCover(i,j)%codes(ilu)
    relhum = Sub(lu)%rh
    !wind
    if(.not.foundws10_met .and. Sub(lu)%z0 > 1.0)&
      u10=Wind_at_h(Grid%u_ref, Grid%z_ref, Z10,&
      Sub(lu)%d, Sub(lu)%z0, Sub(lu)%invL)
    if(relhum>=0.0 .and. u10>=0.0) exit
  end do

  ! precipitation
  if(foundprecip) then
    prec=sum(pr(i,j,:))*60  !mm/h
  else
    prec=surface_precip(i,j)
  end if

! Birch specific emission inhibitors
  pollen_out(1)=pollen_out(1)         & ! out season
    .or.(R(i,j,1)>N_TOT(1)*birch_corr(i,j)) & ! out of pollen
    .or.(relhum>RH_HIGH)              & ! too humid
    .or.(prec>prec_max)               & ! too rainy
    .or.(heatsum(i,j,1)<lim_birch)    & ! too cold
    .or.(t2_nwp(i,j,1)<T_cutoff_birch)& ! too windy
    .or.(heatsum(i,j,1)-birch_h_c(i,j)> dH_d_birch)

! Olive specific emission inhibitors
  pollen_out(2)=pollen_out(2)         & ! out season
    .or.(R(i,j,2)>N_TOT(2))           & ! out of pollen
    .or.(relhum>RH_HIGH)              & ! too humid
    .or.(prec>prec_max)               & ! too rainy
    .or.(heatsum(i,j,2)<lim_olive)    & ! too cold
    .or.(t2_nwp(i,j,1)<T_cutoff_olive)& ! too windy
    .or.(heatsum(i,j,2)-olive_h_c(i,j)> olive_dH(i,j))

! Grass specific emission inhibitors
  pollen_out(4)=pollen_out(4)         & ! out season
    .or.(R(i,j,4)>N_TOT(4))           & ! out of pollen
    .or.(relhum>RH_HIGH)              & ! too humid
    .or.(prec>prec_max)                 ! too rainy

  ! grains to molecular weight: grains/m2/s --> mol/cm3/s
  n2m(:) = 1e-6/Grid%DeltaZ           & ! 1/dZ [1/cm3]
    *grain_wt*AVOG/species_adv(iadv(:))%molwt
!  write(*,*) "n2m_diff ", n2m(:)

!------------------------
! Emission rates: Birch,Olive,Grass
!------------------------
  do g=1,POLLEN_NUM
    if(pollen_out(g))cycle
    ! scale factor meteorological conditions
    select case(POLLEN_GROUP(g))
      case(BIRCH);scale=scale_factor(BIRCH)
      case(OLIVE);scale=scale_factor(OLIVE)
      case(RWEED);cycle
      case(GRASS);scale=scale_factor(GRASS//trim(grass_mode))
    endselect
    R(i,j,g)=R(i,j,g)+N_TOT(g)*scale*dt       ! pollen grains released so far
    rcpoll=N_TOT(g)*scale*pollen_frac(i,j,g)  ! pollen grains production [grains/m2/s]
    EmisNat(inat(g),i,j)      = rcpoll*kgm2h(g) ! [kg/m2/h]
    rcemis (itot(g),KMAX_MID) = rcpoll*n2m(g) ! [mol/cm3/s]
    AreaPOLL(i,j,g)           = rcpoll*3600   ! [grains/m2/h]

    if(debug_ij) write(*,'(a,3(1x,I3),3(1x,es10.3))')&
      POLLEN_GROUP(g),me,i,j,rcemis(itot(g),KMAX_MID),R(i,j,g),AreaPOLL(i,j,g)
  end do
  call write_uset()
contains
!------------------------
! scale factor for different meteorological conditions
!------------------------
function scale_factor(spc) result(scale)
  character(len=*), intent(in)  :: spc
  real :: scale, dHsec
  integer :: g

  scale=f_wind(u10,Grid%wstar)           & ! wind dependence
       *f_cond(relhum,RH_LOW,RH_HIGH)    & ! rh dependence
       *f_cond(prec,PREC_MIN,PREC_MAX)     ! precipitation dependence
  select case(spc)
  case(BIRCH)
    g=1
    scale = scale*birch_corr(i,j) &
      *(t2_nwp(i,j,1)-T_cutoff_birch)/dH_birch &
      *f_in(heatsum(i,j,g),birch_h_c(i,j),PROB_IN_birch) &  ! prob. flowering start
      *f_out(R(i,j,g),N_TOT(g),PROB_OUT_birch)              ! prob. flowering end
  case(OLIVE)
    g=2
    dHsec = olive_dH(i,j)*24*3600   ! Flowering period [degree seconds] olive
    scale = scale &
      *(t2_nwp(i,j,1)-T_cutoff_olive)/dHsec &
      *f_in(heatsum(i,j,g),olive_h_c(i,j),PROB_IN_olive) &  ! prob. flowering start
      *f_out(R(i,j,g),N_TOT(g),PROB_OUT_olive)              ! prob. flowering end
  case(GRASS//'linear')   ! emission mass assuming linear release
    g=4
    scale = scale &
      *f_fade_in (real(daynumber)/grass_start(i,j), &
                  uncert_grass_day/grass_start(i,j))& ! fade-in
      *f_fade_out(real(daynumber)/grass_end(i,j),   &
                  uncert_grass_day/grass_start(i,j))& ! fade-out
      *f_fade_out(R(i,j,g)/N_TOT(g),&
                  1.0-uncert_tot_grass_poll)          ! total-pollen fade-out
    ! Full-emission rate is total pollen divided by the total duration of the season
    scale = scale/(grass_end(i,j)-daynumber+uncert_grass_day)&
                 /(grass_end(i,j)-grass_start(i,j))/86400.0
  case(GRASS//'gamma')    ! assume the modified "taily" Gamma distribution of the season
    g=4
    scale = scale &
      *f_gamma_w_tails((real(daynumber)-grass_start(i,j))  & ! days since season start
                       /(grass_end(i,j)-grass_start(i,j)), & ! season length
                        dt/86400.0                         & ! timestep in days
                       /(grass_end(i,j)-grass_start(i,j)))   ! season length
    scale = scale/dt ! emited fration over dt to emission rate
  case default
    call CheckStop("Unknown pollen type: "//trim(spc))
  end select
end function scale_factor
!------------------------
! Write D2D/USET output at the end of every hour
!------------------------
subroutine write_uset()
  ! indexes for USET/D2D debug output
  integer, save, dimension(POLLEN_NUM) :: &
    n2d_heatsum=-1,n2d_left=-1,n2d_emiss=-1
  integer :: g,n
  if(first_call)then
    do n=1,size(f_2d)
      if(f_2d(n)%class/='USET')cycle
      g=find_index(f_2d(n)%txt,POLLEN_GROUP)
      if(g<0)cycle
      select case(f_2d(n)%subclass)
      case("heatsum")
        call CheckStop(g,[1,size(heatsum,DIM=3)],&
          "USET: '"//trim(f_2d(n)%subclass)//"' out of bounds!")
        f_2d(n)%unit="degree day"
        f_2d(n)%scale=1.0
        n2d_heatsum(g)=n
      case("pollen_left")
        call CheckStop(g,[1,size(R,DIM=3)],&
          "USET: '"//trim(f_2d(n)%subclass)//"' out of bounds!")
        f_2d(n)%unit="1"
        f_2d(n)%scale=1.0
        n2d_left(g)=n
      case("pollen_emiss")
        call CheckStop(g,[1,size(AreaPOLL,DIM=3)],&
          "USET: '"//trim(f_2d(n)%subclass)//"' out of bounds!")
        f_2d(n)%unit="grains/m2/h"
        f_2d(n)%scale=1.0
        n2d_emiss(g)=n
      case default
        cycle
      endselect
    end do
    return
  end if
  do g=1,POLLEN_NUM
    n=n2d_heatsum(g)  ! heatsum
    if(n>0) d_2d(n,i,j,IOU_INST)=heatsum(i,j,g)
    n=n2d_left(g)     ! pollen_left
    if(n>0) d_2d(n,i,j,IOU_INST)=1.0-R(i,j,g)/N_TOT(g)
    n=n2d_emiss(g)    ! pollen_emiss
    if(n>0) d_2d(n,i,j,IOU_INST)=AreaPOLL(i,j,g)
  end do
end subroutine write_uset
end subroutine pollen_flux
function heatsum_calc(t2,T_cutoff) result(ff)
! The temperature needs to be over the cutoff temperature
  real, intent(in) :: t2,T_cutoff
  real             :: ff
! ff = (t2-T_cutoff)*heaviside(t2-T_cutoff)
  ff = DIM(t2,T_cutoff) ! same as max(t2-T_cutoff,0.0)
endfunction heatsum_calc
function f_gamma_w_tails(relTime,relDt) result(ff)
! Returns the pollen prepared for release assuming the modified "taily" Gamma distribution
! of the season. Tails are the reason for many parameters: have to describe the main peak
! via gamma-type distribution, and both elevated tails via add-on corrections.
! formula: rate(x)=exp(-a_1/beta)* sum(scale_i * a_i^power_i), i=1:3
!          where a_i = max(x-timesRel_i,0)
  real, intent(in) :: relTime,relDt ! normalised time, normalised timestep
  real             :: ff
! Fitting parameters
  integer, parameter:: nTerms= 3
  real, parameter::&
    beta = 0.155, &! [main season shape,correction 1,correction 2]
!   timesRel(nTerms)=[ 0.164,-0.26,-0.7 ],&
    timesRel(nTerms)=[ 0.164,-0.3 ,-0.6 ],&
    scales  (nTerms)=[13.1  ,12.6 , 0.25],&
    powers  (nTerms)=[ 1.16 , 2.8 , 1.3 ]
! Local variable
  real :: a1

  ff = 0.0
  a1 = DIM(relTime,timesRel(1)) ! same as max(relTime-timesRel(1),0.0)
  if(a1>beta*10.0)return ! too far from the season peak, as decided by timesRel(1)
  ! Be careful: the rise of the function can be quite steep
  a1=exp(-a1/beta)
! ff=0.5*relDt*(rate(a1,relTime)+rate(a1,relTime+relDt)) ! trapezoidal integration
! WRONG integration to match SILAM code
  ff=0.5*relDt*(rate(a1,relTime)+rate(a1*a1,relTime+relDt))
contains
real function rate(a,x)
  real, intent(in) :: a,x
  rate=a*sum(scales*(x-timesRel)**powers,MASK=(x>timesRel))
end function rate
end function f_gamma_w_tails
function f_wind(u,wstar) result(ff)
! Pollen emission increases with wind speeds up to 5 m/s (u_sat).
! This term cannot be higher than 1.5 (u_max).
  real, intent(in) :: u,wstar
  real             :: ff
  real, parameter  :: u_max=1.5,u_sat=5.0
  ff = u_max - exp(-(u+wstar)/u_sat)
end function f_wind
function f_cond(x,x_min,x_max) result(ff)
! This function is used for both humitidy and rain, as too much
! humidity and rain stops pollen release
  real, intent(in) :: x,x_min,x_max
  real             :: ff
  if(x>x_max) then
    ff=0.0
  elseif(x<x_min)then
    ff=1.0
  else
    ff=(x_max-x)/(x_max-x_min)
  end if
end function f_cond
function f_in(h,h_c,prob) result(ff)
! takes in account the uncertainity that all the trees start
! to release pollen at the same time
  real, intent(in) :: H,H_c,prob
  real             :: ff
  if(h<(1.0-prob)*h_c) then
    ff=0.0
  elseif(h>(1.0+prob)*h_c)then
    ff=1.0
  else
    ff=(h-(1.0-prob)*h_c)/(2.0*prob*h_c)
  end if
end function f_in
function f_out(R,N_tot,prob) result(ff)
! takes in account the uncertainity that all the trees stop
! releasing pollen at the same time
  real, intent(in) :: R,N_tot,prob
  real             :: ff
  if(R<(1.0-prob)*N_tot) then
    ff=1.0
  elseif (R>(1.0+prob)*N_tot) then
    ff=0.0
  else
    ff=1.0-(R-(1.0-prob)*N_tot)/(2.0*prob*N_tot)
  end if
end function f_out
function f_fade_in(value_rel,uncert_rel_) result(ff)
! Computes the linear fade-in function. It is 0 at vaule_rel = 1.-uncert_rel
! Grows to 0.5 at vaule_rel = 1 and grows to 1 at value_rel = 1.+uncert_rel
  real, intent(in) :: value_rel,uncert_rel_
  real             :: uncert_rel,ff

  uncert_rel = max(uncert_rel_,1e-5) ! avoid zero uncertainty
  ff = min(1.0,max(0.0,(value_rel-1.0+uncert_rel)/(2.0*uncert_rel)))
end function f_fade_in
function f_fade_out(value_rel,uncert_rel_) result(ff)
! Computes the linear fade-in function. It is 0 at vaule_rel = 1.-uncert_rel
! Grows to 0.5 at vaule_rel = 1 and grows to 1 at value_rel = 1.+uncert_rel
  real, intent(in) :: value_rel,uncert_rel_
  real             :: uncert_rel,ff

  uncert_rel = max(uncert_rel_,1e-5) ! avoid zero uncertainty
  ff = min(1.0,max(0.0,(1.0+uncert_rel-value_rel)/(2.0*uncert_rel)))
end function f_fade_out
!-------------------------------------------------------------------------!
integer function getRecord(fileName,findDate,fatal)
  character(len=*), intent(in) :: fileName
  type(date), intent(in) :: findDate
  logical, intent(in) :: fatal

  if(MasterProc) getRecord=startRecord()
  call MPI_BCAST(getRecord,1,MPI_INTEGER,MasterPE,MPI_COMM_CALC,IERROR)
  if(getRecord<1)&  ! switch off Pollen_ml when restart file is corrupted
    call MPI_BCAST(USE_POLLEN,1,MPI_LOGICAL,MasterPE,MPI_COMM_CALC,IERROR)
contains
function startRecord() result(nstart)
  integer :: nstart

  logical :: fexist
  integer :: nread
  real(kind=8), parameter :: sec2day=1e0/(24.0*3600.0)
  real, dimension(366), save :: fdays=-1
  real :: ncday(0:1)

  ! ensure file exists
  inquire(file=filename,exist=fexist)
  if(.not.fexist)then
    call CheckStop(fatal,"File not found: "//trim(filename))
    if(MasterProc) write(*,*)&
      "Warning Pollen dump file not found: "//trim(filename)
    call PrintLog("WARNING: Pollen_ml cold start",MasterProc)
    nstart=-1
    return
  end if

  ! ensure file has records
  nread=-1                                ! read all
  fdays(:)=-1.0                           ! times records
  call ReadTimeCDF(filename,fdays,nread)  ! in fname
  if(nread<1)then
    call CheckStop(fatal,"Corrupted file: "//trim(filename))
    if(MasterProc) write(*,*)&
      "Warning Pollen dump file corrupted: "//trim(filename)
    call PrintLog("WARNING: Pollen_ml forced OFF",MasterProc)
    USE_POLLEN=.false.
    nstart=-1
    return
  end if

  ! look for current_date in fdays (time var read from filename)
  call date2nctime(findDate,ncday(1))
  ncday(0)=ncday(1)-dt*sec2day    ! to avoid rounding errors
  ncday(1)=ncday(1)+dt*sec2day    ! look for records in +/- 1dt

  nstart=MINLOC(fdays(:nread),DIM=1,&
    MASK=(fdays(:nread)>=ncday(0)).and.(fdays(:nread)<ncday(1)))

  ! check if we got a match
  if(nstart<1)then  ! ifort compiler needs option: -assume noold_maxminloc
    call CheckStop(fatal,&
      "No records for"//date2string(" YYYY-MM-DD hh:mm ",findDate)//&
      "found in "//trim(filename))
    if(MasterProc) write(*,*)&
      "No records for"//date2string(" YYYY-MM-DD hh:mm ",findDate)//&
      "found in "//trim(filename)
    call PrintLog("WARNING: Pollen_ml cold start",MasterProc)
    nstart=-1
    return
  end if
end function startRecord
end function getRecord
!-------------------------------------------------------------------------!
subroutine pollen_read()
! Read pollen for forecast restart file: heatsum and pollenrest fields
  character(len=*), parameter :: dfmt="('Read: ',A20,' [',L1,'].')"
  character(len=20) :: spc
  logical, save :: first_call = .true.
  logical :: found
  integer :: nstart,g
  character(len=len(template_read)) :: filename
  real,allocatable, dimension(:,:,:) :: data ! Data arrays

  if(.not.checkdates(daynumber,"pollen")) return
  if(.not.first_call) return
  first_call=.false.

  call Config_Pollen()
  filename=date2string(template_read,current_date)
  nstart=getRecord(filename,current_date,.not.FORECAST)
  if(nstart<1) return
  if(MasterProc)&
    write(*,"(3(A,1X),I0)") "Read Pollen dump",trim(filename),"record",nstart

  allocate(data(LIMAX,LJMAX,KMAX_MID))
  do g=1,size(POLLEN_GROUP)
    spc=trim(POLLEN_GROUP(g))
!------------------------
! pollen adv (not written by Nest_ml)
!------------------------
    call GetCDF_modelgrid(trim(spc),filename,data,&
          1,KMAX_MID,nstart,1,needed=.not.FORECAST,found=found)
    if(DEBUG.and.MasterProc) write(*,dfmt)spc,found
    if(found)xn_adv(iadv(g),:,:,:)=data(:,:,:)
!------------------------
! pollen_rest
!------------------------
    spc=trim(POLLEN_GROUP(g))//'_rest'
    call GetCDF_modelgrid(trim(spc),filename,data(:,:,1),&
        1,1,nstart,1,needed=.not.FORECAST,found=found)
    if(DEBUG.and.MasterProc) write(*,dfmt)spc,found
    if(found)R(:,:,g)=N_TOT(g)-data(:,:,1)
!------------------------
! heatsum
!------------------------
    if(g>size(heatsum,DIM=3))cycle
    spc=trim(POLLEN_GROUP(g))//'_heatsum'
    call GetCDF_modelgrid(trim(spc),filename,heatsum(:,:,g),&
        1,1,nstart,1,needed=.not.FORECAST,found=found)
    if(DEBUG.and.MasterProc) write(*,dfmt)spc,found
  end do
  deallocate(data)
end subroutine pollen_read
!-------------------------------------------------------------------------!
subroutine pollen_dump()
! Write pollen for forecast restart file: heatsum and pollenrest fields
  character(len=len(template_write)) :: filename
  character(len=*), parameter :: dfmt="('Write: ',A20,' [',A,'].')"
  character(len=20) :: spc
  logical     :: fexist,create_var
  logical,save:: overwrite=.true. ! append to existing files if false
  integer     :: ncfileID,i,g
  type(Deriv) :: def1
  real,allocatable, dimension(:,:,:) :: data ! Data arrays

  if(.not.checkdates(daynumber,"pollen")) return
  if(.not.compare_date(FORECAST_NDUMP,current_date,&
                       outdate(:FORECAST_NDUMP),wildcard=-1))return
  call CheckStop(allocated(R).NEQV.allocated(heatsum),&
    "Pollen: rest/heatsum allocated error")
  if(.not.(allocated(R).and.allocated(heatsum))) return

  filename = date2string(template_write,current_date)
  if(MasterProc) write(*,*) "Write Pollen dump ",trim(filename)
  inquire(file=filename,exist=fexist)

  def1%avg=.false.                  ! not used
  def1%index=0                      ! not used
  def1%scale=1.0                    ! not used
  def1%iotype=''                    ! not used

  allocate(data(LIMAX,LJMAX,KMAX_MID))
  ncfileID=-1 ! must be <0 as initial value
  do i=1,2                          ! do first one loop to define the fields,
    create_var=(i==1)               ! without writing them (for performance purposes),
    if(all([create_var,fexist,.not.overwrite,&  ! if the file does not exists already
       template_write/=template_write_IC]))cycle! and pollen has its own restart file
    overwrite=overwrite.and.create_var

    do g=1,size(POLLEN_GROUP)
      spc=trim(POLLEN_GROUP(g))
!------------------------
! pollen adv (not written by Nest_ml)
!------------------------
      def1%class='Advected'           ! written
      def1%unit='mix_ratio'           ! written
      def1%name=trim(spc)             ! written
      if(DEBUG.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      if(.not.create_var)data=xn_adv(iadv(g),:,:,:)
      call Out_netCDF(IOU_INST,def1,3,KMAX_MID,data,1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,overwrite=overwrite,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
      overwrite=.false.
!------------------------
! pollen_rest
!------------------------
      def1%class='pollen_out'         ! written
      def1%unit='pollengrains'        ! written
      def1%name=trim(spc)//'_rest'    ! written
      if(DEBUG.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      data(:,:,1)=N_TOT(g)-R(:,:,g)
      call Out_netCDF(IOU_INST,def1,2,1,data(:,:,1),1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
!------------------------
! heatsum
!------------------------
      if(g>size(heatsum,DIM=3))cycle
      def1%class='pollen_out'         ! written
      def1%unit='degreedays'          ! written
      def1%name=trim(spc)//'_heatsum' ! written
      if(DEBUG.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      call Out_netCDF(IOU_INST,def1,2,1,heatsum(:,:,g),1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
    end do
  end do
  if(MasterProc)&
    call CheckNC(nf90_close(ncFileID),"close:"//trim(filename))
  ! ensure time record can be found, fatal error if not
  i=getRecord(filename,current_date,.true.) ! MPI_BCAST inside
  if(MasterProc.and.DEBUG) &
    write(*,"(3(A,1X),I0)") "Found Pollen dump",trim(filename),"record",i
  deallocate(data)
end subroutine pollen_dump
end module Pollen_ml
