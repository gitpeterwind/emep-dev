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
use ChemSpecs,            only: NSPEC_SHL, species
use Chemfields_ml,        only: xn_adv    ! emep model concs.
use Functions_ml,         only: heaviside
use GridValues_ml ,       only: glon, glat, debug_proc, debug_li, debug_lj
use Landuse_ml,           only: LandCover
use LocalVariables_ml,    only: Grid
use MetFields_ml,         only: surface_precip, ws_10m ,rh2m,t2_nwp,&
                                foundws10_met,foundprecip,pr,u_ref,z_bnd,z_mid
use MicroMet_ml,          only: Wind_at_h
use ModelConstants_ml,    only: KMAX_MID, nstep, FORECAST, &
                                METSTEP, MasterProc, IOU_INST, RUNDOMAIN, &
                                dt=>dt_advec, DEBUG=>DEBUG_POLLEN
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
use Io_ml,                only: IO_NML
use mpi,                  only: MPI_COMM_WORLD,MPI_DOUBLE_PRECISION,&
                                MPI_IN_PLACE,MPI_MIN,MPI_MAX!,MPI_ALLREDUCE
! openMPI has no explicit interface for MPI_ALLREDUCE
!-------------------------------------
implicit none
private
public:: pollen_flux,pollen_dump,pollen_read

!** 1) Public (saved) Variables from module:

real, public,save , allocatable, dimension(:,:,:)::&
  AreaPOLL,     & ! emission of pollen 
  heatsum,      & ! heatsum, needs to be remembered for forecast
  R,            & ! Summed pollen release
  Pollen_rest,  & ! what is available pollen pr m2, fr forecast
  Pollen_left     ! amount of pollen left in catkins, relative amount... 0:sr 1:end 

integer, save :: &
  inat_BIRCH=-1,inat_OLIVE=-1,inat_GRASS=-1,&
  itot_BIRCH=-1,itot_OLIVE=-1,itot_GRASS=-1,&
  iadv_BIRCH=-1,iadv_OLIVE=-1,iadv_GRASS=-1

integer, parameter :: &
  max_string_length=200 ! large enough for paths to be set on Pollen_config namelist
character(len=max_string_length), save ::  &  
  birch_frac_nc ="birch_frac.nc",&
  birch_data_nc ="pollen_data.nc",&
  birch_corr_nc ="birch_factor_YYYY.nc",&
  olive_data_nc ="maccoliven_data.nc",&
  grass_field_nc="grass_frac.nc",&
  grass_time_nc ="grass_time.nc",&
  template_read ="POLLEN_IN.nc",&
  template_write="POLLEN_OUT.nc"

type(date), parameter :: &
  date_first_birch=date(-1,3,1,0,0),date_last_birch=date(-1,8,1,0,0),&
  date_first_olive=date(-1,1,1,0,0),date_last_olive=date_last_birch
type(date), save :: & ! will be updated when grass_time.nc is read
  date_first_grass=date(-1,2,7,0,0),date_last_grass=date(-1,9,23,0,0)

integer, save :: day_first=1,day_last=366

real, save, allocatable, dimension(:,:) :: &
  grass_start, & ! Stard day of grass, read in
  grass_len      ! length of gras blooming, read in

logical, parameter :: DEBUG_NC=.false.

!-------------------------------------------------------------------------!
contains
!-------------------------------------------------------------------------!
subroutine Config_Pollen()
  integer :: ios
  logical, save :: first_call = .true.
  NAMELIST /Pollen_config/&
    birch_frac_nc,birch_data_nc,birch_corr_nc,&
    olive_data_nc,grass_field_nc,grass_time_nc,&
    template_read,template_write

  if(.not.first_call)return
  first_call = .false.

  template_read =template_read_IC   ! by default read/write
  template_write=template_write_IC  ! to Next IC/restart file
  rewind(IO_NML)
  read(IO_NML,NML=Pollen_config,iostat=ios)
  call CheckStop(ios,"NML=Pollen_config")  
  if(debug.and.MasterProc)then
    write(*,*) "NAMELIST IS "
    write(*,NML=Pollen_config)
  endif

  inat_BIRCH = find_index(BIRCH,EMIS_BioNat(:))
  inat_OLIVE = find_index(OLIVE,EMIS_BioNat(:))
  inat_GRASS = find_index(GRASS,EMIS_BioNat(:))
  itot_BIRCH = find_index(BIRCH,species(:)%name)
  itot_OLIVE = find_index(OLIVE,species(:)%name)
  itot_GRASS = find_index(GRASS,species(:)%name)
  iadv_BIRCH = itot_BIRCH-NSPEC_SHL
  iadv_OLIVE = itot_OLIVE-NSPEC_SHL
  iadv_GRASS = itot_GRASS-NSPEC_SHL
  if(MasterProc)write(*,"(A,3(' adv#',I3,'=',A))") "Pollen: ",&
    iadv_BIRCH,BIRCH,iadv_OLIVE,OLIVE,iadv_GRASS,GRASS

  allocate(heatsum(LIMAX,LJMAX,2),&     ! Grass does not need heatsum
           Pollen_rest(LIMAX,LJMAX,3),&
           R(LIMAX,LJMAX,3))

  heatsum(:,:,:) = 0.00 
  Pollen_rest(:,:,1) = N_TOT_birch
  Pollen_rest(:,:,2) = N_TOT_olive
  Pollen_rest(:,:,3) = N_TOT_grass
  R(:,:,:) = 0.0 
endsubroutine Config_Pollen
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
  endif
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
  endselect
endfunction checkdates
!-------------------------------------------------------------------------!
subroutine pollen_flux(i,j,debug_flag) 
  implicit none
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug_flag

  real, parameter :: &
    Z10 = 10.0, &  ! 10m height
    UnDef= 0.0
  logical, save :: first_call=.true.
  logical :: &
    debug_ij=.false., out_birch=.false., out_olive=.false., out_grass=.false.
   
  integer :: ii, jj, nlu, ilu, lu, info
  real :: scale,lim_birch,lim_olive,invdz,rcpoll,&
       n2m_birch,n2m_olive,n2m_grass,u10,prec,relhum,&
       dfirst_g,dlast_g,moleccm3s_2_kgm2h

  real, save, allocatable, dimension(:,:) :: &
    birch_frac, & ! fraction of birch, read in
    birch_h_c,  & ! temperature treshold birch, read in
    corr,       & ! correction field for p.emission, read in
    olive_frac, & ! fraction of olive, read in
    olive_h_c,  & ! temperature treshold olive, read in
    grass_frac, & ! fraction of grass, read in
    h_day         ! temperature summed over a day

  integer,save,allocatable, dimension(:,:) :: &
    p_day
 
  ! Read in the different fields
  if(first_call) then 
    if(.not.checkdates(daynumber,"pollen"))return
    call Config_Pollen() 
    if(FORECAST) call pollen_read()

    allocate(Pollen_left(LIMAX,LJMAX,3),p_day(LIMAX,LJMAX))
    Pollen_left(:,:,:) = 1
    p_day(:,:) =current_date%day

    allocate(AreaPOLL(LIMAX,LJMAX,3),h_day(LIMAX,LJMAX))
    AreaPOLL(:,:,:)=0.0
    h_day(:,:)=0.0
      
    allocate(birch_frac(LIMAX,LJMAX),birch_h_c(LIMAX,LJMAX),&
             corr(LIMAX,LJMAX))
    allocate(olive_frac(LIMAX,LJMAX),olive_h_c(LIMAX,LJMAX))
    allocate(grass_frac(LIMAX,LJMAX),grass_start(LIMAX,LJMAX),&
             grass_len(LIMAX,LJMAX))
    call ReadField_CDF(birch_frac_nc,'birch_frac',birch_frac,2, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(birch_data_nc,'h_c',birch_h_c,2, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
! birch: cross_corr default
    call ReadField_CDF(birch_data_nc,'cross',corr,2, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
! birch: cross_corr for specific year
    birch_corr_nc=date2string(birch_corr_nc,current_date,debug=DEBUG_NC.and.MasterProc)
    call ReadField_CDF(birch_corr_nc,'scale_factor',corr,2, &
       interpol='conservative',needed=.false.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(olive_data_nc,'olive_frac',olive_frac,2, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(olive_data_nc,'olive_th',olive_h_c,2, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(grass_field_nc,'grass_frac',grass_frac,2, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(grass_time_nc,'grass_start',grass_start,2, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)
    call ReadField_CDF(grass_time_nc,'grass_length',grass_len,2, &
       interpol='conservative',needed=.true.,debug_flag=DEBUG_NC,UnDef=UnDef)

    ! reduce birch from 60 degrees north:
    forall(ii=1:limax,jj=1:ljmax,glat(ii,jj)>=60.0) &
      birch_frac(ii,jj)=birch_frac(ii,jj)&
                       *max(0.3,1.0-0.005*(glat(ii,jj)-60.0))

    ! start/end of grass pollen season
    dfirst_g=minval(grass_start(:limax,:ljmax),&
                    MASK=(grass_start(:limax,:ljmax)/=UnDef))
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,dfirst_g,1,MPI_DOUBLE_PRECISION,&
                       MPI_MIN,MPI_COMM_WORLD,INFO)
    date_first_grass=date(-1,1,floor(dfirst_g),0,0)
    
    dlast_g=maxval(grass_start(:limax,:ljmax)+grass_len(:limax,:ljmax),&
                   MASK=(grass_start(:limax,:ljmax)/=UnDef))
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,dlast_g ,1,MPI_DOUBLE_PRECISION,&
                       MPI_MAX,MPI_COMM_WORLD,INFO)
    date_last_grass=date(-1,1,ceiling(dlast_g),0,0)

    ! start of pollen forecast
    if(checkdates(daynumber,"pollen",update=.true.).and.MasterProc)&
      write(*,*) "POLLEN setup ",date2string("YYYY-MM-DD hh:mm",current_date)
    first_call = .false. 
  endif !first_call
  debug_ij=all([DEBUG.or.debug_flag,debug_proc,i==debug_li,j==debug_lj])
 
  out_birch=.not.checkdates(daynumber,BIRCH)&
    .or.any([birch_frac(i,j),birch_h_c(i,j)]==UnDef)

  out_olive=.not.checkdates(daynumber,OLIVE)&
    .or.any([olive_frac(i,j),olive_h_c(i,j)]==UnDef)

  out_grass=.not.checkdates(daynumber,GRASS)&
    .or.(daynumber<(grass_start(i,j)-uncert_grass_day))&
    .or.(daynumber>(grass_start(i,j)+uncert_grass_day+grass_len(i,j)))&
    .or.any([grass_frac(i,j),grass_start(i,j)]==UnDef)
  
  if(all([out_olive,out_birch,out_grass]))then
    EmisNat(inat_BIRCH,i,j)      = 0.0
    EmisNat(inat_OLIVE,i,j)      = 0.0
    EmisNat(inat_GRASS,i,j)      = 0.0
    rcemis (itot_BIRCH,KMAX_MID) = 0.0
    rcemis (itot_OLIVE,KMAX_MID) = 0.0
    rcemis (itot_GRASS,KMAX_MID) = 0.0
    return
  endif

 !  Heatsum calculations
 !  Sums up the temperatures that day for each timestep (20 min).
 !  The heatsum are then divided by timesteps, and added using the 
 !  function heatsum_calc. 
  if(nstep==1) then ! Heatsum calculation only each meteorological timestep
    if(p_day(i,j)==current_date%day) then
      h_day(i,j) = h_day(i,j) + t2_nwp(i,j,1)
      if(current_date%hour==21)then
        if(.not.out_birch) &
          heatsum(i,j,1)=heatsum(i,j,1)+heatsum_calc((h_day(i,j)/(24/METSTEP)),T_cutoff_birch)
        if(.not.out_olive) &
           heatsum(i,j,2)=heatsum(i,j,2)+heatsum_calc((h_day(i,j)/(24/METSTEP)),T_cutoff_olive)
      endif
    else 
      p_day(i,j) = current_date%day  
      h_day(i,j) = t2_nwp(i,j,1)
    endif
    ! End of heatsum calculations
  endif ! heatsum calc each timestep
  ! if heatsum is over heatsum threhold for the grid cell, the pollen
  ! emission can start calculating 
  lim_birch = (1-prob_in_birch)*birch_h_c(i,j)
  lim_olive = (1-prob_in_olive)*olive_h_c(i,j)
!  if (i == 15 .and. j == 23 ) write(*,*) "Heatsum ",heatsum(i,j,2),lim_olive,olive_h_c(i,j),i,j,out_birch,out_olive


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
  enddo

  ! precipitation
  if(foundprecip) then 
    prec=sum(pr(i,j,:))*60  !mm/h
  else
    prec=surface_precip(i,j)
  endif

  out_birch=out_birch&
    .or. (heatsum(i,j,1)<lim_birch) &
    .or.(t2_nwp(i,j,1)<T_cutoff_birch) &
    .or.(prec>prec_max)         &
    .or.(relhum > RH_HIGH)      &
    .or.(R(i,j,1)>N_TOT_birch*corr(i,j)) & 
    .or.(heatsum(i,j,1)-birch_h_c(i,j)> dH_d_birch)
  
  out_olive=out_olive&
    .or.(heatsum(i,j,2)<lim_olive)     &
    .or.(t2_nwp(i,j,1)<T_cutoff_olive) &
    .or.(prec>prec_max)         &
    .or.(relhum > RH_HIGH)      &
    .or.(R(i,j,2)>N_TOT_olive)  &
    .or.(heatsum(i,j,2)-olive_h_c(i,j)> dH_d_olive)
  
  out_grass=out_grass&
    .or.(R(i,j,3)>N_TOT_grass)  &  
    .or.(relhum > RH_HIGH)      &
    .or.(prec>prec_max) 

  ! grains to molecular weight: grains/m2/s --> mol/cm3/s 
  invdz     = 1e-6/Grid%DeltaZ         ! 1/dZ [1/cm3]
  n2m_birch = invdz*grain_wt*AVOG/species(itot_BIRCH)%molwt
  n2m_olive = invdz*grain_wt*AVOG/species(itot_OLIVE)%molwt
  n2m_grass = invdz*grain_wt*AVOG/species(itot_GRASS)%molwt
!  write(*,*) "n2m_diff ", n2m_birch, n2m_olive,n2m_grass

  ! EmisNat: molec/cm3/s --> kg/m2/h
  moleccm3s_2_kgm2h = Grid%DeltaZ * 1e6 * 3600.0  &! /cm3/s > /m2/hr
                     /AVOG * 1e-6  ! kg  after *MW

  rcpoll         =0.0
  AreaPOLL(i,j,:)=0.0
  scale          =0.0
!------------------------
! Birch
!------------------------
  if(out_birch)then
    EmisNat(inat_BIRCH,i,j)      = 0.0
    rcemis (itot_BIRCH,KMAX_MID) = 0.0
  else
    !------------------------
    ! scale factor for different meteorological conditions
    !------------------------
    scale=corr(i,j)*(t2_nwp(i,j,1)-T_cutoff_birch)/dH_birch &
      *f_wind(u10,Grid%wstar)              & ! wind dependece
      *f_in(heatsum(i,j,1),birch_h_c(i,j),PROB_IN_birch) & ! probability for flowering to start
      *f_out(R(i,j,1),N_TOT_birch,PROB_OUT_birch)        & ! probability for flowering to end
      *f_cond(relhum,RH_LOW,RH_HIGH)       & ! rh dependence
      *f_cond(prec,PREC_MIN,PREC_MAX)        ! precipitation dependence

    R(i,j,1)=R(i,j,1)+N_TOT_birch*scale*dt              ! R is to check the remains of pollen available
    rcpoll       =N_TOT_birch*scale*birch_frac(i,j) ! The pollen grains production [grains/m2/s]
    AreaPOLL(i,j,1)=rcpoll*3600                 ! [grains/m2/h]

    !------------------------
    ! grains/m2/s --> mol/cm3/s
    !------------------------
    rcpoll = rcpoll * n2m_birch             ! [mol/cm3/s]
    Pollen_rest(i,j,1) = N_TOT_birch - R(i,j,1) ! pollen grains left after release [grains/m2]
    Pollen_left(i,j,1) = Pollen_rest(i,j,1)/N_TOT_birch ! idem [ ]
      
    if(debug_ij) write(6,'(a,3(1x,I3),3(1x,es10.3))')&
      "Result: Birch",me,i,j,rcpoll,R(i,j,1)/N_TOT_birch,AreaPOLL(i,j,1)

     EmisNat(inat_BIRCH,i,j)      = rcpoll*moleccm3s_2_kgm2h*species(itot_BIRCH)%molwt
     rcemis (itot_BIRCH,KMAX_MID) = rcpoll
  endif
!------------------------
! Olive
!------------------------
  if(out_olive)then 
    EmisNat(inat_OLIVE,i,j)      = 0.0
    rcemis (itot_OLIVE,KMAX_MID) = 0.0
  else
    !------------------------
    ! scale factor for different meteorological conditions
    !------------------------
    scale=(t2_nwp(i,j,1)-T_cutoff_olive)/dH_olive &
      *f_wind(u10,Grid%wstar)              & ! wind dependece
      *f_in(heatsum(i,j,2),olive_h_c(i,j),PROB_IN_olive) & ! probability for flowering to start
      *f_out(R(i,j,2),N_TOT_olive,PROB_OUT_olive)        & ! probability for flowering to end
      *f_cond(relhum,RH_LOW,RH_HIGH)       & ! rh dependence
      *f_cond(prec,PREC_MIN,PREC_MAX)        ! precipitation dependence

    R(i,j,2)=R(i,j,2)+N_TOT_olive*scale*dt              ! R is to check the remains of pollen available
    rcpoll       =N_TOT_olive*scale*olive_frac(i,j)/100! The pollen grains production [grains/m2/s]
    AreaPOLL(i,j,2)=rcpoll*3600                 ! [grains/m2/h]

    !------------------------
    ! grains/m2/s --> mol/cm3/s
    !------------------------
    rcpoll = rcpoll * n2m_olive             ! [mol/cm3/s]
    Pollen_rest(i,j,2) = N_TOT_olive - R(i,j,2) ! pollen grains left after release [grains/m2]
    Pollen_left(i,j,2) = Pollen_rest(i,j,2)/N_TOT_olive ! idem [ ]
      
    if(debug_ij) write(6,'(a,3(1x,I3),3(1x,es10.3))')&
      "Result: Olive",me,i,j,rcpoll,R(i,j,2)/N_TOT_olive,AreaPOLL(i,j,2)
     
     EmisNat(inat_OLIVE,i,j)      = rcpoll*moleccm3s_2_kgm2h*species(itot_OLIVE)%molwt
     rcemis (itot_OLIVE,KMAX_MID) = rcpoll
  endif
!------------------------
! Grass
!------------------------
  if(out_grass)then
    EmisNat(inat_GRASS,i,j)      = 0.0
    rcemis (itot_GRASS,KMAX_MID) = 0.0
  else
    !------------------------
    ! scale factor for different meteorological conditions
    !------------------------
    scale=f_wind(u10,Grid%wstar)           & ! wind dependece
      *f_cond(relhum,RH_LOW,RH_HIGH)       & ! rh dependence
      *(1. - prec/PREC_MAX)      &  ! precipitation dependence
      *f_fade_in(real(daynumber)/grass_start(i,j), &
                 uncert_grass_day/grass_start(i,j))  & ! CD-fade-in
      *f_fade_out(real(daynumber)/(grass_start(i,j)+grass_len(i,j)), &
                  uncert_grass_day/grass_start(i,j))&   !CD fade-out
      *f_fade_out(1.- (Pollen_rest(i,j,3)/N_TOT_grass),&
                  (1-uncert_tot_grass_poll))    !total-pollen fade-out
    
     ! Full-emission rate is total pollen divided by the total duration of the season
    scale = scale/(grass_start(i,j)+grass_len(i,j) -daynumber+uncert_grass_day)/(grass_len(i,j)*86400)

    R(i,j,3)=R(i,j,3)+N_TOT_grass*scale*dt              ! R is to check the remains of pollen available
    rcpoll       =N_TOT_grass*scale*grass_frac(i,j)     ! The pollen grains production [grains/m2/s]
    AreaPOLL(i,j,3)=rcpoll*3600                 ! [grains/m2/h]

    !------------------------
    ! grains/m2/s --> mol/cm3/s
    !------------------------
    rcpoll = rcpoll * n2m_grass             ! [mol/cm3/s]
    Pollen_rest(i,j,3) = N_TOT_grass - R(i,j,3) ! pollen grains left after release [grains/m2]
    Pollen_left(i,j,3) = Pollen_rest(i,j,3)/N_TOT_grass ! idem [ ]

    if(debug_ij) write(6,'(a,3(1x,I3),3(1x,es10.3))')&
      "Result: Grass",me,i,j,rcpoll,R(i,j,3)/N_TOT_grass,AreaPOLL(i,j,3)
     
     EmisNat(inat_GRASS,i,j)      = rcpoll*moleccm3s_2_kgm2h*species(itot_GRASS)%molwt
     rcemis (itot_GRASS,KMAX_MID) = rcpoll
  endif
endsubroutine pollen_flux
!-------------------------------------------------------------------------! 
subroutine pollen_read()
! Read pollen for forecast restart file: heatsum and pollenrest fields
  character(len=*), parameter :: dfmt="('Read: ',A20,' [',L1,'].')"
  character(len=20) :: spc
  logical, save :: first_call = .true.
  logical :: fexist,found
  integer :: nread,nstart,g,iadv
  real(kind=8), parameter :: sec2day=1e0/(24.0*3600.0)
  real, dimension(366), save :: fdays=-1
  real :: ncday(0:1),n_tot
  character(len=len(template_read)) :: filename
  real,allocatable, dimension(:,:,:) :: data ! Data arrays

  if(.not.checkdates(daynumber,"pollen")) return
  if(.not.first_call) return
  first_call=.false.

  call Config_Pollen()
  filename=date2string(template_read,current_date)
  inquire(file=filename,exist=fexist)
  if(FORECAST.and..not.fexist)then
    if(MasterProc) write(*,*)"Warning Polen dump file not found: "//trim(filename)
    return
  endif
  call CheckStop(.not.fexist,"Did not find file: "//trim(filename))

  nread=-1                                ! read all
  fdays(:)=-1.0                           ! times records
  call ReadTimeCDF(filename,fdays,nread)  ! in fname
  call date2nctime(current_date,ncday(1))
  ncday(0)=ncday(1)-dt*sec2day    ! to avoid rounding errors
  ncday(1)=ncday(1)+dt*sec2day    ! look for records in +/- 1dt

  nstart=MINLOC(fdays(:nread),DIM=1,&
    MASK=(fdays(:nread)>=ncday(0)).and.(fdays(:nread)<ncday(1)))

  if(.not.FORECAST)&
    call CheckStop(fdays(nstart),ncday(:),&
      "no records for"//date2string(" YYYY-MM-DD hh:mm ",current_date)//&
      "found in "//trim(filename))

  if(MasterProc)&
    write(*,"(3(A,1X),I0)") "Read Pollen dump",trim(filename),"record",nstart

  allocate(data(LIMAX,LJMAX,KMAX_MID))
  do g=1,size(POLLEN_GROUP)
    spc=trim(POLLEN_GROUP(g))
    select case(spc)
      case(BIRCH);iadv=iadv_BIRCH;n_tot=N_TOT_birch
      case(OLIVE);iadv=iadv_OLIVE;n_tot=N_TOT_olive
      case(GRASS);iadv=iadv_GRASS;n_tot=N_TOT_grass
      case default;call CheckStop("Unknown pollen type: "//trim(spc))
    endselect
!------------------------
! pollen adv (not written by Nest_ml)
!------------------------
    call GetCDF_modelgrid(trim(spc),filename,data,&
          1,KMAX_MID,nstart,1,needed=.not.FORECAST,found=found)
    if(DEBUG.and.MasterProc) write(*,dfmt)spc,found
    if(found)xn_adv(iadv,:,:,:)=data(:,:,:)
!------------------------
! pollen_rest
!------------------------
    spc=trim(POLLEN_GROUP(g))//"_rest"
    call GetCDF_modelgrid(trim(spc),filename,Pollen_rest(:,:,g),&
        1,1,nstart,1,needed=.not.FORECAST,found=found)
    if(DEBUG.and.MasterProc) write(*,dfmt)spc,found
    if(found)R(:,:,g)=n_tot-Pollen_rest(:,:,g)
!------------------------
! heatsum
!------------------------
    if(g>size(heatsum,DIM=3))cycle
    spc=trim(POLLEN_GROUP(g))//"_heatsum"
    call GetCDF_modelgrid(trim(spc),filename,heatsum(:,:,g),&
        1,1,nstart,1,needed=.not.FORECAST,found=found)
    if(DEBUG.and.MasterProc) write(*,dfmt)spc,found
  enddo
  deallocate(data)
endsubroutine pollen_read 
!-------------------------------------------------------------------------!
subroutine pollen_dump()
! Write pollen for forecast restart file: heatsum and pollenrest fields
  character(len=len(template_write)) :: filename
  character(len=*), parameter :: dfmt="('Write: ',A20,' [',A,'].')"
  character(len=20) :: spc
  logical     :: fexist,create_var
  logical,save:: overwrite=.true. ! append to existing files if false
  integer     :: ncfileID,i,g,iadv
  type(Deriv) :: def1
  real,allocatable, dimension(:,:,:) :: data ! Data arrays

  if(.not.checkdates(daynumber,"pollen")) return
  if(.not.compare_date(FORECAST_NDUMP,current_date,&
                       outdate(:FORECAST_NDUMP),wildcard=-1))return
  call CheckStop(allocated(Pollen_rest).NEQV.allocated(heatsum),&
    "Pollen: rest/heatsum allocated error")
  if(.not.(allocated(Pollen_rest).and.allocated(heatsum))) return

  filename = date2string(template_write,current_date)
  if(MasterProc) write(*,*) "Write Pollen dump ",trim(filename)
  inquire(file=filename,exist=fexist)
  
  def1%avg=.false.                    ! not used
  def1%index=0                        ! not used
  def1%scale=1.0                      ! not used
  def1%iotype=IOU_INST                ! not used
  
  allocate(data(LIMAX,LJMAX,KMAX_MID))
  ncfileID=-1 ! must be <0 as initial value
  do i=1,2                          ! do first one loop to define the fields,
    create_var=(i==1)               ! without writing them (for performance purposes),
    if(all([create_var,fexist,.not.overwrite,&  ! if the file does not exists already
       template_write/=template_write_IC]))cycle! and pollen has its own restart file
    overwrite=overwrite.and.create_var
    
    do g=1,size(POLLEN_GROUP)
      spc=trim(POLLEN_GROUP(g))
      select case(spc)
        case(BIRCH);iadv=iadv_BIRCH
        case(OLIVE);iadv=iadv_OLIVE
        case(GRASS);iadv=iadv_GRASS
        case default;call CheckStop("Unknown pollen type: "//trim(spc))
      endselect
!------------------------
! pollen adv (not written by Nest_ml)
!------------------------
      def1%class='Advected'           ! written
      def1%unit='mix_ratio'           ! written
      def1%name=trim(spc)             ! written
      if(DEBUG.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      if(.not.create_var)data=xn_adv(iadv,:,:,:)
      call Out_netCDF(IOU_INST,def1,3,KMAX_MID,data,1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,overwrite=overwrite,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
      overwrite=.false.
!------------------------
! pollen_rest
!------------------------
      def1%class='pollen_out'         ! written
      def1%unit="pollengrains"        ! written
      def1%name=trim(spc)//"_rest"    ! written
      if(DEBUG.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      call Out_netCDF(IOU_INST,def1,2,1,Pollen_rest(:,:,g),1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
!------------------------
! heatsum
!------------------------
      if(g>size(heatsum,DIM=3))cycle
      def1%class='pollen_out'         ! written
      def1%unit="degreedays"          ! written
      def1%name=trim(spc)//"_heatsum" ! written
      if(DEBUG.and.MasterProc) write(*,dfmt)def1%name,trim(def1%unit)
      call Out_netCDF(IOU_INST,def1,2,1,heatsum(:,:,g),1.0,CDFtype=CDFtype,&
          out_DOMAIN=out_DOMAIN,create_var_only=create_var,&
          fileName_given=trim(filename),ncFileID_given=ncFileID)
    enddo
  enddo
  if(MasterProc)call CheckNC(nf90_close(ncFileID),"close:"//trim(filename))
  deallocate(data)
endsubroutine pollen_dump   
!-------------------------------------------------------------------------!
function heatsum_calc(t2,T_cutoff) result(ff)
! The temperature needs to be over the cutoff temperature
  real, intent(in) :: t2,T_cutoff
  real             :: ff
  ff = (t2-T_cutoff)*heaviside(t2-T_cutoff)
endfunction heatsum_calc
!-------------------------------------------------------------------------!
function f_wind(u,wstar) result(ff)
! Pollen emission increases with wind speeds up to 5 m/s (u_sat).
! This term cannot be higher than 1.5 (u_max).
  real, intent(in) :: u,wstar
  real             :: ff
  real, parameter  :: u_max=1.5,u_sat=5.0
  ff = u_max - exp(-(u+wstar)/u_sat)
endfunction f_wind
!-------------------------------------------------------------------------!
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
  endif
endfunction f_cond
!-------------------------------------------------------------------------!
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
  endif
endfunction f_in
!-------------------------------------------------------------------------!
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
  endif
endfunction f_out
!-------------------------------------------------------------------------!
function f_fade_in(value_rel,uncert_rel_) result(ff)
! Computes the linear fade-in function. It is 0 at vaule_rel = 1.-uncert_rel
! Grows to 0.5 at vaule_rel = 1 and grows to 1 at value_rel = 1.+uncert_rel
  real, intent(in) :: value_rel,uncert_rel_
  real             :: uncert_rel,ff

  uncert_rel = max(uncert_rel_,1e-5) ! avoid zero uncertainty
  ff = min(1.0,max(0.0,(value_rel-1.0+uncert_rel)/(2.0*uncert_rel)))
endfunction f_fade_in
!-------------------------------------------------------------------------!
function f_fade_out(value_rel,uncert_rel_) result(ff)
! Computes the linear fade-in function. It is 0 at vaule_rel = 1.-uncert_rel
! Grows to 0.5 at vaule_rel = 1 and grows to 1 at value_rel = 1.+uncert_rel
  real, intent(in) :: value_rel,uncert_rel_
  real             :: uncert_rel,ff

  uncert_rel = max(uncert_rel_,1e-5) ! avoid zero uncertainty
  if (value_rel > 1.0 + uncert_rel) then
    ff = 0.0
  elseif (value_rel < 1.0 - uncert_rel) then
    ff = 1.0
  else 
    ff = (1.0+uncert_rel-value_rel)/(2.0*uncert_rel)
  endif
endfunction f_fade_out
!-------------------------------------------------------------------------!
endmodule Pollen_ml
  
  
  
  
  
  
  
  
  
  
  
  
  
