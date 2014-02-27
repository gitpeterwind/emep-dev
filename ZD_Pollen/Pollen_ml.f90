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
  use CheckStop_ml,         only: CheckStop
  use ChemChemicals_ml,     only: species
  use Functions_ml,         only: heaviside
  use GridValues_ml ,       only: glon, glat
  use Landuse_ml,           only: LandCover
  use LocalVariables_ml,    only: Grid, Sub
  use MetFields_ml,         only: surface_precip, ws_10m ,rh2m,t2_nwp,&
                                  foundws10_met,foundprecip,pr,u_ref,z_bnd,z_mid
  use MicroMet_ml,          only: Wind_at_h
  use ModelConstants_ml,    only: KMAX_MID, KMAX_BND, nmax, nstep, FORECAST, &
                                  METSTEP, MasterProc, IOU_INST, RUNDOMAIN, &
                                  dt => dt_advec, IIFULLDOM, JJFULLDOM
  use Nest_ml,              only: outdate
  use NetCDF_ml,            only: ReadField_CDF,Out_netCDF,printCDF,nc_check=>check
  use Par_ml,               only: limax, ljmax, MAXLIMAX,MAXLJMAX,GIMAX,GJMAX, &
                                  me, MSG_READ9, IRUNBEG,JRUNBEG  ! => x, y dimensions
  use OwnDataTypes_ml,      only: Deriv
  use Setup_1dfields_ml,    only: rcemis 
  use SmallUtils_ml,        only: find_index
  use TimeDate_ml,          only: current_date,daynumber,date,day_of_year
  use TimeDate_ExtraUtil_ml,only: date2string
!-------------------------------------
  implicit none
  private
  public:: pollen_flux,pollen_dump,pollen_read

  !** 1) Public (saved) Variables from module:

  real, public,save , allocatable, dimension(:,:)::&
    AreaPOLL,     & ! emission of pollen 
    heatsum,      & ! heatsum, needs to be remembered for forecast
    R,            & ! Summed pollen release
    Pollen_rest,  & ! what is available pollen pr m2, fr forecast
    Pollen_left     ! amount of pollen left in catkins, relative amount... 0:sr 1:end 

 !  logical, public,save, dimension(FORECAST_NDUMP) ::  dump_pollen = .true.
 !  integer                                       ::  fdump = 1
  character(len=*), parameter :: &
    pollen_data="pollen_data.nc",&
    template_pollen_read ="POLLEN_IN.nc",&
    template_pollen_write="POLLEN_OUT.nc"

  type(date), parameter :: &
    date_first=date(-1,3,1,1,0),&
    date_last =date(-1,8,1,1,0)

  integer, save :: day_first=1,day_last=366
!-------------------------------------------------------------------------!
contains
!-------------------------------------------------------------------------!
function checkdates(nday) result(ok)
  integer, intent(in) :: nday
  logical             :: ok
  logical, save :: first_call = .true.
  if(first_call)then
    day_first=day_of_year(current_date%year,date_first%month,date_first%day)
    day_last =day_of_year(current_date%year,date_last%month ,date_last%day)
    first_call = .false.
  endif
  ok=(day_first<=nday).and.(nday<=day_last)
endfunction checkdates
!-------------------------------------------------------------------------!
subroutine pollen_flux(i,j,debug_flag) 
  implicit none
  integer, intent(in) :: i,j    ! coordinates of column
  logical, intent(in) :: debug_flag

  real, parameter :: Z10 = 10.0  ! 10m height

  integer, save :: ipoll,inat_POLL,itot_POLL,poll_day  !DS replaces Pollen_b
  real,    save :: moleccm3s_2_kgm2h 
  logical, save :: first_call = .true.
   
  integer :: ii, jj, nlu, ilu, lu
  real :: scale,lim,invdz,rcpoll,n2m,u10,prec,relhum

  real, save, allocatable, dimension(:,:) :: &
    birch_frac, & ! Fraction of birch, read in
    h_c,        & ! Fraction of birch, read in
    corr,       & ! correction field for p.emission, read in
    h_day         ! Temperature summed over a day
  integer,save,allocatable, dimension(:,:) :: &
    p_day

  if(.not.checkdates(daynumber)) return
! Reas in the different fields

  if(first_call) then 
    inat_POLL = find_index("POLLEN_B",EMIS_BioNat(:))
    itot_POLL = find_index("POLLEN_B",species(:)%name)

    allocate(pollen_left(MAXLIMAX,MAXLJMAX),p_day(MAXLIMAX,MAXLJMAX))
    pollen_left(:,:) = 1
    p_day(:,:) =current_date%day

    if(FORECAST) call pollen_read()
    if(.not.allocated(heatsum))then
      allocate(heatsum(MAXLIMAX,MAXLJMAX))
      heatsum(:,:) = 0.00
    endif
    if(.not.allocated(Pollen_rest))then
      allocate(Pollen_rest(MAXLIMAX,MAXLJMAX))
      Pollen_rest(:,:) = N_TOT
    endif
    if(.not.allocated(R))then
      allocate(R(MAXLIMAX,MAXLJMAX))
      R(:,:) = 0.0
    endif

    allocate(AreaPOLL(MAXLIMAX,MAXLJMAX),h_day(MAXLIMAX,MAXLJMAX))
    AreaPOLL(:,:)=0.0
    h_day(:,:)=0.0
      
    allocate(birch_frac(MAXLIMAX,MAXLJMAX),h_c(MAXLIMAX,MAXLJMAX),corr(MAXLIMAX,MAXLJMAX))
    call ReadField_CDF(pollen_data,'birch',birch_frac,2, &
       interpol = 'conservative',needed=.true.,debug_flag=DEBUG,UnDef=0.0)
    call ReadField_CDF(pollen_data,'h_c',h_c,2, &
       interpol = 'conservative',needed=.true.,debug_flag=DEBUG,UnDef=0.0)
    call ReadField_CDF(pollen_data,'cross',corr,2, &
       interpol = 'conservative',needed=.true.,debug_flag=DEBUG,UnDef=0.0)

    ! reduce birch from 60 degrees north:
    forall(ii=1:MAXLIMAX,jj=1:MAXLIMAX,glat(ii,jj)>=60.0) &
      birch_frac(ii,jj)=birch_frac(ii,jj)&
!!                     *max(0.3,1.0-0.005*(glat(ii,jj)-60.0))
                       *(1-0.7*(glat(ii,jj)-60)/10)
   ! start of pollen forecast
    if(daynumber==day_first.and.MasterProc)&
      write(*,*) "POLLEN setup ",date2string("YYYY-MM-DD hh:mm",current_date)

    first_call = .false. 
  endif  !first_call
 
  if(birch_frac(i,j)==0.0 .or. h_c(i,j)==0.0)then
    EmisNat(inat_POLL,i,j)      = 0.0
    rcemis (itot_POLL,KMAX_MID) = 0.0
    return
  endif

 !  Heatsum calculations
 !  Sums up the temperatures that day for each timestep (20 min).
 !  The heatsum are then divided by timesteps, and added using the 
 !  function heatsum_calc. 
  if(nstep==1) then ! Heatsum calculation only each meteorological timestep
    if(p_day(i,j)==current_date%day) then
      h_day(i,j) = h_day(i,j) + t2_nwp(i,j,1)
      if(current_date%hour==21) &
        heatsum(i,j)=heatsum(i,j)+heatsum_calc((h_day(i,j)/(24/METSTEP)),T_cutoff)
    else 
      p_day(i,j) = current_date%day  
      h_day(i,j) = t2_nwp(i,j,1)
    endif
    ! End of heatsum calculations
  endif ! heatsum calc each timestep
  ! if heatsum is over heatsum threhold for the grid cell, the pollen
  ! emission can start calculating 
  lim = (1-prob_in)*h_c(i,j)
  if(heatsum(i,j)<lim           &
     .or.t2_nwp(i,j,1)<T_cutoff &
     .or. prec>prec_max         &
     .or. R(i,j)>N_TOT*corr(i,j)&
     .or. heatsum(i,j)>h_c(i,j)+dH_d ) then 
    EmisNat(inat_POLL,i,j)      = 0.0
    rcemis (itot_POLL,KMAX_MID) = 0.0
    return
  endif

  ! calculating grains to molecular weight   grains/m2/s --> mol/cm3/s 
  invdz      = 1e-6/Grid%DeltaZ         ! 1/dZ [1/cm3]
  n2m        = invdz*grain_wt*AVOG/species(itot_POLL)%molwt

  ! For EmisNat, need kg/m2/h from molec/cm3/s
  moleccm3s_2_kgm2h = Grid%DeltaZ * 1e6 * 3600.0  &! /cm3/s > /m2/hr
                     /AVOG * 1e-6  ! kg  after *MW

  rcpoll        = 0.0
  AreaPOLL(i,j) = 0.0
  scale         = 0.0

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
    if(.not.foundws10_met)&
      u10=Wind_at_h(Grid%u_ref,Grid%z_ref,Z10,Sub(lu)%d,Sub(lu)%z0,Sub(lu)%invL)

    if(relhum>=0.0 .and. u10>=0.0) exit
  enddo

  ! precipitation
  if(foundprecip) then 
    prec=sum(pr(i,j,:))*60  !mm/h
  else
    prec=surface_precip(i,j)
  endif
!<----------------------------------------------------------------------------------------------- 
!  Scale is calculated from the different meteorological conditions using functions below  
  scale=corr(i,j)*(t2_nwp(i,j,1)-T_cutoff)/dH &
       *f_wind(u10,Grid%wstar)              & ! wind dependece
       *f_in(heatsum(i,j),h_c(i,j),PROB_IN) & ! probability for flowering to start
       *f_out(R(i,j),N_TOT,PROB_OUT)        & ! probability for flowering to end
       *f_cond(relhum,RH_LOW,RH_HIGH)       & ! rh dependence
       *f_cond(prec,PREC_MIN,PREC_MAX)        ! precipitation dependence

    R(i,j)=R(i,j)+N_TOT*scale*dt              ! R is to check the remains of pollen available
    rcpoll       =N_TOT*scale*birch_frac(i,j) ! The pollen grains production [grains/m2/s]
    AreaPOLL(i,j)=rcpoll*3600                 ! [grains/m2/h]
!<-----------------------------------------------------------------------------------------------
  !!! Need to convert grains/m2/s --> mol/cm3/s
    rcpoll = rcpoll * n2m             ! [mol/cm3/s]
    Pollen_rest(i,j) = N_TOT - R(i,j) ! pollen grains left after release [grains/m2]
    pollen_left(i,j) = Pollen_rest(i,j)/N_TOT ! idem [ ]
      
    if(debug_flag.or.DEBUG)&
      write(6,'(a,2f10.5,f10.5)')"Result: ",rcpoll,R(i,j)/N_TOT,AreaPOLL(i,j)

  EmisNat(inat_POLL,i,j)      = rcpoll*moleccm3s_2_kgm2h*species(itot_POLL)%molwt
  rcemis (itot_POLL,KMAX_MID) = rcpoll
endsubroutine pollen_flux
!-------------------------------------------------------------------------! 
subroutine pollen_read()
! Read pollen for forecast restart file: heatsum and pollenrest fields
  use netcdf
  logical, save :: first_call = .true.
  integer :: fid,did,vid
  logical :: fexist
  integer :: nlon,nlat
  character(len=len(template_pollen_read)) :: filename
  real,  allocatable, dimension(:,:) :: h_gl,r_gl

  if(.not.checkdates(daynumber)) return
  if(.not.first_call) return
  first_call=.false.

  filename=date2string(template_pollen_read,current_date)
  inquire(file=filename,exist=fexist)
  if(FORECAST.and..not.fexist)then
    if(MasterProc) write(*,*)"Warning Polen dump file not found: "//trim(filename)
    return
  endif
  call CheckStop(.not.fexist,"Did not find file: "//trim(filename))

  allocate(h_gl(GIMAX,GJMAX),r_gl(GIMAX,GJMAX))
  if(MasterProc)then
    if(DEBUG)write(*,*) "Read Pollen dump ",trim(filename),MAXLIMAX,MAXLJMAX,GIMAX,GJMAX
    call nc_check(nf90_open(path=trim(filename),mode=nf90_nowrite,ncID=fid),"Pollen_dump")
    if(DEBUG)then
      call nc_check(nf90_inq_dimid(ncID=fid,name="lon",dimID=did))
      call nc_check(nf90_inquire_dimension(ncID=fid,dimID=did,len=nlon))
      call nc_check(nf90_inq_dimid(ncID=fid,name="lat",dimID=did))
      call nc_check(nf90_inquire_dimension(ncID=fid,dimID=did,len=nlat))
      write(*,*) "POLLEN TESTER",nlon,nlat
    endif
    call nc_check(nf90_inq_varid(ncID=fid,name="Heatsum"    ,varID=vid),"Heatsum")
    call nc_check(nf90_get_var(ncID=fid,varid=vid,values=h_gl,count=(/GIMAX,GJMAX/)))
    call nc_check(nf90_inq_varid(ncID=fid,name="Pollen_rest",varID=vid),"Pollen_rest")
    call nc_check(nf90_get_var(ncID=fid,varid=vid,values=r_gl,count=(/GIMAX,GJMAX/)))
  endif
  if(allocated(heatsum))    deallocate(heatsum)
  if(allocated(Pollen_rest))deallocate(Pollen_rest)
  if(allocated(R))          deallocate(R)
  allocate(heatsum(MAXLIMAX,MAXLJMAX),Pollen_rest(MAXLIMAX,MAXLJMAX),R(MAXLIMAX,MAXLJMAX))
  call global2local(h_gl,heatsum    ,MSG_READ9,1,GIMAX,GJMAX,1,IRUNBEG,JRUNBEG)
  call global2local(r_gl,Pollen_rest,MSG_READ9,1,GIMAX,GJMAX,1,IRUNBEG,JRUNBEG)
  R(:,:)=N_TOT-Pollen_rest(:,:)
  deallocate(h_gl,r_gl)
! call printCDF("HEAT_GLOB",heatsum,"not")
! call printCDF("REST_GLOB",Pollen_rest,"not") 
endsubroutine pollen_read 
!-------------------------------------------------------------------------!
subroutine pollen_dump ()
! Write pollen for forecast restart file: heatsum and pollenrest fields
  character(len=len(template_pollen_write)) :: filename
  integer     :: i0=60,j0=11,i1=107,j1=58
  type(Deriv) :: def1

! call printCDF("HEAT2",heatsum,"not")
  if(.not.checkdates(daynumber)) return
  if(.not.any((current_date%year   ==outdate%year   .or.outdate%year   ==-1).and.&
              (current_date%month  ==outdate%month  .or.outdate%month  ==-1).and.&
              (current_date%day    ==outdate%day    .or.outdate%day    ==-1).and.&
              (current_date%hour   ==outdate%hour   .or.outdate%hour   ==-1).and.&
              (current_date%seconds==outdate%seconds.or.outdate%seconds==-1)))return
  call CheckStop(allocated(Pollen_rest).NEQV.allocated(heatsum),&
    "Pollen: rest/heatsum allocated error")
  if(.not.(allocated(Pollen_rest).and.allocated(heatsum))) return

  filename = date2string(template_pollen_write,current_date)
  write(*,*) "Pollen filename ",trim(filename)

  i0=RUNDOMAIN(1);j0=RUNDOMAIN(3)
  i1=RUNDOMAIN(2);j1=RUNDOMAIN(4)

  def1%class='pollen_out' !writtenclass='pollen_out' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=1.0        !not used

  def1%name=trim("Pollen_rest")   ! written
  def1%unit=trim("pollengrains")
  call Out_netCDF(IOU_INST,def1,2,1,Pollen_rest,1.0,CDFtype=4,&
        ist=i0,jst=j0,ien=i1,jen=j1,fileName_given=trim(filename)) ! pollen_rest

  def1%name=trim("Heatsum")   ! written
  def1%unit=trim("degreedays")
  call Out_netCDF(IOU_INST,def1,2,1,heatsum,1.0,CDFtype=4,&
        ist=i0,jst=j0,ien=i1,jen=j1,fileName_given=trim(filename)) ! heatsum
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
endmodule Pollen_ml
  
  
  
  
  
  
  
  
  
  
  
  
  
