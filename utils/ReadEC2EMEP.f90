!NBNB neeeds 35GB memory per CPU

!NEW: Read interpolation coefficients for SMI from disk (or write them if they do not exists from before)
!     Simply link the coefficients from anothe rdirectory. The coeff do not change (for a constant IFS resolution).
!NEW: meteofields_class%fscale is used as scaling factor
!Hourly fields are put in a separate file (meteo_h)

!mpif90 -shared-intel -o  ReadEC -Wl,-rpath,/global/apps/netcdf/4.2.1.1/intel/13.0.0/lib  -O3 -r8 -I/global/apps/netcdf/4.2.1.1/intel/13.0.0/include -L/global/apps/netcdf/4.2.1.1/intel/13.0.0/lib ReadEC2EMEP.f90  -lnetcdf -lnetcdff

!mpif90 -shared-intel  -CB -r8  -debug-parameters all -traceback  -ftrapuv -g -fpe0 -O0 -o  ReadEC  -r8   -I/global/apps/netcdf/4.2.1.1/intel/13.0.0/include -L/global/apps/netcdf/4.2.1.1/intel/13.0.0/lib ReadEC2GLOBAL.f90  -lnetcdf -lnetcdff

!qsub -I -lpmem=12GB -lnodes=4 -lwalltime=2:0:0 (Check "top" for Memory!)
!module load mpt
!mpif90 -shared-intel -o  ReadEC  -O3 -r8  -I/sw/sdev/Modules/netcdf/netcdf-4.1.3-intel.11.1.073/include  -L/sw/sdev/Modules/netcdf/netcdf-4.1.3-intel.11.1.073/lib -Wl,-rpath,/sw/sdev/Modules/netcdf/netcdf-4.1.3-intel.11.1.073/lib ReadEC2ECA_12hours.f90  -lnetcdf -lnetcdff
!mpiexec ReadEC

Program callreadCDF
  !example of how to use the routines

  !use ReadCDF_ml
  implicit none
  include 'mpif.h'
  INTEGER STATUS(MPI_STATUS_SIZE),INFO,me, nproc
  integer :: i,j,k,n,month,Nrec,nstart,nfetch,IEMEP_shift,JEMEP_shift,ktest
  integer :: GIMAX,GJMAX,KMAX_MID,nrecords
  integer :: GIMAX_X,GJMAX_X
  real :: period
  character(len=100) pathdir,fileName,fileName_h,sname_req,namefield,validity,pathdirout,fileNameout
  character(len=100)::fileNameagg,fileNamespr,fileNamesur,fileNameuv,meteofields_3Dending_v
  character(len=120)::fileNameash,fileNamepre,fileNamesnw,fileNameEC,fileNameEC_v,command
  character*50 ::units
  integer, parameter::Nrect=31*24*12
  integer timev(Nrect)
  real, allocatable, dimension(:,:,:,:)::varEC,var1,var2,varEC_u,varEC_v
  real, allocatable, dimension(:,:,:)::varsEC,vars1,vars2,psEC,t2mEC,buffer_EC
  real, dimension(10000)::lon_EC,lat_EC,lon_EC_X,lat_EC_X,level
  real ::along,alat,delta,fac,dist1,dist2,dist12
  !OUTPUT GRID SIZE
  integer,parameter :: KMAX1=1,Nhh_8=8,Nhh_24=24
      integer,parameter :: IMAX=1200,JMAX=520,KMAX=37 !Lat Lon
!  integer,parameter :: IMAX=132,JMAX=159,KMAX=37 !EMEP extended
 ! integer,parameter :: IMAX=132*2,JMAX=159*2,KMAX=20 !EMEP extended 25km
  !    integer,parameter :: IMAX=140,JMAX=260,KMAX=20 !UK 
  integer :: Nhh_half,Nhh
  real    :: fi,an,xp,yp,u_ll,v_ll,GRIDWIDTH,alpha
  real    :: xp_u,yp_v
  real, dimension(imax,jmax)    :: ir,jr
  real, dimension(imax+1,jmax+1)    :: gl,gb,angle1
  real, dimension(imax+1,jmax+1)    :: gl_v,gb_v,angle1_v, gl_u,gb_u,angle1_u
  real, dimension(imax,jmax,KMAX)    :: EMEP_3D,EMEP_3D_2
  real, dimension(imax,jmax)    :: P_surf,EMEP_2d
  real, parameter :: PI=3.14159265358979323
  real, parameter :: PS0=1013.25!hPasPS
  real, parameter :: PT=100!hPa
  real,parameter ::   R = 287.0, CP = 1004.0, XKAP = R/CP
  integer :: conv_i(imax),conv_j(jmax),conv_j2(jmax),conv_i1(imax),conv_j1(jmax),conv_i2(imax),conv_i_u(imax),conv_j_u(jmax),conv_i_v(imax),conv_j_v(jmax)
  integer :: conv_i1_u(imax),conv_j1_u(jmax),conv_i1_v(imax),conv_j1_v(jmax),conv_i2_u(imax),conv_j2_u(jmax),conv_i2_v(imax),conv_j2_v(jmax)
  real::lon(imax),lat(jmax),lon_u(imax),lat_v(jmax),di1(imax),dj1(jmax),di1_u(imax),dj1_u(jmax),di1_v(imax),dj1_v(jmax)
  real ::EC_lonstart,EC_latstart,EC_ddeg_lon,EC_ddeg_lat
  real ::EC_ddeg_lon_uv,EC_ddeg_lat_uv
  real ::lonstart,latstart,ddeg,ddeg_lat,ddeg_lon
   integer :: iotyp,IOU_INST=1,IOU_HOUR=5, IOU_YEAR=2 ,IOU_MON=3, IOU_DAY=4  
  integer :: wrec,ndate(4,Nhh_8),ndate_h(4,Nhh_24),nn,nh,ndim,nrecord,nrecord_read
  real ::A(61),B(61),ph(61),Press(61),Alt(61),r1,r2,r3,r4,r5
  real ::Ah(61),Bh(61),Alth(61),sigma_EC(61),sigma_half_EC(61)
  real ::sigma_EMEP(KMAX),xx,EC2EMEP1(KMAX),EC2EMEP2(KMAX)
  real ::EC2EMEP_half1(KMAX),EC2EMEP_half2(KMAX),EC_index_half(KMAX)
  real ::EC_half2EMEP_half1(KMAX),EC_half2EMEP_half2(KMAX),EC_index_half_half(KMAX)
  real ::sigma_bnd_EMEP(KMAX+1)
  integer ::I_EC,J_EC,I_EC1,I_EC2,J_EC1,J_EC2,EC_index(KMAX),k1_EC,k2_EC,deltaindex
  real ::delta_t,ref_latitude
  integer :: yyyy,mm,dd,hh,yyyy1,mm1,dd1
  integer :: nday,nmdays(12)
  real ::di_EC1,di_EC,dj_EC1,dj_EC,d00,d10,d01,d11
  real :: dy,dx,rp,rb
  integer :: N_EC(KMAX),I_EC_N(KMAX,70),I_EC_low,I_EC_high,kk,k_EC,i2d,N2d,i3d,N3d
  real :: X_EC_N(KMAX,70)
  real :: glmin,glmax,rl,om,sum(Nhh_24+1)
  logical :: vector
  integer :: megz,icount,icountready,t2msaved,year,iter,nland,i_EC_cyclic
  real, parameter :: T0 = 273.15   ! zero degrees Celsius in Kelvin 
  real:: fscale
  real :: hyai(KMAX+1),hybi(KMAX+1),hyam(KMAX),hybm(KMAX),xm(imax,jmax)

  integer ::maxnland
  integer*2, allocatable, dimension(:,:) :: N_smi
  integer*2, allocatable, dimension(:,:,:) :: I_smi,J_smi

  type :: meteofields_class
     character*50 ::emep_name
     character*50 ::EC_name
     character*50 ::ending
     character*50 ::unit      !in meteo file (not always the same as EC unit)
     real*8         ::fscale
     character*50 ::treatment
     logical    ::  hourly
     integer      ::ndim
  end type meteofields_class

  type(meteofields_class), dimension(23) :: meteofields_2D =&
       (/ & !emep_name  EC_name  ending unit  fscale treatment hourly ndim
       meteofields_class("surface_pressure", "lnsp"  , "_surfpr.nc",      "hPa" ,1.0,  "log", .false.  ,2  )&
       ,meteofields_class("temperature_2m",   "v2t", "_surf.nc",     "K"  ,1.0, " ",   .false.  ,2  )&
       ,meteofields_class("surface_flux_sensible_heat", "sshf","_asurfhr.nc",        "W/m2" ,1.0/3600,  "accumulated",   .true.  ,2  )&
       ,meteofields_class("surface_flux_latent_heat", "slhf", "_asurfhr.nc",       "W/m2" ,1.0/3600,  "accumulated",   .true.  ,2  )&
       ,meteofields_class("surface_stress", "ewss", "_asurfhr.nc",     "N/m2" ,1.0/3600,  "vector_acc",   .true.  ,2  )&!og 73                            !6
       ,meteofields_class("surface_stress", "nsss","_asurfhr.nc",     "N/m2" ,1.0/3600,  "vector_acc",   .true.  ,2  )&!og 72                           !6
       ,meteofields_class("surface_flux_sensible_heat", "sshf","_asurf.nc",       "W/m2" ,1.0/10800,  "accumulated",   .false.  ,2  )&
       ,meteofields_class("surface_flux_latent_heat", "slhf", "_asurf.nc",     "W/m2" ,1.0/10800,  "accumulated",   .false.  ,2  )&
       ,meteofields_class("surface_stress", "ewss", "_asurf.nc",    "N/m2" ,1.0/10800,  "vector_acc",   .false.  ,2  )&!og 73                            !6
       ,meteofields_class("surface_stress", "nsss","_asurf.nc",     "N/m2" ,1.0/10800,  "vector_acc",   .false.  ,2  )&!og 72                           !6
       ,meteofields_class("snow_depth",  "sd", "_surf.nc",     "m" ,1.0,  " ",   .false.  ,2  )&
       ,meteofields_class("fraction_of_ice","ci", "_surf.nc",     "%" ,100.0,  " ",   .false.  ,2  )&!9
       ,meteofields_class("u10",  "v10u",  "_surf.nc",    "m/s" ,1.0,  " ",   .false.  ,2  )&
       ,meteofields_class("v10",    "v10v", "_surf.nc",   "m/s" ,1.0,  " ",   .false.  ,2  )&!13
       ,meteofields_class("large_scale_precipitations",   "lsp", "_asurf.nc",   "m/s" ,1.0/10800,  "accumulated",   .false.  ,2  )&
       ,meteofields_class("convective_precipitations",   "cp", "_asurf.nc",    "m/s" ,1.0/10800,  "accumulated",   .false.  ,2  )&
       ,meteofields_class("soil_water_content",  "swvl1",   "_surf.nc",    "m3/m3" ,1.0,  " ",   .false.  ,2  )&
       ,meteofields_class("deep_soil_water_content",  "swvl3",   "_add.nc",   "m3/m3" ,1.0,  " ",   .false.  ,2  )&
       ,meteofields_class("sea_surface_temperature","sstk", "_surf.nc","K" ,1.0,  "daily",   .false.  ,2  )&!8
       ,meteofields_class("SMI1",  "swvl1", "_smi1.nc",     "" , 1.0, "SetSea",   .false.  ,2  )&
       ,meteofields_class("SMI3",  "swvl3", "_smi3.nc",     "" , 1.0, "SetSea",   .false.  ,2  )&
       ,meteofields_class("relative_humidity_2m",  "v2d", "_surf.nc",     "%" , 1.0, "dewpoint",   .false.  ,2  )&
       ,meteofields_class("MSLP",   "msl", "_surf.nc",     "Pa" , 1.0, " ",   .false.  ,2  )&
       /)

  type(meteofields_class), dimension(9) :: meteofields_3D =&
       (/ & !emep_name EC_name ending felt_code unit fscale validity ndim
!       meteofields_class("u_wind", "u"  ,"_uv3d.nc",     "m/s" ,1.0,  "stagg",   .false.  ,3  )&!30
!       ,meteofields_class("v_wind","v","_uv3d.nc",         "m/s" ,1.0,  "stagg",   .false.  ,3  )&!
       meteofields_class("u_wind", "u"  ,"_uwin3d.nc",    "m/s" ,1.0,  "NOTstagg",   .false.  ,3  )&!30
       ,meteofields_class("v_wind","v","_vwin3d.nc",       "m/s" ,1.0,  "NOTstagg",   .false.  ,3  )&!NB careful with definition of ending: not used!
!       ,meteofields_class("specific_humidity",  "q","_air3d.nc",   "kg/kg" ,1.0,  " ",   .false.,   3  )&!
!       ,meteofields_class("potential_temperature", "t","_air3d.nc", "K" ,1.0,  "absolute",   .false.,   3  )&!33
!       ,meteofields_class("cloudwater", "clwc", "_air3d.nc",  "kg/kg" ,1.0,  " ",   .false.,   3  )&!
!       ,meteofields_class("3D_cloudcover", "cc", "_air3d.nc",   "%" ,100.0,  " ",   .false.,   3  )&!37
!       ,meteofields_class("total_convective_precipitation_profile", "108.128", "_acc3d.nc",    "kg/m2" ,1.0,  "acc_precip",   .false. ,   3  )&!
!       ,meteofields_class("total_stratiform_precipitation_profile", "109.128", "_acc3d.nc",    "kg/m2" ,1.0,  "acc_precip",   .false. ,   3  )&!
!       ,meteofields_class("turbulent_diffusion_coefficient_for_heat", "110.128",  "_acc3d.nc",  "m2" ,1.0,  "accumulated",   .false. ,   3  )&!??
!       ,meteofields_class("convective_updraft_flux",  "104.128", "_conv.nc",  " kg/m2/s" ,1.0/10800.0,  "acc_bnd",   .false. ,   3  )&!
!       ,meteofields_class("convective_downdraft_flux",   "105.128", "_conv.nc",   "kg/m2/s" ,1.0/10800.0,  "acc_bnd",   .false. ,   3  )&!
!       ,meteofields_class("convective_updraft_detrainment", "106.128", "_conv.nc",   "kg/m3/s" ,1.0/10800.0,  "accumulated",   .false. ,   3  )&!
!       ,meteofields_class("convective_downdraft_detrainment", "107.128", "_conv.nc",    "kg/m3/s" ,1.0/10800.0,  "accumulated",   .false. ,   3  )&!
!       ,meteofields_class("Kz", "110.128",  "_acc3d.nc",  "m2/s" ,1.0/10800.0,  "accumulated",   .false. ,   3  )&!??
       ,meteofields_class("specific_humidity",  "q","_shum3d.nc",      "kg/kg" ,1.0,  " ",   .false.  ,3  )&!
       ,meteofields_class("potential_temperature", "t","_temp3d.nc",    "K" ,1.0,  "absolute",   .false.  ,3  )&!33
       ,meteofields_class("cloudwater", "clwc", "_clwc3d.nc",    "kg/kg" ,1.0,  " ",   .false.  ,3  )&!
       ,meteofields_class("3D_cloudcover", "cc", "_ccov3d.nc",    "%" ,100.0,  " ",   .false.  ,3  )&!37
       ,meteofields_class("total_convective_precipitation_profile", "108.128", "_acpr3d.nc",      "kg/m2" ,1.0,  "acc_precip",   .false.  ,3  )&!
       ,meteofields_class("total_stratiform_precipitation_profile", "109.128", "_aspr3d.nc",      "kg/m2" ,1.0,  "acc_precip",   .false.  ,3  )&!
       ,meteofields_class("etadot", "etadot",  "_etad3d.nc",    "Pa s**-1" ,1.0,  " ",   .false.  ,3  )&!
!       ,meteofields_class("convective_updraft_flux",  "104.128", "_aumf3d.nc",     " kg/m2/s" ,1.0/10800.0,  "acc_bnd",   .false.  ,3  )&!
!       ,meteofields_class("convective_downdraft_flux",   "105.128", "_admf3d.nc",    "kg/m2/s" ,1.0/10800.0,  "acc_bnd",   .false.  ,3  )&!
!       ,meteofields_class("convective_updraft_detrainment", "106.128", "_audr3d.nc",    "kg/m3/s" ,1.0/10800.0,  "accumulated",   .false.  ,3  )&!
!       ,meteofields_class("convective_downdraft_detrainment", "107.128", "_addr3d.nc",     "kg/m3/s" ,1.0/10800.0,  "accumulated",   .false.  ,3  )&!
!       ,meteofields_class("turbulent_diffusion_coefficient_for_heat", "110.128",  "_atdf3d.nc",   "m2" ,1.0,  "accumulated",   .false.  ,3  )&!??
       /)
  meteofields_3Dending_v="_vwin3d.nc"
!  meteofields_3Dending_v="_uv3d.nc"

  !104.128 - Updraught mass flux          [kg/m2]
  !105.128 - Downdraught mass flux        [kg/m2]
  !106.128 - Updraught detrainment rate   [kg/m3]
  !107.128 - Downdraught detrainment rate [kg/m3]
  !108.128 - Total convective precipitation profile  [kg/m2]
  !109.128 - Total stratiform precipitation profile  [kg/m2]
  !110.128 - Turbulent diffusion coefficient for heat  [m2]
  !
  !109.128 kalles også for "large scale precip"
  !110.128 er noe som David ville ha.


  !  pathdir='/fou1/emep/Meteorology/Hemispheric/ERA/'
  !pathdir='/global/work/alvarov/ECDATA/2008/'
year=2015
   ! pathdir='/home/metno/mifaab/work/ECDATA/2012/'
   !pathdir='/project/data1/emep/ECDATA/2012/'
   pathdir='/global/work/mifaab/emep/ECDATA/2015/'
  !pathdirout='/global/work/mifapw/emep/Data/EECCA/metdata_EC/2009/'
  pathdirout='./'

  nproc=1
  call MPI_INIT(info)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, ME, INFO)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, INFO)

 

  !OUTPUT GRID PARAMETERS
  !EXTENDED EMEP with EECCA
  !ref_latitude=60.
  !fi =-32.
  !GRIDWIDTH=50000.
  !an= 6.370e6*(1.0+sin(ref_latitude*PI/180.))/GRIDWIDTH
  !xp=8.0
  !yp=110.0

  !EXTENDED EMEP with EECCA in 25 km resolution
 ! ref_latitude=60.
 ! fi =-32.
 ! GRIDWIDTH=25000.
 ! an= 6.370e6*(1.0+sin(ref_latitude*PI/180.))/GRIDWIDTH
 ! xp=2*8.0-0.5
 ! yp=2*110.0-0.5

  !u grid is shifted by 1/2 gridcell to the right
  !xp_u=xp-0.5
  !v grid is shifted by 1/2 gridcell to the top
  !yp_v=yp-0.5

  !glmin = 0.0
  !glmax = glmin + 360.0

  !do j = 1, JMAX+1
  !   dy  = yp - j  
  !   do i = 1, IMAX+1
  !      dx = i - xp    
  !      rp = sqrt(dx*dx+dy*dy)           ! => distance to pole
  !      rb = 90.0 - 360./PI * atan(rp/AN)  ! => latitude
  !      !        xm(i,j) = 0.5*(1.0+sin(ref_latitude*PI/180.))*(1.0 + (dx*dx+dy*dy)/AN**2 )
  !      if (rp >  1.0e-10)then
  !         rl = fi + 180.0/PI *atan2(dx,dy)
  !      else
  !         rl = fi
  !      endif
  !      if (rl <  glmin)   rl = rl + 360.0
  !      if (rl >  glmax)   rl = rl - 360.0
  !      gl(i,j)=rl                     !     longitude
  !      gb(i,j)=rb                     !     latitude
  !   end do ! i
  !end do ! j
  !do j = 1, JMAX+1
  !   dy  = yp_v - j  
  !   do i = 1, IMAX+1
  !      dx = i - xp    
  !      rp = sqrt(dx*dx+dy*dy)           ! => distance to pole
  !      rb = 90.0 - 360./PI * atan(rp/AN)  ! => latitude
  !      if (rp >  1.0e-10)then
  !         rl = fi + 180.0/PI *atan2(dx,dy)
  !      else
  !         rl = fi
  !|      endif
  !      if (rl <  glmin)   rl = rl + 360.0
  !      if (rl >  glmax)   rl = rl - 360.0
  !      gl_v(i,j)=rl                     !     longitude
  !      gb_v(i,j)=rb                     !     latitude
  !   end do ! i
  !end do ! j
  !do j = 1, JMAX+1
  !   dy  = yp - j  
  !   do i = 1, IMAX+1
  !      dx = i - xp_u    
  !      rp = sqrt(dx*dx+dy*dy)           ! => distance to pole
  !      rb = 90.0 - 360./PI * atan(rp/AN)  ! => latitude
  !      if (rp >  1.0e-10)then
  !         rl = fi + 180.0/PI *atan2(dx,dy)
  !      else
  !         rl = fi
  !      endif
  !      if (rl <  glmin)   rl = rl + 360.0
  !      if (rl >  glmax)   rl = rl - 360.0
  !      gl_u(i,j)=rl                     !     longitude
  !      gb_u(i,j)=rb                     !     latitude
  !   end do ! i
  !end do ! j

  !  stop

  !relative orientation of grids

  !do i=1,IMAX+1
  !   do j=1,JMAX+1
  !      angle1(i,j)=datan2((i-xp),(yp-j))
  !      angle1_u(i,j)=datan2((i-xp_u),(yp-j))
  !      angle1_v(i,j)=datan2((i-xp),(yp_v-j))
  !   enddo
  !enddo
!  write(*,*)angle1(nint(xp),nint(yp))
!  write(*,*)angle1_u(nint(xp),nint(yp))
!  write(*,*)angle1_v(nint(xp),nint(yp))


!GLOBAL
  glmin = 0.0
  glmax = glmin + 360.0
lonstart=-179.0
latstart=-89.5
ddeg=1.0
ddeg_lon = 1.0
ddeg_lat = 1.0

!EMEP01
ddeg=0.1
ddeg_lon = 0.1
ddeg_lat = 0.1
lonstart=-30.0+0.5*ddeg_lon
latstart= 30.0+0.5*ddeg_lat


do i=1,imax
   lon(i)=lonstart+(i-1)*ddeg_lon
   lon_u(i)=lon(i)+0.5*ddeg_lon
enddo

do i=1,jmax
   lat(i)=latstart+(i-1)*ddeg_lat
   lat_v(i)=lat(i)+0.5*ddeg_lat
enddo


  lat_EC(1)=90.0
  lon_EC(1)=0.0

  GIMAX_X=3600
  GJMAX_X=1801
!  GIMAX=1800
!  GJMAX=901
  GIMAX=3600
  GJMAX=1801
  KMAX_MID=37
  write(*,*)GIMAX_X,GJMAX_X,GIMAX,GJMAX,KMAX_MID
  do i=2,gimax
     lon_EC(i)=lon_EC(1)+(i-1)*0.1
  enddo
  do i=2,gjmax
     lat_EC(i)=lat_EC(1)-(i-1)*0.1
  enddo

EC_lonstart=0.0
EC_latstart=90.0
EC_ddeg_lon=0.1
EC_ddeg_lat=-0.1! NB
EC_ddeg_lon_uv=0.1
EC_ddeg_lat_uv=-0.1! NB

do i=1,gimax_x
  lon_EC_X(i)=EC_lonstart+(i-1)*EC_ddeg_lon_uv
enddo
do i=1,gjmax_x
  lat_EC_X(i)=EC_latstart+(i-1)*EC_ddeg_lat_uv
enddo


  !stop

  !check assumed definition of i and j:
  delta=lon_EC(2)-lon_EC(1)
  if(abs(0.1-lon_EC(2)-lon_EC(1))<0.000001)delta=0.1!correct numerical error in long (written as real*4 in file)
  write(*,*)'longitude assumed: start=0 , step=',delta
  write(*,*)'latitude assumed: start=90 , step=',-delta

  do i=1,gimax
     !assumed lonitude for i:
     along=0.+delta*(i-1)
     !     write(*,*)I,along,lon_EC(i)
     if(abs(along-lon_EC(i))>0.01)then
        write(*,*)'WARNING:unexpected longitude definition',i,along,lon_EC(i)
        stop
     endif
  enddo
  !STOP
  do i=1,gjmax
     alat=90.-delta*(i-1)
     if(abs(alat-lat_EC(i))>0.01)then
        write(*,*)'WARNING:unexpected latitude definition',i,alat,lat_EC(i)
        stop
     endif
  enddo
  !  stop
  write(*,*)'longitude and latitude checked'



call findnearestneighbor_cyclic(lon,imax,lon_EC,GIMAX,conv_i1,conv_i2)


199 FORMAT(I4,10F12.5 )
do i=1,imax
   dist1=abs(lon(i)-lon_EC(conv_i1(i)))
   if(dist1>180.0)dist1=abs(dist1-360.0)
   if(dist1>180.0)dist1=abs(dist1-360.0)
   dist2=abs(lon(i)-lon_EC(conv_i2(i)))
   if(dist2>180.0)dist2=abs(dist2-360.0)
   if(dist2>180.0)dist2=abs(dist2-360.0)
   di1(i)=dist2/(dist1+dist2)
!   write(*,199)i,lon(i),lon_EC(conv_i1(i)),lon_EC(conv_i2(i)),di1(i),dist1,dist2,di1(i)*lon_EC(conv_i1(i))+(1.0-di1(i))*lon_EC(conv_i2(i))
enddo
!stop
call findnearestneighbor(lat,jmax,lat_EC,GJMAX,conv_j1,conv_j2)

do i=1,jmax
   dist1=abs(lat(i)-lat_EC(conv_j1(i)))
   dist2=abs(lat(i)-lat_EC(conv_j2(i)))
   dj1(i)=dist2/(dist1+dist2)
!   write(*,199)i,lat(i),lat_EC(conv_j1(i)),lat_EC(conv_j2(i)),dj1(i),dist1,dist2,dj1(i)*lat_EC(conv_j1(i))+(1.0-dj1(i))*lat_EC(conv_j2(i))
enddo
!stop

  if(abs(lon_EC(conv_i1(1))-lonstart)>1.00001*delta/2.and.abs(abs(lon_EC(conv_i1(1))-lonstart)-360.0)>1.00001*delta/2)then
     write(*,*)'error in longitude'
     write(*,*)conv_i1(1),lon_EC(conv_i1(1)),lonstart
     stop
  endif
  if(abs(lat_EC(conv_j1(1))-latstart)>delta/2)then
     write(*,*)'error in latitude'
     write(*,*)conv_j1(1),lat_EC(conv_j1(1)),latstart
     stop
  endif

call findnearestneighbor_cyclic(lon_u,imax,lon_EC_X,GIMAX_X,conv_i1_u,conv_i2_u)
do i=1,imax
   dist1=abs(lon_u(i)-lon_EC_X(conv_i1_u(i)))
   if(dist1>180.0)dist1=abs(dist1-360.0)
   if(dist1>180.0)dist1=abs(dist1-360.0)
   dist2=abs(lon_u(i)-lon_EC_X(conv_i2_u(i)))
   if(dist2>180.0)dist2=abs(dist2-360.0)
   if(dist2>180.0)dist2=abs(dist2-360.0)
   di1_u(i)=dist2/(dist1+dist2)
!   write(*,199)i,lon_u(i),lon_EC_X(conv_i1_u(i)),lon_EC_X(conv_i2_u(i)),di1_u(i),dist1,dist2,di1_u(i)*lon_EC_X(conv_i1_u(i))+(1.0-di1_u(i))*lon_EC_X(conv_i2_u(i))
enddo
!stop
call findnearestneighbor_cyclic(lon,imax,lon_EC_X,GIMAX_X,conv_i1_v,conv_i2_v)
do i=1,imax
   dist1=abs(lon(i)-lon_EC_X(conv_i1_v(i)))
   if(dist1>180.0)dist1=abs(dist1-360.0)
   if(dist1>180.0)dist1=abs(dist1-360.0)
   dist2=abs(lon(i)-lon_EC_X(conv_i2_v(i)))
   if(dist2>180.0)dist2=abs(dist2-360.0)
   if(dist2>180.0)dist2=abs(dist2-360.0)
   di1_v(i)=dist2/(dist1+dist2)
!   write(*,199)i,lon(i),lon_EC_X(conv_i1_v(i)),lon_EC_X(conv_i2_v(i)),di1_v(i),dist1,dist2,di1_v(i)*lon_EC_X(conv_i1_v(i))+(1.0-di1_v(i))*lon_EC_X(conv_i2_v(i))
enddo
!stop
call findnearestneighbor(lat,jmax,lat_EC_X,GJMAX_X,conv_j1_u,conv_j2_u)

do i=1,jmax
   dist1=abs(lat(i)-lat_EC_X(conv_j1_u(i)))
   dist2=abs(lat(i)-lat_EC_X(conv_j2_u(i)))
   dj1_u(i)=dist2/(dist1+dist2)
!   write(*,199)i,lat(i),lat_EC_X(conv_j1_u(i)),lat_EC_X(conv_j2_u(i)),dj1_u(i),dist1,dist2,dj1_u(i)*lat_EC_X(conv_j1_u(i))+(1.0-dj1_u(i))*lat_EC_X(conv_j2_u(i))
enddo
!stop
call findnearestneighbor(lat_v,jmax,lat_EC_X,GJMAX_X,conv_j1_v,conv_j2_v)

do i=1,jmax
   dist1=abs(lat_v(i)-lat_EC_X(conv_j1_v(i)))
   dist2=abs(lat_v(i)-lat_EC_X(conv_j2_v(i)))
   dj1_v(i)=dist2/(dist1+dist2)
!   write(*,199)i,lat_v(i),lat_EC_X(conv_j1_v(i)),lat_EC_X(conv_j2_v(i)),dj1_v(i),dist1,dist2,dj1_v(i)*lat_EC_X(conv_j1_v(i))+(1.0-dj1_v(i))*lat_EC_X(conv_j2_v(i))
enddo
!stop
  do i=1,imax
!give EC index from i index
     conv_i_u(i)=nint((lonstart-EC_lonstart)/EC_ddeg_lon_uv+(i-1+0.5)*ddeg_lon/EC_ddeg_lon_uv+1)
     if(conv_i_u(i)<1)conv_i_u(i)=conv_i_u(i)+GIMAX_X
     if(conv_i_u(i)>GIMAX_X)conv_i_u(i)=conv_i_u(i)-GIMAX_X
     if(conv_i_u(i)<1)then
        write(*,*)'conv_i_u(i)',i,conv_i_u(i)
        stop
     endif
     if(conv_i_u(i)>GIMAX_X)then
        write(*,*)'conv_i_u(i)',i,conv_i_u(i)
        stop
     endif
  enddo
  do j=1,jmax
!give EC index from j index
     conv_j_u(j)=nint((latstart-EC_latstart)/EC_ddeg_lat_uv+(j-1)*ddeg_lat/EC_ddeg_lat_uv+1)
     if(conv_j_u(j)<1)then
        write(*,*)'conv_j_u(i)',j,conv_j_u(j)
        write(*,*)(latstart-EC_latstart),EC_ddeg_lat_uv,(j-1)*ddeg_lat/EC_ddeg_lat_uv
        stop
     endif
     if(conv_j_u(j)>GJMAX_X)then
        write(*,*)'conv_j_u(i)',j,conv_j_u(j)
        stop
     endif
  enddo

  do i=1,imax
!give EC index from i index
     conv_i_v(i)=nint((lonstart-EC_lonstart)/EC_ddeg_lon_uv+(i-1)*ddeg_lon/EC_ddeg_lon_uv+1)
     if(conv_i_v(i)<1)conv_i_v(i)=conv_i_v(i)+GIMAX_X
     if(conv_i_v(i)>GIMAX_X)conv_i_v(i)=conv_i_v(i)-GIMAX_X
     if(conv_i_v(i)<1)stop
     if(conv_i_v(i)>GIMAX_X)stop
  enddo
  do j=1,jmax
!give EC index from j index
     conv_j_v(j)=nint((latstart-EC_latstart)/EC_ddeg_lat_uv+(j-1+0.5)*ddeg_lat/EC_ddeg_lat_uv+1)
     if(conv_j_v(j)<1)stop
     if(conv_j_v(j)>GJMAX_X)stop
  enddo



  ndate(4,1)=3
  ndate(4,2)=6
  ndate(4,3)=9
  ndate(4,4)=12
  ndate(4,5)=15
  ndate(4,6)=18
  ndate(4,7)=21
  ndate(4,8)=24
  do i=1,24
     ndate_h(4,i)=i
  enddo

  !stop


  !vertical level definition for ECMWF fields
  open(32,file='/global/home/mifapw/emep/prog/data/ECvert60.dat',form='formatted')

  !NB: Ah(i) and Bh(i) for half levels! (half level=level boundaries) 
  !NB: ph in hPa, Ah in Pa
  read(32,*)nn,Ah(1),Bh(1),ph(1)

  do i=2,61
     read(32,*)nn,Ah(i),Bh(i),ph(i),Press(i)

77   format(I2,5F)
     if(nn/=i-1)then
        write(*,*)'wrong format',i,nn,Ah(i),Bh(i),ph(i),Press(i)
        stop
     endif

  enddo
  Ah=Ah/100.0 !Pa->hPa

  !Sigma=B + (A+(B-1)*PT)/(PS0-PT) The second term is smallest
  !In a first good approximation we can assume PS0 constant=1013.25 hPa

  do i=1,60
     A(i)=0.5*(Ah(i+1)+Ah(i))
     B(i)=0.5*(Bh(i+1)+Bh(i))
     sigma_EC(i)=B(i) +(A(i)+(B(i)-1)*PT)/(PS0-PT)
     sigma_half_EC(i)=Bh(i) +(Ah(i)+(Bh(i)-1)*PT)/(PS0-PT)
     write(*,*)i,sigma_EC(i),sigma_half_EC(i)
  enddo
  sigma_half_EC(61)=1.

  do i=1,KMAX+1
     hyai(i)=ah(i+(60-KMAX))*100.0!hPa->Pa
     write(*,*)i,hyai(i)
  enddo
  do i=1,KMAX+1
     hybi(i)=bh(i+(60-KMAX))
  enddo

  deltaindex=60-KMAX_MID

  !NB should be defined according to actual integration time !! 
  delta_t=10800. !NB: only correct for 3hours

  iotyp=IOU_INST

  nstart=1 !first record read at 00:00 for accumulated fields
  Nhh=Nhh_8
!  if(me==nproc-1)then
     ! 2D (surface) fields
     allocate(varsEC(GIMAX,GJMAX,Nhh_24+1))
     allocate(vars1(IMAX,JMAX,Nhh_24+1))
     allocate(vars2(IMAX,JMAX,Nhh_24+1))
!     allocate(varsEC(GIMAX,GJMAX,Nhh+1))
!     allocate(vars1(IMAX,JMAX,Nhh+1))
!     allocate(vars2(IMAX,JMAX,Nhh+1))
     allocate(psEC(GIMAX,GJMAX,Nhh+1))
     allocate(t2mEC(GIMAX,GJMAX,Nhh+1))
     allocate(buffer_EC(GIMAX,GJMAX,Nhh+1))
  maxnland=20
     allocate(N_smi(GIMAX,GJMAX))
     allocate(I_smi(GIMAX,GJMAX,maxnland))
     allocate(J_smi(GIMAX,GJMAX,maxnland))
     
     ! 3D  fields
     allocate(varEC(GIMAX,GJMAX,KMAX_MID,Nhh+1))
     allocate(var1(IMAX,JMAX,KMAX,Nhh+1))
     allocate(var2(IMAX,JMAX,KMAX,Nhh+1))
!     allocate(varEC_u(GIMAX_X,GJMAX_X,KMAX_MID,1))
!     allocate(varEC_v(GIMAX_X,GJMAX_X,KMAX_MID,1))
!  endif
  write(*,*)GIMAX,GJMAX,KMAX_MID
!  stop

  if(me==0)then
     OPEN(24,file='InterCoeff',form='unformatted')
     
     i=0;j=0
     read(24,err=222,end=222)I,J,n
222  continue
     if(i/=gimax.or.j/=gjmax.or.n/=maxnland)then
        !needs only to be done once for every EC grid 
        write(*,*)'recalculating interpolation coefficents'
        write(filenameEC,57)trim(pathdir)//'mars',year,1,1,'-0000_smi1.nc'
        call GetCDF("swvl1",fileNameEC,varsEC(1,1,1),GIMAX,GJMAX,1,1,1)
        call  make_smi_interp(I_smi,J_smi,N_smi,varsEC(1,1,1),gimax,gjmax,maxnland)
        write(24)GIMAX,GJMAX,maxnland
        do j=1,GJMAX
           do i=1,GIMAX
              write(24)I_smi(i,j,1:maxnland),J_smi(i,j,1:maxnland),N_smi(i,j)
           enddo
        enddo
     else
        write(*,*)'found interpolation coefficents; reading from disk'
        do j=1,GJMAX
           do i=1,GIMAX
              read(24)I_smi(i,j,1:maxnland),J_smi(i,j,1:maxnland),N_smi(i,j)
           enddo
        enddo
     endif
     close(24)
     do j=1,GJMAX
        do i=1,GIMAX
           do n=1,N_smi(i,j)
              if(I_smi(i,j,n)<1)then
                 write(*,*)i,j,n,N_smi(i,j),I_smi(i,j,n),J_smi(i,j,n)
      stop
              endif
           enddo
        enddo
        enddo
  endif
  n=GIMAX*GJMAX*2
  CALL MPI_BCAST(I_smi,n*maxnland,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
  CALL MPI_BCAST(J_smi,n*maxnland,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
  CALL MPI_BCAST(N_smi,n,MPI_BYTE,0,MPI_COMM_WORLD,INFO) 
  CALL MPI_BARRIER(MPI_COMM_WORLD, INFO)

!  stop
  icount=0
  N3d=size(meteofields_3D)
  if(me==0)write(*,*)'making ',N3d ,' 3D fields '

  do yyyy=year-1,year
!  do yyyy=year,year
     yyyy1=yyyy
     if(4*(yyyy/4)==yyyy)then
        nmdays = (/31,29,31,30,31,30,31,31,30,31,30,31/) 
     else
        nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/) 
     endif
     do mm=1,12
        if(yyyy==year-1.and.mm<12)cycle
        mm1=mm
        do dd=1,nmdays(mm)
           if(yyyy==year-1.and.dd<31)cycle
!           if(mm==10.and.(dd<16))cycle
!           if(mm==10.and.(dd<16.or.dd>19))cycle
           dd1=dd+1 !because the "prevision" can be for the day after 
           if(dd1>nmdays(mm))then
              dd1=1
              mm1=mm+1
              if(mm1==13)then
                 mm1=1
                 yyyy1=yyyy+1
              endif
           endif
           ndate(1,:)=yyyy1
           ndate(2,:)=mm1
           ndate(3,:)=dd1
           ndate_h(1,:)=yyyy1
           ndate_h(2,:)=mm1
           ndate_h(3,:)=dd1

!first part, hours 3 6 9 12


        if(me==mod(icount,nproc))then

56         FORMAT(a5,i4.4,i2.2,i2.2,a3)
57         FORMAT(a,i4.4,i2.2,i2.2,a)
           write(filename,56)'meteo',yyyy1,mm1,dd1,'.nc'
           fileName=trim(pathdirout)//fileName

!           call CreatenetCDFfilePS(fileName,IMAX,JMAX,KMAX,sigma_EMEP,GRIDWIDTH,xp,yp,fi,ref_latitude,hyai,hybi,xm)
           call CreatenetCDFfile(fileName,IMAX,JMAX,KMAX,lonstart,latstart,ddeg_lon,ddeg_lat,hyai,hybi,GRIDWIDTH)
           write(*,*)trim(fileName),' Created'
           write(filename_h,57)'meteo_h',yyyy1,mm1,dd1,'.nc'
           fileName_h=trim(pathdirout)//fileName_h
          call CreatenetCDFfile(fileName_h,IMAX,JMAX,KMAX,lonstart,latstart,ddeg_lon,ddeg_lat,hyai,hybi,GRIDWIDTH)
!           call CreatenetCDFfilePS(fileName_h,IMAX,JMAX,KMAX,sigma_EMEP,GRIDWIDTH,xp,yp,fi,ref_latitude)
!           write(*,*)trim(fileName_h),' Created'

           t2msaved=0
           N2d=22
           do i2d=1,N2d
              sname_req=meteofields_2D(i2d)%EC_name
              validity='instantaneous' !default
! -1200 file contains the data for the day after at 3 6 9 and 12 hours 
              write(filenameEC,57)trim(pathdir)//'mars',yyyy,mm,dd,'-1200'//meteofields_2D(i2d)%ending
              if(meteofields_2D(i2d)%hourly)then
                 Nhh=Nhh_24
              else
                 Nhh=Nhh_8
              endif
              Nhh_half=Nhh/2
              nstart=1
              nfetch=Nhh_half
              !if(meteofields_2D(i2d)%ending=='_surf.nc')then
              !   nstart=2
              !   if(meteofields_2D(i2d)%treatment=="accumulated".or.meteofields_2D(i2d)%treatment=="vector_acc")then
              !      nstart=1
              !      nfetch=Nhh_half+1
              !  endif
              !endif
              if(meteofields_2D(i2d)%treatment=="accumulated".or.meteofields_2D(i2d)%treatment=="vector_acc")then
                 nfetch=Nhh_half+1
              endif

              call GetCDF(sname_req,fileNameEC,varsEC(1,1,1),GIMAX,GJMAX,1,nstart,nfetch)
!Test if the field is constant in the last record as a crude check for corrupt files
              if(maxval(varsEC(:,:,nfetch))==minval(varsEC(:,:,nfetch)))then
                 write(*,*)'WARNING: CONSTANT VALUE ',trim(sname_req)//' ',trim(filenameEC),nstart+nfetch-1
                 stop
              endif
              fscale=meteofields_2D(i2d)%fscale
              do n=1,nfetch
                 varsEC(:,:,n)=varsEC(:,:,n) *fscale
              enddo

              if(meteofields_2D(i2d)%treatment=="accumulated".or.meteofields_2D(i2d)%treatment=="vector_acc")then
!                 if(Nhh==Nhh_24)then
!simulate a 3 hourly metdata set
!                    do n=0,3
!                       varsEC(:,:,3*n+1)=(varsEC(:,:,3*(n+1)+1)-varsEC(:,:,3*n+1))/3.0
!                       varsEC(:,:,3*n+2)=varsEC(:,:,3*n+1)
!                       varsEC(:,:,3*n+3)=varsEC(:,:,3*n+1)
!                    enddo                    
!                 else
                    !convert into differences
                    do n=1,nfetch-1
                       varsEC(:,:,n)=(varsEC(:,:,n+1)-varsEC(:,:,n))
                    enddo
!                 endif
              endif

! -0000 file contains the data for the same day at 15 18 21 and 24 hours 
              write(filenameEC,57)trim(pathdir)//'mars',yyyy1,mm1,dd1,'-0000'//meteofields_2D(i2d)%ending
              nstart=1
              nfetch=Nhh_half
              !if(meteofields_2D(i2d)%ending=='_surf.nc')then
              !   nstart=2
              !   if(meteofields_2D(i2d)%treatment=="accumulated".or.meteofields_2D(i2d)%treatment=="vector_acc")then
              !      nstart=1
              !      nfetch=Nhh_half+1
              !  endif
              !endif
              if(meteofields_2D(i2d)%treatment=="accumulated".or.meteofields_2D(i2d)%treatment=="vector_acc")then
                 nfetch=nfetch+1
              endif

!fill second half of varsEC
              call GetCDF(sname_req,fileNameEC,varsEC(1,1,Nhh_half+1),GIMAX,GJMAX,1,nstart,nfetch)

 !Test if the field is constant in the last record as a crude check for corrupt files
              if(maxval(varsEC(:,:,Nhh_half+nfetch))==minval(varsEC(:,:,Nhh_half+nfetch)))then
                 write(*,*)'WARNING: CONSTANT VALUE 2nd ',trim(sname_req)//' ',trim(filenameEC),nstart+nfetch-1,Nhh_half
                 stop
              endif
              fscale=meteofields_2D(i2d)%fscale
              do n=Nhh_half+1,Nhh_half+nfetch
                 varsEC(:,:,n)=varsEC(:,:,n) *fscale
              enddo

             if(meteofields_2D(i2d)%treatment=="accumulated".or.meteofields_2D(i2d)%treatment=="vector_acc")then
!                 if(Nhh==Nhh_24)then
!simulate a 3 hourly metdata set
!                    do n=4,7
!                       varsEC(:,:,3*n+1)=(varsEC(:,:,3*(n+1)+1)-varsEC(:,:,3*n+1))/3.0
!                       varsEC(:,:,3*n+2)=varsEC(:,:,3*n+1)
!                       varsEC(:,:,3*n+3)=varsEC(:,:,3*n+1)
!                    enddo                    
!                 else
                 !convert into differences
                 do n=Nhh_half+nstart,Nhh_half+nfetch-1
                    varsEC(:,:,n)=(varsEC(:,:,n+1)-varsEC(:,:,n))
                 enddo
!                 endif
                 validity='averaged'
              endif

              if(meteofields_2D(i2d)%treatment=="log")then
                 varsEC=exp(varsEC)*0.01 !Pa->hPa
                 psEC(:,:,1:Nhh)=varsEC(:,:,1:Nhh) !save for future use
              endif
              if(meteofields_2D(i2d)%emep_name=="temperature_2m")then
                 t2mEC(:,:,1:Nhh)=varsEC(:,:,1:Nhh) !save for future use
                 t2msaved=1
              endif
              if(meteofields_2D(i2d)%treatment=="dewpoint")then
                 if(t2msaved/=1)stop
                 varsEC(:,:,1:Nhh)=min(100.0,100*exp(&
                 17.625*(varsEC(:,:,1:Nhh)-T0)/(243.04-T0+varsEC(:,:,1:Nhh))-&
                 17.625*(t2mEC(:,:,1:Nhh)-T0)/(243.04-T0+t2mEC(:,:,1:Nhh)) ))
              endif

              if(meteofields_2D(i2d)%treatment=="SetSea")then
!For undefined sea areas and some islands, define some reasonable values
                 buffer_EC=varsEC!buffer_EC used as help variable so as to interpolate only from the uninterpolated
                 
                 do i=1,gimax
                    do j=1,gjmax
                       sum=0.0
                       if(N_smi(i,j)>0.0)then
                          do n=1,N_smi(i,j)
                             sum(1:Nhh)=sum(1:Nhh)+buffer_EC(I_smi(i,j,n),J_smi(i,j,n),1:Nhh)
                          enddo
                          sum=sum/N_smi(i,j)
                          varsEC(i,j,1:Nhh)=sum(1:Nhh)
                       endif
                    enddo
                 enddo
                 
              endif


                 do n=1,Nhh
                    do i=1,imax
                       i_EC1=conv_i1(i)
                       i_EC2=conv_i2(i)
                       do j=1,jmax

                          j_EC1=conv_j1(j)
                          j_EC2=conv_j2(j)
                          d00 = di1(i)      *dj1(j)
                          d10 = (1.0-di1(i))*dj1(j)
                          d01 = di1(i)      *(1.0-dj1(j))
                          d11 = (1.0-di1(i))*(1.0-dj1(j))

                          vars1(i,j,n)=&
                                varsEC(i_EC1,j_EC1,n)*d00&
                               +varsEC(i_EC2,j_EC1,n)*d10&
                               +varsEC(i_EC1,j_EC2,n)*d01&
                               +varsEC(i_EC2,j_EC2,n)*d11
                      enddo
                    enddo
                 enddo

              if(meteofields_2D(i2d)%treatment=="vector_acc".and.vector==.false.)then
                 vars2=vars1
                 vector=.true.
                 cycle
              endif
              if(vector==.true.)then
                 vars1=sqrt(vars2**2+vars1**2)
              endif
              vector=.false.
 
              namefield=meteofields_2D(i2d)%emep_name
              units=meteofields_2D(i2d)%unit
              if(meteofields_2D(i2d)%hourly)then
 
                call writemeteofield(fileName_h,namefield,vars1 ,IMAX,JMAX,KMAX1,Nhh,ndate_h,validity,units)
              else
                call writemeteofield(fileName,namefield,vars1 ,IMAX,JMAX,KMAX1,Nhh,ndate,validity,units)
             endif
              write(*,*)trim(namefield),'written'

           enddo

           do i3d=1,N3d
              sname_req=meteofields_3D(i3d)%EC_name
              nstart=1
              Nhh=Nhh_8
              Nhh_half=Nhh/2
              nfetch=Nhh_half
              validity='instantaneous'
              if(meteofields_3D(i3d)%treatment(1:3)=='acc')then
                 nfetch=nfetch+1
                 validity='averaged'
              endif
!              if(meteofields_3D(i3d)%treatment(1:3)=='acc')then
!                 nfetch=nfetch+1
!                 validity='averaged'
!              endif
              if(meteofields_3D(i3d)%emep_name/='u_wind'.and.meteofields_3D(i3d)%emep_name/='v_wind')then
!                 if(allocated(varEC_u))then
!                    write(*,*)'wind finished, deallocating'
!                    deallocate(varEC_u)
!                 endif
                 write(filenameEC,57)trim(pathdir)//'mars',yyyy,mm,dd,'-1200'//meteofields_3D(i3d)%ending
                 call GetCDF(sname_req,fileNameEC,varEC(1,1,1,1),GIMAX,GJMAX,KMAX_MID,nstart,nfetch)
!Test if the field is constant in the last record as a crude check for corrupt files
                 ktest=KMAX_MID                
                 if(meteofields_3D(i3d)%emep_name(1:4)=='conv')ktest=KMAX_MID-4!top and bottom may be constant zero                 
                 if(maxval(varEC(:,:,ktest,nfetch))==minval(varEC(:,:,ktest,nfetch)))then
                    write(*,*)'WARNING: CONSTANT VALUE ',trim(sname_req)//' ',trim(filenameEC)
                    stop
                 endif
                 fscale=meteofields_3D(i3d)%fscale
                 do n=1,nfetch
                    varEC(:,:,:,n)=varEC(:,:,:,n)*fscale
                 enddo
                 if(meteofields_3D(i3d)%treatment(1:3)=='acc')then
                    !convert into differences
                    do n=1,nfetch-1
                       varEC(:,:,:,n)=(varEC(:,:,:,n+1)-varEC(:,:,:,n))
                    enddo
                 endif
                 write(filenameEC,57)trim(pathdir)//'mars',yyyy1,mm1,dd1,'-0000'//meteofields_3D(i3d)%ending
                 call GetCDF(sname_req,fileNameEC,varEC(1,1,1,Nhh_half+1),GIMAX,GJMAX,KMAX_MID,nstart,nfetch)
!Test if the field is constant in the last record as a crude check for corrupt files
                 ktest=KMAX_MID                
                 if(meteofields_3D(i3d)%emep_name(1:4)=='conv')ktest=KMAX_MID-4!top and bottom may be constant zero                 
                 if(maxval(varEC(:,:,ktest,Nhh_half+nfetch))==minval(varEC(:,:,ktest,Nhh_half+nfetch)))then
                    write(*,*)'WARNING: CONSTANT VALUE ',trim(sname_req)//' ',trim(filenameEC)
                    stop
                 endif
                fscale=meteofields_3D(i3d)%fscale
                 do n=Nhh_half+1,Nhh_half+nfetch
                    varEC(:,:,:,n)=varEC(:,:,:,n)*fscale
                 enddo
                 if(meteofields_3D(i3d)%treatment(1:3)=='acc')then
                    !convert into differences
                    do n=Nhh_half+1,Nhh_half+nfetch-1
                       varEC(:,:,:,n)=(varEC(:,:,:,n+1)-varEC(:,:,:,n))
                    enddo
                 endif

                    if(meteofields_3D(i3d)%treatment=="absolute")then
                       !convert from absolute to potential temperature
                       do n=1,Nhh
                           do k=1,KMAX_MID
                             varEC(:,:,k,n) = varEC(:,:,k,n)*exp(-XKAP*log((A(k+deltaindex) + B(k+deltaindex)*psEC(:,:,n))*1.e-3))
                          enddo
                       enddo
                    endif
 
                    do n=1,Nhh
                       do k=1,KMAX
                          k1_EC=k!EC_index(k)-deltaindex
!                          k2_EC=k1_EC-1
                          do j=1,jmax
                             j_EC1=conv_j1(j)
                             j_EC2=conv_j2(j)
                             do i=1,imax
                                i_EC1=conv_i1(i)
                                i_EC2=conv_i2(i)
                                d00 = di1(i)      *dj1(j)
                                d10 = (1.0-di1(i))*dj1(j)
                                d01 = di1(i)      *(1.0-dj1(j))
                                d11 = (1.0-di1(i))*(1.0-dj1(j))
                                
                                var1(i,j,k,n)=&
                                     varEC(i_EC1,j_EC1,k1_EC,n)*d00+&
                                      varEC(i_EC2,j_EC1,k1_EC,n)*d10+&
                                      varEC(i_EC1,j_EC2,k1_EC,n)*d01+&
                                      varEC(i_EC2,j_EC2,k1_EC,n)*d11

                             enddo
                          enddo
 
                       enddo
                    enddo

                 if(meteofields_3D(i3d)%emep_name=='total_convective_precipitation_profile')then
                    do n=1,Nhh
                       do k=1,KMAX
                          do j=1,JMAX
                             do i=1,IMAX
                                var2(i,j,k,n)= var1(i,j,k,n)
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
                 namefield=meteofields_3D(i3d)%emep_name

                 if(meteofields_3D(i3d)%emep_name=='total_stratiform_precipitation_profile')then
                    var1= var1+var2
                    namefield='precipitation'
                 endif
                 
                 units=meteofields_3D(i3d)%unit
                 call writemeteofield(fileName,namefield,var1 ,IMAX,JMAX,KMAX,Nhh,ndate,validity,units)
                 write(*,*)trim(namefield),'written'
              

              else
                 !need to take one record at a time because of memory limitations:
                 !(u_og_v*X*Y*K*t*8bytes=)2*3600*1801*37*8*8=30GB
                 write(*,*)'reading wind',GIMAX_X,GJMAX_X,KMAX_MID
!                 if(.not.allocated(varEC_u))then
!                    write(*,*)'wind , reallocating'
!                    allocate(varEC_u(GIMAX_X,GJMAX_X,KMAX_MID,1))
!                 endif

                 Nhh=Nhh_8
                 Nhh_half=Nhh/2
 !                if(meteofields_3D(i3d)%emep_name=='v_wind')cycle!done together with u_wind
                 do nrecord=1,Nhh
                    if(nrecord<=Nhh_half)then
                       write(filenameEC,57)trim(pathdir)//'mars',yyyy,mm,dd,'-1200'//meteofields_3D(i3d)%ending
                       write(filenameEC_v,57)trim(pathdir)//'mars',yyyy,mm,dd,'-1200'//trim(meteofields_3Dending_v)
                       nrecord_read=nrecord
                    else
                       write(filenameEC,57)trim(pathdir)//'mars',yyyy1,mm1,dd1,'-0000'//meteofields_3D(i3d)%ending
                       write(filenameEC_v,57)trim(pathdir)//'mars',yyyy1,mm1,dd1,'-0000'//trim(meteofields_3Dending_v)
                       nrecord_read=nrecord-Nhh_half
                    endif
                    if(meteofields_3D(i3d)%emep_name=='u_wind'.or.meteofields_3D(i3d)%emep_name=='v_wind')then
                       call GetCDF(sname_req,fileNameEC,varEC(1,1,1,1),GIMAX_X,GJMAX_X,KMAX_MID,nrecord_read,1)
                       !Test if the field is constant in the last record as a crude check for corrupt files
                       ktest=KMAX_MID                
                       if(maxval(varEC(:,:,ktest,1))==minval(varEC(:,:,ktest,1)))then
                          write(*,*)'WARNING: CONSTANT VALUE U ',trim(sname_req)//' ',trim(filenameEC)
                          stop
                       endif
!                       call GetCDF('v',fileNameEC_v,varEC_v(1,1,1,1),GIMAX_X,GJMAX_X,KMAX_MID,nrecord_read,1)
!                       !Test if the field is constant in the last record as a crude check for corrupt files
!                       ktest=KMAX_MID                
!                       if(maxval(varEC_v(:,:,ktest,1))==minval(varEC_v(:,:,ktest,1)))then
!                          write(*,*)'WARNING: CONSTANT VALUE V ',trim(sname_req)//' ',trim(filenameEC)
!                          stop
!                       endif
                    else
                       write(*,*)'uv error'
                       stop
                    endif
                    n=1
                    do k=1,KMAX
                       k1_EC=k!EC_index(k)-deltaindex
!                       k2_EC=k1_EC-1
!                       K_EC=k1_EC
                       !write(*,*)n,k,K_EC,EC_index(k)
                       if(EC2EMEP1(k)<EC2EMEP2(k))K_EC=k2_EC
                       do i=1,imax
                          do j=1,jmax
                             !u and v in the nearest defined position (0 order interpolation):
                             !        write(*,*)i,j,nint(gl(i,j)),nint(gb(i,j))
                             !u and v at gridcell boundaries
                             if(meteofields_3D(i3d)%emep_name=='u_wind')then
                                i_EC1=conv_i1_u(i)
                                i_EC2=conv_i2_u(i)
                                j_EC1=conv_j1_u(j)
                                j_EC2=conv_j2_u(j)

                                d00 = di1_u(i)      *dj1_u(j)
                                d10 = (1.0-di1_u(i))*dj1_u(j)
                                d01 = di1_u(i)      *(1.0-dj1_u(j))
                                d11 = (1.0-di1_u(i))*(1.0-dj1_u(j))

                             else
                                i_EC1=conv_i1_v(i)
                                i_EC2=conv_i2_v(i)
                                j_EC1=conv_j1_v(j)
                                j_EC2=conv_j2_v(j)

                                d00 = di1_v(i)      *dj1_v(j)
                                d10 = (1.0-di1_v(i))*dj1_v(j)
                                d01 = di1_v(i)      *(1.0-dj1_v(j))
                                d11 = (1.0-di1_v(i))*(1.0-dj1_v(j))

                             endif
                                var1(i,j,k,nrecord)=&
                                     (varEC(i_EC1,j_EC1,k1_EC,n)*d00+&
                                      varEC(i_EC2,j_EC1,k1_EC,n)*d10+&
                                      varEC(i_EC1,j_EC2,k1_EC,n)*d01+&
                                      varEC(i_EC2,j_EC2,k1_EC,n)*d11)

                             !rotate to EMEP grid
!                             var1(i,j,k,nrecord)=u_ll*dcos(angle1_u(i,j))-v_ll*dsin(angle1_u(i,j))

                             !        write(*,*) EMEP_3D(i,j,k) ,angle1(i,j),EC2EMEP1(k),var0(i_EC,j_EC,k1_EC,1)
                          enddo
                       enddo
                    enddo
                 enddo
                 namefield=meteofields_3D(i3d)%emep_name
                 validity='instantaneous'
                 units=meteofields_3D(i3d)%unit
                 call writemeteofield(fileName,namefield, var1,IMAX,JMAX,KMAX,Nhh,ndate,validity,units)
                 write(*,*)trim(namefield),' finished'

              endif



           enddo

        endif

        icount=icount+1

        enddo
     enddo
  enddo


  CALL MPI_FINALIZE(INFO)


end Program callreadCDF

!module ReadCDF_ml

! implicit none

!public  :: GetCDFhour_1month
!public  :: GetCDF,GetCDFinfo

!contains
!---------------------------------------------------------------------

subroutine GetCDF(varname,fileName,var,I_MAX,J_MAX,K_MAX,nstart,nfetch)
  !
  ! open and reads CDF file
  !
!  varname: field to be fetched (input)
!  fileName:fileName (input)
!  var: array to store field (output)
!  I_MAX,J_MAX,K_MAX, dimensions of "var"(input)
!  nstart: first field read
!  nfetch: number of fields read
  !
  ! The nf90 are functions which return 0 if no error occur.
  ! check is only a subroutine which check wether the function returns zero
  !
  !
!  use typeSizes
  use netcdf
  implicit none

  character (len=*),intent(in) :: fileName 

  character (len = *),intent(in) ::varname
  integer, intent(in) :: nstart,I_MAX,J_MAX,K_MAX
  integer, intent(inout) ::  nfetch
!  real, dimension(I_MAX*J_MAX*K_MAX*20),intent(inout) :: var
  real, dimension(*),intent(inout) :: var
!  real, dimension(:,:,:,:),intent(out) :: var4D
!  real, dimension(132,111,Nrec),intent(out) :: var


  integer :: GIMAX,GJMAX,KMAX_MID,nrecords,period
  integer :: varID,status,ndims,alloc_err
  integer :: n,KMAX,Nrec,ijn,ijkn
  integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,iVarID,jVarID,kVarID,i,j,k
  integer :: var_date(9000),ndate(4)
!  real (kind = FourByteReal), allocatable,dimension(:,:,:,:)  :: values
  real , allocatable,dimension(:,:,:,:)  :: values
  real ::depsum,add_offset,scale_factor
  character*100::attribute,attribute2

!  Nrec=size(var,3)
  Nrec=nfetch

  print *,'  reading ',trim(fileName)
  !open an existing netcdf dataset
  call check(nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID))

  !get global attributes

  !example:
!  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_hour", attribute ))
!  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_date", attribute2 ))
!  print *,'file last modified (yyyymmdd hhmmss.sss) ',attribute2,' ',attribute

  !test if the variable is defined and get varID:
  status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)

  if(status == nf90_noerr) then     
     print *, 'variable exists: ',trim(varname)
  else
    if(trim(varname)=='sstk')then
        status = nf90_inq_varid(ncid = ncFileID, name = 'sst', varID = VarID)
        if(status == nf90_noerr) then     
           print *, 'variable exists: ','sst'
        else
           print *, 'variable does not exist: ',trim(varname),' and sst',nf90_strerror(status)
           stop
           return
        endif
     elseif(trim(varname)=='v2t')then
        status = nf90_inq_varid(ncid = ncFileID, name = 't2m', varID = VarID)
        if(status == nf90_noerr) then     
           print *, 'variable exists: ','t2m'
        else
           print *, 'variable does not exist: ',trim(varname),' and t2m ',nf90_strerror(status)
           stop
           return
        endif
     elseif(trim(varname)=='v10u')then
        status = nf90_inq_varid(ncid = ncFileID, name = 'u10', varID = VarID)
        if(status == nf90_noerr) then     
           print *, 'variable exists: ','u10'
        else
           print *, 'variable does not exist: ',trim(varname),' and u10 ',nf90_strerror(status)
           stop
           return
        endif
     elseif(trim(varname)=='v10v')then
        status = nf90_inq_varid(ncid = ncFileID, name = 'v10', varID = VarID)
        if(status == nf90_noerr) then     
           print *, 'variable exists: ','v10'
        else
           print *, 'variable does not exist: ',trim(varname),' and v10 ',nf90_strerror(status)
           stop
           return
        endif
     elseif(trim(varname)=='v2d')then
        status = nf90_inq_varid(ncid = ncFileID, name = 'd2m', varID = VarID)
        if(status == nf90_noerr) then     
           print *, 'variable exists: ','d2m'
        else
           print *, 'variable does not exist: ',trim(varname),' and d2m ',nf90_strerror(status)
           stop
           return
        endif
     else
        status = nf90_inq_varid(ncid = ncFileID, name = 'p'//trim(varname),varID = VarID)
        if(status == nf90_noerr) then     
           print *, 'variable exists: ','p'//trim(varname)
        else
           status = nf90_inq_varid(ncid = ncFileID, name = 'v~', varID = VarID)
           if(status == nf90_noerr) then     
              print *, 'variable exists: ','v~'
           else
              print *, 'variable does not exist: ',trim(varname),' and v~ or p108.128 ',nf90_strerror(status)
              stop
              return
           endif
        endif
     endif
  endif

  !get dimensions id
  call check(nf90_inq_dimid(ncid = ncFileID, name = "longitude", dimID = idimID))
  call check(nf90_inq_dimid(ncid = ncFileID, name = "latitude", dimID = jdimID))
  if(K_MAX>1)then
     status = nf90_inq_dimid(ncid = ncFileID, name = "levelist", dimID = kdimID)
     if(status /= nf90_noerr) then
        call check(nf90_inq_dimid(ncid = ncFileID, name = "level", dimID = kdimID))
     endif
  endif
 call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

  !get dimensions length
  call check(nf90_inquire_dimension(ncid=ncFileID, dimID=idimID,  len=GIMAX))
  call check(nf90_inquire_dimension(ncid=ncFileID, dimID=jdimID,  len=GJMAX))
  KMAX_MID=1
 if(K_MAX>1)then
  call check(nf90_inquire_dimension(ncid=ncFileID, dimID=kdimID,  len=KMAX_MID))
endif
  call check(nf90_inquire_dimension(ncid=ncFileID, dimID=timedimID,  len=nrecords))

  print *, 'dimensions ',GIMAX,GJMAX,KMAX_MID,nrecords
!  if(GIMAX>size(var,1).or.GJMAX>size(var,2))then
!     write(*,*)'buffer too small',GIMAX,size(var,1),GJMAX,size(var,2)
!     stop
!  endif

  !get variable info
  call check(nf90_inquire_variable(ncFileID, varID, ndims=ndims))
       print *, 'dimensions ',ndims

  if(nstart+nfetch-1>nrecords)then
     write(*,*)'WARNING: did not find all data'
     nfetch=nrecords-nstart+1
     if(nfetch<=0)stop
  endif
  if(nfetch>Nrec)then
     write(*,*)'buffer too small. Increase last dimension'
!     stop
  endif
  if(ndims==3)then
     kmax=1
     !allocate a 2D array 
     allocate(values(GIMAX,GJMAX,nfetch,1), stat=alloc_err)
     if ( alloc_err /= 0 ) then
        print *, 'alloc failed in ReadCDF_ml: ',alloc_err,ndims
        stop
     endif
  elseif(ndims==4)then
     kmax=KMAX_MID
     !allocate a 3D array 
     allocate(values(GIMAX,GJMAX,KMAX_MID,nfetch), stat=alloc_err)
     if ( alloc_err /= 0 ) then
        print *, 'alloc failed in ReadCDF_ml: ',alloc_err,ndims
        stop
     endif

  else
     print *, 'unexpected number of dimensions: ',ndims
     stop
  endif

  !get variable attributes
  !example:
  attribute=''
  if(varname/='108.162')then
  call check(nf90_get_att(ncFileID, VarID, "long_name", attribute))
       print *,'long_name ',attribute
    endif
    
  call check(nf90_get_att(ncFileID, VarID, "add_offset", add_offset))
  call check(nf90_get_att(ncFileID, VarID, "scale_factor", scale_factor))


  !get time variable
  call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeDimID))
  call check(nf90_get_var(ncFileID, timeDimID, var_date,start=(/ nstart /),count=(/ nfetch /)))

  !get variable
  if(ndims==3)then
     !      call check(nf90_get_var(ncFileID, VarID, values,start=(/ 1, 1, nstart /),count=(/ 1, 1, nfetch /)))
     call check(nf90_get_var(ncFileID, VarID, values,start=(/ 1, 1, nstart /),count=(/ GIMAX,GJMAX,nfetch /)))
  elseif(ndims==4)then
     call check(nf90_get_var(ncFileID, VarID, values,start=(/1,1,1,nstart/),count=(/GIMAX,GJMAX,KMAX_MID,nfetch /)))
  endif
!  if(Nrec<nrecords)then
!     write(*,*)'Reading record',nstart,' to ',nstart+nfetch-1
!  endif
!  write(*,*)'date start '
!  call datefromsecondssince1970(ndate,var_date(1))
!  write(*,*)'date end '
!  call datefromsecondssince1970(ndate,var_date(nfetch))
!  period=(var_date(2)-var_date(1))/3600.

  if(ndims==3)then
!     write(*,*)'implemented'
!     stop
     do n=1,nfetch
        do j=1,GJMAX
           do i=1,GIMAX
              ijn=i+(j-1)*I_MAX+(n-1)*I_MAX*J_MAX  
              var(ijn)=values(i,j,n,1)*scale_factor+add_offset
           enddo
        enddo
!        !     if(n<10)write(*,*)n,values(1,1,n,1)
     enddo
  else
!     write(*,*)nfetch,KMAX_MID,I_MAX,J_MAX,K_MAX,scale_factor,add_offset
     do n=1,nfetch
        do k=1,KMAX_MID
        do j=1,GJMAX
           do i=1,GIMAX
              ijkn=i+(j-1)*I_MAX+(k-1)*I_MAX*J_MAX+(n-1)*I_MAX*J_MAX*K_MAX
              var(ijkn)=values(i,j,k,n)*scale_factor+add_offset
!              write(*,*)ijkn,i,j,k,n,ijkn,values(i,j,k,n)
           enddo
        enddo
        enddo
        !     if(n<10)write(*,*)n,values(1,1,n,1)
     enddo
!     stop
  endif
  call check(nf90_close(ncFileID))
  deallocate(values)

end subroutine GetCDF
!_______________________________________________________________________

subroutine check(status)

!  use typeSizes
  use netcdf
  implicit none
  integer, intent ( in) :: status

  if(status /= nf90_noerr) then 
     print *, trim(nf90_strerror(status))
     print *, "error in NetCDF_ml"
     stop
  end if
end subroutine check

subroutine datefromsecondssince1970(ndate,nseconds,printdate)
  !calculate date from seconds that have passed since the start of the year 1970

!  use Dates_ml, only : nmdays
  implicit none

  integer, intent(out) :: ndate(4)
  integer, intent(in) :: nseconds
  integer, intent(in) :: printdate

  integer :: n,nday,nmdays(12),nmdays2(13)
  nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/) 

  nmdays2(1:12)=nmdays
  nmdays2(13)=0
  ndate(1)=1969
  n=0
  do while(n<=nseconds)
     n=n+24*3600*365
     ndate(1)=ndate(1)+1
     if(mod(ndate(1),4)==0)n=n+24*3600
  enddo
  n=n-24*3600*365
  if(mod(ndate(1),4)==0)n=n-24*3600
  if(mod(ndate(1),4)==0)nmdays2(2)=29
  ndate(2)=0
  do while(n<=nseconds)
     ndate(2)=ndate(2)+1
     n=n+24*3600*nmdays2(ndate(2))
  enddo
  n=n-24*3600*nmdays2(ndate(2))
  ndate(3)=0
  do while(n<=nseconds)
     ndate(3)=ndate(3)+1
     n=n+24*3600
  enddo
  n=n-24*3600
  ndate(4)=-1
  do while(n<=nseconds)
     ndate(4)=ndate(4)+1
     n=n+3600
  enddo
  n=n-3600
  !    ndate(5)=nseconds-n
  if(printdate>0)then
  write(*,*)'year: ',ndate(1),', month: ',ndate(2),', day: ',ndate(3),', hour: ',ndate(4),', seconds: ',nseconds-n
  endif
end subroutine datefromsecondssince1970
subroutine datefromhourssince1900(ndate,nhours,printdate)
  !calculate date from hours that have passed since the start of the year 1900
!NB: 1900 is not a leap year

!  use Dates_ml, only : nmdays
  implicit none

  integer, intent(out) :: ndate(4)
  integer, intent(in) :: nhours
  integer, intent(in) :: printdate

  integer :: n,nday,nmdays(12),nmdays2(13)
  nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/) 

  nmdays2(1:12)=nmdays
  nmdays2(13)=0
  ndate(1)=1899
  n=0
  do while(n<=nhours)
     n=n+24*365
     ndate(1)=ndate(1)+1
     if(mod(ndate(1),4)==0)n=n+24
     if(ndate(1)==1900)n=n-24
  enddo
  n=n-24*365
  if(mod(ndate(1),4)==0)n=n-24
  if(ndate(1)==1900)n=n+24
  if(mod(ndate(1),4)==0)nmdays2(2)=29
  if(ndate(1)==1900)nmdays2(2)=28
!  write(*,*)ndate(1),n
  ndate(2)=0
  do while(n<=nhours)
     ndate(2)=ndate(2)+1
     n=n+24*nmdays2(ndate(2))
  enddo
  n=n-24*nmdays2(ndate(2))
!  write(*,*)ndate(2),n
  ndate(3)=0
  do while(n<=nhours)
     ndate(3)=ndate(3)+1
     n=n+24
  enddo
  n=n-24
!  write(*,*)ndate(3),n
  ndate(4)=nhours-n
!  do while(n<=nseconds)
!     ndate(4)=ndate(4)+1
!     n=n+3600
!  enddo
!  n=n-3600
  !    ndate(5)=nseconds-n
  if(printdate>0)then
  write(*,*)'year: ',ndate(1),', month: ',ndate(2),', day: ',ndate(3),', hour: ',ndate(4)
  endif
end subroutine datefromhourssince1900

subroutine GetCDFinfo(fileName,timev,GIMAX,GJMAX,KMAX_MID,nrecords,period,Nrec,long,lat,level)
!
! open and reads CDF file
!
! The nf90 are functions which return 0 if no error occur.
! check is only a subroutine which check wether the function returns zero
!
!
!  use typeSizes
  use netcdf
  implicit none

  character (len=*),intent(in) :: fileName 

 integer, intent(out) :: GIMAX,GJMAX,KMAX_MID,nrecords
 integer, intent(in) :: Nrec
  integer, dimension(Nrec),intent(out) :: timev
  real, intent(out) ::period
  real, dimension(*), intent(out) ::long,lat,level
  integer :: varID,status,ndims,alloc_err
  integer :: n,KMAX,MAXMAX
  integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,iVarID,jVarID,kVarID,i,j,k
  integer :: ndate(4)
  integer, allocatable,dimension(:)  :: values
  real, allocatable,dimension(:)  :: rvalues
  real ::depsum
  character*100::attribute,attribute2
  character*100::name
  integer :: nDimensions,nVariables,nAttributes,xtype

  write(*,*)'  reading ',fileName
!  stop
!open an existing netcdf dataset
  call check(nf90_open(path = trim(fileName), mode = nf90_nowrite, ncid = ncFileID))

! find main properties
 call check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,timeDimID))
  print *, trim(fileName),' properties: '
  print *, 'Nb of dimensions: ',nDimensions
  print *, 'Nb of variables: ',nVariables
!  print *, 'Nb of global attributes: ',nAttributes
  do varid=1,nVariables
     call check(nf90_Inquire_Variable(ncFileID,varid,name,xtype,ndims))
     print *, varid,trim(name)!,'  Nb of dimensions: ',ndims
  enddo

!get global attributes
  call check(nf90_get_att(ncFileID, nf90_global, "Conventions", attribute ))
!write(*,*)'Convention:',attribute

if(attribute(1:2)=='CF')then
write(*,*)'reading longitude :'
     call check(nf90_inq_dimid(ncid = ncFileID, name = "longitude", dimID = idimID))
write(*,*)'reading latitude:'
     call check(nf90_inq_dimid(ncid = ncFileID, name = "latitude", dimID = jdimID))
write(*,*)'reading levelist:'
     status = nf90_inq_dimid(ncid = ncFileID, name = "levelist", dimID = kdimID)
     if(status == nf90_noerr) then
        call check(nf90_inq_dimid(ncid = ncFileID, name = "level", dimID = kdimID))
     endif


else
  !example:
!  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_hour", attribute ))
!  call check(nf90_get_att(ncFileID, nf90_global, "lastmodified_date", attribute2 ))
!  print *,'file last modified (yyyymmdd hhmmss.sss) ',attribute2,' ',attribute

     !get dimensions id
     call check(nf90_inq_dimid(ncid = ncFileID, name = "i", dimID = idimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "j", dimID = jdimID))
     call check(nf90_inq_dimid(ncid = ncFileID, name = "k", dimID = kdimID))
endif
     call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

     !get dimensions length
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=idimID,  len=GIMAX))
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=jdimID,  len=GJMAX))
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=kdimID,  len=KMAX_MID))
     call check(nf90_inquire_dimension(ncid=ncFileID, dimID=timedimID,  len=nrecords))

     print *, 'dimensions ',GIMAX,GJMAX,KMAX_MID,nrecords
     if(nrecords>Nrec)then
        write(*,*)' Nrec too small'
     endif
     allocate(values(nrecords), stat=alloc_err)
     if ( alloc_err /= 0 ) then
        print *, 'alloc failed in ReadCDF_ml: ',alloc_err,nrecords
        stop
     endif


    !get time variable
     call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = timeDimID))
     call check(nf90_get_var(ncFileID, timeDimID, values))

if(attribute(1:2)=='CF')then
     write(*,*)'date start '
     call datefromhourssince1900(ndate,values(1),1)
     write(*,*)'date end '
     call datefromhourssince1900(ndate,values(nrecords),1)
     period=0
     if(nrecords>1)then
        period=(values(2)-values(1))/1.
        if(period<24.1)then
           write(*,*)'Period between 2 first records (hours): ',period
        else
           write(*,*)'Period between 2 first records (days): ',period/24.
        endif
     endif
     allocate(rvalues(GIMAX), stat=alloc_err)
     call check(nf90_get_var(ncFileID, iDimID, rvalues))
     do i=1,gimax
        long(i)=rvalues(i)
     enddo
     deallocate(rvalues)
 
     allocate(rvalues(GJMAX), stat=alloc_err)
     call check(nf90_get_var(ncFileID, jDimID, rvalues))
     do i=1,gjmax
        lat(i)=rvalues(i)
     enddo
     deallocate(rvalues)
 
     allocate(rvalues(KMAX_MID), stat=alloc_err)
     call check(nf90_get_var(ncFileID, kDimID, rvalues))
     do i=1,KMAX_MID
        level(i)=rvalues(i)
     enddo


else
     write(*,*)'date start '
     call datefromsecondssince1970(ndate,values(1),1)
     write(*,*)'date end '
     call datefromsecondssince1970(ndate,values(nrecords),1)
     period=0
     if(nrecords>1)then
        period=(values(2)-values(1))/3600.
        if(period<24.1)then
           write(*,*)'Period between 2 first records (hours): ',period
        else
           write(*,*)'Period between 2 first records (days): ',period/24.
        endif
     endif
endif

     do n=1,min(Nrec,nrecords)
        timev(n)=values(n)
     enddo

     deallocate(values)
if(attribute(1:2)=='CF')then

else
   stop
endif

  call check(nf90_close(ncFileID))

   end subroutine GetCDFinfo
  subroutine lb2ij(imax,jmax,gl,gb,ir2,jr2,fi2,an2,xp2,yp2)
    !-------------------------------------------------------------------! 
    !      calculates coordinates ir2, jr2 (real values) from gl(lat),gb(long) 
    !
    !      input:  xp2,yp2:   coord. of the polar point in grid2
    !              an2:   number of grid-distances from pole to equator in grid2.
    !              fi2:      rotational angle for the grid2 (at i2=0).
    !              i1max,j1max: number of points (grid1) in  x- og y- direction
    !
    !
    !      output: i2(i1,j1): i coordinates in grid2 
    !              j2(i1,j1): j coordinates in grid2 
    !-------------------------------------------------------------------! 

!    use Par_ml,   only : MAXLIMAX, MAXLJMAX

    implicit none


    integer :: imax,jmax,i1, j1
    real    :: fi2,an2,xp2,yp2
    real    :: gl(imax,jmax),gb(imax,jmax)
    real    :: ir2(imax,jmax),jr2(imax,jmax)

    real, parameter :: PI=3.14159265358979323
    real    :: PId4,dr,dr2


    PId4    = PI/4.      
    dr2    = PI/180.0/2.      ! degrees to radians /2
    dr    = PI/180.0      ! degrees to radians 

    do j1 = 1, jmax
       do i1 = 1, imax

          ir2(i1,j1)=xp2+an2*tan(PId4-gb(i1,j1)*dr2)*sin(dr*(gl(i1,j1)-fi2))
          jr2(i1,j1)=yp2-an2*tan(PId4-gb(i1,j1)*dr2)*cos(dr*(gl(i1,j1)-fi2))
!          write(*,*)i1,j1,ir2(i1,j1),jr2(i1,j1),gl(i1,j1),gb(i1,j1)

       end do ! i
    end do ! j

    return
  end subroutine lb2ij

  subroutine ij2lb(imax,jmax,gl,gb,fi,an,xp,yp)
  !-------------------------------------------------------------------! 
  !      calculates l(lat),b(long) (geographical coord.) 
  !      in every grid point. 
  !
  !      input:  xp,yp:   coord. of the polar point.
  !              an:      number of grid-distances from pole to equator.
  !              fi:      rotational angle for the x,y grid (at i=0).
  !              imax,jmax:   number of points in  x- og y- direction
  !              glmin:   gives min.value of geographical lenght
  !                       =>  glmin <= l <= glmin+360.  
  !                           (example glmin = -180. or 0.)
  !                       if "geopos","georek" is used
  !                       then glmin must be the lenght i(1,1) in the
  !                       geographical grid (gl1 to "geopos")
  !      output: gl(ii,jj): longitude glmin <= l <= glmin+360. 
  !              gb(ii,jj): latitude  -90. <= b <= +90. 
  !-------------------------------------------------------------------! 



    implicit none


    integer :: i, j, imax, jmax
    real    :: gl(imax,jmax),gb(imax,jmax)
    real    :: fi, an, xp, yp 
    real    :: om, om2, glmin, glmax,dy, dy2,rp,rb, rl, dx, dr
    real, parameter :: PI=3.14159265358979323


!    fi = -32.0
    glmin = 0.0

    glmax = glmin + 360.0
    dr    = PI/180.0      ! degrees to radians
    om    = 180.0/PI      ! radians to degrees (om=Norwegian omvendt?)
    om2   = om * 2.0

    do j = 1, jmax          
       dy  = yp - j            
       dy2 = dy*dy
       do i = 1, imax       

         dx = i - xp    ! ds - changed
         rp = sqrt(dx*dx+dy2)           ! => distance to pole
         rb = 90.0 - om2 * atan(rp/AN)  ! => latitude
!pw 1/12-04        rl = 0.0
         rl = fi !pw 1/12-04 
         if (rp >  1.0e-10) rl = fi + om*atan2(dx,dy)
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         gl(i,j)=rl                     !     longitude
         gb(i,j)=rb                     !     latitude
!         write(*,*)i,j,gl(i,j),gb(i,j)
       end do ! i
    end do ! j

   return
  end subroutine ij2lb


  subroutine secondssince1970(ndate,nseconds)
    !calculate how many seconds have passed since the start of the year

!    use Dates_ml, only: nmdays,is_leap 
    implicit none

    integer, intent(in) :: ndate(4)
    integer, intent(out) :: nseconds
    integer :: n,nday,nmdays(12)
    nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/) 
    n=ndate(1)
       if(4*(n/4)==n)nmdays(2)=29
    nday=0
    do n=1,ndate(2)-1
       nday=nday+nmdays(n)
    enddo
    nday=nday+ndate(3)

    nseconds=3600*(ndate(4)+24*(nday-1))

!add seconds from each year since 1970
     do n=1970,ndate(1)-1
       nseconds=nseconds+24*3600*365
       if(4*(n/4)==n)nseconds=nseconds+24*3600
     enddo
!     write(*,*)ndate(1),ndate(2),ndate(3),ndate(4),nseconds
  end subroutine secondssince1970


subroutine writemeteofield(fileName,varname, field,IMAX,JMAX,KMAX,Nhh_in,ndate,validity,units)


  use netcdf
  implicit none

  integer,intent(in) :: IMAX,JMAX,KMAX,Nhh_in,ndate(4,Nhh_in)
  character (len=*),intent(in)::fileName,varname,validity,units
  real,intent(inout) :: field(*)
  real :: xmin,xmax,scale_factor,add_offset,ndays
  integer,  parameter :: Int1=1,Int2=2,Int4=3,Real4=4,Real8=5 !CDF typr for outp
  integer :: Wrec,ncFileID,idimID,jdimID,kdimID,VarID,tVarID,xtype,ndims
  integer :: GIMAX_old,GJMAX_old,KMAX_old,OUTtype,ndim,nrecords,nseconds,Nwrite
  character*8 ::lastmodified_date
  character*10 ::lastmodified_hour
  integer:: i,j,k,ijk,ijkh,status,ihh,Nhh
  integer*2, allocatable, dimension(:,:,:)::Ifield_2D
  integer*2, allocatable, dimension(:,:,:,:)::Ifield_3D
  character*100 ::name

  Nhh=Nhh_in

  Wrec=0
  ndim=2
  if(KMAX>1)ndim=3
  Nwrite=IMAX*JMAX*KMAX*Nhh
  
  call check(nf90_open(path=trim(fileName),mode=nf90_write,ncid=ncFileID))
  !test if the defined dimensions are compatible 
  call check(nf90_Inquire_Variable(ncFileID,1,name,xtype,ndims))
  call check(nf90_inq_dimid(ncid = ncFileID, name = trim(name), dimID = idimID))
  call check(nf90_Inquire_Variable(ncFileID,2,name,xtype,ndims))
  call check(nf90_inq_dimid(ncid = ncFileID, name = trim(name), dimID = jdimID))
  call check(nf90_Inquire_Variable(ncFileID,3,name,xtype,ndims))
  call check(nf90_inq_dimid(ncid = ncFileID, name = trim(name), dimID = kdimID))

  call check(nf90_inquire_dimension(ncid=ncFileID,dimID=idimID,len=GIMAX_old))
  call check(nf90_inquire_dimension(ncid=ncFileID,dimID=jdimID,len=GJMAX_old))
  call check(nf90_inquire_dimension(ncid=ncFileID,dimID=kdimID,len=KMAX_old))
     write(6,*)GIMAX_old,IMAX,GJMAX_old,JMAX,KMAX_old,KMAX

  if(GIMAX_old/=IMAX .or. GJMAX_old/=JMAX .or.  KMAX_old<KMAX)then
     write(6,*)'file ', trim(fileName),' has incompatible dimensions'
     write(6,*)GIMAX_old,IMAX,GJMAX_old,JMAX,KMAX_old,KMAX
     stop
  endif

  !test first if the variable is already defined:
  status = nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID)
  if(status == nf90_noerr) then     
     print *, 'variable exists: ',varname
  else
     write(6,*) 'creating variable: ',trim(varname)!,nf90_strerror(status)
      OUTtype=Int2
    call  createnewvariable(ncFileID,varname,ndim,ndate(:,1),OUTtype,validity,units)
  endif


  !Pack the data in "short" type
  !unpacked_value = scale_factor * packed_value + add_offset

!NOTE: a different scaling factor for each level could be more
!accurate for fields which vary much with height. But this may not be 
!compatible with any conventions?

     xmax=maxval(field(1:Nwrite))
     xmin=minval(field(1:Nwrite))

     if(xmax-xmin<1.E-30)then
        add_offset=xmax
        scale_factor=1.0
     else
        add_offset=(xmax+xmin)/2.
        scale_factor=(xmax-xmin)/65532.
     endif

  !get variable id
  call check(nf90_inq_varid(ncid = ncFileID, name = varname, varID = VarID))

  !find the number of records already written
  call check(nf90_get_att(ncFileID, VarID, "numberofrecords",   nrecords))
  write(*,*)'number of dataset already saved: ',nrecords
  !increase the last coordinate by one, to define position of new data
  if(nrecords>0)then
     write(*,*)'WARNING: overwriting fields',nrecords
  endif
  nrecords=1
  write(*,*)'writing on records: ',nrecords,'to ',nrecords+Nhh-1

if(KMAX==1)then
   allocate(Ifield_2D(IMAX,JMAX,Nhh))
   ijkh=0
   do ihh=1,Nhh
!      do j=JMAX,1,-1!NB:REVERSE
      do j=1,JMAX
         do i=1,IMAX
            ijkh=ijkh+1
            Ifield_2D(i,j,ihh)=nint((field(ijkh)-add_offset)/scale_factor)
         enddo
      enddo
   enddo
   call check(nf90_put_var(ncFileID,VarID,Ifield_2D,start=(/1,1,nrecords/)))
   deallocate(Ifield_2D)
else if(ndim==3)then
   allocate(Ifield_3D(IMAX,JMAX,KMAX,Nhh))
   ijkh=0
   do ihh=1,Nhh
      do k=1,KMAX
!         do j=JMAX,1,-1!NB:REVERSE
         do j=1,JMAX
            do i=1,IMAX
               ijkh=ijkh+1
               Ifield_3D(i,j,k,ihh)=nint((field(ijkh)-add_offset)/scale_factor)
            enddo
         enddo
      enddo
   enddo

     call check(nf90_put_var(ncFileID, VarID,&
         Ifield_3D, start = (/ 1, 1, 1,nrecords /)) )
     deallocate(Ifield_3D)

  else
     stop
  endif


     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale_factor ))
     call check(nf90_put_att(ncFileID, varID, "add_offset",  add_offset ))

!update dates
  call Date_And_Time(date=lastmodified_date,time=lastmodified_hour)
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_date", lastmodified_date))
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour))

  call check(nf90_put_att(ncFileID, VarID, "meteo_date_last",ndate(:,Nhh)))

! define time of new records
  nrecords=nrecords-1
  do ihh=1,Nhh
     nrecords=nrecords+1
  call check(nf90_inq_varid(ncid = ncFileID, name = "time", varID = tVarID))
  call dayssince1900(ndate(1,ihh),ndays)
  call check(nf90_put_var(ncFileID, tVarID, ndays, start = (/nrecords/) ) )
!  call secondssince1970(ndate(:,ihh),nseconds)
!  call check(nf90_put_var(ncFileID, tVarID, nseconds, start = (/nrecords/) ) )
  enddo
   call check(nf90_put_att(ncFileID, VarID, "numberofrecords",   nrecords))

  call check(nf90_close(ncFileID))

   end subroutine writemeteofield


subroutine  createnewvariable(ncFileID,varname,ndim,ndate,OUTtype,validity,units)

  !create new netCDF variable

  use netcdf

  implicit none

!  type(Deriv),     intent(in) :: def1 ! definition of fields
  character (len = *),intent(in) ::varname,validity,units
  integer ,intent(in) ::ndim,ncFileID,OUTtype
  integer, intent(in) ::  ndate(4)

  integer :: iDimID,jDimID,kDimID,timeDimID,xtype,ndims
  integer :: varID,nrecords,OUTtypeCDF
  real :: scale,offset
  integer,  parameter :: Int1=1,Int2=2,Int4=3,Real4=4,Real8=5 !CDF typr for output
 character*100 ::name
!  integer :: OUTtypeCDF !NetCDF code for type


  if(OUTtype==Int1)then
     OUTtypeCDF=nf90_byte
  elseif(OUTtype==Int2)then
     OUTtypeCDF=nf90_short
  elseif(OUTtype==Int4)then
     OUTtypeCDF=nf90_int
  elseif(OUTtype==Real4)then
     OUTtypeCDF=nf90_float
  elseif(OUTtype==Real8)then
     OUTtypeCDF=nf90_double
  else
     stop
!     call gc_abort(me,NPROC,"NetCDF_ml: undefined datatype")
  endif


     call check(nf90_redef(ncid = ncFileID))

     !get dimensions id
  !get dimensions id
     call check(nf90_Inquire_Variable(ncFileID,1,name,xtype,ndims))
     call check(nf90_inq_dimid(ncid = ncFileID, name = trim(name), dimID = idimID))
     call check(nf90_Inquire_Variable(ncFileID,2,name,xtype,ndims))
     call check(nf90_inq_dimid(ncid = ncFileID, name = trim(name), dimID = jdimID))
     call check(nf90_Inquire_Variable(ncFileID,3,name,xtype,ndims))
     call check(nf90_inq_dimid(ncid = ncFileID, name = trim(name), dimID = kdimID))

     call check(nf90_inq_dimid(ncid = ncFileID, name = "time", dimID = timeDimID))

     !define new variable
     if(ndim==3)then
        call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = OUTtypeCDF,     &
             dimids = (/ iDimID, jDimID, kDimID , timeDimID/), varID=varID ) )
     elseif(ndim==2)then
        call check(nf90_def_var(ncid = ncFileID, name = varname, xtype = OUTtypeCDF,     &
             dimids = (/ iDimID, jDimID , timeDimID/), varID=varID ) )
     else
         print *, 'createnewvariable: unexpected ndim ',ndim   
     endif
!compress variable
        call check(nf90_def_var_deflate(ncFileid ,varID,shuffle=0,deflate=1 ,deflate_level=4) )

!     FillValue=0.
     !define attributes of new variable
     call check(nf90_put_att(ncFileID, varID, "long_name",  trim(varname) ))
!     call check(nf90_put_att(ncFileID, varID, "coordinates", "lat lon"))
!     call check(nf90_put_att(ncFileID, varID, "grid_mapping", "ECMWF_latlon"))
     nrecords=0
     call check(nf90_put_att(ncFileID, varID, "numberofrecords", nrecords))

     call check(nf90_put_att(ncFileID, varID, "units", units))

     scale=1.
!     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))
     offset=0.
 !    call check(nf90_put_att(ncFileID, varID, "add_offset",  offset ))

  if(OUTtype==Int1)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_byte  ))
      call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))
    call check(nf90_put_att(ncFileID, varID, "add_offset",  offset ))
 elseif(OUTtype==Int2)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_short  ))
     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))
    call check(nf90_put_att(ncFileID, varID, "add_offset",  offset ))
  elseif(OUTtype==Int4)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_int   ))
     call check(nf90_put_att(ncFileID, varID, "scale_factor",  scale ))
    call check(nf90_put_att(ncFileID, varID, "add_offset",  offset ))
  elseif(OUTtype==Real4)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_float  ))
  elseif(OUTtype==Real8)then
     call check(nf90_put_att(ncFileID, varID, "_FillValue", nf90_fill_double  ))
  endif
!     call check(nf90_put_att(ncFileID, varID, "periodlength",   "yearly"))

     call check(nf90_put_att(ncFileID, varID, "meteo_date_first",ndate ))
     call check(nf90_put_att(ncFileID, varID, "meteo_date_last",ndate ))
    call check(nf90_put_att(ncFileID, VarID, "validity", trim(validity)))
  
     call check(nf90_enddef(ncid = ncFileID))

end subroutine  createnewvariable

subroutine CreatenetCDFfile(fileName,GIMAX,GJMAX,KMAX,lonstart,latstart,ddeg_lon,ddeg_lat,hyai,hybi,GRIDWIDTH_M)

!  use typeSizes
  use netcdf
   implicit none

character(len=*),  intent(in)  :: fileName 
integer,  intent(in) :: GIMAX,GJMAX,KMAX
real,  intent(inout) :: hyai(0:KMAX),hybi(0:KMAX),GRIDWIDTH_M
real,  intent(in) :: lonstart,latstart,ddeg_lon,ddeg_lat
character (len=*), parameter :: projection='lon lat'
!character (len=*), parameter :: vert_coord='vertical coordinates = (p-p(top))/(p(surf)-p(top))'
character (len=*), parameter :: vert_coord='atmosphere_hybrid_sigma_pressure_coordinate'
character (len=19) :: projection_params='90.0 -32.0 0.933013' !set later on
integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,VarID,iVarID,jVarID,kVarID
integer :: hyamVarID,hybmVarID,hyaiVarID,hybiVarID,ilevVarID,levVarID,levDimID,ilevDimID
integer :: latVarID,longVarID,xmVarID,xmiVarID,xmjVarID
integer :: i,j,k
!real (kind = FourByteReal),allocatable :: xcoord(:),ycoord(:),kcoord(:)
real , allocatable :: xcoord(:),ycoord(:),kcoord(:)
real,allocatable :: gb_glob(:,:),gl_glob(:,:),xm_i(:,:),xm_j(:,:)

character*8 ::created_date,lastmodified_date
character*10 ::created_hour,lastmodified_hour
real :: glmin,glmax,dy,dx,rp,rb,rl,om
    real, parameter :: PI=3.14159265358979323
    real, parameter :: EARTH_RADIUS= 6.370e6 !EMEP definition
                                ! NB: other can use different values!
    real ::   AN !Distance on the map from pole to equator (No. of cells)
real :: izero,jzero,scale_at_projection_origin
  real :: great_circle_distance
  real :: P0=101325.0!Pa
  real :: P0_Pa=101325.0!Pa
  real :: PS0=1013.250!hPa

!wgs84
real, parameter :: a = 6378137.0
real, parameter :: f = 1.0/298.257223563
real, parameter :: e = sqrt(2*f-f*f)
real:: lat,x
allocate(xm_j(GIMAX,GJMAX))
allocate(xm_i(GIMAX,GJMAX))

write(*,*)'create ',trim(fileName)
write(*,*)'with sizes (IMAX,JMAX,KMAX) ',GIMAX,GJMAX,KMAX

!create file
call check(nf90_create(path = trim(fileName), cmode = nf90_hdf5, ncid = ncFileID))

  ! Define the dimensions
  call check(nf90_def_dim(ncid = ncFileID, name = "lon", len = GIMAX, dimid = iDimID))
  call check(nf90_def_var(ncFileID, "lon", nf90_double, dimids = iDimID, varID = iVarID) )
  call check(nf90_put_att(ncFileID, iVarID, "standard_name", "longitude"))
  call check(nf90_put_att(ncFileID, iVarID, "long_name", "longitude"))
  call check(nf90_put_att(ncFileID, iVarID, "units", "degrees_east"))

  call check(nf90_def_dim(ncid = ncFileID, name = "lat", len = GJMAX, dimid = jDimID))
  call check(nf90_def_var(ncFileID, "lat", nf90_double, dimids = jDimID, varID =jVarID) )
  call check(nf90_put_att(ncFileID, jVarID, "standard_name", "latitude"))
  call check(nf90_put_att(ncFileID, jVarID, "long_name", "latitude"))
  call check(nf90_put_att(ncFileID, jVarID, "units", "degrees_north"))

!  call check(nf90_def_dim(ncid = ncFileID, name = "k", len = KMAX, dimid = kDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "lev", len = KMAX, dimid = levDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "ilev", len = KMAX+1, dimid = ilevDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "time", len = nf90_unlimited, dimid = timeDimID))


! Write global attributes
  call check(nf90_put_att(ncFileID, nf90_global, "Conventions", "CF-1.6" ))
  call check(nf90_put_att(ncFileID, nf90_global, "projection",trim(projection)))
  call check(nf90_put_att(ncFileID, nf90_global, "reference","WGS84"))
!  scale_at_projection_origin=(1.+sin(ref_latitude*PI/180.))/2.
!  write(projection_params,fmt='(''90.0 '',F5.1,F9.6)')fi,scale_at_projection_origin
 ! call check(nf90_put_att(ncFileID, nf90_global, "projection_params",projection_params))
  call check(nf90_put_att(ncFileID, nf90_global, "vert_coord", vert_coord))

 !  call check(nf90_put_att(ncFileID, nf90_global, "xcoordinate_NorthPole",xp ))
!  call check(nf90_put_att(ncFileID, nf90_global, "ycoordinate_NorthPole",yp ))
!  call check(nf90_put_att(ncFileID, nf90_global, "fi",fi ))
!  call check(nf90_put_att(ncFileID, nf90_global, "ref_latitude",ref_latitude))
!GRIDWIDTH_M=EARTH_RADIUS*360.0/GIMAX*PI/180.0
GRIDWIDTH_M=a*ddeg_lat*PI/180.0
  call check(nf90_put_att(ncFileID, nf90_global, "Grid_resolution",GRIDWIDTH_M))


!already defined as dimension
  call check(nf90_def_var(ncFileID, "lev", nf90_double, dimids = levDimID, varID = levVarID) )
  call check(nf90_put_att(ncFileID, levVarID, "standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
!  call check(nf90_put_att(ncFileID, kVarID, "coord_alias", "level"))!  
  call check(nf90_put_att(ncFileID, levVarID, "long_name", "hybrid level at layer midpoints (1000*(A+B))"))
  call check(nf90_put_att(ncFileID, levVarID, "positive", "down"))
  call check(nf90_put_att(ncFileID, levVarID, "formula_terms","a: hyam b: hybm ps: PS p0: P0"))
!p(n,k,j,i) = a(k)*p0 + b(k)*ps(n,j,i)
  call check(nf90_def_var(ncFileID, "P0", nf90_double,  varID = VarID) )
  call check(nf90_put_att(ncFileID, VarID, "units", "Pa"))
  call check(nf90_put_var(ncFileID, VarID, P0_Pa ))

!The hybrid sigma-pressure coordinate for level k is defined as a(k)+b(k) (or ap(k)/p0+b(k), as appropriate). 
  call check(nf90_def_var(ncFileID, "hyam", nf90_double,dimids = levDimID,  varID = hyamVarID) )
  call check(nf90_put_att(ncFileID, hyamVarID, "long_name","hybrid A coefficient at layer midpoints"))
  call check(nf90_def_var(ncFileID, "hybm", nf90_double,dimids = levDimID,  varID = hybmVarID) )
  call check(nf90_put_att(ncFileID, hybmVarID, "long_name","hybrid B coefficient at layer midpoints"))

  call check(nf90_def_var(ncFileID, "ilev", nf90_double, dimids = ilevDimID, varID = ilevVarID) )
  call check(nf90_put_att(ncFileID, ilevVarID, "standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
  call check(nf90_put_att(ncFileID, ilevVarID, "long_name", "hybrid level at layer interfaces (1000*(A+B))."))
  call check(nf90_put_att(ncFileID, ilevVarID, "positive", "down"))
  call check(nf90_put_att(ncFileID, ilevVarID, "formula_terms","a: hyai b: hybi ps: PS p0: P0"))
  call check(nf90_def_var(ncFileID, "hyai", nf90_double, dimids = ilevDimID,  varID = hyaiVarID) )
  call check(nf90_put_att(ncFileID, hyaiVarID, "long_name","hybrid A coefficient at layer interfaces"))
  call check(nf90_def_var(ncFileID, "hybi", nf90_double, dimids = ilevDimID,  varID = hybiVarID) )
  call check(nf90_put_att(ncFileID, hybiVarID, "long_name","hybrid B coefficient at layer interfaces"))
!  call check(nf90_def_var(ncFileID, "k", nf90_double, dimids = kDimID, varID = kVarID) )
!  call check(nf90_put_att(ncFileID, kVarID, "coord_alias", "level"))
!  call check(nf90_put_att(ncFileID, kVarID, "long_name", "vertical sigma coordinates"))
!  call check(nf90_put_att(ncFileID, kVarID, "units", "sigma_level"))
!  call check(nf90_put_att(ncFileID, kVarID, "positive", "down"))


  call check(nf90_def_var(ncFileID, "time", nf90_double, dimids = timeDimID, varID = VarID) )
  call check(nf90_put_att(ncFileID, VarID, "long_name", "time at middle of period"))
  call check(nf90_put_att(ncFileID, VarID, "units", "days since 1900-1-1 00:00:00.0 +00:00"))


!additional attributes

!write time of creation
  call Date_And_Time(date=created_date,time=created_hour)
     write(6,*) 'created_date: ',created_date
     write(6,*) 'created_hour: ',created_hour
  call check(nf90_put_att(ncFileID, nf90_global, "created_date", created_date))
  call check(nf90_put_att(ncFileID, nf90_global, "created_hour", created_hour))
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_date", created_date))
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_hour",created_hour ))

  call check(nf90_def_var(ncFileID, "map_factor_i", nf90_double, dimids = (/ iDimID, jDimID/), varID = xmiVarID) )
  call check(nf90_put_att(ncFileID, xmiVarID, "long_name", "mapping factor in i direction"))
  call check(nf90_put_att(ncFileID, xmiVarID, "units", " "))

  call check(nf90_def_var(ncFileID, "map_factor_j", nf90_double, dimids = (/ iDimID, jDimID/), varID = xmjVarID) )
  call check(nf90_put_att(ncFileID, xmjVarID, "long_name", "mapping factor in j direction"))
  call check(nf90_put_att(ncFileID, xmjVarID, "units", " "))

  ! Leave define mode
  call check(nf90_enddef(ncFileID))


! write values of variables

!  call check(nf90_open(path = trim(fileName), mode = nf90_write, ncid = ncFileID))

! Define horizontal distances in GDV conventions

allocate(xcoord(0:GIMAX+1))
allocate(ycoord(GJMAX+1))
allocate(kcoord(0:KMAX))

  call check(nf90_put_var(ncFileID, hyaiVarID, hyai(0:KMAX)/P0_Pa) )
  call check(nf90_put_var(ncFileID, hybiVarID, hybi(0:KMAX)) )

  do i=1,KMAX
     kcoord(i)=(hyai(i-1)+hyai(i))/2.0/P0_Pa
  enddo
  call check(nf90_put_var(ncFileID, hyamVarID, kcoord(1:KMAX)) )
  do i=1,KMAX
     kcoord(i)=(hybi(i-1)+hybi(i))/2.0
  enddo
  call check(nf90_put_var(ncFileID, hybmVarID, kcoord(1:KMAX)) )

  do i=0,KMAX
     kcoord(i)=1000*(hyai(i)/P0_Pa+hybi(i))
  enddo

  call check(nf90_put_var(ncFileID, ilevVarID, kcoord(0:KMAX)) )
  do i=1,KMAX
     kcoord(i)=500*(hyai(i-1)/P0_Pa+hybi(i-1)+hyai(i)/P0_Pa+hybi(i))
  enddo
  call check(nf90_put_var(ncFileID, levVarID, kcoord(1:KMAX)) )

  do i=0,GIMAX+1
!     xcoord(i)=(i-1)*360.0/GIMAX*PI/180.0
!     xcoord(i)=(i-1)*360.0/GIMAX
     xcoord(i)=(i-1)*ddeg_lon + lonstart
  enddo

  call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAX)) )

  do j=1,GJMAX+1
!     ycoord(j)=-PI/2.0+(j-1)*180.0/(GJMAX-1)*PI/180.0!NB: j reversed in writemeteofield
     ycoord(j)=(j-1)*ddeg_lat + latstart
  enddo

  call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAX)) )


!allocate(gb_glob(GIMAX,GJMAX))
!allocate(gl_glob(GIMAX,GJMAX))



!mapping factor at the "upper" edges 
xcoord=xcoord*PI/180.0
ycoord=ycoord*PI/180.0
!  do j=1,GJMAX
!  do i=1,GIMAX
!     xm(i,j)=great_circle_distance(0.5*(xcoord(i-1)+xcoord(i)),&
!          0.5*(ycoord(j)+ycoord(j+1)),0.5*(xcoord(i)+xcoord(i+1)),&
!          0.5*(ycoord(j)+ycoord(j+1)))
!  enddo
!  enddo
!  xm=GRIDWIDTH_M/(xm*EARTH_RADIUS)
!xm_j is defined at the right edge, but since it is longitude independent it does not matter
lat = latstart
do j=1,gjmax
!   lon = lonstart
!   do i=1,imax
      x=1.0/(1-e*e*(sin(PI*lat/180))**2)
!      write(*,44)lat,(1-e*e)*x*sqrt(x),a*PI/180*(1-e*e)*x*sqrt(x),sqrt(x)*cos(PI*lat/180),a*PI/180*sqrt(x)*cos(PI*lat/180)
      xm_j(:,j)=1.0/((1-e*e)*x*sqrt(x))
!      xm_i(:,j)=1.0/(sqrt(x)*cos(PI*lat/180))
!      lon=lon+ddeg_lon      
!   enddo
   lat=lat+ddeg_lat   
enddo

  call check(nf90_put_var(ncFileID, xmjVarID, xm_j(1:GIMAX,1:GJMAX)) )

!NB: xm_i is defined at the upper edge
lat = latstart+0.5*ddeg_lat 
do j=1,gjmax
!   lon = lonstart
!   do i=1,imax
      x=1.0/(1-e*e*(sin(PI*lat/180))**2)
44 format(f6.1,f16.10,12f15.5)
!      write(*,44)lat,(1-e*e)*x*sqrt(x),a*PI/180*(1-e*e)*x*sqrt(x),sqrt(x)*cos(PI*lat/180),a*PI/180*sqrt(x)*cos(PI*lat/180)
!      xm_j(:,j)=1.0/((1-e*e)*x*sqrt(x))
      xm_i(:,j)=1.0/(sqrt(x)*cos(PI*lat/180))
!      lon=lon+ddeg_lon      
!   enddo
   lat=lat+ddeg_lat   
enddo
!NB: map_factor xm_i are defined with respect to ddeg_lon
  xm_i=xm_i*ddeg_lat/ddeg_lon

  call check(nf90_put_var(ncFileID, xmiVarID, xm_i(1:GIMAX,1:GJMAX)) )

  call check(nf90_close(ncFileID))


deallocate(xcoord)
deallocate(ycoord)
deallocate(kcoord)
deallocate(xm_i,xm_j)
!deallocate(gb_glob)
!deallocate(gl_glob)


end subroutine CreatenetCDFfile

function great_circle_distance(fi1,lambda1,fi2,lambda2) result(dist)

!compute the great circle distance between to points given in 
!spherical coordinates. Sphere has radius 1.
real, intent(in) ::fi1,lambda1,fi2,lambda2
real :: dist

dist=2*asin(sqrt(sin(0.5*(lambda1-lambda2))**2+cos(lambda1)*cos(lambda2)*sin(0.5*(fi1-fi2))**2))

end function great_circle_distance

subroutine CreatenetCDFfilePS(fileName,GIMAX,GJMAX,KMAX,sigma,GRIDWIDTH_M,xp,yp,fi,ref_latitude,hyai,hybi,xm)

!  use typeSizes
  use netcdf
   implicit none

character(len=*),  intent(in)  :: fileName 
integer,  intent(in) :: GIMAX,GJMAX,KMAX
real,  intent(in) :: GRIDWIDTH_M,xp,yp,fi,sigma(KMAX),ref_latitude
real,  intent(inout) :: hyai(0:KMAX),hybi(0:KMAX),xm(GIMAX,GJMAX)

character (len=*), parameter :: projection='Stereographic'
!character (len=*), parameter :: vert_coord='vertical coordinates = (p-p(top))/(p(surf)-p(top))'
character (len=*), parameter :: vert_coord='atmosphere_hybrid_sigma_pressure_coordinate'
character (len=19) :: projection_params='90.0 -32.0 0.933013' !set later on
integer :: ncFileID,iDimID,jDimID,kDimID,timeDimID,VarID,iVarID,jVarID,kVarID
integer :: hyamVarID,hybmVarID,hyaiVarID,hybiVarID,ilevVarID,levVarID,levDimID,ilevDimID
integer :: latVarID,longVarID,xmVarID
integer :: i,j,k
!real (kind = FourByteReal),allocatable :: xcoord(:),ycoord(:),kcoord(:)
real , allocatable :: xcoord(:),ycoord(:),kcoord(:)
real,  allocatable:: gb_glob(:,:),gl_glob(:,:)

character*8 ::created_date,lastmodified_date
character*10 ::created_hour,lastmodified_hour
real :: glmin,glmax,dy,dx,rp,rb,rl,om
    real, parameter :: PI=3.14159265358979323
    real, parameter :: EARTH_RADIUS= 6.370e6 !EMEP definition
                                ! NB: other can use different values!
    real ::   AN !Distance on the map from pole to equator (No. of cells)
  real :: P0=101325.0!Pa
  real :: P0_Pa=101325.0!Pa
  real :: PS0=1013.250!hPa

!AN = EARTH_RADIUS*(1.0+0.5*sqrt(3.0))/GRIDWIDTH_M    ! = 237.7316364 for GRIDWIDTH_M=50 km
AN = EARTH_RADIUS*(1.0+sin(ref_latitude*PI/180.))/GRIDWIDTH_M    ! = 237.7316364 for GRIDWIDTH_M=50 km and ref_latitude=60

write(*,*)'create ',trim(fileName)
write(*,*)'with sizes (IMAX,JMAX,KMAX) ',GIMAX,GJMAX,KMAX

!create file
call check(nf90_create(path = trim(fileName), cmode = nf90_hdf5, ncid = ncFileID))

  ! Define the dimensions
  call check(nf90_def_dim(ncid = ncFileID, name = "i", len = GIMAX, dimid = iDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "j", len = GJMAX, dimid = jDimID))
!  call check(nf90_def_dim(ncid = ncFileID, name = "k", len = KMAX, dimid = kDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "lev", len = KMAX, dimid = levDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "ilev", len = KMAX+1, dimid = ilevDimID))
  call check(nf90_def_dim(ncid = ncFileID, name = "time", len = nf90_unlimited, dimid = timeDimID))


! Write global attributes
  call check(nf90_put_att(ncFileID, nf90_global, "Conventions", "CF-1.6" ))
  call check(nf90_put_att(ncFileID, nf90_global, "projection",projection))
!  write(projection_params,fmt='(''90.0 '',F5.1,F9.6)')fi,(1.+sin(ref_latitude*PI/180.))/2.
!  call check(nf90_put_att(ncFileID, nf90_global, "projection_params",projection_params))
  call check(nf90_put_att(ncFileID, nf90_global, "vert_coord", vert_coord))

  call check(nf90_put_att(ncFileID, nf90_global, "Grid_resolution",GRIDWIDTH_M))
  call check(nf90_put_att(ncFileID, nf90_global, "xcoordinate_NorthPole",xp ))
  call check(nf90_put_att(ncFileID, nf90_global, "ycoordinate_NorthPole",yp ))
  call check(nf90_put_att(ncFileID, nf90_global, "fi",fi ))
  call check(nf90_put_att(ncFileID, nf90_global, "ref_latitude",ref_latitude))

! define coordinate variables
  call check(nf90_def_var(ncFileID, "i", nf90_double, dimids = iDimID, varID = iVarID) )
  call check(nf90_put_att(ncFileID, iVarID, "coord_axis", "x"))
  call check(nf90_put_att(ncFileID, iVarID, "long_name", "EMEP grid x coordinate"))
  call check(nf90_put_att(ncFileID, iVarID, "units", "km"))


  call check(nf90_def_var(ncFileID, "j", nf90_double, dimids = jDimID, varID = jVarID) )
  call check(nf90_put_att(ncFileID, jVarID, "coord_axis", "y"))
  call check(nf90_put_att(ncFileID, jVarID, "long_name", "EMEP grid y coordinate"))
  call check(nf90_put_att(ncFileID, jVarID, "units", "km"))


!already defined as dimension
  call check(nf90_def_var(ncFileID, "lev", nf90_double, dimids = levDimID, varID = levVarID) )
  call check(nf90_put_att(ncFileID, levVarID, "standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
!  call check(nf90_put_att(ncFileID, kVarID, "coord_alias", "level"))!  
  call check(nf90_put_att(ncFileID, levVarID, "long_name", "hybrid level at layer midpoints (1000*(A+B))"))
  call check(nf90_put_att(ncFileID, levVarID, "positive", "down"))
  call check(nf90_put_att(ncFileID, levVarID, "formula_terms","a: hyam b: hybm ps: PS p0: P0"))
!p(n,k,j,i) = a(k)*p0 + b(k)*ps(n,j,i)
  call check(nf90_def_var(ncFileID, "P0", nf90_double,  varID = VarID) )
  call check(nf90_put_att(ncFileID, VarID, "units", "Pa"))
  call check(nf90_put_var(ncFileID, VarID, P0_Pa ))

!The hybrid sigma-pressure coordinate for level k is defined as a(k)+b(k) (or ap(k)/p0+b(k), as appropriate). 
  call check(nf90_def_var(ncFileID, "hyam", nf90_double,dimids = levDimID,  varID = hyamVarID) )
  call check(nf90_put_att(ncFileID, hyamVarID, "long_name","hybrid A coefficient at layer midpoints"))
  call check(nf90_def_var(ncFileID, "hybm", nf90_double,dimids = levDimID,  varID = hybmVarID) )
  call check(nf90_put_att(ncFileID, hybmVarID, "long_name","hybrid B coefficient at layer midpoints"))

  call check(nf90_def_var(ncFileID, "ilev", nf90_double, dimids = ilevDimID, varID = ilevVarID) )
  call check(nf90_put_att(ncFileID, ilevVarID, "standard_name","atmosphere_hybrid_sigma_pressure_coordinate"))
  call check(nf90_put_att(ncFileID, ilevVarID, "long_name", "hybrid level at layer interfaces (1000*(A+B))."))
  call check(nf90_put_att(ncFileID, ilevVarID, "positive", "down"))
  call check(nf90_put_att(ncFileID, ilevVarID, "formula_terms","a: hyai b: hybi ps: PS p0: P0"))
  call check(nf90_def_var(ncFileID, "hyai", nf90_double, dimids = ilevDimID,  varID = hyaiVarID) )
  call check(nf90_put_att(ncFileID, hyaiVarID, "long_name","hybrid A coefficient at layer interfaces"))
  call check(nf90_def_var(ncFileID, "hybi", nf90_double, dimids = ilevDimID,  varID = hybiVarID) )
  call check(nf90_put_att(ncFileID, hybiVarID, "long_name","hybrid B coefficient at layer interfaces"))
!  call check(nf90_def_var(ncFileID, "k", nf90_double, dimids = kDimID, varID = kVarID) )
!  call check(nf90_put_att(ncFileID, kVarID, "coord_alias", "level"))
!  call check(nf90_put_att(ncFileID, kVarID, "long_name", "vertical sigma coordinates"))
!  call check(nf90_put_att(ncFileID, kVarID, "units", "sigma_level"))
!  call check(nf90_put_att(ncFileID, kVarID, "positive", "down"))


!  call check(nf90_def_var(ncFileID, "time", nf90_int, dimids = timeDimID, varID = VarID) )
!  call check(nf90_put_att(ncFileID, VarID, "long_name", "time at middle of period"))
!  call check(nf90_put_att(ncFileID, VarID, "units", "seconds since 1970-1-1 00:00:00.0 +00:00"))
  call check(nf90_def_var(ncFileID, "time", nf90_double, dimids = timeDimID, varID = VarID) )
  call check(nf90_put_att(ncFileID, VarID, "long_name", "time at middle of period"))
  call check(nf90_put_att(ncFileID, VarID, "units", "days since 1900-1-1 00:00:00.0 +00:00"))


  call check(nf90_def_var(ncFileID, "map_factor", nf90_double, dimids = (/ iDimID, jDimID/), varID = xmVarID) )
  call check(nf90_put_att(ncFileID, xmVarID, "long_name", "mapping factor"))

  call check(nf90_put_att(ncFileID, xmVarID, "units", " "))

  call check(nf90_def_var(ncFileID, "lat", nf90_double, dimids = (/ iDimID, jDimID/), varID = latVarID) )
  call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude coordinate"))
  call check(nf90_put_att(ncFileID, latVarID, "standard_name", "latitude"))
  call check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north"))

  call check(nf90_def_var(ncFileID, "lon", nf90_double, dimids = (/ iDimID, jDimID/), varID = longVarID) )
  call check(nf90_put_att(ncFileID, longVarID, "long_name", "longitude coordinate"))
  call check(nf90_put_att(ncFileID, longVarID, "standard_name", "longitude"))
  call check(nf90_put_att(ncFileID, longVarID, "units", "degrees_east"))



!additional attributes

!write time of creation
  call Date_And_Time(date=created_date,time=created_hour)
     write(6,*) 'created_date: ',created_date
     write(6,*) 'created_hour: ',created_hour
  call check(nf90_put_att(ncFileID, nf90_global, "created_date", created_date))
  call check(nf90_put_att(ncFileID, nf90_global, "created_hour", created_hour))
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_date", lastmodified_date))
  call check(nf90_put_att(ncFileID, nf90_global, "lastmodified_hour", lastmodified_hour))

  ! Leave define mode
  call check(nf90_enddef(ncFileID))



! write values of variables

  call check(nf90_open(path = trim(fileName), mode = nf90_write, ncid = ncFileID))

! Define horizontal distances in GDV conventions

allocate(xcoord(GIMAX))
allocate(ycoord(GJMAX))
allocate(kcoord(0:KMAX))

  call check(nf90_put_var(ncFileID, hyaiVarID, hyai(0:KMAX)/P0_Pa) )
  call check(nf90_put_var(ncFileID, hybiVarID, hybi(0:KMAX)) )

  do i=1,KMAX
     kcoord(i)=(hyai(i-1)+hyai(i))/2.0/P0_Pa
  enddo
  call check(nf90_put_var(ncFileID, hyamVarID, kcoord(1:KMAX)) )
  do i=1,KMAX
     kcoord(i)=(hybi(i-1)+hybi(i))/2.0
  enddo
  call check(nf90_put_var(ncFileID, hybmVarID, kcoord(1:KMAX)) )

  do i=0,KMAX
     kcoord(i)=1000*(hyai(i)/P0_Pa+hybi(i))
  enddo

  call check(nf90_put_var(ncFileID, ilevVarID, kcoord(0:KMAX)) )
  do i=1,KMAX
     kcoord(i)=500*(hyai(i-1)/P0_Pa+hybi(i-1)+hyai(i)/P0_Pa+hybi(i))
  enddo
  call check(nf90_put_var(ncFileID, levVarID, kcoord(1:KMAX)) )

  xcoord(1)=(1-xp)*GRIDWIDTH_M/1000.
  do i=2,GIMAX
     xcoord(i)=xcoord(i-1)+GRIDWIDTH_M/1000.
!     print *, i,xcoord(i)
  enddo

  call check(nf90_put_var(ncFileID, iVarID, xcoord(1:GIMAX)) )

  ycoord(1)=(1-yp)*GRIDWIDTH_M/1000.
  do j=2,GJMAX
     ycoord(j)=ycoord(j-1)+GRIDWIDTH_M/1000.
  enddo

  call check(nf90_put_var(ncFileID, jVarID, ycoord(1:GJMAX)) )

!  call check(nf90_put_var(ncFileID, kVarID, sigma(1:KMAX)) )


allocate(gl_glob(gimax,gjmax)) 
allocate(gb_glob(gimax,gjmax)) 

    glmin = -180.0
    glmax = glmin + 360.0
    do j = 1, GJMAX
       dy  = yp - j  
       do i = 1, GIMAX
         dx = i - xp    
         rp = sqrt(dx*dx+dy*dy)           ! => distance to pole
         rb = 90.0 - 360./PI * atan(rp/AN)  ! => latitude
         xm(i,j) = 0.5*(1.0+sin(ref_latitude*PI/180.))*(1.0 + (dx*dx+dy*dy)/AN**2 )

         if (rp >  1.0e-10)then
            rl = fi + 180.0/PI *atan2(dx,dy)
         else
            rl = fi
         endif
         if (rl <  glmin)   rl = rl + 360.0
         if (rl >  glmax)   rl = rl - 360.0
         gl_glob(i,j)=rl                     !     longitude
         gb_glob(i,j)=rb                     !     latitude
!         write(*,*)i,j,gl_glob(i,j)
       end do ! i
    end do ! j

i=8
j=110
write(*,*)'orig',i,j,gl_glob(i,j),gb_glob(i,j)
!stop

  call check(nf90_put_var(ncFileID, latVarID, gb_glob(1:GIMAX,1:GJMAX)) )
  call check(nf90_put_var(ncFileID, longVarID, gl_glob(1:GIMAX,1:GJMAX)) )
  call check(nf90_put_var(ncFileID, xmVarID, xm(1:GIMAX,1:GJMAX)) )

  call check(nf90_close(ncFileID))


deallocate(xcoord)
deallocate(ycoord)
deallocate(kcoord)
deallocate(gl_glob) 
deallocate(gb_glob) 




end subroutine CreatenetCDFfilePS

  subroutine dayssince1900(ndate,ndays)
    !calculate how many days have passed since the start of the year 1900
!NB: 1900 is not a leap year

    implicit none

    integer, intent(in) :: ndate(4)
    real*8, intent(out) :: ndays
    integer :: n,nday,nmdays(12),is_leap
    real*8 ::nseconds
    nmdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    n=ndate(1)
       if(4*(n/4)==n.and.n/=1900)nmdays(2)=29!NB: 1900 is not a leap year

    ndays=0.0
    do n=1,ndate(2)-1
       ndays=ndays+nmdays(n)!entire months since start of last year
    enddo
    ndays=ndays+ndate(3)-1!entire days since start of last month

    ndays=ndays+ndate(4)/24.0!hours since start of last day
!     write(*,*)'days since year start ',ndays

!add days from each entire year since 1900
     do n=1900,ndate(1)-1
       ndays=ndays+365
       if(4*(n/4)==n.and.n/=1900)ndays=ndays+1!NB: 1900 is not a leap year
     enddo
!     write(*,*)ndate(1),ndate(2),ndate(3),ndate(4),ndays


  end subroutine dayssince1900
  subroutine make_smi_interp(I_smi,J_smi,N_smi,varsEC,gimax,gjmax,maxnland)
    implicit none

    integer ::i,j,n,gimax,gjmax,maxnland,iter,i_EC,i_EC_cyclic,j_EC,nland
    real::sum
    real::varsEC(gimax,gjmax)
    integer*2::I_smi(gimax,gjmax,maxnland),J_smi(gimax,gjmax,maxnland),N_smi(gimax,gjmax)

   do i=1,gimax
     write(*,*)i
       do j=1,gjmax
          if(varsEC(i,j)<1.E-4)then
             !either sea, or undefined island, or desert
             if(.not.((j>290.and.j<381.and.(i<313.or.i>1745)).or.(j>321.and.j<367.and.i>1711)).or.varsEC(i,j)>-1.E-4)then
!             if(.not.((j>580.and.j<761.and.(i<625.or.i>3490)).or.(j>642.and.j<733.and.i>3422)).or.varsEC(i,j)>-1.E-4)then
                !not Sahara region
                
                !take nearest values
                do iter=1,gjmax,2
                   sum=0.0
                   nland=0
                   do j_EC=max(1,j-2*iter),min(gjmax,j+2*iter)
                      do i_EC=i-2*iter,i+2*iter!look at broader area in longitude than in latitude    
                         i_EC_cyclic=i_EC!cross the boundaries, cyclic 
                         if(i_EC_cyclic>gimax)i_EC_cyclic=i_EC_cyclic-gimax
                         if(i_EC_cyclic<1)i_EC_cyclic=i_EC_cyclic+gimax
                         if(varsEC(i_EC_cyclic,j_EC)>1.E-4)then
                            nland=nland+1
                            if(nland<=maxnland)then
                               I_smi(i,j,nland)=i_EC_cyclic
                               J_smi(i,j,nland)=j_EC
                            endif
                         endif
                      enddo
                   enddo
                   if(nland>0)then
                      N_smi(i,j)=min(nland,maxnland)
                      go to 643
                   endif
                   
                enddo
                write(*,*)'ERROR: no land found!',i,j,n
                stop
643             continue
             endif
          endif
       enddo
    enddo
    
 
  end subroutine make_smi_interp

subroutine findnearestneighbor_cyclic(val,Nval,coord,Ncoord,index1,index2)
implicit none
integer,intent(in)::Nval,Ncoord
real,intent(in)::val(Nval),coord(Ncoord)
integer,intent(out)::index1(Nval),index2(Nval)
integer::i1,i2,n,i
real ::dist,distmin1,distmin2
do n=1,Nval
   i1=1
   i2=Ncoord
   distmin1=360.0
   distmin2=360.0
   do i=1,Ncoord
      dist=abs(coord(i)-val(n))
      if(dist>180.0)dist=abs(dist-360.0)
      if(dist>180.0)dist=abs(dist-360.0)
      if(dist<distmin2)then
         
         if(dist<distmin1)then
            distmin2=distmin1
            i2=i1
            distmin1=dist
            i1=i
         else
            distmin2=dist
            i2=i
         endif
      endif
   enddo
   index1(n)=i1
   index2(n)=i2
199 FORMAT(3I4,10F12.5 )
!   write(*,199)n,index1(n),index2(n),val(n),coord(index1(n)),coord(index2(n))

enddo
!stop
end subroutine findnearestneighbor_cyclic
subroutine findnearestneighbor(val,Nval,coord,Ncoord,index1,index2)
implicit none
integer,intent(in)::Nval,Ncoord
real,intent(in)::val(Nval),coord(Ncoord)
integer,intent(out)::index1(Nval),index2(Nval)
integer::i1,i2,n,i
real ::dist,distmin1,distmin2
do n=1,Nval
   i1=1
   i2=Ncoord
   distmin1=180.0
   distmin2=180.0
   do i=1,Ncoord
      dist=abs(coord(i)-val(n))
      if(dist<distmin2)then
        if(dist<distmin1)then
            distmin2=distmin1
            i2=i1
            distmin1=dist
            i1=i
         else
            distmin2=dist
            i2=i
         endif
      endif
   enddo
   index1(n)=i1
   index2(n)=i2
199 FORMAT(3I4,10F12.5 )
!   write(*,199)n,index1(n),index2(n),val(n),coord(index1(n)),coord(index2(n))

enddo
!stop

end subroutine findnearestneighbor
subroutine findnearestneighbor_t(val,Nval,coord,Ncoord,index1,index2)
implicit none
integer,intent(in)::Nval,Ncoord
real,intent(in)::val(Nval),coord(Ncoord)
integer,intent(out)::index1(Nval),index2(Nval)
integer::i1,i2,n,i
real ::dist,distmin1,distmin2
do n=Nval,Nval
   i1=1
   i2=Ncoord
   distmin1=180.0
   distmin2=180.0
   do i=1,Ncoord
      dist=abs(coord(i)-val(n))
      write(*,*)i,i1,i2,coord(i),val(n),dist,distmin1,distmin2
      if(dist<distmin2)then
        if(dist<distmin1)then
            distmin2=distmin1
            i2=i1
            distmin1=dist
            i1=i
         else
            distmin2=dist
            i2=i
         endif
      endif
   enddo
   index1(n)=i1
   index2(n)=i2
199 FORMAT(3I4,10F12.5 )
   write(*,199)n,index1(n),index2(n),val(n),coord(index1(n)),coord(index2(n))

enddo
!stop

end subroutine findnearestneighbor_t
