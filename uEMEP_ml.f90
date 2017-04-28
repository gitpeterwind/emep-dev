module uEMEP_ml
!
! all subroutines for uEMEP
!
use CheckStop_ml,     only: CheckStop,StopAll
use Chemfields_ml,    only: xn_adv
use ChemSpecs,        only: NSPEC_ADV, NSPEC_SHL,species_adv
use Country_ml,       only: MAXNLAND,NLAND,Country
use EmisDef_ml,       only: loc_frac, loc_frac_ext, loc_frac_hour, loc_tot_hour, &
                            loc_frac_day, loc_tot_day, loc_frac_month&
                            , loc_tot_month,loc_frac_full,loc_tot_full, NSECTORS,NEMIS_FILE, &
                            EMIS_FILE,nlandcode,landcode,flat_nlandcode,flat_landcode,&
                            sec2tfac_map, sec2hfac_map ,ISNAP_DOM,snapemis,&
                            snapemis_flat,roaddust_emis_pot,KEMISTOP

use EmisGet_ml,       only: nrcemis, iqrc2itot, emis_nsplit,nemis_kprofile, emis_kprofile
use GridValues_ml,    only: dA,dB,xm2, dhs1i, glat, glon, projection
use MetFields_ml,     only: ps,roa
use ModelConstants_ml,only: KMAX_MID, KMAX_BND,USES, USE_uEMEP, uEMEP, IOU_HOUR, IOU_HOUR_INST,&
                            IOU_INST,IOU_YEAR,IOU_MON,IOU_DAY,IOU_HOUR,IOU_HOUR_INST, &
                            KMAX_MID,  MasterProc,dt_advec, RUNDOMAIN
use NetCDF_ml,        only: Real4,Out_netCDF_n
use OwnDataTypes_ml,  only: Deriv
use Par_ml,           only: me, LIMAX, LJMAX,gi0,gj0,li0,li1,lj0,lj1,GIMAX,GJMAX
use PhysicalConstants_ml, only : GRAV, ATWAIR 
use SmallUtils_ml,    only: find_index
use TimeDate_ml,      only: date, current_date,day_of_week
use Timefactors_ml,   only: &
    DegreeDayFactors       & ! degree-days used for SNAP-2
    ,Gridded_SNAP2_Factors, gridfac_HDD & 
    ,fac_min,timefactors   &                  ! subroutine
    ,fac_ehh24x7 ,fac_emm, fac_edd, timefac  ! time-factors

implicit none
private

public  :: init_uEMEP
public  :: out_uEMEP
public  :: av_uEMEP
public  :: uemep_adv_x
public  :: uemep_adv_y
public  :: uemep_adv_k
public  :: uemep_emis

real, private, save ::av_fac_hour,av_fac_day,av_fac_month,av_fac_full

contains
subroutine init_uEMEP
  integer :: i, ix, itot, iqrc, iem

!should be set through config
  
  uEMEP%IOU_wanted = .false.
  uEMEP%IOU_wanted(IOU_YEAR) = .true.
!  uEMEP%IOU_wanted(IOU_DAY) = .true.
  uEMEP%IOU_wanted(IOU_HOUR_INST) = .true.

  uEMEP%Nsec_poll = NSECTORS+1
  uEMEP%dist = 5
  uEMEP%Nvert =7
  uEMEP%DOMAIN = RUNDOMAIN
!  uEMEP%DOMAIN = [280,430,130,340]
!  uEMEP%DOMAIN = [35,90,2,75]!ECCAA 


  iem=find_index(uEMEP%emis ,EMIS_FILE(1:NEMIS_FILE))
  call CheckStop( iem<1, "uEMEP did not find corresponding emission file: "//trim(uEMEP%emis) )
  uEMEP%Nix=emis_nsplit(iem)
  do i=1,uEMEP%Nix
     iqrc=sum(emis_nsplit(1:iem-1)) + i
     itot=iqrc2itot(iqrc)
     ix=itot-NSPEC_SHL
     uEMEP%ix(i)=ix
     uEMEP%mw(i)=species_adv(ix)%molwt
     if(uEMEP%emis=="nox ")then
        ix=find_index("NO2",species_adv(:)%name)
        call CheckStop(ix<0,'Index for NO2 not found')
        uEMEP%mw(i)=species_adv(ix)%molwt
     endif
     if(uEMEP%emis=="sox ")then
        ix=find_index("SO2",species_adv(:)%name)
        call CheckStop(ix<0,'Index for SO2 not found')
        uEMEP%mw(i)=species_adv(ix)%molwt
     endif
  end do

  if(MasterProc)then
     write(*,*)'uEMEP Nsec_poll: ',uEMEP%Nsec_poll
     write(*,*)'uEMEP emission file: ',uEMEP%emis
     write(*,*)'uEMEP number of species in group: ',uEMEP%Nix
     write(*,"(A,30(A,F6.2))")'including:',('; '//trim(species_adv(uEMEP%ix(i))%name)//', mw ',uEMEP%mw(i),i=1,uEMEP%Nix)
  end if
  
  av_fac_hour=0.0
  av_fac_day=0.0
  av_fac_month=0.0
  av_fac_full=0.0
  
  allocate(loc_frac(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
  loc_frac=0.0 !must be initiated to 0 so that outer frame does not contribute.
  allocate(loc_frac_ext(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,0:LIMAX+1,0:LJMAX+1))        
  if(uEMEP%IOU_wanted(IOU_HOUR))then
     allocate(loc_frac_hour(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_frac_hour=0.0
     allocate(loc_tot_hour(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_tot_hour=0.0
  endif  
  if(uEMEP%IOU_wanted(IOU_HOUR_INST))then
     allocate(loc_tot_hour(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_tot_hour=0.0
  endif  
  if(uEMEP%IOU_wanted(IOU_DAY))then
     allocate(loc_frac_day(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_frac_day=0.0
     allocate(loc_tot_day(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_tot_day=0.0
  endif
  if(uEMEP%IOU_wanted(IOU_MON))then
     allocate(loc_frac_month(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_frac_month=0.0
     allocate(loc_tot_month(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_tot_month=0.0
  endif
  if(uEMEP%IOU_wanted(IOU_YEAR))then
     allocate(loc_frac_full(uEMEP%Nsec_poll,-uEMEP%dist:uEMEP%dist,-uEMEP%dist:uEMEP%dist,LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_frac_full=0.0
     allocate(loc_tot_full(LIMAX,LJMAX,KMAX_MID-uEMEP%Nvert+1:KMAX_MID))
     loc_tot_full=0.0
  endif
  
end subroutine init_uEMEP


subroutine out_uEMEP(iotyp)
  integer, intent(in) :: iotyp
  character(len=200) ::filename, varname
  real :: xtot,scale,invtot
  integer ::i,j,k,dx,dy,ix,iix
  integer ::isec_poll,ndim,kmax,CDFtype,dimSizes(10),chunksizes(10)
  integer ::ndim_tot,dimSizes_tot(10),chunksizes_tot(10)
  character (len=20) ::dimNames(10),dimNames_tot(10)
  type(Deriv) :: def1 ! definition of fields
  logical,save :: first_call(10)=.true.

  if(.not. uEMEP%IOU_wanted(iotyp))return

  if(iotyp==IOU_HOUR)then
     fileName='uEMEP_hour.nc'
!     varName='local_fraction'
  else if(iotyp==IOU_HOUR_INST)then
     fileName='uEMEP_hour_inst.nc'
!     varName='local_fraction'
  else if(iotyp==IOU_DAY)then
     fileName='uEMEP_day.nc'
!     varName='local_fraction'
  else if(iotyp==IOU_MON)then
     fileName='uEMEP_month.nc'
!     varName='local_fraction'
  else if(iotyp==IOU_YEAR)then
     fileName='uEMEP_full.nc'
!     varName='local_fraction'
  else
     if(me==0)write(*,*)'IOU not recognized'
  endif
  ndim=6
  ndim_tot=3
  kmax=uEMEP%Nvert
  scale=1.0
  CDFtype=Real4
  dimSizes(1)=uEMEP%Nsec_poll
  dimNames(1)='sector'
  dimSizes(2)=2*uEMEP%dist+1
  dimNames(2)='x_dist'
  dimSizes(3)=2*uEMEP%dist+1
  dimNames(3)='y_dist'
  dimSizes(4)=LIMAX
  dimSizes(5)=LJMAX

  dimSizes_tot(1)=LIMAX
  dimSizes_tot(2)=LJMAX

  select case(projection)
  case('Stereographic')
     dimNames(4)='i'
     dimNames(5)='j'      
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'      
  case('lon lat')
     dimNames(4)='lon'
     dimNames(5)='lat'
     dimNames_tot(1)='lon'
     dimNames_tot(2)='lat'      
  case('Rotated_Spherical')
     dimNames(4)='i'
     dimNames(5)='j'      
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'      
  case('lambert')
     dimNames(4)='i'
     dimNames(5)='j'      
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'      
  case default
     dimNames(4)='i'
     dimNames(5)='j'      
     dimNames_tot(1)='i'
     dimNames_tot(2)='j'      
  end select

  dimSizes(6)=kmax
  dimNames(6)='klevel'
  dimSizes_tot(3)=kmax
  dimNames_tot(3)='klevel'
  def1%class='uEMEP' !written
  def1%avg=.false.      !not used
  def1%index=0          !not used
  def1%scale=1.0      !not used
  def1%name=trim(varName)
  def1%unit=''
  chunksizes=1
  chunksizes(4)=dimSizes(4)
  chunksizes(5)=dimSizes(5)
  chunksizes(6)=dimSizes(6)
  chunksizes_tot=1
  chunksizes_tot(1)=dimSizes_tot(1)
  chunksizes_tot(2)=dimSizes_tot(2)
  chunksizes_tot(3)=dimSizes_tot(3)

  if(first_call(iotyp))then
     def1%name='local_fraction'
     call Out_netCDF_n(iotyp,def1,ndim,kmax,loc_frac_full,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=.true.,create_var_only=.true.,chunksizes=chunksizes)
      def1%name='tot_pollutant'
     call Out_netCDF_n(iotyp,def1,ndim_tot,kmax,loc_tot_full,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=.false.,create_var_only=.true.,chunksizes=chunksizes_tot)
 endif

  if(iotyp==IOU_HOUR .or. iotyp==IOU_HOUR_INST)then
     if(iotyp==IOU_HOUR)then
     do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
        do j=1,ljmax
           do i=1,limax
              invtot=1.0/(loc_tot_hour(i,j,k)+1.E-20)
              do dy=-uEMEP%dist,uEMEP%dist
                 do dx=-uEMEP%dist,uEMEP%dist
                    do isec_poll=1,uEMEP%Nsec_poll
                       loc_frac_hour(isec_poll,dx,dy,i,j,k)=loc_frac_hour(isec_poll,dx,dy,i,j,k)*invtot
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     def1%name='local_fraction'
     call Out_netCDF_n(iotyp,def1,ndim,kmax,loc_frac_hour,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.) 
     loc_tot_hour=loc_tot_hour/av_fac_hour
     av_fac_hour=0.0
     def1%name='tot_pollutant'
     call Out_netCDF_n(iotyp,def1,ndim_tot,kmax,loc_tot_hour,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.) 
     loc_frac_hour=0.0
     loc_tot_hour=0.0
     else
        !IOU_HOUR_INST. No need to average
        
     def1%name='local_fraction'
        call Out_netCDF_n(iotyp,def1,ndim,kmax,loc_frac,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
             fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
      do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
        do j=1,ljmax
           do i=1,limax
              xtot=0.0
              do iix=1,uEMEP%Nix
                 ix=uEMEP%ix(iix)
                 xtot=xtot+(xn_adv(ix,i,j,k)*uEMEP%mw(iix))/ATWAIR&
                      *roa(i,j,k,1)*1.E9 !for ug/m3
                 !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
              end do
              loc_tot_hour(i,j,k)=xtot
           enddo
        enddo
     enddo
     def1%name='tot_pollutant'
        call Out_netCDF_n(iotyp,def1,ndim_tot,kmax,loc_tot_hour,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
             fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
     endif

  else  if(iotyp==IOU_DAY)then

     do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
        do j=1,ljmax
           do i=1,limax
              invtot=1.0/(loc_tot_day(i,j,k)+1.E-20)
              do dy=-uEMEP%dist,uEMEP%dist
                 do dx=-uEMEP%dist,uEMEP%dist
                    do isec_poll=1,uEMEP%Nsec_poll
                       loc_frac_day(isec_poll,dx,dy,i,j,k)=loc_frac_day(isec_poll,dx,dy,i,j,k)*invtot
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     def1%name='local_fraction'
     call Out_netCDF_n(iotyp,def1,ndim,kmax,loc_frac_day,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
     loc_tot_day=loc_tot_day/av_fac_day
     av_fac_day=0.0
     def1%name='tot_pollutant'
     call Out_netCDF_n(iotyp,def1,ndim_tot,kmax,loc_tot_day,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
     loc_frac_day=0.0
     loc_tot_day=0.0

  else  if(iotyp==IOU_MON)then

     do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
        do j=1,ljmax
           do i=1,limax
              invtot=1.0/(loc_tot_month(i,j,k)+1.E-20)
              do dy=-uEMEP%dist,uEMEP%dist
                 do dx=-uEMEP%dist,uEMEP%dist
                    do isec_poll=1,uEMEP%Nsec_poll
                       loc_frac_month(isec_poll,dx,dy,i,j,k)=loc_frac_month(isec_poll,dx,dy,i,j,k)*invtot
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     def1%name='local_fraction'
     call Out_netCDF_n(iotyp,def1,ndim,kmax,loc_frac_month,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
     loc_tot_month=loc_tot_month/av_fac_month
     av_fac_month=0.0
     def1%name='tot_pollutant'
     call Out_netCDF_n(iotyp,def1,ndim_tot,kmax,loc_tot_month,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
          fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
     loc_frac_month=0.0
     loc_tot_month=0.0

  else  if(iotyp==IOU_YEAR)then

     do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
        do j=1,ljmax
           do i=1,limax
              invtot=1.0/(loc_tot_full(i,j,k)+1.E-20)
              do dy=-uEMEP%dist,uEMEP%dist
                 do dx=-uEMEP%dist,uEMEP%dist
                    do isec_poll=1,uEMEP%Nsec_poll
                       loc_frac_full(isec_poll,dx,dy,i,j,k)=loc_frac_full(isec_poll,dx,dy,i,j,k)*invtot
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     def1%name='local_fraction'
     call Out_netCDF_n(iotyp,def1,ndim,kmax,loc_frac_full,scale,CDFtype,dimSizes,dimNames,out_DOMAIN=uEMEP%DOMAIN,&
       fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
     loc_tot_full=loc_tot_full/av_fac_full
     av_fac_full=0.0
     def1%name='tot_pollutant'
     call Out_netCDF_n(iotyp,def1,ndim_tot,kmax,loc_tot_full,scale,CDFtype,dimSizes_tot,dimNames_tot,out_DOMAIN=uEMEP%DOMAIN,&
       fileName_given=trim(fileName),overwrite=.false.,create_var_only=.false.)       
     loc_frac_full=0.0
     loc_tot_full=0.0

  else
     if(me==0)write(*,*)'IOU not recognized'
  endif
  first_call(iotyp)=.false.
end subroutine out_uEMEP

subroutine av_uEMEP(dt,End_of_Day)
  real, intent(in)    :: dt                   ! time-step used in integrations
  logical, intent(in) :: End_of_Day           ! e.g. 6am for EMEP sites
  real :: xtot,scale
  integer ::i,j,k,dx,dy,ix,iix
  integer ::isec_poll,ndim,kmax,CDFtype,dimSizes(10),chunksizes(10),iotyp
  character (len=20) ::dimNames(10)
  type(Deriv) :: def1 ! definition of fields
  logical,save :: first_call=.true.
  
  if(.not. uEMEP%IOU_wanted(IOU_HOUR).and.&
     .not. uEMEP%IOU_wanted(IOU_DAY) .and.&
     .not. uEMEP%IOU_wanted(IOU_MON) .and.&
     .not. uEMEP%IOU_wanted(IOU_YEAR)       )return

  !do the averaging
  do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
     do j=1,ljmax
        do i=1,limax
           xtot=0.0
           do iix=1,uEMEP%Nix
              ix=uEMEP%ix(iix)
              xtot=xtot+(xn_adv(ix,i,j,k)*uEMEP%mw(iix))/ATWAIR&
                   *roa(i,j,k,1)*1.E9 !for ug/m3
              !                   *(dA(k)+dB(k)*ps(i,j,1))/GRAV*1.E6 !for mg/m2
           end do
           if(uEMEP%IOU_wanted(IOU_HOUR))then
              loc_tot_hour(i,j,k)=loc_tot_hour(i,j,k)+xtot
              do dy=-uEMEP%dist,uEMEP%dist
                 do dx=-uEMEP%dist,uEMEP%dist
                    do isec_poll=1,uEMEP%Nsec_poll
                       loc_frac_hour(isec_poll,dx,dy,i,j,k)=loc_frac_hour(isec_poll,dx,dy,i,j,k)+xtot*loc_frac(isec_poll,dx,dy,i,j,k)
                    enddo
                 enddo
              enddo
           else if(uEMEP%IOU_wanted(IOU_DAY))then
              loc_tot_day(i,j,k)=loc_tot_day(i,j,k)+xtot
              do dy=-uEMEP%dist,uEMEP%dist
                 do dx=-uEMEP%dist,uEMEP%dist
                    do isec_poll=1,uEMEP%Nsec_poll
                       loc_frac_day(isec_poll,dx,dy,i,j,k)=loc_frac_day(isec_poll,dx,dy,i,j,k)+xtot*loc_frac(isec_poll,dx,dy,i,j,k)
                    enddo
                 enddo
              enddo
           else if(uEMEP%IOU_wanted(IOU_MON))then
              loc_tot_month(i,j,k)=loc_tot_month(i,j,k)+xtot
              do dy=-uEMEP%dist,uEMEP%dist
                 do dx=-uEMEP%dist,uEMEP%dist
                    do isec_poll=1,uEMEP%Nsec_poll
                       loc_frac_month(isec_poll,dx,dy,i,j,k)=loc_frac_month(isec_poll,dx,dy,i,j,k)+xtot*loc_frac(isec_poll,dx,dy,i,j,k)
                    enddo
                 enddo
              enddo
           else if(uEMEP%IOU_wanted(IOU_YEAR))then
              loc_tot_full(i,j,k)=loc_tot_full(i,j,k)+xtot
              do dy=-uEMEP%dist,uEMEP%dist
                 do dx=-uEMEP%dist,uEMEP%dist
                    do isec_poll=1,uEMEP%Nsec_poll
                       loc_frac_full(isec_poll,dx,dy,i,j,k)=loc_frac_full(isec_poll,dx,dy,i,j,k)+xtot*loc_frac(isec_poll,dx,dy,i,j,k)
                  enddo
                 enddo
              enddo
!              if(i==5.and.j==5.and.k==KMAX_MID)write(*,*)me,xtot,loc_tot_full(i,j,k),loc_frac_full(1,0,0,i,j,k)
           endif
        enddo
     enddo
  enddo
  av_fac_hour=av_fac_hour+1.0
  av_fac_day=av_fac_day+1.0
  av_fac_month=av_fac_month+1.0
  av_fac_full=av_fac_full+1.0

end subroutine av_uEMEP

  subroutine uemep_adv_y(fluxy,i,j,k)
    real, intent(in)::fluxy(NSPEC_ADV,-1:LJMAX+1)
    integer, intent(in)::i,j,k
    real ::x,xn,xx,f_in,inv_tot
    integer ::iix,ix,dx,dy,isec_poll
    xn=0.0
    x=0.0
    xx=0.0
    !positive x or xx means incoming, negative means outgoing
    do iix=1,uEMEP%Nix
       ix=uEMEP%ix(iix)
       xn=xn+xn_adv(ix,i,j,k)*uEMEP%mw(iix)
       x=x-xm2(i,j)*fluxy(ix,j)*uEMEP%mw(iix)!flux through "North" face (Up)
       xx=xx+xm2(i,j)*fluxy(ix,j-1)*uEMEP%mw(iix)!flux through "South" face (Bottom)
    end do
    !NB: here xn already includes the fluxes. Remove them!
    xn=xn-xx-x

    xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux 
    f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
    inv_tot=1.0/(xn+f_in+1.e-20)!incoming dilutes

    xx=max(0.0,xx)*inv_tot!factor due to flux through "South" face (Bottom)
    x =max(0.0,x)*inv_tot!factor due to flux through "North" face (Up)

    do dy=-uEMEP%dist,uEMEP%dist
       do dx=-uEMEP%dist,uEMEP%dist
!        if(k==KMAX_MID.and.i==5.and.j==5.and.)write(*,*)'B ',me,loc_frac(0,1,i,j,k),xn/(xn+f_in+1.e-20)
          do isec_poll=1,uEMEP%Nsec_poll
             loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k) *xn *inv_tot
          enddo
!          if(dx==0 .and. dy==0)cycle!temporary: no return of pollutants
          if(x>0.0.and.dy>-uEMEP%dist)then
             do isec_poll=1,uEMEP%Nsec_poll
                loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)+ loc_frac_ext(isec_poll,dx,dy-1,i,j+1)*x
             enddo
          endif
          if(xx>0.0.and.dy<uEMEP%dist)then
             do isec_poll=1,uEMEP%Nsec_poll
                loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)+ loc_frac_ext(isec_poll,dx,dy+1,i,j-1)*xx
             enddo
          endif
       enddo
    enddo

  end subroutine uemep_adv_y

  subroutine uemep_adv_x(fluxx,i,j,k)
    real, intent(in)::fluxx(NSPEC_ADV,-1:LIMAX+1)
    integer, intent(in)::i,j,k
    real ::x,xn,xx,f_in,inv_tot
    integer ::iix,ix,dx,dy,isec_poll
    xn=0.0
    x=0.0
    xx=0.0
    !positive x or xx means incoming, negative means outgoing
    do iix=1,uEMEP%Nix
       ix=uEMEP%ix(iix)
       xn=xn+xn_adv(ix,i,j,k)*uEMEP%mw(iix)
       x=x-xm2(i,j)*fluxx(ix,i)*uEMEP%mw(iix)!flux through "East" face (Right)
       xx=xx+xm2(i,j)*fluxx(ix,i-1)*uEMEP%mw(iix)!flux through "West" face (Left)
    end do
    !NB: here xn already includes the fluxes. Remove them!
    xn=xn-xx-x

    xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. all outgoing flux 
    f_in=max(0.0,x)+max(0.0,xx)!positive part. all incoming flux
    inv_tot=1.0/(xn+f_in+1.e-20)!incoming dilutes

    x =max(0.0,x)*inv_tot!factor due to flux through "East" face (Right)
    xx=max(0.0,xx)*inv_tot!factor due to flux through "West" face (Left)

    do dy=-uEMEP%dist,uEMEP%dist
       do dx=-uEMEP%dist,uEMEP%dist
          do isec_poll=1,uEMEP%Nsec_poll
             loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k) *xn *inv_tot
          enddo
!          if(dx==0 .and. dy==0)cycle
          if(x>0.0.and.dx>-uEMEP%dist)then
             do isec_poll=1,uEMEP%Nsec_poll
                loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)+ loc_frac_ext(isec_poll,dx-1,dy,i+1,j)*x
             enddo
          endif
          if(xx>0.0.and.dx<uEMEP%dist)then
             do isec_poll=1,uEMEP%Nsec_poll
                loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)+ loc_frac_ext(isec_poll,dx+1,dy,i-1,j)*xx
             enddo
          endif
       enddo
    enddo
 
  end subroutine uemep_adv_x

  subroutine uemep_adv_k(fluxk,i,j)
    real, intent(in)::fluxk(NSPEC_ADV,KMAX_MID)
    integer, intent(in)::i,j
    real ::x,xn,xx,f_in,inv_tot
    integer ::k,iix,ix,dx,dy,isec_poll

    do k = KMAX_MID-uEMEP%Nvert+1,KMAX_MID
       xn=0.0
       x=0.0
       xx=0.0
       do iix=1,uEMEP%Nix
          ix=uEMEP%ix(iix)
          xn=xn+xn_adv(ix,i,j,k)*uEMEP%mw(iix)
          if(k<KMAX_MID)x=x-dhs1i(k+1)*fluxk(ix,k+1)*uEMEP%mw(iix)
          xx=xx+dhs1i(k+1)*fluxk(ix,k)*uEMEP%mw(iix)
       end do
       !NB: here xn already includes the fluxes. Remove them!
       xn=xn-xx-x

       xn=max(0.0,xn+min(0.0,x)+min(0.0,xx))!include negative part. outgoing flux 
       f_in=max(0.0,x)+max(0.0,xx)!positive part. incoming flux
       inv_tot = 1.0/(xn+f_in+1.e-20)

       if(k==KMAX_MID-uEMEP%Nvert+1)then
          !highest level for uemep. Assume zero local fractions coming from above
          xx=0.0
       endif

       x =max(0.0,x)*inv_tot!factor due to flux through bottom face
       xx=max(0.0,xx)*inv_tot!factor due to flux through top face
       if(k<KMAX_MID)then
          if(k>KMAX_MID-uEMEP%Nvert+1)then
             do dy=-uEMEP%dist,uEMEP%dist
                do dx=-uEMEP%dist,uEMEP%dist
                   do isec_poll=1,uEMEP%Nsec_poll
                      loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k) * xn * inv_tot &
                       + loc_frac(isec_poll,dx,dy,i,j,k-1) * xx&
                       + loc_frac(isec_poll,dx,dy,i,j,k+1) * x
                   enddo
                enddo
             enddo
          else
             !k=KMAX_MID-uEMEP%Nvert+1 , assume no local fractions from above
             do dy=-uEMEP%dist,uEMEP%dist
                do dx=-uEMEP%dist,uEMEP%dist
                   do isec_poll=1,uEMEP%Nsec_poll
                      loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)* xn * inv_tot &
                       + loc_frac(isec_poll,dx,dy,i,j,k+1) * x
                   enddo
                enddo
             enddo
          endif
       else
          !k=KMAX_MID , no local fractions from below
          do dy=-uEMEP%dist,uEMEP%dist
             do dx=-uEMEP%dist,uEMEP%dist
                do isec_poll=1,uEMEP%Nsec_poll
                   loc_frac(isec_poll,dx,dy,i,j,k) = loc_frac(isec_poll,dx,dy,i,j,k)* xn * inv_tot &
                       + loc_frac(isec_poll,dx,dy,i,j,k-1) * xx
                enddo
             enddo
          enddo
       endif
    end do
  end subroutine uemep_adv_k


subroutine uEMEP_emis(indate)
!include emission contributions to local fractions

!NB: should replace most of the stuff and use gridrcemis instead!

  implicit none
  type(date), intent(in) :: indate  ! Gives year..seconds
  integer :: i, j, k          ! coordinates, loop variables
  integer :: icc, ncc         ! No. of countries in grid.
  integer :: ficc,fncc        ! No. of countries with
  integer :: iqrc             ! emis indices 
  integer :: isec             ! loop variables: emission sectors
  integer :: iem              ! loop variable over 1..NEMIS_FILE
  integer :: itot             ! index in xn()

  ! Save daytime value between calls, initialise to zero
  integer, save, dimension(MAXNLAND) ::  daytime(1:MAXNLAND) = 0  !  0=night, 1=day
  integer, save, dimension(MAXNLAND) ::  localhour(1:MAXNLAND) = 1  ! 1-24 local hour in the different countries
  integer                         ::  hourloc      !  local hour 
  real, dimension(NRCEMIS)        ::  tmpemis      !  local array for emissions
  real ::  tfac    ! time-factor (tmp variable); dt*h*h for scaling
  real ::  s       ! source term (emis) before splitting
  integer :: iland, iland_timefac  ! country codes, and codes for timefac 
  integer :: daytime_longitude, daytime_iland, hour_longitude, hour_iland,nstart
  integer ::icc_uemep
  integer, save :: wday , wday_loc ! wday = day of the week 1-7
  integer ::ix,iix, neigbor, dx, dy, isec_poll
  real::dt_uemep, xtot, emis_uemep(KMAX_MID,NSECTORS),emis_tot(KMAX_MID)
  logical,save :: first_call=.true.  

  dt_uemep=dt_advec

  wday=day_of_week(indate%year,indate%month,indate%day)
  if(wday==0)wday=7 ! Sunday -> 7
  do iland = 1, NLAND
    daytime(iland) = 0
    hourloc        = indate%hour + Country(iland)%timezone
    localhour(iland) = hourloc  ! here from 0 to 23
    if(hourloc>=7 .and. hourloc<=18) daytime(iland)=1
  end do ! iland
   
  do j = lj0,lj1
    do i = li0,li1
      ncc = nlandcode(i,j)            ! No. of countries in grid
      fncc = flat_nlandcode(i,j) ! No. of countries with flat emissions in grid
      hourloc= mod(nint(indate%hour+24*(1+glon(i,j)/360.0)),24)
      hour_longitude=hourloc
      daytime_longitude=0
      if(hourloc>=7 .and. hourloc<= 18) daytime_longitude=1            
      !*************************************************
      ! First loop over non-flat (one sector) emissions
      !*************************************************
      tmpemis(:)=0.
      icc_uemep=0
      emis_uemep=0.0
      emis_tot=0.0
      do icc = 1, ncc+fncc
        ficc=icc-ncc
        !          iland = landcode(i,j,icc)     ! 1=Albania, etc.
        if(icc<=ncc)then
          iland=find_index(landcode(i,j,icc),Country(:)%icode) !array index
        else
          iland=find_index(flat_landcode(i,j,ficc),Country(:)%icode) 
        end if
        !array index of country that should be used as reference for timefactor
        iland_timefac = find_index(Country(iland)%timefac_index,Country(:)%timefac_index)

        if(Country(iland)%timezone==-100)then
          daytime_iland=daytime_longitude
          hour_iland=hour_longitude + 1   ! add 1 to get 1..24 
        else
          daytime_iland=daytime(iland)
          hour_iland=localhour(iland) + 1
        end if
        !if( hour_iland > 24 ) hour_iland = 1 !DSA12
        wday_loc=wday 
        if(hour_iland>24) then
          hour_iland = hour_iland - 24
          wday_loc=wday + 1
          if(wday_loc==0)wday_loc=7 ! Sunday -> 7
          if(wday_loc>7 )wday_loc=1 
        end if

        do iem = 1, NEMIS_FILE 
          if(trim(EMIS_File(iem))/=trim(uEMEP%emis))cycle
          do isec = 1, NSECTORS       ! Loop over snap codes
            ! Calculate emission rates from snapemis, time-factors, 
            ! and if appropriate any speciation fraction (NEMIS_FRAC)
            iqrc = 0   ! index over emisfrac
            ! kg/m2/s
            
            if(icc<=ncc)then
              tfac = timefac(iland_timefac,sec2tfac_map(isec),iem) &
                   * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc)

              !Degree days - only SNAP-2 
              if(USES%DEGREEDAY_FACTORS .and. &
                   sec2tfac_map(isec)==ISNAP_DOM .and. Gridded_SNAP2_Factors) then
                 ! If INERIS_SNAP2  set, the fac_min will be zero, otherwise
                 ! we make use of a baseload even for SNAP2
                 tfac = ( fac_min(iland,sec2tfac_map(isec),iem) & ! constant baseload
                      + ( 1.0-fac_min(iland,sec2tfac_map(isec),iem) )* gridfac_HDD(i,j) ) &
                      * fac_ehh24x7(sec2tfac_map(isec),hour_iland,wday_loc)
              end if ! =============== HDD 
              
              s = tfac * snapemis(isec,i,j,icc,iem)
            else
              s = snapemis_flat(i,j,ficc,iem)                        
            end if

            do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
              emis_tot(k)=emis_tot(k)+s*emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))*dt_uemep
            end do

            !if(isec==uEMEP%sector .or. uEMEP%sector==0)then
              do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
                emis_uemep(k,isec)=emis_uemep(k,isec)+s*emis_kprofile(KMAX_BND-k,sec2hfac_map(isec))*dt_uemep
              end do
            !end if

          end do ! iem

        end do  ! isec
        !      ==================================================
      end do ! icc  
      
      do k=max(KEMISTOP,KMAX_MID-uEMEP%Nvert+1),KMAX_MID
         if(emis_tot(k)<1.E-20)cycle
         !units kg/m2
         !total pollutant
         xtot=0.0
         do iix=1,uEMEP%Nix
            ix=uEMEP%ix(iix)
            xtot=xtot+(xn_adv(ix,i,j,k)*uEMEP%mw(iix))*(dA(k)+dB(k)*ps(i,j,1))/ATWAIR/GRAV
         end do
         neigbor = 1!local fraction from this i,j
         dx=0 ; dy=0!local fraction from this i,j
         do isec_poll=2,uEMEP%Nsec_poll
            loc_frac(isec_poll,dx,dy,i,j,k)=(loc_frac(isec_poll,dx,dy,i,j,k)*xtot+emis_uemep(k,isec_poll-1))/(xtot+emis_tot(k)+1.e-20)
         enddo

         isec = 0 !sum over all sectors
         isec_poll=1 !sum over all sectors
         loc_frac(isec_poll,dx,dy,i,j,k)=(loc_frac(isec_poll,dx,dy,i,j,k)*xtot+emis_tot(k))/(xtot+emis_tot(k)+1.e-20)

         !local fractions from other cells

         do dy=-uEMEP%dist,uEMEP%dist
            do dx=-uEMEP%dist,uEMEP%dist
               if(dx==0 .and. dy==0)cycle!local fractions from other cells only
               do isec_poll=1,uEMEP%Nsec_poll
                  loc_frac(isec_poll,dx,dy,i,j,k)=(loc_frac(isec_poll,dx,dy,i,j,k)*xtot)/(xtot+emis_tot(k)+1.e-20)
               enddo
            enddo
         enddo

      end do! k

    end do ! i
  end do ! j

  first_call=.false. 

end subroutine uEMEP_emis
end module uEMEP_ml
