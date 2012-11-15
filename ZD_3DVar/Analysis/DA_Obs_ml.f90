module DA_Obs_ml
use Par_ml,           only: me
use MetFields_ml,     only: z_bnd,z_mid,roa,th,q,ps
use TimeDate_ml,      only: current_date
use CheckStop_ml,     only: CheckStop
use Functions_ml,     only: great_circle_distance
use GridValues_ml,    only: glon,glat,coord_in_processor
use Par_ml,           only: limax,ljmax
use Chemfields_ml,    only: xn_adv,cfac,PM25_water,PM25_water_rh50
!se ChemSpecs_shl_ml, only: NSPEC_SHL        ! Maps indices
use ChemSpecs_adv_ml, only: NSPEC_ADV
use ChemChemicals_ml, only: species
use ChemGroups_ml,    only: chemgroups
use ModelConstants_ml,only: KMAX_MID
use SmallUtils_ml,    only: find_index
use TimeDate_ExtraUtil_ml, only: date2string
use Units_ml,         only: Units_Scale,Group_Units
use DA_ml,            only: debug=>DA_DEBUG,dafmt=>da_fmt_msg,damsg=>da_msg
use Util_ml,          only: norm
use spectralcov,      only: nx,ny,nlev,nchem,nxex,nyex,FGSCALE,&
                            nChemObs,nChemNoObs,iChemInv,iChemObs,iChemNoObs
implicit none
integer, save :: nobs!,matidx
! integer, dimension(:), allocatable, save :: pidx
real, dimension(:), allocatable, save :: innov,obsstddev
real, dimension(:,:,:,:), allocatable, save :: H_jac
type H_jac_idx
  integer :: i(4),j(4),l(0:1)
end type
type(H_jac_idx), dimension(:), allocatable, save :: pidx
!-----------------------------------------------------------------------
type obs_data
  logical            :: set=.false.,found=.false.,unitroa=.false.
  integer            :: ispec=-1,ichem=-1,ichemObs=-1
  real               :: min=0.0,max=1e3,error_rel=-1e0,error_rep=-1e0,unitconv=-1e0
  character(len=064) :: name='',unit='',deriv=''
  character(len=256) :: file='',tag=''
end type
integer            :: nobsData
integer, parameter :: nobsDataMax=1000
type(obs_data)     :: obsData(nobsDataMax)
namelist /OBSERVATIONS/nobsData,obsData
!-----------------------------------------------------------------------
real, dimension(:,:,:,:), allocatable :: xn_adv_ex
integer               :: varSpec(NSPEC_ADV),varSpecInv(NSPEC_ADV)
character(len=016)    :: varName(NSPEC_ADV)='',obsVarName(NSPEC_ADV)=''
logical               :: observedVar(NSPEC_ADV)=.false.
contains
subroutine allocate_obs()
!-----------------------------------------------------------------------
! @description
! Allocate the arrays innov and obsstddev. This step is
! required, since the number of observations can vary from
! one time step to the next
! @author AMVB
!-----------------------------------------------------------------------
  integer :: ierr
! no4=nobs*4
! no8=nobs*8
! no12=nobs*12
! matidx=0
  if(.not.allocated(innov))then
    allocate(innov(nobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: INNOV.')
    innov=0.0
  endif
  if(.not.allocated(obsstddev))then
    allocate(obsstddev(nobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: OBSSTDDEV.')
    obsstddev=0.0
  endif
!-----------------------------------------------------------------------
! H_jac is a matrix with components H_{n;i,j,l,k}, where
!   n=1,...,nobs     (observation index)
!   i=1,...,nx       (longitude index)
!   j=1,...,ny       (latitude index)
!   l=1,...,nlev     (altitude index)
!   k=1,...,nchemobs (observed chemical component)
! Due to the 4-point interpolation in the horizontal direction
! of the model results to the observation point, the Jacobian is
! highly sparse wrt i,j,l. For this reason, the matrix is stored using
! the following compact storage scheme:
! H_{n;i,j,l,k}=H_jac(n,p,l,k), where n=1:nobs,p=1:4,l=l(0):l(1)
!   pidx(n)%i(1)=i_n,   pidx(n)%j(1)=j_n,   pidx(n)%l(0)=l0_n,
!   pidx(n)%i(2)=i_n+1, pidx(n)%j(2)=j_n,   pidx(n)%l(1)=l1_n,
!   pidx(n)%i(3)=i_n,   pidx(n)%j(3)=j_n+1,
!   pidx(n)%i(4)=i_n+1, pidx(n)%j(4)=j_n+1,
!-----------------------------------------------------------------------
  if(.not.allocated(H_jac))then
   !allocate(H_jac(nobs,4*nobs,nchem),stat=ierr)
   !allocate(H_jac(nobs,4*nobs,nlev,nchem),stat=ierr)
    allocate(H_jac(nobs,4,nlev,nchemobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: H_JAC.')
    H_jac=0.0
  endif
  if(.not.allocated(pidx))then
   !allocate(pidx(12*nobs),stat=ierr)
    allocate(pidx(nobs),stat=ierr)
    call CheckStop(ierr,'Allocation error: PIDX.')
   !pidx=0
    pidx=H_jac_idx((/0,0,0,0/),(/0,0,0,0/),(/0,0/))
  endif
end subroutine allocate_obs
subroutine deallocate_obs()
!-----------------------------------------------------------------------
! @description
! Deallocate the arrays innov and obsstddev. This step is
! required, since the number of observations can vary from
! one time step to the next
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none
  if(allocated(innov))    deallocate(innov)
  if(allocated(obsstddev))deallocate(obsstddev)
  if(allocated(H_jac))    deallocate(H_jac)
  if(allocated(pidx))     deallocate(pidx)
endsubroutine deallocate_obs
subroutine read_obs(maxobs,flat,flon,falt,y,stddev,ipar,iostat)
!-----------------------------------------------------------------------
! @description
! Read in observations y in standard format
! (Station-No flat flon falt obs-value stddev).
! @author M.Kahnert
!-----------------------------------------------------------------------
  implicit none

  integer maxobs,iostat,ipar(maxobs)
  real y(maxobs),stddev(maxobs)
  real flat(maxobs),flon(maxobs),falt(maxobs)

  integer :: lu1=20
  integer no,iobsData
  character(len=256) file

! #     include "DATAASSIM.INC"
! #     include "MPP.INC"
! #     include "CONST.INC"
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  nobs=0
!-----------------------------------------------------------------------
! open observational data file in standard format, if applicable
!-----------------------------------------------------------------------
  do iobsData=1,nobsData
    if (obsData(iobsData)%file/="") then
      iostat=0
      file=date2string(obsData(iobsData)%file,current_date)
      open(lu1,file=file,form='formatted',status='old',iostat=iostat)
      if(iostat/=0)then
        print dafmt,"Obsdata not opened "//trim(file)
        close(lu1)
        iostat=-1
        cycle
      else
        print dafmt,"Obsdata opened "//trim(file)
      endif
    else
      print dafmt,'No observation file specified'
      return
    endif
!-----------------------------------------------------------------------
! read data
!-----------------------------------------------------------------------
    do while(iostat==0)
      if(nobs>=maxobs)then
        print dafmt,'ERROR reading observations, maxobs too small'
        iostat=-2
        cycle
      endif
      nobs=nobs+1
      ipar(nobs)=iobsData
      read(lu1,*,iostat=iostat)no,flat(nobs),flon(nobs),falt(nobs),y(nobs),stddev(nobs)
      select case(iostat)
      case(0)
        if((.not.coord_in_processor(flon(nobs),flat(nobs))).or.&
           (y(nobs)<obsData(iobsData)%min).or.&
           (y(nobs)>obsData(iobsData)%max))then
         nobs=nobs-1
         cycle
        endif
        stddev(nobs)=max(stddev(nobs),&
          y(nobs)*obsData(iobsData)%error_rel,&
          obsData(iobsData)%error_rep)
        if(stddev(nobs)<=0e0)then
          print dafmt,'WARNING obs stddev <= 0'
          stddev(nobs)=1e-9
        endif
!       if(obsData(iobsData)%scale<0e0)then ! initialize scale
! !         print dafmt,'WARNING obs scale <= 0'
!         obsData(iobsData)%scale=1e0/Units_Scale(obsData(iobsData)%unit,obsData(iobsData)%ispec) !unit conversion factor
!       endif
!       if(obsData(iobsData)%scale/=1e0)then
!         y     (nobs)=y     (nobs)*obsData(iobsData)%scale
!         stddev(nobs)=stddev(nobs)*obsData(iobsData)%scale
!       endif
      case(:-1)     !reached EOF: nothing readed
        nobs=nobs-1
      case(1:)
        print dafmt,'ERROR reading observations'
        nobs=nobs-1
        iostat=-1
      endselect
    enddo
    close(lu1)
  enddo
end subroutine read_obs
subroutine H_op(lat,lon,alt,n,yn,rho0,ipar)
!-----------------------------------------------------------------------
! @description
! H operator mapping from model to observation space;
! Further: interpolate air density to observation point
! @author AMVB
!-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: n,ipar
  real, intent(in)    :: lat,lon,alt
  real, intent(inout) :: yn
  real, intent(out)   :: rho0

  integer i(4),j(4),p,l,ll,k,kk,igrp,g,ispec
  real xn,Hj,H_int(4),T1,QML1,PS1,Hchem(nchem),unitconv
  integer, pointer, dimension(:) :: gspec=>null()      ! group array of indexes
  real,    pointer, dimension(:) :: gunit_conv=>null() ! group array of unit conv. factors
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! no4=4*nobs
! no8=8*nobs
!-----------------------------------------------------------------------
! operator for horizontal four-point interpolation between nearest
! neighbouring grid points:
!-----------------------------------------------------------------------
  call get_interpol_const(lon,lat,i,j,H_int)
!-----------------------------------------------------------------------
! z grid coordinate of observations:
!-----------------------------------------------------------------------
  ll=KMAX_MID
  do while(alt>z_bnd(i(1),j(1),ll).and.ll>=1)
    ll=ll-1
  enddo
  l=ll+nlev-KMAX_MID  ! B might have a reduced number of levels...
  if(debug.and.me==0) then
    print dafmt,'Observation/Model Location:Weight[%],lon[deg],lat[deg],alt[km]'
    print "('#',I0,5(1X,A3,':',F0.1,3(',',F0.2)))",&
    n,'Observation',100.0,lon,lat,alt,&
     ('Model',H_int(p)*100,glon(i(p),j(p)),glat(i(p),j(p)),z_mid(i(p),j(p),ll)*1e-3,p=1,4)
  endif
! interpolated air density:
  rho0=dot_product(H_int(:),&
       (/roa(i(1),j(1),ll,1),roa(i(2),j(2),ll,1),&
         roa(i(3),j(3),ll,1),roa(i(4),j(4),ll,1)/))
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
  if(.not.obsData(ipar)%set)then
    k=find_index(obsData(ipar)%name,obsVarName(:))
    if(k>0)then
      obsData(ipar)%found=.true.
      obsData(ipar)%ichem=k                     ! index observed
      obsData(ipar)%ichemObs=ichemobs(k)        ! index observed/unobserved
      obsData(ipar)%ispec=varSpec(ichemobs(k))  ! index species
      obsData(ipar)%unitconv=Units_Scale(obsData(ipar)%unit,obsData(ipar)%ispec,&
        needroa=obsData(ipar)%unitroa)
    endif
    if(obsData(ipar)%deriv=='')then
                  obsData(ipar)%deriv='mod-lev'
      if(alt==0.0)obsData(ipar)%deriv='surface'
    endif
    write(obsData(ipar)%tag,"(2(A,1X),'(',A,')')")&
      trim(obsData(ipar)%name),trim(obsData(ipar)%deriv),trim(obsData(ipar)%unit)
    obsData(ipar)%set=.true.
  endif
  print dafmt,"Analysis of "//trim(obsData(ipar)%tag)
!-----------------------------------------------------------------------
! analysis of direct observations: model mid-levels
!-----------------------------------------------------------------------
  if(obsData(ipar)%found.and.obsData(ipar)%deriv=='mod-lev')then
    unitconv=obsData(ipar)%unitconv
    if(obsData(ipar)%unitroa)&
      H_int(:)=H_int(:)*(/roa(i(1),j(1),ll,1),roa(i(2),j(2),ll,1),&
                          roa(i(3),j(3),ll,1),roa(i(4),j(4),ll,1)/)
    k=obsData(ipar)%ichem
    kk=ichemobs(k)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    do p=1,4
      H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
        * unitconv                                    ! unit conversion
      xn=xn_adv_ex(i(p),j(p),l,kk)
      yn=yn+H_jac(n,p,l,k)*xn
    enddo
!-----------------------------------------------------------------------
! analysis of direct observations: surface level
!-----------------------------------------------------------------------
  elseif(obsData(ipar)%found.and.obsData(ipar)%deriv=='surface')then
    unitconv=obsData(ipar)%unitconv
    if(obsData(ipar)%unitroa)&
      H_int(:)=H_int(:)*(/roa(i(1),j(1),ll,1),roa(i(2),j(2),ll,1),&
                          roa(i(3),j(3),ll,1),roa(i(4),j(4),ll,1)/)
    ispec=obsData(ipar)%ispec
    k=obsData(ipar)%ichem
    kk=ichemobs(k)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    do p=1,4
      H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
        * unitconv                                  & ! unit conversion
        * cfac(ispec,i(p),j(p))                       ! 50m->3m conversion
      xn=xn_adv_ex(i(p),j(p),l,kk)
      yn=yn+H_jac(n,p,l,k)*xn
    enddo
!-----------------------------------------------------------------------
! analysis of indirect observations: COLUMN (eg NO2tc)
!-----------------------------------------------------------------------
  elseif(obsData(ipar)%found.and.obsData(ipar)%deriv=='Trop.Col.')then
    unitconv=obsData(ipar)%unitconv
    ispec=obsData(ipar)%ispec
    k=obsData(ipar)%ichem
    kk=ichemobs(k)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=(/1,nlev/)
    do ll=1,KMAX_MID-nlev   ! B might have a reduced number of levels...
      do p=1,4
        Hj=H_int(p)                                   & ! interpolation
          * unitconv                                  & ! unit conversion
          * roa(i(p),j(p),ll,1)                       & ! density.
          * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
        xn=xn_adv(ispec,i(p),j(p),ll)*FGSCALE
        yn=yn+Hj*xn
      enddo
    enddo
    do l=1,nlev
      ll=KMAX_MID-nlev+l
      do p=1,4
        H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
          * unitconv                                  & ! unit conversion
          * roa(i(p),j(p),ll,1)                       & ! density.
          * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
        xn=xn_adv_ex(i(p),j(p),l,kk)
        yn=yn+H_jac(n,p,l,k)*xn
      enddo
    enddo
!-----------------------------------------------------------------------
! analysis of indirect observations: model mid-levels
!  PM* (eg PM10), BSC and EXT GROUPs (wavelength is in the %unit)
!-----------------------------------------------------------------------
  elseif((index(obsData(ipar)%name,'PM')>0.or.&
          any(obsData(ipar)%name==(/'BSC','EXT'/)))&
         .and.obsData(ipar)%deriv=='mod-lev')then
!   call CheckStop(.not.obsData(ipar)%unitroa.or.&
!     (index(obsData(ipar)%name,'PM')>0.and.obsData(ipar)%unit(1:2)/='ug'),&
!     "Unsupported obsData unit for: "//trim(obsData(ipar)%tag))
    igrp=find_index(obsData(ipar)%name,chemgroups(:)%name)
    if(igrp<1)then
      if(me==0) print dafmt,"ERROR obsData not yet supported"
      return
    endif
    call Group_Units(igrp,obsData(ipar)%unit,gspec,gunit_conv,debug)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    do g=1,size(gspec)
      kk=varSpecInv(gspec(g))
      if(kk>0.and.observedVar(kk))then
        k=ichemInv(kk)
        do p=1,4
          H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
            * gunit_conv(g)                             & ! unit conversion
            * roa(i(p),j(p),ll,1)                         ! density.
          xn=xn_adv_ex(i(p),j(p),l,kk)
          yn=yn+H_jac(n,p,l,k)*xn
        enddo
      else
        do p=1,4
          Hj=H_int(p)                                   & ! interpolation
            * gunit_conv(g)                             & ! unit conversion
            * roa(i(p),j(p),ll,1)                         ! density.
          xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
          yn=yn+Hj*xn
        enddo
      endif
    enddo
    if(index(obsData(ipar)%name,'PM')>0)then
      do p=1,4
        yn=yn+PM25_water(i(p),j(p),ll)*FGSCALE ! PM water content in ug/m3 at model mid-levels
      enddo
    endif
!-----------------------------------------------------------------------
! analysis of indirect observations: model mid-levels
!  PM* (eg PM10), BSC and EXT GROUPs (wavelength is in the %unit)
!-----------------------------------------------------------------------
  elseif((index(obsData(ipar)%name,'PM')>0.or.&
          any(obsData(ipar)%name==(/'BSC','EXT'/)))&
         .and.obsData(ipar)%deriv=='surface')then
!   call CheckStop(.not.obsData(ipar)%unitroa.or.&
!     (index(obsData(ipar)%name,'PM')>0.and.obsData(ipar)%unit(1:2)/='ug'),&
!     "Unsupported obsData unit for: "//trim(obsData(ipar)%tag))
    igrp=find_index(obsData(ipar)%name,chemgroups(:)%name)
    if(igrp<1)then
      if(me==0) print dafmt,"ERROR obsData not yet supported"
      return
    endif
    print *,igrp,chemgroups(igrp)%name
    call Group_Units(igrp,obsData(ipar)%unit,gspec,gunit_conv,debug)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    do g=1,size(gspec)
      kk=varSpecInv(gspec(g))
      if(kk>0.and.observedVar(kk))then
        k=ichemInv(kk)
        do p=1,4
          H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
            * gunit_conv(g)                             & ! unit conversion
            * roa(i(p),j(p),ll,1)                       & ! density.
            * cfac(gspec(g),i(p),j(p))                    ! 50m->3m conversion
          xn=xn_adv_ex(i(p),j(p),l,kk)
          yn=yn+H_jac(n,p,l,k)*xn
        enddo
      else
        do p=1,4
          Hj=H_int(p)                                   & ! interpolation
            * gunit_conv(g)                             & ! unit conversion
            * roa(i(p),j(p),ll,1)                       & ! density.
            * cfac(gspec(g),i(p),j(p))                    ! 50m->3m conversion
          xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
          yn=yn+Hj*xn
        enddo
      endif
    enddo
    if(index(obsData(ipar)%name,'PM')>0)then
      do p=1,4
        yn=yn+PM25_water_rh50(i(p),j(p))*FGSCALE ! PM water content in ug/m3 at surface level
      enddo
    endif
!-----------------------------------------------------------------------
! analysis of indirect observations: column
!  AOD GROUP (wavelength is in the %unit)
!-----------------------------------------------------------------------
  elseif((obsData(ipar)%name=='AOD').or.&
         (obsData(ipar)%name=='EXT'.and.obsData(ipar)%deriv=='Trop.Col.'))then
    igrp=find_index(obsData(ipar)%name,chemgroups(:)%name)
    if(igrp<1)then
      if(me==0) print dafmt,"ERROR obsData not yet supported"
      return
    endif
    call Group_Units(igrp,obsData(ipar)%unit,gspec,gunit_conv,debug)
    yn=0e0
    pidx(n)%i(:)=i(:)
    pidx(n)%j(:)=j(:)
    pidx(n)%l(:)=l
    pidx(n)%l(:)=(/1,nlev/)
    do g=1,size(gspec)
      kk=varSpecInv(gspec(g))
      if(kk>0.and.observedVar(kk))then
        k=ichemInv(kk)
        do ll=1,KMAX_MID-nlev
          do p=1,4
            Hj=H_int(p)                                   & ! interpolation
              * gunit_conv(g)                             & ! unit conversion
              * roa(i(p),j(p),ll,1)                       & ! density.
              * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
            xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
            yn=yn+Hj*xn
          enddo
        enddo
        do l=1,nlev
          ll=KMAX_MID-nlev+l
          do p=1,4
            H_jac(n,p,l,k)=H_int(p)                       & ! interpolation
              * gunit_conv(g)                             & ! unit conversion
              * roa(i(p),j(p),ll,1)                       & ! density.
              * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
            xn=xn_adv_ex(i(p),j(p),l,kk)
            yn=yn+H_jac(n,p,l,k)*xn
          enddo
        enddo
      else
        do ll=1,KMAX_MID
          do p=1,4
            Hj=H_int(p)                                   & ! interpolation
              * gunit_conv(g)                             & ! unit conversion
              * roa(i(p),j(p),ll,1)                       & ! density.
              * (z_bnd(i(p),j(p),ll)-z_bnd(i(p),j(p),ll+1)) ! level thickness
            xn=xn_adv(gspec(g),i(p),j(p),ll)*FGSCALE
            yn=yn+Hj*xn
          enddo
        enddo
      endif
    enddo
  else
    if(me==0) print dafmt,"ERROR obsData not yet supported"
    return
  endif
end subroutine H_op
! function great_circle_distance(fi1,lambda1,fi2,lambda2) result(dist)
!   !compute the great circle distance between to points given in
!   !spherical coordinates. Sphere has radius 1.
!   real, intent(in) ::fi1,lambda1,fi2,lambda2 !NB: in DEGREES here
!   real :: dist
!   dist=2*asin(sqrt(sind(0.5*(lambda1-lambda2+360.0))**2+&
!         cosd(lambda1+360.0)*cosd(lambda2+360.0)*sind(0.5*(fi1-fi2+360.0))**2))
! end function great_circle_distance
subroutine get_interpol_const(lon,lat,IIij,JJij,Weight)
  !find interpolation constants
  !note that i,j are local
  !find the four closest points
  real,    intent(in)    :: lon,lat
  integer, intent(inout) :: IIij(4),JJij(4)
  real,    intent(out)   :: Weight(4)
  integer :: i,j
  real :: dist(0:4)
  dist=1.0E40
  do j=1,ljmax
    do i=1,limax
      !distance between (lon(i,j),lat(i,j) and (lon,lat)
      dist(0)=great_circle_distance(lon,lat,glon(i,j),glat(i,j))
      if(dist(0)<dist(1))then
        dist(4)=dist(3);dist(3)=dist(2);dist(2)=dist(1);dist(1)=dist(0)
        IIij(4)=IIij(3);IIij(3)=IIij(2);IIij(2)=IIij(1);IIij(1)=i
        JJij(4)=JJij(3);JJij(3)=JJij(2);JJij(2)=JJij(1);JJij(1)=j
      elseif(dist(0)<dist(2))then
        dist(4)=dist(3);dist(3)=dist(2);dist(2)=dist(0)
        IIij(4)=IIij(3);IIij(3)=IIij(2);IIij(2)=i
        JJij(4)=JJij(3);JJij(3)=JJij(2);JJij(2)=j
      elseif(dist(0)<dist(3))then
        dist(4)=dist(3);dist(3)=dist(0)
        IIij(4)=IIij(3);IIij(3)=i
        JJij(4)=JJij(3);JJij(3)=j
      elseif(dist(0)<dist(4))then
        dist(4)=dist(0)
        IIij(4)=i
        JJij(4)=j
      endif
    enddo
  enddo
  dist(0)=(dist(1)+dist(2)+dist(3)+dist(4))
  Weight(1)=1.0-3.0*dist(1)/dist(0)
  dist(0)=(dist(2)+dist(3)+dist(4))
  Weight(2)=(1.0-Weight(1))*(1.0-2.0*dist(2)/dist(0))
  dist(0)=(dist(3)+dist(4))
  Weight(3)=(1.0-Weight(1)-Weight(2))*(1.0-dist(3)/dist(0))
  Weight(4)=1.0-Weight(1)-Weight(2)-Weight(3)
end subroutine get_interpol_const

end module DA_Obs_ml
