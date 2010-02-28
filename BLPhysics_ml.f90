module BLPhysics_ml
 ! Collection of boundary layer met routines. Will move more code from tiphys
 !  here in future. Try to keep 1-D or elemental to allow use in offline codes 
 ! (*No* routines in use, except for testing)
 ! Copiedfrom Unimod.rv18WTA  14th Feb 2010
 ! use GridValues_ml,  only : sigma_bnd

 use ModelConstants_ml,  only : KMAX_MID, KMAX_BND, KWINDTOP, PT
 use PhysicalConstants_ml,  only : KARMAN
 implicit none
 private


! Queries  a
! - should we keep this min value. Very high for stable?
! - max value also seems high. We used to use 2500 m

 real, parameter, public :: PZPBL_MIN=200.   ! Amela's
 real, parameter, public :: PZPBL_MAX=5000.  ! Amela's
 integer, parameter, private :: PIELKE_KZ=1
 integer, parameter, public :: KZ_METHOD=PIELKE_KZ
 ! hb 23.02.2010 Kz from meteo 
 logical, parameter, public :: NWP_Kz=.false.
 logical, parameter, public :: OBRIAN_Kz=.true.

! We don't need to calculate Kz for all layer in future maybe
! Still, for safety  we let this extent to K=1 for now

 real, parameter, public :: KPBL_MAX=1      ! Dave, TMP, 

 real, parameter, private :: EPS=0.01  !prevents div by zero for WS

public :: SeibertRiB_Hmix
public :: SeibertRiB_Hmix_3d
public :: JericevicRiB_Hmix
public :: Venkatram_Hmix
public :: Zilitinkevich_Hmix
public :: TI_BLphysics
public :: fake_zmid
public :: fake_zbnd
!public :: Grisogono
public :: risig1
public :: BrostWyngaardKz

public :: SigmaKz_2_m2s ! hb 23.02.2010 Kz from meteo
 private :: SigmaKz_2_m2s_scalar
 private :: SigmaKz_2_m2s_array

! We put GRAV here to simplify dependenceis

 real, parameter, private :: GRAV = 9.807   ! Gravity, m/s2

! conversion of Kz in sigma coordinates to m2/s, from HF Kz(sigma)=Kz*ro**2*(GRAV/p*)**2
! We can call this conversion routine as either scalar or array
  interface SigmaKz_2_m2s
     module procedure SigmaKz_2_m2s_scalar
     module procedure SigmaKz_2_m2s_array
  end interface SigmaKz_2_m2s

contains

 !----------------------------------------------------------------------------
function BrostWyngaardKz(z,h,ustar,invL) result(Kz)
  real, intent(in) :: z    ! height
  real, intent(in) :: h    ! Boundary layer depth 
  real, intent(in) :: ustar!  u*
  real, intent(in) :: invL !  1/L  
  real :: Kz

     if ( z < h ) then
        Kz = KARMAN * ustar * z * (1-z/h)**1.5 / (1+5*z*invL)
     else
        Kz= 0.0
     end if

end function BrostWyngaardKz

! hb 23.02.2010 Kz from meteo
function SigmaKz_2_m2s_scalar (roa,ps) result(Kz_fac)
  real, intent(in) :: roa
  real, intent(in) :: ps
  real :: fac
  real :: Kz_fac

     fac= (ps - PT)/(GRAV*roa)
     Kz_fac= fac*fac

end function SigmaKz_2_m2s_scalar

function SigmaKz_2_m2s_array (roa,ps) result(Kz_fac)
  real, intent(in), dimension(:,:,:) :: roa
  real, intent(in), dimension(:,:)   :: ps
  real, dimension(size(roa,1),size(roa,2),size(roa,3)) :: Kz_fac
  real :: fac
  integer :: i,j,k

  do k = 1, size(roa,3)
  do j = 1, size(ps,2)
  do i = 1, size(ps,1)
     fac= (ps(i,j) - PT)/(GRAV*roa(i,j,k))
     Kz_fac(i,j,k) = fac*fac
  end do
  end do
  end do

end function SigmaKz_2_m2s_array

 !----------------------------------------------------------------------------

! Two Rib-based mixing height methods.
! Seibert et al., AE, 2000, pp1001-,  eqn (9): 
! Jericevic et al., ACP, 2009,  eqn (17): 
! Boundary layer height is calculated with the bulk Richardson number
! method with critical value of Ric=0.25
! Note, can never give first layer as Hmix. Is
! this intenional? What if we have a really stable BL?

!subroutine SH (Nij,Nk,u,v, zm, theta, pzpbl)
!  integer, intent(in) :: Nij, Nk   ! dims (i x j, k)
!  real, dimension(Nij,Nk), intent(in) :: u,v ! winds
!  real, dimension(Nij,Nk), intent(in) :: zm ! mid-cell height
!  real, dimension(Nij,Nk), intent(in) :: theta !pot. temp
!  real, intent(out) :: pzpbl(Nij)
!  integer :: k, n
!  real, parameter :: Ric = 0.25  ! critical Ric
!  real,dimension(Nij) :: Rib  ! bulk Richardson number
!  real,dimension(Nij) :: Theta1   ! pot temp of lowest cell
!
!! Seibert et al., AE, 2000, pp1001-,  eqn (9): 
!! Although we should use virtual pot temp, not just theta
!
!
!   do n = 1, Nij
!     Theta1(n) = theta(n,KMAX_MID)
!     KLOOP: do k=Nk-1, KWINDTOP, -1 ! KMAX_MID -1
!
!        Rib(n) =      GRAV * zm(n,k) &
!             *(theta(n,k)-Theta1(n) ) / &
!             ( Theta1(n) * ( u(n,k)**2 + v(n,k)**2 )+EPS )
!        if (Rib(n) >= Ric) then
!              pzpbl = zm(n,k)
!              exit KLOOP
!        end if 
!      end do KLOOP
!    enddo !n
!
!   !QUERY  pzpbl=min(pzpbl,pzpbl_max )
!   !QUERY  pzpbl=max(pzpbl, pzpbl_min)
!end subroutine SH

subroutine SeibertRiB_Hmix_3d (u,v, zm, theta, pzpbl)
  real, dimension(:,:,:), intent(in) :: u,v ! winds
  real, dimension(:,:,:), intent(in) :: zm ! mid-cell height
  real, dimension(:,:,:), intent(in) :: theta !pot. temp
  real, intent(out) :: pzpbl(:,:)
  integer :: k, n, kmax
  real, parameter :: Ric = 0.25  ! critical Ric
  real,dimension(size(pzpbl,1),size(pzpbl,2)) :: RiB, Theta1   ! pot temp of lowest cell

! Seibert et al., AE, 2000, pp1001-,  eqn (9): 
! Although we should use virtual pot temp, not just theta

     Theta1(:,:) = theta(:,:,KMAX_MID)
     pzpbl = -999.999
     KLOOP: do k=KMAX_MID-1, KWINDTOP, -1 

        where ( pzpbl < 0 ) 
           Rib(:,:) =      GRAV * zm(:,:,k) &
             *(theta(:,:,k)-Theta1(:,:) ) / &
             ( Theta1(:,:) * ( u(:,:,k)**2 + v(:,:,k)**2 )+EPS )
           where (Rib(:,:) >= Ric) 
              pzpbl(:,:) = zm(:,:,k)
           end where 
         end where 

      end do KLOOP

   !QUERY  pzpbl=min(pzpbl,pzpbl_max )
   !QUERY  pzpbl=max(pzpbl, pzpbl_min)
end subroutine SeibertRiB_Hmix_3d

 !----------------------------------------------------------------------------
subroutine SeibertRiB_Hmix (u,v, zm, theta, pzpbl)
  real, dimension(KWINDTOP:KMAX_MID), intent(in) :: u,v ! winds
  real, dimension(KMAX_MID), intent(in) :: zm ! mid-cell height
  real, dimension(KMAX_MID), intent(in) :: theta !pot. temp
  real, intent(out) :: pzpbl
  integer :: k
  real, parameter :: Ric = 0.25  ! critical Ric
  real :: Rib  ! bulk Richardson number
  real :: Theta1   ! pot temp of lowest cell

! Seibert et al., AE, 2000, pp1001-,  eqn (9): 
! Although we should use virtual pot temp, not just theta

   Theta1 = theta(KMAX_MID)
   do k=KMAX_MID-1, KWINDTOP, -1

       Rib =      GRAV * zm(k) &
             *(theta(k)-Theta1 ) / &
             ( Theta1 * ( u(k)**2 + v(k)**2 )+EPS )
             !print *, k, zm(k), theta(k), sqrt(( u(k)**2 + v(k)**2 )),  RiB
       if(Rib >= Ric) then
              pzpbl = zm(k)
              exit
       endif
    enddo

   !QUERY  pzpbl=min(pzpbl,pzpbl_max )
   !QUERY  pzpbl=max(pzpbl, pzpbl_min)
end subroutine SeibertRiB_Hmix

 !----------------------------------------------------------------------------

subroutine JericevicRiB_Hmix (u,v, zm, theta, zi)
  real, dimension(KWINDTOP:KMAX_MID), intent(in) :: u,v ! winds
  real, dimension(KMAX_MID), intent(in) :: zm ! mid-cell height
  real, dimension(KMAX_MID), intent(in) :: theta !pot. temp
  real, intent(out) :: zi
  integer :: k
  real, parameter :: Ric = 0.25  ! critical Ric
  real :: Rib  ! bulk Richardson number
  real :: Theta1, z1   ! pot temp  and height of lowest cell

! Jericevic et al., ACP, 2009, pp1001-,  eqn (17): 

   Theta1 = theta(KMAX_MID)
   z1     = zm(KMAX_MID)
   zi     = z1  ! start val

   do k=KMAX_MID-1, KWINDTOP, -1

       Rib =   GRAV * ( zm(k) - z1 ) &
             * (theta(k)-Theta1 ) / &
       ( 0.5*(theta(k)+Theta1) * ( u(k)**2 + v(k)**2 )+EPS )
       if(Rib >= Ric) then
              zi = zm(k)
              exit
       endif
    enddo

   !QUERY  pzpbl=min(pzpbl,pzpbl_max )
   !QUERY  pzpbl=max(pzpbl, pzpbl_min)
end subroutine JericevicRiB_Hmix

 !----------------------------------------------------------------------------
function Venkatram_Hmix (ustar) result(zi)
  real, intent(in) :: ustar
  real :: zi
  ! From Venkatram, 1980, super-simple method!

     zi = 2.4e3 * ustar**1.5

end function Venkatram_Hmix
 !----------------------------------------------------------------------------

function Zilitinkevich_Hmix (ustar,invL,lat) result(zi)
 ! Not equator-proof yet!!
   real, intent(in) :: ustar, invL, lat
   real :: f, pi = 4.0*atan(1.0)
   real :: zi
   f = 1.46e-4 * sin(lat*pi/180.0)
   !zi = 0.2 *  ustar/f
   if( invL > f/ustar ) then
     zi = 0.4 * sqrt( ustar/(f*invL) )
   else
     zi = 0.2 *  ustar/f
   end if
end function Zilitinkevich_Hmix

 !----------------------------------------------------------------------------

subroutine TI_BLphysics (u,v, zm, zb, fh, th, exnm, pm, zi, debug_flag)
  real, dimension(:), intent(in) :: u,v ! wind-vels
  real, dimension(:), intent(in) :: zm, zb ! heights of mid and bnd
  real, intent(in) :: fh  ! surface flux of sensible heat, W/m2
  real, dimension(:), intent(in) :: th
  real, dimension(:), intent(in) :: exnm, pm ! exner_mid and p_mid
  real, intent(out) :: zi  ! Final height of ABL
  logical, intent(in) :: debug_flag

  real, dimension(KMAX_BND) :: &
       Ris   &! Richardson number, sigma coords
      ,xksig &! xksig smoothed over three adjacent layers
      ,xksm   ! spacially smoothed Kz in z direction, m2/s.

  integer :: k, km, km1, km2, kp, nh1, nh2
  real :: xl2    ! mixing length squared
  real, save :: zmmin = 200.0 !QUERY, Pielke or EMEP?
  real, save :: zimin = 100.0
  real, save :: zlimax= 3000.0
  real :: Ric, Ric0 !critical Richardson number variables
  real :: dvdz
!  real, parameter :: Ric = 0.25  ! critical Ric
  real, parameter :: DTZ = 3600.0 !time interval for integration of surface 
                                  !heat fluxes in ABL-height calculations, s
  real, parameter :: KZ_MINIMUM = 0.001   ! m2/s
  real, parameter :: KZ_MAXIMUM = 1.0e3 ! m2/s - as old kzmax
  real, parameter :: KZ_SBL_LIMIT = 0.1 ! m2/s - Defines stable BL height
  real :: zis    ! height of the stable ABL, m

  ! For unstable BL
  real :: delq   ! available heat flux for developing the unstable ABL, J/m2
                 ! (heat-input per m2 from the ground during unstable BL)
  real :: thsrf, xdth, dthdzm, dthc, xdthdz
  real, dimension(KMAX_MID) :: thc, dthdz ! Assures th increases with height
  real :: dpidth !heat increasement in accordance with temp. increasement, J/m2
  real :: pidth  !heat used to adjust air temperature, J/m2
  integer :: trc !help variable telling whether or not unstable ABL exists
  integer :: kabl
  real :: ziu    ! height of the unstable ABL, m

 !..the following variables in sigmas-levels:
 !
  do  k=2,KMAX_MID

      km=k-1

     !..wind sheer (first keep as squared)

      dvdz = (u(km)-u(k))**2 + (v(km)-v(k))**2 + EPS 

      Ris(k)=(2*GRAV/(th(km)+th(k))) * (th(km)-th(k))*(zm(km)-zm(k)) &
                  /dvdz

      !........................
      !..mixing length squared:
      !
       xl2=(KARMAN*min(zb(k),zmmin))**2

      !
      !..............................
      !..critical richardsons number:
      !
       Ric0=0.115*((zm(km)-zm(k))*100.0)**0.175
       Ric=max(0.25,Ric0)

       dvdz = sqrt(dvdz)/(zm(km)-zm(k))

      !..................................................................
      !..exchange coefficient (Pielke,...)
      if ( KZ_METHOD == PIELKE_KZ  ) then
         if (Ris(k) > Ric ) then
            xksig(k) = KZ_MINIMUM
         else
            xksig(k) = 1.1 * (Ric-Ris(k)) * xl2 * dvdz /Ric
         end if
      else

         !..exchange coefficient (blackadar, 1979; iversen & nordeng, 1987):
         !
         if(Ris(k) <= 0.0) then
            xksig(k)=xl2*dvdz*sqrt(1.1-87.*Ris(k))
         elseif(Ris(k) <= 0.5*Ric) then
                xksig(k)=xl2*dvdz*(1.1-1.2*Ris(k)/Ric)
         elseif(Ris(k) <= Ric) then
                xksig(k)=xl2*dvdz*(1.-Ris(k)/Ric)
         else
                xksig(k)=0.001
         endif
      end if ! Pielke or Blackadar

  end do ! k

  k=2
  km=1
  kp=3
  xksm(k)=( (zm(km)-zm(k))*xksig(k) + (zm(k)-zm(kp))*xksig(kp) )&
              / ( zm(km) - zm(kp) )

  k=KMAX_MID
  km2=k-2
  km1=k-1
  xksm(k)=( (zm(km2)-zm(km1))*xksig(km1) + (zm(km1)-zm(k))*xksig(k) )&
               / ( zm(km2) - zm(k) )

  do k = 3,KMAX_MID-1
      km1=k-1
      km2=k-2
      kp=k+1
      xksm(k)=(  (zm(km2)-zm(km1))*xksig(km1) + (zm(km1)-zm(k))*xksig(k)&
            + (zm(k)-zm(kp))*xksig(kp) ) / ( zm(km2) - zm(kp) )
  end do ! k

 !............................................................
 !..The height of the stable BL is the lowest level for which:
 !..xksm .le. 1 m2/s (this limit may be changed):
 
  zis=zimin
  nh1 = KMAX_MID
  nh2 = 1

  do k=KMAX_MID,2,-1

     if(xksm(k) >= KZ_SBL_LIMIT .and. nh2 == 1) then
        nh1=k   ! Still unstable
     else
        nh2=0   ! Now stable
     endif
  end do

  k=nh1
  if(zb(nh1) >= zimin) then

      if( abs(xksm(k)-xksm(k-1)) > eps) then

           zis=((xksm(k)-KZ_SBL_LIMIT )*zb(k-1) &
               + (KZ_SBL_LIMIT -xksm(k-1))*zb(k))&
                     /(xksm(k)-xksm(k-1))
      else
          zis=zimin
      endif

   endif


   zi = zis

 !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !---------------------------------------------------------------------
  !....................................
  !..height of unstable boundary layer:
  !
  !
  !..assuring that th is increasing with height.
  !..adjusted th-sounding is assigned to thc-array.
  !..This adjusted th is not meant to be used in
  !..other parts of the model program
  !
  dthdzm = 1.e-4
  thc(KMAX_MID)=th(KMAX_MID)
  do k=KMAX_MID-1,1,-1

      dthc = (th(k)-th(k+1))&
             / (zm(k)-zm(k+1))

      dthdz(k)=max(dthc,dthdzm)

      thc(k)=thc(k+1)+dthdz(k)*(zm(k)-zm(k+1))
if ( debug_flag ) write(6,"(a,i3,3es10.3)") "DEBUG THC ", k, th(k), dthc, thc(k)

  enddo

  !..estimated as the height to which an hour's input
  !..of heat from the ground is vertically distributed,
  !..assuming dry adiabatic adjustment.
  !
  !

  delq=-min(fh,0.0)*DTZ
     thsrf=0.0
     ziu=0.0
if ( debug_flag ) write(6,"(a,3es10.3)") "DEBUG fH, DELQ ", fh, delq

     !.................................

     trc=0                        !..   =0 for stable BL (delq=0):
     if(delq >  0.00001) trc=1    !..trc=1 for unstable BL (delq>0):

     !------------------------------------------------------------
     ! calculating the height of unstable ABL
     !
    kabl = KMAX_MID
    do while( trc == 1)
        kabl = kabl-1
        pidth=0.

        do k=KMAX_MID,kabl,-1
           xdth = thc(kabl)-thc(k)
           dpidth = exnm(k)*xdth*(pm(k+1)-pm(k))/GRAV
           pidth = pidth + dpidth
if ( debug_flag ) write(6,"(a,2i3,6es11.3,i4)") "DEBUG PID ", kabl, k, xdth, exnm(k), pm(k),  dpidth, pidth 
        end do

       if(pidth >= delq.and.trc == 1  ) then

          !  at level kabl or below level kabl and above level kabl+1

          thsrf = thc(kabl)- (thc(kabl)-thc(KMAX_MID))    &
                * (pidth-delq)/pidth

          xdthdz = (thc(kabl)-thc(kabl+1)) / (zm(kabl)-zm(kabl+1))

          ziu = zm(kabl+1) + (thsrf-thc(kabl+1))/xdthdz

          trc=0  
if ( debug_flag ) write(6,"(a,i3,2es10.3,i4)") "DEBUG PICTH ", kabl, delq, pidth, trc 

        endif


if ( debug_flag ) write(6,"(a,i3,es10.3,i5)") "DEBUG mid ", kabl, delq, trc 
        if(kabl <= 4 .and. trc == 1  ) then

           write(6,*)'PBL ziu calculations failed!'

           ziu=zlimax

           trc=0 
        endif

     end do ! while

  zi = max( ziu, zis)
  zi = min( zlimax, zi)


end subroutine TI_BLphysics

!     !AJ-have tested and this can be removed     ustar_nwp(i,j)=max(ustar_nwp(i,j), 0.1)
!
!     !  ustar_nwp(i,j) = max( ustar_nwp(i,j), 1.0e-5 )
!     ! and a consistency between tau and ustar. Here you reset ustar
!     ! and this may affect other routines?
!
!!     Mixing height calculation
!!     Boundary layer height is calculated with the bulk Richardson number
!!      method with critical value of Ric=0.25
!
!            do k=KMAX_MID,2,-1
!               ri(i,j,k) = risig1(u_xmj(i,j,k-1,nr), v_xmi(i,j,k-1,nr), &
!                              th(i,j,KMAX_MID,nr), th(i,j,k-1,nr), &
!                              zm(i,j,k-1))
!               if(ri(i,j,k) >= ric) then
!                  pzpbl(i,j) = zm(i,j,k)
!                  exit
!               endif
!            enddo
!
!            pzpbl(i,j)=min(pzpbl(i,j),pzpbl_max )
!            pzpbl(i,j)=max(pzpbl(i,j), pzpbl_min)
!!AJ You have this minimum pbl height limit in earlier version.
!!AJ I kept it becouse it might influence dry deposition.
!
!! Determination of two input parameters needed in Grisogono aproach
!
!            !DS Kmax(i,j)=const(i,j)*vsm(i,j)*ustar_nwp(i,j)
!            !AJ Yes, this is better
!            Kmax(i,j)=KZ_CONST*pzpbl(i,j)*ustar_nwp(i,j)
!            h(i,j)=max(pzpbl(i,j)/3.,zm(i,j,KMAX_MID-1))
!
! !AJ changed 3 to 4 to , result of tests in 2001 on NO2, SO2 and SO4
! ! Grisogono approach is valid for boundary layer, so z(i,j) is defined as
! ! upper boundary 
! ! above the PBL K(z) calculated with Blackadar approach is used 
!
!!DS - query -what did you mean by nz? z perhaps?
!!AJ -yes z not nz, typing mistake
!
!!DS - faster to write (?):
!!AJ -OK
!               z(i,j)= 1.25 * pzpbl(i,j)
!
!!Calculation of the exponential profile
!
!            do k=KMAX_MID,2,-1
!!AJ               jedan=Kmax(i,j)*zm(i,j,k)*sqrt(exp(1.))/h(i,j)
!!AJ               dva=exp(-0.5*(zm(i,j,k)/h(i,j))**2)
!                  xksig(i,j,k)=Kmax(i,j)*zm(i,j,k)*sqrt(exp(1.))/h(i,j)&
!                              *exp(-0.5*(zm(i,j,k)/h(i,j))**2)+KZ_MINIMUM
!               if(zm(i,j,k)>z(i,j)) exit
!            enddo
!         enddo
!      enddo
!
!      !
!      ! storing of pbl height
!      ! Why use vsm at all, could use pzpbl(i,j) from start?
!      !
!      !AJ changed vsm to pzpbl
!      
!      !smoothing of pbl height
!
!        call smoosp(pzpbl,pzpbl_min,pzpbl_max)
!      !
!      ! transformation of xksig(i,j,k) into skh(i,j,k,nr)
!      ! from HF Kz(sigma)=Kz*ro**2*(GRAV/p*)
!      !
!      do k=2,KMAX_MID
!     !AJ added smoothing of Kz with Shapiro -same as with O_Brien 
!      do i=1,limax
!           do j=1,ljmax
!             help(i,j)=xksig(i,j,k)
!           enddo
!      enddo
!
!      call smoosp(help,KZ_MINIMUM,KZ_MAXIMUM)
!         do i=1,limax
!            do j=1,ljmax
!
!!I have exns=exf1,exnm=exf2,
!               fac = GRAV/(ps(i,j,nr) - PT)
!               fac2 = fac*fac
!               dex12 = th(i,j,k-1,nr)*(exnm(i,j,k)-exns(i,j,k)) &
!                     + th(i,j,k,nr)*(exns(i,j,k)-exnm(i,j,k-1))
!
!               ro = ((ps(i,j,nr)-PT)*sigma_bnd(k)+PT)*CP &
!                  *(exnm(i,j,k)-exnm(i,j,k-1))/(RGAS_KG*exns(i,j,k)*dex12)
!
!               skh(i,j,k,nr) = xksig(i,j,k)*ro*ro*fac2
!            enddo
!         enddo
!      enddo
!   end subroutine Grisogono
 !----------------------------------------------------------------------------
  function risig1(ws,th1,th2,z)    ! Amela, not used in new version
   !calculates the bulk Richardson number
      implicit none
      real, intent(in) :: ws   ! wind-speed
      real, intent(in) :: th1, th2
      real, intent(in) :: z
      real :: risig1
      real :: dvdz

      dvdz = ws*ws+0.001
      risig1=(2.0*GRAV/(th1+th2))*(th2-th1)*z/dvdz
   end function risig1
 !----------------------------------------------------------------------------

subroutine fake_zmid(z)
   
 ! Data for testing. Note that z here is really geopotetntial
 ! heught, but near the ground it doesn't matter so much
 ! from 1 time-period, see DEBUG_Z in Met_ml
   real, dimension(20), intent(out) :: z
   z( 1) =   15244.7131
   z( 2) =   13531.2981
   z( 3) =   12177.9740
   z( 4) =   11004.0286
   z( 5) =    9788.6193
   z( 6) =    8430.8546
   z( 7) =    7013.5860
   z( 8) =    5700.0226
   z( 9) =    4598.2893
   z(10) =    3692.3391
   z(11) =    2934.3594
   z(12) =    2296.2730
   z(13) =    1761.6936
   z(14) =    1317.4565
   z(15) =     951.4678
   z(16) =     655.8842
   z(17) =     424.7443
   z(18) =     254.0498
   z(19) =     137.4284
   z(20) =      45.6146
end subroutine fake_zmid

subroutine fake_zbnd(z)
   
 ! Data for testing. Note that z here is really geopotetntial
 ! heught, but near the ground it doesn't matter so much
 ! from 1 time-period, see DEBUG_Z in Met_ml
   real, dimension(21), intent(out) :: z
  z( 1) =    16276.236
  z( 2) =    14319.868
  z( 3) =    12815.771
  z( 4) =    11598.220
  z( 5) =    10461.131
  z( 6) =     9168.470
  z( 7) =     7747.153
  z( 8) =     6324.928
  z( 9) =     5103.225
  z(10) =     4110.656
  z(11) =     3286.077
  z(12) =     2591.897
  z(13) =     2007.792
  z(14) =     1520.443
  z(15) =     1116.948
  z(16) =      787.332
  z(17) =      525.163
  z(18) =      324.823
  z(19) =      183.602
  z(20) =       91.302
  z(21) =       91.302
end subroutine fake_zbnd


end module BLPhysics_ml

!program testBLM
!  use ModelConstants_ml, only : KMAX_MID
!  use BLPhysics_ml
!  real, dimension(KMAX_MID) :: u,v,zm,theta
!  real :: pzpbl
!  call SeibertRiB_Hmix (u,v, zm, theta, pzpbl)
!end program testBLM
